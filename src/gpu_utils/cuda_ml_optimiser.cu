#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include "src/gpu_utils/cuda_ml_optimiser.h"
#include "src/gpu_utils/cuda_img_operations.h"
#include "src/complex.h"
#include <fstream>
#include <cuda.h>

#define BLOCK_SIZE 128         // This is optimally set as big as possible without its ceil:ed multiple exceeding imagesize by too much.
#define NR_CLASS_MUTEXES 5

static pthread_mutex_t global_mutex2[NR_CLASS_MUTEXES] = { PTHREAD_MUTEX_INITIALIZER };
static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "CUDA ERROR: %s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

class CudaComplex
{
public:
	double real, imag;

	inline
	__device__ __host__ CudaComplex(): real(), imag() {};
	inline
	__device__ __host__ CudaComplex(double real, double imag): real(real), imag(imag) {};
};

class CudaImages
{
public:
	unsigned x,y,xy,num,max_num;
	CudaComplex* start;

	inline
	CudaImages(long unsigned x, long unsigned y, long unsigned max_num):
			x(x), y(y), num(0), max_num(max_num), xy(x*y), start(new CudaComplex[xy*max_num])
	{};

	inline
	CudaComplex* current() { return start + (num*xy); };

	inline
	void increment() { num++; };

	inline
	CudaComplex* operator [](long unsigned i) { return start + (i*((long unsigned) xy)); };

	inline
	long unsigned alloc_size() { return ((long unsigned) num)*((long unsigned) xy); };

	inline
	CudaComplex* data_to_device()
	{
		CudaComplex* d_ptr(0);
		HANDLE_ERROR(cudaMalloc( (void**) &d_ptr, alloc_size() * sizeof(CudaComplex)));
		HANDLE_ERROR(cudaMemcpy( d_ptr, start, alloc_size() * sizeof(CudaComplex), cudaMemcpyHostToDevice));
		return d_ptr;
	}

	inline
	void clear() { delete[] start; }

	inline
	~CudaImages() { delete[] start; }
};

__global__ void cuda_kernel_massive_diff2(	CudaComplex *g_refs, CudaComplex *g_imgs,
									double *g_Minvsigma2, double *g_diff2s,
									const unsigned img_size, const double sum_init,
									bool *g_exp_Mcoarse_significant,
									long int orientation_num,
									long int translation_num,
									long int exp_nr_oversampled_rot,
									long int exp_nr_oversampled_trans)
{
//	int ex = blockIdx.x % orientation_num;
//	int ey = (blockIdx.x - ex) / orientation_num;
	int ex = blockIdx.y * gridDim.x + blockIdx.x;
	//int ey = blockIdx.y;
	int ez = blockIdx.z;

	unsigned long int coarse_rot_idx   = floorf(ex/exp_nr_oversampled_rot);
	unsigned long int coarse_trans_idx = floorf(ez/exp_nr_oversampled_trans);

	// 		Check if it is significant
	//          		AND
	// inside the padded 2D orientation grid
	if(g_exp_Mcoarse_significant + ex + ez*coarse_rot_idx && ex < orientation_num )
	{
		__shared__ double s[BLOCK_SIZE];
		s[threadIdx.x] = 0;

		unsigned pass_num(ceilf((float)img_size/(float)BLOCK_SIZE));
		unsigned long pixel,
		ref_start(ex * img_size),
		img_start(ez * img_size);

		unsigned long ref_pixel_idx;
		unsigned long img_pixel_idx;

		for (unsigned pass = 0; pass < pass_num; pass ++)
		{
			pixel = pass * BLOCK_SIZE + threadIdx.x;

			if (pixel < img_size) //Is inside image
			{
				ref_pixel_idx = ref_start + pixel;
				img_pixel_idx = img_start + pixel;

				double diff_real = g_refs[ref_pixel_idx].real - g_imgs[img_pixel_idx].real;
				double diff_imag = g_refs[ref_pixel_idx].imag - g_imgs[img_pixel_idx].imag;

				s[threadIdx.x] += (diff_real * diff_real + diff_imag * diff_imag) * 0.5 * g_Minvsigma2[pixel];
			}
		}

		// This version should run in             BLOCK_SIZE                  cycles
		// -------------------------------------------------------------------------
	//		if (threadIdx.x == 0)
	//		{
	//			double sum(sum_init);
	//			for (unsigned i = 0; i < BLOCK_SIZE; i++)
	//				sum += s[i];
	//
	//			g_diff2s[ex * translation_num + ey] = sum;
	//		}
		// -------------------------------------------------------------------------

		// This version should run in     BLOCK_SIZE/trads + log2(trads)      cycles
		// ( Runs ~2x as fast as the above one for BLOCK_SIZE=32 )
		// -------------------------------------------------------------------------
		__syncthreads();
		int trads = 32;
		int itr = BLOCK_SIZE/trads;
		if(threadIdx.x<trads)
		{
			for(int i=1; i<itr; i++)
			{
				s[threadIdx.x] += s[i*trads + threadIdx.x];
				//__syncthreads();
			}
		}

		for(int j=(trads/2); j>0; j/=2)
		{
			if(threadIdx.x<j)
			{
				s[threadIdx.x] += s[threadIdx.x+j];
			}
		}
		__syncthreads();
//		if (threadIdx.x*ex == 0)
		{
			g_diff2s[ex * translation_num + ez] = s[0]+sum_init;
		}
		// -------------------------------------------------------------------------
	}
//	else
//	{
//		g_diff2s[ex * translation_num + ey] = 0; //(float)g_exp_Mcoarse_significant[blockIdx.x+blockIdx.y*coarse_rot_idx];
//	}


}

//  Takes a boolean N-by-M matrix and returns pointer pairs to coordinates in two corresponding objects
//__global__ void cuda_kernel_boolToPointers(	bool *matrix,
//												int yLength,
//												int** yPoints)
//{
//	//save the current index of the partial array to a shared location
//	__shared__  long int  length[blockDim.x*BLOCK_SIZE];
//	length[threadIdx.x]=0;
//
//	unsigned yiter(ceilf((float)yLength/(float)BLOCK_SIZE));
//
//	for(i=0; i<yiter; i++)
//	{
//		int pos = ylength*blockIdx.x + i*BLOCK_SIZE + threadIdx.x
//		if(matrix[pos]==1)
//		{
//			yPoints[blockIdx.x][length[blockIdx.x*BLOCK_SIZE+threadidx.x]]=blockIdx.x;
//			length[blockIdx.x*BLOCK_SIZE+threadidx.x]+=1;
//		}
//	}
//
//}

void MlOptimiserCUDA::getAllSquaredDifferences(
		long int my_ori_particle, int exp_current_image_size,
		int exp_ipass, int exp_current_oversampling, int metadata_offset,
		int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
		int exp_itrans_min, int exp_itrans_max, int exp_iclass_min, int exp_iclass_max,
		std::vector<double> &exp_min_diff2,
		std::vector<double> &exp_highres_Xi2_imgs,
		std::vector<MultidimArray<Complex > > &exp_Fimgs,
		std::vector<MultidimArray<double> > &exp_Fctfs,
		MultidimArray<double> &exp_Mweight,
		MultidimArray<bool> &exp_Mcoarse_significant,
		std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
		std::vector<double> &exp_directions_prior, std::vector<double> &exp_psi_prior,
		std::vector<MultidimArray<Complex> > &exp_local_Fimgs_shifted,
		std::vector<MultidimArray<double> > &exp_local_Minvsigma2s,
		std::vector<MultidimArray<double> > &exp_local_Fctfs,
		std::vector<double> &exp_local_sqrtXi2)
{
	// Initialise min_diff and exp_Mweight for this pass
	int exp_nr_particles = mydata.ori_particles[my_ori_particle].particles_id.size();
	long int exp_nr_dir = (do_skip_align || do_skip_rotate) ? 1 : sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
	long int exp_nr_psi = (do_skip_align || do_skip_rotate) ? 1 : sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior);
	long int exp_nr_trans = (do_skip_align) ? 1 : sampling.NrTranslationalSamplings();
	long int exp_nr_oversampled_rot = sampling.oversamplingFactorOrientations(exp_current_oversampling);
	long int exp_nr_oversampled_trans = sampling.oversamplingFactorTranslations(exp_current_oversampling);

	//for scale_correction
	int group_id;
	double myscale;

	//printf("exp_nr_oversampled_rot=%d\n", (unsigned)exp_nr_oversampled_rot);

	exp_Mweight.resize(exp_nr_particles, mymodel.nr_classes * exp_nr_dir * exp_nr_psi * exp_nr_trans * exp_nr_oversampled_rot * exp_nr_oversampled_trans);
	exp_Mweight.initConstant(-999.);
	if (exp_ipass==0)
		exp_Mcoarse_significant.clear();

	exp_min_diff2.clear();
	exp_min_diff2.resize(exp_nr_particles, 99.e99);

	std::vector<MultidimArray<Complex > > dummy;
	precalculateShiftedImagesCtfsAndInvSigma2s(false, my_ori_particle, exp_current_image_size, exp_current_oversampling,
			exp_itrans_min, exp_itrans_max, exp_Fimgs, dummy, exp_Fctfs, exp_local_Fimgs_shifted, dummy,
			exp_local_Fctfs, exp_local_sqrtXi2, exp_local_Minvsigma2s);

	// Loop only from exp_iclass_min to exp_iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++)
	{
		if (mymodel.pdf_class[exp_iclass] > 0.)
		{
			// Local variables
			std::vector< double > oversampled_rot, oversampled_tilt, oversampled_psi;
			std::vector< double > oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
			MultidimArray<Complex > Fref;
			double *Minvsigma2;
			Matrix2D<double> A;

			CudaImages Frefs(exp_local_Minvsigma2s[0].xdim, exp_local_Minvsigma2s[0].ydim,
					(exp_idir_max - exp_idir_min + 1) * (exp_ipsi_max - exp_ipsi_min + 1) * exp_nr_oversampled_rot);

			// Mapping index look-up table
			std::vector< long unsigned > iorientclasses, iover_rots;
			long unsigned orientation_num(0);

			/*=======================================================================================
			                           Generate Reference Projections
			=========================================================================================*/

			//printf("Generate Reference Projections\n");

			Fref.resize(exp_local_Minvsigma2s[0]); //TODO remove this
			Complex* FrefBag = Fref.data; //TODO remove this

			for (long int idir = exp_idir_min, iorient = 0; idir <= exp_idir_max; idir++)
			{
				for (long int ipsi = exp_ipsi_min; ipsi <= exp_ipsi_max; ipsi++, iorient++)
				{
					long int iorientclass = exp_iclass * exp_nr_dir * exp_nr_psi + iorient;

					// Get prior for this direction and skip calculation if prior==0
					double pdf_orientation;
					if (do_skip_align || do_skip_rotate)
					{
						pdf_orientation = mymodel.pdf_class[exp_iclass];
					}
					else if (mymodel.orientational_prior_mode == NOPRIOR)
					{
						pdf_orientation = DIRECT_MULTIDIM_ELEM(mymodel.pdf_direction[exp_iclass], idir);
					}
					else
					{
						pdf_orientation = exp_directions_prior[idir] * exp_psi_prior[ipsi];
					}
					// In the first pass, always proceed
					// In the second pass, check whether one of the translations for this orientation of any of the particles had a significant weight in the first pass
					// if so, proceed with projecting the reference in that direction
					bool do_proceed = (exp_ipass==0) ? true :
						isSignificantAnyParticleAnyTranslation(iorientclass, exp_itrans_min, exp_itrans_max, exp_Mcoarse_significant);
					if (do_proceed && pdf_orientation > 0.)
					{
						// Now get the oversampled (rot, tilt, psi) triplets
						// This will be only the original (rot,tilt,psi) triplet in the first pass (exp_current_oversampling==0)
						sampling.getOrientations(idir, ipsi, exp_current_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
								exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior);
						// Loop over all oversampled orientations (only a single one in the first pass)
						for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
						{
							// Get the Euler matrix
							Euler_angles2matrix(oversampled_rot[iover_rot],
												oversampled_tilt[iover_rot],
												oversampled_psi[iover_rot], A);

							Fref.data = (Complex*) Frefs.current();

							// Project the reference map (into Fref)
							(mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_NOT_INV);

							//TODO REMOVE ONCE YOU KNOW THIS IS ALLWAYS TRUE
							if (Frefs.x != Fref.xdim || Frefs.y != Fref.ydim)
								std::cerr << "!!!!!!! BAD Fref size x:" << Fref.xdim << ":" << Frefs.x << " y:" << Fref.ydim << ":" << Frefs.y << std::endl;

							Frefs.increment();

							orientation_num ++;
							iorientclasses.push_back(iorientclass);
							iover_rots.push_back(iover_rot);
						}
					}
				}
			}
			//printf("Finished generating reference projections\n");

			Fref.data = FrefBag; //TODO remove this

			CudaComplex *d_Frefs = Frefs.data_to_device();

			/*=======================================================================================
			                                  	  Particle Iteration
			=========================================================================================*/

			for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
			{
				/*====================================
				        Generate Translations
				======================================*/

				CudaImages Fimgs(Frefs.x, Frefs.y,
						( exp_itrans_max - exp_itrans_min + 1) * exp_nr_oversampled_trans);

				long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
				long unsigned translation_num(0), ihidden(0);
				std::vector< long unsigned > iover_transes, itranses, ihiddens;

				//printf("Generating translations \n");

				for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++, ihidden++)
				{
					sampling.getTranslations(itrans, exp_current_oversampling,
							oversampled_translations_x, oversampled_translations_y, oversampled_translations_z );

					for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
					{
						/// Now get the shifted image
						// Use a pointer to avoid copying the entire array again in this highly expensive loop
						Complex *myAB;
						if (exp_current_oversampling == 0)
						{
							myAB = (Frefs.y == coarse_size) ? global_fftshifts_ab_coarse[itrans].data
									: global_fftshifts_ab_current[itrans].data;
						}
						else
						{
							int iitrans = itrans * exp_nr_oversampled_trans +  iover_trans;
							myAB = (strict_highres_exp > 0.) ? global_fftshifts_ab2_coarse[iitrans].data
									: global_fftshifts_ab2_current[iitrans].data;
						}
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Fimgs_shifted[ipart])
						{
							double real = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real
									- (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag;
							double imag = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag
									+ (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real;

							//When on gpu, it makes more sense to ctf-correct translated images, rather than anti-ctf-correct ref-projections
							if (do_scale_correction)
							{
								//group_id = mydata.getGroupId(part_id);
								float myscale = mymodel.scale_correction[group_id];
								real /= myscale;
								imag /= myscale;
							}
							if (do_ctf_correction && refs_are_ctf_corrected)
							{
								real /= DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], n);
								imag /= DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], n);
							}
							*(Fimgs.current() + n) = CudaComplex(real, imag);
						}
						Fimgs.increment();
						translation_num ++;

						ihiddens.push_back(ihidden);
						itranses.push_back(itrans);
						iover_transes.push_back(iover_trans);
					}
				}
				//printf("Generating translations finished \n");

				/*===========================================
				   Determine significant comparison indices
				=============================================*/
				//      This section is annoying to test because
				//		it can't complete on first pass, since
				//		the significance has never been set

				long int coarse_num = exp_nr_dir*exp_nr_psi*exp_nr_trans;
				long int significant_num=0;
				std::cerr << "exp_ipass "<< exp_ipass << std::endl;
				if (exp_ipass == 0)
				{
					exp_Mcoarse_significant.resize(coarse_num, 1);
					for (long int i = 0; i < coarse_num; i++)
					{
						DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, i)=1;
//						std::cerr << "exp_Mcoarse_significant("<< i <<") = " <<    DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, i) << std::endl;
//						std::cerr << "exp_Mcoarse_significant("<< i <<") = " << *(&DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, 0)+i*sizeof(bool)) << std::endl;
					}
					significant_num = coarse_num;
				}
				else
				{
					std::vector< long unsigned > transidx, rotidx;
					for (long int i = 0; i < orientation_num; i++)
					{
						long int iover_rot = iover_rots[i];
//						long int iover_rot = i % exp_nr_oversampled_rot
						long int coarse_rot = floor(i/exp_nr_oversampled_rot);
						for (long int j = 0; j < translation_num; j++)
						{
							long int iover_trans = iover_transes[j];
//							long int iover_trans = j % exp_nr_oversampled_trans
							long int coarse_trans = floor(j/exp_nr_oversampled_trans);
							long int ihidden = iorientclasses[i] * exp_nr_trans + ihiddens[j];
							if(DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden)==1)
							{
								 long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden,
										                  exp_current_oversampling, iover_rot, iover_trans);
								 transidx.push_back(i);
								 rotidx.push_back(j);
								 significant_num++;
							}
						}
					}
				}
				std::cerr << "orientation_num "<< orientation_num << std::endl;
				std::cerr << "translation_num "<< translation_num << std::endl;
				std::cerr << "my_nr_significant_coarse_samples "<< DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_NR_SIGN) << std::endl;
				std::cerr << "significant_num "<< significant_num << std::endl;

				/*====================================
				   Initiate Particle Related On GPU
				======================================*/

				//When on gpu, it makes more sense to ctf-correct translated images, rather than anti-ctf-correct ref-projections
				if (do_ctf_correction && refs_are_ctf_corrected)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Fimgs_shifted[ipart])
					{
						DIRECT_MULTIDIM_ELEM(exp_local_Minvsigma2s[ipart], n) *= (DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], n)*DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], n));
					}
				}
				// TODO :    + Assure accuracy with the implemented GPU-based ctf-scaling
				//           + Make setting of myscale robust between here and above.
				//  (scale_correction turns off by default with only one group: ml_optimiser-line 1067,
				//   meaning small-scale test will probably not catch this malfunctioning when/if it breaks.)
				if (do_scale_correction)
				{
					float myscale = mymodel.scale_correction[group_id];
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Fimgs_shifted[ipart])
					{
						DIRECT_MULTIDIM_ELEM(exp_local_Minvsigma2s[ipart], n) *= (myscale*myscale);
					}
				}

				Minvsigma2 = exp_local_Minvsigma2s[ipart].data;
				double *d_Minvsigma2(0);

				CudaComplex *d_Fimgs = Fimgs.data_to_device();

				HANDLE_ERROR(cudaMalloc( (void**) &d_Minvsigma2, Fimgs.xy * sizeof(double)));
				HANDLE_ERROR(cudaMemcpy( d_Minvsigma2, exp_local_Minvsigma2s[ipart].data, Fimgs.xy * sizeof(double), cudaMemcpyHostToDevice));

				double *d_diff2s(0);
				HANDLE_ERROR(cudaMalloc( (void**) &d_diff2s, orientation_num*translation_num * sizeof(double)));
				//HANDLE_ERROR(cudaMemset(d_diff2s, exp_highres_Xi2_imgs[ipart] / 2., orientation_num*translation_num * sizeof(double))); //Initiate diff2 values with zeros


				bool *d_exp_Mcoarse_significant(0);

				HANDLE_ERROR(cudaMalloc( (void**) &d_exp_Mcoarse_significant, coarse_num * sizeof(bool)));
				HANDLE_ERROR(cudaMemcpy( d_exp_Mcoarse_significant, &(exp_Mcoarse_significant.data),  coarse_num * sizeof(bool), cudaMemcpyHostToDevice));
//
//				int *d_rotidx(0);
//				HANDLE_ERROR(cudaMalloc( (void**) &d_rotidx, significant_num * sizeof(int)));
//				HANDLE_ERROR(cudaMemcpy( d_rotidx, rotidx,  significant_num * sizeof(int), cudaMemcpyHostToDevice));

				/*====================================
				    		Kernel Calls
				======================================*/
				short int orient1, orient2;

				if(orientation_num>65535)
				{
					orient1 = ceil(sqrt(orientation_num));
					orient2 = orient1;
				}
				else
				{
					orient1 = orientation_num;
					orient2 = 1;
				}
				dim3 block_dim(orient1,orient2,translation_num);

				cudaEvent_t start, stop;
				float time;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);
				cudaEventRecord(start, 0);
				//printf("Calling kernel with <<(%d,%d), %d>> \n", block_dim.x, block_dim.y, BLOCK_SIZE);
				cuda_kernel_massive_diff2<<<block_dim,BLOCK_SIZE>>>(d_Frefs, d_Fimgs, d_Minvsigma2, d_diff2s,
																	Frefs.xy, exp_highres_Xi2_imgs[ipart] / 2.,
																	d_exp_Mcoarse_significant,
																	orientation_num,
																	translation_num,
																	exp_nr_oversampled_rot,
																	exp_nr_oversampled_trans);

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				cudaEventElapsedTime(&time, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
//				for (long unsigned i = 0; i < orientation_num; i ++)
//				{
//					for (long unsigned j = 0; j < translation_num; j ++)
//					{
//						cuda_diff2_deviceImage( Frefs.xy, (double*) ( d_Frefs + (i * Frefs.xy) ), (double*) ( d_Fimgs + (j * Fimgs.xy) ), d_Minvsigma2, d_diff2s + (i * translation_num + j));
//					}
//				}

				/*====================================
				    	   Retrieve Results
				======================================*/

				HANDLE_ERROR(cudaDeviceSynchronize());
				//printf("Kernel call finished \n");

				std::cerr << "It took "<< time <<" msecs."<< std::endl;
				double* diff2s = new double[orientation_num*translation_num];
				if (exp_ipass == 0)
				{
					exp_Mcoarse_significant.clear();
				}
				HANDLE_ERROR(cudaMemcpy( diff2s, d_diff2s, orientation_num*translation_num*sizeof(double), cudaMemcpyDeviceToHost ));

				/*====================================
				    	Write To Destination
				======================================*/


//				if (exp_current_oversampling > 1)
//				{
//				std::ofstream myfile;
//				std::stringstream sstm;
//				sstm << "diff2s/gpu_part.dat";
//				myfile.open(sstm.str().c_str(), std::ios_base::app);
//				}

				//printf("Writing to destination \n");
				for (long int i = 0; i < orientation_num; i++)
				{
					long int iover_rot = iover_rots[i];

					for (long int j = 0; j < translation_num; j++)
					{
						long int ihidden = iorientclasses[i] * exp_nr_trans + ihiddens[j];
						long int iover_trans = iover_transes[j];

						long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden, exp_current_oversampling,
																							iover_rot, iover_trans);

						double diff2 = diff2s[i * translation_num + j];
						//diff2 += exp_highres_Xi2_imgs[ipart] / 2.;
//
//						if (exp_current_oversampling > 1)
//							myfile << ihidden_over << " " << diff2 << std::endl;

						DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) = diff2;

						// Keep track of minimum of all diff2, only for the last image in this series
						if (diff2 < exp_min_diff2[ipart])
							exp_min_diff2[ipart] = diff2;
					}
				}
				//printf("Writing to destination finished \n");


				cudaFree(d_Fimgs);
				cudaFree(d_diff2s);
				delete [] diff2s;

			} // end loop ipart

			cudaFree(d_Frefs);

		} // end if class significant
	} // end loop iclass
}

void MlOptimiserCUDA::storeWeightedSums(long int my_ori_particle, int exp_current_image_size,
		int exp_current_oversampling, int metadata_offset,
		int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
		int exp_itrans_min, int exp_itrans_max, int exp_iclass_min, int exp_iclass_max,
		std::vector<double> &exp_min_diff2,
		std::vector<double> &exp_highres_Xi2_imgs,
		std::vector<MultidimArray<Complex > > &exp_Fimgs,
		std::vector<MultidimArray<Complex > > &exp_Fimgs_nomask,
		std::vector<MultidimArray<double> > &exp_Fctfs,
		std::vector<MultidimArray<double> > &exp_power_imgs,
		std::vector<Matrix1D<double> > &exp_old_offset,
		std::vector<Matrix1D<double> > &exp_prior,
		MultidimArray<double> &exp_Mweight,
		MultidimArray<bool> &exp_Mcoarse_significant,
		std::vector<double> &exp_significant_weight,
		std::vector<double> &exp_sum_weight,
		std::vector<double> &exp_max_weight,
		std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
		std::vector<double> &exp_directions_prior, std::vector<double> &exp_psi_prior,
		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted_nomask,
		std::vector<MultidimArray<double> > &exp_local_Minvsigma2s,
		std::vector<MultidimArray<double> > &exp_local_Fctfs,
		std::vector<double> &exp_local_sqrtXi2)
{

	int exp_nr_particles = mydata.ori_particles[my_ori_particle].particles_id.size();
	long int exp_nr_dir = (do_skip_align || do_skip_rotate) ? 1 : sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
	long int exp_nr_psi = (do_skip_align || do_skip_rotate) ? 1 : sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior);
	long int exp_nr_trans = (do_skip_align) ? 1 : sampling.NrTranslationalSamplings();
	long int exp_nr_oversampled_rot = sampling.oversamplingFactorOrientations(exp_current_oversampling);
	long int exp_nr_oversampled_trans = sampling.oversamplingFactorTranslations(exp_current_oversampling);

	// Re-do below because now also want unmasked images AND if (stricht_highres_exp >0.) then may need to resize
	precalculateShiftedImagesCtfsAndInvSigma2s(true, my_ori_particle, exp_current_image_size, exp_current_oversampling,
			exp_itrans_min, exp_itrans_max, exp_Fimgs, exp_Fimgs_nomask, exp_Fctfs, exp_local_Fimgs_shifted, exp_local_Fimgs_shifted_nomask,
			exp_local_Fctfs, exp_local_sqrtXi2, exp_local_Minvsigma2s);

	// In doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s() the origin of the exp_local_Minvsigma2s was omitted.
	// Set those back here
	for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
	{
		long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
		int group_id = mydata.getGroupId(part_id);
		DIRECT_MULTIDIM_ELEM(exp_local_Minvsigma2s[ipart], 0) = 1. / (sigma2_fudge * DIRECT_A1D_ELEM(mymodel.sigma2_noise[group_id], 0));
	}

	// Initialise the maximum of all weights to a negative value
	exp_max_weight.clear();
	exp_max_weight.resize(exp_nr_particles, -1.);

	// For norm_correction and scale_correction of all particles of this ori_particle
	std::vector<double> exp_wsum_norm_correction;
	std::vector<MultidimArray<double> > exp_wsum_scale_correction_XA, exp_wsum_scale_correction_AA;
	std::vector<MultidimArray<double> > thr_wsum_signal_product_spectra, thr_wsum_reference_power_spectra;
	exp_wsum_norm_correction.resize(exp_nr_particles, 0.);

	// For scale_correction
	if (do_scale_correction)
	{
		MultidimArray<double> aux;
		aux.initZeros(mymodel.ori_size/2 + 1);
		exp_wsum_scale_correction_XA.resize(exp_nr_particles, aux);
		exp_wsum_scale_correction_AA.resize(exp_nr_particles, aux);
		thr_wsum_signal_product_spectra.resize(mymodel.nr_groups, aux);
		thr_wsum_reference_power_spectra.resize(mymodel.nr_groups, aux);
	}

	std::vector< double> oversampled_rot, oversampled_tilt, oversampled_psi;
	std::vector<double> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
	Matrix2D<double> A;
	MultidimArray<Complex > Fimg, Fref, Frefctf, Fimg_otfshift, Fimg_otfshift_nomask;
	MultidimArray<double> Minvsigma2, Mctf, Fweight;
	double rot, tilt, psi;
	bool have_warned_small_scale = false;
	// Initialising... exp_Fimgs[0] has mymodel.current_size (not coarse_size!)
	Fref.resize(exp_Fimgs[0]);
	Frefctf.resize(exp_Fimgs[0]);
	Fweight.resize(exp_Fimgs[0]);
	Fimg.resize(exp_Fimgs[0]);
	// Initialise Mctf to all-1 for if !do_ctf_corection
	Mctf.resize(exp_Fimgs[0]);
	Mctf.initConstant(1.);
	// Initialise Minvsigma2 to all-1 for if !do_map
	Minvsigma2.resize(exp_Fimgs[0]);
	Minvsigma2.initConstant(1.);
	if (do_shifts_onthefly)
	{
		Fimg_otfshift.resize(Frefctf);
		Fimg_otfshift_nomask.resize(Frefctf);
	}

	// Make local copies of weighted sums (except BPrefs, which are too big)
	// so that there are not too many mutex locks below
	std::vector<MultidimArray<double> > thr_wsum_sigma2_noise, thr_wsum_pdf_direction;
	std::vector<double> thr_wsum_norm_correction, thr_sumw_group, thr_wsum_pdf_class, thr_wsum_prior_offsetx_class, thr_wsum_prior_offsety_class;
	double thr_wsum_sigma2_offset;
	MultidimArray<double> thr_metadata, zeroArray;
	// Wsum_sigma_noise2 is a 1D-spectrum for each group
	zeroArray.initZeros(mymodel.ori_size/2 + 1);
	thr_wsum_sigma2_noise.resize(mymodel.nr_groups, zeroArray);
	// wsum_pdf_direction is a 1D-array (of length sampling.NrDirections()) for each class
	zeroArray.initZeros(sampling.NrDirections());
	thr_wsum_pdf_direction.resize(mymodel.nr_classes, zeroArray);
	// sumw_group is a double for each group
	thr_sumw_group.resize(mymodel.nr_groups, 0.);
	// wsum_pdf_class is a double for each class
	thr_wsum_pdf_class.resize(mymodel.nr_classes, 0.);
	if (mymodel.ref_dim == 2)
	{
		thr_wsum_prior_offsetx_class.resize(mymodel.nr_classes, 0.);
		thr_wsum_prior_offsety_class.resize(mymodel.nr_classes, 0.);
	}
	// wsum_sigma2_offset is just a double
	thr_wsum_sigma2_offset = 0.;


	printf("Entering crazy-ass loop\n");

	// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++)
	{
		for (long int idir = exp_idir_min, iorient = 0; idir <= exp_idir_max; idir++)
		{
			for (long int ipsi = exp_ipsi_min; ipsi <= exp_ipsi_max; ipsi++, iorient++)
			{
				long int iorientclass = exp_iclass * exp_nr_dir * exp_nr_psi + iorient;

				// Only proceed if any of the particles had any significant coarsely sampled translation
				if (isSignificantAnyParticleAnyTranslation(iorientclass, exp_itrans_min, exp_itrans_max, exp_Mcoarse_significant))
				{
					// Now get the oversampled (rot, tilt, psi) triplets
					// This will be only the original (rot,tilt,psi) triplet if (adaptive_oversampling==0)
					sampling.getOrientations(idir, ipsi, adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
							exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior);
					// Loop over all oversampled orientations (only a single one in the first pass)
					for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
					{
						rot = oversampled_rot[iover_rot];
						tilt = oversampled_tilt[iover_rot];
						psi = oversampled_psi[iover_rot];
						// Get the Euler matrix
						Euler_angles2matrix(rot, tilt, psi, A);
						// Project the reference map (into Fref)
						if (!do_skip_maximization)
							(mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_NOT_INV);
						// Inside the loop over all translations and all part_id sum all shift Fimg's and their weights
						// Then outside this loop do the actual backprojection
						Fimg.initZeros();
						Fweight.initZeros();
						/// Now that reference projection has been made loop over all particles inside this ori_particle
						for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
						{
							// This is an attempt to speed up illogically slow updates of wsum_sigma2_offset....
							// It seems to make a big difference!
							double myprior_x, myprior_y, myprior_z, old_offset_z;
							double old_offset_x = XX(exp_old_offset[ipart]);
							double old_offset_y = YY(exp_old_offset[ipart]);
							if (mymodel.ref_dim == 2)
							{
								myprior_x = XX(mymodel.prior_offset_class[exp_iclass]);
								myprior_y = YY(mymodel.prior_offset_class[exp_iclass]);
							}
							else
							{
								myprior_x = XX(exp_prior[ipart]);
								myprior_y = YY(exp_prior[ipart]);
								if (mymodel.data_dim == 3)
								{
									myprior_z = ZZ(exp_prior[ipart]);
									old_offset_z = ZZ(exp_old_offset[ipart]);
								}
							}

							long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
							int group_id = mydata.getGroupId(part_id);
							if (!do_skip_maximization)
							{
								if (do_map)
									Minvsigma2 = exp_local_Minvsigma2s[ipart];
								// else Minvsigma2 was initialised to ones
								// Apply CTF to reference projection
								if (do_ctf_correction)
								{
									Mctf = exp_local_Fctfs[ipart];
									if (refs_are_ctf_corrected)
									{
										FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
										{
											DIRECT_MULTIDIM_ELEM(Frefctf, n) = DIRECT_MULTIDIM_ELEM(Fref, n) * DIRECT_MULTIDIM_ELEM(Mctf, n);
										}
									}
									else
									{
										Frefctf = Fref;
									}
								}
								else
								{
									// initialise because there are multiple particles and Mctf gets selfMultiplied for scale_correction
									Mctf.initConstant(1.);
									Frefctf = Fref;
								}
								if (do_scale_correction)
								{
									double myscale = mymodel.scale_correction[group_id];
									if (myscale > 10000.)
									{
										std::cerr << " rlnMicrographScaleCorrection= " << myscale << " group= " << group_id + 1 << std::endl;
										REPORT_ERROR("ERROR: rlnMicrographScaleCorrection is very high. Did you normalize your data?");
									}
									else if (myscale < 0.001)
									{
										if (!have_warned_small_scale)
										{
											std::cout << " WARNING: ignoring group " << group_id + 1 << " with very small or negative scale (" << myscale <<
													"); Use larger groups for more stable scale estimates." << std::endl;
											have_warned_small_scale = true;
										}
										myscale = 0.001;
									}
									FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
									{
										DIRECT_MULTIDIM_ELEM(Frefctf, n) *= myscale;
									}
									// For CTF-terms in BP
									Mctf *= myscale;
								}
							} // end if !do_skip_maximization

							long int ihidden = iorientclass * exp_nr_trans;
							for (long int itrans = exp_itrans_min, iitrans = 0; itrans <= exp_itrans_max; itrans++, ihidden++)
							{
								sampling.getTranslations(itrans, adaptive_oversampling,
										oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);
								for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++, iitrans++)
								{
									// Only deal with this sampling point if its weight was significant
									long int ihidden_over = ihidden * exp_nr_oversampled_trans * exp_nr_oversampled_rot +
											iover_rot * exp_nr_oversampled_trans + iover_trans;
									double weight = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
									// Only sum weights for non-zero weights
									if (weight >= exp_significant_weight[ipart])
									{
										// Normalise the weight (do this after the comparison with exp_significant_weight!)
										weight /= exp_sum_weight[ipart];
										if (!do_skip_maximization)
										{

											/// Now get the shifted image
											// Use a pointer to avoid copying the entire array again in this highly expensive loop
											Complex *Fimg_shift, *Fimg_shift_nomask;
											if (!do_shifts_onthefly)
											{
												long int ishift = ipart * exp_nr_oversampled_trans * exp_nr_trans + iitrans;
												Fimg_shift = exp_local_Fimgs_shifted[ishift].data;
												Fimg_shift_nomask = exp_local_Fimgs_shifted_nomask[ishift].data;
											}
											else
											{
												Complex* myAB;
												myAB = (adaptive_oversampling == 0 ) ? global_fftshifts_ab_current[iitrans].data : global_fftshifts_ab2_current[iitrans].data;
												FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Fimgs_shifted[ipart])
												{
													double a = (*(myAB + n)).real;
													double b = (*(myAB + n)).imag;
													// Fimg_shift
													double real = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real
															- b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag;
													double imag = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag
															+ b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real;
													DIRECT_MULTIDIM_ELEM(Fimg_otfshift, n) = Complex(real, imag);
													// Fimg_shift_nomask
													real = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).real
															- b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).imag;
													imag = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).imag
															+ b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).real;
													DIRECT_MULTIDIM_ELEM(Fimg_otfshift_nomask, n) = Complex(real, imag);
												}
												Fimg_shift = Fimg_otfshift.data;
												Fimg_shift_nomask = Fimg_otfshift_nomask.data;
											}

											// Store weighted sum of squared differences for sigma2_noise estimation
											// Suggestion Robert Sinkovitz: merge difference and scale steps to make better use of cache
											FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
											{
												int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
												if (ires > -1)
												{
													// Use FT of masked image for noise estimation!
													double diff_real = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real - (*(Fimg_shift + n)).real;
													double diff_imag = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag - (*(Fimg_shift + n)).imag;
													double wdiff2 = weight * (diff_real*diff_real + diff_imag*diff_imag);
													// group-wise sigma2_noise
													DIRECT_MULTIDIM_ELEM(thr_wsum_sigma2_noise[group_id], ires) += wdiff2;
													// For norm_correction
													exp_wsum_norm_correction[ipart] += wdiff2;
												}
											    if (do_scale_correction && DIRECT_A1D_ELEM(mymodel.data_vs_prior_class[exp_iclass], ires) > 3.)
												{
											    	double sumXA, sumA2;
											    	sumXA = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real * (*(Fimg_shift + n)).real;
											    	sumXA += (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag * (*(Fimg_shift + n)).imag;
											    	DIRECT_A1D_ELEM(exp_wsum_scale_correction_XA[ipart], ires) += weight * sumXA;
											    	sumA2 = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real * (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real;
											    	sumA2 += (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag * (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag;
											    	DIRECT_A1D_ELEM(exp_wsum_scale_correction_AA[ipart], ires) += weight * sumA2;
												}
											}

											// Store sum of weights for this group
											thr_sumw_group[group_id] += weight;
											// Store weights for this class and orientation
											thr_wsum_pdf_class[exp_iclass] += weight;

											// The following goes MUCH faster than the original lines below....
											if (mymodel.ref_dim == 2)
											{
												thr_wsum_prior_offsetx_class[exp_iclass] += weight * (old_offset_x + oversampled_translations_x[iover_trans]);
												thr_wsum_prior_offsety_class[exp_iclass] += weight * (old_offset_y + oversampled_translations_y[iover_trans]);
											}
											double diffx = myprior_x - old_offset_x - oversampled_translations_x[iover_trans];
											double diffy = myprior_y - old_offset_y - oversampled_translations_y[iover_trans];
											if (mymodel.data_dim == 3)
											{
												double diffz  = myprior_z - old_offset_z - oversampled_translations_z[iover_trans];
												thr_wsum_sigma2_offset += weight * (diffx*diffx + diffy*diffy + diffz*diffz);
											}
											else
											{
												thr_wsum_sigma2_offset += weight * (diffx*diffx + diffy*diffy);
											}

											// Store weight for this direction of this class
											if (do_skip_align || do_skip_rotate )
											{
												//ignore pdf_direction
											}
											else if (mymodel.orientational_prior_mode == NOPRIOR)
											{
												DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], idir) += weight;
											}
											else
											{
												// In the case of orientational priors, get the original number of the direction back
												long int mydir = exp_pointer_dir_nonzeroprior[idir];
												DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], mydir) += weight;
											}
											// Store sum of weight*SSNR*Fimg in data and sum of weight*SSNR in weight
											// Use the FT of the unmasked image to back-project in order to prevent reconstruction artefacts! SS 25oct11
											FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
											{
												double myctf = DIRECT_MULTIDIM_ELEM(Mctf, n);
												// Note that weightxinvsigma2 already contains the CTF!
												double weightxinvsigma2 = weight * myctf * DIRECT_MULTIDIM_ELEM(Minvsigma2, n);
												// now Fimg stores sum of all shifted w*Fimg
												(DIRECT_MULTIDIM_ELEM(Fimg, n)).real += (*(Fimg_shift_nomask + n)).real * weightxinvsigma2;
												(DIRECT_MULTIDIM_ELEM(Fimg, n)).imag += (*(Fimg_shift_nomask + n)).imag * weightxinvsigma2;
												// now Fweight stores sum of all w
												// Note that CTF needs to be squared in Fweight, weightxinvsigma2 already contained one copy
												DIRECT_MULTIDIM_ELEM(Fweight, n) += weightxinvsigma2 * myctf;
											}
										} // end if !do_skip_maximization

										// Keep track of max_weight and the corresponding optimal hidden variables
										if (weight > exp_max_weight[ipart])
										{
											// Store optimal image parameters
											exp_max_weight[ipart] = weight;

											// TODO: remove, for now to maintain exact numerical version of old threads....
											A = A.inv();
											A = A.inv();
											Euler_matrix2angles(A, rot, tilt, psi);

											DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ROT) = rot;
											DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_TILT) = tilt;
											DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PSI) = psi;
											DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_XOFF) = XX(exp_old_offset[ipart]) + oversampled_translations_x[iover_trans];
											DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_YOFF) = YY(exp_old_offset[ipart]) + oversampled_translations_y[iover_trans];
											if (mymodel.data_dim == 3)
												DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ZOFF) = ZZ(exp_old_offset[ipart]) + oversampled_translations_z[iover_trans];
											DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_CLASS) = (double)exp_iclass + 1;
											DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PMAX) = exp_max_weight[ipart];
										}
									} // end if weight >= exp_significant_weight
								} // end loop iover_trans
							} // end loop itrans
						} // end loop ipart

						if (!do_skip_maximization)
						{
							// Perform the actual back-projection.
							// This is done with the sum of all (in-plane) shifted Fimg's
							// Perform this inside a mutex
							int my_mutex = exp_iclass % NR_CLASS_MUTEXES;
							pthread_mutex_lock(&global_mutex2[my_mutex]);
							(wsum_model.BPref[exp_iclass]).set2DFourierTransform(Fimg, A, IS_NOT_INV, &Fweight);
							pthread_mutex_unlock(&global_mutex2[my_mutex]);
						} // end if !do_skip_maximization
					} // end loop iover_rot
				}// end loop do_proceed
			} // end loop ipsi
		} // end loop idir
	} // end loop iclass

	printf("Exiting crazy-ass loop\n");

	// Extend norm_correction and sigma2_noise estimation to higher resolutions for all particles
	// Also calculate dLL for each particle and store in metadata
	// loop over all particles inside this ori_particle
	double thr_avg_norm_correction = 0.;
	double thr_sum_dLL = 0., thr_sum_Pmax = 0.;
	for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
	{
		long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
		int group_id = mydata.getGroupId(part_id);

		// If the current images were smaller than the original size, fill the rest of wsum_model.sigma2_noise with the power_class spectrum of the images
		for (int ires = mymodel.current_size/2 + 1; ires < mymodel.ori_size/2 + 1; ires++)
		{
			DIRECT_A1D_ELEM(thr_wsum_sigma2_noise[group_id], ires) += DIRECT_A1D_ELEM(exp_power_imgs[ipart], ires);
			// Also extend the weighted sum of the norm_correction
			exp_wsum_norm_correction[ipart] += DIRECT_A1D_ELEM(exp_power_imgs[ipart], ires);
		}

		// Store norm_correction
		// Multiply by old value because the old norm_correction term was already applied to the image
		if (do_norm_correction)
		{
			double old_norm_correction = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_NORM);
			old_norm_correction /= mymodel.avg_norm_correction;
			// The factor two below is because exp_wsum_norm_correctiom is similar to sigma2_noise, which is the variance for the real/imag components
			// The variance of the total image (on which one normalizes) is twice this value!
			double normcorr = old_norm_correction * sqrt(exp_wsum_norm_correction[ipart] * 2.);
			thr_avg_norm_correction += normcorr;
			// Now set the new norm_correction in the relevant position of exp_metadata
			DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_NORM) = normcorr;

			// Print warning for strange norm-correction values
			if (!(iter == 1 && do_firstiter_cc) && DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_NORM) > 10.)
			{
				std::cout << " WARNING: norm_correction= "<< DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_NORM) << " for particle " << part_id << " in group " << group_id + 1 << "; Are your groups large enough?" << std::endl;
			}

		}

		// Store weighted sums for scale_correction
		if (do_scale_correction)
		{
			// Divide XA by the old scale_correction and AA by the square of that, because was incorporated into Fctf
			exp_wsum_scale_correction_XA[ipart] /= mymodel.scale_correction[group_id];
			exp_wsum_scale_correction_AA[ipart] /= mymodel.scale_correction[group_id] * mymodel.scale_correction[group_id];

			thr_wsum_signal_product_spectra[group_id] += exp_wsum_scale_correction_XA[ipart];
			thr_wsum_reference_power_spectra[group_id] += exp_wsum_scale_correction_AA[ipart];
		}

		// Calculate DLL for each particle
		double logsigma2 = 0.;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
		{
			int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
			// Note there is no sqrt in the normalisation term because of the 2-dimensionality of the complex-plane
			// Also exclude origin from logsigma2, as this will not be considered in the P-calculations
			if (ires > 0)
				logsigma2 += log( 2. * PI * DIRECT_A1D_ELEM(mymodel.sigma2_noise[group_id], ires));
		}
		if (exp_sum_weight[ipart]==0)
		{
			std::cerr << " part_id= " << part_id << std::endl;
			std::cerr << " ipart= " << ipart << std::endl;
			std::cerr << " exp_min_diff2[ipart]= " << exp_min_diff2[ipart] << std::endl;
			std::cerr << " logsigma2= " << logsigma2 << std::endl;
			int group_id = mydata.getGroupId(part_id);
			std::cerr << " group_id= " << group_id << std::endl;
			std::cerr << " ml_model.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
			std::cerr << " exp_significant_weight[ipart]= " << exp_significant_weight[ipart] << std::endl;
			std::cerr << " exp_max_weight[ipart]= " << exp_max_weight[ipart] << std::endl;
			std::cerr << " ml_model.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << std::endl;
			REPORT_ERROR("ERROR: exp_sum_weight[ipart]==0");
		}
		double dLL;
		if ((iter==1 && do_firstiter_cc) || do_always_cc)
			dLL = -exp_min_diff2[ipart];
		else
			dLL = log(exp_sum_weight[ipart]) - exp_min_diff2[ipart] - logsigma2;

		// Store dLL of each image in the output array, and keep track of total sum
		DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_DLL) = dLL;
		thr_sum_dLL += dLL;

		// Also store sum of Pmax
		thr_sum_Pmax += DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PMAX);

	}

	// Now, inside a global_mutex, update the other weighted sums among all threads
	if (!do_skip_maximization)
	{
		pthread_mutex_lock(&global_mutex);
		for (int n = 0; n < mymodel.nr_groups; n++)
		{
			wsum_model.sigma2_noise[n] += thr_wsum_sigma2_noise[n];
			wsum_model.sumw_group[n] += thr_sumw_group[n];
			if (do_scale_correction)
			{
				wsum_model.wsum_signal_product_spectra[n] += thr_wsum_signal_product_spectra[n];
				wsum_model.wsum_reference_power_spectra[n] += thr_wsum_reference_power_spectra[n];
			}
		}
		for (int n = 0; n < mymodel.nr_classes; n++)
		{
			wsum_model.pdf_class[n] += thr_wsum_pdf_class[n];
			if (mymodel.ref_dim == 2)
			{
				XX(wsum_model.prior_offset_class[n]) += thr_wsum_prior_offsetx_class[n];
				YY(wsum_model.prior_offset_class[n]) += thr_wsum_prior_offsety_class[n];
			}
			if (!(do_skip_align || do_skip_rotate) )
				wsum_model.pdf_direction[n] += thr_wsum_pdf_direction[n];
		}
		wsum_model.sigma2_offset += thr_wsum_sigma2_offset;
		if (do_norm_correction)
			wsum_model.avg_norm_correction += thr_avg_norm_correction;
		wsum_model.LL += thr_sum_dLL;
		wsum_model.ave_Pmax += thr_sum_Pmax;
		pthread_mutex_unlock(&global_mutex);
	} // end if !do_skip_maximization

	printf("Done doing other stuff\n");
}

//void MlOptimiserCUDA::precalculateModelProjectionsCtfsAndInvSigma2s(bool do_also_unmasked,
//	    int exp_current_image_size, int exp_current_oversampling,
//		std::vector<MultidimArray<Complex > > &model,
//		std::vector<MultidimArray<Complex > > &exp_Fimgs_nomask,
//		std::vector<MultidimArray<double> > &exp_Fctfs,
//		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
//		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted_nomask,
//		std::vector<MultidimArray<double> > &exp_local_Fctfs,
//		std::vector<double> &exp_local_sqrtXi2,
//		std::vector<MultidimArray<double> > &exp_local_Minvsigma2s)
//{
//
//}

void MlOptimiserCUDA::precalculateShiftedImagesCtfsAndInvSigma2s(bool do_also_unmasked,
		long int my_ori_particle, int exp_current_image_size, int exp_current_oversampling,
		int exp_itrans_min, int exp_itrans_max,
		std::vector<MultidimArray<Complex > > &exp_Fimgs,
		std::vector<MultidimArray<Complex > > &exp_Fimgs_nomask,
		std::vector<MultidimArray<double> > &exp_Fctfs,
		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted_nomask,
		std::vector<MultidimArray<double> > &exp_local_Fctfs,
		std::vector<double> &exp_local_sqrtXi2,
		std::vector<MultidimArray<double> > &exp_local_Minvsigma2s)
{

	int exp_nr_particles = mydata.ori_particles[my_ori_particle].particles_id.size();
	int nr_shifts = (do_shifts_onthefly || do_skip_align) ? exp_nr_particles : exp_nr_particles * sampling.NrTranslationalSamplings(exp_current_oversampling);
	// Don't re-do if nothing has changed....
	bool do_ctf_invsig = (exp_local_Fctfs.size() > 0) ? YSIZE(exp_local_Fctfs[0])  != exp_current_image_size : true; // size has changed
	bool do_masked_shifts = (do_ctf_invsig || nr_shifts != exp_local_Fimgs_shifted.size()); // size or nr_shifts has changed

	// Use pre-sized vectors instead of push_backs!!
	exp_local_Fimgs_shifted.resize(nr_shifts);
	if (do_also_unmasked)
		exp_local_Fimgs_shifted_nomask.resize(nr_shifts);
	exp_local_Minvsigma2s.resize(exp_nr_particles);
	exp_local_Fctfs.resize(exp_nr_particles);
	exp_local_sqrtXi2.resize(exp_nr_particles);

	MultidimArray<Complex > Fimg, Fimg_nomask;
	for (int ipart = 0, my_trans_image = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
	{
		long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
		int group_id = mydata.getGroupId(part_id);

		if (do_masked_shifts)
			windowFourierTransform(exp_Fimgs[ipart], Fimg, exp_current_image_size);
		if (do_also_unmasked)
			windowFourierTransform(exp_Fimgs_nomask[ipart], Fimg_nomask, exp_current_image_size);

		if (do_ctf_invsig)
		{
			// Also precalculate the sqrt of the sum of all Xi2
			// Could exp_current_image_size ever be different from mymodel.current_size?
			// Probably therefore do it here rather than in getFourierTransforms
			if ((iter == 1 && do_firstiter_cc) || do_always_cc)
			{
				double sumxi2 = 0.;
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
				{
					sumxi2 += norm(DIRECT_MULTIDIM_ELEM(Fimg, n));
				}
				// Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
				exp_local_sqrtXi2[ipart] = sqrt(sumxi2);
			}

			// Also store downsized Fctfs
			// In the second pass of the adaptive approach this will have no effect,
			// since then exp_current_image_size will be the same as the size of exp_Fctfs
			windowFourierTransform(exp_Fctfs[ipart], exp_local_Fctfs[ipart], exp_current_image_size);

			// Also prepare Minvsigma2
			if (mymodel.data_dim == 3)
				exp_local_Minvsigma2s[ipart].initZeros(ZSIZE(Fimg), YSIZE(Fimg), XSIZE(Fimg));
			else
				exp_local_Minvsigma2s[ipart].initZeros(YSIZE(Fimg), XSIZE(Fimg));

			int *myMresol = (YSIZE(Fimg) == coarse_size) ? Mresol_coarse.data : Mresol_fine.data;
			// With group_id and relevant size of Fimg, calculate inverse of sigma^2 for relevant parts of Mresol
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Minvsigma2s[ipart])
			{
				int ires = *(myMresol + n);
				// Exclude origin (ires==0) from the Probability-calculation
				// This way we are invariant to additive factors
				if (ires > 0)
					DIRECT_MULTIDIM_ELEM(exp_local_Minvsigma2s[ipart], n) = 1. / (sigma2_fudge * DIRECT_A1D_ELEM(mymodel.sigma2_noise[group_id], ires));
			}

		}

		if (do_shifts_onthefly)
		{
			// Store a single, down-sized version of exp_Fimgs[ipart] in exp_local_Fimgs_shifted
			if (do_masked_shifts)
				exp_local_Fimgs_shifted[ipart] = Fimg;
			if (do_also_unmasked)
				exp_local_Fimgs_shifted_nomask[ipart] = Fimg_nomask;
		}
		else
		{
			// Store all translated variants of Fimg
			for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++)
			{
				// First get the non-oversampled translations as defined by the sampling object
				std::vector<double > oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
				sampling.getTranslations(itrans, exp_current_oversampling, oversampled_translations_x,
						oversampled_translations_y, oversampled_translations_z);
				// Then loop over all its oversampled relatives
				for (long int iover_trans = 0; iover_trans < oversampled_translations_x.size(); iover_trans++, my_trans_image++)
				{
					// Shift through phase-shifts in the Fourier transform
					// Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
					if (do_masked_shifts)
					{
						exp_local_Fimgs_shifted[my_trans_image].resize(Fimg);
						if (mymodel.data_dim ==2)
							shiftImageInFourierTransform(Fimg, exp_local_Fimgs_shifted[my_trans_image],
									tab_sin, tab_cos, (double)mymodel.ori_size,
									oversampled_translations_x[iover_trans],
									oversampled_translations_y[iover_trans]);
						else
							shiftImageInFourierTransform(Fimg, exp_local_Fimgs_shifted[my_trans_image],
									tab_sin, tab_cos, (double)mymodel.ori_size,
									oversampled_translations_x[iover_trans],
									oversampled_translations_y[iover_trans],
									oversampled_translations_z[iover_trans]);
					}
					if (do_also_unmasked)
					{
						exp_local_Fimgs_shifted_nomask[my_trans_image].resize(Fimg_nomask);
						if (mymodel.data_dim ==2)
							shiftImageInFourierTransform(Fimg_nomask, exp_local_Fimgs_shifted_nomask[my_trans_image],
								tab_sin, tab_cos, (double)mymodel.ori_size,
								oversampled_translations_x[iover_trans],
								oversampled_translations_y[iover_trans]);
						else
							shiftImageInFourierTransform(Fimg_nomask, exp_local_Fimgs_shifted_nomask[my_trans_image],
								tab_sin, tab_cos, (double)mymodel.ori_size,
								oversampled_translations_x[iover_trans],
								oversampled_translations_y[iover_trans],
								oversampled_translations_z[iover_trans]);
					}
				}
			}
		}
	}
}

bool MlOptimiserCUDA::isSignificantAnyParticleAnyTranslation(long int iorient, int exp_itrans_min, int exp_itrans_max, MultidimArray<bool> &exp_Mcoarse_significant)
{

	long int exp_nr_trans = exp_itrans_max - exp_itrans_min + 1;
	for (long int ipart = 0; ipart < YSIZE(exp_Mcoarse_significant); ipart++)
	{
		long int ihidden = iorient * exp_nr_trans;
		for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++, ihidden++)
		{
			if (DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden))
				return true;
		}
	}
	return false;

}
