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

#define MAX_RESOL_SHARED_MEM 32
#define PIXELWISE_BLOCK_SIZE 128         // This is optimally set as big as possible without its ceil:ed multiple exceeding imagesize by too much.
#define TRANSWISE_BLOCK_SIZE 32			// Most optimal is 32 when doing 84 nr translations
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
		__shared__ double s[PIXELWISE_BLOCK_SIZE];
		s[threadIdx.x] = 0;

		unsigned pass_num(ceilf((float)img_size/(float)PIXELWISE_BLOCK_SIZE));
		unsigned long pixel,
		ref_start(ex * img_size),
		img_start(ez * img_size);

		unsigned long ref_pixel_idx;
		unsigned long img_pixel_idx;

		for (unsigned pass = 0; pass < pass_num; pass ++)
		{
			pixel = pass * PIXELWISE_BLOCK_SIZE + threadIdx.x;

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
		int itr = PIXELWISE_BLOCK_SIZE/trads;
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
				cuda_kernel_massive_diff2<<<block_dim,PIXELWISE_BLOCK_SIZE>>>(d_Frefs, d_Fimgs, d_Minvsigma2, d_diff2s,
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


//				std::ofstream myfile;
//				std::stringstream sstm;
//				sstm << "diff2s/gpu_part.dat";
//				myfile.open(sstm.str().c_str(), std::ios_base::app);
//				myfile << ihidden_over << " " << diff2 << std::endl;

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

__global__ void cuda_kernel_wavg(	CudaComplex *g_refs, CudaComplex *g_imgs, CudaComplex *g_imgs_nomask,
									double* g_weights, double* g_ctfs, double* g_Minvsigma2s,
									double *g_wdiff2s_parts, CudaComplex *g_wavgs, double* g_Fweights,
									const unsigned translation_num, const double weight_norm,
									const double significant_weight)
{
	unsigned otientation = blockIdx.x;
	unsigned pixel = blockIdx.y;
	unsigned image_size = gridDim.y;
	unsigned tid = threadIdx.x;

	//TODO Consider Mresol_fine to speed this kernel

	__shared__ double s_wdiff2[TRANSWISE_BLOCK_SIZE];
	__shared__ double s_weight[TRANSWISE_BLOCK_SIZE];
	__shared__ double s_real[TRANSWISE_BLOCK_SIZE];
	__shared__ double s_imag[TRANSWISE_BLOCK_SIZE];

	s_wdiff2[tid] = 0;
	s_weight[tid] = 0;
	s_real[tid] = 0;
	s_imag[tid] = 0;

	unsigned pass_num(ceilf((float)translation_num/(float)TRANSWISE_BLOCK_SIZE)),translation;

	for (unsigned pass = 0; pass < pass_num; pass ++)
	{
		translation = pass * TRANSWISE_BLOCK_SIZE + tid;

		if (translation < translation_num)
		{
			double weight = g_weights[otientation * translation_num + translation];

			if (weight >= significant_weight)
			{
				weight /= weight_norm;

				unsigned long ref_pixel_idx = otientation * image_size + pixel;
				unsigned long img_pixel_idx = translation * image_size + pixel;

				double diff_real = g_refs[ref_pixel_idx].real - g_imgs[img_pixel_idx].real;
				double diff_imag = g_refs[ref_pixel_idx].imag - g_imgs[img_pixel_idx].imag;
				double wdiff2 = weight * (diff_real*diff_real + diff_imag*diff_imag);

				s_wdiff2[tid] += wdiff2;

				double myctf = g_ctfs[pixel];
				double weightxinvsigma2 = weight * myctf * g_Minvsigma2s[pixel];

				s_real[tid] += g_imgs_nomask[img_pixel_idx].real * weightxinvsigma2;
				s_imag[tid] += g_imgs_nomask[img_pixel_idx].imag * weightxinvsigma2;

				s_weight[tid] += weightxinvsigma2 * myctf;
			}
		}
	}

	__syncthreads();

	//TODO Yeah, lets actually use the GPU! Following four sums should be done as a parallel reduction.

	if (tid == 0)
	{
		double sum(0);
		for (unsigned i = 0; i < TRANSWISE_BLOCK_SIZE; i++)
			sum += s_weight[i];

		g_Fweights[otientation * image_size + pixel] += sum;
	}
	else if (tid == 1)
	{
		double sum(0);
		for (unsigned i = 0; i < TRANSWISE_BLOCK_SIZE; i++)
			sum += s_real[i];

		g_wavgs[otientation * image_size + pixel].real += sum;
	}
	else if (tid == 2)
	{
		double sum(0);
		for (unsigned i = 0; i < TRANSWISE_BLOCK_SIZE; i++)
			sum += s_imag[i];

		g_wavgs[otientation * image_size + pixel].imag += sum;
	}
	else if (tid == 3)
	{
		double sum(0);
		for (unsigned i = 0; i < TRANSWISE_BLOCK_SIZE; i++)
			sum += s_wdiff2[i];

		g_wdiff2s_parts[otientation * image_size + pixel] = sum;
	}
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
	bool have_warned_small_scale = false;
	// Initialising... exp_Fimgs[0] has mymodel.current_size (not coarse_size!)
	Fref.resize(exp_Fimgs[0]);
	Frefctf.resize(exp_Fimgs[0]);
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




	/*
	This should be helpful:

	ihidden_over =
	iclass * exp_nr_dir * exp_nr_psi * exp_nr_oversampled_rot * exp_nr_trans * exp_nr_oversampled_trans +
	idir *                exp_nr_psi * exp_nr_oversampled_rot * exp_nr_trans * exp_nr_oversampled_trans +
	ipsi *                             exp_nr_oversampled_rot * exp_nr_trans * exp_nr_oversampled_trans +
	iover_rot *                                                 exp_nr_trans * exp_nr_oversampled_trans +
	itrans *                                                                   exp_nr_oversampled_trans +
	iover_trans
	*/

	unsigned long class_size = exp_nr_dir * exp_nr_psi * exp_nr_oversampled_rot * exp_nr_trans * exp_nr_oversampled_trans;
//	unsigned long dir_size =                exp_nr_psi * exp_nr_oversampled_rot * exp_nr_trans * exp_nr_oversampled_trans;
//	unsigned long psi_size =                             exp_nr_oversampled_rot * exp_nr_trans * exp_nr_oversampled_trans;
//	unsigned long over_rot_size =                                                 exp_nr_trans * exp_nr_oversampled_trans;
//	unsigned long trans_size =                                                                   exp_nr_oversampled_trans;

	/*
	Inverse map:
	iclass = ihidden_over / class_size
	idir =   ihidden_over / class_size
	ipsi =   ihidden_over / psi_size
	...
	*/

	// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++)
	{



		/*=======================================================================================
		                            REFERENCE PROJECTION GENERATION
		=======================================================================================*/


		CudaImages Frefs(exp_local_Minvsigma2s[0].xdim, exp_local_Minvsigma2s[0].ydim, exp_nr_dir * exp_nr_psi * exp_nr_oversampled_rot);

		std::vector< long unsigned > iorientclasses, idirs, iover_rots;
		std::vector< double > rots, tilts, psis;
		long unsigned orientation_num(0);

		Complex* FrefBag = Fref.data; //TODO fix this

		for (long int idir = exp_idir_min, iorient = 0; idir <= exp_idir_max; idir++)
		{
			for (long int ipsi = exp_ipsi_min; ipsi <= exp_ipsi_max; ipsi++, iorient++)
			{
				long int iorientclass = exp_iclass * exp_nr_dir * exp_nr_psi + iorient;

				sampling.getOrientations(idir, ipsi, adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
						exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior);

				if (isSignificantAnyParticleAnyTranslation(iorientclass, exp_itrans_min, exp_itrans_max, exp_Mcoarse_significant))
				{
					for (long unsigned iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
					{
						double rot = oversampled_rot[iover_rot];
						double tilt = oversampled_tilt[iover_rot];
						double psi = oversampled_psi[iover_rot];
						// Get the Euler matrix
						Euler_angles2matrix(rot, tilt, psi, A);

						rots.push_back(rot);
						tilts.push_back(tilt);
						psis.push_back(psi);

						Fref.data = (Complex*) Frefs.current();
						(mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_NOT_INV);

						Frefs.increment();

						orientation_num ++;
						idirs.push_back(idir);
						iorientclasses.push_back(iorientclass);
						iover_rots.push_back(iover_rot);
					}
				}
			}
		}
		Fref.data = FrefBag; //TODO fix this

		CudaComplex *d_Frefs = Frefs.data_to_device();

		CudaComplex *d_wavgs(0); //Weighted image averages for each projection
		HANDLE_ERROR(cudaMalloc( (void**) &d_wavgs, orientation_num * Frefs.xy * sizeof(CudaComplex)));
		HANDLE_ERROR(cudaMemset(d_wavgs, 0, orientation_num * Frefs.xy * sizeof(CudaComplex))); //Initiate zeros

		double *d_Fweights(0);
		HANDLE_ERROR(cudaMalloc( (void**) &d_Fweights, orientation_num * Frefs.xy * sizeof(double)));
		HANDLE_ERROR(cudaMemset(d_Fweights, 0, orientation_num * Frefs.xy * sizeof(double))); //Initiate zeros


		/*=======================================================================================
										  PARTICLE ITERATION
		=======================================================================================*/

		/// Now that reference projection has been made loop over all particles inside this ori_particle
		for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
		{

			long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
			int group_id = mydata.getGroupId(part_id);

			double myprior_x, myprior_y, myprior_z;
			double old_offset_x = XX(exp_old_offset[ipart]);
			double old_offset_y = YY(exp_old_offset[ipart]);
			double old_offset_z;

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


			/*======================================================
								 TRANSLATIONS
			======================================================*/

			CudaImages Fimgs(Frefs.x, Frefs.y, exp_nr_trans * exp_nr_oversampled_trans);
			CudaImages Fimgs_nomask(Frefs.x, Frefs.y, exp_nr_trans * exp_nr_oversampled_trans);

			long unsigned translation_num(0), ihidden(0);
			std::vector< long unsigned > iover_transes, itranses, ihiddens;

			for (long int itrans = exp_itrans_min, iitrans = 0; itrans <= exp_itrans_max; itrans++, ihidden++)
			{
				sampling.getTranslations(itrans, adaptive_oversampling,
						oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);
				for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++, iitrans++)
				{
					/// Now get the shifted image
					// Use a pointer to avoid copying the entire array again in this highly expensive loop
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
						*(Fimgs.current() + n) = CudaComplex(real, imag);

						// Fimg_shift_nomask
						real = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).real
								- b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).imag;
						imag = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).imag
								+ b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).real;
						*(Fimgs_nomask.current() + n) = CudaComplex(real, imag);
					}

					Fimgs.increment();
					Fimgs_nomask.increment();

					translation_num ++;

					ihiddens.push_back(ihidden);
					itranses.push_back(itrans);
					iover_transes.push_back(iover_trans);
				}
			}

			CudaComplex *d_Fimgs = Fimgs.data_to_device();
			CudaComplex *d_Fimgs_nomask = Fimgs_nomask.data_to_device();

			/*======================================================
					            KERNEL CALL
			======================================================*/

			for (int l = 0; l < translation_num; l++)
				std::cerr << DIRECT_A2D_ELEM(exp_Mweight, ipart, l) << std::endl;
			exit(0);

			double *d_weights(0);
			// TODO Is class_size = exp_Mweight.xdim?
			HANDLE_ERROR(cudaMalloc( (void**) &d_weights, class_size * sizeof(double)));
			HANDLE_ERROR(cudaMemcpy( d_weights, exp_Mweight.data + ipart * class_size,
					class_size * sizeof(double), cudaMemcpyHostToDevice));

			double *d_wdiff2s_parts(0);
			HANDLE_ERROR(cudaMalloc( (void**) &d_wdiff2s_parts,
					orientation_num * Frefs.xy * sizeof(double)));

			double *d_ctfs(0);
			HANDLE_ERROR(cudaMalloc( (void**) &d_ctfs, class_size * sizeof(double)));
			HANDLE_ERROR(cudaMemcpy( d_ctfs, exp_local_Fctfs[ipart].data,
					Frefs.xy * sizeof(double), cudaMemcpyHostToDevice));

			double *d_Minvsigma2s(0);
			HANDLE_ERROR(cudaMalloc( (void**) &d_Minvsigma2s, Fimgs.xy * sizeof(double)));
			HANDLE_ERROR(cudaMemcpy( d_Minvsigma2s, exp_local_Minvsigma2s[ipart].data,
					Fimgs.xy * sizeof(double), cudaMemcpyHostToDevice));

			dim3 block_dim(orientation_num, Frefs.xy);
			std::cerr << "cuda_kernel_wavg<<<" << block_dim.x << "," << block_dim.y << ">>>" << std::endl;

			float time;
			cudaEvent_t start, stop;
			cudaEventCreate(&start);
			cudaEventCreate(&stop);
			cudaEventRecord(start, 0);

			cuda_kernel_wavg<<<block_dim,TRANSWISE_BLOCK_SIZE>>>(
												d_Frefs, d_Fimgs, d_Fimgs_nomask,		//INPUT
												d_weights, d_ctfs, d_Minvsigma2s,		//INPUT
												d_wdiff2s_parts, d_wavgs, d_Fweights,	//OUTPUT
												translation_num, exp_sum_weight[ipart],	//CONTANTS
												exp_significant_weight[ipart]
												);

			HANDLE_ERROR(cudaDeviceSynchronize());
			std::cerr << "cuda_kernel_wavg DONE!" << std::endl;


			cudaEventRecord(stop, 0);
			cudaEventSynchronize(stop);
			cudaEventElapsedTime(&time, start, stop);
			cudaEventDestroy(start);
			cudaEventDestroy(stop);
			std::cerr << "It took "<< time <<" msecs."<< std::endl;


			cudaFree(d_Fimgs);
			cudaFree(d_Fimgs_nomask);

			cudaFree(d_weights);
			cudaFree(d_ctfs);
			cudaFree(d_Minvsigma2s);

			/*======================================================
								COLLECT DATA
			======================================================*/


			//TODO Following reduction should be done on the GPU
			double* wdiff2s_parts = new double[orientation_num*Frefs.xy];
			HANDLE_ERROR(cudaMemcpy( wdiff2s_parts, d_wdiff2s_parts,
					orientation_num*Frefs.xy*sizeof(double), cudaMemcpyDeviceToHost ));
			cudaFree(d_wdiff2s_parts);

			for (long int j = 0; j < Frefs.xy; j++)
			{
				int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, j);
				if (ires > -1)
				{
					double sum = 0;
					for (long int i = 0; i < orientation_num; i++)
						sum += wdiff2s_parts[i * Frefs.xy + j];

					exp_wsum_norm_correction[ipart] += sum;
					thr_wsum_sigma2_noise[group_id].data[ires] += sum;
				}
			}
			delete [] wdiff2s_parts;
			//TODO much in the following double loop can be GPU accelerated
			// exp_nr_dir * exp_nr_psi * exp_nr_oversampled_rot * exp_nr_trans * exp_nr_oversampled_trans
			for (long int i = 0; i < orientation_num; i++)
			{
				long int idir = idirs[i];
				long int iover_rot = iover_rots[i];
				double rot = rots[i];
				double tilt = tilts[i];
				double psi = psis[i];

				for (long int j = 0; j < translation_num; j++)
				{
					long int ihidden = iorientclasses[i] * exp_nr_trans + ihiddens[j];
					long int iover_trans = iover_transes[j];

					long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden, exp_current_oversampling,
																						iover_rot, iover_trans);

					double weight = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
					weight /= exp_sum_weight[ipart];

					// Store sum of weights for this group
					thr_sumw_group[group_id] += weight;
					// Store weights for this class and orientation
					thr_wsum_pdf_class[exp_iclass] += weight;

					double diffx = myprior_x - old_offset_x - oversampled_translations_x[iover_trans];
					double diffy = myprior_y - old_offset_y - oversampled_translations_y[iover_trans];

					if (mymodel.data_dim == 3)
					{
						double diffz  = myprior_z - old_offset_z - oversampled_translations_z[iover_trans];
					}

					// The following goes MUCH faster than the original lines below....
					if (mymodel.ref_dim == 2)
					{
						thr_wsum_prior_offsetx_class[exp_iclass] += weight * (old_offset_x + oversampled_translations_x[iover_trans]);
						thr_wsum_prior_offsety_class[exp_iclass] += weight * (old_offset_y + oversampled_translations_y[iover_trans]);
						thr_wsum_sigma2_offset += weight * (diffx*diffx + diffy*diffy);
					}
					else if (mymodel.data_dim == 3)
					{
						double diffz  = myprior_z - old_offset_z - oversampled_translations_z[iover_trans];
						thr_wsum_sigma2_offset += weight * (diffx*diffx + diffy*diffy + diffz*diffz);
					}

					// Store weight for this direction of this class
					if (mymodel.orientational_prior_mode == NOPRIOR)
					{
						DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], idir) += weight;
					}
					else
					{
						// In the case of orientational priors, get the original number of the direction back
						long int mydir = exp_pointer_dir_nonzeroprior[idir];
						DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], mydir) += weight;
					}

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
				}
			}

		} // end loop ipart

		cudaFree(d_Frefs);

		/*=======================================================================================
										   BACKPROJECTION
		=======================================================================================*/

		CudaComplex* wavgs = new CudaComplex[orientation_num * Frefs.xy];
		HANDLE_ERROR(cudaMemcpy( wavgs, d_wavgs, orientation_num * Frefs.xy*sizeof(CudaComplex), cudaMemcpyDeviceToHost ));
		cudaFree(d_wavgs);

		double* Fweights = new double[orientation_num * Frefs.xy];
		HANDLE_ERROR(cudaMemcpy( Fweights, d_Fweights, orientation_num * Frefs.xy *sizeof(double), cudaMemcpyDeviceToHost ));
		cudaFree(d_Fweights);

		Fweight.resize(exp_Fimgs[0]);
		double* FweightBag = Fweight.data; //TODO fix this

		Fimg.resize(exp_Fimgs[0]);
		Complex* FimgBag = Fimg.data; //TODO fix this

		for (long int i = 0; i < orientation_num; i++)
		{
			Euler_angles2matrix(rots[i], tilts[i], psis[i], A);

			Fimg.data = (Complex*) wavgs + i * Frefs.xy;
			Fweight.data = Fweights + i * Frefs.xy;


			std::ofstream myfile;
			std::string fnm = std::string("out_fweight_gpu.dat");
			myfile.open(fnm.c_str(), std::ios_base::app);
			for (int l = 0; l < Frefs.xy; l++)
				myfile << Fweights[i * Frefs.xy + l] << std::endl;
			myfile.close();
			exit(0);

//			Image<double> tt;
//			tt().resize(exp_current_image_size, exp_current_image_size);
//			FourierTransformer transformer;
//			transformer.inverseFourierTransform(Fimg, tt());
//			CenterFFT(tt(),false);
//			std::string fnm = std::string("out_wavg_gpu.mrc");
//			tt.write(fnm);
//
//			std::cerr << rots[i] << ":" << tilts[i] << ":" << psis[i] << std::endl;
//			exit(0);




			int my_mutex = exp_iclass % NR_CLASS_MUTEXES;
			pthread_mutex_lock(&global_mutex2[my_mutex]);
			(wsum_model.BPref[exp_iclass]).set2DFourierTransform(Fimg, A, IS_NOT_INV, &Fweight);
			pthread_mutex_unlock(&global_mutex2[my_mutex]);
		}

		delete [] wavgs;
		delete [] Fweights;

		Fweight.data = FweightBag; //TODO fix this
		Fimg.data = FimgBag; //TODO fix this
	} // end loop iclass


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
