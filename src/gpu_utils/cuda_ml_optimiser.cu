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

#define BLOCK_SIZE 32
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
									const unsigned img_size, const double sum_init)
{
	__shared__ double s[BLOCK_SIZE];
	s[threadIdx.x] = 0;

	unsigned pass_num(ceilf((float)img_size/(float)BLOCK_SIZE));
	unsigned long p, // image component index
		ref_start(blockIdx.x * img_size),
		img_start(blockIdx.y * img_size);

	for (unsigned pass = 0; pass < pass_num; pass ++)
	{
		p = pass * BLOCK_SIZE + threadIdx.x;

		if (p < img_size) //Is inside image
		{
			unsigned long ref_pixel_idx = ref_start + p;
			unsigned long img_pixel_idx = img_start + p;

		    double diff_real = g_refs[ref_pixel_idx].real - g_imgs[img_pixel_idx].real;
			double diff_imag = g_refs[ref_pixel_idx].imag - g_imgs[img_pixel_idx].imag;

			s[threadIdx.x] += (diff_real * diff_real + diff_imag * diff_imag) * 0.5 * g_Minvsigma2[p];
		}
	}

	// This version should run in             BLOCK_SIZE                  cycles
	// -------------------------------------------------------------------------
	//	if (threadIdx.x == 0)
	//	{
	//		double sum(sum_init);
	//		for (unsigned i = 0; i < BLOCK_SIZE; i++)
	//			sum += s[i];
	//
	//		g_diff2s[blockIdx.x * gridDim.y + blockIdx.y] = sum;
	//	}
	// -------------------------------------------------------------------------



	// This version should run in     BLOCK_SIZE/trads + log2(trads)      cycles
	// (  ~25x faster than above if memory conflicts are avoided )
	// -------------------------------------------------------------------------
	int trads = 32;
	int itr = BLOCK_SIZE/trads;
	if(itr>1)
	{
		for(int i=1; i<itr; i++)
		{
			if((i*trads+threadIdx.x)<BLOCK_SIZE)
			{
				s[threadIdx.x] += s[i*trads + threadIdx.x];
			}
			//__syncthreads();
		}
	}
	//__syncthreads();

	for(int j=(trads/2); j>0; j/=2)
	{
		if(threadIdx.x<j)
		{
			s[threadIdx.x]+=s[threadIdx.x+j];
		}
	}

	if (threadIdx.x*blockIdx.x == 0)
		g_diff2s[blockIdx.x * gridDim.y + blockIdx.y] = s[0]+sum_init;
	// -------------------------------------------------------------------------

}

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
			                            REFERENCE PROJECTION GENERATION
			=======================================================================================*/

			//printf("Generate Reference Projections\n");

			Fref.resize(exp_local_Minvsigma2s[0]);
			Complex* FrefBag = Fref.data; //TODO fix this

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

							Frefs.increment();

							orientation_num ++;
							iorientclasses.push_back(iorientclass);
							iover_rots.push_back(iover_rot);
						}
					}
				}
			}
			//printf("Finished generating reference projections\n");

			Fref.data = FrefBag; //TODO fix this

			CudaComplex *d_Frefs = Frefs.data_to_device();

			/*=======================================================================================
			                                  PARTICLE ITERATION
			=======================================================================================*/

			for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
			{

				/*======================================================
				                     TRANSLATIONS
				======================================================*/

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

				/*======================================================
						INITIATE PARTICLE SPECIFIC VALUES ON GPU
				======================================================*/

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
				//   meaning small-scale test will probably not catch this malfunctioning if it breaks.)
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
				HANDLE_ERROR(cudaMemcpy( d_Minvsigma2, Minvsigma2, Fimgs.xy * sizeof(double), cudaMemcpyHostToDevice));

				double *d_diff2s(0);
				HANDLE_ERROR(cudaMalloc( (void**) &d_diff2s, orientation_num*translation_num * sizeof(double)));
				//HANDLE_ERROR(cudaMemset(d_diff2s, exp_highres_Xi2_imgs[ipart] / 2., orientation_num*translation_num * sizeof(double))); //Initiate diff2 values with zeros

				/*======================================================
									KERNEL CALLS
				======================================================*/
				dim3 block_dim(orientation_num, translation_num);

				//printf("Calling kernel with <<(%d,%d), %d>> \n", block_dim.x, block_dim.y, BLOCK_SIZE);
				cuda_kernel_massive_diff2<<<block_dim,BLOCK_SIZE>>>(d_Frefs, d_Fimgs, d_Minvsigma2, d_diff2s, Frefs.xy, exp_highres_Xi2_imgs[ipart] / 2.);

//				for (long unsigned i = 0; i < orientation_num; i ++)
//				{
//					for (long unsigned j = 0; j < translation_num; j ++)
//					{
//						cuda_diff2_deviceImage( Frefs.xy, (double*) ( d_Frefs + (i * Frefs.xy) ), (double*) ( d_Fimgs + (j * Fimgs.xy) ), d_Minvsigma2, d_diff2s + (i * translation_num + j));
//					}
//				}

				/*======================================================
									RETRIEVE RESULTS
				======================================================*/

				HANDLE_ERROR(cudaDeviceSynchronize());
				//printf("Kernel call finished \n");

				double* diff2s = new double[orientation_num*translation_num];
				HANDLE_ERROR(cudaMemcpy( diff2s, d_diff2s, orientation_num*translation_num*sizeof(double), cudaMemcpyDeviceToHost ));


				/*======================================================
				                    COLLECT DATA
				======================================================*/



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


/*
Kernels Input
	Mresol_fine
	Frefs
	Fimgs
	exp_Mweight


Kernels Output
	thr_wsum_sigma2_noise (dynamic reduction, sad face)
	exp_wsum_norm_correction (sum of wdiff2)
*/
__global__ void cuda_kernel_massive_wdiff2(	CudaComplex *g_refs, CudaComplex *g_imgs, int *g_Mresol_fine, double* exp_Mweight,
									double *g_thr_wsum_sigma2_noise, double *g_exp_wsum_norm_correction,
									const unsigned img_size, const unsigned resol_count, const double exp_sum_weight,
									const unsigned exp_nr_oversampled_rot, const unsigned exp_nr_oversampled_trans)
{
	__shared__ double s_sum[BLOCK_SIZE];
	__shared__ double s_resol[BLOCK_SIZE * resol_count]; //Zero index

	long int ihidden_over = ihidden * exp_nr_oversampled_trans * exp_nr_oversampled_rot +
												iover_rot * exp_nr_oversampled_trans + iover_trans;

	s_sum[threadIdx.x] = 0;

	for (unsigned i = 0; i < resol_count; i ++)
		s_resol[resol_count * i + threadIdx.x ] = 0;

	unsigned pass_num(ceilf((float)img_size/(float)BLOCK_SIZE));
	unsigned long p,
		ref_start(blockIdx.x * img_size),
		img_start(blockIdx.y * img_size);

	for (unsigned pass = 0; pass < pass_num; pass ++)
	{
		p = pass * BLOCK_SIZE + threadIdx.x;

		if (p < img_size) //Is inside image
		{
			int ires = Mresol_fine[p];
			if (ires > -1)
			{
				unsigned long ref_pixel_idx = ref_start + p;
				unsigned long img_pixel_idx = img_start + p;

			    double diff_real = g_refs[ref_pixel_idx].real - g_imgs[img_pixel_idx].real;
				double diff_imag = g_refs[ref_pixel_idx].imag - g_imgs[img_pixel_idx].imag;

				s_sum[threadIdx.x] = weight * (diff_real*diff_real + diff_imag*diff_imag);
				s_resol[resol_count * ires + threadIdx.x] += wdiff2; //Putting neighboring ires adjacent
			}
		}
	}

	__syncthreads();

	if (threadIdx.x == 0)
	{
		double sum(0);
		for (unsigned i = 0; i < BLOCK_SIZE; i++)
			sum += s_sum[i];

		g_exp_wsum_norm_correction[blockIdx.x * gridDim.y + blockIdx.y] = sum;
	}

	if (threadIdx.x == 1)
	{
		for (unsigned ires = 0; ires < resol_count; ires ++)
		{
			double sum(0);
			for (unsigned i = 0; i < BLOCK_SIZE; i++)
				sum += s_resol[resol_count * ires + i];

			 //Putting neighboring ires adjacent
			g_thr_wsum_sigma2_noise[resol_count * ires + (blockIdx.x * gridDim.y + blockIdx.y)] = sum;
		}
	}
}

void MlOptimiserCUDA::storeWeightedSumsCUDA(long int my_ori_particle, int exp_current_image_size,
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




	/*
	This should be helpful:

	ihidden_over =
	iclass * exp_nr_dir * exp_nr_psi * exp_nr_trans * exp_nr_oversampled_rot * exp_nr_oversampled_trans +
	idir *                exp_nr_psi * exp_nr_trans * exp_nr_oversampled_rot * exp_nr_oversampled_trans +
	ipsi *                             exp_nr_trans * exp_nr_oversampled_rot * exp_nr_oversampled_trans +
	itrans *                                          exp_nr_oversampled_rot * exp_nr_oversampled_trans +
	iover_rot *                                                                exp_nr_oversampled_trans +
	iover_trans
	*/

	unsigned long class_size = exp_nr_dir * exp_nr_psi * exp_nr_trans * exp_nr_oversampled_rot * exp_nr_oversampled_trans;
	unsigned long dir_size =                exp_nr_psi * exp_nr_trans * exp_nr_oversampled_rot * exp_nr_oversampled_trans;
	unsigned long psi_size =                             exp_nr_trans * exp_nr_oversampled_rot * exp_nr_oversampled_trans;
	unsigned long trans_size =                                          exp_nr_oversampled_rot * exp_nr_oversampled_trans;
	unsigned long over_rot_size =                                                                exp_nr_oversampled_trans;


	// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++)
	{



		/*=======================================================================================
		                            REFERENCE PROJECTION GENERATION
		=======================================================================================*/


		CudaImages Frefs(exp_local_Minvsigma2s[0].xdim, exp_local_Minvsigma2s[0].ydim,
				(exp_idir_max - exp_idir_min + 1) * (exp_ipsi_max - exp_ipsi_min + 1) * exp_nr_oversampled_rot); //TODO set to exp_nr_something

		Complex* FrefBag = Fref.data; //TODO fix this

		for (long int idir = exp_idir_min; idir <= exp_idir_max; idir++)
		{
			for (long int ipsi = exp_ipsi_min; ipsi <= exp_ipsi_max; ipsi++)
			{
				sampling.getOrientations(idir, ipsi, adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
						exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior);

				for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
				{
					rot = oversampled_rot[iover_rot];
					tilt = oversampled_tilt[iover_rot];
					psi = oversampled_psi[iover_rot];
					// Get the Euler matrix
					Euler_angles2matrix(rot, tilt, psi, A);

					Fref.data = (Complex*) Frefs.current();
					(mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_NOT_INV);

					Frefs.increment();
				}
			}
		}
		Fref.data = FrefBag; //TODO fix this

		CudaComplex *d_Frefs = Frefs.data_to_device();

		/*=======================================================================================
										  PARTICLE ITERATION
		=======================================================================================*/

		// Inside the loop over all translations and all part_id sum all shift Fimg's and their weights
		// Then outside this loop do the actual backprojection
		Fimg.initZeros();
		Fweight.initZeros();
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

			CudaImages Fimgs(Frefs.x, Frefs.y,
					( exp_itrans_max - exp_itrans_min + 1) * exp_nr_oversampled_trans);//TODO set to exp_nr_something
			CudaImages Fimgs_nomask(Frefs.x, Frefs.y,
					( exp_itrans_max - exp_itrans_min + 1) * exp_nr_oversampled_trans);//TODO set to exp_nr_something

			for (long int itrans = exp_itrans_min, iitrans = 0; itrans <= exp_itrans_max; itrans++)
			{
				sampling.getTranslations(itrans, adaptive_oversampling,
						oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);
				for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++, iitrans++)
				{
					/// Now get the shifted image
					// Use a pointer to avoid copying the entire array again in this highly expensive loop
					Complex *Fimg_shift, *Fimg_shift_nomask;
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
					Fimg_shift = Fimg_otfshift.data;
					Fimg_shift_nomask = Fimg_otfshift_nomask.data;

					Fimgs.increment();
					Fimgs_nomask.increment();
				}
			}

			/*======================================================
					INITIATE PARTICLE SPECIFIC VALUES ON GPU
			======================================================*/

			Minvsigma2 = exp_local_Minvsigma2s[ipart].data;

			double *d_Minvsigma2(0);

			CudaComplex *d_Fimgs = Fimgs.data_to_device();

			HANDLE_ERROR(cudaMalloc( (void**) &d_Minvsigma2, Fimgs.xy * sizeof(double)));
			HANDLE_ERROR(cudaMemcpy( d_Minvsigma2, Minvsigma2, Fimgs.xy * sizeof(double), cudaMemcpyHostToDevice));

			double *d_diff2s(0);
			HANDLE_ERROR(cudaMalloc( (void**) &d_diff2s, orientation_num*translation_num * sizeof(double)));
			//HANDLE_ERROR(cudaMemset(d_diff2s, exp_highres_Xi2_imgs[ipart] / 2., orientation_num*translation_num * sizeof(double))); //Initiate diff2 values with zeros

			/*======================================================
								KERNEL CALLS
			======================================================*/

			dim3 block_dim(orientation_num, translation_num);

			//printf("Calling kernel with <<(%d,%d), %d>> \n", block_dim.x, block_dim.y, BLOCK_SIZE);
			cuda_kernel_massive_diff2<<<block_dim,BLOCK_SIZE>>>();


			/*
			Kernels Input
				Mresol_fine
				Frefs
				Fimgs
				exp_Mweight


			Kernels Output
				thr_wsum_sigma2_noise
				exp_wsum_norm_correction (sum of wdiff2)
			*/

			double weight = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
			weight /= exp_sum_weight[ipart];

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
			{
				int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
				if (ires > -1)
				{
					double diff_real = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real - (*(Fimg_shift + n)).real;
					double diff_imag = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag - (*(Fimg_shift + n)).imag;
					double wdiff2 = weight * (diff_real*diff_real + diff_imag*diff_imag);

					DIRECT_MULTIDIM_ELEM(thr_wsum_sigma2_noise[group_id], ires) += wdiff2;

					exp_wsum_norm_correction[ipart] += wdiff2;
				}
			}

			/*
			Kernels Input
				Fimgs_nomask
				exp_Mweight
				Mctf
				Minvsigma2

			Kernels Output
				wFimg
				Fweight
			*/

			// Store sum of weight*SSNR*Fimg in data and sum of weight*SSNR in weight
			// Use the FT of the unmasked image to back-project in order to prevent reconstruction artefacts! SS 25oct11
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wFimg)
			{
				double myctf = DIRECT_MULTIDIM_ELEM(Mctf, n);
				// Note that weightxinvsigma2 already contains the CTF!
				double weightxinvsigma2 = weight * myctf * DIRECT_MULTIDIM_ELEM(Minvsigma2, n);
				// now Fimg stores sum of all shifted w*Fimg
				(DIRECT_MULTIDIM_ELEM(wFimg, n)).real += (*(Fimg_shift_nomask + n)).real * weightxinvsigma2;
				(DIRECT_MULTIDIM_ELEM(wFimg, n)).imag += (*(Fimg_shift_nomask + n)).imag * weightxinvsigma2;
				// now Fweight stores sum of all w
				// Note that CTF needs to be squared in Fweight, weightxinvsigma2 already contained one copy
				DIRECT_MULTIDIM_ELEM(Fweight, n) += weightxinvsigma2 * myctf;
			}






			/*======================================================
								COLLECT DATA
			======================================================*/

			for (long int i = 0; i < orientation_num; i++)
			{
				long int iover_rot = iover_rots[i];

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
				}
			}

			double diffx = myprior_x - old_offset_x - oversampled_translations_x[iover_trans];
			double diffy = myprior_y - old_offset_y - oversampled_translations_y[iover_trans];

			if (mymodel.data_dim == 3)
			{
				double diffz  = myprior_z - old_offset_z - oversampled_translations_z[iover_trans];
			}



			cudaFree(d_Fimgs);
			cudaFree(d_diff2s);
			delete [] diff2s;

		} // end loop ipart

		cudaFree(d_Frefs);
	} // end loop iclass
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
