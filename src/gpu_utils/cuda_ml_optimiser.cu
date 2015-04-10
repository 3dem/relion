#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include "src/gpu_utils/cuda_ml_optimiser.h"
#include "src/gpu_utils/cuda_img_operations.h"
#include "src/gpu_utils/cuda_utils.cuh"
#include "src/complex.h"
#include <fstream>
#include <cuda.h>

#define MAX_RESOL_SHARED_MEM 32
#define BLOCK_SIZE 128         // This is optimally set as big as possible without its ceil:ed multiple exceeding imagesize by too much.
#define NR_CLASS_MUTEXES 5

static pthread_mutex_t global_mutex2[NR_CLASS_MUTEXES] = { PTHREAD_MUTEX_INITIALIZER };
static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

__global__ void cuda_kernel_diff2(	CudaComplex *g_refs, CudaComplex *g_imgs,
									FLOAT *g_Minvsigma2, FLOAT *g_diff2s,
									unsigned img_size, FLOAT sum_init,
									unsigned long significant_num,
									unsigned long translation_num,
									unsigned long *d_rotidx,
									unsigned long *d_transidx)
{
	// blockid
	int ex = blockIdx.y * gridDim.x + blockIdx.x;

	// inside the padded 2D orientation grid
	if( ex < significant_num )
	{
		// index of comparison
		unsigned long int ix=d_rotidx[ex];
		unsigned long int iy=d_transidx[ex];

		__shared__ FLOAT s[BLOCK_SIZE];
		s[threadIdx.x] = 0;

		unsigned pass_num(ceilf((float)img_size/(float)BLOCK_SIZE)), pixel;

		unsigned long ref_start(ix * img_size);
		unsigned long img_start(iy * img_size);
		unsigned long ref_pixel_idx;
		unsigned long img_pixel_idx;

		for (unsigned pass = 0; pass < pass_num; pass ++)
		{
			pixel = pass * BLOCK_SIZE + threadIdx.x;

			if (pixel < img_size) //Is inside image
			{
				ref_pixel_idx = ref_start + pixel;
				img_pixel_idx = img_start + pixel;

				FLOAT diff_real = g_refs[ref_pixel_idx].real - g_imgs[img_pixel_idx].real;
				FLOAT diff_imag = g_refs[ref_pixel_idx].imag - g_imgs[img_pixel_idx].imag;

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
			g_diff2s[ix * translation_num + iy] = s[0]+sum_init;
		}
		// -------------------------------------------------------------------------
	}
}

__global__ void cuda_kernel_cc_diff2(	CudaComplex *g_refs, CudaComplex *g_imgs,
										FLOAT *g_Minvsigma2, FLOAT *g_diff2s,
										unsigned img_size, FLOAT exp_local_sqrtXi2,
										unsigned long significant_num,
										unsigned long translation_num,
										unsigned long *d_rotidx,
										unsigned long *d_transidx)
{
	// blockid
	int ex = blockIdx.y * gridDim.x + blockIdx.x;
	// inside the padded 2D orientation grid
	if( ex < significant_num )
	{
		// index of comparison
		unsigned long int ix=d_rotidx[ex];
		unsigned long int iy=d_transidx[ex];
		__shared__ double    s[BLOCK_SIZE];
		__shared__ double norm[BLOCK_SIZE];
		s[threadIdx.x] = 0;
		unsigned pass_num(ceilf((float)img_size/(float)BLOCK_SIZE));
		unsigned long pixel,
		ref_start(ix * img_size),
		img_start(iy * img_size);
		unsigned long ref_pixel_idx;
		unsigned long img_pixel_idx;
		for (unsigned pass = 0; pass < pass_num; pass ++)
		{
			pixel = pass * BLOCK_SIZE + threadIdx.x;

			if (pixel < img_size) //Is inside image
			{
				ref_pixel_idx = ref_start + pixel;
				img_pixel_idx = img_start + pixel;

				double diff_real = g_refs[ref_pixel_idx].real * g_imgs[img_pixel_idx].real;
				double diff_imag = g_refs[ref_pixel_idx].imag * g_imgs[img_pixel_idx].imag;

				double nR = g_refs[ref_pixel_idx].real*g_refs[ref_pixel_idx].real;
				double nI = g_refs[ref_pixel_idx].imag*g_refs[ref_pixel_idx].imag;

				s[threadIdx.x] -= (diff_real + diff_imag);
				norm[threadIdx.x] += nR+nI;
			}
		}
		// -------------------------------------------------------------------------
		__syncthreads();
		int trads = 32;
		int itr = BLOCK_SIZE/trads;
		if(threadIdx.x<trads)
		{
			for(int i=1; i<itr; i++)
			{
				s[threadIdx.x] += s[i*trads + threadIdx.x];
				norm[threadIdx.x] += norm[i*trads + threadIdx.x];
			}
		}
		for(int j=(trads/2); j>0; j/=2)
		{
			if(threadIdx.x<j)
			{
				s[threadIdx.x] += s[threadIdx.x+j];
				norm[threadIdx.x] += norm[threadIdx.x+j];
			}
		}
		__syncthreads();
		// -------------------------------------------------------------------------
		g_diff2s[ix * translation_num + iy] = s[0]/(sqrt(norm[0])*exp_local_sqrtXi2);
	}
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

	CUDA_CPU_TIC("diff_pre_gpu");

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

	MultidimArray<Complex > Fref;
	Fref.resize(exp_local_Minvsigma2s[0]);

	unsigned image_size = exp_local_Minvsigma2s[0].nzyxdim;

	CUDA_CPU_TOC("diff_pre_gpu");

	// Loop only from exp_iclass_min to exp_iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++)
	{
		if (mymodel.pdf_class[exp_iclass] > 0.)
		{
			// Local variables
			std::vector< double > oversampled_rot, oversampled_tilt, oversampled_psi;
			std::vector< double > oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
			CudaGlobalPtr<FLOAT> gpuMinvsigma2(image_size);
			gpuMinvsigma2.device_alloc();

			Matrix2D<double> A;

			CudaGlobalPtr<CudaComplex> Frefs(image_size * exp_nr_dir * exp_nr_psi * exp_nr_oversampled_rot);

			// Mapping index look-up table
			std::vector< long unsigned > iorientclasses, iover_rots;
			long unsigned orientation_num(0);

			/*=======================================================================================
			                           Generate Reference Projections
			=========================================================================================*/

			CUDA_CPU_TIC("projection_1");

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

							//Fref.data = &Frefs[image_size * orientation_num];
							(mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_NOT_INV);

							for (unsigned i = 0; i < image_size; i++)
							{
								Frefs[image_size * orientation_num + i].real = Fref.data[i].real;
								Frefs[image_size * orientation_num + i].imag = Fref.data[i].imag;
							}

							orientation_num ++;
							iorientclasses.push_back(iorientclass);
							iover_rots.push_back(iover_rot);
						}
					}
				}
			}

			Frefs.size = orientation_num * image_size;
			Frefs.device_alloc();
			Frefs.cp_to_device();
			Frefs.free_host();

			CUDA_CPU_TOC("projection_1");

			/*=======================================================================================
			                                  	  Particle Iteration
			=========================================================================================*/

			for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
			{
				/*====================================
				        Generate Translations
				======================================*/

				CUDA_CPU_TIC("translation_1");

				CudaGlobalPtr<CudaComplex> Fimgs(image_size * exp_nr_trans * exp_nr_oversampled_trans);

				long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
				long unsigned translation_num(0), ihidden(0);
				std::vector< long unsigned > iover_transes, itranses, ihiddens;

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
							myAB = (Fref.ydim == coarse_size) ? global_fftshifts_ab_coarse[itrans].data
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
							FLOAT real = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real
									- (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag;
							FLOAT imag = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag
									+ (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real;

							//When on gpu, it makes more sense to ctf-correct translated images, rather than anti-ctf-correct ref-projections
							if (do_scale_correction)
							{
								//group_id = mydata.getGroupId(part_id);
								FLOAT myscale = mymodel.scale_correction[group_id];
								real /= myscale;
								imag /= myscale;
							}
							if (do_ctf_correction && refs_are_ctf_corrected)
							{
								real /= DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], n);
								imag /= DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], n);
							}
							Fimgs[translation_num * image_size + n].real = real;
							Fimgs[translation_num * image_size + n].imag = imag;
						}
						translation_num ++;

						ihiddens.push_back(ihidden);
						itranses.push_back(itrans);
						iover_transes.push_back(iover_trans);
					}
				}

				CUDA_CPU_TOC("translation_1");

				/*===========================================
				   Determine significant comparison indices
				=============================================*/
				//      This section is annoying to test because
				//		it can't complete on first pass, since
				//		the significance has never been set


				CUDA_CPU_TIC("pair_list_1");

				CudaGlobalPtr<long unsigned> transidx(orientation_num*translation_num), rotidx(orientation_num*translation_num);

				long unsigned coarse_num = exp_nr_dir*exp_nr_psi*exp_nr_trans;
				long unsigned significant_num(0);

				if (exp_ipass == 0)
				{
					exp_Mcoarse_significant.resize(coarse_num, 1);
					for (long unsigned i = 0; i < orientation_num; i++)
					{
						for (long unsigned j = 0; j < translation_num; j++)
						{
							rotidx[significant_num] = i;
							transidx[significant_num] = j;
							significant_num++;
//							DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, i)=1;
//							std::cerr << "exp_Mcoarse_significant("<< i <<") = " <<    DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, i) << std::endl;
//							std::cerr << "exp_Mcoarse_significant("<< i <<") = " << *(&DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, 0)+i*sizeof(bool)) << std::endl;
						}
					}
				}
				else
				{
					for (long unsigned i = 0; i < orientation_num; i++)
					{
						long int iover_rot = iover_rots[i];
//						long int iover_rot = i % exp_nr_oversampled_rot
						long int coarse_rot = floor(i/exp_nr_oversampled_rot);
						for (long unsigned j = 0; j < translation_num; j++)
						{
							long int iover_trans = iover_transes[j];
//							long int iover_trans = j % exp_nr_oversampled_trans
							long int coarse_trans = floor(j/exp_nr_oversampled_trans);
							long int ihidden = iorientclasses[i] * exp_nr_trans + ihiddens[j];
							if(DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden)==1)
							{
								 long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden,
										                  exp_current_oversampling, iover_rot, iover_trans);

								rotidx[significant_num] = i;
								transidx[significant_num] = j;
								significant_num++;
							}
						}
					}
				}

				CUDA_CPU_TOC("pair_list_1");

//				std::cerr << "orientation_num "<< orientation_num << std::endl;
//				std::cerr << "translation_num "<< translation_num << std::endl;
//				std::cerr << "my_nr_significant_coarse_samples "<< DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_NR_SIGN) << std::endl;
//				std::cerr << "significant_num "<< significant_num << std::endl;

				/*====================================
				   Initiate Particle Related On GPU
				======================================*/


				CUDA_CPU_TIC("kernel_init_1");

				// Since we hijack Minvsigma to carry a bit more info into the GPU-kernel
				// we need to make a modified copy, since the global object shouldn't be
				// changed


				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Fimgs_shifted[ipart])
				{
					gpuMinvsigma2[n] = *(exp_local_Minvsigma2s[ipart].data + n );
					//std::cerr <<  *(exp_local_Minvsigma2s[ipart].data + n )<< " ";
				}

				if (do_ctf_correction && refs_are_ctf_corrected)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Fimgs_shifted[ipart])
					{
						gpuMinvsigma2[n] *= (DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], n)*DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], n));
					}
				}
				// TODO :    + Assure accuracy with the implemented GPU-based ctf-scaling
				//           + Make setting of myscale robust between here and above.
				//  (scale_correction turns off by default with only one group: ml_optimiser-line 1067,
				//   meaning small-scale test will probably not catch this malfunctioning when/if it breaks.)
				if (do_scale_correction)
				{
					FLOAT myscale = mymodel.scale_correction[group_id];
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Fimgs_shifted[ipart])
					{
						gpuMinvsigma2[n] *= (myscale*myscale);
					}
				}

				Fimgs.size = translation_num * image_size;
				Fimgs.device_alloc();
				Fimgs.cp_to_device();

				gpuMinvsigma2.cp_to_device();

				CudaGlobalPtr<FLOAT> diff2s(orientation_num*translation_num);
				diff2s.device_alloc();

				rotidx.size = significant_num;
				rotidx.device_alloc();
				rotidx.cp_to_device();

				transidx.size = significant_num;
				transidx.device_alloc();
				transidx.cp_to_device();

				/*====================================
				    		Kernel Calls
				======================================*/
				unsigned orient1, orient2;

				if(significant_num>65535)
				{
					orient1 = ceil(sqrt(significant_num));
					orient2 = orient1;
				}
				else
				{
					orient1 = significant_num;
					orient2 = 1;
				}
				dim3 block_dim(orient1,orient2);

				CUDA_CPU_TOC("kernel_init_1");

				CUDA_GPU_TIC("cuda_kernel_diff2");

				if ((iter == 1 && do_firstiter_cc) || do_always_cc) // do cross-correlation instead of diff
				{
					cuda_kernel_cc_diff2<<<block_dim,BLOCK_SIZE>>>(~Frefs, ~Fimgs, ~gpuMinvsigma2,  ~diff2s,
																	image_size, exp_highres_Xi2_imgs[ipart],
																	significant_num,
																	translation_num,
																	~rotidx,
																	~transidx);
				}
				else
				{
					cuda_kernel_diff2<<<block_dim,BLOCK_SIZE>>>(~Frefs, ~Fimgs, ~gpuMinvsigma2, ~diff2s,
																image_size, exp_highres_Xi2_imgs[ipart] / 2.,
																significant_num,
																translation_num,
																~rotidx,
																~transidx);
				}

				CUDA_GPU_TAC("cuda_kernel_diff2");

				/*====================================
				    	   Retrieve Results
				======================================*/

				HANDLE_ERROR(cudaDeviceSynchronize()); //TODO Apparently this is not required here

				CUDA_GPU_TOC("cuda_kernel_diff2");

				diff2s.cp_to_host();

				if (exp_ipass == 0)
				{
					exp_Mcoarse_significant.clear();
				}

				/*====================================
				    	Write To Destination
				======================================*/


				CUDA_CPU_TIC("collect_data_1");

				for (long unsigned k = 0; k < significant_num; k ++)
				{
					long unsigned i = rotidx[k];
					long unsigned j = transidx[k];
					long int iover_rot = iover_rots[i];

					long int ihidden = iorientclasses[i] * exp_nr_trans + ihiddens[j];
					long int iover_trans = iover_transes[j];

					long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden, exp_current_oversampling,
																						iover_rot, iover_trans);

					double diff2 = diff2s[i * translation_num + j];

					DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) = diff2;

					// Keep track of minimum of all diff2, only for the last image in this series
					if (diff2 < exp_min_diff2[ipart])
						exp_min_diff2[ipart] = diff2;
				}

				CUDA_CPU_TOC("collect_data_1");

			} // end loop ipart

		} // end if class significant
	} // end loop iclass
}

__global__ void cuda_kernel_wavg(	CudaComplex *g_refs, CudaComplex *g_imgs, CudaComplex *g_imgs_nomask,
									FLOAT* g_weights, FLOAT* g_ctfs, FLOAT* g_Minvsigma2s,
									FLOAT *g_wdiff2s_parts, CudaComplex *g_wavgs, FLOAT* g_Fweights,
									unsigned long translation_num, FLOAT weight_norm,
									FLOAT significant_weight, unsigned image_size,
									bool refs_are_ctf_corrected)
{
	unsigned long iorient = blockIdx.y*gridDim.x + blockIdx.x;
	unsigned tid = threadIdx.x;

	//TODO Consider Mresol_fine to speed up this kernel

	unsigned pass_num(ceilf((float)image_size/(float)BLOCK_SIZE)),pixel;
	FLOAT Fweight, wavgs_real, wavgs_imag, wdiff2s_parts;

	for (unsigned pass = 0; pass < pass_num; pass ++)
	{
		wavgs_real = 0;
		wavgs_imag = 0;
		wdiff2s_parts = 0;
		Fweight = 0;

		pixel = pass * BLOCK_SIZE + tid;

		if (pixel < image_size)
		{
			unsigned long orientation_pixel = iorient * image_size + pixel;

			for (unsigned long itrans = 0; itrans < translation_num; itrans++)
			{
				FLOAT weight = g_weights[iorient * translation_num + itrans];

				if (weight >= significant_weight)
				{
					weight /= weight_norm;

					unsigned long img_pixel_idx = itrans * image_size + pixel;

					FLOAT myctf = g_ctfs[pixel];
					CudaComplex ref = g_refs[orientation_pixel];
					if (refs_are_ctf_corrected) //FIXME Create two kernels for the different cases
					{
						ref.real *= myctf;
						ref.imag *= myctf;
					}
					FLOAT diff_real = ref.real - g_imgs[img_pixel_idx].real;
					FLOAT diff_imag = ref.imag - g_imgs[img_pixel_idx].imag;
					wdiff2s_parts += weight * (diff_real*diff_real + diff_imag*diff_imag);
					FLOAT weightxinvsigma2 = weight * myctf * g_Minvsigma2s[pixel];
					wavgs_real += g_imgs_nomask[img_pixel_idx].real * weightxinvsigma2;
					wavgs_imag += g_imgs_nomask[img_pixel_idx].imag * weightxinvsigma2;
					Fweight += weightxinvsigma2 * myctf;
				}
			}

			g_wavgs[orientation_pixel].real = wavgs_real; //TODO should be buffered into shared
			g_wavgs[orientation_pixel].imag = wavgs_imag; //TODO should be buffered into shared
			g_wdiff2s_parts[orientation_pixel] = wdiff2s_parts; //TODO this could be further reduced in here
			g_Fweights[orientation_pixel] = Fweight; //TODO should be buffered into shared
		}
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
	CUDA_CPU_TIC("store_pre_gpu");

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
	MultidimArray<Complex > Fimg, Fimg_otfshift_nomask;
	MultidimArray<double> Fweight, Minvsigma2, Mctf;
	bool have_warned_small_scale = false;

	Fimg.resize(exp_Fimgs[0]);
	Fweight.resize(exp_Fimgs[0]);

	// Initialise Mctf to all-1 for if !do_ctf_corection
	Mctf.resize(exp_Fimgs[0]);
	Mctf.initConstant(1.);
	// Initialise Minvsigma2 to all-1 for if !do_map
	Minvsigma2.resize(exp_Fimgs[0]);
	Minvsigma2.initConstant(1.);

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

	unsigned image_size = exp_local_Minvsigma2s[0].xdim*exp_local_Minvsigma2s[0].ydim;

	CUDA_CPU_TOC("store_pre_gpu");

	// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++)
	{

		/*=======================================================================================
		                            REFERENCE PROJECTION GENERATION
		=======================================================================================*/
		CUDA_CPU_TIC("projection_2");

		CudaGlobalPtr<CudaComplex> Frefs(image_size * exp_nr_dir * exp_nr_psi * exp_nr_oversampled_rot);

		std::vector< long unsigned > iorientclasses, idirs, iover_rots;
		std::vector< double > rots, tilts, psis;
		long unsigned orientation_num(0);

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

						//Fref.data = &Frefs[image_size * orientation_num];
						(mymodel.PPref[exp_iclass]).get2DFourierTransform(Fimg, A, IS_NOT_INV);

						for (unsigned i = 0; i < image_size; i++)
						{
							Frefs[image_size * orientation_num + i].real = Fimg.data[i].real;
							Frefs[image_size * orientation_num + i].imag = Fimg.data[i].imag;
						}

						orientation_num ++;
						idirs.push_back(idir);
						iorientclasses.push_back(iorientclass);
						iover_rots.push_back(iover_rot);
					}
				}
			}
		}

		Frefs.size = orientation_num * image_size;
		Frefs.device_alloc();
		Frefs.cp_to_device();
		Frefs.free_host();

		CudaGlobalPtr<CudaComplex> wavgs(orientation_num * image_size);
		wavgs.device_alloc();
		//wavgs.device_init(0);

		CudaGlobalPtr<FLOAT> Fweights(orientation_num * image_size);
		Fweights.device_alloc();
		//Fweights.device_init(0);

		CUDA_CPU_TOC("projection_2");

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

			CUDA_CPU_TIC("translation_2");

			CudaGlobalPtr<CudaComplex> Fimgs(image_size * exp_nr_trans * exp_nr_oversampled_trans);
			CudaGlobalPtr<CudaComplex> Fimgs_nomask(Fimgs.size);

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
						FLOAT a = (*(myAB + n)).real;
						FLOAT b = (*(myAB + n)).imag;

						// Fimg_shift
						FLOAT real = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real
								- b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag;
						FLOAT imag = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag
								+ b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real;
						Fimgs[translation_num * image_size + n].real = real;
						Fimgs[translation_num * image_size + n].imag = imag;

						// Fimg_shift_nomask
						real = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).real
								- b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).imag;
						imag = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).imag
								+ b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).real;
						Fimgs_nomask[translation_num * image_size + n].real = real;
						Fimgs_nomask[translation_num * image_size + n].imag = imag;
					}

					translation_num ++;

					ihiddens.push_back(ihidden);
					itranses.push_back(itrans);
					iover_transes.push_back(iover_trans);
				}
			}

			Fimgs.size = translation_num * image_size;
			Fimgs.device_alloc();
			Fimgs.cp_to_device();

			Fimgs_nomask.size = translation_num * image_size;
			Fimgs_nomask.device_alloc();
			Fimgs_nomask.cp_to_device();

			CUDA_CPU_TOC("translation_2");

			/*======================================================
					            MAP WEIGHTS
			======================================================*/

			CudaGlobalPtr<FLOAT> sorted_weights(orientation_num * translation_num);

			for (long unsigned i = 0; i < orientation_num; i++)
			{
				long unsigned iover_rot = iover_rots[i];
				for (long unsigned j = 0; j < translation_num; j++)
				{
					long unsigned iover_trans = iover_transes[j];
					long unsigned ihidden = iorientclasses[i] * exp_nr_trans + ihiddens[j];
					long unsigned ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden,
											  exp_current_oversampling, iover_rot, iover_trans);
					sorted_weights[(long unsigned) i * translation_num + j] =
							DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
				}
			}


			/*======================================================
					            KERNEL CALL
			======================================================*/
#ifdef DEBUG_CUDA_MEM
		printf("Before Cpy to Device: ");
		cudaPrintMemInfo();
#endif

			sorted_weights.device_alloc();
			sorted_weights.cp_to_device();
			sorted_weights.free_host();

			CudaGlobalPtr<FLOAT> ctfs(image_size); //TODO Almost same size for all iparts, should be allocated once
			ctfs.device_alloc();

			if (do_ctf_correction)
			{
				for (unsigned i = 0; i < image_size; i++)
					ctfs[i] = (FLOAT) exp_local_Fctfs[ipart].data[i];
				ctfs.cp_to_device();
			}
			else
				ctfs.device_init(1.);

			CudaGlobalPtr<FLOAT> Minvsigma2s(image_size); //TODO Almost same size for all iparts, should be allocated once
			for (unsigned i = 0; i < image_size; i++)
				Minvsigma2s[i] = exp_local_Minvsigma2s[ipart].data[i];

			Minvsigma2s.device_alloc();
			Minvsigma2s.cp_to_device();

			CudaGlobalPtr<FLOAT> wdiff2s_parts(orientation_num * image_size); //TODO Almost same size for all iparts, should be allocated once
			wdiff2s_parts.device_alloc();

			unsigned orient1, orient2;
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
			dim3 block_dim(orient1,orient2);

			CUDA_GPU_TIC("cuda_kernel_wavg");

			cuda_kernel_wavg<<<block_dim,BLOCK_SIZE>>>(
												~Frefs, ~Fimgs, ~Fimgs_nomask,
												~sorted_weights, ~ctfs, ~Minvsigma2s,
												~wdiff2s_parts, ~wavgs, ~Fweights,
												translation_num,
												(FLOAT) exp_sum_weight[ipart],
												(FLOAT) exp_significant_weight[ipart],
												image_size,
												refs_are_ctf_corrected
												);

			CUDA_GPU_TAC("cuda_kernel_wavg");

			HANDLE_ERROR(cudaDeviceSynchronize()); //TODO Apparently this is not required here

			CUDA_GPU_TOC("cuda_kernel_wavg");

			Fimgs.free_device();
			Fimgs_nomask.free_device();

			sorted_weights.free_device();
			ctfs.free_device();
			Minvsigma2s.free_device();

			/*======================================================
								COLLECT DATA
			======================================================*/

			CUDA_CPU_TIC("reduce_wdiff2s");

			//TODO Following reduction should be done on the GPU
			wdiff2s_parts.cp_to_host();
			wdiff2s_parts.free_device();

			for (long int j = 0; j < image_size; j++)
			{
				int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, j);
				if (ires > -1)
				{
					double sum = 0;
					for (long int i = 0; i < orientation_num; i++)
						sum += (double) wdiff2s_parts[i * image_size + j];
					thr_wsum_sigma2_noise[group_id].data[ires] += sum;
					exp_wsum_norm_correction[ipart] += sum;
				}
			}

			wdiff2s_parts.free_host();

			CUDA_CPU_TOC("reduce_wdiff2s");

			CUDA_CPU_TIC("collect_data_2");
			//TODO some in the following double loop can be GPU accelerated
			//TODO should be replaced with loop over pairs of projections and translations (like in the getAllSquaredDifferences-function)

			// exp_nr_dir * exp_nr_psi * exp_nr_oversampled_rot * exp_nr_trans * exp_nr_oversampled_trans
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
								double rot = oversampled_rot[iover_rot];
								double tilt = oversampled_tilt[iover_rot];
								double psi = oversampled_psi[iover_rot];

								// Get the Euler matrix
								Euler_angles2matrix(rot, tilt, psi, A);


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
										if (weight >= exp_significant_weight[ipart])
										{
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

											if (weight > exp_max_weight[ipart])
											{
												// Store optimal image parameters
												exp_max_weight[ipart] = weight;

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
								}
							}
						}
					}
				}
			}

			CUDA_CPU_TOC("collect_data_2");

#ifdef DEBUG_CUDA_MEM
		printf("After Freeing Device Mem: ");
		cudaPrintMemInfo();
#endif

		} // end loop ipart

		Frefs.free_device();

		/*=======================================================================================
										   BACKPROJECTION
		=======================================================================================*/

		CUDA_CPU_TIC("backprojection");

		wavgs.cp_to_host();
		wavgs.free_device();

		Fweights.cp_to_host();
		Fweights.free_device();

#ifdef RELION_TESTING
		std::string fnm = std::string("gpu_out_exp_wsum_norm_correction.txt");
		char *text = &fnm[0];
		freopen(text,"w",stdout);
		for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
		{
			printf("%4.8f \n",exp_wsum_norm_correction[ipart]);
		}
		fclose(stdout);
		//----------
		fnm = std::string("gpu_out_thr_wsum_sigma2_noise.txt");
		text = &fnm[0];
		freopen(text,"w",stdout);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
		{
			printf("%4.8f \n",thr_wsum_sigma2_noise[0].data[n]);
		}
		fclose(stdout);
		//----------
		fnm = std::string("gpu_out_Fweights.txt");
		text = &fnm[0];
		freopen(text,"w",stdout);
		for(int n = 0; n < 1000; n++)
		{
			printf("%4.8f \n",Fweights[n*60+50]);
		}
		fclose(stdout);
#endif

		for (long int i = 0; i < orientation_num; i++)
		{
			Euler_angles2matrix(rots[i], tilts[i], psis[i], A);

			for (unsigned j = 0; j < image_size; j++)
			{
				Fimg.data[j].real = wavgs[i * image_size + j].real;
				Fimg.data[j].imag = wavgs[i * image_size + j].imag;
				Fweight.data[j] = Fweights[i * image_size + j];
			}

			int my_mutex = exp_iclass % NR_CLASS_MUTEXES;
			pthread_mutex_lock(&global_mutex2[my_mutex]);
			(wsum_model.BPref[exp_iclass]).set2DFourierTransform(Fimg, A, IS_NOT_INV, &Fweight);
			pthread_mutex_unlock(&global_mutex2[my_mutex]);
		}

		CUDA_CPU_TOC("backprojection");

	} // end loop iclass

	CUDA_CPU_TIC("store_post_gpu");

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

	CUDA_CPU_TOC("store_post_gpu");
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
