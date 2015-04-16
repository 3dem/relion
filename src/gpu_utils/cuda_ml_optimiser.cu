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
#include "src/parallel.h"

#define MAX_RESOL_SHARED_MEM 32
#define BLOCK_SIZE 128         // This is optimally set as big as possible without its ceil:ed multiple exceeding imagesize by too much.

#define NR_CLASS_MUTEXES 5

static pthread_mutex_t global_mutex2[NR_CLASS_MUTEXES] = { PTHREAD_MUTEX_INITIALIZER };
static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

void MlOptimiserCuda::doThreadExpectationSomeParticles(unsigned thread_id)
{
	size_t first_ipart = 0, last_ipart = 0;
	while (baseMLO->exp_ipart_ThreadTaskDistributor->getTasks(first_ipart, last_ipart))
	{
		for (long unsigned ipart = first_ipart; ipart <= last_ipart; ipart++)
		{
			unsigned my_ori_particle = baseMLO->exp_my_first_ori_particle + ipart;
			SamplingParameters sp;
			sp.nr_particles = baseMLO->mydata.ori_particles[my_ori_particle].particles_id.size();

			OptimisationParamters op(sp.nr_particles, my_ori_particle);

			// In the first iteration, multiple seeds will be generated
			// A single random class is selected for each pool of images, and one does not marginalise over the orientations
			// The optimal orientation is based on signal-product (rather than the signal-intensity sensitive Gaussian)
			// If do_firstiter_cc, then first perform a single iteration with K=1 and cross-correlation criteria, afterwards

			// Decide which classes to integrate over (for random class assignment in 1st iteration)
			sp.iclass_max = baseMLO->mymodel.nr_classes - 1;
			// low-pass filter again and generate the seeds
			if (baseMLO->do_generate_seeds)
			{
				if (baseMLO->do_firstiter_cc && baseMLO->iter == 1)
				{
					// In first (CC) iter, use a single reference (and CC)
					sp.iclass_min = sp.iclass_max = 0;
				}
				else if ( (baseMLO->do_firstiter_cc && baseMLO->iter == 2) ||
						(!baseMLO->do_firstiter_cc && baseMLO->iter == 1))
				{
					// In second CC iter, or first iter without CC: generate the seeds
					// Now select a single random class
					// exp_part_id is already in randomized order (controlled by -seed)
					// WARNING: USING SAME iclass_min AND iclass_max FOR SomeParticles!!
					sp.iclass_min = sp.iclass_max = divide_equally_which_group(baseMLO->mydata.numberOfOriginalParticles(), baseMLO->mymodel.nr_classes, op.my_ori_particle);
				}
			}
			// Global exp_metadata array has metadata of all ori_particles. Where does my_ori_particle start?
			for (long int iori = baseMLO->exp_my_first_ori_particle; iori <= baseMLO->exp_my_last_ori_particle; iori++)
			{
				if (iori == my_ori_particle) break;
				op.metadata_offset += baseMLO->mydata.ori_particles[iori].particles_id.size();
			}

			baseMLO->getFourierTransformsAndCtfs(my_ori_particle, op.metadata_offset, op.Fimgs, op.Fimgs_nomask, op.Fctfs,
					op.old_offset, op.prior, op.power_imgs, op.highres_Xi2_imgs,
					op.pointer_dir_nonzeroprior, op.pointer_psi_nonzeroprior, op.directions_prior, op.psi_prior);

			if (baseMLO->do_realign_movies && baseMLO->movie_frame_running_avg_side > 0)
			{
				baseMLO->calculateRunningAveragesOfMovieFrames(my_ori_particle, op.Fimgs, op.power_imgs, op.highres_Xi2_imgs);
			}

			// To deal with skipped alignments/rotations
			if (baseMLO->do_skip_align)
			{
				sp.itrans_min = sp.itrans_max = sp.idir_min = sp.idir_max = sp.ipsi_min = sp.ipsi_max =
						my_ori_particle - baseMLO->exp_my_first_ori_particle;
			}
			else
			{
				sp.itrans_min = 0;
				sp.itrans_max = baseMLO->sampling.NrTranslationalSamplings() - 1;

				if (baseMLO->do_skip_rotate)
				{
					sp.idir_min = sp.idir_max = sp.ipsi_min = sp.ipsi_max =
							my_ori_particle - baseMLO->exp_my_first_ori_particle;
				}
				else
				{
					sp.idir_min = sp.ipsi_min = 0;
					sp.idir_max = baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior) - 1;
					sp.ipsi_max = baseMLO->sampling.NrPsiSamplings(0, &op.pointer_psi_nonzeroprior ) - 1;
				}
			}

			// Initialise significant weight to minus one, so that all coarse sampling points will be handled in the first pass
			op.significant_weight.resize(sp.nr_particles, -1.);

			// Only perform a second pass when using adaptive oversampling
			int nr_sampling_passes = (baseMLO->adaptive_oversampling > 0) ? 2 : 1;

			for (int ipass = 0; ipass < nr_sampling_passes; ipass++)
			{
				if (baseMLO->strict_highres_exp > 0.)
					// Use smaller images in both passes and keep a maximum on coarse_size, just like in FREALIGN
					sp.current_image_size = baseMLO->coarse_size;
				else if (baseMLO->adaptive_oversampling > 0)
					// Use smaller images in the first pass, larger ones in the second pass
					sp.current_image_size = (ipass == 0) ? baseMLO->coarse_size : baseMLO->mymodel.current_size;
				else
					sp.current_image_size = baseMLO->mymodel.current_size;

				// Use coarse sampling in the first pass, oversampled one the second pass
				sp.current_oversampling = (ipass == 0) ? 0 : baseMLO->adaptive_oversampling;

				sp.nr_dir = (baseMLO->do_skip_align || baseMLO->do_skip_rotate) ? 1 : baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior);
				sp.nr_psi = (baseMLO->do_skip_align || baseMLO->do_skip_rotate) ? 1 : baseMLO->sampling.NrPsiSamplings(0, &op.pointer_psi_nonzeroprior);
				sp.nr_trans = (baseMLO->do_skip_align) ? 1 : baseMLO->sampling.NrTranslationalSamplings();
				sp.nr_oversampled_rot = baseMLO->sampling.oversamplingFactorOrientations(sp.current_oversampling);
				sp.nr_oversampled_trans = baseMLO->sampling.oversamplingFactorTranslations(sp.current_oversampling);

				CUDA_CPU_TIC("getAllSquaredDifferences");
				getAllSquaredDifferences(ipass, op, sp);
				CUDA_CPU_TOC("getAllSquaredDifferences");
				CUDA_CPU_TIC("convertAllSquaredDifferencesToWeights");
				convertAllSquaredDifferencesToWeights(ipass, op, sp);
				CUDA_CPU_TOC("convertAllSquaredDifferencesToWeights");
			}

			// For the reconstruction step use mymodel.current_size!
			sp.current_image_size = baseMLO->mymodel.current_size;

			CUDA_CPU_TIC("storeWeightedSums");
			storeWeightedSums(op, sp);
			CUDA_CPU_TOC("storeWeightedSums");
		}
	}
}



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

__global__ void cuda_kernel_projectAllViews_trilin( CudaComplex *g_model,
													FLOAT *g_eulers,
													CudaComplex *g_Frefs,
													int my_r_max,
													int max_r2,
													int min_r2_nn,
													int image_size,
													int orientation_num,
													int XSIZE_img,
													int YSIZE_img,
													int XSIZE_mdl,
													int YSIZE_mdl,
													int STARTINGY_mdl,
													int STARTINGZ_mdl)
{
	FLOAT fx, fy, fz, xp, yp, zp;
	int x0, x1, y0, y1, z0, z1; //y2;
	long int r2;
	int YXSIZE_mdl = XSIZE_mdl*YSIZE_mdl;
	int pixel;
	bool is_neg_x;
	//FLOAT* A;
	CudaComplex d000, d001, d010, d011, d100, d101, d110, d111;
	CudaComplex dx00, dx01, dx10, dx11, dxy0, dxy1, val;

	// blockid
	int ex = blockIdx.y * gridDim.x + blockIdx.x;

	// inside the padded 2D orientation grid
	if( ex < orientation_num ) // we only need to make
	{
		unsigned pass_num(ceilf(   ((float)image_size) / (float)BLOCK_SIZE  ));

		long ref_pixel = ex*(image_size);

		for (unsigned pass = 0; pass < pass_num; pass++) // finish a reference proj in each block
		{
			pixel = (pass * BLOCK_SIZE) + threadIdx.x;

			if(pixel<image_size)
			{
				int x = pixel % XSIZE_img;
				int y = (int)floorf( (float)pixel / (float)XSIZE_img);

				// Dont search beyond square with side max_r
				if (y > my_r_max)
				{
					if (y >= YSIZE_img - my_r_max)
						y = y - YSIZE_img ;
					else
						x=r2;
				}

				r2 = x*x + y*y;
				if (r2 <= max_r2)
				{
					// Get logical coordinates in the 3D map
	//				if(threadIdx.x==0)
	//				{
	//					&A = g_eulers[blockIdx.x*9];
	//				}
					xp = g_eulers[blockIdx.x*9]   * x + g_eulers[blockIdx.x*9+1] * y;  // FIXME: xp,yp,zp has has accuracy loss
					yp = g_eulers[blockIdx.x*9+3] * x + g_eulers[blockIdx.x*9+4] * y;  // compared to CPU-based projection. This
					zp = g_eulers[blockIdx.x*9+6] * x + g_eulers[blockIdx.x*9+7] * y;  // propagates to dx00, dx10, and so on.
					// Only asymmetric half is stored
					if (xp < 0)
					{
						// Get complex conjugated hermitian symmetry pair
						xp = -xp;
						yp = -yp;
						zp = -zp;
						is_neg_x = true;
					}
					else
					{
						is_neg_x = false;
					}
					//is_neg_x = false; //TODO remove after debugging
					// Trilinear interpolation (with physical coords)
					// Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
					// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
					x0 = floorf(xp);
					fx = xp - x0;
					x1 = x0 + 1;

					y0 = floorf(yp);
					fy = yp - y0;
					y0 -=  STARTINGY_mdl;
					y1 = y0 + 1;

					z0 = floorf(zp);
					fz = zp - z0;
					z0 -= STARTINGZ_mdl;
					z1 = z0 + 1;

	//				  P(z,y,x) = z*YXSIZE+y*XSIZE+x
					d000.real = g_model[z0*YXSIZE_mdl+y0*XSIZE_mdl+x0].real;
					d001.real = g_model[z0*YXSIZE_mdl+y0*XSIZE_mdl+x1].real;
					d010.real = g_model[z0*YXSIZE_mdl+y1*XSIZE_mdl+x0].real;
					d011.real = g_model[z0*YXSIZE_mdl+y1*XSIZE_mdl+x1].real;
					d100.real = g_model[z1*YXSIZE_mdl+y0*XSIZE_mdl+x0].real;
					d101.real = g_model[z1*YXSIZE_mdl+y0*XSIZE_mdl+x1].real;
					d110.real = g_model[z1*YXSIZE_mdl+y1*XSIZE_mdl+x0].real;
					d111.real = g_model[z1*YXSIZE_mdl+y1*XSIZE_mdl+x1].real;

					d000.imag = g_model[z0*YXSIZE_mdl+y0*XSIZE_mdl+x0].imag;
					d001.imag = g_model[z0*YXSIZE_mdl+y0*XSIZE_mdl+x1].imag;
					d010.imag = g_model[z0*YXSIZE_mdl+y1*XSIZE_mdl+x0].imag;
					d011.imag = g_model[z0*YXSIZE_mdl+y1*XSIZE_mdl+x1].imag;
					d100.imag = g_model[z1*YXSIZE_mdl+y0*XSIZE_mdl+x0].imag;
					d101.imag = g_model[z1*YXSIZE_mdl+y0*XSIZE_mdl+x1].imag;
					d110.imag = g_model[z1*YXSIZE_mdl+y1*XSIZE_mdl+x0].imag;
					d111.imag = g_model[z1*YXSIZE_mdl+y1*XSIZE_mdl+x1].imag;
					// Set the interpolated value in the 2D output array
					dx00 = d000 + (d001 - d000)*fx;
					dx01 = d100 + (d101 - d100)*fx;
					dx10 = d010 + (d011 - d010)*fx;
					dx11 = d110 + (d111 - d110)*fx;

					dxy0 = dx00 + (dx10 - dx00)*fy;
					dxy1 = dx01 + (dx11 - dx01)*fy;

					val = dxy0 + (dxy1 - dxy0)*fz;
//					val.real= dx00.real;
//					val.imag= dx00.imag;
					if (is_neg_x)
					{
						val.imag = -val.imag;
					}
				}
				else
				{
					val.real=0;
					val.imag=0;
				}
				g_Frefs[ref_pixel+ pixel].real = val.real;
				g_Frefs[ref_pixel+ pixel].imag = val.imag;

			}
		}
	}
}


void MlOptimiserCuda::getAllSquaredDifferences(unsigned exp_ipass, OptimisationParamters &op, SamplingParameters &sp)
{

	CUDA_CPU_TIC("diff_pre_gpu");

	//for scale_correction
	int group_id;

	//printf("sp.nr_oversampled_rot=%d\n", (unsigned)sp.nr_oversampled_rot);

	op.Mweight.resize(sp.nr_particles, baseMLO->mymodel.nr_classes * sp.nr_dir * sp.nr_psi * sp.nr_trans * sp.nr_oversampled_rot * sp.nr_oversampled_trans);
	op.Mweight.initConstant(-999.);
	if (exp_ipass==0)
		op.Mcoarse_significant.clear();

	op.min_diff2.clear();
	op.min_diff2.resize(sp.nr_particles, 99.e99);

	std::vector<MultidimArray<Complex > > dummy;
	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(false, op.my_ori_particle, sp.current_image_size, sp.current_oversampling,
			sp.itrans_min, sp.itrans_max, op.Fimgs, dummy, op.Fctfs, op.local_Fimgs_shifted, dummy,
			op.local_Fctfs, op.local_sqrtXi2, op.local_Minvsigma2s);

	MultidimArray<Complex > Fref;
	Fref.resize(op.local_Minvsigma2s[0]);

	unsigned image_size = op.local_Minvsigma2s[0].nzyxdim;

	CUDA_CPU_TOC("diff_pre_gpu");

	// Loop only from sp.iclass_min to sp.iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
	{
		if (baseMLO->mymodel.pdf_class[exp_iclass] > 0.)
		{
			// Local variables
			std::vector< double > oversampled_rot, oversampled_tilt, oversampled_psi;
			std::vector< double > oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
			CudaGlobalPtr<FLOAT> gpuMinvsigma2(image_size);
			gpuMinvsigma2.device_alloc();

			Matrix2D<FLOAT> A;

			CudaGlobalPtr<CudaComplex> Frefs;

			// Mapping index look-up table
			std::vector< long unsigned > iorientclasses, iover_rots;
			long unsigned orientation_num(0);

			/*=======================================================================================
			                           Generate Reference Projections
			=========================================================================================*/

			CUDA_CPU_TIC("projection_1");
			bool do_gpu_proj=true;
			if(do_gpu_proj)
			{
				CudaGlobalPtr<FLOAT> eulers(9 * sp.nr_dir * sp.nr_psi * sp.nr_oversampled_rot);

				for (long int idir = sp.idir_min, iorient = 0; idir <= sp.idir_max; idir++)
				{
					for (long int ipsi = sp.ipsi_min; ipsi <= sp.ipsi_max; ipsi++, iorient++)
					{
						long int iorientclass = exp_iclass * sp.nr_dir * sp.nr_psi + iorient;

						// Get prior for this direction and skip calculation if prior==0
						double pdf_orientation;
						if (baseMLO->do_skip_align || baseMLO->do_skip_rotate)
						{
							pdf_orientation = baseMLO->mymodel.pdf_class[exp_iclass];
						}
						else if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
						{
							pdf_orientation = DIRECT_MULTIDIM_ELEM(baseMLO->mymodel.pdf_direction[exp_iclass], idir);
						}
						else
						{
							pdf_orientation = op.directions_prior[idir] * op.psi_prior[ipsi];
						}
						// In the first pass, always proceed
						// In the second pass, check whether one of the translations for this orientation of any of the particles had a significant weight in the first pass
						// if so, proceed with projecting the reference in that direction
						bool do_proceed = (exp_ipass==0) ? true :
								baseMLO->isSignificantAnyParticleAnyTranslation(iorientclass, sp.itrans_min, sp.itrans_max, op.Mcoarse_significant);
						if (do_proceed && pdf_orientation > 0.)
						{
							// Now get the oversampled (rot, tilt, psi) triplets
							// This will be only the original (rot,tilt,psi) triplet in the first pass (sp.current_oversampling==0)
							baseMLO->sampling.getOrientations(idir, ipsi, sp.current_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
									op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);

							// Loop over all oversampled orientations (only a single one in the first pass)
							for (long int iover_rot = 0; iover_rot < sp.nr_oversampled_rot; iover_rot++)
							{
								// Get the Euler matrix
								Euler_angles2matrix(oversampled_rot[iover_rot],
													oversampled_tilt[iover_rot],
													oversampled_psi[iover_rot], A);
	//							std::cerr << "A("<< orientation_num <<")=" << A <<  std::endl;
								if(!IS_NOT_INV)
								{
									A = A.transpose();
								}
								A =  A * (FLOAT) baseMLO->mymodel.PPref[exp_iclass].padding_factor;

								for(unsigned i = 0; i < 9; i++)
									eulers[9 * orientation_num + i] = *(A.mdata + i);

	//							std::cerr << "A("<< orientation_num <<")=" << A <<  std::endl;
								orientation_num ++;
								iorientclasses.push_back(iorientclass);
								iover_rots.push_back(iover_rot);
							}
						}
					}
				}
	//			for(int n=0; n<10; n++)
	//			{
	//				for (int m=0; m<9; m++)
	//				std::cerr << "A("<< n << "," << m <<")=" << eulers[9*n+m] <<  std::endl;
	//			}
				int my_r_max = XMIPP_MIN(baseMLO->mymodel.PPref[exp_iclass].r_max, op.local_Minvsigma2s[0].xdim - 1);
				int max_r2 = my_r_max * my_r_max;
				int min_r2_nn = 0; // r_min_nn * r_min_nn;  //FIXME add nn-algorithm

				CudaGlobalPtr<CudaComplex > model((baseMLO->mymodel.PPref[exp_iclass]).data.nzyxdim);
				for(unsigned i = 0; i < model.size; i++)
				{
					model[i].real = (FLOAT) baseMLO->mymodel.PPref[exp_iclass].data.data[i].real;
					model[i].imag = (FLOAT) baseMLO->mymodel.PPref[exp_iclass].data.data[i].imag;
				}

				model.device_alloc();
				model.cp_to_device();

				eulers.size = orientation_num * 9;
				eulers.device_alloc();
				eulers.cp_to_device();

				Frefs.size = orientation_num * image_size;
				Frefs.device_alloc();

				unsigned int orient1, orient2;
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
				std::cerr << "using block dimensions " << orient1 << "," << orient2 <<  std::endl;

				cuda_kernel_projectAllViews_trilin<<<block_dim,BLOCK_SIZE>>>(~model,
																		~eulers,
																		~Frefs,
																		my_r_max,
																		max_r2,
																		min_r2_nn,
																		image_size,
																		orientation_num,
																		op.local_Minvsigma2s[0].xdim,
																		op.local_Minvsigma2s[0].ydim,
																		baseMLO->mymodel.PPref[exp_iclass].data.xdim,
																		baseMLO->mymodel.PPref[exp_iclass].data.ydim,
																		baseMLO->mymodel.PPref[exp_iclass].data.yinit,
																		baseMLO->mymodel.PPref[exp_iclass].data.zinit);
				eulers.free();
				model.free();

				HANDLE_ERROR(cudaDeviceSynchronize());
			}
			else
			{
				Frefs.size = (image_size * sp.nr_dir * sp.nr_psi * sp.nr_oversampled_rot);
				Frefs.host_alloc();

				for (long int idir = sp.idir_min, iorient = 0; idir <= sp.idir_max; idir++)
				{
					for (long int ipsi = sp.ipsi_min; ipsi <= sp.ipsi_max; ipsi++, iorient++)
					{
						long int iorientclass = exp_iclass * sp.nr_dir * sp.nr_psi + iorient;

						// Get prior for this direction and skip calculation if prior==0
						double pdf_orientation;
						if (baseMLO->do_skip_align || baseMLO->do_skip_rotate)
						{
							pdf_orientation = baseMLO->mymodel.pdf_class[exp_iclass];
						}
						else if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
						{
							pdf_orientation = DIRECT_MULTIDIM_ELEM(baseMLO->mymodel.pdf_direction[exp_iclass], idir);
						}
						else
						{
							pdf_orientation = op.directions_prior[idir] * op.psi_prior[ipsi];
						}
						// In the first pass, always proceed
						// In the second pass, check whether one of the translations for this orientation of any of the particles had a significant weight in the first pass
						// if so, proceed with projecting the reference in that direction
						bool do_proceed = (exp_ipass==0) ? true :
								baseMLO->isSignificantAnyParticleAnyTranslation(iorientclass, sp.itrans_min, sp.itrans_max, op.Mcoarse_significant);
						if (do_proceed && pdf_orientation > 0.)
						{
							// Now get the oversampled (rot, tilt, psi) triplets
							// This will be only the original (rot,tilt,psi) triplet in the first pass (sp.current_oversampling==0)
							baseMLO->sampling.getOrientations(idir, ipsi, sp.current_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
									op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);

							// Loop over all oversampled orientations (only a single one in the first pass)
							for (long int iover_rot = 0; iover_rot < sp.nr_oversampled_rot; iover_rot++)
							{
								// Get the Euler matrix
								Euler_angles2matrix(oversampled_rot[iover_rot],
													oversampled_tilt[iover_rot],
													oversampled_psi[iover_rot], A);

								//Fref.data = &Frefs[image_size * orientation_num];
								(baseMLO->mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_NOT_INV);

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
			}

			CUDA_CPU_TOC("projection_1");

			/*=======================================================================================
			                                  	  Particle Iteration
			=========================================================================================*/

			for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
			{
				/*====================================
				        Generate Translations
				======================================*/

				CUDA_CPU_TIC("translation_1");

				CudaGlobalPtr<CudaComplex> Fimgs(image_size * sp.nr_trans * sp.nr_oversampled_trans);

				long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
				long unsigned translation_num(0), ihidden(0);
				std::vector< long unsigned > iover_transes, itranses, ihiddens;

				for (long int itrans = sp.itrans_min; itrans <= sp.itrans_max; itrans++, ihidden++)
				{
					baseMLO->sampling.getTranslations(itrans, sp.current_oversampling,
							oversampled_translations_x, oversampled_translations_y, oversampled_translations_z );

					for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++)
					{
						/// Now get the shifted image
						// Use a pointer to avoid copying the entire array again in this highly expensive loop
						Complex *myAB;
						if (sp.current_oversampling == 0)
						{
							myAB = (Fref.ydim == baseMLO->coarse_size) ? baseMLO->global_fftshifts_ab_coarse[itrans].data
									: baseMLO->global_fftshifts_ab_current[itrans].data;
						}
						else
						{
							int iitrans = itrans * sp.nr_oversampled_trans +  iover_trans;
							myAB = (baseMLO->strict_highres_exp > 0.) ? baseMLO->global_fftshifts_ab2_coarse[iitrans].data
									: baseMLO->global_fftshifts_ab2_current[iitrans].data;
						}


						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
						{
							FLOAT real = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real
									- (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag;
							FLOAT imag = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag
									+ (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real;

							//When on gpu, it makes more sense to ctf-correct translated images, rather than anti-ctf-correct ref-projections
							if (baseMLO->do_scale_correction)
							{
								//group_id = mydata.getGroupId(part_id);
								FLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
								real /= myscale;
								imag /= myscale;
							}
							if (baseMLO->do_ctf_correction && baseMLO->refs_are_ctf_corrected)
							{
								real /= DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n);
								imag /= DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n);
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

				long unsigned coarse_num = sp.nr_dir*sp.nr_psi*sp.nr_trans;
				long unsigned significant_num(0);

				if (exp_ipass == 0)
				{
					op.Mcoarse_significant.resize(coarse_num, 1);
					for (long unsigned i = 0; i < orientation_num; i++)
					{
						for (long unsigned j = 0; j < translation_num; j++)
						{
							rotidx[significant_num] = i;
							transidx[significant_num] = j;
							significant_num++;
//							DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, i)=1;
//							std::cerr << "op.Mcoarse_significant("<< i <<") = " <<    DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, i) << std::endl;
//							std::cerr << "op.Mcoarse_significant("<< i <<") = " << *(&DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, 0)+i*sizeof(bool)) << std::endl;
						}
					}
				}
				else
				{
					for (long unsigned i = 0; i < orientation_num; i++)
					{
						long int iover_rot = iover_rots[i];
//						long int iover_rot = i % sp.nr_oversampled_rot
						long int coarse_rot = floor(i/sp.nr_oversampled_rot);
						for (long unsigned j = 0; j < translation_num; j++)
						{
							long int iover_trans = iover_transes[j];
//							long int iover_trans = j % sp.nr_oversampled_trans
							long int coarse_trans = floor(j/sp.nr_oversampled_trans);
							long int ihidden = iorientclasses[i] * sp.nr_trans + ihiddens[j];
							if(DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, ihidden)==1)
							{
								 long int ihidden_over = baseMLO->sampling.getPositionOversampledSamplingPoint(ihidden,
										                  sp.current_oversampling, iover_rot, iover_trans);

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


				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
				{
					gpuMinvsigma2[n] = *(op.local_Minvsigma2s[ipart].data + n );
				}

				if (baseMLO->do_ctf_correction && baseMLO->refs_are_ctf_corrected)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
					{
						gpuMinvsigma2[n] *= (DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n)*DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n));
					}
				}
				// TODO :    + Assure accuracy with the implemented GPU-based ctf-scaling
				//           + Make setting of myscale robust between here and above.
				//  (scale_correction turns off by default with only one group: ml_optimiser-line 1067,
				//   meaning small-scale test will probably not catch this malfunctioning when/if it breaks.)
				if (baseMLO->do_scale_correction)
				{
					FLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
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

				if ((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc) // do cross-correlation instead of diff
				{
					cuda_kernel_cc_diff2<<<block_dim,BLOCK_SIZE>>>(~Frefs, ~Fimgs, ~gpuMinvsigma2,  ~diff2s,
																	image_size, op.highres_Xi2_imgs[ipart],
																	significant_num,
																	translation_num,
																	~rotidx,
																	~transidx);
				}
				else
				{
					cuda_kernel_diff2<<<block_dim,BLOCK_SIZE>>>(~Frefs, ~Fimgs, ~gpuMinvsigma2, ~diff2s,
																image_size, op.highres_Xi2_imgs[ipart] / 2.,
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
					op.Mcoarse_significant.clear();
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

					long int ihidden = iorientclasses[i] * sp.nr_trans + ihiddens[j];
					long int iover_trans = iover_transes[j];

					long int ihidden_over = baseMLO->sampling.getPositionOversampledSamplingPoint(ihidden, sp.current_oversampling,
																						iover_rot, iover_trans);

					double diff2 = diff2s[i * translation_num + j];

					DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_over) = diff2;

					// Keep track of minimum of all diff2, only for the last image in this series
					if (diff2 < op.min_diff2[ipart])
						op.min_diff2[ipart] = diff2;
				}

				CUDA_CPU_TOC("collect_data_1");

			} // end loop ipart

		} // end if class significant
	} // end loop iclass
}



void MlOptimiserCuda::convertAllSquaredDifferencesToWeights(unsigned exp_ipass, OptimisationParamters &op, SamplingParameters &sp)
{
	op.sum_weight.clear();
	op.sum_weight.resize(sp.nr_particles, 0.);

	// loop over all particles inside this ori_particle
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		double exp_thisparticle_sumweight = 0.;

		double old_offset_z;
		double old_offset_x = XX(op.old_offset[ipart]);
		double old_offset_y = YY(op.old_offset[ipart]);
		if (baseMLO->mymodel.data_dim == 3)
			old_offset_z = ZZ(op.old_offset[ipart]);

		if ((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc)
		{
			// Binarize the squared differences array to skip marginalisation
			double mymindiff2 = 99.e10;
			long int myminidx = -1;
			// Find the smallest element in this row of op.Mweight
			for (long int i = 0; i < XSIZE(op.Mweight); i++)
			{

				double cc = DIRECT_A2D_ELEM(op.Mweight, ipart, i);
				// ignore non-determined cc
				if (cc == -999.)
					continue;

				// just search for the maximum
				if (cc < mymindiff2)
				{
					mymindiff2 = cc;
					myminidx = i;
				}
			}
			// Set all except for the best hidden variable to zero and the smallest element to 1
			for (long int i = 0; i < XSIZE(op.Mweight); i++)
				DIRECT_A2D_ELEM(op.Mweight, ipart, i)= 0.;

			DIRECT_A2D_ELEM(op.Mweight, ipart, myminidx)= 1.;
			exp_thisparticle_sumweight += 1.;

		}
		else
		{
			// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
			for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			{

				// Make PdfOffset calculation much faster...
				double myprior_x, myprior_y, myprior_z;
				if (baseMLO->mymodel.ref_dim == 2)
				{
					myprior_x = XX(baseMLO->mymodel.prior_offset_class[exp_iclass]);
					myprior_y = YY(baseMLO->mymodel.prior_offset_class[exp_iclass]);
				}
				else
				{
					myprior_x = XX(op.prior[ipart]);
					myprior_y = YY(op.prior[ipart]);
					if (baseMLO->mymodel.data_dim == 3)
						myprior_z = ZZ(op.prior[ipart]);
				}
				for (long int idir = sp.idir_min, iorient = 0; idir <= sp.idir_max; idir++)
				{
					for (long int ipsi = sp.ipsi_min; ipsi <= sp.ipsi_max; ipsi++, iorient++)
					{
						long int iorientclass = exp_iclass * sp.nr_dir * sp.nr_psi + iorient;
						double pdf_orientation;

						// Get prior for this direction
						if (baseMLO->do_skip_align || baseMLO->do_skip_rotate)
						{
							pdf_orientation = baseMLO->mymodel.pdf_class[exp_iclass];
						}
						else if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
						{
							pdf_orientation = DIRECT_MULTIDIM_ELEM(baseMLO->mymodel.pdf_direction[exp_iclass], idir);
						}
						else
						{
							// P(orientation) = P(idir|dir_prior) * P(ipsi|psi_prior)
							// This is the probability of the orientation, given the gathered
							// statistics of all assigned orientations of the dataset, since we
							// are assigning a gaussian prior to all parameters.
							pdf_orientation = op.directions_prior[idir] * op.psi_prior[ipsi];
						}
						// Loop over all translations
						long int ihidden = iorientclass * sp.nr_trans;
						for (long int itrans = sp.itrans_min; itrans <= sp.itrans_max; itrans++, ihidden++)
						{

							// To speed things up, only calculate pdf_offset at the coarse sampling.
							// That should not matter much, and that way one does not need to calculate all the OversampledTranslations
							double offset_x = old_offset_x + baseMLO->sampling.translations_x[itrans];
							double offset_y = old_offset_y + baseMLO->sampling.translations_y[itrans];
							double tdiff2 = (offset_x - myprior_x) * (offset_x - myprior_x) + (offset_y - myprior_y) * (offset_y - myprior_y);
							if (baseMLO->mymodel.data_dim == 3)
							{
								double offset_z = old_offset_z + baseMLO->sampling.translations_z[itrans];
								tdiff2 += (offset_z - myprior_z) * (offset_z - myprior_z);
							}
							// P(offset|sigma2_offset)
							// This is the probability of the offset, given the model offset and variance.
							double pdf_offset;
							if (baseMLO->mymodel.sigma2_offset < 0.0001)
								pdf_offset = ( tdiff2 > 0.) ? 0. : 1.; // FIXME ?? A bit hard, no?
							else
								pdf_offset = exp ( tdiff2 / (-2. * baseMLO->mymodel.sigma2_offset) ) / ( 2. * PI * baseMLO->mymodel.sigma2_offset );

							// TMP DEBUGGING
							if (baseMLO->mymodel.orientational_prior_mode != NOPRIOR && (pdf_offset==0. || pdf_orientation==0.))
							{
								pthread_mutex_lock(&global_mutex);
								std::cerr << " pdf_offset= " << pdf_offset << " pdf_orientation= " << pdf_orientation << std::endl;
								std::cerr << " ipart= " << ipart << " part_id= " << part_id << std::endl;
								std::cerr << " iorient= " << iorient << " idir= " << idir << " ipsi= " << ipsi << std::endl;
								//std::cerr << " sp.nr_psi= " << sp.nr_psi << " exp_nr_dir= " << exp_nr_dir << " sp.nr_trans= " << sp.nr_trans << std::endl;
								for (long int i = 0; i < op.directions_prior.size(); i++)
									std::cerr << " op.directions_prior["<<i<<"]= " << op.directions_prior[i] << std::endl;
								for (long int i = 0; i < op.psi_prior.size(); i++)
									std::cerr << " op.psi_prior["<<i<<"]= " << op.psi_prior[i] << std::endl;
								REPORT_ERROR("ERROR! pdf_offset==0.|| pdf_orientation==0.");
								//pthread_mutex_unlock(&global_mutex);
							}
							if (sp.nr_oversampled_rot == 0)
								REPORT_ERROR("sp.nr_oversampled_rot == 0");
							if (sp.nr_oversampled_trans == 0)
								REPORT_ERROR("sp.nr_oversampled_trans == 0");
							// Now first loop over iover_rot, because that is the order in op.Mweight as well
							long int ihidden_over = ihidden * sp.nr_oversampled_rot * sp.nr_oversampled_trans;
							for (long int iover_rot = 0; iover_rot < sp.nr_oversampled_rot; iover_rot++)
							{
								// Then loop over iover_trans
								for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++, ihidden_over++)
								{
									// Only exponentiate for determined values of op.Mweight
									// (this is always true in the first pass, but not so in the second pass)
									// Only deal with this sampling point if its weight was significant
									if (DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_over) < 0.)
									{
										DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_over) = 0.;
									}
									else
									{
										// Set the weight base to the probability of the parameters given the prior
										double weight = pdf_orientation * pdf_offset;
										double diff2 = DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_over) - op.min_diff2[ipart];
										// next line because of numerical precision of exp-function
										if (diff2 > 700.) weight = 0.;
										// TODO: use tabulated exp function?
										else weight *= exp(-diff2);
										// Store the weight
										DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_over) = weight;

										// Keep track of sum and maximum of all weights for this particle
										// Later add all to exp_thisparticle_sumweight, but inside this loop sum to local thisthread_sumweight first
										exp_thisparticle_sumweight += weight;
									} // end if/else op.Mweight < 0.
								} // end loop iover_trans
							}// end loop iover_rot
						} // end loop itrans
					} // end loop ipsi
				} // end loop idir
			} // end loop exp_iclass
		} // end if iter==1

		//Store parameters for this particle
		op.sum_weight[ipart] = exp_thisparticle_sumweight;

#if defined(DEBUG_CUDA) && defined(__linux__)
		if (exp_thisparticle_sumweight == 0. || std::isnan(exp_thisparticle_sumweight))
		{
			printf("DEBUG_ERROR: zero sum of weights.\n");
			exit( EXIT_FAILURE );
		}
#endif

	} // end loop ipart

	// Initialise op.Mcoarse_significant
	if (exp_ipass==0)
		op.Mcoarse_significant.resize(sp.nr_particles, XSIZE(op.Mweight));

	// Now, for each particle,  find the exp_significant_weight that encompasses adaptive_fraction of op.sum_weight
	op.significant_weight.clear();
	op.significant_weight.resize(sp.nr_particles, 0.);
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		MultidimArray<double> sorted_weight;
		// Get the relevant row for this particle
		op.Mweight.getRow(ipart, sorted_weight);

		// Only select non-zero probabilities to speed up sorting
		long int np = 0;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sorted_weight)
		{
			if (DIRECT_MULTIDIM_ELEM(sorted_weight, n) > 0.)
			{
				DIRECT_MULTIDIM_ELEM(sorted_weight, np) = DIRECT_MULTIDIM_ELEM(sorted_weight, n);
				np++;
			}
		}
		sorted_weight.resize(np);

		// Sort from low to high values
		sorted_weight.sort();

		double frac_weight = 0.;
		double my_significant_weight;
		long int my_nr_significant_coarse_samples = 0;
		for (long int i = XSIZE(sorted_weight) - 1; i >= 0; i--)
		{
			if (exp_ipass==0) my_nr_significant_coarse_samples++;
			my_significant_weight = DIRECT_A1D_ELEM(sorted_weight, i);
			frac_weight += my_significant_weight;
			if (frac_weight > baseMLO->adaptive_fraction * op.sum_weight[ipart])
				break;
		}

		if (exp_ipass==0 && my_nr_significant_coarse_samples == 0)
		{
			std::cerr << " ipart= " << ipart << " adaptive_fraction= " << baseMLO->adaptive_fraction << std::endl;
			std::cerr << " frac-weight= " << frac_weight << std::endl;
			std::cerr << " op.sum_weight[ipart]= " << op.sum_weight[ipart] << std::endl;
			Image<double> It;
			std::cerr << " XSIZE(op.Mweight)= " << XSIZE(op.Mweight) << std::endl;
			It()=op.Mweight;
			It() *= 10000;
			It.write("Mweight2.spi");
			std::cerr << "written Mweight2.spi" << std::endl;
			std::cerr << " np= " << np << std::endl;
			It()=sorted_weight;
			It() *= 10000;
			std::cerr << " XSIZE(sorted_weight)= " << XSIZE(sorted_weight) << std::endl;
			if (XSIZE(sorted_weight) > 0)
			{
				It.write("sorted_weight.spi");
				std::cerr << "written sorted_weight.spi" << std::endl;
			}
			REPORT_ERROR("my_nr_significant_coarse_samples == 0");
		}

		if (exp_ipass==0)
		{
			// Store nr_significant_coarse_samples for this particle
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NR_SIGN) = (double)my_nr_significant_coarse_samples;

			// Keep track of which coarse samplings were significant were significant for this particle
			for (int ihidden = 0; ihidden < XSIZE(op.Mcoarse_significant); ihidden++)
			{
				if (DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden) >= my_significant_weight)
					DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, ihidden) = true;
				else
					DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, ihidden) = false;
			}

		}
		op.significant_weight[ipart] = my_significant_weight;
	} // end loop ipart

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

					FLOAT ctf = g_ctfs[pixel];
					CudaComplex ref = g_refs[orientation_pixel];
					if (refs_are_ctf_corrected) //FIXME Create two kernels for the different cases
					{
						ref.real *= ctf;
						ref.imag *= ctf;
					}
					FLOAT diff_real = ref.real - g_imgs[img_pixel_idx].real;
					FLOAT diff_imag = ref.imag - g_imgs[img_pixel_idx].imag;

					wdiff2s_parts += weight * (diff_real*diff_real + diff_imag*diff_imag);

					FLOAT weightxinvsigma2 = weight * ctf * g_Minvsigma2s[pixel];

					wavgs_real += g_imgs_nomask[img_pixel_idx].real * weightxinvsigma2;
					wavgs_imag += g_imgs_nomask[img_pixel_idx].imag * weightxinvsigma2;

					Fweight += weightxinvsigma2 * ctf;
				}
			}

			g_wavgs[orientation_pixel].real = wavgs_real; //TODO should be buffered into shared
			g_wavgs[orientation_pixel].imag = wavgs_imag; //TODO should be buffered into shared
			g_wdiff2s_parts[orientation_pixel] = wdiff2s_parts; //TODO this could be further reduced in here
			g_Fweights[orientation_pixel] = Fweight; //TODO should be buffered into shared
		}
	}
}



void MlOptimiserCuda::storeWeightedSums(OptimisationParamters &op, SamplingParameters &sp)
{
	CUDA_CPU_TIC("store_pre_gpu");

	// Re-do below because now also want unmasked images AND if (stricht_highres_exp >0.) then may need to resize
	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(true, op.my_ori_particle, sp.current_image_size, sp.current_oversampling,
			sp.itrans_min, sp.itrans_max, op.Fimgs, op.Fimgs_nomask, op.Fctfs, op.local_Fimgs_shifted, op.local_Fimgs_shifted_nomask,
			op.local_Fctfs, op.local_sqrtXi2, op.local_Minvsigma2s);

	// In doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s() the origin of the op.local_Minvsigma2s was omitted.
	// Set those back here
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		int group_id = baseMLO->mydata.getGroupId(part_id);
		DIRECT_MULTIDIM_ELEM(op.local_Minvsigma2s[ipart], 0) = 1. / (baseMLO->sigma2_fudge * DIRECT_A1D_ELEM(baseMLO->mymodel.sigma2_noise[group_id], 0));
	}

	// Initialise the maximum of all weights to a negative value
	op.max_weight.clear();
	op.max_weight.resize(sp.nr_particles, -1.);

	// For norm_correction and scale_correction of all particles of this ori_particle
	std::vector<double> exp_wsum_norm_correction;
	std::vector<MultidimArray<double> > exp_wsum_scale_correction_XA, exp_wsum_scale_correction_AA;
	std::vector<MultidimArray<double> > thr_wsum_signal_product_spectra, thr_wsum_reference_power_spectra;
	exp_wsum_norm_correction.resize(sp.nr_particles, 0.);

	// For scale_correction
	if (baseMLO->do_scale_correction)
	{
		MultidimArray<double> aux;
		aux.initZeros(baseMLO->mymodel.ori_size/2 + 1);
		exp_wsum_scale_correction_XA.resize(sp.nr_particles, aux);
		exp_wsum_scale_correction_AA.resize(sp.nr_particles, aux);
		thr_wsum_signal_product_spectra.resize(baseMLO->mymodel.nr_groups, aux);
		thr_wsum_reference_power_spectra.resize(baseMLO->mymodel.nr_groups, aux);
	}

	std::vector< double> oversampled_rot, oversampled_tilt, oversampled_psi;
	std::vector<double> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
	Matrix2D<double> A;
	MultidimArray<Complex > Fimg, Fimg_otfshift_nomask;
	MultidimArray<double> Fweight, Minvsigma2, Mctf;
	bool have_warned_small_scale = false;

	Fimg.resize(op.Fimgs[0]);
	Fweight.resize(op.Fimgs[0]);

	// Initialise Mctf to all-1 for if !do_ctf_corection
	Mctf.resize(op.Fimgs[0]);
	Mctf.initConstant(1.);
	// Initialise Minvsigma2 to all-1 for if !do_map
	Minvsigma2.resize(op.Fimgs[0]);
	Minvsigma2.initConstant(1.);

	// Make local copies of weighted sums (except BPrefs, which are too big)
	// so that there are not too many mutex locks below
	std::vector<MultidimArray<double> > thr_wsum_sigma2_noise, thr_wsum_pdf_direction;
	std::vector<double> thr_wsum_norm_correction, thr_sumw_group, thr_wsum_pdf_class, thr_wsum_prior_offsetx_class, thr_wsum_prior_offsety_class;
	double thr_wsum_sigma2_offset;
	MultidimArray<double> thr_metadata, zeroArray;
	// Wsum_sigma_noise2 is a 1D-spectrum for each group
	zeroArray.initZeros(baseMLO->mymodel.ori_size/2 + 1);
	thr_wsum_sigma2_noise.resize(baseMLO->mymodel.nr_groups, zeroArray);
	// wsum_pdf_direction is a 1D-array (of length sampling.NrDirections()) for each class
	zeroArray.initZeros(baseMLO->sampling.NrDirections());
	thr_wsum_pdf_direction.resize(baseMLO->mymodel.nr_classes, zeroArray);
	// sumw_group is a double for each group
	thr_sumw_group.resize(baseMLO->mymodel.nr_groups, 0.);
	// wsum_pdf_class is a double for each class
	thr_wsum_pdf_class.resize(baseMLO->mymodel.nr_classes, 0.);
	if (baseMLO->mymodel.ref_dim == 2)
	{
		thr_wsum_prior_offsetx_class.resize(baseMLO->mymodel.nr_classes, 0.);
		thr_wsum_prior_offsety_class.resize(baseMLO->mymodel.nr_classes, 0.);
	}
	// wsum_sigma2_offset is just a double
	thr_wsum_sigma2_offset = 0.;

	unsigned image_size = op.Fimgs[0].nzyxdim;

	CUDA_CPU_TOC("store_pre_gpu");

	// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
	{

		/*=======================================================================================
		                            REFERENCE PROJECTION GENERATION
		=======================================================================================*/
		CUDA_CPU_TIC("projection_2");

		CudaGlobalPtr<CudaComplex> Frefs(image_size * sp.nr_dir * sp.nr_psi * sp.nr_oversampled_rot);

		std::vector< long unsigned > iorientclasses, idirs, iover_rots;
		std::vector< double > rots, tilts, psis;
		long unsigned orientation_num(0);

		for (long int idir = sp.idir_min, iorient = 0; idir <= sp.idir_max; idir++)
		{
			for (long int ipsi = sp.ipsi_min; ipsi <= sp.ipsi_max; ipsi++, iorient++)
			{
				long int iorientclass = exp_iclass * sp.nr_dir * sp.nr_psi + iorient;

				baseMLO->sampling.getOrientations(idir, ipsi, baseMLO->adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
						op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);

				if (baseMLO->isSignificantAnyParticleAnyTranslation(iorientclass, sp.itrans_min, sp.itrans_max, op.Mcoarse_significant))
				{
					for (long unsigned iover_rot = 0; iover_rot < sp.nr_oversampled_rot; iover_rot++)
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
						(baseMLO->mymodel.PPref[exp_iclass]).get2DFourierTransform(Fimg, A, IS_NOT_INV);

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
		for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
		{
			long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
			int group_id = baseMLO->mydata.getGroupId(part_id);

			double myprior_x, myprior_y, myprior_z;
			double old_offset_x = XX(op.old_offset[ipart]);
			double old_offset_y = YY(op.old_offset[ipart]);
			double old_offset_z;

			if (baseMLO->mymodel.ref_dim == 2)
			{
				myprior_x = XX(baseMLO->mymodel.prior_offset_class[exp_iclass]);
				myprior_y = YY(baseMLO->mymodel.prior_offset_class[exp_iclass]);
			}
			else
			{
				myprior_x = XX(op.prior[ipart]);
				myprior_y = YY(op.prior[ipart]);
				if (baseMLO->mymodel.data_dim == 3)
				{
					myprior_z = ZZ(op.prior[ipart]);
					old_offset_z = ZZ(op.old_offset[ipart]);
				}
			}


			/*======================================================
								 TRANSLATIONS
			======================================================*/

			CUDA_CPU_TIC("translation_2");

			CudaGlobalPtr<CudaComplex> Fimgs(image_size * sp.nr_trans * sp.nr_oversampled_trans);
			CudaGlobalPtr<CudaComplex> Fimgs_nomask(Fimgs.size);

			long unsigned translation_num(0), ihidden(0);
			std::vector< long unsigned > iover_transes, itranses, ihiddens;

			for (long int itrans = sp.itrans_min, iitrans = 0; itrans <= sp.itrans_max; itrans++, ihidden++)
			{
				baseMLO->sampling.getTranslations(itrans, baseMLO->adaptive_oversampling,
						oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);
				for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++, iitrans++)
				{
					/// Now get the shifted image
					// Use a pointer to avoid copying the entire array again in this highly expensive loop
					Complex* myAB;
					myAB = (baseMLO->adaptive_oversampling == 0 ) ? baseMLO->global_fftshifts_ab_current[iitrans].data : baseMLO->global_fftshifts_ab2_current[iitrans].data;
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
					{
						FLOAT a = (*(myAB + n)).real;
						FLOAT b = (*(myAB + n)).imag;

						// Fimg_shift
						FLOAT real = a * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real
								- b *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag;
						FLOAT imag = a * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag
								+ b *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real;
						Fimgs[translation_num * image_size + n].real = real;
						Fimgs[translation_num * image_size + n].imag = imag;

						// Fimg_shift_nomask
						real = a * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted_nomask[ipart], n)).real
								- b *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted_nomask[ipart], n)).imag;
						imag = a * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted_nomask[ipart], n)).imag
								+ b *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted_nomask[ipart], n)).real;
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
					            	SCALE
			======================================================*/

			FLOAT part_scale(1.);

			if (baseMLO->do_scale_correction)
			{
				part_scale = baseMLO->mymodel.scale_correction[group_id];
				if (part_scale > 10000.)
				{
					std::cerr << " rlnMicrographScaleCorrection= " << part_scale << " group= " << group_id + 1 << std::endl;
					REPORT_ERROR("ERROR: rlnMicrographScaleCorrection is very high. Did you normalize your data?");
				}
				else if (part_scale < 0.001)
				{
					if (!have_warned_small_scale)
					{
						std::cout << " WARNING: ignoring group " << group_id + 1 << " with very small or negative scale (" << part_scale <<
								"); Use larger groups for more stable scale estimates." << std::endl;
						have_warned_small_scale = true;
					}
					part_scale = 0.001;
				}
			}

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
					long unsigned ihidden = iorientclasses[i] * sp.nr_trans + ihiddens[j];
					long unsigned ihidden_over = baseMLO->sampling.getPositionOversampledSamplingPoint(ihidden,
											  sp.current_oversampling, iover_rot, iover_trans);
					sorted_weights[(long unsigned) i * translation_num + j] =
							DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_over);
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

			CudaGlobalPtr<FLOAT> ctfs(image_size); //TODO Same size for all iparts, should be allocated once
			ctfs.device_alloc();

			if (baseMLO->do_ctf_correction)
			{
				for (unsigned i = 0; i < image_size; i++)
					ctfs[i] = (FLOAT) op.local_Fctfs[ipart].data[i] * part_scale;
			}
			else //TODO should be handled by memset
				for (unsigned i = 0; i < image_size; i++)
					ctfs[i] = part_scale;

			ctfs.cp_to_device();

			CudaGlobalPtr<FLOAT> Minvsigma2s(image_size); //TODO Same size for all iparts, should be allocated once
			Minvsigma2s.device_alloc();

			if (baseMLO->do_map)
				for (unsigned i = 0; i < image_size; i++)
					Minvsigma2s[i] = op.local_Minvsigma2s[ipart].data[i];
			else //TODO should be handled by memset
				for (unsigned i = 0; i < image_size; i++)
					Minvsigma2s[i] = 1;

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
												(FLOAT) op.sum_weight[ipart],
												(FLOAT) op.significant_weight[ipart],
												image_size,
												baseMLO->refs_are_ctf_corrected
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
				int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine, j);
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

			// exp_nr_dir * sp.nr_psi * sp.nr_oversampled_rot * sp.nr_trans * sp.nr_oversampled_trans
			for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			{
				for (long int idir = sp.idir_min, iorient = 0; idir <= sp.idir_max; idir++)
				{
					for (long int ipsi = sp.ipsi_min; ipsi <= sp.ipsi_max; ipsi++, iorient++)
					{
						long int iorientclass = exp_iclass * sp.nr_dir * sp.nr_psi + iorient;

						// Only proceed if any of the particles had any significant coarsely sampled translation
						if (baseMLO->isSignificantAnyParticleAnyTranslation(iorientclass, sp.itrans_min, sp.itrans_max, op.Mcoarse_significant))
						{
							// Now get the oversampled (rot, tilt, psi) triplets
							// This will be only the original (rot,tilt,psi) triplet if (adaptive_oversampling==0)
							baseMLO->sampling.getOrientations(idir, ipsi, baseMLO->adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
									op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);
							// Loop over all oversampled orientations (only a single one in the first pass)
							for (long int iover_rot = 0; iover_rot < sp.nr_oversampled_rot; iover_rot++)
							{
								double rot = oversampled_rot[iover_rot];
								double tilt = oversampled_tilt[iover_rot];
								double psi = oversampled_psi[iover_rot];

								// Get the Euler matrix
								Euler_angles2matrix(rot, tilt, psi, A);


								long int ihidden = iorientclass * sp.nr_trans;
								for (long int itrans = sp.itrans_min, iitrans = 0; itrans <= sp.itrans_max; itrans++, ihidden++)
								{
									baseMLO->sampling.getTranslations(itrans, baseMLO->adaptive_oversampling,
											oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);
									for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++, iitrans++)
									{
										// Only deal with this sampling point if its weight was significant
										long int ihidden_over = ihidden * sp.nr_oversampled_trans * sp.nr_oversampled_rot +
												iover_rot * sp.nr_oversampled_trans + iover_trans;

										double weight = DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_over);
										if (weight >= op.significant_weight[ipart])
										{
											weight /= op.sum_weight[ipart];

											// Store sum of weights for this group
											thr_sumw_group[group_id] += weight;
											// Store weights for this class and orientation
											thr_wsum_pdf_class[exp_iclass] += weight;

											// The following goes MUCH faster than the original lines below....
											if (baseMLO->mymodel.ref_dim == 2)
											{
												thr_wsum_prior_offsetx_class[exp_iclass] += weight * (old_offset_x + oversampled_translations_x[iover_trans]);
												thr_wsum_prior_offsety_class[exp_iclass] += weight * (old_offset_y + oversampled_translations_y[iover_trans]);
											}
											double diffx = myprior_x - old_offset_x - oversampled_translations_x[iover_trans];
											double diffy = myprior_y - old_offset_y - oversampled_translations_y[iover_trans];
											if (baseMLO->mymodel.data_dim == 3)
											{
												double diffz  = myprior_z - old_offset_z - oversampled_translations_z[iover_trans];
												thr_wsum_sigma2_offset += weight * (diffx*diffx + diffy*diffy + diffz*diffz);
											}
											else
											{
												thr_wsum_sigma2_offset += weight * (diffx*diffx + diffy*diffy);
											}

											// Store weight for this direction of this class
											if (baseMLO->do_skip_align || baseMLO->do_skip_rotate )
											{
												//ignore pdf_direction
											}
											else if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
											{
												DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], idir) += weight;
											}
											else
											{
												// In the case of orientational priors, get the original number of the direction back
												long int mydir = op.pointer_dir_nonzeroprior[idir];
												DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], mydir) += weight;
											}

											if (weight > op.max_weight[ipart])
											{
												// Store optimal image parameters
												op.max_weight[ipart] = weight;

												A = A.inv();
												A = A.inv();
												Euler_matrix2angles(A, rot, tilt, psi);

												DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_ROT) = rot;
												DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_TILT) = tilt;
												DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_PSI) = psi;
												DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_XOFF) = XX(op.old_offset[ipart]) + oversampled_translations_x[iover_trans];
												DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_YOFF) = YY(op.old_offset[ipart]) + oversampled_translations_y[iover_trans];
												if (baseMLO->mymodel.data_dim == 3)
													DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_ZOFF) = ZZ(op.old_offset[ipart]) + oversampled_translations_z[iover_trans];
												DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CLASS) = (double)exp_iclass + 1;
												DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_PMAX) = op.max_weight[ipart];
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

		wavgs.cp_to_host();
		wavgs.free_device();

		Fweights.cp_to_host();
		Fweights.free_device();

#ifdef RELION_TESTING
		std::string fnm = std::string("gpu_out_exp_wsum_norm_correction.txt");
		char *text = &fnm[0];
		freopen(text,"w",stdout);
		for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
		{
			printf("%4.8f \n",exp_wsum_norm_correction[ipart]);
		}
		fclose(stdout);
		//----------
		fnm = std::string("gpu_out_thr_wsum_sigma2_noise.txt");
		text = &fnm[0];
		freopen(text,"w",stdout);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(baseMLO->Mresol_fine)
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

		CUDA_CPU_TIC("backprojection");

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
			(baseMLO->wsum_model.BPref[exp_iclass]).set2DFourierTransform(Fimg, A, IS_NOT_INV, &Fweight);
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
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		int group_id = baseMLO->mydata.getGroupId(part_id);

		// If the current images were smaller than the original size, fill the rest of wsum_model.sigma2_noise with the power_class spectrum of the images
		for (int ires = baseMLO->mymodel.current_size/2 + 1; ires < baseMLO->mymodel.ori_size/2 + 1; ires++)
		{
			DIRECT_A1D_ELEM(thr_wsum_sigma2_noise[group_id], ires) += DIRECT_A1D_ELEM(op.power_imgs[ipart], ires);
			// Also extend the weighted sum of the norm_correction
			exp_wsum_norm_correction[ipart] += DIRECT_A1D_ELEM(op.power_imgs[ipart], ires);
		}

		// Store norm_correction
		// Multiply by old value because the old norm_correction term was already applied to the image
		if (baseMLO->do_norm_correction)
		{
			double old_norm_correction = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM);
			old_norm_correction /= baseMLO->mymodel.avg_norm_correction;
			// The factor two below is because exp_wsum_norm_correctiom is similar to sigma2_noise, which is the variance for the real/imag components
			// The variance of the total image (on which one normalizes) is twice this value!
			double normcorr = old_norm_correction * sqrt(exp_wsum_norm_correction[ipart] * 2.);
			thr_avg_norm_correction += normcorr;
			// Now set the new norm_correction in the relevant position of exp_metadata
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM) = normcorr;


			// Print warning for strange norm-correction values
			if (!(baseMLO->iter == 1 && baseMLO->do_firstiter_cc) && DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM) > 10.)
			{
				std::cout << " WARNING: norm_correction= "<< DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM) << " for particle " << part_id << " in group " << group_id + 1 << "; Are your groups large enough?" << std::endl;
			}

		}

		// Store weighted sums for scale_correction
		if (baseMLO->do_scale_correction)
		{
			// Divide XA by the old scale_correction and AA by the square of that, because was incorporated into Fctf
			exp_wsum_scale_correction_XA[ipart] /= baseMLO->mymodel.scale_correction[group_id];
			exp_wsum_scale_correction_AA[ipart] /= baseMLO->mymodel.scale_correction[group_id] * baseMLO->mymodel.scale_correction[group_id];

			thr_wsum_signal_product_spectra[group_id] += exp_wsum_scale_correction_XA[ipart];
			thr_wsum_reference_power_spectra[group_id] += exp_wsum_scale_correction_AA[ipart];
		}

		// Calculate DLL for each particle
		double logsigma2 = 0.;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(baseMLO->Mresol_fine)
		{
			int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine, n);
			// Note there is no sqrt in the normalisation term because of the 2-dimensionality of the complex-plane
			// Also exclude origin from logsigma2, as this will not be considered in the P-calculations
			if (ires > 0)
				logsigma2 += log( 2. * PI * DIRECT_A1D_ELEM(baseMLO->mymodel.sigma2_noise[group_id], ires));
		}
		if (op.sum_weight[ipart]==0)
		{
			std::cerr << " part_id= " << part_id << std::endl;
			std::cerr << " ipart= " << ipart << std::endl;
			std::cerr << " op.min_diff2[ipart]= " << op.min_diff2[ipart] << std::endl;
			std::cerr << " logsigma2= " << logsigma2 << std::endl;
			int group_id = baseMLO->mydata.getGroupId(part_id);
			std::cerr << " group_id= " << group_id << std::endl;
			std::cerr << " ml_model.scale_correction[group_id]= " << baseMLO->mymodel.scale_correction[group_id] << std::endl;
			std::cerr << " exp_significant_weight[ipart]= " << op.significant_weight[ipart] << std::endl;
			std::cerr << " exp_max_weight[ipart]= " << op.max_weight[ipart] << std::endl;
			std::cerr << " ml_model.sigma2_noise[group_id]= " << baseMLO->mymodel.sigma2_noise[group_id] << std::endl;
			REPORT_ERROR("ERROR: op.sum_weight[ipart]==0");
		}
		double dLL;
		if ((baseMLO->iter==1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc)
			dLL = -op.min_diff2[ipart];
		else
			dLL = log(op.sum_weight[ipart]) - op.min_diff2[ipart] - logsigma2;

		// Store dLL of each image in the output array, and keep track of total sum
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_DLL) = dLL;
		thr_sum_dLL += dLL;

		// Also store sum of Pmax
		thr_sum_Pmax += DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_PMAX);

	}

	// Now, inside a global_mutex, update the other weighted sums among all threads
	if (!baseMLO->do_skip_maximization)
	{
		pthread_mutex_lock(&global_mutex);
		for (int n = 0; n < baseMLO->mymodel.nr_groups; n++)
		{
			baseMLO->wsum_model.sigma2_noise[n] += thr_wsum_sigma2_noise[n];
			baseMLO->wsum_model.sumw_group[n] += thr_sumw_group[n];
			if (baseMLO->do_scale_correction)
			{
				baseMLO->wsum_model.wsum_signal_product_spectra[n] += thr_wsum_signal_product_spectra[n];
				baseMLO->wsum_model.wsum_reference_power_spectra[n] += thr_wsum_reference_power_spectra[n];
			}
		}
		for (int n = 0; n < baseMLO->mymodel.nr_classes; n++)
		{
			baseMLO->wsum_model.pdf_class[n] += thr_wsum_pdf_class[n];
			if (baseMLO->mymodel.ref_dim == 2)
			{
				XX(baseMLO->wsum_model.prior_offset_class[n]) += thr_wsum_prior_offsetx_class[n];
				YY(baseMLO->wsum_model.prior_offset_class[n]) += thr_wsum_prior_offsety_class[n];
			}

			if (!(baseMLO->do_skip_align || baseMLO->do_skip_rotate) )
				baseMLO->wsum_model.pdf_direction[n] += thr_wsum_pdf_direction[n];
		}
		baseMLO->wsum_model.sigma2_offset += thr_wsum_sigma2_offset;
		if (baseMLO->do_norm_correction)
			baseMLO->wsum_model.avg_norm_correction += thr_avg_norm_correction;
		baseMLO->wsum_model.LL += thr_sum_dLL;
		baseMLO->wsum_model.ave_Pmax += thr_sum_Pmax;
		pthread_mutex_unlock(&global_mutex);
	} // end if !do_skip_maximization

	CUDA_CPU_TOC("store_post_gpu");
}
