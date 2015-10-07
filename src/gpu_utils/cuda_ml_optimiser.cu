#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_projector.cuh"
#include "src/gpu_utils/cuda_projector_plan.h"
#include "src/gpu_utils/cuda_benchmark_utils.cuh"
#include "src/gpu_utils/cuda_ml_optimiser.h"
#include "src/gpu_utils/cuda_kernels/helper.cuh"
#include "src/gpu_utils/cuda_kernels/diff2.cuh"
#include "src/gpu_utils/cuda_kernels/wavg.cuh"
#include "src/gpu_utils/cuda_helper_functions.cu"
#include "src/gpu_utils/cuda_mem_utils.h"
#include "src/complex.h"
#include <fstream>
#include <cuda_runtime.h>
#include "src/parallel.h"
#include <signal.h>
#include <map>

#ifdef CUDA_FORCESTL
#include "src/gpu_utils/cuda_utils_stl.cuh"
#else
#include "src/gpu_utils/cuda_utils_thrust.cuh"
#endif

static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;


void getFourierTransformsAndCtfs(long int my_ori_particle,
		OptimisationParamters &op,
		SamplingParameters &sp,
		MlOptimiser *baseMLO,
		MlOptimiserCuda *cudaMLO
		)
{
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_FT);
#endif
	FourierTransformer transformer;

	for (int ipart = 0; ipart < baseMLO->mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
	{
		CUDA_CPU_TIC("init");
		FileName fn_img;
		Image<double> img, rec_img;
		MultidimArray<Complex > Fimg, Faux;
		MultidimArray<double> Fctf;

		// Get the right line in the exp_fn_img strings (also exp_fn_recimg and exp_fn_ctfs)
		int istop = 0;
		for (long int ii = baseMLO->exp_my_first_ori_particle; ii < my_ori_particle; ii++)
			istop += baseMLO->mydata.ori_particles[ii].particles_id.size();
		istop += ipart;

		// What is my particle_id?
		long int part_id = baseMLO->mydata.ori_particles[my_ori_particle].particles_id[ipart];
		// Which group do I belong?
		int group_id =baseMLO->mydata.getGroupId(part_id);

		// Get the norm_correction
		double normcorr = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM);

		// Get the optimal origin offsets from the previous iteration
		Matrix1D<double> my_old_offset(2), my_prior(2);
		XX(my_old_offset) = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_XOFF);
		YY(my_old_offset) = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_YOFF);
		XX(my_prior)      = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_XOFF_PRIOR);
		YY(my_prior)      = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_YOFF_PRIOR);
		// Uninitialised priors were set to 999.
		if (XX(my_prior) > 998.99 && XX(my_prior) < 999.01)
			XX(my_prior) = 0.;
		if (YY(my_prior) > 998.99 && YY(my_prior) < 999.01)
			YY(my_prior) = 0.;

		if (baseMLO->mymodel.data_dim == 3)
		{
			my_old_offset.resize(3);
			my_prior.resize(3);
			ZZ(my_old_offset) = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_ZOFF);
			ZZ(my_prior)      = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_ZOFF_PRIOR);
			// Unitialised priors were set to 999.
			if (ZZ(my_prior) > 998.99 && ZZ(my_prior) < 999.01)
				ZZ(my_prior) = 0.;
		}
		CUDA_CPU_TOC("init");

		CUDA_CPU_TIC("nonZeroProb");
		if (baseMLO->mymodel.orientational_prior_mode != NOPRIOR && !(baseMLO->do_skip_align ||baseMLO-> do_skip_rotate))
		{
			// First try if there are some fixed prior angles
			double prior_rot = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_ROT_PRIOR);
			double prior_tilt = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_TILT_PRIOR);
			double prior_psi = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_PSI_PRIOR);

			printf("METADATA_ROT_PRIOR=%f\n",prior_rot);

			// If there were no defined priors (i.e. their values were 999.), then use the "normal" angles
			if (prior_rot > 998.99 && prior_rot < 999.01)
				prior_rot = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_ROT);
			if (prior_tilt > 998.99 && prior_tilt < 999.01)
				prior_tilt = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_TILT);
			if (prior_psi > 998.99 && prior_psi < 999.01)
				prior_psi = DIRECT_A2D_ELEM(baseMLO->exp_metadata,op. metadata_offset + ipart, METADATA_PSI);

			printf("METADATA_ROT_PRIOR=%f\n",prior_rot);

			////////// TODO TODO TODO
			////////// How does this work now: each particle has a different sampling object?!!!
			// Select only those orientations that have non-zero prior probability
			baseMLO->sampling.selectOrientationsWithNonZeroPriorProbability(prior_rot, prior_tilt, prior_psi,
					sqrt(baseMLO->mymodel.sigma2_rot), sqrt(baseMLO->mymodel.sigma2_tilt), sqrt(baseMLO->mymodel.sigma2_psi),
					op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);

			long int nr_orients = baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior) * baseMLO->sampling.NrPsiSamplings(0, &op.pointer_psi_nonzeroprior);
			if (nr_orients == 0)
			{
				std::cerr << " sampling.NrDirections()= " << baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior)
						<< " sampling.NrPsiSamplings()= " << baseMLO->sampling.NrPsiSamplings(0, &op.pointer_psi_nonzeroprior) << std::endl;
				REPORT_ERROR("Zero orientations fall within the local angular search. Increase the sigma-value(s) on the orientations!");
			}

		}
		CUDA_CPU_TOC("nonZeroProb");

		// Get the image and recimg data
		if (baseMLO->do_parallel_disc_io)
		{
			CUDA_CPU_TIC("setXmippOrigin");
			// Read from disc
			FileName fn_img;
			std::istringstream split(baseMLO->exp_fn_img);
			for (int i = 0; i <= istop; i++)
				getline(split, fn_img);

			img.read(fn_img);
			img().setXmippOrigin();
			if (baseMLO->has_converged && baseMLO->do_use_reconstruct_images)
			{
				FileName fn_recimg;
				std::istringstream split2(baseMLO->exp_fn_recimg);
				// Get the right line in the exp_fn_img string
				for (int i = 0; i <= istop; i++)
					getline(split2, fn_recimg);
				rec_img.read(fn_recimg);
				rec_img().setXmippOrigin();
			}
			CUDA_CPU_TOC("setXmippOrigin");
		}
		else
		{
			CUDA_CPU_TIC("setXmippOrigin");
			// Unpack the image from the imagedata
			if (baseMLO->mymodel.data_dim == 3)
			{
				img().resize(baseMLO->mymodel.ori_size, baseMLO->mymodel.ori_size,baseMLO-> mymodel.ori_size);
				// Only allow a single image per call of this function!!! nr_pool needs to be set to 1!!!!
				// This will save memory, as we'll need to store all translated images in memory....
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img())
				{
					DIRECT_A3D_ELEM(img(), k, i, j) = DIRECT_A3D_ELEM(baseMLO->exp_imagedata, k, i, j);
				}
				img().setXmippOrigin();

				if (baseMLO->has_converged && baseMLO->do_use_reconstruct_images)
				{
					rec_img().resize(baseMLO->mymodel.ori_size, baseMLO->mymodel.ori_size,baseMLO-> mymodel.ori_size);
					int offset = (baseMLO->do_ctf_correction) ? 2 * baseMLO->mymodel.ori_size : baseMLO->mymodel.ori_size;
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(rec_img())
					{
						DIRECT_A3D_ELEM(rec_img(), k, i, j) = DIRECT_A3D_ELEM(baseMLO->exp_imagedata, offset + k, i, j);
					}
					rec_img().setXmippOrigin();

				}

			}
			else
			{
				img().resize(baseMLO->mymodel.ori_size, baseMLO->mymodel.ori_size);
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(img())
				{
					DIRECT_A2D_ELEM(img(), i, j) = DIRECT_A3D_ELEM(baseMLO->exp_imagedata, op.metadata_offset + ipart, i, j);
				}
				img().setXmippOrigin();
				if (baseMLO->has_converged && baseMLO->do_use_reconstruct_images)
				{

					////////////// TODO: think this through for no-threads here.....
					rec_img().resize(baseMLO->mymodel.ori_size, baseMLO->mymodel.ori_size);
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(rec_img())
					{
						DIRECT_A2D_ELEM(rec_img(), i, j) = DIRECT_A3D_ELEM(baseMLO->exp_imagedata, baseMLO->exp_nr_images + op.metadata_offset + ipart, i, j);
					}
					rec_img().setXmippOrigin();
				}
			}
			CUDA_CPU_TOC("setXmippOrigin");
		}

		// Apply the norm_correction term
		CUDA_CPU_TIC("normCorr");
		if (baseMLO->do_norm_correction)
		{
			img() *=baseMLO->mymodel.avg_norm_correction / normcorr;
		}
		CUDA_CPU_TOC("normCorr");

		CUDA_CPU_TIC("selfTranslate");
		// Apply (rounded) old offsets first
		my_old_offset.selfROUND();
		selfTranslate(img(), my_old_offset, DONT_WRAP);
		if (baseMLO->has_converged && baseMLO->do_use_reconstruct_images)
			selfTranslate(rec_img(), my_old_offset, DONT_WRAP);

		op.old_offset[ipart] = my_old_offset;
		// Also store priors on translations
		op.prior[ipart] = my_prior;
		CUDA_CPU_TOC("selfTranslate");

		CUDA_CPU_TIC("CenterFFT1");
		// Always store FT of image without mask (to be used for the reconstruction)
		MultidimArray<double> img_aux;
		img_aux = (baseMLO->has_converged && baseMLO->do_use_reconstruct_images) ? rec_img() : img();
		CenterFFT(img_aux, true);
		CUDA_CPU_TOC("CenterFFT1");

		CUDA_CPU_TIC("FourierTransform1");

		CUDA_CPU_TIC("setReal");
		transformer.setReal(img_aux);
		CUDA_CPU_TOC("setReal");
		CUDA_CPU_TIC("Transform");
		transformer.Transform(FFTW_FORWARD);
					if (true)transformer.getFourierCopy(Faux);
					else  transformer.getFourierAlias(Faux);
		CUDA_CPU_TOC("Transform");

//TODO CUFFT is not working properly, waiting with this part
//
//		cudaMLO->inputImageData->setSize(XSIZE(img_aux), YSIZE(img_aux));
//		CUDA_CPU_TIC("Memset1");
//		for (unsigned long i = 0; i < cudaMLO->inputImageData->reals.getSize(); i ++)
//			cudaMLO->inputImageData->reals[i] = (cufftReal) img_aux.data[i];
//		CUDA_CPU_TOC("Memset1");
//		cudaMLO->inputImageData->reals.cp_to_device();
//		cudaMLO->inputImageData->forward();
//		cudaMLO->inputImageData->fouriers.cp_to_host();
//		Faux.resize(ZSIZE(img_aux),YSIZE(img_aux),XSIZE(img_aux)/2+1);
//		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));
//		CUDA_CPU_TIC("Memset2");
//		XFLOAT corrFactor = 1. / cudaMLO->inputImageData->reals.getSize();
//		for (unsigned long i = 0; i < cudaMLO->inputImageData->fouriers.getSize(); i ++)
//		{
//			Faux.data[i].real = (double) cudaMLO->inputImageData->fouriers[i].x * corrFactor;
//			Faux.data[i].imag = (double) cudaMLO->inputImageData->fouriers[i].y * corrFactor;
//		}
//		CUDA_CPU_TOC("Memset2");


		CUDA_CPU_TOC("FourierTransform1");


		CUDA_CPU_TIC("windowFourierTransform1");
		windowFourierTransform(Faux, Fimg, baseMLO->mymodel.current_size);
		CUDA_CPU_TOC("windowFourierTransform1");

		CUDA_CPU_TIC("selfApplyBeamTilt");
		// Here apply the beamtilt correction if necessary
		// This will only be used for reconstruction, not for alignment
		// But beamtilt only affects very high-resolution components anyway...
		//
		double beamtilt_x = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_BEAMTILT_X);
		double beamtilt_y = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_BEAMTILT_Y);
		double Cs = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CTF_CS);
		double V = 1000. * DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CTF_VOLTAGE);
		double lambda = 12.2643247 / sqrt(V * (1. + V * 0.978466e-6));
		if (ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.)
			selfApplyBeamTilt(Fimg, beamtilt_x, beamtilt_y, lambda, Cs,baseMLO->mymodel.pixel_size, baseMLO->mymodel.ori_size);

		op.Fimgs_nomask.at(ipart) = Fimg;

		CUDA_CPU_TOC("selfApplyBeamTilt");

		CUDA_CPU_TIC("zeroMask");
		MultidimArray<double> Mnoise;
		if (!baseMLO->do_zero_mask)
		{
			// Make a noisy background image with the same spectrum as the sigma2_noise

			// Different MPI-distributed subsets may otherwise have different instances of the random noise below,
			// because work is on an on-demand basis and therefore variable with the timing of distinct nodes...
			// Have the seed based on the part_id, so that each particle has a different instant of the noise
			if (baseMLO->do_realign_movies)
				init_random_generator(baseMLO->random_seed + part_id);
			else
				init_random_generator(baseMLO->random_seed + my_ori_particle); // This only serves for exact reproducibility tests with 1.3-code...

			// If we're doing running averages, then the sigma2_noise was already adjusted for the running averages.
			// Undo this adjustment here in order to get the right noise in the individual frames
			MultidimArray<double> power_noise = baseMLO->sigma2_fudge * baseMLO->mymodel.sigma2_noise[group_id];
			if (baseMLO->do_realign_movies)
				power_noise *= (2. * baseMLO->movie_frame_running_avg_side + 1.);

			// Create noisy image for outside the mask
			MultidimArray<Complex > Fnoise;
			Mnoise.resize(img());
			transformer.setReal(Mnoise);
			transformer.getFourierAlias(Fnoise);
			// Fill Fnoise with random numbers, use power spectrum of the noise for its variance
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fnoise)
			{
				int ires = ROUND( sqrt( (double)(kp * kp + ip * ip + jp * jp) ) );
				if (ires >= 0 && ires < XSIZE(Fnoise))
				{
					double sigma = sqrt(DIRECT_A1D_ELEM(power_noise, ires));
					DIRECT_A3D_ELEM(Fnoise, k, i, j).real = rnd_gaus(0., sigma);
					DIRECT_A3D_ELEM(Fnoise, k, i, j).imag = rnd_gaus(0., sigma);
				}
				else
				{
					DIRECT_A3D_ELEM(Fnoise, k, i, j) = 0.;
				}
			}
			// Back to real space Mnoise
			CUDA_CPU_TIC("inverseFourierTransform");
			transformer.inverseFourierTransform();
			CUDA_CPU_TOC("inverseFourierTransform");

			CUDA_CPU_TIC("setXmippOrigin");
			Mnoise.setXmippOrigin();
			CUDA_CPU_TOC("setXmippOrigin");

			CUDA_CPU_TIC("softMaskOutsideMap");
			softMaskOutsideMap(img(), baseMLO->particle_diameter / (2. * baseMLO->mymodel.pixel_size), (double)baseMLO->width_mask_edge, &Mnoise);
			CUDA_CPU_TOC("softMaskOutsideMap");
		}
		else
		{
			CUDA_CPU_TIC("softMaskOutsideMap");

			XFLOAT cosine_width = baseMLO->width_mask_edge;
			XFLOAT radius = (XFLOAT)((double)baseMLO->particle_diameter / (2. *baseMLO-> mymodel.pixel_size));
			if (radius < 0)
				radius = ((double)img.data.xdim)/2.;
			XFLOAT radius_p = radius + cosine_width;


			bool do_softmaskOnGpu = true;
			if(do_softmaskOnGpu)
			{
				CudaGlobalPtr<XFLOAT,false> dev_img(img().nzyxdim);
				dev_img.device_alloc();
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img())
					dev_img[n]=(XFLOAT)img.data.data[n];

				dev_img.cp_to_device();
				dim3 block_dim = 1; //TODO
				cuda_kernel_softMaskOutsideMap<<<block_dim,SOFTMASK_BLOCK_SIZE>>>(	~dev_img,
																					img().nzyxdim,
																					img.data.xdim,
																					img.data.ydim,
																					img.data.zdim,
																					img.data.xdim/2,
																					img.data.ydim/2,
																					img.data.zdim/2,
																					true,
																					radius,
																					radius_p,
																					cosine_width);

				dev_img.cp_to_host();
				DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));

				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img())
				{
					img.data.data[n]=(double)dev_img[n];
				}
			}
			else
				softMaskOutsideMap(img(), radius, (double)cosine_width);

//			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img())
//			{
//				std::cout << img.data.data[n] << std::endl;
//			}
//			exit(0);
			CUDA_CPU_TOC("softMaskOutsideMap");
		}
		CUDA_CPU_TOC("zeroMask");

		CUDA_CPU_TIC("CenterFFT2");
		// Inside Projector and Backprojector the origin of the Fourier Transform is centered!
		CenterFFT(img(), true);
		CUDA_CPU_TOC("CenterFFT2");
		CUDA_CPU_TIC("FourierTransform2");
		// Store the Fourier Transform of the image Fimg
		transformer.FourierTransform(img(), Faux);
		CUDA_CPU_TOC("FourierTransform2");

		CUDA_CPU_TIC("powerClass");
		// Store the power_class spectrum of the whole image (to fill sigma2_noise between current_size and ori_size
		if (baseMLO->mymodel.current_size < baseMLO->mymodel.ori_size)
		{
			MultidimArray<double> spectrum;
			spectrum.initZeros(baseMLO->mymodel.ori_size/2 + 1);
			double highres_Xi2 = 0.;
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux)
			{
				int ires = ROUND( sqrt( (double)(kp*kp + ip*ip + jp*jp) ) );
				// Skip Hermitian pairs in the x==0 column

				if (ires > 0 && ires < baseMLO->mymodel.ori_size/2 + 1 && !(jp==0 && ip < 0) )
				{
					double normFaux = norm(DIRECT_A3D_ELEM(Faux, k, i, j));
					DIRECT_A1D_ELEM(spectrum, ires) += normFaux;
					// Store sumXi2 from current_size until ori_size
					if (ires >= baseMLO->mymodel.current_size/2 + 1)
						highres_Xi2 += normFaux;
				}
			}

			// Let's use .at() here instead of [] to check whether we go outside the vectors bounds
			op.power_imgs.at(ipart) = spectrum;
			op.highres_Xi2_imgs.at(ipart) = highres_Xi2;
		}
		else
		{
			op.highres_Xi2_imgs.at(ipart) = 0.;
		}
		CUDA_CPU_TOC("powerClass");
		// We never need any resolutions higher than current_size
		// So resize the Fourier transforms
		CUDA_CPU_TIC("windowFourierTransform2");
		windowFourierTransform(Faux, Fimg, baseMLO->mymodel.current_size);
		CUDA_CPU_TOC("windowFourierTransform2");
		// Also store its CTF
		CUDA_CPU_TIC("ctfCorr");
		Fctf.resize(Fimg);

		// Now calculate the actual CTF
		if (baseMLO->do_ctf_correction)
		{
			if (baseMLO->mymodel.data_dim == 3)
			{
				Image<double> Ictf;
				if (baseMLO->do_parallel_disc_io)
				{
					// Read CTF-image from disc
					FileName fn_ctf;
					std::istringstream split(baseMLO->exp_fn_ctf);
					// Get the right line in the exp_fn_img string
					for (int i = 0; i <= istop; i++)
						getline(split, fn_ctf);
					Ictf.read(fn_ctf);
				}
				else
				{
					// Unpack the CTF-image from the exp_imagedata array
					Ictf().resize(baseMLO->mymodel.ori_size, baseMLO->mymodel.ori_size, baseMLO->mymodel.ori_size);
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Ictf())
					{
						DIRECT_A3D_ELEM(Ictf(), k, i, j) = DIRECT_A3D_ELEM(baseMLO->exp_imagedata, baseMLO->mymodel.ori_size + k, i, j);
					}
				}
				// Set the CTF-image in Fctf
				Ictf().setXmippOrigin();
				FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf)
				{
					// Use negative kp,ip and jp indices, because the origin in the ctf_img lies half a pixel to the right of the actual center....
					DIRECT_A3D_ELEM(Fctf, k, i, j) = A3D_ELEM(Ictf(), -kp, -ip, -jp);
				}
			}
			else
			{
				CTF ctf;
				ctf.setValues(DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CTF_DEFOCUS_U),
							  DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CTF_DEFOCUS_V),
							  DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CTF_DEFOCUS_ANGLE),
							  DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CTF_VOLTAGE),
							  DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CTF_CS),
							  DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CTF_Q0),
							  DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CTF_BFAC));

				ctf.getFftwImage(Fctf, baseMLO->mymodel.ori_size, baseMLO->mymodel.ori_size, baseMLO->mymodel.pixel_size,
						baseMLO->ctf_phase_flipped, baseMLO->only_flip_phases, baseMLO->intact_ctf_first_peak, true);
			}
		}
		else
		{
			Fctf.initConstant(1.);
		}
		CUDA_CPU_TOC("ctfCorr");
		// Store Fimg and Fctf
		op.Fimgs.at(ipart) = Fimg;
		op.Fctfs.at(ipart) = Fctf;

	} // end loop ipart
	transformer.clear();
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_FT);
#endif
}
void getAllSquaredDifferencesCoarse(
		unsigned exp_ipass,
		OptimisationParamters &op,
		SamplingParameters &sp,
		MlOptimiser *baseMLO,
		MlOptimiserCuda *cudaMLO,
	 	std::vector<cudaStager<unsigned long> > &stagerD2)
{

#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF1);
#endif

	CUDA_CPU_TIC("diff_pre_gpu");

	op.Mweight.resize(sp.nr_particles, baseMLO->mymodel.nr_classes * sp.nr_dir * sp.nr_psi * sp.nr_trans * sp.nr_oversampled_rot * sp.nr_oversampled_trans);
	op.Mweight.initConstant(-999.);

	op.Mcoarse_significant.resize(sp.nr_dir*sp.nr_psi*sp.nr_trans, 1);
	op.Mcoarse_significant.clear();

	op.min_diff2.clear();
	op.min_diff2.resize(sp.nr_particles, 99.e99);

	std::vector<MultidimArray<Complex > > dummy;
	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(false, op.my_ori_particle, sp.current_image_size, sp.current_oversampling,
			sp.itrans_min, sp.itrans_max, op.Fimgs, dummy, op.Fctfs, op.local_Fimgs_shifted, dummy,
			op.local_Fctfs, op.local_sqrtXi2, op.local_Minvsigma2s);

	unsigned image_size = op.local_Minvsigma2s[0].nzyxdim;

	CUDA_CPU_TOC("diff_pre_gpu");

	// Loop only from sp.iclass_min to sp.iclass_max to deal with seed generation in first iteration
	CudaGlobalPtr<XFLOAT> allWeights(cudaMLO->allocator);
	for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		allWeights.setSize(cudaMLO->coarseProjectionPlans[exp_iclass].orientation_num * sp.nr_trans*sp.nr_oversampled_trans * sp.nr_particles + allWeights.getSize());

	allWeights.device_alloc();
	allWeights.host_alloc();
	long int allWeights_pos=0;

	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		long int group_id = baseMLO->mydata.getGroupId(part_id);

		/*====================================
				Generate Translations
		======================================*/

		CUDA_CPU_TIC("translation_1");

		long unsigned translation_num((sp.itrans_max - sp.itrans_min + 1) * sp.nr_oversampled_trans);

		CudaGlobalPtr<XFLOAT> Fimgs_real(cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> Fimgs_imag(cudaMLO->allocator);

		Fimgs_real.device_alloc(image_size * translation_num);
		Fimgs_imag.device_alloc(image_size * translation_num);

		if (baseMLO->do_shifts_onthefly)
		{
			CudaTranslator::Plan transPlan(
					op.local_Fimgs_shifted[ipart].data,
					image_size,
					sp.itrans_min * sp.nr_oversampled_trans,
					( sp.itrans_max + 1) * sp.nr_oversampled_trans,
					cudaMLO->allocator,
					0, //stream
					baseMLO->do_scale_correction ? baseMLO->mymodel.scale_correction[group_id] : 1,
					baseMLO->do_ctf_correction && baseMLO->refs_are_ctf_corrected ? op.local_Fctfs[ipart].data : NULL);

			if (sp.current_oversampling == 0)
			{
				if (op.local_Minvsigma2s[0].ydim == baseMLO->coarse_size)
					cudaMLO->translator_coarse1.translate(transPlan, ~Fimgs_real, ~Fimgs_imag);
				else
					cudaMLO->translator_current1.translate(transPlan, ~Fimgs_real, ~Fimgs_imag);
			}
			else
			{
				if (baseMLO->strict_highres_exp > 0.)
					cudaMLO->translator_coarse2.translate(transPlan, ~Fimgs_real, ~Fimgs_imag);
				else
					cudaMLO->translator_current2.translate(transPlan, ~Fimgs_real, ~Fimgs_imag);
			}
		}
		else
		{
			Fimgs_real.host_alloc();
			Fimgs_imag.host_alloc();

			unsigned long k = 0;
			for (unsigned i = 0; i < op.local_Fimgs_shifted.size(); i ++)
			{
				for (unsigned j = 0; j < op.local_Fimgs_shifted[i].nzyxdim; j ++)
				{
					Fimgs_real[k] = op.local_Fimgs_shifted[i].data[j].real;
					Fimgs_imag[k] = op.local_Fimgs_shifted[i].data[j].imag;
					k++;
				}
			}

			Fimgs_real.cp_to_device();
			Fimgs_imag.cp_to_device();
		}

		CUDA_CPU_TOC("translation_1");


		// To speed up calculation, several image-corrections are grouped into a single pixel-wise "filter", or image-correciton
		CudaGlobalPtr<XFLOAT> corr_img(image_size, cudaMLO->allocator);
		corr_img.device_alloc();

		buildCorrImage(baseMLO,op,corr_img,ipart,group_id);
		corr_img.cp_to_device();


		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			CudaProjectorPlan projectorPlan = cudaMLO->coarseProjectionPlans[exp_iclass];

			if ( projectorPlan.orientation_num > 0 )
			{
				/*====================================
				    	   Kernel Call
				======================================*/

				CudaProjectorKernel projKernel = CudaProjectorKernel::makeKernel(
						cudaMLO->cudaProjectors[exp_iclass],
						op.local_Minvsigma2s[0].xdim,
						op.local_Minvsigma2s[0].ydim,
						op.local_Minvsigma2s[0].xdim-1);

				runDiff2KernelCoarse(
						projKernel,
						~corr_img,
						~Fimgs_real,
						~Fimgs_imag,
						~projectorPlan.eulers,
						&allWeights(allWeights_pos),
						op,
						baseMLO,
						projectorPlan.orientation_num,
						translation_num,
						image_size,
						ipart,
						group_id,
						exp_iclass,
						cudaMLO,
						((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc));

				/*====================================
				    	   Retrieve Results
				======================================*/
				allWeights_pos+=projectorPlan.orientation_num*translation_num;

			} // end if class significant
		} // end loop iclass
		cudaDeviceSynchronize();
		allWeights.cp_to_host();
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));
		op.min_diff2[ipart] = getMinOnDevice(allWeights);
		allWeights_pos=0;

		CUDA_CPU_TIC("diff_coarse_map");
		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			CudaProjectorPlan projectorPlan = cudaMLO->coarseProjectionPlans[exp_iclass];
			if ( projectorPlan.orientation_num > 0 )
			{
				for (unsigned i = 0; i < projectorPlan.orientation_num; i ++)
				{
					unsigned iorientclass = projectorPlan.iorientclasses[i];
					for (unsigned j = 0; j < translation_num; j ++)
						DIRECT_A2D_ELEM(op.Mweight, ipart, iorientclass * translation_num + j) = allWeights[i * translation_num + j + allWeights_pos];
				}
				allWeights_pos+=projectorPlan.orientation_num*translation_num;
			}
		}
		CUDA_CPU_TOC("diff_coarse_map");

	} // end loop ipart
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF1);
#endif
}

void getAllSquaredDifferencesFine(unsigned exp_ipass,
		 	 	 	 	 	 	  OptimisationParamters &op,
		 	 	 	 	 	 	  SamplingParameters &sp,
		 	 	 	 	 	 	  MlOptimiser *baseMLO,
		 	 	 	 	 	 	  MlOptimiserCuda *cudaMLO,
		 	 	 	 	 	 	  std::vector<IndexedDataArray> &FinePassWeights,
		 	 	 	 	 	 	  std::vector<std::vector< IndexedDataArrayMask > > &FPCMasks,
		 	 	 	 	 	 	  std::vector<ProjectionParams> &FineProjectionData,
		 	 	 	 	 	 	  std::vector<cudaStager<unsigned long> > &stagerD2)
{
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2);
#endif

	CUDA_CPU_TIC("diff_pre_gpu");

	op.min_diff2.clear();
	op.min_diff2.resize(sp.nr_particles, 99.e99);
	CUDA_CPU_TIC("precalculateShiftedImagesCtfsAndInvSigma2s");
	std::vector<MultidimArray<Complex > > dummy;
	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(false, op.my_ori_particle, sp.current_image_size, sp.current_oversampling,
			sp.itrans_min, sp.itrans_max, op.Fimgs, dummy, op.Fctfs, op.local_Fimgs_shifted, dummy,
			op.local_Fctfs, op.local_sqrtXi2, op.local_Minvsigma2s);
	CUDA_CPU_TOC("precalculateShiftedImagesCtfsAndInvSigma2s");
	MultidimArray<Complex > Fref;
	Fref.resize(op.local_Minvsigma2s[0]);

	unsigned image_size = op.local_Minvsigma2s[0].nzyxdim;

	CUDA_CPU_TOC("diff_pre_gpu");

	/*=======================================================================================
										  Particle Iteration
	=========================================================================================*/
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		// Reset size without de-allocating: we will append everything significant within
		// the current allocation and then re-allocate the then determined (smaller) volume

		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		long int group_id = baseMLO->mydata.getGroupId(part_id);

		/*====================================
				Generate Translations
		======================================*/

		CUDA_CPU_TIC("translation_2");

		long unsigned translation_num((sp.itrans_max - sp.itrans_min + 1) * sp.nr_oversampled_trans);

		CudaGlobalPtr<XFLOAT> Fimgs_real(cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> Fimgs_imag(cudaMLO->allocator);

		Fimgs_real.device_alloc(image_size * translation_num);
		Fimgs_imag.device_alloc(image_size * translation_num);

		if (baseMLO->do_shifts_onthefly)
		{
			CudaTranslator::Plan transPlan(
					op.local_Fimgs_shifted[ipart].data,
					image_size,
					sp.itrans_min * sp.nr_oversampled_trans,
					( sp.itrans_max + 1) * sp.nr_oversampled_trans,
					cudaMLO->allocator,
					0, //stream
					baseMLO->do_scale_correction ? baseMLO->mymodel.scale_correction[group_id] : 1,
					baseMLO->do_ctf_correction && baseMLO->refs_are_ctf_corrected ? op.local_Fctfs[ipart].data : NULL);

			if (sp.current_oversampling == 0)
			{
				if (op.local_Minvsigma2s[0].ydim == baseMLO->coarse_size)
					cudaMLO->translator_coarse1.translate(transPlan, ~Fimgs_real, ~Fimgs_imag);
				else
					cudaMLO->translator_current1.translate(transPlan, ~Fimgs_real, ~Fimgs_imag);
			}
			else
			{
				if (baseMLO->strict_highres_exp > 0.)
					cudaMLO->translator_coarse2.translate(transPlan, ~Fimgs_real, ~Fimgs_imag);
				else
					cudaMLO->translator_current2.translate(transPlan, ~Fimgs_real, ~Fimgs_imag);
			}
		}
		else
		{
			Fimgs_real.host_alloc();
			Fimgs_imag.host_alloc();

			unsigned long k = 0;
			for (unsigned i = 0; i < op.local_Fimgs_shifted.size(); i ++)
			{
				for (unsigned j = 0; j < op.local_Fimgs_shifted[i].nzyxdim; j ++)
				{
					Fimgs_real[k] = op.local_Fimgs_shifted[i].data[j].real;
					Fimgs_imag[k] = op.local_Fimgs_shifted[i].data[j].imag;
					k++;
				}
			}

			Fimgs_real.cp_to_device();
			Fimgs_imag.cp_to_device();
		}

		CUDA_CPU_TOC("translation_2");


		CUDA_CPU_TIC("kernel_init_1");

		CudaGlobalPtr<XFLOAT> corr_img(image_size, cudaMLO->allocator);
		corr_img.device_alloc();
		buildCorrImage(baseMLO,op,corr_img,ipart,group_id);

		corr_img.cp_to_device();

		CUDA_CPU_TOC("kernel_init_1");
		std::vector< CudaGlobalPtr<XFLOAT> > eulers((sp.iclass_max-sp.iclass_min+1), cudaMLO->allocator);
		cudaStager<XFLOAT> AllEulers(cudaMLO->allocator,9*FineProjectionData[ipart].orientationNumAllClasses);
		AllEulers.prepare_device();
		unsigned long newDataSize(0);
		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			FPCMasks[ipart][exp_iclass].weightNum=0;

			if ((baseMLO->mymodel.pdf_class[exp_iclass] > 0.) && (FineProjectionData[ipart].class_entries[exp_iclass] > 0) )
			{
				// use "slice" constructor with class-specific parameters to retrieve a temporary ProjectionParams with data for this class
				ProjectionParams thisClassProjectionData(	FineProjectionData[ipart],
															FineProjectionData[ipart].class_idx[exp_iclass],
															FineProjectionData[ipart].class_idx[exp_iclass]+FineProjectionData[ipart].class_entries[exp_iclass]);
				// since we retrieved the ProjectionParams for *the whole* class the orientation_num is also equal.

				thisClassProjectionData.orientation_num[0] = FineProjectionData[ipart].class_entries[exp_iclass];
				long unsigned orientation_num  = thisClassProjectionData.orientation_num[0];

				if(orientation_num==0)
					continue;

				CUDA_CPU_TIC("pair_list_1");
				long unsigned significant_num(0);
				long int nr_over_orient = baseMLO->sampling.oversamplingFactorOrientations(sp.current_oversampling);
				long int nr_over_trans = baseMLO->sampling.oversamplingFactorTranslations(sp.current_oversampling);
				// Prepare the mask of the weight-array for this class
				if (FPCMasks[ipart][exp_iclass].weightNum==0)
					FPCMasks[ipart][exp_iclass].firstPos = newDataSize;

				long unsigned ihidden(0);
				std::vector< long unsigned > iover_transes, ihiddens;

				for (long int itrans = sp.itrans_min; itrans <= sp.itrans_max; itrans++, ihidden++)
				{
					for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++)
					{
						ihiddens.push_back(ihidden);
						iover_transes.push_back(iover_trans);
					}
				}

				// Do more significance checks on translations and create jobDivision
				significant_num = makeJobsForDiff2Fine(	op,	sp,												// alot of different type inputs...
														orientation_num, translation_num,
														thisClassProjectionData,
														iover_transes, ihiddens,
														nr_over_orient, nr_over_trans, ipart,
														FinePassWeights[ipart],
														FPCMasks[ipart][exp_iclass]);                // ..and output into index-arrays mask

				// extend size by number of significants found this class
				newDataSize += significant_num;
				FPCMasks[ipart][exp_iclass].weightNum = significant_num;
				FPCMasks[ipart][exp_iclass].lastPos = FPCMasks[ipart][exp_iclass].firstPos + significant_num;
				CUDA_CPU_TOC("pair_list_1");

				CUDA_CPU_TIC("IndexedArrayMemCp2");
//				FPCMasks[ipart][exp_iclass].jobOrigin.cp_to_device();
//				FPCMasks[ipart][exp_iclass].jobExtent.cp_to_device();
				stagerD2[ipart].stage(FPCMasks[ipart][exp_iclass].jobOrigin);
				stagerD2[ipart].stage(FPCMasks[ipart][exp_iclass].jobExtent);
				CUDA_CPU_TOC("IndexedArrayMemCp2");

				CUDA_CPU_TIC("generateEulerMatrices");
				eulers[exp_iclass-sp.iclass_min].setSize(9*FineProjectionData[ipart].class_entries[exp_iclass]);
				eulers[exp_iclass-sp.iclass_min].host_alloc();
				generateEulerMatrices(
						baseMLO->mymodel.PPref[exp_iclass].padding_factor,
						thisClassProjectionData,
						&(eulers[exp_iclass-sp.iclass_min])[0],
						!IS_NOT_INV);
				AllEulers.stage(eulers[exp_iclass-sp.iclass_min]);
				CUDA_CPU_TOC("generateEulerMatrices");
			}
		}
		// copy stagers to device
		stagerD2[ipart].cp_to_device();
		AllEulers.cp_to_device();

		FinePassWeights[ipart].rot_id.cp_to_device(); //FIXME this is not used
		FinePassWeights[ipart].rot_idx.cp_to_device();
		FinePassWeights[ipart].trans_idx.cp_to_device();
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));

		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			if ((baseMLO->mymodel.pdf_class[exp_iclass] > 0.) && (FineProjectionData[ipart].class_entries[exp_iclass] > 0) )
			{
				long unsigned orientation_num  = FineProjectionData[ipart].class_entries[exp_iclass];
				if(orientation_num==0)
					continue;

				long unsigned significant_num(FPCMasks[ipart][exp_iclass].weightNum);
				if(significant_num==0)
					continue;

				CUDA_CPU_TIC("Diff2MakeKernel");
				CudaProjectorKernel projKernel = CudaProjectorKernel::makeKernel(
						cudaMLO->cudaProjectors[exp_iclass],
						op.local_Minvsigma2s[0].xdim,
						op.local_Minvsigma2s[0].ydim,
						op.local_Minvsigma2s[0].xdim-1);
				CUDA_CPU_TOC("Diff2MakeKernel");
				CUDA_CPU_TIC("Diff2CALL");

				// Use the constructed mask to construct a partial class-specific input
				IndexedDataArray thisClassFinePassWeights(FinePassWeights[ipart],FPCMasks[ipart][exp_iclass], cudaMLO->allocator);

				runDiff2KernelFine(
						projKernel,
						~corr_img,
						~Fimgs_real,
						~Fimgs_imag,
						~(eulers[exp_iclass-sp.iclass_min]),
						~thisClassFinePassWeights.rot_id,
						~thisClassFinePassWeights.rot_idx,
						~thisClassFinePassWeights.trans_idx,
						~FPCMasks[ipart][exp_iclass].jobOrigin,
						~FPCMasks[ipart][exp_iclass].jobExtent,
						~thisClassFinePassWeights.weights,
						op,
						baseMLO,
						orientation_num,
						translation_num,
						significant_num,
						image_size,
						ipart,
						exp_iclass,
						cudaMLO,
						FPCMasks[ipart][exp_iclass].jobOrigin.getSize(),
						((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc)
						);

//				DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));
				CUDA_CPU_TOC("Diff2CALL");

			} // end if class significant
		} // end loop iclass
		cudaDeviceSynchronize();
		FinePassWeights[ipart].setDataSize( newDataSize );

		CUDA_CPU_TIC("collect_data_1");
		op.min_diff2[ipart] = std::min(op.min_diff2[ipart],(double)getMinOnDevice(FinePassWeights[ipart].weights));
		CUDA_CPU_TOC("collect_data_1");
//		std::cerr << "  fine pass minweight  =  " << op.min_diff2[ipart] << std::endl;

	}// end loop ipart
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF2);
#endif
}


void convertAllSquaredDifferencesToWeights(unsigned exp_ipass,
											OptimisationParamters &op,
											SamplingParameters &sp,
											MlOptimiser *baseMLO,
											MlOptimiserCuda *cudaMLO,
											std::vector< IndexedDataArray> &PassWeights,
											std::vector< std::vector< IndexedDataArrayMask > > &FPCMasks,
											std::vector<cudaStager<unsigned long> > &stagerD2) // FPCMasks = Fine-Pass Class-Masks
{
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
	{
		if (exp_ipass == 0) baseMLO->timer.tic(baseMLO->TIMING_ESP_WEIGHT1);
		else baseMLO->timer.tic(baseMLO->TIMING_ESP_WEIGHT2);
	}
#endif
	op.sum_weight.clear();
	op.sum_weight.resize(sp.nr_particles, 0.);

	CudaGlobalPtr<XFLOAT> allMweight( &(op.Mweight.data[sp.iclass_min*sp.nr_particles * sp.nr_dir * sp.nr_psi * sp.nr_trans]),(sp.iclass_max-sp.iclass_min+1)*sp.nr_particles * sp.nr_dir * sp.nr_psi * sp.nr_trans, cudaMLO->allocator);
	if(exp_ipass==0) // send all the weights in one go, rather than mess about with sending each weight on it's own -- we'll make a new device pointer for each class instead, which is (almost) free
	{
		allMweight.device_alloc();
		allMweight.cp_to_device();
	}

	// Ready the "prior-containers" for all classes (remake every ipart)
	CudaGlobalPtr<XFLOAT>  pdf_orientation((sp.iclass_max-sp.iclass_min+1) * sp.nr_dir * sp.nr_psi, cudaMLO->allocator);
	CudaGlobalPtr<XFLOAT>  pdf_offset((sp.iclass_max-sp.iclass_min+1)*sp.nr_trans, cudaMLO->allocator);
	pdf_orientation.device_alloc();
	pdf_offset.device_alloc();

	// pdf_orientation is ipart-independent, so we keep it above ipart scope
	CUDA_CPU_TIC("get_orient_priors");
	for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		for (long int idir = sp.idir_min, iorientclass = (exp_iclass-sp.iclass_min) * sp.nr_dir * sp.nr_psi; idir <=sp.idir_max; idir++)
			for (long int ipsi = sp.ipsi_min; ipsi <= sp.ipsi_max; ipsi++, iorientclass++)
				if (baseMLO->do_skip_align || baseMLO->do_skip_rotate)
					pdf_orientation[iorientclass] = baseMLO->mymodel.pdf_class[exp_iclass];
				else if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
					pdf_orientation[iorientclass] = DIRECT_MULTIDIM_ELEM(baseMLO->mymodel.pdf_direction[exp_iclass], idir);
				else
					pdf_orientation[iorientclass] = op.directions_prior[idir] * op.psi_prior[ipsi];
	pdf_orientation.cp_to_device();
	CUDA_CPU_TOC("get_orient_priors");

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

			if(exp_ipass==0)
			{
				int nr_coarse_weights = (sp.iclass_max-sp.iclass_min+1)*sp.nr_particles * sp.nr_dir * sp.nr_psi * sp.nr_trans;
				PassWeights[ipart].weights.setDevPtr(&allMweight(ipart*nr_coarse_weights));
				PassWeights[ipart].weights.setHstPtr(&allMweight[ipart*nr_coarse_weights]);
				PassWeights[ipart].weights.setSize(nr_coarse_weights);
			}
			PassWeights[ipart].weights.h_do_free=false;

			std::pair<int, XFLOAT> min_pair=getArgMinOnDevice(PassWeights[ipart].weights);
			PassWeights[ipart].weights.cp_to_host();
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));

			//Set all device-located weights to zero, and only the smallest one to 1.
			DEBUG_HANDLE_ERROR(cudaMemsetAsync(~(PassWeights[ipart].weights), 0.f, PassWeights[ipart].weights.getSize()*sizeof(XFLOAT),0));

			XFLOAT unity=1;
			DEBUG_HANDLE_ERROR(cudaMemcpyAsync( &(PassWeights[ipart].weights(min_pair.first) ), &unity, sizeof(XFLOAT), cudaMemcpyHostToDevice, 0));

			PassWeights[ipart].weights.cp_to_host();
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));
//
//				// Binarize the squared differences array to skip marginalisation
//				double mymindiff2 = 99.e10;
//				long int myminidx = -1;
//				// Find the smallest element in this row of op.Mweight
//				for (long int i = 0; i < XSIZE(op.Mweight); i++)
//				{
//
//					double cc = DIRECT_A2D_ELEM(op.Mweight, ipart, i);
//					// ignore non-determined cc
//					if (cc == -999.)
//						continue;
//
//					// just search for the maximum
//					if (cc < mymindiff2)
//					{
//						mymindiff2 = cc;
//						myminidx = i;
//					}
//				}
//				// Set all except for the best hidden variable to zero and the smallest element to 1
//				for (long int i = 0; i < XSIZE(op.Mweight); i++)
//					DIRECT_A2D_ELEM(op.Mweight, ipart, i)= 0.;
//
//				DIRECT_A2D_ELEM(op.Mweight, ipart, myminidx)= 1.;

			exp_thisparticle_sumweight += 1.;

		}
		else
		{
			long int sumRedSize=0;
			for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
				sumRedSize+= (exp_ipass==0) ? sp.nr_dir*sp.nr_psi*sp.nr_oversampled_rot/SUMW_BLOCK_SIZE : ceil((float)FPCMasks[ipart][exp_iclass].jobNum / (float)SUMW_BLOCK_SIZE);

			CudaGlobalPtr<XFLOAT> thisparticle_sumweight(sumRedSize, cudaMLO->allocator);
			thisparticle_sumweight.device_alloc();
			long int sumweight_pos=0;

			// loop through making translational priors for all classes this ipart - then copy all at once - then loop through kernel calls ( TODO: group kernel calls into one big kernel)
			CUDA_CPU_TIC("get_offset_priors");
			for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			{
				/*=========================================
						Fetch+generate Translation data
				===========================================*/
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

				for (long int itrans = sp.itrans_min; itrans <= sp.itrans_max; itrans++)
				{
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
					if (baseMLO->mymodel.sigma2_offset < 0.0001)
						pdf_offset[(exp_iclass-sp.iclass_min)*sp.nr_trans + itrans] = ( tdiff2 > 0.) ? 0. : 1.;
					else
						pdf_offset[(exp_iclass-sp.iclass_min)*sp.nr_trans + itrans] = exp ( tdiff2 / (-2. * baseMLO->mymodel.sigma2_offset) ) / ( 2. * PI * baseMLO->mymodel.sigma2_offset );
				}
			}
			pdf_offset.cp_to_device();
			CUDA_CPU_TOC("get_offset_priors");

			for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			{
				CUDA_CPU_TIC("sumweight1");
				long int block_num;

				if(exp_ipass==0)  //use Mweight for now - FIXME use PassWeights.weights (ignore indexArrays)
				{
					CudaGlobalPtr<XFLOAT>  Mweight( &allMweight[(ipart)*(op.Mweight).xdim+
					                                  (exp_iclass-sp.iclass_min) * sp.nr_dir * sp.nr_psi * sp.nr_trans],
													&allMweight((ipart)*(op.Mweight).xdim+
													  (exp_iclass-sp.iclass_min) * sp.nr_dir * sp.nr_psi * sp.nr_trans),
													  sp.nr_dir * sp.nr_psi * sp.nr_trans);
					CudaGlobalPtr<XFLOAT>  pdf_orientation_class(&(pdf_orientation[(exp_iclass-sp.iclass_min)*sp.nr_dir*sp.nr_psi]), &( pdf_orientation((exp_iclass-sp.iclass_min)*sp.nr_dir*sp.nr_psi) ), sp.nr_dir*sp.nr_psi);
					CudaGlobalPtr<XFLOAT>  pdf_offset_class(&(pdf_offset[(exp_iclass-sp.iclass_min)*sp.nr_trans]), &( pdf_offset((exp_iclass-sp.iclass_min)*sp.nr_trans) ), sp.nr_trans);

					block_num = sp.nr_dir*sp.nr_psi/SUMW_BLOCK_SIZE;
					dim3 block_dim(block_num);
//					CUDA_GPU_TIC("cuda_kernel_sumweight");
					cuda_kernel_sumweightCoarse<<<block_dim,SUMW_BLOCK_SIZE,0,cudaMLO->classStreams[exp_iclass]>>>(	~pdf_orientation_class,
																			    ~pdf_offset_class,
																			    ~Mweight,
																			    ~thisparticle_sumweight,
																			    (XFLOAT)op.min_diff2[ipart],
																			    sp.nr_dir*sp.nr_psi,
																			    sp.nr_trans,
																			    sumweight_pos);
//					CUDA_GPU_TAC("cuda_kernel_sumweight");
					sumweight_pos+=block_num;
				}
				else if ((baseMLO->mymodel.pdf_class[exp_iclass] > 0.) && (FPCMasks[ipart][exp_iclass].weightNum > 0) )
				{
					// Use the constructed mask to build a partial (class-specific) input
					// (until now, PassWeights has been an empty placeholder. We now create class-paritals pointing at it, and start to fill it with stuff)
					IndexedDataArray thisClassPassWeights(PassWeights[ipart],FPCMasks[ipart][exp_iclass], cudaMLO->allocator);
					CudaGlobalPtr<XFLOAT>  pdf_orientation_class(&(pdf_orientation[(exp_iclass-sp.iclass_min)*sp.nr_dir*sp.nr_psi]), &( pdf_orientation((exp_iclass-sp.iclass_min)*sp.nr_dir*sp.nr_psi) ), sp.nr_dir*sp.nr_psi);
					CudaGlobalPtr<XFLOAT>  pdf_offset_class(&(pdf_offset[(exp_iclass-sp.iclass_min)*sp.nr_trans]), &( pdf_offset((exp_iclass-sp.iclass_min)*sp.nr_trans) ), sp.nr_trans);

					block_num = ceil((float)FPCMasks[ipart][exp_iclass].jobNum / (float)SUMW_BLOCK_SIZE); //thisClassPassWeights.rot_idx.getSize() / SUM_BLOCK_SIZE;
					dim3 block_dim(block_num);

//					CUDA_GPU_TIC("cuda_kernel_sumweight");
					cuda_kernel_sumweightFine<<<block_dim,SUMW_BLOCK_SIZE,0,cudaMLO->classStreams[exp_iclass]>>>(	~pdf_orientation_class,
																			    ~pdf_offset_class,
																			    ~thisClassPassWeights.weights,
																			    ~thisparticle_sumweight,
																			    (XFLOAT)op.min_diff2[ipart],
																			    sp.nr_oversampled_rot,
																			    sp.nr_oversampled_trans,
																			    ~thisClassPassWeights.rot_id,
																			    ~thisClassPassWeights.trans_idx,
																				~FPCMasks[ipart][exp_iclass].jobOrigin,
																			 	~FPCMasks[ipart][exp_iclass].jobExtent,
																			 	FPCMasks[ipart][exp_iclass].jobNum,
																			 	sumweight_pos);
//					CUDA_GPU_TAC("cuda_kernel_sumweight");
					sumweight_pos+=block_num;
				}
				CUDA_CPU_TOC("sumweight1");
			} // end loop exp_iclass
			cudaDeviceSynchronize();

			if(exp_ipass==0)
				allMweight.cp_to_host();
			else
				PassWeights[ipart].weights.cp_to_host();

			exp_thisparticle_sumweight += getSumOnDevice(thisparticle_sumweight);
		}

		//Store parameters for this particle
		op.sum_weight[ipart] = exp_thisparticle_sumweight;
//		std::cerr << "  sumweight =  " << exp_thisparticle_sumweight << std::endl;

#if defined(DEBUG_CUDA) && defined(__linux__)
		if (exp_thisparticle_sumweight == 0. || std::isnan(exp_thisparticle_sumweight))
		{
			printf("DEBUG_ERROR: zero sum of weights.\n");
			raise(SIGSEGV);
		}
#endif

	} // end loop ipart
	if (exp_ipass==0)
	{
		op.Mcoarse_significant.resize(sp.nr_particles, XSIZE(op.Mweight));
	}

	CUDA_CPU_TIC("convertPostKernel");
	// Now, for each particle,  find the exp_significant_weight that encompasses adaptive_fraction of op.sum_weight

	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];

		double frac_weight = 0.;
		double my_significant_weight;

		if ((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc)
		{
			my_significant_weight = 0.999;
			frac_weight = 1.;
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NR_SIGN) = (double) 1.;
			if (exp_ipass==0) // TODO better memset, 0 => false , 1 => true
				for (int ihidden = 0; ihidden < XSIZE(op.Mcoarse_significant); ihidden++)
					if (DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden) >= my_significant_weight)
						DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, ihidden) = true;
					else
						DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, ihidden) = false;
		}
		else if (exp_ipass!=0)
		{
			CUDA_CPU_TIC("sort");
			CudaGlobalPtr<XFLOAT> sorted_weight_new(PassWeights[ipart].weights.getSize(), cudaMLO->allocator);  // make new sorted weights
			sorted_weight_new.device_alloc();
			sortOnDevice(PassWeights[ipart].weights, sorted_weight_new);
			sorted_weight_new.cp_to_host();							// make host-copy
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));

			CUDA_CPU_TOC("sort");
			for (long int i=sorted_weight_new.getSize()-1; i>=0; i--)
			{
//				if (exp_ipass==0) my_nr_significant_coarse_samples++;
					my_significant_weight = sorted_weight_new[i];
				//std::cerr << "thisweight = " << my_significant_weight << std::endl;
				frac_weight += my_significant_weight;
				if (frac_weight > baseMLO->adaptive_fraction * op.sum_weight[ipart])
					break;
			}
		}
		else
		{
			MultidimArray<XFLOAT> sorted_weight;
			long int my_nr_significant_coarse_samples = 0;
			long int np = 0;

			// Get the relevant row for this particle
			CUDA_CPU_TIC("getRow");
			op.Mweight.getRow(ipart, sorted_weight);
			CUDA_CPU_TOC("getRow");

			// Only select non-zero probabilities to speed up sorting // TODO Remove when mapping is eliminated
			CUDA_CPU_TIC("nonZero");
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sorted_weight)
			{
				if (DIRECT_MULTIDIM_ELEM(sorted_weight, n) > 0.)
				{
//					std::cerr << "weights = " << DIRECT_A2D_ELEM(op.Mweight, ipart, n) << std::endl;
					DIRECT_MULTIDIM_ELEM(sorted_weight, np) = DIRECT_MULTIDIM_ELEM(sorted_weight, n);
					np++;
				}
			}
			sorted_weight.resize(np);
			CUDA_CPU_TOC("nonZero");

			//Do the sorting
			CUDA_CPU_TIC("sort");
			CudaGlobalPtr<XFLOAT> sorted_weight_ptr(sorted_weight.data, np, cudaMLO->allocator);
			CudaGlobalPtr<XFLOAT> sorted_weight_ptr_sorted(sorted_weight.data, np, cudaMLO->allocator);

			sorted_weight_ptr.put_on_device();
			sorted_weight_ptr_sorted.device_alloc();

			sortOnDevice(sorted_weight_ptr, sorted_weight_ptr_sorted);

			sorted_weight_ptr_sorted.cp_to_host();
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));
			CUDA_CPU_TOC("sort");

			for (long int i = XSIZE(sorted_weight) - 1; i >= 0; i--)
			{
				if (exp_ipass==0) my_nr_significant_coarse_samples++;
				my_significant_weight = DIRECT_A1D_ELEM(sorted_weight, i);
				//std::cerr << "thisweight = " << my_significant_weight << std::endl;
				frac_weight += my_significant_weight;
				if (frac_weight > baseMLO->adaptive_fraction * op.sum_weight[ipart])
					break;
			}

			if (my_nr_significant_coarse_samples == 0)
			{
				std::cerr << " ipart= " << ipart << " adaptive_fraction= " << baseMLO->adaptive_fraction << std::endl;
				std::cerr << " frac-weight= " << frac_weight << std::endl;
				std::cerr << " op.sum_weight[ipart]= " << op.sum_weight[ipart] << std::endl;
				Image<XFLOAT> It;
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

		op.significant_weight.clear();
		op.significant_weight.resize(sp.nr_particles, 0.);
		op.significant_weight[ipart] = my_significant_weight;
		//std::cerr << "@sort op.significant_weight[ipart]= " << (XFLOAT)op.significant_weight[ipart] << std::endl;

	} // end loop ipart

	CUDA_CPU_TOC("convertPostKernel");
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
	{
		if (exp_ipass == 0) baseMLO->timer.tic(baseMLO->TIMING_ESP_WEIGHT1);
		else baseMLO->timer.tic(baseMLO->TIMING_ESP_WEIGHT2);
	}
#endif
}

void storeWeightedSums(OptimisationParamters &op, SamplingParameters &sp,
						MlOptimiser *baseMLO,
						MlOptimiserCuda *cudaMLO,
						std::vector<IndexedDataArray> &FinePassWeights,
						std::vector<ProjectionParams> &ProjectionData,
						std::vector<std::vector<IndexedDataArrayMask> > FPCMasks,// FIXME? by ref?
	 	 	 	 	 	std::vector<cudaStager<unsigned long> > &stagerSWS)
{
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_WSUM);
#endif
	CUDA_CPU_TIC("store_init");

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

	std::vector<double> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
	bool have_warned_small_scale = false;

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

	CUDA_CPU_TOC("store_init");

	/*=======================================================================================
	                           COLLECT 2 AND SET METADATA
	=======================================================================================*/

	CUDA_CPU_TIC("collect_data_2");
	int nr_transes = sp.nr_trans*sp.nr_oversampled_trans;
	int nr_fake_classes = (sp.iclass_max-sp.iclass_min+1);
	int oversamples = sp.nr_oversampled_trans * sp.nr_oversampled_rot;
	std::vector<long int> block_nums(sp.nr_particles*nr_fake_classes);

	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		// Allocate space for all classes, so that we can pre-calculate data for all classes, copy in one operation, call kenrels on all classes, and copy back in one operation
		CudaGlobalPtr<XFLOAT>          oo_otrans_x(nr_fake_classes*nr_transes, cudaMLO->allocator); // old_offset_oversampled_trans_x
		CudaGlobalPtr<XFLOAT>          oo_otrans_y(nr_fake_classes*nr_transes, cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> myp_oo_otrans_x2y2z2(nr_fake_classes*nr_transes, cudaMLO->allocator); // my_prior_old_offs....x^2*y^2*z^2
		oo_otrans_x.device_alloc();
		oo_otrans_y.device_alloc();
		myp_oo_otrans_x2y2z2.device_alloc();

		int sumBlockNum =0;
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		int group_id = baseMLO->mydata.getGroupId(part_id);
		CUDA_CPU_TIC("collect_data_2_pre_kernel");
		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			int fake_class = exp_iclass-sp.iclass_min; // if we only have the third class to do, the third class will be the "first" we do, i.e. the "fake" first.
			if ((baseMLO->mymodel.pdf_class[exp_iclass] == 0.) || (ProjectionData[ipart].class_entries[exp_iclass] == 0) )
				continue;

			// Use the constructed mask to construct a partial class-specific input
			IndexedDataArray thisClassFinePassWeights(FinePassWeights[ipart],FPCMasks[ipart][exp_iclass], cudaMLO->allocator);

			// Re-define the job-partition of the indexedArray of weights so that the collect-kernel can work with it.
			block_nums[nr_fake_classes*ipart + fake_class] = makeJobsForCollect(thisClassFinePassWeights, FPCMasks[ipart][exp_iclass], ProjectionData[ipart].orientation_num[exp_iclass]);

			stagerSWS[ipart].stage(FPCMasks[ipart][exp_iclass].jobOrigin);
			stagerSWS[ipart].stage(FPCMasks[ipart][exp_iclass].jobExtent);

			sumBlockNum+=block_nums[nr_fake_classes*ipart + fake_class];

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
								COLLECT 2
			======================================================*/

			//Pregenerate oversampled translation objects for kernel-call
			for (long int itrans = 0, iitrans = 0; itrans < sp.nr_trans; itrans++)
			{
				baseMLO->sampling.getTranslations(itrans, baseMLO->adaptive_oversampling,
						oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);
				for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++, iitrans++)
				{
					oo_otrans_x[fake_class*nr_transes+iitrans] = old_offset_x + oversampled_translations_x[iover_trans];
					oo_otrans_y[fake_class*nr_transes+iitrans] = old_offset_y + oversampled_translations_y[iover_trans];
					double diffx = myprior_x - oo_otrans_x[fake_class*nr_transes+iitrans];
					double diffy = myprior_y - oo_otrans_y[fake_class*nr_transes+iitrans];
					if (baseMLO->mymodel.data_dim == 3)
					{
						double diffz = myprior_z - (old_offset_z + oversampled_translations_z[iover_trans]);
						myp_oo_otrans_x2y2z2[fake_class*nr_transes+iitrans] = diffx*diffx + diffy*diffy + diffz*diffz ;
					}
					else
					{
						myp_oo_otrans_x2y2z2[fake_class*nr_transes+iitrans] = diffx*diffx + diffy*diffy ;
					}
				}
			}
		}

		stagerSWS[ipart].cp_to_device();
		oo_otrans_x.cp_to_device();
		oo_otrans_y.cp_to_device();
		myp_oo_otrans_x2y2z2.cp_to_device();
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));

		CudaGlobalPtr<XFLOAT>                      p_weights(sumBlockNum, cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> p_thr_wsum_prior_offsetx_class(sumBlockNum, cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> p_thr_wsum_prior_offsety_class(sumBlockNum, cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT>       p_thr_wsum_sigma2_offset(sumBlockNum, cudaMLO->allocator);
		p_weights.device_alloc();
		p_thr_wsum_prior_offsetx_class.device_alloc();
		p_thr_wsum_prior_offsety_class.device_alloc();
		p_thr_wsum_sigma2_offset.device_alloc();
		CUDA_CPU_TOC("collect_data_2_pre_kernel");
		int partial_pos=0;
		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			int fake_class = exp_iclass-sp.iclass_min; // if we only have the third class to do, the third class will be the "first" we do, i.e. the "fake" first.
			if ((baseMLO->mymodel.pdf_class[exp_iclass] == 0.) || (ProjectionData[ipart].class_entries[exp_iclass] == 0) )
				continue;

			// Use the constructed mask to construct a partial class-specific input
			IndexedDataArray thisClassFinePassWeights(FinePassWeights[ipart],FPCMasks[ipart][exp_iclass], cudaMLO->allocator);

			int cpos=fake_class*nr_transes;
			int block_num = block_nums[nr_fake_classes*ipart + fake_class];
			dim3 grid_dim_collect2 = block_num;// = splitCudaBlocks(block_num,false);
			cuda_kernel_collect2jobs<<<grid_dim_collect2,SUMW_BLOCK_SIZE>>>(
						&(oo_otrans_x(cpos) ),          // otrans-size -> make const
						&(oo_otrans_y(cpos) ),          // otrans-size -> make const
						&(myp_oo_otrans_x2y2z2(cpos) ), // otrans-size -> make const
						~thisClassFinePassWeights.weights,
					(XFLOAT)op.significant_weight[ipart],
					(XFLOAT)op.sum_weight[ipart],
					sp.nr_trans,
					sp.nr_oversampled_trans,
					sp.nr_oversampled_rot,
					oversamples,
					(baseMLO->do_skip_align || baseMLO->do_skip_rotate ),
						&p_weights(partial_pos),
						&p_thr_wsum_prior_offsetx_class(partial_pos),
						&p_thr_wsum_prior_offsety_class(partial_pos),
						&p_thr_wsum_sigma2_offset(partial_pos),
					~thisClassFinePassWeights.rot_idx,
					~thisClassFinePassWeights.trans_idx,
					~FPCMasks[ipart][exp_iclass].jobOrigin,
					~FPCMasks[ipart][exp_iclass].jobExtent
						);
			partial_pos+=block_num;
		}
		CUDA_CPU_TIC("collect_data_2_post_kernel");
		p_weights.cp_to_host();
		p_thr_wsum_prior_offsetx_class.cp_to_host();
		p_thr_wsum_prior_offsety_class.cp_to_host();
		p_thr_wsum_sigma2_offset.cp_to_host();

		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));
		int iorient = 0;
		partial_pos=0;
		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			int fake_class = exp_iclass-sp.iclass_min; // if we only have the third class to do, the third class will be the "first" we do, i.e. the "fake" first.
			if ((baseMLO->mymodel.pdf_class[exp_iclass] == 0.) || (ProjectionData[ipart].class_entries[exp_iclass] == 0) )
				continue;
			int block_num = block_nums[nr_fake_classes*ipart + fake_class];

			for (long int n = partial_pos; n < partial_pos+block_num; n++)
			{
				iorient= FinePassWeights[ipart].rot_id[FPCMasks[ipart][exp_iclass].jobOrigin[n-partial_pos]];

				long int mydir, idir=floor(iorient/sp.nr_psi);
				if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
					mydir = idir;
				else
					mydir = op.pointer_dir_nonzeroprior[idir];

				// store partials according to indices of the relevant dimension
				DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], mydir) += p_weights[n];
				thr_sumw_group[group_id]                 						+= p_weights[n];
				thr_wsum_pdf_class[exp_iclass]           						+= p_weights[n];
				thr_wsum_sigma2_offset                   						+= p_thr_wsum_sigma2_offset[n];

				if (baseMLO->mymodel.ref_dim == 2)
				{
					thr_wsum_prior_offsetx_class[exp_iclass] += p_thr_wsum_prior_offsetx_class[n];
					thr_wsum_prior_offsety_class[exp_iclass] += p_thr_wsum_prior_offsety_class[n];
				}
			}
			partial_pos+=block_num;
		} // end loop iclass
		CUDA_CPU_TOC("collect_data_2_post_kernel");
	} // end loop ipart

	/*======================================================
	                     SET METADATA
	======================================================*/

	std::vector< double> oversampled_rot, oversampled_tilt, oversampled_psi;
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		CUDA_CPU_TIC("setMetadata");

		CUDA_CPU_TIC("getArgMaxOnDevice");
		std::pair<int, XFLOAT> max_pair = getArgMaxOnDevice(FinePassWeights[ipart].weights);
		CUDA_CPU_TOC("getArgMaxOnDevice");

		Indices max_index;
		max_index.fineIdx = FinePassWeights[ipart].ihidden_overs[max_pair.first];
		op.max_weight[ipart] = max_pair.second;

		//std::cerr << "max val = " << op.max_weight[ipart] << std::endl;
		//std::cerr << "max index = " << max_index.fineIdx << std::endl;
		max_index.fineIndexToFineIndices(sp); // set partial indices corresponding to the found max_index, to be used below

		baseMLO->sampling.getTranslations(max_index.itrans, baseMLO->adaptive_oversampling,
				oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);

		//TODO We already have rot, tilt and psi don't calculated them again
		baseMLO->sampling.getOrientations(max_index.idir, max_index.ipsi, baseMLO->adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
				op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);

		double rot = oversampled_rot[max_index.ioverrot];
		double tilt = oversampled_tilt[max_index.ioverrot];
		double psi = oversampled_psi[max_index.ioverrot];
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_ROT) = rot;
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_TILT) = tilt;
		if (psi>180.)
			psi-=360.;
		else if ( psi<-180.)
			psi+=360.;
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_PSI) = psi;
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_XOFF) = XX(op.old_offset[ipart]) + oversampled_translations_x[max_index.iovertrans];
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_YOFF) = YY(op.old_offset[ipart]) + oversampled_translations_y[max_index.iovertrans];

		if (baseMLO->mymodel.data_dim == 3)
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_ZOFF) = ZZ(op.old_offset[ipart]) + oversampled_translations_z[max_index.iovertrans];
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CLASS) = (double)max_index.iclass + 1;
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_PMAX) = op.max_weight[ipart]/op.sum_weight[ipart];

		CUDA_CPU_TOC("setMetadata");
	}

	CUDA_CPU_TOC("collect_data_2");





	/*=======================================================================================
	                                   MAXIMIZATION
	=======================================================================================*/

	CUDA_CPU_TIC("maximization");

	cudaMLO->clearBackprojectDataBundle();

	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		int group_id = baseMLO->mydata.getGroupId(part_id);

		/*======================================================
		                     TRANSLATIONS
		======================================================*/

		CUDA_CPU_TIC("translation_3");

		long unsigned translation_num((sp.itrans_max - sp.itrans_min + 1) * sp.nr_oversampled_trans);

		CudaGlobalPtr<XFLOAT> Fimgs_real(cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> Fimgs_imag(cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> Fimgs_nomask_real(cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> Fimgs_nomask_imag(cudaMLO->allocator);

		Fimgs_real.device_alloc(image_size * translation_num);
		Fimgs_imag.device_alloc(image_size * translation_num);
		Fimgs_nomask_real.device_alloc(image_size * translation_num);
		Fimgs_nomask_imag.device_alloc(image_size * translation_num);

		if (baseMLO->do_shifts_onthefly)
		{
			CudaTranslator::Plan planMask(
					op.local_Fimgs_shifted[ipart].data,
					image_size,
					sp.itrans_min * sp.nr_oversampled_trans,
					( sp.itrans_max + 1) * sp.nr_oversampled_trans,
					cudaMLO->allocator,
					0 //stream
					);

			CudaTranslator::Plan planNomask(
					op.local_Fimgs_shifted_nomask[ipart].data,
					image_size,
					sp.itrans_min * sp.nr_oversampled_trans,
					( sp.itrans_max + 1) * sp.nr_oversampled_trans,
					cudaMLO->allocator,
					0 //stream
					);

			if (baseMLO->adaptive_oversampling == 0)
			{
				cudaMLO->translator_current1.translate(planMask,   ~Fimgs_real,        ~Fimgs_imag);
				cudaMLO->translator_current1.translate(planNomask, ~Fimgs_nomask_real, ~Fimgs_nomask_imag);
			}
			else
			{
				cudaMLO->translator_current2.translate(planMask,   ~Fimgs_real,        ~Fimgs_imag);
				cudaMLO->translator_current2.translate(planNomask, ~Fimgs_nomask_real, ~Fimgs_nomask_imag);
			}
		}
		else
		{
			Fimgs_real.host_alloc();
			Fimgs_imag.host_alloc();
			Fimgs_nomask_real.host_alloc();
			Fimgs_nomask_imag.host_alloc();

			unsigned long k = 0;
			for (unsigned i = 0; i < op.local_Fimgs_shifted.size(); i ++)
			{
				for (unsigned j = 0; j < op.local_Fimgs_shifted[i].nzyxdim; j ++)
				{
					Fimgs_real[k] = op.local_Fimgs_shifted[i].data[j].real;
					Fimgs_imag[k] = op.local_Fimgs_shifted[i].data[j].imag;
					Fimgs_nomask_real[k] = op.local_Fimgs_shifted_nomask[i].data[j].real;
					Fimgs_nomask_imag[k] = op.local_Fimgs_shifted_nomask[i].data[j].imag;
					k++;
				}
			}

			Fimgs_real.cp_to_device();
			Fimgs_imag.cp_to_device();
			Fimgs_nomask_real.cp_to_device();
			Fimgs_nomask_imag.cp_to_device();
		}

		CUDA_CPU_TOC("translation_3");


		/*======================================================
		                       SCALE
		======================================================*/

		XFLOAT part_scale(1.);

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

		CudaGlobalPtr<XFLOAT> ctfs(image_size, cudaMLO->allocator); //TODO Same size for all iparts, should be allocated once

		if (baseMLO->do_ctf_correction)
		{
			for (unsigned i = 0; i < image_size; i++)
				ctfs[i] = (XFLOAT) op.local_Fctfs[ipart].data[i] * part_scale;
		}
		else //TODO should be handled by memset
			for (unsigned i = 0; i < image_size; i++)
				ctfs[i] = part_scale;

		ctfs.put_on_device();

		/*======================================================
		                       MINVSIGMA
		======================================================*/

		CudaGlobalPtr<XFLOAT> Minvsigma2s(image_size, cudaMLO->allocator); //TODO Same size for all iparts, should be allocated once

		if (baseMLO->do_map)
			for (unsigned i = 0; i < image_size; i++)
				Minvsigma2s[i] = op.local_Minvsigma2s[ipart].data[i];
		else
			for (unsigned i = 0; i < image_size; i++)
				Minvsigma2s[i] = 1;

		Minvsigma2s.put_on_device();

		/*======================================================
		                      CLASS LOOP
		======================================================*/

		CudaGlobalPtr<XFLOAT> wdiff2s_AA(ProjectionData[ipart].orientationNumAllClasses*image_size, 0, cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> wdiff2s_XA(ProjectionData[ipart].orientationNumAllClasses*image_size, 0, cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> wdiff2s_sum(image_size, 0, cudaMLO->allocator);

		wdiff2s_AA.device_alloc();
		wdiff2s_AA.device_init(0.f);
		wdiff2s_XA.device_alloc();
		wdiff2s_XA.device_init(0.f);

		unsigned long AAXA_pos=0;

		wdiff2s_sum.device_alloc();
		wdiff2s_sum.device_init(0.f);

		// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			if((baseMLO->mymodel.pdf_class[exp_iclass] == 0.) || (ProjectionData[ipart].class_entries[exp_iclass] == 0))
				continue;

			// Use the constructed mask to construct a partial class-specific input
			IndexedDataArray thisClassFinePassWeights(FinePassWeights[ipart],FPCMasks[ipart][exp_iclass], cudaMLO->allocator);

			CUDA_CPU_TIC("thisClassProjectionSetupCoarse");
			// use "slice" constructor with class-specific parameters to retrieve a temporary ProjectionParams with data for this class
			ProjectionParams thisClassProjectionData(	ProjectionData[ipart],
														ProjectionData[ipart].class_idx[exp_iclass],
														ProjectionData[ipart].class_idx[exp_iclass]+ProjectionData[ipart].class_entries[exp_iclass]);

			thisClassProjectionData.orientation_num[0] = ProjectionData[ipart].orientation_num[exp_iclass];
			CUDA_CPU_TOC("thisClassProjectionSetupCoarse");

			long unsigned orientation_num(thisClassProjectionData.orientation_num[0]);

			/*======================================================
								PROJECTIONS
			======================================================*/

			BackprojectDataBundle *dataBundle = new BackprojectDataBundle(
					orientation_num * image_size,
					orientation_num * 9,
					0,
					cudaMLO->allocator);

			CUDA_CPU_TIC("generateEulerMatricesProjector");

			generateEulerMatrices(
					baseMLO->mymodel.PPref[exp_iclass].padding_factor,
					thisClassProjectionData,
					&dataBundle->eulers[0],
					!IS_NOT_INV);

			CUDA_CPU_TOC("generateEulerMatricesProjector");

			dataBundle->eulers.device_alloc_end();
			dataBundle->reals.device_alloc_end();
			dataBundle->imags.device_alloc_end();
			dataBundle->weights.device_alloc_end();

			dataBundle->eulers.cp_to_device();


			/*======================================================
								 MAP WEIGHTS
			======================================================*/

			CUDA_CPU_TIC("pre_wavg_map");
			CudaGlobalPtr<XFLOAT> sorted_weights(orientation_num * translation_num, 0, cudaMLO->allocator);

			for (long unsigned i = 0; i < orientation_num*translation_num; i++)
				sorted_weights[i] = -999.;

			for (long unsigned i = 0; i < thisClassFinePassWeights.weights.getSize(); i++)
				sorted_weights[ (thisClassFinePassWeights.rot_idx[i]) * translation_num + thisClassFinePassWeights.trans_idx[i] ]
								= thisClassFinePassWeights.weights[i];

			sorted_weights.put_on_device();
			CUDA_CPU_TOC("pre_wavg_map");

			/*======================================================
								 KERNEL CALL
			======================================================*/

			CudaProjectorKernel projKernel = CudaProjectorKernel::makeKernel(
					cudaMLO->cudaProjectors[exp_iclass],
					op.local_Minvsigma2s[0].xdim,
					op.local_Minvsigma2s[0].ydim,
					op.local_Minvsigma2s[0].xdim-1);

			runWavgKernel(
					projKernel,
					~dataBundle->eulers,
					~Fimgs_real,
					~Fimgs_imag,
					~Fimgs_nomask_real,
					~Fimgs_nomask_imag,
					~sorted_weights,
					~ctfs,
					~Minvsigma2s,
					~wdiff2s_sum,
					&wdiff2s_AA(AAXA_pos),
					&wdiff2s_XA(AAXA_pos),
					~dataBundle->reals,
					~dataBundle->imags,
					~dataBundle->weights,
					op,
					baseMLO,
					orientation_num,
					translation_num,
					image_size,
					ipart,
					group_id,
					exp_iclass,
					part_scale,
					0);

			AAXA_pos+=orientation_num*image_size;

			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));

			/*======================================================
								BACKPROJECTION
			======================================================*/

#ifdef TIMING
			if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
				baseMLO->timer.tic(baseMLO->TIMING_WSUM_BACKPROJ);
#endif

			CUDA_CPU_TIC("backproject");

			cudaMLO->cudaBackprojectors[exp_iclass].backproject(
				~dataBundle->reals,
				~dataBundle->imags,
				~dataBundle->weights,
				~dataBundle->eulers,
				op.local_Minvsigma2s[0].xdim,
				op.local_Minvsigma2s[0].ydim,
				orientation_num);

//			cudaMLO->backprojectDataBundleStack.push(dataBundle);
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaMLO->cudaBackprojectors[exp_iclass].getStream()));
			delete dataBundle;

			CUDA_CPU_TOC("backproject");

#ifdef TIMING
			if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
				baseMLO->timer.toc(baseMLO->TIMING_WSUM_BACKPROJ);
#endif
		} // end loop iclass

		wdiff2s_AA.cp_to_host();
		wdiff2s_XA.cp_to_host();
		wdiff2s_sum.cp_to_host();
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(0));

		AAXA_pos=0;

		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			if((baseMLO->mymodel.pdf_class[exp_iclass] == 0.) || (ProjectionData[ipart].class_entries[exp_iclass] == 0))
				continue;
			for (long int j = 0; j < image_size; j++)
			{
				int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine, j);
				if (baseMLO->do_scale_correction && DIRECT_A1D_ELEM(baseMLO->mymodel.data_vs_prior_class[exp_iclass], ires) > 3.)
				{
					DIRECT_A1D_ELEM(exp_wsum_scale_correction_AA[ipart], ires) += wdiff2s_AA[AAXA_pos+j];
					DIRECT_A1D_ELEM(exp_wsum_scale_correction_XA[ipart], ires) += wdiff2s_XA[AAXA_pos+j];
				}
			}
			AAXA_pos+=ProjectionData[ipart].orientation_num[exp_iclass]*image_size;
		} // end loop iclass
		for (long int j = 0; j < image_size; j++)
		{
			int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine, j);
			if (ires > -1)
			{
				thr_wsum_sigma2_noise[group_id].data[ires] += (double) wdiff2s_sum[j];
				exp_wsum_norm_correction[ipart] += (double) wdiff2s_sum[j]; //TODO could be gpu-reduced
			}
		}
	} // end loop ipart
	CUDA_CPU_TOC("maximization");




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
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_WSUM);
#endif
}

MlOptimiserCuda::MlOptimiserCuda(MlOptimiser *baseMLOptimiser, int dev_id) :
		baseMLO(baseMLOptimiser),
		generateProjectionPlanOnTheFly(false)
{
	unsigned nr_classes = baseMLOptimiser->mymodel.nr_classes;

	/*======================================================
					DEVICE MEM OBJ SETUP
	======================================================*/

	device_id = dev_id;
	int devCount;

	HANDLE_ERROR(cudaGetDeviceCount(&devCount));

	if(dev_id >= devCount)
	{
		std::cerr << " using device_id=" << dev_id << " (device no. " << dev_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		REPORT_ERROR("ERROR: Assigning a thread to a non-existent device (index likely too high)");
	}
	else
		HANDLE_ERROR(cudaSetDevice(dev_id));

	classStreams.resize(nr_classes, 0);
	for (int i = 0; i <= nr_classes; i++)
		HANDLE_ERROR(cudaStreamCreate(&classStreams[i]));

	bpStreams.resize(nr_classes, 0);
	for (int i = 0; i <= nr_classes; i++)
		HANDLE_ERROR(cudaStreamCreateWithPriority(&bpStreams[i], cudaStreamNonBlocking, 1)); //Lower priority stream (1)

	refIs3D = baseMLO->mymodel.ref_dim == 3;

	cudaProjectors.resize(nr_classes);
	cudaBackprojectors.resize(nr_classes);

	//Loop over classes
	for (int iclass = 0; iclass < nr_classes; iclass++)
		cudaBackprojectors[iclass].setStream(bpStreams[iclass]);

	/*======================================================
	                    CUSTOM ALLOCATOR
	======================================================*/

#ifdef CUDA_NO_CUSTOM_ALLOCATION
	printf(" Custom allocator is disabled.\n");
	allocator = new CudaCustomAllocator(0);
#else
	size_t allocationSize(0);

	size_t free, total;
	HANDLE_ERROR(cudaMemGetInfo( &free, &total ));

	if (baseMLO->available_gpu_memory > 0)
		allocationSize = baseMLO->available_gpu_memory * (1000*1000*1000);
	else
		allocationSize = (float)free * .7;

	if (allocationSize > free)
	{
		printf(" WARNING: Required memory per thread, via \"--gpu_memory_per_thread\", not available on device. (Defaulting to less)\n");
		allocationSize = (float)free * .7; //Lets leave some for other processes for now
	}

	printf(" Custom allocator assigned %.2f MB of device memory (on device %d).\n", (float)allocationSize/(1000.*1000.), device_id);

	allocator = new CudaCustomAllocator(allocationSize);
	allocator->setOutOfMemoryHandler(this);
#endif

	//inputImageData = new CufftBundle(0, allocator);
};

void MlOptimiserCuda::resetData()
{
	unsigned nr_classes = baseMLO->mymodel.nr_classes;

	HANDLE_ERROR(cudaSetDevice(device_id));

	/*======================================================
	   PROJECTOR, PROJECTOR PLAN AND BACKPROJECTOR SETUP
	======================================================*/

	//Can we pre-generate projector plan and corresponding euler matrices for all particles
	if (baseMLO->do_skip_align || baseMLO->do_skip_rotate || baseMLO->do_auto_refine || baseMLO->mymodel.orientational_prior_mode != NOPRIOR)
		generateProjectionPlanOnTheFly = true;
	else
		generateProjectionPlanOnTheFly = false;

	coarseProjectionPlans.clear();

#ifdef DEBUG_CUDA
		if (allocator->getNumberOfAllocs() != 0)
		{
			printf("DEBUG_ERROR: Non-zero allocation count encountered in custom allocator between iterations.\n");
			allocator->printState();
			fflush(stdout);
			raise(SIGSEGV);
		}
#endif

	coarseProjectionPlans.resize(nr_classes, allocator);

	//Loop over classes
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		cudaProjectors[iclass].setMdlDim(
				baseMLO->mymodel.PPref[iclass].data.xdim,
				baseMLO->mymodel.PPref[iclass].data.ydim,
				baseMLO->mymodel.PPref[iclass].data.zdim,
				baseMLO->mymodel.PPref[iclass].data.yinit,
				baseMLO->mymodel.PPref[iclass].data.zinit,
				baseMLO->mymodel.PPref[iclass].r_max,
				baseMLO->mymodel.PPref[iclass].padding_factor);

		cudaProjectors[iclass].initMdl(baseMLO->mymodel.PPref[iclass].data.data);

		cudaBackprojectors[iclass].setMdlDim(
				baseMLO->wsum_model.BPref[iclass].data.xdim,
				baseMLO->wsum_model.BPref[iclass].data.ydim,
				baseMLO->wsum_model.BPref[iclass].data.zdim,
				baseMLO->wsum_model.BPref[iclass].data.yinit,
				baseMLO->wsum_model.BPref[iclass].data.zinit,
				baseMLO->wsum_model.BPref[iclass].r_max,
				baseMLO->wsum_model.BPref[iclass].padding_factor);

		cudaBackprojectors[iclass].initMdl();

		//If doing predefined projector plan at all and is this class significant
		if (!generateProjectionPlanOnTheFly && baseMLO->mymodel.pdf_class[iclass] > 0.)
		{
			std::vector<int> exp_pointer_dir_nonzeroprior;
			std::vector<int> exp_pointer_psi_nonzeroprior;
			std::vector<double> exp_directions_prior;
			std::vector<double> exp_psi_prior;

			long unsigned itrans_max = baseMLO->sampling.NrTranslationalSamplings() - 1;
			long unsigned nr_idir = baseMLO->sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
			long unsigned nr_ipsi = baseMLO->sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior );

			coarseProjectionPlans[iclass].setup(
					baseMLO->sampling,
					exp_directions_prior,
					exp_psi_prior,
					exp_pointer_dir_nonzeroprior,
					exp_pointer_psi_nonzeroprior,
					NULL, //Mcoarse_significant
					baseMLO->mymodel.pdf_class,
					baseMLO->mymodel.pdf_direction,
					nr_idir,
					nr_ipsi,
					0, //idir_min
					nr_idir - 1, //idir_max
					0, //ipsi_min
					nr_ipsi - 1, //ipsi_max
					0, //itrans_min
					itrans_max,
					0, //current_oversampling
					1, //nr_oversampled_rot
					iclass,
					true, //coarse
					!IS_NOT_INV,
					baseMLO->do_skip_align,
					baseMLO->do_skip_rotate,
					baseMLO->mymodel.orientational_prior_mode
					);
		}
	}

	/*======================================================
	                  TRANSLATIONS SETUP
	======================================================*/

	if (baseMLO->do_shifts_onthefly)
	{
		if (baseMLO->global_fftshifts_ab_coarse.size() > 0)
			translator_coarse1.setShifters(baseMLO->global_fftshifts_ab_coarse);
		else
			translator_coarse1.clear();

		if (baseMLO->global_fftshifts_ab2_coarse.size() > 0)
			translator_coarse2.setShifters(baseMLO->global_fftshifts_ab2_coarse);
		else
			translator_coarse2.clear();

		if (baseMLO->global_fftshifts_ab_current.size() > 0)
			translator_current1.setShifters(baseMLO->global_fftshifts_ab_current);
		else
			translator_current1.clear();

		if (baseMLO->global_fftshifts_ab2_current.size() > 0)
			translator_current2.setShifters(baseMLO->global_fftshifts_ab2_current);
		else
			translator_current2.clear();
	}
};

void MlOptimiserCuda::doThreadExpectationSomeParticles()
{

//	CUDA_CPU_TOC("interParticle");
	CUDA_CPU_TIC("oneTask");
	DEBUG_HANDLE_ERROR(cudaSetDevice(device_id));
	//std::cerr << " calling on device " << device_id << std::endl;
	//put mweight allocation here
	size_t first_ipart = 0, last_ipart = 0;

	while (baseMLO->exp_ipart_ThreadTaskDistributor->getTasks(first_ipart, last_ipart))
	{
		for (long unsigned ipart = first_ipart; ipart <= last_ipart; ipart++)
		{
			CUDA_CPU_TIC("oneParticle");

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
			CUDA_CPU_TIC("getFourierTransformsAndCtfs");
			getFourierTransformsAndCtfs(my_ori_particle, op, sp, baseMLO, this);
			CUDA_CPU_TOC("getFourierTransformsAndCtfs");
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

			/// -- This is a iframe-indexed vector, each entry of which is a dense data-array. These are replacements to using
			//    Mweight in the sparse (Fine-sampled) pass, coarse is unused but created empty input for convert ( FIXME )
			std::vector <IndexedDataArray> CoarsePassWeights(1, allocator) ,FinePassWeights(sp.nr_particles, allocator);
			// -- This is a iframe-indexed vector, each entry of which is a class-indexed vector of masks, one for each
			//    class in FinePassWeights
			std::vector < std::vector <IndexedDataArrayMask> > FinePassClassMasks(sp.nr_particles, std::vector <IndexedDataArrayMask>(baseMLO->mymodel.nr_classes, allocator));
			// -- This is a iframe-indexed vector, each entry of which is parameters used in the projection-operations *after* the
			//    coarse pass, declared here to keep scope to storeWS
			std::vector < ProjectionParams > FineProjectionData(sp.nr_particles, baseMLO->mymodel.nr_classes);

			std::vector < cudaStager<unsigned long> > stagerD2(sp.nr_particles,allocator), stagerSWS(sp.nr_particles,allocator);

			for (int ipass = 0; ipass < nr_sampling_passes; ipass++)
			{
				CUDA_CPU_TIC("weightPass");

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

				if (ipass == 0)
				{
					 //If particle specific sampling plan required
					if (generateProjectionPlanOnTheFly)
					{
						CUDA_CPU_TIC("generateProjectionSetupCoarse");

						for (int iclass = sp.iclass_min; iclass <= sp.iclass_max; iclass++)
						{
							if (baseMLO->mymodel.pdf_class[iclass] > 0.)
							{
								coarseProjectionPlans[iclass].setup(
										baseMLO->sampling,
										op.directions_prior,
										op.psi_prior,
										op.pointer_dir_nonzeroprior,
										op.pointer_psi_nonzeroprior,
										NULL, //Mcoarse_significant
										baseMLO->mymodel.pdf_class,
										baseMLO->mymodel.pdf_direction,
										sp.nr_dir,
										sp.nr_psi,
										sp.idir_min,
										sp.idir_max,
										sp.ipsi_min,
										sp.ipsi_max,
										sp.itrans_min,
										sp.itrans_max,
										0, //current_oversampling
										1, //nr_oversampled_rot
										iclass,
										true, //coarse
										!IS_NOT_INV,
										baseMLO->do_skip_align,
										baseMLO->do_skip_rotate,
										baseMLO->mymodel.orientational_prior_mode
										);
							}
							else
								coarseProjectionPlans[iclass].clear();
						}
						CUDA_CPU_TOC("generateProjectionSetupCoarse");
					}

					CUDA_CPU_TIC("getAllSquaredDifferencesCoarse");
					getAllSquaredDifferencesCoarse(ipass, op, sp, baseMLO, this, stagerD2);
					CUDA_CPU_TOC("getAllSquaredDifferencesCoarse");

					CUDA_CPU_TIC("convertAllSquaredDifferencesToWeightsCoarse");
					convertAllSquaredDifferencesToWeights(ipass, op, sp, baseMLO, this, CoarsePassWeights, FinePassClassMasks, stagerD2);
					CUDA_CPU_TOC("convertAllSquaredDifferencesToWeightsCoarse");
				}
				else
				{
//					// -- go through all classes and generate projectionsetups for all classes - to be used in getASDF and storeWS below --
//					// the reason to do this globally is subtle - we want the orientation_num of all classes to estimate a largest possible
//					// weight-array, which would be insanely much larger than necessary if we had to assume the worst.
					for (long int iframe = 0; iframe < sp.nr_particles; iframe++)
					{
						FineProjectionData[iframe].orientationNumAllClasses = 0;
						for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
						{
							if(exp_iclass>0)
								FineProjectionData[iframe].class_idx[exp_iclass] = FineProjectionData[iframe].rots.size();
							FineProjectionData[iframe].class_entries[exp_iclass] = 0;

							CUDA_CPU_TIC("generateProjectionSetup");
							FineProjectionData[iframe].orientationNumAllClasses += generateProjectionSetup(
									op,
									sp,
									baseMLO,
									false, //not coarse
									exp_iclass,
									FineProjectionData[iframe]);
							CUDA_CPU_TOC("generateProjectionSetup");

						}
						//set a maximum possible size for all weights (to be reduced by significance-checks)
						FinePassWeights[iframe].setDataSize(FineProjectionData[iframe].orientationNumAllClasses*sp.nr_trans*sp.nr_oversampled_trans);
						FinePassWeights[iframe].dual_alloc_all();
						stagerD2[iframe].size= 2*(FineProjectionData[iframe].orientationNumAllClasses*sp.nr_trans*sp.nr_oversampled_trans);
						stagerD2[iframe].prepare();
					}

//					printf("Allocator used space before 'getAllSquaredDifferencesFine': %.2f MiB\n", (float)allocator->getTotalUsedSpace()/(1024.*1024.));

					CUDA_CPU_TIC("getAllSquaredDifferencesFine");
					getAllSquaredDifferencesFine(ipass, op, sp, baseMLO, this, FinePassWeights, FinePassClassMasks, FineProjectionData, stagerD2);
					CUDA_CPU_TOC("getAllSquaredDifferencesFine");


					CUDA_CPU_TIC("convertAllSquaredDifferencesToWeightsFine");
					convertAllSquaredDifferencesToWeights(ipass, op, sp, baseMLO, this, FinePassWeights, FinePassClassMasks, stagerD2);
					CUDA_CPU_TOC("convertAllSquaredDifferencesToWeightsFine");

				}

				CUDA_CPU_TOC("weightPass");
			}

			// For the reconstruction step use mymodel.current_size!
			sp.current_image_size = baseMLO->mymodel.current_size;
			for (long int iframe = 0; iframe < sp.nr_particles; iframe++)
			{
				stagerSWS[iframe].size= 2*(FineProjectionData[iframe].orientationNumAllClasses);
				stagerSWS[iframe].prepare();
			}
			CUDA_CPU_TIC("storeWeightedSums");
			storeWeightedSums(op, sp, baseMLO, this, FinePassWeights, FineProjectionData, FinePassClassMasks, stagerSWS);
			CUDA_CPU_TOC("storeWeightedSums");


//			CUDA_CPU_TIC("Freefalse");
//			for (long int iframe = 0; iframe < sp.nr_particles; iframe++)
//			{
//				stagerD2[iframe].AllData.h_do_free=false;
//				stagerD2[iframe].AllData.d_do_free=false;
//			}
//			CUDA_CPU_TOC("Freefalse");
			CUDA_CPU_TOC("oneParticle");
		}
	}
	CUDA_CPU_TOC("oneTask");
//	CUDA_CPU_TIC("interParticle");
//	exit(0);
}

