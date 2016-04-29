#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cuda_runtime.h>
#include <signal.h>
#include "src/gpu_utils/cuda_autopicker.h"

#include "src/gpu_utils/cuda_mem_utils.h"
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_benchmark_utils.cuh"
#include "src/gpu_utils/cuda_helper_functions.cuh"
#include "src/gpu_utils/cuda_fft.h"


#ifdef CUDA_FORCESTL
#include "src/gpu_utils/cuda_utils_stl.cuh"
#else
#include "src/gpu_utils/cuda_utils_cub.cuh"
#endif


AutoPickerCuda::AutoPickerCuda(AutoPicker *basePicker, int dev_id) :
	basePckr(basePicker),
	allocator(new CudaCustomAllocator(0, 1)),
	micTransformer(0, allocator),
	cudaTransformer1(0, allocator),
	cudaTransformer2(0, allocator)
{

	cudaProjectors.resize(basePckr->Mrefs.size());
	have_warned_batching=false;
	/*======================================================
	                    DEVICE SETTINGS
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
};

void AutoPickerCuda::run()
{

	int barstep;
	if (basePckr->verb > 0)
	{
		std::cout << " Autopicking ..." << std::endl;
		init_progress_bar(basePckr->fn_micrographs.size());
		barstep = XMIPP_MAX(1,basePckr->fn_micrographs.size() / 60);
	}

	if (!basePckr->do_read_fom_maps)
	{

		for (int iref = 0; iref < (basePckr->Mrefs.size()); iref++)
		{
			CUDA_CPU_TIC("setupProjectors");
			cudaProjectors[iref].setMdlDim(
							basePckr->PPref[iref].data.xdim,
							basePckr->PPref[iref].data.ydim,
							basePckr->PPref[iref].data.zdim,
							basePckr->PPref[iref].data.yinit,
							basePckr->PPref[iref].data.zinit,
							basePckr->PPref[iref].r_max,
							basePckr->PPref[iref].padding_factor);
			cudaProjectors[iref].initMdl(&(basePckr->PPref[iref].data.data[0]));
		}
		CUDA_CPU_TOC("setupProjectors");
	}

	FileName fn_olddir="";
	for (long int imic = 0; imic < basePckr->fn_micrographs.size(); imic++)
	{
		if (basePckr->verb > 0 && imic % barstep == 0)
			progress_bar(imic);


		// Check new-style outputdirectory exists and make it if not!
		FileName fn_dir = basePckr->getOutputRootName(basePckr->fn_micrographs[imic]);
		fn_dir = fn_dir.beforeLastOf("/");
		if (fn_dir != fn_olddir)
		{
			// Make a Particles directory
			int res = system(("mkdir -p " + fn_dir).c_str());
			fn_olddir = fn_dir;
		}
#ifdef TIMING
		basePckr->timer.tic(basePckr->TIMING_A5);
#endif
		autoPickOneMicrograph(basePckr->fn_micrographs[imic]);
	}
#ifdef TIMING
		basePckr->timer.toc(basePckr->TIMING_A5);
#endif
	if (basePckr->verb > 0)
		progress_bar(basePckr->fn_micrographs.size());

	cudaDeviceReset();

}

void AutoPickerCuda::calculateStddevAndMeanUnderMask(CudaGlobalPtr< CUDACOMPLEX > &d_Fmic, CudaGlobalPtr< CUDACOMPLEX > &d_Fmic2, CudaGlobalPtr< CUDACOMPLEX > &d_Fmsk,
		int nr_nonzero_pixels_mask, CudaGlobalPtr< XFLOAT > &d_Mstddev, CudaGlobalPtr< XFLOAT > &d_Mmean,
		size_t x, size_t y, size_t mic_size, size_t workSize)
{
	cudaTransformer2.setSize(workSize,workSize);

	deviceInitValue(d_Mstddev, (XFLOAT)0.);

	RFLOAT normfft = (RFLOAT)(mic_size * mic_size) / (RFLOAT)nr_nonzero_pixels_mask;

	CudaGlobalPtr< CUDACOMPLEX > d_Fcov(d_Fmic.getAllocator());
	d_Fcov.device_alloc(d_Fmic.getSize());

	CUDA_CPU_TIC("PRE-multi_0");
	int Bsize( (int) ceilf(( float)d_Fmic.size/(float)BLOCK_SIZE));
	cuda_kernel_convol_B<<<Bsize,BLOCK_SIZE>>>(   ~d_Fmic,
												  ~d_Fmsk,
												  ~d_Fcov,
												  d_Fmic.getSize());
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
	CUDA_CPU_TOC("PRE-multi_0");

	CUDA_CPU_TIC("PRE-window_0");
	windowFourierTransform2(
			d_Fcov,
			cudaTransformer2.fouriers,
			x, y, 1,
			workSize/2+1, workSize, 1);
	CUDA_CPU_TOC("PRE-window_0");

	CUDA_CPU_TIC("PRE-Transform_0");
	cudaTransformer2.backward();
	CUDA_CPU_TOC("PRE-Transform_0");

	Bsize = ( (int) ceilf(( float)cudaTransformer2.reals.size/(float)BLOCK_SIZE));
	cuda_kernel_multi<<<Bsize,BLOCK_SIZE>>>( cudaTransformer2.reals.d_ptr,
											 cudaTransformer2.reals.d_ptr,
										     (XFLOAT) normfft,
										     cudaTransformer2.reals.size);
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
	CUDA_CPU_TIC("PRE-multi_1");
	cuda_kernel_multi<<<Bsize,BLOCK_SIZE>>>( cudaTransformer2.reals.d_ptr,
											 cudaTransformer2.reals.d_ptr,
											 d_Mstddev.d_ptr,
											 (XFLOAT) -1,
										     cudaTransformer2.reals.size);
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
	CUDA_CPU_TOC("PRE-multi_1");

	CUDA_CPU_TIC("PRE-CenterFFT_0");
	runCenterFFT(cudaTransformer2.reals,
				 (int)cudaTransformer2.xSize,
				 (int)cudaTransformer2.ySize,
				 false,
				 1);
	CUDA_CPU_TOC("PRE-CenterFFT_0");

	cudaTransformer2.reals.cp_on_device(d_Mmean); //TODO remove the need for this

	CUDA_CPU_TIC("PRE-multi_2");
	Bsize = ( (int) ceilf(( float)d_Fmsk.size/(float)BLOCK_SIZE));
	cuda_kernel_convol_A<<<Bsize,BLOCK_SIZE>>>( 	  ~d_Fmsk,
													  ~d_Fmic2,
													  ~d_Fcov,
													  d_Fmsk.size);
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
	CUDA_CPU_TOC("PRE-multi_2");


	CUDA_CPU_TIC("PRE-window_1");
	windowFourierTransform2(
			d_Fcov,
			cudaTransformer2.fouriers,
			x, y, 1,
			workSize/2+1, workSize, 1);
	CUDA_CPU_TOC("PRE-window_1");


	CUDA_CPU_TIC("PRE-Transform_1");
	cudaTransformer2.backward();
	CUDA_CPU_TOC("PRE-Transform_1");

	CUDA_CPU_TIC("PRE-multi_3");
	Bsize = ( (int) ceilf(( float)d_Mstddev.size/(float)BLOCK_SIZE));
	cuda_kernel_finalizeMstddev<<<Bsize,BLOCK_SIZE>>>( 	  d_Mstddev.d_ptr,
														  cudaTransformer2.reals.d_ptr,
														  normfft,
														  d_Mstddev.size);
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
	CUDA_CPU_TOC("PRE-multi_3");

	CUDA_CPU_TIC("PRE-CenterFFT_1");
	runCenterFFT(d_Mstddev,
				 (int)workSize,
				 (int)workSize,
				 false,
				 1);
	CUDA_CPU_TOC("PRE-CenterFFT_1");

}

void AutoPickerCuda::autoPickOneMicrograph(FileName &fn_mic)
{
	Image<RFLOAT> Imic;
	MultidimArray<Complex > Faux, Faux2, Fmic;
	MultidimArray<RFLOAT> Maux, Mstddev, Mccf_best, Mpsi_best, Fctf, Mccf_best_combined;
	MultidimArray<int> Mclass_best_combined;

	CudaGlobalPtr<XFLOAT >  d_Mccf_best(basePckr->workSize*basePckr->workSize, allocator);
	CudaGlobalPtr<XFLOAT >  d_Mpsi_best(basePckr->workSize*basePckr->workSize, allocator);
	d_Mccf_best.device_alloc();
	d_Mpsi_best.device_alloc();

	RFLOAT sum_ref_under_circ_mask, sum_ref2_under_circ_mask;
	int my_skip_side = basePckr->autopick_skip_side + basePckr->particle_size/2;
	CTF ctf;

	int Npsi = 360 / basePckr->psi_sampling;

	int min_distance_pix = ROUND(basePckr->min_particle_distance / basePckr->angpix);
	float scale = (float)basePckr->workSize / (float)basePckr->micrograph_size;

	// Read in the micrograph
#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_A6);
#endif
	CUDA_CPU_TIC("readMicrograph");
	Imic.read(fn_mic);
	CUDA_CPU_TOC("readMicrograph");
	CUDA_CPU_TIC("setXmippOrigin_0");
	Imic().setXmippOrigin();
	CUDA_CPU_TOC("setXmippOrigin_0");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_A6);
#endif

	// Let's just check the square size again....
	RFLOAT my_size, my_xsize, my_ysize;
	my_xsize = XSIZE(Imic());
	my_ysize = YSIZE(Imic());
	my_size = (my_xsize != my_ysize) ? XMIPP_MAX(my_xsize, my_ysize) : my_xsize;

	if (my_size != basePckr->micrograph_size || my_xsize != basePckr->micrograph_xsize || my_ysize != basePckr->micrograph_ysize)
	{
		Imic().printShape();
		std::cerr << " micrograph_size= " << basePckr->micrograph_size << " micrograph_xsize= " << basePckr->micrograph_xsize << " micrograph_ysize= " << basePckr->micrograph_ysize << std::endl;
		REPORT_ERROR("AutoPicker::autoPickOneMicrograph ERROR: No differently sized micrographs are allowed in one run, sorry you will have to run separately for each size...");
	}

	if(!basePckr->do_read_fom_maps)
	{
		CUDA_CPU_TIC("setSize_micTr");
		micTransformer.setSize(basePckr->micrograph_size, basePckr->micrograph_size, 1);
		CUDA_CPU_TOC("setSize_micTr");

		CUDA_CPU_TIC("setSize_cudaTr");
		cudaTransformer1.setSize(basePckr->workSize,basePckr->workSize, Npsi, FFTW_BACKWARD);
		CUDA_CPU_TOC("setSize_cudaTr");
	}
	HANDLE_ERROR(cudaDeviceSynchronize());

	if(cudaTransformer1.batchSize.size()>1 && !have_warned_batching)
	{
		have_warned_batching = true;
		std::cerr << std::endl << "*-----------------------------WARNING------------------------------------------------*"<< std::endl;
		std::cerr 			   << "With the current settings the GPU memory is imposing a soft limit on your performace," << std::endl;
		std::cerr 			   << "since one or more micrographs has to use (at least " << cudaTransformer1.batchSize.size() << ") batches of orientations to "<< std::endl;
		std::cerr              << "achieve the total requested " << Npsi << " orientations. Consider using" << std::endl;
		std::cerr			   << "\t higher --ang" << std::endl;
		std::cerr 			   << "\t harder --shrink" << std::endl;
		std::cerr 			   << "\t higher --lowpass with --shrink 0" << std::endl;
		std::cerr              << "*------------------------------------------------------------------------------------*"<< std::endl;
	}

	// Set mean to zero and stddev to 1 to prevent numerical problems with one-sweep stddev calculations....
    RFLOAT avg0, stddev0, minval0, maxval0;
#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_A7);
#endif
    CUDA_CPU_TIC("computeStats");
	Imic().computeStats(avg0, stddev0, minval0, maxval0);
    CUDA_CPU_TOC("computeStats");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_A7);
#endif
    CUDA_CPU_TIC("middlePassFilter");
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imic())
	{
		// Remove pixel values that are too far away from the mean
		if ( ABS(DIRECT_MULTIDIM_ELEM(Imic(), n) - avg0) / stddev0 > basePckr->outlier_removal_zscore)
			DIRECT_MULTIDIM_ELEM(Imic(), n) = avg0;

		DIRECT_MULTIDIM_ELEM(Imic(), n) = (DIRECT_MULTIDIM_ELEM(Imic(), n) - avg0) / stddev0;
	}
    CUDA_CPU_TOC("middlePassFilter");

	if (basePckr->micrograph_xsize !=basePckr->micrograph_ysize)
	{
		CUDA_CPU_TIC("rewindow");
		// Window non-square micrographs to be a square with the largest side
		rewindow(Imic, basePckr->micrograph_size);
		CUDA_CPU_TOC("rewindow");
		CUDA_CPU_TIC("gaussNoiseOutside");
		// Fill region outside the original window with white Gaussian noise to prevent all-zeros in Mstddev
		FOR_ALL_ELEMENTS_IN_ARRAY2D(Imic())
		{
			if (i < FIRST_XMIPP_INDEX(basePckr->micrograph_ysize)
					|| i > LAST_XMIPP_INDEX(basePckr->micrograph_ysize)
					|| j < FIRST_XMIPP_INDEX(basePckr->micrograph_xsize)
					|| j > LAST_XMIPP_INDEX(basePckr->micrograph_xsize) )
				A2D_ELEM(Imic(), i, j) = rnd_gaus(0.,1.);
		}
		CUDA_CPU_TOC("gaussNoiseOutside");
	}

#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_A8);
#endif
	CUDA_CPU_TIC("CTFread");
	// Read in the CTF information if needed
	if (basePckr->do_ctf)
	{
		// Search for this micrograph in the metadata table
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(basePckr->MDmic)
		{
			FileName fn_tmp;
			basePckr->MDmic.getValue(EMDL_MICROGRAPH_NAME, fn_tmp);
			if (fn_tmp==fn_mic)
			{
				ctf.read(basePckr->MDmic, basePckr->MDmic);
				Fctf.resize(basePckr->workSize,basePckr->workSize/2+1);
				ctf.getFftwImage(Fctf, basePckr->micrograph_size, basePckr->micrograph_size, basePckr->angpix, false, false, basePckr->intact_ctf_first_peak, true);
				break;
			}
		}
	}
	CUDA_CPU_TOC("CTFread");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_A8);
#endif
#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_A9);
#endif
	CUDA_CPU_TIC("mccfResize");
	Mccf_best.resize(basePckr->workSize,basePckr->workSize);
	CUDA_CPU_TOC("mccfResize");
	CUDA_CPU_TIC("mpsifResize");
	Mpsi_best.resize(basePckr->workSize,basePckr->workSize);
	CUDA_CPU_TOC("mpsiResize");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_A9);
#endif
	CudaGlobalPtr< CUDACOMPLEX > d_Fmic(allocator);
	CudaGlobalPtr<XFLOAT > d_Mmean(allocator);
	CudaGlobalPtr<XFLOAT > d_Mstddev(allocator);

#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B1);
#endif
	RFLOAT normfft = (RFLOAT)(basePckr->micrograph_size*basePckr->micrograph_size) / (RFLOAT)basePckr->nr_pixels_circular_mask;;
	if (basePckr->do_read_fom_maps)
	{
		CUDA_CPU_TIC("readFromFomMaps_0");
		FileName fn_tmp=basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_stddevNoise.spi";
		Image<RFLOAT> It;
		It.read(fn_tmp);
		Mstddev = It();
		CUDA_CPU_TOC("readFromFomMaps_0");
	}
	else
	{
		/*
		 * Squared difference FOM:
		 * Sum ( (X-mu)/sig  - A )^2 =
		 *  = Sum((X-mu)/sig)^2 - 2 Sum (A*(X-mu)/sig) + Sum(A)^2
		 *  = (1/sig^2)*Sum(X^2) - (2*mu/sig^2)*Sum(X) + (mu^2/sig^2)*Sum(1) - (2/sig)*Sum(AX) + (2*mu/sig)*Sum(A) + Sum(A^2)
		 *
		 * However, the squared difference with an "empty" ie all-zero reference is:
		 * Sum ( (X-mu)/sig)^2
		 *
		 * The ratio of the probabilities thereby becomes:
		 * P(ref) = 1/sqrt(2pi) * exp (( (X-mu)/sig  - A )^2 / -2 )   // assuming sigma = 1!
		 * P(zero) = 1/sqrt(2pi) * exp (( (X-mu)/sig )^2 / -2 )
		 *
		 * P(ref)/P(zero) = exp(( (X-mu)/sig  - A )^2 / -2) / exp ( ( (X-mu)/sig )^2 / -2)
		 *                = exp( (- (2/sig)*Sum(AX) + (2*mu/sig)*Sum(A) + Sum(A^2)) / - 2 )
		 *
		 *                Therefore, I do not need to calculate (X-mu)/sig beforehand!!!
		 *
		 */

		CUDA_CPU_TIC("Imic_insert");
		for(int i = 0; i< Imic().nzyxdim ; i++)
			micTransformer.reals[i] = (XFLOAT) Imic().data[i];
		micTransformer.reals.cp_to_device();
		CUDA_CPU_TIC("Imic_insert");


		CUDA_CPU_TIC("runCenterFFT_0");
		runCenterFFT(micTransformer.reals, micTransformer.xSize, micTransformer.ySize, true, 1);
		CUDA_CPU_TOC("runCenterFFT_0");


		CUDA_CPU_TIC("FourierTransform_0");
		micTransformer.forward();
		int FMultiBsize = ( (int) ceilf(( float)micTransformer.fouriers.getSize()*2/(float)BLOCK_SIZE));
		cuda_kernel_multi<<<FMultiBsize,BLOCK_SIZE>>>(
				(XFLOAT*)~micTransformer.fouriers,
				(XFLOAT)1/((XFLOAT)(micTransformer.reals.getSize())),
				micTransformer.fouriers.getSize()*2);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
		CUDA_CPU_TOC("FourierTransform_0");

		if (basePckr->highpass > 0.)
		{
			CUDA_CPU_TIC("highpass");
			micTransformer.fouriers.streamSync();
			lowPassFilterMapGPU(	micTransformer.fouriers,
									(size_t)1,
									micTransformer.yFSize,
									micTransformer.xFSize,
									XSIZE(Imic()),
									basePckr->lowpass,
									basePckr->highpass,
									basePckr->angpix,
									2,
									true); //false = lowpass, true=highpass
			micTransformer.fouriers.streamSync();
			micTransformer.backward();
			micTransformer.reals.streamSync();
			CUDA_CPU_TOC("highpass");
		}

		CUDA_CPU_TIC("F_cp");
		CudaGlobalPtr< CUDACOMPLEX > Ftmp(allocator);
		Ftmp.setSize(micTransformer.fouriers.getSize());
		Ftmp.device_alloc();
		micTransformer.fouriers.cp_on_device(Ftmp);
		CUDA_CPU_TOC("F_cp");

		// Also calculate the FFT of the squared micrograph
		CUDA_CPU_TIC("SquareImic");

		cuda_kernel_square<<<FMultiBsize,BLOCK_SIZE>>>(
				~micTransformer.reals,
				micTransformer.reals.getSize());
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
		CUDA_CPU_TOC("SquareImic");

		CUDA_CPU_TIC("FourierTransform_1");

		micTransformer.forward();
		cuda_kernel_multi<<<FMultiBsize,BLOCK_SIZE>>>(
				(XFLOAT*)~micTransformer.fouriers,
				(XFLOAT)1/((XFLOAT)(micTransformer.reals.getSize())),
				micTransformer.fouriers.getSize()*2);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
		CUDA_CPU_TOC("FourierTransform_1");

		// The following calculate mu and sig under the solvent area at every position in the micrograph
		CUDA_CPU_TIC("calculateStddevAndMeanUnderMask");

		d_Mstddev.device_alloc(basePckr->workSize*basePckr->workSize);
		d_Mmean.device_alloc(basePckr->workSize*basePckr->workSize);


		//TODO Do this only once further up in scope
		CudaGlobalPtr< CUDACOMPLEX > d_Fmsk(basePckr->Finvmsk.nzyxdim, allocator);
		for(int i = 0; i< d_Fmsk.size ; i++)
		{
			d_Fmsk[i].x = basePckr->Finvmsk.data[i].real;
			d_Fmsk[i].y = basePckr->Finvmsk.data[i].imag;
		}
		d_Fmsk.put_on_device();
		d_Fmsk.streamSync();

		calculateStddevAndMeanUnderMask(Ftmp, micTransformer.fouriers, d_Fmsk, basePckr->nr_pixels_circular_invmask, d_Mstddev, d_Mmean, micTransformer.xFSize, micTransformer.yFSize, basePckr->micrograph_size, basePckr->workSize);


		//TODO remove this
		d_Mstddev.host_alloc();
		d_Mstddev.cp_to_host();
		d_Mstddev.streamSync();

		Mstddev.resizeNoCp(1, basePckr->workSize, basePckr->workSize);

		//TODO put this in a kernel
		for(int i = 0; i < d_Mstddev.size ; i ++)
		{
			Mstddev.data[i] = d_Mstddev[i];
			if (d_Mstddev[i] > (XFLOAT)1E-10)
				d_Mstddev[i] = 1 / d_Mstddev[i];
			else
				d_Mstddev[i] = 1;
		}

		d_Mstddev.cp_to_device();
		d_Mstddev.streamSync();



		CUDA_CPU_TOC("calculateStddevAndMeanUnderMask");

		// From now on use downsized Fmic, as the cross-correlation with the references can be done at lower resolution
		CUDA_CPU_TIC("windowFourierTransform_0");

		d_Fmic.setSize((basePckr->workSize/2+1)*(basePckr->workSize));
		d_Fmic.device_alloc();
		windowFourierTransform2(
				Ftmp,
				d_Fmic,
				basePckr->micrograph_size/2+1, basePckr->micrograph_size, 1, //Input dimensions
				basePckr->workSize/2+1, basePckr->workSize, 1  //Output dimensions
				);
		CUDA_CPU_TOC("windowFourierTransform_0");

		if (basePckr->do_write_fom_maps)
		{
			CUDA_CPU_TIC("writeToFomMaps");
			// TMP output
			FileName fn_tmp=basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_stddevNoise.spi";
			Image<RFLOAT> It;
			It() = Mstddev;
			It.write(fn_tmp);
			CUDA_CPU_TOC("writeToFomMaps");
		}

	}// end if do_read_fom_maps

	// Now start looking for the peaks of all references
	// Clear the output vector with all peaks
	CUDA_CPU_TIC("initPeaks");
	std::vector<Peak> peaks;
	peaks.clear();
	CUDA_CPU_TOC("initPeaks");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_B1);
#endif

	if (basePckr->autopick_helical_segments)
	{
		if (basePckr->do_read_fom_maps)
		{
			FileName fn_tmp;
			Image<RFLOAT> It_float;
			Image<int> It_int;

			fn_tmp = basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_combinedCCF.spi";
			It_float.read(fn_tmp);
			Mccf_best_combined = It_float();

			fn_tmp = basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_combinedCLASS.spi";
			It_int.read(fn_tmp);
			Mclass_best_combined = It_int();
		}
		else
		{
			Mccf_best_combined.clear();
			Mccf_best_combined.resize(basePckr->workSize, basePckr->workSize);
			Mccf_best_combined.initConstant(-99.e99);
			Mclass_best_combined.clear();
			Mclass_best_combined.resize(basePckr->workSize, basePckr->workSize);
			Mclass_best_combined.initConstant(-1);
		}
	}

	CudaGlobalPtr< XFLOAT >  d_ctf(Fctf.nzyxdim, allocator);
	for(int i = 0; i< d_ctf.size ; i++)
		d_ctf[i]=Fctf.data[i];

	d_ctf.put_on_device();

	for (int iref = 0; iref < basePckr->Mrefs.size(); iref++)
	{

		CUDA_CPU_TIC("OneReference");
		RFLOAT expected_Pratio; // the expectedFOM for this (ctf-corrected) reference
		if (basePckr->do_read_fom_maps)
		{
#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B2);
#endif
			if (!basePckr->autopick_helical_segments)
			{
				CUDA_CPU_TIC("readFromFomMaps");
				FileName fn_tmp;
				Image<RFLOAT> It;

				fn_tmp.compose(basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_ref", iref,"_bestCCF.spi");
				It.read(fn_tmp);
				Mccf_best = It();
				It.MDMainHeader.getValue(EMDL_IMAGE_STATS_MAX, expected_Pratio);  // Retrieve expected_Pratio from the header of the image

				fn_tmp.compose(basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_ref", iref,"_bestPSI.spi");
				It.read(fn_tmp);
				Mpsi_best = It();
				CUDA_CPU_TOC("readFromFomMaps");
			}
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_B2);
#endif

		} //end else if do_read_fom_maps
		else
		{
#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B3);
#endif
			CUDA_CPU_TIC("mccfInit");
			deviceInitValue(d_Mccf_best, (XFLOAT)-LARGE_NUMBER);
			CUDA_CPU_TOC("mccfInit");
			CudaProjectorKernel projKernel = CudaProjectorKernel::makeKernel(
									cudaProjectors[iref],
									(int)basePckr->workSize/2+1,
									(int)basePckr->workSize,
									(int)basePckr->workSize/2+1 -1 );

			int FauxStride = (basePckr->workSize/2+1)*basePckr->workSize;

#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B4);
#endif
	CUDA_CPU_TIC("SingleProjection");
	dim3 blocks((int)ceilf((float)FauxStride/(float)BLOCK_SIZE),1);
	cuda_kernel_rotateAndCtf<<<blocks,BLOCK_SIZE>>>(
													  ~cudaTransformer1.fouriers,
													  ~d_ctf,
													  0,
													  projKernel,
													  0
												);
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
	CUDA_CPU_TOC("SingleProjection");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_B4);
#endif
			/*
			 *    FIRST PSI WAS USED FOR PREP CALCS - THIS IS NOW A DEDICATED SECTION
			 *    -------------------------------------------------------------------
			 */

			CUDA_CPU_TIC("PREP_CALCS");

#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B5);
#endif
			// Sjors 20April2016: The calculation for sum_ref_under_circ_mask, etc below needs to be done on original micrograph_size!
			CUDA_CPU_TIC("windowFourierTransform_FP");
			windowFourierTransform2(cudaTransformer1.fouriers,
									micTransformer.fouriers,
									basePckr->workSize/2+1,        basePckr->workSize,        1, //Input dimensions
									basePckr->micrograph_size/2+1, basePckr->micrograph_size, 1  //Output dimensions
									);
			CUDA_CPU_TOC("windowFourierTransform_FP");

			CUDA_CPU_TIC("inverseFourierTransform_FP");
			micTransformer.backward();
			CUDA_CPU_TOC("inverseFourierTransform_FP");

			CUDA_CPU_TIC("runCenterFFT_FP");
			runCenterFFT(micTransformer.reals,
						 (int)micTransformer.xSize,
						 (int)micTransformer.ySize,
						 false,
						 1);
			CUDA_CPU_TOC("runCenterFFT_FP");

			micTransformer.reals.cp_to_host();

			Maux.resizeNoCp(1,basePckr->micrograph_size, basePckr->micrograph_size);

			micTransformer.reals.streamSync();
			for (int i = 0; i < micTransformer.reals.size ; i ++)
				Maux.data[i] = micTransformer.reals[i];

			CUDA_CPU_TIC("setXmippOrigin_FP_0");
			Maux.setXmippOrigin();
			CUDA_CPU_TOC("setXmippOrigin_FP_0");
			// TODO: check whether I need CenterFFT(Maux, false)
			// Sjors 20apr2016: checked, somehow not needed.

			sum_ref_under_circ_mask = 0.;
			sum_ref2_under_circ_mask = 0.;
			RFLOAT suma2 = 0.;
			RFLOAT sumn = 1.;
			MultidimArray<RFLOAT> Mctfref(basePckr->particle_size, basePckr->particle_size);
			CUDA_CPU_TIC("setXmippOrigin_FP_1");
			Mctfref.setXmippOrigin();
			CUDA_CPU_TOC("setXmippOrigin_FP_1");
			CUDA_CPU_TIC("suma_FP");
			FOR_ALL_ELEMENTS_IN_ARRAY2D(Mctfref) // only loop over smaller Mctfref, but take values from large Maux!
			{
				if (i*i + j*j < basePckr->particle_radius2)
				{
					suma2 += A2D_ELEM(Maux, i, j) * A2D_ELEM(Maux, i, j);
					suma2 += 2. * A2D_ELEM(Maux, i, j) * rnd_gaus(0., 1.);
					sum_ref_under_circ_mask += A2D_ELEM(Maux, i, j);
					sum_ref2_under_circ_mask += A2D_ELEM(Maux, i, j) * A2D_ELEM(Maux, i, j);
					sumn += 1.;
				}
			}
			sum_ref_under_circ_mask /= sumn;
			sum_ref2_under_circ_mask /= sumn;
			expected_Pratio = exp(suma2 / (2. * sumn));

			CUDA_CPU_TOC("suma_FP");
			CUDA_CPU_TOC("PREP_CALCS");

			// for all batches
			int startPsi(0);
			for (int psiIter = 0; psiIter < cudaTransformer1.psiIters; psiIter++) // psi-batches for possible memory-limits
			{

				CUDA_CPU_TIC("Projection");
				dim3 blocks((int)ceilf((float)FauxStride/(float)BLOCK_SIZE),cudaTransformer1.batchSize[psiIter]);
				cuda_kernel_rotateAndCtf<<<blocks,BLOCK_SIZE>>>(
																  ~cudaTransformer1.fouriers,
																  ~d_ctf,
																  DEG2RAD(basePckr->psi_sampling),
																  projKernel,
																  startPsi
															);
				LAUNCH_HANDLE_ERROR(cudaGetLastError());
				CUDA_CPU_TOC("Projection");

				// Now multiply template and micrograph to calculate the cross-correlation
				CUDA_CPU_TIC("convol");
				dim3 blocks2( (int) ceilf(( float)FauxStride/(float)BLOCK_SIZE),cudaTransformer1.batchSize[psiIter]);
				cuda_kernel_batch_convol_A<<<blocks2,BLOCK_SIZE>>>(   cudaTransformer1.fouriers.d_ptr,
																	  d_Fmic.d_ptr,
																	  FauxStride);
				LAUNCH_HANDLE_ERROR(cudaGetLastError());
				CUDA_CPU_TOC("convol");

				CUDA_CPU_TIC("CudaInverseFourierTransform_1");
				cudaTransformer1.backward();
				HANDLE_ERROR(cudaDeviceSynchronize());
				CUDA_CPU_TOC("CudaInverseFourierTransform_1");


				CUDA_CPU_TIC("runCenterFFT_1");
				runCenterFFT(cudaTransformer1.reals,
							(int)cudaTransformer1.xSize,
							 (int)cudaTransformer1.ySize,
							 false,
							 cudaTransformer1.batchSize[psiIter]);
				CUDA_CPU_TOC("runCenterFFT_1");
				// Calculate ratio of prabilities P(ref)/P(zero)
				// Keep track of the best values and their corresponding iref and psi
				// ------------------------------------------------------------------
				// So now we already had precalculated: Mdiff2 = 1/sig*Sum(X^2) - 2/sig*Sum(X) + mu^2/sig*Sum(1)
				// Still to do (per reference): - 2/sig*Sum(AX) + 2*mu/sig*Sum(A) + Sum(A^2)
				CUDA_CPU_TIC("probRatio");
				HANDLE_ERROR(cudaDeviceSynchronize());
				dim3 PR_blocks(ceilf((float)(cudaTransformer1.reals.size/cudaTransformer1.batchSize[psiIter])/(float)PROBRATIO_BLOCK_SIZE));
				cuda_kernel_probRatio<<<PR_blocks,PROBRATIO_BLOCK_SIZE>>>(
						d_Mccf_best.d_ptr,
						d_Mpsi_best.d_ptr,
						cudaTransformer1.reals.d_ptr,
						d_Mmean.d_ptr,
						d_Mstddev.d_ptr,
						cudaTransformer1.reals.size/cudaTransformer1.batchSize[0],
						(XFLOAT) -2*normfft,
						(XFLOAT) 2*sum_ref_under_circ_mask,
						(XFLOAT) sum_ref2_under_circ_mask,
						(XFLOAT) expected_Pratio,
						cudaTransformer1.batchSize[psiIter],
						startPsi,
						Npsi
						);
				LAUNCH_HANDLE_ERROR(cudaGetLastError());
				startPsi += cudaTransformer1.batchSize[psiIter];
				CUDA_CPU_TOC("probRatio");
			    CUDA_CPU_TOC("OneRotation");
			} // end for psi-batches

#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_B6);
#endif
#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B7);
#endif
			CUDA_CPU_TIC("output");
			d_Mccf_best.cp_to_host();
			d_Mpsi_best.cp_to_host();
			d_Mccf_best.streamSync();
			for (int i = 0; i < Mccf_best.nzyxdim; i ++)
			{
				Mccf_best.data[i] = d_Mccf_best[i];
				Mpsi_best.data[i] = d_Mpsi_best[i];
			}
			CUDA_CPU_TOC("output");

			if (basePckr->do_write_fom_maps && !basePckr->autopick_helical_segments)
			{
				CUDA_CPU_TIC("writeFomMaps");
				// TMP output
				FileName fn_tmp;
				Image<RFLOAT> It;
				It() = Mccf_best;
				// Store expected_Pratio in the header of the image..
				It.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, expected_Pratio);  // Store expected_Pratio in the header of the image
				fn_tmp.compose(basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_ref", iref,"_bestCCF.spi");
				It.write(fn_tmp);

				It() = Mpsi_best;
				fn_tmp.compose(basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_ref", iref,"_bestPSI.spi");
				It.write(fn_tmp);
				CUDA_CPU_TOC("writeFomMaps");

			} // end if do_write_fom_maps
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_B7);
#endif
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_B3);
#endif
		} // end if do_read_fom_maps


		//TODO FIX HELICAL SEGMENTS SUPPORT
		if (basePckr->autopick_helical_segments)
		{
			if (!basePckr->do_read_fom_maps)
			{
				// Combine Mccf_best and Mpsi_best from all refs
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mccf_best)
				{
					RFLOAT new_ccf = DIRECT_MULTIDIM_ELEM(Mccf_best, n);
					RFLOAT old_ccf = DIRECT_MULTIDIM_ELEM(Mccf_best_combined, n);
					if (new_ccf > old_ccf)
					{
						DIRECT_MULTIDIM_ELEM(Mccf_best_combined, n) = new_ccf;
						DIRECT_MULTIDIM_ELEM(Mclass_best_combined, n) = iref;
					}
				}
			}
		}
		else
		{
#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B8);
#endif
			// Now that we have Mccf_best and Mpsi_best, get the peaks
			std::vector<Peak> my_ref_peaks;
			CUDA_CPU_TIC("setXmippOriginX3");
			Mstddev.setXmippOrigin();
			Mccf_best.setXmippOrigin();
			Mpsi_best.setXmippOrigin();
			CUDA_CPU_TOC("setXmippOriginX3");

			CUDA_CPU_TIC("peakSearch");
			basePckr->peakSearch(Mccf_best, Mpsi_best, Mstddev, iref, my_skip_side, my_ref_peaks, scale);
			CUDA_CPU_TOC("peakSearch");

			CUDA_CPU_TIC("peakPrune");
			basePckr->prunePeakClusters(my_ref_peaks, min_distance_pix, scale);
			CUDA_CPU_TOC("peakPrune");

			CUDA_CPU_TIC("peakInsert");
			// append the peaks of this reference to all the other peaks
			peaks.insert(peaks.end(), my_ref_peaks.begin(), my_ref_peaks.end());
			CUDA_CPU_TOC("peakInsert");
			CUDA_CPU_TOC("OneReference");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_B8);
#endif

		}
	} // end for iref


	if (basePckr->autopick_helical_segments)
	{
		RFLOAT thres = basePckr->min_fraction_expected_Pratio;
		int peak_r_min = 2;
		std::vector<ccfPeak> ccf_peak_list;
		std::vector<std::vector<ccfPeak> > tube_coord_list, tube_track_list;
		std::vector<RFLOAT> tube_len_list;
		MultidimArray<RFLOAT> Mccfplot;

		Mccf_best_combined.setXmippOrigin();
		Mclass_best_combined.setXmippOrigin();
		basePckr->pickCCFPeaks(Mccf_best_combined, Mclass_best_combined, thres, peak_r_min, (basePckr->particle_diameter / basePckr->angpix),
				ccf_peak_list, Mccfplot, my_skip_side, scale);
		basePckr->extractHelicalTubes(ccf_peak_list, tube_coord_list, tube_len_list, tube_track_list,
				(basePckr->particle_diameter / basePckr->angpix), basePckr->helical_tube_curvature_factor_max,
				(basePckr->min_particle_distance / basePckr->angpix), (basePckr->helical_tube_diameter / basePckr->angpix), scale);
		basePckr->exportHelicalTubes(Mccf_best_combined, Mccfplot, Mclass_best_combined,
					tube_coord_list, tube_track_list, tube_len_list,
					fn_mic, basePckr->fn_out,
					(basePckr->particle_diameter / basePckr->angpix),
					(basePckr->helical_tube_length_min / basePckr->angpix),
					my_skip_side, scale);

		if (basePckr->do_write_fom_maps)
		{
			FileName fn_tmp;
			Image<RFLOAT> It_float;
			Image<int> It_int;

			It_float() = Mccf_best_combined;
			fn_tmp = basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_combinedCCF.spi";
			It_float.write(fn_tmp);

			It_int() = Mclass_best_combined;
			fn_tmp = basePckr->getOutputRootName(fn_mic) + + "_" + basePckr->fn_out + "_combinedCLASS.spi";
			It_int.write(fn_tmp);
		} // end if do_write_fom_maps

		if (basePckr->do_write_fom_maps || basePckr->do_read_fom_maps)
		{
			FileName fn_tmp;
			Image<RFLOAT> It;

			It() = Mccfplot;
			fn_tmp =  basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_combinedPLOT.spi";
			It.write(fn_tmp);
		}
	}
	else
	{
#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B9);
#endif
		//Now that we have done all references, prune the list again...
		CUDA_CPU_TIC("finalPeakPrune");
		basePckr->prunePeakClusters(peaks, min_distance_pix, scale);
		CUDA_CPU_TOC("finalPeakPrune");

		// And remove all too close neighbours
		basePckr->removeTooCloselyNeighbouringPeaks(peaks, min_distance_pix, scale);

		// Write out a STAR file with the coordinates
		MetaDataTable MDout;
		for (int ipeak =0; ipeak < peaks.size(); ipeak++)
		{
			MDout.addObject();
			MDout.setValue(EMDL_IMAGE_COORD_X, (RFLOAT)(peaks[ipeak].x)/scale);
			MDout.setValue(EMDL_IMAGE_COORD_Y, (RFLOAT)(peaks[ipeak].y)/scale);
			MDout.setValue(EMDL_PARTICLE_CLASS, peaks[ipeak].ref + 1); // start counting at 1
			MDout.setValue(EMDL_PARTICLE_AUTOPICK_FOM, peaks[ipeak].fom);
			MDout.setValue(EMDL_ORIENT_PSI, peaks[ipeak].psi);
		}
		FileName fn_tmp = basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + ".star";
		MDout.write(fn_tmp);
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_B9);
#endif
	}

}
