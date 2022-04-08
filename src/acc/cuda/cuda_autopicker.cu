#undef ALTCPU
#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cuda_runtime.h>
#include <signal.h>

#include "src/ml_optimiser.h"
#include "src/acc/acc_ptr.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#include "src/acc/cuda/cuda_mem_utils.h"
#include "src/acc/cuda/cuda_settings.h"
#include "src/acc/cuda/cuda_benchmark_utils.h"
#include "src/acc/cuda/cuda_fft.h"

#include "src/macros.h"
#include "src/error.h"

#ifdef CUDA_FORCESTL
#include "src/acc/cuda/cuda_utils_stl.cuh"
#else
#include "src/acc/cuda/cuda_utils_cub.cuh"
#endif

#include "src/acc/utilities.h"
#include "src/acc/acc_helper_functions.h"

#include "src/acc/cuda/cuda_autopicker.h"


AutoPickerCuda::AutoPickerCuda(AutoPicker *basePicker, const char * timing_fnm) :
	node(NULL),
	basePckr(basePicker),
	allocator(new CudaCustomAllocator(0, 1)),
	micTransformer(0, allocator),
	cudaTransformer1(0, allocator),
#ifdef TIMING_FILES
	timer(timing_fnm),
#endif
	cudaTransformer2(0, allocator)

{
	projectors.resize(basePckr->Mrefs.size());
	have_warned_batching=false;
	/*======================================================
	                    DEVICE SETTINGS
	======================================================*/
	device_id = basePicker->device_id;
	int devCount;
	HANDLE_ERROR(cudaGetDeviceCount(&devCount));

	if(device_id >= devCount)
	{
		//std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		CRITICAL(ERR_GPUID);
	}
	else
		HANDLE_ERROR(cudaSetDevice(device_id));
};

AutoPickerCuda::AutoPickerCuda(AutoPickerMpi *basePicker, const char * timing_fnm) :
	basePckr(basePicker),
	allocator(new CudaCustomAllocator(0, 1)),
	micTransformer(0, allocator),
	cudaTransformer1(0, allocator),
#ifdef TIMING_FILES
	timer(timing_fnm),
#endif
	cudaTransformer2(0, allocator)

{
	node = basePicker->getNode();
	basePicker->verb = (node->isLeader()) ? 1 : 0;

	projectors.resize(basePckr->Mrefs.size());
	have_warned_batching=false;
	/*======================================================
	                    DEVICE SETTINGS
	======================================================*/
	device_id = basePicker->device_id;
	int devCount;
	HANDLE_ERROR(cudaGetDeviceCount(&devCount));

	if(device_id >= devCount)
	{
		//std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		CRITICAL(ERR_GPUID);
	}
	else
		HANDLE_ERROR(cudaSetDevice(device_id));
};

void AutoPickerCuda::run()
{
	long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
	if(node!=NULL)
	{
		// Each node does part of the work
		divide_equally(basePckr->fn_micrographs.size(), node->size, node->rank, my_first_micrograph, my_last_micrograph);
	}
	else
	{
		my_first_micrograph = 0;
		my_last_micrograph = basePckr->fn_micrographs.size() - 1;
	}
	my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;

	int barstep;
	if (basePckr->verb > 0)
	{
		std::cout << " Autopicking ..." << std::endl;
		init_progress_bar(my_nr_micrographs);
		barstep = XMIPP_MAX(1, my_nr_micrographs / 60);
	}

	if (!basePckr->do_read_fom_maps)
	{
		CTIC(timer,"setupProjectors");
		for (int iref = 0; iref < (basePckr->Mrefs.size()); iref++)
		{
			projectors[iref].setMdlDim(
							basePckr->PPref[iref].data.xdim,
							basePckr->PPref[iref].data.ydim,
							basePckr->PPref[iref].data.zdim,
							basePckr->PPref[iref].data.yinit,
							basePckr->PPref[iref].data.zinit,
							basePckr->PPref[iref].r_max,
							basePckr->PPref[iref].padding_factor);
			projectors[iref].initMdl(&(basePckr->PPref[iref].data.data[0]));
		}
		CTOC(timer,"setupProjectors");
	}

	FileName fn_olddir="";

	for (long int imic = my_first_micrograph; imic <= my_last_micrograph; imic++)
	{
		if (basePckr->verb > 0 && imic % barstep == 0)
			progress_bar(imic);


		// Check new-style outputdirectory exists and make it if not!
		FileName fn_dir = basePckr->getOutputRootName(basePckr->fn_micrographs[imic]);
		fn_dir = fn_dir.beforeLastOf("/");
		if (fn_dir != fn_olddir)
		{
			// Make a Particles directory
			mktree(fn_dir);
			fn_olddir = fn_dir;
		}
#ifdef TIMING
		basePckr->timer.tic(basePckr->TIMING_A5);
#endif
		autoPickOneMicrograph(basePckr->fn_micrographs[imic], imic);
	}
#ifdef TIMING
		basePckr->timer.toc(basePckr->TIMING_A5);
#endif
	if (basePckr->verb > 0)
		progress_bar(my_nr_micrographs);

	cudaDeviceReset();

}

void AutoPickerCuda::calculateStddevAndMeanUnderMask(AccPtr< ACCCOMPLEX > &d_Fmic,
        AccPtr< ACCCOMPLEX > &d_Fmic2,
        AccPtr< ACCCOMPLEX > &d_Fmsk,
		int nr_nonzero_pixels_mask, AccPtr< XFLOAT > &d_Mstddev,
        AccPtr< XFLOAT > &d_Mmean,
		size_t x, size_t y, size_t mic_size, size_t workSize)
{
	cudaTransformer2.setSize(workSize,workSize,1);

	deviceInitValue(d_Mstddev, (XFLOAT)0.);

	RFLOAT normfft = (RFLOAT)(mic_size * mic_size) / (RFLOAT)nr_nonzero_pixels_mask;

	AccPtr< ACCCOMPLEX > d_Fcov = d_Fmic.make< ACCCOMPLEX >();
	d_Fcov.deviceAlloc(d_Fmic.getSize());

	CTIC(timer,"PRE-multi_0");
	int Bsize( (int) ceilf(( float)d_Fmic.getSize()/(float)BLOCK_SIZE));
	cuda_kernel_convol_B<<<Bsize,BLOCK_SIZE>>>(   ~d_Fmic,
												  ~d_Fmsk,
												  ~d_Fcov,
												  d_Fmic.getSize());
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
	CTOC(timer,"PRE-multi_0");

	CTIC(timer,"PRE-window_0");
	windowFourierTransform2(
			d_Fcov,
			cudaTransformer2.fouriers,
			x, y, 1,
			workSize/2+1, workSize, 1);
	CTOC(timer,"PRE-window_0");

	CTIC(timer,"PRE-Transform_0");
	cudaTransformer2.backward();
	CTOC(timer,"PRE-Transform_0");

	Bsize = ( (int) ceilf(( float)cudaTransformer2.reals.getSize()/(float)BLOCK_SIZE));
	cuda_kernel_multi<XFLOAT><<<Bsize,BLOCK_SIZE>>>(
			~cudaTransformer2.reals,
			~cudaTransformer2.reals,
			(XFLOAT) normfft,
			cudaTransformer2.reals.getSize());
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
	CTIC(timer,"PRE-multi_1");
	cuda_kernel_multi<XFLOAT><<<Bsize,BLOCK_SIZE>>>(
			~cudaTransformer2.reals,
			~cudaTransformer2.reals,
			~d_Mstddev,
			(XFLOAT) -1,
			cudaTransformer2.reals.getSize());
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
	CTOC(timer,"PRE-multi_1");

	CTIC(timer,"PRE-CenterFFT_0");
	runCenterFFT(cudaTransformer2.reals,
				 (int)cudaTransformer2.xSize,
				 (int)cudaTransformer2.ySize,
				 false,
				 1);
	CTOC(timer,"PRE-CenterFFT_0");

	cudaTransformer2.reals.cpOnAcc(d_Mmean); //TODO remove the need for this

	CTIC(timer,"PRE-multi_2");
	Bsize = ( (int) ceilf(( float)d_Fmsk.getSize()/(float)BLOCK_SIZE));
	cuda_kernel_convol_A<<<Bsize,BLOCK_SIZE>>>( 	  ~d_Fmsk,
													  ~d_Fmic2,
													  ~d_Fcov,
													  d_Fmsk.getSize());
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
	CTOC(timer,"PRE-multi_2");


	CTIC(timer,"PRE-window_1");
	windowFourierTransform2(
			d_Fcov,
			cudaTransformer2.fouriers,
			x, y, 1,
			workSize/2+1, workSize, 1);
	CTOC(timer,"PRE-window_1");


	CTIC(timer,"PRE-Transform_1");
	cudaTransformer2.backward();
	CTOC(timer,"PRE-Transform_1");

	CTIC(timer,"PRE-multi_3");
	Bsize = ( (int) ceilf(( float)d_Mstddev.getSize()/(float)BLOCK_SIZE));
		cuda_kernel_finalizeMstddev<<<Bsize,BLOCK_SIZE>>>(
		~d_Mstddev,
		~cudaTransformer2.reals,
		normfft,
		d_Mstddev.getSize());
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
	CTOC(timer,"PRE-multi_3");

	CTIC(timer,"PRE-CenterFFT_1");
	runCenterFFT(d_Mstddev,
				 (int)workSize,
				 (int)workSize,
				 false,
				 1);
	CTOC(timer,"PRE-CenterFFT_1");

}

void AutoPickerCuda::autoPickOneMicrograph(FileName &fn_mic, long int imic)
{
	Image<RFLOAT> Imic;
	MultidimArray<Complex > Faux, Faux2, Fmic;
	MultidimArray<RFLOAT> Maux, Mstddev, Mmean, Mstddev2, Mavg, Mccf_best, Mpsi_best, Fctf, Mccf_best_combined, Mpsi_best_combined;
	MultidimArray<int> Mclass_best_combined;

	AccPtr<XFLOAT >  d_Mccf_best(basePckr->workSize*basePckr->workSize, allocator);
	AccPtr<XFLOAT >  d_Mpsi_best(basePckr->workSize*basePckr->workSize, allocator);
	d_Mccf_best.deviceAlloc();
	d_Mpsi_best.deviceAlloc();

	// Always use the same random seed
	init_random_generator(basePckr->random_seed + imic);

	RFLOAT sum_ref_under_circ_mask, sum_ref2_under_circ_mask;
	int my_skip_side = basePckr->autopick_skip_side + basePckr->particle_size/2;
	CTF ctf;

	int Npsi = 360 / basePckr->psi_sampling;

	int min_distance_pix = ROUND(basePckr->min_particle_distance / basePckr->angpix);
	XFLOAT scale = (XFLOAT)basePckr->workSize / (XFLOAT)basePckr->micrograph_size;

	// Read in the micrograph
#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_A6);
#endif
	CTIC(timer,"readMicrograph");
	Imic.read(fn_mic);
	CTOC(timer,"readMicrograph");
	CTIC(timer,"setXmippOrigin_0");
	Imic().setXmippOrigin();
	CTOC(timer,"setXmippOrigin_0");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_A6);
#endif

	// Let's just check the square size again....
	RFLOAT my_size, my_xsize, my_ysize;
	my_xsize = XSIZE(Imic());
	my_ysize = YSIZE(Imic());
	my_size = (my_xsize != my_ysize) ? XMIPP_MAX(my_xsize, my_ysize) : my_xsize;
	if (basePckr->extra_padding > 0)
		my_size += 2 * basePckr->extra_padding;

	if (my_size != basePckr->micrograph_size || my_xsize != basePckr->micrograph_xsize || my_ysize != basePckr->micrograph_ysize)
	{
		Imic().printShape();
		std::cerr << " micrograph_size= " << basePckr->micrograph_size << " micrograph_xsize= " << basePckr->micrograph_xsize << " micrograph_ysize= " << basePckr->micrograph_ysize << std::endl;
		REPORT_ERROR("AutoPicker::autoPickOneMicrograph ERROR: No differently sized micrographs are allowed in one run, sorry you will have to run separately for each size...");
	}

	if(!basePckr->do_read_fom_maps)
	{
		CTIC(timer,"setSize_micTr");
		micTransformer.setSize(basePckr->micrograph_size, basePckr->micrograph_size, 1,1);
		CTOC(timer,"setSize_micTr");

		CTIC(timer,"setSize_cudaTr");
		cudaTransformer1.setSize(basePckr->workSize,basePckr->workSize, 1, Npsi, FFTW_BACKWARD);
		CTOC(timer,"setSize_cudaTr");
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
    CTIC(timer,"computeStats");
	Imic().computeStats(avg0, stddev0, minval0, maxval0);
    CTOC(timer,"computeStats");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_A7);
#endif
    CTIC(timer,"middlePassFilter");
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imic())
	{
		// Remove pixel values that are too far away from the mean
		if ( ABS(DIRECT_MULTIDIM_ELEM(Imic(), n) - avg0) / stddev0 > basePckr->outlier_removal_zscore)
			DIRECT_MULTIDIM_ELEM(Imic(), n) = avg0;

		DIRECT_MULTIDIM_ELEM(Imic(), n) = (DIRECT_MULTIDIM_ELEM(Imic(), n) - avg0) / stddev0;
	}
    CTOC(timer,"middlePassFilter");

	if (basePckr->micrograph_xsize != basePckr->micrograph_size || basePckr->micrograph_ysize != basePckr->micrograph_size)
	{
		CTIC(timer,"rewindow");
		// Window non-square micrographs to be a square with the largest side
		rewindow(Imic, basePckr->micrograph_size);
		CTOC(timer,"rewindow");
		CTIC(timer,"gaussNoiseOutside");
		// Fill region outside the original window with white Gaussian noise to prevent all-zeros in Mstddev
		FOR_ALL_ELEMENTS_IN_ARRAY2D(Imic())
		{
			if (i < FIRST_XMIPP_INDEX(basePckr->micrograph_ysize)
					|| i > LAST_XMIPP_INDEX(basePckr->micrograph_ysize)
					|| j < FIRST_XMIPP_INDEX(basePckr->micrograph_xsize)
					|| j > LAST_XMIPP_INDEX(basePckr->micrograph_xsize) )
				A2D_ELEM(Imic(), i, j) = rnd_gaus(0.,1.);
		}
		CTOC(timer,"gaussNoiseOutside");
	}

#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_A8);
#endif
	CTIC(timer,"CTFread");
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
				ctf.readByGroup(basePckr->MDmic, &basePckr->obsModel);
				Fctf.resize(basePckr->workSize,basePckr->workSize/2+1);
				ctf.getFftwImage(Fctf, basePckr->micrograph_size, basePckr->micrograph_size, basePckr->angpix, false, false, basePckr->intact_ctf_first_peak, true);
				break;
			}
		}
	}
	CTOC(timer,"CTFread");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_A8);
#endif
#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_A9);
#endif
	CTIC(timer,"mccfResize");
	Mccf_best.resize(basePckr->workSize,basePckr->workSize);
	CTOC(timer,"mccfResize");
	CTIC(timer,"mpsiResize");
	Mpsi_best.resize(basePckr->workSize,basePckr->workSize);
	CTOC(timer,"mpsiResize");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_A9);
#endif
	AccPtr< ACCCOMPLEX > d_Fmic(allocator);
	AccPtr<XFLOAT > d_Mmean(allocator);
	AccPtr<XFLOAT > d_Mstddev(allocator);

#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B1);
#endif
	RFLOAT normfft = (RFLOAT)(basePckr->micrograph_size*basePckr->micrograph_size) / (RFLOAT)basePckr->nr_pixels_circular_mask;;
	if (basePckr->do_read_fom_maps)
	{
		CTIC(timer,"readFromFomMaps_0");
		FileName fn_tmp=basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_stddevNoise.spi";
		Image<RFLOAT> It;
		It.read(fn_tmp);
		if (basePckr->autopick_helical_segments)
			Mstddev2 = It();
		else
			Mstddev = It();
		fn_tmp=basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_avgNoise.spi";
		It.read(fn_tmp);
		if (basePckr->autopick_helical_segments)
			Mavg = It();
		else
			Mmean = It();
		CTOC(timer,"readFromFomMaps_0");
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

		CTIC(timer,"Imic_insert");
		for(int i = 0; i< Imic().nzyxdim ; i++)
			micTransformer.reals[i] = (XFLOAT) Imic().data[i];
		micTransformer.reals.cpToDevice();
		CTOC(timer,"Imic_insert");


		CTIC(timer,"runCenterFFT_0");
		runCenterFFT(micTransformer.reals, micTransformer.xSize, micTransformer.ySize, true, 1);
		CTOC(timer,"runCenterFFT_0");


		CTIC(timer,"FourierTransform_0");
		micTransformer.forward();
		int FMultiBsize = ( (int) ceilf(( float)micTransformer.fouriers.getSize()*2/(float)BLOCK_SIZE));
		CudaKernels::cuda_kernel_multi<XFLOAT><<<FMultiBsize,BLOCK_SIZE>>>(
				(XFLOAT*)~micTransformer.fouriers,
				(XFLOAT)1/((XFLOAT)(micTransformer.reals.getSize())),
				micTransformer.fouriers.getSize()*2);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
		CTOC(timer,"FourierTransform_0");

		if (basePckr->highpass > 0.)
		{
			CTIC(timer,"highpass");
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
			CTOC(timer,"highpass");
		}

		CTIC(timer,"F_cp");
		AccPtr< ACCCOMPLEX > Ftmp(allocator);
		Ftmp.setSize(micTransformer.fouriers.getSize());
		Ftmp.deviceAlloc();
		micTransformer.fouriers.cpOnAcc(Ftmp);
		CTOC(timer,"F_cp");

		// Also calculate the FFT of the squared micrograph
		CTIC(timer,"SquareImic");

		cuda_kernel_square<<<FMultiBsize,BLOCK_SIZE>>>(
				~micTransformer.reals,
				micTransformer.reals.getSize());
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
		CTOC(timer,"SquareImic");

		CTIC(timer,"FourierTransform_1");

		micTransformer.forward();
		CudaKernels::cuda_kernel_multi<XFLOAT><<<FMultiBsize,BLOCK_SIZE>>>(
				(XFLOAT*)~micTransformer.fouriers,
				(XFLOAT)1/((XFLOAT)(micTransformer.reals.getSize())),
				micTransformer.fouriers.getSize()*2);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
		CTOC(timer,"FourierTransform_1");

		// The following calculate mu and sig under the solvent area at every position in the micrograph
		CTIC(timer,"calculateStddevAndMeanUnderMask");

		d_Mstddev.deviceAlloc(basePckr->workSize*basePckr->workSize);
		d_Mmean.deviceAlloc(basePckr->workSize*basePckr->workSize);

		if (basePckr->autopick_helical_segments)
		{

			AccPtr< ACCCOMPLEX > d_Fmsk2(basePckr->Favgmsk.nzyxdim, allocator);
			AccPtr<XFLOAT > d_Mavg(allocator);
			AccPtr<XFLOAT > d_Mstddev2(allocator);

			d_Fmsk2.deviceAlloc();
			d_Mavg.deviceAlloc(basePckr->workSize*basePckr->workSize);
			d_Mstddev2.deviceAlloc(basePckr->workSize*basePckr->workSize);

			//TODO Do this only once further up in scope
			for(int i = 0; i< d_Fmsk2.getSize() ; i++)
			{
				d_Fmsk2[i].x = basePckr->Favgmsk.data[i].real;
				d_Fmsk2[i].y = basePckr->Favgmsk.data[i].imag;
			}
			d_Fmsk2.cpToDevice();
			d_Fmsk2.streamSync();

			calculateStddevAndMeanUnderMask(Ftmp, micTransformer.fouriers, d_Fmsk2, basePckr->nr_pixels_avg_mask, d_Mstddev2, d_Mavg, micTransformer.xFSize, micTransformer.yFSize, basePckr->micrograph_size, basePckr->workSize);

			d_Mstddev2.hostAlloc();
			d_Mstddev2.cpToHost();
			d_Mstddev2.streamSync();
			Mstddev2.resizeNoCp(1, basePckr->workSize, basePckr->workSize);
			for(int i = 0; i < d_Mstddev2.getSize() ; i ++)
				Mstddev2.data[i] = d_Mstddev2[i];

			d_Mavg.hostAlloc();
			d_Mavg.cpToHost();
			d_Mavg.streamSync();
			Mavg.resizeNoCp(1, basePckr->workSize, basePckr->workSize);
			for(int i = 0; i < d_Mavg.getSize() ; i ++)
				Mavg.data[i] = d_Mavg[i];

		}

		//TODO Do this only once further up in scope
		AccPtr< ACCCOMPLEX > d_Fmsk(basePckr->Finvmsk.nzyxdim, allocator);
		d_Fmsk.deviceAlloc();
		for(int i = 0; i< d_Fmsk.getSize() ; i++)
		{
			d_Fmsk[i].x = basePckr->Finvmsk.data[i].real;
			d_Fmsk[i].y = basePckr->Finvmsk.data[i].imag;
		}
		d_Fmsk.cpToDevice();
		d_Fmsk.streamSync();

		calculateStddevAndMeanUnderMask(Ftmp, micTransformer.fouriers, d_Fmsk, basePckr->nr_pixels_circular_invmask, d_Mstddev, d_Mmean, micTransformer.xFSize, micTransformer.yFSize, basePckr->micrograph_size, basePckr->workSize);


		//TODO remove this
		d_Mstddev.hostAlloc();
		d_Mstddev.cpToHost();
		d_Mstddev.streamSync();

		Mstddev.resizeNoCp(1, basePckr->workSize, basePckr->workSize);

		//TODO put this in a kernel
		for(int i = 0; i < d_Mstddev.getSize() ; i ++)
		{
			Mstddev.data[i] = d_Mstddev[i];
			if (d_Mstddev[i] > (XFLOAT)1E-10)
				d_Mstddev[i] = 1 / d_Mstddev[i];
			else
				d_Mstddev[i] = 1;
		}

		d_Mstddev.cpToDevice();
		d_Mstddev.streamSync();

		d_Mmean.hostAlloc();
		d_Mmean.cpToHost();
		d_Mmean.streamSync();
		Mmean.resizeNoCp(1, basePckr->workSize, basePckr->workSize);
		for(int i = 0; i < d_Mmean.getSize() ; i ++)
			Mmean.data[i] = d_Mmean[i];

		CTOC(timer,"calculateStddevAndMeanUnderMask");

		// From now on use downsized Fmic, as the cross-correlation with the references can be done at lower resolution
		CTIC(timer,"windowFourierTransform_0");

		d_Fmic.setSize((basePckr->workSize/2+1)*(basePckr->workSize));
		d_Fmic.deviceAlloc();
		windowFourierTransform2(
				Ftmp,
				d_Fmic,
				basePckr->micrograph_size/2+1, basePckr->micrograph_size, 1, //Input dimensions
				basePckr->workSize/2+1, basePckr->workSize, 1  //Output dimensions
				);
		CTOC(timer,"windowFourierTransform_0");

		if (basePckr->do_write_fom_maps)
		{
			CTIC(timer,"writeToFomMaps");
			// TMP output
			FileName fn_tmp=basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_stddevNoise.spi";
			Image<RFLOAT> It;
			It() = (basePckr->autopick_helical_segments) ? Mstddev2 : Mstddev;
			It.write(fn_tmp);
			fn_tmp=basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_avgNoise.spi";
			It() = (basePckr->autopick_helical_segments) ? Mavg : Mmean;
			It.write(fn_tmp);
			CTOC(timer,"writeToFomMaps");
		}

	}// end if do_read_fom_maps

	// Now start looking for the peaks of all references
	// Clear the output vector with all peaks
	CTIC(timer,"initPeaks");
	std::vector<Peak> peaks;
	peaks.clear();
	CTOC(timer,"initPeaks");
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

			if (basePckr->do_amyloid)
			{
				fn_tmp = basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_combinedPSI.spi";
				It_float.read(fn_tmp);
				Mpsi_best_combined = It_float();
			}
			else
			{
				fn_tmp = basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_combinedCLASS.spi";
				It_int.read(fn_tmp);
				Mclass_best_combined = It_int();
			}
		}
		else
		{
			Mccf_best_combined.clear();
			Mccf_best_combined.resize(basePckr->workSize, basePckr->workSize);
			Mccf_best_combined.initConstant(-99.e99);
			Mclass_best_combined.clear();
			Mclass_best_combined.resize(basePckr->workSize, basePckr->workSize);
			Mclass_best_combined.initConstant(-1);
			Mpsi_best_combined.clear();
			Mpsi_best_combined.resize(basePckr->workSize, basePckr->workSize);
			Mpsi_best_combined.initConstant(-99.e99);
		}
	}

	AccPtr< XFLOAT >  d_ctf(Fctf.nzyxdim, allocator);
	d_ctf.deviceAlloc();
	if(basePckr->do_ctf)
	{
		for(int i = 0; i< d_ctf.getSize() ; i++)
			d_ctf[i]=Fctf.data[i];
		d_ctf.cpToDevice();
	}

	for (int iref = 0; iref < basePckr->Mrefs.size(); iref++)
	{

		CTIC(timer,"OneReference");
		RFLOAT expected_Pratio; // the expectedFOM for this (ctf-corrected) reference
		if (basePckr->do_read_fom_maps)
		{
#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B2);
#endif
			if (!basePckr->autopick_helical_segments)
			{
				CTIC(timer,"readFromFomMaps");
				FileName fn_tmp;
				Image<RFLOAT> It;

				fn_tmp.compose(basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_ref", iref,"_bestCCF.spi");
				It.read(fn_tmp);
				Mccf_best = It();
				It.MDMainHeader.getValue(EMDL_IMAGE_STATS_MAX, expected_Pratio);  // Retrieve expected_Pratio from the header of the image

				fn_tmp.compose(basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_ref", iref,"_bestPSI.spi");
				It.read(fn_tmp);
				Mpsi_best = It();
				CTOC(timer,"readFromFomMaps");
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
			CTIC(timer,"mccfInit");
			deviceInitValue(d_Mccf_best, (XFLOAT)-LARGE_NUMBER);
			CTOC(timer,"mccfInit");
			AccProjectorKernel projKernel = AccProjectorKernel::makeKernel(
									projectors[iref],
									(int)basePckr->workSize/2+1,
									(int)basePckr->workSize,
									1, // Zdim, always 1 in autopicker.
									(int)basePckr->workSize/2+1 -1 );

			int FauxStride = (basePckr->workSize/2+1)*basePckr->workSize;

#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B4);
#endif
	CTIC(timer,"SingleProjection");
	dim3 blocks((int)ceilf((float)FauxStride/(float)BLOCK_SIZE),1);
	if(basePckr->do_ctf)
	{
		cuda_kernel_rotateAndCtf<<<blocks,BLOCK_SIZE>>>(
													  ~cudaTransformer1.fouriers,
													  ~d_ctf,
													  0,
													  projKernel,
													  0
												);
	}
	else
	{
		cuda_kernel_rotateOnly<<<blocks,BLOCK_SIZE>>>(
													  ~cudaTransformer1.fouriers,
													  0,
													  projKernel,
													  0
												);
	}
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
	CTOC(timer,"SingleProjection");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_B4);
#endif
			/*
			 *    FIRST PSI WAS USED FOR PREP CALCS - THIS IS NOW A DEDICATED SECTION
			 *    -------------------------------------------------------------------
			 */

			CTIC(timer,"PREP_CALCS");

#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B5);
#endif
			// Sjors 20April2016: The calculation for sum_ref_under_circ_mask, etc below needs to be done on original micrograph_size!
			CTIC(timer,"windowFourierTransform_FP");
			windowFourierTransform2(cudaTransformer1.fouriers,
									micTransformer.fouriers,
									basePckr->workSize/2+1,        basePckr->workSize,        1, //Input dimensions
									basePckr->micrograph_size/2+1, basePckr->micrograph_size, 1  //Output dimensions
									);
			CTOC(timer,"windowFourierTransform_FP");

			CTIC(timer,"inverseFourierTransform_FP");
			micTransformer.backward();
			CTOC(timer,"inverseFourierTransform_FP");

			CTIC(timer,"runCenterFFT_FP");
			runCenterFFT(micTransformer.reals,
						 (int)micTransformer.xSize,
						 (int)micTransformer.ySize,
						 false,
						 1);
			CTOC(timer,"runCenterFFT_FP");

			micTransformer.reals.cpToHost();

			Maux.resizeNoCp(1,basePckr->micrograph_size, basePckr->micrograph_size);

			micTransformer.reals.streamSync();
			for (int i = 0; i < micTransformer.reals.getSize() ; i ++)
				Maux.data[i] = micTransformer.reals[i];

			CTIC(timer,"setXmippOrigin_FP_0");
			Maux.setXmippOrigin();
			CTOC(timer,"setXmippOrigin_FP_0");
			// TODO: check whether I need CenterFFT(Maux, false)
			// Sjors 20apr2016: checked, somehow not needed.

			sum_ref_under_circ_mask = 0.;
			sum_ref2_under_circ_mask = 0.;
			RFLOAT suma2 = 0.;
			RFLOAT sumn = 1.;
			MultidimArray<RFLOAT> Mctfref(basePckr->particle_size, basePckr->particle_size);
			CTIC(timer,"setXmippOrigin_FP_1");
			Mctfref.setXmippOrigin();
			CTOC(timer,"setXmippOrigin_FP_1");
			CTIC(timer,"suma_FP");
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

			CTOC(timer,"suma_FP");
			CTOC(timer,"PREP_CALCS");

			// for all batches
			CTIC(timer,"AllPsi");
			int startPsi(0);
			for (int psiIter = 0; psiIter < cudaTransformer1.batchIters; psiIter++) // psi-batches for possible memory-limits
			{

				CTIC(timer,"Projection");
				dim3 blocks((int)ceilf((float)FauxStride/(float)BLOCK_SIZE),cudaTransformer1.batchSize[psiIter]);
				if(basePckr->do_ctf)
				{
					cuda_kernel_rotateAndCtf<<<blocks,BLOCK_SIZE>>>(
															  ~cudaTransformer1.fouriers,
															  ~d_ctf,
															  DEG2RAD(basePckr->psi_sampling),
															  projKernel,
															  startPsi
															);
				}
				else
				{
					cuda_kernel_rotateOnly<<<blocks,BLOCK_SIZE>>>(
															  ~cudaTransformer1.fouriers,
															  DEG2RAD(basePckr->psi_sampling),
															  projKernel,
															  startPsi
															);
				}
				LAUNCH_HANDLE_ERROR(cudaGetLastError());
				CTOC(timer,"Projection");

				// Now multiply template and micrograph to calculate the cross-correlation
				CTIC(timer,"convol");
				dim3 blocks2( (int) ceilf(( float)FauxStride/(float)BLOCK_SIZE),cudaTransformer1.batchSize[psiIter]);
				cuda_kernel_batch_convol_A<<<blocks2,BLOCK_SIZE>>>(
					~cudaTransformer1.fouriers,
					~d_Fmic,
					FauxStride);
				LAUNCH_HANDLE_ERROR(cudaGetLastError());
				CTOC(timer,"convol");

				CTIC(timer,"CudaInverseFourierTransform_1");
				cudaTransformer1.backward();
				HANDLE_ERROR(cudaDeviceSynchronize());
				CTOC(timer,"CudaInverseFourierTransform_1");


				CTIC(timer,"runCenterFFT_1");
				runCenterFFT(cudaTransformer1.reals,
							(int)cudaTransformer1.xSize,
							 (int)cudaTransformer1.ySize,
							 false,
							 cudaTransformer1.batchSize[psiIter]);
				CTOC(timer,"runCenterFFT_1");
				// Calculate ratio of prabilities P(ref)/P(zero)
				// Keep track of the best values and their corresponding iref and psi
				// ------------------------------------------------------------------
				// So now we already had precalculated: Mdiff2 = 1/sig*Sum(X^2) - 2/sig*Sum(X) + mu^2/sig*Sum(1)
				// Still to do (per reference): - 2/sig*Sum(AX) + 2*mu/sig*Sum(A) + Sum(A^2)
				CTIC(timer,"probRatio");
				HANDLE_ERROR(cudaDeviceSynchronize());
				dim3 PR_blocks(ceilf((float)(cudaTransformer1.reals.getSize()/cudaTransformer1.batchSize[psiIter])/(float)PROBRATIO_BLOCK_SIZE));
				cuda_kernel_probRatio<<<PR_blocks,PROBRATIO_BLOCK_SIZE>>>(
						~d_Mccf_best,
						~d_Mpsi_best,
						~cudaTransformer1.reals,
						~d_Mmean,
						~d_Mstddev,
						cudaTransformer1.reals.getSize()/cudaTransformer1.batchSize[0],
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
				CTOC(timer,"probRatio");

			} // end for psi-batches
			CTOC(timer,"AllPsi");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_B6);
#endif
#ifdef TIMING
	basePckr->timer.tic(basePckr->TIMING_B7);
#endif
			CTIC(timer,"output");
			d_Mccf_best.cpToHost();
			d_Mpsi_best.cpToHost();
			d_Mccf_best.streamSync();
			for (int i = 0; i < Mccf_best.nzyxdim; i ++)
			{
				Mccf_best.data[i] = d_Mccf_best[i];
				Mpsi_best.data[i] = d_Mpsi_best[i];
			}
			CTOC(timer,"output");

			if (basePckr->do_write_fom_maps && !basePckr->autopick_helical_segments)
			{
				CTIC(timer,"writeFomMaps");
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
				CTOC(timer,"writeFomMaps");

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
						if (basePckr->do_amyloid)
							DIRECT_MULTIDIM_ELEM(Mpsi_best_combined, n) = DIRECT_MULTIDIM_ELEM(Mpsi_best, n);
						else
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
			CTIC(timer,"setXmippOriginX3");
			Mstddev.setXmippOrigin();
			Mmean.setXmippOrigin();
			Mccf_best.setXmippOrigin();
			Mpsi_best.setXmippOrigin();
			CTOC(timer,"setXmippOriginX3");

			CTIC(timer,"peakSearch");
			basePckr->peakSearch(Mccf_best, Mpsi_best, Mstddev, Mmean, iref, my_skip_side, my_ref_peaks, scale);
			CTOC(timer,"peakSearch");

			CTIC(timer,"peakPrune");
			basePckr->prunePeakClusters(my_ref_peaks, min_distance_pix, scale);
			CTOC(timer,"peakPrune");

			CTIC(timer,"peakInsert");
			// append the peaks of this reference to all the other peaks
			peaks.insert(peaks.end(), my_ref_peaks.begin(), my_ref_peaks.end());
			CTOC(timer,"peakInsert");
			CTOC(timer,"OneReference");
#ifdef TIMING
	basePckr->timer.toc(basePckr->TIMING_B8);
#endif

		}
	} // end for iref

	if (basePckr->autopick_helical_segments)
	{
		if (basePckr->do_write_fom_maps)
		{
			FileName fn_tmp;
			Image<RFLOAT> It_float;
			Image<int> It_int;

			It_float() = Mccf_best_combined;
			fn_tmp = basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_combinedCCF.spi";
			It_float.write(fn_tmp);

			if (basePckr->do_amyloid)
			{
				It_float() = Mpsi_best_combined;
				fn_tmp = basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_combinedPSI.spi";
				It_float.write(fn_tmp);
			}
			else
			{
				It_int() = Mclass_best_combined;
				fn_tmp = basePckr->getOutputRootName(fn_mic) + + "_" + basePckr->fn_out + "_combinedCLASS.spi";
				It_int.write(fn_tmp);
			}
		} // end if do_write_fom_maps

		RFLOAT thres = basePckr->min_fraction_expected_Pratio;
		int peak_r_min = 1;
		std::vector<ccfPeak> ccf_peak_list;
		std::vector<std::vector<ccfPeak> > tube_coord_list, tube_track_list;
		std::vector<RFLOAT> tube_len_list;
		MultidimArray<RFLOAT> Mccfplot;

		Mccf_best_combined.setXmippOrigin();
		Mpsi_best_combined.setXmippOrigin();
		Mstddev2.setXmippOrigin();
		Mavg.setXmippOrigin();
		Mclass_best_combined.setXmippOrigin();
		if (basePckr->do_amyloid)
		{
			basePckr->pickAmyloids(Mccf_best_combined, Mpsi_best_combined, Mstddev2, Mavg, thres, basePckr->amyloid_max_psidiff, fn_mic, basePckr->fn_out,
					(basePckr->helical_tube_diameter / basePckr->angpix), basePckr->autopick_skip_side, scale);
		}
		else
		{
			basePckr->pickCCFPeaks(Mccf_best_combined, Mstddev2, Mavg, Mclass_best_combined, thres, peak_r_min, (basePckr->particle_diameter / basePckr->angpix),
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
		}


		if ((basePckr->do_write_fom_maps || basePckr->do_read_fom_maps) && !basePckr->do_amyloid)
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
		CTIC(timer,"finalPeakPrune");
		basePckr->prunePeakClusters(peaks, min_distance_pix, scale);
		CTOC(timer,"finalPeakPrune");

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
