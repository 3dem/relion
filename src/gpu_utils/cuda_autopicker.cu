#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cuda_runtime.h>
#include <signal.h>
#include "src/projector.h"
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_autopicker.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include "src/gpu_utils/cuda_benchmark_utils.cuh"
#include "src/gpu_utils/cuda_helper_functions.cuh"
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_fft.h"

#include "src/complex.h"

#include "src/image.h"


#ifdef CUDA_FORCESTL
#include "src/gpu_utils/cuda_utils_stl.cuh"
#else
#include "src/gpu_utils/cuda_utils_cub.cuh"
#endif

static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

AutoPickerCuda::AutoPickerCuda(AutoPicker *basePicker, int dev_id) :
	basePckr(basePicker)
{

	cudaProjectors.resize(basePckr->Mrefs.size());
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


	/*======================================================
	                    CUSTOM ALLOCATOR
	======================================================*/

#ifdef CUDA_NO_CUSTOM_ALLOCATION
	printf(" DEBUG: Custom allocator is disabled.\n");
	allocator = new CudaCustomAllocator(0, 1);
#else
	size_t allocationSize(0);

	size_t free, total;
	HANDLE_ERROR(cudaMemGetInfo( &free, &total ));

	if (basePckr->available_gpu_memory > 0)
		allocationSize = basePckr->available_gpu_memory * (1000*1000*1000);
	else
		allocationSize = (float)free * .7;

	if (allocationSize > free)
	{
		printf(" WARNING: Required memory per thread, via \"--gpu_memory_per_thread\", not available on device. (Defaulting to less)\n");
		allocationSize = (float)free * .7; //Lets leave some for other processes for now
	}

	int memAlignmentSize;
	cudaDeviceGetAttribute ( &memAlignmentSize, cudaDevAttrTextureAlignment, dev_id );

	allocator = new CudaCustomAllocator(allocationSize, memAlignmentSize);
#endif
};

void AutoPickerCuda::setupProjectors()
{
	for (int iref = 0; iref < (basePckr->Mrefs.size()); iref++)
	{
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
}

void AutoPickerCuda::run()
{

	int barstep;
	if (basePckr->verb > 0)
	{
		std::cout << " Autopicking ..." << std::endl;
		init_progress_bar(basePckr->fn_micrographs.size());
		barstep = XMIPP_MAX(1,basePckr->fn_micrographs.size() / 60);
	}


	for (long int imic = 0; imic < basePckr->fn_micrographs.size(); imic++)
	{
		if (basePckr->verb > 0 && imic % barstep == 0)
			progress_bar(imic);

		autoPickOneMicrograph(basePckr->fn_micrographs[imic]);
	}

	if (basePckr->verb > 0)
		progress_bar(basePckr->fn_micrographs.size());

	cudaDeviceReset();

}

void AutoPickerCuda::calculateStddevAndMeanUnderMask(const MultidimArray<Complex > &_Fmic, const MultidimArray<Complex > &_Fmic2,
		MultidimArray<Complex > &_Fmsk, int nr_nonzero_pixels_mask, MultidimArray<RFLOAT> &_Mstddev, MultidimArray<RFLOAT> &_Mmean)
{

	CudaGlobalPtr< CUDACOMPLEX >  	d_Fmic(_Fmic.nzyxdim, allocator);
	for(int i = 0; i< d_Fmic.size ; i++)
	{
		d_Fmic[i].x=_Fmic.data[i].real;
		d_Fmic[i].y=_Fmic.data[i].imag;
	}
	d_Fmic.put_on_device();

	CudaGlobalPtr< CUDACOMPLEX >  	d_Fmic2(_Fmic2.nzyxdim, allocator);
	for(int i = 0; i< d_Fmic2.size ; i++)
	{
		d_Fmic2[i].x=_Fmic2.data[i].real;
		d_Fmic2[i].y=_Fmic2.data[i].imag;
	}
	d_Fmic2.put_on_device();

	CudaGlobalPtr< CUDACOMPLEX >  	d_Fmsk(_Fmsk.nzyxdim, allocator);
	for(int i = 0; i< d_Fmsk.size ; i++)
	{
		d_Fmsk[i].x=_Fmsk.data[i].real;
		d_Fmsk[i].y=_Fmsk.data[i].imag;
	}
	d_Fmsk.put_on_device();

	MultidimArray<Complex > Faux, Faux2;
	MultidimArray<RFLOAT> Maux(basePckr->micrograph_size,basePckr->micrograph_size);
	FourierTransformer transformer;

	CudaFFT cudaTransformer(0, allocator);
	cudaTransformer.batchSize = 1;
	cudaTransformer.setSize(basePckr->micrograph_size,basePckr->micrograph_size);

	CudaGlobalPtr< XFLOAT > d_Mstddev(basePckr->micrograph_size*basePckr->micrograph_size, allocator);
	d_Mstddev.device_alloc();
	deviceInitValue(d_Mstddev, (XFLOAT)0.);

	_Mstddev.initZeros(basePckr->micrograph_size, basePckr->micrograph_size);
	RFLOAT normfft = (RFLOAT)(basePckr->micrograph_size * basePckr->micrograph_size) / (RFLOAT)nr_nonzero_pixels_mask;

	// Calculate convolution of micrograph and mask, to get average under mask at all points
	CUDA_CPU_TIC("PRE-resize");
	Faux.resizeNoCp(_Fmic.zdim,_Fmic.ydim,_Fmic.xdim);
	CUDA_CPU_TOC("PRE-resize");
#ifdef DEBUG
	Image<RFLOAT> tt;
#endif

	CUDA_CPU_TIC("PRE-multi_0");
	int Bsize( (int) ceilf(( float)d_Fmic.size/(float)BLOCK_SIZE));
	cuda_kernel_convol_B<<<Bsize,BLOCK_SIZE>>>(   d_Fmic.d_ptr,
												  d_Fmsk.d_ptr,
												  d_Fmic.size);
	CUDA_CPU_TOC("PRE-multi_0");

	CUDA_CPU_TIC("PRE-window_0");
	windowFourierTransform2(d_Fmic,
							cudaTransformer.fouriers,
							_Fmic.xdim, _Fmic.ydim, _Fmic.zdim,
							basePckr->micrograph_size/2+1, basePckr->micrograph_size, 1, 0);
	CUDA_CPU_TOC("PRE-window_0");

	CUDA_CPU_TIC("PRE-Transform_0");
	cudaTransformer.backward();
	CUDA_CPU_TOC("PRE-Transform_0");

	Bsize = ( (int) ceilf(( float)cudaTransformer.reals.size/(float)BLOCK_SIZE));
	cuda_kernel_multi<<<Bsize,BLOCK_SIZE>>>( cudaTransformer.reals.d_ptr,
											 cudaTransformer.reals.d_ptr,
										     (XFLOAT) normfft,
										     cudaTransformer.reals.size);

	CUDA_CPU_TIC("PRE-multi_1");
	cuda_kernel_multi<<<Bsize,BLOCK_SIZE>>>( cudaTransformer.reals.d_ptr,
											 cudaTransformer.reals.d_ptr,
											 d_Mstddev.d_ptr,
											 (XFLOAT) -1,
										     cudaTransformer.reals.size);
	CUDA_CPU_TOC("PRE-multi_1");

	CUDA_CPU_TIC("PRE-CenterFFT_0");
	runCenterFFT(cudaTransformer.reals,
				 (int)cudaTransformer.xSize,
				 (int)cudaTransformer.ySize,
				 false,
				 1);
	CUDA_CPU_TOC("PRE-CenterFFT_0");

	cudaTransformer.reals.host_alloc();
	cudaTransformer.reals.cp_to_host();
	_Mmean.resizeNoCp(1,cudaTransformer.ySize,cudaTransformer.xSize);

	cudaTransformer.reals.streamSync();
	for(int i =0; i< _Mmean.nzyxdim; i++)
		_Mmean.data[i] = cudaTransformer.reals[i];

#ifdef DEBUG
	tt()=Maux;
	CenterFFT(tt(), false);
	tt.write("Mavg_mic.spi");
#endif

	CUDA_CPU_TIC("PRE-multi_2");
	Bsize = ( (int) ceilf(( float)d_Fmsk.size/(float)BLOCK_SIZE));
	cuda_kernel_convol_A<<<Bsize,BLOCK_SIZE>>>( 	  d_Fmsk.d_ptr,
													  d_Fmic2.d_ptr,
													  d_Fmsk.size);
	CUDA_CPU_TOC("PRE-multi_2");


	CUDA_CPU_TIC("PRE-window_1");
	windowFourierTransform2(d_Fmsk,
								cudaTransformer.fouriers,
								_Fmsk.xdim, _Fmsk.ydim, _Fmsk.zdim,
								basePckr->micrograph_size/2+1, basePckr->micrograph_size, 1, 0);
	CUDA_CPU_TOC("PRE-window_1");


	CUDA_CPU_TIC("PRE-Transform_1");
	cudaTransformer.backward();
	CUDA_CPU_TOC("PRE-Transform_1");

	CUDA_CPU_TIC("PRE-multi_3");
	Bsize = ( (int) ceilf(( float)d_Mstddev.size/(float)BLOCK_SIZE));
	cuda_kernel_finalizeMstddev<<<Bsize,BLOCK_SIZE>>>( 	  d_Mstddev.d_ptr,
														  cudaTransformer.reals.d_ptr,
														  normfft,
														  d_Mstddev.size);
	CUDA_CPU_TOC("PRE-multi_3");

	CUDA_CPU_TIC("PRE-CenterFFT_1");
	runCenterFFT(d_Mstddev,
				 (int)cudaTransformer.xSize,
				 (int)cudaTransformer.ySize,
				 false,
				 1);
	CUDA_CPU_TOC("PRE-CenterFFT_1");

	d_Mstddev.cp_to_host();
	d_Mstddev.streamSync();
	for(int i =0; i< d_Mstddev.size; i++)
		_Mstddev.data[i] = d_Mstddev[i];

#ifdef DEBUG
	tt()=_Mstddev;
	tt.write("Msig_mic.spi");
#endif


}

void AutoPickerCuda::autoPickOneMicrograph(FileName &fn_mic)
{
	std::cerr << " cudaAutoPicker being run!" << std::endl;
	Image<RFLOAT> Imic;
	MultidimArray<Complex > Faux, Faux2, Fmic;
	MultidimArray<RFLOAT> Maux, Mstddev, Mmean, Mdiff2, MsumX2, Mccf_best, Mpsi_best, Fctf;
	CudaGlobalPtr<XFLOAT > d_Maux;

	CudaGlobalPtr<XFLOAT >  d_Mccf_best(basePckr->micrograph_size*basePckr->micrograph_size, allocator);
	CudaGlobalPtr<XFLOAT >  d_Mpsi_best(basePckr->micrograph_size*basePckr->micrograph_size, allocator);
	d_Mccf_best.host_alloc();
	d_Mpsi_best.host_alloc();
	d_Mccf_best.device_alloc();
	d_Mpsi_best.device_alloc();

	FourierTransformer transformer;
	RFLOAT sum_ref_under_circ_mask, sum_ref2_under_circ_mask;
	int my_skip_side = basePckr->autopick_skip_side + basePckr->particle_size/2;
	CTF ctf;

	//smaller transformer (separate plan & batch-size) for the first psi prep-stuff
	CudaFFT FPcudaTransformer(0, allocator);
	//General transformer for the rest
	CudaFFT cudaTransformer(0, allocator);

	int min_distance_pix = ROUND(basePckr->min_particle_distance / basePckr->angpix);

	// Read in the micrograph
	CUDA_CPU_TIC("readMicrograph");
	Imic.read(fn_mic);
	CUDA_CPU_TOC("readMicrograph");
	CUDA_CPU_TIC("setXmippOrigin_0");
	Imic().setXmippOrigin();
	CUDA_CPU_TOC("setXmippOrigin_0");

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

	// Set mean to zero and stddev to 1 to prevent numerical problems with one-sweep stddev calculations....
    RFLOAT avg0, stddev0, minval0, maxval0;
    CUDA_CPU_TIC("computeStats");
	Imic().computeStats(avg0, stddev0, minval0, maxval0);
    CUDA_CPU_TOC("computeStats");

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
				Fctf.resize(basePckr->downsize_mic, basePckr->downsize_mic/2 + 1);
				ctf.getFftwImage(Fctf, basePckr->micrograph_size, basePckr->micrograph_size, basePckr->angpix, false, false, basePckr->intact_ctf_first_peak, true);
				break;
			}
		}
	}
	CUDA_CPU_TOC("CTFread");

	CUDA_CPU_TIC("mccfResize");
	Mccf_best.resize(basePckr->micrograph_size, basePckr->micrograph_size);
	CUDA_CPU_TOC("mccfResize");
	CUDA_CPU_TIC("mpsifResize");
	Mpsi_best.resize(basePckr->micrograph_size, basePckr->micrograph_size);
	CUDA_CPU_TOC("mpsiResize");

	RFLOAT normfft = (RFLOAT)(basePckr->micrograph_size * basePckr->micrograph_size) / (RFLOAT)basePckr->nr_pixels_circular_mask;;
	if (basePckr->do_read_fom_maps)
	{
		CUDA_CPU_TIC("readFromFomMaps_0");
		FileName fn_tmp=fn_mic.withoutExtension()+"_"+basePckr->fn_out+"_stddevNoise.spi";
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
		CUDA_CPU_TIC("runCenterFFT_0");
		// Fourier Transform (and downscale) Imic()
		runCenterFFT(Imic(), true, allocator);
		CUDA_CPU_TOC("runCenterFFT_0");


		CUDA_CPU_TIC("FourierTransform_0");
		transformer.FourierTransform(Imic(), Fmic);
		CUDA_CPU_TOC("FourierTransform_0");


		// Also calculate the FFT of the squared micrograph
		CUDA_CPU_TIC("MauxResize");
		Maux.resize(Imic());
		CUDA_CPU_TOC("MauxResize");
		CUDA_CPU_TIC("SquareImic");


		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Maux)
		{
			DIRECT_MULTIDIM_ELEM(Maux, n) = DIRECT_MULTIDIM_ELEM(Imic(), n) * DIRECT_MULTIDIM_ELEM(Imic(), n);
		}


		CUDA_CPU_TOC("SquareImic");
		MultidimArray<Complex > Fmic2;
		CUDA_CPU_TIC("FourierTransform_1");
		transformer.FourierTransform(Maux, Fmic2);
		CUDA_CPU_TOC("FourierTransform_1");



		// The following calculate mu and sig under the solvent area at every position in the micrograph

		CUDA_CPU_TIC("calculateStddevAndMeanUnderMask");
		calculateStddevAndMeanUnderMask(Fmic, Fmic2, basePckr->Finvmsk,basePckr->nr_pixels_circular_invmask, Mstddev, Mmean);
		CUDA_CPU_TOC("calculateStddevAndMeanUnderMask");

		if (basePckr->do_write_fom_maps)
		{
			CUDA_CPU_TIC("writeToFomMaps");
			// TMP output
			FileName fn_tmp=fn_mic.withoutExtension()+"_"+basePckr->fn_out+"_stddevNoise.spi";
			Image<RFLOAT> It;
			It() = Mstddev;
			It.write(fn_tmp);
			CUDA_CPU_TOC("writeToFomMaps");
		}

		// From now on use downsized Fmic, as the cross-correlation with the references can be done at lower resolution
		CUDA_CPU_TIC("windowFourierTransform_0");
		windowFourierTransform(Fmic, Faux, basePckr->downsize_mic);
		CUDA_CPU_TOC("windowFourierTransform_0");
		Fmic = Faux;

	}// end if do_read_fom_maps

	CudaGlobalPtr<XFLOAT >  d_Mmean(Mmean.nzyxdim, allocator);
	CudaGlobalPtr<XFLOAT >  d_Mstddev(Mstddev.nzyxdim, allocator);

	for(int i = 0; i < d_Mmean.size; i++)
	{
		d_Mmean[i] = Mmean.data[i];
		d_Mstddev[i] = Mstddev.data[i];
	}
	d_Mmean.put_on_device();
	d_Mstddev.put_on_device();

	// Now start looking for the peaks of all references
	// Clear the output vector with all peaks
	CUDA_CPU_TIC("initPeaks");
	std::vector<Peak> peaks;
	peaks.clear();
	CUDA_CPU_TOC("initPeaks");

	CudaGlobalPtr< CUDACOMPLEX >  	d_Fmic( Fmic.nzyxdim, allocator);
	for(int i = 0; i< d_Fmic.size ; i++)
	{
		d_Fmic[i].x=Fmic.data[i].real;
		d_Fmic[i].y=Fmic.data[i].imag;
	}
	d_Fmic.put_on_device();

	CUDA_CPU_TIC("setupProjectors");
	setupProjectors();
	CUDA_CPU_TOC("setupProjectors");

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
			CUDA_CPU_TIC("readFromFomMaps");
			FileName fn_tmp;
			Image<RFLOAT> It;
			fn_tmp.compose(fn_mic.withoutExtension()+"_"+basePckr->fn_out+"_ref", iref,"_bestCCF.spi");
			It.read(fn_tmp);
			Mccf_best = It();
			// Retrieve expected_Pratio from the header of the image..
			It.MDMainHeader.getValue(EMDL_IMAGE_STATS_MAX, expected_Pratio);
			fn_tmp.compose(fn_mic.withoutExtension()+"_"+basePckr->fn_out+"_ref", iref,"_bestPSI.spi");
			It.read(fn_tmp);
			Mpsi_best = It();
			CUDA_CPU_TOC("readFromFomMaps");

		} //end else if do_read_fom_maps
		else
		{
			CUDA_CPU_TIC("mccfInit");
			deviceInitValue(d_Mccf_best, (XFLOAT)-LARGE_NUMBER);
			CUDA_CPU_TOC("mccfInit");
			CudaProjectorKernel projKernel = CudaProjectorKernel::makeKernel(
									cudaProjectors[iref],
									(int)Faux.xdim,
									(int)Faux.ydim,
									(int)Faux.xdim-1);

			bool is_first_psi = true;

			int Npsi = 360 / basePckr->psi_sampling;
			int Cpsi = 0;
			int FauxStride = Faux.nzyxdim;

			CudaGlobalPtr<CUDACOMPLEX >  d_FauxNpsi(FauxStride, allocator);

			d_FauxNpsi.size=Npsi*FauxStride;
			d_FauxNpsi.device_alloc();

			CUDA_CPU_TIC("Projection");
			dim3 blocks((int)ceilf((float)FauxStride/(float)BLOCK_SIZE),Npsi);
			cuda_kernel_rotateAndCtf<<<blocks,BLOCK_SIZE>>>(
															  &d_FauxNpsi.d_ptr[Cpsi],
															  d_ctf.d_ptr,
															  DEG2RAD(basePckr->psi_sampling),
															  projKernel
														);
			CUDA_CPU_TOC("Projection");

			d_FauxNpsi.streamSync();


			/*
			 *    FIRST PSI WAS USED FOR PREP CALCS - THIS IS NOW A DEDICATED SECTION
			 *    -------------------------------------------------------------------
			 */

			CUDA_CPU_TIC("PREP_CALCS");

			d_FauxNpsi.cp_to_host(FauxStride);
			for (int i = 0; i < FauxStride ; i ++)
			{
				Faux.data[i].real = d_FauxNpsi[i].x;
				Faux.data[i].imag = d_FauxNpsi[i].y;
			}

			FPcudaTransformer.setSize(Maux.xdim, Maux.ydim);
			CUDA_CPU_TIC("windowFourierTransform_FP");
			windowFourierTransform2(d_FauxNpsi,
									FPcudaTransformer.fouriers,
									Faux.xdim, Faux.ydim, Faux.zdim, //Input dimensions
									basePckr->micrograph_size/2+1, basePckr->micrograph_size, 1,  //Output dimensions
									0 //pos for batched array(s) - which this is not
									);
			CUDA_CPU_TOC("windowFourierTransform_FP");

			CUDA_CPU_TIC("inverseFourierTransform_FP");
			FPcudaTransformer.backward();
			CUDA_CPU_TOC("inverseFourierTransform_FP");

			CUDA_CPU_TIC("runCenterFFT_FP");
			runCenterFFT(FPcudaTransformer.reals,
						 (int)FPcudaTransformer.xSize,
						 (int)FPcudaTransformer.ySize,
						 false,
						 1);
			CUDA_CPU_TOC("runCenterFFT_FP");

			FPcudaTransformer.reals.host_alloc();
			FPcudaTransformer.reals.cp_to_host();

			cudaTransformer.reals.streamSync();
			for (int i = 0; i < FPcudaTransformer.reals.size ; i ++)
				Maux.data[i] = FPcudaTransformer.reals[i];

			CUDA_CPU_TIC("setXmippOrigin_FP_0");
			Maux.setXmippOrigin();
			CUDA_CPU_TOC("setXmippOrigin_FP_0");
			// TODO: check whether I need CenterFFT(Maux, false)

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

			for(int i = 0; i < d_Mstddev.size ; i ++)
			{
				if (d_Mstddev[i] > (XFLOAT)1E-10)
					d_Mstddev[i] = 1 / d_Mstddev[i];
				else
					d_Mstddev[i] = 1;

			}
			d_Mstddev.cp_to_device();
			CUDA_CPU_TOC("suma_FP");
			CUDA_CPU_TOC("PREP_CALCS");

			for (RFLOAT psi = 0.; psi < 360.; psi+=basePckr->psi_sampling, Cpsi += FauxStride)
			{
				CUDA_CPU_TIC("OneRotation");
				cudaTransformer.setSize(Maux.xdim, Maux.ydim);

				// Now multiply template and micrograph to calculate the cross-correlation
				CUDA_CPU_TIC("convol");
				dim3 blocks2( (int) ceilf(( float)FauxStride/(float)BLOCK_SIZE),1); // marker - Npsi y-blocks
				cuda_kernel_batch_convol_A<<<blocks2,BLOCK_SIZE>>>(   &d_FauxNpsi.d_ptr[Cpsi],
															  	  	   d_Fmic.d_ptr,
															  	  	   Faux.nzyxdim);
				CUDA_CPU_TOC("convol");


				CUDA_CPU_TIC("windowFourierTransform_1");
				windowFourierTransform2(
						d_FauxNpsi,
						cudaTransformer.fouriers,
						Faux.xdim, Faux.ydim, Faux.zdim, //Input dimensions
						basePckr->micrograph_size/2+1, basePckr->micrograph_size, 1,  //Output dimensions
						(long int) Cpsi // marker
						);
				CUDA_CPU_TOC("windowFourierTransform_1");

				CUDA_CPU_TIC("CudaInverseFourierTransform_1");
				cudaTransformer.backward();

				CUDA_CPU_TIC("runCenterFFT_1");
				runCenterFFT(cudaTransformer.reals,
							 (int)cudaTransformer.xSize,
							 (int)cudaTransformer.ySize,
							 false,
							 1); // marker - Npsi
				CUDA_CPU_TOC("runCenterFFT_1");

//				cudaTransformer.reals.streamSync();

				CUDA_CPU_TOC("CudaInverseFourierTransform_1");

				// Calculate ratio of prabilities P(ref)/P(zero)
				// Keep track of the best values and their corresponding iref and psi
				// ------------------------------------------------------------------
				// So now we already had precalculated: Mdiff2 = 1/sig*Sum(X^2) - 2/sig*Sum(X) + mu^2/sig*Sum(1)
				// Still to do (per reference): - 2/sig*Sum(AX) + 2*mu/sig*Sum(A) + Sum(A^2)
				CUDA_CPU_TIC("probRatio");

				dim3 PR_blocks(ceilf((float)(cudaTransformer.reals.size/(float)PROBRATIO_BLOCK_SIZE)));
				cuda_kernel_probRatio<<<PR_blocks,PROBRATIO_BLOCK_SIZE>>>(
						d_Mccf_best.d_ptr,
						d_Mpsi_best.d_ptr,
						cudaTransformer.reals.d_ptr,
						d_Mmean.d_ptr,
						d_Mstddev.d_ptr,
						cudaTransformer.reals.size,
						(XFLOAT) -2*normfft,
						(XFLOAT) 2*sum_ref_under_circ_mask,
						(XFLOAT) sum_ref2_under_circ_mask,
						(XFLOAT) expected_Pratio,
						1 // marker Npsi
						);

				CUDA_CPU_TOC("probRatio");
			    is_first_psi = false;
			    CUDA_CPU_TOC("OneRotation");
			} // end for psi

			d_Mccf_best.cp_to_host();
			d_Mpsi_best.cp_to_host();

			for (int i = 0; i < Mccf_best.nzyxdim; i ++)
				Mccf_best.data[i] = d_Mccf_best[i];

			d_Mpsi_best.streamSync();

			if (basePckr->do_write_fom_maps)
			{
				CUDA_CPU_TIC("writeFomMaps");
				// TMP output
				FileName fn_tmp;
				Image<RFLOAT> It;
				It() = Mccf_best;
				// Store expected_Pratio in the header of the image..
				It.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, expected_Pratio);;
				fn_tmp.compose(fn_mic.withoutExtension()+"_"+basePckr->fn_out+"_ref", iref,"_bestCCF.spi");
				It.write(fn_tmp);

				It() = Mpsi_best;
				fn_tmp.compose(fn_mic.withoutExtension()+"_"+basePckr->fn_out+"_ref", iref,"_bestPSI.spi");
				It.write(fn_tmp);
				CUDA_CPU_TOC("writeFomMaps");
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mccf_best)
				{
					std::cerr << DIRECT_MULTIDIM_ELEM(Mccf_best, n) << std::endl;
				}
				exit(0);
			} // end if do_write_fom_maps

		} // end if do_read_fom_maps

		// Now that we have Mccf_best and Mpsi_best, get the peaks
		std::vector<Peak> my_ref_peaks;
		CUDA_CPU_TIC("setXmippOriginX3");
		Mstddev.setXmippOrigin();
		Mccf_best.setXmippOrigin();
		Mpsi_best.setXmippOrigin();
		CUDA_CPU_TOC("setXmippOriginX3");

		CUDA_CPU_TIC("peakSearch");
		basePckr->peakSearch(Mccf_best, Mpsi_best, Mstddev, iref, my_skip_side, my_ref_peaks);
		CUDA_CPU_TOC("peakSearch");

		CUDA_CPU_TIC("peakPrune");
		basePckr->prunePeakClusters(my_ref_peaks, min_distance_pix);
		CUDA_CPU_TOC("peakPrune");

		CUDA_CPU_TIC("peakInsert");
		// append the peaks of this reference to all the other peaks
		peaks.insert(peaks.end(), my_ref_peaks.begin(), my_ref_peaks.end());
		CUDA_CPU_TOC("peakInsert");
		CUDA_CPU_TOC("OneReference");
	} // end for iref


	//Now that we have done all references, prune the list again...
	CUDA_CPU_TIC("finalPeakPrune");
	basePckr->prunePeakClusters(peaks, min_distance_pix);
	CUDA_CPU_TOC("finalPeakPrune");

	// And remove all too close neighbours
	basePckr->removeTooCloselyNeighbouringPeaks(peaks, min_distance_pix);

	// Write out a STAR file with the coordinates
	MetaDataTable MDout;
	for (int ipeak =0; ipeak < peaks.size(); ipeak++)
	{
		MDout.addObject();
		MDout.setValue(EMDL_IMAGE_COORD_X, (RFLOAT)(peaks[ipeak].x));
		MDout.setValue(EMDL_IMAGE_COORD_Y, (RFLOAT)(peaks[ipeak].y));
		MDout.setValue(EMDL_ORIENT_PSI, peaks[ipeak].psi);
		MDout.setValue(EMDL_PARTICLE_CLASS, peaks[ipeak].ref + 1); // start counting at 1
		MDout.setValue(EMDL_PARTICLE_AUTOPICK_FOM, peaks[ipeak].fom);
	}
	FileName fn_tmp = fn_mic.withoutExtension() + "_" + basePckr->fn_out + ".star";
	MDout.write(fn_tmp);

}
