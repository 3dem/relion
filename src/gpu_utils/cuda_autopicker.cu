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

void calculateStddevAndMeanUnderMask2(CudaGlobalPtr< CUDACOMPLEX > &d_Fmic, CudaGlobalPtr< CUDACOMPLEX > &d_Fmic2, CudaGlobalPtr< CUDACOMPLEX > &d_Fmsk,
		int nr_nonzero_pixels_mask, CudaGlobalPtr< XFLOAT > &d_Mstddev, CudaGlobalPtr< XFLOAT > &d_Mmean,
		size_t x, size_t y, size_t workSize)
{
	CudaFFT cudaTransformer(0, d_Fmic.getAllocator());
	cudaTransformer.setSize(workSize,workSize);

	deviceInitValue(d_Mstddev, (XFLOAT)0.);

	RFLOAT normfft = (RFLOAT)(workSize*workSize) / (RFLOAT)nr_nonzero_pixels_mask;

	CudaGlobalPtr< CUDACOMPLEX > d_Fcov(d_Fmic.getAllocator());
	d_Fcov.device_alloc(d_Fmic.getSize());

	CUDA_CPU_TIC("PRE-multi_0");
	int Bsize( (int) ceilf(( float)d_Fmic.size/(float)BLOCK_SIZE));
	cuda_kernel_convol_B<<<Bsize,BLOCK_SIZE>>>(   ~d_Fmic,
												  ~d_Fmsk,
												  ~d_Fcov,
												  d_Fmic.getSize());
	CUDA_CPU_TOC("PRE-multi_0");

	CUDA_CPU_TIC("PRE-window_0");
	windowFourierTransform2(
			d_Fcov,
			cudaTransformer.fouriers,
			x, y, 1,
			workSize/2+1, workSize, 1);
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

	cudaTransformer.reals.cp_on_device(d_Mmean); //TODO remove the need for this

	CUDA_CPU_TIC("PRE-multi_2");
	Bsize = ( (int) ceilf(( float)d_Fmsk.size/(float)BLOCK_SIZE));
	cuda_kernel_convol_A<<<Bsize,BLOCK_SIZE>>>( 	  ~d_Fmsk,
													  ~d_Fmic2,
													  ~d_Fcov,
													  d_Fmsk.size);
	CUDA_CPU_TOC("PRE-multi_2");


	CUDA_CPU_TIC("PRE-window_1");
	windowFourierTransform2(
			d_Fcov,
			cudaTransformer.fouriers,
			x, y, 1,
			workSize/2+1, workSize, 1);
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
				 (int)workSize,
				 (int)workSize,
				 false,
				 1);
	CUDA_CPU_TOC("PRE-CenterFFT_1");
}

void AutoPickerCuda::calculateStddevAndMeanUnderMask(const MultidimArray<Complex > &Fmic, const MultidimArray<Complex > &Fmic2,
		int nr_nonzero_pixels_mask, MultidimArray<RFLOAT> &Mstddev, MultidimArray<RFLOAT> &Mmean)
{

	CudaGlobalPtr< CUDACOMPLEX > d_Fmic(Fmic.nzyxdim, allocator);
	for(int i = 0; i< d_Fmic.size ; i++)
	{
		d_Fmic[i].x=Fmic.data[i].real;
		d_Fmic[i].y=Fmic.data[i].imag;
	}
	d_Fmic.put_on_device();

	CudaGlobalPtr< CUDACOMPLEX > d_Fmic2(Fmic2.nzyxdim, allocator);
	for(int i = 0; i< d_Fmic2.size ; i++)
	{
		d_Fmic2[i].x=Fmic2.data[i].real;
		d_Fmic2[i].y=Fmic2.data[i].imag;
	}
	d_Fmic2.put_on_device();

	CudaGlobalPtr< CUDACOMPLEX > d_Fmsk(basePckr->Finvmsk.nzyxdim, allocator);
	for(int i = 0; i< d_Fmsk.size ; i++)
	{
		d_Fmsk[i].x = basePckr->Finvmsk.data[i].real;
		d_Fmsk[i].y = basePckr->Finvmsk.data[i].imag;
	}
	d_Fmsk.put_on_device();

	CudaGlobalPtr< XFLOAT > d_Mstddev(basePckr->micrograph_size*basePckr->micrograph_size, allocator);
	CudaGlobalPtr< XFLOAT > d_Mmean(basePckr->micrograph_size*basePckr->micrograph_size, allocator);
	d_Mstddev.device_alloc();
	d_Mmean.device_alloc();


	calculateStddevAndMeanUnderMask2(d_Fmic, d_Fmic2, d_Fmsk, nr_nonzero_pixels_mask, d_Mstddev, d_Mmean, Fmic.xdim, Fmic.ydim, basePckr->micrograph_size);

	Mmean.resizeNoCp(1, basePckr->micrograph_size, basePckr->micrograph_size);
	d_Mmean.cp_to_host();
	d_Mmean.streamSync();
	for(int i =0; i< Mmean.nzyxdim; i++)
		Mmean.data[i] = d_Mmean[i];

	Mstddev.resizeNoCp(1, basePckr->micrograph_size, basePckr->micrograph_size);
	d_Mstddev.cp_to_host();
	d_Mstddev.streamSync();
	for(int i =0; i< d_Mstddev.size; i++)
		Mstddev.data[i] = d_Mstddev[i];
}

void AutoPickerCuda::autoPickOneMicrograph(FileName &fn_mic)
{
	Image<RFLOAT> Imic;
	MultidimArray<Complex > Faux, Faux2, Fmic;
	MultidimArray<RFLOAT> Maux, Mstddev, Mmean, Mdiff2, MsumX2, Mccf_best, Mpsi_best, Fctf, Mccf_best_combined;
	MultidimArray<int> Mclass_best_combined;
	CudaGlobalPtr<XFLOAT > d_Maux;

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
	CUDA_CPU_TIC("readMicrograph");
	Imic.read(fn_mic);
	CUDA_CPU_TOC("readMicrograph");
	CUDA_CPU_TIC("setXmippOrigin_0");
	Imic().setXmippOrigin();
	CUDA_CPU_TOC("setXmippOrigin_0");

	size_t downsize_Fmic_x = basePckr->downsize_mic / 2 + 1;
	size_t downsize_Fmic_y = basePckr->downsize_mic;

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

	CudaFFT micTransformer(0, allocator);
	CudaFFT extraMicTransformer(0,allocator);
	CudaFFT FPcudaTransformer(0, allocator);
	CudaFFT cudaTransformer(0, allocator);
	CudaFFT extraCudaTransformer(0, allocator);
	if(!basePckr->do_read_fom_maps)
	{
		CUDA_CPU_TIC("setSize_micTr");
		micTransformer.setSize(Imic().xdim, Imic().ydim,1);
		CUDA_CPU_TOC("setSize_micTr");
		CUDA_CPU_TIC("setSize_micTr");
		extraMicTransformer.setSize(Imic().xdim, Imic().ydim,2);
		CUDA_CPU_TOC("setSize_micTr");

		CUDA_CPU_TIC("setSize_FPudaTr");
		FPcudaTransformer.setSize(basePckr->workSize,basePckr->workSize, 1);
		CUDA_CPU_TOC("setSize_FPcudaTr");
		CUDA_CPU_TIC("setSize_cudaTr");
		cudaTransformer.setSize(basePckr->workSize,basePckr->workSize, Npsi);
		CUDA_CPU_TOC("setSize_cudaTr");
		CUDA_CPU_TIC("setSize_extraCudaTr");
		extraCudaTransformer.setSize(basePckr->workSize,basePckr->workSize, 2);
		CUDA_CPU_TOC("setSize_extraCudaTr");
	}
	HANDLE_ERROR(cudaDeviceSynchronize());

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
				Fctf.resize(downsize_Fmic_y, downsize_Fmic_x);
				ctf.getFftwImage(Fctf, basePckr->micrograph_size, basePckr->micrograph_size, basePckr->angpix, false, false, basePckr->intact_ctf_first_peak, true);
				break;
			}
		}
	}
	CUDA_CPU_TOC("CTFread");

	CUDA_CPU_TIC("mccfResize");
	Mccf_best.resize(basePckr->workSize,basePckr->workSize);
	CUDA_CPU_TOC("mccfResize");
	CUDA_CPU_TIC("mpsifResize");
	Mpsi_best.resize(basePckr->workSize,basePckr->workSize);
	CUDA_CPU_TOC("mpsiResize");

	CudaGlobalPtr< CUDACOMPLEX > d_Fmic(allocator);
	CudaGlobalPtr<XFLOAT > d_Mmean(allocator);
	CudaGlobalPtr<XFLOAT > d_Mstddev(allocator);

	RFLOAT normfft = (RFLOAT)(basePckr->workSize*basePckr->workSize) / (RFLOAT)basePckr->nr_pixels_circular_mask;;
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

		//TODO ADD HIGH PASS FILTER
//		if (highpass > 0.)
//        {
//			lowPassFilterMap(Fmic, XSIZE(Imic()), highpass, angpix, 2, true); // true means highpass instead of lowpass!
//        	transformer.inverseFourierTransform(Fmic, Imic()); // also calculate inverse transform again for squared calculation below
//        }



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
		CUDA_CPU_TOC("FourierTransform_0");

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
		CUDA_CPU_TOC("SquareImic");

		CUDA_CPU_TIC("FourierTransform_1");

		micTransformer.forward();
		cuda_kernel_multi<<<FMultiBsize,BLOCK_SIZE>>>(
				(XFLOAT*)~micTransformer.fouriers,
				(XFLOAT)1/((XFLOAT)(micTransformer.reals.getSize())),
				micTransformer.fouriers.getSize()*2);
		CUDA_CPU_TOC("FourierTransform_1");



		// The following calculate mu and sig under the solvent area at every position in the micrograph

		CUDA_CPU_TIC("calculateStddevAndMeanUnderMask");



//
//		//TODO REMOVED THIS
//		MultidimArray<Complex > _Fmic, _Fmic2;
//		_Fmic.resize(micTransformer.yFSize, micTransformer.xFSize);
//		_Fmic2.resize(micTransformer.yFSize, micTransformer.xFSize);
//
//		Ftmp.host_alloc();
//		Ftmp.cp_to_host();
//		micTransformer.fouriers.cp_to_host();
//		micTransformer.fouriers.streamSync();
//		for (int i = 0; i < Ftmp.getSize(); i ++)
//		{
//			_Fmic.data[i].real = Ftmp[i].x;
//			_Fmic.data[i].imag = Ftmp[i].y;
//			_Fmic2.data[i].real = micTransformer.fouriers[i].x;
//			_Fmic2.data[i].imag = micTransformer.fouriers[i].y;
//		}
//
//
//
//
//
//		calculateStddevAndMeanUnderMask(_Fmic, _Fmic2, basePckr->nr_pixels_circular_invmask, Mstddev, Mmean);
////		basePckr->calculateStddevAndMeanUnderMask(_Fmic, _Fmic2, basePckr->Finvmsk, basePckr->nr_pixels_circular_invmask, Mstddev, Mmean);
//
//
//
//
//
//
//		//TODO REMOVED THIS
//		d_Mmean.setSize(Mmean.nzyxdim);
//		d_Mstddev.setSize(Mstddev.nzyxdim);
//		d_Mmean.host_alloc();
//		d_Mstddev.host_alloc();
//		for(int i = 0; i < d_Mmean.size; i++)
//		{
//			d_Mmean[i] = Mmean.data[i];
//			d_Mstddev[i] = Mstddev.data[i];
//		}
//		d_Mmean.put_on_device();
//		d_Mstddev.put_on_device();






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

		calculateStddevAndMeanUnderMask2(Ftmp, micTransformer.fouriers, d_Fmsk, basePckr->nr_pixels_circular_invmask, d_Mstddev, d_Mmean, micTransformer.xFSize, micTransformer.yFSize, basePckr->workSize);


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

		d_Fmic.setSize(downsize_Fmic_x * downsize_Fmic_y);
		d_Fmic.device_alloc();
		windowFourierTransform2(
				Ftmp,
				d_Fmic,
				Imic().xdim/2+1, Imic().ydim, 1, //Input dimensions
				downsize_Fmic_x, downsize_Fmic_y, 1  //Output dimensions
				);
		CUDA_CPU_TOC("windowFourierTransform_0");

		if (basePckr->do_write_fom_maps)
		{
			CUDA_CPU_TIC("writeToFomMaps");
			// TMP output
			FileName fn_tmp=basePckr->getOutputRootName(fn_mic)+"_"+basePckr->fn_out+"_stddevNoise.spi";
			Image<RFLOAT> It;

//			Mstddev.resizeNoCp(1, basePckr->workSize, basePckr->workSize);
//			for(int i =0; i< d_Mstddev.size; i++)
//				Mstddev.data[i] = d_Mstddev[i];

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


	//TODO FIX HELICAL SEGMENTS SUPPORT
//	if (autopick_helical_segments)
//	{
//		if (do_read_fom_maps)
//		{
//			FileName fn_tmp;
//			Image<RFLOAT> It_float;
//			Image<int> It_int;
//
//			fn_tmp = getOutputRootName(fn_mic)+"_"+fn_out+"_combinedCCF.mrc";
//			It_float.read(fn_tmp);
//			Mccf_best_combined = It_float();
//
//			fn_tmp = getOutputRootName(fn_mic)+"_"+fn_out+"_combinedCLASS.mrc";
//			It_int.read(fn_tmp);
//			Mclass_best_combined = It_int();
//		}
//		else
//		{
//			Mccf_best_combined.clear();
//			Mccf_best_combined.resize(micrograph_size, micrograph_size);
//			Mccf_best_combined.initConstant(-99.e99);
//			Mclass_best_combined.clear();
//			Mclass_best_combined.resize(micrograph_size, micrograph_size);
//			Mclass_best_combined.initConstant(-1);
//		}
//	}



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
									(int)downsize_Fmic_x,
									(int)downsize_Fmic_y,
									(int)downsize_Fmic_x-1);

			int FauxStride = downsize_Fmic_x*downsize_Fmic_y;

			CudaGlobalPtr<CUDACOMPLEX >  d_FauxNpsi(allocator);

			d_FauxNpsi.setSize(Npsi*FauxStride);
			d_FauxNpsi.device_alloc();

			CUDA_CPU_TIC("Projection");
			dim3 blocks((int)ceilf((float)FauxStride/(float)BLOCK_SIZE),Npsi);
			cuda_kernel_rotateAndCtf<<<blocks,BLOCK_SIZE>>>(
															  ~d_FauxNpsi,
															  ~d_ctf,
															  DEG2RAD(basePckr->psi_sampling),
															  projKernel
														);
			CUDA_CPU_TOC("Projection");

			/*
			 *    FIRST PSI WAS USED FOR PREP CALCS - THIS IS NOW A DEDICATED SECTION
			 *    -------------------------------------------------------------------
			 */

			CUDA_CPU_TIC("PREP_CALCS");

//			FPcudaTransformer.setSize(basePckr->workSize,basePckr->workSize);
			CUDA_CPU_TIC("windowFourierTransform_FP");
			windowFourierTransform2(d_FauxNpsi,
									FPcudaTransformer.fouriers,
									downsize_Fmic_x, downsize_Fmic_y, 1, //Input dimensions
									basePckr->workSize/2+1, basePckr->workSize, 1  //Output dimensions
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

			FPcudaTransformer.reals.cp_to_host();

			Maux.resizeNoCp(1,basePckr->workSize,basePckr->workSize)	;

			FPcudaTransformer.reals.streamSync();
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

			CUDA_CPU_TOC("suma_FP");
			CUDA_CPU_TOC("PREP_CALCS");

//			for (RFLOAT psi = 0.; psi < 360.; psi+=basePckr->psi_sampling, Cpsi += FauxStride)
//			{
				CUDA_CPU_TIC("OneRotation");

				// Now multiply template and micrograph to calculate the cross-correlation
				CUDA_CPU_TIC("convol");
				dim3 blocks2( (int) ceilf(( float)FauxStride/(float)BLOCK_SIZE),Npsi);
				cuda_kernel_batch_convol_A<<<blocks2,BLOCK_SIZE>>>(   d_FauxNpsi.d_ptr,
															  	  	  d_Fmic.d_ptr,
															  	  	  FauxStride);
				CUDA_CPU_TOC("convol");
				HANDLE_ERROR(cudaDeviceSynchronize());

				CUDA_CPU_TIC("windowFourierTransform_1");
				windowFourierTransform2(
						d_FauxNpsi,
						cudaTransformer.fouriers,
						downsize_Fmic_x, downsize_Fmic_y, 1, //Input dimensions
						basePckr->workSize/2+1, basePckr->workSize, 1,  //Output dimensions
						Npsi
						);
				CUDA_CPU_TOC("windowFourierTransform_1");
				HANDLE_ERROR(cudaDeviceSynchronize());


				CUDA_CPU_TIC("CudaInverseFourierTransform_1");
				cudaTransformer.backward();
				HANDLE_ERROR(cudaDeviceSynchronize());


				CUDA_CPU_TIC("runCenterFFT_1");
				runCenterFFT(cudaTransformer.reals,
							 (int)cudaTransformer.xSize,
							 (int)cudaTransformer.ySize,
							 false,
							 Npsi);
				CUDA_CPU_TOC("runCenterFFT_1");
				HANDLE_ERROR(cudaDeviceSynchronize());

				CUDA_CPU_TOC("CudaInverseFourierTransform_1");

				// Calculate ratio of prabilities P(ref)/P(zero)
				// Keep track of the best values and their corresponding iref and psi
				// ------------------------------------------------------------------
				// So now we already had precalculated: Mdiff2 = 1/sig*Sum(X^2) - 2/sig*Sum(X) + mu^2/sig*Sum(1)
				// Still to do (per reference): - 2/sig*Sum(AX) + 2*mu/sig*Sum(A) + Sum(A^2)
				CUDA_CPU_TIC("probRatio");
				HANDLE_ERROR(cudaDeviceSynchronize());
				dim3 PR_blocks(ceilf((float)(cudaTransformer.reals.size/Npsi)/(float)PROBRATIO_BLOCK_SIZE));
				cuda_kernel_probRatio<<<PR_blocks,PROBRATIO_BLOCK_SIZE>>>(
						d_Mccf_best.d_ptr,
						d_Mpsi_best.d_ptr,
						cudaTransformer.reals.d_ptr,
						d_Mmean.d_ptr,
						d_Mstddev.d_ptr,
						cudaTransformer.reals.size/Npsi,
						(XFLOAT) -2*normfft,
						(XFLOAT) 2*sum_ref_under_circ_mask,
						(XFLOAT) sum_ref2_under_circ_mask,
						(XFLOAT) expected_Pratio,
						Npsi
						);

				CUDA_CPU_TOC("probRatio");
			    CUDA_CPU_TOC("OneRotation");
//			} // end for psi


//			if(basePckr->workSize!=basePckr->micrograph_size) // if we've been working with a smaller copy, resize it back (in fourier space)
//			{
//				CUDA_CPU_TIC("resize_output");
////				CUDA_CPU_TIC("setSize_EMT");
////				extraMicTransformer.setSize(basePckr->micrograph_size,basePckr->micrograph_size, 2);
////				CUDA_CPU_TOC("setSize_EMT");
////				CUDA_CPU_TIC("setSize_CT");
////				cudaTransformer.setSize(basePckr->workSize,basePckr->workSize,2); //set batchSize to 2 (ccf and psi) to avoid excessive transform calcs
////				CUDA_CPU_TOC("setSize_CT");
//				d_Mccf_best.cp_on_device(extraCudaTransformer.reals);
//				d_Mpsi_best.cp_on_device(&(extraCudaTransformer.reals.d_ptr[d_Mccf_best.size]));
//
//				extraCudaTransformer.forward();
//				int FMultiBsize = ( (int) ceilf(( float)(extraCudaTransformer.fouriers.getSize()*2*2)/(float)BLOCK_SIZE));
//				cuda_kernel_multi<<<FMultiBsize,BLOCK_SIZE>>>(
//						(XFLOAT*)~extraCudaTransformer.fouriers,
//						(XFLOAT)1/((XFLOAT)(extraCudaTransformer.reals.getSize())),
//						extraCudaTransformer.fouriers.getSize()*2*2);
//
//				windowFourierTransform2(    extraCudaTransformer.fouriers,
//											extraMicTransformer.fouriers,
//											basePckr->workSize/2+1, basePckr->workSize, 1,     //Input dimensions
//											basePckr->micrograph_size/2+1, basePckr->micrograph_size, 1,  //Output dimensions
//											2
//											);
//				extraMicTransformer.backward();
//
//				CudaGlobalPtr < RFLOAT > d_aux(allocator); // NOTE - RFLOAT (NOT XFLOAT)
//				d_aux.size = basePckr->micrograph_size*basePckr->micrograph_size*2;
//				d_aux.device_alloc();
//
//				FMultiBsize = ( (int) ceilf(( float)(extraMicTransformer.reals.getSize())/(float)BLOCK_SIZE));
//				cuda_kernel_cast<<<FMultiBsize,BLOCK_SIZE>>>(
//						(XFLOAT*)~extraMicTransformer.reals,
//						(RFLOAT*)~d_aux,
//						d_aux.size);
//
////				cudaCpyDeviceToHost(&d_aux.d_ptr[0],                 &Mccf_best.data[0],Mccf_best.nzyxdim, 0);
////				cudaCpyDeviceToHost(&d_aux.d_ptr[Mccf_best.nzyxdim], &Mccf_best.data[0],Mpsi_best.nzyxdim, 0);
//
//				d_aux.cp_to_host(&Mccf_best.data[0],basePckr->micrograph_size*basePckr->micrograph_size);
//				d_aux.streamSync();
//				d_aux.d_ptr=&d_aux.d_ptr[Mccf_best.nzyxdim];
//				d_aux.cp_to_host(&Mpsi_best.data[0],basePckr->micrograph_size*basePckr->micrograph_size);
//				d_aux.streamSync();
//////				HANDLE_ERROR(cudaDeviceSynchronize());
////				extraMicTransformer.reals.cp_to_host();
////				extraMicTransformer.reals.streamSync();
//////
////				CUDA_CPU_TIC("output");
////				for (int i = 0; i < Mccf_best.nzyxdim; i ++)
////				{
////					Mccf_best.data[i] = extraMicTransformer.reals[i];
////					Mpsi_best.data[i] = extraMicTransformer.reals[Mccf_best.nzyxdim + i];
////				}
////				CUDA_CPU_TOC("output");
//				CUDA_CPU_TOC("resize_output");
//			}
//			else // otherwise just get and prepare for further use.
			{
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
			}

			if (basePckr->do_write_fom_maps)
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

//				for (long int n=0; n<((Mccf_best).nzyxdim); n+=10000)
//				{
//					std::cerr << DIRECT_MULTIDIM_ELEM(Mccf_best, n) << std::endl;
//				}
//				exit(0);
			} // end if do_write_fom_maps

		} // end if do_read_fom_maps


		//TODO FIX HELICAL SEGMENTS SUPPORT
//		if (autopick_helical_segments)
//		{
//			if (!do_read_fom_maps)
//			{
//				// Combine Mccf_best and Mpsi_best from all refs
//				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mccf_best)
//				{
//					RFLOAT new_ccf = DIRECT_MULTIDIM_ELEM(Mccf_best, n);
//					RFLOAT old_ccf = DIRECT_MULTIDIM_ELEM(Mccf_best_combined, n);
//					if (new_ccf > old_ccf)
//					{
//						DIRECT_MULTIDIM_ELEM(Mccf_best_combined, n) = new_ccf;
//						DIRECT_MULTIDIM_ELEM(Mclass_best_combined, n) = iref;
//					}
//				}
//			}
//		}
//		else
		{
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

		}
	} // end for iref


	//TODO FIX HELICAL SEGMENTS SUPPORT
//	if (autopick_helical_segments)
//	{
//		RFLOAT thres = min_fraction_expected_Pratio;
//		int peak_r_min = 2;
//		std::vector<ccfPeak> ccf_peak_list;
//		std::vector<std::vector<ccfPeak> > tube_coord_list, tube_track_list;
//		std::vector<RFLOAT> tube_len_list;
//		MultidimArray<RFLOAT> Mccfplot;
//
//		Mccf_best_combined.setXmippOrigin();
//		Mclass_best_combined.setXmippOrigin();
//		pickCCFPeaks(Mccf_best_combined, Mclass_best_combined, thres, peak_r_min, (particle_diameter / angpix), ccf_peak_list, Mccfplot, micrograph_size, micrograph_minxy_size, my_skip_side);
//		extractHelicalTubes(ccf_peak_list, tube_coord_list, tube_len_list, tube_track_list, (particle_diameter / angpix), helical_tube_curvature_factor_max, (min_particle_distance / angpix), (helical_tube_diameter / angpix));
//		exportHelicalTubes(Mccf_best_combined, Mccfplot, Mclass_best_combined,
//					tube_coord_list, tube_track_list, tube_len_list,
//					fn_mic, fn_out,
//					(particle_diameter / angpix),
//					(helical_tube_length_min / angpix),
//					micrograph_size,
//					micrograph_xsize,
//					micrograph_ysize,
//					my_skip_side);
//
//		if (do_write_fom_maps)
//		{
//			FileName fn_tmp;
//			Image<RFLOAT> It_float;
//			Image<int> It_int;
//
//			It_float() = Mccf_best_combined;
//			fn_tmp = getOutputRootName(fn_mic) + "_" + fn_out + "_combinedCCF.mrc";
//			It_float.write(fn_tmp);
//
//			It_int() = Mclass_best_combined;
//			fn_tmp = getOutputRootName(fn_mic) + + "_" + fn_out + "_combinedCLASS.mrc";
//			It_int.write(fn_tmp);
//		} // end if do_write_fom_maps
//
//		if (do_write_fom_maps || do_read_fom_maps)
//		{
//			FileName fn_tmp;
//			Image<RFLOAT> It;
//
//			It() = Mccfplot;
//			fn_tmp =  getOutputRootName(fn_mic) + "_" + fn_out + "_combinedPLOT.mrc";
//			It.write(fn_tmp);
//		}
//	}
//	else
	{
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
	}

}
