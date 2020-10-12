#ifndef ACC_UTILITIES_IMPL_H_
#define ACC_UTILITIES_IMPL_H_

#include "src/acc/acc_ptr.h"
#include "src/acc/data_types.h"
#include "src/acc/acc_helper_functions.h"
#ifdef _CUDA_ENABLED
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#include "src/acc/cuda/cuda_kernels/wavg.cuh"
#include "src/acc/cuda/cuda_kernels/diff2.cuh"
#include "src/acc/cuda/cuda_fft.h"
#else
#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_kernels/wavg.h"
#include "src/acc/cpu/cpu_kernels/diff2.h"
#endif

void dump_array(char *name, bool *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%d, ", ptr[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_array(char *name, int *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%d, ", ptr[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_array(char *name, size_t *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%zu, ", ptr[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_array(char *name, float *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f, ", ptr[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_complex_array(char *name, ACCCOMPLEX *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f,%f, ", ptr[i].x, ptr[i].y);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_complex_array(char *name, Complex *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f,%f, ", ptr[i].real, ptr[i].imag);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_double_array(char *name, float *ptr, float *ptr2, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f,%f, ", ptr[i], ptr2[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_triple_array(char *name, float *ptr, float *ptr2, float *ptr3, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f,%f,%f, ", ptr[i], ptr2[i], ptr3[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_array(char *name, double *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f, ", ptr[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_double_array(char *name, double *ptr, double *ptr2, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f,%f, ", ptr[i], ptr2[i]);
			count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_triple_array(char *name, double *ptr, double *ptr2, double *ptr3, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f,%f,%f, ", ptr[i], ptr2[i], ptr3[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

namespace AccUtilities
{
template <class MlClass>
void makeNoiseImage(XFLOAT sigmaFudgeFactor,
		MultidimArray<RFLOAT > &sigmaNoiseSpectra,
		long int seed,
		MlClass *accMLO,
		AccPtr<XFLOAT> &RandomImage,
		bool is3D)
{
    // Different MPI-distributed subsets may otherwise have different instances of the random noise below,
    // because work is on an on-demand basis and therefore variable with the timing of distinct nodes...
    // Have the seed based on the part_id, so that each particle has a different instant of the noise
    init_random_generator(seed);

    // Make a holder for the spectral profile and copy to the GPU
    // AccDataTypes::Image<XFLOAT> NoiseSpectra(sigmaNoiseSpectra, ptrFactory);
    AccPtr<XFLOAT> NoiseSpectra = RandomImage.make<XFLOAT>(sigmaNoiseSpectra.nzyxdim);
    NoiseSpectra.allAlloc();

    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigmaNoiseSpectra)
            NoiseSpectra[n] = (XFLOAT)sqrt(sigmaFudgeFactor*sigmaNoiseSpectra.data[n]);

#ifdef _CUDA_ENABLED
    // Set up states to seeda and run randomization on the GPU
    // AccDataTypes::Image<curandState > RandomStates(RND_BLOCK_NUM*RND_BLOCK_SIZE,ptrFactory);
    AccPtr<curandState> RandomStates = RandomImage.make<curandState>(RND_BLOCK_NUM*RND_BLOCK_SIZE);
    RandomStates.deviceAlloc();

    NoiseSpectra.cpToDevice();
    NoiseSpectra.streamSync();
    LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);

    // Initialize randomization by particle ID, like on the CPU-side
    cuda_kernel_initRND<<<RND_BLOCK_NUM,RND_BLOCK_SIZE>>>(
                                     seed,
                                    ~RandomStates);
    LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);

    // Create noise image with the correct spectral profile
    if(is3D)
    {
    	cuda_kernel_RNDnormalDitributionComplexWithPowerModulation3D<<<RND_BLOCK_NUM,RND_BLOCK_SIZE>>>(
                                    ~accMLO->transformer1.fouriers,
                                    ~RandomStates,
									accMLO->transformer1.xFSize,
									accMLO->transformer1.yFSize,
                                    ~NoiseSpectra);
    }
    else
    {
    	cuda_kernel_RNDnormalDitributionComplexWithPowerModulation2D<<<RND_BLOCK_NUM,RND_BLOCK_SIZE>>>(
    	                                    ~accMLO->transformer1.fouriers,
    	                                    ~RandomStates,
    										accMLO->transformer1.xFSize,
    	                                    ~NoiseSpectra);
    }
    LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);

    // Transform to real-space, to get something which look like
    // the particle image without actual signal (a particle)
    accMLO->transformer1.backward();

    // Copy the randomized image to A separate device-array, so that the
    // transformer can be used to set up the actual particle image
    accMLO->transformer1.reals.cpOnDevice(~RandomImage);
    //cudaMLO->transformer1.reals.streamSync();
#else

    // Create noise image with the correct spectral profile
    if(is3D)
    	CpuKernels::RNDnormalDitributionComplexWithPowerModulation3D(accMLO->transformer1.fouriers(), accMLO->transformer1.xFSize, accMLO->transformer1.yFSize, ~NoiseSpectra);
    else
    	CpuKernels::RNDnormalDitributionComplexWithPowerModulation2D(accMLO->transformer1.fouriers(), accMLO->transformer1.xFSize, ~NoiseSpectra);

    // Transform to real-space, to get something which look like
	// the particle image without actual signal (a particle)
    accMLO->transformer1.backward();

	// Copy the randomized image to A separate device-array, so that the
	// transformer can be used to set up the actual particle image
	for(size_t i=0; i<RandomImage.getSize(); i++)
		RandomImage[i]=accMLO->transformer1.reals[i];

#endif
}

static void TranslateAndNormCorrect(MultidimArray<RFLOAT > &img_in,
		AccPtr<XFLOAT> &img_out,
		XFLOAT normcorr,
		RFLOAT xOff,
		RFLOAT yOff,
		RFLOAT zOff,
		bool DATA3D)
{
	//Temporary array because translate is out-of-place
	AccPtr<XFLOAT> temp = img_out.make<XFLOAT>(img_in.nzyxdim);
	temp.allAlloc();

	for (unsigned long i = 0; i < img_in.nzyxdim; i++)
		temp[i] = (XFLOAT) img_in.data[i];

	temp.cpToDevice();
	temp.streamSync();

	// Apply the norm_correction term
	if (normcorr!=1)
	{
#ifdef _CUDA_ENABLED
		int BSZ = ( (int) ceilf(( float)temp.getSize() /(float)BLOCK_SIZE));
		CudaKernels::cuda_kernel_multi<XFLOAT><<<BSZ,BLOCK_SIZE,0,temp.getStream()>>>(temp(),normcorr,temp.getSize());
#else
		CpuKernels::cpu_kernel_multi<XFLOAT>(temp(),normcorr, temp.getSize());
#endif
	}
	//LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);

	if(temp.getAccPtr()==img_out.getAccPtr())
		CRITICAL(ERRUNSAFEOBJECTREUSE);
#ifdef _CUDA_ENABLED
	int BSZ = ( (int) ceilf(( float)temp.getSize() /(float)BLOCK_SIZE));
	if (DATA3D)
		CudaKernels::cuda_kernel_translate3D<XFLOAT><<<BSZ,BLOCK_SIZE,0,temp.getStream()>>>(temp(),img_out(),img_in.zyxdim,img_in.xdim,img_in.ydim,img_in.zdim,xOff,yOff,zOff);
	else
		CudaKernels::cuda_kernel_translate2D<XFLOAT><<<BSZ,BLOCK_SIZE,0,temp.getStream()>>>(temp(),img_out(),img_in.zyxdim,img_in.xdim,img_in.ydim,xOff,yOff);
	//LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);
#else
	if (DATA3D)
		CpuKernels::cpu_translate3D<XFLOAT>(temp(),img_out(),img_in.zyxdim,img_in.xdim,img_in.ydim,img_in.zdim,xOff,yOff,zOff);
	else
		CpuKernels::cpu_translate2D<XFLOAT>(temp(),img_out(),img_in.zyxdim,img_in.xdim,img_in.ydim,xOff,yOff);
#endif
}
template <class MlClass>
void normalizeAndTransformImage(	AccPtr<XFLOAT> &img_in,
									MultidimArray<Complex > &img_out,
									MlClass *accMLO,
									size_t xSize,
									size_t ySize,
									size_t zSize)
{
			img_in.cpOnAcc(accMLO->transformer1.reals);
			runCenterFFT(
					accMLO->transformer1.reals,
					(int)accMLO->transformer1.xSize,
					(int)accMLO->transformer1.ySize,
					(int)accMLO->transformer1.zSize,
					false
					);
			accMLO->transformer1.reals.streamSync();
			accMLO->transformer1.forward();
			accMLO->transformer1.fouriers.streamSync();

			size_t FMultiBsize = ( (int) ceilf(( float)accMLO->transformer1.fouriers.getSize()*2/(float)BLOCK_SIZE));
			AccUtilities::multiply<XFLOAT>(FMultiBsize, BLOCK_SIZE, accMLO->transformer1.fouriers.getStream(),
							(XFLOAT*)~accMLO->transformer1.fouriers,
							(XFLOAT)1/((XFLOAT)(accMLO->transformer1.reals.getSize())),
							accMLO->transformer1.fouriers.getSize()*2);
			//LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);

			AccPtr<ACCCOMPLEX> d_Fimg = img_in.make<ACCCOMPLEX>(xSize * ySize * zSize);
			d_Fimg.allAlloc();
			accMLO->transformer1.fouriers.streamSync();
			windowFourierTransform2(
					accMLO->transformer1.fouriers,
					d_Fimg,
					accMLO->transformer1.xFSize,accMLO->transformer1.yFSize, accMLO->transformer1.zFSize, //Input dimensions
					xSize, ySize, zSize  //Output dimensions
					);
			accMLO->transformer1.fouriers.streamSync();

			d_Fimg.cpToHost();
			d_Fimg.streamSync();
			img_out.initZeros(zSize, ySize, xSize);
			for (unsigned long i = 0; i < img_out.nzyxdim; i ++)
			{
				img_out.data[i].real = (RFLOAT) d_Fimg[i].x;
				img_out.data[i].imag = (RFLOAT) d_Fimg[i].y;
			}

}

static void softMaskBackgroundValue(
	AccDataTypes::Image<XFLOAT> &vol,
	XFLOAT   radius,
	XFLOAT   radius_p,
	XFLOAT   cosine_width,
	AccPtr<XFLOAT> &g_sum,
	AccPtr<XFLOAT> &g_sum_bg)
{
	int block_dim = 128; //TODO: set balanced (hardware-dep?)
#ifdef _CUDA_ENABLED
		cuda_kernel_softMaskBackgroundValue<<<block_dim,SOFTMASK_BLOCK_SIZE,0, vol.getStream()>>>(
				~vol,
				vol.getxyz(),
				vol.getx(),
				vol.gety(),
				vol.getz(),
				vol.getx()/2,
				vol.gety()/2,
				vol.getz()/2,
				radius,
				radius_p,
				cosine_width,
				~g_sum,
				~g_sum_bg);
#else
	CpuKernels::softMaskBackgroundValue(
			block_dim,
			SOFTMASK_BLOCK_SIZE,
			~vol,
			vol.getxyz(),
			vol.getx(),
			vol.gety(),
			vol.getz(),
			vol.getx()/2,
			vol.gety()/2,
			vol.getz()/2,
			radius,
			radius_p,
			cosine_width,
			~g_sum,
			~g_sum_bg);
#endif
}

static void cosineFilter(
		AccDataTypes::Image<XFLOAT> &vol,
		bool do_Mnoise,
		AccDataTypes::Image<XFLOAT> Noise,
		XFLOAT radius,
		XFLOAT radius_p,
		XFLOAT cosine_width,
		XFLOAT sum_bg_total)
{
	int block_dim = 128; //TODO: set balanced (hardware-dep?)
#ifdef _CUDA_ENABLED
	cuda_kernel_cosineFilter<<<block_dim,SOFTMASK_BLOCK_SIZE,0,vol.getStream()>>>(
			~vol,
			vol.getxyz(),
			vol.getx(),
			vol.gety(),
			vol.getz(),
			vol.getx()/2,
			vol.gety()/2,
			vol.getz()/2,
			!do_Mnoise,
			~Noise,
			radius,
			radius_p,
			cosine_width,
			sum_bg_total);
#else
	CpuKernels::cosineFilter(
			block_dim,
			SOFTMASK_BLOCK_SIZE,
			~vol,
			vol.getxyz(),
			vol.getx(),
			vol.gety(),
			vol.getz(),
			vol.getx()/2,
			vol.gety()/2,
			vol.getz()/2,
			!do_Mnoise,
			~Noise,
			radius,
			radius_p,
			cosine_width,
			sum_bg_total);
#endif
}

void initOrientations(AccPtr<RFLOAT> &pdfs, AccPtr<XFLOAT> &pdf_orientation, AccPtr<bool> &pdf_orientation_zeros)
{
#ifdef _CUDA_ENABLED
	int bs = 512;
	int gs = ceil(pdfs.getSize()/(float)(bs));
	cuda_kernel_initOrientations<<<gs, bs, 0, pdfs.getStream()>>>(~pdfs, ~pdf_orientation, ~pdf_orientation_zeros, pdfs.getSize());
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
#else
	for(int iorientclass=0; iorientclass< pdfs.getSize(); iorientclass++)
	{
		if (pdfs[iorientclass] == 0)
		{
			pdf_orientation[iorientclass] = 0.f;
			pdf_orientation_zeros[iorientclass] = true;
		}
		else
		{
			pdf_orientation[iorientclass] = log(pdfs[iorientclass]);
			pdf_orientation_zeros[iorientclass] = false;
		}
	}

#endif
}

void centerFFT_2D(int grid_size, int batch_size, int block_size,
				cudaStream_t stream,
				XFLOAT *img_in,
				size_t image_size,
				int xdim,
				int ydim,
				int xshift,
				int yshift)
{
#ifdef _CUDA_ENABLED
	dim3 blocks(grid_size, batch_size);
	cuda_kernel_centerFFT_2D<<<blocks,block_size,0,stream>>>(
				img_in,
				image_size,
				xdim,
				ydim,
				xshift,
				yshift);
#else
	CpuKernels::centerFFT_2D<XFLOAT>(batch_size, 0, image_size/2,
				img_in,
				image_size,
				xdim,
				ydim,
				xshift,
				yshift);
#endif	
}

void centerFFT_2D(int grid_size, int batch_size, int block_size,
				XFLOAT *img_in,
				size_t image_size,
				int xdim,
				int ydim,
				int xshift,
				int yshift)
{
#ifdef _CUDA_ENABLED
	dim3 blocks(grid_size, batch_size);
	cuda_kernel_centerFFT_2D<<<blocks,block_size>>>(
				img_in,
				image_size,
				xdim,
				ydim,
				xshift,
				yshift);
#else
	CpuKernels::centerFFT_2D<XFLOAT>(batch_size, 0, image_size/2,
			img_in,
			image_size,
			xdim,
			ydim,
			xshift,
			yshift);
#endif
}

void centerFFT_3D(int grid_size, int batch_size, int block_size,
				cudaStream_t stream,
				XFLOAT *img_in,
				size_t image_size,
				int xdim,
				int ydim,
				int zdim,
				int xshift,
				int yshift,
				int zshift)
{
#ifdef _CUDA_ENABLED
	dim3 blocks(grid_size, batch_size);
	cuda_kernel_centerFFT_3D<<<blocks,block_size, 0, stream>>>(
				img_in,
				image_size,
				xdim,
				ydim,
				zdim,
				xshift,
				yshift,
				zshift);
#else
	CpuKernels::centerFFT_3D<XFLOAT>(batch_size, (size_t)0, (size_t)image_size/2,
			img_in,
			image_size,
			xdim,
			ydim,
			zdim,
			xshift,
			yshift,
			zshift);
#endif
}

void kernel_exponentiate_weights_fine(	XFLOAT *g_pdf_orientation,
										bool *g_pdf_orientation_zeros,
										XFLOAT *g_pdf_offset,
										bool *g_pdf_offset_zeros,
										XFLOAT *g_weights,
										XFLOAT min_diff2,
										unsigned long  oversamples_orient,
										unsigned long  oversamples_trans,
										unsigned long *d_rot_id,
										unsigned long *d_trans_idx,
										unsigned long *d_job_idx,
										unsigned long *d_job_num,
										long int job_num,
										cudaStream_t stream)
{
	long block_num = ceil((double)job_num / (double)SUMW_BLOCK_SIZE);

#ifdef _CUDA_ENABLED
cuda_kernel_exponentiate_weights_fine<<<block_num,SUMW_BLOCK_SIZE,0,stream>>>(
		g_pdf_orientation,
		g_pdf_orientation_zeros,
		g_pdf_offset,
		g_pdf_offset_zeros,
		g_weights,
		min_diff2,
		oversamples_orient,
		oversamples_trans,
		d_rot_id,
		d_trans_idx,
		d_job_idx,
		d_job_num,
		job_num);
#else
	CpuKernels::exponentiate_weights_fine(
		g_pdf_orientation,
		g_pdf_orientation_zeros,
		g_pdf_offset,
		g_pdf_offset_zeros,
		g_weights,
		min_diff2,
		oversamples_orient,
		oversamples_trans,
		d_rot_id,
		d_trans_idx,
		d_job_idx,
		d_job_num,
		job_num);
#endif
}

};  // namespace AccUtilities

void run_griddingCorrect(RFLOAT *vol, int interpolator, RFLOAT rrval, RFLOAT r_min_nn,
								size_t iX, size_t iY, size_t iZ)
{
#ifdef CUDA
	dim3 bs(32,4,2);
	dim3 gs(ceil(iX/(float)bs.x), ceil(iY/(float)bs.y), ceil(iZ/(float)bs.z));
	cuda_kernel_griddingCorrect<<<gs,bs>>>(vol, interpolator, rrval, r_min_nn, iX, iY, iZ);
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
#endif
}

void run_padTranslatedMap(
		RFLOAT *d_in, RFLOAT *d_out,
		size_t isX, size_t ieX, size_t isY, size_t ieY, size_t isZ, size_t ieZ, //Input dimensions
		size_t osX, size_t oeX, size_t osY, size_t oeY, size_t osZ, size_t oeZ,  //Output dimensions
		cudaStream_t stream)
{
#ifdef CUDA
	size_t iszX = ieX - isX + 1; 
	size_t iszY = ieY - isY + 1; 
	size_t iszZ = ieZ - isZ + 1; 
	size_t oszX = oeX - osX + 1; 
	size_t oszY = oeY - osY + 1; 
	size_t oszZ = oeZ - osZ + 1; 
   
	if(iszX == oszX && iszY == oszY && iszZ == oszZ)
	{
		cudaCpyDeviceToDevice(d_in, d_out, iszX*iszY*iszZ, stream);
	}
	else
	{
		dim3 block_dim(16,4,2);
		dim3 grid_dim(ceil(oszX / (float) block_dim.x), ceil(oszY / (float) block_dim.y), ceil(oszZ / (float) block_dim.z));
		cuda_kernel_window_transform<RFLOAT><<< grid_dim, block_dim, 0, stream >>>(
				d_in, d_out,
				iszX, iszY, iszZ, //Input dimensions
				isX-osX, isY-osY, isZ-osZ, oszX, oszY, oszZ  //Output dimensions
				);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
	}
#endif
}

void run_CenterFFTbySign(Complex *img_in, int xSize, int ySize, int zSize, cudaStream_t stream)
{
#ifdef CUDA
	dim3 bs(32,4,2);
	dim3 gs(ceil(xSize/(float)bs.x), ceil(ySize/(float)bs.y), ceil(zSize/(float)bs.z));
	if(sizeof(RFLOAT) == sizeof(double))
		cuda_kernel_centerFFTbySign<<<gs,bs, 0, stream>>>(
				(double2*)img_in,
				xSize,
				ySize,
				zSize);
	else
		cuda_kernel_centerFFTbySign<<<gs,bs, 0, stream>>>(
				(float2*)img_in,
				xSize,
				ySize,
				zSize);
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
#endif
}

#endif //ACC_UTILITIES_H_

