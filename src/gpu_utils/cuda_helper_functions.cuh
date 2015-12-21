#ifndef CUDA_HELPER_FUNCTIONS_CUH_
#define CUDA_HELPER_FUNCTIONS_CUH_

#include "src/gpu_utils/cuda_ml_optimiser.h"
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_projector.cuh"
#include "src/gpu_utils/cuda_benchmark_utils.cuh"
#include "src/gpu_utils/cuda_mem_utils.h"
#include "src/gpu_utils/cuda_kernels/helper.cuh"
#include "src/gpu_utils/cuda_kernels/diff2.cuh"
#include "src/gpu_utils/cuda_kernels/wavg.cuh"
#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include "src/complex.h"
#include <fstream>
#include <cuda_runtime.h>
#include "src/parallel.h"
#include <signal.h>

/*
 * This assisting function goes over the orientations determined as significant for this image, and checks
 * which translations should be included in the list of those which differences will be calculated for.
 *
 * Any contiguous translations with a shared orientation are grouped together into a "job" which is supplied
 * to the difference kernel. If there are more contiguous translations than the specified PROJDIFF_CHUNK_SIZE,
 * these are split into separate jobs, to increase paralllelism at the cost of redundant memory reads.
 */
long int makeJobsForDiff2Fine(
		OptimisationParamters &op,  SamplingParameters &sp,
		long int orientation_num, long int translation_num,
		ProjectionParams &FineProjectionData,
		std::vector< long unsigned > &iover_transes,
		std::vector< long unsigned > &ihiddens,
		long int nr_over_orient, long int nr_over_trans, int ipart,
		IndexedDataArray &FPW, // FPW=FinePassWeights
		IndexedDataArrayMask &dataMask);

/*
 * This assisting function goes over the weight-array and groups all weights with shared
 * orientations into 'jobs' which are fed into the collect-kenrel, which reduces all translations
 * with computed differences into a reduced object to be back-projected.
 */
int  makeJobsForCollect(IndexedDataArray &FPW, IndexedDataArrayMask &dataMask, unsigned long NewJobNum); // FPW=FinePassWeights

/*
 * Maps weights to a decoupled indexing of translations and orientations
 */
void mapWeights(
		unsigned long orientation_start,
		XFLOAT *mapped_weights,
		unsigned orientation_num,
		unsigned long idxArr_start,
		unsigned long idxArr_end,
		unsigned translation_num,
		XFLOAT *weights,
		long unsigned *rot_idx,
		long unsigned *trans_idx,
		unsigned long current_oversampling);

void buildCorrImage(MlOptimiser *baseMLO, OptimisationParamters &op, CudaGlobalPtr<XFLOAT> &corr_img, long int ipart, long int group_id);

void generateEulerMatrices(
		XFLOAT padding_factor,
		ProjectionParams &ProjectionData,
		XFLOAT *eulers,
		bool inverse);

long unsigned generateProjectionSetupFine(
		OptimisationParamters &op,
		SamplingParameters &sp,
		MlOptimiser *baseMLO,
		unsigned iclass,
		ProjectionParams &ProjectionData);

void runWavgKernel(
		CudaProjectorKernel &projector,
		XFLOAT *eulers,
		XFLOAT *Fimgs_real,
		XFLOAT *Fimgs_imag,
		XFLOAT *Fimgs_nomask_real,
		XFLOAT *Fimgs_nomask_imag,
		XFLOAT *sorted_weights,
		XFLOAT *ctfs,
		XFLOAT *Minvsigma2s,
		XFLOAT *wdiff2s_parts,
		XFLOAT *wdiff2s_AA,
		XFLOAT *wdiff2s_XA,
		XFLOAT *wavgs_real,
		XFLOAT *wavgs_imag,
		XFLOAT *Fweights,
		OptimisationParamters &op,
		MlOptimiser *baseMLO,
		long unsigned orientation_num,
		long unsigned translation_num,
		unsigned image_size,
		long int ipart,
		int group_id,
		int exp_iclass,
		XFLOAT part_scale,
		cudaStream_t stream);

#define INIT_VALUE_BLOCK_SIZE 512
template< typename T>
__global__ void cuda_kernel_init_value(
		T *data,
		T value,
		size_t size)
{
	size_t idx = blockIdx.x * INIT_VALUE_BLOCK_SIZE + threadIdx.x;
	if (idx < size)
		data[idx] = value;
}

template< typename T>
void deviceInitValue(CudaGlobalPtr<T> &data, T value)
{
	int grid_size = ceil((float)data.getSize()/(float)INIT_VALUE_BLOCK_SIZE);
	cuda_kernel_init_value<T><<< grid_size, INIT_VALUE_BLOCK_SIZE, 0, data.getStream() >>>(
			~data,
			value,
			data.getSize());
}

#define WEIGHT_MAP_BLOCK_SIZE 512
__global__ void cuda_kernel_allweights_to_mweights(
		unsigned long * d_iorient,
		XFLOAT * d_allweights,
		XFLOAT * d_mweights,
		unsigned long orientation_num,
		unsigned long translation_num
		);

void mapAllWeightsToMweights(
		unsigned long * d_iorient, //projectorPlan.iorientclasses
		XFLOAT * d_allweights, //allWeights
		XFLOAT * d_mweights, //Mweight
		unsigned long orientation_num, //projectorPlan.orientation_num
		unsigned long translation_num, //translation_num
		cudaStream_t stream
		);

#define OVER_THRESHOLD_BLOCK_SIZE 512
template< typename T>
__global__ void cuda_kernel_array_over_threshold(
		T *data,
		bool *passed,
		T threshold,
		size_t size)
{
	size_t idx = blockIdx.x * OVER_THRESHOLD_BLOCK_SIZE + threadIdx.x;
	if (idx < size)
	{
		if (data[idx] >= threshold)
			passed[idx] = true;
		else
			passed[idx] = false;
	}
}

template< typename T>
void arrayOverThreshold(CudaGlobalPtr<T> &data, CudaGlobalPtr<bool> &passed, T threshold)
{
	int grid_size = ceil((float)data.getSize()/(float)OVER_THRESHOLD_BLOCK_SIZE);
	cuda_kernel_array_over_threshold<T><<< grid_size, OVER_THRESHOLD_BLOCK_SIZE, 0, data.getStream() >>>(
			~data,
			~passed,
			threshold,
			data.getSize());
}

#define FIND_IN_CUMULATIVE_BLOCK_SIZE 512
template< typename T>
__global__ void cuda_kernel_find_threshold_idx_in_cumulative(
		T *data,
		T threshold,
		size_t size_m1, //data size minus 1
		size_t *idx)
{
	size_t i = blockIdx.x * FIND_IN_CUMULATIVE_BLOCK_SIZE + threadIdx.x;
	if (i < size_m1 && data[i] <= threshold && threshold < data[i+1])
		idx[0] = i+1;
}

size_t findThresholdIdxInCumulativeSum(CudaGlobalPtr<XFLOAT> &data, XFLOAT threshold);

void runDiff2KernelCoarse(
		CudaProjectorKernel &projector,
		XFLOAT *trans_x,
		XFLOAT *trans_y,
		XFLOAT *corr_img,
		XFLOAT *Fimgs_real,
		XFLOAT *Fimgs_imag,
		XFLOAT *d_eulers,
		XFLOAT *diff2s,
		OptimisationParamters &op,
		MlOptimiser *baseMLO,
		long unsigned orientation_num,
		int translation_num,
		int image_size,
		int ipart,
		int group_id,
		int exp_iclass,
		cudaStream_t stream,
		bool do_CC);

void runDiff2KernelFine(
		CudaProjectorKernel &projector,
		XFLOAT *corr_img,
		XFLOAT *Fimgs_real,
		XFLOAT *Fimgs_imag,
		XFLOAT *eulers,
		long unsigned *rot_id,
		long unsigned *rot_idx,
		long unsigned *trans_idx,
		long unsigned *job_idx,
		long unsigned *job_num,
		XFLOAT *diff2s,
		OptimisationParamters &op,
		MlOptimiser *baseMLO,
		long unsigned orientation_num,
		long unsigned translation_num,
		long unsigned significant_num,
		unsigned image_size,
		int ipart,
		int exp_iclass,
		cudaStream_t stream,
		long unsigned job_num_count,
		bool do_CC);

#define WINDOW_FT_BLOCK_SIZE 128
template<bool check_max_r2>
__global__ void cuda_kernel_window_fourier_transform(
		XFLOAT *g_in_real,
		XFLOAT *g_in_imag,
		XFLOAT *g_out_real,
		XFLOAT *g_out_imag,
		unsigned iX, unsigned iY, unsigned iZ, unsigned iYX, //Input dimensions
		unsigned oX, unsigned oY, unsigned oZ, unsigned oYX, //Output dimensions
		unsigned max_idx,
		unsigned max_r2 = 0
		)
{
	unsigned n = threadIdx.x + WINDOW_FT_BLOCK_SIZE * blockIdx.x;
	if (n >= max_idx) return;

	int k, i, kp, ip, jp;

	if (check_max_r2)
	{
		k = n / (iX * iY);
		i = (n % (iX * iY)) / iX;

		kp = k < iX ? k : k - iZ;
		ip = i < iX ? i : i - iY;
		jp = n % iX;

		if (kp*kp + ip*ip + jp*jp > max_r2)
			return;
	}
	else
	{
		k = n / (oX * oY);
		i = (n % (oX * oY)) / oX;

		kp = k < oX ? k : k - oZ;
		ip = i < oX ? i : i - oY;
		jp = n % oX;
	}

	g_out_real[(kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp] = g_in_real[(kp < 0 ? kp + iZ : kp)*iYX + (ip < 0 ? ip + iY : ip)*iX + jp];
	g_out_imag[(kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp] = g_in_imag[(kp < 0 ? kp + iZ : kp)*iYX + (ip < 0 ? ip + iY : ip)*iX + jp];
}

void windowFourierTransform2(
		XFLOAT *d_in_real,
		XFLOAT *d_in_imag,
		XFLOAT *d_out_real,
		XFLOAT *d_out_imag,
		unsigned iX, unsigned iY, unsigned iZ, //Input dimensions
		unsigned oX, unsigned oY, unsigned oZ,  //Output dimensions
		cudaStream_t stream = 0);

#define WINDOW_FT_BLOCK_SIZE 128
template<bool check_max_r2>
__global__ void cuda_kernel_window_fourier_transform(
		CUDACOMPLEX *g_in,
		CUDACOMPLEX *g_out,
		unsigned iX, unsigned iY, unsigned iZ, unsigned iYX, //Input dimensions
		unsigned oX, unsigned oY, unsigned oZ, unsigned oYX, //Output dimensions
		unsigned max_idx,
		unsigned max_r2 = 0
		)
{
	unsigned n = threadIdx.x + WINDOW_FT_BLOCK_SIZE * blockIdx.x;
	if (n >= max_idx) return;

	int k, i, kp, ip, jp;

	if (check_max_r2)
	{
		k = n / (iX * iY);
		i = (n % (iX * iY)) / iX;

		kp = k < iX ? k : k - iZ;
		ip = i < iX ? i : i - iY;
		jp = n % iX;

		if (kp*kp + ip*ip + jp*jp > max_r2)
			return;
	}
	else
	{
		k = n / (oX * oY);
		i = (n % (oX * oY)) / oX;

		kp = k < oX ? k : k - oZ;
		ip = i < oX ? i : i - oY;
		jp = n % oX;
	}

	g_out[(kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp].x = g_in[(kp < 0 ? kp + iZ : kp)*iYX + (ip < 0 ? ip + iY : ip)*iX + jp].x;
	g_out[(kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp].y = g_in[(kp < 0 ? kp + iZ : kp)*iYX + (ip < 0 ? ip + iY : ip)*iX + jp].y;
}

void windowFourierTransform2(
		CudaGlobalPtr<CUDACOMPLEX > &d_in,
		CudaGlobalPtr<CUDACOMPLEX > &d_out,
		unsigned iX, unsigned iY, unsigned iZ, //Input dimensions
		unsigned oX, unsigned oY, unsigned oZ,  //Output dimensions
		cudaStream_t stream = 0);


void selfApplyBeamTilt2(MultidimArray<Complex > &Fimg, RFLOAT beamtilt_x, RFLOAT beamtilt_y,
		RFLOAT wavelength, RFLOAT Cs, RFLOAT angpix, int ori_size);

template <typename T>
void runCenterFFT(MultidimArray< T >& v, bool forward, CudaCustomAllocator *allocator)
{
	CudaGlobalPtr<XFLOAT >  img_in (v.nzyxdim, allocator);   // with original data pointer
//	CudaGlobalPtr<XFLOAT >  img_aux(v.nzyxdim, allocator);   // temporary holder

	for (unsigned i = 0; i < v.nzyxdim; i ++)
		img_in[i] = (XFLOAT) v.data[i];

	img_in.put_on_device();
//	img_aux.device_alloc();

	if ( v.getDim() == 1 )
	{
		std::cerr << "CenterFFT on gpu reverts to cpu for dim!=2 (now dim=1)" <<std::endl;
		// 1D
		MultidimArray< T > aux;
		int l, shift;

		l = XSIZE(v);
		aux.resize(l);
		shift = (int)(l / 2);

		if (!forward)
			shift = -shift;

		// Shift the input in an auxiliar vector
		for (int i = 0; i < l; i++)
		{
			int ip = i + shift;

			if (ip < 0)
				ip += l;
			else if (ip >= l)
				ip -= l;

			aux(ip) = DIRECT_A1D_ELEM(v, i);
		}

		// Copy the vector
		for (int i = 0; i < l; i++)
			DIRECT_A1D_ELEM(v, i) = DIRECT_A1D_ELEM(aux, i);
	}
	else if ( v.getDim() == 2 )
	{
		// 2D
		//std::cerr << "CenterFFT on gpu with dim=2!" <<std::endl;

		long int xshift = (int)(XSIZE(v) / 2);
		long int yshift = (int)(YSIZE(v) / 2);

		if (!forward)
		{
			xshift = -xshift;
			yshift = -yshift;
		}


		dim3 dim(ceilf((float)(v.nzyxdim/(float)(2*CFTT_BLOCK_SIZE))));
		cuda_kernel_centerFFT_2D<<<dim,CFTT_BLOCK_SIZE>>>(img_in.d_ptr,
										  v.nzyxdim,
										  XSIZE(v),
										  YSIZE(v),
										  xshift,
										  yshift);

		img_in.cp_to_host();

//		HANDLE_ERROR(cudaStreamSynchronize(0));

		for (unsigned i = 0; i < v.nzyxdim; i ++)
			v.data[i] = (T) img_in[i];

	}
	else if ( v.getDim() == 3 )
	{
		std::cerr << "CenterFFT on gpu reverts to cpu for dim!=2 (now dim=3)" <<std::endl;
		// 3D
		MultidimArray< T > aux;
		int l, shift;

		// Shift in the X direction
		l = XSIZE(v);
		aux.resize(l);
		shift = (int)(l / 2);

		if (!forward)
			shift = -shift;

		for (int k = 0; k < ZSIZE(v); k++)
			for (int i = 0; i < YSIZE(v); i++)
			{
				// Shift the input in an auxiliar vector
				for (int j = 0; j < l; j++)
				{
					int jp = j + shift;

					if (jp < 0)
						jp += l;
					else if (jp >= l)
						jp -= l;

					aux(jp) = DIRECT_A3D_ELEM(v, k, i, j);
				}

				// Copy the vector
				for (int j = 0; j < l; j++)
					DIRECT_A3D_ELEM(v, k, i, j) = DIRECT_A1D_ELEM(aux, j);
			}

		// Shift in the Y direction
		l = YSIZE(v);
		aux.resize(l);
		shift = (int)(l / 2);

		if (!forward)
			shift = -shift;

		for (int k = 0; k < ZSIZE(v); k++)
			for (int j = 0; j < XSIZE(v); j++)
			{
				// Shift the input in an auxiliar vector
				for (int i = 0; i < l; i++)
				{
					int ip = i + shift;

					if (ip < 0)
						ip += l;
					else if (ip >= l)
						ip -= l;

					aux(ip) = DIRECT_A3D_ELEM(v, k, i, j);
				}

				// Copy the vector
				for (int i = 0; i < l; i++)
					DIRECT_A3D_ELEM(v, k, i, j) = DIRECT_A1D_ELEM(aux, i);
			}

		// Shift in the Z direction
		l = ZSIZE(v);
		aux.resize(l);
		shift = (int)(l / 2);

		if (!forward)
			shift = -shift;

		for (int i = 0; i < YSIZE(v); i++)
			for (int j = 0; j < XSIZE(v); j++)
			{
				// Shift the input in an auxiliar vector
				for (int k = 0; k < l; k++)
				{
					int kp = k + shift;
					if (kp < 0)
						kp += l;
					else if (kp >= l)
						kp -= l;

					aux(kp) = DIRECT_A3D_ELEM(v, k, i, j);
				}

				// Copy the vector
				for (int k = 0; k < l; k++)
					DIRECT_A3D_ELEM(v, k, i, j) = DIRECT_A1D_ELEM(aux, k);
			}
	}
	else
	{
		v.printShape();
		REPORT_ERROR("CenterFFT ERROR: Dimension should be 1, 2 or 3");
	}
}


template <typename T>
void runCenterFFT( CudaGlobalPtr< T > &img_in,
				  int xSize,
				  int ySize,
				  bool forward,
				  CudaCustomAllocator *allocator)
{
//	CudaGlobalPtr<XFLOAT >  img_aux(img_in.h_ptr, img_in.size, allocator);   // temporary holder
//	img_aux.device_alloc();

	int xshift = (xSize / 2);
	int yshift = (ySize / 2);

	if (!forward)
	{
		xshift = -xshift;
		yshift = -yshift;
	}

	dim3 dim(ceilf((float)(img_in.size/(float)(2*CFTT_BLOCK_SIZE))));
	cuda_kernel_centerFFT_2D<<<dim,CFTT_BLOCK_SIZE>>>(img_in.d_ptr,
													  img_in.size,
													  xSize,
													  ySize,
													  xshift,
													  yshift);

//	HANDLE_ERROR(cudaStreamSynchronize(0));
//	img_aux.cp_on_device(img_in.d_ptr); //update input image with centered kernel-output.


}

#endif //CUDA_HELPER_FUNCTIONS_CUH_

