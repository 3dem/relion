#ifndef ACC_UTILITIES_H_
#define ACC_UTILITIES_H_

#include "src/acc/acc_ptr.h"
#include "src/acc/data_types.h"
#ifdef CUDA
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#include "src/acc/cuda/cuda_kernels/wavg.cuh"
#include "src/acc/cuda/cuda_kernels/diff2.cuh"
#else
#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_kernels/wavg.h"
#include "src/acc/cpu/cpu_kernels/diff2.h"
#endif

namespace AccUtilities
{
	
template <typename T>
static void multiply(int block_size, AccDataTypes::Image<T> &ptr, T value)
{
#ifdef CUDA
	int BSZ = ( (int) ceilf(( float)ptr.getSize() /(float)block_size));
	CudaKernels::cuda_kernel_multi<T><<<BSZ,block_size,0,ptr.getStream()>>>(
		ptr(),
		value,
		ptr.getSize());
#else
	CpuKernels::cpu_kernel_multi<T>(
	ptr(),
	value,
	ptr.getSize());
#endif
}
	
template <typename T>
static void multiply(int MultiBsize, int block_size, cudaStream_t stream, T *array, T value, size_t size)
{
#ifdef CUDA
	CudaKernels::cuda_kernel_multi<T><<<MultiBsize,block_size,0,stream>>>(
		array,
		value,
		size);
#else
	CpuKernels::cpu_kernel_multi<T>(
	array,
	value,
	size);
#endif
}

template <typename T>
static void translate(int block_size,
	AccDataTypes::Image<T> &in,
	AccDataTypes::Image<T> &out,
	int dx,
	int dy,
	int dz=0)
{
#ifdef CUDA
int BSZ = ( (int) ceilf(( float)in.getxyz() /(float)block_size));

if (in.is3D())
{
	CudaKernels::cuda_kernel_translate3D<T><<<BSZ,block_size,0,in.getStream()>>>(
		in(),
		out(),
		in.getxyz(),
		in.getx(),
		in.gety(),
		in.getz(),
		dx,
		dy,
		dz);
}
else
{
	CudaKernels::cuda_kernel_translate2D<T><<<BSZ,block_size,0,in.getStream()>>>(
		in(),
		out(),
		in.getxyz(),
		in.getx(),
		in.gety(),
		dx,
		dy);
}
#else
if (in.is3D())
{
	CpuKernels::cpu_translate3D<T>(
		in(),
		out(),
		in.getxyz(),
		in.getx(),
		in.gety(),
		in.getz(),
		dx,
		dy,
		dz);
}
else
{
	CpuKernels::cpu_translate2D<T>(
		in(),
		out(),
		in.getxyz(),
		in.getx(),
		in.gety(),
		dx,
		dy);
}
#endif
}
	
	
template <typename T>
static T getSumOnDevice(AccPtr<T> &ptr)
{
#ifdef CUDA
	return CudaKernels::getSumOnDevice<T>(ptr);
#else
#ifdef DEBUG_CUDA
	if (ptr.size == 0)
		printf("DEBUG_ERROR: getSumOnDevice called with pointer of zero size.\n");
	if (ptr.hPtr == NULL)
		printf("DEBUG_ERROR: getSumOnDevice called with null device pointer.\n");
#endif
	size_t size = ptr.getSize();
	T sum = 0;
	for (int i=0; i<size; i++)
		sum += ptr[i];
	return sum;
#endif
}

template <typename T>
static T getMinOnDevice(AccPtr<T> &ptr)
{
#ifdef CUDA
	return CudaKernels::getMinOnDevice<T>(ptr);
#else
#ifdef DEBUG_CUDA
	if (ptr.size == 0)
		printf("DEBUG_ERROR: getMinOnDevice called with pointer of zero size.\n");
	if (ptr.hPtr == NULL)
		printf("DEBUG_ERROR: getMinOnDevice called with null device pointer.\n");
#endif
	return CpuKernels::getMin<T>(ptr(), ptr.getSize());
#endif
}

template <typename T>
static std::pair<int, T> getArgMinOnDevice(AccPtr<T> &ptr)
{
#ifdef CUDA
	return CudaKernels::getArgMinOnDevice<T>(ptr);
#else
#ifdef DEBUG_CUDA
	if (ptr.size == 0)
		printf("DEBUG_ERROR: getArgMinOnDevice called with pointer of zero size.\n");
	if (ptr.hPtr == NULL)
		printf("DEBUG_ERROR: getArgMinOnDevice called with null device pointer.\n");
#endif
	return CpuKernels::getArgMin<T>(ptr(), ptr.getSize());
#endif
}

template <typename T>
static std::pair<int, T> getArgMaxOnDevice(AccPtr<T> &ptr)
{
#ifdef CUDA
	return CudaKernels::getArgMaxOnDevice<T>(ptr);
#else
#ifdef DEBUG_CUDA
	if (ptr.size == 0)
		printf("DEBUG_ERROR: getArgMaxOnDevice called with pointer of zero size.\n");
	if (ptr.hPtr == NULL)
		printf("DEBUG_ERROR: getArgMaxOnDevice called with null device pointer.\n");
#endif
	return CpuKernels::getArgMax<T>(ptr(), ptr.getSize());
#endif
}

template <typename T>
static int filterGreaterZeroOnDevice(AccPtr<T> &in, AccPtr<T> &out)
{
#ifdef CUDA
	CudaKernels::MoreThanCubOpt<T> moreThanOpt(0.);
	return CudaKernels::filterOnDevice(in, out, moreThanOpt);
#else
	size_t arr_size = in.getSize();
	size_t filt_size = 0;
	size_t outindex = 0;
	// Find how many entries the output array will have
	for(int i=0; i<arr_size; i++)
	{
		if(in[i] > (T)0.0)
			filt_size++;
	}
#ifdef DEBUG_CUDA
	if (filt_size==0)
		ACC_PTR_DEBUG_FATAL("filterGreaterZeroOnDevice - No filtered values greater than 0.\n");
#endif
	out.resizeHost(filt_size);
	// Now populate output array
	for(int i=0; i<arr_size; i++)
		if(in[i] > 0.0) {
			out[outindex] = in[i];
			outindex++;
		}
	return filt_size;
#endif
}

template <typename T>
static void sortOnDevice(AccPtr<T> &in, AccPtr<T> &out)
{
#ifdef CUDA
	CudaKernels::sortOnDevice(in, out);
#else
	//TODO - convert ACCPTR to store data as vector so we don't need to make
	//an extra copies here.  For now, nasty hack
	size_t arr_size = in.getSize();
	std::vector<T> sortVector(in(), in() + in.getSize());
	sort(sortVector.begin(), sortVector.end());
	for (int i=0; i < arr_size; i++)
		out[i] = sortVector[i];
#endif
}

template <typename T>
static void scanOnDevice(AccPtr<T> &in, AccPtr<T> &out)
{
#ifdef CUDA
	CudaKernels::scanOnDevice(in, out);
#else
	T sum = 0.0;
	size_t arr_size = in.getSize();
	for(int i=0; i<arr_size; i++) 
	{
		sum += in[i];
		out[i] = sum;
	}
#endif
}

static void softMaskBackgroundValue(
		int inblock_dim, 
		int inblock_size,
		XFLOAT *vol,
		Image<RFLOAT> &img,
		bool     do_Mnoise,
		XFLOAT   radius,
		XFLOAT   radius_p,
		XFLOAT   cosine_width,
		XFLOAT  *g_sum,
		XFLOAT  *g_sum_bg);


static void cosineFilter(
		int inblock_dim, 
		int inblock_size,
		XFLOAT *vol,
		long int vol_size,
		long int xdim,
		long int ydim,
		long int zdim,
		long int xinit,
		long int yinit,
		long int zinit,
		bool do_Mnoise,
		XFLOAT radius,
		XFLOAT radius_p,
		XFLOAT cosine_width,
		XFLOAT sum_bg_total);

template<bool DATA3D>
void powerClass(int		in_gridSize,
				int			in_blocksize,
				ACCCOMPLEX  *g_image,
				XFLOAT      *g_spectrum,
				int          image_size,
				int          spectrum_size,
				int          xdim,
				int          ydim,
				int          zdim,
				int          res_limit,
				XFLOAT      *g_highres_Xi2)
{
#ifdef CUDA
	dim3 grid_size(in_gridSize);
	cuda_kernel_powerClass<DATA3D><<<grid_size,in_blocksize,0,0>>>(g_image,
		g_spectrum,
		image_size,
		spectrum_size,
		xdim,
		ydim,
		zdim,
		res_limit,
		g_highres_Xi2);
#else
	CpuKernels::powerClass<DATA3D>(in_gridSize,
		g_image,
		g_spectrum,
		image_size,
		spectrum_size,
		xdim,
		ydim,
		zdim,
		res_limit,
		g_highres_Xi2);
#endif
}


template<bool invert>
void acc_make_eulers_2D(int grid_size, int block_size,
		cudaStream_t stream,
		XFLOAT *alphas,
		XFLOAT *eulers,
		unsigned orientation_num)
{
#ifdef CUDA
	cuda_kernel_make_eulers_2D<invert><<<grid_size,block_size,0,stream>>>(
		alphas,
		eulers,
		orientation_num);
#else
	CpuKernels::cpu_kernel_make_eulers_2D<invert>(grid_size, block_size,
		alphas, eulers, orientation_num);
#endif
}

template<bool invert,bool perturb>
void acc_make_eulers_3D(int grid_size, int block_size,
		cudaStream_t stream,
		XFLOAT *alphas,
		XFLOAT *betas,
		XFLOAT *gammas,
		XFLOAT *eulers,
		unsigned orientation_num,
		XFLOAT *R)
{
#ifdef CUDA
	cuda_kernel_make_eulers_3D<invert,perturb><<<grid_size,block_size,0,stream>>>(
		alphas,
		betas,
		gammas,
		eulers,
		orientation_num,
		R);
#else
	CpuKernels::cpu_kernel_make_eulers_3D<invert,perturb>(grid_size, block_size,
		alphas,
		betas,
		gammas,
		eulers,
		orientation_num,
		R);
#endif
}

#ifdef CUDA
	#define INIT_VALUE_BLOCK_SIZE 512
#endif 

template< typename T>
void InitComplexValue(AccPtr<T> &data, XFLOAT value)
{
#ifdef CUDA
	int grid_size = ceil((float)(data.getSize())/(float)INIT_VALUE_BLOCK_SIZE);
	cuda_kernel_init_complex_value<T><<< grid_size, INIT_VALUE_BLOCK_SIZE, 0, data.getStream() >>>(
			~data,
			value,
			data.getSize(), INIT_VALUE_BLOCK_SIZE);
#else
	int Size = data.getSize();
	for(size_t i=0; i<Size; i++)
	{
		data[i].x = value;
		data[i].y = value;
	}
#endif
}

template< typename T>
void InitValue(AccPtr<T> &data, T value)
{
#ifdef CUDA
	int grid_size = ceil((float)data.getSize()/(float)INIT_VALUE_BLOCK_SIZE);
	cuda_kernel_init_value<T><<< grid_size, INIT_VALUE_BLOCK_SIZE, 0, data.getStream() >>>(
			~data,
			value,
			data.getSize(),
			INIT_VALUE_BLOCK_SIZE);
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
#else
	int Size = data.getSize();
	for (size_t i=0; i < Size; i++)
		data.hPtr[i] = value;
#endif
}

template< typename T>
void InitValue(AccPtr<T> &data, T value, size_t Size)
{
#ifdef CUDA
	int grid_size = ceil((float)Size/(float)INIT_VALUE_BLOCK_SIZE);
	cuda_kernel_init_value<T><<< grid_size, INIT_VALUE_BLOCK_SIZE, 0, data.getStream() >>>(
			~data,
			value,
			Size,
			INIT_VALUE_BLOCK_SIZE);
#else
	for (size_t i=0; i < Size; i++)
		data.hPtr[i] = value;
#endif
}

void centerFFT_2D(int grid_size, int batch_size, int block_size,
				cudaStream_t stream,
				XFLOAT *img_in,
				int image_size,
				int xdim,
				int ydim,
				int xshift,
				int yshift);


void centerFFT_2D(int grid_size, int batch_size, int block_size,
				XFLOAT *img_in,
				int image_size,
				int xdim,
				int ydim,
				int xshift,
				int yshift);

void centerFFT_3D(int grid_size, int batch_size, int block_size,
				cudaStream_t stream,
				XFLOAT *img_in,
				int image_size,
				int xdim,
				int ydim,
				int zdim,
				int xshift,
				int yshift,
				int zshift);


template<bool do_highpass>
void frequencyPass(int grid_size, int block_size,
				cudaStream_t stream,
				ACCCOMPLEX *A,
				long int ori_size,
				size_t Xdim,
				size_t Ydim,
				size_t Zdim,
				XFLOAT edge_low,
				XFLOAT edge_width,
				XFLOAT edge_high,
				XFLOAT angpix,
				int image_size)
{
#ifdef CUDA
	dim3 blocks(grid_size);
	cuda_kernel_frequencyPass<do_highpass><<<blocks,block_size, 0, stream>>>(
			A,
			ori_size,
			Xdim,
			Ydim,
			Zdim,
			edge_low,
			edge_width,
			edge_high,
			angpix,
			image_size);
#else
	CpuKernels::kernel_frequencyPass<do_highpass>(grid_size, block_size,
			A,
			ori_size,
			Xdim,
			Ydim,
			Zdim,
			edge_low,
			edge_width,
			edge_high,
			angpix,
			image_size);
#endif
}

template<bool REFCTF, bool REF3D, bool DATA3D, int block_sz>
void kernel_wavg(
		XFLOAT *g_eulers,
		AccProjectorKernel &projector,
		unsigned image_size,
		unsigned long orientation_num,
		XFLOAT *g_img_real,
		XFLOAT *g_img_imag,
		XFLOAT *g_trans_x,
		XFLOAT *g_trans_y,
		XFLOAT *g_trans_z,
		XFLOAT* g_weights,
		XFLOAT* g_ctfs,
		XFLOAT *g_wdiff2s_parts,
		XFLOAT *g_wdiff2s_AA,
		XFLOAT *g_wdiff2s_XA,
		unsigned long translation_num,
		XFLOAT weight_norm,
		XFLOAT significant_weight,
		XFLOAT part_scale,
		cudaStream_t stream)
{
#ifdef CUDA
	//We only want as many blocks as there are chunks of orientations to be treated
	//within the same block (this is done to reduce memory loads in the kernel).
	dim3 block_dim = orientation_num;//ceil((float)orientation_num/(float)REF_GROUP_SIZE);
	
	cuda_kernel_wavg<REFCTF,REF3D,DATA3D,block_sz><<<block_dim,block_sz,(3*block_sz+9)*sizeof(XFLOAT),stream>>>(
		g_eulers,
		projector,
		image_size,
		orientation_num,
		g_img_real,
		g_img_imag,
		g_trans_x,
		g_trans_y,
		g_trans_z,
		g_weights,
		g_ctfs,
		g_wdiff2s_parts,
		g_wdiff2s_AA,
		g_wdiff2s_XA,
		translation_num,
		weight_norm,
		significant_weight,
		part_scale);
#else
	if (DATA3D)
	{
		CpuKernels::wavg_3D<REFCTF>(
			g_eulers,
			projector,
			image_size,
			orientation_num,
			g_img_real,
			g_img_imag,
			g_trans_x,
			g_trans_y,
			g_trans_z,
			g_weights,
			g_ctfs,
			g_wdiff2s_parts,
			g_wdiff2s_AA,
			g_wdiff2s_XA,
			translation_num,
			weight_norm,
			significant_weight,
			part_scale);
	}
	else
	{
		CpuKernels::wavg_ref3D<REFCTF,REF3D>(
			g_eulers,
			projector,
			image_size,
			orientation_num,
			g_img_real,
			g_img_imag,
			g_trans_x,
			g_trans_y,
			g_trans_z,
			g_weights,
			g_ctfs,
			g_wdiff2s_parts,
			g_wdiff2s_AA,
			g_wdiff2s_XA,
			translation_num,
			weight_norm,
			significant_weight,
			part_scale);
	}
#endif
}

template<bool REF3D, bool DATA3D, int block_sz, int eulers_per_block, int prefetch_fraction>
void diff2_coarse(
		int grid_size, int block_size,
		XFLOAT *g_eulers,
		XFLOAT *trans_x,
		XFLOAT *trans_y,
		XFLOAT *trans_z,
		XFLOAT *g_real,
		XFLOAT *g_imag,
		AccProjectorKernel projector,
		XFLOAT *g_corr,
		XFLOAT *g_diff2s,
		int translation_num,
		int image_size,
		cudaStream_t stream
		)
{
#ifdef CUDA
		cuda_kernel_diff2_coarse<REF3D, DATA3D, block_sz, eulers_per_block, prefetch_fraction>
		<<<grid_size,block_size,0,stream>>>(
			g_eulers,
			trans_x,
			trans_y,
			trans_z,
			g_real,
			g_imag,
			projector,
			g_corr,
			g_diff2s,
			translation_num,
			image_size);
#else
	#if 1
		CpuKernels::diff2_coarse<REF3D, DATA3D, block_sz, eulers_per_block, prefetch_fraction>(
			grid_size,
			g_eulers,
			trans_x,
			trans_y,
			trans_z,
			g_real,
			g_imag,
			projector,
			g_corr,
			g_diff2s,
			translation_num,
			image_size
		);
	#else
		if (DATA3D)
			CpuKernels::diff2_coarse_3D<eulers_per_block>(
			grid_size,
			g_eulers,
			trans_x,
			trans_y,
			trans_z,
			g_real,
			g_imag,
			projector,
			g_corr,
			g_diff2s,
			translation_num,
			image_size);
		else
			CpuKernels::diff2_coarse_2D<REF3D, eulers_per_block>(
			grid_size,
			g_eulers,
			trans_x,
			trans_y,
			trans_z,
			g_real,
			g_imag,
			projector,
			g_corr,
			g_diff2s,
			translation_num,
			image_size);
	#endif
#endif
}

template<bool REF3D, bool DATA3D, int block_sz>
void diff2_CC_coarse(
		int grid_size, int block_size,
		XFLOAT *g_eulers,
		XFLOAT *g_imgs_real,
		XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x,
		XFLOAT *g_trans_y,
		XFLOAT *g_trans_z,
		AccProjectorKernel projector,
		XFLOAT *g_corr_img,
		XFLOAT *g_diff2s,
		unsigned translation_num,
		int image_size,
		XFLOAT exp_local_sqrtXi2,
		cudaStream_t stream
		)
{
#ifdef CUDA
	dim3 CCblocks(grid_size,translation_num);
	cuda_kernel_diff2_CC_coarse<REF3D,DATA3D,block_sz>
		<<<CCblocks,block_size,0,stream>>>(
			g_eulers,
			g_imgs_real,
			g_imgs_imag,
			g_trans_x,
			g_trans_y,
			g_trans_z,
			projector,
			g_corr_img,
			g_diff2s,
			translation_num,
			image_size,
			exp_local_sqrtXi2);
#else
	if (DATA3D)
		CpuKernels::diff2_CC_coarse_3D(
			grid_size,
			g_eulers,
			g_imgs_real,
			g_imgs_imag,
			g_trans_x,
			g_trans_y,
			g_trans_z,
			projector,
			g_corr_img,
			g_diff2s,
			translation_num,
			image_size,
			exp_local_sqrtXi2);
	else
		CpuKernels::diff2_CC_coarse_2D<REF3D>(
			grid_size,
			g_eulers,
			g_imgs_real,
			g_imgs_imag,
			g_trans_x,
			g_trans_y,
			projector,
			g_corr_img,
			g_diff2s,
			translation_num,
			image_size,
			exp_local_sqrtXi2);	
#endif	
}

template<bool REF3D, bool DATA3D, int block_sz, int chunk_sz>
void diff2_fine(
		int grid_size, int block_size,
		XFLOAT *g_eulers,
		XFLOAT *g_imgs_real,
		XFLOAT *g_imgs_imag,
		XFLOAT *trans_x,
		XFLOAT *trans_y,
		XFLOAT *trans_z,
		AccProjectorKernel projector,
		XFLOAT *g_corr_img,
		XFLOAT *g_diff2s,
		unsigned image_size,
		XFLOAT sum_init,
		unsigned long orientation_num,
		unsigned long translation_num,
		unsigned long todo_blocks,
		unsigned long *d_rot_idx,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num,
		cudaStream_t stream
		)
{
#ifdef CUDA
		dim3 block_dim = grid_size;
		cuda_kernel_diff2_fine<REF3D,DATA3D, block_sz, chunk_sz>
				<<<block_dim,block_size,0,stream>>>(
					g_eulers,
					g_imgs_real,
					g_imgs_imag,
					trans_x,
					trans_y,
					trans_z,
					projector,
					g_corr_img,    // in these non-CC kernels this is effectively an adjusted MinvSigma2
					g_diff2s,
					image_size,
					sum_init,
					orientation_num,
					translation_num,
					todo_blocks, //significant_num,
					d_rot_idx,
					d_trans_idx,
					d_job_idx,
					d_job_num);
#else
		// TODO - make use of orientation_num, translation_num,todo_blocks on
		// CPU side if CUDA starts to use
	if (DATA3D)
		CpuKernels::diff2_fine_3D(
			grid_size,
			g_eulers,
			g_imgs_real,
			g_imgs_imag,
			trans_x,
			trans_y,
			trans_z,
			projector,
			g_corr_img,    // in these non-CC kernels this is effectively an adjusted MinvSigma2
			g_diff2s,
			image_size,
			sum_init,
			d_rot_idx,
			d_trans_idx,
			d_job_idx,
			d_job_num);
	else
		CpuKernels::diff2_fine_2D<REF3D>(
			grid_size,
			g_eulers,
			g_imgs_real,
			g_imgs_imag,
			trans_x,
			trans_y,
			trans_z,
			projector,
			g_corr_img,    // in these non-CC kernels this is effectively an adjusted MinvSigma2
			g_diff2s,
			image_size,
			sum_init,
			d_rot_idx,
			d_trans_idx,
			d_job_idx,
			d_job_num);
#endif
}

template<bool REF3D, bool DATA3D, int block_sz,int chunk_sz>
void diff2_CC_fine(
		int grid_size, int block_size,
		XFLOAT *g_eulers,
		XFLOAT *g_imgs_real,
		XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x,
		XFLOAT *g_trans_y,
		XFLOAT *g_trans_z,
		AccProjectorKernel &projector,
		XFLOAT *g_corr_img,
		XFLOAT *g_diff2s,
		unsigned image_size,
		XFLOAT sum_init,
		XFLOAT exp_local_sqrtXi2,
		unsigned long orientation_num,
		unsigned long translation_num,
		unsigned long todo_blocks,
		unsigned long *d_rot_idx,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num,
		cudaStream_t stream
		)
{
#ifdef CUDA
	dim3 block_dim = grid_size;
	cuda_kernel_diff2_CC_fine<REF3D,DATA3D,block_sz,chunk_sz>
			<<<block_dim,block_size,0,stream>>>(
				g_eulers,
				g_imgs_real,
				g_imgs_imag,
				g_trans_x,
				g_trans_y,
				g_trans_z,
				projector,
				g_corr_img,
				g_diff2s,
				image_size,
				sum_init,
				exp_local_sqrtXi2,
				orientation_num,
				translation_num,
				todo_blocks,
				d_rot_idx,
				d_trans_idx,
				d_job_idx,
				d_job_num);
#else
		// TODO - Make use of orientation_num, translation_num, todo_blocks on
		// CPU side if CUDA starts to use
	if (DATA3D)
		CpuKernels::diff2_CC_fine_3D(
			grid_size,
			g_eulers,
			g_imgs_real,
			g_imgs_imag,
			g_trans_x,
			g_trans_y,
			g_trans_z,
			projector,
			g_corr_img,
			g_diff2s,
			image_size,
			sum_init,
			exp_local_sqrtXi2,
			d_rot_idx,
			d_trans_idx,
			d_job_idx,
			d_job_num);
	else
		CpuKernels::diff2_CC_fine_2D<REF3D>(
			grid_size,
			g_eulers,
			g_imgs_real,
			g_imgs_imag,
			g_trans_x,
			g_trans_y,
			projector,
			g_corr_img,
			g_diff2s,
			image_size,
			sum_init,
			exp_local_sqrtXi2,
			d_rot_idx,
			d_trans_idx,
			d_job_idx,
			d_job_num);
#endif
}

template<bool failsafe,typename weights_t>
void kernel_exponentiate_weights_coarse(
		int grid_size, 
		int num_classes,
		int block_size,
		XFLOAT *g_pdf_orientation,
		XFLOAT *g_pdf_offset,
		weights_t *g_Mweight,
		XFLOAT avg_diff2,
		XFLOAT min_diff2,
		int nr_coarse_orient,
		int nr_coarse_trans)
{
#ifdef CUDA
		dim3 block_dim(grid_size,num_classes);
		cuda_kernel_exponentiate_weights_coarse<failsafe,weights_t>
		<<<block_dim,block_size,0>>>(
				g_pdf_orientation,
				g_pdf_offset,
				g_Mweight,
				avg_diff2,
				min_diff2,
				nr_coarse_orient,
				nr_coarse_trans);
#else
//						for(int i=0; i<block_num; i++)
//						for(int j=0; j<sp.iclass_max-sp.iclass_min+1; j++)
//							 for(int k=0; k<SUMW_BLOCK_SIZE; k++)
		CpuKernels::exponentiate_weights_coarse<failsafe,weights_t>(
				grid_size, 
				num_classes, 
				block_size,
				g_pdf_orientation,
				g_pdf_offset,
				g_Mweight,
				avg_diff2,
				min_diff2,
				nr_coarse_orient,
				nr_coarse_trans);	
#endif
}

void kernel_exponentiate_weights_fine(	int grid_size, 
										int block_size,
										XFLOAT *g_pdf_orientation,
										XFLOAT *g_pdf_offset,
										XFLOAT *g_weights,
										XFLOAT avg_diff2,
										int oversamples_orient,
										int oversamples_trans,
										unsigned long *d_rot_id,
										unsigned long *d_trans_idx,
										unsigned long *d_job_idx,
										unsigned long *d_job_num,
										long int job_num,
										cudaStream_t stream);

};  // namespace AccUtilities


#endif //ACC_UTILITIES_H_

