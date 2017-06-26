#ifndef ACC_UTILITIES_H_
#define ACC_UTILITIES_H_

#include "src/acc/acc_ptr.h"
#include "src/acc/data_types.h"
#ifdef CUDA
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#else
#include "src/acc/cpu/cpu_kernels/helper.h"
#endif

namespace AccUtilities
{
	
template <typename T>
static void multiply(AccDataTypes::Image<T> &ptr, T value)
{
#ifdef CUDA
	int BSZ = ( (int) ceilf(( float)ptr.getSize() /(float)BLOCK_SIZE));
	CudaKernels::cuda_kernel_multi<T><<<BSZ,BLOCK_SIZE,0,ptr.getStream()>>>(
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
static void multiply(T *array, T value, size_t size, int MultiBsize, cudaStream_t stream)
{
#ifdef CUDA
	CudaKernels::cuda_kernel_multi<T><<<MultiBsize,BLOCK_SIZE,0,stream>>>(
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
static void translate(
	AccDataTypes::Image<T> &in,
	AccDataTypes::Image<T> &out,
	int dx,
	int dy,
	int dz=0)
{
#ifdef CUDA
int BSZ = ( (int) ceilf(( float)in.getxyz() /(float)BLOCK_SIZE));

if (in.is3D())
{
	CudaKernels::cuda_kernel_translate3D<T><<<BSZ,BLOCK_SIZE,0,in.getStream()>>>(
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
	CudaKernels::cuda_kernel_translate2D<T><<<BSZ,BLOCK_SIZE,0,in.getStream()>>>(
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
	size_t size = ptr.getSize();
	T sum;
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
	return CpuKernels::getMin<T>(ptr(), ptr.getSize());
#endif
}

template <typename T>
static std::pair<int, T> getArgMinOnDevice(AccPtr<T> &ptr)
{
#ifdef CUDA
	return CudaKernels::getArgMinOnDevice<T>(ptr);
#else
	return CpuKernels::getArgMin<T>(ptr(), ptr.getSize());
#endif
}

template <typename T>
static std::pair<int, T> getArgMaxOnDevice(AccPtr<T> &ptr)
{
#ifdef CUDA
	return CudaKernels::getArgMaxOnDevice<T>(ptr);
#else
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
	for(int i=0; i<arr_size; i++)
		if(in[i] > 0.0)
			filt_size++;
	out.resizeHost(filt_size);
	for(int i=0; i<arr_size; i++)
		if(in[i] > 0.0)
			out[i] = in[i];
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
		XFLOAT *vol,
		Image<RFLOAT> &img,
		bool     do_Mnoise,
		XFLOAT   radius,
		XFLOAT   radius_p,
		XFLOAT   cosine_width,
		XFLOAT  *g_sum,
		XFLOAT  *g_sum_bg,
		dim3 inblock_dim, 
		int inblock_size)
	{
#ifdef CUDA
		cuda_kernel_softMaskBackgroundValue<<<inblock_dim,inblock_size>>>(	vol,
																	img().nzyxdim,
																	img.data.xdim,
																	img.data.ydim,
																	img.data.zdim,
																	img.data.xdim/2,
																	img.data.ydim/2,
																	img.data.zdim/2, //unused
																	do_Mnoise,
																	radius,
																	radius_p,
																	cosine_width,
																	g_sum,
																	g_sum_bg);
#else
		int block_dim = inblock_dim;
		CpuKernels::softMaskBackgroundValue(block_dim, inblock_size,
			vol,
			img().nzyxdim,
			img.data.xdim,
			img.data.ydim,
			img.data.zdim,
			img.data.xdim/2,
			img.data.ydim/2,
			img.data.zdim/2, //unused
			do_Mnoise,
			radius,
			radius_p,
			cosine_width,
			g_sum,
			g_sum_bg);
#endif
	}
};  // namespace 


#endif //ACC_UTILITIES_H_

