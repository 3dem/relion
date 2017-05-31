#ifndef CUDA_DEVICE_UTILS_CUH_
#define CUDA_DEVICE_UTILS_CUH_

#include <cuda_runtime.h>
#include "src/gpu_utils/cuda_settings.h"

#ifdef CUDA_DOUBLE_PRECISION
__device__ inline double cuda_atomic_add(double* address, double val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do
	{
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
	}
	while (assumed != old); // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
	return __longlong_as_double(old);
}
#else
__device__ inline void cuda_atomic_add(float* address, float value)
{
  atomicAdd(address,value);
}
#endif


/*
 * For the following functions always use fast, low-precision intrinsics
 */

template< typename T1, typename T2 >
static inline
__device__ int floorfracf(T1 a, T2 b)
{
//	return __float2int_rd(__fdividef( (float)a, (float)b ) );
	return (int)(a/b);
}

template< typename T1, typename T2 >
static inline
__device__ int ceilfracf(T1 a, T2 b)
{
//	return __float2int_ru(__fdividef( (float)a, (float)b ) );
	return (int)(a/b + 1);
}

static inline
__device__ XFLOAT no_tex2D(XFLOAT* mdl, XFLOAT xp, XFLOAT yp, int mdlX, int mdlInitY)
{
	int x0 = floorf(xp);
	XFLOAT fx = xp - x0;
	int x1 = x0 + 1;

	int y0 = floorf(yp);
	XFLOAT fy = yp - y0;
	y0 -= mdlInitY;
	int y1 = y0 + 1;

	//-----------------------------
	XFLOAT d00 = mdl[y0*mdlX+x0];
	XFLOAT d01 = mdl[y0*mdlX+x1];
	XFLOAT d10 = mdl[y1*mdlX+x0];
	XFLOAT d11 = mdl[y1*mdlX+x1];
	//-----------------------------
	XFLOAT dx0 = d00 + (d01 - d00)*fx;
	XFLOAT dx1 = d10 + (d11 - d10)*fx;
	//-----------------------------

	return dx0 + (dx1 - dx0)*fy;
}

static inline
__device__ XFLOAT no_tex3D(XFLOAT* mdl, XFLOAT xp, XFLOAT yp, XFLOAT zp, int mdlX, int mdlXY, int mdlInitY, int mdlInitZ)
{
	int x0 = floorf(xp);
	XFLOAT fx = xp - x0;
	int x1 = x0 + 1;

	int y0 = floorf(yp);
	XFLOAT fy = yp - y0;
	y0 -= mdlInitY;
	int y1 = y0 + 1;

	int z0 = floorf(zp);
	XFLOAT fz = zp - z0;
	z0 -= mdlInitZ;
	int z1 = z0 + 1;

	XFLOAT d000 = mdl[z0*mdlXY+y0*mdlX+x0];
	XFLOAT d001 = mdl[z0*mdlXY+y0*mdlX+x1];
	XFLOAT d010 = mdl[z0*mdlXY+y1*mdlX+x0];
	XFLOAT d011 = mdl[z0*mdlXY+y1*mdlX+x1];
	XFLOAT d100 = mdl[z1*mdlXY+y0*mdlX+x0];
	XFLOAT d101 = mdl[z1*mdlXY+y0*mdlX+x1];
	XFLOAT d110 = mdl[z1*mdlXY+y1*mdlX+x0];
	XFLOAT d111 = mdl[z1*mdlXY+y1*mdlX+x1];
	//-----------------------------
	XFLOAT dx00 = d000 + (d001 - d000)*fx;
	XFLOAT dx01 = d100 + (d101 - d100)*fx;
	XFLOAT dx10 = d010 + (d011 - d010)*fx;
	XFLOAT dx11 = d110 + (d111 - d110)*fx;
	//-----------------------------
	XFLOAT dxy0 = dx00 + (dx10 - dx00)*fy;
	XFLOAT dxy1 = dx01 + (dx11 - dx01)*fy;
	//-----------------------------
	return dxy0 + (dxy1 - dxy0)*fz;
}

__device__ __forceinline__ void translatePixel(
		int x,
		int y,
		XFLOAT tx,
		XFLOAT ty,
		XFLOAT &real,
		XFLOAT &imag,
		XFLOAT &tReal,
		XFLOAT &tImag)
{
	XFLOAT s, c;
#ifdef CUDA_DOUBLE_PRECISION
	sincos( x * tx + y * ty , &s, &c );
#else
	sincosf( x * tx + y * ty , &s, &c );
#endif

	tReal = c * real - s * imag;
	tImag = c * imag + s * real;
}

__device__ __forceinline__ void translatePixel(
		int x,
		int y,
		int z,
		XFLOAT tx,
		XFLOAT ty,
		XFLOAT tz,
		XFLOAT &real,
		XFLOAT &imag,
		XFLOAT &tReal,
		XFLOAT &tImag)
{
	XFLOAT s, c;
#ifdef CUDA_DOUBLE_PRECISION
	sincos( x * tx + y * ty + z * tz, &s, &c );
#else
	sincosf( x * tx + y * ty + z * tz, &s, &c );
#endif

	tReal = c * real - s * imag;
	tImag = c * imag + s * real;
}

#endif
