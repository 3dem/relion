#ifndef CUDA_PROJECTOR_H_
#define CUDA_PROJECTOR_H_

#ifdef __CUDACC__
#include <cuda_runtime.h>
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_utils.cuh"
#endif
#include "src/complex.h"

class Cuda3DProjector
{
public:
	int mdlX, mdlY, mdlZ, mdlXYZ, mdlMaxR,
	    mdlInitY, mdlInitZ,
	    imgX, imgY, maxR, maxR2,
	    padding_factor;

#ifdef __CUDACC__

#ifndef CUDA_DOUBLE_PRECISION
	cudaTextureObject_t mdlReal, mdlImag;
	cudaArray *mdlTextureArrayReal, *mdlTextureArrayImag;
#else
	CudaGlobalPtr<FLOAT> mdlReal, mdlImag;
#endif

#endif

	inline
	Cuda3DProjector(int xdim, int ydim, int zdim, int inity, int initz, int max_r, int padding_factor):
		mdlX(xdim), mdlY(ydim), mdlZ(zdim), mdlXYZ(xdim*ydim*zdim), mdlMaxR(max_r),
		mdlInitY(inity), mdlInitZ(initz),
		imgX(0), imgY(0), maxR(0), maxR2(0),
		padding_factor(padding_factor)
	{
#ifdef __CUDACC__
#ifndef CUDA_DOUBLE_PRECISION
		mdlReal = 0;
		mdlImag = 0;
		mdlTextureArrayReal = 0;
		mdlTextureArrayImag = 0;
#else
		mdlReal();
		mdlImag();
#endif
#endif
	};

#ifndef CUDA_DOUBLE_PRECISION
	void setData(float *real, float *imag);
#else
	void setData(double *real, double *imag);
#endif
	void setData(Complex *data);

	inline
	void setImgDim(int x, int y, int r)
	{
		imgX = x;
		imgY = y;

		maxR = mdlMaxR >= r ? r : mdlMaxR;
		maxR2 = maxR*maxR;
	};


#ifdef __CUDACC__

	~Cuda3DProjector()
	{
#ifdef CUDA_DOUBLE_PRECISION
		delete mdlReal;
		delete mdlImag;
#else
		if (mdlReal != 0)  cudaDestroyTextureObject(mdlReal);
		if (mdlImag != 0)  cudaDestroyTextureObject(mdlImag);
		if (mdlTextureArrayReal != 0) cudaFreeArray(mdlTextureArrayReal);
		if (mdlTextureArrayImag != 0) cudaFreeArray(mdlTextureArrayImag);
#endif
	}

	class Kernel
	{
		int mdlX, mdlY, mdlZ,
	        imgX, imgY,
		    mdlInitY, mdlInitZ,
		    padding_factor,
		    maxR, maxR2;

		cudaTextureObject_t mdlReal, mdlImag;

		__device__ inline void project(
				int x, int y,
				FLOAT e0,
				FLOAT e1,
				FLOAT e3,
				FLOAT e4,
				FLOAT e6,
				FLOAT e7,
				FLOAT &real,
				FLOAT &ima);
	};

#endif
};

#endif
