#ifndef CUDA_PROJECTOR_H_
#define CUDA_PROJECTOR_H_

#include "src/complex.h"
#include "src/gpu_utils/cuda_settings.h"

class Cuda3DProjector
{
public:
	int mdlX, mdlY, mdlZ, mdlXYZ, mdlMaxR,
	    mdlInitY, mdlInitZ,
	    imgX, imgY, maxR, maxR2,
	    padding_factor;

#ifndef CUDA_DOUBLE_PRECISION
	void *mdlReal, *mdlImag; //cudaTextureObject_t
	void *texArrayReal, *texArrayImag; //cudaArray_t
#else
	void *mdlReal, *mdlImag; //CudaGlobalPtr<double>
#endif

	Cuda3DProjector():
			mdlX(0), mdlY(0), mdlZ(0),
			mdlXYZ(0), mdlMaxR(0),
			mdlInitY(0), mdlInitZ(0),
			imgX(0), imgY(0), maxR(0), maxR2(0),
			padding_factor(0),
			mdlReal(0), mdlImag(0)
	{
#ifndef CUDA_DOUBLE_PRECISION
		texArrayReal = 0;
		texArrayImag = 0;
#endif
	};

	Cuda3DProjector(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int max_r, int padding_factor):
			mdlX(xdim), mdlY(ydim), mdlZ(zdim),
			mdlXYZ(xdim*ydim*zdim), mdlMaxR(max_r),
			mdlInitY(inity), mdlInitZ(initz),
			imgX(0), imgY(0), maxR(0), maxR2(0),
			padding_factor(padding_factor),
			mdlReal(0), mdlImag(0)
	{
#ifndef CUDA_DOUBLE_PRECISION
		texArrayReal = 0;
		texArrayImag = 0;
#endif
	};

	void setMdlDim(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int maxr, int paddingFactor)
	{
		mdlX = xdim;
		mdlY = ydim;
		mdlZ = zdim;
		mdlXYZ = xdim*ydim*zdim;
		mdlInitY = inity;
		mdlInitZ = initz;
		mdlMaxR = maxr;
		padding_factor = paddingFactor;
	}

	void setMdlData(FLOAT *real, FLOAT *imag);
	void setMdlData(Complex *data);

	void setImgDim(int x, int y, int r)
	{
		imgX = x;
		imgY = y;

		maxR = mdlMaxR >= r ? r : mdlMaxR;
		maxR2 = maxR*maxR;
	};

	~Cuda3DProjector();

//	class Kernel
//	{
//		int mdlX, mdlY, mdlZ,
//	        imgX, imgY,
//		    mdlInitY, mdlInitZ,
//		    padding_factor,
//		    maxR, maxR2;
//
//		cudaTextureObject_t mdlReal, mdlImag;
//
//		__device__ inline void project(
//				int x, int y,
//				FLOAT e0,
//				FLOAT e1,
//				FLOAT e3,
//				FLOAT e4,
//				FLOAT e6,
//				FLOAT e7,
//				FLOAT &real,
//				FLOAT &ima);
//	};

};

#endif
