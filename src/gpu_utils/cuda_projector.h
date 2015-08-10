#ifndef CUDA_PROJECTOR_H_
#define CUDA_PROJECTOR_H_

#include "src/complex.h"
#include "src/gpu_utils/cuda_settings.h"

#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif

class CudaProjector
{
	friend class CudaProjectorKernel;

	int mdlX, mdlY, mdlZ, mdlXYZ, mdlMaxR,
	    mdlInitY, mdlInitZ,
	    padding_factor;

	void *mdlReal, *mdlImag; //cudaTextureObject_t or CudaGlobalPtr<double>

#ifndef CUDA_DOUBLE_PRECISION
	void *texArrayReal, *texArrayImag; //cudaArray_t
#endif

public:
	CudaProjector():
			mdlX(0), mdlY(0), mdlZ(0),
			mdlXYZ(0), mdlMaxR(0),
			mdlInitY(0), mdlInitZ(0),
			padding_factor(0),
			mdlReal(0), mdlImag(0)
	{
#ifndef CUDA_DOUBLE_PRECISION
		texArrayReal = 0;
		texArrayImag = 0;
#endif
	};

	CudaProjector(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int max_r, int padding_factor):
			mdlX(xdim), mdlY(ydim), mdlZ(zdim),
			mdlXYZ(xdim*ydim*zdim), mdlMaxR(max_r),
			mdlInitY(inity), mdlInitZ(initz),
			padding_factor(padding_factor),
			mdlReal(0), mdlImag(0)
	{
#ifndef CUDA_DOUBLE_PRECISION
		texArrayReal = 0;
		texArrayImag = 0;
#endif
	};

	inline
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

	void setMdlData(XFLOAT *real, XFLOAT *imag);
	void setMdlData(Complex *data);

	~CudaProjector();

};

#endif
