#ifndef CUDA_PROJECTOR_H_
#define CUDA_PROJECTOR_H_

#include "src/complex.h"
#include "src/gpu_utils/cuda_settings.h"
#include <cuda_runtime.h>

class CudaProjector
{
	friend class CudaProjectorKernel;

	int mdlX, mdlY, mdlZ, mdlXYZ, mdlMaxR,
	    mdlInitY, mdlInitZ,
	    padding_factor;

#ifndef CUDA_DOUBLE_PRECISION
	cudaArray_t *texArrayReal, *texArrayImag;
	cudaTextureObject_t *mdlReal, *mdlImag;
#else
	double *mdlReal, *mdlImag;
#endif

public:
	CudaProjector():
			mdlX(0), mdlY(0), mdlZ(0),
			mdlXYZ(0), mdlMaxR(0),
			mdlInitY(0), mdlInitZ(0),
			padding_factor(0)
	{
#ifndef CUDA_DOUBLE_PRECISION
		texArrayReal = 0;
		texArrayImag = 0;
		mdlReal = 0;
		mdlImag = 0;
#else
		mdlReal = 0;
		mdlImag = 0;
#endif
	};

	CudaProjector(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int max_r, int padding_factor):
			mdlX(xdim), mdlY(ydim), mdlZ(zdim),
			mdlXYZ(xdim*ydim*zdim), mdlMaxR(max_r),
			mdlInitY(inity), mdlInitZ(initz),
			padding_factor(padding_factor)
	{
#ifndef CUDA_DOUBLE_PRECISION
		texArrayReal = 0;
		texArrayImag = 0;
		mdlReal = 0;
		mdlImag = 0;
#else
		mdlReal = 0;
		mdlImag = 0;
#endif
		if(zdim==1)
			mdlZ=0;
	};

	inline
	void setMdlDim(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int maxr, int paddingFactor)
	{
		mdlX = xdim;
		mdlY = ydim;
		if(zdim==1)
			mdlZ=0;
		else
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
