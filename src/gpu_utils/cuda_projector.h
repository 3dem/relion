#ifndef CUDA_PROJECTOR_H_
#define CUDA_PROJECTOR_H_

#include "src/complex.h"
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include <cuda_runtime.h>

class CudaProjector
{
	friend class CudaProjectorKernel;

	int mdlX, mdlY, mdlZ, mdlXYZ, mdlMaxR,
	    mdlInitY, mdlInitZ,
	    padding_factor;

	size_t allocaton_size;

#ifndef CUDA_NO_TEXTURES

#if(COMPLEXTEXTURE)
	XFLOAT *texArrayComplex2D;
	cudaArray_t *texArrayComplex;
	cudaTextureObject_t *mdlComplex;
#else
	XFLOAT *texArrayReal2D, *texArrayImag2D;
	cudaArray_t *texArrayReal, *texArrayImag;
	cudaTextureObject_t *mdlReal, *mdlImag;
#endif

	size_t pitch2D;
#else
	XFLOAT *mdlReal, *mdlImag;
#endif

public:
	CudaProjector():
			mdlX(0), mdlY(0), mdlZ(0),
			mdlXYZ(0), mdlMaxR(0),
			mdlInitY(0), mdlInitZ(0),
			padding_factor(0),
			allocaton_size(0)
	{
#ifndef CUDA_NO_TEXTURES

#if(COMPLEXTEXTURE)
		texArrayComplex2D = 0;
		texArrayComplex = 0;
		mdlComplex = 0;
#else
		texArrayReal2D = 0;
		texArrayImag2D = 0;
		texArrayReal = 0;
		texArrayImag = 0;
		mdlReal = 0;
		mdlImag = 0;
#endif
		pitch2D = 0;
#else
		mdlReal = 0;
		mdlImag = 0;
#endif
	}

	bool setMdlDim(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int maxr, int paddingFactor);

	void initMdl(XFLOAT *real, XFLOAT *imag);
	void initMdl(Complex *data);

	void clear();

	~CudaProjector()
	{
		clear();
	};

};

#endif
