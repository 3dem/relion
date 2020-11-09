#ifndef ACC_PROJECTOR_H_
#define ACC_PROJECTOR_H_

#include "src/complex.h"
//#include "src/acc/cuda/cuda_settings.h"
//#include "src/acc/cuda/cuda_mem_utils.h"
#include "src/acc/acc_ptr.h"
//#include <cuda_runtime.h>
//#include "src/acc/cuda/cuda_kernels/cuda_device_utils.cuh"
#ifndef _CUDA_ENABLED
#include <complex>
#endif

class AccProjector
{
	friend class AccProjectorKernel;

	int mdlX, mdlY, mdlZ, mdlMaxR,
	    mdlInitY, mdlInitZ;
	XFLOAT padding_factor;
	size_t mdlXYZ;

	size_t allocaton_size;

#ifndef PROJECTOR_NO_TEXTURES

	XFLOAT *texArrayReal2D, *texArrayImag2D;
	cudaArray_t *texArrayReal, *texArrayImag;
	cudaTextureObject_t *mdlReal, *mdlImag;

	size_t pitch2D;
#else
#ifdef _CUDA_ENABLED
	XFLOAT *mdlReal, *mdlImag;
#else
	std::complex<XFLOAT> *mdlComplex;
	int externalFree;
#endif
#endif  // PROJECTOR_NO_TEXTURES

public:
	AccProjector():
			mdlX(0), mdlY(0), mdlZ(0),
			mdlXYZ(0), mdlMaxR(0),
			mdlInitY(0), mdlInitZ(0),
			padding_factor(0),
			allocaton_size(0)
	{
#ifndef PROJECTOR_NO_TEXTURES

		texArrayReal2D = 0;
		texArrayImag2D = 0;
		texArrayReal = 0;
		texArrayImag = 0;
		mdlReal = 0;
		mdlImag = 0;
		pitch2D = 0;
#else
#ifdef _CUDA_ENABLED
		mdlReal = 0;
		mdlImag = 0;
#else
		mdlComplex = 0;
		externalFree = 0;
#endif
#endif
	}

	bool setMdlDim(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int maxr, XFLOAT paddingFactor);

	void initMdl(XFLOAT *real, XFLOAT *imag);
	void initMdl(Complex *data);
#ifndef _CUDA_ENABLED
void initMdl(std::complex<XFLOAT> *data);
#endif

	void clear();

	~AccProjector()
	{
		clear();
	};

};  // AccProjector

#endif
