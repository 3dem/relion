#ifndef ACC_PROJECTOR_H_
#define ACC_PROJECTOR_H_

#include "src/complex.h"
//#include "src/acc/cuda/cuda_settings.h"
//#include "src/acc/cuda/cuda_mem_utils.h"
#include "src/acc/acc_ptr.h"
//#include <cuda_runtime.h>
//#include "src/acc/cuda/cuda_kernels/cuda_device_utils.cuh"
#ifndef CUDA
#include <complex>
#endif

class AccProjector
{
	friend class AccProjectorKernel;

	int mdlX, mdlY, mdlZ, mdlXYZ, mdlMaxR,
	    mdlInitY, mdlInitZ,
	    padding_factor;

	size_t allocaton_size;

#ifndef PROJECTOR_NO_TEXTURES

	XFLOAT *texArrayReal2D, *texArrayImag2D;
	cudaArray_t *texArrayReal, *texArrayImag;
	cudaTextureObject_t *mdlReal, *mdlImag;

	size_t pitch2D;
#else
#ifdef CUDA
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
#ifdef CUDA
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
			int maxr, int paddingFactor);

	void initMdl(XFLOAT *real, XFLOAT *imag);
	void initMdl(Complex *data);
	void initMdl(std::complex<XFLOAT> *data);

	void clear();

	~AccProjector()
	{
		clear();
	};

};  // AccProjector

#endif
