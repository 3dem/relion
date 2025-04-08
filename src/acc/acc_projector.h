#ifndef ACC_PROJECTOR_H_
#define ACC_PROJECTOR_H_

#include "src/complex.h"
//#include "src/acc/cuda/cuda_settings.h"
//#include "src/acc/cuda/cuda_mem_utils.h"
#include "src/acc/acc_ptr.h"
//#include <cuda_runtime.h>
//#include "src/acc/cuda/cuda_kernels/cuda_device_utils.cuh"
#ifdef ALTCPU
#include <complex>
#elif _SYCL_ENABLED
#include "src/acc/sycl/sycl_virtual_dev.h"
using deviceStream_t = virtualSYCL*;
#endif

class AccProjector
{
	friend class AccProjectorKernel;

	int mdlX, mdlY, mdlZ, mdlMaxR,
	    mdlInitY, mdlInitZ;
	XFLOAT padding_factor;
	size_t mdlXYZ;

	size_t allocaton_size;

#ifdef _SYCL_ENABLED
	deviceStream_t devAcc;
#endif

#ifndef PROJECTOR_NO_TEXTURES

	XFLOAT *texArrayReal2D, *texArrayImag2D;
	#ifdef _CUDA_ENABLED
		cudaArray_t *texArrayReal, *texArrayImag;
		cudaTextureObject_t *mdlReal, *mdlImag;
	#elif _HIP_ENABLED
		hipArray_t *texArrayReal, *texArrayImag;
		hipTextureObject_t *mdlReal, *mdlImag;	
	#endif
	size_t pitch2D;
#else
#ifdef ALTCPU
	std::complex<XFLOAT> *mdlComplex;
	int externalFree;
#else
	XFLOAT *mdlComplex;
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
#ifdef _SYCL_ENABLED
		devAcc = nullptr;
#endif
#ifndef PROJECTOR_NO_TEXTURES

		texArrayReal2D = 0;
		texArrayImag2D = 0;
		texArrayReal = 0;
		texArrayImag = 0;
		mdlReal = 0;
		mdlImag = 0;
		pitch2D = 0;
#else
		mdlComplex = 0;
#ifdef ALTCPU
		externalFree = 0;
#endif
#endif
	}

	bool setMdlDim(
#ifdef _SYCL_ENABLED
			deviceStream_t dev,
#endif
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int maxr, XFLOAT paddingFactor);

#if defined(_CUDA_ENABLED) || defined(_HIP_ENABLED)
	void initMdl(XFLOAT *real, XFLOAT *imag);
#endif
#ifdef ALTCPU
	void initMdl(std::complex<XFLOAT> *data);
#else
	void initMdl(Complex *data);
#endif

	void clear();

	~AccProjector()
	{
		clear();
	};

};  // AccProjector

#endif
