#ifndef CUDA_PROJECTOR_CUH_
#define CUDA_PROJECTOR_CUH_

#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_projector.h"
#include <cuda_runtime.h>


class Cuda3DProjectorKernel
{

public:
	int mdlX, mdlY, mdlZ,
		imgX, imgY,
		mdlInitY, mdlInitZ,
		padding_factor,
		maxR, maxR2;

#ifndef CUDA_DOUBLE_PRECISION
	cudaTextureObject_t mdlReal, mdlImag;
#else
	double *mdlReal, *mdlImag;
#endif

	Cuda3DProjectorKernel(
			int mdlX, int mdlY, int mdlZ,
			int imgX, int imgY,
			int mdlInitY, int mdlInitZ,
			int padding_factor,
			int maxR,
#ifndef CUDA_DOUBLE_PRECISION
			cudaTextureObject_t mdlReal, cudaTextureObject_t mdlImag
#else
			double *mdlReal, *mdlImag
#endif
			):
			mdlX(mdlX), mdlY(mdlY), mdlZ(mdlZ),
			imgX(imgX), imgY(imgY),
			mdlInitY(mdlInitY), mdlInitZ(mdlInitZ),
			padding_factor(padding_factor),
			maxR(maxR), maxR2(maxR*maxR),
			mdlReal(mdlReal), mdlImag(mdlImag)
		{};

	__device__ __forceinline__ void project(
			int pixel,
			FLOAT e0,
			FLOAT e1,
			FLOAT e3,
			FLOAT e4,
			FLOAT e6,
			FLOAT e7,
			FLOAT &real,
			FLOAT &imag)
	{
		bool is_neg_x;
		long int r2;

		int x = pixel % imgX;
		int y = (int)floorf( (float)pixel / (float)imgX);

		// Dont search beyond square with side max_r
		if (y > maxR)
		{
			if (y >= imgY - maxR)
				y = y - imgY;
			else
				x=maxR;
		}

		r2 = x*x + y*y;
		if (r2 <= maxR2)
		{
			FLOAT xp = (e0 * x + e1 * y ) * padding_factor;
			FLOAT yp = (e3 * x + e4 * y ) * padding_factor;
			FLOAT zp = (e6 * x + e7 * y ) * padding_factor;

			// Only asymmetric half is stored
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;
				is_neg_x = true;
			}
			else
			{
				is_neg_x = false;
			}
			yp -= mdlInitY;
			zp -= mdlInitZ;

			real = tex3D<FLOAT>(mdlReal, xp + 0.5f, yp + 0.5f, zp + 0.5f);
			imag = tex3D<FLOAT>(mdlImag, xp + 0.5f, yp + 0.5f, zp + 0.5f);

			if (is_neg_x)
			{
				imag = -imag;
			}
		}
		else
		{
			real=0.0f;
			imag=0.0f;
		}
	}




	static Cuda3DProjectorKernel makeKernel(Cuda3DProjector &p, int imgX, int imgY, int imgMaxR)
	{
		int maxR = p.mdlMaxR >= imgMaxR ? imgMaxR : p.mdlMaxR;

		Cuda3DProjectorKernel k(
				p.mdlX, p.mdlY, p.mdlZ,
		        imgX, imgY,
		        p.mdlInitY, p.mdlInitZ,
			    p.padding_factor,
			    maxR,

	#ifndef CUDA_DOUBLE_PRECISION
			    *(cudaTextureObject_t*) p.mdlReal,
			    *(cudaTextureObject_t*) p.mdlImag
	#else
				~(*(CudaGlobalPtr<double>*) p.mdlReal),
				~(*(CudaGlobalPtr<double>*) p.mdlImag)
	#endif

				);

		return k;
	}
};

#endif
