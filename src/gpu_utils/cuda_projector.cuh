#ifndef CUDA_PROJECTOR_CUH_
#define CUDA_PROJECTOR_CUH_

#include "src/gpu_utils/cuda_utils.cuh"
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_projector.h"
#include <cuda_runtime.h>


#ifndef PROJECTOR_SRC_TYPE
#ifndef CUDA_DOUBLE_PRECISION
#define PROJECTOR_SRC_TYPE cudaTextureObject_t
#else
#define PROJECTOR_SRC_TYPE double *
#endif
#endif

class Cuda3DProjectorKernel
{

public:
	int mdlX, mdlXY,
		imgX, imgY,
		mdlInitY, mdlInitZ,
		padding_factor,
		maxR, maxR2;

	PROJECTOR_SRC_TYPE mdlReal;
	PROJECTOR_SRC_TYPE mdlImag;

	Cuda3DProjectorKernel(
			int mdlX, int mdlY,
			int imgX, int imgY,
			int mdlInitY, int mdlInitZ,
			int padding_factor,
			int maxR,
			PROJECTOR_SRC_TYPE mdlReal, PROJECTOR_SRC_TYPE mdlImag
			):
			mdlX(mdlX), mdlXY(mdlX*mdlY),
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
		int r2;

		int x = pixel % imgX;
		int y = (int)floorf( (float)pixel / (float)imgX);

		// Dont search beyond square with side max_r
		if (y > maxR)
		{
			if (y >= imgY - maxR)
				y = y - imgY;
			else
				x = maxR;
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

#ifndef CUDA_DOUBLE_PRECISION

			yp -= mdlInitY;
			zp -= mdlInitZ;

			real = tex3D<FLOAT>(mdlReal, xp + 0.5f, yp + 0.5f, zp + 0.5f);
			imag = tex3D<FLOAT>(mdlImag, xp + 0.5f, yp + 0.5f, zp + 0.5f);

#else

			int x0 = floorf(xp);
			FLOAT fx = xp - x0;
			int x1 = x0 + 1;

			int y0 = floorf(yp);
			FLOAT fy = yp - y0;
			y0 -= mdlInitY;
			int y1 = y0 + 1;

			int z0 = floorf(zp);
			FLOAT fz = zp - z0;
			z0 -= mdlInitZ;
			int z1 = z0 + 1;

			FLOAT d000_real = mdlReal[z0*mdlXY+y0*mdlX+x0];
			FLOAT d001_real = mdlReal[z0*mdlXY+y0*mdlX+x1];
			FLOAT d010_real = mdlReal[z0*mdlXY+y1*mdlX+x0];
			FLOAT d011_real = mdlReal[z0*mdlXY+y1*mdlX+x1];
			FLOAT d100_real = mdlReal[z1*mdlXY+y0*mdlX+x0];
			FLOAT d101_real = mdlReal[z1*mdlXY+y0*mdlX+x1];
			FLOAT d110_real = mdlReal[z1*mdlXY+y1*mdlX+x0];
			FLOAT d111_real = mdlReal[z1*mdlXY+y1*mdlX+x1];

			FLOAT d000_imag = mdlImag[z0*mdlXY+y0*mdlX+x0];
			FLOAT d001_imag = mdlImag[z0*mdlXY+y0*mdlX+x1];
			FLOAT d010_imag = mdlImag[z0*mdlXY+y1*mdlX+x0];
			FLOAT d011_imag = mdlImag[z0*mdlXY+y1*mdlX+x1];
			FLOAT d100_imag = mdlImag[z1*mdlXY+y0*mdlX+x0];
			FLOAT d101_imag = mdlImag[z1*mdlXY+y0*mdlX+x1];
			FLOAT d110_imag = mdlImag[z1*mdlXY+y1*mdlX+x0];
			FLOAT d111_imag = mdlImag[z1*mdlXY+y1*mdlX+x1];

			// Set the interpolated value in the 2D output array
			FLOAT dx00_real = d000_real + (d001_real - d000_real)*fx;
			FLOAT dx01_real = d100_real + (d101_real - d100_real)*fx;
			FLOAT dx10_real = d010_real + (d011_real - d010_real)*fx;
			FLOAT dx11_real = d110_real + (d111_real - d110_real)*fx;

			FLOAT dx00_imag = d000_imag + (d001_imag - d000_imag)*fx;
			FLOAT dx01_imag = d100_imag + (d101_imag - d100_imag)*fx;
			FLOAT dx10_imag = d010_imag + (d011_imag - d010_imag)*fx;
			FLOAT dx11_imag = d110_imag + (d111_imag - d110_imag)*fx;
			//-----------------------------
			FLOAT dxy0_real = dx00_real + (dx10_real - dx00_real)*fy;
			FLOAT dxy1_real = dx01_real + (dx11_real - dx01_real)*fy;

			FLOAT dxy0_imag = dx00_imag + (dx10_imag - dx00_imag)*fy;
			FLOAT dxy1_imag = dx01_imag + (dx11_imag - dx01_imag)*fy;
			//-----------------------------
			real = dxy0_real + (dxy1_real - dxy0_real)*fz;
			imag = dxy0_imag + (dxy1_imag - dxy0_imag)*fz;
			//-----------------------------

#endif

			if (is_neg_x)
			{
				imag = -imag;
			}
		}
		else
		{
			real = 0.0f;
			imag = 0.0f;
		}
	}




	static Cuda3DProjectorKernel makeKernel(Cuda3DProjector &p, int imgX, int imgY, int imgMaxR)
	{
		int maxR = p.mdlMaxR >= imgMaxR ? imgMaxR : p.mdlMaxR;

		Cuda3DProjectorKernel k(
				p.mdlX, p.mdlY,
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
