#ifndef CUDA_PROJECTOR_CUH_
#define CUDA_PROJECTOR_CUH_

#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_projector.h"
#include <cuda_runtime.h>


#ifndef PROJECTOR_SRC_TYPE
#ifndef CUDA_DOUBLE_PRECISION
#define PROJECTOR_PTR_TYPE cudaTextureObject_t
#else
#define PROJECTOR_PTR_TYPE double *
#endif
#endif

class CudaProjectorKernel
{

public:
	int mdlX, mdlXY, mdlZ,
		imgX, imgY,
		mdlInitY, mdlInitZ,
		padding_factor,
		maxR, maxR2;

	PROJECTOR_PTR_TYPE mdlReal;
	PROJECTOR_PTR_TYPE mdlImag;

	CudaProjectorKernel(
			int mdlX, int mdlY, int mdlZ,
			int imgX, int imgY,
			int mdlInitY, int mdlInitZ,
			int padding_factor,
			int maxR,
			PROJECTOR_PTR_TYPE mdlReal, PROJECTOR_PTR_TYPE mdlImag
			):
			mdlX(mdlX), mdlXY(mdlX*mdlY), mdlZ(mdlZ),
			imgX(imgX), imgY(imgY),
			mdlInitY(mdlInitY), mdlInitZ(mdlInitZ),
			padding_factor(padding_factor),
			maxR(maxR), maxR2(maxR*maxR),
			mdlReal(mdlReal), mdlImag(mdlImag)
		{};

	__device__ __forceinline__ void project3Dmodel(
			int pixel,
			XFLOAT e0,
			XFLOAT e1,
			XFLOAT e3,
			XFLOAT e4,
			XFLOAT e6,
			XFLOAT e7,
			XFLOAT &real,
			XFLOAT &imag)
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
			XFLOAT xp = (e0 * x + e1 * y ) * padding_factor;
			XFLOAT yp = (e3 * x + e4 * y ) * padding_factor;
			XFLOAT zp = (e6 * x + e7 * y ) * padding_factor;

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

			real = tex3D<XFLOAT>(mdlReal, xp + 0.5f, yp + 0.5f, zp + 0.5f);
			imag = tex3D<XFLOAT>(mdlImag, xp + 0.5f, yp + 0.5f, zp + 0.5f);

#else

			int x0 = floorf(xp);
			XFLOAT fx = xp - x0;
			int x1 = x0 + 1;

			int y0 = floorf(yp);
			XFLOAT fy = yp - y0;
			y0 -= mdlInitY;
			int y1 = y0 + 1;

			int z0 = floorf(zp);
			XFLOAT fz = zp - z0;
			z0 -= mdlInitZ;
			int z1 = z0 + 1;

			XFLOAT d000_real = mdlReal[z0*mdlXY+y0*mdlX+x0];
			XFLOAT d001_real = mdlReal[z0*mdlXY+y0*mdlX+x1];
			XFLOAT d010_real = mdlReal[z0*mdlXY+y1*mdlX+x0];
			XFLOAT d011_real = mdlReal[z0*mdlXY+y1*mdlX+x1];
			XFLOAT d100_real = mdlReal[z1*mdlXY+y0*mdlX+x0];
			XFLOAT d101_real = mdlReal[z1*mdlXY+y0*mdlX+x1];
			XFLOAT d110_real = mdlReal[z1*mdlXY+y1*mdlX+x0];
			XFLOAT d111_real = mdlReal[z1*mdlXY+y1*mdlX+x1];

			XFLOAT d000_imag = mdlImag[z0*mdlXY+y0*mdlX+x0];
			XFLOAT d001_imag = mdlImag[z0*mdlXY+y0*mdlX+x1];
			XFLOAT d010_imag = mdlImag[z0*mdlXY+y1*mdlX+x0];
			XFLOAT d011_imag = mdlImag[z0*mdlXY+y1*mdlX+x1];
			XFLOAT d100_imag = mdlImag[z1*mdlXY+y0*mdlX+x0];
			XFLOAT d101_imag = mdlImag[z1*mdlXY+y0*mdlX+x1];
			XFLOAT d110_imag = mdlImag[z1*mdlXY+y1*mdlX+x0];
			XFLOAT d111_imag = mdlImag[z1*mdlXY+y1*mdlX+x1];

			// Set the interpolated value in the 2D output array
			XFLOAT dx00_real = d000_real + (d001_real - d000_real)*fx;
			XFLOAT dx01_real = d100_real + (d101_real - d100_real)*fx;
			XFLOAT dx10_real = d010_real + (d011_real - d010_real)*fx;
			XFLOAT dx11_real = d110_real + (d111_real - d110_real)*fx;

			XFLOAT dx00_imag = d000_imag + (d001_imag - d000_imag)*fx;
			XFLOAT dx01_imag = d100_imag + (d101_imag - d100_imag)*fx;
			XFLOAT dx10_imag = d010_imag + (d011_imag - d010_imag)*fx;
			XFLOAT dx11_imag = d110_imag + (d111_imag - d110_imag)*fx;
			//-----------------------------
			XFLOAT dxy0_real = dx00_real + (dx10_real - dx00_real)*fy;
			XFLOAT dxy1_real = dx01_real + (dx11_real - dx01_real)*fy;

			XFLOAT dxy0_imag = dx00_imag + (dx10_imag - dx00_imag)*fy;
			XFLOAT dxy1_imag = dx01_imag + (dx11_imag - dx01_imag)*fy;
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

	__device__ __forceinline__ void project2Dmodel(
				int pixel,
				XFLOAT e0,
				XFLOAT e1,
				XFLOAT e3,
				XFLOAT e4,
				XFLOAT e6,
				XFLOAT e7,
				XFLOAT &real,
				XFLOAT &imag)
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
				XFLOAT xp = (e0 * x + e1 * y ) * padding_factor;
				XFLOAT yp = (e3 * x + e4 * y ) * padding_factor;

				// Only asymmetric half is stored
				if (xp < 0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;
					is_neg_x = true;
				}
				else
				{
					is_neg_x = false;
				}

	#ifndef CUDA_DOUBLE_PRECISION

				yp -= mdlInitY;

				real = tex2D<XFLOAT>(mdlReal, xp + 0.5f, yp + 0.5f);
				imag = tex2D<XFLOAT>(mdlImag, xp + 0.5f, yp + 0.5f);

	#else
				int x0 = floorf(xp);
				XFLOAT fx = xp - x0;
				int x1 = x0 + 1;

				int y0 = floorf(yp);
				XFLOAT fy = yp - y0;
				y0 -= mdlInitY;
				int y1 = y0 + 1;
				//-----------------------------
				XFLOAT d00_real = mdlReal[y0*mdlX+x0];
				XFLOAT d01_real = mdlReal[y0*mdlX+x1];
				XFLOAT d10_real = mdlReal[y1*mdlX+x0];
				XFLOAT d11_real = mdlReal[y1*mdlX+x1];

				XFLOAT d00_imag = mdlImag[y0*mdlX+x0];
				XFLOAT d01_imag = mdlImag[y0*mdlX+x1];
				XFLOAT d10_imag = mdlImag[y1*mdlX+x0];
				XFLOAT d11_imag = mdlImag[y1*mdlX+x1];
				//-----------------------------
				// Set the interpolated value in the 2D output array
				XFLOAT dx0_real = d00_real + (d01_real - d00_real)*fx;
				XFLOAT dx1_real = d10_real + (d11_real - d10_real)*fx;

				XFLOAT dx0_imag = d00_imag + (d01_imag - d00_imag)*fx;
				XFLOAT dx1_imag = d10_imag + (d11_imag - d10_imag)*fx;
				//-----------------------------
				real = dx0_real + (dx1_real - dx0_real)*fy;
				imag = dx0_imag + (dx1_imag - dx0_imag)*fy;
				//-----------------------------
	#endif

//				printf("%f, %f",real,imag);
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


	static CudaProjectorKernel makeKernel(CudaProjector &p, int imgX, int imgY, int imgMaxR)
	{
		int maxR = p.mdlMaxR >= imgMaxR ? imgMaxR : p.mdlMaxR;

		CudaProjectorKernel k(
				p.mdlX, p.mdlY, p.mdlZ,
		        imgX, imgY,
		        p.mdlInitY, p.mdlInitZ,
			    p.padding_factor,
			    maxR,

#ifndef CUDA_DOUBLE_PRECISION
			    *p.mdlReal,
			    *p.mdlImag
#else
				p.mdlReal,
				p.mdlImag
#endif

				);

		return k;
	}
};

#endif
