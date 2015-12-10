#ifndef CUDA_PROJECTOR_CUH_
#define CUDA_PROJECTOR_CUH_

#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_device_utils.cuh"
#include <cuda_runtime.h>


#ifndef CUDA_NO_TEXTURES
#define PROJECTOR_PTR_TYPE cudaTextureObject_t
#else
#define PROJECTOR_PTR_TYPE XFLOAT *
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
	PROJECTOR_PTR_TYPE mdlComplex;

	CudaProjectorKernel(
			int mdlX, int mdlY, int mdlZ,
			int imgX, int imgY,
			int mdlInitY, int mdlInitZ,
			int padding_factor,
			int maxR,
			PROJECTOR_PTR_TYPE mdlComplex
			):
			mdlX(mdlX), mdlXY(mdlX*mdlY), mdlZ(mdlZ),
			imgX(imgX), imgY(imgY),
			mdlInitY(mdlInitY), mdlInitZ(mdlInitZ),
			padding_factor(padding_factor),
			maxR(maxR), maxR2(maxR*maxR),
			mdlComplex(mdlComplex)
		{};

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
			int x,
			int y,
			XFLOAT e0,
			XFLOAT e1,
			XFLOAT e3,
			XFLOAT e4,
			XFLOAT e6,
			XFLOAT e7,
			XFLOAT &real,
			XFLOAT &imag)
	{
		int r2;

		r2 = x*x + y*y;
		if (r2 <= maxR2)
		{
			XFLOAT xp = (e0 * x + e1 * y ) * padding_factor;
			XFLOAT yp = (e3 * x + e4 * y ) * padding_factor;
			XFLOAT zp = (e6 * x + e7 * y ) * padding_factor;

#ifndef CUDA_NO_TEXTURES
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;

				yp -= mdlInitY;
				zp -= mdlInitZ;

				real =    tex3D<XFLOAT>(mdlReal, xp + 0.5f, yp + 0.5f, zp + 0.5f);
				imag =  - tex3D<XFLOAT>(mdlImag, xp + 0.5f, yp + 0.5f, zp + 0.5f);
			}
			else
			{
				yp -= mdlInitY;
				zp -= mdlInitZ;

				real =   tex3D<XFLOAT>(mdlReal, xp + 0.5f, yp + 0.5f, zp + 0.5f);
				imag =   tex3D<XFLOAT>(mdlImag, xp + 0.5f, yp + 0.5f, zp + 0.5f);
			}
#else
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;

				real =   no_tex3D(mdlReal, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
				imag = - no_tex3D(mdlImag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
			}
			else
			{
				real =   no_tex3D(mdlReal, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
				imag =   no_tex3D(mdlImag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
			}
#endif
		}
		else
		{
			real = 0.0f;
			imag = 0.0f;
		}
	}

	__device__ __forceinline__ void project2Dmodel(
				int x,
				int y,
				XFLOAT e0,
				XFLOAT e1,
				XFLOAT e3,
				XFLOAT e4,
				XFLOAT e6,
				XFLOAT e7,
				XFLOAT &real,
				XFLOAT &imag)
		{
			int r2;

			r2 = x*x + y*y;
			if (r2 <= maxR2)
			{
				XFLOAT xp = (e0 * x + e1 * y ) * padding_factor;
				XFLOAT yp = (e3 * x + e4 * y ) * padding_factor;

#ifndef CUDA_NO_TEXTURES
				if (xp < 0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;
					yp -= mdlInitY;

					real =   tex2D<XFLOAT>(mdlReal, xp + 0.5f, yp + 0.5f);
					imag = - tex2D<XFLOAT>(mdlImag, xp + 0.5f, yp + 0.5f);
				}
				else
				{
					yp -= mdlInitY;
					real =   tex2D<XFLOAT>(mdlReal, xp + 0.5f, yp + 0.5f);
					imag =   tex2D<XFLOAT>(mdlImag, xp + 0.5f, yp + 0.5f);
				}
#else
				if (xp < 0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;

					real =   no_tex2D(mdlReal, xp, yp, mdlX, mdlInitY);
					imag = - no_tex2D(mdlImag, xp, yp, mdlX, mdlInitY);
				}
				else
				{
					real =   no_tex2D(mdlReal, xp, yp, mdlX, mdlInitY);
					imag =   no_tex2D(mdlImag, xp, yp, mdlX, mdlInitY);
				}
#endif
			}
			else
			{
				real = 0.0f;
				imag = 0.0f;
			}
		}

	__device__ __forceinline__ void project2Dmodel(
				int x,
				int y,
				XFLOAT e0,
				XFLOAT e1,
				XFLOAT e3,
				XFLOAT e4,
				XFLOAT &real,
				XFLOAT &imag)
		{
			int r2;

			r2 = x*x + y*y;
			if (r2 <= maxR2)
			{
				XFLOAT xp = (e0 * x + e1 * y ) * padding_factor;
				XFLOAT yp = (e3 * x + e4 * y ) * padding_factor;

#ifndef CUDA_NO_TEXTURES
				if (xp < 0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;
					yp -= mdlInitY;

					real =   tex2D<XFLOAT>(mdlReal, xp + 0.5f, yp + 0.5f);
					imag = - tex2D<XFLOAT>(mdlImag, xp + 0.5f, yp + 0.5f);
				}
				else
				{
					yp -= mdlInitY;
					real =   tex2D<XFLOAT>(mdlReal, xp + 0.5f, yp + 0.5f);
					imag =   tex2D<XFLOAT>(mdlImag, xp + 0.5f, yp + 0.5f);
				}
#else
				if (xp < 0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;

					real =   no_tex2D(mdlReal, xp, yp, mdlX, mdlInitY);
					imag = - no_tex2D(mdlImag, xp, yp, mdlX, mdlInitY);
				}
				else
				{
					real =   no_tex2D(mdlReal, xp, yp, mdlX, mdlInitY);
					imag =   no_tex2D(mdlImag, xp, yp, mdlX, mdlInitY);
				}
#endif
			}
			else
			{
				real = 0.0f;
				imag = 0.0f;
			}
		}

#if(COMPLEXTEXTURE)
	__device__ __forceinline__ void project2DComplexModel(
				int x,
				int y,
				XFLOAT e0,
				XFLOAT e1,
				XFLOAT e3,
				XFLOAT e4,
				CUDACOMPLEX &val)
		{
			int r2;

			r2 = x*x + y*y;
			if (r2 <= maxR2)
			{
				XFLOAT xp = (e0 * x + e1 * y ) * padding_factor;
				XFLOAT yp = (e3 * x + e4 * y ) * padding_factor;

				if (xp < 0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp = -xp;
					yp = -yp;
					yp -= mdlInitY;

					val =   tex2D<CUDACOMPLEX>(mdlComplex, xp + 0.5f, yp + 0.5f);
				}
				else
				{
					yp -= mdlInitY;
					val =   tex2D<CUDACOMPLEX>(mdlComplex, xp + 0.5f, yp + 0.5f);
				}
			}
			else
			{
				val.x=(XFLOAT)0;
				val.y=(XFLOAT)0;
			}
		}
#endif

	static CudaProjectorKernel makeKernel(CudaProjector &p, int imgX, int imgY, int imgMaxR)
	{
		int maxR = p.mdlMaxR >= imgMaxR ? imgMaxR : p.mdlMaxR;

		CudaProjectorKernel k(
					p.mdlX, p.mdlY, p.mdlZ,
					imgX, imgY,
					p.mdlInitY, p.mdlInitZ,
					p.padding_factor,
					maxR,
#if(COMPLEXTEXTURE)
					*p.mdlComplex
#else
#ifndef CUDA_NO_TEXTURES
					*p.mdlReal,
					*p.mdlImag
#else
					p.mdlReal,
					p.mdlImag
#endif
#endif
				);
		return k;
	}
};

#endif
