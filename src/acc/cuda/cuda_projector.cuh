#ifndef CUDA_PROJECTOR_CUH_
#define CUDA_PROJECTOR_CUH_

#include <cuda_runtime.h>
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_device_utils.cuh"



#ifndef CUDA_NO_TEXTURES
#define PROJECTOR_PTR_TYPE cudaTextureObject_t
#else
#define PROJECTOR_PTR_TYPE XFLOAT *
#endif

class CudaProjectorKernel
{

public:
	int mdlX, mdlXY, mdlZ,
		imgX, imgY, imgZ,
		mdlInitY, mdlInitZ,
		padding_factor,
		maxR, maxR2;

	PROJECTOR_PTR_TYPE mdlReal;
	PROJECTOR_PTR_TYPE mdlImag;
	PROJECTOR_PTR_TYPE mdlComplex;

	CudaProjectorKernel(
			int mdlX, int mdlY, int mdlZ,
			int imgX, int imgY, int imgZ,
			int mdlInitY, int mdlInitZ,
			int padding_factor,
			int maxR,
			PROJECTOR_PTR_TYPE mdlComplex
			):
			mdlX(mdlX), mdlXY(mdlX*mdlY), mdlZ(mdlZ),
			imgX(imgX), imgY(imgY), imgZ(imgZ),
			mdlInitY(mdlInitY), mdlInitZ(mdlInitZ),
			padding_factor(padding_factor),
			maxR(maxR), maxR2(maxR*maxR),
			mdlComplex(mdlComplex)
		{};

	CudaProjectorKernel(
			int mdlX, int mdlY, int mdlZ,
			int imgX, int imgY, int imgZ,
			int mdlInitY, int mdlInitZ,
			int padding_factor,
			int maxR,
			PROJECTOR_PTR_TYPE mdlReal, PROJECTOR_PTR_TYPE mdlImag
			):
				mdlX(mdlX), mdlXY(mdlX*mdlY), mdlZ(mdlZ),
				imgX(imgX), imgY(imgY), imgZ(imgZ),
				mdlInitY(mdlInitY), mdlInitZ(mdlInitZ),
				padding_factor(padding_factor),
				maxR(maxR), maxR2(maxR*maxR),
				mdlReal(mdlReal), mdlImag(mdlImag)
			{};

	__device__ __forceinline__ void project3Dmodel(
			int x,
			int y,
			int z,
			XFLOAT e0,
			XFLOAT e1,
			XFLOAT e2,
			XFLOAT e3,
			XFLOAT e4,
			XFLOAT e5,
			XFLOAT e6,
			XFLOAT e7,
			XFLOAT e8,
			XFLOAT &real,
			XFLOAT &imag)
	{
		int r2;

		r2 = x*x + y*y + z*z;
		if (r2 <= maxR2)
		{
			XFLOAT xp = (e0 * x + e1 * y + e2 * z ) * padding_factor;
			XFLOAT yp = (e3 * x + e4 * y + e5 * z ) * padding_factor;
			XFLOAT zp = (e6 * x + e7 * y + e8 * z ) * padding_factor;

#ifdef CUDA_NO_TEXTURES
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
#else
	#if(!COMPLEXTEXTURE)
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;

				yp -= mdlInitY;
				zp -= mdlInitZ;

				real =    tex3D<XFLOAT>(mdlReal, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
				imag =  - tex3D<XFLOAT>(mdlImag, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
			}
			else
			{
				yp -= mdlInitY;
				zp -= mdlInitZ;

				real =   tex3D<XFLOAT>(mdlReal, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
				imag =   tex3D<XFLOAT>(mdlImag, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
			}
	#else
			CUDACOMPLEX val;
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;

				yp -= mdlInitY;
				zp -= mdlInitZ;

				val =   tex3D<CUDACOMPLEX>(mdlComplex, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
				val.y = -val.y;
			}
			else
			{
				yp -= mdlInitY;
				zp -= mdlInitZ;

				val =   tex3D<CUDACOMPLEX>(mdlComplex, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
			}
			real=val.x;
			imag=val.y;
	#endif
#endif
		}
		else
		{
			real = 0.0f;
			imag = 0.0f;
		}
	}

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

#ifdef CUDA_NO_TEXTURES
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
#else
	#if(!COMPLEXTEXTURE)
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;

				yp -= mdlInitY;
				zp -= mdlInitZ;

				real =    tex3D<XFLOAT>(mdlReal, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
				imag =  - tex3D<XFLOAT>(mdlImag, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
			}
			else
			{
				yp -= mdlInitY;
				zp -= mdlInitZ;

				real =   tex3D<XFLOAT>(mdlReal, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
				imag =   tex3D<XFLOAT>(mdlImag, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
			}
	#else
			CUDACOMPLEX val;
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;

				yp -= mdlInitY;
				zp -= mdlInitZ;

				val =   tex3D<CUDACOMPLEX>(mdlComplex, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
				val.y = -val.y;
			}
			else
			{
				yp -= mdlInitY;
				zp -= mdlInitZ;

				val =   tex3D<CUDACOMPLEX>(mdlComplex, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
			}
			real=val.x;
			imag=val.y;
	#endif
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
#ifdef CUDA_NO_TEXTURES
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
#else
#if(!COMPLEXTEXTURE)
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				yp -= mdlInitY;

				real =   tex2D<XFLOAT>(mdlReal, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5);
				imag = - tex2D<XFLOAT>(mdlImag, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5);
			}
			else
			{
				yp -= mdlInitY;
				real =   tex2D<XFLOAT>(mdlReal, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5);
				imag =   tex2D<XFLOAT>(mdlImag, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5);
			}
#else
			CUDACOMPLEX val;
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				yp -= mdlInitY;

				val = tex2D<CUDACOMPLEX>(mdlComplex, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5);
				val.y = -val.y;
			}
			else
			{
				yp -= mdlInitY;
				val = tex2D<CUDACOMPLEX>(mdlComplex, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5);
			}
			real=val.x;
			imag=val.y;
#endif
#endif
		}
		else
		{
			real=(XFLOAT)0;
			imag=(XFLOAT)0;
		}
	}

	static CudaProjectorKernel makeKernel(CudaProjector &p, int imgX, int imgY, int imgZ, int imgMaxR)
	{
		int maxR = p.mdlMaxR >= imgMaxR ? imgMaxR : p.mdlMaxR;

		CudaProjectorKernel k(
					p.mdlX, p.mdlY, p.mdlZ,
					imgX, imgY, imgZ,
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
