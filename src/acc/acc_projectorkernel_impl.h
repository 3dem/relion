#ifndef ACC_PROJECTORKERNELIMPL_H_
#define ACC_PROJECTORKERNELIMPL_H_


#ifndef PROJECTOR_NO_TEXTURES
	#ifdef _CUDA_ENABLED
		#define PROJECTOR_PTR_TYPE cudaTextureObject_t
	#elif _HIP_ENABLED
		#define PROJECTOR_PTR_TYPE hipTextureObject_t
	#endif
#else
	#define PROJECTOR_PTR_TYPE XFLOAT *
#endif

class AccProjectorKernel
{

public:
	int mdlX, mdlXY, mdlZ,
		imgX, imgY, imgZ,
		mdlInitY, mdlInitZ,
		maxR, maxR2, maxR2_padded;
	XFLOAT 	padding_factor;

#if defined _CUDA_ENABLED || defined _HIP_ENABLED
	PROJECTOR_PTR_TYPE mdlReal;
	PROJECTOR_PTR_TYPE mdlImag;
#elif _SYCL_ENABLED
	PROJECTOR_PTR_TYPE mdlComplex;
#else
	std::complex<XFLOAT> *mdlComplex;
#endif

#if defined _CUDA_ENABLED || defined _HIP_ENABLED
	AccProjectorKernel(
			int mdlX, int mdlY, int mdlZ,
			int imgX, int imgY, int imgZ,
			int mdlInitY, int mdlInitZ,
			XFLOAT padding_factor,
			int maxR,
			PROJECTOR_PTR_TYPE mdlReal, PROJECTOR_PTR_TYPE mdlImag
			):
				mdlX(mdlX), mdlXY(mdlX*mdlY), mdlZ(mdlZ),
				imgX(imgX), imgY(imgY), imgZ(imgZ),
				mdlInitY(mdlInitY), mdlInitZ(mdlInitZ),
				padding_factor(padding_factor),
				maxR(maxR), maxR2(maxR*maxR), maxR2_padded(maxR*maxR*padding_factor*padding_factor),
				mdlReal(mdlReal), mdlImag(mdlImag)
		{};
#else
	AccProjectorKernel(
			int mdlX, int mdlY, int mdlZ,
			int imgX, int imgY, int imgZ,
			int mdlInitY, int mdlInitZ,
			XFLOAT padding_factor,
			int maxR,
#if _SYCL_ENABLED
			PROJECTOR_PTR_TYPE mdlComplex
#else
			std::complex<XFLOAT> *mdlComplex
#endif
			):
			mdlX(mdlX), mdlXY(mdlX*mdlY), mdlZ(mdlZ),
			imgX(imgX), imgY(imgY), imgZ(imgZ),
			mdlInitY(mdlInitY), mdlInitZ(mdlInitZ),
			padding_factor(padding_factor),
			maxR(maxR), maxR2(maxR*maxR), maxR2_padded(maxR*maxR*padding_factor*padding_factor),
			mdlComplex(mdlComplex)
		{};
#endif

#if defined _CUDA_ENABLED || defined _HIP_ENABLED
	__device__ __forceinline__
#else
	inline
#endif
	void project3Dmodel(
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
		XFLOAT xp = (e0 * x + e1 * y + e2 * z) * padding_factor;
		XFLOAT yp = (e3 * x + e4 * y + e5 * z) * padding_factor;
		XFLOAT zp = (e6 * x + e7 * y + e8 * z) * padding_factor;

		int r2 = xp*xp + yp*yp + zp*zp;

		if (r2 <= maxR2_padded)
		{

#ifdef PROJECTOR_NO_TEXTURES
			bool invers(xp < 0);
			if (invers)
			{
				xp = -xp;
				yp = -yp;
				zp = -zp;
			}

#if defined _CUDA_ENABLED || defined _HIP_ENABLED
			real =   no_tex3D(mdlReal, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
			imag = - no_tex3D(mdlImag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
#elif _SYCL_ENABLED
			syclKernels::no_tex3D(mdlComplex, real, imag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
			imag = -imag;
#else
			CpuKernels::complex3D(mdlComplex, real, imag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
#endif

			if(invers)
			    imag = -imag;


#else
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
#endif
		}
		else
		{
			real = (XFLOAT)0;
			imag = (XFLOAT)0;
		}
	}

#if defined _CUDA_ENABLED || defined _HIP_ENABLED
	__device__ __forceinline__
#else
	inline
#endif
	void project3Dmodel(
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
		XFLOAT xp = (e0 * x + e1 * y ) * padding_factor;
		XFLOAT yp = (e3 * x + e4 * y ) * padding_factor;
		XFLOAT zp = (e6 * x + e7 * y ) * padding_factor;

		int r2 = xp*xp + yp*yp + zp*zp;

		if (r2 <= maxR2_padded)
		{

#ifdef PROJECTOR_NO_TEXTURES
			bool invers(xp < 0);
			if (invers)
			{
				xp = -xp;
				yp = -yp;
				zp = -zp;
			}

	#if defined _CUDA_ENABLED || defined _HIP_ENABLED
			real = no_tex3D(mdlReal, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
			imag = no_tex3D(mdlImag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
	#elif _SYCL_ENABLED
			syclKernels::no_tex3D(mdlComplex, real, imag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
	#else
			CpuKernels::complex3D(mdlComplex, real, imag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
	#endif

			if(invers)
			    imag = -imag;
#else
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
#endif
		}
		else
		{
			real = (XFLOAT)0;
			imag = (XFLOAT)0;
		}
	}

#if defined _CUDA_ENABLED || defined _HIP_ENABLED
__device__ __forceinline__
#else
	inline
#endif
	void project2Dmodel(
				int x,
				int y,
				XFLOAT e0,
				XFLOAT e1,
				XFLOAT e3,
				XFLOAT e4,
				XFLOAT &real,
				XFLOAT &imag)
	{
		XFLOAT xp = (e0 * x + e1 * y ) * padding_factor;
		XFLOAT yp = (e3 * x + e4 * y ) * padding_factor;

		int r2 = xp*xp + yp*yp;

		if (r2 <= maxR2_padded)
		{
#ifdef PROJECTOR_NO_TEXTURES
			bool invers(xp < 0);
			if (invers)
			{
				xp = -xp;
				yp = -yp;
			}

	#if defined _CUDA_ENABLED || defined _HIP_ENABLED
			real = no_tex2D(mdlReal, xp, yp, mdlX, mdlInitY);
			imag = no_tex2D(mdlImag, xp, yp, mdlX, mdlInitY);
	#elif _SYCL_ENABLED
			syclKernels::no_tex2D(mdlComplex, real, imag, xp, yp, mdlX, mdlInitY);
	#else
			CpuKernels::complex2D(mdlComplex, real, imag, xp, yp, mdlX, mdlInitY);
	#endif

			if(invers)
			    imag = -imag;

#else
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
#endif
		}
		else
		{
			real=(XFLOAT)0;
			imag=(XFLOAT)0;
		}
	}

	static AccProjectorKernel makeKernel(AccProjector &p, int imgX, int imgY, int imgZ, int imgMaxR)
	{
		int maxR = p.mdlMaxR >= imgMaxR ? imgMaxR : p.mdlMaxR;

		AccProjectorKernel k(
					p.mdlX, p.mdlY, p.mdlZ,
					imgX, imgY, imgZ,
					p.mdlInitY, p.mdlInitZ,
					p.padding_factor,
					maxR,
#ifndef PROJECTOR_NO_TEXTURES
					*p.mdlReal,
					*p.mdlImag
#else
					p.mdlComplex
#endif
				);
		return k;
	}
};  // class AccProjectorKernel


#endif
