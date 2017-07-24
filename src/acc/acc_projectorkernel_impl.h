#ifndef ACC_PROJECTORKERNELIMPL_H_
#define ACC_PROJECTORKERNELIMPL_H_


#ifndef CUDA_NO_TEXTURES
#define PROJECTOR_PTR_TYPE cudaTextureObject_t
#else
#define PROJECTOR_PTR_TYPE XFLOAT *
#endif

class AccProjectorKernel
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

	AccProjectorKernel(
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

	AccProjectorKernel(
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
			{
#ifndef CUDA			
#ifdef CUDA_NO_TEXTURES
				for(int i=0; i<mdlX * mdlY * mdlZ; i++) {
			        *mdlComplex ++ = *mdlReal ++;
			        *mdlComplex ++ = *mdlImag ++;
			    }
#endif
#endif
			};

#ifdef CUDA
	__device__ __forceinline__ void project3Dmodel(
#else
	void project3Dmodel(
#endif
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
		
        real=(XFLOAT)0;
		imag=(XFLOAT)0;

		r2 = x*x + y*y + z*z;
		if (r2 <= maxR2)
		{
			XFLOAT xp = (e0 * x + e1 * y + e2 * z ) * padding_factor;
			XFLOAT yp = (e3 * x + e4 * y + e5 * z ) * padding_factor;
			XFLOAT zp = (e6 * x + e7 * y + e8 * z ) * padding_factor;

#ifdef CUDA_NO_TEXTURES
			bool invers(xp < 0);
			if (invers)
			{
				xp = -xp;
				yp = -yp;
				zp = -zp;
			}

#ifdef CUDA
			real =   no_tex3D(mdlReal, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
			imag = - no_tex3D(mdlImag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
#else
			CpuKernels::complex3D(mdlComplex, real, imag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);	
#endif
			
			if(invers)
			    imag = -imag;

			
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
			ACCCOMPLEX val;
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;

				yp -= mdlInitY;
				zp -= mdlInitZ;

				val =   tex3D<ACCCOMPLEX>(mdlComplex, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
				val.y = -val.y;
			}
			else
			{
				yp -= mdlInitY;
				zp -= mdlInitZ;

				val =   tex3D<ACCCOMPLEX>(mdlComplex, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
			}
			real=val.x;
			imag=val.y;
	#endif
#endif
		}
		else
		{
			real = (XFLOAT)0;
			imag = (XFLOAT)0;
		}
	}

#ifdef CUDA
	__device__ __forceinline__ void project3Dmodel(
#else
	void project3Dmodel(
#endif
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
		
        real=(XFLOAT)0;
		imag=(XFLOAT)0;

		r2 = x*x + y*y;
		if (r2 <= maxR2)
		{
			XFLOAT xp = (e0 * x + e1 * y ) * padding_factor;
			XFLOAT yp = (e3 * x + e4 * y ) * padding_factor;
			XFLOAT zp = (e6 * x + e7 * y ) * padding_factor;

#ifdef CUDA_NO_TEXTURES
			bool invers(xp < 0);
			if (invers)
			{
				xp = -xp;
				yp = -yp;
				zp = -zp;
			}
			
	#ifdef CUDA
			real = no_tex3D(mdlReal, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
			imag = no_tex3D(mdlImag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
	#else
				CpuKernels::complex3D(mdlComplex, real, imag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
	#endif
			
			if(invers)
			    imag = -imag;
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
			ACCCOMPLEX val;
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;

				yp -= mdlInitY;
				zp -= mdlInitZ;

				val =   tex3D<ACCCOMPLEX>(mdlComplex, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
				val.y = -val.y;
			}
			else
			{
				yp -= mdlInitY;
				zp -= mdlInitZ;

				val =   tex3D<ACCCOMPLEX>(mdlComplex, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5, zp + (XFLOAT)0.5);
			}
			real=val.x;
			imag=val.y;
	#endif
#endif
		}
		else
		{
			real = (XFLOAT)0;
			imag = (XFLOAT)0;
		}
	}

#ifdef CUDA
	__device__ __forceinline__ void project2Dmodel(
#else
	void project2Dmodel(
#endif
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

        real=(XFLOAT)0;
		imag=(XFLOAT)0;
		
		r2 = x*x + y*y;
		if (r2 <= maxR2)
		{
			XFLOAT xp = (e0 * x + e1 * y ) * padding_factor;
			XFLOAT yp = (e3 * x + e4 * y ) * padding_factor;
#ifdef CUDA_NO_TEXTURES
			bool invers(xp < 0);
			if (invers)
			{
				xp = -xp;
				yp = -yp;
			}
			
	#ifdef CUDA
			real = no_tex2D(mdlReal, xp, yp, mdlX, mdlInitY);
			imag = no_tex2D(mdlImag, xp, yp, mdlX, mdlInitY);
	#else
			CpuKernels::complex2D(mdlComplex, real, imag, xp, yp, mdlX, mdlInitY);	
	#endif
			
			if(invers)
			    imag = -imag;
			
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
			ACCCOMPLEX val;
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				yp -= mdlInitY;

				val = tex2D<ACCCOMPLEX>(mdlComplex, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5);
				val.y = -val.y;
			}
			else
			{
				yp -= mdlInitY;
				val = tex2D<ACCCOMPLEX>(mdlComplex, xp + (XFLOAT)0.5, yp + (XFLOAT)0.5);
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

	static AccProjectorKernel makeKernel(AccProjector &p, int imgX, int imgY, int imgZ, int imgMaxR)
	{
		int maxR = p.mdlMaxR >= imgMaxR ? imgMaxR : p.mdlMaxR;

		AccProjectorKernel k(
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
#ifdef CUDA
					p.mdlReal,
					p.mdlImag
#else
					p.mdlComplex
#endif
#endif
#endif
				);
		return k;
	}
};  // class AccProjectorKernel


#endif
