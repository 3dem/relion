#ifndef CPU_UTILITIES_H
#define CPU_UTILITIES_H

#include <src/macros.h>
#include <math.h>
#include "src/acc/cpu/cpu_settings.h"
#include <cassert>

namespace CpuKernels
{
	
#define CHECK_INDEX_DEBUG_FATAL( err ) (HandleCheckIndexPtrDebugFatal( err, __FILE__, __LINE__ ))
static void HandleCheckIndexPtrDebugFatal( const char *err, const char *file, int line )
{
    	fprintf(stderr, "DEBUG ERROR: %s in %s:%d\n", err, file, line );
		fflush(stdout);
		raise(SIGSEGV);
}
template <typename T> 
class checkedArray
{
	private:
		T *underlyingData;
      
	public:
		void initCheckedArray(T *dataToCheck)
		{
			underlyingData = dataToCheck;
		}
		
		T& operator[](size_t idx)
		{
			if (idx > std::numeric_limits<int>::max())
				CHECK_INDEX_DEBUG_FATAL("array index > std::numeric_limits<int>::max()");
			return underlyingData[idx];
		}
		const T& operator[](size_t idx) const
  		{
			if (idx > std::numeric_limits<int>::max())
				CHECK_INDEX_DEBUG_FATAL("const: array index > std::numeric_limits<int>::max()");
			return underlyingData[idx];
		}       
};

/*
 * For the following functions always use fast, low-precision intrinsics
 

template< typename T1, typename T2 >
static inline
int floorfracf(T1 a, T2 b)
{
//	return __float2int_rd(__fdividef( (float)a, (float)b ) );
	return (int)(a/b);
}

template< typename T1, typename T2 >
static inline
int ceilfracf(T1 a, T2 b)
{
//	return __float2int_ru(__fdividef( (float)a, (float)b ) );
	return (int)(a/b + 1);
}
 */
static inline
int floorfracf(int a, int b)
{
        return (int)(a/b);
}

static inline
size_t floorfracf(size_t a, int b)
{
        return (size_t)(a/b);
}

static inline
int floorfracf(int a, size_t b)
{
        return (int)(a/b);
}

static inline
size_t floorfracf(size_t a, size_t b)
{
        return (size_t)(a/b);
}

static inline
int ceilfracf(int a, int b)
{
        return (int)(a/b + 1);
}

static inline
size_t ceilfracf(size_t a, int b)
{
        return (size_t)(a/b + (size_t)1);
}

static inline
int ceilfracf(int a, size_t b)
{
        return (int)(a/b + 1);
}

static inline
size_t ceilfracf(size_t a, size_t b)
{
        return (size_t)(a/b + (size_t)1);
}

#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
static inline
XFLOAT no_tex2D(XFLOAT* mdl, XFLOAT xp, XFLOAT yp, int mdlX, int mdlInitY)
{
	int x0 = floorf(xp);
	XFLOAT fx = xp - x0;
	int x1 = x0 + 1;

	int y0 = floorf(yp);
	XFLOAT fy = yp - y0;
	y0 -= mdlInitY;
	int y1 = y0 + 1;

	//-----------------------------
	XFLOAT d00 = mdl[(size_t)y0*(size_t)mdlX+x0];
	XFLOAT d01 = mdl[(size_t)y0*(size_t)mdlX+x1];
	XFLOAT d10 = mdl[(size_t)y1*(size_t)mdlX+x0];
	XFLOAT d11 = mdl[(size_t)y1*(size_t)mdlX+x1];
	//-----------------------------
	XFLOAT dx0 = d00 + (d01 - d00)*fx;
	XFLOAT dx1 = d10 + (d11 - d10)*fx;
	//-----------------------------

	return dx0 + (dx1 - dx0)*fy;
}

// 2D linear interpolation for complex data that interleaves real and
// imaginary data, rather than storing them in a separate array
#ifdef __INTEL_COMPILER
#pragma omp declare simd uniform(mdlX,mdlInitY)
#else
__attribute__((always_inline))
inline
#endif
static void complex2D(std::complex<XFLOAT> *mdlComplex, XFLOAT &real, XFLOAT &imag,
               XFLOAT xp, XFLOAT yp, int mdlX, int mdlInitY)
{
	int x0 = floorf(xp);
	XFLOAT fx = xp - x0;

	int y0 = floorf(yp);
	XFLOAT fy = yp - y0;
	y0 -= mdlInitY;

    size_t offset1 = ((size_t)y0 * (size_t)mdlX + (size_t)x0);
    size_t offset2 = offset1 + (size_t)1;
    size_t offset3 = offset1 + (size_t)mdlX;
    size_t offset4 = offset3 + (size_t)1;  
    
	//-----------------------------
	XFLOAT d00[2], d01[2], d10[2], d11[2];

	d00[0] = mdlComplex[offset1].real();  d00[1] = mdlComplex[offset1].imag();
	d01[0] = mdlComplex[offset2].real();  d01[1] = mdlComplex[offset2].imag();		    
	d10[0] = mdlComplex[offset3].real();  d10[1] = mdlComplex[offset3].imag();
	d11[0] = mdlComplex[offset4].real();  d11[1] = mdlComplex[offset4].imag();
	
	//-----------------------------
	XFLOAT dx0[2], dx1[2];
	
	dx0[0] = d00[0] + (d01[0] - d00[0]) * fx;
	dx1[0] = d10[0] + (d11[0] - d10[0]) * fx;

	dx0[1] = d00[1] + (d01[1] - d00[1]) * fx;
	dx1[1] = d10[1] + (d11[1] - d10[1]) * fx;
	
	//-----------------------------	
	real = dx0[0] + (dx1[0] - dx0[0])*fy;
	imag = dx0[1] + (dx1[1] - dx0[1])*fy;
}

#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
static inline
XFLOAT no_tex3D(
#ifdef DEBUG_CUDA
				XFLOAT* _mdl, 
#else
				XFLOAT* mdl, 
#endif
				XFLOAT xp, XFLOAT yp, XFLOAT zp, 
				int mdlX, int mdlXY, int mdlInitY, int mdlInitZ)
{
#ifdef DEBUG_CUDA
	checkedArray<XFLOAT> mdl;
	mdl.initCheckedArray(_mdl);
#endif
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

	XFLOAT d000 = mdl[(size_t)z0*(size_t)mdlXY+(size_t)y0*(size_t)mdlX+x0];
	XFLOAT d001 = mdl[(size_t)z0*(size_t)mdlXY+(size_t)y0*(size_t)mdlX+x1];
	XFLOAT d010 = mdl[(size_t)z0*(size_t)mdlXY+(size_t)y1*(size_t)mdlX+x0];
	XFLOAT d011 = mdl[(size_t)z0*(size_t)mdlXY+(size_t)y1*(size_t)mdlX+x1];
	XFLOAT d100 = mdl[(size_t)z1*(size_t)mdlXY+(size_t)y0*(size_t)mdlX+x0];
	XFLOAT d101 = mdl[(size_t)z1*(size_t)mdlXY+(size_t)y0*(size_t)mdlX+x1];
	XFLOAT d110 = mdl[(size_t)z1*(size_t)mdlXY+(size_t)y1*(size_t)mdlX+x0];
	XFLOAT d111 = mdl[(size_t)z1*(size_t)mdlXY+(size_t)y1*(size_t)mdlX+x1];
	//-----------------------------
	XFLOAT dx00 = d000 + (d001 - d000)*fx;
	XFLOAT dx01 = d100 + (d101 - d100)*fx;
	XFLOAT dx10 = d010 + (d011 - d010)*fx;
	XFLOAT dx11 = d110 + (d111 - d110)*fx;
	//-----------------------------
	XFLOAT dxy0 = dx00 + (dx10 - dx00)*fy;
	XFLOAT dxy1 = dx01 + (dx11 - dx01)*fy;
	//-----------------------------
	return dxy0 + (dxy1 - dxy0)*fz;
}

// 3D linear interpolation for complex data that interleaves real and
// imaginary data, rather than storing them in a separate array
#ifdef __INTEL_COMPILER
#pragma omp declare simd uniform(mdlX,mdlXY,mdlInitY,mdlInitZ)
#else
__attribute__((always_inline))
inline
#endif
static void complex3D(
				std::complex<XFLOAT> * mdlComplex, 
				XFLOAT &real, XFLOAT &imag,
				XFLOAT xp, XFLOAT yp, XFLOAT zp, int mdlX, int mdlXY, int mdlInitY, int mdlInitZ
		)
{
	int x0 = floorf(xp);
	XFLOAT fx = xp - x0;

	int y0 = floorf(yp);
	XFLOAT fy = yp - y0;
	y0 -= mdlInitY;

	int z0 = floorf(zp);
	XFLOAT fz = zp - z0;
	z0 -= mdlInitZ;

    size_t offset1 = (size_t)((size_t)z0*(size_t)mdlXY+(size_t)y0*(size_t)mdlX+(size_t)x0);
    size_t offset2 = offset1 + (size_t)1;
    size_t offset3 = offset1 + (size_t)mdlX;
    size_t offset4 = offset3 + (size_t)1;
    size_t offset5 = offset1 + (size_t)mdlXY;
    size_t offset6 = offset2 + (size_t)mdlXY;
    size_t offset7 = offset3 + (size_t)mdlXY;
    size_t offset8 = offset4 + (size_t)mdlXY;    
        
	XFLOAT d000[2], d001[2], d010[2], d011[2];
	XFLOAT d100[2], d101[2], d110[2], d111[2]; 

	d000[0] = mdlComplex[offset1].real(); d000[1] = mdlComplex[offset1].imag();
	d001[0] = mdlComplex[offset2].real(); d001[1] = mdlComplex[offset2].imag();    
	d010[0] = mdlComplex[offset3].real(); d010[1] = mdlComplex[offset3].imag();
	d011[0] = mdlComplex[offset4].real(); d011[1] = mdlComplex[offset4].imag();
	d100[0] = mdlComplex[offset5].real(); d100[1] = mdlComplex[offset5].imag();
	d101[0] = mdlComplex[offset6].real(); d101[1] = mdlComplex[offset6].imag();
	d110[0] = mdlComplex[offset7].real(); d110[1] = mdlComplex[offset7].imag();
	d111[0] = mdlComplex[offset8].real(); d111[1] = mdlComplex[offset8].imag();                    
                                
	//-----------------------------
	XFLOAT dx00[2], dx01[2], dx10[2], dx11[2];
	dx00[0] = d000[0] + (d001[0] - d000[0])*fx;
	dx01[0] = d100[0] + (d101[0] - d100[0])*fx;
	dx10[0] = d010[0] + (d011[0] - d010[0])*fx;
	dx11[0] = d110[0] + (d111[0] - d110[0])*fx;

    dx00[1] = d000[1] + (d001[1] - d000[1])*fx;
    dx01[1] = d100[1] + (d101[1] - d100[1])*fx;
    dx10[1] = d010[1] + (d011[1] - d010[1])*fx;
    dx11[1] = d110[1] + (d111[1] - d110[1])*fx;
	
	//-----------------------------
	XFLOAT dxy0[2], dxy1[2];	
	dxy0[0] = dx00[0] + (dx10[0] - dx00[0])*fy;
	dxy1[0] = dx01[0] + (dx11[0] - dx01[0])*fy;
	
    dxy0[1] = dx00[1] + (dx10[1] - dx00[1])*fy;
    dxy1[1] = dx01[1] + (dx11[1] - dx01[1])*fy;

	//-----------------------------
	
	real = dxy0[0] + (dxy1[0] - dxy0[0])*fz;
	imag = dxy0[1] + (dxy1[1] - dxy0[1])*fz;	
}

} // end of namespace CpuKernels

#endif //CPU_UTILITIES_H
