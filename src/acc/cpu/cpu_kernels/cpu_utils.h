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

// 2D linear interpolation for complex data that interleaves real and
// imaginary data, rather than storing them in a separate array
#if _OPENMP >= 201307	// For OpenMP 4.0 and later
#pragma omp declare simd uniform(mdlX,mdlInitY)
#pragma omp declare simd uniform(mdlX,mdlInitY) simdlen(4)
#pragma omp declare simd uniform(mdlX,mdlInitY) simdlen(8)
#pragma omp declare simd uniform(mdlX,mdlInitY) simdlen(16)
#else
inline
#endif
static void complex2D(std::complex<XFLOAT> *mdlComplex, XFLOAT &real, XFLOAT &imag,
				const XFLOAT xp, const XFLOAT yp, const int mdlX, const int mdlInitY)
{
	const int x0 = floorf(xp);
	const XFLOAT fx = xp - x0;

	const int y0 = floorf(yp);
	const XFLOAT fy = yp - y0;

	const size_t offset1 = ((size_t)(y0-mdlInitY) * (size_t)mdlX + (size_t)x0);
	const size_t offset2 = offset1 + (size_t)1;
	const size_t offset3 = offset1 + (size_t)mdlX;
	const size_t offset4 = offset3 + (size_t)1;

	const XFLOAT d00[2] = { mdlComplex[offset1].real(), mdlComplex[offset1].imag() };
	const XFLOAT d01[2] = { mdlComplex[offset2].real(), mdlComplex[offset2].imag() };
	const XFLOAT d10[2] = { mdlComplex[offset3].real(), mdlComplex[offset3].imag() };
	const XFLOAT d11[2] = { mdlComplex[offset4].real(), mdlComplex[offset4].imag() };

	const XFLOAT dx0[2] = { d00[0] + (d01[0] - d00[0]) * fx, d00[1] + (d01[1] - d00[1]) * fx };
	const XFLOAT dx1[2] = { d10[0] + (d11[0] - d10[0]) * fx, d10[1] + (d11[1] - d10[1]) * fx };

	real = dx0[0] + (dx1[0] - dx0[0])*fy;
	imag = dx0[1] + (dx1[1] - dx0[1])*fy;
}

// 3D linear interpolation for complex data that interleaves real and
// imaginary data, rather than storing them in a separate array
#if _OPENMP >= 201307	// For OpenMP 4.0 and later
#pragma omp declare simd uniform(mdlX,mdlXY,mdlInitY,mdlInitZ)
#pragma omp declare simd uniform(mdlX,mdlXY,mdlInitY,mdlInitZ) simdlen(4)
#pragma omp declare simd uniform(mdlX,mdlXY,mdlInitY,mdlInitZ) simdlen(8)
#pragma omp declare simd uniform(mdlX,mdlXY,mdlInitY,mdlInitZ) simdlen(16)
#else
inline
#endif
static void complex3D(
				std::complex<XFLOAT> * mdlComplex, 
				XFLOAT &real, XFLOAT &imag,
				const XFLOAT xp, const XFLOAT yp, const XFLOAT zp,
				const int mdlX, const int mdlXY, const int mdlInitY, const int mdlInitZ
		)
{
	const int x0 = floorf(xp);
	const XFLOAT fx = xp - x0;

	const int y0 = floorf(yp);
	const XFLOAT fy = yp - y0;

	const int z0 = floorf(zp);
	const XFLOAT fz = zp - z0;

	const size_t offset1 = (size_t)((size_t)(z0-mdlInitZ)*(size_t)mdlXY+(size_t)(y0-mdlInitY)*(size_t)mdlX+(size_t)x0);
	const size_t offset2 = offset1 + (size_t)1;
	const size_t offset3 = offset1 + (size_t)mdlX;
	const size_t offset4 = offset3 + (size_t)1;
	const size_t offset5 = offset1 + (size_t)mdlXY;
	const size_t offset6 = offset2 + (size_t)mdlXY;
	const size_t offset7 = offset3 + (size_t)mdlXY;
	const size_t offset8 = offset4 + (size_t)mdlXY;

	const XFLOAT d000[2] = { mdlComplex[offset1].real(), mdlComplex[offset1].imag() };
	const XFLOAT d001[2] = { mdlComplex[offset2].real(), mdlComplex[offset2].imag() };
	const XFLOAT d010[2] = { mdlComplex[offset3].real(), mdlComplex[offset3].imag() };
	const XFLOAT d011[2] = { mdlComplex[offset4].real(), mdlComplex[offset4].imag() };
	const XFLOAT d100[2] = { mdlComplex[offset5].real(), mdlComplex[offset5].imag() };
	const XFLOAT d101[2] = { mdlComplex[offset6].real(), mdlComplex[offset6].imag() };
	const XFLOAT d110[2] = { mdlComplex[offset7].real(), mdlComplex[offset7].imag() };
	const XFLOAT d111[2] = { mdlComplex[offset8].real(), mdlComplex[offset8].imag() };

	const XFLOAT dx00[2] = { d000[0] + (d001[0] - d000[0])*fx, d000[1] + (d001[1] - d000[1])*fx };
	const XFLOAT dx10[2] = { d010[0] + (d011[0] - d010[0])*fx, d010[1] + (d011[1] - d010[1])*fx };
	const XFLOAT dx01[2] = { d100[0] + (d101[0] - d100[0])*fx, d100[1] + (d101[1] - d100[1])*fx };
	const XFLOAT dx11[2] = { d110[0] + (d111[0] - d110[0])*fx, d110[1] + (d111[1] - d110[1])*fx };

	const XFLOAT dxy0[2] = {  dx00[0] + (dx10[0] - dx00[0])*fy,
								dx00[1] + (dx10[1] - dx00[1])*fy };
	const XFLOAT dxy1[2] = {  dx01[0] + (dx11[0] - dx01[0])*fy,
								dx01[1] + (dx11[1] - dx01[1])*fy };

	real = dxy0[0] + (dxy1[0] - dxy0[0])*fz;
	imag = dxy0[1] + (dxy1[1] - dxy0[1])*fz;
}

} // end of namespace CpuKernels

#endif //CPU_UTILITIES_H
