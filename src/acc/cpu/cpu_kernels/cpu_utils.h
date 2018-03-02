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
 */

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
	XFLOAT d00 = mdl[y0*mdlX+x0];
	XFLOAT d01 = mdl[y0*mdlX+x1];
	XFLOAT d10 = mdl[y1*mdlX+x0];
	XFLOAT d11 = mdl[y1*mdlX+x1];
	//-----------------------------
	XFLOAT dx0 = d00 + (d01 - d00)*fx;
	XFLOAT dx1 = d10 + (d11 - d10)*fx;
	//-----------------------------

	return dx0 + (dx1 - dx0)*fy;
}

// 2D linear interpolation for complex data that interleaves real and
// imaginary data, rather than storing them in a separate array
static inline
void complex2D(std::complex<XFLOAT> *mdlComplex, XFLOAT &real, XFLOAT &imag,
               XFLOAT xp, XFLOAT yp, int mdlX, int mdlInitY
#if 0
,
		std::complex<XFLOAT> &d00, 
		std::complex<XFLOAT> &d01, 
		std::complex<XFLOAT> &d10, 
		std::complex<XFLOAT> &d11,
		std::complex<XFLOAT> &dx0, 
		std::complex<XFLOAT> &dx1,
		std::complex<XFLOAT> &result
#endif
)
{
	int x0 = floorf(xp);
	XFLOAT fx = xp - x0;

	int y0 = floorf(yp);
	XFLOAT fy = yp - y0;
	y0 -= mdlInitY;

    int offset1 = (y0 * mdlX + x0) * 2;
    int offset2 = offset1 + 2;
    int offset3 = offset1 + mdlX * 2;
    int offset4 = offset3 + 2;  
    
	//-----------------------------
//	std::complex<XFLOAT> d00, d01, d10, d11;
	XFLOAT d00[2], d01[2], d10[2], d11[2];

	d00[0] = mdlComplex[offset1].real();  d00[1] = mdlComplex[offset1].imag();
	d01[0] = mdlComplex[offset2].real();  d01[1] = mdlComplex[offset2].imag();		    
	d10[0] = mdlComplex[offset3].real();  d10[1] = mdlComplex[offset3].imag();
	d11[0] = mdlComplex[offset4].real();  d11[1] = mdlComplex[offset4].imag();
	
	//-----------------------------
//	std::complex<XFLOAT> dx0, dx1;
	XFLOAT dx0[2], dx1[2];
	
	dx0[0] = d00[0] + (d01[0] - d00[0]) * fx;
	dx1[0] = d10[0] + (d11[0] - d10[0]) * fx;

	dx0[1] = d00[1] + (d01[1] - d00[1]) * fx;
	dx1[1] = d10[1] + (d11[1] - d10[1]) * fx;
	
	//-----------------------------
//	std::complex<XFLOAT> result;
	
//	result = dx0 + (dx1 - dx0) * fy;
	
//	real = result.real(); //dx0[0] + (dx1[0] - dx0[0])*fy;
//	imag = result.imag(); //dx0[1] + (dx1[1] - dx0[1])*fy;	
	real = dx0[0] + (dx1[0] - dx0[0])*fy;
	imag = dx0[1] + (dx1[1] - dx0[1])*fy;
}

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

	XFLOAT d000 = mdl[z0*mdlXY+y0*mdlX+x0];
	XFLOAT d001 = mdl[z0*mdlXY+y0*mdlX+x1];
	XFLOAT d010 = mdl[z0*mdlXY+y1*mdlX+x0];
	XFLOAT d011 = mdl[z0*mdlXY+y1*mdlX+x1];
	XFLOAT d100 = mdl[z1*mdlXY+y0*mdlX+x0];
	XFLOAT d101 = mdl[z1*mdlXY+y0*mdlX+x1];
	XFLOAT d110 = mdl[z1*mdlXY+y1*mdlX+x0];
	XFLOAT d111 = mdl[z1*mdlXY+y1*mdlX+x1];
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
__attribute__((vector(uniform(mdlX,mdlXY,mdlInitY,mdlInitZ))))
static inline
void complex3D(
#ifdef DEBUG_CUDA
				std::complex<XFLOAT> * _mdlComplex, 
#else
				std::complex<XFLOAT> * mdlComplex, 
#endif
				XFLOAT &real, XFLOAT &imag,
				XFLOAT xp, XFLOAT yp, XFLOAT zp, int mdlX, int mdlXY, int mdlInitY, int mdlInitZ
#if 0
,
		std::complex<XFLOAT> &d000, 
		std::complex<XFLOAT> &d001, 
		std::complex<XFLOAT> &d010, 
		std::complex<XFLOAT> &d011,
		std::complex<XFLOAT> &d100, 
		std::complex<XFLOAT> &d101, 
		std::complex<XFLOAT> &d110, 
		std::complex<XFLOAT> &d111,
		std::complex<XFLOAT> &dx00, 
		std::complex<XFLOAT> &dx01, 
		std::complex<XFLOAT> &dx10, 
		std::complex<XFLOAT> &dx11,
		std::complex<XFLOAT> &dxy0, 
		std::complex<XFLOAT> &dxy1,
		std::complex<XFLOAT> &result
#endif
		)
{
#ifdef DEBUG_CUDA
	checkedArray<std::complex<XFLOAT> > mdlComplex;
	mdlComplex.initCheckedArray(_mdlComplex);
#endif
	int x0 = floorf(xp);
	XFLOAT fx = xp - x0;

	int y0 = floorf(yp);
	XFLOAT fy = yp - y0;
	y0 -= mdlInitY;

	int z0 = floorf(zp);
	XFLOAT fz = zp - z0;
	z0 -= mdlInitZ;

    int offset1 = (z0*mdlXY+y0*mdlX+x0) * 2;
    int offset2 = offset1 + 2;
    int offset3 = offset1 + mdlX * 2;
    int offset4 = offset3 + 2;
    int offset5 = offset1 + mdlXY * 2;
    int offset6 = offset2 + mdlXY * 2;
    int offset7 = offset3 + mdlXY * 2;
    int offset8 = offset4 + mdlXY * 2;    
        
//    std::complex<XFLOAT> d000, d001, d010, d011; //XFLOAT d000[2], d001[2], d010[2], d011[2];
//    std::complex<XFLOAT> d100, d101, d110, d111;//XFLOAT d100[2], d101[2], d110[2], d111[2];   
	XFLOAT d000[2], d001[2], d010[2], d011[2];
	XFLOAT d100[2], d101[2], d110[2], d111[2]; 
//    std::complex<XFLOAT> d100, d101, d110, d111;//XFLOAT d100[2], d101[2], d110[2], d111[2];  
  
	
 //   d000 = mdlComplex[offset1]; //
	d000[0] = mdlComplex[offset1].real(); d000[1] = mdlComplex[offset1].imag();
 //   d001 = mdlComplex[offset2]; // 
	d001[0] = mdlComplex[offset2].real(); d001[1] = mdlComplex[offset2].imag();    
 //   d010 = mdlComplex[offset3]; // 
	d010[0] = mdlComplex[offset3].real(); d010[1] = mdlComplex[offset3].imag();
 //   d011 = mdlComplex[offset4]; // 
	d011[0] = mdlComplex[offset4].real(); d011[1] = mdlComplex[offset4].imag();
 //   d100 = mdlComplex[offset5]; // 
	d100[0] = mdlComplex[offset5].real(); d100[1] = mdlComplex[offset5].imag();
  //  d101 = mdlComplex[offset6]; // 
	d101[0] = mdlComplex[offset6].real(); d101[1] = mdlComplex[offset6].imag();
 //   d110 = mdlComplex[offset7]; // 
	d110[0] = mdlComplex[offset7].real(); d110[1] = mdlComplex[offset7].imag();
 //   d111 = mdlComplex[offset8]; // 
	d111[0] = mdlComplex[offset8].real(); d111[1] = mdlComplex[offset8].imag();                    
                                
	//-----------------------------
//	 std::complex<XFLOAT> dx00, dx01, dx10, dx11; //
	XFLOAT dx00[2], dx01[2], dx10[2], dx11[2];
//    dx00 = d000 + (d001 - d000)*fx; //
	dx00[0] = d000[0] + (d001[0] - d000[0])*fx;
//    dx01 = d100 + (d101 - d100)*fx; //
	dx01[0] = d100[0] + (d101[0] - d100[0])*fx;
//    dx10 = d010 + (d011 - d010)*fx; //
	dx10[0] = d010[0] + (d011[0] - d010[0])*fx;
//    dx11 = d110 + (d111 - d110)*fx; //
	dx11[0] = d110[0] + (d111[0] - d110[0])*fx;

    dx00[1] = d000[1] + (d001[1] - d000[1])*fx;
    dx01[1] = d100[1] + (d101[1] - d100[1])*fx;
    dx10[1] = d010[1] + (d011[1] - d010[1])*fx;
    dx11[1] = d110[1] + (d111[1] - d110[1])*fx;
	
	//-----------------------------
//   std::complex<XFLOAT> dxy0, dxy1; //
	XFLOAT dxy0[2], dxy1[2];	
//	dxy0 = dx00 + (dx10 - dx00)*fy;  //
	dxy0[0] = dx00[0] + (dx10[0] - dx00[0])*fy;
//	dxy1 = dx01 + (dx11 - dx01)*fy;  //
	dxy1[0] = dx01[0] + (dx11[0] - dx01[0])*fy;
	
    dxy0[1] = dx00[1] + (dx10[1] - dx00[1])*fy;
    dxy1[1] = dx01[1] + (dx11[1] - dx01[1])*fy;

	//-----------------------------
//	std::complex<XFLOAT> result;
	
//	result = dxy0 + (dxy1 - dxy0)*fz;
//	real = result.real(); //dxy0[0] + (dxy1[0] - dxy0[0])*fz;
//	imag = result.imag(); //dxy0[1] + (dxy1[1] - dxy0[1])*fz;	
	real = dxy0[0] + (dxy1[0] - dxy0[0])*fz;
	imag = dxy0[1] + (dxy1[1] - dxy0[1])*fz;	
}

} // end of namespace CpuKernels

#endif //CPU_UTILITIES_H
