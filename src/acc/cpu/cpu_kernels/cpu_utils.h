#ifndef CPU_UTILITIES_H
#define CPU_UTILITIES_H

#include <src/macros.h>
#include <math.h>
#include "src/acc/cpu/cpu_settings.h"

namespace CpuKernels
{

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
void complex2D(XFLOAT* mdlComplex, XFLOAT &real, XFLOAT &imag,
               XFLOAT xp, XFLOAT yp, int mdlX, int mdlInitY)
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
	XFLOAT d00[2], d01[2], d10[2], d11[2];
	d00[0] = mdlComplex[offset1], d00[1] = mdlComplex[offset1+1];
	d01[0] = mdlComplex[offset2], d01[1] = mdlComplex[offset2+1];    
	d10[0] = mdlComplex[offset3], d10[1] = mdlComplex[offset3+1];    
	d11[0] = mdlComplex[offset4], d11[1] = mdlComplex[offset4+1];    			    
		
	//-----------------------------
	XFLOAT dx0[2], dx1[2];

    dx0[0] = d00[0] + (d01[0] - d00[0])*fx;
    dx1[0] = d10[0] + (d11[0] - d10[0])*fx;

    dx0[1] = d00[1] + (d01[1] - d00[1])*fx;
    dx1[1] = d10[1] + (d11[1] - d10[1])*fx;

	//-----------------------------

	real = dx0[0] + (dx1[0] - dx0[0])*fy;
	imag = dx0[1] + (dx1[1] - dx0[1])*fy;	
}

static inline
XFLOAT no_tex3D(XFLOAT* mdl, XFLOAT xp, XFLOAT yp, XFLOAT zp, int mdlX, int mdlXY, int mdlInitY, int mdlInitZ)
{
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
static inline
void complex3D(XFLOAT* mdlComplex, XFLOAT &real, XFLOAT &imag,
               XFLOAT xp, XFLOAT yp, XFLOAT zp, int mdlX, int mdlXY, int mdlInitY, int mdlInitZ)
{
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
        
    XFLOAT d000[2], d001[2], d010[2], d011[2];
    XFLOAT d100[2], d101[2], d110[2], d111[2];   
    
    d000[0] = mdlComplex[offset1], d000[1] = mdlComplex[offset1+1];
    d001[0] = mdlComplex[offset2], d001[1] = mdlComplex[offset2+1];    
    d010[0] = mdlComplex[offset3], d010[1] = mdlComplex[offset3+1];
    d011[0] = mdlComplex[offset4], d011[1] = mdlComplex[offset4+1];
    d100[0] = mdlComplex[offset5], d100[1] = mdlComplex[offset5+1];
    d101[0] = mdlComplex[offset6], d101[1] = mdlComplex[offset6+1];
    d110[0] = mdlComplex[offset7], d110[1] = mdlComplex[offset7+1];
    d111[0] = mdlComplex[offset8], d111[1] = mdlComplex[offset8+1];                    
                                
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

template<bool check_max_r2>
void window_fourier_transform(
		XFLOAT *g_in_real,
		XFLOAT *g_in_imag,
		XFLOAT *g_out_real,
		XFLOAT *g_out_imag,
		unsigned iX, unsigned iY, unsigned iZ, unsigned iYX, //Input dimensions
		unsigned oX, unsigned oY, unsigned oZ, unsigned oYX, //Output dimensions
		unsigned max_idx,
		int grid_dim,
		int Npsi,
		int block_size,
		unsigned max_r2 = 0
		)
{
	for(int grd=0; grd<grid_dim; grd++) {
		for(int psi=0; psi<Npsi; psi++) {
			for(int tid=0; tid< block_size; tid++) {
				unsigned n = tid + block_size * grd;
				long int image_offset = oX*oY*oZ*psi;
				if (n >= max_idx) continue;

				int k, i, kp, ip, jp;

				if (check_max_r2)
				{
					k = n / (iX * iY);
					i = (n % (iX * iY)) / iX;

					kp = k < iX ? k : k - iZ;
					ip = i < iX ? i : i - iY;
					jp = n % iX;

					if (kp*kp + ip*ip + jp*jp > max_r2)
						continue;
				}
				else
				{
					k = n / (oX * oY);
					i = (n % (oX * oY)) / oX;

					kp = k < oX ? k : k - oZ;
					ip = i < oX ? i : i - oY;
					jp = n % oX;
				}

				g_out_real[(kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp + image_offset] = g_in_real[(kp < 0 ? kp + iZ : kp)*iYX + (ip < 0 ? ip + iY : ip)*iX + jp + image_offset];
				g_out_imag[(kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp + image_offset] = g_in_imag[(kp < 0 ? kp + iZ : kp)*iYX + (ip < 0 ? ip + iY : ip)*iX + jp + image_offset];
			} // tid
		} // blk
	} // grd
}

#define WINDOW_FT_BLOCK_SIZE 128
template<bool check_max_r2>
void window_fourier_transform(
		ACCCOMPLEX *g_in,
		ACCCOMPLEX *g_out,
		size_t iX, size_t iY, size_t iZ, size_t iYX, //Input dimensions
		size_t oX, size_t oY, size_t oZ, size_t oYX, //Output dimensions
		size_t max_idx,
		int grid_dim,
		int Npsi,
		int block_size,
		size_t max_r2 = 0
		)
{
	for(int grd=0; grd<grid_dim; grd++) {
		for(int psi=0; psi<Npsi; psi++) {
			for(int tid=0; tid< block_size; tid++) {
				size_t n = tid + block_size * grd;
				size_t oOFF = oX*oY*oZ*psi;
				size_t iOFF = iX*iY*iZ*psi;
				if (n >= max_idx) continue;

				long int k, i, kp, ip, jp;

				if (check_max_r2)
				{
					k = n / (iX * iY);
					i = (n % (iX * iY)) / iX;

					kp = k < iX ? k : k - iZ;
					ip = i < iX ? i : i - iY;
					jp = n % iX;

					if (kp*kp + ip*ip + jp*jp > max_r2)
						continue;
				}
				else
				{
					k = n / (oX * oY);
					i = (n % (oX * oY)) / oX;

					kp = k < oX ? k : k - oZ;
					ip = i < oX ? i : i - oY;
					jp = n % oX;
				}

				long int  in_idx = (kp < 0 ? kp + iZ : kp) * iYX + (ip < 0 ? ip + iY : ip)*iX + jp;
				long int out_idx = (kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp;
				g_out[out_idx + oOFF] =  g_in[in_idx + iOFF];
			} // tid
		} // blk
	} // grd
}

} // end of namespace CpuKernels

#endif //CPU_UTILITIES_H
