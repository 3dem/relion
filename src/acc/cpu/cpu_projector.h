#ifndef CPU_PROJECTOR_H_
#define CPU_PROJECTOR_H_

#include "src/complex.h"
#include "src/acc/cpu/cpu_settings.h"

class CpuProjector
{
	friend class CpuProjectorKernel;

	int mdlX, mdlY, mdlZ, mdlXYZ, mdlMaxR,
	    mdlInitY, mdlInitZ,
	    padding_factor;

	size_t allocaton_size;


	XFLOAT *mdlComplex;

public:
	CpuProjector():
			mdlX(0), mdlY(0), mdlZ(0),
			mdlXYZ(0), mdlMaxR(0),
			mdlInitY(0), mdlInitZ(0),
			padding_factor(0),
			allocaton_size(0)
	{
		mdlComplex = 0;
	}

	bool setMdlDim(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int maxr, int paddingFactor);

	void initMdl(XFLOAT *real, XFLOAT *imag);
	void initMdl(Complex *data);

	void clear();

	~CpuProjector()
	{
		clear();
	};

};

static inline 
void complex2D(XFLOAT* mdlComplex, XFLOAT &real, XFLOAT &imag,
               XFLOAT xp, XFLOAT yp, int mdlX, int mdlInitY)
{
#ifdef ACC_DOUBLE_PRECISION
	int x0 = floor(xp);
#else
	int x0 = floorf(xp);
#endif
	
	XFLOAT fx = xp - x0;

#ifdef ACC_DOUBLE_PRECISION
	int y0 = floor(yp);
#else	
    int y0 = floorf(yp);
#endif
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
void complex3D(XFLOAT* mdlComplex, XFLOAT &real, XFLOAT &imag,
               XFLOAT xp, XFLOAT yp, XFLOAT zp, int mdlX, int mdlXY, int mdlInitY, int mdlInitZ)
{
#ifdef ACC_DOUBLE_PRECISION
	int x0 = floor(xp);
#else	
    int x0 = floorf(xp);
#endif
    
	XFLOAT fx = xp - x0;

#ifdef ACC_DOUBLE_PRECISION
	int y0 = floor(yp);
#else	
    int y0 = floorf(yp);
#endif
    
	XFLOAT fy = yp - y0;
	y0 -= mdlInitY;

#ifdef ACC_DOUBLE_PRECISION
	int z0 = floor(zp);
#else
    int z0 = floorf(zp);
#endif	
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


class CpuProjectorKernel
{

public:
	int mdlX, mdlXY, mdlZ,
		imgX, imgY, imgZ, 
		mdlInitY, mdlInitZ,
		padding_factor,
		maxR, maxR2;

	XFLOAT *mdlComplex;

	CpuProjectorKernel(
			int mdlX, int mdlY, int mdlZ,
			int imgX, int imgY, int imgZ,
			int mdlInitY, int mdlInitZ,
			int padding_factor,
			int maxR,
			XFLOAT *mdlComplex
			):
			mdlX(mdlX), mdlXY(mdlX*mdlY), mdlZ(mdlZ),
			imgX(imgX), imgY(imgY), imgZ(imgZ),
			mdlInitY(mdlInitY), mdlInitZ(mdlInitZ),
			padding_factor(padding_factor),
			maxR(maxR), maxR2(maxR*maxR),
			mdlComplex(mdlComplex)
		{};

	CpuProjectorKernel(
			int mdlX, int mdlY, int mdlZ,
			int imgX, int imgY, int imgZ,
			int mdlInitY, int mdlInitZ,
			int padding_factor,
			int maxR,
			XFLOAT *mdlReal, XFLOAT *mdlImag
			):
				mdlX(mdlX), mdlXY(mdlX*mdlY), mdlZ(mdlZ),
				imgX(imgX), imgY(imgY), imgZ(imgZ),
				mdlInitY(mdlInitY), mdlInitZ(mdlInitZ),
				padding_factor(padding_factor),
				maxR(maxR), maxR2(maxR*maxR)
			{
			    for(int i=0; i<mdlX * mdlY * mdlZ; i++) {
			        *mdlComplex ++ = *mdlReal ++;
			        *mdlComplex ++ = *mdlImag ++;			        
			    }
			};

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
        real=(XFLOAT)0;
		imag=(XFLOAT)0;
		
		int r2 = x*x + y*y + z*z;
		if (r2 <= maxR2)
		{
			XFLOAT xp = (e0 * x + e1 * y + e2 * z ) * padding_factor;
			XFLOAT yp = (e3 * x + e4 * y + e5 * z ) * padding_factor;
			XFLOAT zp = (e6 * x + e7 * y + e8 * z ) * padding_factor;

            int neg = 0;
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;
				
				neg = 1;
	        }
	        
			complex3D(mdlComplex, real, imag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);			
			if(neg)
			    imag = -imag;
		}
	}
	
//	 #pragma forceinline
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
        real=(XFLOAT)0;
		imag=(XFLOAT)0;
		
		int r2 = x*x + y*y;
		if (r2 <= maxR2)
		{
			XFLOAT xp = (e0 * x + e1 * y ) * padding_factor;
  		    XFLOAT yp = (e3 * x + e4 * y ) * padding_factor;
			XFLOAT zp = (e6 * x + e7 * y ) * padding_factor;
			
			int neg = 0;
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;
				
				neg = 1;
			}

			complex3D(mdlComplex, real, imag, xp, yp, zp, mdlX, mdlXY, mdlInitY, mdlInitZ);
			if(neg)
			    imag = -imag;
		}
	}

//	#pragma forceinline
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
        real=(XFLOAT)0;
		imag=(XFLOAT)0;

        int r2 = x*x + y*y;
	    if (r2 <= maxR2)
	    {
	        XFLOAT xp = (e0 * x + e1 * y ) * padding_factor;
	        XFLOAT yp = (e3 * x + e4 * y ) * padding_factor;
	        int neg = 0;
	        if (xp < 0)
	        {
	            // Get complex conjugated hermitian symmetry pair
		        xp = -xp;
                yp = -yp;
       
                neg = 1;
            }
               
			complex2D(mdlComplex, real, imag, xp, yp, mdlX, mdlInitY);		
			if(neg)
			    imag = -imag;
	    }
	}

	static CpuProjectorKernel makeKernel(CpuProjector &p, int imgX, int imgY, int imgZ, int imgMaxR)
	{
		int maxR = p.mdlMaxR >= imgMaxR ? imgMaxR : p.mdlMaxR;

		CpuProjectorKernel k(
					p.mdlX, p.mdlY, p.mdlZ,
					imgX, imgY, imgZ, 
					p.mdlInitY, p.mdlInitZ,
					p.padding_factor,
					maxR,
					p.mdlComplex
				);
		return k;
	}
};


#endif
