#include "src/acc/cpu/device_stubs.h"

#include "src/acc/acc_ptr.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/cpu/cpu_benchmark_utils.h"
#include "src/acc/cpu/cpu_helper_functions.h"
#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_kernels/diff2.h"
#include "src/acc/cpu/cpu_kernels/wavg.h"
#include "src/acc/cpu/cpu_kernels/BP.h"
#include "src/acc/utilities.h"
#include "src/acc/data_types.h"

#include "src/acc/acc_helper_functions.h"

#include "src/acc/cpu/cpu_kernels/cpu_utils.h"

namespace CpuKernels
{

/*
 * This draft of a kernel assumes input that has jobs which have a single orientation and sequential translations within each job.
 *
 */
void exponentiate_weights_fine(
		XFLOAT *g_pdf_orientation,
		bool *g_pdf_orientation_zeros,
		XFLOAT *g_pdf_offset,
		bool *g_pdf_offset_zeros,
		XFLOAT *g_weights,
		XFLOAT min_diff2,
		unsigned long  oversamples_orient,
		unsigned long  oversamples_trans,
		unsigned long *d_rot_id,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num,
		long int job_num)
{
	for (long int jobid=0; jobid<job_num; jobid++)
	{
		long int pos = d_job_idx[jobid];
		// index of comparison
		long int ix = d_rot_id   [pos];   // each thread gets its own orient...
		long int iy = d_trans_idx[pos];   // ...and it's starting trans...
		long int in = d_job_num  [jobid]; // ...AND the number of translations to go through

		int c_itrans;
		for (long int itrans=0; itrans < in; itrans++, iy++)
		{
			c_itrans = ( iy - (iy % oversamples_trans))/ oversamples_trans;

			if( g_weights[pos+itrans] < min_diff2 || g_pdf_orientation_zeros[ix] || g_pdf_offset_zeros[c_itrans])
				g_weights[pos+itrans] = std::numeric_limits<XFLOAT>::lowest();
			else
				g_weights[pos+itrans] = g_pdf_orientation[ix] + g_pdf_offset[c_itrans] + min_diff2 - g_weights[pos+itrans];
		}
	}
}
void RNDnormalDitributionComplexWithPowerModulation2D(ACCCOMPLEX* Image, size_t xdim, XFLOAT *spectra)
{
	size_t x,y,size;
	size = xdim*((xdim-1)*2);
	for(size_t i=0; i<size; i++)
    {
		y = ( i / xdim ); // fftshift in one of two dims;
		if(y>=xdim)
			y -= (xdim-1)*2;
		x = i % xdim;

		int ires = (int)(sqrtf(x*x + y*y));

		if(ires<xdim)
		{
			Image[i].x = rnd_gaus(0., spectra[ires]);
			Image[i].y = rnd_gaus(0., spectra[ires]);
		}
		else
		{
			Image[i].x = 0;
			Image[i].y = 0;
		}
    }
}

void RNDnormalDitributionComplexWithPowerModulation3D(ACCCOMPLEX* Image, size_t xdim, size_t ydim, XFLOAT *spectra)
{
	int x,y,z,xydim(xdim*ydim),size;
	size = xdim*((xdim-1)*2);				//assuming square input images (particles)
	for(int i=0; i<size; i++)
    {
	   	z = i / xydim;
        y = ( i - (z*xydim) / xdim );
        x = i % xdim;
        // fftshift in two of three dims;
        if(z>=xdim)
     	   z -= (xdim-1)*2;					//assuming square input images (particles)
        if(y>=xdim)
           y -= (xdim-1)*2;					//assuming square input images (particles)

		int ires = (int)(sqrtf(x*x + y*y + z*z));

		if(ires<xdim)
		{
			Image[i].x = rnd_gaus(0., spectra[ires]);
			Image[i].y = rnd_gaus(0., spectra[ires]);
		}
		else
		{
			Image[i].x = 0;
			Image[i].y = 0;
		}
    }
}

void softMaskBackgroundValue(
	int grid_dim, int block_dim,
	XFLOAT *vol, long int vol_size,
	long int xdim,  long int ydim,  long int zdim,
	long int xinit, long int yinit, long int zinit,
	XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width,
	XFLOAT *g_sum, XFLOAT *g_sum_bg
) {
    const size_t stride = block_dim * grid_dim;
    const size_t xydim = xdim * ydim;

    const auto weight = [radius, radius_p, cosine_width] (XFLOAT r) -> XFLOAT {
        return r < radius   ? 0.0 :
               r > radius_p ? 1.0 :
        #if defined(ACC_DOUBLE_PRECISION)
        0.5  + 0.5  * cos (M_PI * (radius_p - r) / cosine_width);
        #else
        0.5f + 0.5f * cosf(M_PI * (radius_p - r) / cosine_width);
        #endif
    };

	for (int bid = 0; bid < grid_dim;  bid++)
	for (int tid = 0; tid < block_dim; tid++) {
		for (size_t i = bid * blockDim.x + tid; i < vol_size; i += stride) {
			const int z = i / xydim        - zinit;
			const int y = i % xydim / xdim - yinit;
			const int x = i % xydim % xdim - xinit;

			const XFLOAT r = sqrt(XFLOAT(x * x + y * y + z * z));
			const XFLOAT w = weight(r);

			g_sum   [tid] += w;
			g_sum_bg[tid] += w * vol[i];
		}
	}
}
void cosineFilter(
	int block_dim, int block_size,
	XFLOAT *vol, long int vol_size,
	long int xdim,  long int ydim,  long int zdim,
	long int xinit, long int yinit, long int zinit,
	bool do_noise, XFLOAT *noise,
	XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width,
	XFLOAT bg_value
) {
	const size_t stride = block_dim * grid_dim;
    const size_t xydim = xdim * ydim;

    const auto weight = [radius, radius_p, cosine_width] (XFLOAT r) -> XFLOAT {
        return r < radius   ? 0.0 :
               r > radius_p ? 1.0 :
        #if defined(ACC_DOUBLE_PRECISION)
        0.5  + 0.5  * cos (M_PI * (radius_p - r) / cosine_width);
        #else
        0.5f + 0.5f * cosf(M_PI * (radius_p - r) / cosine_width);
        #endif
    };

	for (int bid = 0; bid < grid_dim;  bid++)
	for (int tid = 0; tid < block_dim; tid++) {
		for (size_t i = bid * blockDim.x + tid; i < vol_size; i += stride) {
			const int z = i / xydim        - zinit;
			const int y = i % xydim / xdim - yinit;
			const int x = i % xydim % xdim - xinit;

			const XFLOAT r = sqrt(XFLOAT(x * x + y * y + z * z));
			const XFLOAT w = weight(r);
			const XFLOAT bg = do_noise ? noise[i] : bg_value;

			vol[i] = vol[i] * (1 - w) + bg * w;
		}
	}
}

template <typename T>
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
inline
#endif
void cpu_translate2D(T *	g_image_in,
					T *		g_image_out,
					size_t	image_size,
					int		xdim,
					int		ydim,
					int		dx,
					int		dy)
{
	int x,y,xp,yp;
	size_t new_pixel;

#ifdef DEBUG_CUDA
	if (image_size > (size_t)std::numeric_limits<int>::max())
		ACC_PTR_DEBUG_INFO("cpu_translate2D: image_size > std::numeric_limits<int>::max()");
#endif
	for(size_t pixel=0; pixel<image_size; pixel++)
	{
		x = pixel % xdim;
		y = (pixel-x) / (xdim);

		xp = x + dx;
		yp = y + dy;

		if( yp>=0 && xp>=0 && yp<ydim && xp<xdim)
		{
			new_pixel = yp*xdim + xp;
			if(new_pixel>=0 && new_pixel<image_size) // if displacement is negative, new_pixel could be less than 0
				g_image_out[new_pixel] = g_image_in[pixel];
		}
	}
}

template <typename T>
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
inline
#endif
void cpu_translate3D(T *	g_image_in,
					T*		g_image_out,
					size_t	image_size,
					int		xdim,
					int		ydim,
					int		zdim,
					int		dx,
					int		dy,
					int		dz)
{
	int x,y,z,xp,yp,zp,xy;
	size_t new_voxel;

#ifdef DEBUG_CUDA
	if (image_size > (size_t)std::numeric_limits<int>::max())
		ACC_PTR_DEBUG_INFO("cpu_translate3D: image_size > std::numeric_limits<int>::max()");
#endif
	for(size_t voxel=0; voxel<image_size; voxel++)
	{
		int xydim = xdim*ydim;

		z =  voxel / xydim;
		zp = z + dz;

		xy = voxel % xydim;
		y =  xy / xdim;
		yp = y + dy;

		x =  xy % xdim;
		xp = x + dx;

		if( zp>=0 && yp>=0 && xp>=0 && zp<zdim && yp<ydim && xp<xdim)
		{
			new_voxel = zp*xydim +  yp*xdim + xp;
			if(new_voxel>=0 && new_voxel<image_size) // if displacement is negative, new_pixel could be less than 0
				g_image_out[new_voxel] = g_image_in[voxel];
		}
	}
}

/* TODO - if create optimized CPU version of autopicker
 * All these functions need to be converted to use internal loops rather than
 * block and thread indices to operate like other active functions seen in this file
void probRatio( int       blockIdx_x,
				int       threadIdx_x,
				XFLOAT   *d_Mccf,
				XFLOAT   *d_Mpsi,
				XFLOAT   *d_Maux,
				XFLOAT   *d_Mmean,
				XFLOAT   *d_Mstddev,
				size_t       image_size,
				XFLOAT    normfft,
				XFLOAT    sum_ref_under_circ_mask,
				XFLOAT    sum_ref2_under_circ_mask,
				XFLOAT    expected_Pratio,
				int       NpsiThisBatch,
				int       startPsi,
				int       totalPsis)
{
	|* PLAN TO:
	 *
	 * 1) Pre-filter
	 * 		d_Mstddev[i] = 1 / (2*d_Mstddev[i])   ( if d_Mstddev[pixel] > 1E-10 )
	 * 		d_Mstddev[i] = 1    				  ( else )
	 *
	 * 2) Set
	 * 		sum_ref2_under_circ_mask /= 2.
	 *
	 * 3) Total expression becomes
	 * 		diff2 = ( exp(k) - 1.f ) / (expected_Pratio - 1.f)
	 * 	  where
	 * 	  	k = (normfft * d_Maux[pixel] + d_Mmean[pixel] * sum_ref_under_circ_mask)*d_Mstddev[i] + sum_ref2_under_circ_mask
	 *
	 *|

	size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)PROBRATIO_BLOCK_SIZE;
	if(pixel<image_size)
	{
		XFLOAT Kccf = d_Mccf[pixel];
		XFLOAT Kpsi =(XFLOAT)-1.0;
		for(int psi = 0; psi < NpsiThisBatch; psi++ )
		{
			XFLOAT diff2 = normfft * d_Maux[pixel + image_size*psi];
			diff2 += d_Mmean[pixel] * sum_ref_under_circ_mask;

	//		if (d_Mstddev[pixel] > (XFLOAT)1E-10)
			diff2 *= d_Mstddev[pixel];
			diff2 += sum_ref2_under_circ_mask;

#if defined(ACC_DOUBLE_PRECISION)
			diff2 = exp(-diff2 / 2.); // exponentiate to reflect the Gaussian error model. sigma=1 after normalization, 0.4=1/sqrt(2pi)
#else
			diff2 = expf(-diff2 / 2.f);
#endif

			// Store fraction of (1 - probability-ratio) wrt  (1 - expected Pratio)
			diff2 = (diff2 - (XFLOAT)1.0) / (expected_Pratio - (XFLOAT)1.0);
			if (diff2 > Kccf)
			{
				Kccf = diff2;
				Kpsi = (startPsi + psi)*(360/totalPsis);
			}
		}
		d_Mccf[pixel] = Kccf;
		if (Kpsi >= 0.)
			d_Mpsi[pixel] = Kpsi;
	}
}

void rotateOnly(int              blockIdx_x,
				int              blockIdx_y,
				int              threadIdx_x,
				ACCCOMPLEX     *d_Faux,
				XFLOAT           psi,
				AccProjectorKernel &projector,
				int              startPsi
		       )
{
	int proj = blockIdx_y;
	size_t image_size=(size_t)projector.imgX*(size_t)projector.imgY;
	size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	if(pixel<image_size)
	{
		int y = floorfracf(pixel,projector.imgX);
		int x = pixel % projector.imgX;

		if (y > projector.maxR)
		{
			if (y >= projector.imgY - projector.maxR)
				y = y - projector.imgY;
			else
				x = projector.maxR;
		}

		XFLOAT sa, ca;
#if defined(ACC_DOUBLE_PRECISION)
		sincos((proj+startPsi)*psi, &sa, &ca);
#else
		sincosf((proj+startPsi)*psi, &sa, &ca);
#endif

		ACCCOMPLEX val;

		projector.project2Dmodel(	 x,y,
									 ca,
									-sa,
									 sa,
									 ca,
									 val.x,val.y);

		long int out_pixel = proj*image_size + pixel;

		d_Faux[out_pixel].x =val.x;
		d_Faux[out_pixel].y =val.y;
	}
}

void rotateAndCtf(  int              blockIdx_x,
					int              blockIdx_y,
					int              threadIdx_x,
					ACCCOMPLEX     *d_Faux,
					XFLOAT          *d_ctf,
					XFLOAT           psi,
					AccProjectorKernel &projector,
					int       startPsi
				)
{
	int proj = blockIdx_y;
	size_t image_size=(size_t)projector.imgX*(size_t)projector.imgY;
	size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	if(pixel<image_size)
	{
		int y = floorfracf(pixel,projector.imgX);
		int x = pixel % projector.imgX;

		if (y > projector.maxR)
		{
			if (y >= projector.imgY - projector.maxR)
				y = y - projector.imgY;
			else
				x = projector.maxR;
		}

		XFLOAT sa, ca;
#if defined(ACC_DOUBLE_PRECISION)
		sincos((proj+startPsi)*psi, &sa, &ca);
#else
		sincosf((proj+startPsi)*psi, &sa, &ca);
#endif
		ACCCOMPLEX val;

		projector.project2Dmodel(	 x,y,
									 ca,
									-sa,
									 sa,
									 ca,
									 val.x,val.y);

		long int out_pixel = proj*image_size + pixel;

		d_Faux[out_pixel].x =val.x*d_ctf[pixel];
		d_Faux[out_pixel].y =val.y*d_ctf[pixel];

	}
}


void convol_A(  int           blockIdx_x,
				int           threadIdx_x,
				ACCCOMPLEX  *d_A,
				ACCCOMPLEX  *d_B,
				size_t           image_size)
{
	size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	if(pixel<image_size)
	{
		XFLOAT tr =   d_A[pixel].x;
		XFLOAT ti = - d_A[pixel].y;
		d_A[pixel].x =   tr*d_B[pixel].x - ti*d_B[pixel].y;
		d_A[pixel].y =   ti*d_B[pixel].x + tr*d_B[pixel].y;
	}
}

void convol_A(  int          blockIdx_x,
				int          threadIdx_x,
				ACCCOMPLEX *d_A,
				ACCCOMPLEX *d_B,
				ACCCOMPLEX *d_C,
				size_t          image_size)
{
	size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	if(pixel<image_size)
	{
		XFLOAT tr =   d_A[pixel].x;
		XFLOAT ti = - d_A[pixel].y;
		d_C[pixel].x =   tr*d_B[pixel].x - ti*d_B[pixel].y;
		d_C[pixel].y =   ti*d_B[pixel].x + tr*d_B[pixel].y;
	}
}

void batch_convol_A(int           blockIdx_x,
					int           blockIdx_y,
					int           threadIdx_x,
					ACCCOMPLEX  *d_A,
					ACCCOMPLEX  *d_B,
					size_t           image_size)
{
	size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	int A_off = blockIdx_y * image_size;
	if(pixel<image_size)
	{
		XFLOAT tr =   d_A[pixel + A_off].x;
		XFLOAT ti = - d_A[pixel + A_off].y;
		d_A[pixel + A_off].x =   tr*d_B[pixel].x - ti*d_B[pixel].y;
		d_A[pixel + A_off].y =   ti*d_B[pixel].x + tr*d_B[pixel].y;
	}
}

void batch_convol_A(int           blockIdx_x,
					int           blockIdx_y,
					int           threadIdx_x,
					ACCCOMPLEX  *d_A,
					ACCCOMPLEX  *d_B,
					ACCCOMPLEX  *d_C,
					size_t           image_size)
{
	size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	int A_off = blockIdx_y*image_size;
	if(pixel<image_size)
	{
		XFLOAT tr =   d_A[pixel + A_off].x;
		XFLOAT ti = - d_A[pixel + A_off].y;
		d_C[pixel + A_off].x =   tr*d_B[pixel].x - ti*d_B[pixel].y;
		d_C[pixel + A_off].y =   ti*d_B[pixel].x + tr*d_B[pixel].y;
	}
}

void convol_B(  int          blockIdx_x,
				int          threadIdx_x,
				ACCCOMPLEX *d_A,
				ACCCOMPLEX *d_B,
				size_t          image_size)
{
	size_t pixel =  (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	if(pixel<image_size)
	{
		XFLOAT tr = d_A[pixel].x;
		XFLOAT ti = d_A[pixel].y;
		d_A[pixel].x =   tr*d_B[pixel].x + ti*d_B[pixel].y;
		d_A[pixel].y =   ti*d_B[pixel].x - tr*d_B[pixel].y;
	}
}

void convol_B(  int           blockIdx_x,
				int           threadIdx_x,
				ACCCOMPLEX  *d_A,
				ACCCOMPLEX  *d_B,
				ACCCOMPLEX  *d_C,
				size_t           image_size)
{
	size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	if(pixel<image_size)
	{
		XFLOAT tr = d_A[pixel].x;
		XFLOAT ti = d_A[pixel].y;
		d_C[pixel].x =   tr*d_B[pixel].x + ti*d_B[pixel].y;
		d_C[pixel].y =   ti*d_B[pixel].x - tr*d_B[pixel].y;
	}
}

void batch_convol_B(int          blockIdx_x,
					int          blockIdx_y,
					int          threadIdx_x,
					ACCCOMPLEX *d_A,
					ACCCOMPLEX *d_B,
					size_t          image_size)
{
	long int pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	int A_off = blockIdx_y*image_size;
	if(pixel<image_size)
	{
		XFLOAT tr = d_A[pixel + A_off].x;
		XFLOAT ti = d_A[pixel + A_off].y;
		d_A[pixel + A_off].x =   tr*d_B[pixel].x + ti*d_B[pixel].y;
		d_A[pixel + A_off].y =   ti*d_B[pixel].x - tr*d_B[pixel].y;
	}
}
*/
template <typename T>
void cpu_kernel_multi( T *A,
			T *OUT,
			T  S,
			size_t     image_size)
{
#ifdef DEBUG_CUDA
	if (image_size < 0)
		ACC_PTR_DEBUG_INFO("cpu_kernel_multi:  image_size < 0");
#endif
	for (size_t i = 0; i < image_size; i ++)
		OUT[i] = A[i]*S;
}

template <typename T>
void cpu_kernel_multi( T *A,
			T  S,
			size_t     image_size)
{
#ifdef DEBUG_CUDA
	if (image_size < 0)
		ACC_PTR_DEBUG_INFO("cpu_kernel_multi2:  image_size < 0");
#endif
	for (size_t i = 0; i < image_size; i ++)
		A[i] *= S;
}

template <typename T>
void cpu_kernel_multi( T *A,
			T *B,
			T *OUT,
			T  S,
			size_t     image_size)
{
#ifdef DEBUG_CUDA
	if (image_size < 0)
		ACC_PTR_DEBUG_INFO("cpu_kernel_multi3:  image_size < 0");
#endif
	for (size_t i = 0; i < image_size; i ++)
		OUT[i] = A[i]*B[i]*S;
}

template <typename T>
void cpu_kernel_add(
	T *A,
	T  S,
	size_t size
)
{
#ifdef DEBUG_CUDA
	if (size < 0)
		ACC_PTR_DEBUG_INFO("cpu_kernel_add:  image_size < 0");
#endif
	for (size_t i = 0; i < size; i ++)
		A[i] += S;
}

/*
void batch_multi(   int     blockIdx_x,
					int     blockIdx_y,
					int     threadIdx_x,
					XFLOAT *A,
					XFLOAT *B,
					XFLOAT *OUT,
					XFLOAT  S,
					size_t     image_size)
{
	sise_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	if(pixel<image_size)
		OUT[pixel + blockIdx_y*image_size] = A[pixel + blockIdx_y*image_size]*B[pixel + blockIdx_y*image_size]*S;
}
 */
/* TODO - CPU-optimized autopicker
void finalizeMstddev(   int     blockIdx_x,
						int     threadIdx_x,
						XFLOAT *Mstddev,
						XFLOAT *aux,
						XFLOAT  S,
						size_t     image_size)
{
	int size_t = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	if(pixel<image_size)
	{
		XFLOAT temp = Mstddev[pixel] + S * aux[pixel];
		if(temp > 0)
			Mstddev[pixel] = sqrt(temp);
		else
			Mstddev[pixel] = 0;
	}
}

void square(int     blockIdx_x,
			int     threadIdx_x,
			XFLOAT *A,
			size_t     image_size)
{
	size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	if(pixel<image_size)
		A[pixel] = A[pixel]*A[pixel];
}
*/
template<bool invert>
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
inline
#endif
void cpu_kernel_make_eulers_2D(int grid_size, int block_size,
		XFLOAT *alphas,
		XFLOAT *eulers,
		unsigned long orientation_num)
{
#ifdef DEBUG_CUDA
	if ((size_t)grid_size*(size_t)block_size > (size_t)std::numeric_limits<int>::max())
		ACC_PTR_DEBUG_INFO("cpu_kernel_make_eulers_2D: grid_size*block_size > std::numeric_limits<int>::max()");
#endif
	for(int blockIdx_x=0; blockIdx_x<(int)(grid_size); blockIdx_x++) {
		for(int threadIdx_x=0; threadIdx_x<block_size; threadIdx_x++)  {
			unsigned long oid = (unsigned long)blockIdx_x * (unsigned long)block_size + threadIdx_x; //Orientation id

			if (oid >= orientation_num)
				return;

			XFLOAT ca, sa;
			XFLOAT a = alphas[oid] * (XFLOAT)PI / (XFLOAT)180.0;

#ifdef ACC_DOUBLE_PRECISION
			sincos(a, &sa, &ca);
#else
			sincosf(a, &sa, &ca);
#endif

			if(!invert)
			{
				eulers[9 * oid + 0] = ca;//00
				eulers[9 * oid + 1] = sa;//01
				eulers[9 * oid + 2] = 0 ;//02
				eulers[9 * oid + 3] =-sa;//10
				eulers[9 * oid + 4] = ca;//11
				eulers[9 * oid + 5] = 0 ;//12
				eulers[9 * oid + 6] = 0 ;//20
				eulers[9 * oid + 7] = 0 ;//21
				eulers[9 * oid + 8] = 1 ;//22
			}
			else
			{
				eulers[9 * oid + 0] = ca;//00
				eulers[9 * oid + 1] =-sa;//10
				eulers[9 * oid + 2] = 0 ;//20
				eulers[9 * oid + 3] = sa;//01
				eulers[9 * oid + 4] = ca;//11
				eulers[9 * oid + 5] = 0 ;//21
				eulers[9 * oid + 6] = 0 ;//02
				eulers[9 * oid + 7] = 0 ;//12
				eulers[9 * oid + 8] = 1 ;//22
			}
		} // threadIdx_x
	} // blockIdx_x
}

template<bool invert,bool doL, bool doR>
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
inline
#endif
void cpu_kernel_make_eulers_3D(int grid_size, int block_size,
		XFLOAT *alphas,
		XFLOAT *betas,
		XFLOAT *gammas,
		XFLOAT *eulers,
		unsigned long orientation_num,
		XFLOAT *L,
		XFLOAT *R)
{
#ifdef DEBUG_CUDA
	if ((size_t)grid_size*(size_t)block_size > (size_t)std::numeric_limits<int>::max())
		ACC_PTR_DEBUG_INFO("cpu_kernel_make_eulers_3D: grid_size*block_size > std::numeric_limits<int>::max()");
#endif
	for(int blockIdx_x=0; blockIdx_x<(int)(grid_size); blockIdx_x++) {
        for(int threadIdx_x=0; threadIdx_x<block_size; threadIdx_x++) {
			XFLOAT a(0.f),b(0.f),g(0.f), A[9],B[9];
			XFLOAT ca, sa, cb, sb, cg, sg, cc, cs, sc, ss;

			unsigned long oid = (unsigned long)blockIdx_x * (unsigned long)block_size + threadIdx_x; //Orientation id

			if (oid >= orientation_num)
				return;

			for (int i = 0; i < 9; i ++)
				B[i] = (XFLOAT) 0.f;

			a = alphas[oid] * (XFLOAT)PI / (XFLOAT)180.0;
			b = betas[oid]  * (XFLOAT)PI / (XFLOAT)180.0;
			g = gammas[oid] * (XFLOAT)PI / (XFLOAT)180.0;

#ifdef ACC_DOUBLE_PRECISION
			sincos(a, &sa, &ca);
			sincos(b,  &sb, &cb);
			sincos(g, &sg, &cg);
#else
			sincosf(a, &sa, &ca);
			sincosf(b,  &sb, &cb);
			sincosf(g, &sg, &cg);
#endif

			cc = cb * ca;
			cs = cb * sa;
			sc = sb * ca;
			ss = sb * sa;

			A[0] = ( cg * cc - sg * sa);//00
			A[1] = ( cg * cs + sg * ca);//01
			A[2] = (-cg * sb )         ;//02
			A[3] = (-sg * cc - cg * sa);//10
			A[4] = (-sg * cs + cg * ca);//11
			A[5] = ( sg * sb )         ;//12
			A[6] = ( sc )              ;//20
			A[7] = ( ss )              ;//21
			A[8] = ( cb )              ;//22

			if (doR)
			{
				for (int i = 0; i < 9; i++)
					B[i] = 0.f;

				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						for (int k = 0; k < 3; k++)
							B[i * 3 + j] += A[i * 3 + k] * R[k * 3 + j];
			}
			else
				for (int i = 0; i < 9; i++)
					B[i] = A[i];

			if (doL)
			{
				if (doR)
					for (int i = 0; i < 9; i++)
						A[i] = B[i];

				for (int i = 0; i < 9; i++)
					B[i] = 0.f;

				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						for (int k = 0; k < 3; k++)
							B[i * 3 + j] += L[i * 3 + k] * A[k * 3 + j];
			}

			if(invert)
			{
				if (doL) // this could have anisotropy, so inverse neq transpose!!!
				{
					XFLOAT det;
					det =     B[0] * (B[4] * B[8] - B[7] * B[5])
							- B[1] * (B[3] * B[8] - B[6] * B[5])
							+ B[2] * (B[3] * B[7] - B[6] * B[4]);

					eulers[9 * oid + 0] = (B[4] * B[8] - B[7] * B[5]) / det;
					eulers[9 * oid + 1] = (B[7] * B[2] - B[1] * B[8]) / det;
					eulers[9 * oid + 2] = (B[1] * B[5] - B[4] * B[2]) / det;
					eulers[9 * oid + 3] = (B[5] * B[6] - B[8] * B[3]) / det;
					eulers[9 * oid + 4] = (B[8] * B[0] - B[2] * B[6]) / det;
					eulers[9 * oid + 5] = (B[2] * B[3] - B[5] * B[0]) / det;
					eulers[9 * oid + 6] = (B[3] * B[7] - B[6] * B[4]) / det;
					eulers[9 * oid + 7] = (B[6] * B[1] - B[0] * B[7]) / det;
					eulers[9 * oid + 8] = (B[0] * B[4] - B[3] * B[1]) / det;
				}
				else
				{
					eulers[9 * oid + 0] = B[0];//00
					eulers[9 * oid + 1] = B[3];//01
					eulers[9 * oid + 2] = B[6];//02
					eulers[9 * oid + 3] = B[1];//10
					eulers[9 * oid + 4] = B[4];//11
					eulers[9 * oid + 5] = B[7];//12
					eulers[9 * oid + 6] = B[2];//20
					eulers[9 * oid + 7] = B[5];//21
					eulers[9 * oid + 8] = B[8];//22
				}
			}
			else
			{
				eulers[9 * oid + 0] = B[0];//00
				eulers[9 * oid + 1] = B[1];//10
				eulers[9 * oid + 2] = B[2];//20
				eulers[9 * oid + 3] = B[3];//01
				eulers[9 * oid + 4] = B[4];//11
				eulers[9 * oid + 5] = B[5];//21
				eulers[9 * oid + 6] = B[6];//02
				eulers[9 * oid + 7] = B[7];//12
				eulers[9 * oid + 8] = B[8];//22
			}
		} // threadIdx_x
	} // blockIdx_x
}

} // end of namespace CpuKernels


// -------------------------------  Some explicit template instantiations
template void CpuKernels::cpu_translate2D<XFLOAT>(XFLOAT *,
    XFLOAT*, size_t, int, int, int, int);

template void CpuKernels::cpu_translate3D<XFLOAT>(XFLOAT *,
    XFLOAT *, size_t, int, int, int, int, int, int);

template void CpuKernels::cpu_kernel_multi<XFLOAT>( XFLOAT *, XFLOAT, size_t);
template void CpuKernels::cpu_kernel_add<XFLOAT>( XFLOAT *, XFLOAT, size_t);

template void CpuKernels::cpu_kernel_make_eulers_3D<true, true, true>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void CpuKernels::cpu_kernel_make_eulers_3D<true, true, false>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void CpuKernels::cpu_kernel_make_eulers_3D<true, false,true>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void CpuKernels::cpu_kernel_make_eulers_3D<true, false,false>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void CpuKernels::cpu_kernel_make_eulers_3D<false,true, true>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void CpuKernels::cpu_kernel_make_eulers_3D<false,true, false>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void CpuKernels::cpu_kernel_make_eulers_3D<false,false,true>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void CpuKernels::cpu_kernel_make_eulers_3D<false,false,false>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);

template void CpuKernels::cpu_kernel_make_eulers_2D<true>(int, int,
		XFLOAT *, XFLOAT *, unsigned long);
template void CpuKernels::cpu_kernel_make_eulers_2D<false>(int, int,
		XFLOAT *, XFLOAT *, unsigned long);
// ----------------------------------------------------------------------

