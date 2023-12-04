#include "src/acc/sycl/device_stubs.h"

#include "src/acc/acc_ptr.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/sycl/sycl_benchmark_utils.h"
#include "src/acc/sycl/sycl_helper_functions.h"
#include "src/acc/sycl/sycl_kernels/helper.h"
#include "src/acc/utilities.h"
#include "src/acc/data_types.h"

#include "src/acc/acc_helper_functions.h"

#include "src/acc/sycl/sycl_kernels/sycl_utils.h"

namespace syclKernels
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
				g_weights[pos+itrans] = std::numeric_limits<XFLOAT>::lowest(); //large negative number
			else
				g_weights[pos+itrans] = g_pdf_orientation[ix] + g_pdf_offset[c_itrans] + min_diff2 - g_weights[pos+itrans];
		}
	}
}

void RNDnormalDitributionComplexWithPowerModulation2D(ACCCOMPLEX *Image, size_t xdim, XFLOAT *spectra)
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

void RNDnormalDitributionComplexWithPowerModulation3D(ACCCOMPLEX *Image, size_t xdim, size_t ydim, XFLOAT *spectra)
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

void softMaskBackgroundValue(	int      block_dim,
                                int      block_size,
                                XFLOAT  *vol,
								long int vol_size,
								long int xdim,
								long int ydim,
								long int zdim,
								long int xinit,
								long int yinit,
								long int zinit,
								XFLOAT   radius,
								XFLOAT   radius_p,
								XFLOAT   cosine_width,
								XFLOAT  *g_sum,
								XFLOAT  *g_sum_bg)
{
	for(int bid=0; bid<block_dim; bid++)
	{
		for(int tid=0; tid<block_size; tid++)
		{
	//		vol.setXmippOrigin(); // sets xinit=xdim , also for y z
			XFLOAT r, raisedcos;
			int x,y,z;
			XFLOAT img_pixels;

			size_t texel_pass_num = ceilfracf((size_t)vol_size,(size_t)block_size*(size_t)block_dim);
			size_t texel = (size_t)bid*(size_t)block_size*(size_t)texel_pass_num + tid;

			for (size_t pass = 0; pass < texel_pass_num; pass++, texel+=block_size) // loop through all translations for this orientation
			{
				if(texel<vol_size)
				{
					img_pixels = vol[texel];

					z =   texel / (xdim*ydim) ;
					y = ( texel % (xdim*ydim) ) / xdim ;
					x = ( texel % (xdim*ydim) ) % xdim ;

					z-=zinit;
					y-=yinit;
					x-=xinit;

					r = sqrt(XFLOAT(x*x + y*y + z*z));

					if (r < radius)
						continue;
					else if (r > radius_p)
					{
						g_sum[tid]    += (XFLOAT)1.0;
						g_sum_bg[tid] += img_pixels;
					}
					else
					{
	#if defined(ACC_DOUBLE_PRECISION)
						raisedcos = 0.5 + 0.5  * cos ( (radius_p - r) / cosine_width * M_PI);
	#else
						raisedcos = 0.5 + 0.5  * cosf( (radius_p - r) / cosine_width * M_PI);
	#endif
						g_sum[tid] += raisedcos;
						g_sum_bg[tid] += raisedcos * img_pixels;
					}
				}
			}
		} // tid
	} // bid
}

void cosineFilter(	int      block_dim,
					int      block_size,
					XFLOAT  *vol,
					long int vol_size,
					long int xdim,
					long int ydim,
					long int zdim,
					long int xinit,
					long int yinit,
					long int zinit,
					bool     do_noise,
					XFLOAT *noise,
					XFLOAT   radius,
					XFLOAT   radius_p,
					XFLOAT   cosine_width,
					XFLOAT   bg_value)
{
	for(int bid=0; bid<block_dim; bid++)
	{
		for(int tid=0; tid<block_size; tid++)
		{
//		vol.setXmippOrigin(); // sets xinit=xdim , also for y z
			XFLOAT r, raisedcos, defVal;
			int x,y,z;
			XFLOAT img_pixels;

			size_t texel_pass_num = ceilfracf((size_t)vol_size,(size_t)block_size*(size_t)block_dim);
			size_t texel = (size_t)bid*(size_t)block_size*(size_t)texel_pass_num + tid;

			defVal = bg_value;
			for (size_t pass = 0; pass < texel_pass_num; pass++, texel+=block_size) // loop the available warps enough to complete all translations for this orientation
			{
				if(texel<vol_size)
				{
					img_pixels= vol[texel];

					z =   texel / (xdim*ydim) ;
					y = ( texel % (xdim*ydim) ) / xdim ;
					x = ( texel % (xdim*ydim) ) % xdim ;

					z-=zinit;
					y-=yinit;
					x-=xinit;

					r = sqrt(XFLOAT(x*x + y*y + z*z));

					if(do_noise)
                                		defVal = noise[texel];

					if (r < radius)
						continue;
					else if (r > radius_p)
						img_pixels=defVal;
					else
					{
		#if defined(ACC_DOUBLE_PRECISION)
						raisedcos = 0.5 + 0.5  * cos ( (radius_p - r) / cosine_width * M_PI);
		#else
						raisedcos = 0.5 + 0.5  * cosf( (radius_p - r) / cosine_width * M_PI);
		#endif
						img_pixels= img_pixels*(1-raisedcos) + defVal*raisedcos;

					}
					vol[texel]=img_pixels;
				}
			}
		} // tid
	} // bid
}

template <typename T>
__attribute__((always_inline))
inline
void sycl_translate2D(T		*g_image_in,
					T		*g_image_out,
					size_t	image_size,
					int		xdim,
					int		ydim,
					int		dx,
					int		dy)
{
	int x,y,xp,yp;
	size_t new_pixel;

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
__attribute__((always_inline))
inline
void sycl_translate3D(T		*g_image_in,
					T		*g_image_out,
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

template <typename T>
void sycl_kernel_multi( T *A,
			T *OUT,
			T  S,
			size_t     image_size)
{
	for (size_t i = 0; i < image_size; i ++)
		OUT[i] = A[i]*S;
}

template <typename T>
void sycl_kernel_multi( T *A,
			T  S,
			size_t     image_size)
{
	for (size_t i = 0; i < image_size; i ++)
		A[i] *= S;
}

template <typename T>
void sycl_kernel_multi( T *A,
			T *B,
			T *OUT,
			T  S,
			size_t     image_size)
{
	for (size_t i = 0; i < image_size; i ++)
		OUT[i] = A[i]*B[i]*S;
}

template <typename T>
void sycl_kernel_add(
	T *A,
	T  S,
	size_t size
)
{
	for (size_t i = 0; i < size; i ++)
		A[i] += S;
}

template<bool invert>
__attribute__((always_inline))
inline
void sycl_kernel_make_eulers_2D(int grid_size, int block_size,
		XFLOAT *alphas,
		XFLOAT *eulers,
		unsigned long orientation_num)
{
	for(int blockIdx_x=0; blockIdx_x<(int)(grid_size); blockIdx_x++) {
		for(int threadIdx_x=0; threadIdx_x<block_size; threadIdx_x++) {
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
__attribute__((always_inline))
inline
void sycl_kernel_make_eulers_3D(int grid_size, int block_size,
		XFLOAT *alphas,
		XFLOAT *betas,
		XFLOAT *gammas,
		XFLOAT *eulers,
		unsigned long orientation_num,
		XFLOAT *L,
		XFLOAT *R)
{
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

} // end of namespace syclKernels


// -------------------------------  Some explicit template instantiations
template void syclKernels::sycl_translate2D<XFLOAT>(XFLOAT *,
	XFLOAT *, size_t, int, int, int, int);

template void syclKernels::sycl_translate3D<XFLOAT>(XFLOAT *,
	XFLOAT *, size_t, int, int, int, int, int, int);

template void syclKernels::sycl_kernel_multi<XFLOAT>( XFLOAT *, XFLOAT, size_t);
template void syclKernels::sycl_kernel_add<XFLOAT>( XFLOAT *, XFLOAT, size_t);

template void syclKernels::sycl_kernel_make_eulers_3D<true, true, true>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void syclKernels::sycl_kernel_make_eulers_3D<true, true, false>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void syclKernels::sycl_kernel_make_eulers_3D<true, false,true>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void syclKernels::sycl_kernel_make_eulers_3D<true, false,false>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void syclKernels::sycl_kernel_make_eulers_3D<false,true, true>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void syclKernels::sycl_kernel_make_eulers_3D<false,true, false>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void syclKernels::sycl_kernel_make_eulers_3D<false,false,true>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);
template void syclKernels::sycl_kernel_make_eulers_3D<false,false,false>(int, int,
		XFLOAT *, XFLOAT *, XFLOAT *, XFLOAT *, unsigned long, XFLOAT *, XFLOAT *);

template void syclKernels::sycl_kernel_make_eulers_2D<true>(int, int,
		XFLOAT *, XFLOAT *, unsigned long);
template void syclKernels::sycl_kernel_make_eulers_2D<false>(int, int,
		XFLOAT *, XFLOAT *, unsigned long);
// ----------------------------------------------------------------------

