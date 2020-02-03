#include "src/acc/settings.h"
#include "src/acc/cuda/cuda_kernels/cuda_device_utils.cuh"
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#include "src/acc/cuda/cuda_settings.h"

#include <curand.h>
#include <curand_kernel.h>

/// Needed explicit template instantiations
template __global__ void cuda_kernel_make_eulers_2D<true>(XFLOAT *,
	XFLOAT *, unsigned);
template __global__ void cuda_kernel_make_eulers_2D<false>(XFLOAT *,
	XFLOAT *, unsigned);

template __global__ void cuda_kernel_make_eulers_3D<true, true, true>(XFLOAT *,
		XFLOAT *, XFLOAT *, XFLOAT *, unsigned, XFLOAT *, XFLOAT *);
template __global__ void cuda_kernel_make_eulers_3D<true, true, false>(XFLOAT *,
		XFLOAT *, XFLOAT *, XFLOAT *, unsigned, XFLOAT *, XFLOAT *);
template __global__ void cuda_kernel_make_eulers_3D<true, false,true>(XFLOAT *,
		XFLOAT *, XFLOAT *, XFLOAT *, unsigned, XFLOAT *, XFLOAT *);
template __global__ void cuda_kernel_make_eulers_3D<true, false,false>(XFLOAT *,
		XFLOAT *, XFLOAT *, XFLOAT *, unsigned, XFLOAT *, XFLOAT *);
template __global__ void cuda_kernel_make_eulers_3D<false,true, true>(XFLOAT *,
		XFLOAT *, XFLOAT *, XFLOAT *, unsigned, XFLOAT *, XFLOAT *);
template __global__ void cuda_kernel_make_eulers_3D<false,true, false>(XFLOAT *,
		XFLOAT *, XFLOAT *, XFLOAT *, unsigned, XFLOAT *, XFLOAT *);
template __global__ void cuda_kernel_make_eulers_3D<false,false,true>(XFLOAT *,
		XFLOAT *, XFLOAT *, XFLOAT *, unsigned, XFLOAT *, XFLOAT *);
template __global__ void cuda_kernel_make_eulers_3D<false,false,false>(XFLOAT *,
		XFLOAT *, XFLOAT *, XFLOAT *, unsigned, XFLOAT *, XFLOAT *);

/*
 * This draft of a kernel assumes input that has jobs which have a single orientation and sequential translations within each job.
 *
 */
__global__ void cuda_kernel_exponentiate_weights_fine(
		XFLOAT *g_pdf_orientation,
		bool *g_pdf_orientation_zeros,
		XFLOAT *g_pdf_offset,
		bool *g_pdf_offset_zeros,
		XFLOAT *g_weights,
		XFLOAT min_diff2,
		int oversamples_orient,
		int oversamples_trans,
		unsigned long *d_rot_id,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num,
		long int job_num)
{
	// blockid
	int bid  = blockIdx.x;
	//threadid
	int tid = threadIdx.x;

	long int jobid = bid*SUMW_BLOCK_SIZE+tid;

	if (jobid<job_num)
	{
		long int pos = d_job_idx[jobid];
		// index of comparison
		long int ix = d_rot_id   [pos];   // each thread gets its own orient...
		long int iy = d_trans_idx[pos];   // ...and it's starting trans...
		long int in = d_job_num  [jobid]; // ...AND the number of translations to go through

		int c_itrans;
		for (int itrans=0; itrans < in; itrans++, iy++)
		{
			c_itrans = ( iy - (iy % oversamples_trans))/ oversamples_trans;

			if( g_weights[pos+itrans] < min_diff2 || g_pdf_orientation_zeros[ix] || g_pdf_offset_zeros[c_itrans])
				g_weights[pos+itrans] = -99e99; //large negative number
			else
				g_weights[pos+itrans] = g_pdf_orientation[ix] + g_pdf_offset[c_itrans] + min_diff2 - g_weights[pos+itrans];
		}
	}
}

__global__ void cuda_kernel_initRND(unsigned long seed, curandState *States)
{
       int tid = threadIdx.x;
       int bid = blockIdx.x;

       int id    = bid*RND_BLOCK_SIZE + tid;
       int pixel = bid*RND_BLOCK_SIZE + tid;

       curand_init(seed, pixel, 0, &States[id]);
}

__global__ void cuda_kernel_RNDnormalDitributionComplexWithPowerModulation2D( ACCCOMPLEX *Image,
																		    curandState *States,
																		    long int xdim,
																			XFLOAT * spectra)
{
       int tid = threadIdx.x;
       int bid = blockIdx.x;

       int id    = bid*RND_BLOCK_SIZE + tid;
       int pixel = bid*RND_BLOCK_SIZE + tid;

       //curand_init(1234, pixel, 0, &States[id]);

       int x,y;
       int size = xdim*((xdim-1)*2);   					//assuming square input images (particles)
       int passes = size/(RND_BLOCK_NUM*RND_BLOCK_SIZE) + 1;
       for(int i=0; i<passes; i++)
       {
               if(pixel<size)
               {
                       y = ( pixel / xdim );
                       x = pixel % xdim;

                       // fftshift in one of two dims;
                       if(y>=xdim)
                               y -= (xdim-1)*2;   		//assuming square input images (particles)

                       int ires = rintf(sqrtf(x*x + y*y));
#if defined(ACC_DOUBLE_PRECISION)
                       XFLOAT scale = 0.;
                       if(ires<xdim)
                               scale =  spectra[ires];

                       Image[pixel] = (curand_normal2_double(&States[id]))*scale;
#else
                       XFLOAT scale = 0.f;
                       if(ires<xdim)
                               scale =  spectra[ires];

                       Image[pixel] = (curand_normal2(&States[id]))*scale;
#endif
               }
               pixel += RND_BLOCK_NUM*RND_BLOCK_SIZE;
       }
}
__global__ void cuda_kernel_RNDnormalDitributionComplexWithPowerModulation3D( ACCCOMPLEX *Image,
																		    curandState *States,
																		    long int xdim,
                                                                            long int ydim,
																			XFLOAT * spectra)
{
       int tid = threadIdx.x;
       int bid = blockIdx.x;

       int id    = bid*RND_BLOCK_SIZE + tid;
       int pixel = bid*RND_BLOCK_SIZE + tid;

       //curand_init(1234, pixel, 0, &States[id]);

       int x,y,z,xydim(xdim*ydim);
       int size = xdim*((xdim-1)*2)*((xdim-1)*2);   		//assuming square input images (particles)
       int passes = size/(RND_BLOCK_NUM*RND_BLOCK_SIZE) + 1;
       for(int i=0; i<passes; i++)
       {
               if(pixel<size)
               {
            	   	   z = pixel / xydim;
                       y = ( pixel - (z*xydim) / xdim );
                       x = pixel % xdim;
                       // fftshift in two of three dims;
                       if(z>=xdim)
                    	   z -= (xdim-1)*2;					//assuming square input images (particles)
                       if(y>=xdim)
                           y -= (xdim-1)*2;					//assuming square input images (particles)


                       int ires = rintf(sqrtf(x*x + y*y + z*z));
#if defined(ACC_DOUBLE_PRECISION)
                       XFLOAT scale = 0.;
                       if(ires<xdim)
                               scale =  spectra[ires];

                       Image[pixel] = (curand_normal2_double(&States[id]))*scale;
#else
                       XFLOAT scale = 0.f;
                       if(ires<xdim)
                               scale =  spectra[ires];

                       Image[pixel] = (curand_normal2(&States[id]))*scale;
#endif
               }
               pixel += RND_BLOCK_NUM*RND_BLOCK_SIZE;
       }
}


//__global__ void cuda_kernel_exponentiate_weights_fine2(
//		XFLOAT *g_pdf_orientation,
//		XFLOAT *g_pdf_offset,
//		XFLOAT *g_weights,
//		XFLOAT avg_diff2,
//		int oversamples_orient,
//		int oversamples_trans,
//		unsigned long *d_rot_id,
//		unsigned long *d_trans_idx,
//		unsigned long *d_job_idx,
//		unsigned long *d_job_num,
//		long int job_num)
//{
//	// blockid
//	int bid  = blockIdx.x;
//	//threadid
//	int tid = threadIdx.x;
//
//	long int jobid = bid*SUMW_BLOCK_SIZE+tid;
//
//	if (jobid<job_num)
//	{
//		long int pos = d_job_idx[jobid];
//		// index of comparison
//		long int iy = d_trans_idx[ pos];
//		long int in =  d_job_num[jobid];
//
//		int c_itrans;
//		for (int itrans=0; itrans < in; itrans++, iy++)
//		{
//			XFLOAT a = g_weights[pos+itrans] + avg_diff2;
//
//#if defined(ACC_DOUBLE_PRECISION)
//			if (a < -700.)
//				g_weights[pos+itrans] = 0.;
//			else
//				g_weights[pos+itrans] = exp(a);
//#else
//			if (a < -88.)
//				g_weights[pos+itrans] = 0.f;
//			else
//				g_weights[pos+itrans] = expf(a);
//#endif
//		}
//	}
//}

__global__ void cuda_kernel_softMaskOutsideMap(	XFLOAT *vol,
												long int vol_size,
												long int xdim,
												long int ydim,
												long int zdim,
												long int xinit,
												long int yinit,
												long int zinit,
												bool do_Mnoise,
												XFLOAT radius,
												XFLOAT radius_p,
												XFLOAT cosine_width	)
{

		int tid = threadIdx.x;

//		vol.setXmippOrigin(); // sets xinit=xdim , also for y z
		XFLOAT r, raisedcos;

		__shared__ XFLOAT     img_pixels[SOFTMASK_BLOCK_SIZE];
		__shared__ XFLOAT    partial_sum[SOFTMASK_BLOCK_SIZE];
		__shared__ XFLOAT partial_sum_bg[SOFTMASK_BLOCK_SIZE];

		XFLOAT sum_bg_total =  (XFLOAT)0.0;

		long int texel_pass_num = ceilfracf(vol_size,SOFTMASK_BLOCK_SIZE);
		int texel = tid;

		partial_sum[tid]=(XFLOAT)0.0;
		partial_sum_bg[tid]=(XFLOAT)0.0;
		if (do_Mnoise)
		{
			for (int pass = 0; pass < texel_pass_num; pass++, texel+=SOFTMASK_BLOCK_SIZE) // loop the available warps enough to complete all translations for this orientation
			{
				XFLOAT x,y,z;
				if(texel<vol_size)
				{
					img_pixels[tid]=__ldg(&vol[texel]);

					z = floor( (float) texel                   / (float)((xdim)*(ydim)));
					y = floor( (XFLOAT)(texel-z*(xdim)*(ydim)) / (XFLOAT) xdim );
					x = texel - z*(xdim)*(ydim) - y*xdim;

					z-=zinit;
					y-=yinit;
					x-=xinit;

					r = sqrt(x*x + y*y + z*z);

					if (r < radius)
						continue;
					else if (r > radius_p)
					{
						partial_sum[tid]    += (XFLOAT)1.0;
						partial_sum_bg[tid] += img_pixels[tid];
					}
					else
					{
#if defined(ACC_DOUBLE_PRECISION)
						raisedcos = 0.5 + 0.5  * cospi( (radius_p - r) / cosine_width );
#else
						raisedcos = 0.5f + 0.5f * cospif((radius_p - r) / cosine_width );
#endif
						partial_sum[tid] += raisedcos;
						partial_sum_bg[tid] += raisedcos * img_pixels[tid];
					}
				}
			}
		}

		__syncthreads();
		for(int j=(SOFTMASK_BLOCK_SIZE/2); j>0; j/=2)
		{
			if(tid<j)
			{
				partial_sum[tid] += partial_sum[tid+j];
				partial_sum_bg[tid] += partial_sum_bg[tid+j];
			}
			__syncthreads();
		}

		sum_bg_total  = partial_sum_bg[0] / partial_sum[0];


		texel = tid;
		for (int pass = 0; pass < texel_pass_num; pass++, texel+=SOFTMASK_BLOCK_SIZE) // loop the available warps enough to complete all translations for this orientation
		{
			XFLOAT x,y,z;
			if(texel<vol_size)
			{
				img_pixels[tid]=__ldg(&vol[texel]);

				z =  floor( (float) texel                  / (float)((xdim)*(ydim)));
				y = floor( (XFLOAT)(texel-z*(xdim)*(ydim)) / (XFLOAT)  xdim         );
				x = texel - z*(xdim)*(ydim) - y*xdim;

				z-=zinit;
				y-=yinit;
				x-=xinit;

				r = sqrt(x*x + y*y + z*z);

				if (r < radius)
					continue;
				else if (r > radius_p)
					img_pixels[tid]=sum_bg_total;
				else
				{
#if defined(ACC_DOUBLE_PRECISION)
					raisedcos = 0.5  + 0.5  * cospi( (radius_p - r) / cosine_width );
#else
					raisedcos = 0.5f + 0.5f * cospif((radius_p - r) / cosine_width );
#endif
					img_pixels[tid]= img_pixels[tid]*(1-raisedcos) + sum_bg_total*raisedcos;

				}
				vol[texel]=img_pixels[tid];
			}

		}
}

__global__ void cuda_kernel_softMaskBackgroundValue(	XFLOAT *vol,
														long int vol_size,
														long int xdim,
														long int ydim,
														long int zdim,
														long int xinit,
														long int yinit,
														long int zinit,
														XFLOAT radius,
														XFLOAT radius_p,
														XFLOAT cosine_width,
														XFLOAT *g_sum,
														XFLOAT *g_sum_bg)
{

		int tid = threadIdx.x;
		int bid = blockIdx.x;

//		vol.setXmippOrigin(); // sets xinit=xdim , also for y z
		XFLOAT r, raisedcos;
		int x,y,z;
		__shared__ XFLOAT     img_pixels[SOFTMASK_BLOCK_SIZE];
		__shared__ XFLOAT    partial_sum[SOFTMASK_BLOCK_SIZE];
		__shared__ XFLOAT partial_sum_bg[SOFTMASK_BLOCK_SIZE];

		long int texel_pass_num = ceilfracf(vol_size,SOFTMASK_BLOCK_SIZE*gridDim.x);
		int texel = bid*SOFTMASK_BLOCK_SIZE*texel_pass_num + tid;

		partial_sum[tid]=(XFLOAT)0.0;
		partial_sum_bg[tid]=(XFLOAT)0.0;

		for (int pass = 0; pass < texel_pass_num; pass++, texel+=SOFTMASK_BLOCK_SIZE) // loop the available warps enough to complete all translations for this orientation
		{
			if(texel<vol_size)
			{
				img_pixels[tid]=__ldg(&vol[texel]);

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
					partial_sum[tid]    += (XFLOAT)1.0;
					partial_sum_bg[tid] += img_pixels[tid];
				}
				else
				{
#if defined(ACC_DOUBLE_PRECISION)
					raisedcos = 0.5 + 0.5  * cospi( (radius_p - r) / cosine_width );
#else
					raisedcos = 0.5f + 0.5f * cospif((radius_p - r) / cosine_width );
#endif
					partial_sum[tid] += raisedcos;
					partial_sum_bg[tid] += raisedcos * img_pixels[tid];
				}
			}
		}

		cuda_atomic_add(&g_sum[tid]   , partial_sum[tid]);
		cuda_atomic_add(&g_sum_bg[tid], partial_sum_bg[tid]);
}


__global__ void cuda_kernel_cosineFilter(	XFLOAT *vol,
											long int vol_size,
											long int xdim,
											long int ydim,
											long int zdim,
											long int xinit,
											long int yinit,
											long int zinit,
											bool do_noise,
											XFLOAT *noise,
											XFLOAT radius,
											XFLOAT radius_p,
											XFLOAT cosine_width,
											XFLOAT bg_value)
{

	int tid = threadIdx.x;
	int bid = blockIdx.x;

//		vol.setXmippOrigin(); // sets xinit=xdim , also for y z
	XFLOAT r, raisedcos, defVal;
	int x,y,z;
	__shared__ XFLOAT     img_pixels[SOFTMASK_BLOCK_SIZE];

	long int texel_pass_num = ceilfracf(vol_size,SOFTMASK_BLOCK_SIZE*gridDim.x);
	int texel = bid*SOFTMASK_BLOCK_SIZE*texel_pass_num + tid;

	defVal = bg_value;
	for (int pass = 0; pass < texel_pass_num; pass++, texel+=SOFTMASK_BLOCK_SIZE) // loop the available warps enough to complete all translations for this orientation
	{
		if(texel<vol_size)
		{
			img_pixels[tid]=__ldg(&vol[texel]);

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
				img_pixels[tid]=defVal;
			else
			{
#if defined(ACC_DOUBLE_PRECISION)
				raisedcos = 0.5  + 0.5  * cospi( (radius_p - r) / cosine_width );
#else
				raisedcos = 0.5f + 0.5f * cospif((radius_p - r) / cosine_width );
#endif
				img_pixels[tid]= img_pixels[tid]*(1-raisedcos) + defVal*raisedcos;

			}
			vol[texel]=img_pixels[tid];
		}

	}
}

__global__ void cuda_kernel_centerFFT_2D(XFLOAT *img_in,
										 int image_size,
										 int xdim,
										 int ydim,
										 int xshift,
										 int yshift)
{

	__shared__ XFLOAT buffer[CFTT_BLOCK_SIZE];
	int tid = threadIdx.x;
	int pixel = threadIdx.x + blockIdx.x*CFTT_BLOCK_SIZE;
	long int image_offset = image_size*blockIdx.y;
//	int pixel_pass_num = ceilfracf(image_size, CFTT_BLOCK_SIZE);

//	for (int pass = 0; pass < pixel_pass_num; pass++, pixel+=CFTT_BLOCK_SIZE)
//	{
		if(pixel<(image_size/2))
		{
			int y = floorf((XFLOAT)pixel/(XFLOAT)xdim);
			int x = pixel % xdim;				// also = pixel - y*xdim, but this depends on y having been calculated, i.e. serial evaluation

			int yp = y + yshift;
			if (yp < 0)
				yp += ydim;
			else if (yp >= ydim)
				yp -= ydim;

			int xp = x + xshift;
			if (xp < 0)
				xp += xdim;
			else if (xp >= xdim)
				xp -= xdim;

			int n_pixel = yp*xdim + xp;

			buffer[tid]                    = img_in[image_offset + n_pixel];
			img_in[image_offset + n_pixel] = img_in[image_offset + pixel];
			img_in[image_offset + pixel]   = buffer[tid];
		}
//	}
}

__global__ void cuda_kernel_centerFFT_3D(XFLOAT *img_in,
										 int image_size,
										 int xdim,
										 int ydim,
										 int zdim,
										 int xshift,
										 int yshift,
									 	 int zshift)
{

	__shared__ XFLOAT buffer[CFTT_BLOCK_SIZE];
	int tid = threadIdx.x;
	int pixel = threadIdx.x + blockIdx.x*CFTT_BLOCK_SIZE;
	long int image_offset = image_size*blockIdx.y;

		int xydim = xdim*ydim;
		if(pixel<(image_size/2))
		{
			int z = floorf((XFLOAT)pixel/(XFLOAT)(xydim));
			int xy = pixel % xydim;
			int y = floorf((XFLOAT)xy/(XFLOAT)xdim);
			int x = xy % xdim;

			int xp = x + xshift;
			if (xp < 0)
				xp += xdim;
			else if (xp >= xdim)
				xp -= xdim;

			int yp = y + yshift;
			if (yp < 0)
				yp += ydim;
			else if (yp >= ydim)
				yp -= ydim;

			int zp = z + zshift;
			if (zp < 0)
				zp += zdim;
			else if (zp >= zdim)
				zp -= zdim;

			int n_pixel = zp*xydim + yp*xdim + xp;

			buffer[tid]                    = img_in[image_offset + n_pixel];
			img_in[image_offset + n_pixel] = img_in[image_offset + pixel];
			img_in[image_offset + pixel]   = buffer[tid];
		}
}


__global__ void cuda_kernel_probRatio(  XFLOAT *d_Mccf,
										XFLOAT *d_Mpsi,
										XFLOAT *d_Maux,
										XFLOAT *d_Mmean,
										XFLOAT *d_Mstddev,
										int image_size,
										XFLOAT normfft,
										XFLOAT sum_ref_under_circ_mask,
										XFLOAT sum_ref2_under_circ_mask,
										XFLOAT expected_Pratio,
										int NpsiThisBatch,
										int startPsi,
										int totalPsis)
{
	/* PLAN TO:
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
	 */

	int pixel = threadIdx.x + blockIdx.x*(int)PROBRATIO_BLOCK_SIZE;
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

__global__ void cuda_kernel_rotateOnly(   ACCCOMPLEX *d_Faux,
						  	  	  	  	  XFLOAT psi,
						  	  			  AccProjectorKernel projector,
						  	  			  int startPsi
						  	  			  )
{
	int proj = blockIdx.y;
	int image_size=projector.imgX*projector.imgY;
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
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
		sincos((proj+startPsi)*psi, &sa, &ca);
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

__global__ void cuda_kernel_rotateAndCtf( ACCCOMPLEX *d_Faux,
						  	  	  	  	  XFLOAT *d_ctf,
						  	  	  	  	  XFLOAT psi,
						  	  			  AccProjectorKernel projector,
						  	  			  int startPsi
						  	  			  )
{
	int proj = blockIdx.y;
	int image_size=projector.imgX*projector.imgY;
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
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
		sincos((proj+startPsi)*psi, &sa, &ca);
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


__global__ void cuda_kernel_convol_A( ACCCOMPLEX *d_A,
									 ACCCOMPLEX *d_B,
									 int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
	{
		XFLOAT tr =   d_A[pixel].x;
		XFLOAT ti = - d_A[pixel].y;
		d_A[pixel].x =   tr*d_B[pixel].x - ti*d_B[pixel].y;
		d_A[pixel].y =   ti*d_B[pixel].x + tr*d_B[pixel].y;
	}
}

__global__ void cuda_kernel_convol_A( ACCCOMPLEX *d_A,
									 ACCCOMPLEX *d_B,
									 ACCCOMPLEX *d_C,
									 int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
	{
		XFLOAT tr =   d_A[pixel].x;
		XFLOAT ti = - d_A[pixel].y;
		d_C[pixel].x =   tr*d_B[pixel].x - ti*d_B[pixel].y;
		d_C[pixel].y =   ti*d_B[pixel].x + tr*d_B[pixel].y;
	}
}

__global__ void cuda_kernel_batch_convol_A( ACCCOMPLEX *d_A,
									 	 	ACCCOMPLEX *d_B,
									 	 	int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	int A_off = blockIdx.y*image_size;
	if(pixel<image_size)
	{
		XFLOAT tr =   d_A[pixel + A_off].x;
		XFLOAT ti = - d_A[pixel + A_off].y;
		d_A[pixel + A_off].x =   tr*d_B[pixel].x - ti*d_B[pixel].y;
		d_A[pixel + A_off].y =   ti*d_B[pixel].x + tr*d_B[pixel].y;
	}
}

__global__ void cuda_kernel_batch_convol_A( ACCCOMPLEX *d_A,
									 	 	ACCCOMPLEX *d_B,
									 	 	ACCCOMPLEX *d_C,
									 	 	int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	int A_off = blockIdx.y*image_size;
	if(pixel<image_size)
	{
		XFLOAT tr =   d_A[pixel + A_off].x;
		XFLOAT ti = - d_A[pixel + A_off].y;
		d_C[pixel + A_off].x =   tr*d_B[pixel].x - ti*d_B[pixel].y;
		d_C[pixel + A_off].y =   ti*d_B[pixel].x + tr*d_B[pixel].y;
	}
}

__global__ void cuda_kernel_convol_B(	 ACCCOMPLEX *d_A,
									 	 ACCCOMPLEX *d_B,
									 	 int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
	{
		XFLOAT tr = d_A[pixel].x;
		XFLOAT ti = d_A[pixel].y;
		d_A[pixel].x =   tr*d_B[pixel].x + ti*d_B[pixel].y;
		d_A[pixel].y =   ti*d_B[pixel].x - tr*d_B[pixel].y;
	}
}

__global__ void cuda_kernel_convol_B(	 ACCCOMPLEX *d_A,
									 	 ACCCOMPLEX *d_B,
									 	 ACCCOMPLEX *d_C,
									 	 int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
	{
		XFLOAT tr = d_A[pixel].x;
		XFLOAT ti = d_A[pixel].y;
		d_C[pixel].x =   tr*d_B[pixel].x + ti*d_B[pixel].y;
		d_C[pixel].y =   ti*d_B[pixel].x - tr*d_B[pixel].y;
	}
}

__global__ void cuda_kernel_batch_convol_B(	 ACCCOMPLEX *d_A,
									 	 	 ACCCOMPLEX *d_B,
									 	 	 int image_size)
{
	long int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	int A_off = blockIdx.y*image_size;
	if(pixel<image_size)
	{
		XFLOAT tr = d_A[pixel + A_off].x;
		XFLOAT ti = d_A[pixel + A_off].y;
		d_A[pixel + A_off].x =   tr*d_B[pixel].x + ti*d_B[pixel].y;
		d_A[pixel + A_off].y =   ti*d_B[pixel].x - tr*d_B[pixel].y;
	}
}

__global__ void cuda_kernel_batch_multi( XFLOAT *A,
								   XFLOAT *B,
								   XFLOAT *OUT,
								   XFLOAT S,
		  	  	  	  	  	  	   int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
		OUT[pixel + blockIdx.y*image_size] = A[pixel + blockIdx.y*image_size]*B[pixel + blockIdx.y*image_size]*S;
}

__global__ void cuda_kernel_finalizeMstddev( XFLOAT *Mstddev,
											 XFLOAT *aux,
											 XFLOAT S,
											 int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
	{
		XFLOAT temp = Mstddev[pixel] + S * aux[pixel];
		if(temp > 0)
			Mstddev[pixel] = sqrt(temp);
		else
			Mstddev[pixel] = 0;
	}
}

__global__ void cuda_kernel_square(
		XFLOAT *A,
		int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
		A[pixel] = A[pixel]*A[pixel];
}

template<bool invert>
__global__ void cuda_kernel_make_eulers_2D(
		XFLOAT *alphas,
		XFLOAT *eulers,
		unsigned orientation_num)
{
	unsigned oid = blockIdx.x * BLOCK_SIZE + threadIdx.x; //Orientation id

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
}

template<bool invert, bool doL, bool doR>
__global__ void cuda_kernel_make_eulers_3D(
		XFLOAT *alphas,
		XFLOAT *betas,
		XFLOAT *gammas,
		XFLOAT *eulers,
		unsigned orientation_num,
		XFLOAT *L,
		XFLOAT *R)
{
	XFLOAT a(0.f),b(0.f),g(0.f), A[9],B[9];
	XFLOAT ca, sa, cb, sb, cg, sg, cc, cs, sc, ss;

	unsigned oid = blockIdx.x * BLOCK_SIZE + threadIdx.x; //Orientation id

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
}

__global__ void cuda_kernel_allweights_to_mweights(
		unsigned long * d_iorient,
		XFLOAT * d_allweights,
		XFLOAT * d_mweights,
		unsigned long orientation_num,
		unsigned long translation_num,
        int block_size
		)
{
	size_t idx = blockIdx.x * block_size + threadIdx.x;
	if (idx < orientation_num*translation_num)
		d_mweights[d_iorient[idx/translation_num] * translation_num + idx%translation_num] =
				d_allweights[idx/translation_num * translation_num + idx%translation_num];
                // TODO - isn't this just d_allweights[idx + idx%translation_num]?   Really?
}
