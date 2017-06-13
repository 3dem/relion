#include "src/acc/cpu/settings.h"
#include "src/acc/cpu/utilities.h"
#include "src/acc/cpu/helper.h"


namespace CpuKernels
{

/*
 * This draft of a kernel assumes input that has jobs which have a single orientation and sequential translations within each job.
 * 
 */
void exponentiate_weights_fine(
        int blockIdx_x,
        int threadIdx_x,
		XFLOAT *g_pdf_orientation,
		XFLOAT *g_pdf_offset,
		XFLOAT *g_weights,
		XFLOAT avg_diff2,
		int oversamples_orient,
		int oversamples_trans,
		unsigned long *d_rot_id,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num,
		long int job_num)
{
    // XFLOAT s_weights[SUMW_BLOCK_SIZE];
	XFLOAT s_weights;

	// blockid
	int bid  = blockIdx_x;
	//threadid
	int tid = threadIdx_x;

	long int jobid = bid*SUMW_BLOCK_SIZE+tid;

	if (jobid<job_num)
	{
		long int pos = d_job_idx[jobid];
		// index of comparison
		long int ix =  d_rot_id[   pos];   // each thread gets its own orient...
		long int iy = d_trans_idx[ pos];   // ...and it's starting trans...
		long int in =  d_job_num[jobid];    // ...AND the number of translations to go through

		int c_itrans;//, iorient = bid*SUM_BLOCK_SIZE+tid; //, f_itrans;

		// Because the portion of work is so arbitrarily divided in this kernel,
		// we need to do some brute idex work to get the correct indices.
		for (int itrans=0; itrans < in; itrans++, iy++)
		{
			c_itrans = ( iy - (iy % oversamples_trans))/ oversamples_trans; //floor(x/y) == (x-(x%y))/y  but less sensitive to x>>y and finite precision
//			f_itrans = iy % oversamples_trans;

			XFLOAT prior = g_pdf_orientation[ix] * g_pdf_offset[c_itrans];          	// Same      for all threads - TODO: should be done once for all trans through warp-parallel execution
			XFLOAT diff2 = g_weights[pos+itrans] - avg_diff2;								// Different for all threads
			// next line because of numerical precision of exp-function
	#if defined(ACC_DOUBLE_PRECISION)
				if (diff2 > 700.)
					s_weights = 0.;
				else
					s_weights = prior * exp(-diff2);
	#else
				if (diff2 > 88.)
					s_weights = 0.f;
				else
					s_weights = prior * expf(-diff2);
	#endif
				// TODO: use tabulated exp function? / Sjors  TODO: exp, expf, or __exp in CUDA? /Bjorn
			// Store the weight
			g_weights[pos+itrans] = s_weights; // TODO put in shared mem
		}
	}
}

void SoftMaskBackgroundValue(	int      blockIdx_x,
                                int      threadIdx_x,
                                int      gridDim_x,
                                XFLOAT  *vol,
								long int vol_size,
								long int xdim,
								long int ydim,
								long int zdim,
								long int xinit,
								long int yinit,
								long int zinit,
								bool     do_Mnoise,
								XFLOAT   radius,
								XFLOAT   radius_p,
								XFLOAT   cosine_width,
								XFLOAT  *g_sum,
								XFLOAT  *g_sum_bg)
{

		int tid = threadIdx_x;
		int bid = blockIdx_x;

//		vol.setXmippOrigin(); // sets xinit=xdim , also for y z
		XFLOAT r, raisedcos;
		int x,y,z;
		XFLOAT     img_pixels;

		long int texel_pass_num = ceilfracf(vol_size,SOFTMASK_BLOCK_SIZE*gridDim_x);
		int texel = bid*SOFTMASK_BLOCK_SIZE*texel_pass_num + tid;

		for (int pass = 0; pass < texel_pass_num; pass++, texel+=SOFTMASK_BLOCK_SIZE) // loop the available warps enough to complete all translations for this orientation
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
}


void cosineFilter(	int      blockIdx_x,
					int      threadIdx_x,
					int      gridDim_x,
					XFLOAT  *vol,
					long int vol_size,
					long int xdim,
					long int ydim,
					long int zdim,
					long int xinit,
					long int yinit,
					long int zinit,
					bool     do_Mnoise,
					XFLOAT   radius,
					XFLOAT   radius_p,
					XFLOAT   cosine_width,
					XFLOAT   bg_value)
{

	int tid = threadIdx_x;
	int bid = blockIdx_x;

//		vol.setXmippOrigin(); // sets xinit=xdim , also for y z
	XFLOAT r, raisedcos;
	int x,y,z;
	XFLOAT     img_pixels;

	long int texel_pass_num = ceilfracf(vol_size,SOFTMASK_BLOCK_SIZE*gridDim_x);
	int texel = bid*SOFTMASK_BLOCK_SIZE*texel_pass_num + tid;

	for (int pass = 0; pass < texel_pass_num; pass++, texel+=SOFTMASK_BLOCK_SIZE) // loop the available warps enough to complete all translations for this orientation
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

			if (r < radius)
				continue;
			else if (r > radius_p)
				img_pixels=bg_value;
			else
			{
#if defined(ACC_DOUBLE_PRECISION)
				raisedcos = 0.5 + 0.5  * cos ( (radius_p - r) / cosine_width * M_PI);
#else
				raisedcos = 0.5 + 0.5  * cosf( (radius_p - r) / cosine_width * M_PI);
#endif
				img_pixels= img_pixels*(1-raisedcos) + bg_value*raisedcos;

			}
			vol[texel]=img_pixels;
		}

	}
}
void translate2D(   int      blockIdx_x,
					int      threadIdx_x,
					XFLOAT * g_image_in,
					XFLOAT * g_image_out,
					int      image_size,
					int      xdim,
					int      ydim,
					int      dx,
					int      dy)
{
	int tid = threadIdx_x;
	int bid =  blockIdx_x;

	int x,y,xp,yp;
	int pixel=tid + bid*BLOCK_SIZE;
	int new_pixel;

	if(pixel<image_size)
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

void translate3D(   int      blockIdx_x,	
					int      threadIdx_x,
					XFLOAT * g_image_in,
					XFLOAT * g_image_out,
					int      image_size,
					int      xdim,
					int      ydim,
					int      zdim,
					int      dx,
					int      dy,
					int      dz)
{
	int tid = threadIdx_x;
	int bid =  blockIdx_x;

	int x,y,z,xp,yp,zp,xy;
	int voxel=tid + bid*BLOCK_SIZE;
	int new_voxel;

	int xydim = xdim*ydim;

	if(voxel<image_size)
	{
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

void centerFFT_2D(  int     blockIdx_x,
					int     blockIdx_y,
					int     threadIdx_x,
					XFLOAT *img_in,
					int     image_size,
					int     xdim,
					int     ydim,
					int     xshift,
					int     yshift)
{

	XFLOAT buffer;
	int tid = threadIdx_x;
	int pixel = threadIdx_x + blockIdx_x*CFTT_BLOCK_SIZE;
	long int image_offset = image_size*blockIdx_y;
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

			buffer                         = img_in[image_offset + n_pixel];
			img_in[image_offset + n_pixel] = img_in[image_offset + pixel];
			img_in[image_offset + pixel]   = buffer;
		}
//	}
}

void centerFFT_3D(  int       blockIdx_x,
					int       blockIdx_y,
					int       threadIdx_x,
					XFLOAT   *img_in,
					int       image_size,
					int       xdim,
					int       ydim,
					int       zdim,
					int       xshift,
					int       yshift,
					int       zshift)
{

	XFLOAT buffer;
	int tid = threadIdx_x;
	int pixel = threadIdx_x + blockIdx_x*CFTT_BLOCK_SIZE;
	long int image_offset = image_size*blockIdx_y;

		int xydim = xdim*ydim;
		if(pixel<(image_size/2))
		{
			int z = floorf((XFLOAT)pixel/(XFLOAT)(xydim));
			int xy = pixel % xydim;
			int y = floorf((XFLOAT)xy/(XFLOAT)xdim);
			int x = xy % xdim;

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

			int zp = z + zshift;
			if (zp < 0)
				zp += zdim;
			else if (zp >= zdim)
				zp -= zdim;

			int n_pixel = zp*xydim + yp*xdim + xp;

			buffer                         = img_in[image_offset + n_pixel];
			img_in[image_offset + n_pixel] = img_in[image_offset + pixel];
			img_in[image_offset + pixel]   = buffer;
		}
}


void probRatio( int       blockIdx_x,  
				int       threadIdx_x,
				XFLOAT   *d_Mccf,
				XFLOAT   *d_Mpsi,
				XFLOAT   *d_Maux,
				XFLOAT   *d_Mmean,
				XFLOAT   *d_Mstddev,
				int       image_size,
				XFLOAT    normfft,
				XFLOAT    sum_ref_under_circ_mask,
				XFLOAT    sum_ref2_under_circ_mask,
				XFLOAT    expected_Pratio,
				int       NpsiThisBatch,
				int       startPsi,
				int       totalPsis)
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

	int pixel = threadIdx_x + blockIdx_x*(int)PROBRATIO_BLOCK_SIZE;
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
				ProjectorKernel &projector,
				int              startPsi
		       )
{
	int proj = blockIdx_y;
	int image_size=projector.imgX*projector.imgY;
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
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
					ProjectorKernel &projector,
					int       startPsi
				)
{
	int proj = blockIdx_y;
	int image_size=projector.imgX*projector.imgY;
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
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
				int           image_size)
{
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
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
				int          image_size)
{
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
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
					int           image_size)
{
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
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
					int           image_size)
{
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
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
				int          image_size)
{
	int pixel =  threadIdx_x + blockIdx_x*BLOCK_SIZE;
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
				int           image_size)
{
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
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
					int          image_size)
{
	long int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
	int A_off = blockIdx_y*image_size;
	if(pixel<image_size)
	{
		XFLOAT tr = d_A[pixel + A_off].x;
		XFLOAT ti = d_A[pixel + A_off].y;
		d_A[pixel + A_off].x =   tr*d_B[pixel].x + ti*d_B[pixel].y;
		d_A[pixel + A_off].y =   ti*d_B[pixel].x - tr*d_B[pixel].y;
	}
}

void multi( int     blockIdx_x,
			int     threadIdx_x,
			XFLOAT *A,
			XFLOAT *OUT,
			XFLOAT  S,
			int     image_size)
{
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
	if(pixel<image_size)
		OUT[pixel] = A[pixel]*S;
}

void multi( int     blockIdx_x,
			int     threadIdx_x,
			XFLOAT *A,
			XFLOAT  S,
			int     image_size)
{
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
	if(pixel<image_size)
		A[pixel] = A[pixel]*S;
}

void multi( int     blockIdx_x,
			int     threadIdx_x,
			XFLOAT *A,
			XFLOAT *B,
			XFLOAT *OUT,
			XFLOAT  S,
			int     image_size)
{
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
	if(pixel<image_size)
		OUT[pixel] = A[pixel]*B[pixel]*S;
}

void batch_multi(   int     blockIdx_x,
					int     blockIdx_y,
					int     threadIdx_x,
					XFLOAT *A,
					XFLOAT *B,
					XFLOAT *OUT,
					XFLOAT  S,
					int     image_size)
{
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
	if(pixel<image_size)
		OUT[pixel + blockIdx_y*image_size] = A[pixel + blockIdx_y*image_size]*B[pixel + blockIdx_y*image_size]*S;
}

void finalizeMstddev(   int     blockIdx_x,
						int     threadIdx_x,
						XFLOAT *Mstddev,
						XFLOAT *aux,
						XFLOAT  S,
						int     image_size)
{
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
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
			int     image_size)
{
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
	if(pixel<image_size)
		A[pixel] = A[pixel]*A[pixel];
}

} // end of namespace CpuKernels

