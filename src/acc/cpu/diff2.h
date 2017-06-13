#ifndef DIFF2_KERNELS_H_
#define DIFF2_KERNELS_H_

#include <vector>
#include <iostream>
#include <fstream>

#include "src/acc/cpu/settings.h"
#include "src/acc/cpu/projector.h"
#include "src/acc/cpu/utilities.h"
#include "src/acc/cpu/helper.h"

namespace CpuKernels
{

/*
 *   	DIFFERENCE-BASED KERNELS
 */
template<bool REF3D, bool DATA3D, int block_sz, int eulers_per_block, int prefetch_fraction>
void diff2_coarse(                    
        int     blockIdx_x,
		XFLOAT *g_eulers,
		XFLOAT *trans_x,
		XFLOAT *trans_y,
		XFLOAT *trans_z,		
		XFLOAT *g_real,
		XFLOAT *g_imag,
		ProjectorKernel &projector,
		XFLOAT *g_corr,
		XFLOAT *g_diff2s,
		int translation_num,
		int image_size
		)
{           
    //Prefetch euler matrices
	XFLOAT s_eulers[eulers_per_block * 9];

    for (int i = 0; i < eulers_per_block * 9; i++)
        s_eulers[i] = g_eulers[blockIdx_x * eulers_per_block * 9 + i];

    //Setup variables
	XFLOAT s_ref_real[eulers_per_block][block_sz];
	XFLOAT s_ref_imag[eulers_per_block][block_sz];

	XFLOAT s_real[block_sz];
	XFLOAT s_imag[block_sz];
	XFLOAT s_corr[block_sz]; 
	
	int    x[block_sz], y[block_sz], z[block_sz];

    XFLOAT diff2s[translation_num][eulers_per_block];
    memset(&diff2s[0][0], 0, sizeof(XFLOAT) * translation_num * eulers_per_block);
                    
	//Step through data
	unsigned pass_num(ceilfracf(image_size,block_sz));
	for (unsigned pass = 0; pass < pass_num; pass++) { // finish an entire ref image each block
        int start = pass * block_sz;
        
        // Rotate the reference image per block_sz, saved in cache
        #pragma simd
	    for (int tid=0; tid<block_sz; tid++){
			int pixel = start + tid;
			if(pixel >= image_size)
			    continue;

    		if(DATA3D)
    		{
    			z[tid] =  floorfracf(pixel, projector.imgX*projector.imgY);
    			int xy = pixel % (projector.imgX*projector.imgY);
    			x[tid] =             xy  % projector.imgX;
    			y[tid] = floorfracf( xy,   projector.imgX);
    			if (z[tid] > projector.maxR)
    				z[tid] -= projector.imgZ;
    		}
    		else
    		{
    			x[tid] =            pixel % projector.imgX;
    			y[tid] = floorfracf(pixel, projector.imgX);
    		}
    		if (y[tid] > projector.maxR)
    			y[tid] -= projector.imgY;

            for (int i = 0; i < eulers_per_block; i ++) {
                if(DATA3D) // if DATA3D, then REF3D as well.
                    projector.project3Dmodel(
                            x[tid], y[tid], z[tid],
                            s_eulers[i*9  ],
    						s_eulers[i*9+1],
	    					s_eulers[i*9+2],
	    					s_eulers[i*9+3],
	    					s_eulers[i*9+4],
	    					s_eulers[i*9+5],
	    					s_eulers[i*9+6],
	    					s_eulers[i*9+7],
	    					s_eulers[i*9+8],
                            s_ref_real[i][tid],
                            s_ref_imag[i][tid]);
                else if(REF3D)
    				projector.project3Dmodel(
                            x[tid], y[tid], 
	    					s_eulers[i*9  ],
	    					s_eulers[i*9+1],
	    					s_eulers[i*9+3],
	    					s_eulers[i*9+4],
	    					s_eulers[i*9+6],
	    					s_eulers[i*9+7],
	    					s_ref_real[i][tid],
	    					s_ref_imag[i][tid]);                    
                else
                    projector.project2Dmodel(
                            x[tid], y[tid], 
                            s_eulers[i*9  ],
                            s_eulers[i*9+1],
                            s_eulers[i*9+3],
                            s_eulers[i*9+4],
                            s_ref_real[i][tid],
                            s_ref_imag[i][tid]);
            }

            s_real[tid] = g_real[pixel];
            s_imag[tid] = g_imag[pixel];
            s_corr[tid] = g_corr[pixel] / 2;
        }

        for(int i=0; i<translation_num; i++) {
            XFLOAT tx = trans_x[i];
            XFLOAT ty = trans_y[i];
            XFLOAT tz = trans_z[i];                 
            
            #pragma simd
            for (int tid=0; tid<block_sz; tid++) {
                int pixel = start + tid;
    			if(pixel >= image_size)
	    		    continue;                

                XFLOAT real, imag;
			    if(DATA3D)
				    translatePixel(x[tid], y[tid], z[tid], tx, ty, tz, s_real[tid], s_imag[tid], real, imag);
    			else
	    			translatePixel(x[tid], y[tid],         tx, ty,     s_real[tid], s_imag[tid], real, imag);

                for (int j = 0; j < eulers_per_block; j ++) {
                    XFLOAT diff_real =  s_ref_real[j][tid] - real;
                    XFLOAT diff_imag =  s_ref_imag[j][tid] - imag;
                    
                    diff2s[i][j] += (diff_real * diff_real + diff_imag * diff_imag) * s_corr[tid];
                }             
            } // for tid       
        }  // for each translation
    }  // for each pass
    
    XFLOAT *pData = g_diff2s + blockIdx_x * eulers_per_block * translation_num;
    for(int i=0; i<eulers_per_block; i++) {
        for(int j=0; j<translation_num; j++) {
             *pData += diff2s[j][i];
             pData ++;
        }
    }
}

template<bool REF3D, bool DATA3D, int block_sz,int chunk_sz>
void diff2_fine(
        int     blockIdx_x,
		XFLOAT *g_eulers,
		XFLOAT *g_imgs_real,
		XFLOAT *g_imgs_imag,
		XFLOAT *trans_x,
		XFLOAT *trans_y,
		XFLOAT *trans_z,		
		ProjectorKernel &projector,
		XFLOAT *g_corr_img,
		XFLOAT *g_diff2s,
		unsigned image_size,
		XFLOAT sum_init,
		unsigned long orientation_num,
		unsigned long translation_num,
		unsigned long todo_blocks,
		unsigned long *d_rot_idx,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num
		)
{
    
	unsigned long bid = blockIdx_x;

//    // Specialize BlockReduce for a 1D block of 128 threads on type XFLOAT
//    typedef cub::BlockReduce<XFLOAT, 128> BlockReduce;
//    // Allocate shared memory for BlockReduce
//     typename BlockReduce::TempStorage temp_storage;

	unsigned long pixel;
	XFLOAT ref_real, ref_imag,
		shifted_real, shifted_imag,
		diff_real, diff_imag;


    unsigned trans_num  = (unsigned)d_job_num[bid]; 
    
	XFLOAT s[trans_num][block_sz];    
    memset(&s[0][0], 0, sizeof(XFLOAT) * block_sz * trans_num);
    
    unsigned long int ix = d_rot_idx[d_job_idx[bid]];
	unsigned long int iy;
    unsigned long int iy_part = d_trans_idx[d_job_idx[bid]];    
    
	// index of comparison
	unsigned pass_num(ceilfracf(image_size,block_sz));
	for (unsigned pass = 0; pass < pass_num; pass++) { // finish an entire ref image each block
        #pragma simd 
	    for (int tid=0; tid<block_sz; tid++){
			pixel = (pass * block_sz) + tid;

			if(pixel >= image_size)
			    continue;
			    
            int x,y,z,xy;
			if(DATA3D)
			{
				z =  floorfracf(pixel, projector.imgX*projector.imgY);
				xy = pixel % (projector.imgX*projector.imgY);
				x =             xy  % projector.imgX;
				y = floorfracf( xy,   projector.imgX);
				if (z > projector.maxR)
				{
					if (z >= projector.imgZ - projector.maxR)
						z = z - projector.imgZ;
					else
						x = projector.maxR;
				}
			}
			else
			{
				x =             pixel % projector.imgX;
				y = floorfracf( pixel , projector.imgX);
			}
			if (y > projector.maxR)
			{
				if (y >= projector.imgY - projector.maxR)
					y = y - projector.imgY;
				else
					x = projector.maxR;
			}

			if(DATA3D)
				projector.project3Dmodel(
					x,y,z,
					g_eulers[ix*9  ], g_eulers[ix*9+1], g_eulers[ix*9+2],
					g_eulers[ix*9+3], g_eulers[ix*9+4], g_eulers[ix*9+5],
					g_eulers[ix*9+6], g_eulers[ix*9+7], g_eulers[ix*9+8],
					ref_real, ref_imag);
			else if(REF3D)
				projector.project3Dmodel(
					x,y,
					g_eulers[ix*9  ], g_eulers[ix*9+1],
					g_eulers[ix*9+3], g_eulers[ix*9+4],
					g_eulers[ix*9+6], g_eulers[ix*9+7],
					ref_real, ref_imag);			    
			else
				projector.project2Dmodel(
					x,y,
					g_eulers[ix*9  ], g_eulers[ix*9+1],
					g_eulers[ix*9+3], g_eulers[ix*9+4],
					ref_real, ref_imag);

            XFLOAT imgs_real = g_imgs_real[pixel];
            XFLOAT imgs_imag = g_imgs_imag[pixel];
                
			for (int itrans=0; itrans<trans_num; itrans++) // finish all translations in each partial pass
			{
				iy = iy_part + itrans;

				translatePixel(x, y, trans_x[iy], trans_y[iy], imgs_real, imgs_imag, shifted_real, shifted_imag);

				diff_real =  ref_real - shifted_real;
				diff_imag =  ref_imag - shifted_imag;
				
			    s[itrans][tid] += (diff_real * diff_real + diff_imag * diff_imag) * (XFLOAT)0.5 * g_corr_img[pixel];
			}
		} // for tid
    }  // for pass

    for (int itrans=0; itrans<trans_num; itrans++)
    {
        XFLOAT sum = (XFLOAT)0.0;
        for(int tid=0; tid<block_sz; tid++)
            sum += s[itrans][tid];
		unsigned long int iy = d_job_idx[bid]+itrans;
        g_diff2s[iy] = sum + sum_init;
    }
}


/*
 *   	CROSS-CORRELATION-BASED KERNELS
 */

template<bool REF3D, bool DATA3D, int block_sz>
 void diff2_CC_coarse(
        int blockIdx_x,
		XFLOAT *g_eulers,
		XFLOAT *g_imgs_real,
		XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x,
		XFLOAT *g_trans_y,
		XFLOAT *g_trans_z,	
		ProjectorKernel &projector,
		XFLOAT *g_corr_img,
		unsigned translation_num,
		int      image_size,
		XFLOAT   exp_local_sqrtXi2,
		XFLOAT  *g_diff2
		)
{
	int iorient = blockIdx_x;

	XFLOAT real, imag, ref_real, ref_imag;

	XFLOAT e0,e1,e2,e3,e4,e5,e6,e7,e8;
	e0 = g_eulers[iorient*9  ];
	e1 = g_eulers[iorient*9+1];
	e2 = g_eulers[iorient*9+2];
	e3 = g_eulers[iorient*9+3];
	e4 = g_eulers[iorient*9+4];
	e5 = g_eulers[iorient*9+5];
	e6 = g_eulers[iorient*9+6];
	e7 = g_eulers[iorient*9+7];
	e8 = g_eulers[iorient*9+8];
 
    XFLOAT s_weight[translation_num][block_sz];
    memset(s_weight, 0, sizeof(XFLOAT) * block_sz * translation_num);
    XFLOAT s_norm[translation_num][block_sz];
    memset(s_norm, 0, sizeof(XFLOAT) * block_sz * translation_num);

	unsigned pixel_pass_num( ceilfracf(image_size,block_sz) );
	for (unsigned pass = 0; pass < pixel_pass_num; pass++)
	{	    	
	    #pragma simd    
		for (int tid=0; tid<block_sz; tid++) {
    		unsigned pixel = (pass * block_sz) + tid;

            // TODO: if image_size is divisible by block_sz, the following
            //  if statement can be removed.
	    	if(pixel >= image_size)
		        continue;
		
			int x,y,z,xy;
			if(DATA3D)
			{
				z =  floorfracf(pixel, projector.imgX*projector.imgY);
				xy = pixel % (projector.imgX*projector.imgY);
				x =             xy  % projector.imgX;
				y = floorfracf( xy,   projector.imgX);
				if (z > projector.maxR)
				{
					if (z >= projector.imgZ - projector.maxR)
						z = z - projector.imgZ;
					else
						x = projector.maxR;
				}
			}
			else
			{
				x =             pixel % projector.imgX;
				y = floorfracf( pixel , projector.imgX);
			}
			if (y > projector.maxR)
			{
				if (y >= projector.imgY - projector.maxR)
					y = y - projector.imgY;
				else
					x = projector.maxR;
			}

			if(DATA3D)
				projector.project3Dmodel(
					x,y,z,
					e0,e1,e2,e3,e4,e5,e6,e7,e8,
					ref_real, ref_imag);
			else if(REF3D)
				projector.project3Dmodel(
					x,y,
					e0,e1,e3,e4,e6,e7,
					ref_real, ref_imag);
			else
				projector.project2Dmodel(
					x,y,
					e0,e1,e3,e4,
					ref_real, ref_imag);

            for(int itrans=0; itrans<translation_num; itrans++) {
    			if(DATA3D)
	    			translatePixel(x, y, z, g_trans_x[itrans], g_trans_y[itrans], g_trans_z[itrans], 
	    			               g_imgs_real[pixel], g_imgs_imag[pixel], real, imag);
		    	else
			    	translatePixel(x, y,    g_trans_x[itrans], g_trans_y[itrans],                   
			    	               g_imgs_real[pixel], g_imgs_imag[pixel], real, imag);

			    s_weight[itrans][tid] += (ref_real * real     + ref_imag * imag)      * g_corr_img[pixel];
			    s_norm  [itrans][tid] += (ref_real * ref_real + ref_imag * ref_imag ) * g_corr_img[pixel];
			}
		}
	}
	
	for(int itrans=0; itrans<translation_num; itrans++) {
    	XFLOAT sum_weight = (XFLOAT)0.0;
	    XFLOAT sum_norm   = (XFLOAT)0.0;		
		
    	for(int i=0; i<block_sz; i++){
	        sum_weight += s_weight[itrans][i];
	        sum_norm   += s_norm  [itrans][i];
    	}
		
#ifdef RELION_SINGLE_PRECISION                  
        g_diff2[itrans] = - ( sum_weight / sqrtf(sum_norm));
#else                   
        g_diff2[itrans] = - ( sum_weight / sqrt(sum_norm));
#endif
    }
}

template<bool REF3D, bool DATA3D, int block_sz>
void diff2_CC_fine(
        int blockIdx_x,
		XFLOAT *g_eulers,
		XFLOAT *g_imgs_real,
		XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x,
		XFLOAT *g_trans_y,
		XFLOAT *g_trans_z,		
		ProjectorKernel &projector,
		XFLOAT *g_corr_img,
    	XFLOAT *g_diff2s,
		unsigned image_size,
		XFLOAT sum_init,
		XFLOAT exp_local_sqrtXi2,
		unsigned long todo_blocks,
		unsigned long *d_rot_idx,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num
		)
{                               
	int bid = blockIdx_x;

//    // Specialize BlockReduce for a 1D block of 128 threads on type XFLOAT
//    typedef cub::BlockReduce<XFLOAT, 128> BlockReduce;
//    // Allocate shared memory for BlockReduce
//     typename BlockReduce::TempStorage temp_storage;

	int pixel;
	XFLOAT ref_real, ref_imag, shifted_real, shifted_imag;

	unsigned trans_num   = d_job_num[bid]; //how many transes we have for this rot
	
    XFLOAT  s   [trans_num][block_sz]; 
    XFLOAT  s_cc[trans_num][block_sz];
    memset(&s[0][0],    0, sizeof(XFLOAT) * block_sz * trans_num);
    memset(&s_cc[0][0], 0, sizeof(XFLOAT) * block_sz * trans_num);
                 	
	// index of comparison
	unsigned long int ix = d_rot_idx[d_job_idx[bid]];
	unsigned long int iy;
	unsigned pass_num( ceilfracf(image_size,block_sz) );

	for (unsigned pass = 0; pass < pass_num; pass++) { // finish an entire ref image each block
        int start = pass * block_sz;
        #pragma simd
	    for(int tid=0; tid<block_sz; tid++) {
			pixel = start + tid;

			if(pixel >= image_size)
			    continue;
			    
            int x,y,z,xy;
			if(DATA3D)
			{
				z =  floorfracf(pixel, projector.imgX*projector.imgY);
				xy = pixel % (projector.imgX*projector.imgY);
				x =             xy  % projector.imgX;
				y = floorfracf( xy,   projector.imgX);
				if (z > projector.maxR)
				{
					if (z >= projector.imgZ - projector.maxR)
						z = z - projector.imgZ;
					else
						x = projector.maxR;
				}
			}
			else
			{
				x =             pixel % projector.imgX;
				y = floorfracf( pixel , projector.imgX);
			}

			if (y > projector.maxR)
			{
				if (y >= projector.imgY - projector.maxR)
					y = y - projector.imgY;
				else
					x = projector.maxR;
			}

			if(DATA3D)
				projector.project3Dmodel(
					x,y,z,
					g_eulers[ix*9  ], g_eulers[ix*9+1], g_eulers[ix*9+2],
					g_eulers[ix*9+3], g_eulers[ix*9+4], g_eulers[ix*9+5],
					g_eulers[ix*9+6], g_eulers[ix*9+7], g_eulers[ix*9+8],
					ref_real, ref_imag);
			else if(REF3D)
				projector.project3Dmodel(
					x,y,
					g_eulers[ix*9  ], g_eulers[ix*9+1],
					g_eulers[ix*9+3], g_eulers[ix*9+4],
					g_eulers[ix*9+6], g_eulers[ix*9+7],
					ref_real, ref_imag);
			else
				projector.project2Dmodel(
					x,y,
					g_eulers[ix*9  ], g_eulers[ix*9+1],
					g_eulers[ix*9+3], g_eulers[ix*9+4],
					ref_real, ref_imag);

			for (int itrans=0; itrans<trans_num; itrans++) // finish all translations in each partial pass
			{
				iy = d_trans_idx[d_job_idx[bid]] + itrans;

				translatePixel(x, y, g_trans_x[iy], g_trans_y[iy], g_imgs_real[pixel],
				               g_imgs_imag[pixel], shifted_real, shifted_imag);
				s[itrans][tid]    += (ref_real * shifted_real + ref_imag * shifted_imag) *
				                                 g_corr_img[pixel];
				s_cc[itrans][tid] += (ref_real*ref_real + ref_imag*ref_imag) * g_corr_img[pixel];
			}
		} // loop tid
	} // loop pass
	
    for(int itrans=0; itrans<trans_num; itrans++) {
        XFLOAT sum1 = (XFLOAT)0.0;
        XFLOAT sum2 = (XFLOAT)0.0;
        for(int tid=0; tid<block_sz; tid++){
            sum1 += s   [itrans][tid];
            sum2 += s_cc[itrans][tid]; 
        } 
                     
        unsigned long int iy = d_job_idx[bid] + itrans;
#ifdef RELION_SINGLE_PRECISION         
        g_diff2s[iy] = - sum1 / sqrtf(sum2);
#else
        g_diff2s[iy] = - sum1 / sqrt(sum2);
#endif
        
    }				    	
}

} // end of namespace CpuKernels

#endif /* DIFF2_KERNELS_H_ */
