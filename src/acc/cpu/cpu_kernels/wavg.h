#ifndef WAVG_KERNEL_H_
#define WAVG_KERNEL_H_

#include <vector>
#include <iostream>
#include <fstream>

#include "src/acc/cpu/cpu_settings.h"
#include "src/acc/cpu/cpu_projector.h"
#include "src/acc/cpu/utilities.h"
#include "src/acc/cpu/cpu_kernels/helper.h"

#define RESTRICT

namespace CpuKernels
{

template<bool REFCTF, bool REF3D, bool DATA3D, int block_sz>
void wavg(
        int     blockIdx_x,
		XFLOAT * RESTRICT g_eulers,
		ProjectorKernel &projector,
		unsigned image_size,
		unsigned long orientation_num,
		XFLOAT * RESTRICT g_img_real,
		XFLOAT * RESTRICT g_img_imag,
		XFLOAT * RESTRICT g_trans_x,
		XFLOAT * RESTRICT g_trans_y,
		XFLOAT * RESTRICT g_trans_z,
		XFLOAT * RESTRICT g_weights,
		XFLOAT * RESTRICT g_ctfs,
		XFLOAT * RESTRICT g_wdiff2s_parts,
		XFLOAT * RESTRICT g_wdiff2s_AA,
		XFLOAT * RESTRICT g_wdiff2s_XA,
		unsigned long translation_num,
		XFLOAT weight_norm,
		XFLOAT significant_weight,
		XFLOAT part_scale)
{
	XFLOAT ref_real, ref_imag, img_real, img_imag, trans_real, trans_imag;

	int bid = blockIdx_x; //block ID

	unsigned pass_num(ceilfracf(image_size, block_sz)),pixel;
    XFLOAT weight_norm_inverse = (XFLOAT) 1.0 / weight_norm;
    XFLOAT s_eulers[9];

    for (int i = 0; i < 9; i++)
        s_eulers[i] = g_eulers[bid*9+i];

    XFLOAT trans_x[translation_num], trans_y[translation_num], 
           trans_z[translation_num];
    
    for (unsigned long itrans = 0; itrans < translation_num; itrans++)
    {
        trans_x[itrans] = g_trans_x[itrans];
        trans_y[itrans] = g_trans_y[itrans];
        if(DATA3D)
            trans_z[itrans] = g_trans_z[itrans];
    }

    for (unsigned pass = 0; pass < pass_num; pass++) {// finish a reference proj in each block
        #pragma simd
        for (int tid=0; tid<block_sz; tid++) {
            XFLOAT ss_wdiff2s_parts = (XFLOAT)0.0;
            XFLOAT ss_sumXA = (XFLOAT)0.0;
            XFLOAT ss_sumA2 = (XFLOAT)0.0;

	    	pixel = pass * block_sz + tid;

            if(pixel<image_size)
            {
		    	int x,y,z,xy;
			    if(DATA3D)
			    {
				    z  =  floorfracf(pixel, projector.imgX*projector.imgY);
				    xy = pixel % (projector.imgX*projector.imgY);
			  	    x  =             xy  % projector.imgX;
				    y  = floorfracf( xy,   projector.imgX);
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
					    s_eulers[0], s_eulers[1], s_eulers[2],
					    s_eulers[3], s_eulers[4], s_eulers[5],
					    s_eulers[6], s_eulers[7], s_eulers[8],
					    ref_real, ref_imag);
			    else if(REF3D)
                    projector.project3Dmodel(
                        x,y,
                        s_eulers[0], s_eulers[1],
                        s_eulers[3], s_eulers[4],
                        s_eulers[6], s_eulers[7],
                        ref_real, ref_imag);
                else
                    projector.project2Dmodel(
                            x,y,
                        s_eulers[0], s_eulers[1],
                        s_eulers[3], s_eulers[4],
                        ref_real, ref_imag);

                if (REFCTF)
                {
                    ref_real *= g_ctfs[pixel];
                    ref_imag *= g_ctfs[pixel];
                }
                else
                {
                    ref_real *= part_scale;
                    ref_imag *= part_scale;
                }

                img_real = g_img_real[pixel];
                img_imag = g_img_imag[pixel];

                for (unsigned long itrans = 0; itrans < translation_num; itrans++)
                {
                    XFLOAT weight = g_weights[bid * translation_num + itrans];

                    if (weight >= significant_weight)
                    {
                        weight *= weight_norm_inverse;

				    	if(DATA3D)
						    translatePixel(x, y, z, trans_x[itrans], trans_y[itrans], trans_z[itrans], img_real, img_imag, trans_real, trans_imag);
					    else
						    translatePixel(x, y,    trans_x[itrans], trans_y[itrans],                    img_real, img_imag, trans_real, trans_imag);

				    	XFLOAT diff_real = ref_real - trans_real;
					    XFLOAT diff_imag = ref_imag - trans_imag;

                        ss_wdiff2s_parts += weight * (diff_real*diff_real + diff_imag*diff_imag);
                        ss_sumXA +=  weight * ( ref_real * trans_real + ref_imag * trans_imag);
                        ss_sumA2 +=  weight * ( ref_real*ref_real  +  ref_imag*ref_imag );
                    }   // (weight >= significant_weight)
                }  // for itrans

                // No need to add since pixel is unique:  pixel = pass * block_sz + tid;
				// TODO - atomics needed if threaded.
                g_wdiff2s_XA[pixel] += ss_sumXA;
                g_wdiff2s_AA[pixel] += ss_sumA2;
                g_wdiff2s_parts[pixel] += ss_wdiff2s_parts;
            }  // pixel<image_size
        } // for tid
    } // pass
}

} // end of namespace CpuKernels

#endif /* WAVG_KERNEL_H_ */
