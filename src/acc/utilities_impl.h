#ifndef ACC_UTILITIES_IMPL_H_
#define ACC_UTILITIES_IMPL_H_

#include "src/acc/acc_ptr.h"
#include "src/acc/data_types.h"
#ifdef CUDA
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#include "src/acc/cuda/cuda_kernels/wavg.cuh"
#include "src/acc/cuda/cuda_kernels/diff2.cuh"
#else
#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_kernels/wavg.h"
#include "src/acc/cpu/cpu_kernels/diff2.h"
#endif

namespace AccUtilities
{
	
 static void softMaskBackgroundValue(
		int inblock_dim, 
		int inblock_size,
		XFLOAT *vol,
		Image<RFLOAT> &img,
		bool     do_Mnoise,
		XFLOAT   radius,
		XFLOAT   radius_p,
		XFLOAT   cosine_width,
		XFLOAT  *g_sum,
		XFLOAT  *g_sum_bg)
{
#ifdef CUDA
	dim3 block_dim = inblock_dim;
	cuda_kernel_softMaskBackgroundValue<<<block_dim,inblock_size>>>(	vol,
																img().nzyxdim,
																img.data.xdim,
																img.data.ydim,
																img.data.zdim,
																img.data.xdim/2,
																img.data.ydim/2,
																img.data.zdim/2, //unused
																do_Mnoise,
																radius,
																radius_p,
																cosine_width,
																g_sum,
																g_sum_bg);
#else
	CpuKernels::softMaskBackgroundValue(inblock_dim, inblock_size,
		vol,
		img().nzyxdim,
		img.data.xdim,
		img.data.ydim,
		img.data.zdim,
		img.data.xdim/2,
		img.data.ydim/2,
		img.data.zdim/2, //unused
		do_Mnoise,
		radius,
		radius_p,
		cosine_width,
		g_sum,
		g_sum_bg);
#endif
}

static void cosineFilter(
		int inblock_dim, 
		int inblock_size,
		XFLOAT *vol,
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
		XFLOAT cosine_width,
		XFLOAT sum_bg_total)
{
#ifdef CUDA
	dim3 block_dim = inblock_dim;
	cuda_kernel_cosineFilter<<<block_dim,inblock_size>>>(	vol,
															vol_size,
															xdim,
															ydim,
															zdim,
															xinit,
															yinit,
															zinit, //unused
															do_Mnoise,
															radius,
															radius_p,
															cosine_width,
															sum_bg_total);
#else
	CpuKernels::cosineFilter(inblock_dim, inblock_size,
							vol,
							vol_size,
							xdim,
							ydim,
							zdim,
							xinit,
							yinit,
							zinit, //unused
							do_Mnoise,
							radius,
							radius_p,
							cosine_width,
							sum_bg_total);	
#endif
}

void centerFFT_2D(int grid_size, int batch_size, int block_size,
				cudaStream_t stream,
				XFLOAT *img_in,
				int image_size,
				int xdim,
				int ydim,
				int xshift,
				int yshift)
{
#ifdef CUDA
	dim3 blocks(grid_size, batch_size);
	cuda_kernel_centerFFT_2D<<<blocks,block_size,0,stream>>>(
				img_in,
				image_size,
				xdim,
				ydim,
				xshift,
				yshift);
#else
	CpuKernels::centerFFT_2D(grid_size, batch_size, block_size,
				img_in,
				image_size,
				xdim,
				ydim,
				xshift,
				yshift);
#endif	
}

void centerFFT_2D(int grid_size, int batch_size, int block_size,
				XFLOAT *img_in,
				int image_size,
				int xdim,
				int ydim,
				int xshift,
				int yshift)
{
#ifdef CUDA
	dim3 blocks(grid_size, batch_size);
	cuda_kernel_centerFFT_2D<<<blocks,block_size>>>(
				img_in,
				image_size,
				xdim,
				ydim,
				xshift,
				yshift);
#else
	CpuKernels::centerFFT_2D(grid_size, batch_size, block_size,
				img_in,
				image_size,
				xdim,
				ydim,
				xshift,
				yshift);
#endif
}

void centerFFT_3D(int grid_size, int batch_size, int block_size,
				cudaStream_t stream,
				XFLOAT *img_in,
				int image_size,
				int xdim,
				int ydim,
				int zdim,
				int xshift,
				int yshift,
				int zshift)
{
#ifdef CUDA
	dim3 blocks(grid_size, batch_size);
	cuda_kernel_centerFFT_3D<<<blocks,block_size, 0, stream>>>(
				img_in,
				image_size,
				xdim,
				ydim,
				zdim,
				xshift,
				yshift,
				zshift);
#else
	CpuKernels::centerFFT_3D(grid_size, batch_size, block_size,
				img_in,
				image_size,
				xdim,
				ydim,
				zdim,
				xshift,
				yshift,
				zshift);
#endif
}

void kernel_exponentiate_weights_fine(	int grid_size, 
										int block_size,
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
										long int job_num,
										cudaStream_t stream)
{
#ifdef CUDA
	dim3 block_dim(grid_size);
	cuda_kernel_exponentiate_weights_fine<<<block_dim,block_size,0,stream>>>(
		g_pdf_orientation,
		g_pdf_offset,
		g_weights,
		avg_diff2,
		oversamples_orient,
		oversamples_trans,
		d_rot_id,
		d_trans_idx,
		d_job_idx,
		d_job_num,
		job_num);
#else
	CpuKernels::exponentiate_weights_fine(
		grid_size,
		block_size,
		g_pdf_orientation,
		g_pdf_offset,
		g_weights,
		avg_diff2,
		oversamples_orient,
		oversamples_trans,
		d_rot_id,
		d_trans_idx,
		d_job_idx,
		d_job_num,
		job_num);
#endif
}

};  // namespace AccUtilities


#endif //ACC_UTILITIES_H_

