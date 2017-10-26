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

void dump_array(char *name, bool *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%d, ", ptr[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_array(char *name, int *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%d, ", ptr[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_array(char *name, size_t *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%d, ", ptr[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_array(char *name, float *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f, ", ptr[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_complex_array(char *name, ACCCOMPLEX *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f,%f, ", ptr[i].x, ptr[i].y);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_complex_array(char *name, Complex *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f,%f, ", ptr[i].real, ptr[i].imag);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_double_array(char *name, float *ptr, float *ptr2, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f,%f, ", ptr[i], ptr2[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_triple_array(char *name, float *ptr, float *ptr2, float *ptr3, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f,%f,%f, ", ptr[i], ptr2[i], ptr3[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_array(char *name, double *ptr, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f, ", ptr[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_double_array(char *name, double *ptr, double *ptr2, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f,%f, ", ptr[i], ptr2[i]);
			count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

void dump_triple_array(char *name, double *ptr, double *ptr2, double *ptr3, size_t size)
{
	int count = 0;
	FILE *fp = fopen(name, "w");
	fprintf(fp, "Array size:  %ld\n", size);
	for (size_t i=0; i < size; i++) {
		fprintf(fp, "%f,%f,%f, ", ptr[i], ptr2[i], ptr3[i]);
		count++;
		if (count > 10) {
			fprintf(fp, "\n");
			count = 0;
		}
	}
	fprintf(fp, "\n");
	fflush(fp);	
	fclose(fp);
}

namespace AccUtilities
{
	
static void softMaskBackgroundValue(
	AccDataTypes::Image<XFLOAT> &vol,
	bool     do_Mnoise,
	XFLOAT   radius,
	XFLOAT   radius_p,
	XFLOAT   cosine_width,
	AccPtr<XFLOAT> &g_sum,
	AccPtr<XFLOAT> &g_sum_bg)
{
	int block_dim = 128; //TODO: set balanced (hardware-dep?)
#ifdef CUDA
		cuda_kernel_softMaskBackgroundValue<<<block_dim,SOFTMASK_BLOCK_SIZE,0, vol.getStream()>>>(
				~vol,
				vol.getxyz(),
				vol.getx(),
				vol.gety(),
				vol.getz(),
				vol.getx()/2,
				vol.gety()/2,
				vol.getz()/2,
				do_Mnoise,
				radius,
				radius_p,
				cosine_width,
				~g_sum,
				~g_sum_bg);
#else
	CpuKernels::softMaskBackgroundValue(
			block_dim,
			SOFTMASK_BLOCK_SIZE,
			~vol,
			vol.getxyz(),
			vol.getx(),
			vol.gety(),
			vol.getz(),
			vol.getx()/2,
			vol.gety()/2,
			vol.getz()/2,
			do_Mnoise,
			radius,
			radius_p,
			cosine_width,
			~g_sum,
			~g_sum_bg);
#endif
}

static void cosineFilter(
		AccDataTypes::Image<XFLOAT> &vol,
		bool do_Mnoise,
		XFLOAT radius,
		XFLOAT radius_p,
		XFLOAT cosine_width,
		XFLOAT sum_bg_total)
{
	int block_dim = 128; //TODO: set balanced (hardware-dep?)
#ifdef CUDA
	cuda_kernel_cosineFilter<<<block_dim,SOFTMASK_BLOCK_SIZE,0,vol.getStream()>>>(
			~vol,
			vol.getxyz(),
			vol.getx(),
			vol.gety(),
			vol.getz(),
			vol.getx()/2,
			vol.gety()/2,
			vol.getz()/2,
			do_Mnoise,
			radius,
			radius_p,
			cosine_width,
			sum_bg_total);
#else
	CpuKernels::cosineFilter(
			block_dim,
			SOFTMASK_BLOCK_SIZE,
			~vol,
			vol.getxyz(),
			vol.getx(),
			vol.gety(),
			vol.getz(),
			vol.getx()/2,
			vol.gety()/2,
			vol.getz()/2,
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
	CpuKernels::centerFFT_2D<XFLOAT>(batch_size, 0, image_size/2,
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
	CpuKernels::centerFFT_2D<XFLOAT>(batch_size, 0, image_size/2,
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
	CpuKernels::centerFFT_3D<XFLOAT>(batch_size, 0, image_size/2,
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

void kernel_exponentiate_weights_fine(	XFLOAT *g_pdf_orientation,
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
										long int job_num,
										cudaStream_t stream)
{
	long block_num = ceil((double)job_num / (double)SUMW_BLOCK_SIZE);

#ifdef CUDA
	cuda_kernel_exponentiate_weights_fine<<<block_num,SUMW_BLOCK_SIZE,0,stream>>>(
		g_pdf_orientation,
		g_pdf_orientation_zeros,
		g_pdf_offset,
		g_pdf_offset_zeros,
		g_weights,
		min_diff2,
		oversamples_orient,
		oversamples_trans,
		d_rot_id,
		d_trans_idx,
		d_job_idx,
		d_job_num,
		job_num);
#else
	CpuKernels::exponentiate_weights_fine(
		g_pdf_orientation,
		g_pdf_orientation_zeros,
		g_pdf_offset,
		g_pdf_offset_zeros,
		g_weights,
		min_diff2,
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

