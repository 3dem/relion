#include "src/gpu_utils/cuda_backprojector.h"
#include "src/gpu_utils/cuda_device_utils.cuh"
#include <signal.h>

#define BACKPROJECTION4_BLOCK_SIZE 64
#define BACKPROJECTION4_GROUP_SIZE 16
#define BACKPROJECTION4_PREFETCH_COUNT 3
#define BP_2D_BLOCK_SIZE 128

void CudaBackprojector::setMdlDim(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int max_r, int paddingFactor)
{
	if (xdim != mdlX ||
		ydim != mdlY ||
		zdim != mdlZ ||
		inity != mdlInitY ||
		initz != mdlInitZ ||
		max_r != maxR ||
		paddingFactor != padding_factor)
	{
		mdlX = xdim;
		mdlY = ydim;
		mdlZ = zdim;
		if (mdlZ < 1) mdlZ = 1;
		mdlXYZ = xdim*ydim*zdim;
		mdlInitY = inity;
		mdlInitZ = initz;
		maxR = max_r;
		maxR2 = max_r*max_r;
		padding_factor = paddingFactor;

		clear();

		//Allocate space for model
		HANDLE_ERROR(cudaMalloc( (void**) &d_mdlReal,   mdlXYZ * sizeof(XFLOAT)));
		HANDLE_ERROR(cudaMalloc( (void**) &d_mdlImag,   mdlXYZ * sizeof(XFLOAT)));
		HANDLE_ERROR(cudaMalloc( (void**) &d_mdlWeight, mdlXYZ * sizeof(XFLOAT)));
	}
}

void CudaBackprojector::initMdl()
{
#ifdef CUDA_DEBUG
	if (mdlXYZ == 0)
	{
        printf("Model dimensions must be set with setMdlDim before call to setupMdl.");
		raise(SIGSEGV);
	}
	if (voxelCount != 0)
	{
        printf("DEBUG_ERROR: Duplicated call to model setup");
		raise(SIGSEGV);
	}
#endif

	//Initiate model with zeros
	DEBUG_HANDLE_ERROR(cudaMemset( d_mdlReal,   0, mdlXYZ * sizeof(XFLOAT)));
	DEBUG_HANDLE_ERROR(cudaMemset( d_mdlImag,   0, mdlXYZ * sizeof(XFLOAT)));
	DEBUG_HANDLE_ERROR(cudaMemset( d_mdlWeight, 0, mdlXYZ * sizeof(XFLOAT)));
}

__global__ void cuda_kernel_backproject2D(
		XFLOAT *g_model_real,
		XFLOAT *g_model_imag,
		XFLOAT *g_model_weight,
		XFLOAT *g_wavgs_real,
		XFLOAT *g_wavgs_imag,
		XFLOAT *g_Fweights,
		XFLOAT *g_eulers,
		int max_r,
		int max_r2,
		XFLOAT padding_factor,
		unsigned img_x,
		unsigned img_y,
		unsigned img_xy,
		unsigned mdl_x,
		int mdl_inity)
{
	unsigned tid = threadIdx.x;
	unsigned img = blockIdx.x;

	__shared__ XFLOAT s_eulers[4];

	if (tid == 0)
		s_eulers[0] = g_eulers[img*9+0] * padding_factor;
	else if (tid == 1)
		s_eulers[1] = g_eulers[img*9+1] * padding_factor;
	else if (tid == 2)
		s_eulers[2] = g_eulers[img*9+3] * padding_factor;
	else if (tid == 3)
		s_eulers[3] = g_eulers[img*9+4] * padding_factor;

	__syncthreads();

	int pixel_pass_num(ceilf((float)img_xy/(float)BP_2D_BLOCK_SIZE));

	for (unsigned pass = 0; pass < pixel_pass_num; pass++)
    {
		unsigned pixel = (pass * BP_2D_BLOCK_SIZE) + tid;

		if (pixel >= img_xy)
			continue;

		int x = pixel % img_x;
		int y = (int)floorf( (float)pixel / (float)img_x);

		pixel += img * img_xy;

		// Don't search beyond square with side max_r
		if (y > max_r)
		{
			if (y >= img_y - max_r)
				y -= img_y;
			else
				continue;
		}

		if (x * x + y * y > max_r2)
			continue;

		// Get the weight
		XFLOAT weight = g_Fweights[pixel];

		if (weight > 0.f)
		{
			// Get the relevant value in the input image
			XFLOAT real = g_wavgs_real[pixel];
			XFLOAT imag = g_wavgs_imag[pixel];

			// Get logical coordinates in the 3D map
			XFLOAT xp = (s_eulers[0] * x + s_eulers[1] * y );
			XFLOAT yp = (s_eulers[2] * x + s_eulers[3] * y );

			// Only asymmetric half is stored
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				imag = -imag;
			}

			int x0 = floorf(xp);
			XFLOAT fx = xp - x0;
			int x1 = x0 + 1;

			int y0 = floorf(yp);
			XFLOAT fy = yp - y0;
			y0 -= mdl_inity;
			int y1 = y0 + 1;

			XFLOAT mfx = 1.f - fx;
			XFLOAT mfy = 1.f - fy;

			XFLOAT dd00 = mfy * mfx;
			XFLOAT dd01 = mfy *  fx;
			XFLOAT dd10 =  fy * mfx;
			XFLOAT dd11 =  fy *  fx;

			cuda_atomic_add(&g_model_real  [y0 * mdl_x + x0], dd00 * real);
			cuda_atomic_add(&g_model_imag  [y0 * mdl_x + x0], dd00 * imag);
			cuda_atomic_add(&g_model_weight[y0 * mdl_x + x0], dd00 * weight);

			cuda_atomic_add(&g_model_real  [y0 * mdl_x + x1], dd01 * real);
			cuda_atomic_add(&g_model_imag  [y0 * mdl_x + x1], dd01 * imag);
			cuda_atomic_add(&g_model_weight[y0 * mdl_x + x1], dd01 * weight);

			cuda_atomic_add(&g_model_real  [y1 * mdl_x + x0], dd10 * real);
			cuda_atomic_add(&g_model_imag  [y1 * mdl_x + x0], dd10 * imag);
			cuda_atomic_add(&g_model_weight[y1 * mdl_x + x0], dd10 * weight);

			cuda_atomic_add(&g_model_real  [y1 * mdl_x + x1], dd11 * real);
			cuda_atomic_add(&g_model_imag  [y1 * mdl_x + x1], dd11 * imag);
			cuda_atomic_add(&g_model_weight[y1 * mdl_x + x1], dd11 * weight);
		}
	}
}

__global__ void cuda_kernel_backproject3D_scatter(
		XFLOAT *g_model_real,
		XFLOAT *g_model_imag,
		XFLOAT *g_model_weight,
		XFLOAT *g_wavgs_real,
		XFLOAT *g_wavgs_imag,
		XFLOAT *g_Fweights,
		XFLOAT *g_eulers,
		int max_r,
		int max_r2,
		XFLOAT padding_factor,
		unsigned img_x,
		unsigned img_y,
		unsigned img_xy,
		unsigned mdl_x,
		unsigned mdl_y,
		int mdl_inity,
		int mdl_initz)
{
	unsigned tid = threadIdx.x;
	unsigned img = blockIdx.x;

	__shared__ XFLOAT s_eulers[9];

	if (tid < 9)
		s_eulers[tid] = g_eulers[img*9+tid];

	__syncthreads();

	int pixel_pass_num(ceilf((float)img_xy/(float)BP_2D_BLOCK_SIZE));
	for (unsigned pass = 0; pass < pixel_pass_num; pass++)
    {
		unsigned pixel = (pass * BP_2D_BLOCK_SIZE) + tid;

		if (pixel >= img_xy)
			continue;

		int x = pixel % img_x;
		int y = (int)floorf( (float)pixel / (float)img_x);

		pixel += img * img_xy;

		// Don't search beyond square with side max_r
		if (y > max_r)
		{
			if (y >= img_y - max_r)
				y -= img_y;
		}

		if (x * x + y * y > max_r2)
			continue;

		// Get the weight
		XFLOAT Fweights = g_Fweights[pixel];

		if (Fweights > 0.f)
		{
			// Get the relevant value in the input image
			XFLOAT real = g_wavgs_real[pixel];
			XFLOAT imag = g_wavgs_imag[pixel];

			// Get logical coordinates in the 3D map
			XFLOAT xp = (s_eulers[0] * x + s_eulers[1] * y ) * padding_factor;
			XFLOAT yp = (s_eulers[3] * x + s_eulers[4] * y ) * padding_factor;
			XFLOAT zp = (s_eulers[6] * x + s_eulers[7] * y ) * padding_factor;

			// Only asymmetric half is stored
			if (xp < 0.f)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;
				imag = -imag;
			}

			int x0 = floorf(xp);
			XFLOAT fx = xp - x0;
			int x1 = x0 + 1;

			int y0 = floorf(yp);
			XFLOAT fy = yp - y0;
			y0 -= mdl_inity;
			int y1 = y0 + 1;

			int z0 = floorf(zp);
			XFLOAT fz = zp - z0;
			z0 -= mdl_initz;
			int z1 = z0 + 1;

			XFLOAT mfx = 1.f - fx;
			XFLOAT mfy = 1.f - fy;
			XFLOAT mfz = 1.f - fz;

			XFLOAT dd000 = mfz * mfy * mfx;

			cuda_atomic_add(&g_model_real  [z0 * mdl_x * mdl_y + y0 * mdl_x + x0], dd000 * real);
			cuda_atomic_add(&g_model_imag  [z0 * mdl_x * mdl_y + y0 * mdl_x + x0], dd000 * imag);
			cuda_atomic_add(&g_model_weight[z0 * mdl_x * mdl_y + y0 * mdl_x + x0], dd000 * Fweights);

			XFLOAT dd001 = mfz * mfy *  fx;

			cuda_atomic_add(&g_model_real  [z0 * mdl_x * mdl_y + y0 * mdl_x + x1], dd001 * real);
			cuda_atomic_add(&g_model_imag  [z0 * mdl_x * mdl_y + y0 * mdl_x + x1], dd001 * imag);
			cuda_atomic_add(&g_model_weight[z0 * mdl_x * mdl_y + y0 * mdl_x + x1], dd001 * Fweights);

			XFLOAT dd010 = mfz *  fy * mfx;

			cuda_atomic_add(&g_model_real  [z0 * mdl_x * mdl_y + y1 * mdl_x + x0], dd010 * real);
			cuda_atomic_add(&g_model_imag  [z0 * mdl_x * mdl_y + y1 * mdl_x + x0], dd010 * imag);
			cuda_atomic_add(&g_model_weight[z0 * mdl_x * mdl_y + y1 * mdl_x + x0], dd010 * Fweights);

			XFLOAT dd011 = mfz *  fy *  fx;

			cuda_atomic_add(&g_model_real  [z0 * mdl_x * mdl_y + y1 * mdl_x + x1], dd011 * real);
			cuda_atomic_add(&g_model_imag  [z0 * mdl_x * mdl_y + y1 * mdl_x + x1], dd011 * imag);
			cuda_atomic_add(&g_model_weight[z0 * mdl_x * mdl_y + y1 * mdl_x + x1], dd011 * Fweights);

			XFLOAT dd100 =  fz * mfy * mfx;

			cuda_atomic_add(&g_model_real  [z1 * mdl_x * mdl_y + y0 * mdl_x + x0], dd100 * real);
			cuda_atomic_add(&g_model_imag  [z1 * mdl_x * mdl_y + y0 * mdl_x + x0], dd100 * imag);
			cuda_atomic_add(&g_model_weight[z1 * mdl_x * mdl_y + y0 * mdl_x + x0], dd100 * Fweights);

			XFLOAT dd101 =  fz * mfy *  fx;

			cuda_atomic_add(&g_model_real  [z1 * mdl_x * mdl_y + y0 * mdl_x + x1], dd101 * real);
			cuda_atomic_add(&g_model_imag  [z1 * mdl_x * mdl_y + y0 * mdl_x + x1], dd101 * imag);
			cuda_atomic_add(&g_model_weight[z1 * mdl_x * mdl_y + y0 * mdl_x + x1], dd101 * Fweights);

			XFLOAT dd110 =  fz *  fy * mfx;

			cuda_atomic_add(&g_model_real  [z1 * mdl_x * mdl_y + y1 * mdl_x + x0], dd110 * real);
			cuda_atomic_add(&g_model_imag  [z1 * mdl_x * mdl_y + y1 * mdl_x + x0], dd110 * imag);
			cuda_atomic_add(&g_model_weight[z1 * mdl_x * mdl_y + y1 * mdl_x + x0], dd110 * Fweights);

			XFLOAT dd111 =  fz *  fy *  fx;

			cuda_atomic_add(&g_model_real  [z1 * mdl_x * mdl_y + y1 * mdl_x + x1], dd111 * real);
			cuda_atomic_add(&g_model_imag  [z1 * mdl_x * mdl_y + y1 * mdl_x + x1], dd111 * imag);
			cuda_atomic_add(&g_model_weight[z1 * mdl_x * mdl_y + y1 * mdl_x + x1], dd111 * Fweights);

		}
	}
}


void CudaBackprojector::backproject(
		XFLOAT *d_real,
		XFLOAT *d_imag,
		XFLOAT *d_weight,
		XFLOAT *d_eulers,
		int imgX,
		int imgY,
		unsigned long imageCount)
{

	if(mdlZ==1)
	{
		cuda_kernel_backproject2D<<<imageCount,BP_2D_BLOCK_SIZE,0,stream>>>(
			d_mdlReal,
			d_mdlImag,
			d_mdlWeight,
			d_real,
			d_imag,
			d_weight,
			d_eulers,
			maxR,
			maxR2,
			padding_factor,
			imgX,
			imgY,
			imgX*imgY,
			mdlX,
			mdlInitY);
	}
	else
	{
		cuda_kernel_backproject3D_scatter<<<imageCount,BP_2D_BLOCK_SIZE,0,stream>>>(
				d_mdlReal,
				d_mdlImag,
				d_mdlWeight,
				d_real,
				d_imag,
				d_weight,
				d_eulers,
				maxR,
				maxR2,
				padding_factor,
				imgX,
				imgY,
				imgX*imgY,
				mdlX,
				mdlY,
				mdlInitY,
				mdlInitZ);
	}
}


void CudaBackprojector::getMdlData(XFLOAT *r, XFLOAT *i, XFLOAT * w)
{
	DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream)); //Make sure to wait for remaining kernel executions

	DEBUG_HANDLE_ERROR(cudaMemcpyAsync( r, d_mdlReal,   mdlXYZ * sizeof(XFLOAT), cudaMemcpyDeviceToHost, stream));
	DEBUG_HANDLE_ERROR(cudaMemcpyAsync( i, d_mdlImag,   mdlXYZ * sizeof(XFLOAT), cudaMemcpyDeviceToHost, stream));
	DEBUG_HANDLE_ERROR(cudaMemcpyAsync( w, d_mdlWeight, mdlXYZ * sizeof(XFLOAT), cudaMemcpyDeviceToHost, stream));

	DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream)); //Wait for copy
}

void CudaBackprojector::clear()
{
	if (d_mdlReal != NULL)
	{
		DEBUG_HANDLE_ERROR(cudaFree(d_mdlReal));
		DEBUG_HANDLE_ERROR(cudaFree(d_mdlImag));
		DEBUG_HANDLE_ERROR(cudaFree(d_mdlWeight));

		d_mdlReal = d_mdlImag = d_mdlWeight = NULL;
	}
}

CudaBackprojector::~CudaBackprojector()
{
	clear();
}
