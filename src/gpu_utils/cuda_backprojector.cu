#include "src/gpu_utils/cuda_backprojector.h"
#include "src/gpu_utils/cuda_device_utils.cuh"
#include <signal.h>

#define BACKPROJECTION4_BLOCK_SIZE 64
#define BACKPROJECTION4_GROUP_SIZE 16
#define BACKPROJECTION4_PREFETCH_COUNT 3
#define BP_2D_BLOCK_SIZE 128


void CudaBackprojector::initMdl(int streamPriority)
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

	if(mdlZ==1) //2D case
	{
		//Allocate space for model
		HANDLE_ERROR(cudaMalloc( (void**) &d_mdlReal, mdlXYZ * sizeof(XFLOAT)));
		HANDLE_ERROR(cudaMalloc( (void**) &d_mdlImag, mdlXYZ * sizeof(XFLOAT)));
		HANDLE_ERROR(cudaMalloc( (void**) &d_mdlWeight, mdlXYZ * sizeof(XFLOAT)));

		//Initiate model with zeros
		HANDLE_ERROR(cudaMemset( d_mdlReal, 0, mdlXYZ * sizeof(XFLOAT)));
		HANDLE_ERROR(cudaMemset( d_mdlImag, 0, mdlXYZ * sizeof(XFLOAT)));
		HANDLE_ERROR(cudaMemset( d_mdlWeight, 0, mdlXYZ * sizeof(XFLOAT)));
	}
	else //3D case
	{
		// Generate explicit indicies of voxels

		h_voxelX = new int[mdlXYZ];
		h_voxelY = new int[mdlXYZ];
		h_voxelZ = new int[mdlXYZ];

		for (int x = 0; x < mdlX; x ++)
		{
			for (int y = mdlInitY; y < mdlY + mdlInitY; y++)
			{
				for (int z = mdlInitZ; z < mdlZ + mdlInitZ; z++)
				{
					if (x*x + y*y + z*z <= maxR2 * padding_factor * padding_factor * 1.2f)
					{
						h_voxelX[voxelCount] = x;
						h_voxelY[voxelCount] = y;
						h_voxelZ[voxelCount] = z;
						voxelCount ++;
					}
				}
			}
		}

		//Allocate space for indices
		HANDLE_ERROR(cudaMalloc( (void**) &d_voxelX, voxelCount * sizeof(XFLOAT)));
		HANDLE_ERROR(cudaMalloc( (void**) &d_voxelY, voxelCount * sizeof(XFLOAT)));
		HANDLE_ERROR(cudaMalloc( (void**) &d_voxelZ, voxelCount * sizeof(XFLOAT)));

		//Send over indices
		HANDLE_ERROR(cudaMemcpy( d_voxelX, h_voxelX, voxelCount * sizeof(XFLOAT), cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy( d_voxelY, h_voxelY, voxelCount * sizeof(XFLOAT), cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy( d_voxelZ, h_voxelZ, voxelCount * sizeof(XFLOAT), cudaMemcpyHostToDevice));

		//Allocate space for model
		HANDLE_ERROR(cudaMalloc( (void**) &d_mdlReal, mdlXYZ * sizeof(XFLOAT)));
		HANDLE_ERROR(cudaMalloc( (void**) &d_mdlImag, mdlXYZ * sizeof(XFLOAT)));
		HANDLE_ERROR(cudaMalloc( (void**) &d_mdlWeight, mdlXYZ * sizeof(XFLOAT)));

		//Initiate model with zeros
		HANDLE_ERROR(cudaMemset( d_mdlReal, 0, mdlXYZ * sizeof(XFLOAT)));
		HANDLE_ERROR(cudaMemset( d_mdlImag, 0, mdlXYZ * sizeof(XFLOAT)));
		HANDLE_ERROR(cudaMemset( d_mdlWeight, 0, mdlXYZ * sizeof(XFLOAT)));
	}

	HANDLE_ERROR(cudaStreamCreateWithPriority(&stream, cudaStreamNonBlocking, streamPriority));
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
		s_eulers[1] = g_eulers[img*9+3] * padding_factor;
	else if (tid == 2)
		s_eulers[2] = g_eulers[img*9+1] * padding_factor;
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

__global__ void cuda_kernel_backproject3D_gather(
		int *g_xs,
		int *g_ys,
		int *g_zs,
		XFLOAT *g_model_real,
		XFLOAT *g_model_imag,
		XFLOAT *g_weight,
		XFLOAT *g_wavgs_real,
		XFLOAT *g_wavgs_imag,
		XFLOAT *g_Fweights,
		XFLOAT *g_eulers,
		int max_r2,
		XFLOAT scale2,
		unsigned img_x,
		unsigned img_y,
		unsigned img_xy,
		unsigned long img_count,
		unsigned mdl_x,
		unsigned mdl_y,
		int mdl_inity,
		int mdl_initz,
		int N)
{
	unsigned tid = threadIdx.x;
	unsigned gid = tid / 4; //Group id
	unsigned mid = tid % 4; //Member id
	unsigned local_prefetch_idx = tid * BACKPROJECTION4_PREFETCH_COUNT; //Starting index of the matrix fetch
	unsigned global_idx = blockIdx.x * BACKPROJECTION4_GROUP_SIZE + gid;

	int X(0),Y(0),Z(0);

	if (global_idx < N)
	{
		X = g_xs[global_idx];
		Y = g_ys[global_idx];
		Z = g_zs[global_idx];
	}
	else
		X = mdl_x * 10; // Padding coordinate, place outside images

	int ax(0), ay(0);

	if (mid == 1)
		ax = 1;
	else if (mid == 2)
		ay = 1;
	else if (mid == 3)
	{
		ax = 1;
		ay = 1;
	}

	bool  is_neg_x;
	XFLOAT d, w;
	XFLOAT xp,yp,zp;
	int x,y,idx;

	__shared__ XFLOAT s_e[BACKPROJECTION4_BLOCK_SIZE*BACKPROJECTION4_PREFETCH_COUNT];

	__shared__ XFLOAT s_weight[BACKPROJECTION4_GROUP_SIZE*4];
	__shared__ XFLOAT s_value_real[BACKPROJECTION4_GROUP_SIZE*4];
	__shared__ XFLOAT s_value_imag[BACKPROJECTION4_GROUP_SIZE*4];

	s_weight[tid] = 0.0f;
	s_value_real[tid] = 0.0f;
	s_value_imag[tid] = 0.0f;

	for (int img = 0, b = BACKPROJECTION4_BLOCK_SIZE*BACKPROJECTION4_PREFETCH_COUNT; img < img_count; img ++, b += 9)
	{
		if (b+9 > BACKPROJECTION4_BLOCK_SIZE*BACKPROJECTION4_PREFETCH_COUNT)
		{
			__syncthreads();

			int global_prefetch_idx = img * 9 + local_prefetch_idx;
			if (global_prefetch_idx < img_count*9)
			{
				s_e[local_prefetch_idx+0] = __ldg(&g_eulers[global_prefetch_idx+0]);
				s_e[local_prefetch_idx+1] = __ldg(&g_eulers[global_prefetch_idx+1]);
				s_e[local_prefetch_idx+2] = __ldg(&g_eulers[global_prefetch_idx+2]);
			}

			__syncthreads();
			b = 0;
		}

		zp = (s_e[b+6] * X + s_e[b+7] * Y + s_e[b+8] * Z) / scale2;

		if (fabsf(zp) > 0.87f) continue; //Within the unit cube, sqrt(3)/2=0.866

		yp = (s_e[b+3] * X + s_e[b+4] * Y + s_e[b+5] * Z) / scale2;
		xp = (s_e[b+0] * X + s_e[b+1] * Y + s_e[b+2] * Z) / scale2;

		if (xp < 0.0f)
		{
			yp = -yp;
			xp = -xp;
			is_neg_x = true;
		}
		else
			is_neg_x = false;

		x = (int) floorf(xp) + ax;
		y = (int) floorf(yp) + ay;

		if (x * x + y * y > max_r2) continue;

		if (y < 0 && x == 0)
		{
			is_neg_x = !is_neg_x;
			y = -y;
		}

		xp = (s_e[b+0] * x + s_e[b+3] * y) * scale2;
		yp = (s_e[b+1] * x + s_e[b+4] * y) * scale2;
		zp = (s_e[b+2] * x + s_e[b+5] * y) * scale2;

		if (xp < 0.0f) //Flip sign
		{
			xp = fabsf(X+xp);
			yp = fabsf(Y+yp);
			zp = fabsf(Z+zp);
		}
		else
		{
			xp = fabsf(X-xp);
			yp = fabsf(Y-yp);
			zp = fabsf(Z-zp);
		}

		if (xp < 1.0f && yp < 1.0f && zp < 1.0f)
		{
			if (y < 0) y += img_y;
			idx = img*img_xy + y * img_x + x;
			w = g_Fweights[idx];

			if (w > 0.0f)
			{
				d = (1.0f - xp) * (1.0f - yp) * (1.0f - zp);

				s_weight[tid] += w * d;
				s_value_real[tid] += g_wavgs_real[idx] * d;
				if (is_neg_x) s_value_imag[tid] -= g_wavgs_imag[idx] * d;
				else          s_value_imag[tid] += g_wavgs_imag[idx] * d;
			}
		}
	}

	__syncthreads();

	if (mid == 0)
	{
		XFLOAT sum = s_weight[gid*4 + 0] + s_weight[gid*4 + 1] + s_weight[gid*4 + 2] + s_weight[gid*4 + 3];
		if (sum != 0.0f)
			g_weight[(Z-mdl_initz)*mdl_x*mdl_y + (Y-mdl_inity)*mdl_x + X] += sum;
	}
	else if (mid == 1)
	{
		XFLOAT sum = s_value_real[gid*4 + 0] + s_value_real[gid*4 + 1] + s_value_real[gid*4 + 2] + s_value_real[gid*4 + 3];
		if (sum != 0.0f)
			g_model_real[(Z-mdl_initz)*mdl_x*mdl_y + (Y-mdl_inity)*mdl_x + X] += sum;
	}
	else if (mid == 2)
	{
		XFLOAT sum = s_value_imag[gid*4 + 0] + s_value_imag[gid*4 + 1] + s_value_imag[gid*4 + 2] + s_value_imag[gid*4 + 3];
		if (sum != 0.0f)
			g_model_imag[(Z-mdl_initz)*mdl_x*mdl_y + (Y-mdl_inity)*mdl_x + X] += sum;
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
			XFLOAT xp = (s_eulers[0] * x + s_eulers[3] * y ) * padding_factor;
			XFLOAT yp = (s_eulers[1] * x + s_eulers[4] * y ) * padding_factor;
			XFLOAT zp = (s_eulers[2] * x + s_eulers[5] * y ) * padding_factor;

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
//		int grid_dim = ceil((float)voxelCount / BACKPROJECTION4_GROUP_SIZE);
//		dim3 block_dim( BACKPROJECTION4_GROUP_SIZE * 4 );
//
//		cuda_kernel_backproject3D_gather<<<grid_dim,block_dim,0,stream>>>(
//			d_voxelX,
//			d_voxelY,
//			d_voxelZ,
//			d_mdlReal,
//			d_mdlImag,
//			d_mdlWeight,
//			d_real,
//			d_imag,
//			d_weight,
//			d_eulers,
//			maxR2,
//			padding_factor,
//			imgX,
//			imgY,
//			imgX*imgY,
//			imageCount,
//			mdlX,
//			mdlY,
//			mdlInitY,
//			mdlInitZ,
//			voxelCount);

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
	HANDLE_ERROR(cudaStreamSynchronize(stream)); //Make sure to wait for remaining kernel executions

	HANDLE_ERROR(cudaMemcpyAsync( r, d_mdlReal,   mdlXYZ * sizeof(XFLOAT), cudaMemcpyDeviceToHost, stream));
	HANDLE_ERROR(cudaMemcpyAsync( i, d_mdlImag,   mdlXYZ * sizeof(XFLOAT), cudaMemcpyDeviceToHost, stream));
	HANDLE_ERROR(cudaMemcpyAsync( w, d_mdlWeight, mdlXYZ * sizeof(XFLOAT), cudaMemcpyDeviceToHost, stream));

	HANDLE_ERROR(cudaStreamSynchronize(stream)); //Wait for copy
}


void CudaBackprojector::getMdlData(Complex *data, double * weights)
{
	XFLOAT *r = new XFLOAT[mdlXYZ];
	XFLOAT *i = new XFLOAT[mdlXYZ];
	XFLOAT *w = new XFLOAT[mdlXYZ];

	getMdlData(r, i, w);

	for (unsigned long n = 0; n < mdlXYZ; n++)
	{
		data[n].real = (double) r[n];
		data[n].imag = (double) i[n];
		weights[n] = (double) w[n];
	}

	delete [] r;
	delete [] i;
	delete [] w;
}


CudaBackprojector::~CudaBackprojector()
{
	if (voxelCount != 0)
	{
		HANDLE_ERROR(cudaFree(d_voxelX));
		HANDLE_ERROR(cudaFree(d_voxelY));
		HANDLE_ERROR(cudaFree(d_voxelZ));

		d_voxelX = d_voxelY = d_voxelZ = 0;

		delete [] h_voxelX;
		delete [] h_voxelY;
		delete [] h_voxelZ;

		h_voxelX = h_voxelY = h_voxelZ = 0;

		HANDLE_ERROR(cudaFree(d_mdlReal));
		HANDLE_ERROR(cudaFree(d_mdlImag));
		HANDLE_ERROR(cudaFree(d_mdlWeight));

		d_mdlReal = d_mdlImag = d_mdlWeight = 0;

		HANDLE_ERROR(cudaStreamDestroy(stream));

		stream = 0;
	}
}
