#include "src/gpu_utils/cuda_backprojector.h"
#include "src/gpu_utils/cuda_utils.cuh"
#include <signal.h>


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

	HANDLE_ERROR(cudaMalloc( (void**) &d_voxelX, voxelCount * sizeof(XFLOAT)));
	HANDLE_ERROR(cudaMalloc( (void**) &d_voxelY, voxelCount * sizeof(XFLOAT)));
	HANDLE_ERROR(cudaMalloc( (void**) &d_voxelZ, voxelCount * sizeof(XFLOAT)));

	HANDLE_ERROR(cudaMemcpy( d_voxelX, h_voxelX, voxelCount * sizeof(XFLOAT), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy( d_voxelY, h_voxelY, voxelCount * sizeof(XFLOAT), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy( d_voxelZ, h_voxelZ, voxelCount * sizeof(XFLOAT), cudaMemcpyHostToDevice));

	HANDLE_ERROR(cudaMalloc( (void**) &d_mdlReal, voxelCount * sizeof(XFLOAT)));
	HANDLE_ERROR(cudaMalloc( (void**) &d_mdlImag, voxelCount * sizeof(XFLOAT)));
	HANDLE_ERROR(cudaMalloc( (void**) &d_mdlWeight, voxelCount * sizeof(XFLOAT)));

	HANDLE_ERROR(cudaMemset( d_mdlReal, 0, voxelCount * sizeof(XFLOAT)));
	HANDLE_ERROR(cudaMemset( d_mdlImag, 0, voxelCount * sizeof(XFLOAT)));
	HANDLE_ERROR(cudaMemset( d_mdlWeight, 0, voxelCount * sizeof(XFLOAT)));

	HANDLE_ERROR(cudaStreamCreateWithPriority(&stream, cudaStreamNonBlocking, streamPriority));
}
__global__ void cuda_kernel_backproject2D(
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
	printf("ERROR : No 2D model support yet");
}

__global__ void cuda_kernel_backproject3D(
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
	unsigned pit = (gid * 4 + mid)*BACKPROJECTION4_FETCH_COUNT;
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

	__shared__ XFLOAT s_e[BACKPROJECTION4_BLOCK_SIZE*BACKPROJECTION4_FETCH_COUNT];

	__shared__ XFLOAT s_weight[BACKPROJECTION4_GROUP_SIZE*4];
	__shared__ XFLOAT s_value_real[BACKPROJECTION4_GROUP_SIZE*4];
	__shared__ XFLOAT s_value_imag[BACKPROJECTION4_GROUP_SIZE*4];

	s_weight[tid] = 0.0f;
	s_value_real[tid] = 0.0f;
	s_value_imag[tid] = 0.0f;

	for (int img = 0, b = BACKPROJECTION4_BLOCK_SIZE*BACKPROJECTION4_FETCH_COUNT; img < img_count; img ++, b += 9)
	{
		if (b+9 > BACKPROJECTION4_BLOCK_SIZE*BACKPROJECTION4_FETCH_COUNT)
		{
			__syncthreads();

			int img_9 = img*9+pit;
			if (img_9 < img_count*9)
			{
				s_e[pit+0] = g_eulers[img_9+0];
				s_e[pit+1] = g_eulers[img_9+1];
				s_e[pit+2] = g_eulers[img_9+2];
				s_e[pit+3] = g_eulers[img_9+3];
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
			g_weight[global_idx] += sum;
	}
	else if (mid == 1)
	{
		XFLOAT sum = s_value_real[gid*4 + 0] + s_value_real[gid*4 + 1] + s_value_real[gid*4 + 2] + s_value_real[gid*4 + 3];
		if (sum != 0.0f)
			g_model_real[global_idx] += sum;
	}
	else if (mid == 2)
	{
		XFLOAT sum = s_value_imag[gid*4 + 0] + s_value_imag[gid*4 + 1] + s_value_imag[gid*4 + 2] + s_value_imag[gid*4 + 3];
		if (sum != 0.0f)
			g_model_imag[global_idx] += sum;
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
	int grid_dim = ceil((float)voxelCount / BACKPROJECTION4_GROUP_SIZE);
	dim3 block_dim( BACKPROJECTION4_GROUP_SIZE * 4 );

	if(mdlZ==0)
		cuda_kernel_backproject2D<<<grid_dim,block_dim,0,stream>>>(
			d_voxelX,
			d_voxelY,
			d_voxelZ,
			d_mdlReal,
			d_mdlImag,
			d_mdlWeight,
			d_real,
			d_imag,
			d_weight,
			d_eulers,
			maxR2,
			padding_factor,
			imgX,
			imgY,
			imgX*imgY,
			imageCount,
			mdlX,
			mdlY,
			mdlInitY,
			mdlInitZ,
			voxelCount);
	else
		cuda_kernel_backproject3D<<<grid_dim,block_dim,0,stream>>>(
			d_voxelX,
			d_voxelY,
			d_voxelZ,
			d_mdlReal,
			d_mdlImag,
			d_mdlWeight,
			d_real,
			d_imag,
			d_weight,
			d_eulers,
			maxR2,
			padding_factor,
			imgX,
			imgY,
			imgX*imgY,
			imageCount,
			mdlX,
			mdlY,
			mdlInitY,
			mdlInitZ,
			voxelCount);
}


void CudaBackprojector::getMdlData(XFLOAT *implicitR, XFLOAT *implicitI, XFLOAT * implicitW)
{
	HANDLE_ERROR(cudaStreamSynchronize(stream)); //Make sure to wait for remaining kernel executions

	XFLOAT *explicitR = new XFLOAT[voxelCount];
	XFLOAT *explicitI = new XFLOAT[voxelCount];
	XFLOAT *explicitW = new XFLOAT[voxelCount];

	for (unsigned long i = 0; i < mdlXYZ; i ++)
	{
		implicitR[i] = 0.;
		implicitI[i] = 0.;
		implicitW[i] = 0.;
	}

	HANDLE_ERROR(cudaMemcpy( explicitR, d_mdlReal, voxelCount * sizeof(XFLOAT), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy( explicitI, d_mdlImag, voxelCount * sizeof(XFLOAT), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy( explicitW, d_mdlWeight, voxelCount * sizeof(XFLOAT), cudaMemcpyDeviceToHost));

	for (unsigned long i = 0; i < voxelCount; i ++)
	{
		unsigned long j = (h_voxelZ[i]-mdlInitZ)*mdlX*mdlY + (h_voxelY[i]-mdlInitY)*mdlX + h_voxelX[i];

		implicitR[j] = explicitR[i];
		implicitI[j] = explicitI[i];
		implicitW[j] = explicitW[i];
	}

	delete [] explicitR;
	delete [] explicitI;
	delete [] explicitW;
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
