#include "src/gpu_utils/cuda_backprojector.h"
#include "src/gpu_utils/cuda_utils.cuh"
#include <cuda_runtime.h>
#include <signal.h>


void Cuda3DBackprojector::initMdl()
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

	CudaGlobalPtr<int> *xs = new CudaGlobalPtr<int>(mdlXYZ);
	CudaGlobalPtr<int> *ys = new CudaGlobalPtr<int>(mdlXYZ);
	CudaGlobalPtr<int> *zs = new CudaGlobalPtr<int>(mdlXYZ);

	for (int x = 0; x < mdlX; x ++)
	{
		for (int y = mdlInitY; y < mdlY; y++)
		{
			for (int z = mdlInitZ; z < mdlZ; z++)
			{
				if (x*x + y*y + z*z <= maxR2 * padding_factor * padding_factor * 1.2f)
				{
					(*xs)[voxelCount] = x;
					(*ys)[voxelCount] = y;
					(*zs)[voxelCount] = z;
					voxelCount ++;
				}
			}
		}
	}

	xs->put_on_device(voxelCount);
	ys->put_on_device(voxelCount);
	zs->put_on_device(voxelCount);

	xs->free_host();
	ys->free_host();
	zs->free_host();

	voxelX = (void*) xs;
	voxelY = (void*) ys;
	voxelZ = (void*) zs;


	CudaGlobalPtr<FLOAT> *r = new CudaGlobalPtr<FLOAT>;
	CudaGlobalPtr<FLOAT> *i = new CudaGlobalPtr<FLOAT>;
	CudaGlobalPtr<FLOAT> *w = new CudaGlobalPtr<FLOAT>;

	r->device_alloc(mdlXYZ);
	i->device_alloc(mdlXYZ);
	w->device_alloc(mdlXYZ);

	r->device_init(0);
	i->device_init(0);
	w->device_init(0);

	mdlReal = (void*) r;
	mdlImag = (void*) i;
	mdlWeight = (void*) w;
}

__global__ void cuda_kernel_backproject(
		int *g_xs,
		int *g_ys,
		int *g_zs,
		FLOAT *g_model_real,
		FLOAT *g_model_imag,
		FLOAT *g_weight,
		FLOAT *g_wavgs_real,
		FLOAT *g_wavgs_imag,
		FLOAT *g_Fweights,
		FLOAT *g_eulers,
		int max_r2,
		FLOAT scale2,
		unsigned img_x, unsigned img_y, unsigned img_xy, unsigned long img_count,
		unsigned mdl_x, unsigned mdl_y, int mdl_inity, int mdl_initz,
		int N)
{
	unsigned gid = threadIdx.x / 4; //Group id
	unsigned mid = threadIdx.x % 4; //Member id
	unsigned gm = gid * 4 + mid; //TODO this should be just threadIdx.x
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
	FLOAT d, w;
	FLOAT xp,yp,zp;
	int x,y,idx;

	__shared__ FLOAT s_e[BACKPROJECTION4_BLOCK_SIZE*BACKPROJECTION4_FETCH_COUNT];

	__shared__ FLOAT s_weight[BACKPROJECTION4_GROUP_SIZE*4];
	__shared__ FLOAT s_value_real[BACKPROJECTION4_GROUP_SIZE*4];
	__shared__ FLOAT s_value_imag[BACKPROJECTION4_GROUP_SIZE*4];

	s_weight[gm] = 0.0f;
	s_value_real[gm] = 0.0f;
	s_value_imag[gm] = 0.0f;

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

				s_weight[gm] += w * d;
				s_value_real[gm] += g_wavgs_real[idx] * d;
				if (is_neg_x) s_value_imag[gm] -= g_wavgs_imag[idx] * d;
				else          s_value_imag[gm] += g_wavgs_imag[idx] * d;
			}
		}
	}

	__syncthreads();

	if (mid == 0)
	{
		FLOAT sum = s_weight[gid*4 + 0] + s_weight[gid*4 + 1] + s_weight[gid*4 + 2] + s_weight[gid*4 + 3];
		if (sum != 0.0f)
			g_weight[(Z-mdl_initz)*mdl_x*mdl_y + (Y-mdl_inity)*mdl_x + X] += sum;
	}
	else if (mid == 1)
	{
		FLOAT sum = s_value_real[gid*4 + 0] + s_value_real[gid*4 + 1] + s_value_real[gid*4 + 2] + s_value_real[gid*4 + 3];
		if (sum != 0.0f)
			g_model_real[(Z-mdl_initz)*mdl_x*mdl_y + (Y-mdl_inity)*mdl_x + X] += sum;
	}
	else if (mid == 2)
	{
		FLOAT sum = s_value_imag[gid*4 + 0] + s_value_imag[gid*4 + 1] + s_value_imag[gid*4 + 2] + s_value_imag[gid*4 + 3];
		if (sum != 0.0f)
			g_model_imag[(Z-mdl_initz)*mdl_x*mdl_y + (Y-mdl_inity)*mdl_x + X] += sum;
	}
}

void Cuda3DBackprojector::backproject(
		FLOAT *d_real,
		FLOAT *d_imag,
		FLOAT *d_weight,
		FLOAT *d_eulers,
		int imgX,
		int imgY,
		unsigned long imageCount)
{
	int grid_dim = ceil((float)voxelCount / BACKPROJECTION4_GROUP_SIZE);
	dim3 block_dim( BACKPROJECTION4_GROUP_SIZE * 4 );

	cuda_kernel_backproject<<<grid_dim,block_dim>>>(
			~*((CudaGlobalPtr<int>*) voxelX),
			~*((CudaGlobalPtr<int>*) voxelY),
			~*((CudaGlobalPtr<int>*) voxelZ),
			~*((CudaGlobalPtr<FLOAT>*) mdlReal),
			~*((CudaGlobalPtr<FLOAT>*) mdlImag),
			~*((CudaGlobalPtr<FLOAT>*) mdlWeight),
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


void Cuda3DBackprojector::getMdlData(FLOAT *real, FLOAT *imag, FLOAT * weights)
{
	((CudaGlobalPtr<FLOAT>*) mdlReal)->h_ptr = real;
	((CudaGlobalPtr<FLOAT>*) mdlImag)->h_ptr = imag;
	((CudaGlobalPtr<FLOAT>*) mdlWeight)->h_ptr = weights;

	((CudaGlobalPtr<FLOAT>*) mdlReal)->cp_to_host();
	((CudaGlobalPtr<FLOAT>*) mdlImag)->cp_to_host();
	((CudaGlobalPtr<FLOAT>*) mdlWeight)->cp_to_host();
}


void Cuda3DBackprojector::getMdlData(Complex *data, double * weights)
{
	FLOAT *r = new FLOAT[mdlXYZ];
	FLOAT *i = new FLOAT[mdlXYZ];
	FLOAT *w = new FLOAT[mdlXYZ];

	getMdlData(r, i, w);

	for (unsigned long n = 0; n < mdlXYZ; n++)
	{
		data[n].real = (double) r[n];
		data[n].imag = (double) i[n];
		weights[n] = (double) w[n];
	}

	delete r;
	delete i;
	delete w;
}


Cuda3DBackprojector::~Cuda3DBackprojector()
{
	if (voxelCount != 0)
	{
		delete (CudaGlobalPtr<int>*) voxelX;
		delete (CudaGlobalPtr<int>*) voxelY;
		delete (CudaGlobalPtr<int>*) voxelZ;

		delete (CudaGlobalPtr<FLOAT>*) mdlReal;
		delete (CudaGlobalPtr<FLOAT>*) mdlImag;
		delete (CudaGlobalPtr<FLOAT>*) mdlWeight;
	}
}




