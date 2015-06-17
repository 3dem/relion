#include "src/gpu_utils/cuda_projector.h"
#include <signal.h>

#ifndef CUDA_DOUBLE_PRECISION

void setupTextureArray(
		FLOAT* hostArray,
		cudaArray* deviceArray,
		cudaExtent volumeSize,
		cudaResourceDesc &resDesc,
		cudaTextureDesc &texDesc)
{
	cudaMemcpy3DParms copyParams = {0};
	copyParams.extent   = volumeSize;
	copyParams.kind     = cudaMemcpyHostToDevice;
	copyParams.dstArray = deviceArray;
	copyParams.srcPtr   = make_cudaPitchedPtr(hostArray,volumeSize.width*sizeof(FLOAT), volumeSize.height, volumeSize.depth);
	cudaMemcpy3D(&copyParams);

    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = deviceArray;

    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.filterMode       = cudaFilterModeLinear;
    texDesc.readMode         = cudaReadModeElementType;
    texDesc.normalizedCoords = false;
    for(int n=0; n<3; n++)
    	texDesc.addressMode[n]=cudaAddressModeClamp;

}

void Cuda3DProjector::setData(float *real, float *imag)
{
#ifdef CUDA_DEBUG
	if (mdlReal != 0)
	{
		printf("DEBUG_ERROR: Model data set twice in Cuda3DProjector.\n");
		raise(SIGSEGV);
	}
#endif

	// create channel to describe data type (bits,bits,bits,bits,type)
	cudaChannelFormatDesc channel = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	cudaExtent volumeSize = make_cudaExtent(mdlX, mdlY, mdlZ);

	//allocate device memory for cuda 3D array
	cudaMalloc3DArray(&mdlTextureArrayReal, &channel, volumeSize);
	cudaMalloc3DArray(&mdlTextureArrayImag, &channel, volumeSize);
	struct cudaResourceDesc resDesc_real, resDesc_imag;
	struct cudaTextureDesc texDesc_real, texDesc_imag;

	setupTextureArray(real, mdlTextureArrayReal, volumeSize, resDesc_real, texDesc_real);
	setupTextureArray(imag, mdlTextureArrayImag, volumeSize, resDesc_imag, texDesc_imag);

	cudaCreateTextureObject(&mdlReal, &resDesc_real, &texDesc_real, NULL);
	cudaCreateTextureObject(&mdlImag, &resDesc_imag, &texDesc_imag, NULL);
}

#else

void Cuda3DProjector::setData(double *real, double *imag)
{
#ifdef CUDA_DEBUG
	if (mdlReal.h_ptr != 0)
	{
		printf("DEBUG_ERROR: Model data set twice in CudaProjector.");
		raise(SIGSEGV);
	}
#endif

	mdlReal.h_ptr = real;
	mdlImag.h_ptr = imag;

	mdlReal.size = mdlXYZ;
	mdlImag.size = mdlXYZ;

	mdlReal.device_alloc();
	mdlImag.device_alloc();

	mdlReal.cp_to_device();
	mdlImag.cp_to_device();

	mdlReal.h_ptr = 0;
	mdlImag.h_ptr = 0;
}

#endif


void Cuda3DProjector::setData(Complex *data)
{
	FLOAT *tmpReal = new FLOAT[mdlXYZ];
	FLOAT *tmpImag = new FLOAT[mdlXYZ];

	for (unsigned long i = 0; i < mdlXYZ; i ++)
	{
		tmpReal[i] = (FLOAT) data[i].real;
		tmpImag[i] = (FLOAT) data[i].imag;
	}

	setData(tmpReal, tmpImag);

	delete [] tmpReal;
	delete [] tmpImag;
}




__device__ inline void Cuda3DProjector::Kernel::project(
		int x, int y,
		FLOAT e0,
		FLOAT e1,
		FLOAT e3,
		FLOAT e4,
		FLOAT e6,
		FLOAT e7,
		FLOAT &real,
		FLOAT &imag)
{
	bool is_neg_x;
	long int r2;

	// Dont search beyond square with side max_r
	if (y > maxR)
	{
		if (y >= imgY - maxR)
			y = y - imgY;
		else
			x=maxR;
	}

	r2 = x*x + y*y;
	if (r2 <= maxR2)
	{
		FLOAT xp = (e0 * x + e1 * y ) * padding_factor;
		FLOAT yp = (e3 * x + e4 * y ) * padding_factor;
		FLOAT zp = (e6 * x + e7 * y ) * padding_factor;

		// Only asymmetric half is stored
		if (xp < 0)
		{
			// Get complex conjugated hermitian symmetry pair
			xp = -xp;
			yp = -yp;
			zp = -zp;
			is_neg_x = true;
		}
		else
		{
			is_neg_x = false;
		}
		yp -= mdlInitY;
		zp -= mdlInitZ;

		real=tex3D<FLOAT>(mdlReal,xp+0.5f,yp+0.5f,zp+0.5f);
		imag=tex3D<FLOAT>(mdlImag,xp+0.5f,yp+0.5f,zp+0.5f);

		if (is_neg_x)
		{
			imag = -imag;
		}
	}
	else
	{
		real=0.0f;
		imag=0.0f;
	}
}


