#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_utils.cuh"
#include <cuda_runtime.h>
#include <signal.h>

#ifndef CUDA_DOUBLE_PRECISION

void Cuda3DProjector::setMdlData(float *real, float *imag)
{
#ifdef CUDA_DEBUG
	if (mdlXYZ == 0)
	{
        printf("DEBUG_ERROR: Model dimensions must be set with setMdlDim before call to setMdlData.");
		raise(SIGSEGV);
	}
	if (mdlReal != 0)
	{
        printf("DEBUG_ERROR: Duplicated call to setMdlData.");
		raise(SIGSEGV);
	}
#endif
    mdlReal = (void*) new cudaTextureObject_t();
    mdlImag = (void*) new cudaTextureObject_t();
	texArrayReal = (void*) new cudaArray_t();
	texArrayImag = (void*) new cudaArray_t();

	// create channel to describe data type (bits,bits,bits,bits,type)
	cudaChannelFormatDesc desc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	cudaExtent volumeSize = make_cudaExtent(mdlX, mdlY, mdlZ);

	struct cudaResourceDesc resDesc_real, resDesc_imag;
	struct cudaTextureDesc texDesc_real, texDesc_imag;

	cudaMemcpy3DParms copyParams = {0};
	copyParams.extent   = volumeSize;
	copyParams.kind     = cudaMemcpyHostToDevice;


	HANDLE_ERROR(cudaMalloc3DArray((cudaArray_t*)texArrayReal, &desc, volumeSize));

	copyParams.dstArray = *((cudaArray_t*)texArrayReal);
	copyParams.srcPtr   = make_cudaPitchedPtr((void*)real,mdlX*sizeof(float), mdlY, mdlZ);
	HANDLE_ERROR(cudaMemcpy3D(&copyParams));

	memset(&resDesc_real, 0, sizeof(cudaResourceDesc));
    resDesc_real.resType = cudaResourceTypeArray;
    resDesc_real.res.array.array = copyParams.dstArray;

    memset(&texDesc_real, 0, sizeof(cudaTextureDesc));
    texDesc_real.filterMode       = cudaFilterModeLinear;
    texDesc_real.readMode         = cudaReadModeElementType;
    texDesc_real.normalizedCoords = false;
    for(int n=0; n<3; n++)
    	texDesc_real.addressMode[n]=cudaAddressModeClamp;

	HANDLE_ERROR(cudaCreateTextureObject((cudaTextureObject_t*)mdlReal, &resDesc_real, &texDesc_real, NULL));




	HANDLE_ERROR(cudaMalloc3DArray((cudaArray_t*)texArrayImag, &desc, volumeSize));

	copyParams.dstArray = *((cudaArray_t*)texArrayImag);
	copyParams.srcPtr   = make_cudaPitchedPtr((void*)imag,mdlX*sizeof(float), mdlY, mdlZ);
	HANDLE_ERROR(cudaMemcpy3D(&copyParams));

	memset(&resDesc_imag, 0, sizeof(cudaResourceDesc));
    resDesc_imag.resType = cudaResourceTypeArray;
    resDesc_imag.res.array.array = copyParams.dstArray;

    memset(&texDesc_imag, 0, sizeof(cudaTextureDesc));
    texDesc_imag.filterMode       = cudaFilterModeLinear;
    texDesc_imag.readMode         = cudaReadModeElementType;
    texDesc_imag.normalizedCoords = false;
    for(int n=0; n<3; n++)
    	texDesc_imag.addressMode[n]=cudaAddressModeClamp;

	HANDLE_ERROR(cudaCreateTextureObject((cudaTextureObject_t*)mdlImag, &resDesc_imag, &texDesc_imag, NULL));
}

#else

void Cuda3DProjector::setMdlData(double *real, double *imag)
{
#ifdef CUDA_DEBUG
	if (mdlXYZ == 0)
	{
        printf("DEBUG_ERROR: Model dimensions must be set with setMdlDim before call to setMdlData.");
		raise(SIGSEGV);
	}
	if (mdlReal != 0)
	{
        printf("DEBUG_ERROR: Duplicated call to setMdlData.");
		raise(SIGSEGV);
	}
#endif
	CudaGlobalPtr<double> *r = new CudaGlobalPtr<double>();
	CudaGlobalPtr<double> *i = new CudaGlobalPtr<double>();

	r->h_ptr = real;
	i->h_ptr = imag;

	r->size = mdlXYZ;
	i->size = mdlXYZ;

	r->device_alloc();
	i->device_alloc();

	r->cp_to_device();
	i->cp_to_device();

	r->h_ptr = 0;
	i->h_ptr = 0;

	mdlReal = (void*) r;
	mdlImag = (void*) i;
}
:
#endif


void Cuda3DProjector::setMdlData(Complex *data)
{
	FLOAT *tmpReal = new FLOAT[mdlXYZ];
	FLOAT *tmpImag = new FLOAT[mdlXYZ];

	for (unsigned long i = 0; i < mdlXYZ; i ++)
	{
		tmpReal[i] = (FLOAT) data[i].real;
		tmpImag[i] = (FLOAT) data[i].imag;
	}

	setMdlData(tmpReal, tmpImag);

	delete [] tmpReal;
	delete [] tmpImag;
}


Cuda3DProjector::~Cuda3DProjector()
{
	if (mdlReal != 0)
	{
#ifdef CUDA_DOUBLE_PRECISION
		delete (CudaGlobalPtr<double>*) mdlReal;
		delete (CudaGlobalPtr<double>*) mdlImag;
#else
		cudaDestroyTextureObject(*(cudaTextureObject_t*) mdlReal);
		cudaDestroyTextureObject(*(cudaTextureObject_t*) mdlImag);
		delete (cudaTextureObject_t*) mdlReal;
		delete (cudaTextureObject_t*) mdlImag;

		cudaFreeArray(*((cudaArray_t*) texArrayReal));
		cudaFreeArray(*((cudaArray_t*) texArrayImag));
		delete (cudaArray_t*) texArrayReal;
		delete (cudaArray_t*) texArrayImag;
		texArrayReal = 0;
		texArrayImag = 0;
#endif
		mdlReal = 0;
		mdlImag = 0;
	}
}

