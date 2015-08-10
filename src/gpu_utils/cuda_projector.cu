#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_utils.cuh"
#include <signal.h>

#ifndef CUDA_DOUBLE_PRECISION

void CudaProjector::setMdlData(float *real, float *imag)
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
    mdlReal = new cudaTextureObject_t();
    mdlImag = new cudaTextureObject_t();
	texArrayReal = new cudaArray_t();
	texArrayImag = new cudaArray_t();

	// create channel to describe data type (bits,bits,bits,bits,type)
	cudaChannelFormatDesc desc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	cudaExtent volumeSize = make_cudaExtent(mdlX, mdlY, mdlZ);

	struct cudaResourceDesc resDesc_real, resDesc_imag;
	struct cudaTextureDesc texDesc_real, texDesc_imag;

	cudaMemcpy3DParms copyParams = {0};
	copyParams.extent = volumeSize;
	copyParams.kind   = cudaMemcpyHostToDevice;


	HANDLE_ERROR(cudaMalloc3DArray(texArrayReal, &desc, volumeSize));

	copyParams.dstArray = *texArrayReal;
	copyParams.srcPtr   = make_cudaPitchedPtr(real,mdlX*sizeof(float), mdlY, mdlZ);
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

	HANDLE_ERROR(cudaCreateTextureObject(mdlReal, &resDesc_real, &texDesc_real, NULL));




	HANDLE_ERROR(cudaMalloc3DArray(texArrayImag, &desc, volumeSize));

	copyParams.dstArray = *texArrayImag;
	copyParams.srcPtr   = make_cudaPitchedPtr(imag,mdlX*sizeof(float), mdlY, mdlZ);
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

	HANDLE_ERROR(cudaCreateTextureObject(mdlImag, &resDesc_imag, &texDesc_imag, NULL));
}

#else

void CudaProjector::setMdlData(double *real, double *imag)
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

	HANDLE_ERROR(cudaMalloc( (void**) &mdlReal, mdlXYZ * sizeof(double)));
	HANDLE_ERROR(cudaMalloc( (void**) &mdlImag, mdlXYZ * sizeof(double)));

	HANDLE_ERROR(cudaMemcpy( mdlReal, real, mdlXYZ * sizeof(FLOAT), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy( mdlImag, imag, mdlXYZ * sizeof(FLOAT), cudaMemcpyHostToDevice));

}
#endif


void CudaProjector::setMdlData(Complex *data)
{
	XFLOAT *tmpReal = new XFLOAT[mdlXYZ];
	XFLOAT *tmpImag = new XFLOAT[mdlXYZ];

	for (unsigned long i = 0; i < mdlXYZ; i ++)
	{
		tmpReal[i] = (XFLOAT) data[i].real;
		tmpImag[i] = (XFLOAT) data[i].imag;
	}

	setMdlData(tmpReal, tmpImag);

	delete [] tmpReal;
	delete [] tmpImag;
}


CudaProjector::~CudaProjector()
{
	if (mdlReal != 0)
	{
#ifdef CUDA_DOUBLE_PRECISION
		cudaFree(mdlReal);
		cudaFree(mdlImag);
#else
		cudaDestroyTextureObject(*mdlReal);
		cudaDestroyTextureObject(*mdlImag);
		delete mdlReal;
		delete mdlImag;

		cudaFreeArray(*texArrayReal);
		cudaFreeArray(*texArrayImag);
		delete texArrayReal;
		delete texArrayImag;
		texArrayReal = 0;
		texArrayImag = 0;
#endif
		mdlReal = 0;
		mdlImag = 0;
	}
}

