#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_utils_stl.cuh"
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

	struct cudaResourceDesc resDesc_real, resDesc_imag;
	struct cudaTextureDesc  texDesc;
	// -- Zero all data in objects handlers
	memset(&resDesc_real, 0, sizeof(cudaResourceDesc));
	memset(&resDesc_imag, 0, sizeof(cudaResourceDesc));
	memset(&texDesc, 0, sizeof(cudaTextureDesc));

	//
	if(mdlZ!=0)  // then we are using a 3D reference
	{
		// -- make extents for automatic pitch:ing (aligment) of allocated 3D arrays
		cudaExtent volumeSize = make_cudaExtent(mdlX, mdlY, mdlZ);
		cudaMemcpy3DParms copyParams_real = {0}, copyParams_imag = {0};
		copyParams_real.extent = volumeSize;
		copyParams_real.kind   = cudaMemcpyHostToDevice;
		copyParams_imag.extent = volumeSize;
		copyParams_imag.kind   = cudaMemcpyHostToDevice;

		// -- Allocate and copy data using very celver CUDA memcpy-functions
		HANDLE_ERROR(cudaMalloc3DArray(texArrayReal, &desc, volumeSize));
		copyParams_real.dstArray = *texArrayReal;
		copyParams_real.srcPtr   = make_cudaPitchedPtr(real,mdlX*sizeof(float), mdlY, mdlZ);
		HANDLE_ERROR(cudaMemcpy3D(&copyParams_real));
		// ------------------------------------------
		HANDLE_ERROR(cudaMalloc3DArray(texArrayImag, &desc, volumeSize));
		copyParams_imag.dstArray = *texArrayImag;
		copyParams_imag.srcPtr   = make_cudaPitchedPtr(imag,mdlX*sizeof(float), mdlY, mdlZ);
		HANDLE_ERROR(cudaMemcpy3D(&copyParams_imag));

		// -- Descriptors of the channel(s) in the texture(s)
		resDesc_real.res.array.array = copyParams_real.dstArray;
		resDesc_imag.res.array.array = copyParams_imag.dstArray;
		resDesc_real.resType = cudaResourceTypeArray;
		resDesc_imag.resType = cudaResourceTypeArray;
	}
	else // then we are using a 2D reference
	{
		size_t pitch;

		// -- allocate pitched  (aligned) memory positions, and copy data into them
		HANDLE_ERROR(cudaMallocPitch(&texArrayReal, &pitch, sizeof(float)*mdlX,mdlY));
		HANDLE_ERROR(cudaMemcpy2D(texArrayReal, pitch, real, sizeof(float)*mdlX, sizeof(float)*mdlX, mdlY, cudaMemcpyHostToDevice));
		// ------------------------------------------------
		HANDLE_ERROR(cudaMallocPitch(&texArrayImag, &pitch, sizeof(float)*mdlX,mdlY));
		HANDLE_ERROR(cudaMemcpy2D(texArrayImag, pitch, imag, sizeof(float)*mdlX, sizeof(float)*mdlX, mdlY, cudaMemcpyHostToDevice));


		// -- Descriptors of the channel(s) in the texture(s)
		resDesc_real.resType = cudaResourceTypePitch2D;
		resDesc_real.res.pitch2D.devPtr = texArrayReal;
		resDesc_real.res.pitch2D.pitchInBytes =  pitch;
		resDesc_real.res.pitch2D.width = mdlX;
		resDesc_real.res.pitch2D.height = mdlY;
		resDesc_real.res.pitch2D.desc = desc;
		// -------------------------------------------------
		resDesc_imag.resType = cudaResourceTypePitch2D;
		resDesc_imag.res.pitch2D.devPtr = texArrayImag;
		resDesc_imag.res.pitch2D.pitchInBytes =  pitch;
		resDesc_imag.res.pitch2D.width = mdlX;
		resDesc_imag.res.pitch2D.height = mdlY;
	    resDesc_imag.res.pitch2D.desc = desc;
	}

	// -- Decriptors of the texture(s) and methods used for reading it(them) --
	texDesc.filterMode       = cudaFilterModeLinear;
	texDesc.readMode         = cudaReadModeElementType;
	texDesc.normalizedCoords = false;
	for(int n=0; n<3; n++)
		texDesc.addressMode[n]=cudaAddressModeClamp;

	// -- Create texture object(s)
	HANDLE_ERROR(cudaCreateTextureObject(mdlReal, &resDesc_real, &texDesc, NULL));
	HANDLE_ERROR(cudaCreateTextureObject(mdlImag, &resDesc_imag, &texDesc, NULL));
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

	HANDLE_ERROR(cudaMemcpy( mdlReal, real, mdlXYZ * sizeof(XFLOAT), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy( mdlImag, imag, mdlXYZ * sizeof(XFLOAT), cudaMemcpyHostToDevice));

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

		if(mdlZ!=0)
		{
			cudaFreeArray(*texArrayReal);
			cudaFreeArray(*texArrayImag);
			delete texArrayReal;
			delete texArrayImag;
		}
		else
		{
			cudaFree(texArrayReal);
			cudaFree(texArrayImag);
		}
		texArrayReal = 0;
		texArrayImag = 0;
#endif
		mdlReal = 0;
		mdlImag = 0;
	}
}

