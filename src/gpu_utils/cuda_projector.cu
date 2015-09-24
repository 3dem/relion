#include "src/gpu_utils/cuda_projector.h"
#include <signal.h>


void CudaProjector::setMdlDim(
		int xdim, int ydim, int zdim,
		int inity, int initz,
		int maxr, int paddingFactor)
{
	bool resizeTexure(true);

	if (xdim == mdlX &&
		ydim == mdlY &&
		zdim == mdlZ &&
		inity == mdlInitY &&
		initz == mdlInitZ &&
		maxr == mdlMaxR &&
		paddingFactor == padding_factor)
		resizeTexure = false;

	mdlX = xdim;
	mdlY = ydim;
	if(zdim==1)
		mdlZ=0;
	else
		mdlZ = zdim;
	mdlXYZ = xdim*ydim*zdim;
	mdlInitY = inity;
	mdlInitZ = initz;
	mdlMaxR = maxr;
	padding_factor = paddingFactor;

	if (! resizeTexure) return;

	clear();

#ifndef CUDA_DOUBLE_PRECISION

	mdlReal = new cudaTextureObject_t();
	mdlImag = new cudaTextureObject_t();

	// create channel to describe data type (bits,bits,bits,bits,type)
	cudaChannelFormatDesc desc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

	struct cudaResourceDesc resDesc_real, resDesc_imag;
	struct cudaTextureDesc  texDesc;
	// -- Zero all data in objects handlers
	memset(&resDesc_real, 0, sizeof(cudaResourceDesc));
	memset(&resDesc_imag, 0, sizeof(cudaResourceDesc));
	memset(&texDesc, 0, sizeof(cudaTextureDesc));

	if(mdlZ!=0)  // 3D model
	{
		texArrayReal = new cudaArray_t();
		texArrayImag = new cudaArray_t();

		// -- make extents for automatic pitch:ing (aligment) of allocated 3D arrays
		cudaExtent volumeSize = make_cudaExtent(mdlX, mdlY, mdlZ);

		// -- Allocate and copy data using very celver CUDA memcpy-functions
		HANDLE_ERROR(cudaMalloc3DArray(texArrayReal, &desc, volumeSize));
		HANDLE_ERROR(cudaMalloc3DArray(texArrayImag, &desc, volumeSize));

		// -- Descriptors of the channel(s) in the texture(s)
		resDesc_real.res.array.array = *texArrayReal;
		resDesc_imag.res.array.array = *texArrayImag;
		resDesc_real.resType = cudaResourceTypeArray;
		resDesc_imag.resType = cudaResourceTypeArray;
	}
	else // 2D model
	{
		HANDLE_ERROR(cudaMallocPitch(&texArrayReal2D, &pitch2D, sizeof(float)*mdlX,mdlY));
		HANDLE_ERROR(cudaMallocPitch(&texArrayImag2D, &pitch2D, sizeof(float)*mdlX,mdlY));

		// -- Descriptors of the channel(s) in the texture(s)
		resDesc_real.resType = cudaResourceTypePitch2D;
		resDesc_real.res.pitch2D.devPtr = texArrayReal2D;
		resDesc_real.res.pitch2D.pitchInBytes =  pitch2D;
		resDesc_real.res.pitch2D.width = mdlX;
		resDesc_real.res.pitch2D.height = mdlY;
		resDesc_real.res.pitch2D.desc = desc;
		// -------------------------------------------------
		resDesc_imag.resType = cudaResourceTypePitch2D;
		resDesc_imag.res.pitch2D.devPtr = texArrayImag2D;
		resDesc_imag.res.pitch2D.pitchInBytes =  pitch2D;
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

#else

	DEBUG_HANDLE_ERROR(cudaMalloc( (void**) &mdlReal, mdlXYZ * sizeof(double)));
	DEBUG_HANDLE_ERROR(cudaMalloc( (void**) &mdlImag, mdlXYZ * sizeof(double)));

#endif
}

void CudaProjector::setMdlData(XFLOAT *real, XFLOAT *imag)
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

#ifndef CUDA_DOUBLE_PRECISION

	if(mdlZ!=0)  // 3D model
	{
		// -- make extents for automatic pitch:ing (aligment) of allocated 3D arrays
		cudaExtent volumeSize = make_cudaExtent(mdlX, mdlY, mdlZ);
		cudaMemcpy3DParms copyParams = {0};
		copyParams.extent = volumeSize;
		copyParams.kind   = cudaMemcpyHostToDevice;

		// -- Copy data
		copyParams.dstArray = *texArrayReal;
		copyParams.srcPtr   = make_cudaPitchedPtr(real, mdlX * sizeof(float), mdlY, mdlZ);
		DEBUG_HANDLE_ERROR(cudaMemcpy3D(&copyParams));
		copyParams.dstArray = *texArrayImag;
		copyParams.srcPtr   = make_cudaPitchedPtr(imag, mdlX * sizeof(float), mdlY, mdlZ);
		DEBUG_HANDLE_ERROR(cudaMemcpy3D(&copyParams));
	}
	else // 2D model
	{
		DEBUG_HANDLE_ERROR(cudaMemcpy2D(texArrayReal2D, pitch2D, real, sizeof(float) * mdlX, sizeof(float) * mdlX, mdlY, cudaMemcpyHostToDevice));
		DEBUG_HANDLE_ERROR(cudaMemcpy2D(texArrayImag2D, pitch2D, imag, sizeof(float) * mdlX, sizeof(float) * mdlX, mdlY, cudaMemcpyHostToDevice));
	}

#else

	DEBUG_HANDLE_ERROR(cudaMemcpy( mdlReal, real, mdlXYZ * sizeof(XFLOAT), cudaMemcpyHostToDevice));
	DEBUG_HANDLE_ERROR(cudaMemcpy( mdlImag, imag, mdlXYZ * sizeof(XFLOAT), cudaMemcpyHostToDevice));

#endif

}


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


void CudaProjector::clear()
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

		if(mdlZ!=0) //3D case
		{
			cudaFreeArray(*texArrayReal);
			cudaFreeArray(*texArrayImag);
			delete texArrayReal;
			delete texArrayImag;
		}
		else //2D case
		{
			cudaFree(texArrayReal2D);
			cudaFree(texArrayImag2D);
		}

		texArrayReal = 0;
		texArrayImag = 0;
#endif
		mdlReal = 0;
		mdlImag = 0;
	}
}

