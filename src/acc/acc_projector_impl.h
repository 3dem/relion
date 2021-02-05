#include "src/acc/acc_projector.h"
#include <signal.h>


bool AccProjector::setMdlDim(
		int xdim, int ydim, int zdim,
		int inity, int initz,
		int maxr, XFLOAT paddingFactor)
{
	if(zdim == 1) zdim = 0;

	if (xdim == mdlX &&
		ydim == mdlY &&
		zdim == mdlZ &&
		inity == mdlInitY &&
		initz == mdlInitZ &&
		maxr == mdlMaxR &&
		paddingFactor == padding_factor)
		return false;

	clear();

	mdlX = xdim;
	mdlY = ydim;
	mdlZ = zdim;
	if(zdim == 0)
		mdlXYZ = (size_t)xdim*(size_t)ydim;
	else
		mdlXYZ = (size_t)xdim*(size_t)ydim*(size_t)zdim;
	mdlInitY = inity;
	mdlInitZ = initz;
	mdlMaxR = maxr;
	padding_factor = paddingFactor;

#ifndef PROJECTOR_NO_TEXTURES

	mdlReal = new cudaTextureObject_t();
	mdlImag = new cudaTextureObject_t();

	// create channel to describe data type (bits,bits,bits,bits,type)
	cudaChannelFormatDesc desc;

	desc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

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


		// -- Allocate and copy data using very clever CUDA memcpy-functions
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
		HANDLE_ERROR(cudaMallocPitch(&texArrayReal2D, &pitch2D, sizeof(XFLOAT)*mdlX,mdlY));
		HANDLE_ERROR(cudaMallocPitch(&texArrayImag2D, &pitch2D, sizeof(XFLOAT)*mdlX,mdlY));

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
#ifdef _CUDA_ENABLED
	DEBUG_HANDLE_ERROR(cudaMalloc( (void**) &mdlReal, mdlXYZ * sizeof(XFLOAT)));
	DEBUG_HANDLE_ERROR(cudaMalloc( (void**) &mdlImag, mdlXYZ * sizeof(XFLOAT)));
#else
	mdlComplex = NULL;
#endif
#endif
	return true;
}

void AccProjector::initMdl(XFLOAT *real, XFLOAT *imag)
{
#ifdef DEBUG_CUDA
	if (mdlXYZ == 0)
	{
        printf("DEBUG_ERROR: Model dimensions must be set with setMdlDim before call to setMdlData.");
		CRITICAL(ERR_MDLDIM);
	}
#ifdef _CUDA_ENABLED
	if (mdlReal == NULL)
	{
        printf("DEBUG_ERROR: initMdl called before call to setMdlData.");
		CRITICAL(ERR_MDLSET);
	}
#else
	if (mdlComplex == NULL)
	{
        printf("DEBUG_ERROR: initMdl called before call to setMdlData.");
		CRITICAL(ERR_MDLSET);
	}
#endif
#endif

#ifndef PROJECTOR_NO_TEXTURES
	if(mdlZ!=0)  // 3D model
	{
		// -- make extents for automatic pitching (aligment) of allocated 3D arrays
		cudaMemcpy3DParms copyParams = {0};
		copyParams.extent = make_cudaExtent(mdlX, mdlY, mdlZ);
		copyParams.kind   = cudaMemcpyHostToDevice;

		// -- Copy data
		copyParams.dstArray = *texArrayReal;
		copyParams.srcPtr   = make_cudaPitchedPtr(real, mdlX * sizeof(XFLOAT), mdlY, mdlZ);
		DEBUG_HANDLE_ERROR(cudaMemcpy3D(&copyParams));
		copyParams.dstArray = *texArrayImag;
		copyParams.srcPtr   = make_cudaPitchedPtr(imag, mdlX * sizeof(XFLOAT), mdlY, mdlZ);
		DEBUG_HANDLE_ERROR(cudaMemcpy3D(&copyParams));
	}
	else // 2D model
	{
		DEBUG_HANDLE_ERROR(cudaMemcpy2D(texArrayReal2D, pitch2D, real, sizeof(XFLOAT) * mdlX, sizeof(XFLOAT) * mdlX, mdlY, cudaMemcpyHostToDevice));
		DEBUG_HANDLE_ERROR(cudaMemcpy2D(texArrayImag2D, pitch2D, imag, sizeof(XFLOAT) * mdlX, sizeof(XFLOAT) * mdlX, mdlY, cudaMemcpyHostToDevice));
	}
#else
#ifdef _CUDA_ENABLED
	DEBUG_HANDLE_ERROR(cudaMemcpy( mdlReal, real, mdlXYZ * sizeof(XFLOAT), cudaMemcpyHostToDevice));
	DEBUG_HANDLE_ERROR(cudaMemcpy( mdlImag, imag, mdlXYZ * sizeof(XFLOAT), cudaMemcpyHostToDevice));
#else
	std::complex<XFLOAT> *pData = mdlComplex;
    for(size_t i=0; i<mdlXYZ; i++) {
		std::complex<XFLOAT> arrayval(*real ++, *imag ++);
		pData[i] = arrayval;
    }
#endif
#endif

}

#ifndef _CUDA_ENABLED
void AccProjector::initMdl(std::complex<XFLOAT> *data)
{
	mdlComplex = data;  // No copy needed - everyone shares the complex reference arrays
	externalFree = 1;   // This is shared memory freed outside the projector
}
#endif

void AccProjector::initMdl(Complex *data)
{
	XFLOAT *tmpReal;
	XFLOAT *tmpImag;
	if (posix_memalign((void **)&tmpReal, MEM_ALIGN, mdlXYZ * sizeof(XFLOAT))) CRITICAL(RAMERR);
	if (posix_memalign((void **)&tmpImag, MEM_ALIGN, mdlXYZ * sizeof(XFLOAT))) CRITICAL(RAMERR);


	for (size_t i = 0; i < mdlXYZ; i ++)
	{
		tmpReal[i] = (XFLOAT) data[i].real;
		tmpImag[i] = (XFLOAT) data[i].imag;
	}

	initMdl(tmpReal, tmpImag);

	free(tmpReal);
	free(tmpImag);
}

void AccProjector::clear()
{
#ifdef _CUDA_ENABLED
	if (mdlReal != 0)
	{
#ifndef PROJECTOR_NO_TEXTURES
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
			HANDLE_ERROR(cudaFree(texArrayReal2D));
			HANDLE_ERROR(cudaFree(texArrayImag2D));
		}

		texArrayReal = 0;
		texArrayImag = 0;
#else
		cudaFree(mdlReal);
		cudaFree(mdlImag);
#endif
		mdlReal = 0;
		mdlImag = 0;
	}

	mdlX = 0;
	mdlY = 0;
	mdlZ = 0;
	mdlXYZ = 0;
	mdlInitY = 0;
	mdlInitZ = 0;
	mdlMaxR = 0;
	padding_factor = 0;
	allocaton_size = 0;

#else // ifdef CUDA
	if ((mdlComplex != NULL) && (externalFree == 0))
	{
		delete [] mdlComplex;
		mdlComplex = NULL;
	}
#endif  // ifdef CUDA
}
