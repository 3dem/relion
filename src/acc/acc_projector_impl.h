#include "src/acc/acc_projector.h"
#include <signal.h>


bool AccProjector::setMdlDim(
#ifdef _SYCL_ENABLED
		deviceStream_t dev,
#endif
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
	#ifdef _CUDA_ENABLED
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
	#elif _HIP_ENABLED
		mdlReal = new hipTextureObject_t();
		mdlImag = new hipTextureObject_t();

		// create channel to describe data type (bits,bits,bits,bits,type)
		hipChannelFormatDesc desc;

		desc = hipCreateChannelDesc(32, 0, 0, 0, hipChannelFormatKindFloat);

		struct hipResourceDesc resDesc_real, resDesc_imag;
		struct hipTextureDesc  texDesc;
		// -- Zero all data in objects handlers
		memset(&resDesc_real, 0, sizeof(hipResourceDesc));
		memset(&resDesc_imag, 0, sizeof(hipResourceDesc));
		memset(&texDesc, 0, sizeof(hipTextureDesc));

		if(mdlZ!=0)  // 3D model
		{
			texArrayReal = new hipArray_t();
			texArrayImag = new hipArray_t();

			// -- make extents for automatic pitch:ing (aligment) of allocated 3D arrays
			hipExtent volumeSize = make_hipExtent(mdlX, mdlY, mdlZ);


			// -- Allocate and copy data using very clever HIP memcpy-functions
			HANDLE_ERROR(hipMalloc3DArray(texArrayReal, &desc, volumeSize, hipArrayDefault));
			HANDLE_ERROR(hipMalloc3DArray(texArrayImag, &desc, volumeSize, hipArrayDefault));

			// -- Descriptors of the channel(s) in the texture(s)
			resDesc_real.res.array.array = *texArrayReal;
			resDesc_imag.res.array.array = *texArrayImag;
			resDesc_real.resType = hipResourceTypeArray;
			resDesc_imag.resType = hipResourceTypeArray;
		}
		else // 2D model
		{
			HANDLE_ERROR(hipMallocPitch(reinterpret_cast<void**>(&texArrayReal2D), &pitch2D, sizeof(XFLOAT)*mdlX,mdlY));
			HANDLE_ERROR(hipMallocPitch(reinterpret_cast<void**>(&texArrayImag2D), &pitch2D, sizeof(XFLOAT)*mdlX,mdlY));

			// -- Descriptors of the channel(s) in the texture(s)
			resDesc_real.resType = hipResourceTypePitch2D;
			resDesc_real.res.pitch2D.devPtr = texArrayReal2D;
			resDesc_real.res.pitch2D.pitchInBytes =  pitch2D;
			resDesc_real.res.pitch2D.width = mdlX;
			resDesc_real.res.pitch2D.height = mdlY;
			resDesc_real.res.pitch2D.desc = desc;
			// -------------------------------------------------
			resDesc_imag.resType = hipResourceTypePitch2D;
			resDesc_imag.res.pitch2D.devPtr = texArrayImag2D;
			resDesc_imag.res.pitch2D.pitchInBytes =  pitch2D;
			resDesc_imag.res.pitch2D.width = mdlX;
			resDesc_imag.res.pitch2D.height = mdlY;
			resDesc_imag.res.pitch2D.desc = desc;
		}

		// -- Decriptors of the texture(s) and methods used for reading it(them) --
		if(mdlZ!=0)  // 3D model
			texDesc.filterMode       = hipFilterModePoint;
		else // 2D model
			texDesc.filterMode       = hipFilterModeLinear;
		texDesc.readMode         = hipReadModeElementType;
		texDesc.normalizedCoords = false;

		for(int n=0; n<3; n++)
			texDesc.addressMode[n]=hipAddressModeClamp;

		// -- Create texture object(s)
		HANDLE_ERROR(hipCreateTextureObject(mdlReal, &resDesc_real, &texDesc, NULL));
		HANDLE_ERROR(hipCreateTextureObject(mdlImag, &resDesc_imag, &texDesc, NULL));

	#endif
#else
#ifdef _CUDA_ENABLED
	DEBUG_HANDLE_ERROR(cudaMalloc( (void**) &mdlReal, mdlXYZ * sizeof(XFLOAT)));
	DEBUG_HANDLE_ERROR(cudaMalloc( (void**) &mdlImag, mdlXYZ * sizeof(XFLOAT)));
#elif _HIP_ENABLED
	DEBUG_HANDLE_ERROR(hipMalloc( (void**) &mdlReal, mdlXYZ * sizeof(XFLOAT)));
	DEBUG_HANDLE_ERROR(hipMalloc( (void**) &mdlImag, mdlXYZ * sizeof(XFLOAT)));
#elif _SYCL_ENABLED
	devAcc = dev;
	mdlComplex = (XFLOAT*)devAcc->syclMalloc(2 * mdlXYZ * sizeof(XFLOAT), syclMallocType::device, "mdlComplex");
#else
	mdlComplex = NULL;
#endif
#endif
	return true;
}

#if defined(_CUDA_ENABLED) || defined(_HIP_ENABLED)
void AccProjector::initMdl(XFLOAT *real, XFLOAT *imag)
{
#if defined DEBUG_CUDA || defined DEBUG_HIP
	if (mdlXYZ == 0)
	{
        printf("DEBUG_ERROR: Model dimensions must be set with setMdlDim before call to setMdlData.");
		CRITICAL(ERR_MDLDIM);
	}
#if defined _CUDA_ENABLED || defined _HIP_ENABLED
	if (mdlReal == NULL)
	{
        printf("DEBUG_ERROR: initMdl called before call to setMdlData.");
		CRITICAL(ERR_MDLSET);
	}
#endif
#endif

#ifndef PROJECTOR_NO_TEXTURES
	#ifdef _CUDA_ENABLED
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
	#elif _HIP_ENABLED
	if(mdlZ!=0)  // 3D model
	{
		// -- make extents for automatic pitching (aligment) of allocated 3D arrays
		hipMemcpy3DParms copyParams = {0};
		copyParams.extent = make_hipExtent(mdlX, mdlY, mdlZ);
		copyParams.kind   = hipMemcpyHostToDevice;

		// -- Copy data
		copyParams.dstArray = *texArrayReal;
		copyParams.srcPtr   = make_hipPitchedPtr(real, mdlX * sizeof(XFLOAT), mdlY, mdlZ);
		DEBUG_HANDLE_ERROR(hipMemcpy3D(&copyParams));
		copyParams.dstArray = *texArrayImag;
		copyParams.srcPtr   = make_hipPitchedPtr(imag, mdlX * sizeof(XFLOAT), mdlY, mdlZ);
		DEBUG_HANDLE_ERROR(hipMemcpy3D(&copyParams));
	}
	else // 2D model
	{
		DEBUG_HANDLE_ERROR(hipMemcpy2D(texArrayReal2D, pitch2D, real, sizeof(XFLOAT) * mdlX, sizeof(XFLOAT) * mdlX, mdlY, hipMemcpyHostToDevice));
		DEBUG_HANDLE_ERROR(hipMemcpy2D(texArrayImag2D, pitch2D, imag, sizeof(XFLOAT) * mdlX, sizeof(XFLOAT) * mdlX, mdlY, hipMemcpyHostToDevice));
	}
	#endif
#else
#ifdef _CUDA_ENABLED
	DEBUG_HANDLE_ERROR(cudaMemcpy( mdlReal, real, mdlXYZ * sizeof(XFLOAT), cudaMemcpyHostToDevice));
	DEBUG_HANDLE_ERROR(cudaMemcpy( mdlImag, imag, mdlXYZ * sizeof(XFLOAT), cudaMemcpyHostToDevice));
#elif _HIP_ENABLED
	DEBUG_HANDLE_ERROR(hipMemcpy( mdlReal, real, mdlXYZ * sizeof(XFLOAT), hipMemcpyHostToDevice));
	DEBUG_HANDLE_ERROR(hipMemcpy( mdlImag, imag, mdlXYZ * sizeof(XFLOAT), hipMemcpyHostToDevice));
#endif
#endif
}
#endif

#ifdef ALTCPU
void AccProjector::initMdl(std::complex<XFLOAT> *data)
{
	mdlComplex = data;  // No copy needed - everyone shares the complex reference arrays
	externalFree = 1;   // This is shared memory freed outside the projector
}
#else
void AccProjector::initMdl(Complex *data)
{
#if defined(_CUDA_ENABLED) || defined(_HIP_ENABLED)
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
#elif _SYCL_ENABLED
	XFLOAT *tmpComplex = (XFLOAT*)devAcc->syclMalloc(2 * mdlXYZ * sizeof(XFLOAT), syclMallocType::host);
	if (nullptr == tmpComplex)
	{
		std::string str = "syclMalloc HOST error of size " + std::to_string(2*mdlXYZ * sizeof(XFLOAT)) + ".\n";
		ACC_PTR_DEBUG_FATAL(str.c_str());
		CRITICAL(RAMERR);
	}

	for (size_t i = 0; i < mdlXYZ; i++)
	{
		tmpComplex[2*i  ] = (XFLOAT) data[i].real;
		tmpComplex[2*i+1] = (XFLOAT) data[i].imag;
	}
	devAcc->syclMemcpy(mdlComplex, tmpComplex, 2 * mdlXYZ * sizeof(XFLOAT));
	devAcc->waitAll();
	devAcc->syclFree(tmpComplex);
#endif
}
#endif

void AccProjector::clear()
{
#if defined(_CUDA_ENABLED) || defined(_HIP_ENABLED)
	if (mdlReal != 0)
	{
#ifndef PROJECTOR_NO_TEXTURES
		#ifdef _CUDA_ENABLED
			cudaDestroyTextureObject(*mdlReal);
			cudaDestroyTextureObject(*mdlImag);
		#elif _HIP_ENABLED
			hipDestroyTextureObject(*mdlReal);
			hipDestroyTextureObject(*mdlImag);
		#endif
		delete mdlReal;
		delete mdlImag;

		if(mdlZ!=0) //3D case
		{
		#ifdef _CUDA_ENABLED
			cudaFreeArray(*texArrayReal);
			cudaFreeArray(*texArrayImag);
		#elif _HIP_ENABLED
			hipFreeArray(*texArrayReal);
			hipFreeArray(*texArrayImag);
		#endif
			delete texArrayReal;
			delete texArrayImag;
		}
		else //2D case
		{
		#ifdef _CUDA_ENABLED
			HANDLE_ERROR(cudaFree(texArrayReal2D));
			HANDLE_ERROR(cudaFree(texArrayImag2D));
		#elif _HIP_ENABLED
			HANDLE_ERROR(hipFree(texArrayReal2D));
			HANDLE_ERROR(hipFree(texArrayImag2D));
		#endif
		}

		texArrayReal = 0;
		texArrayImag = 0;
#else
	#ifdef _CUDA_ENABLED
		cudaFree(mdlReal);
		cudaFree(mdlImag);
	#elif _HIP_ENABLED
		hipFree(mdlReal);
		hipFree(mdlImag);
	#endif
#endif
		mdlReal = 0;
		mdlImag = 0;
	}
#elif _SYCL_ENABLED
	if (mdlComplex != NULL)
	{
		devAcc->waitAll();
		devAcc->syclFree(mdlComplex);
		mdlComplex = NULL;
	}
#else
	if ((mdlComplex != NULL) && (externalFree == 0))
	{
		delete [] mdlComplex;
		mdlComplex = NULL;
	}
#endif

	mdlX = 0;
	mdlY = 0;
	mdlZ = 0;
	mdlXYZ = 0;
	mdlInitY = 0;
	mdlInitZ = 0;
	mdlMaxR = 0;
	padding_factor = 0;
	allocaton_size = 0;
}
