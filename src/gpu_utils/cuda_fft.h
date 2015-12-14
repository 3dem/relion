#ifndef CUDA_FFT_H_
#define CUDA_FFT_H_

#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include <cuda_runtime.h>
#include <cufft.h>

#ifdef DEBUG_CUDA
#define HANDLE_CUFFT_ERROR( err ) (CufftHandleError( err, __FILE__, __LINE__ ))
#else
#define HANDLE_CUFFT_ERROR( err ) (err) //Do nothing
#endif
static void CufftHandleError( cufftResult err, const char *file, int line )
{
    if (err != CUFFT_SUCCESS)
    {
        fprintf(stderr, "Cufft error in file '%s' in line %i : %s.\n",
                __FILE__, __LINE__, "error" );
		raise(SIGSEGV);
    }
}

class CudaFFT
{
	bool planSet;
public:
	CudaGlobalPtr<cufftReal> reals;
	CudaGlobalPtr<cufftComplex> fouriers;
	cufftHandle cufftPlanForward, cufftPlanBackward;
	size_t xSize,ySize;

	CudaFFT(cudaStream_t stream, CudaCustomAllocator *allocator):
		reals(stream, allocator),
		fouriers(stream, allocator),
		cufftPlanForward(0),
		cufftPlanBackward(0),
		planSet(false),
		xSize(0), ySize(0)
	{};

	void setSize(size_t x, size_t y)
	{
		if (x == xSize && y == ySize)
			return;

		clear();

		xSize = x;
		ySize = y;

		reals.setSize(x*y);
		reals.device_alloc();

		fouriers.setSize(y*(x/2+1));
		fouriers.device_alloc();

		HANDLE_CUFFT_ERROR( cufftPlan2d(&cufftPlanForward,  x, y, CUFFT_R2C) );
		HANDLE_CUFFT_ERROR( cufftPlan2d(&cufftPlanBackward, x, y, CUFFT_C2R) );

		planSet = true;
	}

	void forward()
	{ HANDLE_CUFFT_ERROR( cufftExecR2C(cufftPlanForward, ~reals, ~fouriers) ); }

	void backward()
	{ HANDLE_CUFFT_ERROR( cufftExecC2R(cufftPlanBackward, ~fouriers, ~reals) ); }

	void backward(CudaGlobalPtr<cufftReal> &dst)
		{ HANDLE_CUFFT_ERROR( cufftExecC2R(cufftPlanBackward, ~fouriers, ~dst) ); }

	void clear()
	{
		if(planSet)
		{
			reals.free_if_set();
			fouriers.free_if_set();
			HANDLE_CUFFT_ERROR(cufftDestroy(cufftPlanForward));
			HANDLE_CUFFT_ERROR(cufftDestroy(cufftPlanBackward));
			planSet = false;
		}
	}

	~CudaFFT()
	{ clear(); }
};

#endif
