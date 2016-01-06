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
#ifdef CUDA_DOUBLE_PRECISION
	CudaGlobalPtr<cufftDoubleReal> reals;
	CudaGlobalPtr<cufftDoubleComplex> fouriers;
#else
	CudaGlobalPtr<cufftReal> reals;
	CudaGlobalPtr<cufftComplex> fouriers;
#endif
	cufftHandle cufftPlanForward, cufftPlanBackward;
	size_t xSize,ySize,xFSize,yFSize;
	int batchSize;

	CudaFFT(cudaStream_t stream, CudaCustomAllocator *allocator):
		reals(stream, allocator),
		fouriers(stream, allocator),
		cufftPlanForward(0),
		cufftPlanBackward(0),
		planSet(false),
		xSize(0), ySize(0),
		xFSize(0), yFSize(0),
		batchSize(1)
	{};

	void setSize(size_t x, size_t y, int batch = 1)
	{
		if (x == xSize && y == ySize && batch == batchSize)
			return;

		clear();

		batchSize = batch;

		xSize = x;
		ySize = y;
		xFSize = x/2 + 1;
		yFSize = y;

		reals.setSize(x*y*batchSize);
		reals.device_alloc();
		reals.host_alloc();

		fouriers.setSize(y*(x/2+1)*batchSize);
		fouriers.device_alloc();
		fouriers.host_alloc();

	    int idist = ySize*xSize;
	    int odist = ySize*(xSize/2+1);

	    int inembed[] = {ySize, xSize};
	    int onembed[] = {ySize, xSize/2+1};

	    int istride = 1;
	    int ostride = 1;

	    int nR[2] = {ySize, xSize};
//	    int nC[2] = {ySize, xSize/2 +1};
#ifdef CUDA_DOUBLE_PRECISION
		HANDLE_CUFFT_ERROR( cufftPlanMany(&cufftPlanForward,  2, nR, inembed, istride, idist, onembed, ostride, odist, CUFFT_D2Z, batchSize));
		HANDLE_CUFFT_ERROR( cufftPlanMany(&cufftPlanBackward, 2, nR, onembed, ostride, odist, inembed, istride, idist, CUFFT_Z2D, batchSize));
//		HANDLE_CUFFT_ERROR( cufftPlan2d(&cufftPlanForward,  x, y, CUFFT_D2Z) );
//		HANDLE_CUFFT_ERROR( cufftPlan2d(&cufftPlanBackward, x, y, CUFFT_Z2D) );

		planSet = true;
	}

	void forward()
	{ HANDLE_CUFFT_ERROR( cufftExecD2Z(cufftPlanForward, ~reals, ~fouriers) ); }

	void backward()
	{ HANDLE_CUFFT_ERROR( cufftExecZ2D(cufftPlanBackward, ~fouriers, ~reals) ); }

	void backward(CudaGlobalPtr<cufftDoubleReal> &dst)
		{ HANDLE_CUFFT_ERROR( cufftExecZ2D(cufftPlanBackward, ~fouriers, ~dst) ); }
#else
//		HANDLE_CUFFT_ERROR( cufftPlan2d(&cufftPlanForward,  x, y, CUFFT_R2C) );
//		HANDLE_CUFFT_ERROR( cufftPlan2d(&cufftPlanBackward, x, y, CUFFT_C2R) );

//		size_t biggness;
//	    HANDLE_CUFFT_ERROR( cufftEstimateMany(2, nR, inembed, istride, idist, onembed, ostride, odist, CUFFT_R2C, batchSize, &biggness));
//		std::cerr  << std::endl << " estimated size for transform = " << biggness << std::endl << std::endl;

	    HANDLE_CUFFT_ERROR( cufftPlanMany(&cufftPlanForward,  2, nR, inembed, istride, idist, onembed, ostride, odist, CUFFT_R2C, batchSize));
		HANDLE_CUFFT_ERROR( cufftPlanMany(&cufftPlanBackward, 2, nR, onembed, ostride, odist, inembed, istride, idist, CUFFT_C2R, batchSize));
//		HANDLE_CUFFT_ERROR( cufftPlanMany(&cufftPlanForward,   2, nR, 0,0,0,0,0,0, CUFFT_R2C, batchSize));
//		HANDLE_CUFFT_ERROR( cufftPlanMany(&cufftPlanBackward,  2, nC, 0,0,0,0,0,0, CUFFT_C2R, batchSize));

		planSet = true;
	}

	void forward()
	{ HANDLE_CUFFT_ERROR( cufftExecR2C(cufftPlanForward, ~reals, ~fouriers) ); }

	void backward()
	{ HANDLE_CUFFT_ERROR( cufftExecC2R(cufftPlanBackward, ~fouriers, ~reals) ); }

	void backward(CudaGlobalPtr<cufftReal> &dst)
		{ HANDLE_CUFFT_ERROR( cufftExecC2R(cufftPlanBackward, ~fouriers, ~dst) ); }
#endif
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
