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
	std::vector< int >  batchSize;
	CudaCustomAllocator *CFallocator;
	int psiSpace, psiIters;

	CudaFFT(cudaStream_t stream, CudaCustomAllocator *allocator):
		reals(stream, allocator),
		fouriers(stream, allocator),
		cufftPlanForward(0),
		cufftPlanBackward(0),
		planSet(false),
		xSize(0), ySize(0),
		xFSize(0), yFSize(0),
		batchSize(1,1),
		CFallocator(allocator)
	{};

	long int estimate(int batch, float fudge = 1.0)
	{
		size_t needed;

	    int idist = ySize*xSize;
	    int odist = ySize*(xSize/2+1);

	    int inembed[] = {ySize, xSize};
	    int onembed[] = {ySize, xSize/2+1};

	    int istride = 1;
	    int ostride = 1;

	    int nR[2] = {ySize, xSize};
	    size_t biggness;

#ifdef CUDA_DOUBLE_PRECISION
		HANDLE_CUFFT_ERROR( cufftEstimateMany(2, nR, inembed, istride, idist, onembed, ostride, odist, CUFFT_D2Z, batch, &biggness));
		needed = biggness;
		HANDLE_CUFFT_ERROR( cufftEstimateMany(2, nR, onembed, ostride, odist, inembed, istride, idist, CUFFT_Z2D, batch, &biggness));
		needed += biggness;
#else
		HANDLE_CUFFT_ERROR( cufftEstimateMany(2, nR, inembed, istride, idist, onembed, ostride, odist, CUFFT_R2C, batch, &biggness));
		needed = biggness;
		HANDLE_CUFFT_ERROR( cufftEstimateMany(2, nR, onembed, ostride, odist, inembed, istride, idist, CUFFT_C2R, batch, &biggness));
		needed += biggness;
#endif
		return (long int)((float)needed*fudge);
	}

	void setSize(size_t x, size_t y, int batch = 1)
	{
		if (x == xSize && y == ySize && batch == batchSize[0] && planSet)
			return;

		clear();

		batchSize.resize(1);
		batchSize[0] = batch;

		xSize = x;
		ySize = y;
		xFSize = x/2 + 1;
		yFSize = y;

		float fudge = 2.0;

		size_t needed, avail, total;
		needed = estimate(batchSize[0],fudge);
		DEBUG_HANDLE_ERROR(cudaMemGetInfo( &avail, &total ));

		double memFrac = (double)needed / (double)avail;

//		std::cout << std::endl << "needed = ";
//		printf("%15li\n", needed);
//		std::cout << "avail  = ";
//		printf("%15li\n", avail);
//		std::cout << "memFrac  = ";
//		printf("%15f\n", memFrac);


		// set the size of the real-array to hold the RESULT of transforms ALL batches.
		reals.setSize(x*y*batch);
		reals.device_alloc();
		reals.host_alloc();

		// set the size of the fouriers-array to hold the RESULT of transforms ALL batches.
		fouriers.setSize(y*(x/2+1)*batch);
		fouriers.device_alloc();
		fouriers.host_alloc();

		// Check if there is enough memory
		//
		//    --- TO HOLD TEMPORARY DATA DURING TRANSFORMS ---
		//
		// If there isn't, find how many there ARE space for and loop through them in batches.

		// batch fudge-factor to avoid running out of temp-space for transform plan

		if(memFrac>1)
		{
			psiIters = CEIL(memFrac);
			psiSpace = CEIL((double) batch / (double)psiIters);
			needed = estimate(psiSpace,fudge);

			while(needed>avail && psiSpace>1)
			{
				psiIters++;
				psiSpace = CEIL((double) batch / (double)psiIters);
				needed = estimate(psiSpace,fudge);
			}

			batchSize.assign(psiIters,psiSpace); // specify psiIters of batches, each with psiSpace orientations
			batchSize[psiIters-1] = psiSpace - (psiSpace*psiIters - batch); // set last to care for remainder.

			if(needed>avail)
				REPORT_ERROR("Not enough memory for even a single orientation.");

			std::cerr << std::endl << "NOTE: Having to use " << psiIters << " batches of orientations ";
			std::cerr << "to achieve the total requested " << batch << " orientations" << std::endl;
//			std::cerr << "( this could affect performance, consider using " << std::endl;
//			std::cerr << "\t higher --ang" << std::endl;
//			std::cerr << "\t harder --shrink" << std::endl;
//			std::cerr << "\t higher --lopass with --shrink 0" << std::endl;

		}
		else
		{
			psiIters = 1;
			psiSpace = batch;
		}

		DEBUG_HANDLE_ERROR(cudaMemGetInfo( &avail, &total ));
		needed = estimate(batchSize[0], fudge);

//		std::cout << "after alloc: " << std::endl << std::endl << "needed = ";
//		printf("%15li\n", needed);
//		std::cout << "avail  = ";
//		printf("%15li\n", avail);

	    int idist = y*x;
	    int odist = y*(x/2+1);

	    int inembed[] = {y, x};
	    int onembed[] = {y, x/2+1};

	    int istride = 1;
	    int ostride = 1;

	    int nR[2] = {y, x};
//	    int nC[2] = {y, x/2 +1};
#ifdef CUDA_DOUBLE_PRECISION
		HANDLE_CUFFT_ERROR( cufftPlanMany(&cufftPlanForward,  2, nR, inembed, istride, idist, onembed, ostride, odist, CUFFT_D2Z, batchSize[0]));
		HANDLE_CUFFT_ERROR( cufftPlanMany(&cufftPlanBackward, 2, nR, onembed, ostride, odist, inembed, istride, idist, CUFFT_Z2D, batchSize[0]));
		HANDLE_CUFFT_ERROR( cufftSetStream(cufftPlanForward, fouriers.getStream()));
		HANDLE_CUFFT_ERROR( cufftSetStream(cufftPlanBackward, reals.getStream()));
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
	    HANDLE_CUFFT_ERROR( cufftPlanMany(&cufftPlanForward,  2, nR, inembed, istride, idist, onembed, ostride, odist, CUFFT_R2C, batchSize[0]));
		HANDLE_CUFFT_ERROR( cufftPlanMany(&cufftPlanBackward, 2, nR, onembed, ostride, odist, inembed, istride, idist, CUFFT_C2R, batchSize[0]));
		HANDLE_CUFFT_ERROR( cufftSetStream(cufftPlanForward, fouriers.getStream()));
		HANDLE_CUFFT_ERROR( cufftSetStream(cufftPlanBackward, reals.getStream()));
//		HANDLE_CUFFT_ERROR( cufftPlanMany(&cufftPlanForward,   2, nR, 0,0,0,0,0,0, CUFFT_R2C, batchSize));
//		HANDLE_CUFFT_ERROR( cufftPlanMany(&cufftPlanBackward,  2, nC, 0,0,0,0,0,0, CUFFT_C2R, batchSize));

		planSet = true;
	}

	void forward()
	{
		if(psiIters>1)
		{
			long int Fstride =xFSize*yFSize;
			long int Rstride =xSize*ySize;
			long int Fpos = 0;
			long int Rpos = 0;
			for (int psiIter = 0; psiIter < psiIters; psiIter++) // psi-batches for possible memory-limits
			{
				HANDLE_CUFFT_ERROR( cufftExecR2C(cufftPlanForward, &reals(Rpos), &fouriers(Fpos)) );
				Fpos += Fstride*batchSize[psiIter];
				Rpos += Rstride*batchSize[psiIter];
			}
		}
		else
		{
			HANDLE_CUFFT_ERROR( cufftExecR2C(cufftPlanForward, ~reals, ~fouriers) );
		}
	}

	void backward()
	{
		if(psiIters>1)
		{
			long int Fstride =xFSize*yFSize;
			long int Rstride =xSize*ySize;
			long int Fpos = 0;
			long int Rpos = 0;
			for (int psiIter = 0; psiIter < psiIters; psiIter++) // psi-batches for possible memory-limits
			{
				HANDLE_CUFFT_ERROR( cufftExecC2R(cufftPlanBackward, &fouriers(Fpos), &reals(Rpos)) );
				Fpos += Fstride*batchSize[psiIter];
				Rpos += Rstride*batchSize[psiIter];
			}
		}
		else
		{
			HANDLE_CUFFT_ERROR( cufftExecC2R(cufftPlanBackward, ~fouriers, ~reals) );
		}
	}

	void backward(CudaGlobalPtr<cufftComplex> &src)
		{ HANDLE_CUFFT_ERROR( cufftExecC2R(cufftPlanBackward, ~src, ~reals) ); }

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
