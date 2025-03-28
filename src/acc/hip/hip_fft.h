/* Portions of this code are under:
   Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
*/
#ifndef HIP_FFT_H_
#define HIP_FFT_H_

#include "src/acc/hip/hip_settings.h"
#include "src/acc/hip/hip_mem_utils.h"
#include <hip/hip_runtime.h>
#include <hipfft/hipfft.h>

#ifdef DEBUG_HIP
#define HANDLE_HIPFFT_ERROR( err ) (HipfftHandleError( err, __FILE__, __LINE__ ))
#else
#define HANDLE_HIPFFT_ERROR( err ) (err) //Do nothing
#endif
static void HipfftHandleError( hipfftResult err, const char *file, int line )
{
    if (err != HIPFFT_SUCCESS)
    {
        fprintf(stderr, "Hipfft error in file '%s' in line %i : %s.\n",
                __FILE__, __LINE__, "error" );
#ifdef DEBUG_HIP
		raise(SIGSEGV);
#else
		CRITICAL(ERRGPUKERN);
#endif
    }
}

class HipFFT
{
	bool planSet;
public:
#ifdef ACC_DOUBLE_PRECISION
	AccPtr<hipfftDoubleReal> reals;
	AccPtr<hipfftDoubleComplex> fouriers;
#else
	AccPtr<hipfftReal> reals;
	AccPtr<hipfftComplex> fouriers;
#endif
	hipfftHandle hipfftPlanForward, hipfftPlanBackward;
	int direction;
	int dimension, idist, odist, istride, ostride;
	int inembed[3];
	int onembed[3];
	size_t xSize,ySize,zSize,xFSize,yFSize,zFSize;
	std::vector< int >  batchSize;
	HipCustomAllocator *CFallocator;
	int batchSpace, batchIters, reqN;

	HipFFT(hipStream_t stream, HipCustomAllocator *allocator, int transformDimension = 2):
		reals(stream, allocator),
		fouriers(stream, allocator),
		hipfftPlanForward(0),
		hipfftPlanBackward(0),
		direction(0),
		dimension((int)transformDimension),
	    idist(0),
	    odist(0),
	    istride(1),
	    ostride(1),
		inembed(),
		onembed(),
		planSet(false),
		xSize(0), ySize(0), zSize(0),
		xFSize(0), yFSize(0), zFSize(0),
		batchSize(1,1),
		reqN(1),
		CFallocator(allocator)
	{};

	void setAllocator(HipCustomAllocator *allocator)
	{
		reals.setAllocator(allocator);
		fouriers.setAllocator(allocator);
		CFallocator = allocator;
	}

	size_t estimate(int batch)
	{
		size_t needed(0);

	    size_t biggness(0);

#ifdef ACC_DOUBLE_PRECISION
	    if(direction<=0)
	    {
			HANDLE_HIPFFT_ERROR( hipfftEstimateMany(dimension, inembed, inembed, istride, idist, onembed, ostride, odist, HIPFFT_D2Z, batch, &biggness));
			needed += biggness;
	    }
		if(direction>=0)
		{
			HANDLE_HIPFFT_ERROR( hipfftEstimateMany(dimension, inembed, onembed, ostride, odist, inembed, istride, idist, HIPFFT_Z2D, batch, &biggness));
			needed += biggness;
		}
#else
		if(direction<=0)
		{
			HANDLE_HIPFFT_ERROR( hipfftEstimateMany(dimension, inembed, inembed, istride, idist, onembed, ostride, odist, HIPFFT_R2C, batch, &biggness));
			needed += biggness;
		}
		if(direction>=0)
		{
			HANDLE_HIPFFT_ERROR( hipfftEstimateMany(dimension, inembed, onembed, ostride, odist, inembed, istride, idist, HIPFFT_C2R, batch, &biggness));
			needed += biggness;
		}
#endif
		size_t res = needed + (size_t)odist*(size_t)batch*sizeof(XFLOAT)*(size_t)2 + (size_t)idist*(size_t)batch*sizeof(XFLOAT);

		return res;
	}

	void setSize(size_t x, size_t y, size_t z, int batch = 1, int setDirection = 0)
	{

		/* Optional direction input restricts transformer to
		 * forwards or backwards tranformation only,
		 * which reduces memory requirements, especially
		 * for large batches of simulatanous transforms.
		 *
		 * FFTW_FORWARDS  === -1
		 * FFTW_BACKWARDS === +1
		 *
		 * The default direction is 0 === forwards AND backwards
		 */

		int checkDim;
		if(z>1)
			checkDim=3;
		else if(y>1)
			checkDim=2;
		else
			checkDim=1;
		if(checkDim != dimension)
			CRITICAL(ERRHIPFFTDIM);

		if( !( (setDirection==-1)||(setDirection==0)||(setDirection==1) ) )
		{
			std::cerr << "*ERROR : Setting a hip transformer direction to non-defined value" << std::endl;
			CRITICAL(ERRHIPFFTDIR);
		}

		direction = setDirection;

		if (x == xSize && y == ySize && z == zSize && batch == reqN && planSet)
			return;

		clear();

		batchSize.resize(1);
		batchSize[0] = batch;
		reqN = batch;

		xSize = x;
		ySize = y;
		zSize = z;

		xFSize = x/2 + 1;
		yFSize = y;
		zFSize = z;

	    idist = zSize*ySize*xSize;
	    odist = zSize*ySize*(xSize/2+1);
	    istride = 1;
	    ostride = 1;

	    if(dimension==3)
	    {
	    	inembed[0] =  zSize;
			inembed[1] =  ySize;
			inembed[2] =  xSize;
			onembed[0] =  zFSize;
			onembed[1] =  yFSize;
			onembed[2] =  xFSize;
	    }
	    else if(dimension==2)
	    {
			inembed[0] =  ySize;
			inembed[1] =  xSize;
			onembed[0] =  yFSize;
			onembed[1] =  xFSize;
	    }
	    else
	    {
			inembed[0] =  xSize;
			onembed[0] =  xFSize;
	    }

		size_t needed, avail, total;
		needed = estimate(batchSize[0]);
		DEBUG_HANDLE_ERROR(hipMemGetInfo( &avail, &total ));

//		std::cout << std::endl << "needed = ";
//		printf("%15zu\n", needed);
//		std::cout << "avail  = ";
//		printf("%15zu\n", avail);

		// Check if there is enough memory
		//
		//    --- TO HOLD TEMPORARY DATA DURING TRANSFORMS ---
		//
		// If there isn't, find how many there ARE space for and loop through them in batches.

		if(needed>avail)
		{
			batchIters = 2;
			batchSpace = CEIL((double) batch / (double)batchIters);
			needed = estimate(batchSpace);

			while(needed>avail && batchSpace>1)
			{
				batchIters++;
				batchSpace = CEIL((double) batch / (double)batchIters);
				needed = estimate(batchSpace);
			}

			if(batchIters>1)
			{
				batchIters = (int)((float)batchIters*1.1 + 1);
				batchSpace = CEIL((double) batch / (double)batchIters);
				needed = estimate(batchSpace);
			}

			batchSize.assign(batchIters,batchSpace); // specify batchIters of batches, each with batchSpace orientations
			batchSize[batchIters-1] = batchSpace - (batchSpace*batchIters - batch); // set last to care for remainder.

			if(needed>avail)
				CRITICAL(ERRFFTMEMLIM);

//			std::cerr << std::endl << "NOTE: Having to use " << batchIters << " batches of orientations ";
//			std::cerr << "to achieve the total requested " << batch << " orientations" << std::endl;
//			std::cerr << "( this could affect performance, consider using " << std::endl;
//			std::cerr << "\t higher --ang" << std::endl;
//			std::cerr << "\t harder --shrink" << std::endl;
//			std::cerr << "\t higher --lopass with --shrink 0" << std::endl;

		}
		else
		{
			batchIters = 1;
			batchSpace = batch;
		}

		reals.setSize(idist*batchSize[0]);
		reals.deviceAlloc();
		reals.hostAlloc();

		fouriers.setSize(odist*batchSize[0]);
		fouriers.deviceAlloc();
		fouriers.hostAlloc();

//		DEBUG_HANDLE_ERROR(hipMemGetInfo( &avail, &total ));
//		needed = estimate(batchSize[0], fudge);

//		std::cout << "after alloc: " << std::endl << std::endl << "needed = ";
//		printf("%15li\n", needed);
//		std::cout << "avail  = ";
//		printf("%15li\n", avail);

#ifdef ACC_DOUBLE_PRECISION
	    if(direction<=0)
	    {
	    	HANDLE_HIPFFT_ERROR( hipfftPlanMany(&hipfftPlanForward,  dimension, inembed, inembed, istride, idist, onembed, ostride, odist, HIPFFT_D2Z, batchSize[0]));
	   		HANDLE_HIPFFT_ERROR( hipfftSetStream(hipfftPlanForward, fouriers.getStream()));
	    }
	    if(direction>=0)
	    {
	    	HANDLE_HIPFFT_ERROR( hipfftPlanMany(&hipfftPlanBackward, dimension, inembed, onembed, ostride, odist, inembed, istride, idist, HIPFFT_Z2D, batchSize[0]));
			HANDLE_HIPFFT_ERROR( hipfftSetStream(hipfftPlanBackward, reals.getStream()));
	    }
		planSet = true;
	}

	void forward()
	{ HANDLE_HIPFFT_ERROR( hipfftExecD2Z(hipfftPlanForward, ~reals, ~fouriers) ); }

	void backward()
	{ HANDLE_HIPFFT_ERROR( hipfftExecZ2D(hipfftPlanBackward, ~fouriers, ~reals) ); }

	void backward(AccPtr<hipfftDoubleReal> &dst)
		{ HANDLE_HIPFFT_ERROR( hipfftExecZ2D(hipfftPlanBackward, ~fouriers, ~dst) ); }
#else
	 	if(direction<=0)
	 	{
	 		HANDLE_HIPFFT_ERROR( hipfftPlanMany(&hipfftPlanForward,  dimension, inembed, inembed, istride, idist, onembed, ostride, odist, HIPFFT_R2C, batchSize[0]));
	 		HANDLE_HIPFFT_ERROR( hipfftSetStream(hipfftPlanForward, fouriers.getStream()));
	 	}
	 	if(direction>=0)
	 	{
	 		HANDLE_HIPFFT_ERROR( hipfftPlanMany(&hipfftPlanBackward, dimension, inembed, onembed, ostride, odist, inembed, istride, idist, HIPFFT_C2R, batchSize[0]));
	 		HANDLE_HIPFFT_ERROR( hipfftSetStream(hipfftPlanBackward, reals.getStream()));
	 	}
		planSet = true;
	}

	void forward()
	{
		if(direction==1)
		{
			std::cout << "trying to execute a forward plan for a hipFFT transformer which is backwards-only" << std::endl;
			CRITICAL(ERRHIPFFTDIRF);
		}
		HANDLE_HIPFFT_ERROR( hipfftExecR2C(hipfftPlanForward, ~reals, ~fouriers) );
	}

	void backward()
	{
		if(direction==-1)
		{
			std::cout << "trying to execute a backwards plan for a hipFFT transformer which is forwards-only" << std::endl;
			CRITICAL(ERRHIPFFTDIRR);
		}
		HANDLE_HIPFFT_ERROR( hipfftExecC2R(hipfftPlanBackward, ~fouriers, ~reals) );
	}

#endif

	void clear()
	{
		if(planSet)
		{
			reals.freeIfSet();
			fouriers.freeIfSet();
			if(direction<=0)
				HANDLE_HIPFFT_ERROR(hipfftDestroy(hipfftPlanForward));
			if(direction>=0)
				HANDLE_HIPFFT_ERROR(hipfftDestroy(hipfftPlanBackward));
			planSet = false;
		}
	}

	~HipFFT()
	{clear();}
};

#endif //HIP_FFT_H_
