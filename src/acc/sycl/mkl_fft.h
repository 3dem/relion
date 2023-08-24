#ifndef MKL_FFT_H_
#define MKL_FFT_H_

#include <vector>
#include <mutex>
#include <fftw3.h>
#include "src/complex.h"

extern std::mutex fft_mutex;

class MklFFT
{
	bool planSet;
public:
	AccPtr<XFLOAT>      reals;
	AccPtr<ACCCOMPLEX> fouriers;

	int    direction;
	int    dimension;
	size_t xSize,ySize,zSize,xFSize,yFSize,zFSize;
	
#ifdef ACC_DOUBLE_PRECISION
    /* fftw Forward plan */
    fftw_plan fPlanForward;

    /* fftw Backward plan */
    fftw_plan fPlanBackward;

#else
    /* fftw Forward plan */
    fftwf_plan fPlanForward;

    /* fftw Backward plan */
    fftwf_plan fPlanBackward;
#endif

	MklFFT(int transformDimension = 2):
	    direction {0}, 
	    dimension {transformDimension},
		planSet {false},
	    xSize {0}, ySize{0}, zSize{0}
	{
        fPlanForward = fPlanBackward = NULL;
	};

	void setSize(size_t x, size_t y, size_t z, int setDirection = 0)
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
			REPORT_ERROR("You are trying to change the dimesion of a MklFFT transformer, which is not allowed");

		if( !( (setDirection==-1)||(setDirection==0)||(setDirection==1) ) )
		{
			std::cerr << "*ERROR : Setting a MklFFT transformer direction to non-defined value" << std::endl;
			return;
		}

		direction = setDirection;

		if( x == xSize && y == ySize && z == zSize && planSet)
			return;
		
		clear();

		xSize = x;
		ySize = y;
		zSize = z;

		xFSize = x/2 + 1;
		yFSize = y;
		zFSize = z;

		if ((xSize * ySize * zSize)==0)
			ACC_PTR_DEBUG_FATAL("Reals array resized to size zero.\n");
//		reals.resizeHostCopy(xSize * ySize * zSize);
		reals.freeHostIfSet();
		reals.setSize(xSize * ySize * zSize);
		reals.hostAlloc();
		
		if ((xFSize * yFSize * zFSize)==0)
			ACC_PTR_DEBUG_FATAL("Fouriers array resized to size zero.\n");
//		fouriers.resizeHostCopy(xFSize * yFSize * zFSize);
		fouriers.freeHostIfSet();
		fouriers.setSize(xFSize * yFSize * zFSize);
		fouriers.hostAlloc();

		int N[3];  
		if(dimension == 1)
			N[0] = xSize;
		else  if(dimension == 2){
			N[0] = ySize;  N[1] = xSize;	    
		}
		else {
			N[0] = zSize;  N[1] = ySize; N[2] = xSize;	    
		}

		{
			std::lock_guard<std::mutex> lock(fft_mutex);
#ifdef ACC_DOUBLE_PRECISION
			fPlanForward = fftw_plan_dft_r2c(dimension, N,  reals(),
										 (fftw_complex*) fouriers(), FFTW_ESTIMATE);
			fPlanBackward = fftw_plan_dft_c2r(dimension, N,
										  (fftw_complex*) fouriers(),  reals(),
										  FFTW_ESTIMATE);

#else
			fPlanForward = fftwf_plan_dft_r2c(dimension, N, reals(),
										 (fftwf_complex*) fouriers(), FFTW_ESTIMATE);
			fPlanBackward = fftwf_plan_dft_c2r(dimension, N,
										  (fftwf_complex*) fouriers(), reals(), FFTW_ESTIMATE);
#endif
			planSet = true;
		}
	}

	void forward()
	{
		if(direction==1)
		{
			std::cout << "trying to execute a forward plan for a MKL FFT transformer which is backwards-only" << std::endl;
			return;
		}
#ifdef ACC_DOUBLE_PRECISION
		fftw_execute_dft_r2c(fPlanForward, reals(), (fftw_complex*) fouriers());
#else
		fftwf_execute_dft_r2c(fPlanForward, reals(),  (fftwf_complex*) fouriers());
#endif
     
	}

	void backward()
	{
	    if(direction==-1)
	    {
	        std::cout << "trying to execute a backwards plan for a MKL FFT transformer which is forwards-only" << std::endl;
	        return;
	    }	     

#ifdef ACC_DOUBLE_PRECISION
        fftw_execute_dft_c2r(fPlanBackward, (fftw_complex*) fouriers(), reals());
#else
        fftwf_execute_dft_c2r(fPlanBackward, (fftwf_complex*) fouriers(), reals());
#endif	
	}

	void clear()
	{
		reals.freeIfSet();
		fouriers.freeIfSet();
	
		if (planSet)
		{
			std::lock_guard<std::mutex> lock(fft_mutex);        
#ifdef ACC_DOUBLE_PRECISION
			fftw_destroy_plan(fPlanForward);
			fftw_destroy_plan(fPlanBackward);
#else
			fftwf_destroy_plan(fPlanForward);
			fftwf_destroy_plan(fPlanBackward);
#endif          
			fPlanForward = fPlanBackward = NULL;
			planSet = false;
		}
	}

	~MklFFT()
	{ clear(); }
};

#endif
