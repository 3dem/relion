#include <math.h>
#include <matrix.h>
#include <mex.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mwSize w, h, d;
	
	if (nrhs < 1) 
	{
		mexErrMsgIdAndTxt("Dynamite::allocateComplexVolumes", "Usage: allocateComplexVolumes(w[,h,d])");
	}
	
	if (nlhs != 2) 
	{
		mexErrMsgTxt("Two output arguments required: data and weight.");
	}
	
	w = (mwSize) mxGetScalar(prhs[0]);
	
	if (nrhs > 1) h = (mwSize) mxGetScalar(prhs[1]);
	else h = w;
	
	if (nrhs > 2) d = (mwSize) mxGetScalar(prhs[2]);
	else d = h;
	
	const mwSize dims[] = {(mwSize)(w/2 + 1), h, d, 2};
	
	plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
}
