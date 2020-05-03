#include <math.h>
#include <matrix.h>
#include <mex.h>

#include <src/jaz/gravis/t4Matrix.h>
#include <img_proc/image.h>
#include <src/jaz/tomo/projection_IO.h>
#include <src/jaz/math/fft.h>


void crashOnArgument(int argumentID);

gravis::d4Matrix anglesToMatrix(
	double tdrot, double tilt, double narot, 
	double cx, double cy,
	double w, double h, double d	);


/*
	IN:
		[0] 	current_data: 3D_array<double> 
		[1] 	current_weight: 3D_array<double> 
		[2] 	Wiener_offset: double
		[3] 	#threads: int   (optional)
		
	OUT:
		[0]	divided_map: 3D_array<double> 
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const int dataInID = 0;
	const int psfInID = 1;
	const int WienerID = 2;
	const int threadNumID = 3;
	
	int num_threads = 1;
	double WienerOffset = 1.0;
	
	if (nrhs > 3) 
	{
		num_threads = mxGetDoubles(prhs[threadNumID])[0];
	}
	
	if (nrhs > 2) 
	{
		WienerOffset = mxGetDoubles(prhs[WienerID])[0];
	}
	
	if (nrhs < 2)
	{
		mexErrMsgIdAndTxt("Dynamite::deconvolveByPsf", 
			"Usage: [divided_map] = deconvolveByPsf(data, PSF, Wiener_offset, #threads)");
	}
	
	if (nlhs != 1) 
	{
		mexErrMsgTxt("Dynamite::deconvolveByPsf: one output argument required");
	}
	
	// image stack:
	
	if (mxGetNumberOfDimensions(prhs[dataInID]) != 3) crashOnArgument(dataInID);
	if (mxGetNumberOfDimensions(prhs[psfInID]) != 3) crashOnArgument(psfInID);
	
	const mwSize* dims3D = mxGetDimensions(prhs[dataInID]);
	const mwSize w = dims3D[0];
	const mwSize h = dims3D[1];
	const mwSize d = dims3D[2];
	
	plhs[0] = mxDuplicateArray(prhs[dataInID]);
	
	Image<double> dataImg(w,h,d,mxGetDoubles(prhs[dataInID]));
	Image<double> psfImg(w,h,d,mxGetDoubles(prhs[psfInID]));
	Image<double> outImg(w,h,d,mxGetDoubles(plhs[0]));
	
	VectorImage<double> dataImgCp(dataImg);
	VectorImage<double> psfImgCp(psfImg);
	
	const int wh = w/2 + 1;
	
	VectorImage<tComplex<double>> data_FS, psf_FS, out_FS(wh,h,d);
	
	FFT::FourierTransform(dataImgCp, data_FS, FFT::Both);
	FFT::FourierTransform(psfImgCp, psf_FS, FFT::Both);
	
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		out_FS(x,y,z) = data_FS(x,y,z) / (psf_FS(x,y,z) + WienerOffset);
	}
	
	FFT::inverseFourierTransform(out_FS, dataImgCp, FFT::Both);
	
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < w; x++)
	{
		outImg(x,y,z) = dataImgCp(x,y,z);
	}
}

void crashOnArgument(int argumentID)
{
	mexPrintf("bad argument #%d\n", argumentID + 1);
	mexErrMsgTxt("aborting backprojection");
}
