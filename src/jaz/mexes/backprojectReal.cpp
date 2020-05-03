#include <math.h>
#include <matrix.h>
#include <mex.h>

#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomo/projection_IO.h>


void crashOnArgument(int argumentID);

gravis::d4Matrix anglesToMatrix(
	double tdrot, double tilt, double narot, 
	double cx, double cy,
	double w, double h, double d	);


/*
	IN:
		[0] 	stack: 3D_array<float> 
		[1]	particle_positions: 3D_array<double> [2 x #F x #P]
		[2]	particle_angles: 3D_array<double> [3 x #F x #P]
		[3] 	current_data: 3D_array<double> 
		[4] 	current_weight: 3D_array<double> 
		[5] 	#threads: int   (optional)
		
	OUT:
		[0]	updated_data: 3D_array<double> 
		[1] 	updated_weight: 3D_array<double> 
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	const int stackID = 0;
	const int posID = 1;
	const int angID = 2;
	const int dataInID = 3;
	const int psfInID = 4;
	const int threadNumID = 5;
	
	int num_threads;
	
	if (nrhs == 6) 
	{
		num_threads = mxGetInt32s(prhs[threadNumID])[0];
	}
	else if (nrhs == 5) 
	{
		num_threads = 1;
	}
	else
	{
		mexErrMsgIdAndTxt("Dynamite::backprojectReal", 
			"Usage: [data, weight] = backprojectReal(stack, pos, ang, data, PSF, #threads)");
	}
	
	if (nlhs != 2) 
	{
		mexErrMsgTxt("Two output arguments required: data and weight.");
	}
	
	// image stack:
	
	if (mxGetNumberOfDimensions(prhs[stackID]) != 3) crashOnArgument(stackID);
	
	const mwSize* dims2D = mxGetDimensions(prhs[stackID]);
	const mwSize w2D = dims2D[0];
	const mwSize h2D = dims2D[1];
	const mwSize fc = dims2D[2];
	
	// positions and angles:
	
	if (mxGetNumberOfDimensions(prhs[posID]) < 2) crashOnArgument(posID);
	if (mxGetNumberOfDimensions(prhs[angID]) < 2) crashOnArgument(angID);
	
	const mwSize* dimsPos = mxGetDimensions(prhs[posID]);
	const mwSize* dimsAng = mxGetDimensions(prhs[posID]);
	
	if (dimsPos[0] != 2) mexErrMsgTxt("particle positions have to be 2D - aborting backprojection");
	if (dimsAng[0] != 3) mexErrMsgTxt("particle angles have to be 3D - aborting backprojection");
	
	if (dimsPos[1] != fc) mexErrMsgTxt("frame number mismatch between positions and image stack - aborting backprojection");
	if (dimsAng[1] != fc) mexErrMsgTxt("frame number mismatch between angles and image stack - aborting backprojection");
	
	const int pc = dimsPos[2];
	if (dimsAng[2] != pc) mexErrMsgTxt("particle number mismatch between angles and positions - aborting backprojection");
	
	const double* posPtr = mxGetDoubles(prhs[posID]);
	const double* angPtr = mxGetDoubles(prhs[angID]);
	
	// 3D volumes:
	
	if (mxGetNumberOfDimensions(prhs[dataInID]) != 3) crashOnArgument(dataInID);
	if (mxGetNumberOfDimensions(prhs[psfInID]) != 3) crashOnArgument(psfInID);
	
	const mwSize* dims3D = mxGetDimensions(prhs[dataInID]);
	const mwSize w = dims3D[0];
	const mwSize h = dims3D[1];
	const mwSize d = dims3D[2];
		
	// prepare projections
	
	std::vector<std::vector<gravis::d4Matrix>> projections(pc);
	
	for (int p = 0; p < pc; p++)
	{
		projections[p] = std::vector<gravis::d4Matrix>(fc);
		
		for (int f = 0; f < fc; f++)
		{
			const int pf = p*fc + f;
			
			projections[p][f] = anglesToMatrix(
				angPtr[3*pf], angPtr[3*pf + 1], angPtr[3*pf + 2],
				posPtr[2*pf], posPtr[2*pf + 1],
				w, h, d);
		}
	}
	
	plhs[0] = mxDuplicateArray(prhs[dataInID]);
	plhs[1] = mxDuplicateArray(prhs[psfInID]);
		
	float* stackPtr = mxGetSingles(prhs[stackID]);
	Image<float> stackImg(w2D, h2D, fc, stackPtr);
	
	double* dataPtr = mxGetDoubles(plhs[0]);
	double* psfPtr = mxGetDoubles(plhs[1]);
	
	Image<double> dataImg(w, h, d, dataPtr);
	Image<double> psfImg(w, h, d, psfPtr);
	
	for (int p = 0; p < pc; p++)
	{
		BackprojectionHelper::backprojectRS(
			stackImg, projections[p], dataImg, num_threads);
		
		BackprojectionHelper::backprojectPsfRS(
			stackImg, projections[p], psfImg, num_threads);
	}
}

void crashOnArgument(int argumentID)
{
	mexPrintf("bad argument #%d\n", argumentID + 1);
	mexErrMsgTxt("aborting backprojection");
}

gravis::d4Matrix anglesToMatrix(
	double tdrot, double tilt, double narot, 
	double cx, double cy,
	double w, double h, double d	)
{
	const double phi = DEG2RAD(tdrot);
	const double theta = DEG2RAD(tilt);
	const double chi = DEG2RAD(narot);
	
	const double sinPhi = sin(phi), cosPhi = cos(phi);
	const double sinTheta = sin(theta), cosTheta = cos(theta);
	const double sinChi = sin(chi), cosChi = cos(chi);
	
	// Center at size/2
	gravis::d4Matrix Tc(
		1, 0, 0, -w/2,
		0, 1, 0, -h/2,
		0, 0, 1, -d/2,
		0, 0, 0, 1);
	
	// The particle rotates clockwise about its z axis.
	gravis::d4Matrix R0(
		 cosPhi, sinPhi, 0, 0,
		-sinPhi, cosPhi, 0, 0,
			  0,      0, 1, 0,
			  0,      0, 0, 1);
	
	// the particle rotates clockwise about its new x axis.
	gravis::d4Matrix R1(
		 1,         0,        0, 0,
		 0,  cosTheta, sinTheta, 0,
		 0, -sinTheta, cosTheta, 0,
		 0,         0,        0, 1);
	
	// the particle rotates clockwise about its new z axis.
	gravis::d4Matrix R2(
		 cosChi, sinChi, 0, 0,
		-sinChi, cosChi, 0, 0,
			  0,      0, 1, 0,
			  0,      0, 0, 1);
	
	// Shift to 3D position of center
	gravis::d4Matrix Ts(
		1, 0, 0, cx,
		0, 1, 0, cy,
		0, 0, 1, 0,
		0, 0, 0, 1);
	
	return Ts * R2 * R1 * R0 * Tc;
}
