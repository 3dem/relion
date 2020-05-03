#include "include/mex_parser.h"
#include <programs/tomo_backproject.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MexParser mp(nlhs, plhs, nrhs, prhs);
	TomoBackprojectProgram tbp;
	
	tbp.stackFn = mp.getString(MexParser::Right,         0, "Tilt sequence image stack (e.g. stack.st:mrc)");
	tbp.projFn = mp.getString(MexParser::Right,          1, "Tilt projections file (*.tlt)");
	tbp.outFn = mp.getString(MexParser::Right,           2, "Output filename");	
	tbp.spacing = mp.getDouble(MexParser::Right,         3, "Binning factor (pixel spacing)", 8.0);		
	tbp.weight = mp.checkOption(MexParser::Right,        4, "Perform weighting in Fourier space (using a Wiener filter)");
	tbp.WienerFract = mp.getDouble(MexParser::Right,     5, "SNR assumed by the Wiener filter", 0.0001);	
	tbp.n_threads = (int) mp.getDouble(MexParser::Right, 6, "Number of threads", 1);
	tbp.zeroDC = mp.checkOption(MexParser::Right,        7, "Zero the DC component of each frame");	
	tbp.taperRad = mp.getDouble(MexParser::Right,        8, "Tapering distance", 0.0);
	tbp.thickness = (int) mp.getDouble(MexParser::Right, 9, "Thickness (read from .proj file by default)", -1);
	
	tbp.x0 = mp.getDouble(MexParser::Right,  10,  "X origin", 1.0);
	tbp.y0 = mp.getDouble(MexParser::Right,  11,  "Y origin", 1.0);
	tbp.z0 = mp.getDouble(MexParser::Right,  12,  "Z origin", 1.0);
	
	tbp.w = (int) mp.getDouble(MexParser::Right,  13,  "Width",  -1.0);
	tbp.h = (int) mp.getDouble(MexParser::Right,  14,  "Height", -1.0);
			
	
	if (!mp.finalize()) return;	
	
	tbp.run();
}
