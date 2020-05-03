#include "include/mex_parser.h"
#include <programs/dark_erase.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MexParser mp(nlhs, plhs, nrhs, prhs);
	
	DarkEraseProgram dep;
	
	dep.stackFn = mp.getString(            MexParser::Right, 0, "Tilt sequence image stack (e.g. stack.st:mrc)");	
	dep.thresh = mp.getDouble(             MexParser::Right, 1, "Threshold value (in sigma)", -2.5);
	dep.rad = mp.getDouble(                MexParser::Right, 2, "Fill radius (in pixels)", 32.0);		
	dep.outFn = mp.getString(              MexParser::Right, 3, "Output filename");
	
	dep.writeNormalized = mp.checkOption(  MexParser::Right, 4, "Only write out a normalized stack to assist in finding an optimal threshold value", false);
	dep.diag = mp.checkOption(             MexParser::Right, 5, "Write out diagnostic data", false);	
	dep.num_threads = (int)mp.getDouble(   MexParser::Right, 6, "Number of threads", 1);	
	
	if (!mp.finalize()) return;
				
	dep.run();
}
