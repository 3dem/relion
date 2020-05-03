#include "include/mex_parser.h"
#include <programs/subtomo.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MexParser mp(nlhs, plhs, nrhs, prhs);
	SubtomoProgram sp;
	
	sp.catFn =         mp.getString(MexParser::Right,   0,  "Catalogue file (.tbl or .cat)");
	sp.tomoListFn =    mp.getString(MexParser::Right,   1,  "Tomogram list filename");
	sp.boxSize = (int) mp.getDouble(MexParser::Right,   2,  "Box size", 100);	
	sp.outTag =        mp.getString(MexParser::Right,   3,  "Output filename pattern");	
	sp.WienerFract =   mp.getDouble(MexParser::Right,   4,  "SNR assumed by the Wiener filter", 0.0001);
	sp.handedness =    mp.getDouble(MexParser::Right,   5,  "Handedness of the tilt geometry", 1.0);	
	sp.optFn =         mp.getString(MexParser::Right,   6,  "Consider a CTF using parameters from the supplied file");	
	sp.num_threads=(int)mp.getDouble(MexParser::Right,  7,  "Number of OMP threads", 6);
	sp.binning =       mp.getDouble(MexParser::Right,   8,  "Binning factor", 1.0);
	sp.do_rotate =     mp.checkOption(MexParser::Right, 9,  "Rotate the particles according to current angles");			
	sp.taper =         mp.getDouble(MexParser::Right,   10, "Taper against the sphere by this number of pixels", 5.0);
	sp.do_center =    !mp.checkOption(MexParser::Right, 11, "Do not subtract the mean from the voxel values");	
	sp.diag =          mp.checkOption(MexParser::Right, 12, "Write out diagnostic information");		
	
	if (!mp.finalize()) return;	
	
	sp.run();
}
