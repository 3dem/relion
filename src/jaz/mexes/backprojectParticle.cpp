#include "include/mex_parser.h"
#include "include/dynamo_helper.h"

#include <programs/backproject.h>

using namespace gravis;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	MexParser mp(nlhs, plhs, nrhs, prhs);
	
	BackprojectProgram bp;
	
	bp.catFn =         mp.getString(MexParser::Right,   0, "Catalogue file (.tbl or .cat)");
	bp.tomoListFn =    mp.getString(MexParser::Right,   1, "Tomogram list filename");
	bp.boxSize = (int) mp.getDouble(MexParser::Right,   2, "Box size", 100);
	bp.outTag =        mp.getString(MexParser::Right,   3, "Output filename pattern");
	bp.WienerFract =   mp.getDouble(MexParser::Right,   4, "SNR assumed by the Wiener filter", 0.0001);
	bp.handedness =    mp.getDouble(MexParser::Right,   5, "Handedness of the tilt geometry", 1.0);
	bp.optFn =         mp.getString(MexParser::Right,   6, "Consider a CTF using parameters from the supplied file");
	bp.num_threads=(int)mp.getDouble(MexParser::Right,  7, "Number of OMP threads", 6);
	bp.realSpace =     mp.checkOption(MexParser::Right, 8, "Perform the backprojection in real space (no CTF)");
	bp.diag =          mp.checkOption(MexParser::Right, 9, "Write out diagnostic information");

	
	if (!mp.finalize()) return;
		
	bp.run();
}
