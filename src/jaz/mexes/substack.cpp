#include "include/mex_parser.h"
#include <programs/substack.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MexParser mp(nlhs, plhs, nrhs, prhs);
	SubstackProgram sp;
	
	sp.catFn =            mp.getString(   MexParser::Right, 0, "Catalogue file (.tbl or .cat)");
	sp.tomoListFn =       mp.getString(   MexParser::Right, 1, "Tomogram list filename");
	sp.boxSize =     (int)mp.getDouble(   MexParser::Right, 2, "Box size", 100);
	sp.binning =          mp.getDouble(   MexParser::Right, 3, "Binning factor", 1.0);
	sp.outTag =           mp.getString(   MexParser::Right, 4, "Output filename pattern");
	sp.num_threads = (int)mp.getDouble(   MexParser::Right, 5, "Number of OMP threads", 6);
	sp.diag =             mp.checkOption( MexParser::Right, 6, "Write out diagnostic information");
	sp.do_center =       !mp.checkOption( MexParser::Right, 7, "Do not subtract the mean from the voxel values");
	sp.diag =             mp.checkOption( MexParser::Right, 8, "Write out diagnostic information");
	
	if (!mp.finalize()) return;	
	
	sp.run();
}
