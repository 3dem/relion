#include "include/mex_parser.h"
#include <programs/find_lattice.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MexParser mp(nlhs, plhs, nrhs, prhs);	
	FindLatticeProgram flp;
	
	flp.tomoFn = mp.getString(        MexParser::Right, 0, "Tomogram filename");
	flp.outFn = mp.getString(         MexParser::Right, 1, "Output filename");		
	flp.spacing_ang = mp.getDouble(   MexParser::Right, 2, "Lattice spacing [Å]");
	flp.angpix = mp.getDouble(        MexParser::Right, 3, "Pixel size [Å]");
	flp.filter_width = mp.getDouble(  MexParser::Right, 4, "Filter width [~Px]", 10.0);
	flp.minValue = mp.getDouble(      MexParser::Right, 5, "Minimal abs. value of bandpass filtered image", 0.333);
	flp.minDensity = mp.getDouble(    MexParser::Right, 6, "Minimal lattice density", 0.0);
	flp.taper = mp.getDouble(         MexParser::Right, 7, "Tapering distance [Px]", 0.0);
	flp.overtones = (int)mp.getDouble(MexParser::Right, 8, "Number of overtones", 0);
	flp.noise = mp.checkOption(       MexParser::Right, 9, "Replace tomogram contents by white noise");
			
	flp.n_threads = (int)mp.getDouble(MexParser::Right, 10, "Number of threads", 1);
	
	if (!mp.finalize()) return;	
	
	flp.run();
}
