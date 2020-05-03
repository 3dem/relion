#include "include/mex_parser.h"
#include <programs/discover_motif.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MexParser mp(nlhs, plhs, nrhs, prhs);
	
    DiscoverMotifProgram dmp;
	
	dmp.pointsFn = mp.getString(           MexParser::Right, 0, "Point cloud (.csv)");
	
	dmp.spacingMin = mp.getDouble(         MexParser::Right, 1, "Min. lattice spacing", 6.5);	
	dmp.spacingMax = mp.getDouble(         MexParser::Right, 2, "Max. lattice spacing", 8.5);
	dmp.neighborCount = (int) mp.getDouble(MexParser::Right, 3, "Number of neighbors", 6);
	dmp.motifCount = (int) mp.getDouble(   MexParser::Right, 4, "Number of motifs to return", 5);
	dmp.cloudBinning =       mp.getDouble( MexParser::Right, 5, "Binning level of the point cloud", 8.0);
	
	dmp.minSc = mp.getDouble(              MexParser::Right, 6, "Minimal score for a point to be used in motif refinement", 3);
	
	dmp.radius = mp.getDouble(             MexParser::Right, 7, "Tolerance [Px]", 0.25);
	dmp.reg = mp.getDouble(                MexParser::Right, 8, "Regularizer constant: inhibits movement of motif points during refinement", 0.1);
	dmp.motIters = (int) mp.getDouble(     MexParser::Right, 9, "Number of motifs to test", 1000);
	dmp.evalIters = (int) mp.getDouble(    MexParser::Right, 10, "Number of points to test against", 10000);

    dmp.refineDiscovery = mp.checkOption(  MexParser::Right, 11, "Refine motif locally during discovery", false);
    dmp.refineDetection = !mp.checkOption( MexParser::Right, 12, "Skip local refinement for motif detection", false);
    
    dmp.num_threads = (int) mp.getDouble(  MexParser::Right, 13, "Number of threads", 1);
    
    dmp.outFn = mp.getString(              MexParser::Right, 14, "Output filename");
    
    dmp.doDistHist = mp.checkOption(       MexParser::Right, 15, "Compute distance density, then quit", false);
    dmp.maxDistHist = mp.getDouble(        MexParser::Right, 16, "Max. distance", 50);
    dmp.dhBins = (int) mp.getDouble(       MexParser::Right, 17, "Number of bins", 256);
    
    dmp.doNcHist = mp.checkOption(         MexParser::Right, 18, "Compute neighbor-count histogram, then quit", false);
    dmp.maxNcHist = (int) mp.getDouble(    MexParser::Right, 19, "Max. neighbor count", 20);

    if (!mp.finalize()) return;

	dmp.run();
}
