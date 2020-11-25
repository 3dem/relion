#include <src/args.h>
#include <src/jaz/tomography/programs/discover_motif.h>


int main(int argc, char *argv[])
{
	IOParser parser;
	
    DiscoverMotifProgram dmp;

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");

	dmp.pointsFn = parser.getOption("--p", "Point cloud (.csv)");
	dmp.spacingMin = textToDouble(parser.getOption("--d0", "Min. lattice spacing", "6.5"));
	dmp.spacingMax = textToDouble(parser.getOption("--d1", "Max. lattice spacing", "8.5"));
	dmp.neighborCount = textToInteger(parser.getOption("--n", "Number of neighbors", "6"));
	dmp.motifCount = textToInteger(parser.getOption("--m", "Number of motifs to return", "5"));
	dmp.cloudBinning = textToDouble(parser.getOption("--bin", "Binning of the point cloud", "8.0"));
	dmp.minSc = textToDouble(parser.getOption("--min_score", "Minimal score for a point to be used in motif refinement", "3"));

	dmp.radius = textToDouble(parser.getOption("--tol", "Tolerance [Px]", "0.25"));
	dmp.reg = textToDouble(parser.getOption("--reg", "Regularizer constant: inhibits movement of motif points during refinement", "0.1"));
	dmp.motIters = textToInteger(parser.getOption("--nm", "Number of motifs to test", "1000"));
	dmp.evalIters = textToInteger(parser.getOption("--ne", "Number of points to test against", "10000"));

	dmp.refineDiscovery = parser.checkOption("--ref-disc", "Refine motif locally during discovery");
	dmp.refineDetection = !parser.checkOption("--no_ref-det", "Skip local refinement for motif detection");

	dmp.num_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));

	dmp.outFn = parser.getOption("--o", "Output filename");

	int dh_section = parser.addSection("Distance histogram");
	dmp.doDistHist = parser.checkOption("--dh", "Compute the distance density, then quit");
	dmp.maxDistHist = textToDouble(parser.getOption("--dhmax", "Max. distance", "50"));
	dmp.dhBins = textToInteger(parser.getOption("--dhb", "Number of bins", "256"));

	int nch_section = parser.addSection("Neighbor-count histogram");
	dmp.doNcHist = parser.checkOption("--nch", "Compute neighbor-count histogram, then quit");
	dmp.maxNcHist = textToInteger(parser.getOption("--nchmax", "Max. neighbor count", "20"));

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	dmp.run();
		
	return 0;
}
