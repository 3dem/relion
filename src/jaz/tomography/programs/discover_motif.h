#ifndef DISCOVER_MOTIF_PROGRAM_H
#define DISCOVER_MOTIF_PROGRAM_H

#include <string>

class DiscoverMotifProgram
{
	public:
		
		DiscoverMotifProgram(){}
		
			int num_threads, dhBins, maxNcHist, neighborCount, motifCount;

			double spacingMin, spacingMax, maxDistHist, minSc, cloudBinning;	
			std::string pointsFn, outFn;
			bool doDistHist, doNcHist, refineDiscovery, refineDetection;

			double radius, reg;
			int motIters, evalIters;
			
		void run();
};

#endif
