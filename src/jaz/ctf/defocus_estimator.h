#ifndef DEFOCUS_ESTIMATOR_H
#define DEFOCUS_ESTIMATOR_H

#include <src/image.h>

class IOParser;
class ReferenceMap;
class ObservationModel;

class DefocusEstimator
{
	public:
		
		DefocusEstimator();
		
		
		void read(IOParser& parser, int argc, char *argv[]);

        void init(
				int verb, int s, int nr_omp_threads,
				bool debug, bool diag,
				std::string outPath,
				ReferenceMap* reference,
				ObservationModel* obsModel);
		
		
		// Fit defocus for all particles on one micrograph
		void processMicrograph(
				long g, MetaDataTable& mdt, 
				const std::vector<Image<Complex>>& obsF);
		
		// Combine all .stars and .eps files
		void merge(const std::vector<MetaDataTable>& mdts, MetaDataTable& mdtOut);
	
		// Write PostScript file with per-particle defocus 
		// plotted onto micrograph in blue-red color scale
		void writeEPS(const MetaDataTable &mdt);
		
		// Has this mdt been processed already?
		bool isFinished(const MetaDataTable& mdt);
		
		
	private:
		
		// cmd. line options (see read()):
		double defocusRange, kmin;	
		bool fitAstigmatism, noGlobAstig, fitPhase, fitCs, globOnly;
		
		// set at init:
		int verb, s, sh, nr_omp_threads;
		bool debug, diag;
		std::string outPath;
		double angpix;
		
		Image<RFLOAT> freqWeight;
		
		ReferenceMap* reference;
		ObservationModel* obsModel;
		
		bool ready;
};

#endif
