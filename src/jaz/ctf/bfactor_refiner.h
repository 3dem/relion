#ifndef BFACTOR_REFINER_H
#define BFACTOR_REFINER_H

#include <src/image.h>
#include <src/jaz/gravis/t2Vector.h>

class IOParser;
class ReferenceMap;
class ObservationModel;

class BFactorRefiner
{
	public:
		
		BFactorRefiner();
		
		
		void read(IOParser& parser, int argc, char *argv[]);

        void init(
				int verb, int s, int nr_omp_threads,
				bool debug, bool diag,
				std::string outPath,
				ReferenceMap* reference,
				ObservationModel* obsModel);
		
		
		// Fit B-factors for all particles on one micrograph
		void processMicrograph(
				long g, MetaDataTable& mdt, 
				const std::vector<Image<Complex>>& obs,
				const std::vector<Image<Complex>>& pred);
		
		// Combine all .stars and .eps files
		std::vector<MetaDataTable> merge(const std::vector<MetaDataTable>& mdts);
	
		// Write PostScript file with per-particle defocus 
		// plotted onto micrograph in blue-red color scale
		void writeEPS(const MetaDataTable &mdt);
		
		// Has this mdt been processed already?
		bool isFinished(const MetaDataTable& mdt);
		
		
	private:
		
		// cmd. line options (see read()):
		double kmin, min_scale, min_B, max_B;	
		bool perMicrograph;
		
		// set at init:
		int verb, s, sh, nr_omp_threads;
		bool debug, diag;
		std::string outPath;
		double angpix;
		
		Image<RFLOAT> freqWeight;
		
		ReferenceMap* reference;
		ObservationModel* obsModel;
		
		bool ready;
		
		static gravis::d2Vector findSigmaKRec1D(
				const std::vector<double>& t_rad, 
				const std::vector<double>& s_rad, 
				double B0, double B1, double min_scale,
				int steps, int depth);
		
		static gravis::d2Vector findBKRec2D(
				const Image<Complex>& obs,
				const Image<Complex>& pred,
				const Image<RFLOAT>& weight,
                double B0, double B1, double min_scale,
                int steps, int depth);
};

#endif
