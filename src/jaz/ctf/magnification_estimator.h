#ifndef MAG_ESTIMATOR_H
#define MAG_ESTIMATOR_H

#include <vector>
#include <string>

#include <src/complex.h>
#include <src/image.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Vector.h>

class IOParser;
class ReferenceMap;
class ObservationModel;
class MetaDataTable;

class MagnificationEstimator
{
	public:
		
		MagnificationEstimator();
		
		void read(IOParser& parser, int argc, char *argv[]);
		
		void init(
				int verb, int s, int nr_omp_threads,
				bool debug, bool diag, std::string outPath,
				ReferenceMap* reference, ObservationModel* obsModel);
		
		// Compute per-pixel information for one micrograph
		void processMicrograph(
				long g, MetaDataTable& mdt, 
				const std::vector<Image<Complex>>& obs,
				const std::vector<Image<Complex>>& pred,
				const std::vector<Volume<gravis::t2Vector<Complex>>>& predGradient);
		
		// Sum up per-pixel information from all micrographs, 
		// then fit beam-tilt model to the per-pixel fit
		void parametricFit(
				std::vector<MetaDataTable>& mdts, 
				MetaDataTable& optOut);
		
		// Has this mdt been processed already?
		bool isFinished(const MetaDataTable& mdt);
		
		
	private:
				
		// cmd. line options (see read())
		double kmin;
		bool adaptAstig, perMgAstig;
		
		// parameters obtained through init()
		int verb, s, sh, nr_omp_threads;
		bool debug, diag, ready;
		std::string outPath;
		double angpix;
		
		ReferenceMap* reference;
		ObservationModel* obsModel;
};

#endif
