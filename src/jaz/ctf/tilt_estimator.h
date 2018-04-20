#ifndef TILT_ESTIMATOR_H
#define TILT_ESTIMATOR_H

#include <src/image.h>

class IOParser;
class ReferenceMap;
class ObservationModel;

class TiltEstimator
{
	public:
		
		TiltEstimator();
		
		void init(
				int verb, int s, int nr_omp_threads,
				bool debug, bool diag, std::string outPath,
				ReferenceMap* reference, ObservationModel* obsModel);

		void processMicrograph(
				long g, MetaDataTable& mdt, 
				const std::vector<Image<Complex>>& obs);
		
		void parametricFit(
				std::vector<MetaDataTable>& mdts, 
				double kmin, double Cs, double lambda,
				bool aniso);
		
	private:
		
		int verb, s, sh, nr_omp_threads;
		bool debug, diag, ready;
		std::string outPath;
		
		ReferenceMap* reference;
		ObservationModel* obsModel;
		
};

#endif
