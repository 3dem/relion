#ifndef ABERRATION_ESTIMATOR_H
#define ABERRATION_ESTIMATOR_H

#include <src/image.h>

class IOParser;
class ReferenceMap;
class ObservationModel;

class AberrationEstimator
{
	public:

		AberrationEstimator();


		void read(IOParser& parser, int argc, char *argv[]);

		void init(
				int verb, int nr_omp_threads,
				bool debug, bool diag, std::string outPath,
				ReferenceMap* reference, ObservationModel* obsModel);

		// Compute per-pixel information for one micrograph
		void processMicrograph(
				long g, MetaDataTable& mdt,
				const std::vector<Image<Complex>>& obs,
				const std::vector<Image<Complex>>& pred);

		// Sum up per-pixel information from all micrographs,
		// then fit beam-tilt model to the per-pixel fit
		void parametricFit(
				const std::vector<MetaDataTable>& mdts,
				MetaDataTable& optOut, std::vector <FileName> &fn_eps);

		// Has this mdt been processed already?
		bool isFinished(const MetaDataTable& mdt);


	private:

		// cmd. line options (see read())
		double kmin;
		int aberr_n_max;
		double xring0, xring1;

		// parameters obtained through init()
		int verb, nr_omp_threads;
		bool debug, diag, ready;
		std::string outPath;

		std::vector<int> s, sh;
		std::vector<double> angpix;

		ReferenceMap* reference;
		ObservationModel* obsModel;

};

#endif
