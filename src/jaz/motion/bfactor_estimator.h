#ifndef B_FACTOR_ESTIMATOR_H
#define B_FACTOR_ESTIMATOR_H

#include <src/image.h>
#include <vector>
#include <string>

class IOParser;
class ObservationModel;
class MicrographHandler;
class ReferenceMap;

class BFactorEstimator
{
	public:
		
		BFactorEstimator();
		
		void read(IOParser& parser, int argc, char *argv[]);

        void init(int verb, int s, int fc, 
				  int nr_omp_threads,
                  std::string outPath, bool debug,
                  ObservationModel* obsModel,
                  MicrographHandler* micrographHandler,
				  ReferenceMap* reference);

        void process(const std::vector<MetaDataTable>& mdts);

        bool doingAnything();
		
		
	protected:

		// read from cmd. line:
		bool doAnything, bfac_diag;
		int nr_helical_asu, f_max, f_batch;
		double helical_rise, helical_twist, 
			fit_minres, perframe_highres;
		FileName fn_sym;
		
		
		// set at init:
		int s, sh, fc;
		int verb, nr_omp_threads;
		std::string outPath;
		bool debug;
		double angpix;

		ObservationModel* obsModel;
		MicrographHandler* micrographHandler;
		ReferenceMap* reference;
		
		bool calculateBfactorSingleFrameReconstruction(
				int frame,
				const MultidimArray<RFLOAT>& fsc_frame,
				const MultidimArray<RFLOAT>& fsc_average,
				double& bfactor, double& offset, double& corr_coeff);
		
		MultidimArray<RFLOAT> maskedFSC( 
				Image<RFLOAT>& I1, 
				Image<RFLOAT>& I2,
				const Image<RFLOAT>& Imask);
		
		void writeStarFileBfactors(
				MultidimArray<RFLOAT>& perframe_bfactors);
		
};

#endif
