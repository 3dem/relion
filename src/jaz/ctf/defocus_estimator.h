#ifndef DEFOCUS_ESTIMATOR_H
#define DEFOCUS_ESTIMATOR_H

class DefocusEstimator
{
	public:
		
		DefocusEstimator();
		
		void read(IOParser& parser, int argc, char *argv[]);

        void init(int verb, int s, int nr_omp_threads,
                  bool debug, bool diag,
                  ReferenceMap* reference,
                  ObservationModel* obsModel);
};

#endif
