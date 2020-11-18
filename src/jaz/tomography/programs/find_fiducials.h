#ifndef FIND_FIDUCIALS_H
#define FIND_FIDUCIALS_H

#include <string>
#include <src/jaz/tomography/optimisation_set.h>


class FindFiducialsProgram
{
	public:

		FindFiducialsProgram(int argc, char *argv[]);

			std::string outDir;
			double thresh, binning_out, binning_in, beadRadius_A, rel_spacing;
			int max_MG, num_threads;
			bool diag;

			OptimisationSet optimisationSet;

		void run();

};

#endif
