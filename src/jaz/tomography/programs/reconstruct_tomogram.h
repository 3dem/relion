#ifndef TOMO_BP_PROGRAM_H
#define TOMO_BP_PROGRAM_H

#include <string>
#include <vector>
#include <src/filename.h>
#include <src/time.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/tomography/tomogram_set.h>

class TomoBackprojectProgram
{
	public:
		
		TomoBackprojectProgram(){}
			
			int n_threads;
			int w, h, d;
			double spacing, angpix_spacing, x0, y0, z0, taperDist, taperFalloff;
			FileName tomoName, outFn;
			bool applyPreWeight, applyWeight, applyCtf, doWiener, zeroDC, FourierCrop;
            bool do_multiple, do_only_unfinished;
			double SNR;

            std::vector<long> tomoIndexTodo;
			OptimisationSet optimisationSet;
			TomogramSet tomogramSet;

		void readParameters(int argc, char *argv[]);
		void initialise();
        void run(int rank = 0, int size = 1);
        void getProjectMatrices(Tomogram &tomogram, MetaDataTable &tomogramTable);
        void reconstructOneTomogram(int tomoIndex);

    private:
        FileName getOutputFileName(int index);
};

#endif
