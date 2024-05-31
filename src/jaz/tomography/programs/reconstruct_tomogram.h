#ifndef TOMO_BP_PROGRAM_H
#define TOMO_BP_PROGRAM_H

#include <string>
#include <vector>
#include <src/filename.h>
#include <src/time.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/optimization/lbfgs.h>

class ReconstructSnrOptimisation : public FastDifferentiableOptimization
{
    public:

        ReconstructSnrOptimisation(const MultidimArray<RFLOAT> &CTF, const MultidimArray<RFLOAT> &SNR, double _lambda, int verb = 0)
        {
            // set CTF^2 vector
            int size = XSIZE(CTF);
            if (XSIZE(SNR) != size) REPORT_ERROR("ReconstructSnrOptimisation ERROR: CTF and SNR arrays are of unequal size!");
            first_peak = -1;
            bool has_reached_one = false;
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(CTF)
            {
                double ctf = DIRECT_MULTIDIM_ELEM(CTF, n);
                if (ctf > 0.99) has_reached_one = true;
                if (has_reached_one && ctf < 0.99 && first_peak < 0) first_peak = n;
            }

            ctf2.resize(size - first_peak);
            snr.resize(size - first_peak);
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(CTF)
            {
                double ctf = DIRECT_MULTIDIM_ELEM(CTF, n);
                if (n >= first_peak)
                {
                    VEC_ELEM(ctf2, n - first_peak) = ctf * ctf;
                    VEC_ELEM(snr, n - first_peak) = DIRECT_MULTIDIM_ELEM(SNR, n);
                }

            }

            // Set D matrix
            D.initIdentity(size - first_peak);
            for (int i=1; i < size - first_peak - 1; i++)
            {
                MAT_ELEM(D, i+1, i) = -1;
            }

            // set lambda
            lambda = _lambda;

        }
        double gradAndValue(const std::vector<double>& x, std::vector<double>& gradDest) const;
        double f(const std::vector<double>& x, void* tempStorage) const;

        int getFirstPeak() {return first_peak;};

    private:
        Matrix1D<double> ctf2, snr;
        Matrix2D<double> D;
        double lambda;
        int first_peak;

};

class TomoBackprojectProgram
{
	public:
		
		TomoBackprojectProgram(){}
			
			int n_threads;
			int w, h, d;
			double spacing, angpix_spacing, x0, y0, z0, taperDist, taperFalloff;
			FileName tomoName, outFn;
			bool applyPreWeight, applyWeight, applyCtf, doWiener, zeroDC, FourierCrop, fourierInversion;
            bool do_multiple, do_only_unfinished;
	     	bool do_even_odd_tomograms, do_2dproj, ctf_intact_first_peak;
            int centre_2dproj, thickness_2dproj;
			double SNR;
            double tiltAngleOffset;
            double BfactorPerElectronDose;
            double lambda;

            std::vector<long> tomoIndexTodo;
			OptimisationSet optimisationSet;
			TomogramSet tomogramSet;

		void readParameters(int argc, char *argv[]);
		void initialise(bool verbose = true);
        void run(int rank = 0, int size = 1);
        void writeOutput(bool do_all_metadata = false);
        void initialiseCtfScaleFactors(int tomoIndex, Tomogram &tomogram);
        MultidimArray<RFLOAT> getCtfCorrectedSNR(const MultidimArray<RFLOAT> &FSC, const MultidimArray<RFLOAT>  &Fctf, double lambda, int verb = 0);

        void reconstructOneTomogram(int tomoIndex, bool doEven, bool doOdd);
        void reconstructOneTomogramFourier(int tomoIndex);
        void setMetaDataAllTomograms();
    private:
        FileName getOutputFileName(int index, bool nameEven, bool nameOdd, bool is_2dproj = false);
};

#endif
