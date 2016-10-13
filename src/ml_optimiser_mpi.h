/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef ML_OPTIMISER_MPI_H_
#define ML_OPTIMISER_MPI_H_
#include "src/mpi.h"
#include "src/ml_optimiser.h"

#define MPITAG_JOB_REQUEST 0
#define MPITAG_JOB_REPLY 1
#define MPITAG_METADATA 2
#define MPITAG_RANDOMSEED 3
#define MPITAG_IMAGE 4
#define MPITAG_PACK 5
#define MPITAG_RFLOAT 6
#define MPITAG_INT 7
#define MPITAG_IDENTIFIER 8

class MlOptimiserMpi: public MlOptimiser
{
	std::vector<int> cudaDeviceShares;

public:
	MpiNode *node;

#ifdef TIMINGMPI
    int MPIR_PACK, MPIR_ALLREDUCE, MPIR_UNPACK, MPIR_EXP, MPIR_MAX, MPIR_BCAST;
#endif

    // Only process unfinished micrographs in movie-refinement on a per-micrograph basis
    bool only_do_unfinished_movies;

    // Original verb
    int ori_verb;

	/** Destructor, calls MPI_Finalize */
    ~MlOptimiserMpi()
    {
        delete node;
    }

    /** Read
     * This could take care of mpi-parallelisation-dependent variables
     */
    void read(int argc, char **argv);

    /** Finalise
     * Free some memory
     */
    void finalise();

    void initialise();

    /** Initialise the work load: divide images equally over all nodes
     * Also initialise the same random seed for all nodes
     */
    void initialiseWorkLoad();

    /** Perform individual power spectra calculation in parallel */
    void calculateSumOfPowerSpectraAndAverageImage(MultidimArray<RFLOAT> &Mavg);

    /** Expectation
     *  This cares care of gathering all weighted sums after the expectation
     */
    void expectation();

    /** After expectation combine all weighted sum arrays across all nodes
     *  Use read/write to temporary files instead of MPI
     */
    void combineAllWeightedSumsViaFile();

    /** Join the sums from two random halves
     *  Use read/write to temporary files instead of MPI
     */
    void combineWeightedSumsTwoRandomHalvesViaFile();

    /** After expectation combine all weighted sum arrays across all nodes
     */
    void combineAllWeightedSums();

    /** Join the sums from two random halves
     */
    void combineWeightedSumsTwoRandomHalves();

    /** Maximization
     * This takes care of the parallel reconstruction of the classes
     */
    void maximization();

    /** Perform unregularized reconstruction
      * With the aim of performing solvent mask corrected FSC inside the auto-refine
      */
    void reconstructUnregularisedMapAndCalculateSolventCorrectedFSC();

    /**
     *  Write temporary data and weight arrays from the backprojector to disc to allow unregularized reconstructions
     */
    void writeTemporaryDataAndWeightArrays();

    /**
     *  Read temporary data and weight arrays from disc and perform unregularized reconstructions
     *  Also write the unregularized reconstructions to disc.
     */
    void readTemporaryDataAndWeightArraysAndReconstruct(int iclass, int ihalf);

    /**
     * Join two independent reconstructions ate the lowest frequencies to avoid convergence in distinct orientations
     */
    void joinTwoHalvesAtLowResolution();

    /** When refining two random halves separately, the master receives both models, calculates FSC and the power of their difference
     *  and sends these curves, together with new tau2_class estimates to all slaves...
     */
    void compareTwoHalves();

    /** Do the real work
     * Expectation is split in image subsets over all nodes, each reconstruction is done on a separate node
     */
    void iterate();


    /**
     * Process movies one micrograph at a time
     * This reduces time and RAM for movie expansion and may thus be very useful for large data sets (with millions of particle movie frames)
     * It does affect the parallelisatione efficiency though, especially when there are few particles per micrograph
     */
    void processMoviesPerMicrograph(int argc, char **argv);


};


#endif /*  ML_OPTIMISER_MPI_H_ */
