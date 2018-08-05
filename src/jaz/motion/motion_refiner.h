/***************************************************************************
 *
 * Author: "Jasenko Zivanov & Sjors H.W. Scheres"
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

#ifndef MOTION_REFINER_H_
#define MOTION_REFINER_H_


#include <src/ctf.h>
#include <src/image.h>
#include <src/fftw.h>
#include <src/micrograph_model.h>
#include <src/jaz/legacy_obs_model.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/parallel_ft.h>

#include <src/jaz/micrograph_handler.h>
#include <src/jaz/reference_map.h>

#include "motion_param_estimator.h"
#include "motion_estimator.h"
#include "frame_recombiner.h"
#include "bfactor_estimator.h"

#include <omp.h>


class MotionRefiner
{
	public:
		
		MotionRefiner();
		
		
		// Read command line arguments
		void read(int argc, char **argv);
		
		// Initialise some general stuff after reading
		void init();
		
		// General Running (Admiral Swimming!)
		void run();
		
		int getVerbosityLevel();
		
		
		// For original particle-polishing-like Bfactors (not used)
		//void calculateSingleFrameReconstruction(int iframe);
		
		// Get output STAR file name for this micrograph
		static FileName getOutputFileNameRoot(std::string outPath, const MetaDataTable& mdt);
		
		
	protected:
		
		// components that do the actual work
		MotionParamEstimator motionParamEstimator;
		MotionEstimator motionEstimator;
		FrameRecombiner frameRecombiner;
		BFactorEstimator bfactorEstimator;
		
		// required components
		LegacyObservationModel obsModel;
		ReferenceMap reference;
		MicrographHandler micrographHandler;
		
		// s: full image size, sh: half-size + 1, fc: frame count
		int s, sh, fc;
		
		// Verbosity
		int verb;
		
		bool debug, findShortestMovie;
		
		double angpix;
		int nr_omp_threads;
		std::string outPath;
		
		std::string starFn, movie_toReplace, movie_replaceBy;
		
		// Allow continuation of crashed jobs
		bool only_do_unfinished;
		
		bool estimateParams,
		     estimateMotion,
		     recombineFrames,
		     estimateOldBfacs,
		     generateStar;
		
		long maxMG, minMG;
		
		MetaDataTable mdt0;
		
		std::vector<MetaDataTable>
			allMdts, // all micrographs (used for B-factor computation)
			chosenMdts, // micrographs between minMG and maxMG
			motionMdts, recombMdts; // unfinished micrographs
		
		
		// combine all EPS files into one logfile.pdf
		void combineEPSAndSTARfiles();
		
		// apply changes to micrograph-filenames implied by
		// movie_path, movie_ending and movie_toReplace/replaceBy
		void adaptMovieNames();
};



#endif
