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

#ifndef CTF_REFINER_H
#define CTF_REFINER_H

#include <src/jaz/obs_model.h>
#include <src/jaz/reference_map.h>
#include <src/image.h>

#include "tilt_estimator.h"
#include "defocus_estimator.h"

class CtfRefiner
{	
	public:
		
		CtfRefiner();
		
		
		void read(int argc, char **argv);		
		void init();		
		void run();		
		void finalise();
		
		
		int getVerbosityLevel();
		
		static FileName getOutputFilenameRoot(
				const MetaDataTable& mdt, std::string outPath);
		
		
	protected:
	
		RFLOAT Cs, lambda, kV;
		
		ObservationModel obsModel;
		ReferenceMap reference;
		
		TiltEstimator tiltEstimator;
		DefocusEstimator defocusEstimator;
	
		// Verbosity
		int verb;
	
		// Allow continuation of crashed jobs
		bool only_do_unfinished;
	
		// Perform per-particle defocus estimation?
		bool do_defocus_fit;
	
		// Perform beamtilt estimation?
		bool do_tilt_fit;
		
		// Estimate anisotropic magnification?
		bool do_anisotropy;
	
		bool debug,     // write out debugging info
		     diag,      // write out diagnostic info
		     clTilt,    // tilt from cmd. line
		     anisoTilt; // use experimental anisotropic tilt model
	
		long maxMG, minMG;
	
		RFLOAT angpix,
			beamtilt_x, beamtilt_y,
			beamtilt_xx, beamtilt_xy, beamtilt_yy;
	
		int nr_omp_threads;
	
		std::string starFn, outPath;
	
		MetaDataTable mdt0;
		std::vector<MetaDataTable> allMdts, unfinishedMdts;	
		
		// s: full image size, sh: half-size + 1, fc: frame count
		int s, sh, fc;
			
		// Fit CTF parameters for all particles on a subset of the micrographs micrograph
		void processSubsetMicrographs(long g_start, long g_end);
};



#endif /* CTF_REFINER_H */
