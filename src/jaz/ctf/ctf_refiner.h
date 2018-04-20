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

class CtfRefiner
{	
	public:
		
		CtfRefiner();
		
		
		void read(int argc, char **argv);
		
		void init();		
		
		void run();
		
		int getVerbosityLevel();
		
		static FileName getOutputFilenameRoot(
				const MetaDataTable& mdt, std::string outPath);
		
		
	protected:
	
		RFLOAT Cs, lambda, kV;
		
		ObservationModel obsModel;
		ReferenceMap reference;
		
		TiltEstimator tiltEstimator;
	
		// Verbosity
		int verb;
	
		// Allow continuation of crashed jobs
		bool only_do_unfinished;
	
		// Perform per-particle defocus estimation?
		bool do_defocus_fit;
	
		// Perform beamtilt estimation?
		bool do_tilt_fit;
	
		bool debug,     // write out debugging info
		     clTilt,    // tilt from cmd. line
		     anisoTilt; // use experimental anisotropic tilt model
	
		long maxMG, minMG;
	
		RFLOAT angpix,
			beamtilt_x, beamtilt_y,
			beamtilt_xx, beamtilt_xy, beamtilt_yy;
	
		int nr_omp_threads;
	
		std::string starFn, outPath;
	
		MetaDataTable mdt0;
		std::vector<MetaDataTable> mdts;
		std::vector<FileName> fn_mics_process, fn_mics_ori;
	
		
		// s: full image size, sh: half-size + 1, fc: frame count
		int s, sh, fc;
	
		// Defocus_fit options
		RFLOAT defocusRange;
	
		// Astigmatism options
		bool fitAstigmatism, noGlobAstig, diag;
	
		// Fit amplitude contrast/phase shift
		bool fitPhase;
	
		// Fit spherical aberration coefficient
		bool fitCs;
	
		// only per-micrograph fits
		bool globOnly;
	
		// Tilt fit options
		RFLOAT kmin,
			testtilt_x, testtilt_y,
			testtilt_xx, testtilt_xy, testtilt_yy;
	
		bool aniso;
	
		Image<Complex> lastXY;
		Image<RFLOAT> lastW;
	
	
		// Get output STAR file name for the gth entry in the mdts
		FileName getOutputFileNameRoot(long int g, bool is_original = false);
	
		// Fir defocus for all particles on one micrograph
		void fitDefocusOneMicrograph(long g, const std::vector<Image<Complex> > &obsF, 
									 const std::vector<Image<Complex> > &preds, int verb = 0);
	
		// Write PostScript file with per-particle defocus plotted onto micrograph in blue-red color scale
		void writePerParticleDefocusEPS(long g);
	
		// Fit CTF parameters for all particles on a subset of the micrographs micrograph
		void processSubsetMicrographs(long g_start, long g_end);
	
		// Read all micrograph metadata tables back in, and combine into onelarge one
		void combineAllDefocusFitAndBeamTiltInformation(
				long g_start, long g_end);
};



#endif /* CTF_REFINER_H */
