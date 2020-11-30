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

#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/reference_map.h>
#include <src/image.h>

#include "tilt_estimator.h"
#include "defocus_estimator.h"
#include "bfactor_refiner.h"
#include "magnification_estimator.h"
#include "aberration_estimator.h"

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

		ObservationModel obsModel;
		ReferenceMap reference;

		TiltEstimator tiltEstimator;
		DefocusEstimator defocusEstimator;
		BFactorRefiner bfactorEstimator;
		AberrationEstimator aberrationEstimator;
		MagnificationEstimator magnificationEstimator;

		// Verbosity
		int verb;

		// Allow continuation of crashed jobs
		bool only_do_unfinished;

		// Whether to estimate defoci, B-factors, antisymmetric aberrations (incl. beam tilt),
		// symmetric aberrations and anisotropic magnification, respectively
		bool do_defocus_fit, do_bfac_fit, do_tilt_fit, do_aberr_fit, do_mag_fit, do_ctf_padding;

		bool debug,     // write out debugging info
		     diag;      // write out diagnostic info

		long maxMG, minMG;

		int nr_omp_threads;

		std::string starFn, outPath;

		MetaDataTable mdt0;
		std::vector<MetaDataTable> allMdts, unfinishedMdts;

		// Fit CTF parameters for all particles on a subset of the micrographs micrograph
		void processSubsetMicrographs(long g_start, long g_end);

		// Combine all .stars and .eps files
		std::vector<MetaDataTable> merge(const std::vector<MetaDataTable>& mdts, std::vector <FileName> &fn_eps);
};



#endif /* CTF_REFINER_H */
