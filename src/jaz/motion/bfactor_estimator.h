/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
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

#ifndef OLD_B_FACTOR_ESTIMATOR_H
#define OLD_B_FACTOR_ESTIMATOR_H

#include <src/image.h>
#include <vector>
#include <string>

class IOParser;
class LegacyObservationModel;
class MicrographHandler;
class ReferenceMap;

class OldBFactorEstimator
{
	public:
		
		OldBFactorEstimator();
		
		void read(IOParser& parser, int argc, char *argv[]);

        void init(int verb, int s, int fc, 
				  int nr_omp_threads,
                  std::string outPath, bool debug,
                  LegacyObservationModel* obsModel,
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

		LegacyObservationModel* obsModel;
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
