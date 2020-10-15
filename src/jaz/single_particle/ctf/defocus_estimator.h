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

#ifndef DEFOCUS_ESTIMATOR_H
#define DEFOCUS_ESTIMATOR_H

#include <src/image.h>

class IOParser;
class ReferenceMap;
class ObservationModel;

class DefocusEstimator
{
	public:
		
		DefocusEstimator();
		
		
		void read(IOParser& parser, int argc, char *argv[]);

        void init(
				int verb, int nr_omp_threads,
				bool debug, bool diag,
				std::string outPath,
				ReferenceMap* reference,
				ObservationModel* obsModel);
		
		
		// Fit defocus for all particles on one micrograph
		void processMicrograph(
				long g, MetaDataTable& mdt, 
				const std::vector<Image<Complex>>& obs,
				const std::vector<Image<Complex>>& pred);
	
		// Write PostScript file with per-particle defocus 
		// plotted onto micrograph in blue-red color scale
		void writeEPS(const MetaDataTable &mdt);
		
		// Has this mdt been processed already?
		bool isFinished(const MetaDataTable& mdt);
		
		
	private:
		
		// cmd. line options (see read()):
		double defocusRange, kmin;	
		int max_iters;
		bool fitAstigmatism, bruteForcePre, bruteForcePost, bruteForceOnly;
		std::string fittingMode;
		
		// set at init:
		int verb, nr_omp_threads;
		bool debug, diag;
		std::string outPath;
		
		std::vector<int> s, sh;
		std::vector<double> angpix;		
		std::vector<Image<RFLOAT> > freqWeights;
		
		ReferenceMap* reference;
		ObservationModel* obsModel;
		
		bool ready;
		
		
		void bruteForceFit(
				long g, MetaDataTable& mdt, 
				const std::vector<Image<Complex>>& obs,
				const std::vector<Image<Complex>>& pred,
				std::string tag);
		
};

#endif
