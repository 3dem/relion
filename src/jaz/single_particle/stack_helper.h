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

#ifndef STACK_HELPER_H
#define STACK_HELPER_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/single_particle/volume.h>
#include <src/jaz/single_particle/gravis/t2Matrix.h>
#include <src/jaz/single_particle/parallel_ft.h>
#include <vector>

class Projector;
class ObservationModel;

class StackHelper
{
	public:
		
		static std::vector<MetaDataTable> splitByMicrographName(const MetaDataTable& mdt);
		
		static MetaDataTable merge(const std::vector<MetaDataTable>& mdts);
		
		static std::vector<MetaDataTable> splitByStack(const MetaDataTable* mdt);
		
		static std::vector<Image<RFLOAT> > loadStack(
				const MetaDataTable* mdt, std::string path = "", int threads = 1);
		
		static std::vector<Image<Complex> > loadStackFS(
				const MetaDataTable& mdt, 
				std::string path = "",
				int threads = 1,
				bool centerParticle = false,
				ObservationModel* obs = 0);
		
		static void saveStack(std::vector<Image<RFLOAT> >& stack, std::string fn);
		
		static std::vector<std::vector<Image<RFLOAT>>> loadMovieStack(
				const MetaDataTable* mdt, std::string moviePath);
	
		// For movies in file
		static std::vector<std::vector<Image<Complex>>> extractMovieStackFS(
				const MetaDataTable* mdt,
				Image<RFLOAT>* gainRef, MultidimArray<bool>* defectMask, std::string movieFn,
				double outPs, double coordsPs, double moviePs, double dataPs,
				int squareSize, int threads,
				bool loadData = true, int firstFrame = 0, int lastFrame = -1,
				RFLOAT hot = -1.0, bool verbose = false, bool saveMemory = false,
				const std::vector<std::vector<gravis::d2Vector>>* offsets_in = 0,
				std::vector<std::vector<gravis::d2Vector>>* offsets_out = 0);
				
		// For movies in memory
		static std::vector<std::vector<Image<Complex>>> extractMovieStackFS(
				const MetaDataTable* mdt, std::vector<MultidimArray<float> > &mgStack,
				double outPs, double coordsPs, double moviePs, double dataPs,
				int squareSize, int threads,
				bool loadData = true, 
				bool verbose = false, 
				const std::vector<std::vector<gravis::d2Vector>>* offsets_in = 0,
				std::vector<std::vector<gravis::d2Vector>>* offsets_out = 0);

		static std::vector<Image<Complex>> FourierTransform(std::vector<Image<RFLOAT> >& stack);
		
		static std::vector<Image<RFLOAT>> inverseFourierTransform(std::vector<Image<Complex> >& stack);
		
		static Image<RFLOAT> toSingleImage(const std::vector<Image<RFLOAT>> stack);
		
		static void varianceNormalize(
					std::vector<Image<Complex>>& movie, 
					bool circleCropped = false);
		
		static std::vector<double> powerSpectrum(
					const std::vector<std::vector<Image<Complex>>>& stack);
		
		static std::vector<double> varSpectrum(
					const std::vector<std::vector<Image<Complex>>>& stack);
		
		static std::vector<double> powerSpectrum(
					const std::vector<std::vector<Image<Complex>>>& obs,
					const std::vector<Image<Complex> >& signal);
};

#endif
