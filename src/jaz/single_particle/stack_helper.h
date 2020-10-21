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
#include <src/jaz/gravis/t2Matrix.h>
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
	
		static std::vector<Image<Complex>> FourierTransform(std::vector<Image<RFLOAT> >& stack);
		
		static std::vector<Image<RFLOAT>> inverseFourierTransform(std::vector<Image<Complex> >& stack);
		
		static Image<RFLOAT> toSingleImage(const std::vector<Image<RFLOAT>> stack);
		
		static void varianceNormalize(
					std::vector<Image<Complex>>& movie, 
					bool circleCropped);
		
		static RFLOAT computePower(
					const RawImage<Complex>& movie, 
					bool circleCropped);
		
		static RFLOAT computePower(
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
