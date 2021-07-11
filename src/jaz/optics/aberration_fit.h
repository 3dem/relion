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

#ifndef ABERRATION_FIT_H
#define ABERRATION_FIT_H

#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/math/tensor2x2.h>

class Tomogram;
class TomoReferenceMap;
class AberrationsCache;


class AberrationBasis
{
	public:

		AberrationBasis(int dims);

		std::vector<double> coefficients;

		virtual void getBasisValues(double x, double y, double* dest) = 0;
		virtual void _offsetCtf(
				double local_Cs, double lambda,
				double rad_azimuth, double defocus_average, double defocus_deviation,
				double K1, double K2, double K3, MetaDataTable& mdt, int particle) = 0;

		void offsetCtf(MetaDataTable& mdt, int particle);
};

class OriginalBasis : public AberrationBasis
{
	public:

		OriginalBasis();

		void getBasisValues(double x, double y, double* dest);
		void _offsetCtf(
				double local_Cs, double lambda,
				double rad_azimuth, double defocus_average, double defocus_deviation,
				double K1, double K2, double K3, MetaDataTable& mdt, int particle);
};

namespace aberration
{
	class EvenData
	{
		public:

			double Axx, Axy, Ayy, bx, by;

			EvenData& operator+=(const EvenData& d);
			EvenData& operator*=(double d);

			static void write(const RawImage<EvenData>& evenData, std::string filename);
			static BufferedImage<EvenData> read(std::string filename);

	};

	class OddData
	{
		public:

			double a;
			dComplex b;

			OddData& operator+=(const OddData& d);
			OddData& operator*=(double d);

			static void write(const RawImage<OddData>& evenData, std::string filename);
			static BufferedImage<OddData> read(std::string filename);
	};

	class EvenSolution
	{
		public:

			BufferedImage<dComplex> optimum;
			BufferedImage<double> phaseShift;
			BufferedImage<Tensor2x2<double>> weight;
	};

	class OddSolution
	{
		public:

			BufferedImage<dComplex> optimum;
			BufferedImage<double> phaseShift;
			BufferedImage<double> weight;
	};
}

class AberrationFit
{
	public:

		static OriginalBasis fitBasic(Image<RFLOAT> phase, Image<RFLOAT> weight, double angpix);
		static Image<RFLOAT> draw(AberrationBasis* fit, double angpix, int s);


		static void considerParticle(
				ParticleIndex part_id,
				const Tomogram& tomogram,
				const TomoReferenceMap& referenceMap,
				const ParticleSet& dataSet,
				const AberrationsCache& aberrationsCache,
				bool flip_value,
				const BufferedImage<float>& freqWeights,
				const BufferedImage<float>& doseWeights,
				const BufferedImage<int>& xRanges,
				int f0, int f1,
				BufferedImage<aberration::EvenData>& even_out,
				BufferedImage<aberration::OddData>& odd_out);


		static aberration::EvenSolution solveEven(
				const BufferedImage<aberration::EvenData>& data);

		static std::vector<double> fitEven(
				const aberration::EvenSolution& solution,
				int n_bands,
				const std::vector<double>& initialCoeffs,
				double pixelSize,
				const std::string& prefix,
				bool writeImages);

		static std::vector<double> solveAndFitEven(
				const BufferedImage<aberration::EvenData>& data,
				int n_bands,
				const std::vector<double>& initialCoeffs,
				double pixelSize,
				const std::string& prefix,
				bool writeImages);


		static aberration::OddSolution solveOdd(
				const BufferedImage<aberration::OddData>& data);

		static std::vector<double> fitOdd(
				const aberration::OddSolution& solution,
				int n_bands,
				const std::vector<double>& initialCoeffs,
				double pixelSize,
				const std::string& prefix,
				bool writeImages);

		static std::vector<double> solveAndFitOdd(
				const BufferedImage<aberration::OddData>& data,
				int n_bands,
				const std::vector<double>& initialCoeffs,
				double pixelSize,
				const std::string& prefix,
				bool writeImages);
};

#endif
