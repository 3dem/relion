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

#ifndef MAGNIFICATION_REFINEMENT_H
#define MAGNIFICATION_REFINEMENT_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/single_particle/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/math/equation2x2.h>
#include <vector>


class MagnificationHelper
{
	public:

		static Matrix2D<RFLOAT> polarToMatrix(
						double scaleMajor = 1.0,
						double scaleMinor = 1.0,
						double angleDeg = 0.0);

		static void matrixToPolar(
						const Matrix2D<RFLOAT>& mat,
						RFLOAT& scaleMajor,
						RFLOAT& scaleMinor,
						RFLOAT& angleDeg);

		static void matrixToPolar(
						const gravis::d2Matrix& mat,
						RFLOAT& scaleMajor,
						RFLOAT& scaleMinor,
						RFLOAT& angleDeg);

		// deprecated: use the other one!
		static void updateScaleFreq(
						const Image<Complex>& prediction,
						const Volume<gravis::t2Vector<Complex>>& predGradient,
						const Image<Complex>& observation,
						CTF& ctf, double angpix,
						Volume<Equation2x2>& eqs,
						bool do_ctf_padding = false);

		static void updateScale(
						const RawImage<fComplex>& prediction,
						const RawImage<gravis::t2Vector<fComplex>>& predGradient,
						const RawImage<fComplex>& observation,
						const RawImage<float>& freqWeights,
						const RawImage<float>& doseWeights,
						CTF& ctf, double angpix,
						RawImage<Equation2x2>& eqs);

		static void updateScaleReal(
						const Image<Complex>& prediction,
						const Image<Complex>& observation,
						const Image<RFLOAT>& snr,
						CTF& ctf, double angpix,
						Volume<Equation2x2>& eqs,
						bool do_ctf_padding = false);

		// deprecated: use the other one!
		static void solvePerPixel(
						const Volume<Equation2x2>& eqs,
						Image<RFLOAT>& vx, Image<RFLOAT>& vy);

		static BufferedImage<double> solvePerPixel(
						const RawImage<Equation2x2>& eqs);

		// deprecated: use the other one!
		static Matrix2D<RFLOAT> solveLinearlyFreq(
						const Volume<Equation2x2>& eqs,
						const Image<RFLOAT>& snr,
						Image<RFLOAT>& vx, Image<RFLOAT>& vy);

		static gravis::d2Matrix solveLinearly(
						const RawImage<Equation2x2>& eqs);

		static void readEQs(std::string path, Volume<Equation2x2>& eqs);
		static void writeEQs(const Volume<Equation2x2>& eqs, std::string path);

		static void updatePowSpec(
						const Image<Complex>& prediction,
						const Image<Complex>& observation,
						CTF& ctf, double angpix,
						Image<RFLOAT>& powSpecPred,
						Image<RFLOAT>& powSpecObs,
						bool do_ctf_padding = false);

		static void adaptAstigmatism(
						const std::vector<Matrix2D<RFLOAT>>& dMs,
						std::vector<MetaDataTable>& partMdts,
						bool perParticle, ObservationModel* obsModel);
};

#endif
