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
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <vector>

#include "equation2x2.h"

class MagnificationHelper
{
	public:
		
		static void updateScaleFreq( 
						const Image<Complex>& prediction,
						const Image<Complex>& observation,
						CTF& ctf, double angpix,
						Volume<Equation2x2>& eqs);
		
		static void updateScaleReal( 
						const Image<Complex>& prediction,
						const Image<Complex>& observation,
						const Image<RFLOAT>& snr,
						CTF& ctf, double angpix,
						Volume<Equation2x2>& eqs);
		
		static void solvePerPixel( 
						const Volume<Equation2x2>& eqs,
						Image<RFLOAT>& vx, Image<RFLOAT>& vy);
		
		static void solveLinearlyFreq( 
						const Volume<Equation2x2>& eqs,
						const Image<RFLOAT>& snr,
						Image<RFLOAT>& mat,
						Image<RFLOAT>& vx, Image<RFLOAT>& vy);
		
		static void readEQs(std::string path, Volume<Equation2x2>& eqs);
		static void writeEQs(const Volume<Equation2x2>& eqs, std::string path);
		
		static void updatePowSpec(
						const Image<Complex>& prediction,
						const Image<Complex>& observation,
						CTF& ctf, double angpix,
						Image<RFLOAT>& powSpecPred,
						Image<RFLOAT>& powSpecObs);
};

#endif
