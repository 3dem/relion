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
