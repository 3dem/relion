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
		
		static Matrix2D<RFLOAT> polarToMatrix(
						double scaleMajor = 1.0, 
						double scaleMinor = 1.0,
						double angleDeg = 0.0);
		
		static void matrixToPolar(
						const Matrix2D<RFLOAT>& mat,
						RFLOAT& scaleMajor, 
						RFLOAT& scaleMinor,
						RFLOAT& angleDeg);
		
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
		
		static Matrix2D<RFLOAT> solveLinearlyFreq( 
						const Volume<Equation2x2>& eqs,
						const Image<RFLOAT>& snr,
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
