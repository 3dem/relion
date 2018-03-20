#ifndef MAGNIFICATION_REFINEMENT_H
#define MAGNIFICATION_REFINEMENT_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <vector>

class Equation2x2
{
    public:

        Equation2x2();

        double Axx, Axy, Ayy, bx, by;
};

class MagnificationRefinement
{
    public:

        static void updateScaleFreq( const Image<Complex>& prediction,
                                     const Image<Complex>& observation,
                                     CTF& ctf, RFLOAT angpix,
                                     Volume<Equation2x2>& eqs);

        static void updateScaleReal( const Image<Complex>& prediction,
                                     const Image<Complex>& observation,
                                     const Image<RFLOAT>& snr,
                                     CTF& ctf, RFLOAT angpix,
                                     Volume<Equation2x2>& eqs);

        static void solvePerPixel( const Volume<Equation2x2>& eqs,
                                   Image<RFLOAT>& vx, Image<RFLOAT>& vy);

        static void solveLinearlyFreq( const Volume<Equation2x2>& eqs,
                                       const Image<RFLOAT>& snr,
                                       Image<RFLOAT>& mat,
                                       Image<RFLOAT>& vx, Image<RFLOAT>& vy);

        static void readEQs(std::string path, Volume<Equation2x2>& eqs);
        static void writeEQs(const Volume<Equation2x2>& eqs, std::string path);

        static void updatePowSpec(const Image<Complex>& prediction,
                              const Image<Complex>& observation,
                              CTF& ctf, RFLOAT angpix,
                              Image<RFLOAT>& powSpecPred,
                              Image<RFLOAT>& powSpecObs);


};

#endif
