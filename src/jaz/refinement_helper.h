#ifndef REFINEMENT_HELPER_H
#define REFINEMENT_HELPER_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <vector>

class Projector;

class RefinementHelper
{
    public:

        static void drawFSC(const MetaDataTable* mdt, std::vector<double>& dest1D,
                            Image<RFLOAT>& dest, double thresh = 0.143);

        static void computeSNR(const MetaDataTable* mdt, Image<RFLOAT>& dest, double eps = 1e-15);
        static void computeSigInvSq(const MetaDataTable* mdt, const std::vector<double>& signalPow,
                                    Image<RFLOAT>& dest, double eps = 1e-15);

        static Image<RFLOAT> correlation(const Image<Complex>& prediction,
                                         const Image<Complex>& observation);

        static Image<RFLOAT> correlation(const std::vector<Image<Complex> >& prediction,
                                         const std::vector<Image<Complex> >& observation);

        static void addToQR(
                const Image<Complex>& prediction, const Image<Complex>& observation,
                Image<Complex>& q, Image<RFLOAT>& r);

        static void addToPQR(
                const Image<Complex>& prediction,
                const Image<Complex>& observation,
                Image<RFLOAT>& p, Image<Complex>& q, Image<RFLOAT>& r);

        static double squaredDiff(
                const Image<Complex>& prediction, const Image<Complex>& observation,
                CTF& ctf, RFLOAT angpix, const Image<RFLOAT>& weight);

        static double squaredDiff(
                const std::vector<Image<Complex>>& predictions,
                const std::vector<Image<Complex>>& observations,
                CTF& ctf, RFLOAT angpix, const Image<RFLOAT>& weight);
};

#endif
