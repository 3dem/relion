#ifndef DEFOCUS_REFINEMENT_H
#define DEFOCUS_REFINEMENT_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <vector>

class AstigmatismOptimization : public Optimization
{
    public:

        AstigmatismOptimization(
                const Image<Complex>& prediction,
                const Image<Complex>& observation,
                const Image<RFLOAT>& weight,
                const CTF& ctf0,
                RFLOAT angpix);

        double f(const std::vector<double>& x) const;

    private:

        const Image<Complex>& prediction, observation;
        const Image<RFLOAT>& weight;
        const CTF& ctf0;
        RFLOAT angpix;
};

class AstigmatismOptimizationAcc : public Optimization
{
    public:

        AstigmatismOptimizationAcc(
                const Image<Complex>& prediction,
                const Image<Complex>& observation,
                const Image<RFLOAT>& weight,
                const CTF& ctf0,
                RFLOAT angpix,
                RFLOAT phiScale = 10);

        double f(const std::vector<double>& x) const;

        double getU(const std::vector<double> &x);
        double getV(const std::vector<double> &x);
        double getPhi(const std::vector<double> &x);

        std::vector<double> getInitialParams();

    private:

        Image<Complex> data;
        const CTF& ctf0;
        RFLOAT angpix;
        RFLOAT phiScale;
};

class DefocusRefinement
{
    public:

        static RFLOAT findDefocus1D(const Image<Complex>& prediction,
            const Image<Complex>& observation,
            const Image<RFLOAT>& weight,
            const CTF& ctf0, RFLOAT angpix,
            RFLOAT* destU, RFLOAT* destV,
            RFLOAT range = 1000.0, int steps = 11,
            int recDepth = 2, RFLOAT recScale = 10.0);

        static void findAstigmatismNM(
            const Image<Complex>& prediction,
            const Image<Complex>& observation,
            const Image<RFLOAT>& weight,
            const CTF& ctf0, RFLOAT angpix,
            RFLOAT* destU, RFLOAT* destV, RFLOAT* destPhi);

        static std::vector<gravis::d2Vector> diagnoseDefocus(const Image<Complex>& prediction,
            const Image<Complex>& observation,
            const Image<RFLOAT>& weight,
            const CTF& ctf0, double angpix,
            double range, int steps, int threads);
};

#endif
