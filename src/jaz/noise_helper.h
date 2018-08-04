#ifndef NOISE_HELPER
#define NOISE_HELPER

#include <vector>
#include "src/projector.h"
#include "src/ctf.h"

class NoiseHelper
{
    public:

        static Image<RFLOAT> predictCCNoise(Projector& prj, double sigma2, double nsamples_ppp, int max_nsamples, int nangles,
                Image<RFLOAT> &dmgWeight, CTF ctf0, double defocusMu, double defocusSigma, double angpix, int thread_num = 1);

        static Image<RFLOAT> visualize(std::vector<double>);

        static std::vector<double> radialAverage(Image<RFLOAT>& map, bool half);
        static Image<RFLOAT> radialMap(std::vector<double>& radAvg, bool centered);
		
        static std::vector<Complex> radialAverage(Image<Complex>& map, bool skipAxes);
        static Image<Complex> radialMap(std::vector<Complex>& radAvg);		
		
		static std::vector<std::pair<double,double>> radialAverageAndStdDevFFTW(Image<RFLOAT>& map);
		

        static std::vector<double> radialWeight(int w, int h, bool half);
        static std::vector<double> fill(Image<RFLOAT>& confusion, double lambda, int iterations);
        static Image<RFLOAT> normalize(const Image<RFLOAT> &confusion);

        static void testVariance(Image<RFLOAT> img);
        static void testColorVariance(Image<RFLOAT> img, std::vector<double> sig2);
        static void testParseval();




};

#endif
