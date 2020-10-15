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
