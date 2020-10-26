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

#ifndef MOTION_PARAM_ESTIMATOR_H
#define MOTION_PARAM_ESTIMATOR_H

#include "alignment_set.h"

#include <src/image.h>
#include <src/time.h>

#include <src/jaz/gravis/t4Vector.h>


//#define TIMING 1

#ifdef TIMING
    #define RCTIC(timer,label) (timer.tic(label))
    #define RCTOC(timer,label) (timer.toc(label))
#else
    #define RCTIC(timer,label)
    #define RCTOC(timer,label)
#endif

class MotionEstimator;
class ReferenceMap;
class ObservationModel;
class ParFourierTransformer;

class MotionParamEstimator
{
    public:

        MotionParamEstimator();


        void read(IOParser& parser, int argc, char *argv[]);

        void init(int verb, int nr_omp_threads, bool debug,
                  std::string outPath, int fc,
                  const std::vector<MetaDataTable>& allMdts,
                  MotionEstimator* motionEstimator,
                  ReferenceMap* reference,
                  ObservationModel* obsModel);

        void run();

        bool anythingToDo();

        // to be used by instances of OptimizationProblem
        void evaluateParams(const std::vector<gravis::d3Vector>& sig_vals,
                            std::vector<double>& TSCs);

        static const double velScale, divScale, accScale;


    protected:

            bool paramsRead, ready;

            AlignmentSet<float> alignmentSet;

            // read from cmd. line:
            bool estim2, estim3;
            int minParticles, maxRange, maxIters, seed, group;
            double sV, sD, sA;
            double iniStep, conv;
			double align_frac, eval_frac;
            double k_cutoff, k_eval;

            // set at init:
            std::vector<MetaDataTable> mdts;

            MotionEstimator* motionEstimator;
            ObservationModel* obsModel;
            ReferenceMap* reference;

            int fc, k_out, verb, nr_omp_threads, s_ref, s, sh;
			bool allGroups;
			
            bool debug;
			std::string outPath;

            #ifdef TIMING
                Timer paramTimer;
                int timeSetup, timeOpt, timeEval;
            #endif


        gravis::d4Vector estimateTwoParamsNM(
                double sig_v_0, double sig_d_0, double sig_acc,
                double inStep, double conv, int maxIters);

        gravis::d4Vector estimateThreeParamsNM(
                double sig_v_0, double sig_d_0, double sig_a_0,
                double inStep, double conv, int maxIters);


        void prepAlignment();
};

#endif
