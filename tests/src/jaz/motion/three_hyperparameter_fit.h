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

#ifndef THREE_HYPERPARAMETER_PROBLEM_H
#define THREE_HYPERPARAMETER_PROBLEM_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/gravis/t3Vector.h>

class MotionParamEstimator;

class ThreeHyperParameterProblem : public Optimization
{
    public:

        ThreeHyperParameterProblem(
            MotionParamEstimator& motionParamEstimator);

        double f(const std::vector<double>& x, void* tempStorage) const;
        void report(int iteration, double cost, const std::vector<double>& x) const;

        static gravis::d3Vector problemToMotion(const std::vector<double>& x);
        static std::vector<double> motionToProblem(gravis::d3Vector vd);

    protected:

        MotionParamEstimator& motionParamEstimator;

        static double accThresh, accEps;
};

#endif
