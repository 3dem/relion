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

#ifndef TWO_HYPERPARAMETER_PROBLEM_H
#define TWO_HYPERPARAMETER_PROBLEM_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/gravis/t2Vector.h>

class MotionParamEstimator;

class TwoHyperParameterProblem : public Optimization
{
    public:

        TwoHyperParameterProblem(
            MotionParamEstimator& motionParamEstimator,
            double s_acc);

        double f(const std::vector<double>& x, void* tempStorage) const;
        void report(int iteration, double cost, const std::vector<double>& x) const;

        static gravis::d2Vector problemToMotion(const std::vector<double>& x);
        static std::vector<double> motionToProblem(gravis::d2Vector vd);

    protected:

        MotionParamEstimator& motionParamEstimator;
        double s_acc;
};

#endif
