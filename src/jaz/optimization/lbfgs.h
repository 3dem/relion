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

#ifndef LBFGS_OPT_H
#define LBFGS_OPT_H

#include <vector>
#include <string>
#include "optimization.h"
#include <src/jaz/lbfgs/lbfgs.h>

class LBFGS
{
    public:

        static std::vector<double> optimize(
            const std::vector<double>& initial,
            const DifferentiableOptimization& opt,
            bool verbose = false,
            int max_iters = 0,
            double epsilon = 1e-5);

        static void test();

    protected:
		
		static std::string translateError(int ret);

        static lbfgsfloatval_t evaluate(
            void *instance,
            const lbfgsfloatval_t *x,
            lbfgsfloatval_t *g,
            const int n,
            const lbfgsfloatval_t step);

        static int progress(
            void *instance,
            const lbfgsfloatval_t *x,
            const lbfgsfloatval_t *g,
            const lbfgsfloatval_t fx,
            const lbfgsfloatval_t xnorm,
            const lbfgsfloatval_t gnorm,
            const lbfgsfloatval_t step,
            int n,
            int k,
            int ls);

        class LibLbfgsAdapter
        {
            public:

                LibLbfgsAdapter(
                    const DifferentiableOptimization& opt,
                    void* tempStorage, int n, bool verbose);

                const DifferentiableOptimization& opt;
                int n;
                bool verbose;
                std::vector<double> x_vec, grad_vec;
                void* tempStorage;


                lbfgsfloatval_t evaluate(
                    const lbfgsfloatval_t *x,
                    lbfgsfloatval_t *g,
                    const int n,
                    const lbfgsfloatval_t step);

                int progress(
                    const lbfgsfloatval_t *x,
                    const lbfgsfloatval_t *g,
                    const lbfgsfloatval_t fx,
                    const lbfgsfloatval_t xnorm,
                    const lbfgsfloatval_t gnorm,
                    const lbfgsfloatval_t step,
                    int n,
                    int k,
                    int ls);
        };
};

#endif
