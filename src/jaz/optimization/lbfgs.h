#ifndef LBFGS_OPT_H
#define LBFGS_OPT_H

#include <vector>
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
