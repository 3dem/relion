#include "lbfgs.h"
#include <src/error.h>

static pthread_mutex_t lib_lbfgs_mutex = PTHREAD_MUTEX_INITIALIZER;

std::vector<double> LBFGS::optimize(
    const std::vector<double> &initial,
    const DifferentiableOptimization &opt,
    bool verbose, int max_iters, double epsilon)
{
    const int N = initial.size();

    lbfgsfloatval_t fx;
    lbfgsfloatval_t* m_x = lbfgs_malloc(N);

    if (m_x == NULL)
    {
        REPORT_ERROR("LBFGS::optimize: Failed to allocate a memory block for variables.\n");
    }

    for (int i = 0; i < N; i++)
    {
        m_x[i] = initial[i];
    }

    LibLbfgsAdapter adapter(opt, N, verbose);

    int ret;

    lbfgs_parameter_t param;

    lbfgs_parameter_init(&param);

    param.max_iterations = max_iters;
    param.epsilon = epsilon;

    pthread_mutex_lock(&lib_lbfgs_mutex);
    {
        ret = lbfgs(N, m_x, &fx, evaluate, progress, &adapter, &param);
    }
    pthread_mutex_unlock(&lib_lbfgs_mutex);

    if (verbose)
    {
        std::cout << "L-BFGS optimization terminated with status code = " << ret << "\n";
        std::cout << "  fx = " << fx << "\n";
    }

    std::vector<double> out(N);

    for (int i = 0; i < N; i++)
    {
        out[i] = m_x[i];
    }

    lbfgs_free(m_x);

    return out;
}

void LBFGS::test()
{
    RosenbrockBanana rb;
    std::vector<double> initial(2);
    initial[0] = 3.0;
    initial[0] = 1.0;

    std::vector<double> x0 = optimize(initial, rb, false);

    std::cout << "should be close to 1, 1: " << x0[0] << ", " << x0[1] << "\n";
}

lbfgsfloatval_t LBFGS::evaluate(
    void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g,
    const int n, const lbfgsfloatval_t step)
{
    LibLbfgsAdapter* adapter = (LibLbfgsAdapter*) instance;

    return adapter->evaluate(x, g, n, step);
}

int LBFGS::progress(
    void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step, int n, int k, int ls)
{
    LibLbfgsAdapter* adapter = (LibLbfgsAdapter*) instance;

    return adapter->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
}

LBFGS::LibLbfgsAdapter::LibLbfgsAdapter(
        const DifferentiableOptimization &opt, int n, bool verbose)
:   opt(opt),
    n(n),
    verbose(verbose),
    x_vec(n),
    grad_vec(n)
{
    tempStorage = opt.allocateTempStorage();
}

lbfgsfloatval_t LBFGS::LibLbfgsAdapter::evaluate(
        const lbfgsfloatval_t *x, lbfgsfloatval_t *g,
        const int n, const lbfgsfloatval_t step)
{
    for (int i = 0; i < n; i++)
    {
        x_vec[i] = x[i];
    }

    double fx = opt.f(x_vec, tempStorage);

    opt.grad(x_vec, grad_vec, tempStorage);

    for (int i = 0; i < n; i++)
    {
        g[i] = grad_vec[i];
    }

    return fx;
}

int LBFGS::LibLbfgsAdapter::progress(
    const lbfgsfloatval_t *x, const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step,
    int n, int k, int ls)
{
    if (verbose)
    {
        std::cout << k << ": " << fx << "\n";
    }
}
