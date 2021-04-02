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

#include "lbfgs.h"
#include <src/error.h>

std::vector<double> LBFGS::optimize(
		const std::vector<double> &initial,
		const DifferentiableOptimization &opt,
		int verbosity, int max_iters,
		double epsilon, double xtol)
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

	void* tempStorage = opt.allocateTempStorage();

	LibLbfgsAdapter adapter(opt, tempStorage, N, verbosity > 0);

	int ret(LBFGSERR_UNKNOWNERROR);

	lbfgs_parameter_t param;

	lbfgs_parameter_init(&param);

	param.max_iterations = max_iters;
	param.epsilon = epsilon;
	param.xtol = xtol;

	#pragma omp critical(LBGFS)
	{
		ret = lbfgs(N, m_x, &fx, evaluate, progress, &adapter, &param);
	}

	if (verbosity > 1)
	{
		std::cout << "\n\nL-BFGS optimization terminated with status code " << ret << ": ";
		std::cout << decodeStatus(ret) << "\n";
		std::cout << "  fx = " << fx << std::endl;
	}

	std::vector<double> out(N);

	for (int i = 0; i < N; i++)
	{
		out[i] = m_x[i];
	}

	lbfgs_free(m_x);

	opt.deallocateTempStorage(tempStorage);

	return out;
}

std::vector<double> LBFGS::optimize(
		const std::vector<double> &initial,
		const FastDifferentiableOptimization &opt,
		int verbosity, int max_iters,
		double epsilon, double xtol)
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

	FastLibLbfgsAdapter adapter(opt, N, verbosity > 0);

	int ret(LBFGSERR_UNKNOWNERROR);

	lbfgs_parameter_t param;

	lbfgs_parameter_init(&param);

	param.max_iterations = max_iters;
	param.epsilon = epsilon;
	param.xtol = xtol;

	//#pragma omp critical(LBGFS)
	{
		ret = lbfgs(N, m_x, &fx, evaluate_fast, progress_fast, &adapter, &param);
	}

	if (verbosity > 1)
	{
		std::cout << "\n\nL-BFGS optimization terminated with status code " << ret << ": ";
		std::cout << decodeStatus(ret) << "\n";
		std::cout << "  fx = " << fx << std::endl;
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

lbfgsfloatval_t LBFGS::evaluate_fast(
		void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g,
		const int n, const lbfgsfloatval_t step)
{
	FastLibLbfgsAdapter* adapter = (FastLibLbfgsAdapter*) instance;

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

int LBFGS::progress_fast(
		void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g,
		const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step, int n, int k, int ls)
{
	FastLibLbfgsAdapter* adapter = (FastLibLbfgsAdapter*) instance;

	return adapter->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
}

// ------------------------------ //

LBFGS::LibLbfgsAdapter::LibLbfgsAdapter(
		const DifferentiableOptimization &opt,
		void *tempStorage, int n, bool verbose)
	:   opt(opt),
	  n(n),
	  verbose(verbose),
	  x_vec(n),
	  grad_vec(n),
	  tempStorage(tempStorage)
{    
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
		for (int i = 0; i < n; i++)
		{
			x_vec[i] = x[i];
		}
		
		opt.report(k, fx, x_vec);
	}
	
	return 0;
}

// ------------------------------ //

LBFGS::FastLibLbfgsAdapter::FastLibLbfgsAdapter(
		const FastDifferentiableOptimization &opt,
		int n, bool verbose)
	:   opt(opt),
	  n(n),
	  verbose(verbose),
	  x_vec(n),
	  grad_vec(n)
{
}

lbfgsfloatval_t LBFGS::FastLibLbfgsAdapter::evaluate(
		const lbfgsfloatval_t *x, lbfgsfloatval_t *g,
		const int n, const lbfgsfloatval_t step)
{
	for (int i = 0; i < n; i++)
	{
		x_vec[i] = x[i];
	}

	double fx = opt.gradAndValue(x_vec, grad_vec);

	for (int i = 0; i < n; i++)
	{
		g[i] = grad_vec[i];
	}

	return fx;
}

int LBFGS::FastLibLbfgsAdapter::progress(
		const lbfgsfloatval_t *x, const lbfgsfloatval_t *g,
		const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm,
		const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step,
		int n, int k, int ls)
{
	if (verbose)
	{
		for (int i = 0; i < n; i++)
		{
			x_vec[i] = x[i];
		}

		opt.report(k, fx, x_vec);
	}

	return 0;
}

// ------------------------------ //

std::string LBFGS::decodeStatus(int ret)
{
	switch(ret)
	{
		case LBFGS_SUCCESS: return "LBFGS_SUCCESS";
		case LBFGS_STOP: return "LBFGS_STOP";
		case LBFGS_ALREADY_MINIMIZED: return "LBFGS_ALREADY_MINIMIZED";
		case LBFGSERR_UNKNOWNERROR: return "LBFGSERR_UNKNOWNERROR";
		case LBFGSERR_LOGICERROR: return "LBFGSERR_LOGICERROR";
		case LBFGSERR_OUTOFMEMORY: return "LBFGSERR_OUTOFMEMORY";
		case LBFGSERR_CANCELED: return "LBFGSERR_CANCELED";
		case LBFGSERR_INVALID_N: return "LBFGSERR_INVALID_N";
		case LBFGSERR_INVALID_N_SSE: return "LBFGSERR_INVALID_N_SSE";
		case LBFGSERR_INVALID_X_SSE: return "LBFGSERR_INVALID_X_SSE";
		case LBFGSERR_INVALID_EPSILON: return "LBFGSERR_INVALID_EPSILON";
		case LBFGSERR_INVALID_TESTPERIOD: return "LBFGSERR_INVALID_TESTPERIOD";
		case LBFGSERR_INVALID_DELTA: return "LBFGSERR_INVALID_DELTA";
		case LBFGSERR_INVALID_LINESEARCH: return "LBFGSERR_INVALID_LINESEARCH";
		case LBFGSERR_INVALID_MINSTEP: return "LBFGSERR_INVALID_MINSTEP";
		case LBFGSERR_INVALID_MAXSTEP: return "LBFGSERR_INVALID_MAXSTEP";
		case LBFGSERR_INVALID_FTOL: return "LBFGSERR_INVALID_FTOL";
		case LBFGSERR_INVALID_WOLFE: return "LBFGSERR_INVALID_WOLFE";
		case LBFGSERR_INVALID_GTOL: return "LBFGSERR_INVALID_GTOL";
		case LBFGSERR_INVALID_XTOL: return "LBFGSERR_INVALID_XTOL";
		case LBFGSERR_INVALID_MAXLINESEARCH: return "LBFGSERR_INVALID_MAXLINESEARCH";
		case LBFGSERR_INVALID_ORTHANTWISE: return "LBFGSERR_INVALID_ORTHANTWISE";
		case LBFGSERR_INVALID_ORTHANTWISE_START: return "LBFGSERR_INVALID_ORTHANTWISE_START";
		case LBFGSERR_INVALID_ORTHANTWISE_END: return "LBFGSERR_INVALID_ORTHANTWISE_END";
		case LBFGSERR_OUTOFINTERVAL: return "LBFGSERR_OUTOFINTERVAL";
		case LBFGSERR_INCORRECT_TMINMAX: return "LBFGSERR_INCORRECT_TMINMAX";
		case LBFGSERR_ROUNDING_ERROR: return "LBFGSERR_ROUNDING_ERROR";
		case LBFGSERR_MINIMUMSTEP: return "LBFGSERR_MINIMUMSTEP";
		case LBFGSERR_MAXIMUMSTEP: return "LBFGSERR_MAXIMUMSTEP";
		case LBFGSERR_MAXIMUMLINESEARCH: return "LBFGSERR_MAXIMUMLINESEARCH";
		case LBFGSERR_MAXIMUMITERATION: return "LBFGSERR_MAXIMUMITERATION";
		case LBFGSERR_WIDTHTOOSMALL: return "LBFGSERR_WIDTHTOOSMALL";
		case LBFGSERR_INVALIDPARAMETERS: return "LBFGSERR_INVALIDPARAMETERS";
		case LBFGSERR_INCREASEGRADIENT: return "LBFGSERR_INCREASEGRADIENT";
	}
	
	return "UNKNOWN_STATUS_CODE";
}
