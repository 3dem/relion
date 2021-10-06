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

#include "gp_motion_fit.h"

#include <src/jaz/single_particle/interpolation.h>
#include <src/jaz/math/Gaussian_process.h>

#include <omp.h>

using namespace gravis;

GpMotionFit::GpMotionFit(
		const std::vector<std::vector<Image<double>>>& correlation,
		double cc_pad,
		double sig_vel_px, double sig_div_px, double sig_acc_px,
		int maxDims,
		const std::vector<d2Vector>& positions,
		const std::vector<d2Vector>& perFrameOffsets,
		int threads, bool expKer)
	:
	  expKer(expKer),
	  pc(correlation.size()),
	  fc(correlation[0].size()),
	  threads(threads),
	  cc_pad(cc_pad),
	  sig_vel_px(sig_vel_px),
	  sig_div_px(sig_div_px),
	  sig_acc_px(sig_acc_px),
	  correlation(correlation),
	  positions(positions),
	  perFrameOffsets(perFrameOffsets)
{
	Matrix2D<double> A(pc,pc);

	const double sv2 = sig_vel_px * sig_vel_px;
	const double sd2 = sig_div_px * sig_div_px;
	const double sd1 = sig_div_px;

	for (int i = 0; i < pc; i++)
	for (int j = i; j < pc; j++)
	{
		const double dd = (positions[i] - positions[j]).norm2();
		const double k = sv2 * (expKer? exp(-sqrt(dd/sd2)) : exp(-0.5*dd/sd1));
		A(i,j) = k;
		A(j,i) = k;
	}

	GaussianProcess::Basis defBasis = GaussianProcess::getBasis(A, maxDims, 1e-10);

	dc = defBasis.eigenvalues.size();

	basis = Matrix2D<RFLOAT>(pc,dc);

	for (int d = 0; d < dc; d++)
	{
		for (int p = 0; p < pc; p++)
		{
			basis(p,d) = (RFLOAT) defBasis.eigenvectors[dc*p + d];
		}
	}

	eigenVals = std::vector<double>(dc);

	for (int d = 0; d < dc; d++)
	{
		eigenVals[d] = (RFLOAT) defBasis.eigenvalues[d];
	}
}


double GpMotionFit::f(const std::vector<double> &x) const
{
	std::vector<std::vector<d2Vector>> pos(pc, std::vector<d2Vector>(fc));

	paramsToPos(x, pos);

	const int pad = 512;
	std::vector<double> e_t(pad*threads, 0.0);

#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
	{
		int t = omp_get_thread_num();

		for (int f = 0; f < fc; f++)
		{
			e_t[pad*t] -= Interpolation::cubicXY(
						correlation[p][f],
						cc_pad * (pos[p][f].x + perFrameOffsets[f].x),
						cc_pad * (pos[p][f].y + perFrameOffsets[f].y),
						0, 0, true);
		}
	}

#pragma omp parallel for num_threads(threads)
	for (int f = 0; f < fc-1; f++)
	{
		int t = omp_get_thread_num();

		for (int d = 0; d < dc; d++)
		{
			const double cx = x[2*(pc + dc*f + d)    ];
			const double cy = x[2*(pc + dc*f + d) + 1];

			e_t[pad*t] += cx*cx + cy*cy;
		}
	}

	if (sig_acc_px > 0.0)
	{
#pragma omp parallel for num_threads(threads)
		for (int f = 0; f < fc-2; f++)
		{
			int t = omp_get_thread_num();

			for (int d = 0; d < dc; d++)
			{
				const double cx0 = x[2*(pc + dc*f + d)    ];
				const double cy0 = x[2*(pc + dc*f + d) + 1];
				const double cx1 = x[2*(pc + dc*(f+1) + d)    ];
				const double cy1 = x[2*(pc + dc*(f+1) + d) + 1];

				const double dcx = cx1 - cx0;
				const double dcy = cy1 - cy0;

				e_t[pad*t] += eigenVals[d]*(dcx*dcx + dcy*dcy) / (sig_acc_px*sig_acc_px);
			}
		}
	}

	double e_tot = 0.0;

	for (int t = 0; t < threads; t++)
	{
		e_tot += e_t[pad*t];
	}

	return e_tot;
}

double GpMotionFit::f(const std::vector<double> &x, void* tempStorage) const
{
	if (tempStorage == 0) return f(x);

	TempStorage* ts = (TempStorage*) tempStorage;

	paramsToPos(x, ts->pos);

	for (int t = 0; t < threads; t++)
	{
		ts->e_t[ts->pad*t] = 0.0;
	}

#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
	{
		int t = omp_get_thread_num();

		for (int f = 0; f < fc; f++)
		{
			const double epf = Interpolation::cubicXY(correlation[p][f],
													  cc_pad * (ts->pos[p][f].x + perFrameOffsets[f].x),
													  cc_pad * (ts->pos[p][f].y + perFrameOffsets[f].y),
													  0, 0, true);
			
			ts->e_t[ts->pad*t] -= epf;
		}
	}

#pragma omp parallel for num_threads(threads)
	for (int f = 0; f < fc-1; f++)
	{
		int t = omp_get_thread_num();

		for (int d = 0; d < dc; d++)
		{
			const double cx = x[2*(pc + dc*f + d)    ];
			const double cy = x[2*(pc + dc*f + d) + 1];

			ts->e_t[ts->pad*t] += cx*cx + cy*cy;
		}
	}

	if (sig_acc_px > 0.0)
	{
#pragma omp parallel for num_threads(threads)
		for (int f = 0; f < fc-2; f++)
		{
			int t = omp_get_thread_num();

			for (int d = 0; d < dc; d++)
			{
				const double cx0 = x[2*(pc + dc*f + d)    ];
				const double cy0 = x[2*(pc + dc*f + d) + 1];
				const double cx1 = x[2*(pc + dc*(f+1) + d)    ];
				const double cy1 = x[2*(pc + dc*(f+1) + d) + 1];

				const double dcx = cx1 - cx0;
				const double dcy = cy1 - cy0;

				ts->e_t[ts->pad*t] += eigenVals[d]*(dcx*dcx + dcy*dcy) / (sig_acc_px*sig_acc_px);
			}
		}
	}

	double e_tot = 0.0;

	for (int t = 0; t < threads; t++)
	{
		e_tot += ts->e_t[ts->pad*t];
	}

	return e_tot;
}

void GpMotionFit::grad(const std::vector<double> &x,
					   std::vector<double> &gradDest) const
{
	std::vector<std::vector<d2Vector>> pos(pc, std::vector<d2Vector>(fc));
	paramsToPos(x, pos);

	int pad = 512;
	std::vector<std::vector<d2Vector>> ccg_pf(pc, std::vector<d2Vector>(fc + pad));

#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
	{
		for (int f = 0; f < fc; f++)
		{
			d2Vector vr = Interpolation::cubicXYgrad(
						correlation[p][f],
						cc_pad * (pos[p][f].x + perFrameOffsets[f].x),
						cc_pad * (pos[p][f].y + perFrameOffsets[f].y),
						0, 0, true);

			ccg_pf[p][f] = d2Vector(vr.x, vr.y);
		}
	}

	std::vector<std::vector<double>> gradDestT(threads, std::vector<double>(gradDest.size()+pad, 0.0));

#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
		for (int f = 0; f < fc; f++)
		{
			int t = omp_get_thread_num();

			gradDestT[t][2*p  ] -= ccg_pf[p][f].x;
			gradDestT[t][2*p+1] -= ccg_pf[p][f].y;
		}

#pragma omp parallel for num_threads(threads)
	for (int d = 0; d < dc; d++)
		for (int p = 0; p < pc; p++)
		{
			int t = omp_get_thread_num();

			d2Vector g(0.0, 0.0);

			const double bpd = basis(p,d);

			for (int f = fc-2; f >= 0; f--)
			{
				g.x += bpd * ccg_pf[p][f+1].x;
				g.y += bpd * ccg_pf[p][f+1].y;

				gradDestT[t][2*(pc + dc*f + d)  ] -= g.x;
				gradDestT[t][2*(pc + dc*f + d)+1] -= g.y;
			}
		}

#pragma omp parallel for num_threads(threads)
	for (int f = 0; f < fc-1; f++)
		for (int d = 0; d < dc; d++)
		{
			int t = omp_get_thread_num();

			gradDestT[t][2*(pc + dc*f + d)  ] += 2.0 * x[2*(pc + dc*f + d)  ];
			gradDestT[t][2*(pc + dc*f + d)+1] += 2.0 * x[2*(pc + dc*f + d)+1];
		}

	if (sig_acc_px > 0.0)
	{
		const double sa2 = sig_acc_px*sig_acc_px;

#pragma omp parallel for num_threads(threads)
		for (int f = 0; f < fc-2; f++)
			for (int d = 0; d < dc; d++)
			{
				int t = omp_get_thread_num();

				const double cx0 = x[2*(pc + dc*f + d)    ];
				const double cy0 = x[2*(pc + dc*f + d) + 1];
				const double cx1 = x[2*(pc + dc*(f+1) + d)    ];
				const double cy1 = x[2*(pc + dc*(f+1) + d) + 1];

				const double dcx = cx1 - cx0;
				const double dcy = cy1 - cy0;

				//e_tot += eigenVals[d]*(dcx*dcx + dcy*dcy) / (sig_acc_px*sig_acc_px);

				gradDestT[t][2*(pc + dc*f + d)  ] -= 2.0 * eigenVals[d] * dcx / sa2;
				gradDestT[t][2*(pc + dc*f + d)+1] -= 2.0 * eigenVals[d] * dcy / sa2;
				gradDestT[t][2*(pc + dc*(f+1) + d)  ] += 2.0 * eigenVals[d] * dcx / sa2;
				gradDestT[t][2*(pc + dc*(f+1) + d)+1] += 2.0 * eigenVals[d] * dcy / sa2;
			}
	}

	for (int i = 0; i < gradDest.size(); i++)
	{
		gradDest[i] = 0.0;
	}

	for (int t = 0; t < threads; t++)
		for (int i = 0; i < gradDest.size(); i++)
		{
			gradDest[i] += gradDestT[t][i];
		}
}

void GpMotionFit::grad(const std::vector<double> &x,
					   std::vector<double> &gradDest,
					   void* tempStorage) const
{
	if (tempStorage == 0) return grad(x, gradDest);

	TempStorage* ts = (TempStorage*) tempStorage;

	paramsToPos(x, ts->pos);

#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
	{
		for (int f = 0; f < fc; f++)
		{
			d2Vector vr = Interpolation::cubicXYgrad(
						correlation[p][f],
						cc_pad * (ts->pos[p][f].x + perFrameOffsets[f].x),
						cc_pad * (ts->pos[p][f].y + perFrameOffsets[f].y),
						0, 0, true);

			ts->ccg_pf[p][f] = d2Vector(vr.x, vr.y);
		}
	}

#pragma omp parallel for num_threads(threads)
	for (int t = 0; t < threads; t++)
		for (int i = 0; i < gradDest.size(); i++)
		{
			ts->gradDestT[t][i] = 0.0;
		}

#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
		for (int f = 0; f < fc; f++)
		{
			int t = omp_get_thread_num();

			ts->gradDestT[t][2*p  ] -= ts->ccg_pf[p][f].x;
			ts->gradDestT[t][2*p+1] -= ts->ccg_pf[p][f].y;
		}

#pragma omp parallel for num_threads(threads)
	for (int d = 0; d < dc; d++)
		for (int p = 0; p < pc; p++)
		{
			int t = omp_get_thread_num();

			d2Vector g(0.0, 0.0);

			const double bpd = basis(p,d);

			for (int f = fc-2; f >= 0; f--)
			{
				g.x += bpd * ts->ccg_pf[p][f+1].x;
				g.y += bpd * ts->ccg_pf[p][f+1].y;

				ts->gradDestT[t][2*(pc + dc*f + d)  ] -= g.x;
				ts->gradDestT[t][2*(pc + dc*f + d)+1] -= g.y;
			}
		}

#pragma omp parallel for num_threads(threads)
	for (int f = 0; f < fc-1; f++)
		for (int d = 0; d < dc; d++)
		{
			int t = omp_get_thread_num();

			ts->gradDestT[t][2*(pc + dc*f + d)  ] += 2.0 * x[2*(pc + dc*f + d)  ];
			ts->gradDestT[t][2*(pc + dc*f + d)+1] += 2.0 * x[2*(pc + dc*f + d)+1];
		}

	if (sig_acc_px > 0.0)
	{
		const double sa2 = sig_acc_px*sig_acc_px;

#pragma omp parallel for num_threads(threads)
		for (int f = 0; f < fc-2; f++)
			for (int d = 0; d < dc; d++)
			{
				int t = omp_get_thread_num();

				const double cx0 = x[2*(pc + dc*f + d)    ];
				const double cy0 = x[2*(pc + dc*f + d) + 1];
				const double cx1 = x[2*(pc + dc*(f+1) + d)    ];
				const double cy1 = x[2*(pc + dc*(f+1) + d) + 1];

				const double dcx = cx1 - cx0;
				const double dcy = cy1 - cy0;

				//e_tot += eigenVals[d]*(dcx*dcx + dcy*dcy) / (sig_acc_px*sig_acc_px);

				ts->gradDestT[t][2*(pc + dc*f + d)  ] -= 2.0 * eigenVals[d] * dcx / sa2;
				ts->gradDestT[t][2*(pc + dc*f + d)+1] -= 2.0 * eigenVals[d] * dcy / sa2;
				ts->gradDestT[t][2*(pc + dc*(f+1) + d)  ] += 2.0 * eigenVals[d] * dcx / sa2;
				ts->gradDestT[t][2*(pc + dc*(f+1) + d)+1] += 2.0 * eigenVals[d] * dcy / sa2;
			}
	}

	for (int i = 0; i < gradDest.size(); i++)
	{
		gradDest[i] = 0.0;
	}

	for (int t = 0; t < threads; t++)
		for (int i = 0; i < gradDest.size(); i++)
		{
			gradDest[i] += ts->gradDestT[t][i];
		}
}

void *GpMotionFit::allocateTempStorage() const
{
	TempStorage* ts = new TempStorage;

	const int pad = 512;
	const int parCt = 2*(pc + dc*(fc-1));

	ts->pad = pad;
	ts->pos = std::vector<std::vector<d2Vector>>(pc, std::vector<d2Vector>(fc + pad));
	ts->ccg_pf = std::vector<std::vector<d2Vector>>(pc, std::vector<d2Vector>(fc + pad));
	ts->gradDestT = std::vector<std::vector<double>>(threads, std::vector<double>(parCt + pad, 0.0));
	ts->e_t = std::vector<double>(pad*threads, 0.0);

	return ts;
}

void GpMotionFit::deallocateTempStorage(void* ts) const
{
	if (ts) delete (TempStorage*) ts;
}

void GpMotionFit::paramsToPos(
		const std::vector<double>& x,
		std::vector<std::vector<d2Vector>>& pos) const
{
#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
	{
		d2Vector pp(x[2*p], x[2*p+1]);

		for (int f = 0; f < fc; f++)
		{
			pos[p][f] = pp;

			if (f < fc-1)
			{
				d2Vector vel(0.0, 0.0);

				for (int d = 0; d < dc; d++)
				{
					const double cx = x[2*(pc + dc*f + d)    ];
					const double cy = x[2*(pc + dc*f + d) + 1];

					vel.x += cx * basis(p,d);
					vel.y += cy * basis(p,d);
				}

				pp += vel;
			}
		}
	}
}

void GpMotionFit::posToParams(
		const std::vector<std::vector<d2Vector>>& pos,
		std::vector<double>& x) const
{
	x.resize(2*(pc + dc*(fc-1)));

	for (int p = 0; p < pc; p++)
	{
		x[2*p]   = pos[p][0].x;
		x[2*p+1] = pos[p][0].y;
	}

	for (int f = 0; f < fc-1; f++)
		for (int d = 0; d < dc; d++)
		{
			d2Vector c(0.0, 0.0);

			for (int p = 0; p < pc; p++)
			{
				d2Vector v = pos[p][f+1] - pos[p][f];

				c.x += v.x * basis(p,d);
				c.y += v.y * basis(p,d);
			}

			x[2*(pc + dc*f + d)  ] = c.x/eigenVals[d];
			x[2*(pc + dc*f + d)+1] = c.y/eigenVals[d];
		}
}
