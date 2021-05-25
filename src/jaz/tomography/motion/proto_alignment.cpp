#include "proto_alignment.h"
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/ctf.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/util/log.h>
#include <src/jaz/math/Tait_Bryan_angles.h>
#include <omp.h>
#include <iomanip>

using namespace gravis;


ProtoAlignment::ProtoAlignment(
	const std::vector<BufferedImage<double>>& CCs,
	const std::vector<gravis::d4Matrix>& frameProj, 
	const ParticleSet& dataSet,
	const std::vector<ParticleIndex>& partIndices,
	const std::vector<BufferedImage<fComplex>>& referenceFS,
	bool constParticles,
	bool constAngles,
	bool constShifts,
	bool doAnisotropy,
	bool perTiltAnisotropy,
	int maxRange, 
	d3Vector tomoCentre,
	int progressBarOffset,
	int num_threads,
	double paddingFactor)
	:	
	  CCs(CCs),
	  frameProj(frameProj),
	  dataSet(dataSet),
	  partIndices(partIndices),
	  referenceFS(referenceFS),
	  constParticles(constParticles),
	  constAngles(constAngles),
	  constShifts(constShifts),
	  doAnisotropy(doAnisotropy),
	  perTiltAnisotropy(perTiltAnisotropy),
	  devMode(false),
	  fc(frameProj.size()),
	  pc(partIndices.size()),
	  maxRange(maxRange),
	  tomoCentre(tomoCentre),
	  progressBarOffset(progressBarOffset),
	  num_threads(num_threads),
	  paddingFactor(paddingFactor)
{	
	initialPos.resize(pc);
	
	for (int p = 0; p < pc; p++)
	{
		initialPos[p] = dataSet.getPosition(partIndices[p]);
	}
	
	minusCentre = d4Matrix(
			1, 0, 0, -tomoCentre.x, 
			0, 1, 0, -tomoCentre.y, 
			0, 0, 1, -tomoCentre.z, 
			0, 0, 0, 1 );
	
	plusCentre = d4Matrix(
			1, 0, 0, tomoCentre.x, 
			0, 1, 0, tomoCentre.y, 
			0, 0, 1, tomoCentre.z, 
			0, 0, 0, 1 );
}

double ProtoAlignment::f(const std::vector<double> &x, void *tempStorage) const
{
	const int data_pad = 512;
	std::vector<double> cost_par(num_threads * data_pad, 0.0);

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		const int t = omp_get_thread_num();
		
		double phi, theta, psi, dx, dy, skew, y_scale;
		readTiltParameters(x, f, phi, theta, psi, dx, dy, skew, y_scale);
		
		const d4Matrix Q = TaitBryan::anglesToMatrix4(phi, theta, psi);
		
		d4Matrix P = plusCentre * Q * minusCentre * frameProj[f];
		
		P(0,3) += dx;
		P(1,3) += dy;

		const d4Matrix A(
			1.0, skew, 0.0, 0.0,
			0.0, y_scale, 0.0, 0.0,
			0.0, 0.0, 1.0, 0.0,
			0.0, 0.0, 0.0, 1.0 );

		const d4Matrix AP = A * P;

		const int part_off = getParticleDataOffset();

				
		for (int p = 0; p < pc; p++)
		{
			const d3Vector off = constParticles? 
						d3Vector(0.0, 0.0, 0.0) : 
						d3Vector(x[part_off + 3*p], x[part_off + 3*p+1], x[part_off + 3*p+2]);
			
			const d4Vector pos4(initialPos[p] + off);

			const d4Vector p1 = AP * pos4;
			const d4Vector p0 = frameProj[f] * d4Vector(initialPos[p]);
			
			const d4Vector dp  = p1 - p0;
			
			const double dx_img = (dp.x + maxRange) * paddingFactor;
			const double dy_img = (dp.y + maxRange) * paddingFactor;
			
			double val = Interpolation::cubicXY_clip(
						CCs[p], dx_img, dy_img, f);
			
			cost_par[t * data_pad] -= val;
		}
	}
		
	double cost(0.0);
	
	for (int t = 0; t < num_threads; t++)
	{
		cost += cost_par[t * data_pad] / (fc * (double) pc);
	}
	
	return cost;
}

void ProtoAlignment::grad(const std::vector<double> &x, std::vector<double> &gradDest, void *tempStorage) const
{
	const int data_pad = 512;
	const int xs = x.size();
	const int thread_stride = xs + data_pad;
	
	std::vector<double> grad_par(thread_stride * num_threads, 0.0);
			
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		const int t = omp_get_thread_num();
		
		double phi, theta, psi, dx, dy, skew, y_scale;
		readTiltParameters(x, f, phi, theta, psi, dx, dy, skew, y_scale);
		
		const d4Matrix Q = TaitBryan::anglesToMatrix4(phi, theta, psi);	
		t4Vector<gravis::d3Matrix> dQ = TaitBryan::anglesToMatrixAndDerivatives(phi, theta, psi);
		
		d4Matrix Q_phi(dQ[0]);
		d4Matrix Q_theta(dQ[1]);
		d4Matrix Q_psi(dQ[2]);
		
		Q_phi(3,3) = 0.0;
		Q_theta(3,3) = 0.0;
		Q_psi(3,3) = 0.0;
		
		const d4Matrix centProj = minusCentre * frameProj[f];
		
		d4Matrix P = plusCentre * Q * centProj;
		
		d4Matrix P_phi   = plusCentre * Q_phi   * centProj;
		d4Matrix P_theta = plusCentre * Q_theta * centProj;
		d4Matrix P_psi   = plusCentre * Q_psi   * centProj;
		
		P(0,3) += dx;
		P(1,3) += dy;

		d4Matrix A(
			1.0, skew, 0.0, 0.0,
			0.0, y_scale, 0.0, 0.0,
			0.0, 0.0, 1.0, 0.0,
			0.0, 0.0, 0.0, 1.0 );

		d4Matrix A_skew(
			0.0, 1.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0 );

		d4Matrix A_y_scale(
			0.0, 0.0, 0.0, 0.0,
			0.0, 1.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0,
			0.0, 0.0, 0.0, 0.0 );

		const d4Matrix AP = A * P;


		const int part_off = getParticleDataOffset();
		
		for (int p = 0; p < pc; p++)
		{
			const d3Vector off = constParticles? 
						d3Vector(0.0, 0.0, 0.0) : 
						d3Vector(x[part_off + 3*p], x[part_off + 3*p+1], x[part_off + 3*p+2]);
			
			const d4Vector pos4(initialPos[p] + off);

			const d4Vector p1_0 = P * pos4;
			const d4Vector p0 = frameProj[f] * d4Vector(initialPos[p]);

			const d4Vector p1 = A * p1_0;
			
			const d4Vector dp  = p1 - p0;

			const double dx_img = (dp.x + maxRange) * paddingFactor;
			const double dy_img = (dp.y + maxRange) * paddingFactor;
			
			const d4Vector dp_phi   = A * P_phi   * pos4;
			const d4Vector dp_theta = A * P_theta * pos4;
			const d4Vector dp_psi   = A * P_psi   * pos4;
			const d4Vector dp_dx      = d4Vector(1.0, 0.0, 0.0, 0.0);
			const d4Vector dp_dy      = d4Vector(skew, y_scale, 0.0, 0.0);
			const d4Vector dp_skew    = A_skew    * p1_0;
			const d4Vector dp_y_scale = A_y_scale * p1_0;
			
			d2Vector g = -((double)paddingFactor) * Interpolation::cubicXYGrad_clip(
						CCs[p], dx_img, dy_img, f);

			int offset = getTiltDataOffset(f);
			const int t0 = t * thread_stride;

			if (!constAngles)
			{
				grad_par[t0 + offset    ]  +=  dp_phi.x   * g.x  +  dp_phi.y   * g.y;
				grad_par[t0 + offset + 1]  +=  dp_theta.x * g.x  +  dp_theta.y * g.y;
				grad_par[t0 + offset + 2]  +=  dp_psi.x   * g.x  +  dp_psi.y   * g.y;

				offset += 3;
			}

			if (!constShifts)
			{
				grad_par[t0 + offset    ]  +=  dp_dx.x * g.x + dp_dx.y * g.y;
				grad_par[t0 + offset + 1]  +=  dp_dy.x * g.x + dp_dy.y * g.y;

				offset += 2;
			}

			if (doAnisotropy)
			{
				const double skew_grad    = dp_skew.x    * g.x  +  dp_skew.y    * g.y;
				const double y_scale_grad = dp_y_scale.x * g.x  +  dp_y_scale.y * g.y;


				if (perTiltAnisotropy)
				{
					grad_par[t0 + offset    ]  +=  skew_grad;
					grad_par[t0 + offset + 1]  +=  y_scale_grad;
				}
				else
				{
					grad_par[t0    ]  +=  skew_grad;
					grad_par[t0 + 1]  +=  y_scale_grad;
				}
			}
			
			if (!constParticles)
			{
				grad_par[t0 + part_off + 3*p    ]  +=  AP(0,0) * g.x  +  AP(1,0) * g.y;
				grad_par[t0 + part_off + 3*p + 1]  +=  AP(0,1) * g.x  +  AP(1,1) * g.y;
				grad_par[t0 + part_off + 3*p + 2]  +=  AP(0,2) * g.x  +  AP(1,2) * g.y;
			}
		}
	}

	for (int i = 0; i < xs; i++)
	{
		gradDest[i] = 0.0;
	}
	
	for (int t = 0; t < num_threads; t++)
	for (int i = 0; i < xs; i++)
	{
		gradDest[i] += grad_par[t*thread_stride + i] / (fc * (double) pc);
	}
}

std::vector<d4Matrix> ProtoAlignment::getProjections(
		const std::vector<double> &x) const
{
	std::vector<d4Matrix> out(fc);
	
	for (int f = 0; f < fc; f++)
	{
		double phi, theta, psi, dx, dy, skew, y_scale;
		readTiltParameters(x, f, phi, theta, psi, dx, dy, skew, y_scale);
		
		const d4Matrix Q = TaitBryan::anglesToMatrix4(phi, theta, psi);

		d4Matrix A(
			1.0, skew, 0.0, 0.0,
			0.0, y_scale, 0.0, 0.0,
			0.0, 0.0, 1.0, 0.0,
			0.0, 0.0, 0.0, 1.0 );

		out[f] = plusCentre * Q * minusCentre * frameProj[f];
		
		out[f](0,3) += dx;
		out[f](1,3) += dy;

		out[f] = A * out[f];
	}
		
	return out;
}

void ProtoAlignment::shiftParticles(
		const std::vector<double> &x,
		const std::vector<ParticleIndex>& partIndices, 
		ParticleSet& target) const
{
	if (constParticles) return;
	
	for (int p = 0; p < pc; p++)
	{
		const int part_off = getParticleDataOffset();

		const d3Vector origin = initialPos[p] + d3Vector(
			x[part_off + 3*p], x[part_off + 3*p+1], x[part_off + 3*p+2]);
		
		target.moveParticleTo(partIndices[p], origin);
	}
}

std::vector<d3Vector> ProtoAlignment::getParticlePositions(
		const std::vector<double>& x) const
{
	std::vector<d3Vector> out(pc);

	for (int p = 0; p < pc; p++)
	{
		const int part_off = getParticleDataOffset();

		out[p] = initialPos[p] + d3Vector(
			x[part_off + 3*p], x[part_off + 3*p+1], x[part_off + 3*p+2]);
	}

	return out;
}

int ProtoAlignment::getParamCount()
{
	int out = getParticleDataOffset();
	if (!constParticles) out += 3*pc;
	
	return out;
}

std::vector<BufferedImage<double>> ProtoAlignment::drawShiftedCCs(const std::vector<double> &x) const
{
	const int d = CCs[0].xdim;
	std::vector<BufferedImage<double>> out(pc);
	
	for (int p = 0; p < pc; p++)
	{
		out[p] = BufferedImage<double>(d,d,fc);
	}

	
	BufferedImage<dComplex> CCsFS(d/2+1,d,fc);

	const int part_off = getParticleDataOffset();
		
	for (int p = 0; p < pc; p++)
	{
		const d3Vector off = constParticles? 
					d3Vector(0.0, 0.0, 0.0) : 
					d3Vector(x[part_off + 3*p], x[part_off + 3*p+1], x[part_off + 3*p+2]);
		
		const d4Vector pos4(initialPos[p] + off);
		
		std::vector<d2Vector> posInNewImg(fc);
		
		for (int f = 0; f < fc; f++)
		{
			double phi, theta, psi, dx, dy, skew, y_scale;
			readTiltParameters(x, f, phi, theta, psi, dx, dy, skew, y_scale);
			
			const d4Matrix Q = TaitBryan::anglesToMatrix4(phi, theta, psi);
			
			d4Matrix P = plusCentre * Q * minusCentre * frameProj[f];
			
			P(0,3) += dx;
			P(1,3) += dy;

			d4Matrix A(
				1.0, skew, 0.0, 0.0,
				0.0, y_scale, 0.0, 0.0,
				0.0, 0.0, 1.0, 0.0,
				0.0, 0.0, 0.0, 1.0 );
			
			const d4Vector dp = A * P * pos4 - frameProj[f] * d4Vector(initialPos[p]);
			
			posInNewImg[f] = d2Vector(
						dp.x * paddingFactor,
						dp.y * paddingFactor);
		}
			
		NewStackHelper::FourierTransformStack(CCs[p], CCsFS, true, num_threads);
		NewStackHelper::shiftStack(CCsFS, posInNewImg, CCsFS, true, num_threads);
		NewStackHelper::inverseFourierTransformStack(CCsFS, out[p], false, num_threads);
	}
	
	return out;
}

void ProtoAlignment::report(int iteration, double cost, const std::vector<double> &x) const
{
	if (devMode)
	{
		int prec = (int)(log(iteration) / log(10));
		int step = pow(10, prec);

		if (prec == 0 || iteration % step == 0)
		{
			std::cout << iteration << " \t " << std::setw(11) << std::setprecision(10) << cost << std::endl;
		}
	}
	else
	{
		Log::updateProgress(progressBarOffset + iteration);
	}
}

std::vector<double> ProtoAlignment::createInitial()
{
	std::vector<double> out(getParamCount(), 0.0);

	if (doAnisotropy)
	{
		if (perTiltAnisotropy)
		{
			int y_scale_off = 1;

			if (!constAngles) y_scale_off += 3;
			if (!constShifts) y_scale_off += 2;

			for (int f = 0; f < fc; f++)
			{
				out[getTiltDataOffset(f) + y_scale_off] = 1.0;
			}
		}
		else
		{
			out[1] = 1.0;
		}
	}

	return out;
}
