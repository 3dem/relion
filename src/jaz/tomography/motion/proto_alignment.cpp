#include "proto_alignment.h"
#include <src/jaz/tomography/data_set.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/ctf.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/util/log.h>
#include <src/jaz/math/Tait_Bryan_angles.h>
#include <omp.h>

using namespace gravis;


ProtoAlignment::ProtoAlignment(
	const std::vector<BufferedImage<double>>& CCs,
	const std::vector<gravis::d4Matrix>& frameProj, 
	DataSet* dataSet, 
	const std::vector<int>& partIndices,
	const std::vector<BufferedImage<fComplex>>& referenceFS,
	bool constParticles,
	bool constAngles,
	bool constShifts,
	int maxRange, 
	d3Vector tomoCentre,
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
	  fc(frameProj.size()),
	  pc(partIndices.size()),
	  maxRange(maxRange),
	  tomoCentre(tomoCentre),
	  num_threads(num_threads),
	  paddingFactor(paddingFactor)
{	
	initialPos.resize(pc);
	
	for (int p = 0; p < pc; p++)
	{
		initialPos[p] = dataSet->getPosition(partIndices[p]);
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
	
	const int fs = getFrameStride();
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		const int t = omp_get_thread_num();
		
		double phi, theta, psi, dx, dy;		
		readParams(x, fs*f, phi, theta, psi, dx, dy);
		
		const d4Matrix Q = TaitBryan::anglesToMatrix4(phi, theta, psi);		
		
		d4Matrix P = plusCentre * Q * minusCentre * frameProj[f];
		
		P(0,3) += dx;
		P(1,3) += dy;
				
		for (int p = 0; p < pc; p++)
		{
			const d3Vector off = constParticles? 
						d3Vector(0.0, 0.0, 0.0) : 
						d3Vector(x[fs*fc + 3*p], x[fs*fc + 3*p+1], x[fs*fc + 3*p+2]);
			
			const d4Vector pos4(initialPos[p] + off);
			
			const d4Vector dp  = P * pos4 - frameProj[f] * d4Vector(initialPos[p]);
			
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
	
	const int fs = getFrameStride();
	
	std::vector<double> grad_par(thread_stride * num_threads, 0.0);
			
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		const int t = omp_get_thread_num();
		
		double phi, theta, psi, dx, dy;		
		readParams(x, fs*f, phi, theta, psi, dx, dy);
		
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
		
		
		for (int p = 0; p < pc; p++)
		{
			const d3Vector off = constParticles? 
						d3Vector(0.0, 0.0, 0.0) : 
						d3Vector(x[fs*fc + 3*p], x[fs*fc + 3*p+1], x[fs*fc + 3*p+2]);
			
			const d4Vector pos4(initialPos[p] + off);
			
			const d4Vector dp  = P * pos4 - frameProj[f] * d4Vector(initialPos[p]);
						
			const double dx_img = (dp.x + maxRange) * paddingFactor;
			const double dy_img = (dp.y + maxRange) * paddingFactor;
			
			const d4Vector dp_phi   = P_phi   * pos4;
			const d4Vector dp_theta = P_theta * pos4;
			const d4Vector dp_psi   = P_psi   * pos4;
			
			d2Vector g = -((double)paddingFactor) * Interpolation::cubicXYGrad_clip(
						CCs[p], dx_img, dy_img, f);
			
			if (constAngles)
			{
				if (!constShifts)
				{
					grad_par[t*thread_stride + fs*f    ]  +=  g.x;
					grad_par[t*thread_stride + fs*f + 1]  +=  g.y;
				}
			}
			else
			{
				if (constShifts)
				{
					grad_par[t*thread_stride + fs*f    ]  +=  dp_phi.x   * g.x  +  dp_phi.y   * g.y;
					grad_par[t*thread_stride + fs*f + 1]  +=  dp_theta.x * g.x  +  dp_theta.y * g.y;
					grad_par[t*thread_stride + fs*f + 2]  +=  dp_psi.x   * g.x  +  dp_psi.y   * g.y;
				}
				else
				{
					grad_par[t*thread_stride + fs*f    ]  +=  dp_phi.x   * g.x  +  dp_phi.y   * g.y;
					grad_par[t*thread_stride + fs*f + 1]  +=  dp_theta.x * g.x  +  dp_theta.y * g.y;
					grad_par[t*thread_stride + fs*f + 2]  +=  dp_psi.x   * g.x  +  dp_psi.y   * g.y;
					grad_par[t*thread_stride + fs*f + 3]  +=  g.x;
					grad_par[t*thread_stride + fs*f + 4]  +=  g.y;
				}
			}
			
			if (!constParticles)
			{
				grad_par[t*thread_stride + fs*fc + 3*p    ]  +=  P(0,0) * g.x  +  P(1,0) * g.y;
				grad_par[t*thread_stride + fs*fc + 3*p + 1]  +=  P(0,1) * g.x  +  P(1,1) * g.y;
				grad_par[t*thread_stride + fs*fc + 3*p + 2]  +=  P(0,2) * g.x  +  P(1,2) * g.y;
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

std::vector<d4Matrix> ProtoAlignment::getProjections(const std::vector<double> &x) const
{
	std::vector<d4Matrix> out(fc);
	
	const int fs = getFrameStride();
	
	for (int f = 0; f < fc; f++)
	{	
		double phi, theta, psi, dx, dy;		
		readParams(x, fs*f, phi, theta, psi, dx, dy);
		
		const d4Matrix Q = TaitBryan::anglesToMatrix4(phi, theta, psi);		
		
		out[f] = plusCentre * Q * minusCentre * frameProj[f];
		
		out[f](0,3) += dx;
		out[f](1,3) += dy;
	}
		
	return out;
}

void ProtoAlignment::shiftParticles(
		const std::vector<double> &x,
		const std::vector<int>& partIndices, 
		DataSet *target) const
{
	if (constParticles) return;
	
	const int fs = getFrameStride();
	
	for (int p = 0; p < pc; p++)
	{
		const d3Vector origin = initialPos[p] + d3Vector(
					x[fs*fc + 3*p], x[fs*fc + 3*p+1], x[fs*fc + 3*p+2]);
		
		target->moveParticleTo(partIndices[p], origin);
	}
}

int ProtoAlignment::getParamCount()
{
	const int fs = getFrameStride();
	
	int out = fs * fc;
	
	if (!constParticles) out += 3*pc;
	
	return out;
}

void ProtoAlignment::protocol(const std::vector<double> &x, int frame0, int frame1) const
{
	double sum(0.0);
	
	const int fs = getFrameStride();
	
	for (int f = frame0; f <= frame1; f++)
	{
		double phi, theta, psi, dx, dy;		
		readParams(x, fs*f, phi, theta, psi, dx, dy);
		
		const d4Matrix Q = TaitBryan::anglesToMatrix4(phi, theta, psi);		
		
		d4Matrix P = plusCentre * Q * minusCentre * frameProj[f];
		
		P(0,3) += dx;
		P(1,3) += dy;
		
		std::cout << "frame " << f << ": " << std::endl;
		std::cout << "P = \n" << P << std::endl;
				
		for (int p = 0; p < pc; p++)
		{
			const d3Vector off = constParticles? 
						d3Vector(0.0, 0.0, 0.0) : 
						d3Vector(x[fs*fc + 3*p], x[fs*fc + 3*p+1], x[fs*fc + 3*p+2]);
			
			const d4Vector pos4(initialPos[p] + off);
			
			std::cout << "   " << p << ": " << pos4 << std::endl;
			std::cout << "      " << initialPos[p] << " + " << off << std::endl;
			
			const d4Vector dp  = P * pos4 - frameProj[f] * d4Vector(initialPos[p]);
			
			std::cout << "      " << dp << std::endl;
			
			double val = Interpolation::cubicXY_clip(
						CCs[p], dp.x + maxRange, dp.y + maxRange, f);
			
			sum -= val / (fc * (double) pc);
			
			std::cout << "      " << val << " -> " << sum << std::endl;
		}
		
		std::cout << "sum = " << sum << std::endl;
	}
}

std::vector<BufferedImage<double>> ProtoAlignment::drawShiftedCCs(const std::vector<double> &x) const
{
	const int d = CCs[0].xdim;
	std::vector<BufferedImage<double>> out(pc);
	
	for (int p = 0; p < pc; p++)
	{
		out[p] = BufferedImage<double>(d,d,fc);
	}
	
	const int fs = getFrameStride();
	
	BufferedImage<dComplex> CCsFS(d/2+1,d,fc);
	
		
	for (int p = 0; p < pc; p++)
	{
		const d3Vector off = constParticles? 
					d3Vector(0.0, 0.0, 0.0) : 
					d3Vector(x[fs*fc + 3*p], x[fs*fc + 3*p+1], x[fs*fc + 3*p+2]);
		
		const d4Vector pos4(initialPos[p] + off);
		
		std::vector<d2Vector> posInNewImg(fc);
		
		for (int f = 0; f < fc; f++)
		{
			double phi, theta, psi, dx, dy;		
			readParams(x, fs*f, phi, theta, psi, dx, dy);
			
			const d4Matrix Q = TaitBryan::anglesToMatrix4(phi, theta, psi);		
			
			d4Matrix P = plusCentre * Q * minusCentre * frameProj[f];
			
			P(0,3) += dx;
			P(1,3) += dy;
			
			const d4Vector dp = P * pos4 - frameProj[f] * d4Vector(initialPos[p]);
			
			posInNewImg[f] = d2Vector(
						dp.x * paddingFactor,
						dp.y * paddingFactor);
		}
			
		StackHelper::FourierTransformStack(CCs[p], CCsFS, true, num_threads);
		StackHelper::shiftStack(CCsFS, posInNewImg, CCsFS, true, num_threads);
		StackHelper::inverseFourierTransformStack(CCsFS, out[p], false, num_threads);
	}
	
	return out;
}

void ProtoAlignment::report(int iteration, double cost, const std::vector<double> &x) const
{
	Log::updateProgress(iteration);
}
