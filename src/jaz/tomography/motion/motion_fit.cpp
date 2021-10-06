#include "motion_fit.h"
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/ctf.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/math/Tait_Bryan_angles.h>
#include <src/jaz/math/Gaussian_process.h>
#include <src/jaz/util/zio.h>
#include <omp.h>

using namespace gravis;


MotionFit::MotionFit(
	const std::vector<BufferedImage<double>>& CCs,
	const std::vector<gravis::d4Matrix>& frameProj, 
	ParticleSet& dataSet,
	const std::vector<ParticleIndex>& partIndices,
	MotionParameters motionParameters,
	Settings settings,
	gravis::d3Vector tomoCentre,
	double frameDose,
	double pixelSize,
	double paddingFactor,
	int progressBarOffset,
	int num_threads,
	bool verbose)
	:	
	  CCs(CCs),
	  frameProj(frameProj),
	  dataSet(dataSet),
	  partIndices(partIndices),
	  motionParameters(motionParameters),
	  settings(settings),
	  tomoCentre(tomoCentre),
	  frameDose(frameDose),
	  pixelSize(pixelSize),
	  paddingFactor(paddingFactor),
	  progressBarOffset(progressBarOffset),
	  num_threads(num_threads),
	  verbose(verbose),
	  fc(frameProj.size()),
	  pc(partIndices.size()),
	  maxRange(CCs[0].xdim / (2 * paddingFactor))
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
	
	double sig_vel_px, sig_div_px; 
	
	
	if (settings.params_scaled_by_dose)
	{
		sig_vel_px = frameDose * motionParameters.sig_vel / pixelSize;
		sig_div_px = frameDose * motionParameters.sig_div / pixelSize;
	}
	else
	{
		sig_vel_px = motionParameters.sig_vel / pixelSize;		
		sig_div_px = motionParameters.sig_div / pixelSize;
	}
	
	if (verbose)
	{
		Log::beginSection("Effective motion parameters:");

		Log::print("σ_vel = "+ZIO::itoa(sig_vel_px)+" px/frame");
		Log::print("σ_div = "+ZIO::itoa(sig_div_px)+" px/frame");

		Log::endSection();
	}
	
	GpKernel* kernel(0);
	
	if (settings.sqExpKernel)
	{
		kernel = new SquareExponentialKernel(sig_vel_px, sig_div_px);
	}
	else
	{
		kernel = new ExponentialKernel(sig_vel_px, sig_div_px);
	}
	
	Matrix2D<double> C = GaussianProcess::computeCovariance(initialPos, kernel);
	
	delete kernel;
	
	if (verbose)
	{
		Log::print("Decomposing covariance matrix");
	}
	
	GaussianProcess::Basis defBasis = GaussianProcess::getBasis(C, settings.maxEDs);
	
	deformationBasis = defBasis.eigenvectors;
	deformationLambda = defBasis.eigenvalues;
	
	bc = deformationLambda.size();
	
	if (verbose)
	{
		if (bc == pc)
		{
			Log::print("Keeping all " + ZIO::itoa(bc) + " eigendeformations");
		}
		else
		{
			Log::print("Keeping " + ZIO::itoa(bc) + " out of " + ZIO::itoa(pc) + " eigendeformations");
		}
	}
}

/*
	Parameter Layout:
	
	0:                  ([phi, theta, psi], [dx, dy]) * fc           frame alignment: fs * fc
	fs * fc:             [dx0, dy0, dz0]  *  pc;                     static part. shifts: 3 * pc
	fs * fc + 3 * pc:    [b0x, b0y, b0z][b1x, ..., bBz] * (fc-1);    motion: 3 * bc * (fc - 1)
	
	fs * fc  +  3 * pc  +  3 * bc * (fc - 1)   total
*/

double MotionFit::f(const std::vector<double> &x, void *tempStorage) const
{
	const int data_pad = 512;
	std::vector<double> cost_par(num_threads * data_pad, 0.0);
	
	const int fs = getFrameStride();
	
	std::vector<d4Matrix> P(fc);
	
	for (int f = 0; f < fc; f++)
	{
		double phi, theta, psi, dx, dy;		
		readParams(x, fs*f, phi, theta, psi, dx, dy);
		
		const d4Matrix Q = TaitBryan::anglesToMatrix4(phi, theta, psi);		
		
		P[f] = plusCentre * Q * minusCentre * frameProj[f];
		
		P[f](0,3) += dx;
		P[f](1,3) += dy;
	}
	
	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{
		const int t = omp_get_thread_num();
		
		d3Vector shift = settings.constParticles? 
					d3Vector(0.0, 0.0, 0.0) : 
					d3Vector(x[fs*fc + 3*p], x[fs*fc + 3*p+1], x[fs*fc + 3*p+2]);
		
		for (int f = 0; f < fc; f++)
		{
			const d4Vector pos4(initialPos[p] + shift);
			
			const d4Vector dp  = P[f] * pos4 - frameProj[f] * d4Vector(initialPos[p]);
			
			const double dx_img = (dp.x + maxRange) * paddingFactor;
			const double dy_img = (dp.y + maxRange) * paddingFactor;
			
			double val = -Interpolation::cubicXY_clip(CCs[p], dx_img, dy_img, f);
			
			cost_par[t * data_pad] += val;
			
			if (!settings.constParticles && f < fc-1)
			{
				shift += getPosChange(x, p, f, fs * fc + 3 * pc);
			}
		}
	}
		
	double cost(0.0);
	
	for (int t = 0; t < num_threads; t++)
	{
		cost += cost_par[t * data_pad];
	}
	
	if (!settings.constParticles)
	{
		for (int m = 0; m < fc-1; m++)
		{
			for (int b = 0; b < bc; b++)
			{
				const int i0 = fs*fc + 3*(pc + m*bc + b);
				
				cost += x[i0]*x[i0] + x[i0+1]*x[i0+1] + x[i0+2]*x[i0+2];
			}
		}
	}
	
	return cost;
}

void MotionFit::grad(const std::vector<double> &x, std::vector<double> &gradDest, void *tempStorage) const
{
	const int fs = getFrameStride();
	const int xs = x.size();
	const int data_pad = 512;
	const int step_grad = xs + data_pad;
	const int step_frame = fc + data_pad;

	std::vector<d4Matrix> P(fc), P_phi(fc), P_theta(fc), P_psi(fc);

	for (int f = 0; f < fc; f++)
	{
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

		P[f] = plusCentre * Q * centProj;

		P_phi[f]   = plusCentre * Q_phi   * centProj;
		P_theta[f] = plusCentre * Q_theta * centProj;
		P_psi[f]   = plusCentre * Q_psi   * centProj;

		P[f](0,3) += dx;
		P[f](1,3) += dy;
	}


	std::vector<double> grad_par(step_grad * num_threads, 0.0);
	std::vector<d3Vector> dC_dPos(step_frame * num_threads, d3Vector(0.0, 0.0, 0.0));


	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{

		const int th = omp_get_thread_num();

		d3Vector shift = settings.constParticles?
					d3Vector(0.0, 0.0, 0.0) :
					d3Vector(x[fs*fc + 3*p], x[fs*fc + 3*p+1], x[fs*fc + 3*p+2]);

		for (int f = 0; f < fc; f++)
		{
			const d4Vector pos4(initialPos[p] + shift);

			const d4Vector dp  = P[f] * pos4 - frameProj[f] * d4Vector(initialPos[p]);

			const double dx_img = (dp.x + maxRange) * paddingFactor;
			const double dy_img = (dp.y + maxRange) * paddingFactor;

			const d4Vector dp_phi   = P_phi[f]   * pos4;
			const d4Vector dp_theta = P_theta[f] * pos4;
			const d4Vector dp_psi   = P_psi[f]   * pos4;


			d2Vector g = -((double)paddingFactor) * Interpolation::cubicXYGrad_clip(
						CCs[p], dx_img, dy_img, f);


			if (settings.constAngles)
			{
				if (!settings.constShifts)
				{
					grad_par[th*step_grad + fs*f    ]  +=  g.x;
					grad_par[th*step_grad + fs*f + 1]  +=  g.y;
				}
			}
			else
			{
				if (settings.constShifts)
				{
					grad_par[th*step_grad + fs*f    ]  +=  dp_phi.x   * g.x  +  dp_phi.y   * g.y;
					grad_par[th*step_grad + fs*f + 1]  +=  dp_theta.x * g.x  +  dp_theta.y * g.y;
					grad_par[th*step_grad + fs*f + 2]  +=  dp_psi.x   * g.x  +  dp_psi.y   * g.y;
				}
				else
				{
					grad_par[th*step_grad + fs*f    ]  +=  dp_phi.x   * g.x  +  dp_phi.y   * g.y;
					grad_par[th*step_grad + fs*f + 1]  +=  dp_theta.x * g.x  +  dp_theta.y * g.y;
					grad_par[th*step_grad + fs*f + 2]  +=  dp_psi.x   * g.x  +  dp_psi.y   * g.y;
					grad_par[th*step_grad + fs*f + 3]  +=  g.x;
					grad_par[th*step_grad + fs*f + 4]  +=  g.y;
				}
			}

			if (!settings.constParticles)
			{
				const d3Vector dC_dPos_f(
					P[f](0,0) * g.x  +  P[f](1,0) * g.y,
					P[f](0,1) * g.x  +  P[f](1,1) * g.y,
					P[f](0,2) * g.x  +  P[f](1,2) * g.y);

				dC_dPos[th*step_frame + f] = dC_dPos_f;

				grad_par[th*step_grad + fs*fc + 3*p    ]  +=  dC_dPos_f.x;
				grad_par[th*step_grad + fs*fc + 3*p + 1]  +=  dC_dPos_f.y;
				grad_par[th*step_grad + fs*fc + 3*p + 2]  +=  dC_dPos_f.z;

				if (f < fc-1)
				{
					shift += getPosChange(x, p, f, fs * fc + 3 * pc);
				}
			}
		}

		if (!settings.constParticles)
		{
			for (int m = 0; m < fc-1; m++)
			{
				for (int b = 0; b < bc; b++)
				{
					d3Vector dC_dXm(0.0, 0.0, 0.0);
					const double def = deformationBasis[p*bc + b];

					for (int f = m+1; f < fc; f++)
					{
						dC_dXm += def * dC_dPos[th*step_frame + f];
					}

					const int i0 = th*step_grad + fs*fc + 3*(pc + m*bc + b);

					grad_par[i0    ] += dC_dXm.x;
					grad_par[i0 + 1] += dC_dXm.y;
					grad_par[i0 + 2] += dC_dXm.z;
				}
			}
		}
	}

	for (int i = 0; i < xs; i++)
	{
		gradDest[i] = 0.0;
	}

	for (int th = 0; th < num_threads; th++)
	for (int i = 0; i < xs; i++)
	{
		gradDest[i] += grad_par[th*step_grad + i];
	}

	if (!settings.constParticles)
	{
		for (int m = 0; m < fc-1; m++)
		{
			for (int b = 0; b < bc; b++)
			{
				const int i0 = fs*fc + 3*(pc + m*bc + b);

				gradDest[i0] += 2.0 * (x[i0] + x[i0+1] + x[i0+2]);
			}
		}
	}
}

double MotionFit::gradAndValue(const std::vector<double> &x, std::vector<double> &gradDest) const
{
	const int fs = getFrameStride();
	const int xs = x.size();
	const int data_pad = 512;
	const int step_grad = xs + data_pad;
	const int step_frame = fc + data_pad;

	std::vector<d4Matrix> P(fc), P_phi(fc), P_theta(fc), P_psi(fc);

	for (int f = 0; f < fc; f++)
	{
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

		P[f] = plusCentre * Q * centProj;

		P_phi[f]   = plusCentre * Q_phi   * centProj;
		P_theta[f] = plusCentre * Q_theta * centProj;
		P_psi[f]   = plusCentre * Q_psi   * centProj;

		P[f](0,3) += dx;
		P[f](1,3) += dy;
	}


	std::vector<double> grad_par(step_grad * num_threads, 0.0);
	std::vector<double> val_par(data_pad * num_threads, 0.0);
	std::vector<d3Vector> dC_dPos(step_frame * num_threads, d3Vector(0.0, 0.0, 0.0));


	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{

		const int th = omp_get_thread_num();

		d3Vector shift = settings.constParticles?
					d3Vector(0.0, 0.0, 0.0) :
					d3Vector(x[fs*fc + 3*p], x[fs*fc + 3*p+1], x[fs*fc + 3*p+2]);

		for (int f = 0; f < fc; f++)
		{
			const d4Vector pos4(initialPos[p] + shift);

			const d4Vector dp  = P[f] * pos4 - frameProj[f] * d4Vector(initialPos[p]);


			const double dx_img = (dp.x + maxRange) * paddingFactor;
			const double dy_img = (dp.y + maxRange) * paddingFactor;

			const d4Vector dp_phi   = P_phi[f]   * pos4;
			const d4Vector dp_theta = P_theta[f] * pos4;
			const d4Vector dp_psi   = P_psi[f]   * pos4;


			const d3Vector g = -((double)paddingFactor) *
				Interpolation::cubicXYGradAndValue_clip(CCs[p], dx_img, dy_img, f);

			val_par[th*data_pad] += g.z;


			if (settings.constAngles)
			{
				if (!settings.constShifts)
				{
					grad_par[th*step_grad + fs*f    ]  +=  g.x;
					grad_par[th*step_grad + fs*f + 1]  +=  g.y;
				}
			}
			else
			{
				if (settings.constShifts)
				{
					grad_par[th*step_grad + fs*f    ]  +=  dp_phi.x   * g.x  +  dp_phi.y   * g.y;
					grad_par[th*step_grad + fs*f + 1]  +=  dp_theta.x * g.x  +  dp_theta.y * g.y;
					grad_par[th*step_grad + fs*f + 2]  +=  dp_psi.x   * g.x  +  dp_psi.y   * g.y;
				}
				else
				{
					grad_par[th*step_grad + fs*f    ]  +=  dp_phi.x   * g.x  +  dp_phi.y   * g.y;
					grad_par[th*step_grad + fs*f + 1]  +=  dp_theta.x * g.x  +  dp_theta.y * g.y;
					grad_par[th*step_grad + fs*f + 2]  +=  dp_psi.x   * g.x  +  dp_psi.y   * g.y;
					grad_par[th*step_grad + fs*f + 3]  +=  g.x;
					grad_par[th*step_grad + fs*f + 4]  +=  g.y;
				}
			}

			if (!settings.constParticles)
			{
				const d3Vector dC_dPos_f(
					P[f](0,0) * g.x  +  P[f](1,0) * g.y,
					P[f](0,1) * g.x  +  P[f](1,1) * g.y,
					P[f](0,2) * g.x  +  P[f](1,2) * g.y);

				dC_dPos[th*step_frame + f] = dC_dPos_f;

				grad_par[th*step_grad + fs*fc + 3*p    ]  +=  dC_dPos_f.x;
				grad_par[th*step_grad + fs*fc + 3*p + 1]  +=  dC_dPos_f.y;
				grad_par[th*step_grad + fs*fc + 3*p + 2]  +=  dC_dPos_f.z;

				if (f < fc-1)
				{
					shift += getPosChange(x, p, f, fs * fc + 3 * pc);
				}
			}
		}

		if (!settings.constParticles)
		{
			for (int m = 0; m < fc-1; m++)
			{
				for (int b = 0; b < bc; b++)
				{
					d3Vector dC_dXm(0.0, 0.0, 0.0);
					const double def = deformationBasis[p*bc + b];

					for (int f = m+1; f < fc; f++)
					{
						dC_dXm += def * dC_dPos[th*step_frame + f];
					}

					const int i0 = th*step_grad + fs*fc + 3*(pc + m*bc + b);

					grad_par[i0    ] += dC_dXm.x;
					grad_par[i0 + 1] += dC_dXm.y;
					grad_par[i0 + 2] += dC_dXm.z;
				}
			}
		}
	}

	for (int i = 0; i < xs; i++)
	{
		gradDest[i] = 0.0;
	}

	double cost = 0.0;

	for (int th = 0; th < num_threads; th++)
	{
		cost += val_par[th*data_pad];
	}

	for (int th = 0; th < num_threads; th++)
	for (int i = 0; i < xs; i++)
	{
		gradDest[i] += grad_par[th*step_grad + i];
	}

	if (!settings.constParticles)
	{
		double reg_cost = 0.0;

		for (int m = 0; m < fc-1; m++)
		{
			for (int b = 0; b < bc; b++)
			{
				const int i0 = fs*fc + 3*(pc + m*bc + b);

				gradDest[i0  ] += 2.0 * x[i0  ];
				gradDest[i0+1] += 2.0 * x[i0+1];
				gradDest[i0+2] += 2.0 * x[i0+2];

				reg_cost += x[i0]*x[i0] + x[i0+1]*x[i0+1] + x[i0+2]*x[i0+2];
			}
		}

		cost += reg_cost;
	}

	return cost;
}

std::vector<d4Matrix> MotionFit::getProjections(
		const std::vector<double> &x,
		const std::vector<int>& frameSequence) const
{
	std::vector<d4Matrix> out(fc);
	
	const int fs = getFrameStride();
	
	for (int f = 0; f < fc; f++)
	{	
		double phi, theta, psi, dx, dy;		
		readParams(x, fs*f, phi, theta, psi, dx, dy);
		
		const d4Matrix Q = TaitBryan::anglesToMatrix4(phi, theta, psi);		
		
		const int fa = frameSequence[f];
		
		out[fa] = plusCentre * Q * minusCentre * frameProj[f];
		
		out[fa](0,3) += dx;
		out[fa](1,3) += dy;
	}
		
	return out;
}

std::vector<d3Vector> MotionFit::getParticlePositions(
		const std::vector<double>& x) const
{
	std::vector<d3Vector> out(pc);
	
	const int fs = getFrameStride();
	
	for (int p = 0; p < pc; p++)
	{
		out[p] = initialPos[p] + d3Vector(
			x[fs*fc + 3*p],
			x[fs*fc + 3*p+1],
			x[fs*fc + 3*p+2]);
	}

	return out;
}

Trajectory MotionFit::getTrajectory(
		const std::vector<double> &x, int p,
		const std::vector<int>& frameSequence) const
{
	Trajectory out(fc);
	
	if (settings.constParticles) return out;
	
	const int fs = getFrameStride();

	d3Vector shift(0.0, 0.0, 0.0);
	
	for (int f = 0; f < fc; f++)
	{
		const int fa = frameSequence[f];
		
		out.shifts_Ang[fa] = pixelSize * shift;
		
		if (f < fc-1)
		{
			shift += getPosChange(x, p, f, fs*fc + 3*pc);
		}
	}
	
	return out;
}

std::vector<Trajectory> MotionFit::exportTrajectories(
		const std::vector<double>& x, 
		const ParticleSet& dataSet,
		const std::vector<int>& frameSequence) const
{
	std::vector<Trajectory> out(pc);

	for (int p = 0; p < pc; p++)
	{
		const int pp = partIndices[p].value;

		out[p] = dataSet.motionTrajectories[pp] + getTrajectory(x, p, frameSequence);
	}

	return out;
}

int MotionFit::getParamCount()
{
	const int fs = getFrameStride();
	
	int out = fs * fc;
	
	if (!settings.constParticles) out += 3 * (pc + bc * (fc - 1));
	
	return out;
}

std::vector<BufferedImage<double>> MotionFit::drawShiftedCCs(const std::vector<double> &x) const
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
		d3Vector shift = settings.constParticles? 
					d3Vector(0.0, 0.0, 0.0) : 
					d3Vector(x[fs*fc + 3*p], x[fs*fc + 3*p+1], x[fs*fc + 3*p+2]);
		
		const int offset = fs * fc + 3 * pc;
		
		
		std::vector<d2Vector> posInNewImg(fc);
		
		for (int f = 0; f < fc; f++)
		{
			double phi, theta, psi, dx, dy;		
			readParams(x, fs*f, phi, theta, psi, dx, dy);
			
			const d4Matrix Q = TaitBryan::anglesToMatrix4(phi, theta, psi);		
			
			d4Matrix P = plusCentre * Q * minusCentre * frameProj[f];
			
			P(0,3) += dx;
			P(1,3) += dy;
			
			const d4Vector pos4(initialPos[p] + shift);			
			const d4Vector dp = P * pos4 - frameProj[f] * d4Vector(initialPos[p]);
			
			posInNewImg[f] = d2Vector(
						dp.x * paddingFactor,
						dp.y * paddingFactor);
			
			if (!settings.constParticles && f < fc-1)
			{
				shift += getPosChange(x, p, f, offset);
			}
		}
			
		NewStackHelper::FourierTransformStack(CCs[p], CCsFS, true, num_threads);
		NewStackHelper::shiftStack(CCsFS, posInNewImg, CCsFS, true, num_threads);
		NewStackHelper::inverseFourierTransformStack(CCsFS, out[p], false, num_threads);
	}
	
	return out;
}

Mesh MotionFit::visualiseTrajectories3D(const std::vector<double> &x, double scale)
{
	Mesh out;
	
	std::vector<int> timeSeq(fc);
	
	for (int f = 0; f < fc; f++)
	{
		timeSeq[f] = f;
	}
	
	for (int p = 0; p < pc; p++)
	{
		Trajectory track = getTrajectory(x, p, timeSeq);
		
		for (int f = 0; f < fc-1; f++)
		{
			const d3Vector a = initialPos[p] + scale * track.shifts_Ang[f] / pixelSize;
			const d3Vector b = initialPos[p] + scale * track.shifts_Ang[f+1] / pixelSize;
			
			const double c = (f == 0)? 1.0 : 0.66 - 0.33*(f%2);
						
			MeshBuilder::addColouredBar(a,b, 0.1, 7, dRGB(c,c,c), out);
		}
	}
	
	return out;
}

void MotionFit::visualiseTrajectories2D(
		const std::vector<double>& x,
		double scale,
		const std::string& tomo_name,
		const std::string& file_name_root)
{
	std::vector<int> timeSeq(fc);

	for (int f = 0; f < fc; f++)
	{
		timeSeq[f] = f;
	}

	std::vector<std::string> plot_names {"XY", "XZ", "YZ"};
	std::vector<i2Vector> dim_indices{i2Vector(0,1), i2Vector(0,2), i2Vector(1,2)};

	for (int dim = 0; dim < 3; dim++)
	{
		CPlot2D plot2D(tomo_name + " Motion " + plot_names[dim]);
		plot2D.SetXAxisSize(600);
		plot2D.SetYAxisSize(600);
		plot2D.SetDrawLegend(false);
		plot2D.SetFlipY(true);

		for (int p = 0; p < pc; p++)
		{
			Trajectory track = getTrajectory(x, p, timeSeq);

			// Mark start of each track
			CDataSet start;
			start.SetDrawMarker(true);
			start.SetMarkerSize(8);
			start.SetDatasetColor(0.2,0.5,1.0);

			const d3Vector a = initialPos[p] + scale * track.shifts_Ang[0] / pixelSize;
			const i2Vector di = dim_indices[dim];

			CDataPoint point3(a[di[0]],a[di[1]]);
			start.AddDataPoint(point3);
			plot2D.AddDataSet(start);
		}

		for (int p = 0; p < pc; p++)
		{
			Trajectory track = getTrajectory(x, p, timeSeq);

			CDataSet curve;
			curve.SetDrawMarker(false);
			curve.SetDatasetColor(0.0,0.0,0.0);
			curve.SetLineWidth(0.5);

			for (int f = 0; f < fc; f++)
			{
				const d3Vector a = initialPos[p] + scale * track.shifts_Ang[f] / pixelSize;
				const i2Vector di = dim_indices[dim];

				CDataPoint point(a[di[0]],a[di[1]]);

				curve.AddDataPoint(point);
			}

			plot2D.AddDataSet(curve);
		}

		std::string label_x = plot_names[dim].substr(0,1) +
				" (in pixels; trajectory scaled by " + ZIO::itoa(scale) + ")";

		std::string label_y = plot_names[dim].substr(1,1);

		plot2D.SetXAxisTitle(label_x);
		plot2D.SetYAxisTitle(label_y);

		FileName fn_eps = file_name_root + "_" + plot_names[dim] + ".eps";

		plot2D.OutputPostScriptPlot(fn_eps);
	}
}

void MotionFit::report(int iteration, double cost, const std::vector<double> &x) const
{
	if (verbose)
	{
		Log::updateProgress(progressBarOffset + iteration);
	}
}

void MotionFit::analyseGradient(const std::vector<double>& x, int particle, int frame, double epsilon)
{
	/*
		compute d(pos_pf.xy) / d(params)  ->  vec<d2> dPos_dx
		shift each param and note the change in pos_pf.xy  ->  vec<d2>
		
	*/
	
	
	const int fs = getFrameStride();
	const int xs = x.size();
	const int data_pad = 512;	
	const int step_grad = xs + data_pad;
	const int step_frame = fc + data_pad;
	
	std::vector<d4Matrix> P(fc), P_phi(fc), P_theta(fc), P_psi(fc);
	
	for (int f = 0; f < fc; f++)
	{
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
		
		P[f] = plusCentre * Q * centProj;
		
		P_phi[f]   = plusCentre * Q_phi   * centProj;
		P_theta[f] = plusCentre * Q_theta * centProj;
		P_psi[f]   = plusCentre * Q_psi   * centProj;
		
		P[f](0,3) += dx;
		P[f](1,3) += dy;
	}
	
	
	std::vector<d2Vector> grad_par(step_grad, d2Vector(0.0, 0.0));
	std::vector<d3Vector> dx_dPos(step_frame, d3Vector(0.0, 0.0, 0.0));
	std::vector<d3Vector> dy_dPos(step_frame, d3Vector(0.0, 0.0, 0.0));
	
	
	const int p = particle;		
	const int th = 0;
	
	d3Vector shift = settings.constParticles? 
				d3Vector(0.0, 0.0, 0.0) : 
				d3Vector(x[fs*fc + 3*p], x[fs*fc + 3*p+1], x[fs*fc + 3*p+2]);

	
	for (int f = 0; f < fc; f++)
	{
		const d4Vector pos4(initialPos[p] + shift);
		
		const d4Vector dp_phi   = P_phi[f]   * pos4;
		const d4Vector dp_theta = P_theta[f] * pos4;
		const d4Vector dp_psi   = P_psi[f]   * pos4;
		
		
		if (f == frame)
		{
			if (settings.constAngles)
			{
				if (!settings.constShifts)
				{
					grad_par[th*step_grad + fs*f    ].x  +=  1;
					grad_par[th*step_grad + fs*f + 1].y  +=  1;
				}
			}
			else
			{
				if (settings.constShifts)
				{
					grad_par[th*step_grad + fs*f    ]  +=  d2Vector(dp_phi.x, dp_phi.y);
					grad_par[th*step_grad + fs*f + 1]  +=  d2Vector(dp_theta.x, dp_theta.y);
					grad_par[th*step_grad + fs*f + 2]  +=  d2Vector(dp_psi.x, dp_psi.y);
				}
				else
				{
					grad_par[th*step_grad + fs*f    ]  +=  d2Vector(dp_phi.x, dp_phi.y);
					grad_par[th*step_grad + fs*f + 1]  +=  d2Vector(dp_theta.x, dp_theta.y);
					grad_par[th*step_grad + fs*f + 2]  +=  d2Vector(dp_psi.x, dp_psi.y);
					grad_par[th*step_grad + fs*f + 3].x  +=  1;
					grad_par[th*step_grad + fs*f + 4].y  +=  1;
				}
			}
			
			if (!settings.constParticles)
			{
				const d3Vector dx_dPos_f(
					P[f](0,0),
					P[f](0,1),
					P[f](0,2));
						
				const d3Vector dy_dPos_f(
					P[f](1,0),
					P[f](1,1),
					P[f](1,2));
						
				dx_dPos[th*step_frame + f] = dx_dPos_f;
				dy_dPos[th*step_frame + f] = dy_dPos_f;
						
				grad_par[th*step_grad + fs*fc + 3*p    ]  +=  d2Vector(dx_dPos_f.x, dy_dPos_f.x);
				grad_par[th*step_grad + fs*fc + 3*p + 1]  +=  d2Vector(dx_dPos_f.y, dy_dPos_f.y);
				grad_par[th*step_grad + fs*fc + 3*p + 2]  +=  d2Vector(dx_dPos_f.z, dy_dPos_f.z);
			}
		}
		
		if (f < fc-1)
		{
			//shift += getPosChange(x, p, f, fs * fc + 3 * pc);
			
			{
				d3Vector out(0.0, 0.0, 0.0);
				
				for (int b = 0; b < bc; b++)
				{
					const int i0 = fs*fc + 3*(pc + f*bc + b);
					const double def = deformationBasis[p*bc + b];
					
					for (int i = 0; i < 3; i++)
					{
						out[i] += x[i0+i] * def;
					}
				}
				
				shift += out;
			}
		}
	}
		
	if (!settings.constParticles)
	{
		for (int m = 0; m < fc-1; m++)
		{
			for (int b = 0; b < bc; b++)
			{
				d3Vector dx_dXm(0.0, 0.0, 0.0);
				d3Vector dy_dXm(0.0, 0.0, 0.0);	
				
				const double def = deformationBasis[p*bc + b];
				
				for (int f = m+1; f < fc; f++)
				{
					if (f == frame)
					{
						dx_dXm += def * dx_dPos[th*step_frame + f];
						dy_dXm += def * dy_dPos[th*step_frame + f];
					}
				}
									
				const int i0 = th*step_grad + fs*fc + 3*(pc + m*bc + b);
				
				grad_par[i0    ] += d2Vector(dx_dXm.x, dy_dXm.x);
				grad_par[i0 + 1] += d2Vector(dx_dXm.y, dy_dXm.y);
				grad_par[i0 + 2] += d2Vector(dx_dXm.z, dy_dXm.z);					
			}
		}
	}
	

	std::vector<d2Vector> dPos_dx(xs, d2Vector(0.0, 0.0));
	
	for (int i = 0; i < xs; i++)
	{
		dPos_dx[i] += grad_par[th*step_grad + i];
	}
	
	
	d2Vector pos0;
	
	shift = settings.constParticles? 
				d3Vector(0.0, 0.0, 0.0) : 
				d3Vector(x[fs*fc + 3*p], x[fs*fc + 3*p+1], x[fs*fc + 3*p+2]);
	
	for (int f = 0; f <= frame; f++)
	{
		if (f == frame)
		{
			const d4Vector pos4(initialPos[p] + shift);
			d4Vector dp4 = P[f] * pos4 - frameProj[f] * d4Vector(initialPos[p]);
			
			pos0 = d2Vector(dp4.x, dp4.y);
		}
		else if (!settings.constParticles)
		{
			shift += getPosChange(x, p, f, fs * fc + 3 * pc);
		}
	}
		
	for (int param_id = 0; param_id < xs; param_id++)
	{
		std::vector<double> xoff = x;
		xoff[param_id] += epsilon;
		
		d3Vector shift = settings.constParticles? 
					d3Vector(0.0, 0.0, 0.0) : 
					d3Vector(xoff[fs*fc + 3*p], xoff[fs*fc + 3*p+1], xoff[fs*fc + 3*p+2]);
		
		for (int f = 0; f <= frame; f++)
		{
			if (f == frame)
			{
				const d4Vector pos4(initialPos[p] + shift);				
				const d4Vector dp4  = P[f] * pos4 - frameProj[f] * d4Vector(initialPos[p]);
				const d2Vector pos1(dp4.x, dp4.y);
				
				std::cout << param_id << ": " << dPos_dx[param_id] << " vs. " << (pos1 - pos0)/epsilon << std::endl;
			}
			else if (!settings.constParticles)
			{
				shift += getPosChange(xoff, p, f, fs * fc + 3 * pc);
			}
		}
	}
}
