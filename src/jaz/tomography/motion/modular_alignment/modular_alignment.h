#ifndef MODULAR_ALIGNMENT_H
#define MODULAR_ALIGNMENT_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/jaz/util/log.h>
#include <src/jaz/tomography/motion/trajectory.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/ctf.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/math/Tait_Bryan_angles.h>
#include <src/jaz/math/Gaussian_process.h>
#include <src/jaz/util/zio.h>
#include <omp.h>

class CTF;

struct ModularAlignmentSettings
{
	bool constParticles, constAngles, constShifts,
		 params_scaled_by_dose, sqExpKernel;

	int maxEDs;
};

struct ModularAlignmentMotionParameters
{
	double sig_vel, sig_div;
};

template<class MotionModel, class DeformationModel2D>
class ModularAlignment : public FastDifferentiableOptimization
{
	public:

		ModularAlignment(
				const std::vector<BufferedImage<double>>& CCs,
				const std::vector<gravis::d4Matrix>& frameProj,
				ParticleSet& dataSet,
				const std::vector<ParticleIndex>& partIndices,
				ModularAlignmentMotionParameters motionParameters,
				ModularAlignmentSettings settings,
				gravis::d3Vector tomoCentre,
				double frameDose,
				double pixelSize,
				double paddingFactor,
				int progressBarOffset,
				int num_threads,
				bool verbose);
		

			MotionModel motionModel;
			DeformationModel2D deformationModel2D;

			const std::vector<BufferedImage<double>>& CCs;  // one frame stack for each particle
			const std::vector<gravis::d4Matrix>& frameProj; // initial projection matrices
			ParticleSet& dataSet;
			const std::vector<ParticleIndex>& partIndices;

			ModularAlignmentMotionParameters motionParameters;
			ModularAlignmentSettings settings;

			gravis::d3Vector tomoCentre;
			double frameDose, pixelSize, paddingFactor;
			int progressBarOffset, num_threads;

			bool verbose;

			int fc, pc, bc, maxRange;


			std::vector<gravis::d3Vector> initialPos;

			gravis::d4Matrix minusCentre, plusCentre;
			
			

		double gradAndValue(const std::vector<double>& x, std::vector<double>& gradDest) const;
			
		std::vector<gravis::d4Matrix> getProjections(
				const std::vector<double>& x,
				const std::vector<int>& frameSequence) const;

		std::vector<gravis::d3Vector> getParticlePositions(
				const std::vector<double>& x) const;

		Trajectory getTrajectory(
				const std::vector<double>& x,
				int p,
				const std::vector<int>& frameSequence) const;

		std::vector<Trajectory> exportTrajectories(
				const std::vector<double>& x,
				const ParticleSet& dataSet,
				const std::vector<int>& frameSequence) const;

		int getParamCount();

		std::vector<BufferedImage<double>> drawShiftedCCs(const std::vector<double>& x) const;
		Mesh visualiseTrajectories(const std::vector<double>& x, double scale);

		void report(int iteration, double cost, const std::vector<double>& x) const;


	protected:

		inline int getFrameStride() const
		{
			int out = 0;

			if (!settings.constAngles) out += 3;
			if (!settings.constShifts) out += 2;

			return out;
		}

		inline void readViewParams(
			const std::vector<double>& x, int f,
			double& phi, double& theta, double& psi,
			double& dx, double& dy) const;
};

template<class MotionModel, class DeformationModel2D>
ModularAlignment<MotionModel, DeformationModel2D>::ModularAlignment(
		const std::vector<BufferedImage<double>>& CCs,
		const std::vector<gravis::d4Matrix>& frameProj,
		ParticleSet& dataSet,
		const std::vector<ParticleIndex>& partIndices,
		ModularAlignmentMotionParameters motionParameters,
		ModularAlignmentSettings settings,
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
	
	minusCentre = gravis::d4Matrix(
			1, 0, 0, -tomoCentre.x, 
			0, 1, 0, -tomoCentre.y, 
			0, 0, 1, -tomoCentre.z, 
			0, 0, 0, 1 );
	
	plusCentre = gravis::d4Matrix(
			1, 0, 0, tomoCentre.x, 
			0, 1, 0, tomoCentre.y, 
			0, 0, 1, tomoCentre.z, 
			0, 0, 0, 1 );
}

template<class MotionModel, class DeformationModel2D>
double ModularAlignment<MotionModel, DeformationModel2D>::gradAndValue(
        const std::vector<double>& x, std::vector<double>& gradDest) const
{
	const int fs = getFrameStride();
	const int xs = x.size();
	const int data_pad = 512;
	const int step_grad = xs + data_pad;
	const int step_frame = fc + data_pad;

	std::vector<gravis::d4Matrix> P(fc), P_phi(fc), P_theta(fc), P_psi(fc);

	for (int f = 0; f < fc; f++)
	{
		double phi, theta, psi, dx, dy;
		readViewParams(x, f, phi, theta, psi, dx, dy);

		const gravis::d4Matrix Q = 
		        TaitBryan::anglesToMatrix4(phi, theta, psi);
		
		gravis::t4Vector<gravis::d3Matrix> dQ = 
		        TaitBryan::anglesToMatrixAndDerivatives(phi, theta, psi);

		gravis::d4Matrix Q_phi(dQ[0]);
		gravis::d4Matrix Q_theta(dQ[1]);
		gravis::d4Matrix Q_psi(dQ[2]);

		Q_phi(3,3) = 0.0;
		Q_theta(3,3) = 0.0;
		Q_psi(3,3) = 0.0;

		const gravis::d4Matrix centProj = minusCentre * frameProj[f];

		P[f] = plusCentre * Q * centProj;

		P_phi[f]   = plusCentre * Q_phi   * centProj;
		P_theta[f] = plusCentre * Q_theta * centProj;
		P_psi[f]   = plusCentre * Q_psi   * centProj;

		P[f](0,3) += dx;
		P[f](1,3) += dy;
	}


	std::vector<double> grad_par(step_grad * num_threads, 0.0);
	std::vector<double> val_par(data_pad * num_threads, 0.0);
	std::vector<gravis::d3Vector> dC_dPos(step_frame * num_threads, gravis::d3Vector(0.0, 0.0, 0.0));


	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{

		const int th = omp_get_thread_num();

		gravis::d3Vector shift = settings.constParticles?
			gravis::d3Vector(0.0, 0.0, 0.0) :
			gravis::d3Vector(x[fs*fc + 3*p], x[fs*fc + 3*p+1], x[fs*fc + 3*p+2]);

		for (int f = 0; f < fc; f++)
		{
			const gravis::d4Vector pos4(initialPos[p] + shift);

			const gravis::d4Vector dp  = P[f] * pos4 
			        - frameProj[f] * gravis::d4Vector(initialPos[p]);

			const double dx_img = (dp.x + maxRange) * paddingFactor;
			const double dy_img = (dp.y + maxRange) * paddingFactor;

			const gravis::d4Vector dp_phi   = P_phi[f]   * pos4;
			const gravis::d4Vector dp_theta = P_theta[f] * pos4;
			const gravis::d4Vector dp_psi   = P_psi[f]   * pos4;


			const gravis::d3Vector g = -((double)paddingFactor) *
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
				const gravis::d3Vector dC_dPos_f(
					P[f](0,0) * g.x  +  P[f](1,0) * g.y,
					P[f](0,1) * g.x  +  P[f](1,1) * g.y,
					P[f](0,2) * g.x  +  P[f](1,2) * g.y);

				dC_dPos[th*step_frame + f] = dC_dPos_f;

				grad_par[th*step_grad + fs*fc + 3*p    ]  +=  dC_dPos_f.x;
				grad_par[th*step_grad + fs*fc + 3*p + 1]  +=  dC_dPos_f.y;
				grad_par[th*step_grad + fs*fc + 3*p + 2]  +=  dC_dPos_f.z;

				if (f < fc-1)
				{
					shift += motionModel.getPosChange(x, p, f, fs * fc + 3 * pc);
				}
			}
		}

		if (!settings.constParticles)
		{
			motionModel.updateCostGradient(
				dC_dPos, th*step_frame, p, 
				grad_par, th*step_grad + fs*fc + 3*pc);
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
		const int offset = fs*fc + 3*pc;
		
		cost += deformationModel2D.computePriorCostAndGradient(
			x, offset, fc, gradDest);
	}

	return cost;
}

template<class MotionModel, class DeformationModel2D>
inline void ModularAlignment<MotionModel, DeformationModel2D>::readViewParams(
		const std::vector<double>& x, int f,
		double& phi, double& theta, double& psi,
		double& dx, double& dy) const
{
	const int fs = getFrameStride();
	int offset = f * fs;
	
	if (settings.constAngles)
	{
		phi   = 0.0;
		theta = 0.0;
		psi   = 0.0;
	}
	else
	{
		phi   = x[offset  ];
		theta = x[offset+1];
		psi   = x[offset+2];

		offset += 3;
	}

	if (settings.constShifts)
	{
		dx    = 0.0;
		dy    = 0.0;
	}
	else
	{
		dx    = x[offset  ];
		dy    = x[offset+1];

		offset += 2;
	}
}





#endif
