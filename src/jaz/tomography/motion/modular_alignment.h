#ifndef MODULAR_ALIGNMENT_H
#define MODULAR_ALIGNMENT_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/jaz/util/log.h>
#include "trajectory.h"

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

template<class MotionModel, class DeformationModel>
class ModularAlignment : public FastDifferentiableOptimization
{
	public:



		ModularAlignment(
				const std::vector<BufferedImage<double>>& CCs,
				const std::vector<gravis::d4Matrix>& frameProj,
				ParticleSet& dataSet,
				const std::vector<ParticleIndex>& partIndices,
				const std::vector<BufferedImage<fComplex>>& referenceFS,
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
			DeformationModel deformationModel;

			const std::vector<BufferedImage<double>>& CCs; // one frame stack for each particle
			const std::vector<gravis::d4Matrix>& frameProj;
			ParticleSet& dataSet;
			const std::vector<ParticleIndex>& partIndices;
			const std::vector<BufferedImage<fComplex>>& referenceFS;

			ModularAlignmentMotionParameters motionParameters;
			ModularAlignmentSettings settings;

			gravis::d3Vector tomoCentre;
			double frameDose, pixelSize, paddingFactor;
			int progressBarOffset, num_threads;

			bool verbose;

			int fc, pc, bc, maxRange;


			std::vector<gravis::d3Vector> initialPos;

			gravis::d4Matrix minusCentre, plusCentre;



		double f(const std::vector<double>& x, void* tempStorage) const;
		void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;
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
				const std::vector<double>& x, int offset,
				double& phi, double& theta, double& psi,
				double& dx, double& dy) const
		{
			if (settings.constAngles)
			{
				if (settings.constShifts)
				{
					phi   = 0.0;
					theta = 0.0;
					psi   = 0.0;

					dx    = 0.0;
					dy    = 0.0;
				}
				else
				{
					phi   = 0.0;
					theta = 0.0;
					psi   = 0.0;

					dx    = x[offset  ];
					dy    = x[offset+1];
				}
			}
			else
			{
				if (settings.constShifts)
				{
					phi   = x[offset  ];
					theta = x[offset+1];
					psi   = x[offset+2];

					dx    = 0.0;
					dy    = 0.0;
				}
				else
				{
					phi   = x[offset  ];
					theta = x[offset+1];
					psi   = x[offset+2];
					dx    = x[offset+3];
					dy    = x[offset+4];
				}
			}
		}
};

class StaticMotionModel
{
	public:

		StaticMotionModel();

		inline gravis::d3Vector getPosChange(
			const std::vector<double>& x,
			int particle,
			int mode,
			int offset) const
		{
			return gravis::d3Vector(0.0, 0.0, 0.0);
		}
};

class LinearMotionModel
{
	public:

		LinearMotionModel(
			ModularAlignmentMotionParameters motionParameters,
			ModularAlignmentSettings settings,
			std::vector<gravis::d3Vector> initialPos,
			double frameDose,
			double pixelSize,
			bool verbose);


			std::vector<double> deformationBasis, deformationLambda;
			int bc;


		inline gravis::d3Vector getPosChange(
			const std::vector<double>& x,
			int particle,
			int mode,
			int offset) const
		{
			gravis::d3Vector out(0.0, 0.0, 0.0);

			for (int b = 0; b < bc; b++)
			{
				const int i0 = offset + 3*(mode*bc + b);
				const double def = deformationBasis[particle*bc + b];

				for (int i = 0; i < 3; i++)
				{
					out[i] += x[i0+i] * def;
				}
			}

			return out;
		}

};

#endif
