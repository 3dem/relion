#ifndef GP_MOTION_MODEL_H
#define GP_MOTION_MODEL_H

#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/tomography/particle_set.h>

class Tomogram;

class GPMotionModel
{
	public:
		
			struct Parameters
			{
				double sig_vel, sig_div;
				bool params_scaled_by_dose, sqExpKernel;
				int maxEDs;
			};

			
		GPMotionModel(
				const ParticleSet& dataSet,
				const std::vector<ParticleIndex>& partIndices,
				const Tomogram& tomogram,
				Parameters parameters,
				bool verbose);

		
			std::vector<gravis::d3Vector> initialPos;
			std::vector<double> deformationBasis, deformationLambda;
			int pc, bc;


		inline void updatePosition(
				const double* x,
				int particle_index,
				gravis::d3Vector& position) const;

		inline void updateDataTermGradient(
				const gravis::d3Vector *dC_dPos,
				int particle_index,
				int fc,
				double *target) const;

		inline double computePriorValueAndGradient(
				const double *x,
				int fc,
				double* gradDest) const;

		inline int getParameterCount() const;


};

inline void GPMotionModel::updatePosition(
    const double* x,
	int particle_index,
	gravis::d3Vector& position) const
{
	gravis::d3Vector d(0.0, 0.0, 0.0);

	for (int b = 0; b < bc; b++)
	{
		const double def = deformationBasis[particle_index*bc + b];
		
		d.x += x[3*b    ] * def;
		d.y += x[3*b + 1] * def;
		d.z += x[3*b + 2] * def;
	}

	position += d;
}

inline void GPMotionModel::updateDataTermGradient(
	const gravis::d3Vector* dC_dPos,
	int particle_index,
	int fc,
	double* target) const
{
	for (int m = 0; m < fc - 1; m++)
	{
		for (int b = 0; b < bc; b++)
		{
			gravis::d3Vector dC_dXm(0.0, 0.0, 0.0);
			const double def = deformationBasis[particle_index*bc + b];

			for (int f = m+1; f < fc; f++)
			{
				dC_dXm += def * dC_dPos[f];
			}

			const int i0 = 3*(m*bc + b);

			target[i0    ] += dC_dXm.x;
			target[i0 + 1] += dC_dXm.y;
			target[i0 + 2] += dC_dXm.z;
		}
	}
}

inline double GPMotionModel::computePriorValueAndGradient(
	const double* x,
	int fc,
	double* gradDest) const
{
	double cost = 0.0;
	
	for (int m = 0; m < fc-1; m++)
	{
		for (int b = 0; b < bc; b++)
		{
			const int i0 = 3*(m*bc + b);
			
			gradDest[i0  ] += 2.0 * x[i0  ];
			gradDest[i0+1] += 2.0 * x[i0+1];
			gradDest[i0+2] += 2.0 * x[i0+2];
			
			cost += x[i0]*x[i0] + x[i0+1]*x[i0+1] + x[i0+2]*x[i0+2];
		}
	}
	
	return cost;
}

inline int GPMotionModel::getParameterCount() const
{
	return 3 * bc;
}

#endif
