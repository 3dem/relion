#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/tomography/particle_set.h>

class Tomogram;

class GPMotionModel
{
	public:
		
			struct Settings
			{
				bool constParticles, constAngles, constShifts,
					 params_scaled_by_dose, sqExpKernel;
				
				int maxEDs;
			};
			
			struct MotionParameters
			{
				double sig_vel, sig_div;
			};

			
		GPMotionModel(const ParticleSet& dataSet,
			const std::vector<ParticleIndex>& partIndices,
		    const Tomogram& tomogram,
		    MotionParameters motionParameters,    
			Settings settings, 
		    bool verbose);

		
			std::vector<gravis::d3Vector> initialPos;
			std::vector<double> deformationBasis, deformationLambda;
			int pc, bc;


		inline gravis::d3Vector getPosChange(
			const std::vector<double>& x,
			int particle,
			int mode,
			int offset) const;
		
		inline void updateCostGradient(
			const std::vector<gravis::d3Vector>& dC_dPos,
			int offset_in,
			int particle_index,
			std::vector<double>& target,
			int offset_out) const;
		
		inline double computePriorCostAndGradient(
		    const std::vector<double>& x,
			int offset,
			int fc,
			std::vector<double>& gradDest) const;
		        
		        

};


inline gravis::d3Vector GPMotionModel::getPosChange(
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

inline void GPMotionModel::updateCostGradient(
	const std::vector<gravis::d3Vector>& dC_dPos,
	int offset_in,
	int particle_index,
	std::vector<double>& target,
	int offset_out) const
{
	const int fc = dC_dPos.size();
	
	for (int m = 0; m < fc-1; m++)
	{
		for (int b = 0; b < bc; b++)
		{
			gravis::d3Vector dC_dXm(0.0, 0.0, 0.0);
			const double def = deformationBasis[particle_index*bc + b];

			for (int f = m+1; f < fc; f++)
			{
				dC_dXm += def * dC_dPos[offset_in + f];
			}

			const int i0 = offset_out + 3*(m*bc + b);

			target[i0    ] += dC_dXm.x;
			target[i0 + 1] += dC_dXm.y;
			target[i0 + 2] += dC_dXm.z;
		}
	}
}

inline double GPMotionModel::computePriorCostAndGradient(
    const std::vector<double>& x,
	int offset, 
	int fc,
	std::vector<double>& gradDest) const
{
	double cost = 0.0;
	
	for (int m = 0; m < fc-1; m++)
	{
		for (int b = 0; b < bc; b++)
		{
			const int i0 = offset + 3*(m*bc + b);
			
			gradDest[i0  ] += 2.0 * x[i0  ];
			gradDest[i0+1] += 2.0 * x[i0+1];
			gradDest[i0+2] += 2.0 * x[i0+2];
			
			cost += x[i0]*x[i0] + x[i0+1]*x[i0+1] + x[i0+2]*x[i0+2];
		}
	}
	
	return cost;
}
