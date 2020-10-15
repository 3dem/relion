#ifndef PROTO_ALIGNMENT_H
#define PROTO_ALIGNMENT_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/image/buffered_image.h>

class ParticleSet;
class CTF;

class ProtoAlignment : public DifferentiableOptimization
{
	public:
		
		ProtoAlignment(
				const std::vector<BufferedImage<double>>& CCs,
				const std::vector<gravis::d4Matrix>& frameProj, 
				const ParticleSet& dataSet,
				const std::vector<int>& partIndices,
				const std::vector<BufferedImage<fComplex>>& referenceFS,
				bool constParticles,
				bool constAngles,
				bool constShifts,
				int maxRange, 
				gravis::d3Vector tomoCentre,
				int num_threads,
				double paddingFactor);
		
		
			std::vector<gravis::d4Matrix> frameProj; // make a reference again
			const ParticleSet& dataSet;
			const std::vector<int>& partIndices;
			const std::vector<BufferedImage<fComplex>>& referenceFS;
			
			bool constParticles, constAngles, constShifts;
			int fc, pc, maxRange;	
			gravis::d3Vector tomoCentre;
			int num_threads;
			double paddingFactor;
			
			const std::vector<BufferedImage<double>>& CCs; // one frame stack for each particle
			std::vector<gravis::d3Vector> initialPos;
			
			gravis::d4Matrix minusCentre, plusCentre;
			
			
			
		double f(const std::vector<double>& x, void* tempStorage) const;	
		void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;
				
		std::vector<gravis::d4Matrix> getProjections(const std::vector<double>& x) const;	
		
		void shiftParticles(
				const std::vector<double>& x,
				const std::vector<int>& partIndices,
				ParticleSet& target) const;
		
		int getParamCount();
		
		void protocol(const std::vector<double>& x, int frame0, int frame1) const;
	
		std::vector<BufferedImage<double>> drawShiftedCCs(const std::vector<double>& x) const;
		
		void report(int iteration, double cost, const std::vector<double>& x) const;
	
		
	protected:
	
		inline int getFrameStride() const
		{
			int out = 0;
			
			if (!constAngles) out += 3;
			if (!constShifts) out += 2;
			
			return out;
		}
		
		inline void readParams(
				const std::vector<double>& x, int offset,
				double& phi, double& theta, double& psi,
				double& dx, double& dy) const
		{
			if (constAngles)
			{
				if (constShifts)
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
				if (constShifts)
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

#endif
