#ifndef PROTO_ALIGNMENT_H
#define PROTO_ALIGNMENT_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/image/buffered_image.h>

class CTF;

class ProtoAlignment : public DifferentiableOptimization
{
	public:
		
		ProtoAlignment(
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
				gravis::d3Vector tomoCentre,
				int progressBarOffset,
				int num_threads,
				double paddingFactor);
		
		
			std::vector<gravis::d4Matrix> frameProj; // make a reference again
			const ParticleSet& dataSet;
			const std::vector<ParticleIndex>& partIndices;
			const std::vector<BufferedImage<fComplex>>& referenceFS;
			
			bool	constParticles, constAngles, constShifts,
					doAnisotropy, perTiltAnisotropy,
					devMode;

			int fc, pc, maxRange;
			gravis::d3Vector tomoCentre;
			int progressBarOffset, num_threads;
			double paddingFactor;

			mutable int lastIterationNumber;
			
			const std::vector<BufferedImage<double>>& CCs; // one frame stack for each particle
			std::vector<gravis::d3Vector> initialPos;
			
			gravis::d4Matrix minusCentre, plusCentre;
			
			
			
		double f(const std::vector<double>& x, void* tempStorage) const;
		void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;
				
		std::vector<gravis::d4Matrix> getProjections(
				const std::vector<double>& x) const;

		void shiftParticles(
				const std::vector<double>& x,
				const std::vector<ParticleIndex>& partIndices,
				ParticleSet& target) const;

		std::vector<gravis::d3Vector> getParticlePositions(
				const std::vector<double>& x) const;
		
		int getParamCount();
		
		void protocol(const std::vector<double>& x, int frame0, int frame1) const;
	
		std::vector<BufferedImage<double>> drawShiftedCCs(const std::vector<double>& x) const;
		
		void report(int iteration, double cost, const std::vector<double>& x) const;

		std::vector<double> createInitial();

		void visualiseShifts(
				const std::vector<double> &x,
				const std::string &tomo_name,
				const std::string &file_name_root) const;
		
	protected:
	
		inline int getFrameStride() const
		{
			int out = 0;
			
			if (!constAngles) out += 3;
			if (!constShifts) out += 2;
			if (doAnisotropy && perTiltAnisotropy) out += 2;
			
			return out;
		}

		inline int getTiltDataOffset(int f) const
		{
			if (!doAnisotropy || perTiltAnisotropy)
			{
				return f * getFrameStride();
			}
			else
			{
				return f * getFrameStride() + 2;
			}
		}

		inline int getParticleDataOffset() const
		{
			return getTiltDataOffset(fc);
		}
		
		inline void readTiltParameters(
				const std::vector<double>& x, int f,
				double& phi, double& theta, double& psi,
				double& dx, double& dy,
				double& skew, double& y_scale) const
		{
			int offset = getTiltDataOffset(f);

			if (constAngles)
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

			if (constShifts)
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

			if (doAnisotropy)
			{
				if (perTiltAnisotropy)
				{
					skew    = x[offset  ];
					y_scale = x[offset+1];
				}
				else
				{
					skew    = x[0];
					y_scale = x[1];
				}
			}
			else
			{
				skew    = 0.0;
				y_scale = 1.0;
			}
		}
			   
			   
};

#endif
