#ifndef ABERRATION_FIT_PROGRAM_H
#define ABERRATION_FIT_PROGRAM_H

#include "refinement.h"
#include <src/jaz/math/tensor2x2.h>


class AberrationFitProgram : public RefinementProgram
{
	public:
		
		enum Granularity 
		{
			PerTomogram, Global
		};
		
		class EvenData
		{
			public:
				
				double Axx, Axy, Ayy, bx, by;
				
				EvenData& operator+=(const EvenData& d);				
		};
		
		class OddData
		{
			public:
				
				double a;
				dComplex b;
				
				OddData& operator+=(const OddData& d);
		};
			 
		
		AberrationFitProgram(int argc, char *argv[]);
		
		
			Granularity granularity;
			int n_even, n_odd;
		
		void readParams(IOParser& parser);
		void run();
				
		
		static void considerParticle(
				int part_id, 
				const Tomogram& tomogram, 
				const TomoReferenceMap& referenceMap, 
				const ParticleSet* dataSet,
				bool flip_value,
				const BufferedImage<float>& frqWeight,
				int f0, int f1,
				BufferedImage<EvenData>& even_out, 
				BufferedImage<OddData>& odd_out);
		
		static std::vector<double> solveEven(
				const BufferedImage<EvenData>& data,
				int n_bands,
				double pixelSize, 
				std::string prefix,
				bool writeImages);
		
		static std::vector<double> solveOdd(
				const BufferedImage<OddData>& data,
				int n_bands,
				double pixelSize, 
				std::string prefix,
				bool writeImages);
};

#endif
