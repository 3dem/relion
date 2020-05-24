#ifndef PREDICTION_H
#define PREDICTION_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t4Matrix.h>

#include "reference_map.h"

class ParticleSet;
class CTF;
class TomoList;
class Tomogram;


class Prediction
{
	public:
		
		enum HalfSet
		{
			OwnHalf = 0,
			OppositeHalf = 1
		};
		
		enum Modulation
		{
			Unmodulated = 0,
			AmplitudeModulated = 1,
			PhaseModulated = 2,
			AmplitudeAndPhaseModulated = 3
		};
		
		static BufferedImage<fComplex> predictFS(
				int particle_id, 
				const ParticleSet* dataSet,
				gravis::d4Matrix proj,
				int s, const CTF& centreCtf,
				gravis::d3Vector tomo_centre,
				double handedness, 
				double pixelSize,
				const std::vector<BufferedImage<fComplex>>& referenceFS,
				HalfSet halfSet = OwnHalf,
				Modulation modulation = AmplitudeAndPhaseModulated);
		
		static BufferedImage<fComplex> predictFS(
				int particle_id, 
				const ParticleSet* dataSet,
				gravis::d4Matrix proj,
				int s, 
				const std::vector<BufferedImage<fComplex>>& referenceFS,
				HalfSet halfSet = OwnHalf);
		
		static BufferedImage<float> predictRS(
				int particle_id, 
				const ParticleSet* dataSet,
				gravis::d4Matrix proj,
				int s, 
				const std::vector<BufferedImage<fComplex>>& referenceFS,
				HalfSet halfSet = OwnHalf);
		
		static std::vector<BufferedImage<double>> computeCroppedCCs(
				const ParticleSet* dataSet,
				const std::vector<int>& partIndices,
				const Tomogram& tomogram,
				const TomoReferenceMap& referenceMap,
				const BufferedImage<float>& frqWghts,
				const std::vector<int>& sequence,
				int maxRange,				
				bool flip_value,
				int num_threads,
				double paddingFactor,
				HalfSet halfSet = OwnHalf);
				
				
};

#endif
