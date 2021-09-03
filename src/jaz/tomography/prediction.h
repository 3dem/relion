#ifndef PREDICTION_H
#define PREDICTION_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/optics/aberrations_cache.h>
#include <src/jaz/tomography/particle_set.h>

#include "reference_map.h"

class CTF;
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

		enum DoseWeight
		{
			NotDoseWeighted = 0,
			DoseWeighted = 1
		};

		enum CtfScale
		{
			CtfUnscaled = 0,
			CtfScaled = 1
		};
		

		static BufferedImage<fComplex> predictModulated(
				ParticleIndex particle_id,
				const ParticleSet& dataSet,
				gravis::d4Matrix proj,
				int s, const CTF& ctf, double pixelSize,
				const AberrationsCache& aberrationsCache,
				const std::vector<BufferedImage<fComplex>>& referenceFS,
				HalfSet halfSet = OwnHalf,
				Modulation modulation = AmplitudeAndPhaseModulated,
				const RawImage<float>* doseWeight = 0,
				CtfScale ctfScale = CtfScaled,
				const int* xRanges = 0);
		
		static BufferedImage<fComplex> predictFS(
				ParticleIndex particle_id,
				const ParticleSet& dataSet,
				gravis::d4Matrix proj,
				int s, 
				const std::vector<BufferedImage<fComplex>>& referenceFS,
				HalfSet halfSet = OwnHalf,
				const int* xRanges = 0);

		static std::vector<BufferedImage<double>> computeCroppedCCs(
				const ParticleSet& dataSet,
				const std::vector<ParticleIndex>& partIndices,
				const Tomogram& tomogram,
				const AberrationsCache& aberrationsCache,
				const TomoReferenceMap& referenceMap,
				const BufferedImage<float>& freqWeights,
				const BufferedImage<float>& doseWeights,
				const std::vector<int>& sequence,
				const BufferedImage<int>& xRanges,
				int maxRange,
				bool flip_value,
				int num_threads,
				double paddingFactor,
				HalfSet halfSet = OwnHalf,
				bool verbose = true);

		static void predictMicrograph(
				int frame_index,
				const ParticleSet& dataSet,
				const std::vector<ParticleIndex>& partIndices,
				const Tomogram& tomogram,
				const AberrationsCache& aberrationsCache,
				const TomoReferenceMap& referenceMap,
				const RawImage<float>* doseWeights,
				RawImage<float>& target_slice,
				HalfSet halfSet = OwnHalf,
				Modulation modulation = AmplitudeAndPhaseModulated,
				CtfScale ctfScale = CtfScaled);

				
				
};

#endif
