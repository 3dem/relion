#ifndef LOCAL_PARTICLE_REFINEMENT_H
#define LOCAL_PARTICLE_REFINEMENT_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/reference_map.h>
#include <src/jaz/optics/aberrations_cache.h>


class LocalParticleRefinement : public FastDifferentiableOptimization
{
	public:

		LocalParticleRefinement(
				ParticleIndex particle_id,
				const ParticleSet& particleSet,
				const Tomogram& tomogram,
				const TomoReferenceMap& reference,
				const BufferedImage<float>& freqWeight,
				const BufferedImage<float>& doseWeight,
				const AberrationsCache& aberrationsCache,
				double dose_cutoff,
				int minFrame,
				int maxFrame);


			const ParticleIndex particle_id;
			const ParticleSet& particleSet;
			const Tomogram& tomogram;
			const TomoReferenceMap& reference;
			const BufferedImage<float>& freqWeight;
			const BufferedImage<float>& doseWeight;
			const AberrationsCache& aberrationsCache;

			int minFrame, maxFrame;

			BufferedImage<fComplex> observations;
			std::vector<gravis::d3Matrix> Pt;
			std::vector<CTF> CTFs;
			gravis::d3Vector position;
			std::vector<int> max_radius;
			BufferedImage<float> precomputedCTFs;
			std::vector<bool> isVisible;


		double f(const std::vector<double>& x, void* tempStorage) const;
		void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;

		double gradAndValue(const std::vector<double>& x, std::vector<double>& gradDest) const;

		static void applyChange(
				const std::vector<double>& x,
				ParticleSet& target,
				ParticleIndex particle_id,
				double pixel_size);
};

#endif
