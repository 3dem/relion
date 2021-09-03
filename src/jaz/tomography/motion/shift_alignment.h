#ifndef SHIFT_ALIGNMENT_H
#define SHIFT_ALIGNMENT_H

#include <vector>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/reference_map.h>
#include <src/jaz/optics/aberrations_cache.h>


class ShiftAlignment
{
	public:

		static std::vector<gravis::d2Vector> alignGlobally(
				const Tomogram& tomogram,
				const std::vector<ParticleIndex>& particles,
				const ParticleSet& particleSet,
				const TomoReferenceMap& referenceMap,
				const RawImage<float>& doseWeights,
				const AberrationsCache& aberrationsCache,
				bool do_whiten,
				int num_threads,
				bool diag,
				const std::string& tag,
				const std::string& outDir);

		static std::vector<gravis::d2Vector> alignPerParticle(
				const Tomogram& tomogram,
				const std::vector<BufferedImage<double>>& CCs,
				double padding,
				int range,
				int verbosity,
				bool diag,
				bool per_tomogram_progress,
				const std::string& tag,
				const std::string& outDir);

		static void visualiseShifts(
				const std::vector<gravis::d2Vector>& shifts,
				const std::vector<int>& sequence,
				const std::string& tomo_name,
				const std::string& file_name_root);

};

#endif
