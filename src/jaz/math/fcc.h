#ifndef RELION_TOMO_FCC_H
#define RELION_TOMO_FCC_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/particle_set.h>


class FCC
{
	public:
		
		static BufferedImage<double> compute(
				const ParticleSet& dataSet,
				const std::vector<ParticleIndex>& partIndices,
				const Tomogram& tomogram,
				const std::vector<BufferedImage<fComplex>>& referenceFS,
				bool flip_value,
				int num_threads);
		
		static BufferedImage<double> compute3(
				const ParticleSet& dataSet,
				const std::vector<ParticleIndex>& partIndices,
				const Tomogram& tomogram,
				const std::vector<BufferedImage<fComplex>>& referenceFS,
				bool flip_value,
				int num_threads);
		
		static BufferedImage<double> divide(
				const BufferedImage<double>& fcc3);

		static BufferedImage<double> sumOverTime(
				const RawImage<double>& fcc3);
};

#endif
