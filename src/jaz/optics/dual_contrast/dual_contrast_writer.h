#ifndef DUAL_CONTRAST_WRITER_H
#define DUAL_CONTRAST_WRITER_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/cutting.h>
#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/raw_image.h>
#include <src/jaz/gravis/tRGB.h>
#include <src/jaz/gravis/tImage.h>


class DualContrastWriter
{
	public:

		template <typename T>
		static void writeAxialSlices(
				const RawImage<T>& phase_map,
				const RawImage<T>& amplitude_map,
				std::string filename_prefix,
				double output_scale = 0.2);

};

template <typename T>
void DualContrastWriter :: writeAxialSlices(
		const RawImage<T>& phase_map,
		const RawImage<T>& amplitude_map,
		std::string filename_prefix,
		double output_scale)
{
	std::vector<std::string> slice_names {"yz", "zx", "xy"};

	const double mean_phase = Normalization::computeMean(phase_map);
	const double mean_amplitude = Normalization::computeMean(amplitude_map);

	const double phase_variance = Normalization::computeVariance(phase_map, mean_phase);
	const double amplitude_variance = Normalization::computeVariance(amplitude_map, mean_amplitude);

	const double phase_scale = 1.0 / sqrt(phase_variance);
	const double amplitude_scale = 1.0 / sqrt(amplitude_variance);

	const int s = phase_map.xdim;

	for (int dim = 0; dim < 3; dim++)
	{
		BufferedImage<double> slice_amp = Cutting::extractAxialSlice(
					amplitude_map, dim, s/2);

		BufferedImage<double> slice_phase = Cutting::extractAxialSlice(
					phase_map, dim, s/2);

		gravis::tImage<gravis::dRGB> pngOut(s,s);
		pngOut.fill(gravis::dRGB(0.f));

		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			double p = phase_scale * (slice_phase(x,y) - mean_phase);
			double a = amplitude_scale * (slice_amp(x,y) - mean_amplitude);

			pngOut(x,y) = output_scale * gravis::dRGB(p-a,p,p+a) + 0.5 * gravis::dRGB(1,1,1);
			pngOut(x,y).clamp();
		}

		pngOut.writePNG(filename_prefix + slice_names[dim] + ".png");
	}
}


#endif
