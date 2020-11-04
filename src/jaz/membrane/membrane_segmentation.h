#ifndef MEMBRANE_SEGMENTATION_H
#define MEMBRANE_SEGMENTATION_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t3Vector.h>

class MembraneSegmentation
{
	public:

		static BufferedImage<float> constructMembraneKernel(
				int w, int h, int d,
				double falloff,
				double kernel_width,
				double spacing,
				double ratio,
				double depth,
				double angle = 0);

		static BufferedImage<float> correlateWithMembrane(
				BufferedImage<float>& map,
				double falloff,
				double kernel_width,
				double spacing,
				double ratio,
				double depth,
				double angle = 0);

		static BufferedImage<float> correlateWithMembraneMultiAngle(
				BufferedImage<float>& map,
				double falloff,
				double kernel_width,
				double spacing,
				double ratio,
				double depth,
				double max_tilt,
				int tilt_steps);
		
		static BufferedImage<float> determineMembraniness(
				const RawImage<float>& tomo,
				RawImage<Tensor3x3<float>>& J,
				double sigma_diff,
				double lambda_fin,
				double thresh_edge,
				int num_threads,
				int itersAcross = 30,
				int itersAlong = 100,
				float stepSize = 0.05f);
		
		static BufferedImage<float> makeRegionalTerm();
		
		static BufferedImage<float> findSilhouette3D(
				const RawImage<float>& membrane,
				RawImage<Tensor3x3<float>>& J,
				gravis::d4Matrix view,
				double kappa);
		
		static BufferedImage<float> softmaxMembraneDist(
				const RawImage<float>& membrane,
				float lambda);


};

#endif
