#ifndef DIFFUSION_TENSORS_H
#define DIFFUSION_TENSORS_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/math/tensor3x3.h>

class DiffusionTensors
{
	public:
		
		static BufferedImage<Tensor3x3<float>> plateDiffusion(
				const RawImage<Tensor3x3<float>>& J,
				double sigma, double alpha = 0.0);
		
		static BufferedImage<Tensor3x3<float>> weightedPlateDiffusion(
				const RawImage<Tensor3x3<float>>& J,
				const RawImage<float>& edges,
				double sigma, double gamma, double alpha = 0.0);
		
		static BufferedImage<Tensor3x3<float>> membraneDiffusion(
				const RawImage<Tensor3x3<float>>& J,
				const RawImage<float>& edges);
		
		static BufferedImage<Tensor3x3<float>> antiPlateDiffusion(
				const RawImage<Tensor3x3<float>>& J,
				double alpha = 0.0);
};

#endif
