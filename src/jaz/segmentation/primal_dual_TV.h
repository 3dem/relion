#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/math/tensor3x3.h>

class PrimalDualTV
{
	public:
		
		static BufferedImage<float> anisotropicTV(
				RawImage<float> regionalCost,
				RawImage<Tensor3x3<float>> diffusionTensor,
				int maxIterations = 100,
				float sigma = 0.1f, 
				float tau = 0.1f,
				float nu = 1.8f);
};

#endif
