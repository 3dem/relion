#include "structure_tensor.h"
#include <src/jaz/segmentation/diffusion.h>
#include <src/jaz/segmentation/diffusion_tensors.h>

BufferedImage<Tensor3x3<float>> StructureTensor :: computeNonLinear(
	const RawImage<float>& src, 
	double rho, 
	double sigma0,
	double taper,
	int diffusion_iterations, 
	int outer_iterations,
	bool planar,
	double lambda,
	double stepSize,
	int num_threads)
{
	BufferedImage<Tensor3x3<float>> J = compute3D(src, rho, sigma0, taper);	
	J = StructureTensor::forwardAverage(J);
	
	for (int i = 0; i < outer_iterations; i++)
	{
		std::cout << i << " / " << outer_iterations << std::endl;
		
		BufferedImage<Tensor3x3<float>> D = DiffusionTensors::plateDiffusion(J, lambda);		
		std::vector<BufferedImage<float>> JmVec = StructureTensor::split(J);
		
		for (int i = 0; i < 6; i++)
		{
			JmVec[i] = Diffusion::diffuse(
				JmVec[i], D, stepSize, diffusion_iterations, num_threads);
		}
		
		J = StructureTensor::join3D(JmVec);
	}
	
	return J;
}
