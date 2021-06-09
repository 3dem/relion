#include "Fourier_2D_deformation_model.h"

Fourier2DDeformationModel::Fourier2DDeformationModel()
{
	
}

Fourier2DDeformationModel::Fourier2DDeformationModel(
		const Deformation2D::Parameters& parameters,
		const gravis::i2Vector& imageSize)
	:
	  imageSize(imageSize),
	  gridSize(parameters.grid_width, parameters.grid_height),
	  gridSpacing(
		imageSize.x / (double)(parameters.grid_width - 1),
		imageSize.y / (double)(parameters.grid_height - 1))
{
}
