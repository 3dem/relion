#include "Fourier_2D_deformation_model.h"

Fourier2DDeformationModel::Fourier2DDeformationModel()
{
	
}

Fourier2DDeformationModel::Fourier2DDeformationModel(
		const Deformation2D::Parameters& parameters,
		const gravis::i2Vector& imageSize)
	:
	  imageSize(imageSize),
	  gridSize(parameters.grid_width, parameters.grid_height)
{
	spatialFrequencies.reserve((gridSize.x/2 + 1) * gridSize.y);

	for (int y = 0; y < gridSize.y; y++)
	{
		if (y < gridSize.y/2)
		{
			for (int x = (y > 0? 0 : 1); x < (gridSize.x/2 + 1); x++)
			{
				spatialFrequencies.push_back(
					gravis::d2Vector(
						x * PI / imageSize.x,
						y * PI / imageSize.y));
			}
		}
		else
		{
			for (int x = 1; x < (gridSize.x/2 + 1); x++)
			{
				spatialFrequencies.push_back(
					gravis::d2Vector(
						x * PI / imageSize.x,
						(y - gridSize.y) * PI / imageSize.y));
			}
		}
	}
}
