#include "spline_2D_deformation_model.h"

Spline2DDeformationModel::Spline2DDeformationModel(
		const gravis::i2Vector& imageSize, 
		const gravis::i2Vector& gridSize)
	:
	  imageSize(imageSize),
	  gridSize(gridSize),
	  gridSpacing(
		imageSize.x / (double)(gridSize.x - 1), 
		imageSize.y / (double)(gridSize.y - 1))
{
}
