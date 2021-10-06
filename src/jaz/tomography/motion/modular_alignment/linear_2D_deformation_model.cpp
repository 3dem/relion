#include "linear_2D_deformation_model.h"

Linear2DDeformationModel::Linear2DDeformationModel()
{
}

Linear2DDeformationModel::Linear2DDeformationModel(
		const Deformation2D::Parameters& parameters, 
		const gravis::i2Vector& imageSize)
:	imageSize(imageSize),
	imageCentre(0.5 * imageSize.x, 0.5 * imageSize.y)
{
}
