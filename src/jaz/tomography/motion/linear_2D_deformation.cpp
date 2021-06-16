#include "linear_2D_deformation.h"

using namespace gravis;

Linear2DDeformation::Linear2DDeformation()
{
	
}

Linear2DDeformation::Linear2DDeformation(
		i2Vector imageSize, 
		const double* coefficients)
{
	Deformation2D::Parameters parameters;
	
	parameters.grid_width  = 1;
	parameters.grid_height = 1;
	
	deformationModel = Linear2DDeformationModel(parameters, imageSize);
	
	this->coefficients = d3Vector(coefficients[0], coefficients[1], coefficients[2]);
}

d2Vector Linear2DDeformation::apply(const d2Vector pl) const
{
	return pl + deformationModel.computeShift(pl, &coefficients[0]);
}

const double* Linear2DDeformation::getCoefficients() const
{
	return &coefficients[0];
}
