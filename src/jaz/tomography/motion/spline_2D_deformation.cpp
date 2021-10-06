#include "spline_2D_deformation.h"

using namespace gravis;

Spline2DDeformation::Spline2DDeformation()
{
	
}

Spline2DDeformation::Spline2DDeformation(
		i2Vector imageSize, 
		i2Vector gridSize, 
		const double* coefficients)
{
	Deformation2D::Parameters parameters;
	
	parameters.grid_width  = gridSize.x;
	parameters.grid_height = gridSize.y;
	
	deformationModel = Spline2DDeformationModel(
		parameters, imageSize);
	
	RawImage<d4Vector> asImage(gridSize.x, gridSize.y, 2, (d4Vector*)coefficients);
	this->coefficients = asImage;
}

d2Vector Spline2DDeformation::apply(const d2Vector pl) const
{
	return pl + deformationModel.computeShift(pl, (double*)coefficients.data);
}

const double* Spline2DDeformation::getCoefficients() const
{
	return (double*) &coefficients[0][0];
}
