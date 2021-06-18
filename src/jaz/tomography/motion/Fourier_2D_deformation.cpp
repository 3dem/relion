#include "Fourier_2D_deformation.h"

using namespace gravis;

Fourier2DDeformation::Fourier2DDeformation()
{
	
}

Fourier2DDeformation::Fourier2DDeformation(
		i2Vector imageSize, 
		i2Vector gridSize, 
		const double* coefficients)
{
	Deformation2D::Parameters parameters;
	
	parameters.grid_width  = gridSize.x;
	parameters.grid_height = gridSize.y;
	
	deformationModel = Fourier2DDeformationModel(parameters, imageSize);
	
	RawImage<dComplex> asImage(gridSize.x/2 + 1, gridSize.y, 2, (dComplex*)coefficients);
	this->coefficients = asImage;
}

d2Vector Fourier2DDeformation::apply(const d2Vector pl) const
{
	return pl + deformationModel.computeShift(pl, (double*)coefficients.data);
}

const double* Fourier2DDeformation::getCoefficients() const
{
	return (double*) &coefficients[0].real;
}
