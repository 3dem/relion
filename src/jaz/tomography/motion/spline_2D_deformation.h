#ifndef SPLINE_2D_DEFORMATION_H
#define SPLINE_2D_DEFORMATION_H

#include <src/jaz/gravis/t2Vector.h>
#include "modular_alignment/spline_2D_deformation_model.h"

class Spline2DDeformation
{
	public:
		
		Spline2DDeformation();
		
		Spline2DDeformation(
			gravis::i2Vector imageSize,
			gravis::i2Vector gridSize,
			const double* coefficients);
		
		
			Spline2DDeformationModel deformationModel;
			BufferedImage<gravis::d4Vector> coefficients;
			
		
		gravis::d2Vector apply(const gravis::d2Vector pl) const;
};

#endif
