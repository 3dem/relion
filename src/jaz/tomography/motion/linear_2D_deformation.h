#ifndef LINEAR_2D_DEFORMATION_H
#define LINEAR_2D_DEFORMATION_H

#include <src/jaz/gravis/t3Vector.h>
#include "modular_alignment/linear_2D_deformation_model.h"
#include "2D_deformation.h"

class Linear2DDeformation : public Deformation2D
{
	public:
		
		Linear2DDeformation();
		
		Linear2DDeformation(
			gravis::i2Vector imageSize,
			const double* coefficients);
		
		
			Linear2DDeformationModel deformationModel;
			gravis::d3Vector coefficients;
			
		
		gravis::d2Vector apply(const gravis::d2Vector pl) const;
		const double* getCoefficients() const;
};

#endif
