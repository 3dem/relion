#ifndef FOURIER_2D_DEFORMATION_H
#define FOURIER_2D_DEFORMATION_H

#include <src/jaz/gravis/t2Vector.h>
#include "modular_alignment/Fourier_2D_deformation_model.h"
#include "2D_deformation.h"

class Fourier2DDeformation : public Deformation2D
{
	public:
		
		Fourier2DDeformation();
		
		Fourier2DDeformation(
			gravis::i2Vector imageSize,
			gravis::i2Vector gridSize,
			const double* coefficients);
		
		
			Fourier2DDeformationModel deformationModel;
			BufferedImage<dComplex> coefficients;
			
		
		gravis::d2Vector apply(const gravis::d2Vector pl) const;
		const double* getCoefficients() const;
};

#endif
