#ifndef DEFORMATION_2D_H
#define DEFORMATION_2D_H

#include <src/jaz/gravis/t2Vector.h>

class Deformation2D
{
	public:
		
		struct Parameters
		{
			int grid_width, grid_height;
		};
		
		virtual gravis::d2Vector apply(const gravis::d2Vector pl) const = 0;
		virtual const double* getCoefficients() const = 0;
};

#endif
