#ifndef NO_DEFORMATION_MODEL_H
#define NO_DEFORMATION_MODEL_H

#include <src/jaz/gravis/t2Vector.h>

class No2DDeformationModel
{
	public:

		inline int getParameterCount() const;
		
		inline void computeShiftAndGradient(
				const gravis::d2Vector& pl,
				const double* parameters,
				gravis::d2Vector& def,
				gravis::d2Vector& def_x,
				gravis::d2Vector& def_y) const;
		
		inline gravis::d2Vector transformImageGradient(
				const gravis::d2Vector& g0,
				const gravis::d2Vector& mx,
				const gravis::d2Vector& my) const;
		
		inline void updateDataTermGradient(
				const gravis::d2Vector& pl,
				const gravis::d2Vector& g0,
				const double* parameters,
				double* target) const;
		
		inline double computePriorValueAndGradient(
				double lambda,
				const double* parameters,
				double* gradDest) const;
		
};

inline int No2DDeformationModel::getParameterCount() const
{
	return 0;
}

inline void No2DDeformationModel::computeShiftAndGradient(
		const gravis::d2Vector &pl,
		const double* parameters,
		gravis::d2Vector &def,
		gravis::d2Vector &def_x,
		gravis::d2Vector &def_y) const
{
	def = gravis::d2Vector(0.0, 0.0);
	def_x = gravis::d2Vector(0.0, 0.0);
	def_y = gravis::d2Vector(0.0, 0.0);
}

inline gravis::d2Vector No2DDeformationModel::transformImageGradient(
		const gravis::d2Vector &g0,
		const gravis::d2Vector &mx,
		const gravis::d2Vector &my) const
{
	return g0;
}

inline void No2DDeformationModel::updateDataTermGradient(
		const gravis::d2Vector &pl,
		const gravis::d2Vector &g0,
		const double *parameters,
		double *target) const
{
	
}

double No2DDeformationModel::computePriorValueAndGradient(
		double lambda, 
		const double* parameters, 
		double* gradDest) const
{
	return 0.0;
}

#endif
