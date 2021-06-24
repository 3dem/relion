#ifndef LINEAR_2D_DEFORMATION_MODEL_H
#define LINEAR_2D_DEFORMATION_MODEL_H

#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/tomography/motion/2D_deformation.h>


class Linear2DDeformationModel
{
	public:
		
		
		Linear2DDeformationModel();
		
		Linear2DDeformationModel(
				const Deformation2D::Parameters& parameters,
				const gravis::i2Vector& imageSize);
		
		
			gravis::i2Vector imageSize;
			gravis::d2Vector imageCentre;
		
		
		inline int getParameterCount() const;
		
		inline gravis::d2Vector computeShift(
				const gravis::d2Vector& pl,
				const double* parameters) const;
		
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



inline int Linear2DDeformationModel::getParameterCount() const
{
	return 3;
}

inline gravis::d2Vector Linear2DDeformationModel::computeShift(
		const gravis::d2Vector& pl,
		const double* parameters) const
{	
	const double axx = parameters[0];
	const double axy = parameters[1];
	const double ayy = parameters[2];
	
	const gravis::d2Vector r = pl - imageCentre;
	
	return gravis::d2Vector(
		axx * r.x,
		axy * r.x + ayy * r.y);
}

inline void Linear2DDeformationModel::computeShiftAndGradient(
		const gravis::d2Vector& pl,
		const double* parameters,
		gravis::d2Vector& def,
		gravis::d2Vector& def_x,
		gravis::d2Vector& def_y) const
{
	const double axx = parameters[0];
	const double axy = parameters[1];
	const double ayy = parameters[2];
	
	const gravis::d2Vector r = pl - imageCentre;
	
	def = gravis::d2Vector(
		axx * r.x,
		axy * r.x + ayy * r.y);
	
	def_x = gravis::d2Vector(
		axx,
		axy);
	
	def_y = gravis::d2Vector(
		0.0,
		ayy);
}

inline gravis::d2Vector Linear2DDeformationModel::transformImageGradient(
		const gravis::d2Vector &g0,
		const gravis::d2Vector &mx,
		const gravis::d2Vector &my) const
{
	return gravis::d2Vector (
				(mx.x + 1.0) * g0.x  +        mx.y  * g0.y,
				       my.x  * g0.x  + (my.y + 1.0) * g0.y );
}

inline void Linear2DDeformationModel::updateDataTermGradient(
		const gravis::d2Vector &pl,
		const gravis::d2Vector &g0,
		const double *parameters,
		double *target) const
{
	const gravis::d2Vector r = pl - imageCentre;
	
	target[0] += r.x * g0.x;
	target[1] += r.x * g0.y;
	target[2] += r.y * g0.y;
}

inline double Linear2DDeformationModel::computePriorValueAndGradient(
		double lambda,
		const double* parameters,
		double* gradDest) const
{
	const double axx = parameters[0];
	const double axy = parameters[1];
	const double ayy = parameters[2];
	
	const double out = lambda * (axx * axx + axy * axy + ayy * ayy);
	
	gradDest[0] += 2.0 * lambda * axx;
	gradDest[1] += 2.0 * lambda * axy;
	gradDest[2] += 2.0 * lambda * ayy;
	
	return out;
}

#endif
	
