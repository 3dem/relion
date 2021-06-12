#ifndef FOURIER_DEFORMATION_MODEL_H
#define FOURIER_DEFORMATION_MODEL_H

#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/motion/2D_deformation.h>


class Fourier2DDeformationModel
{
	public:
		
		Fourier2DDeformationModel();
		
		Fourier2DDeformationModel(
				const Deformation2D::Parameters& parameters,
				const gravis::i2Vector& imageSize);
		
		
			gravis::i2Vector imageSize, gridSize;
			std::vector<gravis::d2Vector> spatialFrequencies;
		
		
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



inline int Fourier2DDeformationModel::getParameterCount() const
{
	return 4 * spatialFrequencies.size();
}

inline gravis::d2Vector Fourier2DDeformationModel::computeShift(
		const gravis::d2Vector& pl,
		const double* parameters) const
{
	const RawImage<dComplex> P(spatialFrequencies.size(), 2, 1, (dComplex*)parameters);
	
	gravis::d2Vector def(0.0, 0.0);
		
	for (int dim = 0; dim < 2; dim++)
	{	
		for (int i = 0; i < spatialFrequencies.size(); i++)
		{
			const double t = spatialFrequencies[i].dot(pl);
			
			def[dim] += P(i,dim).real * cos(t) + P(i,dim).imag * sin(t);
		}
	}
	
	return def;
}

inline void Fourier2DDeformationModel::computeShiftAndGradient(
		const gravis::d2Vector& pl,
		const double* parameters,
		gravis::d2Vector& def,
		gravis::d2Vector& def_x,
		gravis::d2Vector& def_y) const
{
	const RawImage<dComplex> P(spatialFrequencies.size(), 2, 1, (dComplex*)parameters);
	
	for (int dim = 0; dim < 2; dim++)
	{	
		def[dim]   = 0.0;
		def_x[dim] = 0.0;
		def_y[dim] = 0.0;

		for (int i = 0; i < spatialFrequencies.size(); i++)
		{
			const double t = spatialFrequencies[i].dot(pl);
			
			const dComplex z = P(i,dim);
			
			def[dim] += z.real * cos(t) + z.imag * sin(t);
			
			const double def_t = z.real * (-sin(t)) + z.imag * cos(t);
			
			def_x[dim] += def_t * spatialFrequencies[i].x;
			def_y[dim] += def_t * spatialFrequencies[i].y;
		}
	}
}

inline gravis::d2Vector Fourier2DDeformationModel::transformImageGradient(
		const gravis::d2Vector &g0,
		const gravis::d2Vector &mx,
		const gravis::d2Vector &my) const
{
	return gravis::d2Vector (
				(mx.x + 1.0) * g0.x  +        mx.y  * g0.y,
				       my.x  * g0.x  + (my.y + 1.0) * g0.y );
}

inline void Fourier2DDeformationModel::updateDataTermGradient(
		const gravis::d2Vector &pl,
		const gravis::d2Vector &g0,
		const double *parameters,
		double *target) const
{
	RawImage<dComplex> grad(spatialFrequencies.size(), 2, 1, (dComplex*)target);
	
	for (int dim = 0; dim < 2; dim++)
	{
		for (int i = 0; i < spatialFrequencies.size(); i++)
		{
			const double t = spatialFrequencies[i].dot(pl);
			
			grad(i,dim).real += cos(t) * g0[dim];
			grad(i,dim).imag += sin(t) * g0[dim];
		}
	}
}

inline double Fourier2DDeformationModel::computePriorValueAndGradient(
		double lambda,
		const double* parameters,
		double* gradDest) const
{
	const RawImage<dComplex> data(spatialFrequencies.size(), 2, 1, (dComplex*)parameters);
	RawImage<dComplex> grad(spatialFrequencies.size(), 2, 1, (dComplex*)gradDest);
	
	double out = 0.0;
	
	for (int dim = 0; dim < 2; dim++)
	{
		for (int i = 0; i < spatialFrequencies.size(); i++)
		{
			const double e = lambda * spatialFrequencies[i].norm2();
			
			out += e * data(i,dim).norm();
			
			grad(i,dim).real += 2.0 * e * data(i,dim).real;
			grad(i,dim).imag += 2.0 * e * data(i,dim).imag;
		}
	}
	
	return out;
}

#endif
