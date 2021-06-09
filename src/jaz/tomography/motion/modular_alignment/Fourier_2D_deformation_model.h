#ifndef FOURIER_DEFORMATION_MODEL_H
#define FOURIER_DEFORMATION_MODEL_H

#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/buffered_image.h>


class Fourier2DDeformationModel
{
	public:
		
		struct Parameters
		{
			int grid_width, grid_height;
		};
		
		Fourier2DDeformationModel();
		
		Fourier2DDeformationModel(
				const Parameters& parameters,
				const gravis::i2Vector& imageSize);
		
		
			gravis::i2Vector imageSize, gridSize;
			gravis::d2Vector gridSpacing;
		
		
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
		
		inline void updateCostGradient(
				const gravis::d2Vector& pl,
				const gravis::d2Vector& g0,
				const double* parameters,
				double* target) const;
};



inline int Fourier2DDeformationModel::getParameterCount() const
{
	return 4 * (gridSize.x/2 + 1) * gridSize.y;
}

inline gravis::d2Vector Fourier2DDeformationModel::computeShift(
		const gravis::d2Vector& pl,
		const double* parameters) const
{
	const RawImage<dComplex> data((gridSize.x/2 + 1), gridSize.y, 2, (dComplex*)parameters);
	
	gravis::d2Vector def;
		
	for (int dim = 0; dim < 2; dim++)
	{	
		for (int y = 0; y < gridSize.y; y++)
		for (int x = (y < gridSize.y && y > 0? 0 : 1); x < gridSize.x; x++)
		{
			const gravis::d2Vector d(
				x * PI / imageSize.x,
				y * PI / imageSize.y);
			
			const double t = d.dot(pl);
			
			def[dim] += data(x,y,dim).real * cos(t) + data(x,y,dim).imag * sin(t);
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
	const RawImage<dComplex> data((gridSize.x/2 + 1), gridSize.y, 2, (dComplex*)parameters);
		
	for (int dim = 0; dim < 2; dim++)
	{	
		for (int y = 0; y < gridSize.y; y++)
		for (int x = (y < gridSize.y && y > 0? 0 : 1); x < gridSize.x; x++)
		{
			const gravis::d2Vector d(
				x * PI / imageSize.x,
				y * PI / imageSize.y);
			
			const double t = d.dot(pl);
			
			const dComplex z = data(x,y,dim);
			
			def[dim] += z.real * cos(t) + z.imag * sin(t);
			
			const double def_t = z.real * (-sin(t)) + z.imag * cos(t);
			
			def_x[dim] += def_t * d.x;
			def_y[dim] += def_t * d.y;
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

inline void Fourier2DDeformationModel::updateCostGradient(
		const gravis::d2Vector &pl,
		const gravis::d2Vector &g0,
		const double *parameters,
		double *target) const
{
	RawImage<dComplex> grad((gridSize.x/2 + 1), gridSize.y, 2, (dComplex*)target);
	
	for (int dim = 0; dim < 2; dim++)
	{	
		for (int y = 0; y < gridSize.y; y++)
		for (int x = (y < gridSize.y && y > 0? 0 : 1); x < gridSize.x; x++)
		{
			const gravis::d2Vector d(
				x * PI / imageSize.x,
				y * PI / imageSize.y);
			
			const double t = d.dot(pl);
			
			grad(x,y,dim).real += cos(t) * g0[dim];
			grad(x,y,dim).imag += sin(t) * g0[dim];
		}
	}
}

#endif
