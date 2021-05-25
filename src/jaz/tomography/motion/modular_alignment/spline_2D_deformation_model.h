#ifndef SPLINE_DEFORMATION_MODEL_H
#define SPLINE_DEFORMATION_MODEL_H

#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/image/raw_image.h>


class Spline2DDeformationModel
{
	public:
		
		struct DataPoint
		{
			double value, slope_x, slope_y, twist;
		};

		struct Parameters
		{
			int grid_width, grid_height;
		};
		
		Spline2DDeformationModel(
				const Parameters& parameters,
				const gravis::i2Vector& imageSize);
		
		
			gravis::i2Vector imageSize, gridSize;
			gravis::d2Vector gridSpacing;
		
		
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
		
		inline void updateCostGradient(
				const gravis::d2Vector& pl,
				const gravis::d2Vector& g0,
				const double* parameters,
				double* target) const;
		
};

inline int Spline2DDeformationModel::getParameterCount() const
{
	return 8 * gridSize.x * gridSize.y;
}

inline void Spline2DDeformationModel::computeShiftAndGradient(
		const gravis::d2Vector& pl,
		const double* parameters,
		gravis::d2Vector& def,
		gravis::d2Vector& def_x,
		gravis::d2Vector& def_y) const
{
	gravis::d2Vector pl_grid(pl.x / gridSpacing.x, pl.y / gridSpacing.y);
	
	const double eps = 1e-10;
	
	for (int dim = 0; dim < 2; dim++)
	{
		if (pl_grid[dim] < 0.0)
		{
			pl_grid[dim] = 0.0;
		}
		else if (pl_grid[dim] > gridSize[dim] - 1 - eps)
		{
			pl_grid[dim] = gridSize[dim] - 1 - eps;
		}
	}
	
	const gravis::i2Vector cell((int)pl_grid.x, (int)pl_grid.y);
	
	const gravis::d2Vector frc(
			pl_grid.x - cell.x,
			pl_grid.y - cell.y);
	
	const double x  = frc.x;
	const double x2 = x * x;
	const double x3 = x2 * x;
	
	const double y  = frc.y;
	const double y2 = y * y;
	const double y3 = y2 * y;
	
	const RawImage<DataPoint> data(gridSize.x, gridSize.y, 2, (DataPoint*)parameters);
		
	for (int dim = 0; dim < 2; dim++)
	{		
		DataPoint d00 = data(cell.x,     cell.y    , dim);
		DataPoint d01 = data(cell.x,     cell.y + 1, dim);
		DataPoint d10 = data(cell.x + 1, cell.y    , dim);
		DataPoint d11 = data(cell.x + 1, cell.y + 1, dim);
		
		const gravis::d4Matrix F
			( d00.value,   d01.value,   d00.slope_y, d01.slope_y, 
			  d10.value,   d11.value,   d10.slope_y, d11.slope_y, 
			  d00.slope_x, d01.slope_x, d00.twist,   d01.twist, 
			  d10.slope_x, d11.slope_x, d10.twist,   d11.twist );
		
		const gravis::d4Vector vx(
			1.0      -  3.0 * x2  +  2.0 * x3,
			            3.0 * x2  -  2.0 * x3,
			      x  -  2.0 * x2  +        x3,
			                - x2  +        x3);
		
		const gravis::d4Vector vy(
			1.0      -  3.0 * y2  +  2.0 * y3,
			            3.0 * y2  -  2.0 * y3,
			      y  -  2.0 * y2  +        y3,
			                - y2  +        y3);
		
		const gravis::d4Vector dx(
			      -  6.0 * x  +  6.0 * x2,
			         6.0 * x  -  6.0 * x2,
			 1.0  -  4.0 * x  +  3.0 * x2,
			        -2.0 * x  +  3.0 * x2);
		
		const gravis::d4Vector dy(
			      -  6.0 * y  +  6.0 * y2,
			         6.0 * y  -  6.0 * y2,
			 1.0  -  4.0 * y  +  3.0 * y2,
			        -2.0 * y  +  3.0 * y2);
		
		def[dim]   = vx.dot(F * vy);
		
		def_x[dim] = dx.dot(F * vy) / gridSpacing.x;
		def_y[dim] = vx.dot(F * dy) / gridSpacing.y;
	}
}

inline gravis::d2Vector Spline2DDeformationModel::transformImageGradient(
		const gravis::d2Vector &g0,
		const gravis::d2Vector &mx,
		const gravis::d2Vector &my) const
{
	return gravis::d2Vector (
				(mx.x + 1.0) * g0.x  +        mx.y  * g0.y,
				       my.x  * g0.x  + (my.y + 1.0) * g0.y );
}

inline void Spline2DDeformationModel::updateCostGradient(
		const gravis::d2Vector &pl,
		const gravis::d2Vector &g0,
		const double *parameters,
		double *target) const
{
	gravis::d2Vector pl_grid(pl.x / gridSpacing.x, pl.y / gridSpacing.x);
	
	const double eps = 1e-10;
	
	for (int dim = 0; dim < 2; dim++)
	{
		if (pl_grid[dim] < 0.0)
		{
			pl_grid[dim] = 0.0;
		}
		else if (pl_grid[dim] > gridSize[dim] - 1 - eps)
		{
			pl_grid[dim] = gridSize[dim] - 1 - eps;
		}
	}
	
	const gravis::i2Vector cell((int)pl_grid.x, (int)pl_grid.y);
	
	const gravis::d2Vector frc(
			pl_grid.x - cell.x,
			pl_grid.y - cell.y);
	
	const double x  = frc.x;
	const double x2 = x * x;
	const double x3 = x2 * x;
	
	const double y  = frc.y;
	const double y2 = y * y;
	const double y3 = y2 * y;
	
	const gravis::d4Vector vx(
		1.0      -  3.0 * x2  +  2.0 * x3,
					3.0 * x2  -  2.0 * x3,
			  x  -  2.0 * x2  +        x3,
						- x2  +        x3);
	
	const gravis::d4Vector vy(
		1.0      -  3.0 * y2  +  2.0 * y3,
					3.0 * y2  -  2.0 * y3,
			  y  -  2.0 * y2  +        y3,
						- y2  +        y3);
	
	RawImage<DataPoint> grad(gridSize.x, gridSize.y, 2, (DataPoint*)target);
	
	for (int dim = 0; dim < 2; dim++)
	{
		grad(cell.x,     cell.y,     dim).value   += vx[0] * vy[0] * g0[dim];
		grad(cell.x,     cell.y,     dim).slope_x += vx[2] * vy[0] * g0[dim];
		grad(cell.x,     cell.y,     dim).slope_y += vx[0] * vy[2] * g0[dim];
		grad(cell.x,     cell.y,     dim).twist   += vx[2] * vy[2] * g0[dim];
		
		grad(cell.x + 1, cell.y,     dim).value   += vx[1] * vy[0] * g0[dim];
		grad(cell.x + 1, cell.y,     dim).slope_x += vx[3] * vy[0] * g0[dim];
		grad(cell.x + 1, cell.y,     dim).slope_y += vx[1] * vy[2] * g0[dim];
		grad(cell.x + 1, cell.y,     dim).twist   += vx[2] * vy[2] * g0[dim];
		
		grad(cell.x,     cell.y + 1, dim).value   += vx[0] * vy[1] * g0[dim];
		grad(cell.x,     cell.y + 1, dim).slope_x += vx[2] * vy[1] * g0[dim];
		grad(cell.x,     cell.y + 1, dim).slope_y += vx[0] * vy[3] * g0[dim];
		grad(cell.x,     cell.y + 1, dim).twist   += vx[3] * vy[3] * g0[dim];
		
		grad(cell.x + 1, cell.y + 1, dim).value   += vx[1] * vy[1] * g0[dim];
		grad(cell.x + 1, cell.y + 1, dim).slope_x += vx[3] * vy[1] * g0[dim];
		grad(cell.x + 1, cell.y + 1, dim).slope_y += vx[1] * vy[3] * g0[dim];
		grad(cell.x + 1, cell.y + 1, dim).twist   += vx[3] * vy[3] * g0[dim];
	}	
}

#endif
