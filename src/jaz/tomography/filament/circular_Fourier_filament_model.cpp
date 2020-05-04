#include "circular_Fourier_filament_model.h"
#include "coordinate_conventions.h"
#include <omp.h>

using namespace gravis;


CircularFourierFilamentModel::CircularFourierFilamentModel(
		double length, 
		const std::vector<d4Matrix>& proj,
		double padding )

:	length(length),
	padding(padding),
	proj(proj),
	totalLength(padding * length)
{	
}

void CircularFourierFilamentModel::apply(
		const std::vector<double>& params, 
		const std::vector<std::vector<float>>& encodedCoords, 
		std::vector<std::vector<float>>& signedDist_out,
		int num_threads) const
{
	const int fc = encodedCoords.size();
	
	/* 6 (real) dims for shift in 3D, 2 for radius: */	
	const int paramsDims = 8;
	const int paramsComplexDims = paramsDims / 2;
	const int paramFreqs = params.size() / paramsDims;
	
	signedDist_out.resize(fc);
	
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		const size_t pc = encodedCoords[f].size() / COORD_DIM;
		
		signedDist_out[f].resize(pc);
		
		for (size_t p = 0; p < pc; p++)			
		{
			const double t  = encodedCoords[f][COORD_DIM * p + COORD_T];
			const double q0 = encodedCoords[f][COORD_DIM * p + COORD_Q];
			
			const d2Vector d( encodedCoords[f][COORD_DIM * p + COORD_X],
			                  encodedCoords[f][COORD_DIM * p + COORD_Y]);
			
			const double phi = 1.0 * PI * t / totalLength;
			
			d4Vector x(0.0, 0.0, 0.0, 0.0);
			
			for (int i = 0; i < paramFreqs; i++)
			for (int j = 0; j < paramsComplexDims; j++)
			{
				x[j] += params[2 * (i*paramsComplexDims + j)    ] * sin( (i+1) * phi )
				      + params[2 * (i*paramsComplexDims + j) + 1] * cos( (i+1) * phi );
			}
			
			d4Vector d_3D(x[0], x[1], x[2], 0.0);
			d4Vector d_2D = proj[f] * d_3D;
			
			const double shift = d_2D.x * d.x + d_2D.y * d.y;
			const double grow = q0 > 0? -x[3] : x[3];
			
			signedDist_out[f][p] = q0 + shift + grow;
		}
	}
}

void CircularFourierFilamentModel::computeGradient(
		const std::vector<double>& params, 
		const std::vector<std::vector<float>>& encodedCoords,
		const std::vector<std::vector<float>>& slopes, 
		std::vector<double>& grad_out,
		int num_threads) const
{
	const int fc = encodedCoords.size();
	
	const int totParams = params.size();
	const int paramsDims = 8;
	const int paramsComplexDims = paramsDims / 2;
	const int paramFreqs = totParams / paramsDims;
	
	std::vector<std::vector<double>> grad_par(num_threads);
	
	for (int th = 0; th < num_threads; th++)
	{
		grad_par[th] = std::vector<double>(totParams, 0.0);
	}
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		int th = omp_get_thread_num();
		
		const size_t pc = encodedCoords[f].size() / COORD_DIM;
		
		for (size_t p = 0; p < pc; p++)			
		{
			const double t  = encodedCoords[f][COORD_DIM * p + COORD_T];
			const double q0 = encodedCoords[f][COORD_DIM * p + COORD_Q];
			
			const d2Vector d( encodedCoords[f][COORD_DIM * p + COORD_X],
			                  encodedCoords[f][COORD_DIM * p + COORD_Y]);
			
			const double phi = 1.0 * PI * t / totalLength;
			
			const d2Vector dd_2D_dx(proj[f](0,0), proj[f](1,0));
			const d2Vector dd_2D_dy(proj[f](0,1), proj[f](1,1));
			const d2Vector dd_2D_dz(proj[f](0,2), proj[f](1,2));
					
			const double dshift_dx = dd_2D_dx.x * d.x + dd_2D_dx.y * d.y;
			const double dshift_dy = dd_2D_dy.x * d.x + dd_2D_dy.y * d.y;
			const double dshift_dz = dd_2D_dz.x * d.x + dd_2D_dz.y * d.y;
			
			const double dgrow_dr = q0 > 0? -1.0 : 1.0;
						
			for (int i = 0; i < paramFreqs; i++)
			{
				const double slopeRe = cos( (i+1) * phi ) * slopes[f][p];
				const double slopeIm = sin( (i+1) * phi ) * slopes[f][p];
				
				grad_par[th][2*(i*paramsComplexDims + 0) + 0] += dshift_dx * slopeIm;
				grad_par[th][2*(i*paramsComplexDims + 0) + 1] += dshift_dx * slopeRe;
				
				grad_par[th][2*(i*paramsComplexDims + 1) + 0] += dshift_dy * slopeIm;
				grad_par[th][2*(i*paramsComplexDims + 1) + 1] += dshift_dy * slopeRe;
				
				grad_par[th][2*(i*paramsComplexDims + 2) + 0] += dshift_dz * slopeIm;
				grad_par[th][2*(i*paramsComplexDims + 2) + 1] += dshift_dz * slopeRe;
				
				grad_par[th][2*(i*paramsComplexDims + 3) + 0] += dgrow_dr  * slopeIm;				
				grad_par[th][2*(i*paramsComplexDims + 3) + 1] += dgrow_dr  * slopeRe;
			}
		}
	}
	
	grad_out = std::vector<double>(totParams, 0.0);
	
	for (int th = 0; th < num_threads; th++)
	for (int p = 0; p < totParams; p++)
	{
		grad_out[p] += grad_par[th][p];
	}
}
