#include "filament_model.h"
#include "coordinate_conventions.h"
#include <omp.h>


using namespace gravis;


void DummyFilamentModel::apply(
		const std::vector<double>& params, 
		const std::vector<std::vector<float>>& encodedCoords, 
		std::vector<std::vector<float>>& signedDist_out,
		int num_threads) const
{
	const int fc = encodedCoords.size();
	
	signedDist_out.resize(fc);			
	
	for (int f = 0; f < fc; f++)
	{
		const size_t pc = encodedCoords[f].size() / COORD_DIM;
		
		signedDist_out[f].resize(pc);
		
		for (size_t p = 0; p < pc; p++)
		{
			signedDist_out[f][p] = encodedCoords[f][COORD_DIM * p + COORD_Q];
		}
	}
}

void DummyFilamentModel::computeGradient(
		const std::vector<double>& params,
		const std::vector<std::vector<float>>& encodedCoords,
		const std::vector<std::vector<float>>& slopes,
		std::vector<double>& grad_out,
		int num_threads) const
{
	grad_out.resize(0);
}
