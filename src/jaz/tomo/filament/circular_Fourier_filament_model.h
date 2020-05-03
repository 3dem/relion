#ifndef CIRCULAR_FOURIER_FILAMENT_MODEL_H
#define CIRCULAR_FOURIER_FILAMENT_MODEL_H

#include <src/jaz/image/buffered_image.h>
#include "filament_model.h"

class CircularFourierFilamentModel : public FilamentModel
{
	public:
		
		CircularFourierFilamentModel(
				double length,
				const std::vector<gravis::d4Matrix>& proj,
				double padding = 2.0);
		
		
			double length, padding, totalLength;
			const std::vector<gravis::d4Matrix>& proj;
		
		
		void apply(
				const std::vector<double>& params,
				const std::vector<std::vector<float>>& encodedCoords,
				std::vector<std::vector<float>>& signedDist_out,
				int num_threads) const;
		
		void computeGradient(
				const std::vector<double>& params,
				const std::vector<std::vector<float>>& encodedCoords,
				const std::vector<std::vector<float>>& slopes,
				std::vector<double>& grad_out,
				int num_threads) const;
};


#endif
