#ifndef FILAMENT_MODEL_H
#define FILAMENT_MODEL_H

#include <src/jaz/image/buffered_image.h>

class FilamentModel
{
	public:
		
		// constructs a displacement field from a set of parameters
		
		virtual void apply(
				const std::vector<double>& params,
				const std::vector<std::vector<float>>& encodedCoords,
				std::vector<std::vector<float>>& signedDist_out,
				int num_threads) const = 0;
		
		virtual void computeGradient(
				const std::vector<double>& params,
				const std::vector<std::vector<float>>& encodedCoords,
				const std::vector<std::vector<float>>& slopes,
				std::vector<double>& grad_out,
				int num_threads) const = 0;
};

class DummyFilamentModel : public FilamentModel
{
	public:
		
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
