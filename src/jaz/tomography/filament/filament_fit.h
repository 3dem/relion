#ifndef FILAMENT_FIT_H
#define FILAMENT_FIT_H


#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/optimization/optimization.h>
#include <map>

#include <src/jaz/math/spline.h>
#include <src/jaz/image/buffered_image.h>

#include "filament_model.h"
#include "filament_mapping.h"
#include "pixel_set.h"

class Filament;

class FilamentFit : public DifferentiableOptimization
{
	public:
		
		
		
		FilamentFit(
				const Filament& filament, 
				const FilamentMapping& mapping,
				const std::vector<gravis::d4Matrix>& proj,
				const RawImage<float> binnedStack, 
				const RawImage<float> outlier_mask, 
				const FilamentModel* model,
				double binning, 
				int num_threads);
		
		
			const Filament& filament;
			const FilamentMapping& mapping;	
			const std::vector<gravis::d4Matrix>& proj;		
			const RawImage<float> binnedStack, outlier_mask;
			
			const FilamentModel* model;
			double maxDistBinned, binning;
			
			PixelSet optimisationPixels, allPixels;
			
			int num_threads;
			
			
		double f(
				const std::vector<double>& x, 
				void* tempStorage) const;
		
		void grad(
				const std::vector<double>& x, 
				std::vector<double>& gradDest, 
				void* tempStorage) const;		
		
		void report(int iteration, double cost, const std::vector<double>& x) const; 
		
		
		BufferedImage<float> visualise(
				const std::vector<double>& x, 
				bool subtract = true,
				BufferedImage<float>* original = 0,
				float scale = 1.f);
		
		BufferedImage<float> computeCostByOffset(
				const std::vector<double>& x,
				double minOffset,
				double maxOffset,
				int samples, 
				double sigma);
		
		
	protected:
				
		void encode(const BufferedImage<float>& mask, PixelSet& out);
		
		void averageValues(
				const std::vector<float>& phases,
				const std::vector<float>& pixelValues,  
				std::vector<float>& out) const;
		
		void expandValues(
				const std::vector<float>& phases, 
				const std::vector<float>& phaseValues, 
				std::vector<float>& out) const;
		
		void getSlopes(
				const std::vector<float>& phases, 
				const std::vector<float>& phaseValues, 
				std::vector<float>& slopes_out) const;		
		
		BufferedImage<float> costByOffset(
				const std::vector<float>& phases,       // pc
				const std::vector<float>& phaseValues,  // phc
				const std::vector<float>& arcCoords,    // pc
				const std::vector<float>& pixelValues,  // pc
				double minOffset,
				double maxOffset,
				int samples,
				double arcLength,
				double arcScale, 
				double sigma) const;
		
		double computeZeroLevel(const std::vector<float>& phaseValues) const;
		void correctZeroLevel(std::vector<float>& phaseValues) const;
};

#endif
