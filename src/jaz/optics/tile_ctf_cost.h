#ifndef TILE_CTF_COST_H
#define TILE_CTF_COST_H

#include <src/jaz/image/buffered_image.h>
#include "spectral_ctf_cost.h"


class TileCtfCost : public DifferentiableOptimization
{
	public:
		
		/*typedef enum 
		{
			PerFrame  = 0, 
			PerTile   = 1 
		} 
		Mode;
		
		typedef enum 
		{
			Defocus     = 0, 
			Astigmatism = 1, 
			Envelope    = 2, 
			Noise       = 3,
			ParamNum    = 4 
		} 
		Parameter;*/
		
		
		TileCtfCost(RawImage<float>& tiles, int tw, int th,
					RawImage<double>& tileOffZ,
					RawImage<float>& background, 
					double pixelSize, double voltage, double Cs, double Q0, 
					double hand,
					int num_threads, int r0, int r1);
		
				
			int tw, th, num_threads;
			double hand;
			std::vector<SpectralCtfCost> tileCosts;
			BufferedImage<double> tileOffZ;
			
			
		double f(const std::vector<double>& x, void* tempStorage) const;
		void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;
		
		std::vector<double> getInitialParams(double defocus);
		std::vector<double> getInitialParams(const std::vector<double>& global);
		
		BufferedImage<float> renderTable(const std::vector<double>& x);
		BufferedImage<float> renderCentralFit(const std::vector<double>& x);
		BufferedImage<float> renderAlignedSpectrum(const std::vector<double>& x);
		BufferedImage<float> renderAvgSpectrum();
		
		std::vector<double> radialAverageAligned(const std::vector<double>& x);
		std::vector<double> radialAverage(const std::vector<double>& x);
		std::vector<double> radialAverageFit(const std::vector<double>& x);
			
		
};

#endif
