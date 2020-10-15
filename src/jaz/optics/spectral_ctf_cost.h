#ifndef SPECTRAL_CTF_FIT_H
#define SPECTRAL_CTF_FIT_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/image/buffered_image.h>
#include <src/ctf.h>


/*
  TODO: add class SpectralTiltCtfFit that shuffles parameters to 
  individual SpectralCtfCost instances, allowing for:
	- (k, B) and gamma per tile or per frame
	- alpha and beta per frame or global
*/
  
class SpectralCtfCost : public DifferentiableOptimization
{
	public:
		
		SpectralCtfCost();
		
		SpectralCtfCost(
				RawImage<float> spectrum,
				RawImage<float> background,
				double pixelSize, double voltage, double Cs, double Q0,
				int num_threads, int k0, int k1);
		
			RawImage<float> spectrum, background;
			int num_threads, k0, k1;
			double pixelSize, voltage, Cs, Q0;
			double lambda;
			
			double R1, R2, R3_5; // correspond to K1..K5 in Relion's CTF, but in pixels
		
			
		double f(const std::vector<double>& x, void* tempStorage) const;
		void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;

		void report(int iteration, double cost, const std::vector<double>& x) const;
		
		std::vector<double> getInitialParams(double defocus);
		//CTF getCtf(const std::vector<double>& x);
		
		BufferedImage<float> render(const std::vector<double>& x);
		void addAlignedSpectrum(double z0, const std::vector<double>& x, 
								BufferedImage<double>& accum, BufferedImage<double>& weight);
		
		void addAlignedRadialAverage(double z0, const std::vector<double>& x, 
							   std::vector<double>& accum, std::vector<double>& weight);
		
		void addRadialAverage(std::vector<double>& accum, std::vector<double>& weight);
		
		void addRadialAverageFit(const std::vector<double>& x, 
							   std::vector<double>& accum, std::vector<double>& weight);
		
		double offsetDefocusParam(double x0, double deltaZ_Ang) const;
				
};

#endif
