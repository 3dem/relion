
#ifndef DAMAGE_H
#define DAMAGE_H

#include <string>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/image/buffered_image.h>
#include <vector>

class Damage
{
    public:

		static void applyWeight(
				RawImage<float>& stack,
				double pixelSize,
				const std::vector<double> &doses,
				int num_threads);

		static void applyWeight(
				RawImage<fComplex>& stack,
				double pixelSize,
				const std::vector<double> &doses,
				int num_threads);
		
		static std::vector<double> criticalDamageVector(double pixelSize, int boxSize);
		
		static std::vector<double> weightVector(double dose, double pixelSize, int boxSize);
		
		static BufferedImage<float> weightImage(double dose, double pixelSize, int boxSize);
		
		static BufferedImage<float> weightStack_GG(const std::vector<double>& doses,
										 double pixelSize, int boxSize);
		
		static double getWeight(double dose, double k);
		
		static std::vector<double> estimateDoseWeights(
				const BufferedImage<float>& tiltSeries,
				double pixelSize, 
				int tileSize,
				int firstFrame = -1, 
				double k0_Ang = 25.0, 
				bool diag = false);
		
		static gravis::d2Vector findDoseWeight(
			int f, int bestFrame, int tileSize,
			const std::vector<std::vector<double>>& powSum,
			std::vector<double>& Nk, double k0,
			double N0, double N1, double dN);
		
		
		
		static std::pair<std::vector<double>, std::vector<double>> 
			fitBkFactors(
				const BufferedImage<double>& fcc3, 
				double boxSize, 
				double pixelSize, 
				int k0, int k1, 
				BufferedImage<double>& plot,
				bool useL2);
				
		
		static BufferedImage<double> renderFit(
				const std::vector<double>& B, 
				const std::vector<double>& scale, 
				double boxSize, double pixelSize,
				int kc);
		
		static BufferedImage<double> renderFitNrm(
				const std::vector<double>& B, 
				const std::vector<double>& scale, 
				double boxSize, double pixelSize,
				const BufferedImage<double>& FCC3);
		
		
		
		static gravis::d2Vector findBkRec_L2(
				const BufferedImage<double>& fcc3, int f,
				double boxSize, double pixelSize,
				int k0, int k1,
				double B0, double B1,
				int steps, int depth, double q);
		
		static gravis::d2Vector findBkRec_FCC(
				const BufferedImage<double>& fcc3, int f,
				double boxSize, double pixelSize,
				int k0, int k1,
				double B0, double B1,
				int steps, int depth, double q);
		
		static double B_to_sigma(double B, double angbox);
		
		static double sigma_to_B(double sigma, double angbox);
		
		
		static BufferedImage<float> computeBfactorWeights(
				const std::vector<double>& bFacs, 
				const std::vector<double>& aFacs,
				int boxSize, int pixelSize, 
				bool normalize = false);
		
		static void renormalise(
				std::vector<std::vector<double>>& B_t,
				std::vector<std::vector<double>>& k_t);
};

#endif
