
#ifndef TOMO_CTF_FIND_H
#define TOMO_CTF_FIND_H

#include <string>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/cutting.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/error.h>

class TomoCtfFind
{
    public:

        TomoCtfFind(const BufferedImage<float>& tiltSeries, 
					const std::vector<gravis::d4Matrix>& projections,
					bool diag,
					double pixelSize, double voltage, double Cs, double Q0 = 0.07,
					int tileSize = 512, double r0_ang = 50.0, double r1_ang = 10.0,
					int firstFrame = -1);
		
		
            const BufferedImage<float>& tiltSeries;
			const std::vector<gravis::d4Matrix>& projections;
			
			bool diag;
            double pixelSize, voltage, Cs, Q0, lambda;
			double r0_ang, r1_ang;
			int tileSize, firstFrame;
			
		
			
		void initialize(int num_threads);
		
		void findInitialCtf(double z0, double z1, double dz);
		
		std::vector<double> findDefoci(double z0, double z1, double dz);
		void findAstigmatism();
		
		double evaluateAstigmatic(double deltaF, double a0, double a1, int frame);
		double evaluate1D(double deltaF, int frame, double hand);
				
		
	protected:
		
		
			BufferedImage<float> tiles;
			std::vector<BufferedImage<float>> avgTiles, backgrounds;
			BufferedImage<float> avgTile, avgTileBgSub, avgBackground;
			
			int tw, th, fc;
			double tdx, tdy;
			int r0_pix, r1_pix, r0ref_pix, r1ref_pix;
			
			int num_threads;
			
			BufferedImage<double> tileOffZ;
			std::vector<double> frameDose;
			std::vector<std::vector<double>> damageWeight;
			BufferedImage<std::vector<float>> radAvg, radVar, radAvgNrm;
			std::vector<double> radAvgAllNrm;
			std::vector<double> frqWgh;
			
			std::vector<double> globalParams;
		
		
		void extractTiles();
		void averageTiles();
		void subtractBackground();
		void averageRadially();
		
		BufferedImage<float> estimateBackground(
				const RawImage<float>& img_fftwHalf, 
				double sigma_blur_px, 
				double sigma_cent_px, 
				double sigma_out_px, 
				bool debug);
};

class AstigmatismOptimization : public Optimization
{
    public:
		
		AstigmatismOptimization();
		
		double f(const std::vector<double>& x, void* tempStorage) const;
		void report(int iteration, double cost, const std::vector<double>& x) const {}
};

#endif
