#ifndef CTF_EQUIPHASE_FIT_H
#define CTF_EQUIPHASE_FIT_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/image/buffered_image.h>
#include <src/ctf.h>

class CtfEquiphaseFit : public Optimization
{
	public:
		
		CtfEquiphaseFit(
				const BufferedImage<float>& spectrum,
				double pixelSize, double voltage, double Cs, double Q0,
				double mapScale,
				int num_threads, int k0, int k1);
		
			const BufferedImage<float>& spectrum;
			int num_threads, k0, k1;
			double pixelSize, voltage, Cs, Q0;
			double mapScale;
			CTF exampleCtf;
		
			
		double f(const std::vector<double>& x, void* tempStorage) const;
		void report(int iteration, double cost, const std::vector<double>& x) const;
		
		CTF paramsToCtf(const std::vector<double>& x) const;
		
		
		void averageAlongIso(
				const CTF& ctf, 
				std::vector<double>& accum, 
				std::vector<double>& weight,
				int boundary = 2) const;
		
		double compareWithExpansion(
				const CTF& ctf,
				const std::vector<double>& accum,
				int boundary = 2) const;
		
		BufferedImage<double> computeExpansion(
				const CTF& ctf, 
				const std::vector<double>& accum,
				int boundary = 2) const;
};

#endif
