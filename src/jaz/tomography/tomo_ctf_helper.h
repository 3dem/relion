#ifndef TOMO_CTF_HELPER_H
#define TOMO_CTF_HELPER_H

#include <src/ctf.h>
#include <src/jaz/gravis/t4Matrix.h>


class TomoCtfHelper
{
	public:
		
		static std::vector<CTF> loadCtffind4(
			std::string path, 
			int imageCount,
			double voltage = 300.0, 
			double Cs = 2.2,
			double Q0 = 0.1, 
			double Bfac = 0.0,
			double scale = 1.0);
		
		static std::vector<CTF> loadCtfplotter(
			std::string path, 
			int imageCount,
			double voltage = 300.0, 
			double Cs = 2.2,
			double Q0 = 0.1, 
			double Bfac = 0.0,
			double scale = 1.0);
				
		static CTF setFromFile(
			std::stringstream& line, 
			double voltage, double Cs, double Q0, double Bfac, double scale);
		
		static void writeToFile(
			int frameIndex, 
			const CTF& ctf,
			std::ofstream& file);
		
		static void read(
			std::string doseFn, std::string ctfFn, 
			int s, int fc, 
			double voltage, double Cs, double pixelSize,
			std::string stackFn, 
			std::vector<CTF>& centreCTFs_out, 
			BufferedImage<float>& doseWeights_out,
			std::vector<double>& cumulativeDose_out);
};

#endif
