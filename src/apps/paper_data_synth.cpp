#include <vector>
#include <src/strings.h>
#include <src/filename.h>
#include <src/jaz/single_particle/distribution_helper.h>
#include <src/jaz/gravis/t2Vector.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	const int w = 200;
	const int n = 100;
	const double mu = 25;
	const double sigR = 5;
	const double sigN = 2.0;
	
	std::vector<double> ref(w);
	std::vector<std::vector<double>> arrays(n);
	
	std::string out = "illustration/";
	mktree(out);	
	
	std::ofstream refOut(out+"ref.dat");
			
	for (int i = 0; i < w; i++)
	{
		const double x = (i-mu)/sigR;
		ref[i] = exp(-x*x);
		
		refOut << i << " " << ref[i] << "\n";
	}
	
	std::vector<double> avg = std::vector<double>(w, 0.0);
	
	std::ofstream maxOut(out+"maxima.dat");
	
	d2Vector maxAvg(0.0, 0.0);
		
	for (int k = 0; k < n; k++)
	{
		std::ofstream refOut(out+"noise_"+integerToString(k, 3, '0')+".dat");
		
		arrays[k] = std::vector<double>(w);
		
		int imax = 0;
		double vmax = -100.0;
				
		for (int i = 0; i < w; i++)
		{
			arrays[k][i] = ref[i] + DistributionHelper::sampleGauss(0, sigN);
			
			refOut << i << " " << arrays[k][i] << "\n";
			
			if (vmax < arrays[k][i])
			{
				vmax = arrays[k][i];
				imax = i;
			}
			
			avg[i] += arrays[k][i] / (double)n;
		}
		
		refOut << "\n" << imax << " " << vmax << "\n";
		
		maxOut << imax << " " << vmax << "\n\n";
		
		maxAvg.x += imax / (double)n;
		maxAvg.y += vmax / n;
	}	
	
	std::ofstream avgOut(out+"avg.dat");	
	
	int imaxA = 0;
	double vmaxA = -100.0;
	
	for (int i = 0; i < w; i++)
	{		
		avgOut << i << " " << avg[i] << "\n";
		
		if (vmaxA < avg[i])
		{
			vmaxA = avg[i];
			imaxA = i;
		}
	}
	
	avgOut << "\n" << imaxA << " " << vmaxA << "\n";
	avgOut << "\n" << maxAvg.x << " " << maxAvg.y << "\n";
	
	return RELION_EXIT_SUCCESS;	
}
