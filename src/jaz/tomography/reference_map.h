#ifndef TOMO_REFERENCE_MAP_H
#define TOMO_REFERENCE_MAP_H

#include <src/jaz/image/buffered_image.h>
#include <vector>

class OptimisationSet;

class TomoReferenceMap
{
	public:

		static void presharpen(BufferedImage<float>& map_RS, double padding);
		

		TomoReferenceMap();

		void read(IOParser& parser);
		void read(const OptimisationSet& optimisationSet);
		
			std::string mapFilenames[2], maskFilename, fscFilename;
			double fscThresholdWidth, freqCutoff_A, pixelSize;
			int lastShell;
			bool flatWeight;

			std::vector<BufferedImage<float>> image_real;
			BufferedImage<float> mask;
			std::vector<BufferedImage<fComplex>> image_FS;
			std::vector<double> SNR_weight;


		void load(int boxSize = -1, int verbosity = 1);

		template <typename T>
		void contributeWeight(RawImage<T> freqWeights) const;
			
		int getBoxSize() const;

};

template <typename T>
void TomoReferenceMap::contributeWeight(
	RawImage<T> freqWeights) const
{
	const int s0 = image_FS[0].ydim;
	const int sh0 = image_FS[0].xdim;

	const int s1 = freqWeights.ydim;
	const int sh1 = freqWeights.xdim;

	const double scale = s1 / (double) s0;

	for (int y = 0; y < s1; y++)
	for (int x = 0; x < sh1; x++)
	{
		double xx = x;
		double yy = y < s1/2? y : y - s1;

		double r = sqrt(xx*xx + yy*yy) / scale;

		int r0 = (int) r;
		int r1 = r0 + 1;

		if (r1 >= sh0 || r1 > lastShell + fscThresholdWidth * scale)
		{
			freqWeights(x,y) *= 0.0;
		}
		else
		{
			const double f = r - r0;

			const double g0 = SNR_weight[r0];
			const double g1 = SNR_weight[r1];

			const double g = (1 - f) * g0 + f * g1;

			double env = 1.0;

			if (r > lastShell &&
				r < lastShell + fscThresholdWidth * scale)
			{
				const double t = (r - lastShell) / (fscThresholdWidth * scale);

				env = (cos(PI * t) + 1.0) / 2.0;
			}

			freqWeights(x,y) *= env * g;
		}
	}
}

#endif
