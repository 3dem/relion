#include "phaseline_average.h"

using namespace gravis;

std::vector<double> PhaseLineAverage::averageNN(
		const RawImage<float>& data, 
		const RawImage<float>& phase,
		int bins)
{
	std::vector<double> avg(bins, 0.0);
	std::vector<double> cnt(bins, 0.0);
	
	const int w = phase.xdim;
	const int h = phase.ydim;
	const int d = phase.zdim;
	
	if (data.xdim != w || data.ydim != h || data.zdim != d)
	{
		REPORT_ERROR_STR("PhaseLineAverage::averageNN: incompatible image sizes (data: "
						 << data.getSizeString() << ", phase: " << phase.getSizeString() );
	}
	
	f2Vector bounds = findBounds(phase);
	float rng = bounds[1] - bounds[0];
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		float phi0 = phase(x,y,z);
		float idxf = (bins - 1) * (phi0 - bounds[0]) / rng;
		
		int idx = (int)(idxf + 0.5f);
		
		avg[idx] += data(x,y,z);
		cnt[idx] += 1.0;
	}
	
	for (int i = 0; i < bins; i++)
	{
		if (cnt[i] > 0.0)
		{
			avg[i] /= cnt[i];
		}
	}
	
	return avg;
}

BufferedImage<float> PhaseLineAverage::expandLIN(
		const RawImage<float>& data, 
		const RawImage<float>& phase, 
		const std::vector<double>& avg)
{
	const int w = phase.xdim;
	const int h = phase.ydim;
	const int d = phase.zdim;
	const int bins = avg.size();
	
	if (data.xdim != w || data.ydim != h || data.zdim != d)
	{
		REPORT_ERROR_STR("PhaseLineAverage::expand: incompatible image sizes (data: "
						 << data.getSizeString() << ", phase: " << phase.getSizeString() );
	}
	
	BufferedImage<float> out(w,h,d);
	
	f2Vector bounds = findBounds(phase);
	float rng = bounds[1] - bounds[0];
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		float phi0 = phase(x,y,z);
		float idxf = (bins - 1) * (phi0 - bounds[0]) / rng;
		
		int idx0 = (int)idxf;		
		int idx1 = idx0 + 1;
		float fr = idxf - idx0;
		
		if (idx1 >= bins)
		{
			out(x,y,z) = avg[bins-1];
		}
		else
		{
			out(x,y,z) = fr * avg[idx1] + (1-fr) * avg[idx0];
		}
	}
	
	return out;
}

f2Vector PhaseLineAverage::findBounds(
		const RawImage<float>& phase)
{
	const int w = phase.xdim;
	const int h = phase.ydim;
	const int d = phase.zdim;
	
	float minVal = std::numeric_limits<float>::max();
	float maxVal = -std::numeric_limits<float>::max();
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const float phi = phase(x,y,z);
		
		if (phi > maxVal) maxVal = phi;
		if (phi < minVal) minVal = phi;
	}
	
	return f2Vector(minVal, maxVal);
}
