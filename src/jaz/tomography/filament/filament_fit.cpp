#include "filament_fit.h"
#include "filament_model.h"
#include "filament.h"
#include "coordinate_conventions.h"
#include <src/jaz/util/zio.h>
#include <omp.h>


using namespace gravis;


FilamentFit::FilamentFit(
		const Filament& filament, 
		const FilamentMapping& mapping,
		const std::vector<gravis::d4Matrix>& proj,
		const RawImage<float> binnedStack, 
		const RawImage<float> outlier_mask, 
		const FilamentModel* model,
		double binning, 
		int num_threads)

: filament(filament),
  mapping(mapping),
  proj(proj),
  binnedStack(binnedStack),
  outlier_mask(outlier_mask),
  model(model),
  maxDistBinned(filament.maxRadius / binning),
  binning(binning),
  num_threads(num_threads)
{
	BufferedImage<float> optimisationMask = mapping.mask * outlier_mask;
	
	encode(optimisationMask, optimisationPixels);
	encode(mapping.mask, allPixels);
}

double FilamentFit::f(
		const std::vector<double> &x, 
		void *tempStorage) const
{
	const int fc = mapping.index.zdim;
	const int diam = (int)(2 * maxDistBinned);
	
	std::vector<std::vector<float>> phases(fc);	
	
	model->apply(x, optimisationPixels.splineCoords, phases, num_threads);
	
	const int data_padding = 512;
	std::vector<double> sums(data_padding * num_threads, 0.0);
	
	// DEBUG!!!
	
	/*#pragma omp parallel for num_threads(num_threads)		
	for (int f = 0; f < fc; f++)*/
	for (int f = 20; f < 21; f++)
	{
		if (f != 20) continue;
		
		int th = omp_get_thread_num();
		
		std::vector<float> average(diam);
		
		const size_t pc = phases[f].size();
		
		std::vector<float> prediction(pc);
		
		averageValues(phases[f], optimisationPixels.pixelValue[f], average);
		expandValues(phases[f], average, prediction);
				
		double sum_f(0.0);
		
		for (size_t p = 0; p < pc; p++)
		{
			double d = prediction[p] - optimisationPixels.pixelValue[f][p];
			
			sum_f += d * d;
		}
		
		sums[data_padding * th] = sum_f;
	}
	
	double sum(0.0);
	
	for (int th = 0; th < num_threads; th++)
	{
		sum += sums[data_padding * th];
	}
	
	return sum / optimisationPixels.totalNumber;
}

void FilamentFit::grad(
		const std::vector<double> &x, 
		std::vector<double> &gradDest, 
		void *tempStorage) const
{
	const int fc = mapping.index.zdim;
	const int diam = (int)(2 * maxDistBinned);
	
	std::vector<std::vector<float>> phases(fc), slopes(fc);
	
	model->apply(x, optimisationPixels.splineCoords, phases, num_threads);
	
	#pragma omp parallel for num_threads(num_threads)		
	for (int f = 0; f < fc; f++)
	{
		const size_t pc = phases[f].size();
		
		std::vector<float> prediction(pc);
		
		slopes[f].resize(pc);
		std::vector<float> average(diam);
		
		averageValues(phases[f], optimisationPixels.pixelValue[f], average);
		expandValues(phases[f], average, prediction);		
		getSlopes(phases[f], average, slopes[f]);
		
		
		// DEBUG!!!
		if (f == 20)
		{		
			for (size_t p = 0; p < pc; p++)
			{
				// E = (v(x) - v0)^2
				// => dE/dx = 2 (v(x) - v0) dv(x)/dx
				slopes[f][p] = 2.0 * (prediction[p] - optimisationPixels.pixelValue[f][p]) * slopes[f][p];
			}
		}
		else
		{
			for (size_t p = 0; p < pc; p++)
			{
				slopes[f][p] = 0.0;
			}
		}
	}
	
	model->computeGradient(x, optimisationPixels.splineCoords, slopes, gradDest, num_threads);
	
	for (int i = 0; i < gradDest.size(); i++)
	{
		gradDest[i] /= optimisationPixels.totalNumber;
	}
}

void FilamentFit::report(int iteration, double cost, const std::vector<double> &x) const
{
	std::cout << "         " << iteration << ": " << cost << std::endl;
}

BufferedImage<float> FilamentFit::visualise(
		const std::vector<double> &x, 
		bool subtract,
		BufferedImage<float>* original,
		float scale)
{
	const int fc = mapping.index.zdim;
	const int diam = (int)(2 * maxDistBinned);
	
	std::vector<std::vector<float>> estPhases(fc), allPhases(fc);
	std::vector<float> average(diam);
	
	model->apply(x, optimisationPixels.splineCoords, estPhases, num_threads);
	model->apply(x, allPixels.splineCoords, allPhases, num_threads);
		
	BufferedImage<float> out = original == 0? binnedStack : *original;
	
	const bool constantCorrection = false;

	for (int f = 0; f < fc; f++)
	{
		const size_t pc_all = allPhases[f].size();
		
		std::vector<float> prediction(pc_all);
		
		averageValues(estPhases[f], optimisationPixels.pixelValue[f], average);
		
		double zero_level;
		
		if (constantCorrection)
		{
			zero_level = computeZeroLevel(average);
		}
		else
		{	
			correctZeroLevel(average);
		}	
		
		expandValues(allPhases[f], average, prediction);
		
		if (subtract)
		{
			if (constantCorrection)
			{
				for (size_t p = 0; p < pc_all; p++)
				{
					out[(long long int)allPixels.pixelIndex[f][p]] -= scale * (prediction[p] - zero_level);
				}
			}
			else
			{
				for (size_t p = 0; p < pc_all; p++)
				{
					out[(long long int)allPixels.pixelIndex[f][p]] -= scale * prediction[p];
				}
			}
		}
		else
		{
			for (size_t p = 0; p < pc_all; p++)
			{
				out[(size_t)allPixels.pixelIndex[f][p]] = prediction[p];
			}
		}
	}
	
	return out;
}

BufferedImage<float> FilamentFit::computeCostByOffset(
		const std::vector<double>& x, 
		double minOffset, double maxOffset, int samples, double sigma)
{
	const int fc = mapping.index.zdim;
	const int diam = (int)(2 * maxDistBinned);
	
	std::vector<std::vector<float>> phases(fc);	
	
	model->apply(x, optimisationPixels.splineCoords, phases, num_threads);
		
	BufferedImage<float> out((int)(filament.arcLen / binning), 2 * samples, fc);
	out.fill(0.f);
			
	// DEBUG!!!
	
	/*#pragma omp parallel for num_threads(num_threads)		
	for (int f = 0; f < fc; f++)*/
	for (int f = 20; f < 21; f++)
	{
		if (f != 20) continue;
				
		std::vector<float> average(diam);
		
		averageValues(phases[f], optimisationPixels.pixelValue[f], average);
		
		BufferedImage<float> slice = costByOffset(
			phases[f], 
			average, 
			optimisationPixels.splineCoords[f],
			optimisationPixels.pixelValue[f],
			minOffset, maxOffset, samples,
			filament.arcLen,
			1.0 / binning,
			sigma);
		
		out.copySliceFrom(f, slice, 0);
	}
	
	return out;
}

BufferedImage<float> FilamentFit::costByOffset(
	const std::vector<float>& phases,       // pc
	const std::vector<float>& phaseValues,  // phc
	const std::vector<float>& arcCoords,    // pc
	const std::vector<float>& pixelValues,  // pc
	double minOffset,
	double maxOffset,
	int samples,
	double arcLength,
	double arcScale, 
	double sigma) const
{
	const int w = (int) (arcScale * arcLength);
	const int h = 2 * samples;
	const size_t pc = phases.size();
	
	const double y_scale = (maxOffset - minOffset) / (double) (samples - 1);
	const double x_scale = arcScale * arcLength / (double) w;
		
	std::vector<float> 
			offsetPhases(pc), 
			offsetValues(pc);
	
	BufferedImage<float> out(w,h), count(w,h);
	out.fill(0.f);
	count.fill(0.f);
	
	for (int off = 0; off < samples; off++)
	{
		const double offset = minOffset + off * y_scale;
		
		for (size_t p = 0; p < pc; p++)
		{
			const double q = arcCoords[COORD_DIM * p + COORD_Q];
			offsetPhases[p] = phases[p] + ((q > 0)? offset : -offset) * binning;
		}
		
		expandValues(offsetPhases, phaseValues, offsetValues);
				
		for (size_t p = 0; p < pc; p++)
		{
			const double d = offsetValues[p] - pixelValues[p];
			const double q = arcCoords[COORD_DIM * p + COORD_Q];
			
			int x = (int) (arcScale * arcCoords[COORD_DIM * p + COORD_T]);
			
			if (x < 0) x = 0;
			if (x >= w) x = w-1;
			
			const int y = (q > 0)? off : off + samples;
			
			out(x,y)   += d*d;
			count(x,y) += 1.f;
		}
	}
	
	if (sigma > 0.0)
	{
		BufferedImage<float> out_smooth(w,h), count_smooth(w,h);
		
		const int r = (int)(3.0 * sigma);
		const double beta = 1.0 / (2.0 * sigma * sigma);
		
		std::vector<double> kernel(2*r + 1);
		
		for (int i = 0; i < 2*r + 1; i++)
		{
			const double x = i - r;			
			kernel[i] = exp( -beta*x*x );
		}
		
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			double sumU(0.0), sumC(0.0), weight(0.0);
			
			for (int xx = x - r; xx <= x + r; xx++)
			{
				if (xx < 0) continue;
				if (xx >= w) break;
				
				const double k = kernel[xx - x + r];
				
				sumU += k * out(xx,y);
				sumC += k * count(xx,y);
				weight += k;
			}
			
			out_smooth(x,y) = sumU / weight;
			count_smooth(x,y) = sumC / weight;
		}
		
		out = out_smooth;
		count = count_smooth;
	}
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		if (count(x,y) > 0.f) 
		{
			out(x,y) /= count(x,y);
		}
		else
		{
			out(x,y) = 0.f;
		}
	}
	
	for (int x = 0; x < w; x++)
	for (int d = 0; d < 2; d++)
	{
		double minVal = std::numeric_limits<double>::max();
		
		for (int y = d * samples; y < (d + 1) * samples; y++)
		{
			if (out(x,y) < minVal) 
			{
				minVal = out(x,y);
			}
		}
		
		for (int y = d * samples; y < (d + 1) * samples; y++)
		{
			out(x,y) -= minVal;
		}
	}
	
	return out;
}
	
void FilamentFit::encode(const BufferedImage<float>& mask, PixelSet& out)
{
	const int w  = mapping.index.xdim;
	const int h  = mapping.index.ydim;
	const int fc = mapping.index.zdim;
	
	out.splineCoords.resize(fc);
	out.pixelValue.resize(fc);
	out.pixelIndex.resize(fc);	
	
	out.totalNumber = 0;
	
	for (long long int f = 0; f < fc; f++)
	{
		long long int count = 0;
		
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			if (mask(x,y,f) > 0.5f)
			{
				count++;
			}
		}
		
		out.totalNumber += count;
		
		out.splineCoords[f].resize(COORD_DIM * count);
		out.pixelValue[f].resize(count);
		out.pixelIndex[f].resize(count);
		
		long long int index = 0;
		
		Spline<gravis::d2Vector,double> spline2D = Filament::projectSpline(filament.spline, proj[f]);
		
		for (long long int y = 0; y < h; y++)
		for (long long int x = 0; x < w; x++)
		{
			if (mask(x,y,f) > 0.5f)
			{
				out.splineCoords[f][COORD_DIM * index + COORD_T] = mapping.index(x,y,f);
				out.splineCoords[f][COORD_DIM * index + COORD_Q] = mapping.signed_dist(x,y,f);
				
				d2Vector d = (d2Vector(x,y) - spline2D.getValue(mapping.index(x,y,f))).normalize();
						
				out.splineCoords[f][COORD_DIM * index + COORD_X] = (float) d.x;
				out.splineCoords[f][COORD_DIM * index + COORD_Y] = (float) d.y;
				
				out.pixelValue[f][index] = binnedStack(x,y,f);
				
				out.pixelIndex[f][index] = f*w*h + y*w + x;
						
				index++;
			}
		}
	}
}

void FilamentFit::averageValues(
		const std::vector<float>& phases,
		const std::vector<float>& pixelValues, 
		std::vector<float>& out) const
{
	const int pc = phases.size();
	const int cc = out.size();
	
	const double mid = cc / 2.0;
	
	for (size_t c = 0; c < cc; c++)
	{
		out[c] = 0.f;
	}
	
	std::vector<float> count(cc, 0.f);
	
	for (size_t p = 0; p < pc; p++)
	{
		const float phi = phases[p] / binning;		
		int ci = (int) std::round(phi + mid);
		
		if (ci < 0) ci = 0;
		if (ci >= cc) ci = cc-1;
		
		out[ci] += pixelValues[p];
		count[ci] += 1.f;
	}
	
	for (size_t c = 0; c < cc; c++)
	{
		if (count[c] > 0.f) 
		{
			out[c] /= count[c];
		}
	}
}

void FilamentFit::expandValues(
		const std::vector<float>& phases, 
		const std::vector<float>& phaseValues, 
		std::vector<float>& out) const
{
	const int pc = phases.size();
	const int os = phaseValues.size();
	
	const double mid = os / 2.0;
	
	for (size_t p = 0; p < pc; p++)
	{
		const float phi = phases[p] / binning;
		
		int p0 = (int) std::floor(phi + mid);
		int p1 = (int) std::ceil(phi + mid);
		
		if (p0 < 0) p0 = 0;
		if (p0 >= os) p0 = os - 1;
		
		if (p1 < 0) p1 = 0;
		if (p1 >= os) p1 = os - 1;
		
		float dphi = phi + mid - p0;
		
		if (dphi < 0) dphi = 0; 
		else if (dphi > 1) dphi = 1;
		
		out[p] = (1.f - dphi) * phaseValues[p0] + dphi * phaseValues[p1];
	}
}

void FilamentFit::getSlopes(
		const std::vector<float>& phases, 
		const std::vector<float>& phaseValues, 
		std::vector<float>& slopes_out) const
{
	const int pc = phases.size();
	const int os = phaseValues.size();
	
	slopes_out.resize(pc);
	
	const double mid = os / 2.0;
	
	for (size_t p = 0; p < pc; p++)
	{
		const float phi = phases[p] / binning;
		
		int p0 = (int) std::floor(phi + mid);
		int p1 = (int) std::ceil(phi + mid);
		
		if (p0 < 0) p0 = 0;
		if (p0 >= os) p0 = os - 1;
		
		if (p1 < 0) p1 = 0;
		if (p1 >= os) p1 = os - 1;
		
		if (p0 > 0 && p1 < os - 1) 
		{
			slopes_out[p] = phaseValues[p1] - phaseValues[p0];
		}
		else 
		{
			slopes_out[p] = 0.f; 
		}
	}	
}

double FilamentFit::computeZeroLevel(const std::vector<float> &phaseValues) const
{
	const int d = phaseValues.size();
	const int d0 = 5 * d / 100;
	const int d1 = 95 * d / 100;
	
	double sum(0.0), cnt(0.0);
	
	for (int i = 0; i < d; i++)
	{
		if (i < d0 || i > d1)
		{
			sum += phaseValues[i];
			cnt += 1.0;
		}
	}
	
	return sum / cnt;
}

void FilamentFit::correctZeroLevel(std::vector<float> &phaseValues) const
{
	const int d = phaseValues.size();
	const int d0 = 5 * d / 100;
	const int d1 = 95 * d / 100;
	
	double sum0(0.0), cnt0(0.0), sum1(0.0), cnt1(0.0);
	
	for (int i = 0; i < d; i++)
	{
		if (i < d0)
		{
			sum0 += phaseValues[i];
			cnt0 += 1.0;
		}
		else if (i > d1)
		{
			sum1 += phaseValues[i];
			cnt1 += 1.0;
		}
	}
	
	double v0 = sum0 / cnt0;
	double v1 = sum1 / cnt1;
	
	for (int i = 0; i < d; i++)
	{
		phaseValues[i] -= v0 + (v1 - v0) * i / (double) (d-1);
		
		/*if (i < d0)
		{
			phaseValues[i] -= v0;
		}
		else if (i < d1)
		{
			phaseValues[i] -= v0 + (v1 - v0) * (i - d0) / (double) (d1 - d0);
		}
		else
		{
			phaseValues[i] -= v1;
		}*/
	}
}


