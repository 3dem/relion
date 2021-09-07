#include "fcc.h"
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/projection/fwd_projection.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/tomogram.h>

#include <src/jaz/image/centering.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/radial_avg.h>
#include <src/jaz/optics/aberrations_cache.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>

#include <src/time.h>

using namespace gravis;


BufferedImage<double> FCC::compute(
	const ParticleSet& dataSet,
	const std::vector<ParticleIndex>& partIndices,
	const Tomogram& tomogram,
	const std::vector<BufferedImage<fComplex>>& referenceFS, 
	bool flip_value, 
	int num_threads)
{
	BufferedImage<double> FCC3 = compute3(
		dataSet, partIndices, tomogram, referenceFS, flip_value, num_threads);

	return divide(FCC3);
}

BufferedImage<double> FCC::compute3(
	const ParticleSet& dataSet,
	const std::vector<ParticleIndex>& partIndices,
	const Tomogram& tomogram,
	const std::vector<BufferedImage<fComplex>>& referenceFS,
	bool flip_value,
	int num_threads)
{
	const int s = referenceFS[0].ydim;
	const int sh = s/2 + 1;
		
	const int pc = partIndices.size();
	const int fc = tomogram.frameCount;
	
	std::vector<BufferedImage<double>> FCC_by_thread(num_threads);
	
	for (int th = 0; th < num_threads; th++)
	{
		FCC_by_thread[th] = BufferedImage<double>(sh, fc, 3);
		FCC_by_thread[th].fill(0.0);
	}

	AberrationsCache aberrationsCache(dataSet.optTable, s, dataSet.getOriginalPixelSize(0));
	BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s, 1.0);

	
	Log::beginProgress("Computing Fourier-cylinder correlations", pc/num_threads);
		
	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{
		const int th = omp_get_thread_num();
		
		if (th == 0)
		{
			Log::updateProgress(p);
		}
		
		const ParticleIndex part_id = partIndices[p];
		
		const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(part_id, fc, tomogram.optics.pixelSize);
		const std::vector<bool> isVisible = tomogram.determineVisiblity(traj, s/2.0);
		
		d4Matrix projCut;
				
		BufferedImage<fComplex> observation(sh,s);
		
		for (int f = 0; f < fc; f++)
		{
			if (!isVisible[f]) continue;
			
			TomoExtraction::extractFrameAt3D_Fourier(
					tomogram.stack, f, s, 1.0, tomogram, traj[f],
					observation, projCut, 1, true);

			RawImage<float> doseSlice = doseWeights.getSliceRef(f);

			BufferedImage<fComplex> prediction = Prediction::predictModulated(
					part_id, dataSet, projCut, s,
					tomogram.getCtf(f, dataSet.getPosition(part_id)),
					tomogram.optics.pixelSize,
					aberrationsCache,
					referenceFS,
					Prediction::OppositeHalf,
					Prediction::AmplitudeAndPhaseModulated,
					&doseSlice,
					Prediction::CtfUnscaled);

			const float scale = flip_value? -1.f : 1.f;
			
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				const fComplex z0 = observation(x,y);
				const fComplex z1 = scale * prediction(x,y);
				
				const double yy = y < s/2? y : y - s;
				const double r = sqrt(x*x + yy*yy);
				
				int ri  = (int) (r + 0.5);
				
				if (ri < sh)
				{
					FCC_by_thread[th](ri,f,0) += (double)(z0.real * z1.real + z0.imag * z1.imag);
					FCC_by_thread[th](ri,f,1) += (double)(z0.real * z0.real + z0.imag * z0.imag);
					FCC_by_thread[th](ri,f,2) += (double)(z1.real * z1.real + z1.imag * z1.imag);
				}
			}
		}
	}
	
	Log::endProgress();
	
	BufferedImage<double> FCC(sh, fc, 3);
	FCC.fill(0.0);
	
	for (int th = 0; th < num_threads; th++)
	{
		FCC += FCC_by_thread[th];
	}

	return FCC;	
}

BufferedImage<double> FCC::divide(const BufferedImage<double>& fcc3)
{
	const int sh = fcc3.xdim;
	const int fc = fcc3.ydim;
	
	BufferedImage<double> out(sh,fc);
	
	for (int f = 0; f < fc; f++)
	for (int x = 0; x < sh; x++)
	{
		const double n = fcc3(x,f,0);
		const double d1 = fcc3(x,f,1);
		const double d2 = fcc3(x,f,2);
		
		const double d = sqrt(d1 * d2);
		
		out(x,f) = d > 0.0? n / d : 0.0;
	}
	
	return out;
}

BufferedImage<double> FCC::sumOverTime(const RawImage<double>& fcc3)
{
	const int sh = fcc3.xdim;
	const int fc = fcc3.ydim;

	BufferedImage<double> sum(sh,1,3);

	sum.fill(0.0);

	for (int x = 0; x < sh; x++)
	{
		for (int f = 0; f < fc; f++)
		{
			for (int d = 0; d < 3; d++)
			{
				sum(x,0,d) += fcc3(x,f,d);
			}
		}
	}

	return sum;
}
