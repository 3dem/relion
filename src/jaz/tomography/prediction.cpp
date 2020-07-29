#include "prediction.h"
#include "particle_set.h"
#include "projection/fwd_projection.h"
#include "reconstruction.h"
#include "tomo_ctf_helper.h"
#include "tomolist.h"
#include "tomogram.h"
#include "reference_map.h"

#include <src/jaz/image/centering.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/radial_avg.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>

#include <src/time.h>


#define TEST_CC_SCALE 0
#define TEST_SPECTRAL_POWER 0

using namespace gravis;


BufferedImage<fComplex> Prediction::predictModulated(
		int particle_id, const ParticleSet* dataSet, d4Matrix proj, int s,
		const CTF& ctf, double pixelSize,
		const std::vector<BufferedImage<fComplex>>& referenceFS,
		HalfSet halfSet,
		Modulation modulation)
{
	BufferedImage<fComplex> prediction = predictFS(
				particle_id, dataSet, proj, s, referenceFS, halfSet);
	
	if (  modulation == AmplitudeModulated
	   || modulation == AmplitudeAndPhaseModulated)
	{
		const int sh = s / 2 + 1;
		BufferedImage<float> ctfImg(sh, s);
		
		ctf.draw(s, s, pixelSize, &ctfImg[0]);
		
		prediction *= ctfImg;
	}
	
	if (  modulation == PhaseModulated
	   || modulation == AmplitudeAndPhaseModulated)
	{
		// @TODO: modulate phase
	}
	
	return prediction;
}

BufferedImage<fComplex> Prediction::predictFS(
		int particle_id, const ParticleSet* dataSet, d4Matrix proj, int s,
		const std::vector<BufferedImage<fComplex>>& referenceFS,
		HalfSet halfSet)
{
	const int sh = s/2 + 1;
			
	BufferedImage<float> predictionReal = predictRS(particle_id, dataSet, proj, s, referenceFS, halfSet);
	
	BufferedImage<fComplex> prediction(sh,s);
	FFT::FourierTransform(predictionReal, prediction, FFT::Both);
	
	return prediction;
}

BufferedImage<float> Prediction::predictRS(
		int particle_id, const ParticleSet* dataSet, d4Matrix proj, int s,
		const std::vector<BufferedImage<fComplex>>& referenceFS,
		HalfSet halfSet)
{
	const int sh = s/2 + 1;
	
	const d4Matrix particleToTomo = dataSet->getMatrix4x4(particle_id, s, s, s);
	const d4Matrix projPart = proj * particleToTomo;	
	
	const int hs0 = dataSet->getHalfSet(particle_id);
	const int hs = (halfSet == OppositeHalf)? 1 - hs0: hs0;
	
	BufferedImage<fComplex> prediction(sh,s), psf(sh,s);

	ForwardProjection::forwardProject(
			referenceFS[hs], {projPart}, prediction, psf, 1);
	
	BufferedImage<float> predictionReal(s,s);
	
	Reconstruction::correctStack(prediction, psf, predictionReal, false, 1);
	
	return predictionReal;
}

std::vector<BufferedImage<double> > Prediction::computeCroppedCCs(
		const ParticleSet* dataSet,
		const std::vector<int>& partIndices,
		const Tomogram& tomogram,
		const TomoReferenceMap& referenceMap,
		const BufferedImage<float>& frqWghts,
		const std::vector<int>& sequence,
		int maxRange,				
		bool flip_value,
		int num_threads,
		double paddingFactor,
		HalfSet halfSet)
{
	const int s = referenceMap.getBoxSize();
	const int sh = s/2 + 1;
	
	const int s_act = s * paddingFactor;
	const int sh_act = s_act / 2 + 1;
	
	const int pc = partIndices.size();
	const int fc = tomogram.frameCount;
	const int diam = (int)(2 * maxRange * paddingFactor);
	const int border = (int)(s * paddingFactor - diam) / 2;
	
	std::vector<BufferedImage<double>> CCs(pc);
	
	for (int p = 0; p < pc; p++)
	{
		CCs[p] = BufferedImage<double>(diam, diam, fc);
	}
	
	
	#if TEST_SPECTRAL_POWER
	
		const int data_str = sh + 2048;
		const int buf_size = data_str * num_threads;
		
		std::vector<double> 
				sumPowObs(buf_size,0.0), 
				sumPowPred(buf_size,0.0), 
				sumPowCCunw(buf_size,0.0), 
				sumPowCC(buf_size,0.0), 
				sumWgh(buf_size,0.0);
		
		Image<float> 
				predSumFS(sh,s,num_threads),
				obsSumFS(sh,s,num_threads),
				ccSumFS(sh,s,num_threads);
		
		predSumFS.fill(0.f);
		obsSumFS.fill(0.f);
		ccSumFS.fill(0.f);
		
	#endif	
		
		
	Log::beginProgress("Computing cross correlations", pc/num_threads);
	
	#pragma omp parallel for num_threads(num_threads)		
	for (int p = 0; p < pc; p++)
	{
		const int th = omp_get_thread_num();
		
		if (th == 0)
		{
			Log::updateProgress(p);
		}
		
		const int part_id = partIndices[p];	
		
		const std::vector<d3Vector> traj = dataSet->getTrajectoryInPixels(part_id, fc, tomogram.optics.pixelSize);
		
		d4Matrix projCut;	
		
		BufferedImage<fComplex> observation(sh,s);
		
		for (int ft = 0; ft < fc; ft++)
		{
			const int f = sequence[ft];
			
			TomoExtraction::extractFrameAt3D_Fourier(
					tomogram.stack, f, s, 1.0, tomogram.projectionMatrices[f], traj[f],
					observation, projCut, 1, false, true);
						
			BufferedImage<fComplex> prediction = Prediction::predictModulated(
					part_id, dataSet, projCut, s, 
					tomogram.getCtf(f, dataSet->getPosition(part_id)),
					tomogram.optics.pixelSize,
					referenceMap.image_FS, halfSet);
					
			BufferedImage<fComplex> ccFS(sh,s);
			
			const float scale = flip_value? -1.f : 1.f;
			
			observation(0,0) = fComplex(0.f, 0.f);
			prediction(0,0) = fComplex(0.f, 0.f);
			
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				ccFS(x,y) = scale * frqWghts(x,y,f) * observation(x,y) * prediction(x,y).conj();
				
				#if TEST_SPECTRAL_POWER
				
				if (ft == 0)
				{
					const double yy = y < s/2? y : y - s;
					const int ri = (int)sqrt(x*x + yy*yy);
					
					const double nobs = observation(x,y).norm();
					const double npred = prediction(x,y).norm();
					const double ncc = ccFS(x,y).norm() / (scale*scale);
					
					obsSumFS(x,y,th) += nobs;
					predSumFS(x,y,th) += npred;
					ccSumFS(x,y,th) += ncc;
					
					if (ri < sh)
					{	
						const int idx = th * data_str + ri;
						
						sumPowObs[idx] += nobs;
						sumPowPred[idx] += npred;
						sumPowCCunw[idx] += npred*nobs;
						sumPowCC[idx] += ncc;
						sumWgh[idx] += 1.0;
					}
				}
					
				#endif
			}
			
			BufferedImage<fComplex> ccFS_padded = Padding::padCorner2D_half(ccFS, sh_act, s_act);
			
			BufferedImage<float> ccRS_full;
			
			/* Note: we don't use "FFT::Both" here, because the values of the CC 
			   are already normalised in Fourier space due to whitening:		*/
			FFT::inverseFourierTransform(ccFS_padded, ccRS_full);
			
			CCs[p].copySliceFrom(ft, 
				Centering::fftwFullToHumanFull(
					Padding::unpadCorner2D_full(ccRS_full, border)));
			
			
			#if TEST_CC_SCALE
			
			if (th == 0)
			{
				{
					Image<fComplex> prediction2 = prediction;
							
					for (int y = 0; y < s;  y++)
					for (int x = 0; x < sh; x++)
					{
						prediction2(x,y) = scale * frqWghts(x,y,f) * prediction(x,y);
					}
					
					Image<float> predRS, obsRS;
					
					FFT::inverseFourierTransform(observation, obsRS, FFT::Both);
					FFT::inverseFourierTransform(prediction2, predRS, FFT::Both);
					
					double CC_real = 0.0;
					
					for (int y = 0; y < s; y++)
					for (int x = 0; x < s; x++)
					{
						CC_real += obsRS(x,y) * predRS(x,y);
					}
											
					std::cout 
						<< "real, cplx, ratio:   " << CC_real 
						<< "  " << ccRS_full(0,0)
						<< "  " << (CC_real / ccRS_full(0,0)) << '\n';
					
					obsRS.write("debug/obsRS_"+ZIO::itoa(p)+"_"+ZIO::itoa(f)+".mrc");
					predRS.write("debug/predRS_"+ZIO::itoa(p)+"_"+ZIO::itoa(f)+".mrc");
				}
				
				{				
					dComplex CC_manual(0.0, 0.0);
					
					for (int yy = 0; yy < s;  yy++)
					for (int xx = 0; xx < s; xx++)
					{
						if (xx < sh)
						{
							int x = xx;
							int y = yy;
									
							dComplex zd(ccFS(x,y).real, ccFS(x,y).imag);
							CC_manual += zd;
						}
						else
						{
							int x = s - xx;
							int y = (s - yy) % s;
							
							dComplex zd(ccFS(x,y).real, -ccFS(x,y).imag);
							CC_manual += zd;
						}
					}
					
					
					std::cout 
						<< "real, cplx, ratio:   " << CC_manual 
						<< "  " << ccRS_full(0,0)
						<< "  " << (CC_manual.real / ccRS_full(0,0)) << '\n';
				}
			}
			
			#endif
			
		}
	}
	
	Log::endProgress();
	
	#if TEST_SPECTRAL_POWER
			
		std::ofstream 
				ofsObs("debug/spec_obs.dat"),
				ofsPred("debug/spec_pred.dat"),
				ofsCCunw("debug/spec_CCunw.dat"),
				ofsCC("debug/spec_CC.dat");
							 
		for (int i = 0; i < sh; i++)
		{
			for (int th = 1; th < num_threads; th++)
			{
				sumPowObs[i] += sumPowObs[i + th * data_str];
				sumPowPred[i] += sumPowPred[i + th * data_str];
				sumPowCCunw[i] += sumPowCCunw[i + th * data_str];
				sumPowCC[i] += sumPowCC[i + th * data_str];
				sumWgh[i] += sumWgh[i + th * data_str];
			}
			
			const double wg = sumWgh[i];
			
			if (wg > 0.0)
			{
				sumPowObs[i] /= wg;
				sumPowPred[i] /= wg;
				sumPowCCunw[i] /= wg;
				sumPowCC[i] /= wg;
			}
			
			const double toA = 1.0 / (s * tomogram.optics.pixelSize);
			
			ofsObs << toA * i << ' ' << sumPowObs[i] << '\n';
			ofsPred << toA * i << ' ' << sumPowPred[i] << '\n';
			ofsCCunw << toA * i << ' ' << sumPowCCunw[i] << '\n';
			ofsCC << toA * i << ' ' << sumPowCC[i] << '\n';
		}
		
		for (int th = 1; th < num_threads; th++)
		{
			obsSumFS.getSliceRef(0) += obsSumFS.getSliceRef(th);
			predSumFS.getSliceRef(0) += predSumFS.getSliceRef(th);
			ccSumFS.getSliceRef(0) += ccSumFS.getSliceRef(th);
		}
		
		obsSumFS.getSliceRef(0).write("debug/obsSumFS.mrc");
		predSumFS.getSliceRef(0).write("debug/predSumFS.mrc");
		ccSumFS.getSliceRef(0).write("debug/ccSumFS.mrc");
		
	#endif	
	
	return CCs;
}
