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
#include <src/jaz/optics/damage.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>

#include <src/time.h>


#define TEST_CC_SCALE 0
#define TEST_SPECTRAL_POWER 0

using namespace gravis;


BufferedImage<fComplex> Prediction::predictModulated(
		ParticleIndex particle_id, const ParticleSet& dataSet, d4Matrix proj, int s,
		const CTF& ctf, double pixelSize,
		const AberrationsCache& aberrationsCache,
		const std::vector<BufferedImage<fComplex>>& referenceFS,
		HalfSet halfSet,
		Modulation modulation,
		const RawImage<float>* doseWeight,
		CtfScale ctfScale,
		const int* xRanges)
{
	BufferedImage<fComplex> prediction = predictFS(
				particle_id, dataSet, proj, s, referenceFS, halfSet, xRanges);

	const int og = dataSet.getOpticsGroup(particle_id);
	
	if (  modulation == AmplitudeModulated
	   || modulation == AmplitudeAndPhaseModulated)
	{
		const int sh = s / 2 + 1;
		BufferedImage<float> ctfImg(sh, s);

		const BufferedImage<double>* gammaOffset =
			aberrationsCache.hasSymmetrical? &aberrationsCache.symmetrical[og] : 0;

		if (ctfScale == CtfUnscaled)
		{
			CTF ctf0 = ctf;
			ctf0.scale = 1.0;
			ctf0.draw(s, s, pixelSize, gammaOffset, &ctfImg[0]);
		}
		else
		{
			ctf.draw(s, s, pixelSize, gammaOffset, &ctfImg[0]);
		}
		
		prediction *= ctfImg;
	}
	
	if (aberrationsCache.hasAntisymmetrical &&
		 (modulation == PhaseModulated || modulation == AmplitudeAndPhaseModulated) )
	{
		if (aberrationsCache.phaseShift[og].ydim != s)
		{
			REPORT_ERROR_STR(
				"Prediction::predictModulated: wrong cached phase-shift size. Box size: "
				<< s << ", cache size: " << aberrationsCache.phaseShift[og].ydim);
		}

		prediction *= aberrationsCache.phaseShift[og];
	}

	if (doseWeight != 0)
	{
		prediction *= (*doseWeight);
	}
	
	return prediction;
}

BufferedImage<fComplex> Prediction::predictFS(
		ParticleIndex particle_id, const ParticleSet& dataSet, d4Matrix proj, int s,
		const std::vector<BufferedImage<fComplex>>& referenceFS,
		HalfSet halfSet,
		const int* xRanges)
{
	const int sh = s/2 + 1;

	const d4Matrix particleToTomo = dataSet.getMatrix4x4(particle_id, s, s, s);
	const d4Matrix projPart = proj * particleToTomo;

	const int hs0 = dataSet.getHalfSet(particle_id);
	const int hs = (halfSet == OppositeHalf)? 1 - hs0: hs0;

	BufferedImage<fComplex> prediction(sh,s);

	if (xRanges != 0)
	{
		ForwardProjection::forwardProjectWithinRange(
			xRanges,
			referenceFS[hs], {projPart}, prediction, 1);
	}
	else
	{
		ForwardProjection::forwardProject(
			referenceFS[hs], {projPart}, prediction, 1);
	}

	return prediction;
}

std::vector<BufferedImage<double> > Prediction::computeCroppedCCs(
		const ParticleSet& dataSet,
		const std::vector<ParticleIndex>& partIndices,
		const Tomogram& tomogram,
		const AberrationsCache& aberrationsCache,
		const TomoReferenceMap& referenceMap,
		const BufferedImage<float>& freqWeights,
		const BufferedImage<float>& doseWeights,
		const std::vector<int>& sequence,
		const BufferedImage<int>& xRanges,
		int maxRange,
		bool flip_value,
		int num_threads,
		double paddingFactor,
		HalfSet halfSet,
		bool verbose)
{
	const int s = referenceMap.getBoxSize();
	const int sh = s/2 + 1;
	
	const int s_act = s * paddingFactor;
	const int sh_act = s_act / 2 + 1;

	/* pad the CCs by 3 pixels of zeros on each side, so that under cubic interpolation,
	   the gradient vanishes at the boundary (1 pixel away from the edge): */
	const bool pad_by_3 = true;
	
	const int pc = partIndices.size();
	const int fc = tomogram.frameCount;
	const int diam = (int)(2 * maxRange * paddingFactor) + (pad_by_3? 6 : 0);
	const int border = (int)(s * paddingFactor - diam) / 2;
	
	std::vector<BufferedImage<double>> CCs(pc);
	
	for (int p = 0; p < pc; p++)
	{
		CCs[p] = BufferedImage<double>(diam, diam, fc);
	}

	if (verbose)
	{
		Log::beginProgress("Computing cross correlations", pc/num_threads);
	}
	
	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{
		const int th = omp_get_thread_num();
		
		if (verbose && th == 0)
		{
			Log::updateProgress(p);
		}
		
		const ParticleIndex part_id = partIndices[p];
		
		const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(
					part_id, fc, tomogram.optics.pixelSize);
		
		d4Matrix projCut;
		
		BufferedImage<fComplex> observation(sh,s);


		for (int ft = 0; ft < fc; ft++)
		{
			const int f = sequence[ft];

			const RawImage<float> doseSlice = doseWeights.getConstSliceRef(f);

			if (!tomogram.isVisible(traj[f], f, s/2.0))
			{
				CCs[p].getSliceRef(f).fill(0.0);

				continue;
			}
			
			TomoExtraction::extractFrameAt3D_Fourier(
					tomogram.stack, f, s, 1.0, tomogram,
					traj[f], observation, projCut, 1, true);

			BufferedImage<fComplex> prediction = Prediction::predictModulated(
					part_id, dataSet, projCut, s, 
					tomogram.getCtf(f, dataSet.getPosition(part_id)),
					tomogram.optics.pixelSize,
					aberrationsCache,
					referenceMap.image_FS, halfSet,
					AmplitudeAndPhaseModulated,
					&doseSlice,
					CtfScaled,
					&xRanges(0,f));
					
			BufferedImage<fComplex> ccFS(sh,s);
			
			const float scale = flip_value? -1.f : 1.f;
			
			observation(0,0) = fComplex(0.f, 0.f);
			prediction(0,0) = fComplex(0.f, 0.f);
			
			for (int y = 0; y < s;  y++)
			{
				for (int x = 0; x < xRanges(y,f); x++)
				{
					ccFS(x,y) = scale * freqWeights(x,y,f) * observation(x,y) * prediction(x,y).conj();
				}
				for (int x = xRanges(y,f); x < sh; x++)
				{
					ccFS(x,y) = fComplex(0.f, 0.f);
				}
			}
			
			BufferedImage<fComplex> ccFS_padded = Padding::padCorner2D_half(ccFS, sh_act, s_act);
			
			BufferedImage<float> ccRS_full;
			
			/* Note: we don't use "FFT::Both" here, because the values of the CC 
			   are already normalised in Fourier space due to whitening:		*/
			FFT::inverseFourierTransform(ccFS_padded, ccRS_full);
			
			CCs[p].copySliceFrom(ft, 
				Centering::fftwFullToHumanFull(
					Padding::unpadCorner2D_full(ccRS_full, border)));

			if (pad_by_3)
			{
				for (int y = 0; y < 3; y++)
				{
					for (int x = 0; x < diam; x++)
					{
						CCs[p](x,y,ft) = 0.0;
					}
				}

				for (int y = 3; y < diam - 3; y++)
				{
					for (int x = 0; x < 3; x++)
					{
						CCs[p](x,y,ft) = 0.0;
					}

					for (int x = diam - 3; x < diam; x++)
					{
						CCs[p](x,y,ft) = 0.0;
					}
				}

				for (int y = diam - 3; y < diam; y++)
				{
					for (int x = 0; x < diam; x++)
					{
						CCs[p](x,y,ft) = 0.0;
					}
				}
			}

			#if TEST_CC_SCALE
			
			if (th == 0)
			{
				{
					Image<fComplex> prediction2 = prediction;
							
					for (int y = 0; y < s;  y++)
					for (int x = 0; x < sh; x++)
					{
						prediction2(x,y) = scale * freqWeights(x,y,f) * prediction(x,y);
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
	
	if (verbose)
	{
		Log::endProgress();
	}
	
	return CCs;
}

void Prediction::predictMicrograph(
		int frame_index,
		const ParticleSet &dataSet,
		const std::vector<ParticleIndex> &partIndices,
		const Tomogram &tomogram,
		const AberrationsCache &aberrationsCache,
		const TomoReferenceMap &referenceMap,
		const RawImage<float>* doseWeights,
		RawImage<float> &target_slice,
		HalfSet halfSet,
		Modulation modulation,
		CtfScale ctfScale)
{
	const int s = referenceMap.getBoxSize();

	const int pc = partIndices.size();
	const int fc = tomogram.frameCount;

	const int f = frame_index;
	const int w = target_slice.xdim;
	const int h = target_slice.ydim;

	BufferedImage<float> prediction_RS(s,s);

	target_slice.fill(0.f);


	for (int p = 0; p < pc; p++)
	{
		const ParticleIndex part_id = partIndices[p];

		const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(
					part_id, fc, tomogram.optics.pixelSize);

		BufferedImage<fComplex> prediction_FS = Prediction::predictModulated(
				part_id, dataSet, tomogram.projectionMatrices[f], s,
				tomogram.getCtf(f, dataSet.getPosition(part_id)),
				tomogram.optics.pixelSize,
				aberrationsCache,
				referenceMap.image_FS,
				halfSet,
				modulation,
				doseWeights,
				ctfScale);

		const d2Vector q = tomogram.projectPoint(traj[f], f);
		const i2Vector c(round(q.x) - s/2, round(q.y) - s/2);
		const d2Vector d(q.x - s/2 - c.x, q.y - s/2 - c.y);

		NewStackHelper::shiftStack(prediction_FS, {-d}, prediction_FS, true, 1);

		FFT::inverseFourierTransform(prediction_FS, prediction_RS, FFT::Both);

		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			int xx = x + c.x;
			int yy = y + c.y;

			if (xx >= 0 && xx < w && yy >= 0 && yy < h)
			{
				target_slice(xx,yy) -= prediction_RS(x,y);
			}
		}
	}
}
