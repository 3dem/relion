#include "shift_alignment.h"
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/color_helper.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>
#include <src/CPlot2D.h>

using namespace gravis;


std::vector<gravis::d2Vector> ShiftAlignment::alignGlobally(
		const Tomogram& tomogram,
		const std::vector<ParticleIndex>& particles,
		const ParticleSet& particleSet,
		const TomoReferenceMap& referenceMap,
		const RawImage<float>& doseWeights,
		const AberrationsCache& aberrationsCache,
		bool do_whiten,
		int num_threads,
		bool diag,
		const std::string& tag,
		const std::string& outDir)
{
	const int w0  = tomogram.stack.xdim;
	const int wh0 = w0 / 2 + 1;
	const int h0  = tomogram.stack.ydim;
	const int fc = tomogram.frameCount;

	BufferedImage<float> allCCs;

	if (diag)
	{
		allCCs = BufferedImage<float>(w0,h0,fc);
	}

	std::vector<gravis::d2Vector> out(fc);

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		const RawImage<float> obs_slice = tomogram.stack.getConstSliceRef(f);
		const RawImage<float> dose_slice = doseWeights.getConstSliceRef(f);

		BufferedImage<float> pred_slice(w0,h0);

		Prediction::predictMicrograph(
			f, particleSet, particles, tomogram,
			aberrationsCache, referenceMap, &dose_slice,
			pred_slice,
			Prediction::OwnHalf,
			Prediction::AmplitudeAndPhaseModulated,
			Prediction::CtfUnscaled);


		BufferedImage<fComplex> obs_slice_hat, pred_slice_hat;

		FFT::FourierTransform(obs_slice, obs_slice_hat);
		FFT::FourierTransform(pred_slice, pred_slice_hat);

		BufferedImage<fComplex> CC_hat(wh0, h0);

		if (do_whiten)
		{
			std::vector<double> pow_spec = PowerSpectrum::fromFftwHalf(obs_slice_hat);

			for (int y = 0; y < h0;  y++)
			for (int x = 0; x < wh0; x++)
			{
				const double sp = PowerSpectrum::interpolate(x, y, w0, h0, pow_spec);

				CC_hat(x,y) = obs_slice_hat(x,y) * pred_slice_hat(x,y).conj() / sp;
			}
		}
		else
		{
			for (int y = 0; y < h0;  y++)
			for (int x = 0; x < wh0; x++)
			{
				CC_hat(x,y) = obs_slice_hat(x,y) * pred_slice_hat(x,y).conj();
			}
		}

		referenceMap.contributeWeight<fComplex>(CC_hat);

		BufferedImage<float> CC;
		FFT::inverseFourierTransform(CC_hat, CC);

		CC = Centering::fftwFullToHumanFull(CC);

		if (diag)
		{
			allCCs.getSliceRef(f).copyFrom(CC);
		}

		d2Vector opt = Interpolation::quadraticMaxXY(CC) - d2Vector(w0/2, h0/2);

		out[f] = opt;
	}

	if (diag)
	{
		allCCs.write(outDir + "allCCs_" + tag + ".mrc");
	}

	return out;
}

std::vector<d2Vector> ShiftAlignment::alignPerParticle(
		const Tomogram& tomogram,
		const std::vector<BufferedImage<double>>& CCs,
		double padding,
		int range,
		int verbosity,
		bool diag,
		bool per_tomogram_progress,
		const std::string& tag,
		const std::string& outDir)
{
	const int diam = CCs[0].xdim;
	const int pc = CCs.size();

	if (pc == 0)
	{
		return std::vector<d2Vector>(0);
	}

	const int fc = CCs[0].zdim;

	BufferedImage<float> CCsum(diam, diam, fc);

	CCsum.fill(0.f);

	if (verbosity > 0 && per_tomogram_progress)
	{
		Log::beginProgress("Adding up cross-correlations", pc);
	}

	for (int p = 0; p < pc; p++)
	{
		if (verbosity > 0 && per_tomogram_progress)
		{
			Log::updateProgress(p);
		}

		CCsum += CCs[p];
	}

	if (verbosity > 0 && per_tomogram_progress)
	{
		Log::endProgress();
	}

	if (diag)
	{
		CCsum.write(outDir + "CCsum_" + tag + ".mrc");
	}

	d2Vector origin(padding*range + 3, padding*range + 3);

	std::vector<gravis::d2Vector> out(fc);

	for (int f = 0; f < fc; f++)
	{
		d2Vector opt = (Interpolation::quadraticMaxXY(CCsum.getSliceRef(f)) - origin) / padding;
		const int ff = tomogram.frameSequence[f];
		out[ff] = opt;
	}

	return out;
}

void ShiftAlignment::visualiseShifts(
	const std::vector<d2Vector> &shifts,
	const std::vector<int> &sequence,
	const std::string &tomo_name,
	const std::string &file_name_root)
{
	const int fc = shifts.size();

	CPlot2D plot2D(tomo_name + ": frame shifts");
	plot2D.SetXAxisSize(600);
	plot2D.SetYAxisSize(600);
	plot2D.SetDrawLegend(false);
	plot2D.SetFlipY(true);
	plot2D.SetDrawXAxisGridLines(false);
	plot2D.SetDrawYAxisGridLines(false);

	std::vector<CDataSet> points_by_frame(fc);

	for (int ft = 0; ft < fc; ft++)
	{
		gravis::dRGB c = ColorHelper::signedToRedBlue(ft/(double)fc);

		const int f = sequence[ft];

		points_by_frame[f].SetDrawMarker(true);
		points_by_frame[f].SetDrawLine(false);
		points_by_frame[f].SetDatasetColor(c.r,c.g,c.b);
		points_by_frame[f].SetMarkerSize(3);

		const d2Vector d = shifts[f];

		points_by_frame[f].AddDataPoint(CDataPoint(d.x,d.y));
	}

	for (int ft = fc-1; ft >= 0; ft--)
	{
		const int f = sequence[ft];
		plot2D.AddDataSet(points_by_frame[f]);
	}

	FileName fn_eps = file_name_root + ".eps";

	plot2D.OutputPostScriptPlot(fn_eps);
}
