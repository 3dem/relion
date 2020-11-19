#include "template_picker.h"
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/symmetry.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/fiducials.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/optics/aberrations_cache.h>
#include <src/jaz/math/Euler_angles_relion.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/time.h>
#include <iostream>

using namespace gravis;


void TemplatePickerProgram::readBasicParameters(IOParser& parser, int argc, char *argv[])
{
	optimisation_set.read(
		parser,
		true,           // optimisation set
		false,  false,  // particles
		true,   true,   // tomograms
		false,  false,  // trajectories
		false,  false,  // manifolds
		false,  false); // reference

	template_filename = parser.getOption("--template", "Template file name");
	fiducials_radius_A = textToDouble(parser.getOption("--frad", "Fiducial marker radius [Å]", "100"));
	num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "8"));

	out_dir = parser.getOption("--o", "Output directory");
}

void TemplatePickerProgram::readParameters(int argc, char *argv[])
{
	IOParser parser;

	parser.setCommandLine(argc, argv);

	readBasicParameters(parser, argc, argv);

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	out_dir = ZIO::prepareTomoOutputDirectory(out_dir, argc, argv);
}

void TemplatePickerProgram::initialise()
{
	template_map_RS.read(template_filename);

	TomoReferenceMap::presharpen(template_map_RS, 1.0);

	const double taper_edge_width = 5.0;
	Reconstruction::taper(template_map_RS, taper_edge_width, true, 1);

	FFT::FourierTransform(template_map_RS, template_map_FS, FFT::Both);

	Centering::shiftInSitu(template_map_FS);

	tomogramSet = TomogramSet(optimisation_set.tomograms);

	if (!tomogramSet.globalTable.containsLabel(EMDL_TOMO_FIDUCIALS_STARFILE))
	{
		Log::warn("No fiducial markers present: you are advised to run relion_tomo_find_fiducials first.");
	}
}

void TemplatePickerProgram::run()
{
	initialise();


	Tomogram tomogram = tomogramSet.loadTomogram(0, true);

	const int w = tomogram.stack.xdim;
	const int wh = w / 2 + 1;
	const int h = tomogram.stack.ydim;
	const int fc = tomogram.frameCount;
	const double pixel_size = tomogram.optics.pixelSize;
	const double ba = pixel_size;

	bool has_fiducials = tomogram.hasFiducials();

	Log::print("Filtering");


	const double segmentation_binning = 1;

	/*Tomogram tomogram_binned = tomogram0.FourierCrop(segmentation_binning, num_threads);*/

	std::vector<d3Vector> fiducials(0);


	if (has_fiducials)
	{
		Log::print("Erasing fiducial markers");

		const double fiducials_radius = fiducials_radius_A / pixel_size;

		fiducials = Fiducials::read(tomogram.fiducialsFilename, tomogram.optics.pixelSize);

		Fiducials::erase(
			fiducials,
			fiducials_radius / segmentation_binning,
			tomogram,
			num_threads);
	}

	tomogram.stack.write(out_dir+"stack_no_fid.mrc");

	BufferedImage<fComplex> framesFS(wh,h,fc);

	NewStackHelper::FourierTransformStack(tomogram.stack, framesFS, true);

	const int s_powSpec = 512;

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		BufferedImage<double> powSpec = PowerSpectrum::periodogramAverage2D(
			tomogram.stack, s_powSpec, s_powSpec, 2.0, f, false);

		std::vector<double> powSpec1D = RadialAvg::fftwHalf_2D_lin(powSpec);

		std::vector<float> frqWghts1D(powSpec1D.size());

		for (int i = 0; i < powSpec1D.size(); i++)
		{
			// divide by the spectral power (a^2), so that the resulting CC
			// becomes a product of (REF / a) * (DATA / a)
			frqWghts1D[i] = (float)(1.0 / powSpec1D[i]);
		}

		const double dose = tomogram.cumulativeDose[f];

		for (int yy = 0; yy < h; yy++)
		for (int xx = 0; xx < wh; xx++)
		{
			const double x = xx;
			const double y = yy < h/2? yy : yy - h;
			const double ru = sqrt(x*x/(w*w) + y*y/(h*h));
			const double rd = s_powSpec * ru;

			const int r0 = (int)(rd);
			const int r1 = r0 + 1;

			if (r1 >= powSpec1D.size())
			{
				framesFS(xx,yy,f) *= frqWghts1D[powSpec1D.size()-1];
			}
			else
			{
				const double t = rd - r0;

				framesFS(xx,yy,f) *= (1 - t) * frqWghts1D[r0] + t * frqWghts1D[r1];
			}

			framesFS(xx,yy,f) *= Damage::getWeight(dose, ru / ba);
		}
	}

	BufferedImage<float> framesFilteredRS(w,h,fc);
	NewStackHelper::inverseFourierTransformStack(framesFS, framesFilteredRS, true);
	framesFilteredRS.write(out_dir+"stack_filtered.mrc");

	BufferedImage<float> maskRS(w,h);
	const double maskRad = template_map_RS.xdim / 2.0;
	const double maskFalloff = 20;

	double maskSum = 0.0;

	for (int yy = 0; yy < h; yy++)
	for (int xx = 0; xx < w; xx++)
	{
		const double x = xx < w/2? xx : xx - w;
		const double y = yy < h/2? yy : yy - h;

		const double r = sqrt(x*x + y*y);

		if (r > maskRad)
		{
			maskRS(xx,yy) = 0.f;
		}
		else if (r > maskRad - maskFalloff)
		{
			const double dr = (r - maskRad + maskFalloff) / maskFalloff;
			maskRS(xx,yy) = 0.5f * (cos(PI * dr) + 1.f);
		}
		else
		{
			maskRS(xx,yy) = 1.f;
		}

		maskSum += maskRS(xx,yy);
	}

	maskRS /= maskSum;

	maskRS.write(out_dir + "maskRS.mrc");

	BufferedImage<fComplex> maskFS;

	FFT::FourierTransform(maskRS, maskFS, FFT::Both);

	BufferedImage<float> framesSqRS(w,h,fc);

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		for (int yy = 0; yy < h; yy++)
		for (int xx = 0; xx < w; xx++)
		{
			framesSqRS(xx,yy,f) = framesFilteredRS(xx,yy,f) * framesFilteredRS(xx,yy,f);
		}
	}

	framesSqRS.write(out_dir + "framesSqRS.mrc");

	BufferedImage<fComplex> framesSqFS(wh,h,fc);

	NewStackHelper::FourierTransformStack(framesSqRS, framesSqFS, true);

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		for (int yy = 0; yy < h; yy++)
		for (int xx = 0; xx < wh; xx++)
		{
			framesSqFS(xx,yy,f) = framesSqFS(xx,yy,f) * maskFS(xx,yy).conj();
		}
	}

	NewStackHelper::inverseFourierTransformStack(framesSqFS, framesSqRS, true);

	framesSqRS.write(out_dir + "framesSqRS_filt.mrc");

	pick(0, DEG2RAD(45), 0, tomogram, framesFS, framesSqRS, num_threads);
}

void TemplatePickerProgram::pick(
		double rot, double tilt, double psi,
		const Tomogram& tomogram,
		const BufferedImage<fComplex>& framesFS,
		const BufferedImage<float>& maskedFramesSqRS,
		int num_threads)
{
	const int s = template_map_FS.ydim;
	const int sh = s / 2 + 1;
	const int w = tomogram.stack.xdim;
	const int wh = w / 2 + 1;
	const int h = tomogram.stack.ydim;
	const int fc = tomogram.frameCount;
	const int ba = tomogram.optics.pixelSize * s;

	const d4Matrix P_ref = Euler::anglesToMatrix4(rot, tilt, psi);


	std::vector<BufferedImage<fComplex>> prediction2D_small_FS(
				num_threads, BufferedImage<fComplex>(sh,s));

	std::vector<BufferedImage<float>> prediction2D_small_RS(
				num_threads, BufferedImage<float>(s,s));


	std::vector<BufferedImage<fComplex>> predictions2D_large_FS(
				num_threads, BufferedImage<fComplex>(wh,h));

	std::vector<BufferedImage<float>> predictions2D_large_RS(
				num_threads, BufferedImage<float>(w,h));


	std::vector<BufferedImage<fComplex>> CC_FS(
				num_threads, BufferedImage<fComplex>(wh,h));

	std::vector<BufferedImage<float>> CC_RS_temp(
				num_threads, BufferedImage<float>(w,h));

	BufferedImage<float> CC_RS(w, h, fc);


	for (int th = 0; th < num_threads; th++)
	{
		predictions2D_large_RS[th].fill(0.f);
	}

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		const int th = omp_get_thread_num();

		const d4Matrix P = tomogram.projectionMatrices[f] * P_ref;

		ForwardProjection::forwardProject(template_map_FS, {P}, prediction2D_small_FS[th], 1);

		CTF ctf = tomogram.centralCTFs[f];

		for (int y = 0; y < s; y++)
		for (int x = 0; x < sh; x++)
		{
			const double xa = x / ba;
			const double ya = (y < s/2? y : y - s) / ba;

			prediction2D_small_FS[th](x,y) *= -ctf.getCTF(xa,ya);
		}

		FFT::inverseFourierTransform(prediction2D_small_FS[th], prediction2D_small_RS[th], FFT::Both);

		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			const int xx = x < s/2? x : w + x - s;
			const int yy = y < s/2? y : w + y - s;

			predictions2D_large_RS[th](xx,yy) = prediction2D_small_RS[th](x,y);
		}

		FFT::FourierTransform(predictions2D_large_RS[th], predictions2D_large_FS[th], FFT::Both);


		for (int y = 0; y < h;  y++)
		for (int x = 0; x < wh; x++)
		{
			const double mod = (1 - 2*(x%2)) * (1 - 2*(y%2));

			CC_FS[th](x,y) = mod * framesFS(x,y,f) * predictions2D_large_FS[th](x,y).conj();
		}

		FFT::inverseFourierTransform(CC_FS[th], CC_RS_temp[th], FFT::Both);

		CC_RS.getSliceRef(f).copyFrom(CC_RS_temp[th]);

		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			CC_RS(x,y,f) -= 0.5f * maskedFramesSqRS(x,y,f);
		}
	}

	CC_RS.write(out_dir+"DEBUG_CC_RS_nrm.mrc");


	{
		const double binning = 8;
		const double taper_dist = template_map_RS.xdim / binning;

		const d3Vector origin(0.0);
		const d3Vector spacing(1.0);
		const d3Vector diagonal = d3Vector(tomogram.w0, tomogram.h0, tomogram.d0);

		const int w2D = CC_RS.xdim;
		const int h2D = CC_RS.ydim;

		const int wt2D = tomogram.stack.xdim;
		const int ht2D = tomogram.stack.ydim;
		const int fc = tomogram.stack.zdim;



		if (w2D != wt2D || h2D != ht2D)
		{
			REPORT_ERROR_STR("Detection::findLocalMaxima: images are of incorrect size: "
							 << w2D << "x" << h2D << " vs. " << wt2D << "x" << ht2D);
		}

		const int w2Db = w2D / binning;
		const int h2Db = h2D / binning;

		BufferedImage<float> binnedSimilarity(w2Db, h2Db, fc);

		#pragma omp parallel for num_threads(num_threads)
		for (int f = 0; f < fc; f++)
		{
			BufferedImage<float> binned = Resampling::downsampleMax_2D_full(
						CC_RS.getConstSliceRef(f), w2Db, h2Db);

			binnedSimilarity.getSliceRef(f).copyFrom(binned);
		}

		Tomogram binnedTomogram = tomogram.FourierCrop(binning, num_threads, false);

		binnedSimilarity.write(out_dir+"DEBUG_CC_RS_binned_max.mrc");


		const int w3D = diagonal.x;
		const int h3D = diagonal.y;
		const int d3D = diagonal.z;

		const int w3Db = w3D / binning;
		const int h3Db = h3D / binning;
		const int d3Db = d3D / binning;


		BufferedImage<float> coarseVol(w3Db, h3Db, d3Db), coarseMask(w3Db, h3Db, d3Db);

		RealSpaceBackprojection::backprojectRaw(
			binnedTomogram.projectionMatrices, binnedSimilarity,
			coarseVol, coarseMask,
			origin, spacing * binning, num_threads,
			RealSpaceBackprojection::Linear,
			taper_dist, taper_dist, 3.0);

		coarseVol.write(out_dir+"DEBUG_coarseVol.mrc");
	}
}
