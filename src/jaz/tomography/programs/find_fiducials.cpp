#include "find_fiducials.h"
#include <src/args.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/projection/real_backprojection.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/image/filter.h>
#include <src/jaz/image/detection.h>
#include <src/jaz/image/similarity.h>
#include <src/jaz/mesh/mesh.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/tomography/fiducials.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>

using namespace gravis;


FindFiducialsProgram::FindFiducialsProgram(int argc, char *argv[])
{
	IOParser parser;

	parser.setCommandLine(argc, argv);

	optimisationSet.read(
		parser,
		true,            // optimisation set
		false,  false,   // particles
		true,   true,    // tomograms
		false,  false,   // trajectories
		false,  false,   // manifolds
		false,  false);  // reference

	int gen_section = parser.addSection("General options");

	outDir = parser.getOption("--o", "Output directory");
	thresh = textToDouble(parser.getOption("--d", "Detection threshold", "5"));
	beadRadius_A = textToDouble(parser.getOption("--r", "Bead radius [Å]", "50"));
	rel_spacing = textToDouble(parser.getOption("--sp", "Minimal bead spacing (as a multiple of radius, L^inf distance)", "1.5"));
	binning_in = textToDouble(parser.getOption("--bin0", "Search binning level", "4"));
	binning_out = textToDouble(parser.getOption("--bin1", "CC binning level", "4"));
	diag = parser.checkOption("--diag", "Write out diagnostic information");
	max_MG = textToInteger(parser.getOption("--max_MG", "Last tilt series to consider", "-1"));

	num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);
}

void FindFiducialsProgram::run()
{
	if (diag)
	{
		ZIO::makeDir(outDir+"Diagnostic");
	}

	TomogramSet tomogramSet(optimisationSet.tomograms);
	int tc = tomogramSet.size();
	
	if (max_MG >= 0 && max_MG < tc)
	{
		tc = max_MG + 1;
	}

	BufferedImage<float> visualisation(0,0,0);

	double vis_pixel_size = 1;

	for (int t = 0; t < tc; t++)
	{
		Log::beginSection("Tomogram " + ZIO::itoa(t));

		Tomogram tomogram0 = tomogramSet.loadTomogram(t, true);

		Tomogram tomogram = tomogram0.FourierCrop(binning_out, num_threads);
		const int fc = tomogram.frameCount;

		const int w = tomogram.stack.xdim;
		const int h = tomogram.stack.ydim;

		vis_pixel_size = tomogram.optics.pixelSize;

		const double beadRadius_px = beadRadius_A / tomogram.optics.pixelSize;

		BufferedImage<float> fidKernel = Detection::smallCircleKernel<float>(
					beadRadius_px, w, h);

		BufferedImage<tComplex<float>> fidKernelFS;

		FFT::FourierTransform(fidKernel, fidKernelFS);

		BufferedImage<float> fidCC(w, h, fc);

		#pragma omp parallel for num_threads(num_threads)
		for (int f = 0; f < fc; f++)
		{
			BufferedImage<float> slice = tomogram.stack.getSliceRef(f);
			BufferedImage<float> sliceHP = ImageFilter::highpassStackGaussPadded(slice, beadRadius_px/4);
			BufferedImage<float> sliceBP = ImageFilter::Gauss2D(sliceHP, 0, beadRadius_px/2, true);
			BufferedImage<float> CC2D = Similarity::CC_2D(fidKernel, sliceBP);

			fidCC.getSliceRef(f).copyFrom(CC2D);
		}

		const d3Vector origin(0.0);
		const d3Vector spacing(binning_out);
		const d3Vector diagonal = d3Vector(tomogram.w0, tomogram.h0, tomogram.d0) / binning_out;


		std::vector<gravis::d3Vector> detections = Detection::findLocalMaxima(
			tomogram, fidCC, origin, spacing, diagonal,
			(float)thresh, 10000, rel_spacing * beadRadius_px,
			num_threads, binning_in,
			diag? outDir+"Diagnostic/" + tomogram.name + "_" : "");

		Log::print(ZIO::itoa(detections.size()) + " beads found.");

		if (diag)
		{
			Mesh mesh;

			for (int i = 0; i < detections.size(); i++)
			{
				MeshBuilder::addOctahedron(
					detections[i] * tomogram0.optics.pixelSize,
					beadRadius_A / tomogram0.optics.pixelSize,
					mesh);
			}

			mesh.writePly(outDir+"Diagnostic/fiducials_"+tomogram0.name+".ply");
		}

		std::string fidFn = Fiducials::write(
			detections,
			tomogram0.optics.pixelSize,
			tomogram0.name,
			outDir);
		
		tomogramSet.setFiducialsFile(t, fidFn);


		if (visualisation.xdim == 0)
		{
			visualisation.resize(w,h,tc);
		}

		const int best_frame = tomogram0.getLeastDoseFrame();

		visualisation.getSliceRef(t).copyFrom(tomogram.stack.getSliceRef(best_frame));

		const float mu = Normalization::computeMean(visualisation.getSliceRef(t));
		const float var = Normalization::computeVariance(visualisation.getSliceRef(t), mu);
		const float val = mu + 6 * sqrt(var);


		for (int i = 0; i < detections.size(); i++)
		{
			d3Vector pw = detections[i];
			d4Vector pi = tomogram.projectionMatrices[best_frame] * d4Vector(pw);

			const int d = 21;
			const int q = 1;
			const int x0 = std::round(pi.x);
			const int y0 = std::round(pi.y);

			for (int i = 0; i < d; i++)
			for (int j = -q; j <= q; j++)
			{
				const int xx = x0 + i - d/2;
				const int yy = y0 + i - d/2;
				const int x1 = x0 + j;
				const int y1 = y0 + j;

				if (xx >= 0 && xx < w && y1 >= 0 && y1 < h)
				{
					visualisation(xx,y1,t) = val;
				}

				if (x1 >= 0 && x1 < w && yy >= 0 && yy < h)
				{
					visualisation(x1,yy,t) = val;
				}
			}
		}
		
		Log::endSection();
	}
	
	tomogramSet.write(outDir + "tomograms.star");

	optimisationSet.tomograms = outDir + "tomograms.star";
	optimisationSet.write(outDir + "optimisation_set.star");

	visualisation.write(outDir+"detections.mrc", vis_pixel_size);
}
