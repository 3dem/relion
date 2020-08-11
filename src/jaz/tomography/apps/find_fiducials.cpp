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
#include <src/jaz/tomography/fiducials.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>

using namespace gravis;


int main(int argc, char *argv[])
{
	std::string tomoSetFn, outDir;
	double thresh, binning_out, binning_in, beadRadius_A;
	int tomoIndex, num_threads;
	bool diag;
	
	IOParser parser;

	try
	{	
		parser.setCommandLine(argc, argv);

		int gen_section = parser.addSection("General refinement options");

		outDir = parser.getOption("--o", "Output directory");
		tomoSetFn = parser.getOption("--t", "Tomogram set", "tomograms.star");
		thresh = textToDouble(parser.getOption("--d", "Detection threshold", "5"));
		beadRadius_A = textToDouble(parser.getOption("--r", "Bead radius [Å]", "100"));
		binning_in = textToDouble(parser.getOption("--bin0", "Search binning level", "4"));
		binning_out = textToDouble(parser.getOption("--bin1", "CC binning level", "4"));
		diag = parser.checkOption("--diag", "Write out diagnostic information");
		
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		
		Log::readParams(parser);

		if (parser.checkForErrors())
		{
			parser.writeUsage(std::cout);
			exit(1);
		}
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	outDir = ZIO::makeOutputDir(outDir);
	
	TomogramSet tomogramSet(tomoSetFn);	
	const int tc = tomogramSet.size();
	

	for (int t = 0; t < tc; t++)
	{
		Log::beginSection("Tomogram " + ZIO::itoa(t));
		            
		Tomogram tomogram0 = tomogramSet.loadTomogram(t, true);

		Tomogram tomogram = tomogram0.FourierCrop(binning_out, num_threads);
		const int fc = tomogram.frameCount;

		const int w = tomogram.stack.xdim;
		const int h = tomogram.stack.ydim;

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
			BufferedImage<float> sliceHP = ImageFilter::highpassStackGaussPadded(slice, 2 * beadRadius_px);
			BufferedImage<float> sliceBP = ImageFilter::Gauss2D(sliceHP, 0, beadRadius_px / 2, true);
			BufferedImage<float> CC2D = Similarity::CC_2D(fidKernel, sliceBP);

			fidCC.getSliceRef(f).copyFrom(CC2D);
		}

		const d3Vector origin(0.0);
		const d3Vector spacing(binning_out);
		const d3Vector diagonal = d3Vector(tomogram.w0, tomogram.h0, tomogram.d0) / binning_out;


		std::vector<gravis::d3Vector> detections = Detection::findLocalMaxima(
			tomogram, fidCC, origin, spacing, diagonal,
			(float)thresh, 10000, beadRadius_px, num_threads, binning_in, "debug_");

		Log::print(ZIO::itoa(detections.size()) + " blobs found.");
		

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

			mesh.writePly(outDir+"fiducials_"+tomogram0.name+".ply");
		}

		std::string fidFn = Fiducials::write(
			detections,
			tomogram0.optics.pixelSize,
			tomogram0.name,
			outDir);
		
		tomogramSet.setFiducialsFile(t, fidFn);
		
		Log::endSection();
	}
	
	tomogramSet.write(outDir+"tomograms.star");
}
