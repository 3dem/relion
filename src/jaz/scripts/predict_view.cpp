
#include <unistd.h>
#include <string.h>
#include <fstream>

#include <src/args.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/new_ft.h>
#include <src/jaz/single_particle/noise_helper.h>
#include <src/jaz/single_particle/fftw_helper.h>
#include <src/jaz/single_particle/reference_map.h>
#include <src/jaz/single_particle/new_reference_map.h>
#include <src/jaz/image/stack_helper.h>
#include <src/jaz/image/filter.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/util/zio.h>

#include <omp.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	std::string starFn, outPath;
	double minFreqPx, lowPass;
	bool oppositeHalf, predictCTF, writeObs;
	int minMG, maxMG, threads;

	//ReferenceMap reference;
	NewReferenceMap reference;

	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);

		parser.addSection("General options");

		starFn = parser.getOption("--i", "Input particle *.star file");
		reference.read(parser, argc, argv);

		minFreqPx = textToDouble(parser.getOption("--min_freq", "Min. image frequency [px]", "0"));

		oppositeHalf = parser.checkOption("--opposite_half", "Correlate with opposite half-set");
		predictCTF = parser.checkOption("--predict_CTF", "Modulate prediction by CTF");
		writeObs = parser.checkOption("--write_obs", "Write out the observed images as well");

		lowPass = textToDouble(parser.getOption("--lowpass", "Blur by a Gaussian with this sigma [px]", "-1"));

		minMG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
		maxMG = textToInteger(parser.getOption("--max_MG", "Last micrograph index (default is to process all)", "-1"));

		threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));
		outPath = parser.getOption("--o", "Output path");

		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		return RELION_EXIT_FAILURE;
	}

	if (outPath[outPath.length()-1] != '/')
	{
		outPath = outPath + "/";
	}

	std::string command = "mkdir -p " + outPath;
	int res = system(command.c_str());


	ObservationModel obsModel;
	MetaDataTable mdt0;

	ObservationModel::loadSafely(starFn, obsModel, mdt0);

	std::vector<MetaDataTable> allMdts = StackHelper::splitByMicrographName(mdt0);

	reference.load(1, false);

	const int s = reference.s;
	const int sh = s/2 + 1;

	if (maxMG < 0) maxMG = allMdts.size() - 1;

	for (int m = 0; m <= maxMG; m++)
	{
		int opticsGroup;
		allMdts[m].getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, 0);
		opticsGroup--;

		std::vector<Image<Complex>> pred = reference.predictAll(
			allMdts[m], obsModel,
		//	oppositeHalf? ReferenceMap::Opposite : ReferenceMap::Own,
			oppositeHalf? NewReferenceMap::Opposite : NewReferenceMap::Own,
			threads, predictCTF, true, false);

		const int pc = pred.size();

		BufferedImage<RFLOAT> output(s,s,pc);
		BufferedImage<Complex> output_FS(sh,s,pc);

		for (int p = 0; p < pc; p++)
		{
			BufferedImage<Complex> slice_FS(pred[p]);

			output_FS.getSliceRef(p).copyFrom(slice_FS);

			BufferedImage<RFLOAT> slice_RS;

			Centering::shiftInSitu(slice_FS);

			FFT::inverseFourierTransform(slice_FS, slice_RS);

			for (int y = 0; y < s; y++)
			for (int x = 0; x < s; x++)
			{
				output(x,y,p) = slice_RS(x,y);
			}
		}

		if (lowPass > 0)
		{
			output = ImageFilter::GaussStack(output, lowPass, true);
		}

		output.write(
			outPath +  "micrograph_" + ZIO::itoa(m) + "_prediction.mrc",
			obsModel.getPixelSize(opticsGroup));

		output_FS.writeVtk(
			outPath +  "micrograph_" + ZIO::itoa(m) + "_prediction_FS.vtk",
			d3Vector(0.0), (1.0 / obsModel.getPixelSize(opticsGroup)) * d3Vector(1));

		if (writeObs)
		{
			std::vector<Image<RFLOAT>> stack = StackHelper::loadStack(&allMdts[m], "", threads);

			std::cout << stack.size() << std::endl;

			BufferedImage<RFLOAT> stackOut(s,s,pc);

			for (int p = 0; p < pc; p++)
			{
				stackOut.getSliceRef(p).copyFrom(stack[p]);
			}

			if (lowPass > 0)
			{
				stackOut = ImageFilter::GaussStack(stackOut, lowPass, true);
			}

			stackOut.write(
				outPath +  "micrograph_" + ZIO::itoa(m) + "_observation.mrc",
				obsModel.getPixelSize(opticsGroup));
		}
	}

	return RELION_EXIT_SUCCESS;
}
