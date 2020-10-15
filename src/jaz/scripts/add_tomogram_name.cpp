#include <src/args.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/tomography/tomogram_set.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	IOParser parser;

	std::string tableFn, inputFn, outputFn;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		inputFn = parser.getOption("--i", "Input particle set");
		tableFn = parser.getOption("--m", "Table mapping micrograph names (first column) to tomogram names (second column)");
		outputFn = parser.getOption("--o", "Output particle set");

		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	std::map<std::string,std::string> mic2tomo;

	std::ifstream ifs(tableFn);

	if (!ifs)
	{
		REPORT_ERROR_STR("Unable to read "+tableFn);
	}

	std::string line;

	while (std::getline(ifs, line))
	{
		std::stringstream sts(line);
		std::string mic, tomo;

		sts >> mic;
		sts >> tomo;

		std::cout << mic << " -> " << tomo << std::endl;

		mic2tomo[mic] = tomo;
	}

	MetaDataTable opticsTable, particleTable;

	opticsTable.read(inputFn, "optics");
	particleTable.read(inputFn, "particles");

	particleTable.addLabel(EMDL_TOMO_NAME);
	particleTable.addLabel(EMDL_TOMO_TILT_SERIES_NAME);

	for (int p = 0; p < particleTable.numberOfObjects(); p++)
	{
		std::string mic = particleTable.getString(EMDL_MICROGRAPH_NAME, p);

		if (mic2tomo.find(mic) == mic2tomo.end())
		{
			REPORT_ERROR_STR("Micrograph name '" << mic << "' not found in "
							 << tableFn);
		}

		std::string tomo = mic2tomo[mic];

		particleTable.setValue(EMDL_TOMO_NAME, tomo, p);
		particleTable.setValue(EMDL_TOMO_TILT_SERIES_NAME, mic, p);
	}

	std::ofstream ofs(outputFn);

	opticsTable.write(ofs);
	particleTable.write(ofs);

	return 0;
}
