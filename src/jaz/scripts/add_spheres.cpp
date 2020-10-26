#include <src/args.h>
#include <src/jaz/tomography/manifold/CMM_loader.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/manifold/manifold_set.h>
#include <src/jaz/tomography/manifold/sphere.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	IOParser parser;

	std::string listFn, outPath, tomoSetFn;
	double spheres_binning;
	bool diag;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		tomoSetFn = parser.getOption("--t", "Tomogram set filename", "tomograms.star");
		listFn = parser.getOption("--i", "File containing a list of tomogram-name/spheres-file pairs");

		spheres_binning = textToDouble(parser.getOption("--sbin", "Binning factor of the sphere coordinates"));

		diag = parser.checkOption("--diag", "Write out diagnostic information");

		outPath = parser.getOption("--o", "Output filename pattern");

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

	outPath = ZIO::makeOutputDir(outPath);

	if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag");
	}

	ManifoldSet manifold_set;

	std::ifstream list(listFn);

	if (!list)
	{
		REPORT_ERROR_STR("Unable to read "+listFn);
	}

	std::map<std::string, std::string> tomoToSpheres;

	std::string line;

	while (std::getline(list, line))
	{
		std::stringstream sts;
		sts << line;

		std::string tomogram_name, spheresFn;
		sts >> tomogram_name;
		sts >> spheresFn;

		tomoToSpheres[tomogram_name] = spheresFn;
	}




	for (std::map<std::string, std::string>::iterator it = tomoToSpheres.begin();
		 it != tomoToSpheres.end(); it++)
	{
		const std::string tomogram_name = it->first;
		const std::string spheresFn = it->second;

		std::cout << tomogram_name << ": " << spheresFn << std::endl;

		Log::beginSection("Tomogram " + tomogram_name);

		std::vector<gravis::d4Vector> spheres = CMM_Loader::readSpheres(spheresFn, spheres_binning);

		TomogramManifoldSet tomogram_manifold_set;

		for (int m = 0; m < spheres.size(); m++)
		{
			std::vector<double> coefficients{
				spheres[m].x,
				spheres[m].y,
				spheres[m].z,
				spheres[m].w};

			tomogram_manifold_set.add(new Sphere(coefficients, m));
		}

		manifold_set.add(tomogram_name, tomogram_manifold_set);

		Log::endSection();
	}


	manifold_set.write(outPath + "manifolds.star");


	return 0;
}
