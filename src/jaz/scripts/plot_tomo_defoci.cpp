#include <src/args.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/tomography/tomogram_set.h>


using namespace gravis;


int main(int argc, char *argv[])
{
	std::string inFn, outDir;


	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		inFn = parser.getOption("--t", "Input tomogram set");
		outDir = parser.getOption("--o", "Output directory");


		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	TomogramSet tomogramSet(inFn);

	outDir = ZIO::makeOutputDir(outDir);

	std::ofstream file(outDir + "defoci.dat");

	const int tc = tomogramSet.size();

	for (int dim = 0; dim < 2; dim++)
	{
		double x_coord = 0;

		for (int t = 0; t < tc; t++)
		{
			Tomogram tomogram = tomogramSet.loadTomogram(t, false);

			const int fc = tomogram.frameCount;

			for (int f = 0; f < fc; f++)
			{
				if (dim == 0)
				{
					file << x_coord << ' ' << tomogram.centralCTFs[f].DeltafU << '\n';
				}
				else
				{
					file << x_coord << ' ' << tomogram.centralCTFs[f].DeltafV << '\n';
				}

				x_coord += 1.0 / (double) (fc + 1);
			}

			x_coord += 1.0 / (double) (fc + 1);
		}

		file << '\n';
	}

	return 0;
}

