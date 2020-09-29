#include "sample_manifold.h"
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/manifold/manifold_set.h>
#include <src/jaz/mesh/mesh_builder.h>

using namespace gravis;


void SampleManifoldProgram::readParameters(int argc, char *argv[])
{
	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		tomogram_set_filename = parser.getOption("--t", "Tomogram set", "tomograms.star");
		manifolds_filename = parser.getOption("--m", "Manifolds set", "manifolds.star");

		spacing = textToDouble(parser.getOption("--spacing", "Particle spacing [A]"));
		depth = textToDouble(parser.getOption("--depth", "Depth below surface [A]"));
		const double max_tilt_deg = textToDouble(parser.getOption("--max_tilt", "Maximum tilt angle [degrees]"));
		max_tilt = DEG2RAD(max_tilt_deg);
		avoid_missing_wedge = parser.checkOption("--nmw", "Do not sample particles from the missing wedges");
		avoid_present_wedge = parser.checkOption("--npw", "Do not sample particles from the present wedges");

		output_path = parser.getOption("--o", "Output filename pattern");

		Log::readParams(parser);

		parser.checkForErrors();

	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
}

void SampleManifoldProgram::run()
{
	output_path = ZIO::makeOutputDir(output_path);

	TomogramSet tomogram_set(tomogram_set_filename);
	ManifoldSet manifold_set(manifolds_filename);

	const int tomogram_count = 1;//tomogram_set.size();

	for (int t = 0; t < tomogram_count; t++)
	{
		Tomogram tomogram = tomogram_set.loadTomogram(t, false);

		const double pixel_size = tomogram.optics.pixelSize;


		std::map<int, const Manifold*> manifolds_map
				= manifold_set.getManifoldsInTomogram(tomogram.name);

		Mesh diagnostic;

		for (std::map<int, const Manifold*>::iterator it = manifolds_map.begin();
			 it != manifolds_map.end(); it++)
		{
			const int manifold_index = it->first;
			const Manifold* manifold = it->second;

			std::vector<RigidAlignment> locations = manifold->sampleParticles(
						spacing / pixel_size,
						depth / pixel_size,
						-max_tilt,  max_tilt,
						!avoid_present_wedge,
						!avoid_missing_wedge,
						tomogram.projectionMatrices);

			for (int i = 0; i < locations.size(); i++)
			{
				RigidAlignment ra = locations[i];

				Matrix2D<double> A;
				Euler_angles2matrix(ra.rot, ra.tilt, ra.psi, A, false);
				const d3Vector n(A(2,0), A(2,1), A(2,2));

				MeshBuilder::addCone(
					pixel_size * ra.position,
					pixel_size * ra.position + spacing * n,
					0.5 * spacing, 20, diagnostic);
			}
		}

		diagnostic.writeObj(output_path + tomogram.name + ".obj");
	}
}
