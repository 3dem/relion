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

	parser.setCommandLine(argc, argv);

	optimisationSet.read(
		parser,
		true,            // optimisation set
		false,  false,   // particles
		true,   true,    // tomograms
		false,  false,   // trajectories
		true,   true,    // manifolds
		false,  false);  // reference

	int gen_section = parser.addSection("General options");

	spacing = textToDouble(parser.getOption("--spacing", "Particle spacing [A]"));
	depth = textToDouble(parser.getOption("--depth", "Depth below surface [A]"));
	const double max_tilt_deg = textToDouble(parser.getOption("--max_tilt", "Maximum tilt angle [degrees]"));
	max_tilt = DEG2RAD(max_tilt_deg);
	avoid_missing_wedge = parser.checkOption("--nmw", "Do not sample particles from the missing wedges");
	avoid_present_wedge = parser.checkOption("--npw", "Do not sample particles from the present wedges");
	store_tilt_series = parser.checkOption("--ts", "Store the name of the tilt series in the star file");

	output_path = parser.getOption("--o", "Output directory");

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	output_path = ZIO::prepareTomoOutputDirectory(output_path, argc, argv);
}

void SampleManifoldProgram::run()
{
	TomogramSet tomogram_set(optimisationSet.tomograms);
	ManifoldSet manifold_set(optimisationSet.manifolds);

	MetaDataTable optics_table, particles_table;

	optics_table.setName("optics");
	particles_table.setName("particles");

	Tomogram tomogram0 = tomogram_set.loadTomogram(0, false);

	const std::string optics_group_name = "optics_group_1";

	optics_table.addObject();
	optics_table.setValue(EMDL_IMAGE_OPTICS_GROUP, 1);
	optics_table.setValue(EMDL_IMAGE_OPTICS_GROUP_NAME, optics_group_name);
	optics_table.setValue(EMDL_CTF_CS, tomogram0.optics.Cs);
	optics_table.setValue(EMDL_CTF_VOLTAGE, tomogram0.optics.voltage);
	optics_table.setValue(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, tomogram0.optics.pixelSize);


	const int tomogram_count = tomogram_set.size();

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

				Matrix2D<RFLOAT> A;
				Euler_angles2matrix(ra.rot, ra.tilt, ra.psi, A, false);
				const d3Vector n(A(0,2), A(1,2), A(2,2));


				particles_table.addObject();
				const int j = particles_table.numberOfObjects() - 1;


				particles_table.setValue(EMDL_TOMO_PARTICLE_NAME,  tomogram.name + "/" + ZIO::itoa(j), j);

				particles_table.setValue(EMDL_TOMO_NAME, tomogram.name, j);
				particles_table.setValue(EMDL_TOMO_MANIFOLD_INDEX, manifold_index, j);

				if (store_tilt_series)
				{
					particles_table.setValue(EMDL_TOMO_TILT_SERIES_NAME, tomogram.tiltSeriesFilename, j);
				}

				particles_table.setValue(EMDL_IMAGE_COORD_X, ra.position.x, j);
				particles_table.setValue(EMDL_IMAGE_COORD_Y, ra.position.y, j);
				particles_table.setValue(EMDL_IMAGE_COORD_Z, ra.position.z, j);

				particles_table.setValue(EMDL_TOMO_SUBTOMOGRAM_ROT, ra.rot, j);
				particles_table.setValue(EMDL_TOMO_SUBTOMOGRAM_TILT, ra.tilt, j);
				particles_table.setValue(EMDL_TOMO_SUBTOMOGRAM_PSI, ra.psi, j);

				particles_table.setValue(EMDL_ORIENT_TILT_PRIOR, 0.0, j);

				particles_table.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, 0.0, j);
				particles_table.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, 0.0, j);
				particles_table.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, 0.0, j);

				particles_table.setValue(EMDL_ORIENT_ORIGIN_X_PRIOR_ANGSTROM, 0.0, j);
				particles_table.setValue(EMDL_ORIENT_ORIGIN_Y_PRIOR_ANGSTROM, 0.0, j);
				particles_table.setValue(EMDL_ORIENT_ORIGIN_Z_PRIOR_ANGSTROM, 0.0, j);

				particles_table.setValue(EMDL_ORIENT_ROT, 0.0, j);
				particles_table.setValue(EMDL_ORIENT_TILT, 0.0, j);
				particles_table.setValue(EMDL_ORIENT_PSI, 0.0, j);

				particles_table.setValue(EMDL_MLMODEL_GROUP_NO, t, j);
				particles_table.setValue(EMDL_PARTICLE_CLASS, 1, j);
				particles_table.setValue(EMDL_IMAGE_OPTICS_GROUP, 1, j);
				particles_table.setValue(EMDL_PARTICLE_RANDOM_SUBSET, j % 2 + 1, j);



				MeshBuilder::addCone(
					pixel_size * ra.position,
					pixel_size * ra.position + spacing * n,
					0.5 * spacing, 20, diagnostic);
			}
		}

		diagnostic.writeObj(output_path + tomogram.name + ".obj");
	}


	std::ofstream ofs(output_path + "particles.star");

	optics_table.write(ofs);
	particles_table.write(ofs);

	optimisationSet.particles = output_path + "particles.star";
	optimisationSet.write(output_path + "optimisation_set.star");
}
