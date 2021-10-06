#include <src/args.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/particle_set.h>


int main(int argc, char *argv[])
{
	try
	{

		std::string outDir;
		OptimisationSet optimisation_set;

		IOParser parser;
		parser.setCommandLine(argc, argv);

		optimisation_set.read(
			parser,
			true,           // optimisation set
			true,   true,   // particles
			true,   true,   // tomograms
			true,   false,  // trajectories
			false,  false,  // manifolds
			true,   true);  // reference

		int gen_section = parser.addSection("General optics splitting options");

		outDir = parser.getOption("--o", "Output directory");

		Log::readParams(parser);

		if (parser.checkForErrors())
		{
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
		}

		outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);

		TomogramSet tomogram_set(optimisation_set.tomograms, true);
		ParticleSet particle_set(optimisation_set.particles, optimisation_set.trajectories, true);

		std::vector<std::vector<ParticleIndex>> particles = particle_set.splitByTomogram(tomogram_set, true);

		if (particle_set.numberOfOpticsGroups() > 1)
		{
			Log::warn("This data set already contains multiple optics groups.");
		}

		MetaDataTable old_optics_table = particle_set.optTable;

		particle_set.optTable.clear();
		particle_set.optTable.setName(old_optics_table.getName());

		const int tc = particles.size();

		for (int t = 0; t < tc; t++)
		{
			particle_set.optTable.addObject();
			particle_set.optTable.setObject(old_optics_table.getObject(0));

			particle_set.optTable.setValue(EMDL_IMAGE_OPTICS_GROUP, t+1, t);
			particle_set.optTable.setValue(EMDL_IMAGE_OPTICS_GROUP_NAME, tomogram_set.getTomogramName(t), t);

			const int pc = particles[t].size();

			for (int p = 0; p < pc; p++)
			{
				particle_set.setOpticsGroup(particles[t][p], t);
			}
		}

		tomogram_set.write(outDir + "tomograms.star");
		optimisation_set.tomograms = outDir+"tomograms.star";

		particle_set.write(outDir + "particles.star");
		optimisation_set.particles = outDir+"particles.star";

		optimisation_set.write(outDir + "optimisation_set.star");

	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}
