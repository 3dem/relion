#include <src/args.h>
#include <src/jaz/tomography/imod_import.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>
#include <src/ctf.h>

using namespace gravis;

struct ParticleData
{
	d3Vector coordinates, angles;
	std::string tomoName;
	bool hasAngles;
};

void addParticles(
	MetaDataTable& table,
	const std::string& knownTomoName,
	std::vector<ParticleData>& outputData,
	bool& hasAngles)
{
	hasAngles =
			table.containsLabel(EMDL_ORIENT_ROT)
		 && table.containsLabel(EMDL_ORIENT_TILT)
		 && table.containsLabel(EMDL_ORIENT_PSI);

	outputData.resize(table.numberOfObjects());

	for (int i = 0; i < table.numberOfObjects(); i++)
	{
		ParticleData& d = outputData[i];

		d.coordinates.x = table.getDouble(EMDL_IMAGE_COORD_X, i);
		d.coordinates.y = table.getDouble(EMDL_IMAGE_COORD_Y, i);
		d.coordinates.z = table.getDouble(EMDL_IMAGE_COORD_Z, i);

		if (hasAngles)
		{
			d.angles[0] = table.getDouble(EMDL_ORIENT_ROT, i);
			d.angles[1] = table.getDouble(EMDL_ORIENT_TILT, i);
			d.angles[2] = table.getDouble(EMDL_ORIENT_PSI, i);
		}
		else
		{
			d.angles[0] = 0.0;
			d.angles[1] = 0.0;
			d.angles[2] = 0.0;
		}

		if (knownTomoName == "")
		{
			d.tomoName = table.getString(EMDL_TOMO_NAME, i);
		}
		else
		{
			d.tomoName = knownTomoName;
		}
	}
}

int main(int argc, char *argv[])
{
	try
	{
		IOParser parser;

		std::string inStarFn, tomoFn, outDir;
		bool flipZ;

		try
		{
			parser.setCommandLine(argc, argv);

			int gen_section = parser.addSection("General options");

			inStarFn = parser.getOption("--i", "Input STAR file containing either a set of particle coordinates or a set of names of files containing those");
			tomoFn = parser.getOption("--t", "Input tomogram set");
			outDir = parser.getOption("--o", "Output directory");
			flipZ = parser.checkOption("--flipZ", "Flip the Z coordinate");

			Log::readParams(parser);

			if (parser.checkForErrors())
			{
				exit(RELION_EXIT_FAILURE);
			}
		}
		catch (RelionError XE)
		{
			parser.writeUsage(std::cout);
			std::cerr << XE;
			exit(RELION_EXIT_FAILURE);
		}

		outDir = ZIO::makeOutputDir(outDir);


		std::vector<ParticleData> allParticleData;

		MetaDataTable inputTable, outParticles;
		inputTable.read(inStarFn);
		outParticles.setName("particles");

		if (!inputTable.containsLabel(EMDL_TOMO_NAME))
		{
			REPORT_ERROR_STR("Input star file " << inStarFn
				<< " does not contain tomogram names (rlnTomoName)");
		}

		if (inputTable.containsLabel(EMDL_IMAGE_COORD_X) &&
			inputTable.containsLabel(EMDL_IMAGE_COORD_Y) &&
			inputTable.containsLabel(EMDL_IMAGE_COORD_Z)    )
		{
			bool hasAngles;

			addParticles(inputTable, "", allParticleData, hasAngles);
			outParticles.append(inputTable);
		}
		else if (inputTable.containsLabel(EMDL_TOMO_IMPORT_PARTICLE_FILE))
		{
			int numberWithAngles = 0;

			for (int i = 0; i < inputTable.numberOfObjects(); i++)
			{
				const std::string subTableName = inputTable.getString(
							EMDL_TOMO_IMPORT_PARTICLE_FILE, i);

				const std::string tomoName = inputTable.getString(
							EMDL_TOMO_NAME, i);

				MetaDataTable inputSubTable;
				inputSubTable.read(subTableName);

				bool hasAngles;
				std::vector<ParticleData> particleData;

				addParticles(inputSubTable, tomoName, particleData, hasAngles);
				outParticles.append(inputSubTable);

				if (hasAngles)
				{
					numberWithAngles++;
				}

				allParticleData.reserve(allParticleData.size() + particleData.size());

				for (int p = 0; p < particleData.size(); p++)
				{
					allParticleData.push_back(particleData[p]);
				}
			}

			if (numberWithAngles != 0 && numberWithAngles != inputTable.numberOfObjects())
			{
				Log::warn("Only " + ZIO::itoa(numberWithAngles) + " out of "
					+ ZIO::itoa(inputTable.numberOfObjects()) + " partial particle files contain angles.");
			}
		}
		else
		{
			REPORT_ERROR_STR("Input star file " << inStarFn
				<< " contains neither particle coordinates (rlnCoordinateX/Y/Z)"
				<< " nor a list of STAR files that do (rlnTomoImportParticleFile)");
		}

		// split tomograms into optics groups

		std::map<std::string,std::vector<std::string>> opticsGroupName_to_tomoName;
		std::map<std::string, int> tomoSizeZ;

		TomogramSet tomogramSet(tomoFn);

		const int tc = tomogramSet.size();

		for (int t = 0; t < tc; t++)
		{
			const std::string opticsGroupName = tomogramSet.getOpticsGroupName(t);
			const std::string tomogramName = tomogramSet.getTomogramName(t);

			opticsGroupName_to_tomoName[opticsGroupName].push_back(tomogramName);

			if (flipZ)
			{
				int zSize;
				tomogramSet.globalTable.getValueSafely(EMDL_TOMO_SIZE_Z, zSize, t);
				tomoSizeZ[tomogramName] = zSize - 1;
			}
		}


		std::map<std::string, int> tomoName_to_optGroup;
		int optGroup = 0;

		MetaDataTable opticsTable;
		opticsTable.setName("optics");

		for (std::map<std::string,std::vector<std::string>>::iterator it = opticsGroupName_to_tomoName.begin();
			 it != opticsGroupName_to_tomoName.end(); it++)
		{
			const std::string opticsGroupName = it->first;
			const std::vector<std::string> allTomoNames = it->second;

			for (int t = 0; t < allTomoNames.size(); t++)
			{
				tomoName_to_optGroup[allTomoNames[t]] = optGroup;
			}

			Tomogram tomogram = tomogramSet.loadTomogram(
						tomogramSet.getTomogramIndex(allTomoNames[0]), false);

			opticsTable.addObject();

			opticsTable.setValue(EMDL_IMAGE_OPTICS_GROUP, optGroup + 1, optGroup);
			opticsTable.setValue(EMDL_IMAGE_OPTICS_GROUP_NAME, opticsGroupName, optGroup);
			opticsTable.setValue(EMDL_CTF_CS, tomogram.optics.Cs, optGroup);
			opticsTable.setValue(EMDL_CTF_VOLTAGE, tomogram.optics.voltage, optGroup);
			opticsTable.setValue(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, tomogram.optics.pixelSize, optGroup);

			// _rlnMicrographBinning, _rlnCtfDataAreCtfPremultiplied, _rlnImageDimensionality,
			// _rlnImageSize and _rlnImagePixelSize need to be specified by subtomo

			optGroup++;
		}



		int particle_id = 1;
		std::string last_tomo_name = "";

		for (int p = 0; p < allParticleData.size(); p++)
		{

			const ParticleData& d = allParticleData[p];

			if (tomoName_to_optGroup.find(d.tomoName) == tomoName_to_optGroup.end())
			{
				REPORT_ERROR_STR("No tomogram named '" << d.tomoName << "' found in " << tomoFn);
			}

			const int opticsGroup = tomoName_to_optGroup[d.tomoName];

			if (d.tomoName != last_tomo_name)
			{
				particle_id = 1;
				last_tomo_name = d.tomoName;
			}

			outParticles.setValue(EMDL_TOMO_NAME, d.tomoName, p);
			outParticles.setValue(EMDL_TOMO_PARTICLE_NAME, d.tomoName + "/"
															   + ZIO::itoa(particle_id), p);
			int subset;
			if (!inputTable.getValue(EMDL_PARTICLE_RANDOM_SUBSET, subset, p))
			{
				outParticles.setValue(EMDL_PARTICLE_RANDOM_SUBSET, p%2 + 1, p);
			}

			outParticles.setValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup + 1, p);

			outParticles.setValue(EMDL_IMAGE_COORD_X, d.coordinates.x, p);
			outParticles.setValue(EMDL_IMAGE_COORD_Y, d.coordinates.y, p);
			if (flipZ)
			{
				outParticles.setValue(EMDL_IMAGE_COORD_Z, tomoSizeZ[d.tomoName] - d.coordinates.z, p);
			}
			else
			{
				outParticles.setValue(EMDL_IMAGE_COORD_Z, d.coordinates.z, p);
			}

			outParticles.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, 0.0, p);
			outParticles.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, 0.0, p);
			outParticles.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, 0.0, p);

			outParticles.setValue(EMDL_ORIENT_ROT,  d.angles[0], p);
			outParticles.setValue(EMDL_ORIENT_TILT, d.angles[1], p);
			outParticles.setValue(EMDL_ORIENT_PSI,  d.angles[2], p);

			particle_id++;
		}

		std::ofstream ofs(outDir + "particles.star");
		opticsTable.write(ofs);
		outParticles.write(ofs);

		OptimisationSet optimisationSet;
		optimisationSet.particles = outDir + "particles.star";
		optimisationSet.tomograms = tomoFn;

		optimisationSet.write(outDir + "optimisation_set.star");

	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}
