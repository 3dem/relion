#include <src/args.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/particle_set.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	std::string tomoSetFn, outFn;

	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		tomoSetFn = parser.getOption("--t", "Tomogram set");
		outFn = parser.getOption("--o", "Output file name");

		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	std::ifstream ifs(tomoSetFn);

	if (!ifs)
	{
		REPORT_ERROR_STR("update_tomogram_set: Unable to read " << tomoSetFn);
	}


	MetaDataTable globalTable;
	std::vector<MetaDataTable> tomogramTables;

	globalTable.readStar(ifs, "global");

	const int tc = globalTable.numberOfObjects();

	tomogramTables.resize(tc);

	std::vector<MetaDataTable> allTables = MetaDataTable::readAll(ifs, tc+1);

	bool namesAreOld = false;

	for (int t = 0; t < tc; t++)
	{
		const std::string expectedOldName = "tomo_" + ZIO::itoa(t);
		const std::string expectedNewName = globalTable.getString(EMDL_TOMO_NAME, t);
		const std::string name = allTables[t+1].getName();

		if (name == expectedOldName)
		{
			namesAreOld = true;
		}
		else if (name != expectedNewName)
		{
			REPORT_ERROR_STR("TomogramSet::TomogramSet: file is corrupted " << tomoSetFn);
		}

		tomogramTables[t] = allTables[t+1];
		tomogramTables[t].setName(expectedNewName);
	}

	if (!namesAreOld)
	{
		Log::warn("Tomogram set " + tomoSetFn + " was already up to date.");
	}

	globalTable.setName("global");

	if (outFn.find_last_of('/') != std::string::npos)
	{
		std::string path = outFn.substr(0, outFn.find_last_of('/'));
		mktree(path);
	}

	std::ofstream ofs(outFn);

	if (!ofs)
	{
		REPORT_ERROR("update_tomogram_set: unable to write to " + outFn);
	}

	globalTable.write(ofs);

	for (int t = 0; t < tc; t++)
	{
		tomogramTables[t].write(ofs);
	}

	return 0;
}
