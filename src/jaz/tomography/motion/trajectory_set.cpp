#include "trajectory_set.h"
#include <src/metadata_table.h>
#include <src/jaz/util/zio.h>
#include <fstream>

using namespace gravis;


TrajectorySet::TrajectorySet()
{
	
}

TrajectorySet::TrajectorySet(const std::vector<std::vector<d3Vector>>& shiftsInPix, double pixelSize)
{
	const int pc = shiftsInPix.size();
	shifts.resize(pc);
	
	for (int p = 0; p < pc; p++)
	{
		const int fc = shiftsInPix[p].size();
		shifts[p].resize(fc);
		
		for (int f = 0; f < fc; f++)
		{
			shifts[p][f] = shiftsInPix[p][f] * pixelSize;
		}
	}
}

TrajectorySet::TrajectorySet(std::string filename)
{
	std::ifstream ifs(filename);

	if (!ifs)
	{
		REPORT_ERROR("TrajectorySet::TrajectorySet: unable to read " + filename + ".");
	}

	MetaDataTable mdt;

	mdt.readStar(ifs, "general");

	int pc;

	if (!mdt.getValue(EMDL_PARTICLE_NUMBER, pc))
	{
		REPORT_ERROR("TrajectorySet::TrajectorySet: missing particle number in "+filename+".");
	}

	std::vector<std::vector<d3Vector>> out(pc);

	for (int p = 0; p < pc; p++)
	{
		mdt.readStar(ifs, ZIO::itoa(p));

		int fc = mdt.numberOfObjects();

		out[p] = std::vector<d3Vector>(fc);

		for (int f = 0; f < fc; f++)
		{
			mdt.getValueSafely(EMDL_ORIENT_ORIGIN_X_ANGSTROM, out[p][f].x, f);
			mdt.getValueSafely(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, out[p][f].y, f);
			mdt.getValueSafely(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, out[p][f].z, f);
		}
	}
}

void TrajectorySet::write(std::string filename)
{
	const int pc = shifts.size();

	std::string path = filename.substr(0, filename.find_last_of('/'));
	mktree(path);

	std::ofstream ofs(filename);
	MetaDataTable mdt;

	mdt.setName("general");
	mdt.setIsList(true);
	mdt.addObject();
	mdt.setValue(EMDL_PARTICLE_NUMBER, pc);

	mdt.write(ofs);
	mdt.clear();

	for (int p = 0; p < pc; p++)
	{
		mdt.setName(ZIO::itoa(p));
		
		const int fc = shifts[p].size();

		for (int f = 0; f < fc; f++)
		{
			mdt.addObject();
			
			mdt.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, shifts[p][f].x);
			mdt.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, shifts[p][f].y);
			mdt.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, shifts[p][f].z);
		}

		mdt.write(ofs);
		mdt.clear();
	}
}
