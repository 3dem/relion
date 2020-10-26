#include "trajectory.h"
#include <src/metadata_table.h>
#include <src/jaz/util/zio.h>
#include <fstream>

using namespace gravis;


std::vector<Trajectory> Trajectory::read(std::string filename)
{
	std::ifstream ifs(filename);

	if (!ifs)
	{
		REPORT_ERROR("Trajectory::read: unable to read " + filename + ".");
	}

	MetaDataTable mdt0;

	mdt0.readStar(ifs, "general");

	int pc;

	if (!mdt0.getValue(EMDL_PARTICLE_NUMBER, pc))
	{
		REPORT_ERROR("Trajectory::read: missing particle number in "+filename+".");
	}

	std::vector<Trajectory> out(pc);
	
	
	std::vector<MetaDataTable> mdts = MetaDataTable::readAll(ifs, pc+1);

	for (int p = 0; p < pc; p++)
	{
		MetaDataTable& mdt = mdts[p+1];

		int fc = mdt.numberOfObjects();
		
		out[p] = Trajectory(fc);

		for (int f = 0; f < fc; f++)
		{
			mdt.getValueSafely(EMDL_ORIENT_ORIGIN_X_ANGSTROM, out[p].shifts_Ang[f].x, f);
			mdt.getValueSafely(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, out[p].shifts_Ang[f].y, f);
			mdt.getValueSafely(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, out[p].shifts_Ang[f].z, f);
		}
	}
	
	return out;
}

void Trajectory::write(const std::vector<Trajectory>& shifts, std::string filename)
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
		
		const int fc = shifts[p].shifts_Ang.size();

		for (int f = 0; f < fc; f++)
		{
			mdt.addObject();
			
			mdt.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, shifts[p].shifts_Ang[f].x);
			mdt.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, shifts[p].shifts_Ang[f].y);
			mdt.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, shifts[p].shifts_Ang[f].z);
		}

		mdt.write(ofs);
		mdt.clear();
	}
}

Trajectory::Trajectory()
{	
}

Trajectory::Trajectory(int fc)
	:	shifts_Ang(fc, d3Vector(0.0, 0.0, 0.0))
{	
}

std::vector<d3Vector> Trajectory::getShiftsInPix(d3Vector origin, double pixelSize) const
{
	std::vector<d3Vector> out(shifts_Ang.size());

	const int fc = shifts_Ang.size();

	for (int f = 0; f < fc; f++)
	{
		out[f] = origin + shifts_Ang[f] / pixelSize;
	}
	
	return out;
}
