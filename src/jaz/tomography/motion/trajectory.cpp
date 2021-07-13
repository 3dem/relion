#include "trajectory.h"
#include <src/metadata_table.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <fstream>

using namespace gravis;


std::vector<Trajectory> Trajectory::read(std::string filename, ParticleSet& particleSet)
{
	std::ifstream ifs(filename);

	if (!ifs)
	{
		REPORT_ERROR("Trajectory::read: unable to read " + filename + ".");
	}

	MetaDataTable mdt0;

	mdt0.readStar(ifs, "general");

	const int pc = particleSet.getTotalParticleNumber();

	{
		int trc;

		if (!mdt0.getValue(EMDL_PARTICLE_NUMBER, trc))
		{
			REPORT_ERROR("Trajectory::read: missing particle number in "+filename+".");
		}

		if (trc < pc)
		{
			REPORT_ERROR("Trajectory::read: insufficient number of particles in "+filename+".");
		}
	}

	std::vector<Trajectory> out(pc);
	
	std::vector<MetaDataTable> mdts = MetaDataTable::readAll(ifs, pc+1);

	if (mdts[1].getName() == "0")
	{
		Log::warn("The particle trajectories in "+filename+
			" do not appear to have names. Please do not edit the particle table under any circumstances.");

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
	}
	else
	{
		const long int tableCount = mdts.size();

		std::map<std::string, long int> name_to_table;

		for (int t = 1; t < tableCount; t++)
		{
			const std::string name = mdts[t].getName();

			name_to_table[name] = t;
		}

		for (int pp = 0; pp < pc; pp++)
		{
			const std::string name = particleSet.getName(ParticleIndex(pp));

			if (name_to_table.find(name) == name_to_table.end())
			{
				REPORT_ERROR_STR("Trajectory::read: no trajectory found for particle '"
					<< name << "' in " << filename);
			}

			MetaDataTable& mdt = mdts[name_to_table[name]];

			int fc = mdt.numberOfObjects();

			out[pp] = Trajectory(fc);

			for (int f = 0; f < fc; f++)
			{
				mdt.getValueSafely(EMDL_ORIENT_ORIGIN_X_ANGSTROM, out[pp].shifts_Ang[f].x, f);
				mdt.getValueSafely(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, out[pp].shifts_Ang[f].y, f);
				mdt.getValueSafely(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, out[pp].shifts_Ang[f].z, f);
			}
		}
	}
	
	return out;
}

void Trajectory::write(
		const std::vector<std::vector<Trajectory>>& shifts,
		const ParticleSet& particleSet,
		const std::vector<std::vector<ParticleIndex>>& particles,
		std::string filename)
{
	const int tc = particles.size();
	int pc = 0;

	for (int t = 0; t < tc; t++)
	{
		pc += shifts[t].size();
	}

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

	for (int t = 0; t < tc; t++)
	for (int pp = 0; pp < particles[t].size(); pp++)
	{
		mdt.setName(particleSet.getName(particles[t][pp]));

		const int fc = shifts[t][pp].shifts_Ang.size();

		for (int f = 0; f < fc; f++)
		{
			mdt.addObject();
			
			mdt.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, shifts[t][pp].shifts_Ang[f].x);
			mdt.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, shifts[t][pp].shifts_Ang[f].y);
			mdt.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, shifts[t][pp].shifts_Ang[f].z);
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

Trajectory operator + (const Trajectory& t1, const Trajectory& t2)
{
	Trajectory out = t1;
	out += t2;
	return out;
}
