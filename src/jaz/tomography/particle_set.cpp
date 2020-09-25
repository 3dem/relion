#include "particle_set.h"
#include "tomogram_set.h"
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>
#include <src/euler.h>
#include <sstream>
#include <set>
#include <map>

using namespace gravis;


ParticleSet* ParticleSet::load(std::string filename, std::string motionFilename)
{
	return new ParticleSet(filename, motionFilename);
}

ParticleSet::ParticleSet()
{}

ParticleSet::ParticleSet(std::string filename, std::string motionFilename)
{
	optTable.read(filename, "optics");
	partTable.read(filename, "particles");
	
	if (!optTable.labelExists(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE))
	{
		REPORT_ERROR("ParticleSet::ParticleSet: "
					 + EMDL::label2Str(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE)
					 + " missing from optics MetaDataTable.\n");
	}

	hasMotion = motionFilename != "";

	if (hasMotion)
	{
		motionTrajectories = Trajectory::read(motionFilename);
	}
}

std::vector<std::vector<int> > ParticleSet::splitByTomogram(const TomogramSet& tomogramSet) const
{
	std::vector<std::vector<int>> out(0);

	if (!partTable.labelExists(EMDL_TOMO_NAME))
	{
		REPORT_ERROR("ParticleSet::splitByTomogram: "
					 + EMDL::label2Str(EMDL_TOMO_NAME)
					 + " missing from MetaDataTable. Please run add_tomo_name on this particle set.\n");
	}

	const int tc = tomogramSet.size();

	out.resize(tc);

	std::set<std::string> unknown_tomo_names;
	std::vector<bool> any_particles_found(tc, false);

	std::map<std::string, int> name_to_index;

	for (int t = 0; t < tc; t++)
	{
		const std::string name = tomogramSet.globalTable.getString(EMDL_TOMO_NAME, t);
		name_to_index[name] = t;
	}

	std::vector<int> pc_t(tc, 0);

	const long pc = partTable.numberOfObjects();

	for (int i = 0; i < pc; i++)
	{
		const std::string name = partTable.getString(EMDL_TOMO_NAME, i);

		if (name_to_index.find(name) == name_to_index.end())
		{
			unknown_tomo_names.insert(name);
		}
		else
		{
			const int t = name_to_index[name];
			pc_t[t]++;

			any_particles_found[t] = true;
		}
	}

	bool any_empty = false;
	bool all_empty = true;

	for (int t = 0; t < tc; t++)
	{
		out[t] = std::vector<int>(0);
		out[t].reserve(pc_t[t]);

		if (any_particles_found[t])
		{
			all_empty = false;
		}
		else
		{
			any_empty = true;
		}
	}

	if (any_empty)
	{
		std::stringstream empty_names_list;
		int empty_so_far = 0;

		for (int t = 0; t < tc; t++)
		{
			if (!any_particles_found[t])
			{
				if (empty_so_far > 0)
				{
					empty_names_list << ", ";
				}

				empty_names_list << tomogramSet.globalTable.getString(EMDL_TOMO_NAME, t);
				empty_so_far++;
			}
		}

		Log::warn("No particles were found for the following tomograms: "+empty_names_list.str());
	}

	if (!unknown_tomo_names.empty())
	{
		std::stringstream empty_names_list;
		int empty_so_far = 0;

		for (std::set<std::string>::iterator it = unknown_tomo_names.begin();
			 it != unknown_tomo_names.end(); it++)
		{
			if (empty_so_far > 0)
			{
				empty_names_list << ", ";
			}

			empty_names_list << *it;
			empty_so_far++;
		}

		Log::warn("Particles were found belonging to the following undefined tomograms: "+empty_names_list.str());
	}

	if (all_empty)
	{
		REPORT_ERROR_STR(
			"ParticleSet::splitByTomogram: no particles found that belong to any of the defined tomograms. "
			<< "Please compare the tomogram names in the particle file to those in the tomogram file.");
	}

	for (int i = 0; i < pc; i++)
	{
		const std::string name = partTable.getString(EMDL_TOMO_NAME, i);

		if (name_to_index.find(name) != name_to_index.end())
		{
			const int t = name_to_index[name];
			out[t].push_back(i);
		}
	}

	return out;
}

int ParticleSet::getTotalParticleNumber() const
{
	return partTable.numberOfObjects();
}

d3Vector ParticleSet::getPosition(long int particle_id) const
{
	d3Vector pos, off;
	
	partTable.getValueSafely(EMDL_IMAGE_COORD_X, pos.x, particle_id);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Y, pos.y, particle_id);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Z, pos.z, particle_id);
	
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_X_ANGSTROM, off.x, particle_id);
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, off.y, particle_id);
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, off.z, particle_id);
	
	const int og = getOpticsGroup(particle_id);
	
	const double originalPixelSize = optTable.getDouble(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, og);
	
	d3Vector out = pos - off / originalPixelSize;
	
	out.x += 1.0;
	out.y += 1.0;
	out.z += 1.0;
	
	return out;
}

d3Matrix ParticleSet::getSubtomogramMatrix(long particle_id) const
{
	if (partTable.labelExists(EMDL_TOMO_SUBTOMOGRAM_ROT))
	{
		double phi, theta, psi;
		
		partTable.getValueSafely(EMDL_TOMO_SUBTOMOGRAM_ROT,  phi,   particle_id);
		partTable.getValueSafely(EMDL_TOMO_SUBTOMOGRAM_TILT, theta, particle_id);
		partTable.getValueSafely(EMDL_TOMO_SUBTOMOGRAM_PSI,  psi,   particle_id);
		
		Matrix2D<double> A;
		Euler_angles2matrix(phi, theta, psi, A, false);	
		return convert(A);
	}
	else
	{
		return d3Matrix(
			1, 0, 0,
			0, 1, 0,
			0, 0, 1 );
	}
}

d3Matrix ParticleSet::getParticleMatrix(long particle_id) const
{
	double phi, theta, psi;
	
	partTable.getValueSafely(EMDL_ORIENT_ROT,  phi,   particle_id);
	partTable.getValueSafely(EMDL_ORIENT_TILT, theta, particle_id);
	partTable.getValueSafely(EMDL_ORIENT_PSI,  psi,   particle_id);
	
	Matrix2D<double> A;
	Euler_angles2matrix(phi, theta, psi, A, false);	
	return convert(A);
	
}

d3Matrix ParticleSet::getMatrix3x3(long int particle_id) const
{
	const d3Matrix A_particle = getParticleMatrix(particle_id);
	const d3Matrix A_subtomogram = getSubtomogramMatrix(particle_id);
	
	return A_subtomogram * A_particle;
}

d4Matrix ParticleSet::getMatrix4x4(long int particle_id, double w, double h, double d) const
{
	d3Matrix A = getMatrix3x3(particle_id);	
	d3Vector pos = getPosition(particle_id);
	
	int cx = ((int)w) / 2;
	int cy = ((int)h) / 2;
	int cz = ((int)d) / 2;
	
	gravis::d4Matrix Tc(
		1, 0, 0, -cx,
		0, 1, 0, -cy,
		0, 0, 1, -cz,
		0, 0, 0, 1);
	
	d4Matrix R(
		A(0,0), A(0,1), A(0,2), pos.x, 
		A(1,0), A(1,1), A(1,2), pos.y, 
		A(2,0), A(2,1), A(2,2), pos.z, 
		0.0,    0.0,    0.0,    1.0   );
	
	d4Matrix Ts(
		1, 0, 0, pos.x,
		0, 1, 0, pos.y,
		0, 0, 1, pos.z,
		0, 0, 0, 1);
	
	return Ts * R * Tc;
}

std::string ParticleSet::getName(long int particle_id) const
{
	std::stringstream sts;
	sts << particle_id;
	
	return sts.str();
}

int ParticleSet::getHalfSet(long int particle_id) const
{
	int s;
	partTable.getValueSafely(EMDL_PARTICLE_RANDOM_SUBSET, s, particle_id);
	return s - 1;
}


void ParticleSet::moveParticleTo(long int particle_id, gravis::d3Vector pos)
{
	partTable.setValue(EMDL_IMAGE_COORD_X, pos.x - 1.0, particle_id);
	partTable.setValue(EMDL_IMAGE_COORD_Y, pos.y - 1.0, particle_id);
	partTable.setValue(EMDL_IMAGE_COORD_Z, pos.z - 1.0, particle_id);
	
	partTable.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, 0.0, particle_id);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, 0.0, particle_id);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, 0.0, particle_id);
}

void ParticleSet::shiftParticleBy(long int particle_id, gravis::d3Vector shift)
{
	double x, y, z;
	
	partTable.getValueSafely(EMDL_IMAGE_COORD_X, x, particle_id);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Y, y, particle_id);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Z, z, particle_id);
	
	partTable.setValue(EMDL_IMAGE_COORD_X, x + shift.x, particle_id);
	partTable.setValue(EMDL_IMAGE_COORD_Y, y + shift.y, particle_id);
	partTable.setValue(EMDL_IMAGE_COORD_Z, z + shift.z, particle_id);
}

void ParticleSet::write(std::string fn) const
{
	std::ofstream ofs(fn);
	
	optTable.write(ofs);
	partTable.write(ofs);
}

void ParticleSet::setImageFileNames(std::string data, std::string weight, long int particle_id)
{
	partTable.setValue(EMDL_IMAGE_NAME, data, particle_id);
	partTable.setValue(EMDL_CTF_IMAGE, weight, particle_id);
}

d3Vector ParticleSet::getParticleOffset(long particle_id) const
{
	d3Vector out;
	
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_X_ANGSTROM, out.x, particle_id);
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, out.y, particle_id);
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, out.z, particle_id);
	
	return out;
}

void ParticleSet::setParticleOffset(long particle_id, const d3Vector& v)
{
	partTable.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, v.x, particle_id);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, v.y, particle_id);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, v.z, particle_id);
}

d3Vector ParticleSet::getParticleCoord(long particle_id) const
{
	d3Vector out;
	
	partTable.getValueSafely(EMDL_IMAGE_COORD_X, out.x, particle_id);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Y, out.y, particle_id);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Z, out.z, particle_id);
	
	return out;
}

void ParticleSet::setParticleCoord(long particle_id, const d3Vector& v)
{
	partTable.setValue(EMDL_IMAGE_COORD_X, v.x, particle_id);
	partTable.setValue(EMDL_IMAGE_COORD_Y, v.y, particle_id);
	partTable.setValue(EMDL_IMAGE_COORD_Z, v.z, particle_id);
}

int ParticleSet::getOpticsGroup(long particle_id) const
{
	if (!partTable.labelExists(EMDL_IMAGE_OPTICS_GROUP))
	{
		REPORT_ERROR("ParticleSet::getPixelSize: pixel size (rlnImagePixelSize) missing from optics table");
	}
	
	int out;
	partTable.getValueSafely(EMDL_IMAGE_OPTICS_GROUP, out, particle_id);
	return out - 1;
}

int ParticleSet::numberOfOpticsGroups() const
{
	return optTable.numberOfObjects();
}

double ParticleSet::getBinnedPixelSize(int opticsGroup) const
{
	if (!optTable.labelExists(EMDL_IMAGE_PIXEL_SIZE))
	{
		REPORT_ERROR("ParticleSet::getBinnedPixelSize: pixel size (rlnImagePixelSize) missing from optics table");
	}
	
	double out;
	optTable.getValueSafely(EMDL_IMAGE_PIXEL_SIZE, out, opticsGroup);
	return out;
}

double ParticleSet::getOriginalPixelSize(int opticsGroup) const
{
	if (!optTable.labelExists(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE))
	{
		REPORT_ERROR("ParticleSet::getOriginalPixelSize: pixel size (rlnMicrographOriginalPixelSize) missing from optics table");
	}
	
	double out;
	optTable.getValueSafely(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, out, opticsGroup);
	return out;
}

std::vector<d3Vector> ParticleSet::getTrajectoryInPixels(long particle_id, int fc, double pixelSize) const
{
	const d3Vector p0 = getPosition(particle_id);

	if (hasMotion)
	{
		return motionTrajectories[particle_id].getShiftsInPix(p0, pixelSize);
	}
	else
	{
		return std::vector<d3Vector>(fc, p0);
	}
}

void ParticleSet::checkTrajectoryLengths(int p0, int np, int fc, std::string caller) const
{
	if (hasMotion)
	{
		for (int p = p0; p < p0 + np; p++)
		{
			if (motionTrajectories[p].shifts_Ang.size() != fc)
			{
				REPORT_ERROR_STR(caller << ": bad trajectory lengths; expected " << fc << " frames, found "
								 << motionTrajectories[p].shifts_Ang.size());
			}
		}
	}
}

d3Matrix ParticleSet::convert(const Matrix2D<double> &A)
{
	return d3Matrix(
				A(0,0), A(0,1), A(0,2), 
				A(1,0), A(1,1), A(1,2), 
				A(2,0), A(2,1), A(2,2) );
}
