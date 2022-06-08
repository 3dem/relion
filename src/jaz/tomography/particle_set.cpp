#include "particle_set.h"
#include "tomogram_set.h"
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/math/Euler_angles_relion.h>
#include <sstream>
#include <set>
#include <map>

using namespace gravis;


ParticleSet::ParticleSet()
{}

ParticleSet::ParticleSet(std::string filename, std::string motionFilename, bool verbose)
{
	optTable.read(filename, "optics");
	partTable.read(filename, "particles");
	
	if (!optTable.containsLabel(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE))
	{
		REPORT_ERROR("ParticleSet::ParticleSet: "
					 + EMDL::label2Str(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE)
					 + " missing from optics MetaDataTable.\n");
	}

	if (!partTable.containsLabel(EMDL_TOMO_PARTICLE_NAME))
	{
		if (verbose)
		{
			Log::warn("The particles in "+filename+
				" do not have names (rlnTomoParticleName). They are being added now.");
		}

		std::map<std::string, int> tomoParticleCount;

		for (int p = 0; p < partTable.numberOfObjects(); p++)
		{
			const std::string tomoName = partTable.getString(EMDL_TOMO_NAME, p);

			if (tomoParticleCount.find(tomoName) == tomoParticleCount.end())
			{
				tomoParticleCount[tomoName] = 1;
			}

			const int id = tomoParticleCount[tomoName];
			tomoParticleCount[tomoName]++;

			partTable.setValue(EMDL_TOMO_PARTICLE_NAME, tomoName + "/" + ZIO::itoa(id), p);
		}
	}
    else
    {
        // If we do have particle names, make sure the input particles are sorted on their name, as this is implicitly assumed for the motion trajectories
        // This will sort 1, 11, 12, 2, 3, 4, 5, ...
        // But that's OK. The important thing is that all particles from each tomogram are together
        // SHWS 8jun2022: because checkTrajectories has now been repaired, I don't think this is necessary anymore....
        //partTable.newSort(EMDL_TOMO_PARTICLE_NAME);

    }

	hasMotion = motionFilename != "";

	if (hasMotion)
	{
		motionTrajectories = Trajectory::read(motionFilename, *this);
	}
}

void ParticleSet::reserve(int particleNumber)
{
	motionTrajectories.reserve(particleNumber);
}

ParticleIndex ParticleSet::addParticle(const ParticleSet &particleSet, ParticleIndex index)
{
	partTable.addObject(particleSet.partTable.getObject(index.value));

	if (particleSet.hasMotion)
	{
		if (!hasMotion && partTable.numberOfObjects() > 0)
		{
			REPORT_ERROR("ParticleSet::addParticle: trying to add a particle with motion to a particle set without.");
		}

		hasMotion = true;

		if (partTable.numberOfObjects() != motionTrajectories.size() + 1)
		{
			REPORT_ERROR("ParticleSet::addParticle: number of particles is not the same as the number of trajectories.");
		}

		motionTrajectories.push_back(particleSet.motionTrajectories[index.value]);
	}

	return ParticleIndex(partTable.numberOfObjects() - 1);
}

void ParticleSet::clearParticles()
{
	partTable.clear();
	partTable.setName("particles");
	motionTrajectories.clear();
}

std::vector<std::vector<ParticleIndex> > ParticleSet::splitByTomogram(const TomogramSet& tomogramSet, bool verbose) const
{
	std::vector<std::vector<ParticleIndex>> out(0);

	if (!partTable.containsLabel(EMDL_TOMO_NAME))
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
		out[t] = std::vector<ParticleIndex>(0);
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

		if (verbose)
		{
			Log::warn("No particles were found for the following tomograms: "+empty_names_list.str());
		}
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

		if (verbose)
		{
			Log::warn("Particles were found belonging to the following undefined tomograms: "+empty_names_list.str());
		}
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
			out[t].push_back(ParticleIndex(i));
		}
	}

	return out;
}

int ParticleSet::getTotalParticleNumber() const
{
	return partTable.numberOfObjects();
}

d3Vector ParticleSet::getPosition(ParticleIndex particle_id) const
{
	const int og = getOpticsGroup(particle_id);
	
	const double originalPixelSize = optTable.getDouble(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, og);

	const d3Matrix A_subtomogram = getSubtomogramMatrix(particle_id);

	d3Vector out = getParticleCoord(particle_id) - (A_subtomogram * getParticleOffset(particle_id)) / originalPixelSize;
	
	out.x += 1.0;
	out.y += 1.0;
	out.z += 1.0;
	
	return out;
}

d3Matrix ParticleSet::getSubtomogramMatrix(ParticleIndex particle_id) const
{
	if (partTable.containsLabel(EMDL_TOMO_SUBTOMOGRAM_ROT))
	{
		const double phi  =  partTable.getAngleInRad(EMDL_TOMO_SUBTOMOGRAM_ROT,  particle_id.value);
		const double theta = partTable.getAngleInRad(EMDL_TOMO_SUBTOMOGRAM_TILT, particle_id.value);
		const double psi  =  partTable.getAngleInRad(EMDL_TOMO_SUBTOMOGRAM_PSI,  particle_id.value);

		return Euler::anglesToMatrix3(phi, theta, psi);
	}
	else
	{
		return d3Matrix(
			1, 0, 0,
			0, 1, 0,
			0, 0, 1 );
	}
}

d3Matrix ParticleSet::getParticleMatrix(ParticleIndex particle_id) const
{
	if (partTable.containsLabel(EMDL_ORIENT_ROT))
	{
		const double phi  =  partTable.getAngleInRad(EMDL_ORIENT_ROT,  particle_id.value);
		const double theta = partTable.getAngleInRad(EMDL_ORIENT_TILT, particle_id.value);
		const double psi  =  partTable.getAngleInRad(EMDL_ORIENT_PSI,  particle_id.value);

		return Euler::anglesToMatrix3(phi, theta, psi);
	}
	else
	{
		return d3Matrix(
				1, 0, 0,
				0, 1, 0,
				0, 0, 1 );
	}
}

d3Matrix ParticleSet::getMatrix3x3(ParticleIndex particle_id) const
{
	const d3Matrix A_particle = getParticleMatrix(particle_id);
	const d3Matrix A_subtomogram = getSubtomogramMatrix(particle_id);
	
	return A_subtomogram * A_particle;
}


// This maps coordinates from particle space to tomogram space.
d4Matrix ParticleSet::getMatrix4x4(ParticleIndex particle_id, double w, double h, double d) const
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
		A(0,0), A(0,1), A(0,2), 0.0,
		A(1,0), A(1,1), A(1,2), 0.0,
		A(2,0), A(2,1), A(2,2), 0.0,
		0.0,    0.0,    0.0,    1.0   );
	
	d4Matrix Ts(
		1, 0, 0, pos.x,
		0, 1, 0, pos.y,
		0, 0, 1, pos.z,
		0, 0, 0, 1);
	
	return Ts * R * Tc;
}

t4Vector<d3Matrix> ParticleSet::getMatrixDerivativesOverParticleAngles(
		ParticleIndex particle_id) const
{
	const double phi  =  partTable.getAngleInRad(EMDL_ORIENT_ROT,  particle_id.value);
	const double theta = partTable.getAngleInRad(EMDL_ORIENT_TILT, particle_id.value);
	const double psi  =  partTable.getAngleInRad(EMDL_ORIENT_PSI,  particle_id.value);

	t4Vector<d3Matrix> dA_particle =
		Euler::anglesToMatrixAndDerivatives(phi, theta, psi);

	const d3Matrix A_subtomogram = getSubtomogramMatrix(particle_id);

	for (int i = 0; i < 4; i++)
	{
		dA_particle[i] = A_subtomogram * dA_particle[i];
	}

	return dA_particle;
}

std::string ParticleSet::getName(ParticleIndex particle_id) const
{
	return partTable.getString(EMDL_TOMO_PARTICLE_NAME, particle_id.value);
}

int ParticleSet::getHalfSet(ParticleIndex particle_id) const
{
	int s;
	partTable.getValueSafely(EMDL_PARTICLE_RANDOM_SUBSET, s, particle_id.value);
	return s - 1;
}


void ParticleSet::moveParticleTo(ParticleIndex particle_id, gravis::d3Vector pos)
{
	partTable.setValue(EMDL_IMAGE_COORD_X, pos.x - 1.0, particle_id.value);
	partTable.setValue(EMDL_IMAGE_COORD_Y, pos.y - 1.0, particle_id.value);
	partTable.setValue(EMDL_IMAGE_COORD_Z, pos.z - 1.0, particle_id.value);
	
	partTable.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, 0.0, particle_id.value);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, 0.0, particle_id.value);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, 0.0, particle_id.value);
}

void ParticleSet::shiftParticleBy(ParticleIndex particle_id, gravis::d3Vector shift)
{
	double x, y, z;
	
	partTable.getValueSafely(EMDL_IMAGE_COORD_X, x, particle_id.value);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Y, y, particle_id.value);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Z, z, particle_id.value);
	
	partTable.setValue(EMDL_IMAGE_COORD_X, x + shift.x, particle_id.value);
	partTable.setValue(EMDL_IMAGE_COORD_Y, y + shift.y, particle_id.value);
	partTable.setValue(EMDL_IMAGE_COORD_Z, z + shift.z, particle_id.value);
}

void ParticleSet::write(const std::string& filename) const
{
	std::ofstream ofs(filename);

	optTable.write(ofs);
	partTable.write(ofs);
}

void ParticleSet::writeTrajectories(const std::string &filename) const
{
	const int pc = motionTrajectories.size();

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

	for (int pp = 0; pp < pc; pp++)
	{
		mdt.setName(getName(ParticleIndex(pp)));

		const int fc = motionTrajectories[pp].shifts_Ang.size();

		for (int f = 0; f < fc; f++)
		{
			mdt.addObject();

			mdt.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, motionTrajectories[pp].shifts_Ang[f].x);
			mdt.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, motionTrajectories[pp].shifts_Ang[f].y);
			mdt.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, motionTrajectories[pp].shifts_Ang[f].z);
		}

		mdt.write(ofs);
		mdt.clear();
	}
}

void ParticleSet::setImageFileNames(std::string data, std::string weight, ParticleIndex particle_id)
{
	partTable.setValue(EMDL_IMAGE_NAME, data, particle_id.value);
	partTable.setValue(EMDL_CTF_IMAGE, weight, particle_id.value);
}

d3Vector ParticleSet::getParticleOffset(ParticleIndex particle_id) const
{
	d3Vector out;
	
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_X_ANGSTROM, out.x, particle_id.value);
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, out.y, particle_id.value);
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, out.z, particle_id.value);
	
	return out;
}

void ParticleSet::setParticleOffset(ParticleIndex particle_id, const d3Vector& v)
{
	partTable.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, v.x, particle_id.value);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, v.y, particle_id.value);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, v.z, particle_id.value);
}

d3Vector ParticleSet::getParticleCoord(ParticleIndex particle_id) const
{
	d3Vector out;
	
	partTable.getValueSafely(EMDL_IMAGE_COORD_X, out.x, particle_id.value);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Y, out.y, particle_id.value);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Z, out.z, particle_id.value);
	
	return out;
}

void ParticleSet::setParticleCoord(ParticleIndex particle_id, const d3Vector& v)
{
	partTable.setValue(EMDL_IMAGE_COORD_X, v.x, particle_id.value);
	partTable.setValue(EMDL_IMAGE_COORD_Y, v.y, particle_id.value);
	partTable.setValue(EMDL_IMAGE_COORD_Z, v.z, particle_id.value);
}

int ParticleSet::getOpticsGroup(ParticleIndex particle_id) const
{
	if (!partTable.containsLabel(EMDL_IMAGE_OPTICS_GROUP))
	{
		REPORT_ERROR("ParticleSet::getPixelSize: pixel size (rlnImagePixelSize) missing from optics table");
	}
	
	int out;
	partTable.getValueSafely(EMDL_IMAGE_OPTICS_GROUP, out, particle_id.value);
	return out - 1;
}

void ParticleSet::setOpticsGroup(ParticleIndex particle_id, int zeroBasedId)
{
	partTable.setValue(EMDL_IMAGE_OPTICS_GROUP, zeroBasedId+1, particle_id.value);
}

int ParticleSet::numberOfOpticsGroups() const
{
	return optTable.numberOfObjects();
}

double ParticleSet::getBinnedPixelSize(int opticsGroup) const
{
	if (!optTable.containsLabel(EMDL_IMAGE_PIXEL_SIZE))
	{
		REPORT_ERROR("ParticleSet::getBinnedPixelSize: pixel size (rlnImagePixelSize) missing from optics table");
	}
	
	double out;
	optTable.getValueSafely(EMDL_IMAGE_PIXEL_SIZE, out, opticsGroup);
	return out;
}

double ParticleSet::getOriginalPixelSize(int opticsGroup) const
{
	if (!optTable.containsLabel(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE))
	{
		REPORT_ERROR("ParticleSet::getOriginalPixelSize: tilt series pixel size (rlnTomoTiltSeriesPixelSize) missing from optics table");
	}
	
	double out;
	optTable.getValueSafely(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, out, opticsGroup);
	return out;
}

std::vector<d3Vector> ParticleSet::getTrajectoryInPixels(ParticleIndex particle_id, int fc, double pixelSize, bool from_original_coordinate) const
{
	const d3Vector p0 = (from_original_coordinate) ? getParticleCoord(particle_id) : getPosition(particle_id);

	if (hasMotion)
	{
		return motionTrajectories[particle_id.value].getShiftsInPix(p0, pixelSize);
	}
	else
	{
		return std::vector<d3Vector>(fc, p0);
	}
}

void ParticleSet::checkTrajectoryLengths(const std::vector<ParticleIndex> &tomogram_particles, int fc, std::string caller) const
{
	if (hasMotion)
	{
		for (int i = 0; i < tomogram_particles.size(); i++)
		{
            long int p = tomogram_particles[i].value;
            if (motionTrajectories[p].shifts_Ang.size() != fc)
			{
				std::cerr << " i= " << i << " p= " << p << " fc= " << fc << " name= " << getName(tomogram_particles[i]) << std::endl;
                REPORT_ERROR_STR(caller << ": bad trajectory lengths; expected " << fc << " frames, found "
								 << motionTrajectories[p].shifts_Ang.size());
			}
		}
	}
}

std::vector<std::vector<int>> ParticleSet::splitEvenly(
		const std::vector<std::vector<ParticleIndex>>& particlesByTomogram,
		int segment_count)
{
	const int tc = particlesByTomogram.size();
	const int sc = segment_count;

	std::vector<int> tomo_weight(tc);
	std::vector<std::set<int>> segments(sc);
	std::vector<int> segment_weight(sc, 0);

	int total_weight = 0;

	for (int t = 0; t < tc; t++)
	{
		tomo_weight[t] = particlesByTomogram[t].size();
		total_weight += tomo_weight[t];

		segments[0].insert(t);
	}


	while (true)
	{
		for (int s = 0; s < sc; s++)
		{
			segment_weight[s] = 0;

			for (std::set<int>::iterator it = segments[s].begin();
				 it != segments[s].end(); it++)
			{
				segment_weight[s] += tomo_weight[*it];
			}
		}

		std::vector<int> order = IndexSort<int>::sortIndices(segment_weight);

		int heaviest_segment = order[order.size()-1];
		int lightest_segment = order[0];

		int lightest_tomo = 0;
		int lightest_tomo_weight = total_weight;

		for (std::set<int>::iterator it = segments[heaviest_segment].begin();
			 it != segments[heaviest_segment].end(); it++)
		{
			if (tomo_weight[*it] < lightest_tomo_weight)
			{
				lightest_tomo = *it;
				lightest_tomo_weight = tomo_weight[*it];
			}
		}

		if (segment_weight[lightest_segment] + tomo_weight[lightest_tomo]
				< segment_weight[heaviest_segment])
		{
			segments[lightest_segment].insert(lightest_tomo);
			segments[heaviest_segment].erase(lightest_tomo);

		}
		else
		{
			break;
		}
	}

	std::vector<std::vector<int>> out(sc);

	for (int s = 0; s < sc; s++)
	{
		for (std::set<int>::iterator it = segments[s].begin();
			 it != segments[s].end(); it++)
		{
			out[s].push_back(*it);
		}
	}

	return out;
}

std::vector<int> ParticleSet::enumerate(
		const std::vector<std::vector<ParticleIndex>>& particlesByTomogram)
{
	const int tc = particlesByTomogram.size();

	std::vector<int> indices(tc);

	for (int t = 0; t < tc; t++)
	{
		indices[t] = t;
	}

	return indices;
}

std::vector<int> ParticleSet::enumerateNonEmpty(const std::vector<std::vector<ParticleIndex> > &particlesByTomogram)
{
	const int tc0 = particlesByTomogram.size();
	int count = 0;

	for (int t = 0; t < tc0; t++)
	{
		if (particlesByTomogram[t].size() > 0)
		{
			count++;
		}
	}

	std::vector<int> out(count);
	int index = 0;

	for (int t = 0; t < tc0; t++)
	{
		if (particlesByTomogram[t].size() > 0)
		{
			out[index] = t;
			index++;
		}
	}

	return out;
}


