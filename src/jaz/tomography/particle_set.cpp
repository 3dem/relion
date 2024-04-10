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

ParticleSet::ParticleSet(std::string filename, std::string motionFilename, bool verbose, const TomogramSet *tomogramSet)
{
    if (!read(filename, motionFilename, verbose, tomogramSet))
        REPORT_ERROR("ERROR: there are no particles in " + filename);
}

bool ParticleSet::read(std::string filename, std::string motionFilename, bool verbose, const TomogramSet *tomogramSet)
{

    if (genTable.read(filename,"general"))
    {
        genTable.getValueSafely(EMDL_TOMO_SUBTOMOGRAM_STACK2D, is_stack2d);
    }
    else
    {
        is_stack2d = false;
        genTable.setIsList(true);
        genTable.addObject();
        genTable.setValue(EMDL_TOMO_SUBTOMOGRAM_STACK2D, is_stack2d);
    }

    optTable.read(filename, "optics");
    partTable.read(filename, "particles");

    // subtomo can call to this function without a good particleSet yet, when coming straight out of subtomogram particle picking
    // In that case, initialise optics groups from the tomogramSet
    if (optTable.numberOfObjects() == 0)
    {

        // If the tomogramSet globalTable contains opticsGroupName, then use those groups, otherwise each tomogram is its own optics_group
        std::map<std::string, std::string> tomoname_to_opticsgroupname;
        size_t nr_tomos = tomogramSet->globalTable.numberOfObjects();
        for (size_t t = 0; t < nr_tomos; t++)
        {
            std::string tomo_name, optics_group_name;
            tomogramSet->globalTable.getValue(EMDL_TOMO_NAME, tomo_name, t);
            if (tomogramSet->globalTable.containsLabel(EMDL_IMAGE_OPTICS_GROUP_NAME))
                tomogramSet->globalTable.getValue(EMDL_IMAGE_OPTICS_GROUP_NAME, optics_group_name, t);
            else
                optics_group_name = tomo_name;
            tomoname_to_opticsgroupname.insert(std::make_pair(tomo_name, optics_group_name));
        }

        // Check which tomo_names are present in the partTable
        std::vector<std::string> present_optics_groupnames;
        std::string my_prev_name="";
        std::map<std::string, int> opticsgroupname_to_opticsgroup;
        std::map<std::string, std::string> presentopticsgroupname_to_firsttomoname;
        FOR_ALL_OBJECTS_IN_METADATA_TABLE(partTable)
        {
            std::string myname, myopticsgroupname;
            partTable.getValue(EMDL_TOMO_NAME, myname);
            myopticsgroupname = tomoname_to_opticsgroupname[myname];
            if (myopticsgroupname != my_prev_name)
            {
                bool is_new = true;
                for (size_t i = 0; i < present_optics_groupnames.size(); i++)
                {
                    if (myopticsgroupname == present_optics_groupnames[i]) is_new = false;
                }
                if (is_new)
                {
                    present_optics_groupnames.push_back(myopticsgroupname);
                    opticsgroupname_to_opticsgroup.insert(std::make_pair(myopticsgroupname, present_optics_groupnames.size()));
                    presentopticsgroupname_to_firsttomoname.insert(std::make_pair(myopticsgroupname, myname));
                }
            }
            my_prev_name = myopticsgroupname;
        }

        // construct optics table with those tomo_names that are present in the partTable
        optTable.setName("optics");
        if (!tomogramSet->globalTable.containsLabel(EMDL_CTF_VOLTAGE))
            REPORT_ERROR("ERROR: tomogramSet->globalTable does not contain rlnVoltage label");
         if (!tomogramSet->globalTable.containsLabel(EMDL_CTF_Q0))
            REPORT_ERROR("ERROR: tomogramSet->globalTable does not contain rlnAmplitudeContrast label");
       if (!tomogramSet->globalTable.containsLabel(EMDL_CTF_CS))
            REPORT_ERROR("ERROR: tomogramSet->globalTable does not contain rlnSphericalAberration label");
        if (!tomogramSet->globalTable.containsLabel(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE))
            REPORT_ERROR("ERROR: tomogramSet->globalTable does not contain rlnTomoTiltSeriesPixelSize label");
        for (size_t i = 0; i < present_optics_groupnames.size(); i++)
        {
            int idx = tomogramSet->getTomogramIndex(presentopticsgroupname_to_firsttomoname[present_optics_groupnames[i] ]);
            double Q0, Cs, kV, tiltSeriesPixelSize;
            tomogramSet->globalTable.getValue(EMDL_CTF_VOLTAGE, kV, idx);
            tomogramSet->globalTable.getValue(EMDL_CTF_CS, Cs, idx);
            tomogramSet->globalTable.getValue(EMDL_CTF_Q0, Q0, idx);
            tomogramSet->globalTable.getValue(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, tiltSeriesPixelSize, idx);

            optTable.addObject();
            optTable.setValue(EMDL_CTF_VOLTAGE, kV);
            optTable.setValue(EMDL_CTF_CS, Cs);
            optTable.setValue(EMDL_CTF_Q0, Q0);
            optTable.setValue(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, tiltSeriesPixelSize);
            optTable.setValue(EMDL_IMAGE_OPTICS_GROUP, opticsgroupname_to_opticsgroup[ present_optics_groupnames[i] ]);
            optTable.setValue(EMDL_IMAGE_OPTICS_GROUP_NAME, present_optics_groupnames[i]);
        }

        // Now also set optics groups in partTable
        FOR_ALL_OBJECTS_IN_METADATA_TABLE(partTable)
        {
            std::string myname;
            partTable.getValue(EMDL_TOMO_NAME, myname);
            partTable.setValue(EMDL_IMAGE_OPTICS_GROUP, opticsgroupname_to_opticsgroup[tomoname_to_opticsgroupname[myname] ]);
        }

    }

    
    long int pc = partTable.numberOfObjects();
    if (pc == 0) return false;

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

    return true;
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

d3Vector ParticleSet::getPosition(ParticleIndex particle_id, const gravis::d3Vector &tomo_centre, bool apply_origin_shifts) const
{
	const int og = getOpticsGroup(particle_id);
	
	const double tiltSeriesPixelSize = getTiltSeriesPixelSize(og);

    d3Vector out = getParticleCoordDecenteredPixel(particle_id, tomo_centre, tiltSeriesPixelSize);

    if (apply_origin_shifts)
    {
        const d3Matrix A_subtomogram = getSubtomogramMatrix(particle_id);

        out -= (A_subtomogram * getParticleOffset(particle_id)) / tiltSeriesPixelSize;
    }

    // SHWS & ABurt 19Jul2022: let's no longer do this in relion-4.1 (also 25nov22)
    // SHWS 17jan24: put '+1' back in because relion5 doesn't go as high resolution as tutorial data set as relion4....
    // SHWS 2feb24: we can't do this anymore now...
    //out.x += 1.0;
	//out.y += 1.0;
	//out.z += 1.0;


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
d4Matrix ParticleSet::getMatrix4x4(ParticleIndex particle_id, const gravis::d3Vector &tomo_centre, double w, double h, double d) const
{
        d3Matrix A = getMatrix3x3(particle_id);
        d3Vector pos = getPosition(particle_id, tomo_centre, true);

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
	if (!hasHalfSets())
        REPORT_ERROR("ERROR: function getHalfSet was called without having halfsets in the particle star file.");

    int s;
	partTable.getValue(EMDL_PARTICLE_RANDOM_SUBSET, s, particle_id.value);
	return s - 1;
}

bool ParticleSet::hasHalfSets() const
{
    return partTable.containsLabel(EMDL_PARTICLE_RANDOM_SUBSET);
}

void ParticleSet::write(const std::string& filename)
{
	std::ofstream ofs(filename);

    genTable.setName("general");
    genTable.setValue(EMDL_TOMO_SUBTOMOGRAM_STACK2D, is_stack2d);
    genTable.write(ofs);
	optTable.write(ofs);

    // Remove spurious uncentered coordinates in pixels from relion-4 files
    if (partTable.containsLabel(EMDL_IMAGE_COORD_X)) partTable.deactivateLabel(EMDL_IMAGE_COORD_X);
    if (partTable.containsLabel(EMDL_IMAGE_COORD_Y)) partTable.deactivateLabel(EMDL_IMAGE_COORD_Y);
    if (partTable.containsLabel(EMDL_IMAGE_COORD_Z)) partTable.deactivateLabel(EMDL_IMAGE_COORD_Z);

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
    if (weight != "")
	    partTable.setValue(EMDL_CTF_IMAGE, weight, particle_id.value);
}

d3Vector ParticleSet::getParticleOffset(ParticleIndex particle_id) const
{
	d3Vector out;
	
	if (partTable.containsLabel(EMDL_ORIENT_ORIGIN_X_ANGSTROM))
        partTable.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, out.x, particle_id.value);
	else
        out.x = 0.;
	if (partTable.containsLabel(EMDL_ORIENT_ORIGIN_Y_ANGSTROM))
        partTable.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, out.y, particle_id.value);
	else
        out.y = 0.;
    if (partTable.containsLabel(EMDL_ORIENT_ORIGIN_Z_ANGSTROM))
        partTable.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, out.z, particle_id.value);
    else
        out.z = 0.;

	return out;
}

void ParticleSet::setParticleOffset(ParticleIndex particle_id, const d3Vector& v)
{
	partTable.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, v.x, particle_id.value);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, v.y, particle_id.value);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, v.z, particle_id.value);
}

d3Vector ParticleSet::getParticleCoordDecenteredPixel(ParticleIndex particle_id, const gravis::d3Vector &tomo_centre, RFLOAT tiltSeriesPixelSize) const
{
	d3Vector out;

    if (partTable.containsLabel(EMDL_IMAGE_CENT_COORD_X_ANGST) &&
        partTable.containsLabel(EMDL_IMAGE_CENT_COORD_Y_ANGST) &&
        partTable.containsLabel(EMDL_IMAGE_CENT_COORD_Z_ANGST))
    {
        partTable.getValue(EMDL_IMAGE_CENT_COORD_X_ANGST, out.x, particle_id.value);
        partTable.getValue(EMDL_IMAGE_CENT_COORD_Y_ANGST, out.y, particle_id.value);
        partTable.getValue(EMDL_IMAGE_CENT_COORD_Z_ANGST, out.z, particle_id.value);

        // Inside Jasenko's code, all coordinates are in decentered pixels, convert now from centered Angstroms (centre_tomo is in pixels)
        out /= tiltSeriesPixelSize;
        out += tomo_centre;
    }
    else if (partTable.containsLabel(EMDL_IMAGE_COORD_X) &&
             partTable.containsLabel(EMDL_IMAGE_COORD_Y) &&
             partTable.containsLabel(EMDL_IMAGE_COORD_Z))
    {
        // Maintain backwards compatibility with relion-4
        partTable.getValue(EMDL_IMAGE_COORD_X, out.x, particle_id.value);
        partTable.getValue(EMDL_IMAGE_COORD_Y, out.y, particle_id.value);
        partTable.getValue(EMDL_IMAGE_COORD_Z, out.z, particle_id.value);

    }
    else
        REPORT_ERROR("Cannot find particle coordinates (rlnCenteredCoordinateX/Y/ZAngst) in particle star file");

	return out;
}

void ParticleSet::setParticleCoordDecenteredPixel(ParticleIndex particle_id, d3Vector v, const gravis::d3Vector &tomo_centre, RFLOAT tiltSeriesPixelSize)
{

    // Inside Jasenko's code all coordinates are in decentered pixels, convert now to centered Angstroms (centre_tomo is in pixels)
    v -= tomo_centre;
    v *= tiltSeriesPixelSize;

	partTable.setValue(EMDL_IMAGE_CENT_COORD_X_ANGST, v.x, particle_id.value);
	partTable.setValue(EMDL_IMAGE_CENT_COORD_Y_ANGST, v.y, particle_id.value);
	partTable.setValue(EMDL_IMAGE_CENT_COORD_Z_ANGST, v.z, particle_id.value);
}

int ParticleSet::getOpticsGroup(ParticleIndex particle_id) const
{
	if (!partTable.containsLabel(EMDL_IMAGE_OPTICS_GROUP))
	{
		REPORT_ERROR("ParticleSet::getOpticsGroup: optics group (rlnOpticsGroup) is missing from optics table");
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

double ParticleSet::getTiltSeriesPixelSize(int opticsGroup) const
{
	if (!optTable.containsLabel(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE))
	{
		REPORT_ERROR("ParticleSet::getTiltSeriesPixelSize: tilt series pixel size (rlnTomoTiltSeriesPixelSize) missing from optics table");
	}
	
	double out;
	optTable.getValueSafely(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, out, opticsGroup);
	return out;
}

std::vector<int> ParticleSet::getVisibleFrames(ParticleIndex particle_id) const
{
   	if (!partTable.containsLabel(EMDL_TOMO_VISIBLE_FRAMES))
    {
        REPORT_ERROR("ParticleSet::getVisibileFrames: frame visibility vector (rlnTomoVisibleFrames) missing from particles table");
    }

       std::vector<int> out;
       partTable.getValue(EMDL_TOMO_VISIBLE_FRAMES, out, particle_id.value);
       return out;

}

std::vector<d3Vector> ParticleSet::getTrajectoryInPixels(ParticleIndex particle_id, int fc, const gravis::d3Vector &tomo_centre, double pixelSize, bool from_original_coordinate) const
{
    const d3Vector p0 = getPosition(particle_id, tomo_centre, !from_original_coordinate);

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


