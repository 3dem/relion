/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
#ifndef EXP_MODEL_H_
#define EXP_MODEL_H_
#include <fstream>
#include "src/matrix2d.h"
#include "src/image.h"
#include "src/multidim_array.h"
#include "src/metadata_table.h"
#include "src/time.h"

/// Reserve large vectors with some reasonable estimate
// Larger numbers will still be OK, but memory management might suffer
#define MAX_NR_PARTICLES_PER_MICROGRAPH 1000
#define MAX_NR_MICROGRAPHS 2000
#define MAX_NR_FRAMES_PER_MOVIE 100

////////////// Hierarchical metadata model for tilt series

class ExpParticle
{
public:
	// Particle id
	long int id;

	// ID of the micrograph that this particle comes from
	long int micrograph_id;

	// ID of the group that this particle comes from
	long int group_id;

	// Random subset this particle belongs to
	int random_subset;

	// Pre-read array of the image in RAM
	MultidimArray<float> img;

	// Empty Constructor
	ExpParticle()
	{
		clear();
	}

	// Destructor needed for work with vectors
	~ExpParticle()
	{
		clear();
	}

	// Initialise
	void clear()
	{
		id = micrograph_id = group_id = -1;
		random_subset = 0;
		img.clear();
	}

};

class ExpOriginalParticle
{
public:
	// Name of this particle (by this name it will be recognised upon reading)
	std::string name;

	// Random subset this original_particle belongs to
	int random_subset;

	// All the id's of the particles that were derived from this original particle
	std::vector<long int> particles_id;

	// Order of those particles in the original particle (extracted from mic_name)
	std::vector<int> particles_order;

	// Empty Constructor
	ExpOriginalParticle()
	{
		clear();
	}

	// Destructor needed for work with vectors
	~ExpOriginalParticle()
	{
		clear();
	}

	// Initialise
	void clear()
	{
		name="undefined";
		particles_id.clear();
		particles_order.clear();
		particles_id.reserve(MAX_NR_FRAMES_PER_MOVIE);
		particles_order.reserve(MAX_NR_FRAMES_PER_MOVIE);
	}

	void addParticle(long int _particle_id, int _random_subset, int _order);

};


class ExpMicrograph
{
public:
	// ID of this micrograph, i.e. which number in the MDmic am I?
	long int id;

	// Name of this micrograph (by this name it will be recognised upon reading)
	std::string name;

	// All the particles that were recorded on this micrograph
	std::vector<long int> particle_ids;

	// All the original particles that were recorded on this average micrograph
	std::vector<long int> ori_particle_ids;

	// Empty Constructor
	ExpMicrograph()
	{
		clear();
	}

	// Destructor needed for work with vectors
	~ExpMicrograph()
	{
		clear();
	}

	// Copy constructor needed for work with vectors
	ExpMicrograph(ExpMicrograph const& copy)
	{
		id = copy.id;
		name = copy.name;
		particle_ids = copy.particle_ids;
		ori_particle_ids = copy.ori_particle_ids;

	}

	// Define assignment operator in terms of the copy constructor
	ExpMicrograph& operator=(ExpMicrograph const& copy)
	{
		id = copy.id;
		name = copy.name;
		particle_ids = copy.particle_ids;
		ori_particle_ids = copy.ori_particle_ids;
		return *this;
	}

	// Initialise
	void clear()
	{
		id = -1;
		name="";
		particle_ids.clear();
		particle_ids.reserve(MAX_NR_PARTICLES_PER_MICROGRAPH);
		ori_particle_ids.clear();
	}

};

class ExpGroup
{
public:
	// ID of this group
	long int id;

	// Name of this group (by this name it will be recognised upon reading)
	std::string name;

	// Empty Constructor
	ExpGroup()
	{
		clear();
	}

	// Destructor needed for work with vectors
	~ExpGroup()
	{
		clear();
	}

	// Copy constructor needed for work with vectors
	ExpGroup(ExpGroup const& copy)
	{
		id = copy.id;
		name = copy.name;
	}

	// Define assignment operator in terms of the copy constructor
	ExpGroup& operator=(ExpGroup const& copy)
	{
		id = copy.id;
		name = copy.name;
		return *this;
	}

	// Initialise
	void clear()
	{
		id = -1;
		name="";
	}

};


class Experiment
{
public:
	// All groups in the experiment
	std::vector<ExpGroup> groups;

	// All micrographs in the experiment
	std::vector<ExpMicrograph> micrographs;

	// All particles in the experiment
	std::vector<ExpParticle> particles;

	// All original particles in the experiment
	std::vector<ExpOriginalParticle> ori_particles;

	// Number of particles in random subsets 1 and 2;
	long int nr_ori_particles_subset1, nr_ori_particles_subset2;

	// Experiment-related metadata
    MetaDataTable MDexp;

    // One large MetaDataTable for all images
    MetaDataTable MDimg;

    // Number of bodies in multi-body refinement
    int nr_bodies;

    // Vector with MetaDataTables for orientations of different bodies in the multi-body refinement
    std::vector<MetaDataTable> MDbodies;

    // One large MetaDataTable for all micrographs
    MetaDataTable MDmic;

    // Directory on scratch disk to copy particles to
    FileName fn_scratch;

    // Number of particles saved on the scratchdir
    long int nr_parts_on_scratch;

    // Number of Gb on scratch disk before copying particles
    long int free_space_Gb;

    // Is this sub-tomograms?
    bool is_3D;

	// Empty Constructor
	Experiment()
	{
		clear();
	}

	~Experiment()
	{
		clear();
	}

	void clear()
	{
		groups.clear();
		groups.reserve(MAX_NR_MICROGRAPHS);
		micrographs.clear();
		micrographs.reserve(MAX_NR_MICROGRAPHS);
		particles.clear(); // reserve upon reading
		ori_particles.clear(); // TODO: reserve upon reading
		nr_ori_particles_subset1 = nr_ori_particles_subset2 = 0;
		nr_bodies = 1;
		fn_scratch = "";
		nr_parts_on_scratch = 0;
		free_space_Gb = 10;
		is_3D = false;
		MDexp.clear();
		MDexp.setIsList(true);
		MDimg.clear();
		MDimg.setIsList(false);
		MDbodies.clear();
		MDmic.clear();
		MDmic.setIsList(false);
		MDimg.setName("images");
		MDmic.setName("micrographs");
		MDexp.setName("experiment");
	}

	// Calculate the total number of particles in this experiment
	long int numberOfParticles(int random_subset = 0);

	// Calculate the total number of particles in this experiment
	long int numberOfOriginalParticles(int random_subset = 0);

	// Calculate the total number of micrographs in this experiment
	long int numberOfMicrographs();

	// Calculate the total number of groups in this experiment
	long int numberOfGroups();

	// Get the random_subset for this particle
	int getRandomSubset(long int part_id);

	// Get the micrograph_id for the N'th image for this particle
	long int getMicrographId(long int part_id);

	// Get the group_id for the N'th image for this particle
	long int getGroupId(long int part_id);

	// Get the metadata-row for this image in a separate MetaDataTable
	MetaDataTable getMetaDataImage(long int part_id);

	// Add a particle
	long int addParticle(long int group_id, long int micrograph_id, int random_subset = 0);

	// Add an original particle
	long int addOriginalParticle(std::string part_name, int random_subset = 0);

	// Add a group
	long int addGroup(std::string mic_name);

	// Add a micrograph
	long int addMicrograph(std::string mic_name);

	// for separate refinement of random halves of the data
	void divideOriginalParticlesInRandomHalves(int seed, bool do_helical_refine = false);

	// Randomise the order of the original_particles
	void randomiseOriginalParticlesOrder(int seed, bool do_split_random_halves = false);

	// calculate maximum number of images for a particle (possibly within a range of particles)
	int maxNumberOfImagesPerOriginalParticle(long int first_particle_id = -1, long int last_particle_id = -1);

	// Make sure the particles inside each orriginal_particle are in the right order
	// After they have been ordered, get rid of the particles_order vector inside the ori_particles
	void orderParticlesInOriginalParticles();

	// Add a given number of new bodies (for multi-body refinement) to the Experiment,
	// by copying the relevant entries from MDimg into MDbodies
	void initialiseBodies(int _nr_bodies);

	// Get the image name for a given particle_id
	bool getImageNameOnScratch(long int particle_id, FileName &fn_img, bool is_ctf_image = false);

	// For parallel executions, lock the scratch directory with a unique code, so we won't copy the same data many times to the same position
	// This determines the lockname and removes the lock if it exists
	FileName initialiseScratchLock(FileName _fn_scratch, FileName _fn_out);

	// Returns true if particles need to be copied, and creates a lock file.
	// Returns false if the particles do not need to be copied. In that case, only the number of particles on the scratch disk needs to be counted
	// Also checks how much free space there is on the scratch dir
	bool prepareScratchDirectory(FileName _fn_scratch, FileName fn_lock = "");

	// Wipe the generic scratch directory clean
	void deleteDataOnScratch();

	// Copy particles from their original position to a scratch directory
	// Monitor when the scratch disk gets to have fewer than free_scratch_Gb space,
	// in that case, stop copying, and keep reading particles from where they were...
	void copyParticlesToScratch(int verb, bool do_copy = true, bool also_do_ctf_image = false, long int free_scratch_Gb = 10);


	// Print help message for possible command-line options
	void usage();

	// Read from file
	void read(FileName fn_in, bool do_ignore_original_particle_name = false,
			bool do_ignore_group_name = false, bool do_preread_images = false,
			bool need_tiltpsipriors_for_helical_refine = false);

	// Write
	void write(FileName fn_root);



};

#endif /* METADATA_MODEL_H_ */
