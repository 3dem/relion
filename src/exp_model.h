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
	MultidimArray<RFLOAT> img;

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

// This class describes which OriginalParticles in the data set belong to the same frame-average micrograph
class AverageMicrograph
{
public:
	// ID of this average micrograph, i.e. which number in the MDmic am I?
	long int id;

	// Name of this average micrograph
	std::string name;

	// All the original particles that were recorded on this average micrograph
	std::vector<long int> ori_particles_id;

	// Empty Constructor
	AverageMicrograph()
	{
		clear();
	}

	// Destructor needed for work with vectors
	~AverageMicrograph()
	{
		clear();
	}

	// Copy constructor needed for work with vectors
	AverageMicrograph(AverageMicrograph const& copy)
	{
		id = copy.id;
		name = copy.name;
		ori_particles_id = copy.ori_particles_id;
	}

	// Define assignment operator in terms of the copy constructor
	AverageMicrograph& operator=(AverageMicrograph const& copy)
	{
		id = copy.id;
		name = copy.name;
		ori_particles_id = copy.ori_particles_id;
		return *this;
	}

	// Initialise
	void clear()
	{
		id = -1;
		name="";
		ori_particles_id.clear();
		ori_particles_id.reserve(MAX_NR_PARTICLES_PER_MICROGRAPH);
	}
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

	}

	// Define assignment operator in terms of the copy constructor
	ExpMicrograph& operator=(ExpMicrograph const& copy)
	{
		id = copy.id;
		name = copy.name;
		particle_ids = copy.particle_ids;
		return *this;
	}

	// Initialise
	void clear()
	{
		id = -1;
		name="";
		particle_ids.clear();
		particle_ids.reserve(MAX_NR_PARTICLES_PER_MICROGRAPH);
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

	// All average micrographs in this experiment (only used for movie-processing, i.e. by the particle_polisher
	std::vector<AverageMicrograph> average_micrographs;

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
		micrographs.clear();
		micrographs.reserve(MAX_NR_MICROGRAPHS);
		groups.reserve(MAX_NR_MICROGRAPHS);
		particles.clear(); // reserve upon reading
		ori_particles.clear(); // TODO: reserve upon reading
		MDexp.clear();
		MDexp.setIsList(true);
		MDimg.clear();
		MDimg.setIsList(false);
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

	// Add an AverageMicrograph
	long int addAverageMicrograph(std::string avg_mic_name);

	// for separate refinement of random halves of the data
	void divideOriginalParticlesInRandomHalves(int seed);

	// Randomise the order of the original_particles
	void randomiseOriginalParticlesOrder(int seed, bool do_split_random_halves = false);

	// calculate maximum number of images for a particle (possibly within a range of particles)
	int maxNumberOfImagesPerOriginalParticle(long int first_particle_id = -1, long int last_particle_id = -1);

	// Given the STAR file of a set of movieframes, expand the current Experiment to contain all movie frames
	// rlnParticleName entries in the movie-frame Experiment should coincide with rlnImageName entries in the current Experiment
	// the entries rlnAngleRot, rlnAngleTilt, rlnAnglePsi, rlnOriginX and rlnOriginY will be taken from the current Experiment and
	// copied into the new moevieframe Experiment. In addition, these values will be used to center the corresponding Priors
	void expandToMovieFrames(FileName fn_data_movie, int verb = 0);

	// Make sure the particles inside each orriginal_particle are in the right order
	// After they have been ordered, get rid of the particles_order vector inside the ori_particles
	void orderParticlesInOriginalParticles();

	// Add a given number of new bodies (for multi-body refinement) to the Experiment,
	// by copying the relevant entries from MDimg into MDbodies
	void initialiseBodies(int _nr_bodies);

	// Print help message for possible command-line options
	void usage();

	// Read from file
	// May23,2015 - Shaoda, Helical refinement
	void read(FileName fn_in, bool do_ignore_original_particle_name = false, bool do_ignore_group_name = false, bool do_preread_images = false, bool need_tiltpsipriors_for_helical_refine = false);

	// Write
	void write(FileName fn_root);

};

#endif /* METADATA_MODEL_H_ */
