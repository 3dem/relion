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
#include "src/ctf.h"
#include <src/jaz/single_particle/obs_model.h>

/// Reserve large vectors with some reasonable estimate
// Larger numbers will still be OK, but memory management might suffer
#define MAX_NR_GROUPS 2000

////////////// Hierarchical metadata model

class ExpImage
{
public:
	// Position of the image in the original input STAR file
	long int id;

	// To which particle does this image belong
	long int particle_id;

	// This is the Nth image in this optics_group, for writing to scratch disk: filenames
	long int optics_group_id;

	// Name of this image (by this name it will be recognised upon reading)
	std::string name;

	// ID of the group that this image comes from
	long int group_id;

	// The optics group for this image
	int optics_group;

	// Pre-read array of the image in RAM
	MultidimArray<float> img;

	// Empty Constructor
	ExpImage() {}

	// Destructor needed for work with vectors
	~ExpImage() {}

	// Copy constructor needed for work with vectors
	ExpImage(ExpImage const& copy)
	{
		id = copy.id;
		particle_id = copy.particle_id;
		optics_group_id = copy.optics_group_id;
		name = copy.name;
		group_id = copy.group_id;
		optics_group = copy.optics_group;
		img = copy.img;

	}

	// Define assignment operator in terms of the copy constructor
	ExpImage& operator=(ExpImage const& copy)
	{
		id = copy.id;
		particle_id = copy.particle_id;
		optics_group_id = copy.optics_group_id;
		name = copy.name;
		group_id = copy.group_id;
		optics_group = copy.optics_group;
		img = copy.img;
		return *this;
	}
};

class ExpParticle
{
public:

	// Name of this particle (by this name all the images inside it will be grouped)
	std::string name;

	// Random subset this particle belongs to
	int random_subset;

	// Vector of all the images for this particle
	std::vector<ExpImage> images;

	// Empty Constructor
	ExpParticle() {}

	// Destructor needed for work with vectors
	~ExpParticle() {}

	// Copy constructor needed for work with vectors
	ExpParticle(ExpParticle const& copy)
	{
		name = copy.name;
		random_subset = copy.random_subset;
		images = copy.images;
	}

	// Define assignment operator in terms of the copy constructor
	ExpParticle& operator=(ExpParticle const& copy)
	{
		name = copy.name;
		random_subset = copy.random_subset;
		images = copy.images;
		return *this;
	}

	int numberOfImages()
	{
		return images.size();
	}
};

class ExpGroup
{
public:
	// ID of this group
	long int id;

	// The optics_group for this group
	int optics_group;

	// Name of this group (by this name it will be recognised upon reading)
	std::string name;

	// Empty Constructor
	ExpGroup() {}

	// Destructor needed for work with vectors
	~ExpGroup() {}

	// Copy constructor needed for work with vectors
	ExpGroup(ExpGroup const& copy)
	{
		id = copy.id;
		optics_group = copy.optics_group;
		name = copy.name;
	}

	// Define assignment operator in terms of the copy constructor
	ExpGroup& operator=(ExpGroup const& copy)
	{
		id = copy.id;
		optics_group = copy.optics_group;
		name = copy.name;
		return *this;
	}
};

class Experiment
{
public:
	// All groups in the experiment
	std::vector<ExpGroup> groups;

	// All particles in the experiment
	std::vector<ExpParticle> particles;

	// Indices of the sorted particles
	std::vector<long int> sorted_idx;

	// Number of particles in random subsets 1 and 2;
	long int nr_particles_subset1, nr_particles_subset2;

	// Number of images per optics group
	std::vector<long int> nr_images_per_optics_group;

	// One large MetaDataTable for all images
	MetaDataTable MDimg;

	// Number of bodies in multi-body refinement
	int nr_bodies;

	// Vector with MetaDataTables for orientations of different bodies in the multi-body refinement
	std::vector<MetaDataTable> MDbodies;

	// Observation model holding the data for all optics groups
	ObservationModel obsModel;

	// Directory on scratch disk to copy particles to
	FileName fn_scratch;

	// Number of particles saved on the scratchdir, one for each optics_group
	std::vector<long int> nr_parts_on_scratch;

	// Number of Gb on scratch disk before copying particles
	RFLOAT free_space_Gb;

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
		groups.reserve(MAX_NR_GROUPS);
		particles.clear(); // reserve upon reading
		sorted_idx.clear();
		nr_particles_subset1 = nr_particles_subset2 = 0;
		nr_bodies = 1;
		fn_scratch = "";
		nr_parts_on_scratch.clear();
		free_space_Gb = 10;
		is_3D = false;
		MDimg.clear();
		MDimg.setIsList(false);
		MDbodies.clear();
		MDimg.setName("images");
	}

	// Calculate the total number of particles in this experiment
	long int numberOfParticles(int random_subset = 0);

	// Get the total number of images in a given particle
	long int numberOfImagesInParticle(long int part_id);

	// Calculate the total number of groups in this experiment
	long int numberOfGroups();

	// Calculate the total number of optics groups in this experiment
	int numberOfOpticsGroups();

	// Is any of the optics groups CTF-premultiplied?
	bool hasCtfPremultiplied();

	// Get the pixel size for this optics group
	RFLOAT getOpticsPixelSize(int optics_group);

	// Get the original image size for this optics group
	int getOpticsImageSize(int optics_group);

	// Get the random_subset for this particle
	int getRandomSubset(long int part_id);

	// Get the group_id for the N'th image for this particle
	long int getGroupId(long int part_id, int img_id);

	// Get the optics group to which the N'th image for this particle belongs
	int getOpticsGroup(long int part_id, int img_id);

	// Get the original position in the input STAR file for the N'th image for this particle
	int getOriginalImageId(long int part_id, int img_id);

	// Get the pixel size for the N-th image of this particle
	RFLOAT getImagePixelSize(long int part_id, int img_id);

	// Get the vector of number of images per group_id
	void getNumberOfImagesPerGroup(std::vector<long int> &nr_particles_per_group, int random_subset = 0);

	// Get the vector of number of images per group_id
	void getNumberOfImagesPerOpticsGroup(std::vector<long int> &nr_particles_per_group, int random_subset = 0);

	// Get the metadata-row for this image in a separate MetaDataTable
	MetaDataTable getMetaDataImage(long int part_id, int img_id);

	// Which micrograph (or tomogram) doe this particle image comes from?
	FileName getMicrographName(long int ori_img_id);
	FileName getMicrographName(long int part_id, int img_id);

	// Add a particle
	long int addParticle(std::string part_name, int random_subset = 0);

 	// Add an image to the given particle
	int addImageToParticle(long int part_id, std::string img_name, long int ori_img_id, long int group_id,
	                       int optics_group, bool unique);

	// Add a group
	long int addGroup(std::string mic_name, int optics_group);

	// for separate refinement of random halves of the data
	void divideParticlesInRandomHalves(int seed, bool do_helical_refine = false);

	// Randomise the order of the particles
	void randomiseParticlesOrder(int seed, bool do_split_random_halves = false, int subsets_size = -1);

	// Make sure the images inside each particle are in the right order
	void orderImagesInParticles();

	// Add a given number of new bodies (for multi-body refinement) to the Experiment,
	// by copying the relevant entries from MDimg into MDbodies
	void initialiseBodies(int _nr_bodies);

	// Get the image name for a given part_id
	bool getImageNameOnScratch(long int part_id, int img_id, FileName &fn_img, bool is_ctf_image = false);

	// For parallel executions, lock the scratch directory with a unique code, so we won't copy the same data many times to the same position
	// This determines the lockname and removes the lock if it exists
	FileName initialiseScratchLock(FileName _fn_scratch, FileName _fn_out);

	// Returns true if particles need to be copied, and creates a lock file.
	// Returns false if the particles do not need to be copied. In that case, only the number of particles on the scratch disk needs to be counted
	// Also checks how much free space there is on the scratch dir
	bool prepareScratchDirectory(FileName _fn_scratch, FileName fn_lock = "");

	void setScratchDirectory(FileName _fn_scratch, bool do_reuse_scratch, int verb=0);

	// Wipe the generic scratch directory clean
	void deleteDataOnScratch();

	// Copy particles from their original position to a scratch directory
	// Monitor when the scratch disk gets to have fewer than free_scratch_Gb space,
	// in that case, stop copying, and keep reading particles from where they were...
	void copyParticlesToScratch(int verb, bool do_copy = true, bool also_do_ctf_image = false, RFLOAT free_scratch_Gb = 10);

	// Read from file
	void read(
		FileName fn_in,
		bool do_ignore_particle_name = false,
		bool do_ignore_group_name = false, bool do_preread_images = false,
		bool need_tiltpsipriors_for_helical_refine = false, int verb = 0);

	// Write
	void write(FileName fn_root);


private:

	struct compareOpticsGroupsParticles
	{
	    const std::vector<ExpParticle>& particles;
	    compareOpticsGroupsParticles(const std::vector<ExpParticle>& particles) : particles(particles) { }
	    bool operator()(const long int i, const long int j) { return particles[i].images[0].optics_group < particles[j].images[0].optics_group;}
	};

	struct compareRandomSubsetParticles
	{
	    const std::vector<ExpParticle>& particles;
	    compareRandomSubsetParticles(const std::vector<ExpParticle>& particles) : particles(particles) { }
	    bool operator()(const long int i, const long int j) { return particles[i].random_subset < particles[j].random_subset;}
	};


};

#endif /* METADATA_MODEL_H_ */
