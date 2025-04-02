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
#include "src/exp_model.h"
#include <sys/statvfs.h>
using namespace gravis;

long int Experiment::numberOfParticles(int random_subset)
{
	if (random_subset == 0)
		return particles.size();
	else if (random_subset == 1)
		return nr_particles_subset1;
	else if (random_subset == 2)
		return nr_particles_subset2;
	else
		REPORT_ERROR("ERROR: Experiment::numberOfParticles invalid random_subset: " + integerToString(random_subset));
}

// Get the total number of images in a given particle
long int Experiment::numberOfImagesInParticle(long int part_id)
{
	return particles[part_id].images.size();
}

long int Experiment::numberOfGroups()
{
	return groups.size();
}

int Experiment::numberOfOpticsGroups()
{
	return obsModel.numberOfOpticsGroups();
}

bool Experiment::hasCtfPremultiplied()
{
	for (int og = 0; og < numberOfOpticsGroups(); og++)
		if (obsModel.getCtfPremultiplied(og)) return true;

	return false;
}

bool Experiment::hasCtfCorrected()
{
    for (int og = 0; og < numberOfOpticsGroups(); og++)
        if (obsModel.getCtfCorrected(og)) return true;

    return false;
}

RFLOAT Experiment::getOpticsPixelSize(int optics_group)
{
	return obsModel.getPixelSize(optics_group);
}

int Experiment::getOpticsImageSize(int optics_group)
{
	return obsModel.getBoxSize(optics_group);
}

long int Experiment::getGroupId(long int part_id)
{
	return particles[part_id].group_id;
}

int Experiment::getOpticsGroup(long part_id)
{
	return particles[part_id].optics_group;
}

int Experiment::getRandomSubset(long int part_id)
{
	return particles[part_id].random_subset;
}

RFLOAT Experiment::getImagePixelSize(long int part_id)
{
	int optics_group = particles[part_id].optics_group;
	return obsModel.getPixelSize(optics_group);
}

Matrix2D<RFLOAT> Experiment::getRotationMatrix(long int part_id, int img_id)
{
    return particles[part_id].images[img_id].Aproj;
}

void Experiment::getTranslationInTiltSeries(long int part_id, int img_id,
                                                        RFLOAT shift3d_x, RFLOAT shift3d_y, RFLOAT shift3d_z,
                                                        RFLOAT &shift2d_x, RFLOAT &shift2d_y, RFLOAT &shift2d_z)
{
    Matrix2D<RFLOAT> Aproj = particles[part_id].images[img_id].Aproj;
    shift2d_x = Aproj(0,0) * shift3d_x + Aproj(0,1) * shift3d_y + Aproj(0,2) * shift3d_z;
    shift2d_y = Aproj(1,0) * shift3d_x + Aproj(1,1) * shift3d_y + Aproj(1,2) * shift3d_z;
    shift2d_z = 0.;
}

void Experiment::getNumberOfParticlesPerGroup(std::vector<long int> &nr_particles_per_group, int random_subset)
{
	nr_particles_per_group.resize(groups.size());
	for (long int i = 0; i < nr_particles_per_group.size(); i++)
		nr_particles_per_group[i] = 0.;

	for (long int part_id = 0; part_id < particles.size(); part_id++)
	{
		if (random_subset == 0 || particles[part_id].random_subset == random_subset)
		{
			nr_particles_per_group[particles[part_id].group_id] += 1;
		}
	}

}

void Experiment::getNumberOfParticlesPerOpticsGroup(std::vector<long int> &nr_particles_per_optics_group, int random_subset)
{
	nr_particles_per_optics_group.resize(obsModel.numberOfOpticsGroups());
	for (long int i = 0; i < nr_particles_per_optics_group.size(); i++)
		nr_particles_per_optics_group[i] = 0.;

	for (long int part_id = 0; part_id < particles.size(); part_id++)
	{
		if (random_subset == 0 || particles[part_id].random_subset == random_subset)
		{
			nr_particles_per_optics_group[particles[part_id].optics_group] += 1;
		}
	}
}

MetaDataTable Experiment::getMetaDataParticle(long int part_id)
{
	MetaDataTable result;
	result.addObject(MDimg.getObject(part_id));
	return result;
}

FileName Experiment::getMicrographName(long int part_id)
{
	FileName micname="";
	if (is_3D || is_tomo)
	{
		MDimg.getValue(EMDL_TOMO_NAME, micname, part_id);
	}
	else
	{
		MDimg.getValue(EMDL_MICROGRAPH_NAME, micname, part_id);
	}

	return micname;
}

void Experiment::addParticle(std::string img_name, int optics_group, long int group_id, int random_subset, int tomogram_id)
{

    if (optics_group >= obsModel.numberOfOpticsGroups())
        REPORT_ERROR("Experiment::addImageToParticle: optics_group out of range");

    if (group_id >= groups.size())
        REPORT_ERROR("Experiment::addImageToParticle: group_id out of range");

    ExpParticle particle;
	particle.id = particles.size();
    particle.name = img_name;
    particle.tomogram_id = tomogram_id;
	particle.random_subset = random_subset;
    particle.optics_group = optics_group;
    particle.group_id = group_id;

    nr_particles_per_optics_group[optics_group]++;
    particle.optics_group_id = nr_particles_per_optics_group[optics_group] - 1;

	// Push back this particle in the particles vector and its sorted index in sorted_idx
	sorted_idx.push_back(particles.size());
	particles.push_back(particle);

    return;
}

void Experiment::addImageToParticle(long int part_id, d4Matrix *Aproj, CTF *ctf, float dose, double BfactorPerElectronDose)
{

    Matrix2D<RFLOAT> A(3,3);
    if (Aproj == NULL)
    {
        A.initIdentity();
    }
    else
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                A(i, j) = (*Aproj)(i, j);
    }

	ExpImage img;
    if (ctf == NULL)
    {
        img.defU = img.defV = img.defAngle = img.phase_shift = 0.;
        img.scale = 1.;
    }
    else
    {
        img.defU = ctf->DeltafU;
        img.defV = ctf->DeltafV;
        img.defAngle = ctf->azimuthal_angle;
        img.scale = ctf->scale;
        img.phase_shift = ctf->phase_shift;
    }

	img.particle_id = part_id;
	img.Aproj = A;

    if (BfactorPerElectronDose > 0.)
    {
        img.dose = -999.;
        img.bfactor = BfactorPerElectronDose * dose;
    }
    else
    {
        img.dose = dose;
        img.bfactor = 0.;
    }

	// Push back this particle in the particles vector
	particles[part_id].images.push_back(img);

	return;
}

long int Experiment::addGroup(std::string group_name, int _optics_group)
{
	// Add new group to this Experiment
	ExpGroup group;
	group.id = groups.size(); // start counting groups at 0!
	group.optics_group = _optics_group;
	group.name = group_name;

	// Push back this group
	groups.push_back(group);

	// Return the id in the groups vector
	return group.id;
}

void Experiment::divideParticlesInRandomHalves(int seed, bool do_helical_refine)
{
	// Only do this if the random_subset of all original_particles is zero
	bool all_are_zero = true;
	bool some_are_zero = false;
	nr_particles_subset1 = 0;
	nr_particles_subset2 = 0;
	for (long int i = 0; i < particles.size(); i++)
	{
		int random_subset = particles[i].random_subset;
		if (random_subset != 0)
		{
			all_are_zero = false;
			// Keep track of how many particles there are in each subset
			if (random_subset == 1)
				nr_particles_subset1++;
			else if (random_subset == 2)
				nr_particles_subset2++;
			else
				REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: invalid number for random subset (i.e. not 1 or 2): " + integerToString(random_subset));
		}
		else
			some_are_zero = true;

		if (!all_are_zero && some_are_zero)
			REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: some random subset values are zero and others are not. They should all be zero, or all bigger than zero!");
	}

	if (all_are_zero)
	{
		// Only randomise them if the random_subset values were not read in from the STAR file
		srand(seed);
		if (do_helical_refine)
		{
			std::string mic_name, img_name;
			int nr_swaps, nr_segments_subset1, nr_segments_subset2, helical_tube_id;
			std::map<std::string, int> map_mics;
			std::map<std::string, int>::const_iterator ii_map;
			std::vector<std::pair<std::string, int> > vec_mics;

			bool divide_according_to_helical_tube_id = false;
			if (MDimg.containsLabel(EMDL_PARTICLE_HELICAL_TUBE_ID))
				divide_according_to_helical_tube_id = true;

			// Count micrograph names
			map_mics.clear();
			for (long int part_id = 0; part_id < particles.size(); part_id++)
			{
				// Get name of micrograph of the first image in this particle
				mic_name = getMicrographName(part_id);
				if (divide_according_to_helical_tube_id)
				{
					MDimg.getValue(EMDL_PARTICLE_HELICAL_TUBE_ID, helical_tube_id, part_id);
					if (helical_tube_id < 1)
						REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: Helical tube ID should be positive integer!");
					mic_name += std::string("_TUBEID_");
					mic_name += std::string(integerToString(helical_tube_id));
				}
				if ((map_mics.insert(std::make_pair(mic_name, 1))).second == false)
					map_mics[mic_name]++;
			}

			vec_mics.clear();
			for (ii_map = map_mics.begin(); ii_map != map_mics.end(); ii_map++)
				vec_mics.push_back(*ii_map);

			// NEW RANDOMISATION (better than the old one)
			nr_swaps = 0;
			for (int ptr_a = 0; ptr_a < (vec_mics.size() - 1); ptr_a++)
			{
				std::pair<std::string, int> tmp;
				int ptr_b = ROUND(rnd_unif(ptr_a, vec_mics.size() - 1));
				if ( (ptr_b <= ptr_a) || (ptr_b >= vec_mics.size()) )
					continue;
				nr_swaps++;
				tmp = vec_mics[ptr_a];
				vec_mics[ptr_a] = vec_mics[ptr_b];
				vec_mics[ptr_b] = tmp;
			}

			// Divide micrographs into halves
			map_mics.clear();
			nr_segments_subset1 = nr_segments_subset2 = 0;
			for (int ii = 0; ii < vec_mics.size(); ii++)
			{
				if (nr_segments_subset1 < nr_segments_subset2)
				{
					nr_segments_subset1 += vec_mics[ii].second;
					vec_mics[ii].second = 1;
				}
				else
				{
					nr_segments_subset2 += vec_mics[ii].second;
					vec_mics[ii].second = 2;
				}
				map_mics.insert(vec_mics[ii]);
			}

			for (long int part_id = 0; part_id < particles.size(); part_id++)
			{
				// Get name of micrograph of the first image in this particle
				mic_name = getMicrographName(part_id);
				if (divide_according_to_helical_tube_id)
				{
					MDimg.getValue(EMDL_PARTICLE_HELICAL_TUBE_ID, helical_tube_id, part_id);
					if (helical_tube_id < 1)
						REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: Helical tube ID should be positive integer!");
					mic_name += std::string("_TUBEID_");
					mic_name += std::string(integerToString(helical_tube_id));
				}
				particles[part_id].random_subset = map_mics[mic_name];
			}
		}
		else
		{
			for (long int part_id = 0; part_id < particles.size(); part_id++)
			{
				int random_subset = rand() % 2 + 1;
				particles[part_id].random_subset = random_subset; // randomly 1 or 2
			}
		}

		// Now that random subsets have been assigned, count the number of particles in each subset and set new labels in entire MDimg
		for (long int part_id = 0; part_id < particles.size(); part_id++)
		{
			int random_subset = getRandomSubset(part_id);

			if (random_subset == 1)
				nr_particles_subset1++;
			else if (random_subset == 2)
				nr_particles_subset2++;
			else
				REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: invalid number for random subset (i.e. not 1 or 2): " + integerToString(random_subset));

            MDimg.setValue(EMDL_PARTICLE_RANDOM_SUBSET, random_subset, part_id);
		}
	}

	if (nr_particles_subset2 == 0 || nr_particles_subset1 == 0)
		REPORT_ERROR("ERROR: one of your half sets has no segments. Is rlnRandomSubset set to 1 or 2 in your particles STAR file? Or in case you're doing helical, half-sets are always per-filament, so provide at least 2 filaments.");

	std::stable_sort(sorted_idx.begin(), sorted_idx.end(), compareRandomSubsetParticles(particles));

}

void Experiment::randomiseParticlesOrder(int seed, bool do_split_random_halves, int subsets_size)
{
	//This static flag is for only randomize once
	static bool randomised = false;
	const bool doing_subset = 0 < subsets_size && subsets_size < numberOfParticles() ;
	if (!randomised || doing_subset)
	{
		srand(seed);

		if (do_split_random_halves)
		{
			std::stable_sort(sorted_idx.begin(), sorted_idx.end(), compareRandomSubsetParticles(particles));

			// sanity check
			long int nr_half1 = 0, nr_half2 = 0;
			for (long int i = 0; i < particles.size(); i++)
			{
				const int random_subset = particles[i].random_subset;
				if (random_subset == 1)
					nr_half1++;
				else if (random_subset == 2)
					nr_half2++;
				else
					REPORT_ERROR("ERROR Experiment::randomiseParticlesOrder: invalid number for random subset (i.e. not 1 or 2): " + integerToString(random_subset));
			}

			if (nr_half1 != nr_particles_subset1)
				REPORT_ERROR("ERROR Experiment::randomiseParticlesOrder: invalid half1 size:" + integerToString(nr_half1) + " != " + integerToString(nr_particles_subset1));
			if (nr_half2 != nr_particles_subset2)
				REPORT_ERROR("ERROR Experiment::randomiseParticlesOrder: invalid half2 size:" + integerToString(nr_half2) + " != " + integerToString(nr_particles_subset2));

			// Randomise the two particle lists
			std::random_shuffle(sorted_idx.begin(), sorted_idx.begin() + nr_half1);
			std::random_shuffle(sorted_idx.begin() + nr_half1, sorted_idx.end());

			// Make sure the particles are sorted on their optics_group.
			// Otherwise CudaFFT re-calculation of plans every time image size changes slows down things a lot!
			long max_nr1 = doing_subset ? subsets_size : nr_half1;
			long max_nr2 = doing_subset ? subsets_size : nr_half2;
			std::stable_sort(sorted_idx.begin(), sorted_idx.begin() + max_nr1, compareOpticsGroupsParticles(particles));
			std::stable_sort(sorted_idx.begin() + nr_half1, sorted_idx.begin() + nr_half1 + max_nr2, compareOpticsGroupsParticles(particles));
		}
		else
		{
			// Just randomise the entire vector
			std::random_shuffle(sorted_idx.begin(), sorted_idx.end());

			// Make sure the particles are sorted on their optics_group.
			// Otherwise CudaFFT re-calculation of plans every time image size changes slows down things a lot!
			long max_nr = doing_subset ? subsets_size : numberOfParticles();
 			std::stable_sort(sorted_idx.begin(), sorted_idx.begin() + max_nr, compareOpticsGroupsParticles(particles));
		}

		randomised = true;
	}
}

void Experiment::initialiseBodies(int _nr_bodies)
{
	if (_nr_bodies < 2)
	{
		return;
	}
	else
	{
		nr_bodies = _nr_bodies;
		MetaDataTable MDbody;
		MDbody.setIsList(false);
		bool is_3d = (MDimg.containsLabel(EMDL_ORIENT_ORIGIN_Z));
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDimg)
		{
			MDbody.addObject();
			RFLOAT norm, zero=0., ninety=90.;
			MDimg.getValue(EMDL_IMAGE_NORM_CORRECTION, norm);
			MDbody.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, zero);
			MDbody.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, zero);
			MDbody.setValue(EMDL_ORIENT_ROT, zero);
			MDbody.setValue(EMDL_ORIENT_TILT, ninety);
			MDbody.setValue(EMDL_ORIENT_PSI, zero);
			MDbody.setValue(EMDL_IMAGE_NORM_CORRECTION, norm);
			if (is_3d)
			{
				MDbody.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, zero);
			}
		}
		// Now just fill all bodies with that MDbody
		MDbodies.resize(nr_bodies, MDbody);
		for (int ibody = 0; ibody < nr_bodies; ibody++)
		{
			std::string tablename = "images_body_" + integerToString(ibody+1);
			MDbodies[ibody].setName(tablename);
		}
	}
}

bool Experiment::getImageNameOnScratch(long int part_id, FileName &fn_img, bool is_ctf_image)
{
	int optics_group = getOpticsGroup(part_id);
    //TODO: move optics_group_id to particle, not image!!!!
	long int my_id = particles[part_id].optics_group_id;

    //REPORT_ERROR("DEBUG: STILL NEED TO ACCOUNT FOR MULTIPLE IMAGES PER PARTICLE HEREE!!!!!! UNFINISHED CODE....");

#ifdef DEBUG_SCRATCH
	std::cerr << "part_id = " << part_id << " img_id = " << img_id << " my_id = " << my_id << " nr_parts_on_scratch[" << optics_group << "] = " << nr_parts_on_scratch[optics_group] << std::endl;
#endif

	if (fn_scratch != "" && my_id < nr_parts_on_scratch[optics_group])
	{
		if (is_3D)
		{
			if (is_ctf_image)
				fn_img = fn_scratch + "opticsgroup" + integerToString(optics_group+1) + "_particle_ctf" + integerToString(my_id+1)+".mrc";
			else
				fn_img = fn_scratch + "opticsgroup" + integerToString(optics_group+1) + "_particle" + integerToString(my_id+1)+".mrc";
		}
        else if (is_tomo)
        {
            fn_img = fn_scratch + "opticsgroup" + integerToString(optics_group+1) + "_particle" + integerToString(my_id+1)+".mrcs";
        }
		else
		{
			// Write different optics groups into different stacks, as sizes might be different
			FileName fn_tmp = fn_scratch + "opticsgroup" + integerToString(optics_group+1) + "_particles.mrcs";
			fn_img.compose(my_id+1, fn_tmp);
		}

#ifdef DEBUG_SCRATCH
		std::cerr << "getImageNameOnScratch: " << particles[part_id].name << " is cached at " << fn_img << std::endl;
#endif
		return true;
	}
	else
	{
		return false;
	}
}

void Experiment::setScratchDirectory(FileName _fn_scratch, bool do_reuse_scratch, int verb)
{
	// Make sure fn_scratch ends with a slash
	if (_fn_scratch[_fn_scratch.length()-1] != '/')
		_fn_scratch += '/';
	fn_scratch = _fn_scratch + "relion_volatile/";

	if (do_reuse_scratch)
	{
		nr_parts_on_scratch.resize(numberOfOpticsGroups(), 0);
		for (int optics_group = 0; optics_group < numberOfOpticsGroups(); optics_group++)
		{
			if (is_3D)
			{
				FileName fn_tmp = fn_scratch + "opticsgroup" + integerToString(optics_group+1) + "_particle*.mrc";
				std::vector<FileName> fn_all;
				fn_tmp.globFiles(fn_all, true);
				nr_parts_on_scratch[optics_group] = fn_all.size();

			}
			else
			{
				FileName fn_tmp = fn_scratch + "opticsgroup" + integerToString(optics_group+1) + "_particles.mrcs";
				if (exists(fn_tmp))
				{
					Image<RFLOAT> Itmp;
					Itmp.read(fn_tmp, false);
					nr_parts_on_scratch[optics_group] = NSIZE(Itmp());
				}
#ifdef DEBUG_SCRATCH
				if (verb > 0)
					std::cerr << " optics_group= " << (optics_group + 1) << " nr_parts_on_scratch[optics_group]= " << nr_parts_on_scratch[optics_group] << std::endl;
#endif
			}
		}
	}
}

FileName Experiment::initialiseScratchLock(FileName _fn_scratch, FileName _fn_out)
{
	// Get a unique lockname for this run
	int uniqnr = rand() % 100000;
	FileName fn_uniq = _fn_out;
	fn_uniq.replaceAllSubstrings("/", "_");
	fn_uniq += "_lock" + integerToString(uniqnr);
	FileName fn_lock = fn_scratch + fn_uniq;

	if (exists(fn_lock))
		remove(fn_lock.c_str());

	return fn_lock;
}

bool Experiment::prepareScratchDirectory(FileName _fn_scratch, FileName fn_lock)
{
	if (fn_lock != "" && exists(fn_lock))
	{
		// Still measure how much free space there is
		struct statvfs vfs;
		statvfs(_fn_scratch.c_str(), &vfs);
		char nodename[64] = "undefined";
		gethostname(nodename,sizeof(nodename));
		std::string myhost(nodename);
		free_space_Gb = (RFLOAT)vfs.f_bsize * vfs.f_bfree / (1024 * 1024 * 1024);

		return false;
	}
	else
	{
		// Wipe the directory clean and make a new one
		std::string command;
		deleteDataOnScratch();

		// Make the scratch directory with write permissions
		command = "install -d -m 0777 " + fn_scratch;
		if (system(command.c_str()))
			REPORT_ERROR("ERROR: cannot execute: " + command);

		// Touch the lock file
		if(fn_lock != "")
		{
			touch(fn_lock);
			command = "chmod 0777 " + fn_lock;
			if (system(command.c_str()))
				REPORT_ERROR("ERROR: cannot execute: " + command);
		}

		// Measure how much free space there is
		struct statvfs vfs;
		statvfs(_fn_scratch.c_str(), &vfs);
		char nodename[64] = "undefined";
		gethostname(nodename,sizeof(nodename));
		std::string myhost(nodename);
		free_space_Gb = (RFLOAT)vfs.f_bsize * vfs.f_bfree / (1024 * 1024 * 1024);
		std::cout << " + On host " << myhost << ": free scratch space = " << free_space_Gb << " Gb." << std::endl;

		return true;
	}
}

void Experiment::deleteDataOnScratch()
{
	// Wipe the scratch directory
	if (fn_scratch != "" && exists(fn_scratch))
	{
		std::string command = " rm -rf " + fn_scratch;
		if (system(command.c_str()))
			REPORT_ERROR("ERROR: cannot execute: " + command);
	}
}

void Experiment::copyParticlesToScratch(int verb, bool do_copy, bool also_do_ctf_image, RFLOAT keep_free_scratch_Gb)
{

    // This function relies on prepareScratchDirectory() being called before!

	long int nr_part = particles.size();
	int barstep;
	if (verb > 0 && do_copy)
	{
		std::cout << " Copying particles to scratch directory: " << fn_scratch << std::endl;
		init_progress_bar(nr_part);
		barstep = XMIPP_MAX(1, nr_part / 60);
	}

	long int one_part_space, used_space = 0.;
	long int max_space = (free_space_Gb - keep_free_scratch_Gb) * 1024 * 1024 * 1024; // in bytes
#ifdef DEBUG_SCRATCH
	std::cerr << " free_space_Gb = " << free_space_Gb << " GB, keep_free_scratch_Gb = " << keep_free_scratch_Gb << " GB.\n";
	std::cerr << " Max space RELION can use = " << max_space << " bytes" << std::endl;
#endif
	// Loop over all particles and copy them one-by-one
	FileName fn_open_stack = "";
	fImageHandler hFile;
	long int total_nr_parts_on_scratch = 0;
	nr_parts_on_scratch.resize(numberOfOpticsGroups(), 0);

	const int check_abort_frequency=100;

	FileName prev_img_name = "/Unlikely$filename$?*!";
	int prev_optics_group = -999;
    for (long int part_id = 0; part_id < particles.size(); part_id++)
	{
		// TODO: think about MPI_Abort here....
		if (part_id % check_abort_frequency == 0 && pipeline_control_check_abort_job())
			exit(RELION_EXIT_ABORTED);

        long int imgno;
        FileName fn_ctf, fn_stack, fn_new;
        Image<RFLOAT> img;
        FileName fn_img = particles[part_id].name;
        int optics_group = particles[part_id].optics_group;

        // Get the size of the first particle
		if (nr_parts_on_scratch[optics_group] == 0)
		{
			Image<RFLOAT> tmp;
			tmp.read(fn_img, false); // false means: only read the header!
			one_part_space = NZYXSIZE(tmp())*sizeof(float); // MRC images are stored in floats!
			bool myis3D = (ZSIZE(tmp()) > 1);
			if (myis3D != is_3D)
				REPORT_ERROR("BUG: inconsistent is_3D values!");
			// add MRC header size for subtomograms (in 3D or as 2D stack), which are stored as 1 MRC file each
			if (is_3D || is_tomo) one_part_space += 1024;
            if (is_3D)
			{
				also_do_ctf_image = MDimg.containsLabel(EMDL_CTF_IMAGE);
				if (also_do_ctf_image)
					one_part_space *= 2;
			}
#ifdef DEBUG_SCRATCH
			std::cerr << "one_part_space[" << optics_group << "] = " << one_part_space << std::endl;
#endif
		}

        bool is_duplicate = (prev_img_name == fn_img && prev_optics_group == optics_group);
        // Read in the particle image, and write out on scratch
        if (do_copy && !is_duplicate)
        {
#ifdef DEBUG_SCRATCH
            std::cerr << "used_space = " << used_space << std::endl;
#endif
            // Now we have the particle in memory
            // See how much space it occupies
            used_space += one_part_space;
            // If there is no more space, exit the loop over all objects to stop copying files and change filenames in MDimg
            if (used_space > max_space)
            {
                char nodename[64] = "undefined";
                gethostname(nodename,sizeof(nodename));
                std::string myhost(nodename);
                std::cerr << " Warning: scratch space full on " << myhost << ". Remaining " << nr_part - total_nr_parts_on_scratch << " particles will be read from where they were."<< std::endl;
                break;
            }

            if (is_3D)
            {
                // For subtomograms, write individual .mrc files,possibly also CTF images
                img.read(fn_img);
                fn_new = fn_scratch + "opticsgroup" + integerToString(optics_group+1) + "_particle" + integerToString(nr_parts_on_scratch[optics_group]+1)+".mrc";
                img.write(fn_new);
                if (also_do_ctf_image)
                {
                    FileName fn_ctf;
                    MDimg.getValue(EMDL_CTF_IMAGE, fn_ctf, part_id);
                    img.read(fn_ctf);
                    fn_new = fn_scratch + "opticsgroup" + integerToString(optics_group+1) + "_particle_ctf" + integerToString(nr_parts_on_scratch[optics_group]+1)+".mrc";
                    img.write(fn_new);
                }
            }
            else if (is_tomo)
            {
                 // For subtomograms as 2D stacks, write individual .mrcs files
                img.read(fn_img);
                fn_new = fn_scratch + "opticsgroup" + integerToString(optics_group+1) + "_particle" + integerToString(nr_parts_on_scratch[optics_group]+1)+".mrcs";
                img.write(fn_new);
            }
            else
            {
                // Only open/close new stacks, so check if this is a new stack
                fn_img.decompose(imgno, fn_stack);
                if (fn_stack != fn_open_stack)
                {
                    // Manual closing isn't necessary: if still open, then openFile will first close the filehandler
                    // Also closing the last one isn't necessary, as destructor will do this.
                    //if (fn_open_stack != "")
                    //	hFile.closeFile();
                    hFile.openFile(fn_stack, WRITE_READONLY);
                    fn_open_stack = fn_stack;
                }
                img.readFromOpenFile(fn_img, hFile, -1, false);

                fn_new.compose(nr_parts_on_scratch[optics_group]+1, fn_scratch + "opticsgroup" + integerToString(optics_group+1) + "_particles.mrcs");
                if (nr_parts_on_scratch[optics_group] == 0)
                    img.write(fn_new, -1, false, WRITE_OVERWRITE);
                else
                    img.write(fn_new, -1, true, WRITE_APPEND);

#ifdef DEBUG_SCRATCH
                std::cerr << "Cached " << fn_img << " to " << fn_new << std::endl;
#endif
            }
        }

        // Update the counter and progress bar
        if (!is_duplicate)
            nr_parts_on_scratch[optics_group]++;
        total_nr_parts_on_scratch++;

        prev_img_name = fn_img;
        prev_optics_group = optics_group;

        if (verb > 0 && total_nr_parts_on_scratch % barstep == 0)
            progress_bar(total_nr_parts_on_scratch);

	} // end loop part_id

	if (verb)
	{
		progress_bar(nr_part);
		for (int i = 0; i < nr_parts_on_scratch.size(); i++)
		{
			std::cout << " For optics_group " << (i + 1) << ", there are " << nr_parts_on_scratch[i] << " particles on the scratch disk." << std::endl;
		}
	}

	if (do_copy && total_nr_parts_on_scratch>1)
	{
		std::string command = " chmod -R 777 " + fn_scratch + "/";
		if (system(command.c_str()))
			REPORT_ERROR("ERROR in executing: " + command);
	}
}


// Read from file
bool Experiment::read(FileName fn_exp, FileName fn_tomo, FileName fn_motion,
                      bool do_ignore_particle_name, bool do_ignore_group_name, bool do_preread_images,
                      bool need_tiltpsipriors_for_helical_refine, bool set_offset_priors_to_offsets, int verb)
{

//#define DEBUG_READ
#ifdef DEBUG_READ
	std::cerr << "Entering Experiment::read" << std::endl;
	Timer timer;
	int tall = timer.setNew("ALL");
	int tread = timer.setNew("read");
	int tsort = timer.setNew("sort");
	int tfill = timer.setNew("fill");
	int tgroup = timer.setNew("find group");
	int tdef = timer.setNew("set defaults");
	int tend = timer.setNew("ending");
	char c;
	timer.tic(tall);
	timer.tic(tread);
#endif

    bool remove_priors_again = false;

	is_tomo = false;

    // Only open stacks once and then read multiple images
	fImageHandler hFile;
	long int dump;
	FileName fn_stack, fn_open_stack="";

	// Initialize by emptying everything
	clear();
	long int group_id = 0, part_id = 0;

	if (!fn_exp.isStarFile())
	{
		REPORT_ERROR("ERROR: Input particles should be provides in a STAR file.");
	}
	else
	{

        // MDimg and MDopt have to be read at the same time, so that the optics groups can be
		// renamed in case they are non-contiguous or not sorted
		ObservationModel::loadSafely(fn_exp, obsModel, MDimg, "particles", verb);
        is_tomo = obsModel.isTomoStack2D;

		// The below is useful for filenames on scratch disk (related to avoiding copying duplicate particles when doing symmetry expansion)
        nr_particles_per_optics_group.resize(obsModel.numberOfOpticsGroups(), 0);

        std::vector<std::vector<ParticleIndex>> particles_idx;
        if (is_tomo)
        {

            if (fn_tomo == "") REPORT_ERROR("ERROR: you need to provide --tomograms when refining 2D stacks of tilt series images");

            // For now read in particle table twice: once into MDimg and once into particleSet...
            // TODO: work with pointer to avoid duplication?
            // TODO: ObservationModel::loadSafely(fn_exp, obsModel, MDimg, "particles", verb, true, true);
            tomogramSet.read(fn_tomo, false);
            particleSet.read(fn_exp, fn_motion, false);

            // For checking tomogram sanity below
            particles_idx = particleSet.splitByTomogram(tomogramSet, false);

        }

		// Set is_3D from MDopt
		int mydim=2;
		obsModel.opticsMdt.getValue(EMDL_IMAGE_DIMENSIONALITY, mydim, 0);
		is_3D = (mydim == 3);

#ifdef DEBUG_READ
		std::cerr << "Done reading MDimg" << std::endl;
		timer.toc(tread);
		timer.tic(tsort);
		//std::cerr << "Press any key to continue..." << std::endl;
		//std::cin >> c;
#endif

		// Sort input particles on micrograph (or tomogram) name
		EMDLabel my_sort_column = (is_3D) ? EMDL_TOMO_NAME : EMDL_MICROGRAPH_NAME;
		if (MDimg.containsLabel(my_sort_column))  MDimg.newSort(my_sort_column);

#ifdef DEBUG_READ
		std::cerr << "Done sorting MDimg" << std::endl;
		std::cerr << " MDimg.numberOfObjects()= " << MDimg.numberOfObjects() << std::endl;
		timer.toc(tsort);
		timer.tic(tfill);
		long nr_read = 0;
#endif
		// allocate 1 block of memory
		particles.reserve(MDimg.numberOfObjects());

        Tomogram tomogram;
        FileName prev_tomo_name = "";

		// Now Loop over all objects in the metadata file and fill the logical tree of the experiment
		for (long int part_id = 0; part_id < MDimg.numberOfObjects(); part_id++)
		{
			// Get the optics group of this particle
			int optics_group = obsModel.getOpticsGroup(MDimg, part_id);

#ifdef DEBUG_READ
			timer.tic(tgroup);
#endif

			// For example in particle_polishing the groups are not needed...
			if (!do_ignore_group_name)
			{
				FileName group_name="";
				// Check whether there is a group label, if not use a group for each micrograph
				if (MDimg.containsLabel(EMDL_MLMODEL_GROUP_NAME))
				{
					MDimg.getValue(EMDL_MLMODEL_GROUP_NAME, group_name, part_id);
				}
				else
				{
					FileName mic_name = getMicrographName(part_id);
                    FileName fn_pre, fn_jobnr;
					decomposePipelineFileName(mic_name, fn_pre, fn_jobnr, group_name);
				}

				// If this group did not exist yet, add it to the experiment
				group_id = -1;
				for (long int i = groups.size() - 1; i >= 0; i--) // search backwards to find match faster
				{
					if (groups[i].name == group_name)
					{
						group_id = groups[i].id;
						break;
					}
				}
				if (group_id < 0)
				{
					group_id = addGroup(group_name, optics_group);
				}
			}
			else
			{
				// All images belong to the same micrograph and group
				if (part_id == 0) group_id = addGroup("group", 0);
                else group_id = 0;
			}

            // The group number is only set upon reading: it is not read from the STAR file itself,
            // there the only thing that matters is the order of the micrograph_names
            // Write igroup+1, to start numbering at one instead of at zero
            MDimg.setValue(EMDL_MLMODEL_GROUP_NO, group_id + 1, part_id);

#ifdef DEBUG_READ
				timer.toc(tgroup);
#endif

			// If there is an EMDL_PARTICLE_RANDOM_SUBSET entry in the input STAR-file, then set the random_subset, otherwise use default (0)
			int my_random_subset;
			if (!MDimg.getValue(EMDL_PARTICLE_RANDOM_SUBSET, my_random_subset, part_id))
			{
				my_random_subset = 0;
			}


            FileName img_name;
            MDimg.getValue(EMDL_IMAGE_NAME, img_name, part_id);

			// Now get all the images for this particle (only one for SPA, multiple for STA)
            if (is_tomo)
            {

                // Find the tomogram this particle belongs to
                std::string tomo_name = MDimg.getString(EMDL_TOMO_NAME, part_id);
                int tomo_id = tomogramSet.getTomogramIndex(tomo_name);

                if (tomo_name != prev_tomo_name)
                {
                    tomogram = tomogramSet.loadTomogram(tomo_id, false);
                    prev_tomo_name = tomo_name;

                    // Tomogram sanity checks
                    tomogram.validateParticleOptics(particles_idx[tomo_id], particleSet);
                    particleSet.checkTrajectoryLengths(particles_idx[tomo_id], tomogram.frameCount, "exp_model.read");
                }

                // Add this particle to the Experiment, with its tomogram
                addParticle(img_name, optics_group, group_id, my_random_subset, tomo_id);

                // Pre-orientation of this particle in the tomogram
                ParticleIndex id(part_id);
                d3Matrix A = particleSet.getSubtomogramMatrix(id);
                const d3Vector pos = particleSet.getPosition(id, tomogram.centre, true);


                // Add only the visible images for this particle
                std::vector<int> isVisible = particleSet.getVisibleFrames(id);
                for (int f = 0; f < tomogram.frameCount; f++)
                {

                    if (isVisible[f] > 0)
                    {
                        float dose = tomogram.getCumulativeDose(f);

                        d4Matrix P = tomogram.projectionMatrices[f] * d4Matrix(A);

                        CTF ctf = tomogram.getCtf(f, pos);

                        addImageToParticle(part_id, &P, &ctf, dose, tomogram.BfactorPerElectronDose);
                    }
                }
            }
            else
            {
                // Add this particle to the Experiment
                addParticle(img_name, optics_group, group_id, my_random_subset);

                // Create a new image in this particle
                addImageToParticle(part_id);
            }
#ifdef DEBUG_READ
                timer.tic(tori);
#endif
            if (do_preread_images)
            {
                Image<float> img;
                if (is_tomo || is_3D)
                {
                    img.read(img_name);
                    particles[part_id].img = img();
                }
                else
                {
                    img_name.decompose(dump, fn_stack);
                    if (fn_stack != fn_open_stack)
                    {
                        hFile.openFile(fn_stack, WRITE_READONLY);
                        fn_open_stack = fn_stack;
                    }
                    img.readFromOpenFile(img_name, hFile, -1, false);
                    img().setXmippOrigin();
                    particles[part_id].img = img();
    			}
            }

#ifdef DEBUG_READ
			timer.toc(tori);
#endif

#ifdef DEBUG_READ
			nr_read++;
#endif
		} // end loop over all objects in MDimg (part_id)

#ifdef DEBUG_READ
		timer.toc(tfill);
		timer.tic(tdef);
		std::cerr << " nr_read= " << nr_read << " particles.size()= " << particles.size() << " groups.size()= " << groups.size() << std::endl;
#endif

		// Check for the presence of multiple bodies (for multi-body refinement)
		bool is_done = false;
		nr_bodies = 0;
		while (!is_done)
		{
			std::string tablename = "images_body_" + integerToString(nr_bodies+1);
                        MetaDataTable MDimgin;
			if (MDimgin.read(fn_exp, tablename) > 0)
			{
				nr_bodies++;
				MDbodies.push_back(MDimgin);
			}
			else
			{
				is_done = true;
			}
		}
		// Even if we don't do multi-body refinement, then nr_bodies is still 1
		nr_bodies = XMIPP_MAX(nr_bodies, 1);
	}

#ifdef DEBUG_READ
	std::cerr << "Done filling MDimg" << std::endl;
	//std::cerr << "Press any key to continue..." << std::endl;
	//std::cin >> c;
#endif

	// Make sure some things are always set in the MDimg
	bool have_rot  = MDimg.containsLabel(EMDL_ORIENT_ROT);
	bool have_tilt = MDimg.containsLabel(EMDL_ORIENT_TILT);
	bool have_psi  = MDimg.containsLabel(EMDL_ORIENT_PSI);
	bool have_xoff = MDimg.containsLabel(EMDL_ORIENT_ORIGIN_X_ANGSTROM);
	bool have_yoff = MDimg.containsLabel(EMDL_ORIENT_ORIGIN_Y_ANGSTROM);
	bool have_zoff = MDimg.containsLabel(EMDL_ORIENT_ORIGIN_Z_ANGSTROM);
	bool have_zcoord = MDimg.containsLabel(EMDL_IMAGE_COORD_Z);
	bool have_clas = MDimg.containsLabel(EMDL_PARTICLE_CLASS);
	bool have_norm = MDimg.containsLabel(EMDL_IMAGE_NORM_CORRECTION);

	// Jan20,2016 - Helical reconstruction
	bool have_tilt_prior = MDimg.containsLabel(EMDL_ORIENT_TILT_PRIOR);
	bool have_psi_prior = MDimg.containsLabel(EMDL_ORIENT_PSI_PRIOR);
	bool have_tiltpsi = (have_tilt) && (have_psi);
	bool have_tiltpsi_prior = (have_tilt_prior) && (have_psi_prior);
	if (need_tiltpsipriors_for_helical_refine)
	{
		if (!have_tiltpsi_prior)
		{
			if (!have_tiltpsi)
				REPORT_ERROR("exp_model.cpp: Experiment::read(): Tilt and psi priors of helical segments are missing!");
		}
	}
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDimg)
	{
		RFLOAT dzero=0., done=1.;
		int izero = 0;
		if (!have_rot)
			MDimg.setValue(EMDL_ORIENT_ROT, dzero);
		if (!have_tilt)
			MDimg.setValue(EMDL_ORIENT_TILT, dzero);
		if (!have_psi)
			MDimg.setValue(EMDL_ORIENT_PSI, dzero);
		if (!have_xoff)
			MDimg.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, dzero);
		if (!have_yoff)
			MDimg.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, dzero);
		if ( (!have_zoff) && (have_zcoord) )
			MDimg.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, dzero);
		if (!have_clas)
			MDimg.setValue(EMDL_PARTICLE_CLASS, izero);
		if (!have_norm)
			MDimg.setValue(EMDL_IMAGE_NORM_CORRECTION, done);
		if (need_tiltpsipriors_for_helical_refine && have_tiltpsi_prior) // If doing 3D helical reconstruction and PRIORs exist
		{
			RFLOAT tilt = 0., psi = 0.;
			if (have_tiltpsi)
				MDimg.getValue(EMDL_ORIENT_TILT, tilt);
			// If ANGLEs do not exist or they are all set to 0 (from a Class2D job), copy values of PRIORs to ANGLEs
			if ( (!have_tiltpsi) || ((have_tiltpsi) && (ABS(tilt) < 0.001)) )
			{
				MDimg.getValue(EMDL_ORIENT_TILT_PRIOR, tilt);
				MDimg.getValue(EMDL_ORIENT_PSI_PRIOR, psi);
				MDimg.setValue(EMDL_ORIENT_TILT, tilt);
				MDimg.setValue(EMDL_ORIENT_PSI, psi);
			}
		}

        if (set_offset_priors_to_offsets)
        {
            RFLOAT off;
            if (!MDimg.containsLabel(EMDL_ORIENT_ORIGIN_X_PRIOR_ANGSTROM))
            {
                MDimg.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, off);
                MDimg.setValue(EMDL_ORIENT_ORIGIN_X_PRIOR_ANGSTROM, off);
                remove_priors_again = true;
            }
            if (!MDimg.containsLabel(EMDL_ORIENT_ORIGIN_Y_PRIOR_ANGSTROM))
            {
                MDimg.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, off);
                MDimg.setValue(EMDL_ORIENT_ORIGIN_Y_PRIOR_ANGSTROM, off);
            }
            if (have_zcoord && !MDimg.containsLabel(EMDL_ORIENT_ORIGIN_Z_PRIOR_ANGSTROM))
            {
                MDimg.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, off);
                MDimg.setValue(EMDL_ORIENT_ORIGIN_Z_PRIOR_ANGSTROM, off);
            }
        }
	}

    // Make sure the partTable in particleSet is kept the same as MDimg in Experiment, for reading from it in ml_optimiser.cpp
    // TODO! Make use of pointers to avoid duplication of entire MDimg here...
    if (is_tomo) particleSet.partTable = MDimg;

    // Keep track whether priors that were added above should be removed again later...
    return remove_priors_again;

#ifdef DEBUG_READ
	timer.toc(tdef);
	std::cerr << "Done setting defaults MDimg" << std::endl;
	timer.toc(tall);
	timer.printTimes(false);
	//std::cerr << "Writing out debug_data.star" << std::endl;
	//write("debug");
	//exit(0);
#endif
}

// Write to file
void Experiment::write(FileName fn_out, bool remove_offset_priors)
{
    std::ofstream ofs(fn_out);

    if (is_tomo) obsModel.generalMdt.write(ofs);

    obsModel.opticsMdt.write(ofs);
    if (remove_offset_priors)
    {
        if (MDimg.containsLabel(EMDL_ORIENT_ORIGIN_X_PRIOR_ANGSTROM))
            MDimg.deactivateLabel(EMDL_ORIENT_ORIGIN_X_PRIOR_ANGSTROM);
        if (MDimg.containsLabel(EMDL_ORIENT_ORIGIN_Y_PRIOR_ANGSTROM))
            MDimg.deactivateLabel(EMDL_ORIENT_ORIGIN_Y_PRIOR_ANGSTROM);
        if (MDimg.containsLabel(EMDL_ORIENT_ORIGIN_Z_PRIOR_ANGSTROM))
            MDimg.deactivateLabel(EMDL_ORIENT_ORIGIN_Z_PRIOR_ANGSTROM);
    }
    MDimg.write(ofs);

    if (nr_bodies > 1)
    {
        for (int ibody = 0; ibody < nr_bodies; ibody++)
        {
            MDbodies[ibody].write(ofs);
        }
    }

}
