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

void ExpOriginalParticle::addParticle(long int _particle_id, int _random_subset, int _order)
{
	// Keep random_subsets equal in each original particle
	if (random_subset != _random_subset)
	{
		std::cerr << " random_subset= " << random_subset << " _random_subset= " << _random_subset << std::endl;
		std::cerr << " particles_id.size()= " << particles_id.size() << " name= "<<name << std::endl;
		REPORT_ERROR("ExpOriginalParticle:addParticle: incompatible random subsets between particle and its original particle.");
	}
	particles_id.push_back(_particle_id);
	particles_order.push_back(_order);
}

long int Experiment::numberOfParticles(int random_subset)
{
	if (random_subset == 0)
		return particles.size();
	else
	{
		long int result = 0;
		for (long int i = 0; i < ori_particles.size(); i++)
		{
			if (ori_particles[i].random_subset == random_subset)
			{
				result += ori_particles[i].particles_id.size();
			}
		}
		return result;
	}
}

long int Experiment::numberOfOriginalParticles(int random_subset)
{
	if (random_subset == 0)
		return ori_particles.size();
	else if (random_subset == 1)
		return nr_ori_particles_subset1;
	else if (random_subset == 2)
		return nr_ori_particles_subset2;
	else
		REPORT_ERROR("ERROR: Experiment::numberOfOriginalParticles invalid random_subset: " + integerToString(random_subset));
}


long int Experiment::numberOfMicrographs()
{
	return micrographs.size();
}

long int Experiment::numberOfGroups()
{
	return groups.size();
}

long int Experiment::getMicrographId(long int part_id)
{
#ifdef DEBUG_CHECKSIZES
	if (part_id >= particles.size())
	{
		std::cerr<< "part_id= "<<part_id<<" particles.size()= "<< particles.size() <<std::endl;
		REPORT_ERROR("part_id >= particles.size()");
	}
#endif
	return (particles[part_id]).micrograph_id;
}

long int Experiment::getGroupId(long int part_id)
{
#ifdef DEBUG_CHECKSIZES
	if (part_id >= particles.size())
	{
		std::cerr<< "part_id= "<<part_id<<" particles.size()= "<< particles.size() <<std::endl;
		REPORT_ERROR("part_id >= particles.size()");
	}
#endif
	return (particles[part_id]).group_id;
}

int Experiment::getRandomSubset(long int part_id)
{
	return particles[part_id].random_subset;
}

MetaDataTable Experiment::getMetaDataImage(long int part_id)
{
	MetaDataTable result;
	result.addObject(MDimg.getObject(part_id));
	return result;
}

long int Experiment::addParticle(long int group_id,
		long int micrograph_id, int random_subset)
{

	if (group_id >= groups.size())
		REPORT_ERROR("Experiment::addImage: group_id out of range");

	if (micrograph_id >= micrographs.size())
		REPORT_ERROR("Experiment::addImage: micrograph_id out of range");

	ExpParticle particle;
	particle.id = particles.size();
	particle.group_id = group_id;
	particle.micrograph_id = micrograph_id;
	particle.random_subset = random_subset;
	// Push back this particle in the particles vector
	particles.push_back(particle);
	(micrographs[micrograph_id].particle_ids).push_back(particle.id);

	// Return the id in the particles vector
	return particle.id;

}

long int Experiment::addOriginalParticle(std::string part_name, int _random_subset)
{

	ExpOriginalParticle ori_particle;
	ori_particle.random_subset = _random_subset;
	ori_particle.name = part_name;
	long int id = ori_particles.size();
	ori_particles.push_back(ori_particle);

	// Return the id in the ori_particles vector
	return id;

}

long int Experiment::addGroup(std::string group_name)
{
	// Add new group to this Experiment
	ExpGroup group;
	group.id = groups.size(); // start counting groups at 0!
	group.name = group_name;

	// Push back this micrograph
	groups.push_back(group);

	// Return the id in the micrographs vector
	return group.id;

}

long int Experiment::addMicrograph(std::string mic_name)
{
	// Add new micrograph to this Experiment
	ExpMicrograph micrograph;
	micrograph.id = micrographs.size();
	micrograph.name = mic_name;

	// Push back this micrograph
	micrographs.push_back(micrograph);

	// Return the id in the micrographs vector
	return micrograph.id;

}

void Experiment::divideOriginalParticlesInRandomHalves(int seed, bool do_helical_refine)
{

	// Only do this if the random_subset of all original_particles is zero
	bool all_are_zero = true;
	bool some_are_zero = false;
	nr_ori_particles_subset1 = 0;
	nr_ori_particles_subset2 = 0;
	for (long int i = 0; i < ori_particles.size(); i++)
	{
		int random_subset = ori_particles[i].random_subset;
		if (random_subset != 0)
		{
			all_are_zero = false;
			// Keep track of how many particles there are in each subset
			if (random_subset == 1)
				nr_ori_particles_subset1++;
			else if (random_subset == 2)
				nr_ori_particles_subset2++;
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
			bool divide_according_to_helical_tube_id;
			std::map<std::string, int> map_mics;
			std::map<std::string, int>::const_iterator ii_map;
			std::vector<std::pair<std::string, int> > vec_mics;

			for (long int i = 0; i < ori_particles.size(); i++)
			{
				if ( ((ori_particles[i]).particles_id.size() != 1) || ((ori_particles[i]).particles_order.size() != 1) )
					REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: cannot divide helical segments into random halves with tilt series or movie frames!");
			}

			if ( (!MDimg.containsLabel(EMDL_IMAGE_NAME)) || (!MDimg.containsLabel(EMDL_MICROGRAPH_NAME)) )
				REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: Input MetadataTable should contain rlnImageName and rlnMicrographName!");

			divide_according_to_helical_tube_id = false;
			if (MDimg.containsLabel(EMDL_PARTICLE_HELICAL_TUBE_ID))
				divide_according_to_helical_tube_id = true;

			// Count micrograph names
			map_mics.clear();
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDimg)
			{
				// == OLD ==
				//MDimg.getValue(EMDL_IMAGE_NAME, img_name);
				//mic_name = img_name.substr(img_name.find("@") + 1);
				// == NEW == compatible with 3D subtomograms
				MDimg.getValue(EMDL_MICROGRAPH_NAME, img_name);
				mic_name = img_name;

				if (divide_according_to_helical_tube_id)
				{
					MDimg.getValue(EMDL_PARTICLE_HELICAL_TUBE_ID, helical_tube_id);
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

//#define OLD_RANDOMISATION
#ifdef OLD_RANDOMISATION
			// Perform swaps (OLD RANDOMISATION - not perfect)
			nr_swaps = ROUND(rnd_unif(vec_mics.size(), 2. * vec_mics.size()));
			for (int ii = 0; ii < nr_swaps; ii++)
			{
				int ptr_a, ptr_b;
				std::pair<std::string, int> tmp;
				ptr_a = ROUND(rnd_unif(0, vec_mics.size()));
				ptr_b = ROUND(rnd_unif(0, vec_mics.size()));
				if ( (ptr_a == ptr_b) || (ptr_a < 0) || (ptr_b < 0) || (ptr_a >= vec_mics.size()) || (ptr_b >= vec_mics.size()) )
					continue;
				tmp = vec_mics[ptr_a];
				vec_mics[ptr_a] = vec_mics[ptr_b];
				vec_mics[ptr_b] = tmp;

				// DEBUG
				//std::cout << " Swap mic_id= " << ptr_a << " with mic_id= " << ptr_b << "." << std::endl;
			}
#else
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

				// DEBUG
				//std::cout << " Swap mic_id= " << ptr_a << " with mic_id= " << ptr_b << "." << std::endl;
			}
#endif
			// DEBUG
			//if (divide_according_to_helical_tube_id)
			//	std::cout << " Helical tubes= " << vec_mics.size() << ", nr_swaps= " << nr_swaps << std::endl;
			//else
			//	std::cout << " Micrographs= " << vec_mics.size() << ", nr_swaps= " << nr_swaps << std::endl;

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

			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDimg)
			{
				MDimg.getValue(EMDL_MICROGRAPH_NAME, img_name);
				mic_name = img_name;

				if (divide_according_to_helical_tube_id)
				{
					MDimg.getValue(EMDL_PARTICLE_HELICAL_TUBE_ID, helical_tube_id);
					if (helical_tube_id < 1)
						REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: Helical tube ID should be positive integer!");
					mic_name += std::string("_TUBEID_");
					mic_name += std::string(integerToString(helical_tube_id));
				}
				MDimg.setValue(EMDL_PARTICLE_RANDOM_SUBSET, map_mics[mic_name]);
			}

			nr_ori_particles_subset1 = nr_segments_subset1;
			nr_ori_particles_subset2 = nr_segments_subset2;

			// DEBUG
			//std::cout << " Helical segments in two half sets = " << nr_ori_particles_subset1 << ", " << nr_ori_particles_subset2 << std::endl;

			// Loop over all particles in each ori_particle and set their random_subset
			for (long int i = 0; i < ori_particles.size(); i++)
			{
				int random_subset;
				long int part_id = (ori_particles[i]).particles_id[0];
				MDimg.getValue(EMDL_PARTICLE_RANDOM_SUBSET, random_subset, part_id);
				ori_particles[i].random_subset = random_subset;
				particles[part_id].random_subset = random_subset;
			}
		}
		else
		{
			for (long int i = 0; i < ori_particles.size(); i++)
			{
				int random_subset = rand() % 2 + 1;
				ori_particles[i].random_subset = random_subset; // randomly 1 or 2
				if (random_subset == 1)
					nr_ori_particles_subset1++;
				else if (random_subset == 2)
					nr_ori_particles_subset2++;
				else
					REPORT_ERROR("ERROR Experiment::divideParticlesInRandomHalves: invalid number for random subset (i.e. not 1 or 2): " + integerToString(random_subset));

				// Loop over all particles in each ori_particle and set their random_subset
				for (long int j = 0; j < ori_particles[i].particles_id.size(); j++)
				{
					long int part_id = (ori_particles[i]).particles_id[j];
					{
						particles[part_id].random_subset = random_subset;
						MDimg.setValue(EMDL_PARTICLE_RANDOM_SUBSET, random_subset, part_id);
					}
				}
			}
		}
	}


	// Now re-order such that half1 is in first half, and half2 is in second half of the particle list (for MPI_parallelisattion)

	std::vector<long int> ori_particle_list1, ori_particle_list2;
	// Fill the two particle lists
	for (long int i = 0; i < ori_particles.size(); i++)
	{
		int random_subset = ori_particles[i].random_subset;
		if (random_subset == 1)
			ori_particle_list1.push_back(i);
		else if (random_subset == 2)
			ori_particle_list2.push_back(i);
		else
			REPORT_ERROR("ERROR: invalid number for random subset (i.e. not 1 or 2): " + integerToString(random_subset));
	}

	// Just a silly check for the sizes of the ori_particle_lists (to be sure)
	if (ori_particle_list1.size() != nr_ori_particles_subset1)
		REPORT_ERROR("ERROR: invalid ori_particle_list1 size:" + integerToString(ori_particle_list1.size()) + " != " + integerToString(nr_ori_particles_subset1));
	if (ori_particle_list2.size() != nr_ori_particles_subset2)
		REPORT_ERROR("ERROR: invalid ori_particle_list2 size:" + integerToString(ori_particle_list2.size()) + " != " + integerToString(nr_ori_particles_subset2));

	// First fill new_ori_particles with the first subset, then with the second
	std::vector<ExpOriginalParticle> new_ori_particles;

	for (long int i = 0; i < ori_particle_list1.size(); i++)
		new_ori_particles.push_back(ori_particles[ori_particle_list1[i]]);
	for (long int i = 0; i < ori_particle_list2.size(); i++)
		new_ori_particles.push_back(ori_particles[ori_particle_list2[i]]);

	ori_particles=new_ori_particles;

	if (nr_ori_particles_subset2 == 0 || nr_ori_particles_subset1 == 0)
		REPORT_ERROR("ERROR: one of your half sets has no segments. Helical half-sets are always per-filament. Provide at least 2 filaments.");


}

void Experiment::randomiseOriginalParticlesOrder(int seed, bool do_split_random_halves)
{
	//This static flag is for only randomize once
	static bool randomised = false;
	if (!randomised)
	{

		srand(seed);
		std::vector<ExpOriginalParticle> new_ori_particles;

		if (do_split_random_halves)
		{
			std::vector<long int> ori_particle_list1, ori_particle_list2;
			ori_particle_list1.clear();
			ori_particle_list2.clear();
			// Fill the two particle lists
			for (long int i = 0; i < ori_particles.size(); i++)
			{
				int random_subset = ori_particles[i].random_subset;
				if (random_subset == 1)
					ori_particle_list1.push_back(i);
				else if (random_subset == 2)
					ori_particle_list2.push_back(i);
				else
					REPORT_ERROR("ERROR Experiment::randomiseParticlesOrder: invalid number for random subset (i.e. not 1 or 2): " + integerToString(random_subset));
			}

			// Just a silly check for the sizes of the ori_particle_lists (to be sure)
			if (ori_particle_list1.size() != nr_ori_particles_subset1)
				REPORT_ERROR("ERROR Experiment::randomiseParticlesOrder: invalid ori_particle_list1 size:" + integerToString(ori_particle_list1.size()) + " != " + integerToString(nr_ori_particles_subset1));
			if (ori_particle_list2.size() != nr_ori_particles_subset2)
				REPORT_ERROR("ERROR Experiment::randomiseParticlesOrder: invalid ori_particle_list2 size:" + integerToString(ori_particle_list2.size()) + " != " + integerToString(nr_ori_particles_subset2));

			// Randomise the two particle lists
			std::random_shuffle(ori_particle_list1.begin(), ori_particle_list1.end());
			std::random_shuffle(ori_particle_list2.begin(), ori_particle_list2.end());

			// First fill new_ori_particles with the first subset, then with the second
			for (long int i = 0; i < ori_particle_list1.size(); i++)
				new_ori_particles.push_back(ori_particles[ori_particle_list1[i]]);
			for (long int i = 0; i < ori_particle_list2.size(); i++)
				new_ori_particles.push_back(ori_particles[ori_particle_list2[i]]);

		}
		else
		{

			// First fill in order
			std::vector<long int> ori_particle_list;
			ori_particle_list.resize(ori_particles.size());
			for (long int i = 0; i < ori_particle_list.size(); i++)
				ori_particle_list[i] = i;

			// Randomise
			std::random_shuffle(ori_particle_list.begin(), ori_particle_list.end());

			// Refill new_ori_particles
			for (long int i = 0; i < ori_particle_list.size(); i++)
				new_ori_particles.push_back(ori_particles[ori_particle_list[i]]);
		}

		ori_particles=new_ori_particles;
		randomised = true;

	}
}


void Experiment::orderParticlesInOriginalParticles()
{
	// If the orders are negative (-1) then dont sort anything
	if (ori_particles[0].particles_order[0] < 0)
		return;

	for (long int i = 0; i < ori_particles.size(); i++)
	{
		int nframe = ori_particles[i].particles_order.size();

		std::vector<std::pair<long int, long int> > vp;
        vp.reserve(nframe);
        for (long int j = 0; j < nframe; j++)
        	vp.push_back(std::make_pair(ori_particles[i].particles_order[j], j));
        // Sort on the first elements of the pairs
        std::sort(vp.begin(), vp.end());

        // tmp copy of particles_id
        std::vector<long int> _particles_id = ori_particles[i].particles_id;
        for (int j = 0; j < nframe; j++)
			ori_particles[i].particles_id[j] = _particles_id[vp[j].second];

		// We now no longer need the particles_order vector, clear it to save memory
		ori_particles[i].particles_order.clear();
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
			RFLOAT norm, zero=0.;
			MDimg.getValue(EMDL_IMAGE_NORM_CORRECTION, norm);
			MDbody.setValue(EMDL_ORIENT_ORIGIN_X, zero);
			MDbody.setValue(EMDL_ORIENT_ORIGIN_Y, zero);
			MDbody.setValue(EMDL_ORIENT_ROT, zero);
			MDbody.setValue(EMDL_ORIENT_TILT, zero);
			MDbody.setValue(EMDL_ORIENT_PSI, zero);
			MDbody.setValue(EMDL_IMAGE_NORM_CORRECTION, norm);
			if (is_3d)
			{
				MDbody.setValue(EMDL_ORIENT_ORIGIN_Z, zero);
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
	if (fn_scratch != "" && part_id < nr_parts_on_scratch)
	{
		if (is_3D)
		{
			if (is_ctf_image)
				fn_img = fn_scratch + "particle_ctf" + integerToString(part_id+1, 5)+".mrc";
			else
				fn_img = fn_scratch + "particle" + integerToString(part_id+1, 5)+".mrc";
		}
		else
		{

			if (ori_particles[part_id].particles_id.size() > 1)
			{
				std::cerr << " part_id= " << part_id << " ori_particles[part_id].particles_id.size()= " << ori_particles[part_id].particles_id.size() << " ori_particles[part_id].name= " << ori_particles[part_id].name << std::endl;
				REPORT_ERROR("BUG: getImageNameOnScratch cannot work with movies!");
			}
			fn_img.compose(part_id+1, fn_scratch + "particles.mrcs");
		}
		return true;
	}
	else
	{
		return false;
	}

}

FileName Experiment::initialiseScratchLock(FileName _fn_scratch, FileName _fn_out)
{
    // Make sure fn_scratch ends with a slash
	if (_fn_scratch[_fn_scratch.length()-1] != '/')
		_fn_scratch += '/';
	fn_scratch = _fn_scratch + "relion_volatile/";

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

    // Make sure fn_scratch ends with a slash
	if (_fn_scratch[_fn_scratch.length()-1] != '/')
		_fn_scratch += '/';
	fn_scratch = _fn_scratch + "relion_volatile/";

	if (fn_lock != "" && exists(fn_lock))
	{
		// Still measure how much free space there is
		struct statvfs vfs;
		statvfs(_fn_scratch.c_str(), &vfs);
		long int free_Gb = vfs.f_bsize*vfs.f_bfree/(1024*1024*1024);
	    char nodename[64] = "undefined";
	    gethostname(nodename,sizeof(nodename));
	    std::string myhost(nodename);
	    free_space_Gb = vfs.f_bsize*vfs.f_bfree/(1024*1024*1024);

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
		long int free_Gb = vfs.f_bsize*vfs.f_bfree/(1024*1024*1024);
	    char nodename[64] = "undefined";
	    gethostname(nodename,sizeof(nodename));
	    std::string myhost(nodename);
	    free_space_Gb = vfs.f_bsize*vfs.f_bfree/(1024*1024*1024);
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

void Experiment::copyParticlesToScratch(int verb, bool do_copy, bool also_do_ctf_image, long int keep_free_scratch_Gb)
{

	// This function relies on prepareScratchDirectory() being called before!

	long int nr_part = MDimg.numberOfObjects();
	int barstep;
	if (verb > 0)
	{
		std::cout << " Copying particles to scratch directory: " << fn_scratch << std::endl;
		init_progress_bar(nr_part);
		barstep = XMIPP_MAX(1, nr_part / 60);
	}

	long int one_part_space, used_space = 0.;
	long int max_space = (free_space_Gb - keep_free_scratch_Gb)*1024*1024*1024; // in bytes

	// Loop over all particles and copy them one-by-one
	FileName fn_open_stack = "";
	fImageHandler hFile;
	nr_parts_on_scratch = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDimg)
	{
		long int imgno;
		FileName fn_img, fn_ctf, fn_stack, fn_new;
		Image<RFLOAT> img;
		MDimg.getValue(EMDL_IMAGE_NAME, fn_img);

		// Get the size of the first particle
		if (nr_parts_on_scratch == 0)
		{
			Image<RFLOAT> tmp;
			tmp.read(fn_img, false); // false means: only read the header!
			one_part_space = ZYXSIZE(tmp())*sizeof(float); // MRC images are stored in floats!
			bool myis3D = (ZSIZE(tmp()) > 1);
			if (myis3D != is_3D)
				REPORT_ERROR("BUG: inconsistent is_3D values!");
			// add MRC header size for subtomograms, which are stored as 1 MRC file each
			if (is_3D)
			{
				one_part_space += 1024;
				also_do_ctf_image = MDimg.containsLabel(EMDL_CTF_IMAGE);
				if (also_do_ctf_image)
					one_part_space *= 2;
			}
		}

		// Now we have the particle in memory
		// See how much space it occupies
		used_space += one_part_space;
		// If there is no more space, exit the loop over all objects to stop copying files and change filenames in MDimg
		if (used_space > max_space)
		{
			char nodename[64] = "undefined";
			gethostname(nodename,sizeof(nodename));
			std::string myhost(nodename);
			std::cerr << " Warning: scratch space full on " << myhost << ". Remaining " << nr_part - nr_parts_on_scratch << " particles will be read from where they were."<< std::endl;
			break;
		}

		// Read in the particle image, and write out on scratch
		if (do_copy)
		{
			if (is_3D)
			{
				// For subtomograms, write individual .mrc files,possibly also CTF images
				img.read(fn_img);
				fn_new = fn_scratch + "particle" + integerToString(nr_parts_on_scratch+1, 5)+".mrc";
				img.write(fn_new);
				if (also_do_ctf_image)
				{
					FileName fn_ctf;
					MDimg.getValue(EMDL_CTF_IMAGE, fn_ctf);
					img.read(fn_ctf);
					fn_new = fn_scratch + "particle_ctf" + integerToString(nr_parts_on_scratch+1, 5)+".mrc";
					img.write(fn_new);
				}
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

				fn_new.compose(nr_parts_on_scratch+1, fn_scratch + "particles.mrcs");
				if (nr_parts_on_scratch == 0)
					img.write(fn_new, -1, false, WRITE_OVERWRITE);
				else
					img.write(fn_new, -1, true, WRITE_APPEND);
			}
		}

		// Update the counter and progress bar
		nr_parts_on_scratch++;

		if (verb > 0 && nr_parts_on_scratch % barstep == 0)
			progress_bar(nr_parts_on_scratch);

	}

	if (verb > 0)
		progress_bar(nr_part);

	if (do_copy && nr_parts_on_scratch>1)
	{
		std::string command = " chmod 777 " + fn_scratch + "particle*";
		if (system(command.c_str()))
			REPORT_ERROR("ERROR in executing: " + command);
	}

}

void Experiment::usage()
{
	std::cout
	<< "  -i                     : Starfile with input images\n"
	;
}

// Read from file
void Experiment::read(FileName fn_exp, bool do_ignore_original_particle_name,
		bool do_ignore_group_name, bool do_preread_images,
		bool need_tiltpsipriors_for_helical_refine)
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
	int tori = timer.setNew("find ori_particle");
	int tdef = timer.setNew("set defaults");
	int tend = timer.setNew("ending");
	char c;
	timer.tic(tall);
	timer.tic(tread);
#endif

	// Only open stacks once and then read multiple images
	fImageHandler hFile;
	long int dump;
	FileName fn_stack, fn_open_stack="";

	// Initialize by emptying everything
	clear();
	long int group_id, mic_id, part_id;

	if (!fn_exp.isStarFile())
	{
		// Read images from stack. Ignore all metadata, just use filenames
		// Add a single Micrograph
		group_id = addGroup("group");
		mic_id = addMicrograph("micrograph");

		// Check that a MRC stack ends in .mrcs, not .mrc (which will be read as a MRC 3D map!)
		if (fn_exp.contains(".mrc") && !fn_exp.contains(".mrcs"))
			REPORT_ERROR("Experiment::read: ERROR: MRC stacks of 2D images should be have extension .mrcs, not .mrc!");

		// Read in header-only information to get the NSIZE of the stack
		Image<RFLOAT> img;
		img.read(fn_exp, false); // false means skip data, only read header

		// allocate 1 block of memory
		particles.reserve(NSIZE(img()));
		ori_particles.reserve(NSIZE(img()));
		for (long int n = 0; n <  NSIZE(img()); n++)
		{
			FileName fn_img;
			fn_img.compose(n+1, fn_exp); // fn_img = integerToString(n) + "@" + fn_exp;
			// Add the particle to my_area = 0
			part_id = addParticle(group_id, mic_id);
			MDimg.addObject();
			if (do_preread_images)
			{
				Image<float> img;
				fn_img.decompose(dump, fn_stack);
				if (fn_stack != fn_open_stack)
				{
					hFile.openFile(fn_stack, WRITE_READONLY);
					fn_open_stack = fn_stack;
				}
				img.readFromOpenFile(fn_img, hFile, -1, false);
				img().setXmippOrigin();
				particles[part_id].img = img();
			}
			// Also add OriginalParticle
			(ori_particles[addOriginalParticle("particle")]).addParticle(part_id, 0, -1);
			// Set the filename and other metadata parameters
			MDimg.setValue(EMDL_IMAGE_NAME, fn_img, part_id);
		}

	}
	else
	{
		// Just read first data block
		MDimg.read(fn_exp);

		// If this for movie-processing, then the STAR file might be a list of STAR files with original movie particles.
		// If that is the case: then just append all into a new MDimg.
		if (MDimg.containsLabel(EMDL_STARFILE_MOVIE_PARTICLES))
		{
			MetaDataTable MDjoined, MDone;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDimg)
			{
				FileName fn_star;
				MDimg.getValue(EMDL_STARFILE_MOVIE_PARTICLES, fn_star);
				MDone.read(fn_star);
				MDjoined.append(MDone);
			}
			MDimg = MDjoined;
		}

#ifdef DEBUG_READ
		std::cerr << "Done reading MDimg" << std::endl;
		timer.toc(tread);
		timer.tic(tsort);
		//std::cerr << "Press any key to continue..." << std::endl;
		//std::cin >> c;
#endif

		// Sort input particles on micrographname
		bool is_mic_a_movie=false, star_contains_micname;
		star_contains_micname = MDimg.containsLabel(EMDL_MICROGRAPH_NAME);
		if (star_contains_micname)
		{
			// See if the micrograph names contain an "@", i.e. whether they are movies and we are inside polishing or so.
			FileName fn_mic;
			MDimg.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
			if (fn_mic.contains("@"))
			{
				is_mic_a_movie = true;
				MDimg.newSort(EMDL_MICROGRAPH_NAME, false, true); // sort on part AFTER "@"
			}
			else
			{
				is_mic_a_movie = false;
				MDimg.newSort(EMDL_MICROGRAPH_NAME); // just sort on fn_mic
			}

			if (do_ignore_group_name)
				group_id = addGroup("group");
		}
		else
		{
			// If there is no EMDL_MICROGRAPH_NAME, then just use a single group and micrograph
			group_id = addGroup("group");
			mic_id = addMicrograph("micrograph");
		}
#ifdef DEBUG_READ
	std::cerr << "Done sorting MDimg" << std::endl;
	std::cerr << " MDimg.numberOfObjects()= " << MDimg.numberOfObjects() << std::endl;
	timer.toc(tsort);
	timer.tic(tfill);
	long nr_read = 0;
#endif
		// allocate 1 block of memory
		particles.reserve(MDimg.numberOfObjects());

		// Now Loop over all objects in the metadata file and fill the logical tree of the experiment
		long int last_oripart_idx = -1;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDimg)
		{
			// Add new micrographs or get mic_id for existing micrograph
			FileName mic_name=""; // Filename instead of string because will decompose below
			if (star_contains_micname)
			{
				long int idx = micrographs.size();
				std::string last_mic_name = (idx > 0) ? micrographs[idx-1].name : "";

				MDimg.getValue(EMDL_MICROGRAPH_NAME, mic_name);

				// All frames of a movie belong to the same micrograph
				if (is_mic_a_movie)
					mic_name = mic_name.substr(mic_name.find("@")+1);

				mic_id = -1;
				if (last_mic_name == mic_name)
				{
					// This particle belongs to the previous micrograph
					mic_id = micrographs[idx - 1].id;
				}
				else
				{
					// A new micrograph
					last_oripart_idx = ori_particles.size();
				}

				// Make a new micrograph
				if (mic_id < 0)
					mic_id = addMicrograph(mic_name);

#ifdef DEBUG_READ
				timer.tic(tgroup);
#endif

				// For example in particle_polishing the groups are not needed...
				if (!do_ignore_group_name)
				{
					std::string group_name;
					// Check whether there is a group label, if not use a group for each micrograph
					if (MDimg.containsLabel(EMDL_MLMODEL_GROUP_NAME))
					{
						MDimg.getValue(EMDL_MLMODEL_GROUP_NAME, group_name);
					}
					else
					{
						FileName fn_pre, fn_jobnr, fn_post;
						decomposePipelineFileName(mic_name, fn_pre, fn_jobnr, fn_post);
						group_name = fn_post;
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
						group_id = addGroup(group_name);
				}

#ifdef DEBUG_READ
				timer.toc(tgroup);
#endif

			}
			else
			{
				// All images belong to the same micrograph
				mic_id = 0;
				group_id = 0;
			}

			// If there is an EMDL_PARTICLE_RANDOM_SUBSET entry in the input STAR-file, then set the random_subset, otherwise use default (0)
			int my_random_subset;
			if (!MDimg.getValue(EMDL_PARTICLE_RANDOM_SUBSET, my_random_subset))
				my_random_subset = 0;

			// Create a new particle
			part_id = addParticle(group_id, mic_id, my_random_subset);

#ifdef DEBUG_READ
			timer.tic(tori);
#endif

			if (do_preread_images)
			{
				FileName fn_img;
				MDimg.getValue(EMDL_IMAGE_NAME, fn_img);
				Image<float> img;
				fn_img.decompose(dump, fn_stack);
				if (fn_stack != fn_open_stack)
				{
					hFile.openFile(fn_stack, WRITE_READONLY);
					fn_open_stack = fn_stack;
				}
				img.readFromOpenFile(fn_img, hFile, -1, false);
				img().setXmippOrigin();
				particles[part_id].img = img();
			}

			// Add this particle to an existing OriginalParticle, or create a new OriginalParticle
			std::string ori_part_name;
			long int ori_part_id = -1;

			if (MDimg.containsLabel(EMDL_PARTICLE_ORI_NAME))
				MDimg.getValue(EMDL_PARTICLE_ORI_NAME, ori_part_name);
			else
				MDimg.getValue(EMDL_IMAGE_NAME, ori_part_name);

			if (MDimg.containsLabel(EMDL_PARTICLE_ORI_NAME) && !do_ignore_original_particle_name)
			{
				// Only search ori_particles for the last (original) micrograph
				for (long int i = last_oripart_idx; i < ori_particles.size(); i++)
				{
					if (ori_particles[i].name == ori_part_name)
					{
						ori_part_id = i;
						break;
					}
				}
			}

			// If no OriginalParticles with this name was found,
			// or if no EMDL_PARTICLE_ORI_NAME in the input file, or if do_ignore_original_particle_name
			// then add a new ori_particle
			if (ori_part_id < 0)
			{
				ori_part_id = addOriginalParticle(ori_part_name, my_random_subset);
				// Also add this original_particle to an original_micrograph (only for movies)
				if (is_mic_a_movie)
				{
					micrographs[mic_id].ori_particle_ids.push_back(ori_part_id);
				}
			}
#ifdef DEBUG_READ
			timer.toc(tori);
#endif


			// Add this particle to the OriginalParticle
			std::string fnt;
			long int my_order;
			mic_name.decompose(my_order, fnt);
			(ori_particles[ori_part_id]).addParticle(part_id, my_random_subset, my_order);

			// The group number is only set upon reading: it is not read from the STAR file itself,
			// there the only thing that matters is the order of the micrograph_names
			// Write igroup+1, to start numbering at one instead of at zero
			MDimg.setValue(EMDL_MLMODEL_GROUP_NO, group_id + 1, part_id);

#ifdef DEBUG_READ
			nr_read++;
#endif
		} // end loop over all objects in MDimg

#ifdef DEBUG_READ
		timer.toc(tfill);
		timer.tic(tdef);
		std::cerr << " MDimg.lastObject()= " << MDimg.lastObject() << std::endl;
		std::cerr << " nr_read= " << nr_read << " particles.size()= " << particles.size() << " ori_particles.size()= " << ori_particles.size()  << " micrographs.size()= " << micrographs.size() << " groups.size()= " << groups.size() << std::endl;
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
	bool have_xoff = MDimg.containsLabel(EMDL_ORIENT_ORIGIN_X);
	bool have_yoff = MDimg.containsLabel(EMDL_ORIENT_ORIGIN_Y);
	bool have_zoff = MDimg.containsLabel(EMDL_ORIENT_ORIGIN_Z);
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
			MDimg.setValue(EMDL_ORIENT_ORIGIN_X, dzero);
		if (!have_yoff)
			MDimg.setValue(EMDL_ORIENT_ORIGIN_Y, dzero);
		if ( (!have_zoff) && (have_zcoord) )
			MDimg.setValue(EMDL_ORIENT_ORIGIN_Z, dzero);
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
	}

#ifdef DEBUG_READ
	timer.toc(tdef);
	std::cerr << "Done setting defaults MDimg" << std::endl;
	timer.tic(tend);
	//std::cerr << "Press any key to continue..." << std::endl;
	//std::cin >> c;
#endif

	// Also set the image_size (use the last image for that, still in fn_img)
	FileName fn_img;
	Image<RFLOAT> img;
	MDimg.getValue(EMDL_IMAGE_NAME, fn_img, MDimg.firstObject());
	if (fn_img != "")
	{
		img.read(fn_img, false); //false means read only header, skip real data
		is_3D = (ZSIZE(img()) > 1);
		int image_size = XSIZE(img());
		if (image_size != YSIZE(img()))
			REPORT_ERROR("Experiment::read: xsize != ysize: only squared images allowed");
		// Add a single object to MDexp
		MDexp.addObject();
		MDexp.setValue(EMDL_IMAGE_SIZE, image_size);
		if (ZSIZE(img()) > 1)
		{
			if (image_size != ZSIZE(img()))
				REPORT_ERROR("Experiment::read: xsize != zsize: only cubed images allowed");
			MDexp.setValue(EMDL_IMAGE_DIMENSIONALITY, 3);
		}
		else
		{
			MDexp.setValue(EMDL_IMAGE_DIMENSIONALITY, 2);
		}
	}
	else
	{
		REPORT_ERROR("There are no images read in: please check your input file...");
	}

	// Order the particles in each ori_particle (only useful for realignment of movie frames)
	orderParticlesInOriginalParticles();

#ifdef DEBUG_READ
	timer.toc(tend);
	timer.toc(tall);
	timer.printTimes(false);
	//std::cerr << "Writing out debug_data.star" << std::endl;
	//write("debug");
	//exit(0);
#endif
}

// Write to file
void Experiment::write(FileName fn_root)
{

	std::ofstream  fh;
	FileName fn_tmp = fn_root+"_data.star";
    fh.open((fn_tmp).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR( (std::string)"Experiment::write: Cannot write file: " + fn_tmp);

    // Always write MDimg
    MDimg.write(fh);

    if (nr_bodies > 1)
    {
		for (int ibody = 0; ibody < nr_bodies; ibody++)
		{
			MDbodies[ibody].write(fh);
		}
    }

	fh.close();

}
