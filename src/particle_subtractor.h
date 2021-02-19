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

#ifndef PARTICLE_SUBTRACTOR_H_
#define PARTICLE_SUBTRACTOR_H_

#include "src/ml_optimiser.h"
#include "src/fftw.h"
#include "src/time.h"
#include "src/mask.h"
#include "src/funcs.h"


class ParticleSubtractor
{
public:
	// I/O Parser
	IOParser parser;

	// The ml_optimiser model, including both the particles and the references to subtract
	MlOptimiser opt;

	// FileName for the optimiser.star file, a possibly more restricted particle subset, mask and output files
	FileName fn_opt, fn_sel, fn_msk, fn_out, fn_revert;

	//Output metadatatable for each rank
	MetaDataTable MDimg_out;

	// For conventional 3D classifications/ refinements: center the subtracted particles?
	bool do_center;

	// Calculate average spectral SNRs?
	bool do_ssnr;

	// Write in half-precision 16 bit floating point numbers (MRC mode 12)
	bool write_float16;

	// Running sums of power of signal and noise for SSNR calculation (keep public for MPI access)
	MultidimArray<RFLOAT> sum_count, sum_S2, sum_N2;

	// verbosity
	int verb;

public:
	// Read command line arguments
	void read(int argc, char **argv);

	// Print usage instructions
	void usage();

	// Division of labour for parallelisation
	void divideLabour(int _rank, int _size, long int &my_first, long int &my_last);

	// Initialise some stuff after reading
	void initialise(int rank = 0, int size = 1);

	// Run over my subset of particles
	void run();

	// Revert particle subtraction
	void revert();

	// Write out a STAR file for this rank
	void saveStarFile(int myrank=0);

	// Write out the final STAR file
	void combineStarFile(int myrank=0);

	// Get name of a single subtracted particle
	FileName getParticleName(long int imgno, int myrank, int optics_group=-1);

	// subtract one particle
	void subtractOneParticle(long int part_id, long int imgno, long int counter);

private:
	// Pre-calculated rotation matrix for (0,90,0) rotation, and its transpose, for multi-body orientations
	Matrix2D<RFLOAT> A_rot90, A_rot90T;

	// Which orientations/masks to use for the subtraction, i.e. relative to which body?
	int subtract_body;

	// Output size of the subtracted particles
	int boxsize;

	// Recenter the subtracted particles at the centre-of-mass of the input mask?
	bool do_recenter_on_mask;

	// Ignore class in the particle STAR file
	bool ignore_class;

	// User-provided 3D-center (in pixels), which will be projected to provide the 2D centre of the subtracted particle images
	Matrix1D<RFLOAT> new_center;

	// For division of labour
	long int my_first_part_id, my_last_part_id;

	// For MPI parallelisation
	int rank, size;

	// image to particle mapping
	std::vector<long int> nr_particles_in_optics_group;
	std::map<long int, FileName> imgno_to_filename;
};

#endif /* PARTICLE_SUBTRACTOR_H_ */
