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
#ifndef SRC_RECONSTRUCTOR_H_
#define SRC_RECONSTRUCTOR_H_

#include <src/backprojector.h>
#include <src/funcs.h>
#include <src/ctf.h>
#include <src/args.h>
#include <src/error.h>
#include <src/euler.h>
#include <src/time.h>
#include <src/ml_model.h>
#include <src/jaz/single_particle/obs_model.h>

class Reconstructor
{
public:
	// I/O Parser
	IOParser parser;

	FileName fn_out, fn_sel, fn_img, fn_sym, fn_sub, fn_fsc, fn_debug, fn_noise, image_path;

	MetaDataTable DF;
	ObservationModel obsModel;
	MlModel model;

	int r_max, r_min_nn, blob_order, ref_dim, interpolator, iter,
	    debug_ori_size, debug_size,
	    ctf_dim, nr_helical_asu, newbox, width_mask_edge, nr_sectors, subset, chosen_class,
	    data_dim, output_boxsize, verb;

	RFLOAT blob_radius, blob_alpha, angular_error, shift_error, angpix, maxres,
	       helical_rise, helical_twist;

	bool do_ctf, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak,
	     do_fom_weighting, do_3d_rot, do_reconstruct_ctf, do_ewald, skip_weighting, skip_mask, do_debug,
	     do_ignore_optics, skip_subtomo_correction, normalised_subtomo, ctf3d_squared;


	bool skip_gridding, do_reconstruct_ctf2, do_reconstruct_meas, is_reverse, read_weights, do_external_reconstruct;

	float padding_factor, mask_diameter;

	// All backprojectors needed for parallel reconstruction
	BackProjector backprojector;

	// A single projector is needed for parallel reconstruction
	Projector projector;

public:
	/** Empty constructor
	 *
	 * A default Projector is created.
	 *
	 * @code
	 * Projector PPref;
	 * @endcode
	 */
	Reconstructor() { }

	// Read command line arguments
	void read(int argc, char **argv);

	// Print usage instructions
	void usage();

	// Initialise some stuff after reading
	void initialise();

	// Execute
	void run();

	// Reconstruct with debug arrays
	void readDebugArrays();

	// Loop over all particles to be back-projected
	void backproject(int rank = 0, int size = 1);

	// For parallelisation purposes
	void backprojectOneParticle(long int ipart);

	// perform the gridding reconstruction
	void reconstruct();

	void applyCTFPandCTFQ(MultidimArray<Complex> &Fin, CTF &ctf, FourierTransformer &transformer,
	                      MultidimArray<Complex> &outP, MultidimArray<Complex> &outQ, bool skip_mask=false);
};

#endif /* SRC_RECONSTRUCTOR_H_ */
