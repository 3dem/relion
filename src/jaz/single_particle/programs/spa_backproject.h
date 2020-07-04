/***************************************************************************
 * Author: "Jasenko Zivanov"
 *
 * Derived from Sjors' SpaBackproject:
 *
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
#ifndef SPA_BACKPROJECT_PROG_H
#define SPA_BACKPROJECT_PROG_H

#include <src/backprojector.h>
#include <src/funcs.h>
#include <src/ctf.h>
#include <src/args.h>
#include <src/error.h>
#include <src/euler.h>
#include <src/time.h>
#include <src/ml_model.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/optics/dual_contrast/dual_contrast_voxel.h>


class SpaBackproject
{
		
	public:
		
		IOParser parser;
		
		FileName fn_out, fn_sel, fn_img, fn_sym, fn_sub, fn_fsc, fn_noise, image_path;
		
		MetaDataTable DF;
		ObservationModel obsModel;
		MlModel model;
		
		int 
			r_max, r_min_nn, blob_order, ref_dim, interpolator, iter,
			ctf_dim, nr_helical_asu, newbox, width_mask_edge, nr_sectors, subset, chosen_class,
			output_boxsize, verb, 
			num_threads_in, num_threads_out, num_threads_total;
		
		RFLOAT 
			blob_radius, blob_alpha, angular_error, shift_error, angpix, maxres,
			helical_rise, helical_twist,

			SNR, dual_contrast_lambda;
		
		bool 
			do_ctf, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak,				 
			do_fom_weighting, do_reconstruct_ctf, do_ewald, 
			skip_weighting, skip_mask, do_debug,
			skip_gridding, do_reconstruct_ctf2, do_reconstruct_meas, 
			is_reverse, read_weights, do_external_reconstruct,
			use_legacy_fwd_mapping, use_new_fwd_mapping,
			explicit_spreading_function;
		
		float 
			padding_factor, mask_diameter;
		
		

		std::vector<BackProjector> backprojectors;
		
		Projector projector;
		
		
		
		struct AccumulationVolume
		{
			BufferedImage<Complex> data;
			BufferedImage<RFLOAT> weight;
		};

		std::vector<AccumulationVolume> accumulation_volumes;
		std::vector<BufferedImage<RFLOAT>> multiplicities, spreading_functions;
		std::vector<BufferedImage<DualContrastVoxel<RFLOAT>>> dual_contrast_accumulation_volumes;
		
		bool compute_multiplicity, do_dual_contrast, do_isotropic_Wiener;
		int padded_box_size;
		
		
	public:
		
		SpaBackproject() { }
		
		
		void read(int argc, char **argv);
		void usage();
		void initialise();
		void determineOutputBoxSize();
		void run();
		void readDebugArrays();
		void backprojectAllParticles();
		void backprojectOneParticle(long int p, int thread_id);
		void reconstructLegacy();
		void reconstructNew();
		void reconstructDualContrast();
		
		void applyCTFPandCTFQ(
				MultidimArray<Complex> &Fin, 
				CTF &ctf, 
				FourierTransformer &transformer,
				MultidimArray<Complex> &outP, 
				MultidimArray<Complex> &outQ, 
				bool skip_mask=false);
};

#endif
