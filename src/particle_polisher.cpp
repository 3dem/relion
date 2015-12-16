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

#include "src/particle_polisher.h"
//#define DEBUG_TILT

void ParticlePolisher::read(int argc, char **argv)
{

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "STAR file with the aligned movie frames, e.g. run1_ct25_data.star");
	fn_out = parser.getOption("--o", "Output rootname", "shiny");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "1"));
	running_average_width = textToInteger(parser.getOption("--movie_frames_running_avg", "Number of movie frames in each running average", "5"));
	do_start_all_over = parser.checkOption("--dont_read_old_files", "Do not read intermediate results from disc, but re-do all calculations from scratch");

	int fit_section = parser.addSection("Beam-induced movement fitting options");
	sigma_neighbour_distance = textToFloat(parser.getOption("--sigma_nb", "Standard deviation for a Gaussian weight on the particle distance", "100."));
	if (parser.checkOption("--log_fit",       "Fit line on a logarithmic time-scale      (default=linear fit)"))
		fitting_mode = LOGARITHMIC_FIT;
	else if (parser.checkOption("--sqrt_fit", "Fit line on a square-root time-scale      (default=linear fit)"))
		fitting_mode = SQRT_FIT;
	else if (parser.checkOption("--no_fit",   "Do not fit any function through movements (default=linear fit)"))
		fitting_mode = NO_FIT;
	else
		fitting_mode = LINEAR_FIT;

	int post_section = parser.addSection("Per-frame B-factor estimation options");
	do_weighting = !parser.checkOption("--skip_bfactor_weighting", "Don't perform B-factor weighting");
	frame_running_average = textToInteger(parser.getOption("--bfactor_running_avg", "Number of movie frames to join in B-factor determination", "1"));
	perframe_highres = textToFloat(parser.getOption("--perframe_highres", "Highest resolution (in A) to include in the per-frame reconstructions", "-1."));
	fit_minres = textToFloat(parser.getOption("--autob_lowres", "Lowest resolution (in A) to include in fitting of the B-factor", "20."));
	fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
	// Sep24,2015 - Shaoda, Helical reconstruction
	nr_helical_asu = textToInteger(parser.getOption("--h_asu", "Number of new helical asymmetric units (asu) per box (1 means no helical symmetry is present)", "1"));
	helical_twist = textToFloat(parser.getOption("--h_twist", "Helical twist (in degrees, positive values for right-handedness)", "0."));
	helical_rise = textToFloat(parser.getOption("--h_rise", "Helical rise (in Angstroms)", "0."));

	fn_mask = parser.getOption("--mask", "Postprocessing mask for B-factor determination of per-frame reconstructions (1=protein, 0=solvent, all values in range [0,1])", "");

	int ctf_section = parser.addSection("CTF options");
   	do_ctf = !parser.checkOption("--no_ctf", "Don't apply CTF correction");
	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");
	ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Images have been phase flipped");
	only_flip_phases = parser.checkOption("--only_flip_phases", "Do not correct CTF-amplitudes, only flip phases");
	defocus_shift_max = textToFloat(parser.getOption("--defocus_shift_max", "Maximum shift in defocus (in A) to search for each particle", "0"));
	defocus_shift_step = textToFloat(parser.getOption("--defocus_shift_step", "Maximum shift in defocus (in A) to search for each particle", "100"));

	int norm_section = parser.addSection("Normalisation options");
	do_normalise = !parser.checkOption("--skip_normalise", "Do not normalise the polsihed particles?");
	bg_radius =  textToInteger(parser.getOption("--bg_radius", "Radius of the circular mask that will be used to define the background area (in pixels)", "-1"));
	do_ramp = !parser.checkOption("--no_ramp", "Just subtract the background mean in the normalisation, instead of subtracting a fitted ramping background. ");
	white_dust_stddev = textToFloat(parser.getOption("--white_dust", "Sigma-values above which white dust will be removed (negative value means no dust removal)","-1"));
	black_dust_stddev = textToFloat(parser.getOption("--black_dust", "Sigma-values above which black dust will be removed (negative value means no dust removal)","-1"));

	int beamtilt_section = parser.addSection("Beamtilt refinement options");
	beamtilt_max = textToFloat(parser.getOption("--beamtilt_max", "Maximum beamtilt (in mrad) to search", "0."));
	beamtilt_step = textToFloat(parser.getOption("--beamtilt_step", "Step-size for beamtilt searches (in mrad)", "0.2"));
	minres_beamtilt = textToFloat(parser.getOption("--minres_beamtilt", "Lowest resolution to include in beam-tilt correction (in A)", "6"));

	int out_section = parser.addSection("Polished particles output options");
	first_frame = textToInteger(parser.getOption("--first_frame", "First frame to include in the polished particles", "1"));
	last_frame = textToInteger(parser.getOption("--last_frame", "First frame to include in the polished particles (default is all)", "-1"));
	step_frame  = textToInteger(parser.getOption("--avg_movie_frames", "Give the same value as used in particle extraction (usually just 1)", "1"));

	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void ParticlePolisher::usage()
{
	parser.writeUsage(std::cerr);
}

void ParticlePolisher::initialise()
{
#ifndef DEBUG_TILT
    if (verb > 0)
    	std::cout << " + Reading the input STAR file ... " << std::endl;
    exp_model.read(fn_in, false, true); // false means NOT do_ignore_particle_name here, true means DO_ignore_group_name...

	// Pre-size the fitted_movements array (for parallelised output...)
	long int total_nr_images = exp_model.numberOfParticles();
	fitted_movements.resize(total_nr_images, 2);

	last_frame = (last_frame < 0) ? exp_model.ori_particles[0].particles_id.size() : last_frame;
	perframe_bfactors.initZeros( 3 * (last_frame - first_frame + 1)); // now store bfactor, offset and corr_coeff
#endif

	// Get the pixel size from the input STAR file
	if (exp_model.MDimg.containsLabel(EMDL_CTF_MAGNIFICATION) && exp_model.MDimg.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
	{
		RFLOAT mag, dstep;
		exp_model.MDimg.getValue(EMDL_CTF_MAGNIFICATION, mag);
		exp_model.MDimg.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
		angpix = 10000. * dstep / mag;
		if (verb > 0)
			std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << std::endl;
	}

    if (do_weighting)
    {

    	if (fn_mask == "")
            REPORT_ERROR("When doing B-factor weighting, you have to provide a mask!");

    	// Read in the Imask
    	Imask.read(fn_mask);
    	ori_size = XSIZE(Imask());
    }

	if (do_normalise && bg_radius < 0)
		REPORT_ERROR("ERROR: please provide a radius for a circle that defines the background area when normalising...");

	// Sep24,2015 - Shaoda, Helical reconstruction
	if (nr_helical_asu > 1)
	{
		if ( (fabs(helical_twist) < 0.01) || (fabs(helical_twist) > 179.99) || (angpix < 0.001) || ((helical_rise / angpix) < 0.001))
			REPORT_ERROR("ERROR: Invalid helical twist or rise!");
	}
}

// Fit the beam-induced translations for all average micrographs
void ParticlePolisher::fitMovementsAllMicrographs()
{

	// Loop over all average micrographs
	int barstep;
	int my_nr_micrographs = exp_model.average_micrographs.size();
	if (verb > 0)
	{
		std::cout << " + Fitting straight paths for beam-induced movements in all micrographs ... " << std::endl;
		init_progress_bar(my_nr_micrographs);
		barstep = XMIPP_MAX(1, my_nr_micrographs/ 60);
	}

    for (long int i = 0; i < my_nr_micrographs; i++)
	{
    	if (verb > 0 && i % barstep == 0)
			progress_bar(i);

		fitMovementsOneMicrograph(i);
	}

    // Set the fitted movements in the xoff and yoff columns of the exp_model.MDimg
    for (long int ipart = 0; ipart < exp_model.numberOfParticles(); ipart++)

	{
		long int part_id = exp_model.particles[ipart].id;
		RFLOAT xoff = DIRECT_A2D_ELEM(fitted_movements, part_id, 0);
		RFLOAT yoff = DIRECT_A2D_ELEM(fitted_movements, part_id, 1);
		exp_model.MDimg.setValue(EMDL_ORIENT_ORIGIN_X, xoff, part_id);
		exp_model.MDimg.setValue(EMDL_ORIENT_ORIGIN_Y, yoff, part_id);
	}

    if (verb > 0)
	{
		progress_bar(my_nr_micrographs);
	}

	// Write out the STAR file with all the fitted movements
	FileName fn_tmp = fn_in.withoutExtension() + "_" + fn_out + ".star";
	exp_model.MDimg.write(fn_tmp);
	std::cout << " + Written out all fitted movements in STAR file: " << fn_tmp << std::endl;


}



void ParticlePolisher::fitMovementsOneMicrograph(long int imic)
{

	std::vector<RFLOAT> x_pick, y_pick, x_off_prior, y_off_prior, x_start, x_end, dummy; // X and Y-coordinates for the average particles in the micrograph
	std::vector< std::vector<RFLOAT> > x_off, y_off; // X and Y shifts w.r.t. the average for each movie frame
	RFLOAT x_pick_p, y_pick_p, x_off_p, y_off_p, x_off_prior_p, y_off_prior_p;


	// Loop over all original_particles in this average_micrograph
	// And fill the x_pick, y_;pick, x_off and y_off vectors
	for (long int ipar = 0; ipar < (exp_model.average_micrographs[imic]).ori_particles_id.size(); ipar++)
	{
		long int ori_part_id = exp_model.average_micrographs[imic].ori_particles_id[ipar];

		x_off.push_back(dummy);
		y_off.push_back(dummy);
		bool is_first = true;
		for (long int i_frame = 0; i_frame < exp_model.ori_particles[ori_part_id].particles_id.size(); i_frame++ )
		{
			long int part_id = exp_model.ori_particles[ori_part_id].particles_id[i_frame];

			exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, x_off_p, part_id);
			exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, y_off_p, part_id);
			exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_X_PRIOR, x_off_prior_p, part_id);
			exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y_PRIOR, y_off_prior_p, part_id);

			// Store the value of the picked coordinates for the first movie frame
			if (is_first)
			{
				is_first = false;
				exp_model.MDimg.getValue(EMDL_IMAGE_COORD_X, x_pick_p, part_id);
				exp_model.MDimg.getValue(EMDL_IMAGE_COORD_Y, y_pick_p, part_id);
				x_pick.push_back(x_pick_p);
				y_pick.push_back(y_pick_p);
				x_off_prior.push_back(x_off_prior_p);
				y_off_prior.push_back(y_off_prior_p);
			}
			else
			{
				if (ABS(x_off_prior_p - x_off_prior[ipar]) > 0.01 || ABS(y_off_prior_p - y_off_prior[ipar]) > 0.01)
					REPORT_ERROR("ParticlePolisher::processFramesOneMicrograph ERROR: unequal priors on x and y for different movie frames!");
			}

			// Store the offsets for all movie frames, relative to the prior (to get the movements of the frames)
			x_off[ipar].push_back(x_off_p);
			y_off[ipar].push_back(y_off_p);
		}
	}

	// Now do the actual fitting
	RFLOAT gauss_const = 1. / sqrt(2 * PI * sigma_neighbour_distance * sigma_neighbour_distance);
	RFLOAT min2sigma2 = - 2. * sigma_neighbour_distance * sigma_neighbour_distance;
	// Loop over all ori_particles
	for (long int ipar = 0; ipar < (exp_model.average_micrographs[imic]).ori_particles_id.size(); ipar++)
	{
		long int ori_part_id = exp_model.average_micrographs[imic].ori_particles_id[ipar];

		// Sjors 14sep2015: bug reported by Kailu Yang
		RFLOAT my_pick_x = x_pick[ipar] - x_off_prior[ipar];
		RFLOAT my_pick_y = y_pick[ipar] - y_off_prior[ipar];

		fit_point2D      onepoint;
		std::vector<fit_point2D> points_x, points_y;

		// Loop over all other ori_particles on this micrograph and determine weight for contribution to this ori_particle
		for (long int ii = 0; ii < x_pick.size(); ii++)
		{

			// Sjors 14sep2015: bug reported by Kailu Yang
			RFLOAT nb_pick_x = x_pick[ii] - x_off_prior[ii]; // add prior to center at average position
			RFLOAT nb_pick_y = y_pick[ii] - y_off_prior[ii]; // add prior to center at average position
			RFLOAT dist2 = (nb_pick_x - my_pick_x) * (nb_pick_x - my_pick_x) + (nb_pick_y - my_pick_y) * (nb_pick_y - my_pick_y);
			RFLOAT weight;
			if (ABS(min2sigma2) < 0.01)
				weight = (dist2 < 0.01) ? 1. : 0.;
			else
				weight = gauss_const * exp(dist2 / min2sigma2);

			if (weight > 0.00001) // ignore very small weights
			{
				int first_frame_fit = running_average_width / 2;
				int last_frame_fit = x_off[ii].size() - (running_average_width / 2);
				// only fit using those frames that contribute fully to the running average, i.e. exclude the first and last few frames...
				for (long int i_frame = first_frame_fit; i_frame < last_frame_fit; i_frame++)
				{
					if (fitting_mode == LINEAR_FIT)
						onepoint.x = i_frame + 1; // linear movement with time.... that may not be true!
					else if (fitting_mode == LOGARITHMIC_FIT)
						onepoint.x = log(i_frame + 1); // logarithmic movement with time: better?!
					else if (fitting_mode == SQRT_FIT)
						onepoint.x = sqrt(i_frame + 1); // logarithmic movement with time: better?!
					onepoint.w = weight;
					onepoint.y = x_off[ii][i_frame] - x_off_prior[ii]; // subtract prior to center movement of average position
					points_x.push_back(onepoint);
					onepoint.y = y_off[ii][i_frame] - y_off_prior[ii]; // subtract prior to center movement of average position
					points_y.push_back(onepoint);
				}
			}
		}

		RFLOAT slope_x, intercept_x, ccf_x;
		RFLOAT slope_y, intercept_y, ccf_y;
		fitStraightLine(points_x, slope_x, intercept_x, ccf_x);
		fitStraightLine(points_y, slope_y, intercept_y, ccf_y);

		// Now set the interpolated values
		// Do this for all frames! (whereas fitting was done with excluding the running_average_width halves from both ends)
		// Rather than setting them straight back into the exp_model.MDimg, use the fitted_movements array
		// In the parallel version, this will be required for gathering the results from all nodes
		for (long int i_frame = 0; i_frame < x_off[ipar].size(); i_frame++)
		{
			if (fitting_mode == LINEAR_FIT)
			{
				x_off_p = slope_x * (i_frame + 1) + intercept_x;
				y_off_p = slope_y * (i_frame + 1) + intercept_y;
			}
			else if (fitting_mode == LOGARITHMIC_FIT)
			{
				x_off_p = slope_x * log(i_frame + 1) + intercept_x;
				y_off_p = slope_y * log(i_frame + 1) + intercept_y;
			}
			else if (fitting_mode == SQRT_FIT)
			{
				x_off_p = slope_x * sqrt(i_frame + 1) + intercept_x;
				y_off_p = slope_y * sqrt(i_frame + 1) + intercept_y;
			}

			long int part_id = exp_model.ori_particles[ori_part_id].particles_id[i_frame];
			A2D_ELEM(fitted_movements, part_id, 0) = x_off_p + x_off_prior[ipar];
			A2D_ELEM(fitted_movements, part_id, 1) = y_off_p + y_off_prior[ipar];
		}
	}

}

void ParticlePolisher::calculateAllSingleFrameReconstructionsAndBfactors()
{

	FileName fn_star = fn_in.withoutExtension() + "_" + fn_out + "_bfactors.star";
	if (!do_start_all_over && readStarFileBfactors(fn_star))
	{
		if (verb > 0)
			std::cout << " + " << fn_star << " already exists: skipping calculation of per-frame B-factors." <<std::endl;
		return;
	}

	RFLOAT bfactor, offset, corr_coeff;

	// Loop over all frames to be included in the reconstruction
	int my_nr_frames = last_frame - first_frame + 1;
	int barstep;
	if (verb > 0)
	{
		std::cout << " + Calculating per-frame reconstructions ... " << std::endl;
		init_progress_bar(my_nr_frames);
		barstep = XMIPP_MAX(1, my_nr_frames/ 60);
	}

    for (long int i = first_frame; i <= last_frame; i++)
	{

    	calculateSingleFrameReconstruction(i, 1);
    	calculateSingleFrameReconstruction(i, 2);

    	if (verb > 0 && i % barstep == 0)
    		progress_bar(i);

	}

    // Also calculate the average of all single-frames for both halves
    calculateAverageAllSingleFrameReconstructions(1);
    calculateAverageAllSingleFrameReconstructions(2);

    // Now calculate the FSC between the average of all single-frame reconstructions
    calculateBfactorSingleFrameReconstruction(-1, bfactor, offset, corr_coeff);

	if (verb > 0)
	{
		std::cout << " + Calculating per-frame B-factors and offsets ... " << std::endl;
		init_progress_bar(my_nr_frames);
		barstep = XMIPP_MAX(1, my_nr_frames/ 60);
	}
    for (long int i = first_frame; i <= last_frame; i++)
	{

       	calculateBfactorSingleFrameReconstruction(i, bfactor, offset, corr_coeff);
       	int iframe = i - first_frame;
       	DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 0) = bfactor;
       	DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 1) = offset;
       	DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 2) = corr_coeff;

       	if (verb > 0 && i % barstep == 0)
    		progress_bar(i);
	}


    if (verb > 0)
	{
		progress_bar(my_nr_frames);
	}

    // Write STAR file
    writeStarFileBfactors(fn_star);


    // Also write a STAR file with the relative contributions of each frame to all frequencies
    fn_star = fn_in.withoutExtension() + "_" + fn_out + "_relweights.star";
    writeStarFileRelativeWeights(fn_star);


}

bool ParticlePolisher::readStarFileBfactors(FileName fn_star)
{
	if (exists(fn_star))
	{
		MetaDataTable MD;
		MD.read(fn_star);
		int i = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			MD.getValue(EMDL_POSTPROCESS_BFACTOR, DIRECT_A1D_ELEM(perframe_bfactors, i * 3 + 0) );
			MD.getValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, DIRECT_A1D_ELEM(perframe_bfactors, i * 3 + 1) );
			i++;
		}
		return true;
	}
	else
	{
		return false;
	}

}

void ParticlePolisher::writeStarFileBfactors(FileName fn_star)
{
	MetaDataTable MDout;
	MDout.setName("perframe_bfactors");

	for (int i = 0; i < XSIZE(perframe_bfactors)/3; i++ )
	{
		MDout.addObject();
		MDout.setValue(EMDL_IMAGE_FRAME_NR, first_frame + i);
		MDout.setValue(EMDL_POSTPROCESS_BFACTOR, DIRECT_A1D_ELEM(perframe_bfactors, i * 3 + 0) );
		MDout.setValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, DIRECT_A1D_ELEM(perframe_bfactors, i * 3 + 1) );
	}

	MDout.write(fn_star);
}

void ParticlePolisher::writeStarFileRelativeWeights(FileName fn_star)
{
	std::ofstream  fh;
    fh.open((fn_star).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR( (std::string)"ParticlePolisher::writeStarFileRelativeWeights: Cannot write file: " + fn_star);

    // First pre-calculate the sum of all weights at every frequency
    MultidimArray<RFLOAT> sumweight_per_shell(ori_size/2), cumulative_relweight_per_shell(ori_size/2);
    sumweight_per_shell.initZeros();
    for (int iframe = 0; iframe < XSIZE(perframe_bfactors)/3; iframe++ )
   	{
    	RFLOAT bfactor = DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 0);
    	RFLOAT offset = DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 1);
		for (int i = 1; i < ori_size/2; i++) // ignore origin
		{
			RFLOAT res = (RFLOAT)i / (ori_size * angpix); // resolution in 1/A
			RFLOAT w = exp( (bfactor / 4.)  * res * res  + offset);

			DIRECT_A1D_ELEM(sumweight_per_shell, i) += w;
		}
   	}

    // Now calculate the relative weights and their cumulative curves
    cumulative_relweight_per_shell.initZeros();
    for (int iframe = 0; iframe < XSIZE(perframe_bfactors)/3; iframe++ )
    {

    	RFLOAT bfactor = DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 0);
    	RFLOAT offset = DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 1);

    	MetaDataTable MDout;
		std::string fn_table = "relative_weights_frame_" + integerToString(first_frame + iframe);
		MDout.setName(fn_table);

		for (int i = 1; i < ori_size/2; i++) // ignore origin
		{
			RFLOAT res = (RFLOAT)i / (ori_size * angpix); // resolution in 1/A
			RFLOAT w = exp( (bfactor / 4.)  * res * res  + offset);

			w /= DIRECT_A1D_ELEM(sumweight_per_shell, i); // from absolute to relative weight
			DIRECT_A1D_ELEM(cumulative_relweight_per_shell, i) += w;

			MDout.addObject();
			MDout.setValue(EMDL_RESOLUTION, res); //res in 1/A
			MDout.setValue(EMDL_PERFRAME_RELATIVE_WEIGHT, w);
			MDout.setValue(EMDL_PERFRAME_CUMULATIVE_WEIGHT, DIRECT_A1D_ELEM(cumulative_relweight_per_shell, i));

		}
		MDout.write(fh);

	}

    fh.close();
}

void ParticlePolisher::calculateAverageAllSingleFrameReconstructions(int this_half)
{

	FileName fn_sum;
	fn_sum = fn_in.withoutExtension() + "_" + fn_out + "_avgframes_half" + integerToString(this_half) + "_class001_unfil.mrc";

	if (!do_start_all_over && exists(fn_sum))
	{
		if (verb > 0)
			std::cout << std::endl << " + " << fn_sum << " already exists: skipping calculation average of all frames." << std::endl;
		return;
	}

	Image<RFLOAT> Isum, Ione;
	for (long int this_frame = first_frame; this_frame <= last_frame; this_frame++)
	{
    	FileName fn_vol;
    	fn_vol.compose(fn_in.withoutExtension() + "_" + fn_out + "_frame", this_frame, "", 3);
    	fn_vol += "_half" + integerToString(this_half) + "_class001_unfil.mrc";

    	if (this_frame == first_frame)
    	{
    		Isum.read(fn_vol);
    	}
    	else
    	{
    		Ione.read(fn_vol);
    		Isum() += Ione();
    	}
	}

	// Write the average map to disc
	Isum.write(fn_sum);

}

void ParticlePolisher::calculateSingleFrameReconstruction(int this_frame, int this_half)
{

	FileName fn_vol;
	fn_vol.compose(fn_in.withoutExtension() + "_" + fn_out + "_frame", this_frame, "", 3);
	fn_vol += "_half" + integerToString(this_half) + "_class001_unfil.mrc";
	if (!do_start_all_over && exists(fn_vol))
	{
		if (verb > 0)
			std::cout << std::endl << " + " << fn_vol << " already exists: skipping per-frame reconstruction." << std::endl;
		return;
	}



	int image_size, current_size;
	// get image size from metadatatable
	exp_model.MDexp.getValue(EMDL_IMAGE_SIZE, image_size);
	// Set current_size
	if (perframe_highres > 0.)
		current_size = 2 * ROUND(image_size * angpix / perframe_highres);
	else
		current_size = image_size;
	BackProjector backprojector(image_size, 3, fn_sym);
	backprojector.initZeros(current_size);

	CTF ctf;
	std::string dum;
	Matrix2D<RFLOAT> A3D;
	MultidimArray<Complex > Faux, F2D, F2Dp, Fsub;
	MultidimArray<RFLOAT> Fweight, Fctf;
	Image<RFLOAT> img, vol;
	FourierTransformer transformer;
	RFLOAT xtrans, ytrans;
	RFLOAT rot, tilt, psi;
	int i_half;
	long int i_frame;
	FileName fn_img, fn_mic;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(exp_model.MDimg)
	{

		exp_model.MDimg.getValue(EMDL_PARTICLE_RANDOM_SUBSET, i_half);
		exp_model.MDimg.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
		fn_mic.decompose(i_frame, dum);
		// If we used --avg_movie_frame upon extraction, now count less
		// e.g. instead of counting 4, 8, 12 count 1, 2, 3
		i_frame /= step_frame;

		// TODO!! Make a running average window here
		if (ABS(i_frame - this_frame) <= frame_running_average/2 && i_half == this_half)
		{
			exp_model.MDimg.getValue(EMDL_IMAGE_NAME, fn_img);
			img.read(fn_img);
			CenterFFT(img(), true);
			transformer.FourierTransform(img(), F2Dp);
			windowFourierTransform(F2Dp, F2D, current_size);

			// Use the prior-angles, as these were determined from the average particles
			// The individual-frame-determined angles would be too noisy....
			exp_model.MDimg.getValue(EMDL_ORIENT_ROT_PRIOR, rot);
			exp_model.MDimg.getValue(EMDL_ORIENT_TILT_PRIOR, tilt);
			exp_model.MDimg.getValue(EMDL_ORIENT_PSI_PRIOR, psi);
			Euler_angles2matrix(rot, tilt, psi, A3D);

			// Translations
			xtrans=ytrans=0.;
			exp_model.MDimg.getValue( EMDL_ORIENT_ORIGIN_X, xtrans);
			exp_model.MDimg.getValue( EMDL_ORIENT_ORIGIN_Y, ytrans);
			if (ABS(xtrans) > 0. || ABS(ytrans) > 0. )
				shiftImageInFourierTransform(F2D, F2D, image_size, xtrans, ytrans );

			// CTF
			Fctf.resize(F2D);
			Fctf.initConstant(1.);
			if (do_ctf)
			{
				ctf.read(exp_model.MDimg, exp_model.MDimg);
				ctf.getFftwImage(Fctf, image_size, image_size, angpix, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
				{
					DIRECT_MULTIDIM_ELEM(F2D, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
				}
			}

			backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
		}

	}

	backprojector.symmetrise(nr_helical_asu, helical_twist, helical_rise / angpix);

	// Now do the reconstruction
	MultidimArray<RFLOAT> dummy;
	backprojector.reconstruct(vol(), 10, false, 1., dummy, dummy, dummy, dummy);

	vol.write(fn_vol);

}

void ParticlePolisher::calculateBfactorSingleFrameReconstruction(int this_frame, RFLOAT &bfactor, RFLOAT &offset, RFLOAT &corr_coeff)
{

	if (NZYXSIZE(Imask()) == 0)
		REPORT_ERROR("ParticlePolisher::calculateBfactorSingleFrameReconstruction BUG: first read Imask, then call this function.");

	if (this_frame > 0 && XSIZE(fsc_average) == 0)
		REPORT_ERROR("ParticlePolisher::calculateBfactorSingleFrameReconstruction BUG: first call this function with this_frame < 0 to calculate FSC_average.");


		// Read in the 2 half-reconstructions and calculate their FSC
	FileName fn_root_half;
	// Make sure that the first call to this function is with this_frame < 0!!!
	if (this_frame < 0)
		fn_root_half = fn_in.withoutExtension() + "_" + fn_out+"_avgframes";
	else
		fn_root_half.compose(fn_in.withoutExtension() + "_" + fn_out+"_frame",this_frame,"", 3);

	FileName fn_half1, fn_half2;
	Image<RFLOAT> I1, I2;
	MultidimArray<RFLOAT> fsc_frame;
	I1.read(fn_root_half + "_half1_class001_unfil.mrc");
	I2.read(fn_root_half + "_half2_class001_unfil.mrc");

	// Mask both maps
	I1() *= Imask();
	I2() *= Imask();
	getFSC(I1(), I2(), fsc_frame);

	if (this_frame < 0)
	{
		fsc_average = fsc_frame;
		return;
	}
	else
	{
		// Now use relative ratio of signal amplitudes w.r.t. the average of all single-frame reconstructions
		// SSNR= tau^2/sigma^2 = FSC / (1 - FSC)
		// tau_frame / tau_avg = tau_f / tau_a = sqrt (FSC_f / (1 - FSC_f)) / sqrt (FSC_a / (1 - FSC_a))
		// = sqrt( {FSC_f / (1 - FSC_f)} * {(1 - FSC_a) / FSC_a} )
		// = sqrt( (FSC_f - FSC_f * FSC_a) / (FSC_a - FSC_f * FSC_a)  )
		// Then make a Guinier plot of that: ln(tau) versus 1/d^2 and
		// fit a line through all points where FSC_f < 1 && FSC_f > 0.143
		// Store the bfactor (4*slope) AND the offset of that line


		MetaDataTable MDout;
		MDout.setName("relative_guinier");

		fit_point2D      onepoint;
		std::vector<fit_point2D> guinier;
		for (int i = 1; i < XSIZE(fsc_frame); i++) // ignore origin
		{
			RFLOAT res = (RFLOAT)i / (XSIZE(I1()) * angpix); // resolution in 1/A
			RFLOAT resang = 1. / res;
			RFLOAT res2 = res*res;

			RFLOAT fsc_f = DIRECT_A1D_ELEM(fsc_frame, i);
			RFLOAT fsc_a = DIRECT_A1D_ELEM(fsc_average, i);

			if (resang < fit_minres && resang > perframe_highres && fsc_f < 1 && fsc_a < 1 && fsc_f > 0.143 && fsc_a > 0.143)
			{
				// ln(tau_f / tau_a)
				// I could have calculated the same using: ln(tau_f / tau_a)   = ln(tau_f) - ln(tau_a)
				// where tau_f = sqrt (FSC_f / (1 - FSC_f)) and tau_a = sqrt (FSC_a / (1 - FSC_a))
				// This is numerically identical
				RFLOAT logreltau = log( sqrt( (fsc_f - fsc_f * fsc_a) / (fsc_a - fsc_f * fsc_a) ) );

				onepoint.x = res2;
				onepoint.y = logreltau;
				onepoint.w = 1.;
				guinier.push_back(onepoint);

				MDout.addObject();
				MDout.setValue(EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, res2);
				MDout.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_IN, logreltau);
			}
		}
		MDout.write(fn_root_half + "_guinier.star");

		// Now do the fit
		fitStraightLine(guinier, bfactor, offset, corr_coeff);
		// this is the B-factor relative to the average from all single-frame reconstructions!
		// in this case: positive values mean BETTER than average, and thus HIGHER WEIGHTS!
		bfactor *= 4.;
	}

}

void ParticlePolisher::polishParticlesAllMicrographs()
{

	if (!do_start_all_over && exists(fn_out + ".star"))
	{
		if (verb > 0)
			std::cout << std::endl << " + " << fn_out << ".star already exists: skipping polishing of the particles." << std::endl;

		return;
	}

	// Loop over all average micrographs
	int barstep;
	int my_nr_micrographs = exp_model.average_micrographs.size();
	if (verb > 0)
	{
		std::cout << " + Write out polished particles for all " << my_nr_micrographs << " micrographs ... " << std::endl;
		init_progress_bar(my_nr_micrographs);
		barstep = XMIPP_MAX(1, my_nr_micrographs/ 60);
	}

    for (long int i = 0; i < my_nr_micrographs; i++)
	{
    	if (verb > 0 && i % barstep == 0)
			progress_bar(i);

		polishParticlesOneMicrograph(i);
	}

    writeStarFilePolishedParticles();

    if (verb > 0)
	{
		progress_bar(my_nr_micrographs);
	}

}

void ParticlePolisher::writeStarFilePolishedParticles()
{

	// Also write the STAR file with the MetaData of the polished particles
	// Loop over all original_particles in this average_micrograph

	MDshiny.clear();
	for (long int imic = 0; imic < exp_model.average_micrographs.size(); imic++)
	{
		for (long int ipar = 0; ipar < exp_model.average_micrographs[imic].ori_particles_id.size(); ipar++)
		{
			long int ori_part_id = exp_model.average_micrographs[imic].ori_particles_id[ipar];
			long int part_id = exp_model.ori_particles[ori_part_id].particles_id[first_frame - 1];

			// Get the corresponding line from the input STAR file
			MDshiny.addObject(exp_model.MDimg.getObject(part_id));

			// The orientations and offsets in this line of the metadatatable come only from a single frame!!!
			// Use the ones stored in the prior instead.
			RFLOAT xoff, yoff, rot, tilt, psi;
			MDshiny.getValue(EMDL_ORIENT_ORIGIN_X_PRIOR, xoff);
			MDshiny.getValue(EMDL_ORIENT_ORIGIN_Y_PRIOR, yoff);
			MDshiny.getValue(EMDL_ORIENT_ROT_PRIOR, rot);
			MDshiny.getValue(EMDL_ORIENT_TILT_PRIOR, tilt);
			MDshiny.getValue(EMDL_ORIENT_PSI_PRIOR, psi);
			MDshiny.setValue(EMDL_ORIENT_ORIGIN_X, xoff);
			MDshiny.setValue(EMDL_ORIENT_ORIGIN_Y, yoff);
			MDshiny.setValue(EMDL_ORIENT_ROT, rot);
			MDshiny.setValue(EMDL_ORIENT_TILT, tilt);
			MDshiny.setValue(EMDL_ORIENT_PSI, psi);

			// also get the "@" out of the micrograph name!
			FileName fn_tmp;
			MDshiny.getValue(EMDL_MICROGRAPH_NAME, fn_tmp);
			fn_tmp = fn_tmp.afterFirstOf("@");
			MDshiny.setValue(EMDL_MICROGRAPH_NAME, fn_tmp);

			// Also change this particle's image_name
			FileName fn_part, fn_img;
			exp_model.MDimg.getValue(EMDL_PARTICLE_ORI_NAME, fn_part, part_id);
			fn_part = fn_part.withoutExtension();
			std::string mic_name;
			long int nr;
			fn_part.decompose(nr, mic_name);
			fn_part = mic_name + "_" + fn_out + ".mrcs";
			fn_img.compose(ipar + 1, fn_part);

			MDshiny.setValue(EMDL_IMAGE_NAME, fn_img);
		}
	}

	// Deactivate PRIORs so they are not written out
	MDshiny.deactivateLabel(EMDL_ORIENT_ORIGIN_X_PRIOR);
	MDshiny.deactivateLabel(EMDL_ORIENT_ORIGIN_Y_PRIOR);
	MDshiny.deactivateLabel(EMDL_ORIENT_ROT_PRIOR);
	MDshiny.deactivateLabel(EMDL_ORIENT_TILT_PRIOR);
	MDshiny.deactivateLabel(EMDL_ORIENT_PSI_PRIOR);

	// Write output metadatatable
	MDshiny.write(fn_out + ".star");
	std::cout << " + Written out all polished particles and their corresponding STAR file: " << fn_out << ".star" << std::endl;

}


void ParticlePolisher::polishParticlesOneMicrograph(long int imic)
{

	// Then read in all individual movie frames, apply frame x,y-movements as phase shifts and calculate polished (shiny) particles
	// as average of the re-aligned frames
	RFLOAT x_off_p, y_off_p, x_off_prior_p, y_off_prior_p;
	FileName fn_img, fn_part;
	Image<RFLOAT> img;
	FourierTransformer transformer;
	MultidimArray<Complex > Fimg, Fwsum;
	MultidimArray<RFLOAT> Fsumw;
	RFLOAT xtrans, ytrans;
	RFLOAT all_minval = 99999., all_maxval = -99999., all_avg = 0., all_stddev = 0.;

	// Loop over all original_particles in this average_micrograph
	for (long int ipar = 0; ipar < exp_model.average_micrographs[imic].ori_particles_id.size(); ipar++)
	{
		long int ori_part_id = exp_model.average_micrographs[imic].ori_particles_id[ipar];

		// Loop over all frames for motion corrections and possibly dose-dependent weighting
		for (long int i_frame = first_frame; i_frame <= last_frame; i_frame++ )
		{
			long int part_id = exp_model.ori_particles[ori_part_id].particles_id[i_frame - 1]; // start counting frames at 0, not 1

			exp_model.MDimg.getValue(EMDL_IMAGE_NAME, fn_img, part_id);
			exp_model.MDimg.getValue(EMDL_PARTICLE_ORI_NAME, fn_part, part_id);
			exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_X_PRIOR, x_off_prior_p, part_id);
			exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y_PRIOR, y_off_prior_p, part_id);
			if (fitting_mode != NO_FIT)
			{
				x_off_p = A2D_ELEM(fitted_movements, part_id, 0);
				y_off_p = A2D_ELEM(fitted_movements, part_id, 1);
			}
			else
			{
				exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, x_off_p, part_id);
				exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, y_off_p, part_id);
			}

			img.read(fn_img);
			RFLOAT ori_size= XSIZE(img());

			// Get the image shifts relative to the prior
			xtrans = x_off_p - x_off_prior_p;
			ytrans = y_off_p - y_off_prior_p;

#ifdef DEBUG
			if (fn_part =="000001@Particles/Micrographs/Falcon_2012_06_13-01_05_13_0_rh15particles.mrcs")
			{
				RFLOAT xp, yp;
				exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, xp, part_id);
				exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, yp, part_id);
				std::cerr << " x_off_prior_p= " << x_off_prior_p << " x_off_p= " << x_off_p << " xp= " << xp << std::endl;
				std::cerr << " y_off_prior_p= " << y_off_prior_p << " y_off_p= " << y_off_p << " yp= " << yp << std::endl;
				std::cerr << " iframe= " << i_frame << " XX(trans)= " << XX(trans) << " YY(trans)= " << YY(trans) << std::endl;
			}
#endif

			// Apply the phase shifts for this translation in Fourier space
			transformer.FourierTransform(img(), Fimg);

			if (i_frame == first_frame)
			{
				Fwsum.initZeros(Fimg);
				Fsumw.initZeros(Fimg);
			}

			shiftImageInFourierTransform(Fimg, Fimg, ori_size, xtrans, ytrans);

			// Apply (positive!!) B-factor weighting and store weighted sums
			int iframe = i_frame - first_frame;
			RFLOAT bfactor = (do_weighting) ? DIRECT_A1D_ELEM(perframe_bfactors, 3*iframe + 0) : 0.;
			RFLOAT offset  = (do_weighting) ? DIRECT_A1D_ELEM(perframe_bfactors, 3*iframe + 1) : 0.;
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(Fimg)
			{
				RFLOAT res = sqrt((RFLOAT)ip * ip + jp * jp) / (XSIZE(img()) * angpix); // get resolution in 1/Angstrom
				if (res <= 1. / (angpix * 2.) ) // Apply B-factor weighting until Nyquist
				{
					RFLOAT w = exp( (bfactor / 4)  * res * res  + offset);
					DIRECT_A2D_ELEM(Fwsum, i, j) += w * DIRECT_A2D_ELEM(Fimg, i, j);
					DIRECT_A2D_ELEM(Fsumw, i, j) += w;
				}
			}
		}

		// Calculate the weighted average image
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fwsum)
		{
			if (DIRECT_MULTIDIM_ELEM(Fsumw, n) > 0.)
				DIRECT_MULTIDIM_ELEM(Fwsum, n) /= DIRECT_MULTIDIM_ELEM(Fsumw, n);
		}
		transformer.inverseFourierTransform(Fwsum, img());

		if (do_normalise)
			normalise(img, bg_radius, white_dust_stddev, black_dust_stddev, do_ramp);

		// Calculate statistics for the (normalised?) particle
		RFLOAT minval, maxval, avgval, stddev;
		img().computeStats(avgval, stddev, minval, maxval);

		// Keep track of overall statistics for all particles in this field-of-view
		all_minval = XMIPP_MIN(minval, all_minval);
		all_maxval = XMIPP_MAX(maxval, all_maxval);
		all_avg	+= avgval;
		all_stddev += stddev*stddev;

		// write the new average (i.e. the shiny, or polished particle)
		fn_part = fn_part.withoutExtension();
		std::string mic_name;
		long int nr;
		fn_part.decompose(nr, mic_name);
		fn_part = mic_name + "_" + fn_out + ".mrcs";
		fn_img.compose(ipar + 1, fn_part);

		// When last particle, also write the correct header
		if (ipar == exp_model.average_micrographs[imic].ori_particles_id.size() - 1)
		{
			all_avg /= exp_model.average_micrographs[imic].ori_particles_id.size();
			all_stddev = sqrt(all_stddev/exp_model.average_micrographs[imic].ori_particles_id.size());
			img.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, all_minval);
			img.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, all_maxval);
			img.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, all_avg);
			img.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, all_stddev);
		}

		if (ipar == 0)
			img.write(fn_img, -1, false, WRITE_OVERWRITE);
		else
			img.write(fn_img, -1, true, WRITE_APPEND);

	}

}

void ParticlePolisher::reconstructShinyParticlesAndFscWeight(int ipass)
{

	if (verb > 0)
		std::cout << "+ Reconstructing two halves of shiny particles ..." << std::endl;

	// Re-read the shiny particles' MetaDataTable into exp_model
	exp_model.read(fn_out + ".star", true);

	// Do the reconstructions for both halves
	reconstructShinyParticlesOneHalf(1);
	reconstructShinyParticlesOneHalf(2);


	FileName fn_post = (ipass == 1) ? "_post" : "_post2";
	if (!do_start_all_over && exists(fn_in.withoutExtension() + "_" + fn_out + fn_post + "_masked.mrc")
					       && exists(fn_in.withoutExtension() + "_" + fn_out + fn_post + ".star") )
	{
		if (verb > 0)
			std::cout << std::endl << " + " << fn_in.withoutExtension() << "_" << fn_out << fn_post << "_masked.mrc already exists: re-reading map into memory." << std::endl;

		if (verb > 0)
			std::cout << std::endl << " + " << fn_in.withoutExtension() << "_" << fn_out << fn_post << ".star already exists: re-reading resolution from it." << std::endl;

		MetaDataTable MD;
		MD.read(fn_in.withoutExtension() + "_" + fn_out + fn_post + ".star", "general");
		MD.getValue(EMDL_POSTPROCESS_FINAL_RESOLUTION, maxres_model);
	}
	else
	{
		// Re-read the two halves to calculate FSCs
		Postprocessing prm;

		prm.clear();
		prm.fn_in = fn_in.withoutExtension() + "_" + fn_out;
		prm.fn_out = prm.fn_in + fn_post;
		prm.angpix = angpix;
		prm.do_auto_mask = false;
		prm.fn_mask = fn_mask;
		prm.do_auto_bfac = false;
		prm.do_fsc_weighting = true;
		prm.verb = 0;
		prm.run();

		maxres_model = prm.global_resol;
	}

	if (verb > 0)
		std::cout << " + Setting Fourier transforms of the two shiny half-reconstructions ..." << std::endl;

	MultidimArray<RFLOAT> dum;
	Image<RFLOAT> refvol;
	FileName fn_vol;
	fn_vol = fn_in.withoutExtension() + "_" + fn_out + "_half1_class001_unfil.mrc";
	refvol.read(fn_vol);
	PPrefvol_half1.ori_size = XSIZE(refvol());
	PPrefvol_half1.padding_factor = 2;
	PPrefvol_half1.interpolator = TRILINEAR;
	PPrefvol_half1.r_min_nn = 10;
	PPrefvol_half1.data_dim = 2;
	PPrefvol_half1.computeFourierTransformMap(refvol(), dum);
	fn_vol = fn_in.withoutExtension() + "_" + fn_out + "_half2_class001_unfil.mrc";
	refvol.read(fn_vol);
	PPrefvol_half2.ori_size = XSIZE(refvol());
	PPrefvol_half2.padding_factor = 2;
	PPrefvol_half2.interpolator = TRILINEAR;
	PPrefvol_half2.r_min_nn = 10;
	PPrefvol_half2.data_dim = 2;
	PPrefvol_half2.computeFourierTransformMap(refvol(), dum);

}

void ParticlePolisher::reconstructShinyParticlesOneHalf(int this_half)
{

	FileName fn_vol = fn_in.withoutExtension() + "_" + fn_out + "_half" + integerToString(this_half) + "_class001_unfil.mrc";
	if (!do_start_all_over && exists(fn_vol))
	{
		if (verb > 0)
			std::cout << std::endl << " + " << fn_vol << " already exists: skipping shiny-particle reconstruction." << std::endl;
		return;
	}

	// get image size, angpix (from metadatatable), fn_sym
	int image_size;
	exp_model.MDexp.getValue(EMDL_IMAGE_SIZE, image_size);
	BackProjector backprojector(image_size, 3, fn_sym);
	backprojector.initZeros();
	Projector projector(image_size);


	CTF ctf;
	std::string dum;
	Matrix2D<RFLOAT> A3D;
	MultidimArray<Complex > Faux, F2D, Fsub;
	MultidimArray<RFLOAT> Fweight, Fctf;
	Image<RFLOAT> img, vol;
	FourierTransformer transformer;
	RFLOAT xtrans, ytrans;
	RFLOAT rot, tilt, psi;
	int i_half;
	long int i_frame;
	FileName fn_img;

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(exp_model.MDimg)
	{

		exp_model.MDimg.getValue(EMDL_PARTICLE_RANDOM_SUBSET, i_half);
		if (i_half == this_half)
		{
			exp_model.MDimg.getValue(EMDL_IMAGE_NAME, fn_img);
			img.read(fn_img);
			CenterFFT(img(), true);
			transformer.FourierTransform(img(), F2D);

			// Use the prior-angles, as these were determined from the average particles
			// The individual-frame-determined angles would be too noisy....
			exp_model.MDimg.getValue(EMDL_ORIENT_ROT, rot);
			exp_model.MDimg.getValue(EMDL_ORIENT_TILT, tilt);
			exp_model.MDimg.getValue(EMDL_ORIENT_PSI, psi);
			Euler_angles2matrix(rot, tilt, psi, A3D);

			// Translations
			xtrans = ytrans = 0.;
			exp_model.MDimg.getValue( EMDL_ORIENT_ORIGIN_X, xtrans);
			exp_model.MDimg.getValue( EMDL_ORIENT_ORIGIN_Y, ytrans);
			if (ABS(xtrans) > 0. || ABS(ytrans) > 0. )
				shiftImageInFourierTransform(F2D, F2D, image_size, xtrans, ytrans );

			// CTF
			Fctf.resize(F2D);
			Fctf.initConstant(1.);
			if (do_ctf)
			{
				ctf.read(exp_model.MDimg, exp_model.MDimg);
				ctf.getFftwImage(Fctf, image_size, image_size, angpix, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
				{
					DIRECT_MULTIDIM_ELEM(F2D, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
				}
			}

			backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
		}

	}


	backprojector.symmetrise(nr_helical_asu, helical_twist, helical_rise / angpix);

	// Now do the reconstruction
	MultidimArray<RFLOAT> dummy;
	backprojector.reconstruct(vol(), 10, false, 1., dummy, dummy, dummy, dummy);

	vol.write(fn_vol);

}



void ParticlePolisher::optimiseBeamTiltAndDefocus()
{

	// This function assumes the shiny particles are in exp_mdel.MDimg!!
	if (beamtilt_max <= 0. && defocus_shift_max <= 0.)
		return;

	if (minres_beamtilt < maxres_model)
	{
		if (verb > 0)
			std::cout << " Skipping beamtilt correction, as the resolution of the shiny reconstruction  does not go beyond minres_beamtilt of " << minres_beamtilt << " Ang." << std::endl;
		return;
	}

	getBeamTiltGroups();

	initialiseSquaredDifferenceVectors();

	// Loop over all average micrographs
	int barstep;
	int my_nr_micrographs = exp_model.micrographs.size();
	if (verb > 0)
	{
		std::cout << " + Optimising beamtilts in all micrographs ... " << std::endl;
		init_progress_bar(my_nr_micrographs);
		barstep = XMIPP_MAX(1, my_nr_micrographs/ 60);
	}

    for (long int i = 0; i < my_nr_micrographs; i++)
	{
    	if (verb > 0 && i % barstep == 0)
			progress_bar(i);

    	optimiseBeamTiltAndDefocusOneMicrograph(i);
	}

    if (verb > 0)
	{
		progress_bar(my_nr_micrographs);
	}

    // Now get the final optimised beamtilts!
    applyOptimisedBeamTiltsAndDefocus();

    // Write the new MDTable to disc
	if (verb > 0)
		exp_model.MDimg.write(fn_out + ".star");

}

void ParticlePolisher::getBeamTiltGroups()
{
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(exp_model.MDimg)
	{

		// Get the name of the beamtilt group (micrograph name if no groups are defined...)
		FileName fn_group;
		if (exp_model.MDimg.containsLabel(EMDL_IMAGE_BEAMTILT_GROUP))
			exp_model.MDimg.getValue(EMDL_IMAGE_BEAMTILT_GROUP, fn_group);
		else
			exp_model.MDimg.getValue(EMDL_MICROGRAPH_NAME, fn_group);

		bool is_unique = true;
		for (int igroup = 0; igroup < fn_beamtilt_groups.size(); igroup++)
		{
			if (fn_beamtilt_groups[igroup] == fn_group)
			{
				is_unique = false;
				break;
			}
		}
		if (is_unique)
		{
			fn_beamtilt_groups.push_back(fn_group);
		}
	}

}

void ParticlePolisher::initialiseSquaredDifferenceVectors()
{

	nr_sampled_beam_tilts= 0;
	if (beamtilt_max > 0.)
	{
		int n_steps = CEIL(beamtilt_max / beamtilt_step);
		for (RFLOAT tilt_y = -n_steps*beamtilt_step; tilt_y <= n_steps*beamtilt_step; tilt_y += beamtilt_step)
		{
			for (RFLOAT tilt_x = -n_steps*beamtilt_step; tilt_x <= n_steps*beamtilt_step; tilt_x += beamtilt_step)
			{
				if (sqrt(tilt_y*tilt_y + tilt_x*tilt_x) <= beamtilt_max)
				{
					nr_sampled_beam_tilts++;
				}
			}
		}
		// Store squared differences for all beam tilts and for all data sets
		diff2_beamtilt.initZeros(fn_beamtilt_groups.size(), nr_sampled_beam_tilts);
	}

	if (defocus_shift_max > 0.)
	{
		defocus_shift_allmics.initZeros(exp_model.micrographs.size());
	}


}

void ParticlePolisher::applyOptimisedBeamTiltsAndDefocus()
{
	if (beamtilt_max > 0.)
	{
		// Use in two different loops
		int n_steps = CEIL(beamtilt_max / beamtilt_step);

		best_beamtilts.clear();
		Matrix1D<RFLOAT> my_tilts(2);
		for (int igroup = 0; igroup < fn_beamtilt_groups.size(); igroup++)
		{
			RFLOAT mindiff2 = LARGE_NUMBER;
			RFLOAT best_tilt_x, best_tilt_y;
			if (verb > 0)
				std::cout << " + Beamtilt group " << fn_beamtilt_groups[igroup] << std::endl;
			int n_tilt= 0;
			for (RFLOAT tilt_y = -n_steps*beamtilt_step; tilt_y <= n_steps*beamtilt_step; tilt_y += beamtilt_step)
			{
				for (RFLOAT tilt_x = -n_steps*beamtilt_step; tilt_x <= n_steps*beamtilt_step; tilt_x += beamtilt_step)
				{
					if (sqrt(tilt_y*tilt_y + tilt_x*tilt_x) <= beamtilt_max)
					{
						RFLOAT diff2 = DIRECT_A2D_ELEM(diff2_beamtilt, igroup, n_tilt);
						if (verb > 1)
							std::cout << " + tilt_x = " << tilt_x << " + tilt_y = " << tilt_y << " diff2= " << diff2 << std::endl;
#ifdef DEBUG_TILT

						if (verb > 0)
							std::cerr << " igroup= " << igroup << " n_tilt= " << n_tilt
							<< " DIRECT_A2D_ELEM(diff2_beamtilt, igroup, n_tilt)= " << DIRECT_A2D_ELEM(diff2_beamtilt, igroup, n_tilt)
							<< " diff2= " << diff2 << " mindiff2= " << mindiff2
							<< std::endl;
#endif
						if (diff2 <= mindiff2)
						{
							best_tilt_x = tilt_x;
							best_tilt_y = tilt_y;
							mindiff2 = diff2;
						}
						n_tilt ++;
					}
				}
			}
			if (verb > 0)
				std::cout << " + Best tilt_x = " << best_tilt_x << " best tilt_y = " << best_tilt_y << " mindiff2= " << mindiff2 << std::endl;
			XX(my_tilts) = best_tilt_x;
			YY(my_tilts) = best_tilt_y;
			best_beamtilts.push_back(my_tilts);
		}
	}


	// Now set beamtilts in the MetaDataTable
    int i_group;
	for (long int imic = 0; imic < exp_model.micrographs.size(); imic++)
	{

		RFLOAT best_defocus_shift = (defocus_shift_max > 0.) ? DIRECT_A1D_ELEM(defocus_shift_allmics, imic) : 0.;
		for (long int ipart = 0; ipart < exp_model.micrographs[imic].particle_ids.size(); ipart++)
    	{
			long int part_id = exp_model.micrographs[imic].particle_ids[ipart];

    		// Set the optimised beamtilts in the MetadataTable
    		if (beamtilt_max > 0.)
    		{
				// First get which beamtilt group this micrograph comes from
				if (ipart == 0)
				{
					FileName fn_group;
					if (exp_model.MDimg.containsLabel(EMDL_IMAGE_BEAMTILT_GROUP))
						exp_model.MDimg.getValue(EMDL_IMAGE_BEAMTILT_GROUP, fn_group, part_id);
					else
						exp_model.MDimg.getValue(EMDL_MICROGRAPH_NAME, fn_group, part_id);
					bool found = false;
					for (int igroup = 0; igroup < fn_beamtilt_groups.size(); igroup++)
					{
						if (fn_group == fn_beamtilt_groups[igroup])
						{
							i_group = igroup;
							found = true;
							break;
						}
					}
					if (!found)
						REPORT_ERROR("ParticlePolisher::optimiseBeamTiltOneMicrograph ERROR: could not find data set name for " + fn_group);

				}

				// Then set the corresponding beamtilt values for each particles
				exp_model.MDimg.setValue(EMDL_IMAGE_BEAMTILT_X, XX(best_beamtilts[i_group]), part_id);
				exp_model.MDimg.setValue(EMDL_IMAGE_BEAMTILT_Y, YY(best_beamtilts[i_group]), part_id);

    		}

			// Also set the optimised defocus values for each micrograph
			if (defocus_shift_max > 0.)
			{

				if (ipart == 0 && verb > 0)
					std::cout << " + Micrograph " << exp_model.micrographs[imic].name << " defocus_shift= " << best_defocus_shift << std::endl;

				RFLOAT ori_defocusU, ori_defocusV;
				exp_model.MDimg.getValue(EMDL_CTF_DEFOCUSU, ori_defocusU, part_id);
				exp_model.MDimg.getValue(EMDL_CTF_DEFOCUSV, ori_defocusV, part_id);
				exp_model.MDimg.setValue(EMDL_CTF_DEFOCUSU, ori_defocusU + best_defocus_shift, part_id);
				exp_model.MDimg.setValue(EMDL_CTF_DEFOCUSV, ori_defocusV + best_defocus_shift, part_id);
			}
    	}
	}



}

void ParticlePolisher::optimiseBeamTiltAndDefocusOneMicrograph(int imic)
{

	// get image size, angpix (from metadatatable), fn_sym
	int image_size;
	exp_model.MDexp.getValue(EMDL_IMAGE_SIZE, image_size);

	CTF ctf;
	Matrix2D<RFLOAT> A3D;
	MultidimArray<Complex > F2D, F2Dtilt, Fref;
	MultidimArray<RFLOAT> Fctf;
	Image<RFLOAT> img;
	FourierTransformer transformer;
	RFLOAT xtrans, ytrans;
	RFLOAT rot, tilt, psi;
	FileName fn_img;
	int i_group = -1, my_half = 0;

	RFLOAT xsize = angpix * PPrefvol_half1.ori_size;

	// Weighted squared-differences for all defocusses
	int nr_sampled_defocus_shifts;
	MultidimArray<RFLOAT> wdiff2_defocus;
	if (defocus_shift_max > 0.)
	{
		nr_sampled_defocus_shifts = CEIL(defocus_shift_max / defocus_shift_step);
		wdiff2_defocus.initZeros(2*nr_sampled_defocus_shifts + 1);
	}

	for (long int ipart = 0; ipart < exp_model.micrographs[imic].particle_ids.size(); ipart++)
	{
		long int part_id = exp_model.micrographs[imic].particle_ids[ipart];

		// Get which beamtilt group this micrograph comes from
		if (ipart == 0)
		{
			FileName fn_group;
			if (exp_model.MDimg.containsLabel(EMDL_IMAGE_BEAMTILT_GROUP))
				exp_model.MDimg.getValue(EMDL_IMAGE_BEAMTILT_GROUP, fn_group, part_id);
			else
				exp_model.MDimg.getValue(EMDL_MICROGRAPH_NAME, fn_group, part_id);
			bool found = false;
			for (int igroup = 0; igroup < fn_beamtilt_groups.size(); igroup++)
			{
				if (fn_group == fn_beamtilt_groups[igroup])
				{
					i_group = igroup;
					found = true;
					break;
				}
			}
			if (!found)
				REPORT_ERROR("ParticlePolisher::optimiseBeamTiltOneMicrograph ERROR: could not find beamtilt group name for " + fn_group);

		}


		exp_model.MDimg.getValue(EMDL_IMAGE_NAME, fn_img, part_id);
		img.read(fn_img);
		CenterFFT(img(), true);
		transformer.FourierTransform(img(), F2D);

		// Which half do I belong to?
		exp_model.MDimg.getValue(EMDL_PARTICLE_RANDOM_SUBSET, my_half, part_id);

		// Use the prior-angles, as these were determined from the average particles
		// The individual-frame-determined angles would be too noisy....
		exp_model.MDimg.getValue(EMDL_ORIENT_ROT, rot, part_id);
		exp_model.MDimg.getValue(EMDL_ORIENT_TILT, tilt, part_id);
		exp_model.MDimg.getValue(EMDL_ORIENT_PSI, psi, part_id);
		Euler_angles2matrix(rot, tilt, psi, A3D);

		// Get the reference projection
		Fref.resize(F2D);
		if (my_half == 1)
			PPrefvol_half1.get2DFourierTransform(Fref, A3D, IS_NOT_INV);
		else if (my_half == 2)
			PPrefvol_half2.get2DFourierTransform(Fref, A3D, IS_NOT_INV);
		else
			REPORT_ERROR("ERROR unrecognised random subset, not 1 or 2...");

		// Shift the experimental image
		xtrans = ytrans = 0.;
		exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, xtrans, part_id);
		exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, ytrans, part_id);
		if (ABS(xtrans) > 0. || ABS(ytrans) > 0. )
			shiftImageInFourierTransform(F2D, F2D, image_size, xtrans, ytrans );

		// apply CTF to the reference
		if (do_ctf)
		{

			Fctf.resize(F2D);
			ctf.read(exp_model.MDimg, exp_model.MDimg, part_id);

			// Store the original values of defocusU and defocusV
			RFLOAT ori_defocusU = ctf.DeltafU;
			RFLOAT ori_defocusV = ctf.DeltafV;
			RFLOAT best_defocus_shift;
			RFLOAT mindiff_thispart = LARGE_NUMBER;
			// Optimise per-particle CTF defocus
			// For now only non-anisotropically
			if (defocus_shift_max > 0.)
			{

				int n_defocus = 0;
				for (RFLOAT defocus_shift = -nr_sampled_defocus_shifts*defocus_shift_step;
						defocus_shift <= nr_sampled_defocus_shifts*defocus_shift_step;
						defocus_shift += defocus_shift_step)
				{
					// Get modified CTF and apply to reference
					ctf.DeltafU = ori_defocusU + defocus_shift;
					ctf.DeltafV = ori_defocusV + defocus_shift;
					ctf.initialise();
					ctf.getFftwImage(Fctf, image_size, image_size, angpix, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);

					// Calculate squared difference with the image
					RFLOAT diff2 = 0.;
					FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(F2D)
					{
						RFLOAT res = xsize/sqrt((RFLOAT)(ip * ip + jp * jp)); // get resolution in 1/pixel
						if (res <= minres_beamtilt && res >= maxres_model)
				    	{
							RFLOAT diff_real = DIRECT_A2D_ELEM(Fctf, i, j) * (DIRECT_A2D_ELEM(Fref, i, j)).real - (DIRECT_A2D_ELEM(F2D, i, j)).real;
							RFLOAT diff_imag = DIRECT_A2D_ELEM(Fctf, i, j) * (DIRECT_A2D_ELEM(Fref, i, j)).imag - (DIRECT_A2D_ELEM(F2D, i, j)).imag;
							diff2 += (diff_real * diff_real + diff_imag * diff_imag);
				    	}
					}
					// Store the accumulated squared differences...
					if (diff2 < mindiff_thispart)
					{
						mindiff_thispart = diff2;
						best_defocus_shift = defocus_shift;
					}
					//std::cerr << " iimg= " << iimg << " defocus_shift= " << defocus_shift << " diff2= " << diff2 << std::endl;

					DIRECT_A1D_ELEM(wdiff2_defocus, n_defocus) += diff2;
					n_defocus++;
				}

				// Set best defocus for each individual particle...
				// TODO!!! This only works in sequential version!!!
				//exp_model.MDimg.setValue(EMDL_CTF_DEFOCUSU, ori_defocusU + best_defocus_shift, part_id);
				//exp_model.MDimg.setValue(EMDL_CTF_DEFOCUSV, ori_defocusV + best_defocus_shift, part_id);

				// Re-set the original defocus values in the ctf object
				ctf.DeltafU = ori_defocusU;
				ctf.DeltafV = ori_defocusV;
				ctf.initialise();
			}


		}

		if (beamtilt_max > 0.)
		{

			if (do_ctf)
			{
				// After per-particle assessment, just re-read the stored CTF and apply to the reference projection
				ctf.getFftwImage(Fctf, image_size, image_size, angpix, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
				{
					DIRECT_MULTIDIM_ELEM(Fref, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
				}
			}

			F2Dtilt.resize(F2D);
			// Loop over all beam tilts
			int n_tilts = 0;
			int n_steps = CEIL(beamtilt_max / beamtilt_step);
			for (RFLOAT tilt_y = -n_steps*beamtilt_step; tilt_y <= n_steps*beamtilt_step; tilt_y += beamtilt_step)
			{
				for (RFLOAT tilt_x = -n_steps*beamtilt_step; tilt_x <= n_steps*beamtilt_step; tilt_x += beamtilt_step)
				{
					if (sqrt(tilt_y*tilt_y + tilt_x*tilt_x) <= beamtilt_max)
					{
						// Now calculate the squared differences, taking the beamtilt into account
						applyBeamTilt(F2D, F2Dtilt, tilt_x, tilt_y, ctf.lambda, ctf.Cs, angpix, image_size);

						RFLOAT diff2 = 0.;
						RFLOAT ndiff = 0.;
						FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(F2Dtilt)
						{
							// Store amplitude-weighted phase difference and sum of amplitude-weights
							RFLOAT res = xsize/sqrt((RFLOAT)(ip * ip + jp * jp)); // get resolution in 1/pixel
							if (res <= minres_beamtilt && res >= maxres_model)
							{

								RFLOAT diff_real = (DIRECT_A2D_ELEM(Fref, i, j)).real - (DIRECT_A2D_ELEM(F2Dtilt, i, j)).real;
								RFLOAT diff_imag = (DIRECT_A2D_ELEM(Fref, i, j)).imag - (DIRECT_A2D_ELEM(F2Dtilt, i, j)).imag;
								diff2 += (diff_real * diff_real + diff_imag * diff_imag);
								ndiff += 1.;
							}
						}
						diff2 /= ndiff;

						// Store the accumulate weighted differences...
						DIRECT_A2D_ELEM(diff2_beamtilt, i_group, n_tilts) += diff2;

						n_tilts++;
						if (n_tilts > nr_sampled_beam_tilts)
						{
							std::cerr << " n_tilts= " << n_tilts << " nr_sampled_beam_tilts= " << nr_sampled_beam_tilts << " beamtilt_max= " << beamtilt_max << std::endl;
							std::cerr << " tilt_x= " << tilt_x << " tilt_y= " << tilt_y << " beamtilt_step= " << beamtilt_step << std::endl;
							REPORT_ERROR("BUG: too large n_tilts....");
						}

					} // end if sqrt(tilt_y*tilt_y + tilt_x*tilt_x) <= beamtilt_max
				} // end for tilt_x
			} // end for tilt_y
		} // end if beamtilt_max > 0

	} // end for over all particles in this micrograph



	// Set optimal defocus shift averaged over all particles in this micrograph:
	if (defocus_shift_max > 0.)
	{
		RFLOAT mindiff2=LARGE_NUMBER, best_defocus_shift;
		int n_defocus = 0;
		for (RFLOAT defocus_shift = -nr_sampled_defocus_shifts*defocus_shift_step;
				defocus_shift <= nr_sampled_defocus_shifts*defocus_shift_step;
				defocus_shift += defocus_shift_step)
		{
			RFLOAT diff2 = DIRECT_A1D_ELEM(wdiff2_defocus, n_defocus);
			//std::cerr << std::setprecision(10) << " imic= " << imic << " defocus_shift= " << defocus_shift << " diff2= " << diff2 << std::endl;
			if (diff2 < mindiff2)
			{
				mindiff2 = diff2;
				best_defocus_shift = defocus_shift;
			}
			n_defocus++;
		}
		DIRECT_A1D_ELEM(defocus_shift_allmics, imic) = best_defocus_shift;

		//exp_model.MDimg.setValue(EMDL_CTF_DEFOCUSU, ori_defocusU + best_defocus_shift, part_id);
		//exp_model.MDimg.setValue(EMDL_CTF_DEFOCUSV, ori_defocusV + best_defocus_shift, part_id);
		// Set best defocus for all particles on this micrograph
		// TODO: this only works in sequential version!!!
		for (long int ipart = 0; ipart < exp_model.micrographs[imic].particle_ids.size(); ipart++)
		{
			long int part_id = exp_model.micrographs[imic].particle_ids[ipart];
			RFLOAT ori_defocusU, ori_defocusV;
			exp_model.MDimg.getValue(EMDL_CTF_DEFOCUSU, ori_defocusU, part_id);
			exp_model.MDimg.getValue(EMDL_CTF_DEFOCUSV, ori_defocusV, part_id);


			exp_model.MDimg.setValue(EMDL_CTF_DEFOCUSU, ori_defocusU + best_defocus_shift, part_id);
			exp_model.MDimg.setValue(EMDL_CTF_DEFOCUSV, ori_defocusV + best_defocus_shift, part_id);
		}

	}
}

void ParticlePolisher::run()
{

	// Fit straight lines through all beam-induced translations
	if (fitting_mode != NO_FIT)
		fitMovementsAllMicrographs();

	// Perform single-frame reconstructions and estimate dose-dependent B-factors
        if (do_weighting)
	     calculateAllSingleFrameReconstructionsAndBfactors();

	// Write out the intermediately polished particles
	polishParticlesAllMicrographs();

	// Now reconstruct with all polished particles: two independent halves, FSC-weighting of the sum of the two...
	reconstructShinyParticlesAndFscWeight(1);

	// Optimise beam-tilt and defocus per beamtilt group and/or micrograph
	optimiseBeamTiltAndDefocus();

	// Reconstruct again two halves to see whether the beamtilt and/or defocus optimisation has helped
	if (beamtilt_max > 0. || defocus_shift_max > 0.)
		reconstructShinyParticlesAndFscWeight(2);

	if (verb > 0)
		std::cout << " done!" << std::endl;

}


