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
	fn_out = parser.getOption("--o", "Output directory", "Polish/");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "-1"));
	only_do_unfinished = parser.checkOption("--only_do_unfinished", "Skip those steps for which output files already exist.");

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
	is_helix = parser.checkOption("--helix", "Polish helical segments?");
	helical_nr_asu = textToInteger(parser.getOption("--helical_nr_asu", "Number of new helical asymmetric units (asu) per box (1 means no helical symmetry is present)", "1"));
	helical_twist = textToFloat(parser.getOption("--helical_twist", "Helical twist (in degrees, positive values for right-handedness)", "0."));
	helical_rise = textToFloat(parser.getOption("--helical_rise", "Helical rise (in Angstroms)", "0."));
	fn_mask = parser.getOption("--mask", "Postprocessing mask for B-factor determination of per-frame reconstructions (1=protein, 0=solvent, all values in range [0,1])", "");

	int ctf_section = parser.addSection("CTF options");
   	do_ctf = !parser.checkOption("--no_ctf", "Don't apply CTF correction");
	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");
	ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Images have been phase flipped");
	only_flip_phases = parser.checkOption("--only_flip_phases", "Do not correct CTF-amplitudes, only flip phases");
	//defocus_shift_max = textToFloat(parser.getOption("--defocus_shift_max", "Maximum shift in defocus (in A) to search for each particle", "0"));
	//defocus_shift_step = textToFloat(parser.getOption("--defocus_shift_step", "Maximum shift in defocus (in A) to search for each particle", "100"));

	int norm_section = parser.addSection("Normalisation options");
	do_normalise = !parser.checkOption("--skip_normalise", "Do not normalise the polsihed particles?");
	bg_radius =  textToInteger(parser.getOption("--bg_radius", "Radius of the circular mask that will be used to define the background area (in pixels)", "-1"));
	do_ramp = !parser.checkOption("--no_ramp", "Just subtract the background mean in the normalisation, instead of subtracting a fitted ramping background. ");
	white_dust_stddev = textToFloat(parser.getOption("--white_dust", "Sigma-values above which white dust will be removed (negative value means no dust removal)","-1"));
	black_dust_stddev = textToFloat(parser.getOption("--black_dust", "Sigma-values above which black dust will be removed (negative value means no dust removal)","-1"));

	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	// Make sure fn_out ends with a slash
	if (fn_out[fn_out.length()-1] != '/')
		fn_out += "/";

}

void ParticlePolisher::usage()
{
	parser.writeUsage(std::cout);
}

void ParticlePolisher::generateMicrographList()
{

    if (only_do_unfinished && exists(fn_out+"micrograph_list.star"))
    {
    	std::cout << " + Reading in pre-existing micrograph_list.star file "<<std::endl;

    	MetaDataTable MDmics;
		MDmics.read(fn_out+"micrograph_list.star");
		fn_mics.clear();
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDmics)
		{
			FileName fn_mic;
			MDmics.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
			fn_mics.push_back(fn_mic);
		}
    }
    else
    {

		if (verb > 0)
			std::cout << " + Reading the input STAR file ... " << std::endl;

		MetaDataTable MDin;
		MDin.read(fn_in);

		// list of input STAR files
		std::vector<FileName> fn_stars;
		bool input_is_movie_data = false;
		if (MDin.containsLabel(EMDL_STARFILE_MOVIE_PARTICLES))
		{
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
			{
				FileName fnt;
				MDin.getValue(EMDL_STARFILE_MOVIE_PARTICLES, fnt);
				fn_stars.push_back(fnt);
			}
		}
		else
		{
			// The input is already a movie-data file
			fn_stars.push_back(fn_in);
			input_is_movie_data = true;
		}

		// Now break up all STAR files into pieces, one for each micrograph and store in fn_mics
		for (int istar = 0; istar < fn_stars.size(); istar++)
		{

			// No need to re-read MDin if the input was already movie data!
			if (!input_is_movie_data)
				MDin.read(fn_stars[istar]);

			MDin.newSort(EMDL_MICROGRAPH_NAME, false, true); // false=no reverse, true= do sort only on string after "@"

			FileName fn_old="", fn_curr, fn_pre, fn_jobnr, fn_post;
			MetaDataTable MDonemic;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
			{
				MDin.getValue(EMDL_MICROGRAPH_NAME, fn_curr);
				fn_curr=fn_curr.substr(fn_curr.find("@")+1);
				decomposePipelineFileName(fn_curr, fn_pre, fn_jobnr, fn_curr);

				if (fn_curr != fn_old && fn_old != "")
				{
					FileName fn_star = fn_out + fn_old.withoutExtension()+"_input.star";
					fn_mics.push_back(fn_star);
					FileName fn_dir = fn_star.beforeLastOf("/");
					if (!exists(fn_dir))
						int res = system(("mkdir -p " + fn_dir).c_str());
					MDonemic.write(fn_star);
					MDonemic.clear();
				}

				MDonemic.addObject(MDin.getObject());

				// Reset the old micrograph name
				fn_old=fn_curr;

			} // end loop all objects in input STAR file

			// Also write the last MDonemic
			FileName fn_star = fn_out + fn_old.withoutExtension()+"_input.star";
			fn_mics.push_back(fn_star);
			FileName fn_dir = fn_star.beforeLastOf("/");
			if (!exists(fn_dir))
				int res = system(("mkdir -p " + fn_dir).c_str());
			MDonemic.write(fn_star);
			MDonemic.clear();

		} // end loop over all fn_stars

		// Write fn_mics into a STAR file to disk (this is used by MPI version to read back in for all slaves, and may be generally useful for tracing errors)
		MetaDataTable MDmics;
		for (long int imic = 0; imic < fn_mics.size(); imic++)
		{
			MDmics.addObject();
			MDmics.setValue(EMDL_MICROGRAPH_NAME, fn_mics[imic]);
		}
		MDmics.write(fn_out+"micrograph_list.star");
    }
}


void ParticlePolisher::initialise()
{

	// Read in first micrograph STAR file to get some general informatopm
    Experiment exp_model;
	exp_model.read(fn_mics[0]);

	// Get the vector with which movie-frame numbers are stored here
	movie_frame_numbers.clear();
	for (int ipart = 0; ipart < exp_model.ori_particles[0].particles_id.size(); ipart++)
	{
		FileName fn_mic, dum;
		long int i_frame;
		long int part_id = exp_model.ori_particles[0].particles_id[ipart];
		exp_model.MDimg.getValue(EMDL_MICROGRAPH_NAME, fn_mic, part_id);
		fn_mic.decompose(i_frame, dum);
		movie_frame_numbers.push_back(i_frame);
		// Also get the image size
		if (ipart == 0)
		{
			FileName fn_img;
			exp_model.MDimg.getValue(EMDL_IMAGE_NAME, fn_img, part_id);
			Image<RFLOAT> It;
			It.read(fn_img, false); // false means only read header
			ori_size = XSIZE(It());
		}
	}

	// Initialise array for storing bfactor, offset and corr_coeff
	perframe_bfactors.initZeros( 3 * movie_frame_numbers.size());

	// Get the pixel size from the input STAR file
	if (angpix < 0.)
	{
		if (exp_model.MDimg.containsLabel(EMDL_CTF_MAGNIFICATION) && exp_model.MDimg.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
		{
			RFLOAT mag, dstep;
			exp_model.MDimg.goToObject(0);
			exp_model.MDimg.getValue(EMDL_CTF_MAGNIFICATION, mag);
			exp_model.MDimg.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
			angpix = 10000. * dstep / mag;
			if (verb > 0)
				std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << std::endl;
		}
		else
		{
			REPORT_ERROR("ParticlePolisher::initialise ERROR: The input STAR files do not contain information about the pixel size. Please provide --angpix");
		}
	}

	if (do_weighting)
    {

    	if (fn_mask == "")
            REPORT_ERROR("When doing B-factor weighting, you have to provide a mask!");

    	// Read in the Imask
    	Imask.read(fn_mask);
    	if (ori_size != XSIZE(Imask()))
    	{
    		std::cerr << " image box size= " << ori_size << " mask box size= " << XSIZE(Imask()) << std::endl;
    		REPORT_ERROR("ParticlePolisher:: ERROR: size of the mask is not equal to the size of the movie-particles.");
    	}
    }

	if (do_normalise && bg_radius < 0)
	{
		bg_radius = ROUND(0.375 * ori_size);

		if (verb > 0)
			std::cout << " + Using default particle diameter of 75% of the box size, so bg_radius= " << bg_radius << std::endl;
	}

	// Sep24,2015 - Shaoda, Helical reconstruction
	if ( (is_helix) && (helical_nr_asu > 1) )
	{
		if ( (fabs(helical_twist) > 360.) || (angpix < 0.001) || ((helical_rise / angpix) < 0.001))
			REPORT_ERROR("ParticlePolisher::initialise ERROR: Invalid helical twist or rise!");
	}

	fn_olddir = "";
}

// Fit the beam-induced translations for all average micrographs
void ParticlePolisher::fitMovementsAllMicrographs()
{

	// Loop over all average micrographs
	int barstep;
	long int my_nr_micrographs = fn_mics.size();
	if (verb > 0)
	{
		std::cout << " + Fitting straight paths for beam-induced movements in all " << my_nr_micrographs << " micrographs ... " << std::endl;
		init_progress_bar(my_nr_micrographs);
		barstep = XMIPP_MAX(1, my_nr_micrographs/ 60);
	}

    for (long int i = 0; i < my_nr_micrographs; i++)
	{
    	if (verb > 0 && i % barstep == 0)
			progress_bar(i);

		fitMovementsOneMicrograph(i);
	}

    if (verb > 0)
	{
		progress_bar(my_nr_micrographs);
	}


}



void ParticlePolisher::fitMovementsOneMicrograph(long int imic)
{

	FileName fn_fit = fn_mics[imic].withoutExtension() + "_fit.star";
	if (only_do_unfinished && exists(fn_fit))
		return;



	// Also write out a postscript file with the fits
	FileName fn_eps = fn_mics[imic].withoutExtension() + "_fit.eps";
	CPlot2D *plot2D=new CPlot2D(fn_mics[imic]);
	plot2D->SetXAxisSize(600);
	plot2D->SetYAxisSize(600);
	plot2D->SetDrawLegend(false);

	std::vector<RFLOAT> x_pick, y_pick, x_off_prior, y_off_prior, x_start, x_end, dummy; // X and Y-coordinates for the average particles in the micrograph
	std::vector< std::vector<RFLOAT> > x_off, y_off; // X and Y shifts w.r.t. the average for each movie frame
	RFLOAT x_pick_p, y_pick_p, x_off_p, y_off_p, x_off_prior_p, y_off_prior_p;

	Experiment exp_model;
	exp_model.read(fn_mics[imic]);

	// Running average window with used for the determination of the frame movements, now taken from _data.star
	int running_average_width;
	if (!exp_model.MDimg.containsLabel(EMDL_PARTICLE_MOVIE_RUNNING_AVG))
		REPORT_ERROR("ParticlePolisher::fitMovementsOneMicrograph ERROR: input STAR file does not contain rlnMovieFramesRunningAverage label");

	// Just testing we've read a single micrograph only!
	if (exp_model.micrographs.size()>1)
		REPORT_ERROR("BUG: exp_model.micrographs.size()= " + integerToString(exp_model.micrographs.size()));

	// Loop over all original_particles in this average_micrograph
	// And fill the x_pick, y_;pick, x_off and y_off vectors
	for (long int ipar = 0; ipar < (exp_model.micrographs[0]).ori_particle_ids.size(); ipar++)
	{
		long int ori_part_id = exp_model.micrographs[0].ori_particle_ids[ipar];

		CDataSet dataSet;
		dataSet.SetDrawMarker(false);
		dataSet.SetDatasetColor(0.0,0.0,0.0);
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
			exp_model.MDimg.getValue(EMDL_PARTICLE_MOVIE_RUNNING_AVG, running_average_width, part_id);
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
			CDataPoint point(x_pick_p + 25*x_off_p, y_pick_p + 25*y_off_p);
			dataSet.AddDataPoint(point);
		}
		plot2D->AddDataSet(dataSet);
	}

	// Now do the actual fitting
	RFLOAT gauss_const = 1. / sqrt(2 * PI * sigma_neighbour_distance * sigma_neighbour_distance);
	RFLOAT min2sigma2 = - 2. * sigma_neighbour_distance * sigma_neighbour_distance;
	// Loop over all ori_particles
	for (long int ipar = 0; ipar < (exp_model.micrographs[0]).ori_particle_ids.size(); ipar++)
	{
		long int ori_part_id = exp_model.micrographs[0].ori_particle_ids[ipar];

		CDataSet dataSet;
		dataSet.SetDrawMarker(false);
		dataSet.SetDatasetColor(1.0,0.0,0.0);

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
			exp_model.MDimg.setValue(EMDL_ORIENT_ORIGIN_X, x_off_p + x_off_prior[ipar], part_id);
			exp_model.MDimg.setValue(EMDL_ORIENT_ORIGIN_Y, y_off_p + y_off_prior[ipar], part_id);

			if (i_frame == 0)
			{
				exp_model.MDimg.getValue(EMDL_IMAGE_COORD_X, x_pick_p, part_id);
				exp_model.MDimg.getValue(EMDL_IMAGE_COORD_Y, y_pick_p, part_id);
			}
			CDataPoint point(x_pick_p + 25*(x_off_p + x_off_prior[ipar]), y_pick_p + 25*(y_off_p + y_off_prior[ipar]));
			dataSet.AddDataPoint(point);
		}
		plot2D->AddDataSet(dataSet);
	}

	// Write the STAR file with the fitted coordinates, make directory if it doesn't exist
	FileName fn_dir = fn_fit.beforeLastOf("/");
	if (!exists(fn_dir))
		int res = system(("mkdir -p " + fn_dir).c_str());

	exp_model.MDimg.write(fn_fit);

	// Write the movement plot as well
	plot2D->SetXAxisTitle("X-coordinate");
	plot2D->SetYAxisTitle("Y-coordinate");
	plot2D->OutputPostScriptPlot(fn_eps);
}

void ParticlePolisher::calculateAllSingleFrameReconstructionsAndBfactors()
{

	FileName fn_star = fn_out + "bfactors.star";
	if (only_do_unfinished && readStarFileBfactors(fn_star))
	{
		if (verb > 0)
			std::cout << " + " << fn_star << " already exists: skipping calculation of per-frame B-factors." <<std::endl;
		return;
	}

	RFLOAT bfactor, offset, corr_coeff;

	// Loop over all frames to be included in the reconstruction
	int my_nr_frames = movie_frame_numbers.size();
	int barstep;
	if (verb > 0)
	{
		std::cout << " + Calculating per-frame reconstructions ... " << std::endl;
		init_progress_bar(my_nr_frames);
		barstep = XMIPP_MAX(1, my_nr_frames/ 60);
	}

	for (long int iframe = 0; iframe < movie_frame_numbers.size(); iframe++)
	{

		calculateSingleFrameReconstruction(iframe, 1);
		calculateSingleFrameReconstruction(iframe, 2);

		if (verb > 0 && iframe % barstep == 0)
			progress_bar(iframe);

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
    for (long int iframe = 0; iframe < my_nr_frames; iframe++)
	{

		calculateBfactorSingleFrameReconstruction(iframe, bfactor, offset, corr_coeff);
		DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 0) = bfactor;
		DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 1) = offset;
		DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 2) = corr_coeff;

		if (verb > 0 && iframe % barstep == 0)
			progress_bar(iframe);
	}


    if (verb > 0)
	{
		progress_bar(my_nr_frames);
	}

    // Write STAR file
    writeStarFileBfactors(fn_star);


    // Also write a STAR file with the relative contributions of each frame to all frequencies
    fn_star = fn_out + "relweights.star";
    writeStarFileRelativeWeights(fn_star);


}

bool ParticlePolisher::readStarFileBfactors(FileName fn_star)
{
	if (exists(fn_star))
	{
		MetaDataTable MD;
		MD.read(fn_star);
		int iframe = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			int itest;
			MD.getValue(EMDL_IMAGE_FRAME_NR, itest);
			if (itest != movie_frame_numbers[iframe])
				REPORT_ERROR("ParticlePolisher::readStarFileBfactors BUG: itest= " + integerToString(itest) + " != " + integerToString(movie_frame_numbers[iframe]));
			MD.getValue(EMDL_POSTPROCESS_BFACTOR, DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 0) );
			MD.getValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 1) );
			iframe++;
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

	for (int iframe = 0; iframe < movie_frame_numbers.size(); iframe++ )
	{
		MDout.addObject();
		MDout.setValue(EMDL_IMAGE_FRAME_NR, movie_frame_numbers[iframe]);
		MDout.setValue(EMDL_POSTPROCESS_BFACTOR, DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 0) );
		MDout.setValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 1) );
	}

	MDout.write(fn_star);

	CPlot2D *plot2D=new CPlot2D("Polishing B-factors");
	plot2D->SetXAxisSize(600);
	plot2D->SetYAxisSize(400);
	plot2D->SetDrawLegend(false);
	plot2D->SetXAxisTitle("movie frame");
	plot2D->SetYAxisTitle("B-factor");
	MDout.addToCPlot2D(plot2D, EMDL_IMAGE_FRAME_NR, EMDL_POSTPROCESS_BFACTOR);
	plot2D->OutputPostScriptPlot(fn_out + "bfactors.eps");

	CPlot2D *plot2Db=new CPlot2D("Polishing scale-factors");
	plot2Db->SetXAxisSize(600);
	plot2Db->SetYAxisSize(400);
	plot2Db->SetDrawLegend(false);
	plot2Db->SetXAxisTitle("movie frame");
	plot2Db->SetYAxisTitle("Scale-factor");
	MDout.addToCPlot2D(plot2Db, EMDL_IMAGE_FRAME_NR, EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT);
	plot2Db->OutputPostScriptPlot(fn_out + "scalefactors.eps");


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
    for (int iframe = 0; iframe < movie_frame_numbers.size(); iframe++ )
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
	for (int iframe = 0; iframe < movie_frame_numbers.size(); iframe++ )
	{
		RFLOAT bfactor = DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 0);
		RFLOAT offset = DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 1);

		MetaDataTable MDout;
		std::string fn_table = "relative_weights_frame_" + integerToString(movie_frame_numbers[iframe]);
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
	fn_sum = fn_out + "avgframes_half" + integerToString(this_half) + "_class001_unfil.mrc";

	if (only_do_unfinished && exists(fn_sum))
	{
		if (verb > 0)
			std::cout << std::endl << " + " << fn_sum << " already exists: skipping calculation average of all frames." << std::endl;
		return;
	}

	Image<RFLOAT> Isum, Ione;
	for (long int iframe = 0; iframe < movie_frame_numbers.size(); iframe++)
	{

		int this_frame = movie_frame_numbers[iframe];
		FileName fn_vol;
		fn_vol.compose(fn_out + "frame", this_frame, "", 3);
		fn_vol += "_half" + integerToString(this_half) + "_class001_unfil.mrc";

		if (iframe == 0)
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

void ParticlePolisher::calculateSingleFrameReconstruction(int iframe, int this_half)
{

	// convert iframe to actual frame number in the STAR files!
	int this_frame = movie_frame_numbers[iframe];

	FileName fn_vol;
	fn_vol.compose(fn_out + "frame", this_frame, "", 3);
	fn_vol += "_half" + integerToString(this_half) + "_class001_unfil.mrc";
	if (only_do_unfinished && exists(fn_vol))
	{
		if (verb > 0)
			std::cout << std::endl << " + " << fn_vol << " already exists: skipping per-frame reconstruction." << std::endl;
		return;
	}

	// Set current_size
	int current_size;
	if (perframe_highres > 0.)
		current_size = 2 * ROUND(ori_size * angpix / perframe_highres);
	else
		current_size = ori_size;
	BackProjector backprojector(ori_size, 3, fn_sym);
	backprojector.initZeros(current_size);

	// Loop over all individual micrographs
	for (long int imic = 0; imic < fn_mics.size(); imic++)
	{
		FileName fn_fit = (fitting_mode == NO_FIT) ? fn_mics[imic] : fn_mics[imic].withoutExtension() + "_fit.star";
		Experiment exp_model;
		exp_model.read(fn_fit);

		CTF ctf;
		std::string dum;
		Matrix2D<RFLOAT> A3D;
		MultidimArray<Complex > Faux, F2D, F2Dp, Fsub;
		MultidimArray<RFLOAT> Fweight, Fctf;
		Image<RFLOAT> img;
		RFLOAT xtrans, ytrans;
		RFLOAT rot, tilt, psi;
		int i_half;
		long int my_frame;
		FileName fn_img, fn_mic;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(exp_model.MDimg)
		{

			exp_model.MDimg.getValue(EMDL_PARTICLE_RANDOM_SUBSET, i_half);
			exp_model.MDimg.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
			fn_mic.decompose(my_frame, dum);
			if (ABS(my_frame - this_frame) <= frame_running_average/2 && i_half == this_half)
			{
				exp_model.MDimg.getValue(EMDL_IMAGE_NAME, fn_img);
				img.read(fn_img);
				CenterFFT(img(), true);
				FourierTransformer transformer;
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
				exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, xtrans);
				exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, ytrans);
				if (ABS(xtrans) > 0. || ABS(ytrans) > 0. )
					shiftImageInFourierTransform(F2D, F2D, ori_size, xtrans, ytrans );

				// CTF
				Fctf.resize(F2D);
				Fctf.initConstant(1.);
				if (do_ctf)
				{
					ctf.read(exp_model.MDimg, exp_model.MDimg);
					ctf.getFftwImage(Fctf, ori_size, ori_size, angpix, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
					{
						DIRECT_MULTIDIM_ELEM(F2D, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
						DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
				}

				backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
			}

		} // end loop over all movie frames in exp_model
	} // end loop over all micrographs in fn_mics

	backprojector.symmetrise(helical_nr_asu, helical_twist, helical_rise / angpix);

	// Now do the reconstruction
	MultidimArray<RFLOAT> dummy;
	Image<RFLOAT> vol;
	backprojector.reconstruct(vol(), 10, false, 1., dummy, dummy, dummy, dummy, dummy);

	vol.write(fn_vol);

}

void ParticlePolisher::calculateBfactorSingleFrameReconstruction(int iframe, RFLOAT &bfactor, RFLOAT &offset, RFLOAT &corr_coeff)
{

	if (NZYXSIZE(Imask()) == 0)
		REPORT_ERROR("ParticlePolisher::calculateBfactorSingleFrameReconstruction BUG: first read Imask, then call this function.");

	if (iframe > 0 && XSIZE(fsc_average) == 0)
		REPORT_ERROR("ParticlePolisher::calculateBfactorSingleFrameReconstruction BUG: first call this function with this_frame < 0 to calculate FSC_average.");


		// Read in the 2 half-reconstructions and calculate their FSC
	FileName fn_root_half;
	// Make sure that the first call to this function is with this_frame < 0!!!
	int this_frame;
	if (iframe < 0)
	{
		fn_root_half = fn_out+"avgframes";
	}
	else
	{
		this_frame = movie_frame_numbers[iframe];
		fn_root_half.compose(fn_out+"frame",this_frame,"", 3);
	}

	FileName fn_half1, fn_half2;
	Image<RFLOAT> I1, I2;
	MultidimArray<RFLOAT> fsc_frame;
	I1.read(fn_root_half + "_half1_class001_unfil.mrc");
	I2.read(fn_root_half + "_half2_class001_unfil.mrc");

	// Mask both maps
	I1() *= Imask();
	I2() *= Imask();
	getFSC(I1(), I2(), fsc_frame);

	if (iframe < 0)
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

		// Check if any points were included in the Guinier plot
		if (guinier.size() < 3)
		{
			std::cerr << " WARNING: insufficient number of points in the Guinier plot of movie frame: " << this_frame << std::endl;
			std::cerr << " Consider lowering the lowres-limit, or average over multiple frames in the B-factor estimation." << std::endl;
		}

		// Now do the fit
		fitStraightLine(guinier, bfactor, offset, corr_coeff);
		// this is the B-factor relative to the average from all single-frame reconstructions!
		// in this case: positive values mean BETTER than average, and thus HIGHER WEIGHTS!
		bfactor *= 4.;

		CDataSet dataSet;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDout)
		{
			RFLOAT res2;
			MDout.getValue(EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, res2);
			CDataPoint point(res2, offset + (res2 * bfactor / 4.));
			dataSet.AddDataPoint(point);
		}
		dataSet.SetDatasetColor(1., 0., 0.);
		dataSet.SetDatasetTitle("Fitted straight line");
		CPlot2D *plot2D=new CPlot2D("Guinier plot frame " + integerToString(iframe+1));
		plot2D->SetXAxisSize(600);
		plot2D->SetYAxisSize(400);
		plot2D->SetXAxisTitle("resolution (1/A^2)");
		plot2D->SetYAxisTitle("ln(amplitudes)");
		MDout.addToCPlot2D(plot2D, EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, EMDL_POSTPROCESS_GUINIER_VALUE_IN);
		plot2D->AddDataSet(dataSet);
		plot2D->OutputPostScriptPlot(fn_out + "frame_"+integerToString(iframe+1, 3, '0')+"_guinier.eps");
	}

}

void ParticlePolisher::polishParticlesAllMicrographs()
{

	if (only_do_unfinished && exists(fn_out + "shiny.star"))
	{
		if (verb > 0)
			std::cout << std::endl << " + " << fn_out << "shiny.star already exists: skipping polishing of the particles." << std::endl;

		return;
	}

	// Loop over all average micrographs
	int barstep;
	long int my_nr_micrographs = fn_mics.size();
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

void ParticlePolisher::changeParticleStackName(FileName &fn_part)
{

	long int nr;
	FileName fn_stack;
	fn_part.decompose(nr, fn_stack);
	FileName fn_pre, fn_jobnr, fn_post;
	decomposePipelineFileName(fn_stack, fn_pre, fn_jobnr, fn_post);
	fn_part = fn_out + fn_post;
}

void ParticlePolisher::writeStarFilePolishedParticles()
{

	// Also write the STAR file with the MetaData of all polished particles
	// Loop over all original_particles in this micrograph

	MDshiny.clear();
	for (long int imic = 0; imic < fn_mics.size(); imic++)
	{
		FileName fn_fit = (fitting_mode == NO_FIT) ? fn_mics[imic] : fn_mics[imic].withoutExtension() + "_fit.star";
		Experiment exp_model;
		exp_model.read(fn_fit);

		// Just testing we've read a single micrograph
		if (exp_model.micrographs.size()>1)
			REPORT_ERROR("BUG: exp_model.micrographs.size()= " + integerToString(exp_model.micrographs.size()));

		for (long int ipar = 0; ipar < exp_model.micrographs[0].ori_particle_ids.size(); ipar++)
		{
			long int ori_part_id = exp_model.micrographs[0].ori_particle_ids[ipar];
			long int part_id = exp_model.ori_particles[ori_part_id].particles_id[0];

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
			changeParticleStackName(fn_part);
			fn_img.compose(ipar + 1, fn_part);
			MDshiny.setValue(EMDL_IMAGE_NAME, fn_img);
		}
	}

	// Deactivate PRIORs so they are not written out
	MDshiny.deactivateLabel(EMDL_ORIENT_ORIGIN_X_PRIOR);
	MDshiny.deactivateLabel(EMDL_ORIENT_ORIGIN_Y_PRIOR);
	MDshiny.deactivateLabel(EMDL_ORIENT_ROT_PRIOR);
	if (!is_helix)
	{
		MDshiny.deactivateLabel(EMDL_ORIENT_TILT_PRIOR);
		MDshiny.deactivateLabel(EMDL_ORIENT_PSI_PRIOR);
	}

	// Write output metadatatable
	MDshiny.write(fn_out + "shiny.star");
	std::cout << " + Written out all polished particles and their corresponding STAR file: " << fn_out << "shiny.star" << std::endl;

}

void ParticlePolisher::polishParticlesOneMicrograph(long int imic)
{

	FileName fn_fit = (fitting_mode == NO_FIT) ? fn_mics[imic] : fn_mics[imic].withoutExtension() + "_fit.star";
	Experiment exp_model;
	exp_model.read(fn_fit);

	// Just testing we've read a single micrograph
	if (exp_model.micrographs.size()>1)
		REPORT_ERROR("BUG: exp_model.micrographs.size()= " + integerToString(exp_model.micrographs.size()));

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

	// Loop over all original_particles in this micrograph
	for (long int ipar = 0; ipar < exp_model.micrographs[0].ori_particle_ids.size(); ipar++)
	{
		long int ori_part_id = exp_model.micrographs[0].ori_particle_ids[ipar];

		// Loop over all frames for motion corrections and possibly dose-dependent weighting
		for (long int iframe = 0; iframe < exp_model.ori_particles[ori_part_id].particles_id.size(); iframe++)
		{
			long int part_id = exp_model.ori_particles[ori_part_id].particles_id[iframe];

			exp_model.MDimg.getValue(EMDL_IMAGE_NAME, fn_img, part_id);
			exp_model.MDimg.getValue(EMDL_PARTICLE_ORI_NAME, fn_part, part_id);
			exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_X_PRIOR, x_off_prior_p, part_id);
			exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y_PRIOR, y_off_prior_p, part_id);
			exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, x_off_p, part_id);
			exp_model.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, y_off_p, part_id);

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
				std::cerr << " iframe= " << iframe << " XX(trans)= " << XX(trans) << " YY(trans)= " << YY(trans) << std::endl;
			}
#endif

			// Apply the phase shifts for this translation in Fourier space
			transformer.FourierTransform(img(), Fimg);

			if (iframe == 0)
			{
				Fwsum.initZeros(Fimg);
				Fsumw.initZeros(Fimg);
			}

			shiftImageInFourierTransform(Fimg, Fimg, ori_size, xtrans, ytrans);

			// Apply (positive!!) B-factor weighting and store weighted sums
			RFLOAT bfactor = (do_weighting) ? DIRECT_A1D_ELEM(perframe_bfactors, 3 * iframe + 0) : 0.;
			RFLOAT offset  = (do_weighting) ? DIRECT_A1D_ELEM(perframe_bfactors, 3 * iframe + 1) : 0.;
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
		changeParticleStackName(fn_part);
		fn_img.compose(ipar + 1, fn_part);

		// Only make directory if needed
		if (ipar == 0)
		{
			FileName fn_dir = fn_part.beforeLastOf("/");
			if (fn_dir != fn_olddir)
			{
				// Make a Particles directory
				int res = system(("mkdir -p " + fn_dir).c_str());
				fn_olddir = fn_dir;
			}
		}

		// When last particle, also write the correct header
		if (ipar == exp_model.micrographs[0].ori_particle_ids.size() - 1)
		{
			all_avg /= exp_model.micrographs[0].ori_particle_ids.size();
			all_stddev = sqrt(all_stddev/exp_model.micrographs[0].ori_particle_ids.size());
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
	Experiment exp_model;
	exp_model.read(fn_out + "shiny.star", true);

	// Do the reconstructions for both halves
	reconstructShinyParticlesOneHalf(1, exp_model);
	reconstructShinyParticlesOneHalf(2, exp_model);


	FileName fn_post = (ipass == 1) ? "shiny_post" : "shiny_post2";
	if (only_do_unfinished && exists(fn_out + fn_post + "_masked.mrc")
					       && exists(fn_out + fn_post + ".star") )
	{
		if (verb > 0)
			std::cout << std::endl << " + " << fn_out << fn_post << "_masked.mrc already exists: re-reading map into memory." << std::endl;

		if (verb > 0)
			std::cout << std::endl << " + " << fn_out << fn_post << ".star already exists: re-reading resolution from it." << std::endl;

		MetaDataTable MD;
		MD.read(fn_out + fn_post + ".star", "general");
		MD.getValue(EMDL_POSTPROCESS_FINAL_RESOLUTION, maxres_model);
	}
	else
	{
		// Re-read the two halves to calculate FSCs
		Postprocessing prm;

		prm.clear();
		prm.fn_in = fn_out + "shiny";
		prm.fn_out = fn_out + fn_post;
		prm.angpix = angpix;
		prm.do_auto_mask = false;
		prm.fn_mask = fn_mask;
		prm.do_auto_bfac = false;
		prm.do_fsc_weighting = true;
		prm.verb = 0;
		prm.run();

		maxres_model = prm.global_resol;
	}
	std::cout << " Resolution of reconstructions from shiny particles: " << maxres_model << std::endl;
	std::cout << " But you probably want to re-run at least a 3D auto-refinement with the shiny particles." << std::endl;


}

void ParticlePolisher::reconstructShinyParticlesOneHalf(int this_half, Experiment &exp_model)
{


	FileName fn_vol = fn_out + "shiny_half" + integerToString(this_half) + "_class001_unfil.mrc";
	if (only_do_unfinished && exists(fn_vol))
	{
		if (verb > 0)
			std::cout << std::endl << " + " << fn_vol << " already exists: skipping shiny-particle reconstruction." << std::endl;
		return;
	}

	// get image size, angpix (from metadatatable), fn_sym
	BackProjector backprojector(ori_size, 3, fn_sym);
	backprojector.initZeros();
	Projector projector(ori_size);

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
				shiftImageInFourierTransform(F2D, F2D, ori_size, xtrans, ytrans );

			// CTF
			Fctf.resize(F2D);
			Fctf.initConstant(1.);
			if (do_ctf)
			{
				ctf.read(exp_model.MDimg, exp_model.MDimg);
				ctf.getFftwImage(Fctf, ori_size, ori_size, angpix, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
				{
					DIRECT_MULTIDIM_ELEM(F2D, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
				}
			}

			backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
		}

	}


	backprojector.symmetrise(helical_nr_asu, helical_twist, helical_rise / angpix);

	// Now do the reconstruction
	MultidimArray<RFLOAT> dummy;
	backprojector.reconstruct(vol(), 10, false, 1., dummy, dummy, dummy, dummy, dummy);

	vol.write(fn_vol);

}

void ParticlePolisher::generateLogFilePDF()
{

	if (!exists(fn_out + "logfile.pdf"))
	{
		std::vector<FileName> fn_eps;

		fn_eps.push_back(fn_out + "bfactors.eps");
		fn_eps.push_back(fn_out + "scalefactors.eps");
		fn_eps.push_back(fn_out + "frame_???_guinier.eps");

		FileName fn_prev="";
		for (long int i = 0; i < fn_mics.size(); i++)
		{
			if (fn_prev != fn_mics[i].beforeLastOf("/"))
			{
				fn_prev = fn_mics[i].beforeLastOf("/");
				fn_eps.push_back(fn_prev+"/*.eps");
			}
		}

		joinMultipleEPSIntoSinglePDF(fn_out + "logfile.pdf ", fn_eps);

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

	// Make a logfile in pdf format
	generateLogFilePDF();

	// Write out the intermediately polished particles
	polishParticlesAllMicrographs();

	// Now reconstruct with all polished particles: two independent halves, FSC-weighting of the sum of the two...
	reconstructShinyParticlesAndFscWeight(1);

	if (verb > 0)
		std::cout << " done!" << std::endl;

}


