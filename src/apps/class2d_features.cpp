/*
 * class2D_features.cpp
 *
 *  Created on: 20 Oct 2019
 *      Author: ldong
 */

/*
 * fourier_features.cpp
 *
 *  Created on: 24 Jul 2018
 *      Author: ldong
 */
/*
 * calculatePvsLBP() comes from Xmipp source code,
 * https://github.com/I2PC/xmipp/blob/devel/libraries/reconstruction/classify_extract_features.cpp
 * following is the license.
 */
/***************************************************************************
 *
 * Authors:    Tomas Majtner            tmajtner@cnb.csic.es (2017)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <math.h>
#include <src/args.h>
#include <src/image.h>
#include <src/fftw.h>
#include <src/time.h>
#include <src/ml_model.h>
#include <src/ml_optimiser.h>
#include <src/exp_model.h>
#include <src/ctf.h>
#include <stdio.h>


// This contains 3 moments for an image
class moments {
public:
    RFLOAT mean, stddev, skew, kurt;
    moments(): mean(0), stddev(0), skew(0), kurt(0)
    {
    }
};

// This defines all information about a single class
class class_features {
public:
	// Class-wise features
    FileName name;
    Image<RFLOAT> img;
    long class_index;
    int is_selected, resol_limit;
    RFLOAT class_distribution, accuracy_rotation, accuracy_translation, estimated_resolution, particle_nr;
    RFLOAT class_score, edge_signal, scattered_signal, weighted_resolution;
    RFLOAT lowpass_filtered_img_avg, lowpass_filtered_img_stddev, lowpass_filtered_img_minval, lowpass_filtered_img_maxval;
    std::vector<RFLOAT> resolutions, lbp, lbp_p, lbp_s;
    moments circular_mask_moments, ring_moments, inner_circle_moments, fft_moments, protein_moments, solvent_moments;

    // Job-wise features
    RFLOAT PixelSize, SigmaOffSets, AveragePmax, ParticleDiameter, HighresLimitExpectation, job_score;
    int NrClasses, CurrentIteration, OriginalImageSize;

    class_features(): name(""),
    		class_index(0),
			is_selected(0),
			particle_nr(0),
    		class_distribution(0),
			accuracy_rotation(0),
			accuracy_translation(0),
			estimated_resolution(999.0),
			weighted_resolution(999.0),
			PixelSize(0),
			OriginalImageSize(0),
			SigmaOffSets(0),
			AveragePmax(0),
			ParticleDiameter(0),
			HighresLimitExpectation(0),
			NrClasses(0),
			CurrentIteration(0),
			class_score(-1),
			job_score(-1),
			resol_limit(-1),
			edge_signal(-1.),
			scattered_signal(-1.),
			lowpass_filtered_img_avg(0),
			lowpass_filtered_img_stddev(0),
			lowpass_filtered_img_minval(0),
			lowpass_filtered_img_maxval(0)
    {
    }
};


class liyi_class_features {

public:
	IOParser parser;
	FileName fn_out, fn_optimiser, fn_model, fn_select, fn_job_score, fn_cf;

	RFLOAT minRes, job_score;
	RFLOAT radius_ratio, radius;
	RFLOAT circular_mask_radius, uniform_angpix = 4.0;
	RFLOAT binary_threshold, lowpass;
    int debug;

    // Save some time by limting calculations
	int only_use_this_class;
	bool do_skip_angular_errors, do_skip_protein_vs_solvent, do_skip_LBP;

	MlOptimiser myopt;
	MetaDataTable MD_optimiser, MD_select;
	std::vector<class_features> features_all_classes;

	bool do_save_masks, save_masks_only;


	liyi_class_features(): job_score(-1),
		minRes(-1)
	{
	}

	/** ========================== I/O operations  =========================== */

//read class averages
	void read(int argc, char **argv){

			parser.setCommandLine(argc, argv);

			// TODO: optional input files, eg. job score file
	    	int general_section = parser.addSection("General options");
	    	fn_optimiser = parser.getOption("--opt", "Input optimiser.star file", "");
			fn_out = parser.getOption("--o", "Name for output class features", "cf_lbp_rescaled.star");
	    	fn_select = parser.getOption("--select", "Input class_averages.star from the Selection job or backup_selection.star", "");
	    	fn_job_score = parser.getOption("--fn_score", "Input job score file", "");

	    	int part_section = parser.addSection("Partial calculations");
			fn_cf = parser.getOption("--cf_file", "Input class feature star file", "");
	    	only_use_this_class = textToInteger(parser.getOption("--only_class_nr", "Class number of the class of interest", "-1"));
	    	do_skip_angular_errors = parser.checkOption("--skip_angular_errors", "Skip angular error calculation");
	    	do_skip_protein_vs_solvent = parser.checkOption("--skip_pvs", "Skip protein-solvent mask calculation");
	    	do_skip_LBP = parser.checkOption("--skip_lbp", "Skip LBP calculation");

	    	int expert_section = parser.addSection("Expert options");
	    	radius_ratio = textToFloat(parser.getOption("--radius_ratio", "Ratio of inner radius of the interested ring area in proportion to the current circular mask radius", "0.95"));
			radius = textToFloat(parser.getOption("--radius", "Inner radius of the interested ring area to the current circular mask radius", "0"));
			binary_threshold = textToFloat(parser.getOption("--binary_threshold", "Threshold for generating binary masks.", "0."));
			lowpass = textToFloat(parser.getOption("--lowpass", "Low-pass filter frequency (in A)", "25."));
			do_save_masks = parser.checkOption("--save_masks", "Save the protein and solvent masks ");
		    debug = textToInteger(parser.getOption("--debug", "Debug level", "0"));

	    	// Check for errors in the command-line option
	    	if (radius > 0 && radius < 1)
	    		REPORT_ERROR("Incorrect value(s) for inner_radius / outer_radius. Try inner_radius_ratio / outer_radius_ratio.");
	    	if (parser.checkForErrors())
	    		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	/** ========================== Initialisation  =========================== */

	void initialise(){

		if (do_skip_angular_errors || do_skip_LBP || do_skip_protein_vs_solvent)
		{
			if (fn_cf == "") REPORT_ERROR("ERROR: you need to provide a class feature input file if you wish to skip some calculations!");
		}

		// Read in the MD_optimiser table from the STAR file, get model.star and data.star
		if (fn_optimiser != "")
		{
			myopt.read(fn_optimiser); // true means skip_groups_and_pdf_direction from mlmodel; only read 1000 particles...
			if (debug>0) std::cerr << "Done with reading optimiser ..." << std::endl;

			// A bit ugly, but we need fn_model...
			MetaDataTable MDopt;
			MDopt.read(fn_optimiser);
			MDopt.getValue(EMDL_OPTIMISER_MODEL_STARFILE, fn_model);

			if (myopt.intact_ctf_first_peak)
			{
				if (debug>0) std::cerr << "Doing first peak CTF correction ..." << std::endl;

				// Calculate avg. defocus
				RFLOAT def_avg = 0, def_u, def_v;
				for (long int part_id = 0; part_id < myopt.mydata.MDimg.numberOfObjects(); part_id++)
				{
					myopt.mydata.MDimg.getValue(EMDL_CTF_DEFOCUSU, def_u, part_id);
					myopt.mydata.MDimg.getValue(EMDL_CTF_DEFOCUSV, def_v, part_id);
					def_avg += def_u + def_v;
				}
				def_avg /= (2. * myopt.mydata.MDimg.numberOfObjects());

				// some people may have used a too small , or even zero amplitude contrast, lets forbid that...
				//q0 = XMIPP_MAX(q0, 0.07);
				CTF avgctf;
				avgctf.setValuesByGroup(&myopt.mydata.obsModel, 0, def_u, def_v, 0.);

				// Loop over all classes in myopt.mymodel.Iref
				for (long iref =0; iref < myopt.mymodel.Iref.size(); iref++)
				{
					correctCtfUntilFirstPeak(myopt.mymodel.Iref[iref], avgctf);
				}
			}
			myopt.mymodel.setFourierTransformMaps(false);

			// Calculate mask radius
			circular_mask_radius = (myopt.particle_diameter * 0.5) / myopt.mymodel.pixel_size;

		}

		// Open selected_class table
		if (fn_select != "")
		{
			MD_select.read(fn_select);
		}

		// Read in class features from a previous run if fn_cf is provided
		readClassFeatures();

		// Get job score
		if (fn_job_score != "")
		{
		    std::ifstream in(fn_job_score, std::ios_base::in);
		    if (in.fail()) REPORT_ERROR( (std::string) "ERROR: File " + fn_job_score + " does not exists" );
		    std::string line;
		    getline(in, line, '\n');
		    job_score = textToFloat(line);
//		    std::cout << " Read job score = " << job_score << std::endl;
		    in.close();
		}
		else
		{
			job_score = -1.;
		}

		if (radius_ratio > 0 && radius > 0)
		{
			std::cout << "WARNING: Should not provide radius ratio and radius at the same time. Ignoring the radius ratio..." << std::endl;
		}
	}

	/** ========================== Moments Calculators ===========================  */

	// Calculate moments for one class
	void meanMomentsCalculator (MultidimArray<RFLOAT> &img, RFLOAT inner_radius, RFLOAT outer_radius, moments &mmts){
		int pix_num = 0;
		RFLOAT inner_radius_square, outer_radius_square;

		inner_radius_square = pow((inner_radius), 2);
		outer_radius_square = pow((outer_radius), 2);

		// Calculating image means
		RFLOAT sum = 0;
		for (long int i= STARTINGY(img) ; i <= FINISHINGY(img); i++)
		{
			for (long int j= STARTINGX(img) ; j <= FINISHINGX(img); j++)
			{
				if (i*i+j*j >= inner_radius_square && i*i+j*j <= outer_radius_square)
				{
					pix_num++;
					sum += A2D_ELEM(img, i, j);
				}
			}
		}
		mmts.mean = sum / pix_num;

		// Calculating image moments
		RFLOAT square_sum = 0, cube_sum = 0, quad_sum = 0;
		for (long int i= STARTINGY(img) ; i <= FINISHINGY(img); i++)
		{
			for (long int j= STARTINGX(img) ; j <= FINISHINGX(img); j++)
			{
				if (i*i+j*j >= inner_radius_square && i*i+j*j <= outer_radius_square)
				{
					square_sum += pow((A2D_ELEM(img, i, j)-mmts.mean), 2);
					cube_sum += pow((A2D_ELEM(img, i, j)- mmts.mean), 3);
					quad_sum += pow((A2D_ELEM(img, i, j)- mmts.mean), 4);
				}
			}
		}

		mmts.stddev = sqrt(square_sum / pix_num);
		if (square_sum == 0)
		{
			mmts.skew = 0;
			mmts.kurt = 0;
		} else
		{
			mmts.skew = (cube_sum * sqrt (pix_num))/(pow(square_sum, 1.5));
			mmts.kurt = (quad_sum * pix_num)/pow(square_sum, 2);
		}

	}

	/** =================================================== Getting Other Class Metadata ======================================================  */

	// Find estimated resolution from MD_model_classes_ table for one class
	void findResolution(class_features &cf)
	{

		std::string class_name = "model_class_"+integerToString(cf.class_index);
		MetaDataTable MD_model_this_class;
		MD_model_this_class.read(fn_model, class_name);
		bool isGood = true;

		// Test for numberOfObjects() function
//		std::cout <<"Number of spectral of class " << cf.class_index <<": "<< MD_model_this_class.numberOfObjects() << std::endl;
		for (int i=1; i<MD_model_this_class.numberOfObjects();i++)
		{
			RFLOAT ssnrMap;
			MD_model_this_class.getValue(EMDL_MLMODEL_DATA_VS_PRIOR_REF, ssnrMap, i);
			if (ssnrMap < 1)
			{
				MD_model_this_class.getValue(EMDL_RESOLUTION_ANGSTROM, cf.estimated_resolution, i-1);
				isGood = false;
				break;
			}
		}
		if (isGood)
		{
			cf.estimated_resolution = myopt.mymodel.pixel_size*2;
		}
	}

	// Calculate accuracy rotation and translation for each of the non-empty classes
	void calculateExpectedAngularErrors(int iclass, class_features &cf)
	{
		// Set current_image_size to the coarse_size to calculate expected angular errors
		int current_image_size = myopt.mymodel.current_size;

		// Separate angular error estimate for each of the classes
		RFLOAT acc_rot = 999., acc_trans = 999.;

		// P(X | X_1) / P(X | X_2) = exp ( |F_1 - F_ 2|^2 / (-2 sigma2) )
		// exp(-4.60517) = 0.01
		RFLOAT pvalue = 4.60517;

		// Randomise particle orders only the first time
		//if (iclass == 0) myopt.mydata.randomiseParticlesOrder(0, false, false);

		// calculate acc rot and trans for large classes (particle number > 100)
		// for small classes, set acc rot to 5 degrees and acc trans to 8 pixels
		if (cf.particle_nr <= 100. )
		{
			cf.accuracy_rotation = 5.;
			cf.accuracy_translation = 8.;
		}
		else
		{
			RFLOAT acc_rot_class = 0.;
			RFLOAT acc_trans_class = 0.;
			// Particles are already in random order, so just move from 0 to n_trials
			//LOOP OVER 100 RANDOM PARTICLES HERE
			int n_trials = 100;
			for (long int part_id = 0; part_id <= n_trials; part_id++)
			{
				int group_id = myopt.mydata.getGroupId(part_id, 0);
				RFLOAT my_pixel_size = myopt.mydata.getImagePixelSize(part_id, 0);
				const int optics_group = myopt.mydata.getOpticsGroup(part_id, 0);
				int my_image_size = myopt.mydata.getOpticsImageSize(optics_group);
				bool ctf_premultiplied = myopt.mydata.obsModel.getCtfPremultiplied(optics_group);

				MultidimArray<RFLOAT> Fctf;
				// Get CTF for this particle
				if (myopt.do_ctf_correction)
				{
					Fctf.resize(current_image_size, current_image_size/ 2 + 1);

					// Get parameters that change per-particle from the exp_metadata
					CTF ctf;
					RFLOAT def_u, def_v, def_angle, voltage, cs, q0;
					myopt.mydata.MDimg.getValue(EMDL_CTF_DEFOCUSU, def_u, part_id);                 //??
					myopt.mydata.MDimg.getValue(EMDL_CTF_DEFOCUSV, def_v, part_id);
					myopt.mydata.MDimg.getValue(EMDL_CTF_DEFOCUS_ANGLE, def_angle, part_id);
					ctf.setValuesByGroup(&myopt.mydata.obsModel, optics_group, def_u, def_v, def_angle);
					ctf.getFftwImage(Fctf, my_image_size, my_image_size, myopt.mymodel.pixel_size,
							myopt.ctf_phase_flipped, myopt.only_flip_phases, myopt.intact_ctf_first_peak, true, myopt.do_ctf_padding);
				}
				// Search 2 times: ang and off
				for (int imode = 0; imode < 2; imode++)
				{
					RFLOAT ang_error = 0.;
					RFLOAT sh_error = 0.;
					RFLOAT ang_step;
					RFLOAT sh_step;
					RFLOAT my_snr = 0.;

					// Search for ang_error and sh_error where there are at least 3-sigma differences!
					// 13feb12: change for explicit probability at P=0.01
					while (my_snr <= pvalue)
					{
						// Gradually increase the step size
						if (ang_error < 0.2)
						  ang_step = 0.05;
						else if (ang_error < 1.)
						  ang_step = 0.1;
						else if (ang_error < 2.)
						  ang_step = 0.2;
						else if (ang_error < 5.)
						  ang_step = 0.5;
						else if (ang_error < 10.)
						  ang_step = 1.0;
						else if (ang_error < 20.)
						  ang_step = 2;
						else
						  ang_step = 5.0;

						if (sh_error < 1.)
						  sh_step = 0.1;
						else if (sh_error < 2.)
						  sh_step = 0.2;
						else if (sh_error < 5.)
						  sh_step = 0.5;
						else if (sh_error < 10.)
						  sh_step = 1.0;
						else
						  sh_step = 2.0;

						ang_error += ang_step;
						sh_error += sh_step;

						// Prevent an endless while by putting boundaries on ang_error and sh_error
						if ( (imode == 0 && ang_error > 30.) || (imode == 1 && sh_error > 10.) )
						  break;

						MultidimArray<Complex > F1, F2;
						Matrix2D<RFLOAT> A1, A2;

						// INSTEAD of below: get rot tilt and psi from data.star file for this 1/100 random particles
						RFLOAT rot1;
						RFLOAT tilt1;
						RFLOAT psi1;
						RFLOAT xoff1 = 0.;
						RFLOAT yoff1 = 0.;
						RFLOAT zoff1 = 0.;
						myopt.mydata.MDimg.getValue(EMDL_ORIENT_ROT, rot1, part_id);
						myopt.mydata.MDimg.getValue(EMDL_ORIENT_TILT, tilt1, part_id);
						myopt.mydata.MDimg.getValue(EMDL_ORIENT_PSI, psi1, part_id);

						F1.initZeros(current_image_size, current_image_size/ 2 + 1);

						// Get the FT of the first image
						Euler_angles2matrix(rot1, tilt1, psi1, A1, false);
						A1 = myopt.mydata.obsModel.applyAnisoMag(A1, optics_group);
						A1 = myopt.mydata.obsModel.applyScaleDifference(A1, optics_group, myopt.mymodel.ori_size, myopt.mymodel.pixel_size);
						(myopt.mymodel.PPref[iclass]).get2DFourierTransform(F1, A1);    //?
						// Apply the angular or shift error
						RFLOAT rot2 = rot1;
						RFLOAT tilt2 = tilt1;
						RFLOAT psi2 = psi1;
						RFLOAT xshift = xoff1;
						RFLOAT yshift = yoff1;
						RFLOAT zshift = zoff1;

						// Perturb psi or xoff , depending on the mode
						if (imode == 0)
						{
							if (myopt.mymodel.ref_dim == 3)
							{
								// Randomly change rot, tilt or psi
								RFLOAT ran = rnd_unif();
								if (ran < 0.3333)
								  rot2 = rot1 + ang_error;
								else if (ran < 0.6667)
								  tilt2 = tilt1 + ang_error;
								else
								  psi2  = psi1 + ang_error;
							}
							else
							{
								psi2  = psi1 + ang_error;
							}
						}
						else
						{
							// Randomly change xoff or yoff
							RFLOAT ran = rnd_unif();
							if (myopt.mymodel.data_dim == 3)
							{
								if (ran < 0.3333)
								  xshift = xoff1 + sh_error;
								else if (ran < 0.6667)
								  yshift = yoff1 + sh_error;
								else
								  zshift = zoff1 + sh_error;
							}
							else
							{
							   if (ran < 0.5)
								  xshift = xoff1 + sh_error;
								else
								  yshift = yoff1 + sh_error;
							}
						}

						// Get the FT of the second image
						if (myopt.mymodel.data_dim == 2)
						  F2.initZeros(current_image_size, current_image_size/ 2 + 1);
						else
						  F2.initZeros(current_image_size, current_image_size, current_image_size/ 2 + 1);

						if (imode == 0)
						{
							// Get new rotated version of reference
							Euler_angles2matrix(rot2, tilt2, psi2, A2, false);
							A2 = myopt.mydata.obsModel.applyAnisoMag(A2, optics_group);
							A2 = myopt.mydata.obsModel.applyScaleDifference(A2, optics_group, myopt.mymodel.ori_size, myopt.mymodel.pixel_size);
							(myopt.mymodel.PPref[iclass]).get2DFourierTransform(F2, A2);    //??
						}
						else
						{
							// Get shifted version
							shiftImageInFourierTransform(F1, F2, (RFLOAT) myopt.mymodel.ori_size, -xshift, -yshift, -zshift);
						}

						// Apply CTF to F1 and F2 if necessary
						if (myopt.do_ctf_correction)
						{
							FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
							{
								DIRECT_MULTIDIM_ELEM(F1, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
								DIRECT_MULTIDIM_ELEM(F2, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
							}

							if (myopt.mydata.hasCtfPremultiplied())
							{
								FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
									{
								  DIRECT_MULTIDIM_ELEM(F1, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
								  DIRECT_MULTIDIM_ELEM(F2, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
								}
							}
						}

						my_snr = 0.;
						FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(F1)
						{
							int ires = ROUND(sqrt((RFLOAT)(ip*ip + jp*jp)));
							if (ires < myopt.mymodel.ori_size / 2 + 1 && !(jp==0 && ip < 0))
							{
								my_snr += norm(DIRECT_A2D_ELEM(F1, i ,j) - DIRECT_A2D_ELEM(F2, i, j)) / (2 * myopt.mymodel.sigma2_noise[group_id](ires) );
							}
						}

					} // end while my_snr >= pvalue

					if (imode == 0)
					  acc_rot_class += ang_error;
					else if (imode == 1)
					  acc_trans_class += sh_error;
				} // end for imode
			} // end for part_id

			cf.accuracy_rotation = acc_rot_class / (RFLOAT)n_trials;
			cf.accuracy_translation = acc_trans_class / (RFLOAT)n_trials;

		} // end for if large class condition

	}


	// Read in job score txt file (if there is one) and calculate class score for one class
	void calculateClassScore(class_features &cf, RFLOAT minRes){
		RFLOAT weight;

		if (job_score != -1)
		{
			weight = pow((minRes/cf.estimated_resolution), 2);
	//		std::cout << "weight: " << weight << std::endl;
			switch (cf.is_selected)
			{
			case 1:
				cf.class_score = (0.75+weight*0.25)*job_score;
				break;
			case 2:
				cf.class_score = (0.25+weight*0.25)*job_score;
				break;
			case 3:
				cf.class_score = (0+weight*0.25)*job_score;
				break;
			case 4:
				cf.class_score = (0+weight*0.25)*job_score;
				std::cout << "Class " << cf.class_index << "is labelled cyan!" << std::endl;
				break;
			case 5:
				cf.class_score = (0.5+weight*0.25)*job_score;
				break;
			case 6:
				cf.class_score = 0.0;
				std::cout << "Class " << cf.class_index << "is labelled yellow!" << std::endl;
				break;
			case 0:
				cf.class_score = 0.0;
				break;
			default:
				std::cout << "Illegal selection label..." << std::endl;
			}
		}
		else
		{
			cf.class_score = -1;
		}
//		std::cout << "class score: " << cf.class_score << std::endl;
	}


	void makeFilteredMasks(MultidimArray<RFLOAT> img, MultidimArray<RFLOAT> &lpf, MultidimArray<RFLOAT> &p_mask, MultidimArray<RFLOAT> &s_mask,
			RFLOAT &scattered_signal, long &protein_area, long &solvent_area)
	{
		MultidimArray<bool> visited;
		std::stack<std::pair<long int, long int>> white_stack;
		std::vector<std::vector<std::pair<long int, long int>>> all_islands;
		long int scattered_signal_nr = 0;
		int dx[4] = {-1, 1, 0, 0};
		int dy[4] = {0, 0, -1, 1};

		lpf = img;
		lpf.setXmippOrigin();
		visited.clear();
		visited.resize(lpf);
		visited.initConstant(false);
		p_mask.setXmippOrigin();
		s_mask.setXmippOrigin();
		p_mask.initZeros(lpf);
		s_mask.initZeros(lpf);
		protein_area = 0;
		long circular_area = 0;

		// A hyper-parameter to adjust: definition of central area: 0.7 of radius (~ half of the area)
		//		RFLOAT sq_ctr_r = 0.49*circular_mask_radius*circular_mask_radius;
		//		RFLOAT central_area = 3.14*0.49*circular_mask_radius*circular_mask_radius;

		lowPassFilterMap(lpf, lowpass, myopt.mymodel.pixel_size);

		for (long int i= STARTINGY(lpf) ; i <= FINISHINGY(lpf); i++)
		{
			for (long int j= STARTINGX(lpf) ; j <= FINISHINGX(lpf); j++)
			{
				if (i*i+j*j <= 0.49*circular_mask_radius*circular_mask_radius)
				{
					if (A2D_ELEM(visited, i, j) == false)
					{
						// Mark
						A2D_ELEM(visited, i, j) = true;
						if (A2D_ELEM(lpf, i, j) > binary_threshold)
						{                  // Find an initial white pixel and identify a new islands
							std::vector<std::pair<long int, long int>> island;
							long int inside = 1;
							white_stack.push(std::make_pair(i, j));
							island.push_back(std::make_pair(i, j));
							while (! white_stack.empty())
							{                             // For each island
								std::pair<long int, long int> center = white_stack.top();
								white_stack.pop();
								for (int a=0; a<4; a++)
								{                               // Explore neighbours
									long y = center.first + dy[a];
									long x = center.second + dx[a];
									// This boundary is added in addition to the particle diameter boundary encircled
									// in case some particle diameters are set even larger than the box sizes by accident
									if (y >= STARTINGY(lpf) && y<= FINISHINGY(lpf) && x >= STARTINGX(lpf) && x<= FINISHINGX(lpf))
									{
										if (y*y+x*x<=circular_mask_radius*circular_mask_radius && (A2D_ELEM(visited, y, x) == false))
										{
											A2D_ELEM(visited, y, x) = true;
											if (A2D_ELEM(lpf, y, x)> binary_threshold)
											{    // White neighbours
												white_stack.push(std::make_pair(y, x));
												island.push_back(std::make_pair(y, x));
												if (y*y+x*x <= 0.49*circular_mask_radius*circular_mask_radius)
												{
													inside++;                               // Count white pixel number inside the central area
												}
											}
										}
									}
								} // end of for loop for looking at 8 neighbours
							}
							if (inside > 0.2*3.14*0.49*circular_mask_radius*circular_mask_radius)
							{
								all_islands.push_back(island);
								protein_area += island.size();
							}
							else
							{
								scattered_signal_nr += island.size();
							}
						}
					}
				}
			}
		}
		// Calculate two member variables: scattered_signal and solvent_area
		scattered_signal = RFLOAT(scattered_signal_nr)/ RFLOAT((scattered_signal_nr + protein_area));
		// Make the masks
		for (long int i= STARTINGY(lpf) ; i <= FINISHINGY(lpf); i++)
		{
			for (long int j= STARTINGX(lpf) ; j <= FINISHINGX(lpf); j++)
			{
				if (i*i+j*j <= circular_mask_radius*circular_mask_radius)
				{
					A2D_ELEM(s_mask, i, j) = 1.;
					circular_area++;
				}
			}
		}
		solvent_area = circular_area - protein_area;
		for(int p=0; p<all_islands.size(); p++)
		{
			for (int q=0; q<all_islands[p].size(); q++)
			{
				long int i = all_islands[p][q].first;
				long int j = all_islands[p][q].second;
				A2D_ELEM(p_mask, i, j) = 1.;
				A2D_ELEM(s_mask, i, j) = 0.;
			}
		}
	}

	// protein mask filtered by DFS
	void proteinVsSolventFeatures(class_features &cf)
	{
		FileName fp_out, fs_out, flpf_out, f_original;

		MultidimArray<RFLOAT> lpf, p_mask, s_mask;

		// Make filtered masks and collect the class feature "scattered_signal"
		long protein_area, solvent_area;
		makeFilteredMasks(cf.img(), lpf, p_mask, s_mask, cf.scattered_signal, protein_area, solvent_area);

		if (do_save_masks)
		{
			// Create folders to save protein and solvent region masks
			char foldername[256];
			snprintf(foldername, 255, "pvs_masks_lp%.0fth%.3f", lowpass, binary_threshold);
			FileName filtered_mask_folder = foldername;
			if (!exists(filtered_mask_folder))
			{
				std::string filtered_mask_fn_command = "mkdir -p " + filtered_mask_folder;
				system(filtered_mask_fn_command.c_str());
			}

			Image<RFLOAT> p_out, s_out, lpf_out;
			lpf_out() = lpf;
			p_out() = p_mask;
			s_out() = s_mask;

			flpf_out = filtered_mask_folder+"/class"+integerToString(cf.class_index)+"_lowpassfiltered.mrc";
			fp_out = filtered_mask_folder+"/class"+integerToString(cf.class_index)+"_p_mask.mrc";
			fs_out = filtered_mask_folder+"/class"+integerToString(cf.class_index)+"_s_mask.mrc";
			f_original = filtered_mask_folder+"/class"+integerToString(cf.class_index)+"_original.mrc";

			lpf_out.write(flpf_out);
			p_out.write(fp_out);
			s_out.write(fs_out);
			cf.img.write(f_original);
		}

		// Store the mean, stddev, minval and maxval of the lowpassed image as features
		lpf.computeStats(cf.lowpass_filtered_img_avg, cf.lowpass_filtered_img_stddev,
				cf.lowpass_filtered_img_minval, cf.lowpass_filtered_img_maxval);

		// Apply masks and calculate features
		RFLOAT p_sum = 0, s_sum = 0;
		RFLOAT p_mean = 0, s_mean = 0;
		// Mean
		for (long int n = 0; n< YXSIZE(p_mask); n++)
		{
			if (p_mask.data[n] > 0.5)
			{
				p_sum += cf.img().data[n];
			}
			else if (s_mask.data[n] > 0.5)
			{
				s_sum += cf.img().data[n];
			}
		}
		p_mean = p_sum / protein_area;
		s_mean = s_sum / solvent_area;
		cf.protein_moments.mean = p_mean;
		cf.solvent_moments.mean = s_mean;

		// Moments
		RFLOAT p_square_sum = 0, p_cube_sum = 0, p_quad_sum = 0;
		RFLOAT s_square_sum = 0, s_cube_sum = 0, s_quad_sum = 0;
		for (long int n = 0; n< YXSIZE(p_mask); n++)
		{
			if (p_mask.data[n] > 0.5)
			{
				p_square_sum += pow((cf.img().data[n] - p_mean), 2);
				p_cube_sum += pow((cf.img().data[n] - p_mean), 3);
				p_quad_sum += pow((cf.img().data[n] - p_mean), 4);
			}
			else if (s_mask.data[n] > 0.5)
			{
				s_square_sum += pow((cf.img().data[n]- s_mean), 2);
				s_cube_sum += pow((cf.img().data[n] - s_mean), 3);
				s_quad_sum += pow((cf.img().data[n] - s_mean), 4);
			}
//			std::cout << cf.img().data[n] << "\t" << p_mean << std::endl;
		}

		cf.protein_moments.stddev = sqrt(p_square_sum / protein_area);
		cf.solvent_moments.stddev = sqrt(s_square_sum / solvent_area);
		if (p_square_sum == 0)
		{
			cf.protein_moments.skew = 0;
			cf.protein_moments.kurt = 0;
		}
		else
		{
			cf.protein_moments.skew = (p_cube_sum * sqrt (protein_area))/(pow(p_square_sum, 1.5));
			cf.protein_moments.kurt = (p_quad_sum * protein_area)/pow(p_square_sum, 2);
		}
		if (s_square_sum == 0)
		{
			cf.solvent_moments.skew = 0;
			cf.solvent_moments.kurt = 0;
		}
		else
		{
			cf.solvent_moments.skew = (s_cube_sum * sqrt (solvent_area))/(pow(s_square_sum, 1.5));
			cf.solvent_moments.kurt = (s_quad_sum * solvent_area)/pow(s_square_sum, 2);
		}

		// Edge signal
		long int edge_pix = 0, edge_white = 0;
		for (long int i= STARTINGY(p_mask) ; i <= FINISHINGY(p_mask); i++)
		{
			for (long int j= STARTINGX(p_mask) ; j <= FINISHINGX(p_mask); j++)
			{
//				long int r = round(sqrt((RFLOAT)(i * i + j * j)));
				if (round(sqrt((RFLOAT)(i * i + j * j))) == round(circular_mask_radius))
				{
					edge_pix++;
					if (A2D_ELEM(p_mask, i, j) > 0)
					{
						edge_white++;
					}
				}
			}
		}
		cf.edge_signal = RFLOAT(edge_white) / RFLOAT(edge_pix);

	}

	void calculatePvsLBP(MultidimArray<RFLOAT> Iin, std::vector<double> &lbp, std::vector<double> &lbp_p, std::vector<double> &lbp_s)
	{
		FileName img_name;
		MultidimArray<RFLOAT> I, lpf, p_mask, s_mask;
		std::vector<double> min_idxs, min_idxs_sort;

		// Re-scale the image to have uniform pixel size of 4 angstrom
		int newsize = ROUND(XSIZE(Iin) * (myopt.mymodel.pixel_size / uniform_angpix));
		newsize -= newsize%2; //make even in case it is not already
		I = Iin;
		resizeMap(I, newsize);

		// Reset pixel_size and mask diameter for the mask making
		//pixel_size = uniform_angpix;
		circular_mask_radius = myopt.particle_diameter / uniform_angpix * 0.5;

		// Make filtered masks
		RFLOAT dummy;
		long idummy;
		makeFilteredMasks(I, lpf, p_mask, s_mask, dummy, idummy, idummy);

		unsigned char code;
		double center;
		double lbp_hist[256] = {}, lbp_hist_p[256] = {}, lbp_hist_s[256] = {};

		// Make the map from 256 possible original rotation-variant LBP values to 36 rotation-invariant LBP values
		for (int i = 0; i < 256; i++)
		{
			code = i;
			int code_min = (int) code;
			for (int ii = 0; ii < 7; ii++)
			{
				unsigned char c = code & 1;
				code >>= 1;
				code |= (c << 7);
				if ((int) code < code_min)
					code_min = (int) code;
			}
			min_idxs.push_back(code_min);
		}
		min_idxs_sort = min_idxs;
		std::sort(min_idxs_sort.begin(), min_idxs_sort.end());
		std::unique(min_idxs_sort.begin(), min_idxs_sort.end());

		// Calculate rotation invariant LBP(8, 1) value for one pixel
		double sum=0., sum_p=0., sum_s=0.;
		for (int y = 1; y < (YSIZE(I)-1); y++)
		{
			for (int x = 1; x < (XSIZE(I)-1); x++)
			{
				code = 0;
				// Generating original LBP value (rotation variant) for a pixel
				center = DIRECT_A2D_ELEM(I,y,x);
				code |= (DIRECT_A2D_ELEM(I,y-1,x-1) > center) << 7;
				code |= (DIRECT_A2D_ELEM(I,y-1,x  ) > center) << 6;
				code |= (DIRECT_A2D_ELEM(I,y-1,x+1) > center) << 5;
				code |= (DIRECT_A2D_ELEM(I,y,  x+1) > center) << 4;
				code |= (DIRECT_A2D_ELEM(I,y+1,x+1) > center) << 3;
				code |= (DIRECT_A2D_ELEM(I,y+1,x  ) > center) << 2;
				code |= (DIRECT_A2D_ELEM(I,y+1,x-1) > center) << 1;
				code |= (DIRECT_A2D_ELEM(I,y  ,x-1) > center) << 0;
				// Map to rotation invariant value
				int idx = min_idxs[(int) code];
				if (DIRECT_A2D_ELEM(p_mask,y,x)> 0.5)
				{
					lbp_hist[idx]+=1.;
					sum += 1.;
					lbp_hist_p[idx]+=1.;
					sum_p += 1.;
				} else if (DIRECT_A2D_ELEM(s_mask,y,x)> 0.5)
				{
					lbp_hist[idx]+=1.;
					sum += 1.;
					lbp_hist_s[idx]+=1.;
					sum_s += 1.;
				}
			}
		}
		// Trim to include only the 36 rotation invariant LBP values
		 for (int i = 0; i < 36; i++)
		{
			int idx = min_idxs_sort[i];
			if (sum>0.) lbp_hist[idx] /= sum;
			if (sum_p>0.) lbp_hist_p[idx] /= sum_p;
			if (sum_s>0.) lbp_hist_s[idx] /= sum_s;
			lbp.push_back(lbp_hist[idx]);
			lbp_p.push_back(lbp_hist_p[idx]);
			lbp_s.push_back(lbp_hist_s[idx]);
		}
	}

/** ===========================================================Correct CTF until first peak ====================================================== */
	void correctCtfUntilFirstPeak(MultidimArray<RFLOAT> &in, CTF ctf)
	{

		FourierTransformer transformer;
		MultidimArray<Complex > Faux;

		RFLOAT xs = (RFLOAT)XSIZE(in) * myopt.mymodel.pixel_size;
		RFLOAT ys = (RFLOAT)YSIZE(in) * myopt.mymodel.pixel_size;
	    transformer.FourierTransform(in, Faux, false);

	    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(Faux)
	    {
			RFLOAT x = (RFLOAT)jp / xs;
			RFLOAT y = (RFLOAT)ip / ys;

			// ctf should become part of the pvs_features class, perhaps call it avgctf
			RFLOAT ctf_val = ctf.getCTF(x, y, false, false, false, false, 0., true);
			if (ctf_val > 0.)
	    		DIRECT_A2D_ELEM(Faux, i, j) /= ctf_val;
	    }
	    transformer.inverseFourierTransform(Faux, in);
	}

	/** ======================================================== Collecting All Features ========================================================= */

	// Get features for non-empty classes
	void getFeatures ()
	{

		minRes = 999.0;
		features_all_classes.clear();

		// Collect other features for one class
		int start_class = 0;
		int end_class = myopt.mymodel.nr_classes;
		if (only_use_this_class > 0)
		{
			start_class = only_use_this_class-1;
			end_class = only_use_this_class;

		}
		for (int iclass = start_class; iclass < end_class; iclass++)
		{
			if (debug>0) std::cerr << " dealing with class: " << iclass+1 << std::endl;
			class_features features_this_class;

			// Get class distribution and excluding empty classes
			features_this_class.class_distribution = myopt.mymodel.pdf_class[iclass];
			if (features_this_class.class_distribution > 0)
			{

				features_this_class.name = myopt.mymodel.ref_names[iclass];
				features_this_class.img() = myopt.mymodel.Iref[iclass];

				// Get class indexes from reference image name
				// TODO: why do we need this? remove?
				std::string root;
				features_this_class.name.decompose(features_this_class.class_index, root);

				// Get number of particles in the class from data.star file
				features_this_class.particle_nr = features_this_class.class_distribution * myopt.mydata.numberOfParticles(0);
				if (debug > 0) std::cerr << " features_this_class.particle_nr= " << features_this_class.particle_nr << std::endl;

				// Get selection label (if training data)
				if (MD_select.numberOfObjects() > 0)
				{
					MD_select.getValue(EMDL_SELECTED, features_this_class.is_selected, iclass);
				}
				else
				{
					features_this_class.is_selected = 1;
				}
				if (debug > 0) std::cerr << " features_this_class.is_selected= " << features_this_class.is_selected << std::endl;

				// Get estimated resolution (regardless of whether it is already in model_classes table or not)
				if (myopt.mymodel.estimated_resolution[iclass] > 0.)
				{
					features_this_class.estimated_resolution = myopt.mymodel.estimated_resolution[iclass];
				}
				else
				{
					// TODO: this still relies on mlmodel!!!
					findResolution(features_this_class);
				}
				if (debug > 0) std::cerr << " features_this_class.estimated_resolution= " << features_this_class.estimated_resolution << std::endl;

				// Calculate particle number-weighted resolution
				features_this_class.weighted_resolution = (1. / (features_this_class.estimated_resolution*features_this_class.estimated_resolution)) / log(features_this_class.particle_nr);
				if (debug > 0) std::cerr << " features_this_class.weighted_resolution= " << features_this_class.weighted_resolution << std::endl;

				// Calculate moments for class average image
				RFLOAT image_radius = XSIZE(features_this_class.img())/2.;
				circular_mask_radius = std::min(image_radius, myopt.particle_diameter / (myopt.mymodel.pixel_size*2));
				// Determining radius to use
				if (radius_ratio > 0 && radius <= 0) radius = radius_ratio * circular_mask_radius;
				if (radius > 0)
				{
					meanMomentsCalculator(features_this_class.img(), radius, circular_mask_radius, features_this_class.ring_moments);
					meanMomentsCalculator(features_this_class.img(), 0, radius, features_this_class.inner_circle_moments);
				}
				if (debug > 0) std::cerr << " done with moments" << std::endl;

				// Find job-wise best resolution among selected (red) classes in preparation for class score calculation called in the write_output function
				if (features_this_class.is_selected == 1 && features_this_class.estimated_resolution < minRes)
				{
					minRes = features_this_class.estimated_resolution;
				}

				if (!do_skip_angular_errors)
				{
					// Get class accuracy rotation and translation from model.star if present
					features_this_class.accuracy_rotation = myopt.mymodel.acc_rot[iclass];
					features_this_class.accuracy_translation = myopt.mymodel.acc_trans[iclass];
					if (debug>0) std::cerr << " myopt.mymodel.acc_rot[iclass]= " << myopt.mymodel.acc_rot[iclass] << " myopt.mymodel.acc_trans[iclass]= " << myopt.mymodel.acc_trans[iclass] << std::endl;
					if (features_this_class.accuracy_rotation > 99. || features_this_class.accuracy_translation > 99.)
					{
						calculateExpectedAngularErrors(iclass, features_this_class);
					}
					if (debug > 0) std::cerr << " done with angular errors" << std::endl;
				}

				if (!do_skip_protein_vs_solvent)
				{
					// Calculate protein and solvent region moments
					proteinVsSolventFeatures(features_this_class);
					if (debug > 0) std::cerr << " done with pvs" << std::endl;
				}

				if (!do_skip_LBP)
				{
					// Calculate whole image LBP and protein and solvent area LBP
					calculatePvsLBP(features_this_class.img(), features_this_class.lbp, features_this_class.lbp_p, features_this_class.lbp_s);
					if (debug > 0) std::cerr << " done with lbp" << std::endl;
				}

				features_all_classes.push_back(features_this_class);
			}

		} // end iterating all classes

	}

	/** ================================================== Writing Output ========================================================= */


	// TODO: Liyi: make a read
	void readClassFeatures()
	{
		if (fn_cf == "") return;

		MetaDataTable MD_class_features;
		MD_class_features.read(fn_cf);
		features_all_classes.resize(MD_class_features.numberOfObjects());
		for (int i=0; i<features_all_classes.size();i++)
		{

			MD_class_features.getValue(EMDL_MLMODEL_REF_IMAGE, features_all_classes[i].name);
			MD_class_features.getValue(EMDL_CLASS_FEAT_CLASS_INDEX, features_all_classes[i].class_index);
			MD_class_features.getValue(EMDL_CLASS_FEAT_IS_SELECTED,features_all_classes[i].is_selected);
			MD_class_features.getValue(EMDL_MLMODEL_PDF_CLASS, features_all_classes[i].class_distribution);
			MD_class_features.getValue(EMDL_MLMODEL_ACCURACY_ROT, features_all_classes[i].accuracy_rotation);
			MD_class_features.getValue(EMDL_MLMODEL_ACCURACY_TRANS, features_all_classes[i].accuracy_translation);
			MD_class_features.getValue(EMDL_MLMODEL_ESTIM_RESOL_REF, features_all_classes[i].estimated_resolution);
			MD_class_features.getValue(EMDL_CLASS_FEAT_WEIGHTED_RESOLUTION, features_all_classes[i].weighted_resolution);
			MD_class_features.getValue(EMDL_CLASS_FEAT_PARTICLE_NR, features_all_classes[i].particle_nr);

			// Job-wise features
			MD_class_features.getValue(EMDL_MLMODEL_PIXEL_SIZE, myopt.mymodel.pixel_size);
			MD_class_features.getValue(EMDL_MLMODEL_ORIGINAL_SIZE, myopt.mymodel.ori_size);//??
			MD_class_features.getValue(EMDL_MLMODEL_NR_CLASSES, myopt.mymodel.nr_classes);
			MD_class_features.getValue(EMDL_MLMODEL_SIGMA_OFFSET, myopt.mymodel.sigma2_offset);
			MD_class_features.getValue(EMDL_MLMODEL_AVE_PMAX, myopt.mymodel.ave_Pmax);
			MD_class_features.getValue(EMDL_OPTIMISER_DO_CORRECT_CTF, myopt.do_ctf_correction);
			MD_class_features.getValue(EMDL_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, myopt.intact_ctf_first_peak);
			MD_class_features.getValue(EMDL_OPTIMISER_ITERATION_NO, myopt.iter); //??
			MD_class_features.getValue(EMDL_OPTIMISER_PARTICLE_DIAMETER, myopt.particle_diameter);
			MD_class_features.getValue(EMDL_OPTIMISER_HIGHRES_LIMIT_EXP, myopt.strict_highres_exp);//??
			MD_class_features.getValue(EMDL_CLASS_FEAT_JOB_SCORE, job_score);

			// Class score
			MD_class_features.getValue(EMDL_CLASS_FEAT_CLASS_SCORE, features_all_classes[i].class_score);

			// Moments for the ring, inner circle, and outer circle
			if (radius > 0)
			{
				MD_class_features.getValue(EMDL_CLASS_FEAT_RING_MEAN, features_all_classes[i].ring_moments.mean);
				MD_class_features.getValue(EMDL_CLASS_FEAT_RING_STDDEV, features_all_classes[i].ring_moments.stddev);
				MD_class_features.getValue(EMDL_CLASS_FEAT_RING_SKEW, features_all_classes[i].ring_moments.skew);
				MD_class_features.getValue(EMDL_CLASS_FEAT_RING_KURT, features_all_classes[i].ring_moments.kurt);
			}

			// Protein and solvent region moments
			MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_MEAN, features_all_classes[i].protein_moments.mean);
			MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_STDDEV, features_all_classes[i].protein_moments.stddev);
			MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_SKEW, features_all_classes[i].protein_moments.skew);
			MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_KURT, features_all_classes[i].protein_moments.kurt);
			MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_MEAN, features_all_classes[i].solvent_moments.mean);
			MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_STDDEV, features_all_classes[i].solvent_moments.stddev);
			MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_SKEW, features_all_classes[i].solvent_moments.skew);
			MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_KURT, features_all_classes[i].solvent_moments.kurt);
			MD_class_features.getValue(EMDL_CLASS_FEAT_SCATTERED_SIGNAL, features_all_classes[i].scattered_signal);
			MD_class_features.getValue(EMDL_CLASS_FEAT_EDGE_SIGNAL, features_all_classes[i].edge_signal);

		}

	}


	// Generate star file: write feature value of each of the features for all classes in the job out in the format of a star file
	void writeClassFeatures()
	{
		std::cout << "writing star file..." << std::endl;
		MetaDataTable MD_class_features;
		MD_class_features.setName("class_features");
		for (int i=0; i<features_all_classes.size();i++)
		{
			// First calculate class score based on best resolution among red classes of the job, job score, estimated resolution, and selection label
			calculateClassScore(features_all_classes[i], minRes);

			MD_class_features.addObject();
			MD_class_features.setValue(EMDL_MLMODEL_REF_IMAGE, features_all_classes[i].name);
			MD_class_features.setValue(EMDL_CLASS_FEAT_CLASS_INDEX, features_all_classes[i].class_index);
			MD_class_features.setValue(EMDL_CLASS_FEAT_IS_SELECTED,features_all_classes[i].is_selected);
			MD_class_features.setValue(EMDL_MLMODEL_PDF_CLASS, features_all_classes[i].class_distribution);
			MD_class_features.setValue(EMDL_MLMODEL_ACCURACY_ROT, features_all_classes[i].accuracy_rotation);
			MD_class_features.setValue(EMDL_MLMODEL_ACCURACY_TRANS, features_all_classes[i].accuracy_translation);
			MD_class_features.setValue(EMDL_MLMODEL_ESTIM_RESOL_REF, features_all_classes[i].estimated_resolution);
			MD_class_features.setValue(EMDL_CLASS_FEAT_WEIGHTED_RESOLUTION, features_all_classes[i].weighted_resolution);
			MD_class_features.setValue(EMDL_CLASS_FEAT_PARTICLE_NR, features_all_classes[i].particle_nr);

			// Job-wise features
			MD_class_features.setValue(EMDL_MLMODEL_PIXEL_SIZE, myopt.mymodel.pixel_size);
			MD_class_features.setValue(EMDL_MLMODEL_ORIGINAL_SIZE, myopt.mymodel.ori_size);//??
			MD_class_features.setValue(EMDL_MLMODEL_NR_CLASSES, myopt.mymodel.nr_classes);
			MD_class_features.setValue(EMDL_MLMODEL_SIGMA_OFFSET, myopt.mymodel.sigma2_offset);
			MD_class_features.setValue(EMDL_MLMODEL_AVE_PMAX, myopt.mymodel.ave_Pmax);
			MD_class_features.setValue(EMDL_OPTIMISER_DO_CORRECT_CTF, myopt.do_ctf_correction);
			MD_class_features.setValue(EMDL_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, myopt.intact_ctf_first_peak);
			MD_class_features.setValue(EMDL_OPTIMISER_ITERATION_NO, myopt.iter); //??
			MD_class_features.setValue(EMDL_OPTIMISER_PARTICLE_DIAMETER, myopt.particle_diameter);
			MD_class_features.setValue(EMDL_OPTIMISER_HIGHRES_LIMIT_EXP, myopt.strict_highres_exp);//??
			MD_class_features.setValue(EMDL_CLASS_FEAT_JOB_SCORE, job_score);

			// Class score
			MD_class_features.setValue(EMDL_CLASS_FEAT_CLASS_SCORE, features_all_classes[i].class_score);

			// Moments for the ring, inner circle, and outer circle
			if (radius > 0)
			{
				MD_class_features.setValue(EMDL_CLASS_FEAT_RING_MEAN, features_all_classes[i].ring_moments.mean);
				MD_class_features.setValue(EMDL_CLASS_FEAT_RING_STDDEV, features_all_classes[i].ring_moments.stddev);
				MD_class_features.setValue(EMDL_CLASS_FEAT_RING_SKEW, features_all_classes[i].ring_moments.skew);
				MD_class_features.setValue(EMDL_CLASS_FEAT_RING_KURT, features_all_classes[i].ring_moments.kurt);
			}

			// Protein and solvent region moments
			MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_MEAN, features_all_classes[i].protein_moments.mean);
			MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_STDDEV, features_all_classes[i].protein_moments.stddev);
			MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_SKEW, features_all_classes[i].protein_moments.skew);
			MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_KURT, features_all_classes[i].protein_moments.kurt);
			MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_MEAN, features_all_classes[i].solvent_moments.mean);
			MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_STDDEV, features_all_classes[i].solvent_moments.stddev);
			MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_SKEW, features_all_classes[i].solvent_moments.skew);
			MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_KURT, features_all_classes[i].solvent_moments.kurt);
			MD_class_features.setValue(EMDL_CLASS_FEAT_SCATTERED_SIGNAL, features_all_classes[i].scattered_signal);
			MD_class_features.setValue(EMDL_CLASS_FEAT_EDGE_SIGNAL, features_all_classes[i].edge_signal);


		}
		MD_class_features.write(fn_out);
	}


};



int main(int argc, char *argv[]){

	liyi_class_features class_features;

	try
	{
		class_features.read(argc, argv);
		class_features.initialise();
		class_features.getFeatures();
		class_features.writeClassFeatures();
	}
	catch (RelionError XE)
	{
		std::cerr << XE;
		exit(1);
	}
	return 0;
}














