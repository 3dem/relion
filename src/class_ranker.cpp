 /***************************************************************************
 *
 * Author: "Liyi Dong and Sjors H.W. Scheres"
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

#include "src/class_ranker.h"


//
// Calculates n! (uses double arithmetic to avoid overflow)
//
double ZernikeMomentsExtractor::factorial(long n)
{
	if(n < 0)
		return(0.0) ;
	if(n == 0)
		return(1.0) ;
	else
		return(n * factorial(n-1)) ;
}

//
// Equations for Zernike moments from Prokop and Reeves, CVGIP: Graphical
//   Models and Image Processing, v 54 n 5, Sep 92, pp 438-460.
//

// Function to calculate the Zernike polynomial Rnl(r).
//  Note: 0 <= r <= 1
double ZernikeMomentsExtractor::zernikeR(int n, int l, double r)
{
	int m ;
	double sum = 0.0 ;

	if( ((n-l) % 2) != 0 )
	{
		std::cerr << "zernikeR(): improper values of n,l\n" ;
		return(0) ;
	}

	for(m = 0; m <= (n-l)/2; m++)
	{
		sum += (pow((double)-1.0,(double)m)) * ( factorial(n-m) ) /
				( factorial(m) * (factorial((n - 2*m + l) / 2)) *
					(factorial((n - 2*m - l) / 2)) ) *
					( pow((double)r, (double)(n - 2*m)) );
	}

	return(sum) ;
}

Complex ZernikeMomentsExtractor::zernikeZ(MultidimArray<double> img, int n, int l, double r_max)
{
  double rho ;		// radius of pixel from COM
  double theta ;    // angle of pixel
  Complex integral(0.,0.) ;

  FOR_ALL_ELEMENTS_IN_ARRAY2D(img)
  {

	  if(r_max > 0.0)
		  rho = sqrt((double)(i*i + j*j)) / r_max;
	  else
		  rho = 0.0;

	  if(rho <= 1.0)
	  {
		  theta = (i == 0 && j == 0) ? 0.0 :  atan2(i, j);
		  Complex aux(cos(l*theta), sin(l*theta));
		  integral += zernikeR(n,l,rho) * A2D_ELEM(img, i, j) * rho * conj(aux);
	  }

  }

  return(integral * (n+1)/PI) ;

}


std::vector<double> ZernikeMomentsExtractor::getZernikeMoments(MultidimArray<double> img, long z_order, double radius, bool verb)
{
	if (z_order > 20 || z_order < 0)
		REPORT_ERROR("BUG: zernike(): You choice of z_order is invalid; choose a value between 0 and 20");


	std::vector<double> zfeatures;

	// Normalise images to be intensity from [0,1]
	double minval, maxval, range;
	MultidimArray<int> mask;
	mask.resize(img);
	mask.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY2D(mask)
	{
		if ((double)(i*i+j*j) <= radius*radius)
			A2D_ELEM(mask, i, j) = 1;
		else
			A2D_ELEM(mask, i, j) = 0;
	}
	img.computeDoubleMinMax(minval, maxval, &mask);
	range = maxval -minval;
	if (range > 0.)
	{
		img = img - minval;
		img /= range;
	}
	else
	{
		long z_num_features = (long)( ((z_order + 4) *
					   (z_order + 1) - 2 *
					   (long)(((long)z_order + 1) / (long)2) ) / 4 ) ;
		zfeatures.resize(z_num_features, 0.);
		return zfeatures;
	}

	// Calculate Zernike moments
	for (int n = 0; n <= z_order; n++)
	{
		for (int l = 0; l <= n; l++)
		{
			if ((n-l) % 2 == 0)
			{
				zfeatures.push_back(abs(zernikeZ(img, n, l, radius))) ;
			}
		}
	}

	if (verb)
	{
		for (int i=0; i < zfeatures.size(); i++)
		{
			std::cerr << " i= " << i << " zfeatures[i]= " << zfeatures[i] << std::endl;
		}
	}

	return zfeatures;
}


double HaralickExtractor::Entropy(MultidimArray<double> arr)
{
	double result = 0.0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(arr)
	{
		result += DIRECT_MULTIDIM_ELEM(arr, n) * std::log(DIRECT_MULTIDIM_ELEM(arr, n) + HARALICK_EPS);
	}
	return -1 * result;
}

/*calculates probsum, probdiff, margprobx and y at once*/
void HaralickExtractor::fast_init()
{
	if (NZYXSIZE(matcooc) == 0) return;

	margprobx.initZeros(XSIZE(matcooc));
	margproby.initZeros(YSIZE(matcooc));
	probsum.initZeros(2*XSIZE(matcooc));
	probdiff.initZeros(XSIZE(matcooc));

	double local;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(matcooc)
	{
		local = DIRECT_A2D_ELEM(matcooc, i, j);
		DIRECT_A1D_ELEM(margprobx, j) += local;
		DIRECT_A1D_ELEM(margproby, i) += local;
		DIRECT_A1D_ELEM(probsum, i + j) += local;
		DIRECT_A1D_ELEM(probdiff, abs(i - j)) += local;

	}
	hx = Entropy(margprobx);
	hy = Entropy(margproby);
	margprobx.computeAvgStddev(meanx, stddevx);
	margproby.computeAvgStddev(meany, stddevy);
	//Everything set up
	initial = true;
}

/*0 => energy, 1 => entropy, 2=> inverse difference */
/*3 => correlation, 4=> info measure 1, 5 => info measure 2*/
std::vector<double> HaralickExtractor::cooc_feats()
{
	std::vector<double> ans(7, 0.0);
	double hxy1 = 0.0;
	double hxy2 = 0.0;
	double local, xy;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(matcooc)
	{
		local = DIRECT_A2D_ELEM(matcooc, i, j);
		ans[0] += local * local;
		ans[1] -= local * log(local + HARALICK_EPS);
		ans[2] += local * (1 / (1 + (i - j) * (i - j)));
		ans[3] += (i * j * local) - (meanx * meany);
		ans[6] += (i-meany)*(j-meanx)*local;
		xy = DIRECT_A1D_ELEM(margprobx, j) * DIRECT_A1D_ELEM(margproby, i);
		hxy1 -= local * log(xy + HARALICK_EPS);
		hxy2 -= xy * log(xy + HARALICK_EPS);
	}

	ans[3] = ans[3] / (stddevx * stddevy);
	ans[4] = (ans[1] - hxy1) / std::max(hx, hy);
	//std::cerr << " hxy1= " << hxy1 << " hxy2= " << hxy2 << " ans[1]= " << ans[1] << "  exp(-2 *(hxy2 - ans[1]))= " <<  exp(-2 *(hxy2 - ans[1])) << " arg= " << -2 *(hxy2 - ans[1]) << std::endl;
	double arg = -2. * (hxy1 - ans[1]);
	if (arg >= 0.)
		ans[5] = 0.;
	else if (arg < -50.)
		ans[5] = 1.;
	else
		ans[5] = sqrt(1 - exp(arg));
	return ans;
}

/*0 => contrast, 1 => diff entropy, 2 => diffvariance */
/*3 => sum average, 4 => sum entropy, 5 => sum variance */
std::vector<double> HaralickExtractor::margprobs_feats()
{
	std::vector<double> ans(6, 0.0);
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(probdiff)
	{
		ans[0] += i * i * DIRECT_A1D_ELEM(probdiff, i);
		ans[1] += -1 * DIRECT_A1D_ELEM(probdiff, i) * log(DIRECT_A1D_ELEM(probdiff, i) + HARALICK_EPS);
	}
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(probdiff)
	{
		ans[2] += (i - ans[1]) * (i - ans[1]) * DIRECT_A1D_ELEM(probdiff, i);
	}

	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(probsum)
	{
		ans[3] += i *  DIRECT_A1D_ELEM(probsum, i);
		ans[4] += -1 * DIRECT_A1D_ELEM(probsum, i) * log(DIRECT_A1D_ELEM(probsum, i) + HARALICK_EPS);
	}
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(probsum)
	{
		ans[5] += (i - ans[4]) * (i - ans[4]) * DIRECT_A1D_ELEM(probsum, i);
	}
	return ans;
}


MultidimArray<double> HaralickExtractor::fast_feats(bool verbose)
{
	MultidimArray<double> result;
	result.initZeros(13);
	if (NZYXSIZE(matcooc) ==0) return result;
	if (!initial) fast_init();
	std::vector<double> margfeats = margprobs_feats();
	std::vector<double> coocfeats = cooc_feats();
	for (int i = 0; i < 7; i++)
	{
		result(i) = coocfeats[i];
	}
	for (int i = 0; i < 6; i++)
	{
		result(7 + i) = margfeats[i];
	}
	return result;
}

MultidimArray<RFLOAT> HaralickExtractor::MatCooc(MultidimArray<int> img, int N,
		int deltax, int deltay, MultidimArray<int> *mask)
{
	int target, next;
	int newi, newj;
	MultidimArray<RFLOAT> ans;
	ans.initZeros(N + 1, N + 1);
	RFLOAT counts = 0.;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(img)
	{

		if (mask == NULL || DIRECT_A2D_ELEM(*mask, i, j) > 0)
		{
			newi = i + deltay;
			newj = j + deltax;

			// Stay inside the image
			if (newi < YSIZE(img) && newj < XSIZE(img) && newj >= 0 && newi >= 0 &&
					(mask == NULL || DIRECT_A2D_ELEM(*mask, i, j) > 0) 	)
			{
				target = DIRECT_A2D_ELEM(img, i, j);
				next  = DIRECT_A2D_ELEM(img, newi, newj);
				DIRECT_A2D_ELEM(ans, target, next) += 1.0;
				// this is not in original code from Abello, but that's how I understand it should be done...
				DIRECT_A2D_ELEM(ans, next, target) += 1.0;
				counts += 2.;
			}
		}
	}

	ans /= counts;

	return ans;
}

std::vector<double> HaralickExtractor::getHaralickFeatures(MultidimArray<RFLOAT> img,
		MultidimArray<int> *mask, bool verbose)
{
	std::vector<double> ans;
	ans.resize(13, 0.);

	MultidimArray<double> avg;
	avg.initZeros(13);

	// Check there are non-zero elements in the mask
	if (mask != NULL && (*mask).sum() == 0) return ans;

	// Convert greyscale image to integer image with much fewer (32) grey-scale values
	MultidimArray<int> imgint;
	imgint.resize(img);
	double minval, maxval, range;
	img.computeDoubleMinMax(minval, maxval, mask);
	range = maxval -minval;
	if (range > 0.)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img)
		{
			if (mask == NULL || DIRECT_MULTIDIM_ELEM(*mask, n) > 0)
			{
				DIRECT_MULTIDIM_ELEM(imgint, n) = floor( ((DIRECT_MULTIDIM_ELEM(img, n) - minval) * 31.0) / range);
			}
		}

		// Average over Haralick features in horizontal, vertical, and 2 diagonal directions
		/*
		x o
		O o
		o o
		*/
		matcooc = MatCooc(imgint, 31, 0, 1, mask);
		fast_init(); //initialize internal variables
		avg += fast_feats();
		matcooc = MatCooc(imgint, 31, 1, 0, mask);
		fast_init(); //initialize internal variables
		avg += fast_feats();
		matcooc = MatCooc(imgint, 31, 1, 1, mask);
		fast_init(); //initialize internal variables
		avg += fast_feats();
		matcooc = MatCooc(imgint, 31, 1, -1, mask);
		fast_init(); //initialize internal variables
		avg += fast_feats();

		for (int i = 0; i < ans.size(); i++) ans[i] = avg(i) / 4.;

		if (verbose)
		{
			std::cout << " - Energy: " << ans[0] << std::endl;
			std::cout << " - Entropy: " << ans[1] << std::endl;
			std::cout << " - Inverse Difference Moment: " << ans[2] << std::endl;
			std::cout << " - Correlation: " << ans[3] << std::endl;
			std::cout << " - Info Measure of Correlation 1: " << ans[4] << std::endl;
			std::cout << " - Info Measure of Correlation 2: " << ans[5] << std::endl;
			std::cout << " - Sum of squares variance: " << ans[6] << std::endl;
			std::cout << " - Contrast: " << ans[7] << std::endl;
			std::cout << " - Difference Entropy: " << ans[8] << std::endl;
			std::cout << " - Difference Variance: " << ans[9] << std::endl;
			std::cout << " - Sum Average: " << ans[10] << std::endl;
			std::cout << " - Sum Entropy: " << ans[11] << std::endl;
			std::cout << " - Sum Variance: " << ans[12] << std::endl;
		}
	}

	return ans;

}




void ClassRanker::read(int argc, char **argv, int rank)
{
	parser.setCommandLine(argc, argv);

	// TODO: optional input files, eg. job score file
	int general_section = parser.addSection("General options");
	fn_optimiser = parser.getOption("--opt", "Input optimiser.star file", "");
	fn_out = parser.getOption("--o", "Directory name for output files", "./");
	fn_ext = parser.getOption("--ext", "Extension for root filenames of output optimiser.star and model.star", "ranked");
	do_select = parser.checkOption("--auto_select", "Perform auto-selection of particles based on below thresholds for the score");
	select_min_score = textToFloat(parser.getOption("--min_score", "Minimum selected score to be included in class selection", "0.5"));
	select_max_score = textToFloat(parser.getOption("--max_score", "Maximum selected score to be included in class selection", "999."));
	do_relative_threshold = parser.checkOption("--relative_thresholds", "If true, interpret the above min and max_scores as fractions of the maximum score of all predicted classes in the input");

	int part_section = parser.addSection("Network training options (only used in development!)");
	do_ranking  = !parser.checkOption("--train", "Only write output files for training purposes (don't rank classes)");
	fn_select = parser.getOption("--select", "Input class_averages.star from the Selection job or backup_selection.star", "backup_selection.star");
	fn_job_score = parser.getOption("--fn_score", "Input job score file", "job_score.txt");
	fn_cf = parser.getOption("--cf_file", "Input class feature star file", "");
	only_use_this_class = textToInteger(parser.getOption("--only_class_nr", "Class number of the class of interest", "-1"));
	do_skip_angular_errors = parser.checkOption("--skip_angular_errors", "Skip angular error calculation");
	do_granularity_features  = parser.checkOption("--do_granularity_features", "Calculate granularity features");
	do_save_masks = parser.checkOption("--save_masks", "Save the automatically generated 2D solvent masks for all references");
	do_save_mask_c = parser.checkOption("--save_mask_c", "Write out images of protein mask circumferences.");
	fn_mask_dir = parser.getOption("--mask_folder_name", "Folder name for saving all masks.", "protein_solvent_masks");

	int subimg_section = parser.addSection("Extract subimage for deep convolutional neural network analysis");
	do_subimages = parser.checkOption("--extract_subimages", "Extract subimages within protein mask for each class");
	nr_subimages = textToInteger(parser.getOption("--nr_subimages", "Number of subimage to extract (randonly)", "25"));
	subimage_boxsize = textToInteger(parser.getOption("--subimage_boxsize", "Boxsize (in pixels) for subimages", "16"));

	int expert_section = parser.addSection("Expert options");
	radius_ratio = textToFloat(parser.getOption("--radius_ratio", "Ratio of inner radius of the interested ring area in proportion to the current circular mask radius", "0.95"));
	radius = textToFloat(parser.getOption("--radius", "Inner radius of the interested ring area to the current circular mask radius", "-1"));
	lowpass = textToFloat(parser.getOption("--lowpass", "Image lowpass filter threshold for generating binary masks.", "25"));
	binary_threshold = textToFloat(parser.getOption("--binary_threshold", "Threshold for generating binary masks.", "0."));
    debug = textToInteger(parser.getOption("--debug", "Debug level", "0"));
	verb = textToInteger(parser.getOption("--verb", "Verbosity level", "1"));
	fn_features = parser.getOption("--fn_features", "Filename for output features star file", "features.star");
	fn_sel_parts = parser.getOption("--fn_sel_parts", "Filename for output star file with selected particles", "particles.star");
	fn_sel_classavgs = parser.getOption("--fn_sel_classavgs", "Filename for output star file with selected class averages", "class_averages.star");
	fn_root = parser.getOption("--fn_root", "rootname for output model.star and optimiser.star files", "run");
	fn_subimages = parser.getOption("--output_subimages", "Name of the mrcs file for saving output subimages.", "subimages.mrcs");
	fn_subimage_star = parser.getOption("--output_subimage_star", "Name of the subimage star file.", "subimages.star");

	fn_torch_model = parser.getOption("--fn_torch_model", "Filename for the serialized Torch model.", ""); // Default should be compile-time defined

	// Check for errors in the command-line option
	if (radius > 0 && radius < 1)
		REPORT_ERROR("Incorrect value(s) for inner_radius / outer_radius. Try inner_radius_ratio / outer_radius_ratio.");
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	// Initialise verb for non-parallel execution
	verb = 1;

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}

void ClassRanker::usage()
{
	parser.writeUsage(std::cout);
}

void ClassRanker::initialise()
{

	if (verb > 0) std::cout << " Initialising... " << std::endl;

	// Make sure output rootname ends with a "/" and make directory
	if (fn_out[fn_out.length()-1] != '/') fn_out += "/";
	if (!exists(fn_out))
	{
		int res = system(("mkdir -p " + fn_out).c_str());
	}

	if (do_skip_angular_errors)
	{
		if (fn_cf == "") REPORT_ERROR("ERROR: you need to provide a class feature input file if you wish to skip some calculations!");
	}

	// Read in the MD_optimiser table from the STAR file, get model.star and data.star
	if (fn_optimiser != "")
	{

		MD_optimiser.read(fn_optimiser, "optimiser_general");
		MD_optimiser.getValue(EMDL_OPTIMISER_MODEL_STARFILE, fn_model);
		MD_optimiser.getValue(EMDL_OPTIMISER_DATA_STARFILE, fn_data);

		MD_optimiser.getValue(EMDL_OPTIMISER_DO_CORRECT_CTF, do_ctf_correction);
		MD_optimiser.getValue(EMDL_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, intact_ctf_first_peak);
		MD_optimiser.getValue(EMDL_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED, ctf_phase_flipped);
		MD_optimiser.getValue(EMDL_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES, only_flip_phases);
		MD_optimiser.getValue(EMDL_OPTIMISER_PARTICLE_DIAMETER, particle_diameter);

		if (!do_ctf_correction)
		{
			std::cerr << " Skipping this job because it hasn't been done with CTF correction ..." << std::endl;
			exit(1);
		}

		//Sjors 06022020: go back to just reading MD_optimiser for speed
		mymodel.read(fn_model, true); // true means: read only one group!
		if (debug>0) std::cerr << "Done with reading model.star ..." << std::endl;

		//myopt.read(fn_optimiser); // true means skip_groups_and_pdf_direction from mlmodel; only read 1000 particles...
		//if (debug>0) std::cerr << "Done with reading optimiser ..." << std::endl;

		// Collect features for all specified classes
		start_class = 0;
		end_class = mymodel.nr_classes;
		if (only_use_this_class > 0)
		{
			start_class = only_use_this_class-1;
			end_class = only_use_this_class;
		}

		// Check if any angular accuracies are missing
		bool haveAllAccuracies = true;

		// Pre-calculate Fourier Transforms of all classes, but only if angular accuracies need to be computed
		if (!do_skip_angular_errors)
		{
			for (int iclass = start_class; iclass < end_class; iclass++)
			{
				// Only consider features with non-zero class distributions
				if (mymodel.pdf_class[iclass] > 0)
				{

					if (mymodel.acc_rot[iclass] > 99. || mymodel.acc_rot[iclass] > 99.)
					{
						haveAllAccuracies = false;
						break;
					}
				}
			}
			if (!haveAllAccuracies)
			{
				// Set FTs of the references
				mymodel.setFourierTransformMaps(false);
			}
		}

		if (intact_ctf_first_peak || do_ranking || (!do_skip_angular_errors && !haveAllAccuracies))
		{
			// Read in particles (otherwise wait until haveAllAccuracies or performRanking, as Liyi sometimes doesn't need mydata)
			mydata.read(fn_data, true, true); // true true means: ignore particle_name and group name!
			total_nr_particles = mydata.numberOfParticles(0);
			if (debug>0) std::cerr << "Done with reading data.star ..." << std::endl;

		}
		else
		{
			MetaDataTable MDtmp;
			total_nr_particles = MDtmp.read(fn_data, "", NULL, "", true); // true means do_only_count
		}

		if (intact_ctf_first_peak)
		{
			if (verb > 0) std::cout << " Doing first peak CTF correction ..." << std::endl;

			// Calculate avg. defocus
			RFLOAT def_avg = 0, def_u, def_v;
			for (long int part_id = 0; part_id < mydata.MDimg.numberOfObjects(); part_id++)
			{
				mydata.MDimg.getValue(EMDL_CTF_DEFOCUSU, def_u, part_id);
				mydata.MDimg.getValue(EMDL_CTF_DEFOCUSV, def_v, part_id);
				def_avg += def_u + def_v;
			}
			def_avg /= (2. * mydata.MDimg.numberOfObjects());

			// some people may have used a too small , or even zero amplitude contrast, lets forbid that...
			//q0 = XMIPP_MAX(q0, 0.07);
			CTF avgctf;
			avgctf.setValuesByGroup(&mydata.obsModel, 0, def_u, def_v, 0.);

			// Loop over all classes in myopt.mymodel.Iref
			for (int iref =0; iref < mymodel.Iref.size(); iref++)
			{
				correctCtfUntilFirstPeak(mymodel.Iref[iref], avgctf);
			}
		}

	}

	if (!do_select)
	{
		// These options are for generating training data

		// Open selected_class table
		if (fn_select != "")
		{
			if (verb > 0) std::cout << " Reading in selection file: " << fn_select << std::endl;
			MD_select.read(fn_select);
		}

		// Read in class features from a previous run if fn_cf is provided
		if (fn_cf != "")
		{

			if (verb > 0) std::cout << " Reading in previously calculated feature star file: " << fn_cf << std::endl;
			readFeatures();
		}

		// Get job score
		if (fn_job_score != "")
		{
			std::ifstream in(fn_job_score, std::ios_base::in);
			if (in.fail()) REPORT_ERROR( (std::string) "ERROR: File " + fn_job_score + " does not exists" );
			std::string line;
			getline(in, line, '\n');
			job_score = textToFloat(line);
			in.close();
			if (verb > 0) std::cout << " Read in job_score of: " << job_score << std::endl;
		}
		else
		{
			job_score = -1.;
		}
	}

	if (radius_ratio > 0 && radius > 0)
	{
		std::cout << "WARNING: Should not provide radius ratio and radius at the same time. Ignoring the radius ratio..." << std::endl;
	}
}


void ClassRanker::run()
{

	initialise();
	getFeatures();
	if (do_ranking) performRanking();
	writeFeatures();

	if (verb > 0) std::cout << "Done!" << std::endl;

	return;
}

int ClassRanker::getClassIndex(FileName &name)
{
	// Get class indexes from reference image name
	long int result;
	std::string root;
	name.decompose(result, root);
	return result;

}

// SHWS 15072020: extract subimages within protein mask
MultidimArray<RFLOAT> ClassRanker::getSubimages(MultidimArray<RFLOAT> &img, int boxsize, int nr_images, MultidimArray<int> *mask)
{

	// TODO: make the if statement for only making subimages for classes with non-zero protein area inside this function

	MultidimArray<RFLOAT> subimages(nr_images, 1, boxsize, boxsize);
	subimages.setXmippOrigin();
	int imgsize = XSIZE(img);

	for (int my_image = 0; my_image < nr_images; my_image++)
	{
		int xpos, ypos;

		bool got_one = false;
		while (!got_one)
		{
			// uniformly sample x and y and see if it is within the mask
			xpos = FLOOR(rnd_unif(FIRST_XMIPP_INDEX(imgsize), LAST_XMIPP_INDEX(imgsize)));
			ypos = FLOOR(rnd_unif(FIRST_XMIPP_INDEX(imgsize), LAST_XMIPP_INDEX(imgsize)));
			if ((mask == NULL) || (A2D_ELEM(*mask, ypos, xpos) > 0) )
			{
				got_one = true;
			}
		}

		// Extract subimage here
		int x0 = xpos + FIRST_XMIPP_INDEX(boxsize);
		int xF = xpos + LAST_XMIPP_INDEX(boxsize);
		int y0 = ypos + FIRST_XMIPP_INDEX(boxsize);
		int yF = ypos + LAST_XMIPP_INDEX(boxsize);

		MultidimArray<RFLOAT> subimg;
		img.window(subimg, y0, x0, yF, xF, 0.);
		subimages.setImage(my_image, subimg);

	}
	return subimages;
}

// Calculate moments for one class
moments ClassRanker::calculateMoments(MultidimArray<RFLOAT> &img, RFLOAT inner_radius, RFLOAT outer_radius, MultidimArray<int> *mask)
{
	moments result;

	int pix_num = 0;
	RFLOAT inner_radius_square, outer_radius_square;

	if (mask != NULL  && (*mask).sum() == 0)
	{
		result.sum = result.mean = result.stddev = result.skew = result.kurt = 0.;
		return result;
	}

	inner_radius_square = pow((inner_radius), 2);
	outer_radius_square = pow((outer_radius), 2);

	// Calculating image means
	RFLOAT sum = 0;
	for (long int i= STARTINGY(img) ; i <= FINISHINGY(img); i++)
	{
		for (long int j= STARTINGX(img) ; j <= FINISHINGX(img); j++)
		{
			if (i*i+j*j >= inner_radius_square && i*i+j*j <= outer_radius_square && ((mask == NULL) || (A2D_ELEM(*mask, i, j) > 0) ) )
			{
				pix_num++;
				sum += A2D_ELEM(img, i, j);
			}
		}
	}
	result.sum = sum;
	result.mean = sum / pix_num;

	// Calculating other image moments
	RFLOAT square_sum = 0, cube_sum = 0, quad_sum = 0;
	for (long int i= STARTINGY(img) ; i <= FINISHINGY(img); i++)
	{
		for (long int j= STARTINGX(img) ; j <= FINISHINGX(img); j++)
		{
			if (i*i+j*j >= inner_radius_square && i*i+j*j <= outer_radius_square && ((mask == NULL) || (A2D_ELEM(*mask, i, j) > 0) ) )
			{
				square_sum += pow((A2D_ELEM(img, i, j)-result.mean), 2);
				cube_sum += pow((A2D_ELEM(img, i, j)- result.mean), 3);
				quad_sum += pow((A2D_ELEM(img, i, j)- result.mean), 4);
			}
		}
	}

	result.stddev = sqrt(square_sum / pix_num);
	if (square_sum == 0)
	{
		result.skew = 0;
		result.kurt = 0;
	} else
	{
		result.skew = (cube_sum * sqrt (pix_num))/(pow(square_sum, 1.5));
		result.kurt = (quad_sum * pix_num)/pow(square_sum, 2);
	}

	return result;
}


void ClassRanker::calculatePvsLBP(MultidimArray<RFLOAT> I, MultidimArray<int> &p_mask, MultidimArray<int> &s_mask, classFeatures &cf)
{

	unsigned char code;
	double center;
	double lbp_hist[256] = {}, lbp_hist_p[256] = {}, lbp_hist_s[256] = {};

	// Make the map from 256 possible original rotation-variant LBP values to 36 rotation-invariant LBP values
	std::vector<double> min_idxs, min_idxs_sort;
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
		cf.lbp.push_back(lbp_hist[idx]);
		cf.lbp_p.push_back(lbp_hist_p[idx]);
		cf.lbp_s.push_back(lbp_hist_s[idx]);
	}
}

std::vector<RFLOAT> ClassRanker::calculateGranulo(const MultidimArray<double> &I)
{

	std::vector<RFLOAT> result;

	Image<double> G;
    G().resize(I);
    double m, M;
    I.computeDoubleMinMax(m, M);

    if (XSIZE(I) < 15 || YSIZE(I) < 15)
    {
        std::cerr << "ERROR: Input image must be at least 15x15px"
                  << "to extract granulo features!\n";
        exit(1);
    }

    for (int N = 1; N < 7; N++)
    {
        // creating circular structuring element
        int size = N*2 + 1;
        bool struct_elem[size][size];
        for (int y = 0; y < size; y++)
        {
            for (int x = 0; x < size; x++)
                struct_elem[x][y] = ((x-N)*(x-N) + (y-N)*(y-N)) <= N*N;
        }

        // morphological erosion
        double sum = 0.0;
        for (int y = 0; y < YSIZE(I); y++)
        {
            for (int x = 0; x < XSIZE(I); x++)
            {
                double struct_min = M;

                for (int yy = y-N; yy <= y+N; yy++)
                {
                    if (yy < 0 || yy >= YSIZE(I)) continue;

                    for (int xx = x-N; xx <= x+N; xx++)
                    {
                        if (xx < 0 || xx >= XSIZE(I)) continue;

                        if (struct_elem[xx+N-x][yy+N-y] &&
                            DIRECT_A2D_ELEM(I, yy, xx) < struct_min)
                            struct_min = DIRECT_A2D_ELEM(I, yy, xx);
                    }
                }
                DIRECT_A2D_ELEM(G(), y, x) = struct_min;
            }
        }

        // morphological dilation (dilation after erosion = opening)
        for (int y = 0; y < YSIZE(I); y++)
        {
            for (int x = 0; x < XSIZE(I); x++)
            {
                double struct_max = m;

                for (int yy = y-N; yy <= y+N; yy++)
                {
                    if (yy < 0 || yy >= YSIZE(I)) continue;

                    for (int xx = x-N; xx <= x+N; xx++)
                    {
                        if (xx < 0 || xx >= XSIZE(I)) continue;

                        if (struct_elem[xx+N-x][yy+N-y] &&
                            DIRECT_A2D_ELEM(G(), yy, xx) > struct_max)
                            struct_max = DIRECT_A2D_ELEM(G(), yy, xx);
                    }
                }

                sum += struct_max;
            }
        }

        result.push_back(sum);
    }

    return result;
}


/** =================================================== Getting Other Class Metadata ======================================================  */

// Find estimated resolution from MD_model_classes_ table for one class
RFLOAT ClassRanker::findResolution(classFeatures &cf)
{

	RFLOAT result;

	std::string class_name = "model_class_" + integerToString(cf.class_index);
	MetaDataTable MD_model_this_class;
	MD_model_this_class.read(fn_model, class_name);
	bool isGood = true;

	// Test for numberOfObjects() function
	for (int i=1; i<MD_model_this_class.numberOfObjects();i++)
	{
		RFLOAT ssnrMap;
		MD_model_this_class.getValue(EMDL_MLMODEL_DATA_VS_PRIOR_REF, ssnrMap, i);
		if (ssnrMap < 1)
		{
			MD_model_this_class.getValue(EMDL_RESOLUTION_ANGSTROM, result, i-1);
			isGood = false;
			break;
		}
	}
	if (isGood)
	{
		result = mymodel.pixel_size*2;
	}

	return result;
}


// Calculate accuracy rotation and translation for each of the non-empty classes
void ClassRanker::calculateExpectedAngularErrors(int iclass, classFeatures &cf)
{
	// Set current_image_size to the coarse_size to calculate expected angular errors
	int current_image_size = mymodel.current_size;

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
			// SHWS 6Feb2020: just work with noise spectrum from group 0 to save time!
			int group_id = 0; // mydata.getGroupId(part_id, 0);
			RFLOAT my_pixel_size = mydata.getImagePixelSize(part_id, 0);
			const int optics_group = mydata.getOpticsGroup(part_id, 0);
			int my_image_size = (mydata.obsModel.hasBoxSizes) ? mydata.getOpticsImageSize(optics_group) : mymodel.ori_size;
			bool ctf_premultiplied = mydata.obsModel.getCtfPremultiplied(optics_group);

			MultidimArray<RFLOAT> Fctf;
			// Get CTF for this particle (do_ctf_correction is always true for this program)!
			if (do_ctf_correction)
			{
				Fctf.resize(current_image_size, current_image_size/ 2 + 1);

				// Get parameters that change per-particle from the exp_metadata
				CTF ctf;
				RFLOAT def_u, def_v, def_angle, voltage, cs, q0;
				mydata.MDimg.getValue(EMDL_CTF_DEFOCUSU, def_u, part_id);                 //??
				mydata.MDimg.getValue(EMDL_CTF_DEFOCUSV, def_v, part_id);
				mydata.MDimg.getValue(EMDL_CTF_DEFOCUS_ANGLE, def_angle, part_id);
				ctf.setValuesByGroup(&mydata.obsModel, optics_group, def_u, def_v, def_angle);
				ctf.getFftwImage(Fctf, my_image_size, my_image_size, mymodel.pixel_size,
						ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true, false);
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
					mydata.MDimg.getValue(EMDL_ORIENT_ROT, rot1, part_id);
					mydata.MDimg.getValue(EMDL_ORIENT_TILT, tilt1, part_id);
					mydata.MDimg.getValue(EMDL_ORIENT_PSI, psi1, part_id);

					F1.initZeros(current_image_size, current_image_size/ 2 + 1);

					// Get the FT of the first image
					Euler_angles2matrix(rot1, tilt1, psi1, A1, false);
					// Older versions of RELION dont have the metadata to do this, so disable then
					if (mydata.obsModel.hasBoxSizes)
					{
						A1 = mydata.obsModel.applyAnisoMag(A1, optics_group);
						A1 = mydata.obsModel.applyScaleDifference(A1, optics_group, mymodel.ori_size, mymodel.pixel_size);
					}
					(mymodel.PPref[iclass]).get2DFourierTransform(F1, A1);    //?
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
						if (mymodel.ref_dim == 3)
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
						if (mymodel.data_dim == 3)
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
					if (mymodel.data_dim == 2)
					  F2.initZeros(current_image_size, current_image_size/ 2 + 1);
					else
					  F2.initZeros(current_image_size, current_image_size, current_image_size/ 2 + 1);

					if (imode == 0)
					{
						// Get new rotated version of reference
						Euler_angles2matrix(rot2, tilt2, psi2, A2, false);
						// Older versions of RELION dont have the metadata to do this, so disable then
						if (mydata.obsModel.hasBoxSizes)
						{
							A2 = mydata.obsModel.applyAnisoMag(A2, optics_group);
							A2 = mydata.obsModel.applyScaleDifference(A2, optics_group, mymodel.ori_size, mymodel.pixel_size);
						}
						(mymodel.PPref[iclass]).get2DFourierTransform(F2, A2);    //??
					}
					else
					{
						// Get shifted version
						shiftImageInFourierTransform(F1, F2, (RFLOAT) mymodel.ori_size, -xshift, -yshift, -zshift);
					}

					// Apply CTF to F1 and F2 if necessary
					if (do_ctf_correction)
					{
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
						{
							DIRECT_MULTIDIM_ELEM(F1, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
							DIRECT_MULTIDIM_ELEM(F2, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
						}

						if (mydata.hasCtfPremultiplied())
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
						if (ires < mymodel.ori_size / 2 + 1 && !(jp==0 && ip < 0))
						{
							my_snr += norm(DIRECT_A2D_ELEM(F1, i ,j) - DIRECT_A2D_ELEM(F2, i, j)) / (2 * mymodel.sigma2_noise[group_id](ires) );
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
RFLOAT ClassRanker::getClassScoreFromJobScore(classFeatures &cf, RFLOAT minRes)
{
	RFLOAT weight, result;

	if (job_score != -1)
	{
		weight = pow((minRes/cf.estimated_resolution), 2);
		if (debug > 0) std::cout << " job_score: " << job_score << " weight: " << weight << std::endl;
		switch (cf.is_selected)
		{
		case 1:
			result = (0.75+weight*0.25)*job_score;
			break;
		case 2:
			result = (0.25+weight*0.25)*job_score;
			break;
		case 3:
			result = (0+weight*0.25)*job_score;
			break;
		case 4:
			result = (0+weight*0.25)*job_score;
			if (verb > 0) std::cout << "Class " << cf.name << " is labelled cyan!" << std::endl;
			break;
		case 5:
			result = (0.5+weight*0.25)*job_score;
			break;
		case 6:
			result = 0.0;
//			if (verb > 0) std::cout << "Class " << cf.name << " is labelled yellow!" << std::endl;
			break;
		case 0:
			result = 0.0;
			break;
		default:
			if (verb > 0)
			{
				//std::cout << "Class " << cf.name << ": Illegal selection label " << cf.is_selected << std::endl;
				std::cout << "Class " << cf.class_index << ": Illegal selection label " << cf.is_selected << std::endl;
			}

		}
	}
	else
	{
		result = -1;
	}
	if (debug > 0) std::cout << "class score: " << result << std::endl;

	return result;

}


void ClassRanker::makeSolventMasks(classFeatures cf, MultidimArray<RFLOAT> img, MultidimArray<RFLOAT> &lpf, MultidimArray<int> &p_mask, MultidimArray<int> &s_mask,
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
	p_mask.initZeros(lpf);
	s_mask.initZeros(lpf);
	p_mask.setXmippOrigin();
	s_mask.setXmippOrigin();

	protein_area = 0;
	long circular_area = 0;

	binary_threshold = 0.05*cf.lowpass_filtered_img_stddev;

	// A hyper-parameter to adjust: definition of central area: 0.7 of radius (~ half of the area)
	lowPassFilterMap(lpf, lowpass, uniform_angpix);

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
											   // Find a new white pixel that was never visited before and use it to identify a new island
					if (A2D_ELEM(lpf, i, j) > binary_threshold)
					{
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
						} // finish finding all pixels for the new island

						//  Debug
//						std::cerr << "Class " << cf.class_index << " inside: " << inside << std::endl;

						if (inside > 0.4*3.14*0.49*circular_mask_radius*circular_mask_radius)
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
	scattered_signal = (RFLOAT(scattered_signal_nr))/ (RFLOAT(scattered_signal_nr + protein_area + 1));
	// Make the masks
	for (long int i= STARTINGY(lpf); i <= FINISHINGY(lpf); i++)
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
			A2D_ELEM(p_mask, i, j) = 1;
			A2D_ELEM(s_mask, i, j) = 0;
		}
	}
}


void ClassRanker::saveMasks(Image<RFLOAT> &img, MultidimArray<RFLOAT> &lpf, MultidimArray<int> &p_mask, MultidimArray<int> &s_mask,
							classFeatures &cf)
{

	// Create folders to save protein and solvent region masks
	char foldername[256];
	snprintf(foldername, 255, "%s", fn_mask_dir.c_str());
	FileName filtered_mask_folder = foldername;
	if (!exists(filtered_mask_folder))
	{
		std::string filtered_mask_fn_command = "mkdir -p " + fn_out+filtered_mask_folder;
		int result = system(filtered_mask_fn_command.c_str());
	}

	Image<RFLOAT> p_out, s_out, lpf_out;
	lpf_out() = lpf;
	p_out().resize(YSIZE(p_mask), XSIZE(p_mask));
	s_out().resize(YSIZE(p_mask), XSIZE(p_mask));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(p_mask)
	{
		DIRECT_MULTIDIM_ELEM(p_out(), n) = (double)DIRECT_MULTIDIM_ELEM(p_mask, n);
		DIRECT_MULTIDIM_ELEM(s_out(), n) = (double)DIRECT_MULTIDIM_ELEM(s_mask, n);
	}

	FileName fp_out, fs_out, flpf_out, fimg_out;
	fimg_out = fn_out+filtered_mask_folder+"/class"+integerToString(cf.class_index)+"_original.mrc";
	flpf_out = fn_out+filtered_mask_folder+"/class"+integerToString(cf.class_index)+"_lowpassfiltered.mrc";
	fp_out = fn_out+filtered_mask_folder+"/class"+integerToString(cf.class_index)+"_p_mask.mrc";
	fs_out = fn_out+filtered_mask_folder+"/class"+integerToString(cf.class_index)+"_s_mask.mrc";

	img.write(fimg_out);
	lpf_out.write(flpf_out);
	p_out.write(fp_out);
	s_out.write(fs_out);

}

void ClassRanker::maskCircumference(MultidimArray<int> p_mask, RFLOAT &protein_C, classFeatures cf, bool do_save_mask_c)
{
	// Sanity check
	Image<RFLOAT> c_out;
	c_out().initZeros(p_mask);
	c_out().setXmippOrigin();

	// Neighbours' relative coordinate. Notice that
	int dx[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
	int dy[8] = {0, 1, 1, 1, 0, -1, -1, -1};

	for (long int i= STARTINGY(p_mask) ; i <= FINISHINGY(p_mask); i++)
	{
		for (long int j= STARTINGX(p_mask) ; j <= FINISHINGX(p_mask); j++)
		{								// inside the protein area
			if (A2D_ELEM(p_mask, i, j) > 0.5)
			{
				for (int a=0; a<8; a++)
				{                               // Neighbours
					long y = i + dy[a];
					long x = j + dx[a];
					if (A2D_ELEM(p_mask, y, x)<0.5)     // If any neighbour is black
					{
						A2D_ELEM(c_out(), i, j) = 1;
						protein_C++;
						break;
					}
				}
			}
		}
	}
	if (do_save_mask_c)
	{
		std::string mask_c_fn_command = "mkdir -p mask_C";
		int result = system(mask_c_fn_command.c_str());

		FileName fc_out;
		fc_out = "./mask_C/class_"+integerToString(cf.class_index)+"_mask_c.mrc";
		c_out.write(fc_out);
	}
}

void ClassRanker::correctCtfUntilFirstPeak(MultidimArray<RFLOAT> &in, CTF ctf)
{

	FourierTransformer transformer;
	MultidimArray<Complex > Faux;

	RFLOAT xs = (RFLOAT)XSIZE(in) * mymodel.pixel_size;
	RFLOAT ys = (RFLOAT)YSIZE(in) * mymodel.pixel_size;
    transformer.FourierTransform(in, Faux, false);

    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(Faux)
    {
		RFLOAT x = (RFLOAT)jp / xs;
		RFLOAT y = (RFLOAT)ip / ys;

		// ctf should become part of the pvs_features class, perhaps call it avgctf
		RFLOAT ctf_val = ctf.getCTF(x, y, false, false, false, false, 0., true);
		// Don't use any value smaller than 7%. Some people accidentally use 0% amplitude contrast, which would be bad!
		if (ctf_val < 0.07) ctf_val = 0.07;
		DIRECT_A2D_ELEM(Faux, i, j) /= ctf_val;
    }
    transformer.inverseFourierTransform(Faux, in);
}

// Data local normalisation (wrote into a separate function as not sure if this can be done in a nicer way...)
void ClassRanker::localNormalisation(std::vector<classFeatures> &cf_all)
{
	// protein_sum
	RFLOAT ps_sum = 0.0, ps_mean = 0.0, ps_stddev = 0.0;
	for (int i = 0; i < cf_all.size(); i++) ps_sum += cf_all[i].protein_moments.sum;
	ps_mean = ps_sum / cf_all.size();
	for (int i = 0; i < cf_all.size(); i++) ps_stddev += pow((cf_all[i].protein_moments.sum - ps_mean), 2);
	ps_stddev = sqrt(ps_stddev / cf_all.size());
	for (int i = 0; i < cf_all.size(); i++) cf_all[i].normalized_features.protein_sum = (cf_all[i].protein_moments.sum - ps_mean) / ps_stddev;

	// solvent_sum
	RFLOAT ss_sum = 0.0, ss_mean = 0.0, ss_stddev = 0.0;
	for (int i = 0; i < cf_all.size(); i++) ss_sum += cf_all[i].solvent_moments.sum;
	ss_mean = ss_sum / cf_all.size();
	for (int i = 0; i < cf_all.size(); i++) ss_stddev += pow((cf_all[i].solvent_moments.sum - ss_mean), 2);
	ss_stddev = sqrt(ss_stddev / cf_all.size());
	for (int i = 0; i < cf_all.size(); i++) cf_all[i].normalized_features.solvent_sum = (cf_all[i].solvent_moments.sum - ss_mean) / ss_stddev;

	// relative_signal_intensity
	RFLOAT sum = 0.0, mean = 0.0, stddev = 0.0;
	for (int i = 0; i < cf_all.size(); i++) sum += cf_all[i].relative_signal_intensity;
	mean = sum / cf_all.size();
	for (int i = 0; i < cf_all.size(); i++) stddev += pow((cf_all[i].relative_signal_intensity - mean), 2);
	stddev = sqrt(stddev / cf_all.size());
	for (int i = 0; i < cf_all.size(); i++) cf_all[i].normalized_features.relative_signal_intensity = (cf_all[i].relative_signal_intensity - mean) / stddev;
}

/** ======================================================== Collecting All Features ========================================================= */

// Get features for non-empty classes
void ClassRanker::getFeatures()
{

	minRes = 999.0;
	features_all_classes.clear();
	features_all_classes.reserve(mymodel.nr_classes);

	if (verb > 0)
	{
		std::cout << " Calculating features for each class ..." << std::endl;
		init_progress_bar(end_class-start_class);
	}

	int ith_nonzero_class = 0;
	for (int iclass = start_class; iclass < end_class; iclass++)
	{
		if (debug > 0) std::cerr << " dealing with class: " << iclass+1 << std::endl;
		classFeatures features_this_class;

		// Get class distribution and particle number in the class and exclude classes with less than 10 particles (so that particle-number
// weighted resolution is sensible)
		features_this_class.class_distribution = mymodel.pdf_class[iclass];
		features_this_class.particle_nr = features_this_class.class_distribution * total_nr_particles;
		if (features_this_class.particle_nr > 10)
		{
			features_this_class.name = mymodel.ref_names[iclass];
			features_this_class.class_index = getClassIndex(features_this_class.name);
			Image<RFLOAT> img;
			img() = mymodel.Iref[iclass];

			// Get selection label (if training data)
			if (MD_select.numberOfObjects() > 0)
			{
				MD_select.getValue(EMDL_SELECTED, features_this_class.is_selected, iclass);
			}
			else
			{
				features_this_class.is_selected = 1;
			}

			// Get estimated resolution (regardless of whether it is already in model_classes table or not)
			if (mymodel.estimated_resolution[iclass] > 0.)
			{
				features_this_class.estimated_resolution = mymodel.estimated_resolution[iclass];
			}
			else
			{
				// TODO: this still relies on mlmodel!!!
				features_this_class.estimated_resolution = findResolution(features_this_class);
			}

			// Calculate particle number-weighted resolution
			features_this_class.weighted_resolution = (1. / (features_this_class.estimated_resolution*features_this_class.estimated_resolution)) / log(features_this_class.particle_nr);

			// Calculate image size weighted resolution
			features_this_class.relative_resolution = features_this_class.estimated_resolution / (mymodel.ori_size * mymodel.pixel_size);

			// Find job-wise best resolution among selected (red) classes in preparation for class score calculation called in the write_output function
			if (features_this_class.is_selected == 1 && features_this_class.estimated_resolution < minRes)
			{
				minRes = features_this_class.estimated_resolution;
			}

			if (do_skip_angular_errors)
			{
				features_this_class.accuracy_rotation = (preread_features_all_classes[ith_nonzero_class]).accuracy_rotation;
				features_this_class.accuracy_translation = (preread_features_all_classes[ith_nonzero_class]).accuracy_translation;
			}
			else
			{
				// Calculate class accuracy rotation and translation from model.star if present
				features_this_class.accuracy_rotation = mymodel.acc_rot[iclass];
				features_this_class.accuracy_translation = mymodel.acc_trans[iclass];
				if (debug>0) std::cerr << " mymodel.acc_rot[iclass]= " << mymodel.acc_rot[iclass] << " mymodel.acc_trans[iclass]= " << mymodel.acc_trans[iclass] << std::endl;
				if (features_this_class.accuracy_rotation > 99. || features_this_class.accuracy_translation > 99.)
				{
					calculateExpectedAngularErrors(iclass, features_this_class);
				}
				if (debug > 0) std::cerr << " done with angular errors" << std::endl;
			}

			// Now that we are going to calculate image-based features,
			// re-scale the image to have uniform pixel size of 4 angstrom
			int newsize = ROUND(XSIZE(img()) * (mymodel.pixel_size / uniform_angpix));
			newsize -= newsize%2; //make even in case it is not already
			resizeMap(img(), newsize);
			img().setXmippOrigin();

			// Determining radius to use
			circular_mask_radius = particle_diameter / (uniform_angpix * 2.);
			circular_mask_radius = std::min( (XSIZE(img())/2.) , circular_mask_radius);
			if (radius_ratio > 0 && radius <= 0) radius = radius_ratio * circular_mask_radius;
			// Calculate moments in ring area
			if (radius > 0)
			{
				features_this_class.ring_moments = calculateMoments(img(), radius, circular_mask_radius);
//				features_this_class.inner_circle_moments = calculateMoments(img(), 0, radius); // no longer written out
			}
			if (debug > 0) std::cerr << " done with ring moments" << std::endl;

			// Store the mean, stddev, minval and maxval of the lowpassed image as features
			MultidimArray<RFLOAT> lpf;
			lpf = img();
			lowPassFilterMap(lpf, lowpass, uniform_angpix);
			lpf.computeStats(features_this_class.lowpass_filtered_img_avg, features_this_class.lowpass_filtered_img_stddev,
					features_this_class.lowpass_filtered_img_minval, features_this_class.lowpass_filtered_img_maxval);

		 	// Make filtered masks
//			MultidimArray<RFLOAT> lpf;
			MultidimArray<int> p_mask, s_mask;
			long protein_area=0, solvent_area=0;
			makeSolventMasks(features_this_class, img(), lpf, p_mask, s_mask, features_this_class.scattered_signal, protein_area, solvent_area);
			// Protein and solvent area
			if (protein_area > 1) features_this_class.protein_area = 1;
			if (solvent_area > 0.08*3.14*circular_mask_radius*circular_mask_radius) features_this_class.solvent_area = 1;
			if (do_save_masks) saveMasks(img, lpf, p_mask, s_mask, features_this_class);

			// Circumference to area ratio
			RFLOAT protein_C = 0.;
			if (features_this_class.protein_area > 0.5)
			{
				maskCircumference(p_mask, protein_C, features_this_class, do_save_mask_c);
				features_this_class.CAR = protein_C / (2*sqrt(3.14*protein_area));
				// Debug
//				std::cerr << "Class " << features_this_class.class_index << ": protein area: " << protein_area << " mask circumference: " << protein_C << std::endl;
			}
			// Store entropy features on overall, protein and solvent region
			features_this_class.solvent_entropy = img().entropy(&s_mask);
			features_this_class.protein_entropy = img().entropy(&p_mask);
			features_this_class.total_entropy = img().entropy();

			// Moments for the protein and solvent area
			features_this_class.protein_moments = calculateMoments(img(), 0., circular_mask_radius, &p_mask);
			features_this_class.solvent_moments = calculateMoments(img(), 0., circular_mask_radius, &s_mask);

			// Signal intensity in the protein area relative to the solvent area
			features_this_class.relative_signal_intensity = features_this_class.protein_moments.sum - features_this_class.solvent_moments.mean*protein_area;

			// Fraction of white pixels in the protein mask on the edge
			long int edge_pix = 0, edge_white = 0;
			FOR_ALL_ELEMENTS_IN_ARRAY2D(p_mask)
			{
				if (round(sqrt(RFLOAT(i * i + j * j))) == round(circular_mask_radius))
				{
					edge_pix++;
					if (A2D_ELEM(p_mask, i, j) == 1) edge_white++;
				}
			}
			features_this_class.edge_signal = RFLOAT(edge_white) / RFLOAT(edge_pix);
			if (debug > 0) std::cerr << " done with edge signal" << std::endl;

			if (do_granularity_features)
			{
				// Calculate whole image LBP and protein and solvent area LBP
				calculatePvsLBP(img(), p_mask, s_mask, features_this_class);
				if (debug > 0) std::cerr << " done with lbp" << std::endl;

				// Calculate Haralick features
				if (debug>0) std::cerr << "Haralick features for protein area:" << std::endl;
				features_this_class.haralick_p = haralick_extractor.getHaralickFeatures(img(), &p_mask, debug>0);
				if (debug>0) std::cerr << "Haralick features for solvent area:" << std::endl;
				features_this_class.haralick_s = haralick_extractor.getHaralickFeatures(img(), &s_mask, debug>0);
				if (debug > 0) std::cerr << " done with haralick" << std::endl;

				// Calculate Zernike moments
				features_this_class.zernike_moments = zernike_extractor.getZernikeMoments(img(), 7, circular_mask_radius, debug>0);
				if (debug> 0 ) std::cerr << " done with Zernike moments" << std::endl;

				// Calculate granulo feature
				features_this_class.granulo = calculateGranulo(img());
			}
//			std::cout << "protein_area: " << features_this_class.protein_area << std::endl;
//			std::cout << "solvent_area: " << features_this_class.solvent_area << std::endl;

			// Collect selected features and apply global normalisation to part of them
			features_this_class.normalized_features.accuracy_rotation = (features_this_class.accuracy_rotation - global_mean[0]) / global_stddev[0];
			features_this_class.normalized_features.accuracy_translation = (features_this_class.accuracy_translation - global_mean[1]) / global_stddev[1];
			features_this_class.normalized_features.weighted_resolution = (features_this_class.weighted_resolution - global_mean[2]) / global_stddev[2];
			features_this_class.normalized_features.relative_resolution = (features_this_class.relative_resolution - global_mean[3]) / global_stddev[3];
			features_this_class.normalized_features.ring_mean = (features_this_class.ring_moments.mean - global_mean[4]) / global_stddev[4];
			features_this_class.normalized_features.ring_stddev = (features_this_class.ring_moments.stddev - global_mean[5]) / global_stddev[5];
			features_this_class.normalized_features.protein_stddev = (features_this_class.protein_moments.stddev - global_mean[6]) / global_stddev[6];
			features_this_class.normalized_features.solvent_mean = (features_this_class.solvent_moments.mean - global_mean[7]) / global_stddev[7];
			features_this_class.normalized_features.solvent_stddev = (features_this_class.solvent_moments.stddev - global_mean[8]) / global_stddev[8];
			features_this_class.normalized_features.scattered_signal = (features_this_class.scattered_signal - global_mean[9]) / global_stddev[9];
			features_this_class.normalized_features.edge_signal = (features_this_class.edge_signal - global_mean[10]) / global_stddev[10];
			features_this_class.normalized_features.lowpass_filtered_img_avg = (features_this_class.lowpass_filtered_img_avg - global_mean[11]) / global_stddev[11];
			features_this_class.normalized_features.lowpass_filtered_img_stddev = (features_this_class.lowpass_filtered_img_stddev - global_mean[12]) / global_stddev[12];
			features_this_class.normalized_features.lowpass_filtered_img_minval = (features_this_class.lowpass_filtered_img_minval - global_mean[13]) / global_stddev[13];
			features_this_class.normalized_features.lowpass_filtered_img_maxval = (features_this_class.lowpass_filtered_img_maxval - global_mean[14]) / global_stddev[14];
			features_this_class.normalized_features.granulo0 = (features_this_class.granulo[0] - global_granulo_mean[0]) / global_granulo_stddev[0];
			features_this_class.normalized_features.granulo1 = (features_this_class.granulo[1] - global_granulo_mean[1]) / global_granulo_stddev[1];
			features_this_class.normalized_features.granulo2 = (features_this_class.granulo[2] - global_granulo_mean[2]) / global_granulo_stddev[2];
			features_this_class.normalized_features.granulo3 = (features_this_class.granulo[3] - global_granulo_mean[3]) / global_granulo_stddev[3];
			features_this_class.normalized_features.granulo4 = (features_this_class.granulo[4] - global_granulo_mean[4]) / global_granulo_stddev[4];
			features_this_class.normalized_features.granulo5 = (features_this_class.granulo[5] - global_granulo_mean[5]) / global_granulo_stddev[5];

			// SHWS 15072020: new try small subimages with fixed boxsize at uniform_angpix for image-based CNN
			if (do_subimages && protein_area>0)
			{
				features_this_class.subimages = getSubimages(img(), subimage_boxsize, nr_subimages, &p_mask);
				if (debug> 0 ) std::cerr << " done with getSubimages" << std::endl;
			}

			// Push back the features of this class in the vector for all classes
			features_all_classes.push_back(features_this_class);

			ith_nonzero_class++;

			if (verb > 0)
				progress_bar(iclass-start_class+1);

		} // end if non-zero class

	} // end iterating all classes

	// Apply local normalisation for protein_sum, solvent_sum, and relative_signal_intensity
	ClassRanker::localNormalisation(features_all_classes);

	// If training, auto-labelled class score will be calculated and written out in writeFeatures()

	if (verb > 0)
		progress_bar(end_class-start_class);

}

// TODO: Liyi: make a read
void ClassRanker::readFeatures()
{

	preread_features_all_classes.clear();

	if (fn_cf == "") return;

	MetaDataTable MD_class_features;
	MD_class_features.read(fn_cf);

	int i =0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_class_features)
	{

		classFeatures this_class_feature;

		MD_class_features.getValue(EMDL_MLMODEL_REF_IMAGE, this_class_feature.name);
		MD_class_features.getValue(EMDL_CLASS_FEAT_CLASS_SCORE, this_class_feature.class_score);
		MD_class_features.getValue(EMDL_CLASS_FEAT_IS_SELECTED,this_class_feature.is_selected);
		MD_class_features.getValue(EMDL_CLASS_FEAT_CLASS_INDEX, this_class_feature.class_index);
		MD_class_features.getValue(EMDL_MLMODEL_PDF_CLASS, this_class_feature.class_distribution);
		MD_class_features.getValue(EMDL_MLMODEL_ACCURACY_ROT, this_class_feature.accuracy_rotation);
		MD_class_features.getValue(EMDL_MLMODEL_ACCURACY_TRANS, this_class_feature.accuracy_translation);
		MD_class_features.getValue(EMDL_MLMODEL_ESTIM_RESOL_REF, this_class_feature.estimated_resolution);
		MD_class_features.getValue(EMDL_CLASS_FEAT_WEIGHTED_RESOLUTION, this_class_feature.weighted_resolution);
		MD_class_features.getValue(EMDL_CLASS_FEAT_RELATIVE_RESOLUTION, this_class_feature.relative_resolution);


		// Moments for the ring
		if (radius > 0)
		{
			MD_class_features.getValue(EMDL_CLASS_FEAT_RING_MEAN, this_class_feature.ring_moments.mean);
			MD_class_features.getValue(EMDL_CLASS_FEAT_RING_STDDEV, this_class_feature.ring_moments.stddev);
			MD_class_features.getValue(EMDL_CLASS_FEAT_RING_SKEW, this_class_feature.ring_moments.skew);
			MD_class_features.getValue(EMDL_CLASS_FEAT_RING_KURT, this_class_feature.ring_moments.kurt);
		}

		// Protein and solvent region moments
		MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_AREA, this_class_feature.protein_area);
		MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_SUM, this_class_feature.protein_moments.sum);
		MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_MEAN, this_class_feature.protein_moments.mean);
		MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_STDDEV, this_class_feature.protein_moments.stddev);
		MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_SKEW, this_class_feature.protein_moments.skew);
		MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_KURT, this_class_feature.protein_moments.kurt);
		MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_AREA, this_class_feature.solvent_area);
		MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_SUM, this_class_feature.solvent_moments.sum);
		MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_MEAN, this_class_feature.solvent_moments.mean);
		MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_STDDEV, this_class_feature.solvent_moments.stddev);
		MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_SKEW, this_class_feature.solvent_moments.skew);
		MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_KURT, this_class_feature.solvent_moments.kurt);
		MD_class_features.getValue(EMDL_CLASS_FEAT_RELATIVE_SIGNAL_INT, this_class_feature.relative_signal_intensity);
		MD_class_features.getValue(EMDL_CLASS_FEAT_SCATTERED_SIGNAL, this_class_feature.scattered_signal);
		MD_class_features.getValue(EMDL_CLASS_FEAT_EDGE_SIGNAL, this_class_feature.edge_signal);
		MD_class_features.getValue(EMDL_CLASS_FEAT_CAR, this_class_feature.CAR);

        // Lowpass filtered image features
        MD_class_features.getValue(EMDL_CLASS_FEAT_LOWPASS_FILTERED_IMAGE_MEAN, this_class_feature.lowpass_filtered_img_avg);
        MD_class_features.getValue(EMDL_CLASS_FEAT_LOWPASS_FILTERED_IMAGE_STDDEV, this_class_feature.lowpass_filtered_img_stddev);
        MD_class_features.getValue(EMDL_CLASS_FEAT_LOWPASS_FILTERED_IMAGE_MIN, this_class_feature.lowpass_filtered_img_minval);
        MD_class_features.getValue(EMDL_CLASS_FEAT_LOWPASS_FILTERED_IMAGE_MAX, this_class_feature.lowpass_filtered_img_maxval);

        if (do_granularity_features)
        {
			// Protein and solvent region LBP's
			MD_class_features.getValue(EMDL_CLASS_FEAT_LBP, this_class_feature.lbp);
			MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_LBP, this_class_feature.lbp_p);
			MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_LBP, this_class_feature.lbp_s);

			// Protein and solvent region entropy
			MD_class_features.getValue(EMDL_CLASS_FEAT_TOTAL_ENTROPY, this_class_feature.total_entropy);
			MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_ENTROPY, this_class_feature.protein_entropy);
			MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_ENTROPY, this_class_feature.solvent_entropy);

			MD_class_features.getValue(EMDL_CLASS_FEAT_PROTEIN_HARALICK, this_class_feature.haralick_p);
			MD_class_features.getValue(EMDL_CLASS_FEAT_SOLVENT_HARALICK, this_class_feature.haralick_s);

			// Zernike moments
			MD_class_features.getValue(EMDL_CLASS_FEAT_ZERNIKE_MOMENTS, this_class_feature.zernike_moments);

			// Granulo
			MD_class_features.getValue(EMDL_CLASS_FEAT_GRANULO, this_class_feature.granulo);
        }

        preread_features_all_classes.push_back(this_class_feature);
		i++;
	}

}


float ClassRanker::deployTorchModel(FileName &model_path, std::vector<float> &features) {

#ifdef _TORCH_ENABLED
	torch::jit::script::Module module;

	// Deserialize model
	try {
		module = torch::jit::load(model_path);
	}
	catch (const c10::Error& e) {
		REPORT_ERROR("Error loading Torch model.");
	}

	// Create an inputs with batch size=1.
	std::vector<torch::jit::IValue> inputs;
	torch::Tensor t = torch::from_blob(features.data(), {1, (int)features.size()});
	inputs.push_back(t);

	// Execute the model and turn its output into a tensor.
	at::Tensor output = module.forward(inputs).toTensor();

	// Extract score value from tensor and return
	return output[0][0].item<float>();

#else //_TORCH_ENABLED
	REPORT_ERROR("RELION was not compiled with Torch support.\nConfigure with -DTORCH=ON and re-compile.");
#endif //_TORCH_ENABLED
}

void ClassRanker::performRanking()
{
	if (mydata.numberOfParticles() == 0)
	{
		// Read in particles if we haven't done this already
		mydata.read(fn_data, true, true); // true true means: ignore particle_name and group name!
		if (debug>0) std::cerr << "Done with reading data.star ..." << std::endl;
	}

	// Initialise all scores to -999 (including empty classes!)
	std::vector<RFLOAT> predicted_scores(mymodel.nr_classes, -999.);

	// to preserve original particle order in the data.star file
	if (do_select)
	{
		// Store original image order
		long int nr_parts = mydata.MDimg.numberOfObjects();
		for (long int j = 0; j < nr_parts; j++)
		{
			mydata.MDimg.setValue(EMDL_SORTED_IDX, j, j);
		}
	}

	MetaDataTable MDselected_particles, MDselected_classavgs;

	long int nr_sel_parts = 0;
	long int nr_sel_classavgs = 0;
	RFLOAT highscore = 0.0;
	std::vector<int> selected_classes;
	std::vector<float> scores(features_all_classes.size());
	float max_score = -999.;
	for (int i = 0; i < features_all_classes.size(); i++)
	{
		std::vector<float> feature_vector = features_all_classes[i].toVector();
		scores[i] = (RFLOAT) deployTorchModel(fn_torch_model, feature_vector);
		if (scores[i] > max_score) max_score = scores[i];
	}

	RFLOAT my_min = select_min_score;
	RFLOAT my_max = select_max_score;
	if (do_relative_threshold)
	{
		my_min = select_min_score * max_score;
		my_max = select_max_score * max_score;
	}

	for (int i = 0; i < features_all_classes.size(); i++)
	{
		if (do_select && scores[i] >= my_min && scores[i] <= my_max)
		{
			nr_sel_classavgs++;
			selected_classes.push_back(features_all_classes[i].class_index);

			MDselected_classavgs.addObject();
			MDselected_classavgs.setValue(EMDL_MLMODEL_REF_IMAGE, features_all_classes[i].name);
			MDselected_classavgs.setValue(EMDL_CLASS_PREDICTED_SCORE, scores[i]);
			MDselected_classavgs.setValue(EMDL_MLMODEL_PDF_CLASS, features_all_classes[i].class_distribution);
			MDselected_classavgs.setValue(EMDL_MLMODEL_ACCURACY_ROT, features_all_classes[i].accuracy_rotation);
			MDselected_classavgs.setValue(EMDL_MLMODEL_ACCURACY_TRANS_ANGSTROM, features_all_classes[i].accuracy_translation);
			MDselected_classavgs.setValue(EMDL_MLMODEL_ESTIM_RESOL_REF, features_all_classes[i].estimated_resolution);
		}

		// Set myscore in the vector that now runs over ALL classes (including empty ones)
		int iclass = features_all_classes[i].class_index - 1; // class counting in STAR files starts at 1!
		predicted_scores.at(iclass) = scores[i];
	}

	// Write optimiser.star and model.star in the output directory.
	FileName fn_opt_out, fn_model_out;
	fn_opt_out = fn_out + fn_root + "_optimiser.star";
	fn_model_out = fn_out + fn_root + "_model.star";
	MD_optimiser.setValue(EMDL_OPTIMISER_MODEL_STARFILE, fn_model_out);
	MD_optimiser.write(fn_opt_out);

	MetaDataTable MDclass;
	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
	{
		MDclass.addObject();
		MDclass.setValue(EMDL_MLMODEL_REF_IMAGE, mymodel.ref_names[iclass]);
		MDclass.setValue(EMDL_CLASS_FEAT_CLASS_SCORE, predicted_scores[iclass]);
		MDclass.setValue(EMDL_MLMODEL_PDF_CLASS, mymodel.pdf_class[iclass]);
		MDclass.setValue(EMDL_MLMODEL_ACCURACY_ROT, mymodel.acc_rot[iclass]);
		MDclass.setValue(EMDL_MLMODEL_ACCURACY_TRANS_ANGSTROM, mymodel.acc_trans[iclass]);
		MDclass.setValue(EMDL_MLMODEL_ESTIM_RESOL_REF, mymodel.estimated_resolution[iclass]);
		MDclass.setValue(EMDL_MLMODEL_FOURIER_COVERAGE_TOTAL_REF, mymodel.total_fourier_coverage[iclass]);
		if (mymodel.ref_dim==2)
		{
			MDclass.setValue(EMDL_MLMODEL_PRIOR_OFFX_CLASS, XX(mymodel.prior_offset_class[iclass]));
			MDclass.setValue(EMDL_MLMODEL_PRIOR_OFFY_CLASS, YY(mymodel.prior_offset_class[iclass]));
		}
	}
	MDclass.write(fn_model_out);
	if (verb > 0)
	{
		std::cout << " Written model with all predicted scores to: " << fn_opt_out << std::endl;
	}

	if (do_select)
	{

		// Select all particles in the data.star that have classes inside the selected_classes vector
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(mydata.MDimg)
		{
			int classnr;
			mydata.MDimg.getValue(EMDL_PARTICLE_CLASS, classnr);
			if (std::find(selected_classes.begin(), selected_classes.end(), classnr) != selected_classes.end())
			{
				MDselected_particles.addObject(mydata.MDimg.getObject());
				nr_sel_parts++;
			}
		}
		// Maintain the original image ordering and obsModel in mydata too
		MDselected_particles.sort(EMDL_SORTED_IDX);
		mydata.MDimg = MDselected_particles;
		mydata.write(fn_out+fn_sel_parts);

		// Also write out class_averages.star with the selected classes
		MDselected_classavgs.write(fn_out+fn_sel_classavgs);

		if (verb > 0)
		{
			std::cout << " Written " << nr_sel_parts << " selected particles to: " << fn_out << fn_sel_parts << std::endl;
			std::cout << " Written " << nr_sel_classavgs << " selected classes to: " << fn_out << fn_sel_classavgs << std::endl;
		}

	}
}


// Generate star file: write feature value of each of the features for all classes in the job out in the format of a star file
// Notice: not writing out normalized features as they are for c++ only execution. Hence the written out features still need to be normalized before used for training.
void ClassRanker::writeFeatures()
{
	MetaDataTable MD_class_features;
	MD_class_features.setName("class_features");
	for (int i=0; i<features_all_classes.size();i++)
	{
        // First calculate class scores if this is for training purpose (as the calculation depends on minRes, which can only be obtained at the end of the loop of getFeatures())
		if (!do_ranking)
        {
        	// for training: get a target class score from the jobscore and the resolution
    		// First calculate class score based on best resolution among red classes of the job, job score, estimated resolution, and the selection label
			if (features_all_classes[i].particle_nr > 10)
			{
				features_all_classes[i].class_score = getClassScoreFromJobScore(features_all_classes[i], minRes);
			}
        }

        MD_class_features.addObject();
		MD_class_features.setValue(EMDL_MLMODEL_REF_IMAGE, features_all_classes[i].name);
		if (!do_ranking)
		{
			MD_class_features.setValue(EMDL_CLASS_FEAT_CLASS_SCORE, features_all_classes[i].class_score);
		}
		MD_class_features.setValue(EMDL_CLASS_FEAT_IS_SELECTED,features_all_classes[i].is_selected);
		MD_class_features.setValue(EMDL_CLASS_FEAT_CLASS_INDEX, features_all_classes[i].class_index);
		MD_class_features.setValue(EMDL_MLMODEL_PDF_CLASS, features_all_classes[i].class_distribution);
		MD_class_features.setValue(EMDL_MLMODEL_ACCURACY_ROT, features_all_classes[i].accuracy_rotation);
		MD_class_features.setValue(EMDL_MLMODEL_ACCURACY_TRANS, features_all_classes[i].accuracy_translation);

		MD_class_features.setValue(EMDL_MLMODEL_ESTIM_RESOL_REF, features_all_classes[i].estimated_resolution);
		MD_class_features.setValue(EMDL_CLASS_FEAT_WEIGHTED_RESOLUTION, features_all_classes[i].weighted_resolution);
		MD_class_features.setValue(EMDL_CLASS_FEAT_RELATIVE_RESOLUTION, features_all_classes[i].relative_resolution);
		MD_class_features.setValue(EMDL_CLASS_FEAT_PARTICLE_NR, features_all_classes[i].particle_nr);

		// Moments for the ring
		if (radius > 0)
		{
			MD_class_features.setValue(EMDL_CLASS_FEAT_RING_MEAN, features_all_classes[i].ring_moments.mean);
			MD_class_features.setValue(EMDL_CLASS_FEAT_RING_STDDEV, features_all_classes[i].ring_moments.stddev);
			MD_class_features.setValue(EMDL_CLASS_FEAT_RING_SKEW, features_all_classes[i].ring_moments.skew);
			MD_class_features.setValue(EMDL_CLASS_FEAT_RING_KURT, features_all_classes[i].ring_moments.kurt);
		}

		// Protein and solvent region moments
		MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_AREA, features_all_classes[i].protein_area);
		MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_SUM, features_all_classes[i].protein_moments.sum);
		MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_MEAN, features_all_classes[i].protein_moments.mean);
		MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_STDDEV, features_all_classes[i].protein_moments.stddev);
		MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_SKEW, features_all_classes[i].protein_moments.skew);
		MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_KURT, features_all_classes[i].protein_moments.kurt);
		MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_AREA, features_all_classes[i].solvent_area);
		MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_SUM, features_all_classes[i].solvent_moments.sum);
		MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_MEAN, features_all_classes[i].solvent_moments.mean);
		MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_STDDEV, features_all_classes[i].solvent_moments.stddev);
		MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_SKEW, features_all_classes[i].solvent_moments.skew);
		MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_KURT, features_all_classes[i].solvent_moments.kurt);
		MD_class_features.setValue(EMDL_CLASS_FEAT_RELATIVE_SIGNAL_INT, features_all_classes[i].relative_signal_intensity);
		MD_class_features.setValue(EMDL_CLASS_FEAT_SCATTERED_SIGNAL, features_all_classes[i].scattered_signal);
		MD_class_features.setValue(EMDL_CLASS_FEAT_EDGE_SIGNAL, features_all_classes[i].edge_signal);
		MD_class_features.setValue(EMDL_CLASS_FEAT_CAR, features_all_classes[i].CAR);

        // Lowpass filtered image features
		MD_class_features.setValue(EMDL_CLASS_FEAT_LOWPASS_FILTERED_IMAGE_MEAN, features_all_classes[i].lowpass_filtered_img_avg);
		MD_class_features.setValue(EMDL_CLASS_FEAT_LOWPASS_FILTERED_IMAGE_STDDEV, features_all_classes[i].lowpass_filtered_img_stddev);
		MD_class_features.setValue(EMDL_CLASS_FEAT_LOWPASS_FILTERED_IMAGE_MIN, features_all_classes[i].lowpass_filtered_img_minval);
		MD_class_features.setValue(EMDL_CLASS_FEAT_LOWPASS_FILTERED_IMAGE_MAX, features_all_classes[i].lowpass_filtered_img_maxval);

		if (do_granularity_features)
		{
			// Protein and solvent region LBP's
			MD_class_features.setValue(EMDL_CLASS_FEAT_LBP, features_all_classes[i].lbp);
			MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_LBP, features_all_classes[i].lbp_p);
			MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_LBP, features_all_classes[i].lbp_s);

			// Protein and solvent region entropy
			MD_class_features.setValue(EMDL_CLASS_FEAT_TOTAL_ENTROPY, features_all_classes[i].total_entropy);
			MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_ENTROPY, features_all_classes[i].protein_entropy);
			MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_ENTROPY, features_all_classes[i].solvent_entropy);

			MD_class_features.setValue(EMDL_CLASS_FEAT_PROTEIN_HARALICK, features_all_classes[i].haralick_p);
			MD_class_features.setValue(EMDL_CLASS_FEAT_SOLVENT_HARALICK, features_all_classes[i].haralick_s);

			// Zernike moments
			MD_class_features.setValue(EMDL_CLASS_FEAT_ZERNIKE_MOMENTS, features_all_classes[i].zernike_moments);
			MD_class_features.setValue(EMDL_CLASS_FEAT_GRANULO, features_all_classes[i].granulo);
		}

	}
	MD_class_features.write(fn_out+fn_features);
	if (verb > 0) std::cout << " Written features to star file: " << fn_out << fn_features << std::endl;

	if (do_subimages)
	{
		MetaDataTable MD_subimages;
		MD_subimages.setName("subimages");
		long int nr_stack = 0;
		for (int f=0; f<features_all_classes.size();f++)
		{
			if (NZYXSIZE(features_all_classes[f].subimages) > 1)
			{
				int nr_images = NSIZE(features_all_classes[f].subimages);
				for (long n=0; n < nr_images; n++)
				{
					FileName fnt;
					fnt.compose(nr_stack+1, fn_out + fn_subimages);
					MD_subimages.addObject();
					MD_subimages.setValue(EMDL_IMAGE_NAME, fnt);
					MD_subimages.setValue(EMDL_MLMODEL_REF_IMAGE, features_all_classes[f].name);
					MD_subimages.setValue(EMDL_CLASS_FEAT_CLASS_SCORE, features_all_classes[f].class_score);

					Image<RFLOAT> img;
					features_all_classes[f].subimages.getImage(n, img());
					if (img().getSize() > 1)
					{
						img.write(fn_out + fn_subimages, -1, false, (nr_stack == 0) ? WRITE_OVERWRITE : WRITE_APPEND);
						nr_stack++;
					}
				}
			}
		}
		MD_subimages.write(fn_out + fn_subimage_star);
		if (verb > 0) std::cout << " Written subimages to: " << fn_out + fn_subimage_star << std::endl;

	}


}
