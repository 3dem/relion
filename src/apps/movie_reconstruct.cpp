/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres" "Takanori Nakane"
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

#include <src/backprojector.h>
#include <src/funcs.h>
#include <src/ctf.h>
#include <src/args.h>
#include <src/error.h>
#include <src/euler.h>
#include <src/time.h>
#include <src/ml_model.h>
#include <src/jaz/obs_model.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/motion/motion_helper.h>
#include <src/jaz/img_proc/filter_helper.h>

class MovieReconstructor
{
public:
	// I/O Parser
	IOParser parser;

	FileName fn_out, fn_sym, fn_sel, fn_gain, traj_path, movie_path;

	MetaDataTable DF;
	ObservationModel obsModel;

	int r_max, r_min_nn, blob_order, ref_dim, interpolator, iter,
	    debug_ori_size, debug_size,
	    ctf_dim, nr_helical_asu, width_mask_edge, nr_sectors, subset, chosen_class,
	    data_dim, output_boxsize, movie_boxsize, verb, frame;

	RFLOAT blob_radius, blob_alpha, angular_error, shift_error, angpix, maxres,
	       coord_angpix, movie_angpix, helical_rise, helical_twist;
	std::vector<double> data_angpixes;

	bool do_ctf, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak,
	     do_ewald, skip_weighting, skip_mask;

	bool skip_gridding, is_reverse, read_weights, do_external_reconstruct;

	float padding_factor, mask_diameter;

	// All backprojectors needed for parallel reconstruction
	BackProjector backprojector;

public:
	MovieReconstructor() { }

	// Read command line arguments
	void read(int argc, char **argv);

	// Initialise some stuff after reading
	void initialise();

	// Execute
	void run();

	// Loop over all particles to be back-projected
	void backproject(int rank = 0, int size = 1);

	// For parallelisation purposes
	void backprojectOneParticle(MetaDataTable &mdt, long int ipart, MultidimArray<Complex> &F2D);

	// perform the gridding reconstruction
	void reconstruct();

	void applyCTFPandCTFQ(MultidimArray<Complex> &Fin, CTF &ctf, FourierTransformer &transformer,
	                      MultidimArray<Complex> &outP, MultidimArray<Complex> &outQ, bool skip_mask=false);
};

int main(int argc, char *argv[])
{
	MovieReconstructor app;

	try
	{
		app.read(argc, argv);
		app.initialise();
		app.run();

	}
	catch (RelionError XE)
	{
		std::cerr << XE;
		return RELION_EXIT_FAILURE;
	}
	return RELION_EXIT_SUCCESS;
}

void MovieReconstructor::run()
{
	backproject(0, 1);

	reconstruct();
}

void MovieReconstructor::read(int argc, char **argv)
{
	parser.setCommandLine(argc, argv);

	int general_section = parser.addSection("General options");
	fn_sel = parser.getOption("--i", "Input STAR file with the projection images and their orientations", "");
	fn_out = parser.getOption("--o", "Name for output reconstruction","relion.mrc");
	fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
	maxres = textToFloat(parser.getOption("--maxres", "Maximum resolution (in Angstrom) to consider in Fourier space (default Nyquist)", "-1"));
	padding_factor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));
	traj_path = parser.getOption("--traj_path", "Trajectory path prefix", "");
	movie_path = parser.getOption("--movie_path", "Movie path prefix", "");
	fn_gain = parser.getOption("--gain", "Gain reference (in the right direction)", "");
	subset = textToInteger(parser.getOption("--subset", "Subset of images to consider (1: only reconstruct half1; 2: only half2; other: reconstruct all)", "-1"));
	movie_angpix = textToFloat(parser.getOption("--movie_angpix", "Pixel size in the movie", "-1"));
	if (movie_angpix < 0)
		REPORT_ERROR("For this program, you have to explicitly specify the movie pixel size (--movie_angpix).");
	coord_angpix = textToFloat(parser.getOption("--coord_angpix", "Pixel size of particle coordinates", "-1"));
	if (coord_angpix < 0)
		REPORT_ERROR("For this program, you have to explicitly specify the coordinate pixel size (--coord_angpix).");
	frame = textToInteger(parser.getOption("--frame", "Movie frame to reconstruct", "1"));
	movie_boxsize = textToInteger(parser.getOption("--window", "Box size to extract from raw movies", "-1"));
	if (movie_boxsize < 0 || movie_boxsize % 2 != 0)
		REPORT_ERROR("You have to specify the extraction box size (--window) as an even number.");
	output_boxsize = textToInteger(parser.getOption("--scale", "Box size after down-sampling", "-1"));
	if (output_boxsize < 0 || output_boxsize % 2 != 0)
		REPORT_ERROR("You have to specify the reconstruction box size (--scale) as an even number.");

	int ctf_section = parser.addSection("CTF options");
	do_ctf = parser.checkOption("--ctf", "Apply CTF correction");
	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");
	ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Images have been phase flipped");
	only_flip_phases = parser.checkOption("--only_flip_phases", "Do not correct CTF-amplitudes, only flip phases");

	int ewald_section = parser.addSection("Ewald-sphere correction options");
	do_ewald = parser.checkOption("--ewald", "Correct for Ewald-sphere curvature (developmental)");
	mask_diameter  = textToFloat(parser.getOption("--mask_diameter", "Diameter (in A) of mask for Ewald-sphere curvature correction", "-1."));
	width_mask_edge = textToInteger(parser.getOption("--width_mask_edge", "Width (in pixels) of the soft edge on the mask", "3"));
	is_reverse = parser.checkOption("--reverse_curvature", "Try curvature the other way around");
	nr_sectors = textToInteger(parser.getOption("--sectors", "Number of sectors for Ewald sphere correction", "2"));
	skip_mask = parser.checkOption("--skip_mask", "Do not apply real space mask during Ewald sphere correction");
	skip_weighting = parser.checkOption("--skip_weighting", "Do not apply weighting during Ewald sphere correction");
	if (verb > 0 && do_ewald && mask_diameter < 0 && !(skip_mask && skip_weighting))
		REPORT_ERROR("To apply Ewald sphere correction (--ewald), you have to specify the mask diameter(--mask_diameter).");

	int helical_section = parser.addSection("Helical options");
	nr_helical_asu = textToInteger(parser.getOption("--nr_helical_asu", "Number of helical asymmetrical units", "1"));
	helical_rise = textToFloat(parser.getOption("--helical_rise", "Helical rise (in Angstroms)", "0."));
	helical_twist = textToFloat(parser.getOption("--helical_twist", "Helical twist (in degrees, + for right-handedness)", "0."));

	int expert_section = parser.addSection("Expert options");
	if (parser.checkOption("--NN", "Use nearest-neighbour instead of linear interpolation before gridding correction"))
		interpolator = NEAREST_NEIGHBOUR;
	else
		interpolator = TRILINEAR;
	blob_radius = textToFloat(parser.getOption("--blob_r", "Radius of blob for gridding interpolation", "1.9"));
	blob_order = textToInteger(parser.getOption("--blob_m", "Order of blob for gridding interpolation", "0"));
	blob_alpha = textToFloat(parser.getOption("--blob_a", "Alpha-value of blob for gridding interpolation", "15"));
	iter = textToInteger(parser.getOption("--iter", "Number of gridding-correction iterations", "10"));
	ref_dim = 3;
	skip_gridding = parser.checkOption("--skip_gridding", "Skip gridding part of the reconstruction");
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

	// Hidden
	r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}

void MovieReconstructor::initialise()
{
	angpix = movie_angpix * movie_boxsize / output_boxsize;
	std::cout << "Movie box size = " << movie_boxsize << " px at " << movie_angpix << " A/px" << std::endl;
	std::cout << "Reconstruction box size = " << output_boxsize << " px at " << angpix << " A/px" << std::endl;
	std::cout << "Coordinate pixel size = " << coord_angpix << " A/px" << std::endl;

	// Read MetaData file, which should have the image names and their angles!
	ObservationModel::loadSafely(fn_sel, obsModel, DF, "particles", 0, false);
	std::cout << "Read " << DF.numberOfObjects() << " particles." << std::endl;
	data_angpixes = obsModel.getPixelSizes();

	if (verb > 0 && (subset == 1 || subset == 2) && !DF.containsLabel(EMDL_PARTICLE_RANDOM_SUBSET))
	{
		REPORT_ERROR("The rlnRandomSubset column is missing in the input STAR file.");
	}

	if (verb > 0 && (chosen_class >= 0) && !DF.containsLabel(EMDL_PARTICLE_CLASS))
	{
		REPORT_ERROR("The rlnClassNumber column is missing in the input STAR file.");
	}

	randomize_random_generator();

	if (do_ewald) do_ctf = true;
	data_dim = 2;

	if (maxres < 0.)
		r_max = -1;
	else
		r_max = CEIL(output_boxsize * angpix / maxres);
}

void MovieReconstructor::backproject(int rank, int size)
{
	backprojector = BackProjector(output_boxsize, ref_dim, fn_sym, interpolator,
					padding_factor, r_min_nn, blob_order,
					blob_radius, blob_alpha, data_dim, skip_gridding);
	backprojector.initZeros(2 * r_max);

	std::vector<MetaDataTable> mdts = StackHelper::splitByMicrographName(DF);
	
	const int nr_movies= mdts.size();
	if (verb > 0)
	{
		std::cout << " + Back-projecting all images ..." << std::endl;
		time_config();
		init_progress_bar(nr_movies);

	}

	FileName fn_img, fn_mic, fn_stack, fn_traj, fn_movie;
	long int stack_id;
	MetaDataTable trajectories;
	FourierTransformer transformer;
	Image<float> Iframe, Igain;
	Image<RFLOAT> Iparticle;
	Image<Complex> Fparticle;
	if (fn_gain != "")
		Igain.read(fn_gain);

	int frame_no = frame; // 1-indexed
	for (int imov = 0; imov < nr_movies; imov++)
	{	
		mdts[imov].getValue(EMDL_MICROGRAPH_NAME, fn_mic, 0);
		// TODO: Actually we should look at MotionCorr STAR file to find the right movie and gain filename
		//       For the time being, let's assume this is sorted out by Python script.
		fn_traj = traj_path + "/" + fn_mic.withoutExtension().getBaseName() + "_tracks.star";
		fn_movie = movie_path + "/" + fn_mic.withoutExtension().getBaseName() + ".tiff"; 
#ifdef DEBUG
		std::cout << "fn_mic = " << fn_mic << "\n\tfn_traj = " << fn_traj << "\n\tfn_movie = " << fn_movie << std::endl;
#endif
		// Read trajectories. Both particle ID and frame ID are 0-indexed in this array.
		std::vector<std::vector<gravis::d2Vector>> trajectories = MotionHelper::readTracksInPix(fn_traj, movie_angpix);

		// TODO: loop over relevant frames
		// TODO: support EER
		FileName fn_frame;
		fn_frame.compose(frame_no, fn_movie);
		Iframe.read(fn_frame);
		const int w0 = XSIZE(Iframe());
		const int h0 = YSIZE(Iframe());

		// Apply gain correction
		// Probably we can ignore defect correction, because we are not re-aligning.
		if (fn_gain == "")
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iframe())
				DIRECT_MULTIDIM_ELEM(Iframe(), n) *= -1;
		else
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iframe())
				DIRECT_MULTIDIM_ELEM(Iframe(), n) *= -DIRECT_MULTIDIM_ELEM(Igain(), n);

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(mdts[imov])
		{
#ifndef DEBUG
			progress_bar(imov);
#endif

			int random_subset = 0;
			mdts[imov].getValue(EMDL_PARTICLE_RANDOM_SUBSET, random_subset);

			// TODO: Use two threads for half sets to avoid reading movies twice??
			if (subset >= 1 && subset <= 2 && random_subset != subset)
				continue;

			const int opticsGroup = obsModel.getOpticsGroup(mdts[imov]); // 0-indexed
			const RFLOAT data_angpix = data_angpixes[opticsGroup];
			mdts[imov].getValue(EMDL_IMAGE_NAME, fn_img);
			fn_img.decompose(stack_id, fn_stack);
#ifdef DEBUG
			std::cout << "\tstack_id = " << stack_id << " fn_stack = " << fn_stack << std::endl;
#endif
			if (stack_id > trajectories.size())
				REPORT_ERROR("Missing trajectory!");

			RFLOAT coord_x, coord_y, origin_x, origin_y, traj_x, traj_y;
			mdts[imov].getValue(EMDL_IMAGE_COORD_X, coord_x); // in micrograph pixel
			mdts[imov].getValue(EMDL_IMAGE_COORD_Y, coord_y);
			mdts[imov].getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, origin_x); // in Angstrom
			mdts[imov].getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, origin_y);
#ifdef DEBUG
			std::cout << "\t\tcoord_mic_px = (" << coord_x << ", " << coord_y << ")";
			std::cout << " origin_angst = (" << origin_x << ", " << origin_y << ")";
			std::cout << " traj_movie_px = (" << trajectories[stack_id - 1][frame_no - 1].x <<  ", " << trajectories[stack_id - 1][frame_no - 1].y << ")" << std::endl;
#endif

			// Below might look overly complicated but is necessary to have the same rounding behaviour as Extract & Polish.
			Iparticle().initZeros(movie_boxsize, movie_boxsize);

			// ORIGINAL CODE (does not work when angpix changes)
			// pixel coordinate of the top left corner of the extractio box after down-sampling
			double xpO = (int)(coord_x * coord_angpix / angpix) - output_boxsize / 2;
			double ypO = (int)(coord_y * coord_angpix / angpix) - output_boxsize / 2;
			// pixel coordinate in the movie
			int x0 = (int)round(xpO * angpix / movie_angpix);
			int y0 = (int)round(ypO * angpix / movie_angpix);

			// Revised code: use data_angpix
			double xpOD = (int)(coord_x * coord_angpix / data_angpix);
			double ypOD = (int)(coord_y * coord_angpix / data_angpix);
			int x0D = (int)round(xpOD * data_angpix / movie_angpix) - movie_boxsize / 2;
			int y0D = (int)round(ypOD * data_angpix / movie_angpix) - movie_boxsize / 2;

			// pixel coordinate in the movie: cleaner but not compatible with existing files...
			int x0N = (int)round(coord_x * coord_angpix / movie_angpix) - movie_boxsize / 2;
			int y0N = (int)round(coord_y * coord_angpix / movie_angpix) - movie_boxsize / 2;
#ifdef DEBUG
			std::cout << "DEBUG: xpOD = " << xpO << " ypOD = " << ypO << std::endl;
			std::cout << "DEBUG: xpO  = " << xpO << " ypO  = " << ypO << std::endl;
			std::cout << "DEBUG: x0 = " << x0 << " y0 = " << y0 << " data_angpix = " << data_angpix << " angpix = " << angpix << std::endl;
			std::cout << "DEBUG: x0D = " << x0D << " y0D = " << y0D << " data_angpix = " << data_angpix << " angpix = " << angpix << std::endl;
			std::cout << "DEBUG: x0N = " << x0N << " y0N = " << y0N << std::endl;
#endif
			// Use revised code
			x0 = x0D; y0 = y0D;

			double dxM = trajectories[stack_id - 1][frame_no - 1].x;
			double dyM = trajectories[stack_id - 1][frame_no - 1].y;

			int dxI = (int)round(dxM);
			int dyI = (int)round(dyM);

			x0 += dxI;
			y0 += dyI;

			for (long int y = 0; y < movie_boxsize; y++)
			for (long int x = 0; x < movie_boxsize; x++)
			{
				int xx = x0 + x;
				int yy = y0 + y;

				if (xx < 0) xx = 0;
				else if (xx >= w0) xx = w0 - 1;

				if (yy < 0) yy = 0;
				else if (yy >= h0) yy = h0 - 1;

				DIRECT_NZYX_ELEM(Iparticle(), 0, 0, y, x) = DIRECT_NZYX_ELEM(Iframe(), 0, 0, yy, xx);
			}

			// Residual shifts in Angstrom. They don't contain OriginX/Y. Note the NEGATIVE sign.
			double dxR = - (dxM - dxI) * movie_angpix;
			double dyR = - (dyM - dyI) * movie_angpix;

			// Further shifts by OriginX/Y. Note that OriginX/Y are applied as they are 
			// (defined as "how much shift" we have to move particles).
			dxR += origin_x;
			dyR += origin_y;

			Iparticle().setXmippOrigin();
			transformer.FourierTransform(Iparticle(), Fparticle());
			if (output_boxsize != movie_boxsize) // TODO: FIXME: Something wrong. Output FSC is worse.
				Fparticle = FilterHelper::cropCorner2D(Fparticle, output_boxsize / 2 + 1, output_boxsize);
			shiftImageInFourierTransform(Fparticle(), Fparticle(), output_boxsize, dxR / angpix, dyR / angpix);
			CenterFFTbySign(Fparticle());

			backprojectOneParticle(mdts[imov], current_object, Fparticle());
		}
	}

	if (verb > 0)
		progress_bar(nr_movies);
}

void MovieReconstructor::backprojectOneParticle(MetaDataTable &mdt, long int p, MultidimArray<Complex> &F2D)
{
	RFLOAT rot, tilt, psi, fom, r_ewald_sphere;
	Matrix2D<RFLOAT> A3D;
	MultidimArray<RFLOAT> Fctf;
	Matrix1D<RFLOAT> trans(2);
	FourierTransformer transformer;

	// Rotations
	mdt.getValue(EMDL_ORIENT_ROT, rot, p);
	mdt.getValue(EMDL_ORIENT_TILT, tilt, p);
	mdt.getValue(EMDL_ORIENT_PSI, psi, p);
	Euler_angles2matrix(rot, tilt, psi, A3D);

	// If we are considering Ewald sphere curvature, the mag. matrix
	// has to be provided to the backprojector explicitly
	// (to avoid creating an Ewald ellipsoid)
	const bool ctf_premultiplied = false;
	const int opticsGroup = obsModel.getOpticsGroup(mdt, p);
	if (obsModel.getPixelSize(opticsGroup) != angpix)
	{
		std::cout << "before: " << obsModel.getPixelSize(opticsGroup) << " after " << angpix << std::endl;
		obsModel.setPixelSize(opticsGroup, angpix);
	}
	if (obsModel.getBoxSize(opticsGroup) != output_boxsize)
		obsModel.setBoxSize(opticsGroup, output_boxsize);
	//ctf_premultiplied = obsModel.getCtfPremultiplied(opticsGroup);
	
	if (do_ewald && ctf_premultiplied)
		REPORT_ERROR("We cannot perform Ewald sphere correction on CTF premultiplied particles.");
	Matrix2D<RFLOAT> magMat;
	if (!do_ewald)
	{
		A3D = obsModel.applyAnisoMag(A3D, opticsGroup);
	}

	// We don't need this, since we are backprojecting as is.
	/*
	std::cout << "before: " << A3D << std::endl;
	A3D = obsModel.applyScaleDifference(A3D, opticsGroup, output_boxsize, angpix);
	std::cout << "after: " << A3D << std::endl;
	*/

	MultidimArray<Complex> F2DP, F2DQ;
	FileName fn_img;

	Fctf.resize(F2D);
	Fctf.initConstant(1.);

	// Apply CTF if necessary
	if (do_ctf)
	{
		{
			CTF ctf;
			
			ctf.readByGroup(mdt, &obsModel, p);

			ctf.getFftwImage(Fctf, output_boxsize, output_boxsize, angpix,
			                 ctf_phase_flipped, only_flip_phases,
			                 intact_ctf_first_peak, true);

			obsModel.demodulatePhase(mdt, p, F2D); // This internally uses angpix!!
			obsModel.divideByMtf(mdt, p, F2D);

			// Ewald-sphere curvature correction
			if (do_ewald)
			{
				applyCTFPandCTFQ(F2D, ctf, transformer, F2DP, F2DQ, skip_mask);

				if (!skip_weighting)
				{
					// Also calculate W, store again in Fctf
					ctf.applyWeightEwaldSphereCurvature_noAniso(Fctf, output_boxsize, output_boxsize, angpix, mask_diameter);
				}

				// Also calculate the radius of the Ewald sphere (in pixels)
				r_ewald_sphere = output_boxsize * angpix / ctf.lambda;
			}
		}
	}

	if (true) // not subtract
	{
		if (do_ewald)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
			{
				DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
			}
		}
		else if (do_ctf) // "Normal" reconstruction, multiply X by CTF, and W by CTF^2
		{
			if (!ctf_premultiplied)
			{
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
				{
					DIRECT_MULTIDIM_ELEM(F2D, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
				}
			}
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fctf)
			{
				DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
			}
		}

		DIRECT_A2D_ELEM(F2D, 0, 0) = 0.0;

		if (do_ewald)
		{
			Matrix2D<RFLOAT> magMat;

			if (obsModel.hasMagMatrices)
			{
				magMat = obsModel.getMagMatrix(opticsGroup);
			}
			else
			{
				magMat = Matrix2D<RFLOAT>(2,2);
				magMat.initIdentity();
			}

			backprojector.set2DFourierTransform(F2DP, A3D, &Fctf, r_ewald_sphere, true, &magMat);
			backprojector.set2DFourierTransform(F2DQ, A3D, &Fctf, r_ewald_sphere, false, &magMat);
		}
		else
		{
			backprojector.set2DFourierTransform(F2D, A3D, &Fctf);
		}
	}
}

void MovieReconstructor::reconstruct()
{
	bool do_map = false;
	bool do_use_fsc = false;
	MultidimArray<RFLOAT> fsc, dummy;
	Image<RFLOAT> vol;
	fsc.resize(output_boxsize/2+1);

	if (verb > 0)
		std::cout << " + Starting the reconstruction ..." << std::endl;

	backprojector.symmetrise(nr_helical_asu, helical_twist, helical_rise/angpix);

	MultidimArray<RFLOAT> tau2;
	backprojector.reconstruct(vol(), iter, do_map, tau2);

	vol.setSamplingRateInHeader(angpix);
	vol.write(fn_out);
	if (verb > 0)
		std::cout << " + Done! Written output map in: "<<fn_out<<std::endl;
}



void MovieReconstructor::applyCTFPandCTFQ(MultidimArray<Complex> &Fin, CTF &ctf, FourierTransformer &transformer,
                                          MultidimArray<Complex> &outP, MultidimArray<Complex> &outQ, bool skip_mask)
{
	//FourierTransformer transformer;
	outP.resize(Fin);
	outQ.resize(Fin);
	float angle_step = 180./nr_sectors;
	for (float angle = 0.; angle < 180.;  angle +=angle_step)
	{
		MultidimArray<Complex> CTFP(Fin), Fapp(Fin);
		MultidimArray<RFLOAT> Iapp(YSIZE(Fin), YSIZE(Fin));
		// Two passes: one for CTFP, one for CTFQ
		for (int ipass = 0; ipass < 2; ipass++)
		{
			bool is_my_positive = (ipass == 1) ? is_reverse : !is_reverse;

			// Get CTFP and multiply the Fapp with it
			ctf.getCTFPImage(CTFP, YSIZE(Fin), YSIZE(Fin), angpix, is_my_positive, angle);

			Fapp = Fin * CTFP; // element-wise complex multiplication!

			if (!skip_mask)
			{
				// inverse transform and mask out the particle....
				CenterFFTbySign(Fapp);
				transformer.inverseFourierTransform(Fapp, Iapp);

				softMaskOutsideMap(Iapp, ROUND(mask_diameter/(angpix*2.)), (RFLOAT)width_mask_edge);

				// Re-box to a smaller size if necessary....
				if (output_boxsize < YSIZE(Fin))
				{
					Iapp.setXmippOrigin();
					Iapp.window(FIRST_XMIPP_INDEX(output_boxsize), FIRST_XMIPP_INDEX(output_boxsize),
					            LAST_XMIPP_INDEX(output_boxsize),  LAST_XMIPP_INDEX(output_boxsize));

				}
				// Back into Fourier-space
				transformer.FourierTransform(Iapp, Fapp, false); // false means: leave Fapp in the transformer
				CenterFFTbySign(Fapp);
			}

			// First time round: resize the output arrays
			if (ipass == 0 && fabs(angle) < XMIPP_EQUAL_ACCURACY)
			{
				outP.resize(Fapp);
				outQ.resize(Fapp);
			}

			// Now set back the right parts into outP (first pass) or outQ (second pass)
			float anglemin = angle + 90. - (0.5*angle_step);
			float anglemax = angle + 90. + (0.5*angle_step);

			// angles larger than 180
			bool is_reverse = false;
			if (anglemin >= 180.)
			{
				anglemin -= 180.;
				anglemax -= 180.;
				is_reverse = true;
			}
			MultidimArray<Complex> *myCTFPorQ, *myCTFPorQb;
			if (is_reverse)
			{
				myCTFPorQ  = (ipass == 0) ? &outQ : &outP;
				myCTFPorQb = (ipass == 0) ? &outP : &outQ;
			}
			else
			{
				myCTFPorQ  = (ipass == 0) ? &outP : &outQ;
				myCTFPorQb = (ipass == 0) ? &outQ : &outP;
			}

			// Deal with sectors with the Y-axis in the middle of the sector...
			bool do_wrap_max = false;
			if (anglemin < 180. && anglemax > 180.)
			{
				anglemax -= 180.;
				do_wrap_max = true;
			}

			// use radians instead of degrees
			anglemin = DEG2RAD(anglemin);
			anglemax = DEG2RAD(anglemax);
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(CTFP)
			{
				RFLOAT x = (RFLOAT)jp;
				RFLOAT y = (RFLOAT)ip;
				RFLOAT myangle = (x*x+y*y > 0) ? acos(y/sqrt(x*x+y*y)) : 0; // dot-product with Y-axis: (0,1)
				// Only take the relevant sector now...
				if (do_wrap_max)
				{
					if (myangle >= anglemin)
						DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					else if (myangle < anglemax)
						DIRECT_A2D_ELEM(*myCTFPorQb, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
				}
				else
				{
					if (myangle >= anglemin && myangle < anglemax)
						DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
				}
			}
		}
	}
}


