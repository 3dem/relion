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
#include "spa_backproject.h"
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/image/symmetry.h>
#include <omp.h>

using namespace gravis;


void SpaBackproject::read(int argc, char **argv)
{
	parser.setCommandLine(argc, argv);

	int general_section = parser.addSection("General options");
	
	fn_sel = parser.getOption("--i", "Input STAR file with the projection images and their orientations", "");
	fn_out = parser.getOption("--o", "Name prefix for output reconstructions");
	fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
	maxres = textToFloat(parser.getOption("--maxres", "Maximum resolution (in Angstrom) to consider in Fourier space (default Nyquist)", "-1"));
	padding_factor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));
	image_path = parser.getOption("--img", "Optional: image path prefix", "");
	subset = textToInteger(parser.getOption("--subset", "Subset of images to consider (1: only reconstruct half1; 2: only half2; other: reconstruct all)", "-1"));
	chosen_class = textToInteger(parser.getOption("--class", "Consider only this class (-1: use all classes)", "-1"));
	angpix  = textToFloat(parser.getOption("--angpix", "Pixel size in the reconstruction (take from first optics group by default)", "-1"));

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
	newbox = textToInteger(parser.getOption("--newbox", "Box size of reconstruction after Ewald sphere correction", "-1"));
	nr_sectors = textToInteger(parser.getOption("--sectors", "Number of sectors for Ewald sphere correction", "2"));
	skip_mask = parser.checkOption("--skip_mask", "Do not apply real space mask during Ewald sphere correction");
	skip_weighting = parser.checkOption("--skip_weighting", "Do not apply weighting during Ewald sphere correction");

	if (verb > 0 && do_ewald && mask_diameter < 0 && !(skip_mask && skip_weighting))
	{
		REPORT_ERROR("To apply Ewald sphere correction (--ewald), you have to specify the mask diameter(--mask_diameter).");
	}

	int helical_section = parser.addSection("Helical options");
	
	nr_helical_asu = textToInteger(parser.getOption("--nr_helical_asu", "Number of helical asymmetrical units", "1"));
	helical_rise = textToFloat(parser.getOption("--helical_rise", "Helical rise (in Angstroms)", "0."));
	helical_twist = textToFloat(parser.getOption("--helical_twist", "Helical twist (in degrees, + for right-handedness)", "0."));

	int expert_section = parser.addSection("Expert options");
	
	fn_sub = parser.getOption("--subtract","Subtract projections of this map from the images used for reconstruction", "");
	if (parser.checkOption("--NN", "Use nearest-neighbour instead of linear interpolation before gridding correction"))
	{
		interpolator = NEAREST_NEIGHBOUR;
	}
	else
	{
		interpolator = TRILINEAR;
	}
	
	blob_radius = textToFloat(parser.getOption("--blob_r", "Radius of blob for gridding interpolation", "1.9"));
	blob_order = textToInteger(parser.getOption("--blob_m", "Order of blob for gridding interpolation", "0"));
	blob_alpha = textToFloat(parser.getOption("--blob_a", "Alpha-value of blob for gridding interpolation", "15"));
	iter = textToInteger(parser.getOption("--iter", "Number of gridding-correction iterations", "10"));
	ref_dim = textToInteger(parser.getOption("--refdim", "Dimension of the reconstruction (2D or 3D)", "3"));
	angular_error = textToFloat(parser.getOption("--angular_error", "Apply random deviations with this standard deviation (in degrees) to each of the 3 Euler angles", "0."));
	shift_error = textToFloat(parser.getOption("--shift_error", "Apply random deviations with this standard deviation (in Angstrom) to each of the 2 translations", "0."));
	do_fom_weighting = parser.checkOption("--fom_weighting", "Weight particles according to their figure-of-merit (_rlnParticleFigureOfMerit)");
	fn_fsc = parser.getOption("--fsc", "FSC-curve for regularized reconstruction", "");
	ctf_dim  = textToInteger(parser.getOption("--reconstruct_ctf", "Perform a 3D reconstruction from 2D CTF-images, with the given size in pixels", "-1"));
	do_reconstruct_ctf2 = parser.checkOption("--ctf2", "Reconstruct CTF^2 and then take the sqrt of that");
	skip_gridding = parser.checkOption("--skip_gridding", "Skip gridding part of the reconstruction");
	fn_noise = parser.getOption("--reconstruct_noise","Reconstruct noise using sigma2 values in this model STAR file", "");
	read_weights = parser.checkOption("--read_weights", "Developmental: read freq. weight files");
	do_debug = parser.checkOption("--write_debug_output", "Write out arrays with data and weight terms prior to reconstruct");
	do_external_reconstruct = parser.checkOption("--external_reconstruct", "Write out BP denominator and numerator for external_reconstruct program");
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));
	
	num_threads_in = textToInteger(parser.getOption("--j_in", "Number of inner threads", "1"));
	num_threads_out = textToInteger(parser.getOption("--j_out", "Number of outer threads", "1"));
	num_threads_total = num_threads_in * num_threads_out;
	
	use_fwd_mapping = parser.checkOption("--fwd_mapping", "Use the legacy forward-mapping algorithm");

	// Hidden
	r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));
		
	Log::readParams(parser);

	// Check for errors in the command-line option
	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
}

void SpaBackproject::usage()
{
	parser.writeUsage(std::cout);
}

void SpaBackproject::initialise()
{
	do_reconstruct_ctf = (ctf_dim > 0);
	
	if (do_reconstruct_ctf)
	{
		do_ctf = false;
		padding_factor = 1.;
	}

	ObservationModel::loadSafely(fn_sel, obsModel, DF, "particles", 0, false);
	
	if (obsModel.opticsMdt.numberOfObjects() == 0)
	{
		REPORT_ERROR("Optics table has zero entries");
	}

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

	if (fn_noise != "")
	{
		model.read(fn_noise);
	}

	// Get dimension of the images
	if (do_reconstruct_ctf)
	{
		output_boxsize = ctf_dim;
	}
	else
	{
		(DF).firstObject();
		DF.getValue(EMDL_IMAGE_NAME, fn_img);

		if (image_path != "")
		{
			fn_img = image_path + "/" + fn_img.substr(fn_img.find_last_of("/")+1);
		}

		Image<RFLOAT> img0;
		img0.read(fn_img, false);
		output_boxsize = (int)XSIZE(img0());
		
		// Allow reconstructing smaller box than the input images (which should have large boxes!)
		if (newbox > 0)
		{
			output_boxsize = newbox;
		}

		const int data_dim = obsModel.opticsMdt.getInt(EMDL_IMAGE_DIMENSIONALITY, 0);
		
		if (data_dim != 2)
		{
			REPORT_ERROR("Only 2D images are supported.");
		}
	}

	if (angpix < 0.)
	{
		angpix = obsModel.getPixelSize(0);
		Log::print("Reconstructing at the pixel size of the first optics group: " + ZIO::itoa(angpix) + " Å");
	}

	if (maxres < 0.)
	{
		r_max = -1;
	}
	else
	{
		r_max = CEIL(output_boxsize * angpix / maxres);
	}
}

void SpaBackproject::run()
{
	initialise();
	backprojectAllParticles();
	
	if (use_fwd_mapping)
	{
		reconstructForward();
	}
	else
	{
		reconstructBackward();
	}
}

void SpaBackproject::backprojectAllParticles()
{
	if (fn_sub != "")
	{
		projector = Projector(output_boxsize, interpolator, padding_factor, r_min_nn);
		Image<RFLOAT> sub;
		sub.read(fn_sub);
		MultidimArray<RFLOAT> dummy;
		projector.computeFourierTransformMap(sub(), dummy, 2 * r_max);
	}

	if (use_fwd_mapping)
	{
		backprojectors = std::vector<BackProjector>(
					2 * num_threads_out, 
					BackProjector(
						output_boxsize, ref_dim, fn_sym, interpolator,
						padding_factor, r_min_nn, blob_order,
						blob_radius, blob_alpha, 2, skip_gridding));
		
		for (int i = 0; i < backprojectors.size(); i++)
		{
			backprojectors[i].initZeros(2 * r_max);
		}
	}
	else
	{
		padded_box_size = (int)(padding_factor * output_boxsize);
		
		const int s = padded_box_size;
		const int sh = s/2 + 1;
		
		accumulation_volumes = std::vector<AccumulationVolume>(2 * num_threads_out);
		
		for (int i = 0; i < accumulation_volumes.size(); i++)
		{
			accumulation_volumes[i].data = BufferedImage<Complex>(sh,s,s);
			accumulation_volumes[i].weight = BufferedImage<RFLOAT>(sh,s,s);
			accumulation_volumes[i].multiplicity = BufferedImage<RFLOAT>(sh,s,s);
			accumulation_volumes[i].spreading_function = BufferedImage<RFLOAT>(sh,s,s);
			
			accumulation_volumes[i].data.fill(Complex(0.0, 0.0));
			accumulation_volumes[i].weight.fill(0.0);
			accumulation_volumes[i].multiplicity.fill(0.0);
			accumulation_volumes[i].spreading_function.fill(0.0);
		}
	}

	long int nr_parts = DF.numberOfObjects();
	
	if (verb > 0)
	{
		time_config();
		Log::beginProgress("Back-projecting all images", nr_parts / num_threads_out);
	}
	
	#pragma omp parallel for num_threads(num_threads_out)
	for (long int p = 0; p < nr_parts; p++)
	{
		const int th = omp_get_thread_num();
		
		backprojectOneParticle(p, th);

		if (th == 0 && verb > 0)
		{
			Log::updateProgress(p);
		}
	}

	if (verb > 0)
	{
		Log::endProgress();
	}
}

void SpaBackproject::backprojectOneParticle(long int p, int thread_id)
{
	if (chosen_class >= 0 && chosen_class != DF.getInt(EMDL_PARTICLE_CLASS, p))
	{
		return;
	}
	
	const int subset = DF.getIntMinusOne(EMDL_PARTICLE_RANDOM_SUBSET, p);
	
	const int data_id = 2*thread_id + subset;
	
	BackProjector* backprojector;
	AccumulationVolume* accumulationVolume;
	
	if (use_fwd_mapping)
	{
		backprojector = &backprojectors[data_id];
	}
	else
	{
		accumulationVolume = &accumulation_volumes[data_id];
	}
	
	Matrix2D<RFLOAT> A3D;
	
	// Rotations
	{
		RFLOAT rot, tilt, psi;

		if (ref_dim == 2)
		{
			rot = tilt = 0.;
		}
		else
		{
			DF.getValue(EMDL_ORIENT_ROT, rot, p);
			DF.getValue(EMDL_ORIENT_TILT, tilt, p);
		}
	
		psi = 0.;
		DF.getValue(EMDL_ORIENT_PSI, psi, p);
	
		if (angular_error > 0.)
		{
			rot += rnd_gaus(0., angular_error);
			tilt += rnd_gaus(0., angular_error);
			psi += rnd_gaus(0., angular_error);
		}
	
		Euler_angles2matrix(rot, tilt, psi, A3D);
	}

	// If we are considering Ewald sphere curvature, the mag. matrix
	// has to be provided to the backprojector explicitly
	// (to avoid creating an Ewald ellipsoid)	
	
	const int opticsGroup = obsModel.getOpticsGroup(DF, p);
	const int myBoxSize = obsModel.getBoxSize(opticsGroup);
	const RFLOAT myPixelSize = obsModel.getPixelSize(opticsGroup);
	const bool ctf_premultiplied = obsModel.getCtfPremultiplied(opticsGroup);
	
	if (do_ewald && ctf_premultiplied)
	{
		REPORT_ERROR("We cannot perform Ewald sphere correction on CTF premultiplied particles.");
	}
	
	if (!do_ewald)
	{
		A3D = obsModel.applyAnisoMag(A3D, opticsGroup);
	}
	
	A3D = obsModel.applyScaleDifference(A3D, opticsGroup, output_boxsize, angpix);
	
	
	
	MultidimArray<Complex> F2D, F2DP, F2DQ;
	FileName fn_img;
	Image<RFLOAT> img;
	FourierTransformer transformer;

	if (!do_reconstruct_ctf && fn_noise == "")
	{
		RFLOAT tx = DF.getRfloat(EMDL_ORIENT_ORIGIN_X_ANGSTROM, p);
		RFLOAT ty = DF.getRfloat(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, p);
		
		if (shift_error > 0.)
		{
			tx += rnd_gaus(0., shift_error);
			ty += rnd_gaus(0., shift_error);
		}
	
		// As of v3.1, shifts are in Angstroms in the STAR files, convert back to pixels here
		tx /= myPixelSize;
		ty /= myPixelSize;
		
		DF.getValue(EMDL_IMAGE_NAME, fn_img, p);
		img.read(fn_img);
		img().setXmippOrigin();
		transformer.FourierTransform(img(), F2D);
		CenterFFTbySign(F2D);

		shiftImageInFourierTransform(F2D, F2D, img.data.xdim, tx, ty, 0);
	}
	else
	{
		F2D.resize(myBoxSize, myBoxSize / 2 + 1);
	}

	if (fn_noise != "")
	{
		FileName fn_group;
		
		if (DF.containsLabel(EMDL_MLMODEL_GROUP_NAME))
		{
			DF.getValue(EMDL_MLMODEL_GROUP_NAME, fn_group);
		}
		else if (DF.containsLabel(EMDL_MICROGRAPH_NAME))
		{
			DF.getValue(EMDL_MICROGRAPH_NAME, fn_group);
		}
		else
		{
			REPORT_ERROR("ERROR: cannot find rlnGroupName or rlnMicrographName in the input --i file...");
		}

		int my_mic_id = -1;
		
		for (int mic_id = 0; mic_id < model.group_names.size(); mic_id++)
		{
			if (fn_group == model.group_names[mic_id])
			{
				my_mic_id = mic_id;
				break;
			}
		}

		if (my_mic_id < 0) 
		{
			REPORT_ERROR("ERROR: cannot find " + fn_group + " in the input model file...");
		}

		RFLOAT normcorr = 1.;
		
		if (DF.containsLabel(EMDL_IMAGE_NORM_CORRECTION)) 
		{
			DF.getValue(EMDL_IMAGE_NORM_CORRECTION, normcorr);
		}

		// Make coloured-noise image
		for (long int k = 0, kp = 0; k<ZSIZE(F2D); k++, kp = (k < XSIZE(F2D)) ? k : k - ZSIZE(F2D))
		for (long int i = 0, ip = 0; i<YSIZE(F2D); i++, ip = (i < XSIZE(F2D)) ? i : i - YSIZE(F2D))
		for (long int j = 0, jp = 0; j<XSIZE(F2D); j++, jp = j)
		{
			int ires = ROUND(sqrt((RFLOAT)(kp*kp + ip*ip + jp*jp)));
			ires = XMIPP_MIN(ires, myBoxSize/2); // at freqs higher than Nyquist: use last sigma2 value

			RFLOAT sigma = sqrt(DIRECT_A1D_ELEM(model.sigma2_noise[my_mic_id], ires));
			DIRECT_A3D_ELEM(F2D, k, i, j).real += rnd_gaus(0., sigma);
			DIRECT_A3D_ELEM(F2D, k, i, j).imag += rnd_gaus(0., sigma);
		}
	}
	
	MultidimArray<RFLOAT> Fctf;
	Fctf.resize(F2D);
	Fctf.initConstant(1.);
	
	RawImage<RFLOAT> ctfImage(Fctf);
	RawImage<Complex> dataImage(F2D);
		
	RFLOAT r_ewald_sphere = std::numeric_limits<RFLOAT>::max();
	
	// Apply CTF if necessary
	if (do_ctf || do_reconstruct_ctf)
	{
		CTF ctf;		
		ctf.readByGroup(DF, &obsModel, p);
		ctf.getFftwImage(
				Fctf, myBoxSize, myBoxSize, myPixelSize,
				ctf_phase_flipped, only_flip_phases,
				intact_ctf_first_peak, true);

		obsModel.demodulatePhase(DF, p, F2D);
		obsModel.divideByMtf(DF, p, F2D);

		// Ewald-sphere curvature correction
		if (do_ewald)
		{
			applyCTFPandCTFQ(F2D, ctf, transformer, F2DP, F2DQ, skip_mask);

			if (!skip_weighting)
			{
				// Also calculate W, store again in Fctf
				ctf.applyWeightEwaldSphereCurvature_noAniso(Fctf, myBoxSize, myBoxSize, myPixelSize, mask_diameter);
			}

			// Also calculate the radius of the Ewald sphere (in pixels)
			r_ewald_sphere = myBoxSize * myPixelSize / ctf.lambda;
		}
	}

	// Subtract reference projection
	if (fn_sub != "")
	{
		MultidimArray<Complex> Fsub;		
		Fsub.resize(F2D);
		projector.get2DFourierTransform(Fsub, A3D);
		
		RawImage<Complex> clutterImage(Fsub);
		
		if (do_ctf)
		{
			clutterImage *= ctfImage;
		}

		dataImage -= clutterImage;
		
		if (use_fwd_mapping)
		{
			backprojector->set2DFourierTransform(F2D, A3D);
		}
		else
		{
			REPORT_ERROR("Backward mapping not implemented yet.");
		}
	}
	else
	{
		if (do_reconstruct_ctf)
		{
			dataImage.copyFrom(ctfImage);
			
			if (do_reconstruct_ctf2)
			{
				dataImage *= ctfImage;
			}
				
			ctfImage.fill(1.f);
		}
		else if (do_ewald)
		{
			ctfImage *= ctfImage;
		}
		// "Normal" reconstruction, multiply X by CTF, and W by CTF^2
		else if (do_ctf)
		{
			if (!ctf_premultiplied)
			{
				dataImage *= ctfImage;
			}
			
			ctfImage *= ctfImage;
		}

		// Do the following after squaring the CTFs!
		if (do_fom_weighting)
		{
			const RFLOAT fom = DF.getRfloat(EMDL_PARTICLE_FOM, p);
			
			dataImage *= fom;
			ctfImage *= fom;
		}

		if (read_weights)
		{
			std::string name, fullName;

			DF.getValue(EMDL_IMAGE_NAME, fullName, 0);
			name = fullName.substr(fullName.find("@")+1);

			if (image_path != "")
			{
				name = image_path + "/" + name.substr(name.find_last_of("/")+1);
			}

			std::string wghName = name;
			wghName = wghName.substr(0, wghName.find_last_of('.')) + "_weight.mrc";

			Image<RFLOAT> wgh;
			wgh.read(wghName);
			RawImage<RFLOAT> weightImage(wgh);

			if (!weightImage.hasEqualSize(ctfImage))
			{
				REPORT_ERROR(wghName + " and " + name + " are of unequal size.\n");
			}

			ctfImage *= weightImage;
		}

		dataImage(0,0) = 0.0;

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

			if (use_fwd_mapping)
			{
				backprojector->set2DFourierTransform(F2DP, A3D, &Fctf, r_ewald_sphere, true, &magMat);
				backprojector->set2DFourierTransform(F2DQ, A3D, &Fctf, r_ewald_sphere, false, &magMat);
			}
			else
			{
				REPORT_ERROR("Ewald's sphere curvature is currently only supported with forward mapping.");
			}
		}
		else
		{
			if (use_fwd_mapping)
			{
				backprojector->set2DFourierTransform(F2D, A3D, &Fctf);
			}
			else
			{
				d4Matrix proj(
						A3D(0,0), A3D(0,1), A3D(0,2), 0,
						A3D(1,0), A3D(1,1), A3D(1,2), 0,
						A3D(2,0), A3D(2,1), A3D(2,2), 0,
						       0,        0,        0, 1 );
				
				FourierBackprojection::backprojectSlice_noSF(
					dataImage, ctfImage, proj,
					accumulationVolume->data,
					accumulationVolume->weight,
					accumulationVolume->multiplicity,
					padding_factor,
					num_threads_in);
						
				FourierBackprojection::backprojectSpreadingFunction(
					proj,
					accumulationVolume->spreading_function,
					padding_factor);
			}
		}
	}
}

void SpaBackproject::reconstructForward()
{
	const bool do_use_fsc = fn_fsc != "";
	bool do_MAP = do_use_fsc;
	
	Image<RFLOAT> vol_xmipp;
		
	if (verb > 0)
	{
		Log::print("Merging volumes");
	}
	
	for (int subset = 0; subset < 2; subset++)	
	for (int th = 1; th < num_threads_out; th++)
	{
		backprojectors[subset].data += backprojectors[2*th + subset].data;
		backprojectors[subset].weight += backprojectors[2*th + subset].weight;
	}
	
	
	for (int subset = 0; subset < 2; subset++)
	{
		if (verb > 0)
		{
			Log::beginSection("Subset " + ZIO::itoa(subset+1));
		}
		
		BackProjector& backprojector = backprojectors[subset];
		
		if (verb > 0)
		{
			Log::print("Applying symmetries");
		}
		
		backprojector.symmetrise(
					nr_helical_asu, helical_twist, helical_rise/angpix, 
					num_threads_total);
	
		const long int s = ctf_dim;
		
		if (do_reconstruct_ctf)
		{
			if (verb > 0)
			{
				Log::print("Dividing CTF volume");
			}
			
			vol_xmipp().initZeros(s, s, s);
			vol_xmipp().setXmippOrigin();
			
			RawImage<RFLOAT> vol(vol_xmipp);
			RawImage<Complex> backproj_data(backprojector.data);
			RawImage<RFLOAT> backproj_weight(backprojector.weight);
					
			const long int mid_out = s/2;
			const long int mid_in = backprojector.data.ydim / 2;
			
			
			for (long int z = 0; z < s; z++)
			for (long int y = 0; y < s; y++)
			for (long int x = 0; x < s; x++)
			{
				int xp, yp, zp;
	
				if (x < mid_out)
				{
					xp = mid_out - x;
					
					yp = mid_in - (y - mid_out);
					zp = mid_in - (z - mid_out);
				}
				else
				{
					xp = x - mid_out;
					
					yp = mid_in + (y - mid_out);
					zp = mid_in + (z - mid_out);
				}
				
				if ( xp >= 0 && xp < backproj_data.xdim &&
					 yp >= 0 && yp < backproj_data.ydim &&
					 zp >= 0 && zp < backproj_data.zdim)
				{
					const RFLOAT wg = backproj_weight(xp,yp,zp);
					RFLOAT val = 0.f;
					
					if (wg > 0.)
					{
						val = backproj_data(xp,yp,zp) / wg;
					}
					
					if (do_reconstruct_ctf2)
					{
						val = sqrt(val);
					}
					
					vol(x,y,z) = val;
				}
			}
		}
		else
		{
			if (do_debug)
			{
				if (verb > 0)
				{
					Log::print("Writing debugging data");
				}
				
				Image<RFLOAT> It;
				FileName fn_tmp = fn_out.withoutExtension() + "_half_" + ZIO::itoa(subset+1);
				It().resize(backprojector.data);
				
				for (long int n = 0; n < It().nzyxdim; n++)
				{
					It().data[n] = backprojector.data.data[n].real;
				}
				
				It.write(fn_tmp+"_data_real.mrc");
				
				for (long int n = 0; n < It().nzyxdim; n++)
				{
					It().data[n] = backprojector.data.data[n].imag;
				}
				
				It.write(fn_tmp+"_data_imag.mrc");
				It() = backprojector.weight;
				It.write(fn_tmp+"_weight.mrc");
			}
	
			MultidimArray<RFLOAT> tau2;
			
			if (do_use_fsc) 
			{
				if (verb > 0)
				{
					Log::print("Reading FSC");
				}
				
				do_MAP = true;
				MetaDataTable MDfsc;
				MDfsc.read(fn_fsc);
				
				MultidimArray<RFLOAT> fsc;
				fsc.resize(output_boxsize/2+1);
				
				for (int i = 0; i < MDfsc.numberOfObjects(); i++)
				{
					int idx = MDfsc.getInt(EMDL_SPECTRAL_IDX, i);
					RFLOAT val = MDfsc.getRfloat(EMDL_MLMODEL_FSC_HALVES_REF, i);
					
					fsc(idx) = val;
				}
				
				MultidimArray<RFLOAT> dummy;
				
				backprojector.updateSSNRarrays(1., tau2, dummy, dummy, dummy, fsc, do_use_fsc, true);
			}
	
			if (do_external_reconstruct)
			{
				if (verb > 0)
				{
					Log::print("Preparing for external reconstruction");
				}
				
				FileName fn_root = fn_out.withoutExtension() + "_half_" + ZIO::itoa(subset+1);			
				MultidimArray<RFLOAT> dummy;
				
				backprojector.externalReconstruct(vol_xmipp(),
						fn_root,
						tau2, dummy, dummy, dummy, false, 1., 1);
			}
			else
			{
				if (verb > 0)
				{
					Log::print("Reconstructing");
				}
				
				backprojector.reconstruct(vol_xmipp(), iter, do_MAP, tau2);
			}
		}
		
		std::string my_fn_out = fn_out + "_half_" + ZIO::itoa(subset+1) + ".mrc";
		
		vol_xmipp.setSamplingRateInHeader(angpix);
		vol_xmipp.write(my_fn_out);
		
		if (verb > 0)
		{
			Log::print("Done! Written output map in: " + my_fn_out);
		}
		
		Log::endSection();
	}
}

void SpaBackproject::reconstructBackward()
{
	const int outCount = accumulation_volumes.size();
			
	if (outCount > 2)
	{		
		Log::print("Merging volumes");
	
		for (int i = 2; i < outCount; i++)
		{
			accumulation_volumes[i%2].data += accumulation_volumes[i].data;
			accumulation_volumes[i%2].weight += accumulation_volumes[i].weight;
			accumulation_volumes[i%2].multiplicity += accumulation_volumes[i].multiplicity;
			accumulation_volumes[i%2].spreading_function += accumulation_volumes[i].spreading_function;
		}
	}
	
	std::vector<BufferedImage<Complex>> dataImgFS = {
		accumulation_volumes[0].data, 
		accumulation_volumes[1].data};
	
	std::vector<BufferedImage<RFLOAT>> psfImgFS = {
		accumulation_volumes[0].spreading_function, 
		accumulation_volumes[1].spreading_function};
	
	std::vector<BufferedImage<RFLOAT>> ctfImgFS = {
		accumulation_volumes[0].weight, 
		accumulation_volumes[1].weight};
	
	if (fn_sym != "C1")
	{		
		Log::print("Applying symmetries");
		
		for (int half = 0; half < 2; half++)
		{
			dataImgFS[half] = Symmetry::symmetrise_FS_complex(
				dataImgFS[half], fn_sym, num_threads_total);
			
			psfImgFS[half] = Symmetry::symmetrise_FS_real(
				psfImgFS[half], fn_sym, num_threads_total);
			
			ctfImgFS[half] = Symmetry::symmetrise_FS_real(
				ctfImgFS[half], fn_sym, num_threads_total);
		}
	}
	
	const int s = dataImgFS[0].ydim;	
	
	std::vector<BufferedImage<double>> dataImgRS(2), dataImgDivRS(2);
	
	BufferedImage<dComplex> dataImgFS_both = dataImgFS[0] + dataImgFS[1];
	BufferedImage<double> psfImgFS_both = psfImgFS[0] + psfImgFS[1];
	BufferedImage<double> ctfImgFS_both = ctfImgFS[0] + ctfImgFS[1];
	
	Log::beginSection("Reconstructing");
	
	const double WienerFract = 0.01;
	
	for (int half = 0; half < 2; half++)
	{		
		Log::print("Half " + ZIO::itoa(half));
		
		dataImgRS[half] = BufferedImage<double>(s,s,s);
		dataImgDivRS[half] = BufferedImage<double>(s,s,s);
		
		Reconstruction::griddingCorrect3D(
					dataImgFS[half], psfImgFS[half], dataImgRS[half],
					true, num_threads_total);
		
		Reconstruction::ctfCorrect3D(
					dataImgRS[half], ctfImgFS[half], dataImgDivRS[half],
					1.0 / WienerFract, num_threads_total);
		
		dataImgDivRS[half].write(fn_out+"_half"+ZIO::itoa(half+1)+".mrc", angpix);		
		dataImgRS[half].write(fn_out+"_data_half"+ZIO::itoa(half+1)+".mrc", angpix);
		
		Centering::fftwHalfToHumanFull(ctfImgFS[half]).write(
					fn_out+"_weight_half"+ZIO::itoa(half+1)+".mrc", angpix);
	}
	
	Log::endSection();
	
	Reconstruction::griddingCorrect3D(
		dataImgFS_both, psfImgFS_both, dataImgRS[0], true, num_threads_total);
	
	Reconstruction::ctfCorrect3D(
		dataImgRS[0], ctfImgFS_both, dataImgDivRS[0], 1.0 / WienerFract, num_threads_total);
	
	dataImgDivRS[0].write(fn_out+"_merged.mrc", angpix);	
	dataImgRS[0].write(fn_out+"_data_merged.mrc", angpix);
	
	Centering::fftwHalfToHumanFull(ctfImgFS[0]).write(fn_out+"_weight_merged.mrc", angpix);
}

void SpaBackproject::applyCTFPandCTFQ(
		MultidimArray<Complex> &Fin, CTF &ctf, FourierTransformer &transformer,
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
				if (newbox > 0 && newbox < YSIZE(Fin))
				{
					Iapp.setXmippOrigin();
					Iapp.window(FIRST_XMIPP_INDEX(newbox), FIRST_XMIPP_INDEX(newbox),
					            LAST_XMIPP_INDEX(newbox),  LAST_XMIPP_INDEX(newbox));
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
					{
						DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					}
					else if (myangle < anglemax)
					{
						DIRECT_A2D_ELEM(*myCTFPorQb, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					}
				}
				else
				{
					if (myangle >= anglemin && myangle < anglemax)
					{
						DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
					}
				}
			}
		}
	}
}
