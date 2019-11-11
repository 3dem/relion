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

#include "src/particle_subtractor.h"

void ParticleSubtractor::read(int argc, char **argv)
{
	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	fn_opt = parser.getOption("--i", "Name of optimiser.star file from refinement/classification to use for subtraction", "");
	fn_out = parser.getOption("--o", "Output directory name", "Subtract/");
	fn_msk = parser.getOption("--mask", "Name of the 3D mask with all density that should be kept, i.e. not subtracted", "");
	fn_sel = parser.getOption("--data", "Name of particle STAR file, in case not all particles from optimiser are to be used", "");
	fn_revert = parser.getOption("--revert", "Name of particle STAR file to revert. When this is provided, all other options are ignored.", "");

	int center_section = parser.addSection("Centering options");
	do_recenter_on_mask = parser.checkOption("--recenter_on_mask", "Use this flag to center the subtracted particles on projections of the centre-of-mass of the input mask");
	new_center.resize(3);
	XX(new_center) = textToInteger(parser.getOption("--center_x", "X-coordinate of 3D coordinate, which will be projected to center the subtracted particles.", "9999"));
	YY(new_center) = textToInteger(parser.getOption("--center_y", "Y-coordinate of 3D coordinate, which will be projected to center the subtracted particles.", "9999"));
	ZZ(new_center) = textToInteger(parser.getOption("--center_z", "Z-coordinate of 3D coordinate, which will be projected to center the subtracted particles.", "9999"));
	boxsize = textToInteger(parser.getOption("--new_box", "Output size of the subtracted particles", "-1"));

	verb = 1;
	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	if ((fn_opt == "" && fn_revert == "") || (fn_opt != "" && fn_revert != ""))
		REPORT_ERROR("Please specify only one of --i OR --revert");
}

void ParticleSubtractor::usage()
{
	parser.writeUsage(std::cout);
}

void ParticleSubtractor::divideLabour(int _rank, int _size, long int &my_first, long int &my_last)
{
	if (opt.do_split_random_halves)
	{
		if (size < 2) REPORT_ERROR("You need to run this program with at least 2 MPI processes for subtraction of 3D auto-refine jobs");
		opt.mydata.divideParticlesInRandomHalves(0, opt.do_helical_refine);

		int my_halfset = (_rank % 2 == 1) ? 1 : 2;
		int mysize = (my_halfset == 1) ? _size/2 : _size/2 + _size%2;
		divide_equally(opt.mydata.numberOfParticles(my_halfset), mysize, _rank/2, my_first, my_last);
		if (my_halfset == 2)
		{
    		my_first += opt.mydata.numberOfParticles(1);
    		my_last += opt.mydata.numberOfParticles(1);
		}
	}
	else
	{
		divide_equally(opt.mydata.numberOfParticles(), _size, _rank, my_first, my_last);
	}
	// std::cerr << "_size= " <<_size << " _rank= " << _rank << " my_first= " << my_first << " my_last= " << my_last << std::endl;
}

// Initialise some stuff after reading
void ParticleSubtractor::initialise(int _rank, int _size)
{
	rank = _rank;
	size = _size;

	if (rank > 0) verb = 0;

	if (fn_revert != "")
		return;

	// Make directory for output particles
	if (fn_out[fn_out.length()-1] != '/') fn_out += "/";
	if (verb > 0)
	{
		int res = system(("mkdir -p " + fn_out + "Particles").c_str());
	}

	opt.read(fn_opt, rank, true); // true means: prevent prereading all particle images

	// Overwrite the particles STAR file with a smaller subset
	if (fn_sel != "")
	{
		opt.mydata.clear();
		bool is_helical_segment = (opt.do_helical_refine) || ((opt.mymodel.ref_dim == 2) && (opt.helical_tube_outer_diameter > 0.));
		opt.mydata.read(fn_sel, false, false, false, is_helical_segment);
	}

	divideLabour(rank, size, my_first_part_id, my_last_part_id);

	if (verb > 0) std::cout << " + Reading in mask ... " << std::endl;
	// Mask stuff
	Image<RFLOAT> Imask;
	Imask.read(fn_msk);
	Imask().setXmippOrigin();

	RFLOAT minval, maxval;
	Imask().computeDoubleMinMax(minval, maxval);
	if (minval < 0. || maxval > 1.)
	{
		REPORT_ERROR("ERROR: the keep_inside mask has values outside the range [0,1]");
	}

	if (do_recenter_on_mask)
	{
		Imask().centerOfMass(new_center);
		do_center = true;
	}
	else if ((int)XX(new_center) != 9999 && (int)YY(new_center) != 9999 && (int)ZZ(new_center) != 9999)
	{
		do_center = true;
	}
	else
	{
		new_center.initZeros();
		do_center = false;
	}

	if (verb > 0 && (do_center || opt.fn_body_masks != "None"))
	{
		std::cout << " + The subtracted particles will be re-centred on projections of 3D-coordinate: (" << new_center(0) << " , " << new_center(1) << " , " << new_center(2) << ")" << std::endl;
	}

	if (opt.fn_body_masks != "None")
	{
		if (verb > 0)
		{
			std::cout << " + Initialising multi-body subtraction ..." << std::endl;
		}

		// This creates a rotation matrix for (rot,tilt,psi) = (0,90,0)
		// It will be used to make all Abody orientation matrices relative to (0,90,0) instead of the more logical (0,0,0)
		// This is useful, as psi-priors are ill-defined around tilt=0, as rot becomes the same as -psi!!
		rotation3DMatrix(-90., 'Y', A_rot90, false);
		A_rot90T = A_rot90.transpose();

		// Find out which body has the biggest overlap with the keepmask, use these orientations
		RFLOAT best_overlap = 0.;
		subtract_body = -1;
		for (int ibody = 0; ibody < opt.mymodel.nr_bodies; ibody++)
		{
			if (!Imask().sameShape(opt.mymodel.masks_bodies[ibody]))
			{
				Imask().printShape();
				opt.mymodel.masks_bodies[ibody].printShape();
				REPORT_ERROR("ERROR: input mask is not of same shape as body masks.");
			}

			RFLOAT overlap = 0.;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imask())
				overlap += DIRECT_MULTIDIM_ELEM(opt.mymodel.masks_bodies[ibody], n) * DIRECT_MULTIDIM_ELEM(Imask(), n);

			if (overlap > best_overlap)
			{
				best_overlap = overlap;
				subtract_body = ibody;
			}
		}

		if (subtract_body < 0) REPORT_ERROR("ERROR: input mask does not overlap with any of the bodies....");

		// Apply the inverse of the keepmask to all the mask_bodies
		for (int obody = 0; obody < opt.mymodel.nr_bodies; obody++)
		{
			int ii = DIRECT_A2D_ELEM(opt.mymodel.pointer_body_overlap, subtract_body, obody);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imask())
			{
				DIRECT_MULTIDIM_ELEM(opt.mymodel.masks_bodies[ii], n) *= (1. - DIRECT_MULTIDIM_ELEM(Imask(), n));
			}
		}
	}
	else
	{
		// For normal refinement/classification: just apply the inverse of keepmask to the references
		for (int iclass = 0; iclass < opt.mymodel.nr_classes; iclass++)
		{
			if (!Imask().sameShape(opt.mymodel.Iref[iclass]))
			{
				Imask().printShape();
				opt.mymodel.Iref[iclass].printShape();
				REPORT_ERROR("ERROR: input mask is not of same shape as reference inside the optimiser.");
			}

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imask())
			{
				DIRECT_MULTIDIM_ELEM(opt.mymodel.Iref[iclass], n) *= (1. - DIRECT_MULTIDIM_ELEM(Imask(), n));
			}
		}
	}

	if (verb > 0)
	{
		std::cout << " + Calculating Fourier transforms of the maps ..." << std::endl;
	}

	// Now set up the Projectors inside the model
	opt.mymodel.setFourierTransformMaps(false); // false means ignore tau2_class

	// ensure even boxsize of subtracted images
	if (boxsize > 0)
		boxsize -= boxsize%2;
}


void ParticleSubtractor::revert()
{
	ObservationModel obsModel;
	MetaDataTable MD;

	ObservationModel::loadSafely(fn_revert, obsModel, MD);

	if (!MD.containsLabel(EMDL_IMAGE_ORI_NAME))
		REPORT_ERROR("The input STAR file does not contain the rlnImageOriginalName column.");

	if (!MD.containsLabel(EMDL_IMAGE_NAME))
		REPORT_ERROR("The input STAR file does not contain the rlnImageName column");

	// Swap image names
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		FileName f1, f2;
		MD.getValue(EMDL_IMAGE_ORI_NAME, f1);
		MD.getValue(EMDL_IMAGE_NAME, f2);
		MD.setValue(EMDL_IMAGE_ORI_NAME, f2);
		MD.setValue(EMDL_IMAGE_NAME, f1);
	}

	// Fix box size
	// TODO: BUG: What happens if the input particles have several box sizes?
	std::vector<bool> fixed_box_size(obsModel.numberOfOpticsGroups(), false);
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		const int og = obsModel.getOpticsGroup(MD);
		if (fixed_box_size[og])
			continue;

		FileName img_name, fn_img;
		long int dummy;
		MD.getValue(EMDL_IMAGE_NAME, img_name);
		img_name.decompose(dummy, fn_img);

		if (!exists(fn_img))
			REPORT_ERROR("Failed to read " + fn_img + " to determine the box size.");
		Image<RFLOAT> Ihead;
		Ihead.read(img_name, false, -1, false, true);
		if (XSIZE(Ihead()) != YSIZE(Ihead()))
			REPORT_ERROR("Particle " + img_name + " is not square.");
		obsModel.setBoxSize(og, XSIZE(Ihead()));
		obsModel.opticsMdt.setValue(EMDL_IMAGE_SIZE, XSIZE(Ihead()), og);

		fixed_box_size[og] = true;
	}

	for (int i = 0; i < obsModel.numberOfOpticsGroups(); i++)
	{
		if (!fixed_box_size[i])
			std::cerr << "WARNING: could not determine the box size of optics group " << obsModel.getGroupName(i) << std::endl;
		else
			std::cout << "The box size of the optics group " << obsModel.getGroupName(i) << " was updated to " << obsModel.getBoxSize(i) << " px" << std::endl;
	}

	obsModel.save(MD, fn_out + "original.star");
	std::cout << "Writen " << (fn_out + "original.star") << std::endl;
}

void ParticleSubtractor::run()
{
	int my_halfset = 0;

	long int nr_parts = my_last_part_id - my_first_part_id + 1;
	long int barstep = XMIPP_MAX(1, nr_parts/120);
	if (verb > 0)
	{
		std::cout << " + Subtracting all particles ..." << std::endl;
		time_config();
		init_progress_bar(nr_parts);
	}

	long int imgno = 0;
	// For multibody: store, xoff, yoff, zoff, rot, tilt and psi; for normal, only store xoff, yoff, zoff
	if (opt.fn_body_masks != "None")
	{
		orients.resize(nr_parts, 6);
	}
	else if (do_center)
	{
		orients.resize(nr_parts, 3);
	}

	for (long int part_id = my_first_part_id, imgno = 0; part_id <= my_last_part_id; part_id++, imgno++)
	{

		if (imgno % barstep == 0)
		{
			if (pipeline_control_check_abort_job())
				exit(RELION_EXIT_ABORTED);
		}

		subtractOneParticle(part_id, imgno, orients);

		if (imgno % barstep == 0 && verb > 0) progress_bar(imgno);
	}

	if (verb > 0) progress_bar(nr_parts);
}

void ParticleSubtractor::setLinesInStarFile(int myrank)
{
	long int my_first, my_last;
	divideLabour(myrank, size, my_first, my_last);

	long int imgno = 0;
	for (long int part_id = my_first; part_id <= my_last; part_id++)
	{
		if (do_center || opt.fn_body_masks != "None")
		{
			opt.mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, DIRECT_A2D_ELEM(orients, imgno, 0), part_id);
			opt.mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, DIRECT_A2D_ELEM(orients, imgno, 1), part_id);
			if (opt.mymodel.data_dim == 3) opt.mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, DIRECT_A2D_ELEM(orients, imgno, 2), part_id);
		}
		if (opt.fn_body_masks != "None")
		{
			opt.mydata.MDimg.setValue(EMDL_ORIENT_ROT, DIRECT_A2D_ELEM(orients, imgno, 3), part_id);
			opt.mydata.MDimg.setValue(EMDL_ORIENT_TILT, DIRECT_A2D_ELEM(orients, imgno, 4), part_id);
			opt.mydata.MDimg.setValue(EMDL_ORIENT_PSI, DIRECT_A2D_ELEM(orients, imgno, 5), part_id);
		}

		// Store the original particle name, and also set the subtracted name
		FileName fn_img;
		opt.mydata.MDimg.getValue(EMDL_IMAGE_NAME, fn_img, part_id);
		opt.mydata.MDimg.setValue(EMDL_IMAGE_ORI_NAME, fn_img, part_id);
		fn_img = getParticleName(imgno, myrank);
		opt.mydata.MDimg.setValue(EMDL_IMAGE_NAME, fn_img, part_id);

		imgno++;
	}

	// Remove origin prior columns if present, as we have re-centered.
	if (do_center || opt.fn_body_masks != "None")
	{
		if (opt.mydata.MDimg.containsLabel(EMDL_ORIENT_ORIGIN_X_PRIOR_ANGSTROM))
		{
			opt.mydata.MDimg.deactivateLabel(EMDL_ORIENT_ORIGIN_X_PRIOR_ANGSTROM);
		}
		if (opt.mydata.MDimg.containsLabel(EMDL_ORIENT_ORIGIN_Y_PRIOR_ANGSTROM))
		{
			opt.mydata.MDimg.deactivateLabel(EMDL_ORIENT_ORIGIN_Y_PRIOR_ANGSTROM);
		}
	}
}

void ParticleSubtractor::saveStarFile()
{
	// Reset image size in optics table, if the images were rewindowed in a different box
	if (boxsize > 0)
	{
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(opt.mydata.obsModel.opticsMdt)
		{
			opt.mydata.obsModel.opticsMdt.setValue(EMDL_IMAGE_SIZE, boxsize);
		}
	}

	opt.mydata.obsModel.save(opt.mydata.MDimg, fn_out + "particles_subtracted.star");
	std::cout << " + Saved STAR file with " << opt.mydata.MDimg.numberOfObjects()
	          << " subtracted particles in " << fn_out <<"particles_subtracted.star" << std::endl;
}

FileName ParticleSubtractor::getParticleName(long int imgno, int myrank)
{
	// Now write out the image
	FileName fn_img;
	if (opt.mymodel.data_dim == 3)
	{
		fn_img.compose(fn_out+"Particles/subtracted", imgno+1,"mrc");
	}
	else
	{
		FileName fn_stack;
		if (size > 1)  fn_stack.compose(fn_out+"Particles/subtracted_", myrank + 1, "mrcs");
		else fn_stack = fn_out+"Particles/subtracted.mrcs";
		fn_img.compose(imgno+1,fn_stack);
	}

	return fn_img;
}

void ParticleSubtractor::subtractOneParticle(long int part_id, long int imgno, MultidimArray<RFLOAT> &orients)
{
	// Read the particle image
	Image<RFLOAT> img;
	FileName fn_img;
	opt.mydata.MDimg.getValue(EMDL_IMAGE_NAME, fn_img, part_id);
	img.read(fn_img);
	img().setXmippOrigin();
	int optics_group = opt.mydata.getOpticsGroup(part_id, 0);

	// Get the consensus class, orientational parameters and norm (if present)
	RFLOAT my_pixel_size = opt.mydata.getImagePixelSize(part_id, 0);
	Matrix1D<RFLOAT> my_old_offset(3), my_residual_offset(3), centering_offset(3);
	Matrix2D<RFLOAT> Aori;
	RFLOAT rot, tilt, psi, xoff, yoff, zoff, norm, scale;
	int myclass=0;
	if (opt.mydata.MDimg.containsLabel(EMDL_PARTICLE_CLASS))
	{
		opt.mydata.MDimg.getValue(EMDL_PARTICLE_CLASS, myclass, part_id);
		myclass--; // start counting at zero instead of one
	}
	opt.mydata.MDimg.getValue(EMDL_ORIENT_ROT, rot, part_id);
	opt.mydata.MDimg.getValue(EMDL_ORIENT_TILT, tilt, part_id);
	opt.mydata.MDimg.getValue(EMDL_ORIENT_PSI, psi, part_id);
	opt.mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, XX(my_old_offset), part_id);
	opt.mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, YY(my_old_offset), part_id);
	if (opt.mymodel.data_dim == 3) opt.mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, ZZ(my_old_offset), part_id);
	// As of v3.1, offsets are in Angstrom: convert back to pixels!
	my_old_offset /= my_pixel_size;

	// Apply the norm_correction term
	if (!opt.mydata.MDimg.getValue(EMDL_IMAGE_NORM_CORRECTION, norm, part_id)) norm = 1.;
	if (opt.do_norm_correction) img() *= opt.mymodel.avg_norm_correction / norm;

	Matrix1D<RFLOAT> my_projected_com(3), my_refined_ibody_offset(3);
	if (opt.fn_body_masks != "None")
	{
		// 17May2017: Shift image to the projected COM for this body!
		// Aori is the original transformation matrix of the consensus refinement
		Euler_angles2matrix(rot, tilt, psi, Aori, false);
		my_projected_com = Aori * opt.mymodel.com_bodies[subtract_body];

		// Subtract the projected COM offset, to position this body in the center
		my_old_offset -= my_projected_com;
		my_residual_offset = my_old_offset;
		// Apply the old_offset (rounded to avoid interpolation errors)
		my_old_offset.selfROUND();
		selfTranslate(img(), my_old_offset, WRAP);
		// keep track of the differences between the rounded and the original offsets
		my_residual_offset -= my_old_offset;

	}

	// Now that the particle is centered (for multibody), get the FourierTransform of the particle
	MultidimArray<Complex> Faux, Fimg;
	MultidimArray<RFLOAT> Fctf;
	FourierTransformer transformer;
	CenterFFT(img(), true);
	transformer.FourierTransform(img(), Fimg);
	Fctf.resize(Fimg);

	if (opt.do_ctf_correction)
	{
		if (opt.mymodel.data_dim == 3)
		{
			Image<RFLOAT> Ictf;
			FileName fn_ctf;
			opt.mydata.MDimg.getValue(EMDL_CTF_IMAGE, fn_ctf, part_id);
			Ictf.read(fn_ctf);

			// If there is a redundant half, get rid of it
			if (XSIZE(Ictf()) == YSIZE(Ictf()))
			{
				Ictf().setXmippOrigin();
				FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf)
				{
					// Use negative kp,ip and jp indices, because the origin in the ctf_img lies half a pixel to the right of the actual center....
					DIRECT_A3D_ELEM(Fctf, k, i, j) = A3D_ELEM(Ictf(), -kp, -ip, -jp);
				}
			}
			// otherwise, just window the CTF to the current resolution
			else if (XSIZE(Ictf()) == YSIZE(Ictf()) / 2 + 1)
			{
				windowFourierTransform(Ictf(), Fctf, YSIZE(Fctf));
			}
			// if dimensions are neither cubical nor FFTW, stop
			else
			{
				REPORT_ERROR("3D CTF volume must be either cubical or adhere to FFTW format!");
			}
		}
		else
		{
			CTF ctf;
			ctf.readByGroup(opt.mydata.MDimg, &opt.mydata.obsModel, part_id);
			ctf.getFftwImage(Fctf, opt.mymodel.ori_size, opt.mymodel.ori_size, opt.mymodel.pixel_size,
					opt.ctf_phase_flipped, false, opt.intact_ctf_first_peak, true);
		}
	}
	else
	{
		Fctf.initConstant(1.);
	}

	MultidimArray<Complex> Fsubtract;
	Fsubtract.initZeros(Fimg);

	if (opt.fn_body_masks != "None")
	{
		// For multi-body refinement
		Matrix2D<RFLOAT> Aresi_subtract;
		for (int obody = 0; obody < opt.mymodel.nr_bodies; obody++)
		{
			// Unlike getFourierTransformsAndCtfs, no check for ibody==obody: also subtract rest of subtract_body!

			Matrix1D<RFLOAT> body_offset(3);
			RFLOAT body_rot, body_tilt, body_psi;
			opt.mydata.MDbodies[obody].getValue(EMDL_ORIENT_ROT, body_rot, part_id);
			opt.mydata.MDbodies[obody].getValue(EMDL_ORIENT_TILT, body_tilt, part_id);
			opt.mydata.MDbodies[obody].getValue(EMDL_ORIENT_PSI,  body_psi, part_id);
			opt.mydata.MDbodies[obody].getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, XX(body_offset), part_id);
			opt.mydata.MDbodies[obody].getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, YY(body_offset), part_id);
			if (opt.mymodel.data_dim == 3)
			{
				opt.mydata.MDbodies[obody].getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, ZZ(body_offset), part_id);
			}

			// As of v3.1, offsets are in Angstrom: convert back to pixels!
			body_offset /= my_pixel_size;

			Matrix2D<RFLOAT> Aresi,  Abody;
			// Aresi is the residual orientation for this obody
			Euler_angles2matrix(body_rot, body_tilt, body_psi, Aresi, false);
			if (obody == subtract_body) Aresi_subtract = Aresi;
			// The real orientation to be applied is the obody transformation applied and the original one
			Abody = Aori * (opt.mymodel.orient_bodies[obody]).transpose() * A_rot90 * Aresi * opt.mymodel.orient_bodies[obody];
			// Apply anisotropic mag and scaling
			Abody = opt.mydata.obsModel.applyAnisoMag(Abody, optics_group);
			Abody = opt.mydata.obsModel.applyScaleDifference(Abody, optics_group, opt.mymodel.ori_size, opt.mymodel.pixel_size);

			// Get the FT of the projection in the right direction
			MultidimArray<Complex> FTo;
			FTo.initZeros(Fimg);
			// The following line gets the correct pointer to account for overlap in the bodies
			int oobody = DIRECT_A2D_ELEM(opt.mymodel.pointer_body_overlap, subtract_body, obody);
			opt.mymodel.PPref[oobody].get2DFourierTransform(FTo, Abody);

			// Body is centered at its own COM: move it back to its place in the original particle image

			// Projected COM for this body (using Aori, just like above for ibody and my_projected_com!!!)
			Matrix1D<RFLOAT> other_projected_com(3);
			other_projected_com = Aori * (opt.mymodel.com_bodies[obody]);

			// Subtract refined obody-displacement
			other_projected_com -= body_offset;

			// Subtract the projected COM already applied to this image for ibody
			other_projected_com -= my_projected_com;

			shiftImageInFourierTransform(FTo, Faux, (RFLOAT)opt.mymodel.ori_size,
					XX(other_projected_com), YY(other_projected_com), ZZ(other_projected_com));

			// Sum the Fourier transforms of all the obodies
			Fsubtract += Faux;
		} // end for obody

		// Set orientations back into the original RELION system of coordinates
		Matrix2D<RFLOAT> Abody;

		// Write out the rot,tilt,psi as the combination of Aori and Aresi!! So get rid of the rotations around the tilt=90 axes,
		Abody = Aori * (opt.mymodel.orient_bodies[subtract_body]).transpose() * A_rot90 * Aresi_subtract * opt.mymodel.orient_bodies[subtract_body];
		Euler_matrix2angles(Abody, rot, tilt, psi);

		// Store the optimal orientations in the orients array
		DIRECT_A2D_ELEM(orients, imgno, 3) = rot;
		DIRECT_A2D_ELEM(orients, imgno, 4) = tilt;
		DIRECT_A2D_ELEM(orients, imgno, 5) = psi;

		// Also get refined offset for this body
		opt.mydata.MDbodies[subtract_body].getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, XX(my_refined_ibody_offset), part_id);
		opt.mydata.MDbodies[subtract_body].getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, YY(my_refined_ibody_offset), part_id);
		if (opt.mymodel.data_dim == 3)
		{
			opt.mydata.MDbodies[subtract_body].getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, ZZ(my_refined_ibody_offset), part_id);
		}
		// As of v3.1, offsets are in Angstrom: convert back to pixels!
		my_refined_ibody_offset /= my_pixel_size;

		// re-center to new_center
		my_residual_offset += my_refined_ibody_offset;
		my_residual_offset += Abody * (opt.mymodel.com_bodies[subtract_body] - new_center);

	}
	else
	{
		// Normal 3D classification/refinement: get the projection in rot,tilt,psi for the corresponding class
		Matrix2D<RFLOAT> A3D;
		Euler_angles2matrix(rot, tilt, psi, A3D);

		// Apply anisotropic mag and scaling
		A3D = opt.mydata.obsModel.applyAnisoMag(A3D, optics_group);
		A3D = opt.mydata.obsModel.applyScaleDifference(A3D, optics_group, opt.mymodel.ori_size, opt.mymodel.pixel_size);
		opt.mymodel.PPref[myclass].get2DFourierTransform(Fsubtract, A3D);

		// Shift in opposite direction as offsets in the STAR file
		shiftImageInFourierTransform(Fsubtract, Fsubtract, (RFLOAT)opt.mymodel.ori_size,
				-XX(my_old_offset), -YY(my_old_offset), -ZZ(my_old_offset));

		if (do_center)
		{
			// Re-center the output particle to a new centre...
			my_residual_offset = my_old_offset - A3D * (new_center);
		}
	}


	// Apply the CTF to the to-be-subtracted projection
	if (opt.do_ctf_correction)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fsubtract)
		{
			DIRECT_MULTIDIM_ELEM(Fsubtract, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
		}

		// Also do phase modulation, for beam tilt correction and other asymmetric aberrations
		opt.mydata.obsModel.demodulatePhase(optics_group, Fsubtract, true); // true means do_modulate_instead
		opt.mydata.obsModel.divideByMtf(optics_group, Fsubtract, true); // true means do_multiply_instead
	}

	// Do the actual subtraction
	Fimg -= Fsubtract;

	// And go finally back to real-space
	transformer.inverseFourierTransform(Fimg, img());
	CenterFFT(img(), false);

	if (do_center || opt.fn_body_masks != "None")
	{
		// Recenter the particles
		centering_offset = my_residual_offset;
		centering_offset.selfROUND();
		my_residual_offset -= centering_offset;
		selfTranslate(img(), centering_offset, WRAP);

		// Set the non-integer difference between the rounded centering offset and the actual offsets in the STAR file
		DIRECT_A2D_ELEM(orients, imgno, 0) = my_pixel_size * XX(my_residual_offset);
		DIRECT_A2D_ELEM(orients, imgno, 1) = my_pixel_size * YY(my_residual_offset);
		if (opt.mymodel.data_dim == 3) DIRECT_A2D_ELEM(orients, imgno, 2) = my_pixel_size * ZZ(my_residual_offset);
	}

	// Rebox the image
	if (boxsize > 0)
	{
		if (img().getDim() == 2)
		{
			img().window(FIRST_XMIPP_INDEX(boxsize), FIRST_XMIPP_INDEX(boxsize),
					   LAST_XMIPP_INDEX(boxsize),  LAST_XMIPP_INDEX(boxsize));
		}
		else if (img().getDim() == 3)
		{
			img().window(FIRST_XMIPP_INDEX(boxsize), FIRST_XMIPP_INDEX(boxsize), FIRST_XMIPP_INDEX(boxsize),
					   LAST_XMIPP_INDEX(boxsize),  LAST_XMIPP_INDEX(boxsize),  LAST_XMIPP_INDEX(boxsize));
		}
	}

	// Now write out the image
	fn_img = getParticleName(imgno, rank);
	img.setSamplingRateInHeader(my_pixel_size);
	if (opt.mymodel.data_dim == 3)
	{
		img.write(fn_img);
	}
	else
	{
		if (imgno == 0)
			img.write(fn_img, -1, false, WRITE_OVERWRITE);
		else
			img.write(fn_img, -1, false, WRITE_APPEND);
	}
}
