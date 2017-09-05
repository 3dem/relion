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

#include "flex_analyser.h"


void FlexAnalyser::read(int argc, char **argv)
{
	parser.setCommandLine(argc, argv);

	int gen_section = parser.addSection("General options");
	fn_data = parser.getOption("--data", "The _data.star file with the orientations to be analysed", "");
	fn_model = parser.getOption("--model", " The corresponding _model.star file with the refined model", "");
	fn_bodies = parser.getOption("--bodies", "The corresponding star file with the definition of the bodies", "");
	fn_out = parser.getOption("--o", "Output rootname", "analyse");

	int model_section = parser.addSection("3D model options");
	do_3dmodels = parser.checkOption("--3dmodels", "Generate a 3D model for each experimental particles");
	size_3dmodels = textToInteger(parser.getOption("--size_3dmodels", "Output size of the 3D models (default is same as input particles)", "-1"));

	int pca_section = parser.addSection("PCA options");
	do_PCA_3dmodels = parser.checkOption("--PCA_3dmodels", "Perform a principal components analysis on input 3D models");
	do_PCA_orient = parser.checkOption("--PCA_orient", "Perform a principal components analysis on the multibody orientations");

	int subtract_section = parser.addSection("Subtract options");
	do_subtract = parser.checkOption("--subtract", "Generate subtracted experimental particles");
	fn_keepmask = parser.getOption("--keep_inside", "Subtract everything except the density inside this mask", "");
	boxsize = textToInteger(parser.getOption("--boxsize", "Boxsize (in pixels) of the subtracted particles (default is keep original)", "-1"));

	// TODO: implement scale correction!
	//do_scale = parser.checkOption("--scale", "Apply scale correction in the subtraction");
	do_norm = parser.checkOption("--norm", "Apply normalisation correction in the subtraction");

	int ctf_section = parser.addSection("CTF options");
	do_ctf = parser.checkOption("--ctf", "Apply CTF to reference projections");
	ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Have the data been CTF phase-flipped?");
	ctf_premultiplied = parser.checkOption("--ctf_multiplied", "Have the data been premultiplied with their CTF?");
	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Ignore CTFs until their first peak?");

	// Initialise verb for non-parallel execution
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void FlexAnalyser::initialise()
{
	if (do_3dmodels || do_subtract)
	{
		if (fn_data == "")
			REPORT_ERROR("ERROR: please provide the --data argument!");
		else
			data.read(fn_data);

		if (fn_model == "")
			REPORT_ERROR("ERROR: please provide the --model argument!");
		else
			model.read(fn_model);

		if (fn_bodies != "")
			model.initialiseBodies(fn_bodies, fn_out);
		else
			model.nr_bodies = 1;

		if (model.nr_bodies != data.nr_bodies)
			REPORT_ERROR("ERROR: Unequal number of bodies in bodies.star and data.star files!");

		if (do_3dmodels && model.nr_bodies == 1)
			REPORT_ERROR("ERROR: --3dmodels option is only valid for multibody refinements.");

		if (model.nr_bodies > 1)
		{
			// This creates a rotation matrix for (rot,tilt,psi) = (0,90,0)
			// It will be used to make all Abody orientation matrices relative to (0,90,0) instead of the more logical (0,0,0)
			// This is useful, as psi-priors are ill-defined around tilt=0, as rot becomes the same as -psi!!
			rotation3DMatrix(-90., 'Y', A_rot90, false);
			A_rot90T = A_rot90.transpose();
		}

		// ensure even boxsize of subtracted images
		if (boxsize > 0)
			boxsize -= boxsize%2;

	}


	// TODO setup PCA stuff as well.... for 3D models and/or multibody orientations ...

}

void FlexAnalyser::run()
{

	if (do_subtract)
		setupSubtractionMasksAndProjectors();
	else if (do_3dmodels)
		setup3DModels();

	// Loop through all particles
	loopThroughParticles();

}

void FlexAnalyser::setupSubtractionMasksAndProjectors()
{
	Image<RFLOAT> Imask;
	Imask.read(fn_keepmask);
	Imask().setXmippOrigin();

	RFLOAT minval, maxval;
	Imask().computeDoubleMinMax(minval, maxval);
	if (minval < 0. || maxval > 1.)
		REPORT_ERROR("ERROR: the keep_inside mask has values outside the range [0,1]");

	if (model.nr_bodies > 1)
	{
		// Find out which body has the biggest overlap with the keepmask, use these orientations
		RFLOAT best_overlap = 0.;
		subtract_body = -1;
		for (int ibody = 0; ibody < model.nr_bodies; ibody++)
		{
			if (!Imask().sameShape(model.masks_bodies[ibody]))
			{
				Imask().printShape();
				model.masks_bodies[ibody].printShape();
				REPORT_ERROR("ERROR: input keep_inside mask is not of same shape as body masks.");
			}

			RFLOAT overlap = 0.;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imask())
				overlap += DIRECT_MULTIDIM_ELEM(model.masks_bodies[ibody], n) * DIRECT_MULTIDIM_ELEM(Imask(), n);

			if (overlap > best_overlap)
			{
				best_overlap = overlap;
				subtract_body = ibody;
			}
		}

		if (subtract_body < 0)
			REPORT_ERROR("ERROR: the input keep_inside mask does not overlap with any of the bodies....");

		// Apply the inverse of the keepmask to all the mask_bodies
		for (int obody = 0; obody < model.nr_bodies; obody++)
		{
			int ii = DIRECT_A2D_ELEM(model.pointer_body_overlap, subtract_body, obody);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imask())
			{
				DIRECT_MULTIDIM_ELEM(model.masks_bodies[ii], n) *= (1. - DIRECT_MULTIDIM_ELEM(Imask(), n));
			}
		}
	}
	else
	{

		subtract_body = 0;
		// Apply mask to all classes (could be more than one)
		for (int iclass = 0; iclass < model.nr_classes; iclass++)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imask())
			{
				model.Iref[iclass] *= (1. - DIRECT_MULTIDIM_ELEM(Imask(), n));
			}
		}

		// misuse the model.com_bodies and model.max_radius_mask_bodies vectors for centering the subtracted images
		// find center-of-mass for rotations around it
		int mydim = Imask().getDim();
		Matrix1D<RFLOAT> com(mydim);
		Imask().centerOfMass(com);
		// find maximum radius of mask around it's COM
		int max_d2 = 0.;
		FOR_ALL_ELEMENTS_IN_ARRAY3D(Imask())
		{
			if (A3D_ELEM(Imask(), k, i, j) > 0.05)
			{
				int d2 = (k - ZZ(com)) * (k - ZZ(com)) + (i - YY(com)) * (i - YY(com)) + (j - XX(com)) * (j - XX(com));
				max_d2 = XMIPP_MAX(max_d2, d2);
			}
		}
		model.max_radius_mask_bodies.push_back(CEIL(sqrt((RFLOAT)max_d2)));

		com.selfROUND();
		model.com_bodies.push_back(com);
	}

	// Now set up the Projectors inside the model
	model.setFourierTransformMaps(false); // false means ignore tau2_class

}

void FlexAnalyser::setup3DModels()
{
	for (int ibody = 0; ibody < model.nr_bodies; ibody++)
	{
		// Premultiply the map with the mask (otherwise need to do this again for every particle
		model.Iref[ibody] *= model.masks_bodies[ibody];
		// Place each body with its center-of-mass in the center of the box, as that's where the rotations are around
		selfTranslate(model.Iref[ibody], -model.com_bodies[ibody], DONT_WRAP);
		// And do the same for the masks
		selfTranslate(model.masks_bodies[ibody], -model.com_bodies[ibody], DONT_WRAP);
}

}

void FlexAnalyser::loopThroughParticles(int rank, int size)
{
	long int total_nr_particles = data.numberOfOriginalParticles();
	// Allow parallelisation
	long int my_first_ori_particle = 0, my_last_ori_particle = total_nr_particles-1;
	if (size > 1)
		divide_equally(total_nr_particles, size, rank, my_first_ori_particle, my_last_ori_particle);

	DFo.clear();
	DFo.setIsList(false);

	long int todo_particles = my_last_ori_particle-my_first_ori_particle+1;
	long int update_interval = todo_particles / 60;

	if (verb > 0)
		init_progress_bar(todo_particles);
	long int imgno = 0;
	for (long int ori_particle = my_first_ori_particle; ori_particle <= my_last_ori_particle; ori_particle++)
	{
		if (do_subtract)
		{
			subtractOneParticle(ori_particle, imgno, rank, size);
		}
		else if (do_3dmodels)
		{
			make3DModelOneParticle(ori_particle, imgno, rank, size);
		}
		else if (do_PCA_orient)
		{
			REPORT_ERROR("do_PCA_orient not implemented yet....");
		}

        if (imgno%update_interval==0 && verb > 0)
        	progress_bar(imgno);
        imgno++;
    }
	if (verb > 0)
		progress_bar(todo_particles);

	FileName fn_star, fn_addon = (do_subtract) ? "subtracted" : "3dmodels";
	if (size > 1)
		fn_star.compose(fn_out + "_rank", rank + 1, "_" + fn_addon + ".star");
	else
		fn_star = fn_out + "_" + fn_addon + ".star";
	DFo.write(fn_star);

}

void FlexAnalyser::subtractOneParticle(long int ori_particle, long int imgno, int rank, int size)
{
	// don't allow multiple particles per ori_particle!!!!
	if (data.ori_particles[ori_particle].particles_id.size() > 1)
		REPORT_ERROR("BUG: no movie particles allowed here...");

	long int part_id = data.ori_particles[ori_particle].particles_id[0];

	Image<RFLOAT> img;
	FileName fn_img;
	data.MDimg.getValue(EMDL_IMAGE_NAME, fn_img, part_id);
	img.read(fn_img);
	img().setXmippOrigin();

	// Get the consensus class, orientational parameters and norm (if present)
	Matrix1D<RFLOAT> my_old_offset(3), my_residual_offset(3);
	Matrix2D<RFLOAT> Aori;
	RFLOAT rot, tilt, psi, xoff, yoff, zoff, norm, scale;
	int myclass=0;
	if (data.MDimg.containsLabel(EMDL_PARTICLE_CLASS))
		data.MDimg.getValue(EMDL_PARTICLE_CLASS, myclass, part_id);
	data.MDimg.getValue(EMDL_ORIENT_ROT, rot, part_id);
	data.MDimg.getValue(EMDL_ORIENT_TILT, tilt, part_id);
	data.MDimg.getValue(EMDL_ORIENT_PSI, psi, part_id);
	data.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, XX(my_old_offset), part_id);
	data.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, YY(my_old_offset), part_id);
	if (model.data_dim == 3)
		data.MDimg.getValue(EMDL_ORIENT_ORIGIN_Z, ZZ(my_old_offset), part_id);
	if (!data.MDimg.getValue(EMDL_IMAGE_NORM_CORRECTION, norm, part_id))
		norm = 1.;


	// 17May2017: Shift image to the projected COM for this body!
	// Aori is the original transformation matrix of the consensus refinement
	Matrix1D<RFLOAT> my_projected_com(3), my_refined_ibody_offset(3);
	Euler_angles2matrix(rot, tilt, psi, Aori, false);
	my_projected_com = Aori * model.com_bodies[subtract_body];

	// Subtract the projected COM offset, to position this body in the center
	my_old_offset -= my_projected_com;
	my_residual_offset = my_old_offset;
	if (model.nr_bodies > 1)
	{
		// Also get refined offset for this body
		data.MDbodies[subtract_body].getValue(EMDL_ORIENT_ORIGIN_X, XX(my_refined_ibody_offset), part_id);
		data.MDbodies[subtract_body].getValue(EMDL_ORIENT_ORIGIN_Y, YY(my_refined_ibody_offset), part_id);
		if (model.data_dim == 3)
			data.MDbodies[subtract_body].getValue(EMDL_ORIENT_ORIGIN_Z, ZZ(my_refined_ibody_offset), part_id);
	}
	else
		// Without bodies, there is no refined offset, so set to zero
		my_refined_ibody_offset.initZeros();

	// Apply the norm_correction term
	if (do_norm)
		img() *= model.avg_norm_correction / norm;

	// Apply the old_offset (rounded to avoid interpolation errors)
	my_old_offset.selfROUND();
	selfTranslate(img(), my_old_offset, DONT_WRAP);
	// keep track of the differences between the rounded and the original offsets
	my_residual_offset -= my_old_offset;

	// Get the FourierTransform of the particle
	MultidimArray<Complex> Faux, Fimg;
	MultidimArray<RFLOAT> Fctf;
	FourierTransformer transformer;
	CenterFFT(img(), true);
	transformer.FourierTransform(img(), Faux);
	windowFourierTransform(Faux, Fimg, model.current_size);
	Fctf.resize(Fimg);

	if (do_ctf)
	{
		if (model.data_dim == 3)
		{
			Image<RFLOAT> Ictf;
			FileName fn_ctf;
			data.MDimg.getValue(EMDL_CTF_IMAGE, fn_ctf, part_id);
			Ictf.read(fn_ctf);
			Ictf().setXmippOrigin();
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf)
			{
				// Use negative kp,ip and jp indices, because the origin in the ctf_img lies half a pixel to the right of the actual center....
				DIRECT_A3D_ELEM(Fctf, k, i, j) = A3D_ELEM(Ictf(), -kp, -ip, -jp);
			}
		}
		else
		{
			CTF ctf;
			ctf.read(data.MDimg, data.MDimg, part_id);
			ctf.getFftwImage(Fctf, model.ori_size, model.ori_size, model.pixel_size,
					ctf_phase_flipped, false, intact_ctf_first_peak, true);
		}
	}
	else
		Fctf.initConstant(1.);

	MultidimArray<Complex> Fsubtract;
	Fsubtract.initZeros(Fimg);

	// For multi-body refinement
	Matrix2D<RFLOAT> Aresi_subtract;
	if (model.nr_bodies > 1)
	{
		for (int obody = 0; obody < model.nr_bodies; obody++)
		{
			// Unlike getFourierTransformsAndCtfs, no check for ibody==obody: also subtract rest of subtract_body!

			Matrix1D<RFLOAT> body_offset(3);
			RFLOAT body_rot, body_tilt, body_psi;
			data.MDbodies[obody].getValue(EMDL_ORIENT_ROT, body_rot, part_id);
			data.MDbodies[obody].getValue(EMDL_ORIENT_TILT, body_tilt, part_id);
			data.MDbodies[obody].getValue(EMDL_ORIENT_PSI,  body_psi, part_id);
			data.MDbodies[obody].getValue(EMDL_ORIENT_ORIGIN_X, XX(body_offset), part_id);
			data.MDbodies[obody].getValue(EMDL_ORIENT_ORIGIN_Y, YY(body_offset), part_id);
			if (model.data_dim == 3)
				data.MDbodies[obody].getValue(EMDL_ORIENT_ORIGIN_Z, ZZ(body_offset), part_id);

			Matrix2D<RFLOAT> Aresi,  Abody;
			// Aresi is the residual orientation for this obody
			Euler_angles2matrix(body_rot, body_tilt, body_psi, Aresi, false);
			if (obody == subtract_body)
				Aresi_subtract = Aresi;
			// The real orientation to be applied is the obody transformation applied and the original one
			Abody = Aori * (model.orient_bodies[obody]).transpose() * A_rot90 * Aresi * model.orient_bodies[obody];

			// Get the FT of the projection in the right direction
			MultidimArray<Complex> FTo;
			FTo.initZeros(Fimg);
			// The following line gets the correct pointer to account for overlap in the bodies
			int oobody = DIRECT_A2D_ELEM(model.pointer_body_overlap, subtract_body, obody);
			model.PPref[oobody].get2DFourierTransform(FTo, Abody, IS_NOT_INV);

			// Body is centered at its own COM: move it back to its place in the original particle image

			// Projected COM for this body (using Aori, just like above for ibody and my_projected_com!!!)
			Matrix1D<RFLOAT> other_projected_com(3);
			other_projected_com = Aori * (model.com_bodies[obody]);

			// Subtract refined obody-displacement
			other_projected_com -= body_offset;

			// Subtract the projected COM already applied to this image for ibody
			other_projected_com -= my_projected_com;

			shiftImageInFourierTransform(FTo, Faux, (RFLOAT)model.ori_size,
					XX(other_projected_com), YY(other_projected_com), ZZ(other_projected_com));

			// Sum the Fourier transforms of all the obodies
			Fsubtract += Faux;
		} // end for obody
	} // end if nr_bodies > 1
	else
	{
		Faux.initZeros(Fimg);
		// The original PPref is in the original position in the volume,
		// but Fimg was selfTranslated according to my_old_offset (which has placed it at the COM of the keep_inside mask)
		model.PPref[myclass].get2DFourierTransform(Faux, Aori, IS_NOT_INV);
		shiftImageInFourierTransform(Faux, Fsubtract, (RFLOAT)model.ori_size,
				XX(my_old_offset), YY(my_old_offset), ZZ(my_old_offset));

	}

	// Now that we have all the summed projections of the obodies, apply CTF, masks etc
	// Apply the CTF to this reference projection
	if (do_ctf)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fsubtract)
		{
			DIRECT_MULTIDIM_ELEM(Fsubtract, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
		}
	}

	// Do the actual subtraction
	Fimg -= Fsubtract;

	// And go finally back to real-space
	windowFourierTransform(Fimg, Faux, model.ori_size);
	transformer.inverseFourierTransform(Faux, img());
	CenterFFT(img(), false);

    // Set the entry in the output STAR file
	DFo.addObject();
    DFo.setObject(data.MDimg.getObject(part_id));

    // Set orientations back into the original RELION system of coordinates
	if (model.nr_bodies > 1)
	{
		Matrix2D<RFLOAT> Abody;
		// Write out the rot,tilt,psi as the combination of Aori and Aresi!! So get rid of the rotations around the tilt=90 axes,
		Abody = Aori * (model.orient_bodies[subtract_body]).transpose() * A_rot90 * Aresi_subtract * model.orient_bodies[subtract_body];
		Euler_matrix2angles(Abody, rot, tilt, psi);
		DFo.setValue(EMDL_ORIENT_ROT, rot, part_id);
		DFo.setValue(EMDL_ORIENT_TILT, tilt, part_id);
		DFo.setValue(EMDL_ORIENT_PSI, psi, part_id);
		my_residual_offset += my_refined_ibody_offset;
	}

	// Set the difference between the rounded my_old_offset and the actual offsets, plus (for multibody) my_refined_ibody_offset
	DFo.setValue(EMDL_ORIENT_ORIGIN_X, XX(my_residual_offset), part_id);
	DFo.setValue(EMDL_ORIENT_ORIGIN_Y, YY(my_residual_offset), part_id);
	if (model.data_dim == 3)
		DFo.setValue(EMDL_ORIENT_ORIGIN_Z, ZZ(my_residual_offset), part_id);


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


	// Now write it all out
	if (model.data_dim == 3)
    {
        fn_img.compose(fn_out, imgno+1,"mrc");
        img.write(fn_img);
    }
    else
    {
    	// Write this particle to the stack on disc
        // First particle: write stack in overwrite mode, from then on just append to it
        FileName fn_stack;
        if (size > 1)
        	fn_stack.compose(fn_out + "_rank", rank+1, "_subtracted.mrcs");
        else
        	fn_stack = fn_out + "_subtracted.mrcs";
    	fn_img.compose(imgno+1,fn_stack);
        if (imgno == 0)
            img.write(fn_img, -1, false, WRITE_OVERWRITE);
        else
            img.write(fn_img, -1, false, WRITE_APPEND);
    }
    DFo.setValue(EMDL_IMAGE_NAME,fn_img);

}

void FlexAnalyser::make3DModelOneParticle(long int ori_particle, long int imgno, int rank, int size)
{
	// don't allow multiple particles per ori_particle!!!!
	if (data.ori_particles[ori_particle].particles_id.size() > 1)
		REPORT_ERROR("BUG: no movie particles allowed here...");

	long int part_id = data.ori_particles[ori_particle].particles_id[0];

	// Get the consensus class, orientational parameters and norm (if present)
	Matrix1D<RFLOAT> my_old_offset(3);
	Matrix2D<RFLOAT> Aori;
	RFLOAT rot, tilt, psi, xoff, yoff, zoff;
	data.MDimg.getValue(EMDL_ORIENT_ROT, rot, part_id);
	data.MDimg.getValue(EMDL_ORIENT_TILT, tilt, part_id);
	data.MDimg.getValue(EMDL_ORIENT_PSI, psi, part_id);
	data.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, XX(my_old_offset), part_id);
	data.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, YY(my_old_offset), part_id);
	if (model.data_dim == 3)
		data.MDimg.getValue(EMDL_ORIENT_ORIGIN_Z, ZZ(my_old_offset), part_id);
	Euler_angles2matrix(rot, tilt, psi, Aori, false);

	Image<RFLOAT> img;
	MultidimArray<RFLOAT> sumw;
	img().initZeros(model.Iref[0]);
	sumw.initZeros(model.Iref[0]);

	for (int ibody = 0; ibody < model.nr_bodies; ibody++)
	{
	    MultidimArray<RFLOAT> Mbody, Mmask;
		Matrix1D<RFLOAT> body_offset(3), body_offset_3d(3);
		RFLOAT body_rot, body_tilt, body_psi;
		data.MDbodies[ibody].getValue(EMDL_ORIENT_ROT, body_rot, part_id);
		data.MDbodies[ibody].getValue(EMDL_ORIENT_TILT, body_tilt, part_id);
		data.MDbodies[ibody].getValue(EMDL_ORIENT_PSI,  body_psi, part_id);
		data.MDbodies[ibody].getValue(EMDL_ORIENT_ORIGIN_X, XX(body_offset), part_id);
		data.MDbodies[ibody].getValue(EMDL_ORIENT_ORIGIN_Y, YY(body_offset), part_id);
		if (model.data_dim == 3)
			data.MDbodies[ibody].getValue(EMDL_ORIENT_ORIGIN_Z, ZZ(body_offset), part_id);

		Matrix2D<RFLOAT> Aresi,  Abody;
		// Aresi is the residual orientation for this ibody
		Euler_angles2matrix(body_rot, body_tilt, body_psi, Aresi);
		// Only apply the residual orientation now!!!
		Abody = (model.orient_bodies[ibody]).transpose() * A_rot90 * Aresi * model.orient_bodies[ibody];

		// Now we have to get back from the 2D refined body_offset to some 3D translation of the body (with one direction non-defined)
		// We will need the original projection direction, Aori for that!!
		// Because one direction is ill-defined, this may not be such a good idea?
		// But anyway, this should bring it closer to truth than not doing anything at all...
		Aori = Aori * Abody;
		body_offset_3d = Aori.inv() * (-body_offset);

		// Also put back at the centre-of-mass of this body
		body_offset_3d += model.com_bodies[ibody];
		Abody.resize(4,4);

	    MAT_ELEM(Abody, 0, 3) = XX(body_offset_3d);
	    MAT_ELEM(Abody, 1, 3) = YY(body_offset_3d);
	    MAT_ELEM(Abody, 2, 3) = ZZ(body_offset_3d);
	    MAT_ELEM(Abody, 3, 3) = 1.;

		Mbody.resize(model.Iref[ibody]);
		Mmask.resize(model.masks_bodies[ibody]);
		applyGeometry(model.Iref[ibody], Mbody, Abody, IS_NOT_INV, DONT_WRAP);
	    applyGeometry(model.masks_bodies[ibody], Mmask, Abody, IS_NOT_INV, DONT_WRAP);

		img() += Mbody;
		sumw += Mmask;
	} // end for ibody

	// Divide the img by sumw to deal with overlapping bodies: just take average
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img())
	{
		if (DIRECT_MULTIDIM_ELEM(sumw, n) > 1.)
			DIRECT_MULTIDIM_ELEM(img(), n) /= DIRECT_MULTIDIM_ELEM(sumw, n);
	}
	// Write the image to disk
	FileName fn_img;
    fn_img.compose(fn_out, imgno+1,"mrc");
    img.setSamplingRateInHeader(model.pixel_size);
    img.write(fn_img);

    DFo.addObject();
    DFo.setValue(EMDL_MLMODEL_REF_IMAGE, fn_img);
    data.MDimg.getValue(EMDL_IMAGE_NAME, fn_img, part_id);
    DFo.setValue(EMDL_IMAGE_NAME, fn_img);
}
