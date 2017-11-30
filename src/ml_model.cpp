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

#include "src/ml_model.h"

#ifdef MDL_TIMING
    Timer mdl_timer;
	int TIMING_MDL_1 = 		proj_timer.setNew("MDL_1");
#define TIMING_TOC(id) mdl_timer.toc(id)
#else
#define TIMING_TIC(id)
#define TIMING_TOC(id)
#endif

void MlModel::initialise(bool _do_sgd)
{

	// Auxiliary vector with relevant size in Fourier space
	MultidimArray<RFLOAT > aux;
    aux.initZeros(ori_size / 2 + 1);

	// Now resize all relevant vectors
    Iref.resize(nr_classes * nr_bodies);
    masks_bodies.resize(nr_bodies);
    com_bodies.resize(nr_bodies);
    pdf_class.resize(nr_classes, 1./(RFLOAT)nr_classes);
    pdf_direction.resize(nr_classes);
    group_names.resize(nr_groups, "");
    sigma2_noise.resize(nr_groups, aux);
    nr_particles_group.resize(nr_groups);
    tau2_class.resize(nr_classes * nr_bodies, aux);
    fsc_halves_class.resize(aux);
    sigma2_class.resize(nr_classes * nr_bodies, aux);
    data_vs_prior_class.resize(nr_classes * nr_bodies, aux);
    fourier_coverage_class.resize(nr_classes * nr_bodies, aux);
    // TODO handle these two correctly.
    bfactor_correction.resize(nr_groups, 0.);
    scale_correction.resize(nr_groups, 1.);

	acc_rot.resize(nr_classes * nr_bodies, 0);
	acc_trans.resize(nr_classes * nr_bodies, 0);
	estimated_resolution.resize(nr_classes * nr_bodies, 0);
	total_fourier_coverage.resize(nr_classes * nr_bodies, 0);

	helical_twist.resize(nr_classes, 0);
	helical_rise.resize(nr_classes, 0);

	if (ref_dim==2)
	{
		Matrix1D<RFLOAT> empty(2);
		prior_offset_class.resize(nr_classes * nr_bodies, empty);
	}
	// These arrays will be resized when they are filled
	orientability_contrib.resize(nr_classes * nr_bodies);

	Projector ref(ori_size, interpolator, padding_factor, r_min_nn, data_dim);
    PPref.clear();
    PPrefRank.clear();
    // Now fill the entire vector with instances of "ref"
    PPref.resize(nr_classes * nr_bodies, ref);

    do_sgd = _do_sgd;
    if (do_sgd)
    	Igrad.resize(nr_classes);

}

// Reading from a file
void MlModel::read(FileName fn_in)
{

	// Clear current model
    clear();

    // Open input file
    std::ifstream in(fn_in.data(), std::ios_base::in);
    if (in.fail())
        REPORT_ERROR( (std::string) "MlModel::readStar: File " + fn_in + " cannot be read." );

    MetaDataTable MDclass, MDgroup, MDlog, MDsigma, MDbodies;

    // Read general stuff
    MDlog.readStar(in, "model_general");

	if (!MDlog.getValue(EMDL_MLMODEL_DIMENSIONALITY, ref_dim) ||
		!MDlog.getValue(EMDL_MLMODEL_ORIGINAL_SIZE, ori_size) ||
		!MDlog.getValue(EMDL_MLMODEL_CURRENT_RESOLUTION, current_resolution) ||
		!MDlog.getValue(EMDL_MLMODEL_CURRENT_SIZE, current_size) ||
		!MDlog.getValue(EMDL_MLMODEL_PADDING_FACTOR, padding_factor) ||
		!MDlog.getValue(EMDL_MLMODEL_INTERPOLATOR, interpolator) ||
		!MDlog.getValue(EMDL_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION, r_min_nn) ||
		!MDlog.getValue(EMDL_MLMODEL_PIXEL_SIZE, pixel_size) ||
		!MDlog.getValue(EMDL_MLMODEL_NR_CLASSES, nr_classes) ||
		!MDlog.getValue(EMDL_MLMODEL_NR_GROUPS, nr_groups) ||
		!MDlog.getValue(EMDL_MLMODEL_TAU2_FUDGE_FACTOR, tau2_fudge_factor) ||
		!MDlog.getValue(EMDL_MLMODEL_NORM_CORRECTION_AVG, avg_norm_correction) ||
		!MDlog.getValue(EMDL_MLMODEL_SIGMA_OFFSET, sigma2_offset) ||
		!MDlog.getValue(EMDL_MLMODEL_PRIOR_MODE, orientational_prior_mode) ||
		!MDlog.getValue(EMDL_MLMODEL_SIGMA_ROT, sigma2_rot) ||
		!MDlog.getValue(EMDL_MLMODEL_SIGMA_TILT, sigma2_tilt) ||
		!MDlog.getValue(EMDL_MLMODEL_SIGMA_PSI, sigma2_psi) ||
		!MDlog.getValue(EMDL_MLMODEL_LL, LL) ||
		!MDlog.getValue(EMDL_MLMODEL_AVE_PMAX, ave_Pmax) )
		REPORT_ERROR("MlModel::readStar: incorrect model_general table");

	// Retain compability with model files written by Relion prior to 1.4
	if (!MDlog.getValue(EMDL_MLMODEL_DIMENSIONALITY_DATA, data_dim))
		data_dim = 2;
	if (!MDlog.getValue(EMDL_MLMODEL_NR_BODIES, nr_bodies))
		nr_bodies = 1;
	if (!MDlog.getValue(EMDL_MLMODEL_IS_HELIX, is_helix))
		is_helix = false;
	if (is_helix)
	{
		if (nr_bodies != 1)
			REPORT_ERROR("MlModel::readStar: incorrect nr_bodies for helix");
		if (ref_dim == 2)
			REPORT_ERROR("MlModel::readStar: incorrect ref_dim for helix");
	}
	if (!MDlog.getValue(EMDL_MLMODEL_HELICAL_NR_ASU, helical_nr_asu))
		helical_nr_asu = 1;
	if (!MDlog.getValue(EMDL_MLMODEL_HELICAL_TWIST_MIN, helical_twist_min))
		helical_twist_min = 0.;
	if (!MDlog.getValue(EMDL_MLMODEL_HELICAL_TWIST_MAX, helical_twist_max))
		helical_twist_max = 0.;
	if (!MDlog.getValue(EMDL_MLMODEL_HELICAL_TWIST_INITIAL_STEP, helical_twist_inistep))
		helical_twist_inistep = 0.;
	if (!MDlog.getValue(EMDL_MLMODEL_HELICAL_RISE_MIN, helical_rise_min))
		helical_rise_min = 0.;
	if (!MDlog.getValue(EMDL_MLMODEL_HELICAL_RISE_MAX, helical_rise_max))
		helical_rise_max = 0.;
	if (!MDlog.getValue(EMDL_MLMODEL_HELICAL_RISE_INITIAL_STEP, helical_rise_inistep))
		helical_rise_inistep = 0.;

    // Treat classes or bodies (for multi-body refinement) in the same way...
    int nr_classes_bodies = (nr_bodies > 1) ? nr_bodies : nr_classes;

    if (nr_classes > 1 && nr_bodies > 1)
    	REPORT_ERROR("MlModel::readStar: nr_classes and nr_bodies cannot be both larger than one.");

    // Take inverse again of current resolution:
    current_resolution = 1. / current_resolution;

    sigma2_offset *= sigma2_offset;
	sigma2_rot *= sigma2_rot;
	sigma2_tilt *= sigma2_tilt;
	sigma2_psi *= sigma2_psi;

	// Resize vectors
	initialise();

	// Read classes
	FileName fn_tmp, fn_tmp2;
	Image<RFLOAT> img;
	if (nr_bodies > 1)
		MDclass.readStar(in, "model_bodies");
	else
		MDclass.readStar(in, "model_classes");
	int iclass = 0;
	do_sgd = false;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDclass)
	{
		if (!MDclass.getValue(EMDL_MLMODEL_REF_IMAGE, fn_tmp) ||
			!MDclass.getValue(EMDL_MLMODEL_ACCURACY_ROT, acc_rot[iclass]) ||
			!MDclass.getValue(EMDL_MLMODEL_ACCURACY_TRANS, acc_trans[iclass]) 	)
			REPORT_ERROR("MlModel::readStar: incorrect model_classes/bodies table");
		// backwards compatible
		if (!MDclass.getValue(EMDL_MLMODEL_ESTIM_RESOL_REF, estimated_resolution[iclass]))
			estimated_resolution[iclass] = 0.;
		if (!MDclass.getValue(EMDL_MLMODEL_FOURIER_COVERAGE_TOTAL_REF, total_fourier_coverage[iclass]))
			total_fourier_coverage[iclass] = 0.;
		if (ref_dim==2)
			if (!MDclass.getValue(EMDL_MLMODEL_PRIOR_OFFX_CLASS, XX(prior_offset_class[iclass])) ||
				!MDclass.getValue(EMDL_MLMODEL_PRIOR_OFFY_CLASS, YY(prior_offset_class[iclass])) )
				REPORT_ERROR("MlModel::readStar: incorrect model_classes/bodies table: no offset priors for 2D classes");
		if (nr_bodies == 1)
		{
			if (!MDclass.getValue(EMDL_MLMODEL_PDF_CLASS, pdf_class[iclass]) )
				REPORT_ERROR("MlModel::readStar: incorrect model_classes table: no pdf_class");
		}
		else
		{
			// Read in mask for this body
			if (!MDclass.getValue(EMDL_MASK_NAME, fn_tmp2) )
				REPORT_ERROR("MlModel::readStar: incorrect model_classes table: no body mask name");
			Image<RFLOAT> It;
			It.read(fn_tmp2);
			It().setXmippOrigin();
			masks_bodies[iclass] = It();
		}
		if (is_helix)
		{
			if (!MDclass.getValue(EMDL_MLMODEL_HELICAL_RISE, helical_rise[iclass]) ||
			    !MDclass.getValue(EMDL_MLMODEL_HELICAL_TWIST, helical_twist[iclass]) )
				REPORT_ERROR("MlModel::readStar: incorrect helical parameters");
		}

		// Read in actual reference image
		img.read(fn_tmp);
		Iref[iclass] = img();

		// Check to see whether there is a SGD-gradient entry as well
		if (MDclass.getValue(EMDL_MLMODEL_SGD_GRADIENT_IMAGE, fn_tmp))
		{
			do_sgd=true;
			if (iclass == 0)
				Igrad.resize(nr_classes);
			img.read(fn_tmp);
			Igrad[iclass] = img();
		}
		iclass++;
	}

	// Read group stuff
	MDgroup.readStar(in, "model_groups");
	long int igroup;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDgroup)
	{
        if (!MDgroup.getValue(EMDL_MLMODEL_GROUP_NO, igroup))
                REPORT_ERROR("MlModel::readStar: incorrect model_groups table");
        //Start counting of groups at 1, not at 0....
        if (!MDgroup.getValue(EMDL_MLMODEL_GROUP_SCALE_CORRECTION, scale_correction[igroup-1]) ||
                !MDgroup.getValue(EMDL_MLMODEL_GROUP_NR_PARTICLES, nr_particles_group[igroup-1]) ||
                !MDgroup.getValue(EMDL_MLMODEL_GROUP_NAME, group_names[igroup-1]))
                REPORT_ERROR("MlModel::readStar: incorrect model_groups table");
	}

	// Read SSNR, noise reduction, tau2_class spectra for each class
	for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
	{
		if (nr_bodies > 1)
			MDsigma.readStar(in, "model_body_" + integerToString(iclass + 1));
		else
			MDsigma.readStar(in, "model_class_" + integerToString(iclass + 1));
		int idx;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsigma)
		{
			if (!MDsigma.getValue(EMDL_SPECTRAL_IDX, idx))
				REPORT_ERROR("MlModel::readStar: incorrect table model_class/body_"+integerToString(iclass));
			if (!MDsigma.getValue(EMDL_MLMODEL_DATA_VS_PRIOR_REF, data_vs_prior_class[iclass](idx)) ||
			    !MDsigma.getValue(EMDL_MLMODEL_TAU2_REF, tau2_class[iclass](idx)) ||
			    !MDsigma.getValue(EMDL_MLMODEL_FSC_HALVES_REF, fsc_halves_class(idx)) ||
			    !MDsigma.getValue(EMDL_MLMODEL_SIGMA2_REF, sigma2_class[iclass](idx)))
				REPORT_ERROR("MlModel::readStar: incorrect table model_class/body_"+integerToString(iclass));
			// backwards compatible with STAR files without Fourier coverage
			if (!MDsigma.getValue(EMDL_MLMODEL_FOURIER_COVERAGE_REF, fourier_coverage_class[iclass](idx)))
				fourier_coverage_class[iclass](idx) = 0.;
		}
	}

	// Read sigma models for each group
	for (int igroup = 0; igroup < nr_groups; igroup++)
	{
		if (nr_particles_group[igroup] > 0)
		{
			MDsigma.readStar(in, "model_group_" + integerToString(igroup + 1));
			int idx;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsigma)
			{
				if (!MDsigma.getValue(EMDL_SPECTRAL_IDX, idx))
					REPORT_ERROR("MlModel::readStar: incorrect table model_group_"+integerToString(igroup));
				if (!MDsigma.getValue(EMDL_MLMODEL_SIGMA2_NOISE, sigma2_noise[igroup](idx)))
					REPORT_ERROR("MlModel::readStar: incorrect table model_group_"+integerToString(igroup));
			}
		}
		else
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
			{
				DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n) = 0.;
			}
		}
	}

	// Read pdf_direction models for each class
	if (ref_dim == 3)
	{
		for (int iclass = 0; iclass < nr_classes; iclass++)
		{
			MDclass.readStar(in, "model_pdf_orient_class_" + integerToString(iclass + 1));
			pdf_direction[iclass].clear();
			RFLOAT aux;
			std::vector<RFLOAT> vaux;
			vaux.clear();
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDclass)
			{
				if (!MDclass.getValue(EMDL_MLMODEL_PDF_ORIENT, aux))
					REPORT_ERROR("MlModel::readStar: incorrect table model_pdf_orient_class"+integerToString(iclass));
				vaux.push_back(aux);
			}
			pdf_direction[iclass].resize(vaux.size());
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(pdf_direction[iclass])
			{
				DIRECT_A1D_ELEM(pdf_direction[iclass], i) = vaux[i];
			}
			nr_directions = vaux.size();
		}
	}
	else
	{
		// For 2D case, just fill pdf_direction with ones.
		for (int iclass = 0; iclass < nr_classes; iclass++)
		{
			pdf_direction[iclass].clear();
			pdf_direction[iclass].resize(1);
			DIRECT_A1D_ELEM(pdf_direction[iclass], 0) = 1.;
		}
		nr_directions = 1;
	}

	// Close file handler
	in.close();

}

void MlModel::write(FileName fn_out, HealpixSampling &sampling, bool do_write_bild)
{

	MetaDataTable MDclass, MDgroup, MDlog, MDsigma, MDbodies;
    FileName fn_tmp, fn_tmp2;
    RFLOAT aux;
    std::ofstream  fh;

    // Treat classes or bodies (for multi-body refinement) in the same way...
    int nr_classes_bodies = (nr_bodies > 1) ? nr_bodies : nr_classes;
    // A. Write images
    if (ref_dim == 2)
    {
    	Image<RFLOAT> img(XSIZE(Iref[0]), YSIZE(Iref[0]), 1, nr_classes_bodies);
    	for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
    	{
    		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iref[iclass])
			{
    			DIRECT_NZYX_ELEM(img(), iclass, 0, i, j) = DIRECT_A2D_ELEM(Iref[iclass], i, j);
			}
    	}
    	if (nr_bodies > 1)
    		img.write(fn_out + "_bodies.mrcs");
    	else
    		img.write(fn_out + "_classes.mrcs");

    	if (do_sgd)
    	{
        	for (int iclass = 0; iclass < nr_classes; iclass++)
        	{
        		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Igrad[iclass])
    			{
        			DIRECT_NZYX_ELEM(img(), iclass, 0, i, j) = DIRECT_A2D_ELEM(Igrad[iclass], i, j);
    			}
        	}
       		img.write(fn_out + "_gradients.mrcs");
    	}
    }
    else
    {
    	Image<RFLOAT> img;
    	// Set correct voxel size in the header
		img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, pixel_size);
		img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, pixel_size);
		img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, pixel_size);
    	for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
    	{
    		if (nr_bodies > 1)
    			fn_tmp.compose(fn_out+"_body", iclass+1, "mrc", 3);
    		else
    			fn_tmp.compose(fn_out+"_class", iclass+1, "mrc", 3);

    		img() = Iref[iclass];
    		img.write(fn_tmp);
    	}
    	if (do_sgd)
    	{
			for (int iclass = 0; iclass < nr_classes; iclass++)
			{
				fn_tmp.compose(fn_out+"_grad", iclass+1, "mrc", 3);

				img() = Igrad[iclass];
				img.write(fn_tmp);
			}
    	}

    	if (do_write_bild)
    	{
			// Also write out bild files with the orientational distribution of each class
			// Also write out angular distributions
    		// Don't do this for bodies, only for classes!
			for (int iclass = 0; iclass < nr_classes; iclass++)
			{
				FileName fn_bild;
				fn_bild.compose(fn_out+"_class",iclass+1,"", 3);
				fn_bild += "_angdist.bild";
				RFLOAT offset = ori_size * pixel_size / 2.;
				sampling.writeBildFileOrientationalDistribution(pdf_direction[iclass], fn_bild, offset, offset);
			}
    	}
	}

    // B. Write STAR file with metadata
    fn_tmp = fn_out + "_model.star";
    fh.open((fn_tmp).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR( (std::string)"MlModel::write: Cannot write file: " + fn_tmp);

	// Write the output STAR file
	MDlog.setIsList(true);
	MDlog.addObject();
	MDlog.setName("model_general");
	MDlog.setValue(EMDL_MLMODEL_DIMENSIONALITY, ref_dim);
	MDlog.setValue(EMDL_MLMODEL_DIMENSIONALITY_DATA, data_dim);
	MDlog.setValue(EMDL_MLMODEL_ORIGINAL_SIZE, ori_size);
	MDlog.setValue(EMDL_MLMODEL_CURRENT_RESOLUTION, 1./current_resolution);
	MDlog.setValue(EMDL_MLMODEL_CURRENT_SIZE, current_size);
	MDlog.setValue(EMDL_MLMODEL_PADDING_FACTOR, padding_factor);
	MDlog.setValue(EMDL_MLMODEL_IS_HELIX, is_helix);
	if (is_helix)
	{
		MDlog.setValue(EMDL_MLMODEL_HELICAL_NR_ASU, helical_nr_asu);
		MDlog.setValue(EMDL_MLMODEL_HELICAL_TWIST_MIN, helical_twist_min);
		MDlog.setValue(EMDL_MLMODEL_HELICAL_TWIST_MAX, helical_twist_max);
		MDlog.setValue(EMDL_MLMODEL_HELICAL_TWIST_INITIAL_STEP, helical_twist_inistep);
		MDlog.setValue(EMDL_MLMODEL_HELICAL_RISE_MIN, helical_rise_min);
		MDlog.setValue(EMDL_MLMODEL_HELICAL_RISE_MAX, helical_rise_max);
		MDlog.setValue(EMDL_MLMODEL_HELICAL_RISE_INITIAL_STEP, helical_rise_inistep);
	}
	MDlog.setValue(EMDL_MLMODEL_INTERPOLATOR, interpolator);
	MDlog.setValue(EMDL_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION, r_min_nn);
	MDlog.setValue(EMDL_MLMODEL_PIXEL_SIZE, pixel_size);
	MDlog.setValue(EMDL_MLMODEL_NR_CLASSES, nr_classes);
	MDlog.setValue(EMDL_MLMODEL_NR_BODIES, nr_bodies);
	MDlog.setValue(EMDL_MLMODEL_NR_GROUPS, nr_groups);
	MDlog.setValue(EMDL_MLMODEL_TAU2_FUDGE_FACTOR, tau2_fudge_factor);
	MDlog.setValue(EMDL_MLMODEL_NORM_CORRECTION_AVG, avg_norm_correction);
	MDlog.setValue(EMDL_MLMODEL_SIGMA_OFFSET, sqrt(sigma2_offset));
	MDlog.setValue(EMDL_MLMODEL_PRIOR_MODE, orientational_prior_mode);
	MDlog.setValue(EMDL_MLMODEL_SIGMA_ROT, sqrt(sigma2_rot));
	MDlog.setValue(EMDL_MLMODEL_SIGMA_TILT, sqrt(sigma2_tilt));
	MDlog.setValue(EMDL_MLMODEL_SIGMA_PSI, sqrt(sigma2_psi));
	MDlog.setValue(EMDL_MLMODEL_LL, LL);
	MDlog.setValue(EMDL_MLMODEL_AVE_PMAX, ave_Pmax);
	MDlog.write(fh);

	// Calculate resolutions and total Fourier coverages for each class
	calculateTotalFourierCoverage();

    // Write metadata and images for all classes
	FileName fn_root;
	fn_root = fn_out.beforeFirstOf("_it");
	if (nr_bodies > 1)
		MDclass.setName("model_bodies");
	else
		MDclass.setName("model_classes");
	for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
	{
		MDclass.addObject();
		Image<RFLOAT> Itmp;
		if (ref_dim==2)
		{
			if (nr_bodies > 1)
			{
				fn_tmp = fn_out + "_bodies.mrcs";
				fn_tmp2.compose(fn_root+"_body",iclass+1,"", 3); // class number from 1 to K!
				fn_tmp2 += "_mask.mrc";
			}
			else
			{
				fn_tmp = fn_out + "_classes.mrcs";
			}
			fn_tmp.compose(iclass+1, fn_tmp); // fn_tmp = integerToString(iclass) + "@" + fn_tmp;
		}
		else
		{
			if (nr_bodies > 1)
			{
				fn_tmp.compose(fn_out+"_body",iclass+1,"mrc", 3); // class number from 1 to K!
				fn_tmp2.compose(fn_root+"_body",iclass+1,"", 3); // class number from 1 to K!
				fn_tmp2 += "_mask.mrc";
			}
			else
				fn_tmp.compose(fn_out+"_class",iclass+1,"mrc", 3); // class number from 1 to K!
		}
		MDclass.setValue(EMDL_MLMODEL_REF_IMAGE, fn_tmp);
		if (do_sgd)
		{
			if (ref_dim==2)
				fn_tmp.compose(iclass+1, fn_out + "_gradients.mrcs");
			else
				fn_tmp.compose(fn_out+"_grad",iclass+1,"mrc", 3);
			MDclass.setValue(EMDL_MLMODEL_SGD_GRADIENT_IMAGE, fn_tmp);
		}

		// Also set he maskname for multi-body refinement
		if (nr_bodies > 1)
			MDclass.setValue(EMDL_MASK_NAME, fn_tmp2);

		// For multiple bodies: only star PDF_CLASS in the first one!
		int myclass = (nr_bodies > 1) ? 0 : iclass; // for multi-body: just set iclass=0
		MDclass.setValue(EMDL_MLMODEL_PDF_CLASS, pdf_class[myclass]);
		MDclass.setValue(EMDL_MLMODEL_ACCURACY_ROT, acc_rot[iclass]);
		MDclass.setValue(EMDL_MLMODEL_ACCURACY_TRANS, acc_trans[iclass]);
		MDclass.setValue(EMDL_MLMODEL_ESTIM_RESOL_REF, estimated_resolution[iclass]);
		MDclass.setValue(EMDL_MLMODEL_FOURIER_COVERAGE_TOTAL_REF, total_fourier_coverage[iclass]);

		if (ref_dim==2)
		{
			MDclass.setValue(EMDL_MLMODEL_PRIOR_OFFX_CLASS, XX(prior_offset_class[iclass]));
			MDclass.setValue(EMDL_MLMODEL_PRIOR_OFFY_CLASS, YY(prior_offset_class[iclass]));
		}

		if (is_helix)
		{
			MDclass.setValue(EMDL_MLMODEL_HELICAL_RISE, helical_rise[iclass]);
			MDclass.setValue(EMDL_MLMODEL_HELICAL_TWIST, helical_twist[iclass]);
		}
	}
	MDclass.write(fh);

	// Write radial_average of tau2_class and data_vs_prior_class for each reference
	for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
	{
		MDsigma.clear();
		if (nr_bodies > 1)
			MDsigma.setName("model_body_"+integerToString(iclass+1));
		else
			MDsigma.setName("model_class_"+integerToString(iclass+1));
		for (int ii = 0; ii < XSIZE(tau2_class[iclass]); ii++)
		{
			MDsigma.addObject();
			MDsigma.setValue(EMDL_SPECTRAL_IDX, ii);
			MDsigma.setValue(EMDL_RESOLUTION, getResolution(ii));
			MDsigma.setValue(EMDL_RESOLUTION_ANGSTROM, getResolutionAngstrom(ii));
			MDsigma.setValue(EMDL_MLMODEL_DATA_VS_PRIOR_REF, data_vs_prior_class[iclass](ii));
			MDsigma.setValue(EMDL_MLMODEL_FSC_HALVES_REF, fsc_halves_class(ii));
			MDsigma.setValue(EMDL_MLMODEL_FOURIER_COVERAGE_REF, fourier_coverage_class[iclass](ii));
			MDsigma.setValue(EMDL_MLMODEL_SIGMA2_REF, sigma2_class[iclass](ii));
			MDsigma.setValue(EMDL_MLMODEL_TAU2_REF, tau2_class[iclass](ii));
			// Only write orientabilities if they have been determined
			if (XSIZE(orientability_contrib[iclass]) == XSIZE(tau2_class[iclass]))
				MDsigma.setValue(EMDL_MLMODEL_ORIENTABILITY_CONTRIBUTION, orientability_contrib[iclass](ii));
		}
		MDsigma.write(fh);
	}

    // Write scale-correction for all groups
    MDgroup.setName("model_groups");
    for (long int igroup = 0; igroup < nr_groups; igroup++)
    {
		MDgroup.addObject();
		//Start counting of groups at 1, not at 0....
		MDgroup.setValue(EMDL_MLMODEL_GROUP_NO, igroup+1);
		MDgroup.setValue(EMDL_MLMODEL_GROUP_NAME, group_names[igroup]);
		MDgroup.setValue(EMDL_MLMODEL_GROUP_NR_PARTICLES, nr_particles_group[igroup]);
		MDgroup.setValue(EMDL_MLMODEL_GROUP_SCALE_CORRECTION, scale_correction[igroup]);
    }
    MDgroup.write(fh);

	// Write sigma models for each group
	for (int igroup = 0; igroup < nr_groups; igroup++)
	{
		MDsigma.clear();
		MDsigma.setName("model_group_"+integerToString(igroup+1));
		for (int ii = 0; ii < XSIZE(sigma2_noise[igroup]); ii++)
		{
			MDsigma.addObject();
			// Some points in sigma2_noise arrays are never used...
			aux = sigma2_noise[igroup](ii);
			if (aux > 0.)
			{
				MDsigma.setValue(EMDL_SPECTRAL_IDX, ii);
				MDsigma.setValue(EMDL_RESOLUTION, getResolution(ii));
				MDsigma.setValue(EMDL_MLMODEL_SIGMA2_NOISE, aux);
			}
		}
		MDsigma.write(fh);
	}

	// Write pdf_direction models for each class
	if (ref_dim == 3)
	{
		for (int iclass = 0; iclass < nr_classes; iclass++)
		{
			MDclass.clear();
			MDclass.setName("model_pdf_orient_class_"+integerToString(iclass+1));
			for (int ii=0; ii < XSIZE(pdf_direction[iclass]); ii++)
			{
				MDclass.addObject();
				MDclass.setValue(EMDL_MLMODEL_PDF_ORIENT, pdf_direction[iclass](ii));
			}
			MDclass.write(fh);
		}
	}

}




void  MlModel::readTauSpectrum(FileName fn_tau, int verb)
{
	MetaDataTable MDtau;
	RFLOAT val;
	int idx;
	MDtau.read(fn_tau);
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDtau)
	{
		MDtau.getValue(EMDL_SPECTRAL_IDX, idx);
		MDtau.getValue(EMDL_MLMODEL_TAU2_REF, val);
		if (idx < XSIZE(tau2_class[0]))
			tau2_class[0](idx) = tau2_fudge_factor * val;
	}
	if (idx < XSIZE(tau2_class[0]) - 1)
	{
		if (verb > 0) std::cerr<< " Warning: provided tau2-spectrum has fewer entries ("<<idx+1<<") than needed ("<<XSIZE(tau2_class[0])<<"). Set rest to zero..."<<std::endl;
	}
	// Use the same spectrum for all classes
	for (int iclass = 0; iclass < nr_classes; iclass++)
		tau2_class[iclass] =  tau2_class[0];

}

// Reading images from disc
void MlModel::readImages(FileName fn_ref, bool _is_3d_model, int _ori_size, Experiment &_mydata,
			bool &do_average_unaligned, bool &do_generate_seeds, bool &refs_are_ctf_corrected, bool _do_sgd)
{

	// Set some stuff
	nr_groups = _mydata.groups.size();
	ori_size = _ori_size;
	RFLOAT avg_norm_correction = 1.;

	// Data dimensionality
	_mydata.MDexp.getValue(EMDL_IMAGE_DIMENSIONALITY, data_dim);

	// Read references into memory
	Image<RFLOAT> img;
	FileName fn_tmp;
	if (fn_ref != "None")
	{
		// Read the references into memory
		do_average_unaligned = false;
		// If this is a STAR file, ignore nr_classes and read all references from this file
		if (fn_ref.isStarFile())
		{
			MetaDataTable MDref;
			MDref.read(fn_ref,"model_classes");
			if(!MDref.getValue(EMDL_MLMODEL_REF_IMAGE, fn_tmp)) // if we did not find the meta-data label _rlnReferenceImage in a directed search, try more generally
				MDref.read(fn_ref);
			if(!MDref.getValue(EMDL_MLMODEL_REF_IMAGE, fn_tmp)) // if we still did not find the meta-data label _rlnReferenceImage, report an error
				REPORT_ERROR("When specifying a .star-file as --ref input, you need to have the _rlnReferenceImage field");

			do_generate_seeds = false;
			// ignore nr_classes from the command line, use number of entries in STAR file
			nr_classes = 0;
			Iref.clear();
			Igrad.clear();
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDref)
			{
				MDref.getValue(EMDL_MLMODEL_REF_IMAGE, fn_tmp);
				img.read(fn_tmp);
				ref_dim = img().getDim();
				if (ori_size != XSIZE(img()) || ori_size != YSIZE(img()))
				{
					std::cerr << " ori_size= " << ori_size << " XSIZE(img())= " << XSIZE(img()) << std::endl;
					REPORT_ERROR("MlOptimiser::read: size of reference images is not the same as the experimental images!");
				}
				Iref.push_back(img());
				if (_do_sgd)
				{
					img() *= 0.;
					Igrad.push_back(img());
				}
				nr_classes++;
			}
		}
		// For a single image, read this image as reference and set it in all nr_classes Irefs
		else
		{
			img.read(fn_ref);
			img().setXmippOrigin();
			ref_dim = img().getDim();
			if (ori_size != XSIZE(img()) || ori_size != YSIZE(img()))
			{
				std::cerr << " ori_size= " << ori_size << " XSIZE(img())= " << XSIZE(img()) << std::endl;
				REPORT_ERROR("MlOptimiser::read: size of reference image is not the same as the experimental images!");
			}
			Iref.clear();
			Igrad.clear();
			if (nr_bodies > 1)
			{
				for (int ibody = 0; ibody < nr_bodies; ibody++)
				{
					Iref.push_back(img());
					if (masks_bodies.size() <= ibody)
						REPORT_ERROR("BUG: masks_bodies.size() < ibody. Did you initialise the body masks before reading the references?");
					Iref[ibody] *= masks_bodies[ibody];
				}
			}
			else
			{
				for (int iclass = 0; iclass < nr_classes; iclass++)
				{
					Iref.push_back(img());
					if (_do_sgd)
					{
						img() *= 0.;
						Igrad.push_back(img());
					}
				}
			}
			if (nr_classes > 1)
				do_generate_seeds = true;
			else
				do_generate_seeds = false;
		}
	}
	else
	{
		// If no -ref is given, calculate average of all unaligned images later on.
		do_average_unaligned = true;
		do_generate_seeds = false; // after SGD introduction, this is now done in the estimation of initial sigma2 step!
		refs_are_ctf_corrected = true;
		if (_is_3d_model || data_dim == 3)
		{
			ref_dim = 3;
			img().initZeros(ori_size, ori_size, ori_size);
		}
		else
		{
			ref_dim = 2;
			img().initZeros(ori_size, ori_size);
		}
		Iref.clear();
		Igrad.clear();
		for (int iclass = 0; iclass < nr_classes; iclass++)
		{
			Iref.push_back(img());
			if (_do_sgd)
				Igrad.push_back(img());
		}
	}

	initialise(_do_sgd);

	// Now set the group names from the Experiment groups list
	for (int i=0; i< nr_groups; i++)
		group_names[i] = _mydata.groups[i].name;

}

void MlModel::reassignGroupsForMovies(Experiment &mydata, std::string &movie_name)
{

	std::vector<long int> rename_ids(mydata.groups.size());
	for (long int igr = 0; igr < mydata.groups.size(); igr++)
	{
		FileName data_name = (mydata.groups[igr].name);
		data_name = data_name.beforeLastOf("_"+movie_name);
		long int rename_id = -1;
		for (long int id = 0; id < group_names.size(); id++)
		{
			if (data_name == group_names[id].withoutExtension())
			{
				rename_id = id;
				break;
			}
		}

		if (rename_id < 0)
			REPORT_ERROR("MlModel::adjustGroupsForMovies ERROR: cannot find " + data_name + " among the groups of the model!");

		rename_ids[igr] = rename_id;
	}

	// Now change the group_ids of all particles!
	for (long int ipart = 0; ipart < mydata.particles.size(); ipart++)
	{
		long int old_id = mydata.particles[ipart].group_id;
		mydata.particles[ipart].group_id = rename_ids[old_id];
	}

}

void MlModel::initialisePdfDirection(int newsize)
{

	// If the pdf_direction were already filled (size!=0), and newsize=oldsize then leave them as they were
	// If they were still empty, or if the size changes, then initialise them with an even distribution
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		int oldsize = MULTIDIM_SIZE(pdf_direction[iclass]);
		if (oldsize == 0 || oldsize != newsize)
		{
			pdf_direction[iclass].resize(newsize);
			pdf_direction[iclass].initConstant(1./((RFLOAT) nr_classes * newsize));
		}
	}
	nr_directions = newsize;

}

void MlModel::initialiseBodyMasks(FileName fn_masks, FileName fn_root_out)
{
	MetaDataTable MD;
	MD.read(fn_masks);
	if (!MD.containsLabel(EMDL_MASK_NAME))
		REPORT_ERROR("ERROR MlModel::initialiseBodyMasks: body-mask STAR file does not contain rlnBodyMaskName label.");

	nr_bodies = 0;
	masks_bodies.resize(MD.numberOfObjects());
	com_bodies.resize(MD.numberOfObjects());
	FileName fn_mask;
	Image<RFLOAT> Imask;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		MD.getValue(EMDL_MASK_NAME, fn_mask);
		Imask.read(fn_mask);
		Imask().setXmippOrigin();
		masks_bodies[nr_bodies] = Imask();
		// find center-of-mass for rotations around it
		int mydim = Imask().getDim();
		Matrix1D<RFLOAT> com(mydim);
		Imask().centerOfMass(com);
		com_bodies[nr_bodies].resize(3);
		XX(com_bodies[nr_bodies]) = ROUND(XX(com)); // ROUND so no interpolation artifacts in selfTranslate(Iref)
		YY(com_bodies[nr_bodies]) = ROUND(YY(com));
		if (mydim == 3)
			ZZ(com_bodies[nr_bodies]) = ROUND(ZZ(com));
		else
			ZZ(com_bodies[nr_bodies]) = 0.;
		// Also write the mask with the standard name to disk

		fn_mask.compose(fn_root_out + "_body", nr_bodies + 1, "", 3); // body number from 1 to K!
		fn_mask += "_mask.mrc";

		Imask.write(fn_mask);
		// update counter at the end!
		nr_bodies++;
	}

}



void MlModel::setFourierTransformMaps(bool update_tau2_spectra, int nr_threads, bool do_gpu)
{
	bool do_heavy(true);
	int nr_classes_bodies = nr_classes * nr_bodies; // also set multiple bodies!
	for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
    {

		MultidimArray<RFLOAT> Irefp;
		//19may2015: if multi-body refinement: place each body with its center-of-mass in the center
		if (nr_bodies > 1)
		{
			translate(Iref[iclass], Irefp, -com_bodies[iclass], DONT_WRAP);
		}
		else
		{
			Irefp = Iref[iclass];
		}

		if(PPrefRank.size() > 1)
			do_heavy = PPrefRank[iclass];

        if (update_tau2_spectra)
        {
        	PPref[iclass].computeFourierTransformMap(Irefp, tau2_class[iclass], current_size, nr_threads, true, do_heavy);
        }
        else
        {
        	MultidimArray<RFLOAT> dummy;
        	PPref[iclass].computeFourierTransformMap(Irefp, dummy, current_size, nr_threads, true, do_heavy);
        }
    }

}

void MlModel::initialiseDataVersusPrior(bool fix_tau)
{

    // Get total number of particles
	RFLOAT nr_particles = 0.;
	for (int igroup = 0; igroup < nr_particles_group.size(); igroup++)
		nr_particles += (RFLOAT)nr_particles_group[igroup];

	// Calculate average sigma2_noise over all image groups
	MultidimArray<RFLOAT> avg_sigma2_noise;
	avg_sigma2_noise.initZeros(sigma2_noise[0]);
	for (int igroup = 0; igroup < nr_particles_group.size(); igroup++)
	{
		avg_sigma2_noise += (RFLOAT)(nr_particles_group[igroup]) * sigma2_noise[igroup];
	}
	avg_sigma2_noise /= nr_particles;

	// Get the FT of all reference structures
    // The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
    // And spectrum is squared, so ori_size*ori_size in the 3D case!
	RFLOAT normfft = (ref_dim == 3 && data_dim == 2) ? (RFLOAT)(ori_size * ori_size) : 1.;

    for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		// Initialise output arrays to correct size
		tau2_class[iclass].resize(sigma2_noise[0]);

		// Get the power spectrum of the reference
		MultidimArray<RFLOAT> spectrum(sigma2_noise[0]);
		getSpectrum(Iref[iclass], spectrum, POWER_SPECTRUM);

		// Factor two because of two-dimensionality of the complex plane
		// (just like sigma2_noise estimates, the power spectra should be divided by 2)
		spectrum *= normfft / 2.;

		// Update the tau2_class spectrum for this reference
		// This is only for writing out in the it000000_model.star file
		if (!fix_tau)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(tau2_class[iclass])
			{
				DIRECT_A1D_ELEM(tau2_class[iclass], i) = tau2_fudge_factor * DIRECT_A1D_ELEM(spectrum, i);
			}
		}

		// Calculate data_vs_prior_class as spectral_nr_observations_per_class/sigma2_noise vs 1/tau2_class
		data_vs_prior_class[iclass].resize(sigma2_noise[0]);
		fsc_halves_class.initZeros(sigma2_noise[0]);
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(tau2_class[iclass])
		{
			RFLOAT evidence = nr_particles * pdf_class[iclass] / DIRECT_A1D_ELEM(avg_sigma2_noise, i);
			// empirical accounting for ratio of pixels in 3D shells compared to 2D shells
			if (ref_dim == 3 && i > 0)
				evidence /= (2. * (RFLOAT)i);
			RFLOAT prior = 1. /  DIRECT_A1D_ELEM(tau2_class[iclass], i);
			RFLOAT myssnr = evidence / prior;
			DIRECT_A1D_ELEM(data_vs_prior_class[iclass], i ) = myssnr;
			// Also initialise FSC-halves here (...)
			//DIRECT_A1D_ELEM(fsc_halves_class[iclass], i ) = myssnr / (myssnr + 1);
		}
	} // end loop iclass

}

void MlModel::initialiseHelicalParametersLists(RFLOAT _helical_twist, RFLOAT _helical_rise)
{
    if (nr_classes < 1)
    	REPORT_ERROR("MlModel.cpp::initialiseHelicalParametersLists  nr_classes is smaller than 1");
    helical_twist.resize(nr_classes);
    helical_rise.resize(nr_classes);
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
    	helical_twist[iclass] = _helical_twist;
    	helical_rise[iclass] = _helical_rise;
    }
}

void MlModel::calculateTotalFourierCoverage()
{
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		int maxres = 0;
		for (int ires = 0; ires < XSIZE(data_vs_prior_class[iclass]); ires++)
		{
			if (DIRECT_A1D_ELEM(data_vs_prior_class[iclass], ires) < 1.)
				break;
			maxres = ires;
		}

		estimated_resolution[iclass] = 1./getResolution(maxres);
		total_fourier_coverage[iclass] = 0.;
		RFLOAT count = 0;
		for (long int k=FIRST_XMIPP_INDEX(maxres+2); k<=FIRST_XMIPP_INDEX(maxres+2) + maxres+1; k++) \
		   for (long int i=FIRST_XMIPP_INDEX(maxres+2); i<=FIRST_XMIPP_INDEX(maxres+2) + maxres+1; i++) \
			   for (long int j=FIRST_XMIPP_INDEX(maxres+2); j<=FIRST_XMIPP_INDEX(maxres+2) + maxres+1; j++) \
			   {
				   int r = sqrt(RFLOAT(k*k+i*i+j*j));
				   if (r <= maxres)
				   {
					   total_fourier_coverage[iclass] += DIRECT_A1D_ELEM(fourier_coverage_class[iclass], r);
					   count += 1.;
				   }
			   }
		total_fourier_coverage[iclass] /= count;
	}

}


/////////// MlWsumModel
void MlWsumModel::initialise(MlModel &_model, FileName fn_sym, bool asymmetric_padding, bool _skip_gridding)
{
	nr_classes = _model.nr_classes;
	nr_bodies = _model.nr_bodies;
    nr_groups = _model.nr_groups;
    nr_directions = _model.nr_directions;
    ref_dim = _model.ref_dim;
    data_dim = _model.data_dim;
    ori_size = _model.ori_size;
    pdf_class = _model.pdf_class;
    if (ref_dim == 2)
    	prior_offset_class = _model.prior_offset_class;
    pdf_direction = _model.pdf_direction;
    sigma2_offset = _model.sigma2_offset;
    sigma2_noise = _model.sigma2_noise;
    sigma2_rot = _model.sigma2_rot;
    sigma2_tilt = _model.sigma2_tilt;
    sigma2_psi = _model.sigma2_psi;
    interpolator = _model.interpolator;
    r_min_nn = _model.r_min_nn;
	is_helix = _model.is_helix;
	helical_nr_asu = _model.helical_nr_asu;
	helical_twist_min = _model.helical_twist_min;
	helical_twist_max = _model.helical_twist_max;
	helical_twist_inistep = _model.helical_twist_inistep;
	helical_rise_min = _model.helical_rise_min;
	helical_rise_max = _model.helical_rise_max;
	helical_rise_inistep = _model.helical_rise_inistep;

    padding_factor = _model.padding_factor;
    if (asymmetric_padding)
    	padding_factor ++;

    // Don't need forward projectors in MlWsumModel!
    PPref.clear();
    // Don't need scale_correction and bfactor_correction, keep wsum_signal_product_spectra and wsum_reference_power_spectra instead
    scale_correction.clear();
    bfactor_correction.clear();
    tau2_class.clear();
    data_vs_prior_class.clear();
	acc_rot.clear();
	acc_trans.clear();
	estimated_resolution.clear();
	total_fourier_coverage.clear();
    orientability_contrib.clear();

	helical_twist.resize(nr_classes);
	helical_rise.resize(nr_classes);
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		helical_twist[iclass] = _model.helical_twist[iclass];
		helical_rise[iclass] = _model.helical_rise[iclass];
	}

    MultidimArray<RFLOAT> aux(ori_size / 2 + 1);
    wsum_signal_product_spectra.resize(nr_groups, aux);
    wsum_reference_power_spectra.resize(nr_groups, aux);

    // Resize MlWsumModel-specific vectors
    BackProjector BP(ori_size, ref_dim, fn_sym, interpolator, padding_factor, r_min_nn,
    		         ML_BLOB_ORDER, ML_BLOB_RADIUS, ML_BLOB_ALPHA, data_dim, _skip_gridding);
    BPref.clear();
    BPref.resize(nr_classes * nr_bodies, BP); // also set multiple bodies
    sumw_group.resize(nr_groups);

}

void MlWsumModel::initZeros()
{

    LL = 0.;
    ave_Pmax = 0.;
    sigma2_offset = 0.;
    avg_norm_correction = 0.;
    sigma2_rot = 0.;
    sigma2_tilt = 0.;
    sigma2_psi = 0.;

    // Set all weighted sums to zero

    for (int iclass = 0; iclass < nr_classes * nr_bodies; iclass++)
    {
    	BPref[iclass].initZeros(current_size);
    }

    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        pdf_class[iclass] = 0.;
        // Assume pdf_direction is already of the right size...
        pdf_direction[iclass].initZeros();
        if (ref_dim == 2)
        	prior_offset_class[iclass].initZeros();
    }

    // Initialise sigma2_noise spectra and sumw_group
    for (int igroup = 0; igroup < nr_groups; igroup++)
    {
        sumw_group[igroup] = 0.;
        sigma2_noise[igroup].initZeros();
        wsum_signal_product_spectra[igroup].initZeros();
        wsum_reference_power_spectra[igroup].initZeros();
    }
}

//#define DEBUG_PACK
#ifdef DEBUG_PACK
#define MAX_PACK_SIZE     100000
#else
// Approximately 1024*1024*1024/8/2 ~ 0.5 Gb
#define MAX_PACK_SIZE 671010000
#endif

void MlWsumModel::pack(MultidimArray<RFLOAT> &packed)
{
    // for LL & avePmax & sigma2_offset & avg_norm_correction & sigma2_rot & sigma2_tilt & sigma2_psi
    unsigned long long packed_size = 0;
    int spectral_size = (ori_size / 2) + 1;

    packed_size += 7 ;
    // for all group-related stuff
    packed_size += nr_groups * spectral_size;
    packed_size += nr_groups * spectral_size;
    packed_size += nr_groups * spectral_size;
    // for sumw_group
    packed_size += nr_groups;
    // for all class-related stuff
    // data is complex: multiply by two!
    packed_size += nr_classes * nr_bodies * 2 * BPref[0].getSize();
    packed_size += nr_classes * nr_bodies * BPref[0].getSize();
    packed_size += nr_classes * nr_directions;
    // for pdf_class
    packed_size += nr_classes;
    // for priors for each class
    if (ref_dim==2)
    	packed_size += nr_classes*2;

    // Get memory for the packed array
    packed.clear();
    packed.resize(packed_size);

    // Start packing
    unsigned long long idx = 0;

    DIRECT_MULTIDIM_ELEM(packed, idx++) = LL;
    DIRECT_MULTIDIM_ELEM(packed, idx++) = ave_Pmax;
    DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_offset;
    DIRECT_MULTIDIM_ELEM(packed, idx++) = avg_norm_correction;
    DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_rot;
    DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_tilt;
    DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_psi;

    for (int igroup = 0; igroup < nr_groups; igroup++)
    {
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
        {
            DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n);
        }
    	sigma2_noise[igroup].clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_signal_product_spectra[igroup])
        {
        	DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(wsum_signal_product_spectra[igroup], n);
        }
        wsum_signal_product_spectra[igroup].clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_reference_power_spectra[igroup])
        {
            DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(wsum_reference_power_spectra[igroup], n);
        }
        wsum_reference_power_spectra[igroup].clear();

        DIRECT_MULTIDIM_ELEM(packed, idx++) = sumw_group[igroup];

    }
    for (int iclass = 0; iclass < nr_classes * nr_bodies; iclass++)
    {

    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].data)
        {
            DIRECT_MULTIDIM_ELEM(packed, idx++) = (DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).real;
            DIRECT_MULTIDIM_ELEM(packed, idx++) = (DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).imag;
        }
    	BPref[iclass].data.clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].weight)
        {
            DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(BPref[iclass].weight, n);
        }
        BPref[iclass].weight.clear();
    }
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdf_direction[iclass])
        {
            DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], n);
        }
        pdf_direction[iclass].clear();

        DIRECT_MULTIDIM_ELEM(packed, idx++) = pdf_class[iclass];

        if (ref_dim==2)
        {
        	DIRECT_MULTIDIM_ELEM(packed, idx++) = XX(prior_offset_class[iclass]);
        	DIRECT_MULTIDIM_ELEM(packed, idx++) = YY(prior_offset_class[iclass]);
        }
    }
#ifdef DEBUG_PACK
    std::cerr << " idx= " << idx << " packed_size= " << packed_size << std::endl;
#endif

    // Just to check whether we went outside our memory...
    if (idx != packed_size)
    {
       	std::cerr << "idx= " << idx << "packed_size= " << packed_size << std::endl;
        REPORT_ERROR("MlWsumModel::pack: idx != packed_size");
    }

}
void MlWsumModel::unpack(MultidimArray<RFLOAT> &packed)
{
    int spectral_size = (ori_size / 2) + 1;

    unsigned long long idx = 0;

    LL = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ave_Pmax = DIRECT_MULTIDIM_ELEM(packed, idx++);
    sigma2_offset = DIRECT_MULTIDIM_ELEM(packed, idx++);
    avg_norm_correction = DIRECT_MULTIDIM_ELEM(packed, idx++);
    sigma2_rot = DIRECT_MULTIDIM_ELEM(packed, idx++);
    sigma2_tilt = DIRECT_MULTIDIM_ELEM(packed, idx++);
    sigma2_psi = DIRECT_MULTIDIM_ELEM(packed, idx++);

    for (int igroup = 0; igroup < nr_groups; igroup++)
    {
    	sigma2_noise[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
        {
        	DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
        wsum_signal_product_spectra[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_signal_product_spectra[igroup])
        {
        	DIRECT_MULTIDIM_ELEM(wsum_signal_product_spectra[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
        wsum_reference_power_spectra[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_reference_power_spectra[igroup])
        {
        	DIRECT_MULTIDIM_ELEM(wsum_reference_power_spectra[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
        sumw_group[igroup] = DIRECT_MULTIDIM_ELEM(packed, idx++);
    }

    for (int iclass = 0; iclass < nr_classes * nr_bodies; iclass++)
    {
    	BPref[iclass].initialiseDataAndWeight(current_size);
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].data)
        {
    		(DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).real = DIRECT_MULTIDIM_ELEM(packed, idx++);
    		(DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).imag = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].weight)
        {
    		DIRECT_MULTIDIM_ELEM(BPref[iclass].weight, n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
    }
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
    	pdf_direction[iclass].resize(nr_directions);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdf_direction[iclass])
        {
        	DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
        pdf_class[iclass] = DIRECT_MULTIDIM_ELEM(packed, idx++);

        if (ref_dim==2)
        {
        	XX(prior_offset_class[iclass]) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        	YY(prior_offset_class[iclass]) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
    }

    unsigned long long packed_size = MULTIDIM_SIZE(packed);
    packed.clear();

    // Just to check whether we went outside our memory...
    if (idx != packed_size)
    {
       	std::cerr << "idx= " << idx << " packed_size= " << packed_size << std::endl;
        REPORT_ERROR("MlWsumModel::unpack: idx != idx_stop-idx_start");
    }
}


void MlWsumModel::pack(MultidimArray<RFLOAT> &packed, int &piece, int &nr_pieces, bool do_clear)
{


    // Determine size of the packed array
    int nr_groups = sigma2_noise.size();
    int nr_classes_bodies = BPref.size();
    int nr_classes = pdf_class.size();
    int spectral_size = (ori_size / 2) + 1;
    unsigned long long packed_size = 0;
    unsigned long long idx_start, idx_stop;

	// for LL & avePmax & sigma2_offset & avg_norm_correction & sigma2_rot & sigma2_tilt & sigma2_psi
    packed_size += 7 ;
    // for all group-related stuff
    packed_size += nr_groups * spectral_size;
    packed_size += nr_groups * spectral_size;
    packed_size += nr_groups * spectral_size;
    // for sumw_group
    packed_size += nr_groups;
    // for all class-related stuff
    // data is complex: multiply by two!
    packed_size += nr_classes_bodies * 2 * BPref[0].getSize();
    packed_size += nr_classes_bodies * BPref[0].getSize();
    packed_size += nr_classes * nr_directions;
    // for pdf_class
    packed_size += nr_classes;
    // for priors for each class
    if (ref_dim==2)
    	packed_size += nr_classes*2;

    if (piece < 0 && nr_pieces < 0)
    {
    	// Special case: prevent making multiple pieces if input piece and nr_pieces are both negative
        idx_start = 0;
        idx_stop = packed_size;
    }
    else if (packed_size > MAX_PACK_SIZE)
    {
        idx_start = (unsigned long long)piece * MAX_PACK_SIZE;
        idx_stop = XMIPP_MIN(idx_start + MAX_PACK_SIZE, packed_size);
        nr_pieces = CEIL((RFLOAT)packed_size/(RFLOAT)MAX_PACK_SIZE);
    }
    else
    {
        idx_start = 0;
        idx_stop = packed_size;
        nr_pieces = 1;
    }

    // increment piece so that pack will be called again
    piece++;
#ifdef DEBUG_PACK
    std::cerr << " PACK: idx_start= " << idx_start << " idx_stop= " << idx_stop << " piece= " << piece << " nr_pieces= " << nr_pieces <<" packed_size= "<<packed_size<< std::endl;
    std::cerr << " nr_classes= " << nr_classes << " nr_groups= " << nr_groups << " packed_size= " << packed_size << std::endl;
    std::cerr << " MULTIDIM_SIZE(sigma2_noise[0])= " << MULTIDIM_SIZE(sigma2_noise[0]) << " MULTIDIM_SIZE(wsum_signal_product_spectra[0])= " << MULTIDIM_SIZE(wsum_signal_product_spectra[0]) << " MULTIDIM_SIZE(wsum_reference_power_spectra[0])= " << MULTIDIM_SIZE(wsum_reference_power_spectra[0]) << std::endl;
    std::cerr << " sigma2_noise.size()= " << sigma2_noise.size() << " wsum_signal_product_spectra.size()= " << wsum_signal_product_spectra.size() << " wsum_signal_product_spectra.size()= " << wsum_signal_product_spectra.size() << std::endl;
    std::cerr << " MULTIDIM_SIZE(pdf_direction[0])= " << MULTIDIM_SIZE(pdf_direction[0]) << " pdf_direction.size()= " << pdf_direction.size()<<std::endl;
#endif

    // Get memory for the packed array
    packed.clear();
    packed.resize(idx_stop - idx_start);

    unsigned long long idx = 0;
    unsigned long long ori_idx = 0;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = LL;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = ave_Pmax;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_offset;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = avg_norm_correction;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_rot;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_tilt;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = sigma2_psi;
    ori_idx++;

    for (int igroup = 0; igroup < nr_groups; igroup++)
    {
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n);
            ori_idx++;
        }
    	if (idx == ori_idx && do_clear)
            sigma2_noise[igroup].clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_signal_product_spectra[igroup])
        {
        	if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(wsum_signal_product_spectra[igroup], n);
        	ori_idx++;
        }
        if (idx == ori_idx && do_clear)
            wsum_signal_product_spectra[igroup].clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_reference_power_spectra[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(wsum_reference_power_spectra[igroup], n);
            ori_idx++;
        }
        if (idx == ori_idx && do_clear)
            wsum_reference_power_spectra[igroup].clear();

        if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = sumw_group[igroup];
        ori_idx++;

    }
    for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
    {
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].data)
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = (DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).real;
            ori_idx++;
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = (DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).imag;
            ori_idx++;
        }
        // Only clear after the whole array has been packed... i.e. not when we reached the pack_size halfway through
        if (idx == ori_idx && do_clear)
            BPref[iclass].data.clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].weight)
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(BPref[iclass].weight, n);
            ori_idx++;
        }
        if (idx == ori_idx && do_clear)
            BPref[iclass].weight.clear();
    }

    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdf_direction[iclass])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], n);
            ori_idx++;
        }
        if (idx == ori_idx && do_clear)
        	pdf_direction[iclass].clear();

        if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = pdf_class[iclass];
        ori_idx++;

        if (ref_dim==2)
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = XX(prior_offset_class[iclass]);
            ori_idx++;
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = YY(prior_offset_class[iclass]);
            ori_idx++;
        }
    }
#ifdef DEBUG_PACK
    std::cerr << " idx= " << idx << " packed_size= " << packed_size << std::endl;
#endif

    // Just to check whether we went outside our memory...
    //std::cerr << " PACK piece= " << piece-1 << " nr_pieces= " << nr_pieces << " ori_idx= " << ori_idx<< " packed_size= " << packed_size << std::endl;
    //std::cerr << " PACK idx= " << idx << " idx_stop-idx_start= " << idx_stop-idx_start << " idx_start= " << idx_start << " idx_stop= " << idx_stop    << std::endl;
    if (idx != idx_stop-idx_start)
    {
       	std::cerr << "idx= " << idx << "ori_idx= " << ori_idx << " idx_start= " << idx_start << " idx_stop= " << idx_stop << " packed_size= " << packed_size << std::endl;
        REPORT_ERROR("MlWsumModel::pack: idx != idx_stop-idx_start");

    }

}

void MlWsumModel::unpack(MultidimArray<RFLOAT> &packed, int piece, bool do_clear)
{


    int nr_groups = sigma2_noise.size();
    int nr_classes_bodies = BPref.size();
    int nr_classes = pdf_class.size();
    int spectral_size = (ori_size / 2) + 1;
    unsigned long long idx_start;
    unsigned long long idx_stop;
    if (piece < 0)
    {
    	// Special case: prevent making multiple pieces if input piece is negative
        idx_start = 0;
        idx_stop  = MULTIDIM_SIZE(packed);
    }
    else
    {
    	idx_start = (unsigned long long)piece * MAX_PACK_SIZE;
    	idx_stop  = idx_start + (unsigned long long)MULTIDIM_SIZE(packed);
    }
    unsigned long long ori_idx = 0;
    unsigned long long idx = 0;
#ifdef DEBUG_PACK
    std::cerr << " UNPACK piece= " << piece << " idx_start= " << idx_start << " idx_stop= " << idx_stop << std::endl;
#endif

    if (ori_idx >= idx_start && ori_idx < idx_stop) LL = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) ave_Pmax = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) sigma2_offset = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) avg_norm_correction = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) sigma2_rot = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) sigma2_tilt = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) sigma2_psi = DIRECT_MULTIDIM_ELEM(packed, idx++);
    ori_idx++;

    for (int igroup = 0; igroup < nr_groups; igroup++)
    {

    	if (idx == ori_idx)
    		sigma2_noise[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
            	DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
            ori_idx++;
        }

        if (idx == ori_idx)
    		wsum_signal_product_spectra[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_signal_product_spectra[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
            	DIRECT_MULTIDIM_ELEM(wsum_signal_product_spectra[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
            ori_idx++;
        }

        if (idx == ori_idx)
    		wsum_reference_power_spectra[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(wsum_reference_power_spectra[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
            	DIRECT_MULTIDIM_ELEM(wsum_reference_power_spectra[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
            ori_idx++;
        }

        if (ori_idx >= idx_start && ori_idx < idx_stop)
        	sumw_group[igroup] = DIRECT_MULTIDIM_ELEM(packed, idx++);
        ori_idx++;

    }

    for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
    {
    	if (idx == ori_idx)
    		BPref[iclass].initialiseDataAndWeight(current_size);
    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].data)
        {
        	if (ori_idx >= idx_start && ori_idx < idx_stop)
				(DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).real = DIRECT_MULTIDIM_ELEM(packed, idx++);
        	ori_idx++;

        	if (ori_idx >= idx_start && ori_idx < idx_stop)
            	(DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).imag = DIRECT_MULTIDIM_ELEM(packed, idx++);
        	ori_idx++;
            //DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n) = Complex(re, im);
        }

    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].weight)
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
				DIRECT_MULTIDIM_ELEM(BPref[iclass].weight, n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
            ori_idx++;
        }
    }

    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
    	if (idx == ori_idx)
    		pdf_direction[iclass].resize(nr_directions);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdf_direction[iclass])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
				DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
            ori_idx++;
        }

        if (ori_idx >= idx_start && ori_idx < idx_stop)
        	pdf_class[iclass] = DIRECT_MULTIDIM_ELEM(packed, idx++);
        ori_idx++;

        if (ref_dim == 2)
        {
			if (ori_idx >= idx_start && ori_idx < idx_stop)
				XX(prior_offset_class[iclass]) = DIRECT_MULTIDIM_ELEM(packed, idx++);
			ori_idx++;
			if (ori_idx >= idx_start && ori_idx < idx_stop)
				YY(prior_offset_class[iclass]) = DIRECT_MULTIDIM_ELEM(packed, idx++);
			ori_idx++;
        }
    }


    unsigned long long packed_size = MULTIDIM_SIZE(packed);
    // Free memory
    if (do_clear)
        packed.clear();

    // Just to check whether we went outside our memory...
    //std::cerr << " UNPACK piece= " << piece << " idx= " << idx << " idx_stop-idx_start= " << idx_stop-idx_start << " idx_start= " << idx_start << " idx_stop= " << idx_stop    << std::endl;
    if (idx != idx_stop-idx_start)
    {
       	std::cerr << "idx= " << idx << "ori_idx= " << ori_idx << " idx_start= " << idx_start << " idx_stop= " << idx_stop << " packed_size= " << packed_size << std::endl;
        REPORT_ERROR("MlWsumModel::unpack: idx != idx_stop-idx_start");
    }


}



