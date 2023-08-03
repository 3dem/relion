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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include "src/ml_model.h"

#define MOM2_INIT_CONSTANT 1

#ifdef MDL_TIMING
	Timer mdl_timer;
	int TIMING_MDL_1 = proj_timer.setNew("MDL_1");
#define TIMING_TOC(id) mdl_timer.toc(id)
#else
#define TIMING_TIC(id)
#define TIMING_TOC(id)
#endif

void MlModel::initialise(bool _do_grad, bool _pseudo_halfsets)
{

	// Auxiliary vector with relevant size in Fourier space
	MultidimArray<RFLOAT > aux;
	aux.initZeros(ori_size / 2 + 1);

	// Now resize all relevant vectors
	Iref.resize(nr_classes * nr_bodies);
	masks_bodies.resize(nr_bodies);
	com_bodies.resize(nr_bodies);
	rotate_direction_bodies.resize(nr_bodies);
	orient_bodies.resize(nr_bodies);
	sigma_tilt_bodies.resize(nr_bodies, 0.);
	sigma_psi_bodies.resize(nr_bodies, 0.);
	sigma_offset_bodies.resize(nr_bodies, 0.);
	keep_fixed_bodies.resize(nr_bodies, 0);
	pointer_body_overlap.resize(nr_bodies, nr_bodies);
	max_radius_mask_bodies.resize(nr_bodies, -1);
	pdf_class.resize(nr_classes, 1./(RFLOAT)nr_classes);
	class_age.resize(nr_classes, 0);
	pdf_direction.resize(nr_classes * nr_bodies);
	group_names.resize(nr_groups, "");
	sigma2_noise.resize(nr_optics_groups);
	nr_particles_per_group.resize(nr_groups);
	tau2_class.resize(nr_classes * nr_bodies, aux);
	fsc_halves_class.resize(nr_classes * nr_bodies, aux);
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
	if(nr_classes != 1 && nr_bodies !=1)
		REPORT_ERROR("MlModel::initialise() - nr_bodies or nr_classes must be 1");
	PPref.resize(nr_classes * nr_bodies, ref);

	do_grad = _do_grad;
	pseudo_halfsets = _pseudo_halfsets;
	if (_do_grad) {
		if (_pseudo_halfsets)
			Igrad1.resize(2 * nr_classes);
		else
			Igrad1.resize(nr_classes);
		Igrad2.resize(nr_classes);
	}

	ref_names.resize(nr_classes * nr_bodies);
}

// Reading from a file
void MlModel::read(FileName fn_in, int nr_optics_groups_from_mydata, bool _do_grad, bool _pseudo_halfsets)
{

	// Clear current model
	clear();

	// Open input file
	std::ifstream in(fn_in.data(), std::ios_base::in);
	if (in.fail())
		REPORT_ERROR( (std::string) "MlModel::readStar: File " + fn_in + " cannot be read." );

	MetaDataTable MDclass, MDgroup, MDopticsgroup, MDlog, MDsigma, MDbodies;

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
	    !MDlog.getValue(EMDL_MLMODEL_PRIOR_MODE, orientational_prior_mode) ||
	    !MDlog.getValue(EMDL_MLMODEL_SIGMA_ROT, sigma2_rot) ||
	    !MDlog.getValue(EMDL_MLMODEL_SIGMA_TILT, sigma2_tilt) ||
	    !MDlog.getValue(EMDL_MLMODEL_SIGMA_PSI, sigma2_psi) ||
	    !MDlog.getValue(EMDL_MLMODEL_LL, LL) ||
	    !MDlog.getValue(EMDL_MLMODEL_AVE_PMAX, ave_Pmax) )
		REPORT_ERROR("MlModel::readStar: incorrect model_general table");

	if (!MDlog.getValue(EMDL_MLMODEL_SIGMA_OFFSET_ANGSTROM, sigma2_offset))
	{
		if (MDlog.getValue(EMDL_MLMODEL_SIGMA_OFFSET, sigma2_offset))
		{
			sigma2_offset *= pixel_size;
		}
		else
		{
			REPORT_ERROR("MlModel::readStar: incorrect model_general table: cannot find sigma_offset");
		}
	}

	bool is_pre4 = false;
	if (!MDlog.getValue(EMDL_MLMODEL_NR_OPTICS_GROUPS, nr_optics_groups))
	{
		// Just set to one here, and then confirm later in ml_optimiser where also mydata is available whether this was correct
		nr_optics_groups = nr_optics_groups_from_mydata;
		is_pre4 = true;
	}
	else
	{
		if (nr_optics_groups != nr_optics_groups_from_mydata)
		{
			REPORT_ERROR("ERROR: nr_optics_group from data.star does not match that in model.star!");
		}
	}

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
	initialise(_do_grad, _pseudo_halfsets);

	// Read classes
	Image<RFLOAT> img;
	if (nr_bodies > 1)
		MDclass.readStar(in, "model_bodies");
	else
		MDclass.readStar(in, "model_classes");

	int iclass = 0;
	do_grad = false;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDclass)
	{
		if (!MDclass.getValue(EMDL_MLMODEL_ACCURACY_TRANS_ANGSTROM, acc_trans[iclass]))
		{
			if (MDclass.getValue(EMDL_MLMODEL_ACCURACY_TRANS, acc_trans[iclass]))
			{
				acc_trans[iclass] *= pixel_size;
			}
			else
			{
				REPORT_ERROR("MlModel::readStar: incorrect model_classes/bodies table: no acc_trans");
			}
		}

		if (!MDclass.getValue(EMDL_MLMODEL_REF_IMAGE, ref_names[iclass]) ||
		    !MDclass.getValue(EMDL_MLMODEL_ACCURACY_ROT, acc_rot[iclass]) )
			REPORT_ERROR("MlModel::readStar: incorrect model_classes/bodies table: no ref_image or acc_rot");
		// backwards compatible
		if (!MDclass.getValue(EMDL_MLMODEL_ESTIM_RESOL_REF, estimated_resolution[iclass]))
			estimated_resolution[iclass] = 0.;
		if (!MDclass.getValue(EMDL_MLMODEL_FOURIER_COVERAGE_TOTAL_REF, total_fourier_coverage[iclass]))
			total_fourier_coverage[iclass] = 0.;
		if (ref_dim==2)
			if (!MDclass.getValue(EMDL_MLMODEL_PRIOR_OFFX_CLASS, XX(prior_offset_class[iclass])) ||
			    !MDclass.getValue(EMDL_MLMODEL_PRIOR_OFFY_CLASS, YY(prior_offset_class[iclass])) )
				REPORT_ERROR("MlModel::readStar: incorrect model_classes/bodies table: no offset priors for 2D classes");
		if (iclass == 0 || nr_bodies == 1) // there is only one pdf_class for multibody, but multiple for classification!
			if (!MDclass.getValue(EMDL_MLMODEL_PDF_CLASS, pdf_class[iclass]) )
				REPORT_ERROR("MlModel::readStar: incorrect model_classes table: no pdf_class");
		if (is_helix)
		{
			if (!MDclass.getValue(EMDL_MLMODEL_HELICAL_RISE, helical_rise[iclass]) ||
			    !MDclass.getValue(EMDL_MLMODEL_HELICAL_TWIST, helical_twist[iclass]) )
				REPORT_ERROR("MlModel::readStar: incorrect helical parameters");
		}
		if (nr_bodies > 1)
		{
			if (MDclass.containsLabel(EMDL_BODY_KEEP_FIXED))
				MDclass.getValue(EMDL_BODY_KEEP_FIXED, keep_fixed_bodies[iclass]);
			else
				keep_fixed_bodies[iclass] = 0;
		}

		// Read in actual reference image
		img.read(ref_names[iclass]);
		img().setXmippOrigin();
		Iref[iclass] = img();

		FileName fn_tmp;
		// Check to see whether there are gradient tracking entry as well
		if (MDclass.getValue(EMDL_MLMODEL_GRADIENT_MOMENT1_IMAGE, fn_tmp))
		{
			bool is_2d(Iref[0].zdim == 1);
			Image<RFLOAT> img;
			do_grad = true;
			if (iclass == 0)
			{
				if (pseudo_halfsets)
					Igrad1.resize(2 * nr_classes);
				else
					Igrad1.resize(nr_classes);
			}
			img.read(fn_tmp);

			Igrad1[iclass].resize(
			        is_2d ? 1: Iref[0].zdim * padding_factor,
					Iref[0].ydim * padding_factor,
					Iref[0].xdim * padding_factor/2+1
			);

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Igrad1[iclass]) {
				DIRECT_MULTIDIM_ELEM(Igrad1[iclass], n).real = DIRECT_MULTIDIM_ELEM(img(), n * 2 + 0);
				DIRECT_MULTIDIM_ELEM(Igrad1[iclass], n).imag = DIRECT_MULTIDIM_ELEM(img(), n * 2 + 1);
			}

			if (pseudo_halfsets)
			{
				FileName fileName;

				long int img_idx;
				fn_tmp.decompose(img_idx, fileName);

				if (img_idx >= 0)
					fileName.compose(img_idx + nr_classes, fileName);
				else
				{
					fileName = fn_tmp.withoutExtension();
					img_idx = atoi( fileName.substr(fileName.length() - 3, 3).c_str() );
					fileName = fileName.substr(0, fileName.length() - 3);
					fileName.compose(fileName, img_idx + nr_classes, "mrc", 3);
				}
				img.read(fileName);

				Igrad1[iclass + nr_classes].resize(
						is_2d ? 1 : Iref[0].zdim * padding_factor,
						Iref[0].ydim * padding_factor,
						Iref[0].xdim * padding_factor / 2 + 1
				);

				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Igrad1[iclass + nr_classes]) {
					DIRECT_MULTIDIM_ELEM(Igrad1[iclass + nr_classes], n).real = DIRECT_MULTIDIM_ELEM(img(), n * 2 + 0);
					DIRECT_MULTIDIM_ELEM(Igrad1[iclass + nr_classes], n).imag = DIRECT_MULTIDIM_ELEM(img(), n * 2 + 1);
				}
			}
		}

		if (MDclass.getValue(EMDL_MLMODEL_GRADIENT_MOMENT2_IMAGE, fn_tmp))
		{
			Image<RFLOAT> img;
			do_grad = true;
			if (iclass == 0)
				Igrad2.resize(nr_classes);
			img.read(fn_tmp);

			Igrad2[iclass].resize(
					Iref[0].zdim == 1 ? 1: Iref[0].zdim * padding_factor,
					Iref[0].ydim * padding_factor,
					Iref[0].xdim * padding_factor/2+1
			);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Igrad2[iclass]) {
				DIRECT_MULTIDIM_ELEM(Igrad2[iclass], n).real = DIRECT_MULTIDIM_ELEM(img(), n * 2 + 0);
				DIRECT_MULTIDIM_ELEM(Igrad2[iclass], n).imag = DIRECT_MULTIDIM_ELEM(img(), n * 2 + 1);
			}
		}

		iclass++;
	}

	// Read group stuff
	MDgroup.readStar(in, "model_groups");
	long int igroup;
	long long nr_particles = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDgroup)
	{
		if (!MDgroup.getValue(EMDL_MLMODEL_GROUP_NO, igroup))
		{
			REPORT_ERROR("MlModel::readStar: incorrect model_groups table");
		}
		//Start counting of groups at 1, not at 0....
		if (!MDgroup.getValue(EMDL_MLMODEL_GROUP_SCALE_CORRECTION, scale_correction[igroup - 1]) ||
			!MDgroup.getValue(EMDL_MLMODEL_GROUP_NR_PARTICLES, nr_particles_per_group[igroup - 1]) ||
		    !MDgroup.getValue(EMDL_MLMODEL_GROUP_NAME, group_names[igroup-1]))
			REPORT_ERROR("MlModel::readStar: incorrect model_groups table");
		nr_particles += nr_particles_per_group[igroup - 1];
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
				REPORT_ERROR("MlModel::readStar: incorrect table model_class/body_"+integerToString(iclass + 1));
			if (!MDsigma.getValue(EMDL_MLMODEL_DATA_VS_PRIOR_REF, data_vs_prior_class[iclass](idx)) ||
			    !MDsigma.getValue(EMDL_MLMODEL_TAU2_REF, tau2_class[iclass](idx)) ||
			    !MDsigma.getValue(EMDL_MLMODEL_FSC_HALVES_REF, fsc_halves_class[iclass](idx)) ||
			    !MDsigma.getValue(EMDL_MLMODEL_SIGMA2_REF, sigma2_class[iclass](idx)))
				REPORT_ERROR("MlModel::readStar: incorrect table model_class/body_"+integerToString(iclass + 1));
			// backwards compatible with STAR files without Fourier coverage
			if (!MDsigma.getValue(EMDL_MLMODEL_FOURIER_COVERAGE_REF, fourier_coverage_class[iclass](idx)))
				fourier_coverage_class[iclass](idx) = 0.;
		}
	}

	// Read sigma2_noise models for each optics group
	for (int igroup = 0; igroup < nr_optics_groups; igroup++)
	{
		// Allow sigma2_noise with different sizes!
		sigma2_noise[igroup].resize(ori_size/2 + 1);

		if (is_pre4)
		{
			// If this is a pre-relion-4.0 model.star file, then just read the first sigma2_noise spectrum into all optics_groups
			MDsigma.readStar(in, "model_group_1");
		}
		else
		{
			MDsigma.readStar(in, "model_optics_group_" + integerToString(igroup + 1));
		}
		int idx;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsigma)
		{
			if (!MDsigma.getValue(EMDL_SPECTRAL_IDX, idx))
				REPORT_ERROR("MlModel::readStar: incorrect table model_group_" + integerToString(igroup + 1));
			if (!MDsigma.getValue(EMDL_MLMODEL_SIGMA2_NOISE, sigma2_noise[igroup](idx)))
				REPORT_ERROR("MlModel::readStar: incorrect table model_group_" + integerToString(igroup + 1));
		}

		// For empty optics groups (not so likely, but who knows...)
		if (MDsigma.numberOfObjects() == 0)
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
		for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
		{
			if (nr_bodies > 1)
				MDclass.readStar(in, "model_pdf_orient_body_" + integerToString(iclass + 1));
			else
				MDclass.readStar(in, "model_pdf_orient_class_" + integerToString(iclass + 1));
			pdf_direction[iclass].clear();
			RFLOAT aux;
			std::vector<RFLOAT> vaux;
			vaux.clear();
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDclass)
			{
				if (!MDclass.getValue(EMDL_MLMODEL_PDF_ORIENT, aux))
					REPORT_ERROR("MlModel::readStar: incorrect table model_pdf_orient_class_" + integerToString(iclass + 1));
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
		for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
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

void MlModel::write(FileName fn_out, HealpixSampling &sampling, bool do_write_bild, bool only_write_images)
{

	MetaDataTable MDclass, MDgroup, MDopticsgroup, MDlog, MDsigma, MDbodies;
	FileName fn_tmp, fn_tmp2, fn_mom1, fn_mom2;
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
		img.setSamplingRateInHeader(pixel_size);
		if (nr_bodies > 1)
			img.write(fn_out + "_bodies.mrcs");
		else
			img.write(fn_out + "_classes.mrcs");

		if (do_grad)
		{
			int nr_grads = pseudo_halfsets ? nr_classes_bodies * 2 : nr_classes_bodies;
			Image<RFLOAT> img1(XSIZE(Igrad1[0])*2, YSIZE(Igrad1[0]), 1, nr_grads);
			for (int iclass = 0; iclass < nr_classes; iclass++)
			{
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Igrad1[iclass])
				{
					DIRECT_NZYX_ELEM(img1(), iclass, 0, i, j*2+0) = DIRECT_A2D_ELEM(Igrad1[iclass], i, j).real;
					DIRECT_NZYX_ELEM(img1(), iclass, 0, i, j*2+1) = DIRECT_A2D_ELEM(Igrad1[iclass], i, j).imag;

				}
				if (pseudo_halfsets)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Igrad1[iclass])
					{
						DIRECT_NZYX_ELEM(img1(), iclass + nr_classes, 0, i, j*2+0) = DIRECT_A2D_ELEM(Igrad1[iclass + nr_classes], i, j).real;
						DIRECT_NZYX_ELEM(img1(), iclass + nr_classes, 0, i, j*2+1) = DIRECT_A2D_ELEM(Igrad1[iclass + nr_classes], i, j).imag;
					}
				}
			}
			img1.write(fn_out + "_1moment.mrcs");

			Image<RFLOAT> img2(XSIZE(Igrad1[0])*2, YSIZE(Igrad1[0]), 1, nr_classes_bodies);
			for (int iclass = 0; iclass < nr_classes; iclass++)
			{
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Igrad2[iclass])
					{
						DIRECT_NZYX_ELEM(img2(), iclass, 0, i, j*2+0) = DIRECT_A2D_ELEM(Igrad2[iclass], i, j).real;
						DIRECT_NZYX_ELEM(img2(), iclass, 0, i, j*2+1) = DIRECT_A2D_ELEM(Igrad2[iclass], i, j).imag;
					}

			}
			img2.write(fn_out + "_2moment.mrcs");
		}
	}
	else
	{
		Image<RFLOAT> img;
		// Set correct voxel size in the header
		for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
		{
			img() = Iref[iclass];
			img.setSamplingRateInHeader(pixel_size);
			if (nr_bodies > 1)
			{
				fn_tmp.compose(fn_out+"_body", iclass+1, "mrc", 3);
				// apply the body mask for output to the user
				// No! That interferes with a clean continuation of multibody refinement, as ref will be masked 2x then!
				// img() *= masks_bodies[iclass];
			}
			else
				fn_tmp.compose(fn_out+"_class", iclass+1, "mrc", 3);

			img.write(fn_tmp);
		}

		if (do_grad)
		{
			for (int iclass = 0; iclass < nr_classes; iclass++)
			{
				fn_tmp.compose(fn_out+"_1moment", iclass+1, "mrc", 3);

				Image<RFLOAT> img(XSIZE(Igrad1[0])*2, YSIZE(Igrad1[0]), ZSIZE(Igrad1[0]));
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Igrad1[iclass]) {
					DIRECT_A3D_ELEM(img(), k, i, j*2+0) = DIRECT_A3D_ELEM(Igrad1[iclass], k, i, j).real;
					DIRECT_A3D_ELEM(img(), k, i, j*2+1) = DIRECT_A3D_ELEM(Igrad1[iclass], k, i, j).imag;
				}
				img.write(fn_tmp);

				if (pseudo_halfsets)
				{
					fn_tmp.compose(fn_out+"_1moment", iclass+1+nr_classes, "mrc", 3);

					Image<RFLOAT> img(XSIZE(Igrad1[0])*2, YSIZE(Igrad1[0]), ZSIZE(Igrad1[0]));
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Igrad1[iclass+nr_classes]) {
								DIRECT_A3D_ELEM(img(), k, i, j*2+0) = DIRECT_A3D_ELEM(Igrad1[iclass+nr_classes], k, i, j).real;
								DIRECT_A3D_ELEM(img(), k, i, j*2+1) = DIRECT_A3D_ELEM(Igrad1[iclass+nr_classes], k, i, j).imag;
							}
					img.write(fn_tmp);
				}

				fn_tmp.compose(fn_out+"_2moment", iclass+1, "mrc", 3);
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Igrad2[iclass]) {
							DIRECT_A3D_ELEM(img(), k, i, j*2+0) = DIRECT_A3D_ELEM(Igrad2[iclass], k, i, j).real;
							DIRECT_A3D_ELEM(img(), k, i, j*2+1) = DIRECT_A3D_ELEM(Igrad2[iclass], k, i, j).imag;
						}
				img.write(fn_tmp);
			}
		}

		if (do_write_bild)
		{
			// Also write out bild files with the orientational distribution of each class
			// Also write out angular distributions
			// Don't do this for bodies, only for classes!
			for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
			{
				FileName fn_bild;
				if (nr_bodies > 1)
					fn_bild.compose(fn_out+"_body",iclass+1,"", 3);
				else
					fn_bild.compose(fn_out+"_class",iclass+1,"", 3);
				fn_bild += "_angdist.bild";
				RFLOAT offset = ori_size * pixel_size / 2.;
				if (nr_bodies > 1)
				{
					// 14jul2017: rotations are all relative to (rot,tilt)=(0,90) to prevent problems with psi-prior around  tilt=0!
					sampling.writeBildFileOrientationalDistribution(pdf_direction[iclass], fn_bild, offset, offset,
					                                                &orient_bodies[iclass], &com_bodies[iclass]);
				}
				else
				{
					sampling.writeBildFileOrientationalDistribution(pdf_direction[iclass], fn_bild, offset, offset);
				}
			}
		}

	}

	if (only_write_images)
		return;

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
	MDlog.setValue(EMDL_MLMODEL_NR_OPTICS_GROUPS, nr_optics_groups);
	MDlog.setValue(EMDL_MLMODEL_TAU2_FUDGE_FACTOR, tau2_fudge_factor);
	MDlog.setValue(EMDL_MLMODEL_NORM_CORRECTION_AVG, avg_norm_correction);
	MDlog.setValue(EMDL_MLMODEL_SIGMA_OFFSET_ANGSTROM, sqrt(sigma2_offset));
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

			fn_mom1.compose(iclass + 1, fn_out + "_1moment.mrcs");
			fn_mom2.compose(iclass + 1, fn_out + "_2moment.mrcs");
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

			fn_mom1.compose(fn_out + "_1moment", iclass + 1, "mrc", 3);
			fn_mom2.compose(fn_out + "_2moment", iclass + 1, "mrc", 3);
		}
		MDclass.setValue(EMDL_MLMODEL_REF_IMAGE, fn_tmp);

		if (do_grad) {
			MDclass.setValue(EMDL_MLMODEL_GRADIENT_MOMENT1_IMAGE, fn_mom1);
			MDclass.setValue(EMDL_MLMODEL_GRADIENT_MOMENT2_IMAGE, fn_mom2);
		}

		// For multiple bodies: only star PDF_CLASS in the first one!
		int myclass = (nr_bodies > 1) ? 0 : iclass; // for multi-body: just set iclass=0

		MDclass.setValue(EMDL_MLMODEL_PDF_CLASS, pdf_class[myclass]);
		MDclass.setValue(EMDL_MLMODEL_ACCURACY_ROT, acc_rot[iclass]);
		MDclass.setValue(EMDL_MLMODEL_ACCURACY_TRANS_ANGSTROM, acc_trans[iclass]);
		MDclass.setValue(EMDL_MLMODEL_ESTIM_RESOL_REF, estimated_resolution[iclass]);
		MDclass.setValue(EMDL_MLMODEL_FOURIER_COVERAGE_TOTAL_REF, total_fourier_coverage[iclass]);
		if (nr_bodies > 1)
		{
			MDclass.setValue(EMDL_BODY_ROTATE_DIRECTION_X, XX(rotate_direction_bodies[iclass]));
			MDclass.setValue(EMDL_BODY_ROTATE_DIRECTION_Y, YY(rotate_direction_bodies[iclass]));
			MDclass.setValue(EMDL_BODY_ROTATE_DIRECTION_Z, ZZ(rotate_direction_bodies[iclass]));
			MDclass.setValue(EMDL_BODY_KEEP_FIXED, keep_fixed_bodies[iclass]);
		}

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
			MDsigma.setValue(EMDL_MLMODEL_FSC_HALVES_REF, fsc_halves_class[iclass](ii));
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
		MDgroup.setValue(EMDL_MLMODEL_GROUP_NR_PARTICLES, nr_particles_per_group[igroup]);
		MDgroup.setValue(EMDL_MLMODEL_GROUP_SCALE_CORRECTION, scale_correction[igroup]);
	}
	MDgroup.write(fh);

	// Write sigma models for each optics group
	for (int igroup = 0; igroup < nr_optics_groups; igroup++)
	{
		MDsigma.clear();
		MDsigma.setName("model_optics_group_"+integerToString(igroup+1));
		for (int ii = 0; ii < XSIZE(sigma2_noise[igroup]); ii++)
		{
			// Some points in sigma2_noise arrays are never used...
			aux = sigma2_noise[igroup](ii);
			if (aux > 0.)
			{
				MDsigma.addObject();
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
		for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
		{
			MDclass.clear();
			if (nr_bodies > 1)
				MDclass.setName("model_pdf_orient_body_"+integerToString(iclass+1));
			else
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
void MlModel::initialiseFromImages(
	FileName fn_ref, bool _is_3d_model, Experiment &_mydata,
	bool &do_average_unaligned, bool &do_generate_seeds, bool &refs_are_ctf_corrected,
	RFLOAT _ref_angpix, bool _do_grad, bool _pseudo_halfsets, bool _do_trust_ref_size, bool verb)
{


	// Data dimensionality
	if (!_mydata.obsModel.opticsMdt.containsLabel(EMDL_IMAGE_DIMENSIONALITY))
	{
		if (verb > 0) std::cerr << " WARNING: input particles STAR file does not have a column for image dimensionality, assuming 2D images ..." << std::endl;
		data_dim = 2;
	}
	else
	{
		_mydata.obsModel.opticsMdt.getValue(EMDL_IMAGE_DIMENSIONALITY, data_dim, 0);
	}

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
			Igrad1.clear();
			Igrad2.clear();

			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDref)
			{
				MDref.getValue(EMDL_MLMODEL_REF_IMAGE, fn_tmp);
				img.read(fn_tmp);
				img().setXmippOrigin();
				if (_ref_angpix > 0.)
				{
					pixel_size = _ref_angpix;
				}
				else
				{
					RFLOAT header_pixel_size;
					if (nr_classes == 0)
					{
						img.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_X, header_pixel_size);
						pixel_size = header_pixel_size;
					}
					else
					{
						img.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_X, header_pixel_size);
						if (fabs(header_pixel_size - pixel_size) > 0.001)
						{
							REPORT_ERROR("MlModel::readImages ERROR: different models have different pixel sizes in their headers!");
						}
					}
				}

				ori_size = XSIZE(img());
				ref_dim = img().getDim();
				Iref.push_back(img());

				if (_do_grad)
				{
					MultidimArray<Complex> zeros(
							Iref[0].zdim == 1 ? 1: Iref[0].zdim * padding_factor,
							img().ydim * padding_factor,
							img().xdim  * padding_factor / 2 + 1
					);
					zeros.initZeros();
					MultidimArray<Complex> constv(zeros);
					constv.initConstant(Complex(MOM2_INIT_CONSTANT, MOM2_INIT_CONSTANT));

					Igrad1.push_back(zeros);
					if (_pseudo_halfsets)
						Igrad1.push_back(zeros);

					Igrad2.push_back(constv); // Mom2 init value
				}
				nr_classes++;
			}
		}
		// For a single image, read this image as reference and set it in all nr_classes Irefs
		else
		{
			img.read(fn_ref);
			img().setXmippOrigin();
			if (_ref_angpix > 0.)
			{
				pixel_size = _ref_angpix;
			}
			else
			{
				RFLOAT header_pixel_size;
				img.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_X, header_pixel_size);
				if (header_pixel_size <= 0)
				{
					std::cerr << " header_pixel_size = " << header_pixel_size << std::endl;
					REPORT_ERROR("MlModel::initialiseFromImages: Pixel size of reference image is not set!");
				}
				pixel_size = header_pixel_size;
			}
			ori_size = XSIZE(img());
			ref_dim = img().getDim();
			if (ori_size != XSIZE(img()) || ori_size != YSIZE(img()))
			{
				std::cerr << " ori_size= " << ori_size << " XSIZE(img())= " << XSIZE(img()) << std::endl;
				REPORT_ERROR("MlOptimiser::read: size of reference image is not the same as the experimental images!");
			}
			Iref.clear();
			Igrad1.clear();
			Igrad2.clear();
			if (nr_bodies > 1)
			{
				for (int ibody = 0; ibody < nr_bodies; ibody++)
				{
					Iref.push_back(img());
					if (masks_bodies.size() <= ibody)
						REPORT_ERROR("BUG: masks_bodies.size() < ibody. Did you initialise the body masks before reading the references?");
				}
			}
			else
			{
				MultidimArray<Complex> zeros(
						img().zdim * padding_factor,
						img().ydim * padding_factor,
						img().xdim  * padding_factor / 2 + 1
				);
				zeros.initZeros();
				MultidimArray<Complex> constv(zeros);
				constv.initConstant(Complex(MOM2_INIT_CONSTANT, MOM2_INIT_CONSTANT));

				for (int iclass = 0; iclass < nr_classes; iclass++)
				{
					Iref.push_back(img());

					if (_do_grad) {
						Igrad1.push_back(zeros);
						if (_pseudo_halfsets)
							Igrad1.push_back(zeros);

						Igrad2.push_back(constv); // Mom2 init value
					}
				}
			}
			if (nr_classes > 1)
				do_generate_seeds = true;
			else
				do_generate_seeds = false;
		}

	}

	// Make sure that the model has the same box and pixel size as (the first optics group of) the data
	RFLOAT pixel_size_first_optics_group = _mydata.getOpticsPixelSize(0);
	int box_size_first_optics_group = _mydata.getOpticsImageSize(0);

	if (fn_ref != "None")
	{

		if (fabs(pixel_size - pixel_size_first_optics_group) > 0.001 ||
		    ori_size != box_size_first_optics_group)
		{

			std::string mesg = "";
			if (fabs(pixel_size - pixel_size_first_optics_group) > 0.001)
			{
				mesg = " The reference pixel size is " + floatToString(pixel_size)
				     + " A/px, but the pixel size of the first optics group of the data is "
					 + floatToString(pixel_size_first_optics_group) + " A/px! \n";
			}
			if (ori_size != box_size_first_optics_group)
			{
				mesg += " The reference box size is " + integerToString(ori_size)
				     + " px, but the box size of the first optics group of the data is "
					 + integerToString(box_size_first_optics_group) + " px!\n";
			}

			if (!_do_trust_ref_size)
				REPORT_ERROR("ERROR " + mesg + "\nIf you want to re-scale and/or re-box input particles into the pixel size and the box size of the reference, re-run the program with the --trust_ref_size option.");
			else if (verb)
				std::cerr << " WARNING " << mesg;
		}

	}
	else
	{
		pixel_size = pixel_size_first_optics_group;
		ori_size = box_size_first_optics_group;

		// Calculate average of all unaligned images later on.
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
		img().setXmippOrigin();
		Iref.clear();
		Igrad1.clear();
		Igrad2.clear();
		for (int iclass = 0; iclass < nr_classes; iclass++)
		{
			Iref.push_back(img());

			if (_do_grad) {

				MultidimArray<Complex> zeros(
						Iref[0].zdim == 1 ? 1: Iref[0].zdim * padding_factor,
						img().ydim * padding_factor,
						img().xdim  * padding_factor / 2 + 1
				);
				zeros.initZeros();
				MultidimArray<Complex> constv(zeros);
				constv.initConstant(Complex(MOM2_INIT_CONSTANT, MOM2_INIT_CONSTANT));

				Igrad1.push_back(zeros);
				if (_pseudo_halfsets)
					Igrad1.push_back(zeros);

				Igrad2.push_back(constv); // Mom2 init value
			}
		}
	}


	// Set some group stuff
	nr_optics_groups = _mydata.numberOfOpticsGroups();
	MultidimArray<RFLOAT> aux;
	aux.initZeros(ori_size/2 + 1);
	sigma2_noise.resize(nr_optics_groups, aux);
	nr_groups = _mydata.groups.size();

	initialise(_do_grad, _pseudo_halfsets);

	for (int i=0; i< nr_groups; i++)
		group_names[i] = _mydata.groups[i].name;

}

void MlModel::initialisePdfDirection(long long int newsize)
{

	// If the pdf_direction were already filled (size!=0), and newsize=oldsize then leave them as they were
	// If they were still empty, or if the size changes, then initialise them with an even distribution
	for (int iclass = 0; iclass < nr_classes * nr_bodies; iclass++)
	{
		long long int oldsize = MULTIDIM_SIZE(pdf_direction[iclass]);
		if (oldsize == 0 || oldsize != newsize)
		{
			pdf_direction[iclass].resize(newsize);
			pdf_direction[iclass].initConstant(1./((RFLOAT) nr_classes * newsize));
		}
	}
	nr_directions = newsize;

}

void MlModel::initialiseBodies(FileName fn_masks, FileName fn_root_out, bool also_initialise_rest, int rank)
{
	MetaDataTable MD;
	MD.read(fn_masks);
	if (!MD.containsLabel(EMDL_BODY_MASK_NAME))
		REPORT_ERROR("ERROR MlModel::initialiseBodyMasks: body-mask STAR file does not contain rlnBodyMaskName label.");

	nr_bodies = 0;
	masks_bodies.resize(MD.numberOfObjects());
	com_bodies.resize(MD.numberOfObjects());
	rotate_direction_bodies.resize(MD.numberOfObjects());
	orient_bodies.resize(MD.numberOfObjects());
	sigma_tilt_bodies.resize(MD.numberOfObjects());
	sigma_psi_bodies.resize(MD.numberOfObjects());
	sigma_offset_bodies.resize(MD.numberOfObjects());
	keep_fixed_bodies.resize(MD.numberOfObjects());
	max_radius_mask_bodies.resize(MD.numberOfObjects());
	FileName fn_mask;
	Image<RFLOAT> Imask;
	std::vector<int> relatives_to;
	Matrix1D<RFLOAT> one_direction(3);
	bool has_rotate_directions = false;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		MD.getValue(EMDL_BODY_MASK_NAME, fn_mask);
		Imask.read(fn_mask);
		RFLOAT minval, maxval;
		Imask().computeDoubleMinMax(minval, maxval);
		if (minval < 0. || maxval > 1.)
			REPORT_ERROR("ERROR: the mask " + fn_mask + " has values outside the range [0,1]");

		Imask().setXmippOrigin();
		masks_bodies[nr_bodies] = Imask();
		Imask.setSamplingRateInHeader(pixel_size);
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
		max_radius_mask_bodies[nr_bodies] = CEIL(pixel_size * sqrt((RFLOAT)max_d2));

		// Get which body to rotate relative to
		int relative_to = -1;
		if (MD.containsLabel(EMDL_BODY_ROTATE_RELATIVE_TO))
		{
			MD.getValue(EMDL_BODY_ROTATE_RELATIVE_TO, relative_to);
			relative_to--;// numbering in STAR file starts with 1
		}
		relatives_to.push_back(relative_to);

		if (MD.containsLabel(EMDL_BODY_ROTATE_DIRECTION_X) &&
		    MD.containsLabel(EMDL_BODY_ROTATE_DIRECTION_Y) &&
		    MD.containsLabel(EMDL_BODY_ROTATE_DIRECTION_Z))
		{
			has_rotate_directions = true;
			MD.getValue(EMDL_BODY_ROTATE_DIRECTION_X, XX(one_direction));
			MD.getValue(EMDL_BODY_ROTATE_DIRECTION_Y, YY(one_direction));
			MD.getValue(EMDL_BODY_ROTATE_DIRECTION_Z, ZZ(one_direction));
			rotate_direction_bodies.push_back(one_direction);
		}

		RFLOAT val;
		if (MD.containsLabel(EMDL_BODY_SIGMA_ANG))
		{
			MD.getValue(EMDL_BODY_SIGMA_ANG, val);
			sigma_tilt_bodies[nr_bodies] = val;
			sigma_psi_bodies[nr_bodies] = val;
		}
		else
		{
			if (!(MD.containsLabel(EMDL_BODY_SIGMA_TILT) && MD.containsLabel(EMDL_BODY_SIGMA_PSI)) )
				REPORT_ERROR("ERROR: either provide rlnBodySigmaAngles OR provide rlnBodySigmaTilt and rlnBodySigmaPsi in the body STAR file.");
			MD.getValue(EMDL_BODY_SIGMA_TILT, val);
			sigma_tilt_bodies[nr_bodies] = val;
			MD.getValue(EMDL_BODY_SIGMA_PSI, val);
			sigma_psi_bodies[nr_bodies] = val;
		}

		if (MD.getValue(EMDL_BODY_SIGMA_OFFSET_ANGSTROM, val))
		{
			sigma_offset_bodies[nr_bodies] = val;
		}
		else if (MD.getValue(EMDL_BODY_SIGMA_OFFSET, val))
		{
			val *= pixel_size;
		}
		else
		{
			REPORT_ERROR("ERROR: the body STAR file should contain a rlnBodySigmaOffsetAngst column for the prior on the offsets for each body");
		}

		// Also write the mask with the standard name to disk
		fn_mask.compose(fn_root_out + "_body", nr_bodies + 1, "", 3); // body number from 1 to K!
		fn_mask += "_mask.mrc";

		if (rank == 0)
			Imask.write(fn_mask);

		// update counter at the end!
		nr_bodies++;
	}

	// Now that we have the COMs, also get the orientation matrix and the direction of rotation for each body
	for (int ibody = 0; ibody < nr_bodies; ibody++)
	{
		if (relatives_to[ibody] >= 0)
		{
			// If another body was given in the input STAR file, rotate this body wrt the COM of the other body
			rotate_direction_bodies[ibody] = com_bodies[relatives_to[ibody]];
			rotate_direction_bodies[ibody] -= com_bodies[ibody];
		}
		else if (has_rotate_directions)
		{
			// If the rotation vector is specified directly, just use this one
		}
		else
		{
			// if no relative-bodies, nor explicit rotation directions are specified in the STAR file, then rotate relative to (0,0,0)
			rotate_direction_bodies[ibody].initZeros();
			rotate_direction_bodies[ibody] -= com_bodies[ibody];
		}

		rotate_direction_bodies[ibody].selfNormalize();
		alignWithZ(-rotate_direction_bodies[ibody], orient_bodies[ibody], false);
	}

	if (also_initialise_rest)
	{
		if (Iref.size() != 1)
			REPORT_ERROR("BUG: at this point, there should only be a single reference!");

		for (int ibody = 1; ibody < nr_bodies; ibody++)
		{

			Iref.push_back(Iref[0]);
			tau2_class.push_back(tau2_class[0]);
			fsc_halves_class.push_back(fsc_halves_class[0]);
			sigma2_class.push_back(sigma2_class[0]);
			data_vs_prior_class.push_back(data_vs_prior_class[0]);
			fourier_coverage_class.push_back(fourier_coverage_class[0]);
			acc_rot.push_back(acc_rot[0]);
			acc_trans.push_back(acc_trans[0]);
			estimated_resolution.push_back(estimated_resolution[0]);
			total_fourier_coverage.push_back(total_fourier_coverage[0]);
			if (ref_dim==2)
				prior_offset_class.push_back(prior_offset_class[0]);
			orientability_contrib.push_back(orientability_contrib[0]);
			PPref.push_back(PPref[0]);
			pdf_direction.push_back(pdf_direction[0]);

			// If all sigmas are zero, ignore this body in the refinement
			if (sigma_tilt_bodies[ibody] < 0.001 &&
			    sigma_psi_bodies[ibody] < 0.001 &&
			    sigma_offset_bodies[ibody] < 0.001)
				keep_fixed_bodies[ibody] = 1;
			else
				keep_fixed_bodies[ibody] = 0;
		}

		// If provided a specific reference, re-set the corresponding Iref entry
		if (MD.containsLabel(EMDL_BODY_REFERENCE_NAME))
		{
			int ibody = 0;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
			{
				FileName fn_ref;
				MD.getValue(EMDL_BODY_REFERENCE_NAME, fn_ref);
				if (fn_ref != "None")
				{
					Image<RFLOAT> img;
					img.read(fn_ref);
					img().setXmippOrigin();
					Iref[ibody] = img();
				}
				ibody++;
			}
		}
	}

	// Find the overlap of the bodies, and extend the Iref, PPref and masks_bodies vectors
	pointer_body_overlap.resize(nr_bodies, nr_bodies);
	pointer_body_overlap_inv.resize(nr_bodies);

//#define DEBUG_OVERLAP
	if (norm_body_mask_overlap)
	{
		MultidimArray<RFLOAT> sum_mask = masks_bodies[0];
		for (int ibody = 1; ibody < nr_bodies; ibody++)
			sum_mask += masks_bodies[ibody];

		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(sum_mask)
			if (DIRECT_A1D_ELEM(sum_mask, i) > 1.)
				for (int ibody = 0; ibody < nr_bodies; ibody++)
					DIRECT_A1D_ELEM(masks_bodies[ibody], i) /= DIRECT_A1D_ELEM(sum_mask, i);

		for (int ibody = 0; ibody < nr_bodies; ibody++)
		{
			for (int obody = 0; obody < nr_bodies; obody++)
				DIRECT_A2D_ELEM(pointer_body_overlap, ibody, obody) = obody;
			pointer_body_overlap_inv[ibody] = ibody;
		}

#ifdef DEBUG_OVERLAP
		for (int ibody = 0; ibody < nr_bodies; ibody++)
		{
			Image<RFLOAT> It;
			It()= masks_bodies[ibody];
			fnt = "mask_ibody"+integerToString(ibody)+".spi";
			It.write(fnt);
			std::cerr << " PPref.size()= " << PPref.size() << std::endl;
		}
#endif
	}
	else
	{
		for (int ibody = 0; ibody < nr_bodies; ibody++)
		{
#ifdef DEBUG_OVERLAP
			Image<RFLOAT> It;
			FileName fnt;
			It()= masks_bodies[ibody];
			fnt = "mask_ibody"+integerToString(ibody)+".spi";
			It.write(fnt);
#endif
			for (int obody = 0; obody < nr_bodies; obody++)
			{
				if (ibody == obody)
				{
					DIRECT_A2D_ELEM(pointer_body_overlap, ibody, obody) = obody;
					pointer_body_overlap_inv[obody] = obody;
				}
				else
				{
					// Sum all the previously done obody masks to see whether there is also overlap with any of them
					MultidimArray<RFLOAT> overlap_mask = masks_bodies[ibody];
					for (int oldobody = 0; oldobody < obody; oldobody++)
					{
						if (oldobody != ibody)
						{
							int ii = DIRECT_A2D_ELEM(pointer_body_overlap, ibody, oldobody);
							overlap_mask += masks_bodies[ii];
						}
					}
					// Calculate the overlap between the sum of ibody and all the old obodies until now
					overlap_mask *= masks_bodies[obody]; // element-wise multiplication
					// If there is overlap, generate another PPref
					if (overlap_mask.sum() > 0.)
					{
						// Calculate the mask that has the overlap subtracted from the obody mask
						overlap_mask = masks_bodies[obody] - overlap_mask;
						// set the right pointer in the 2D matrix
						DIRECT_A2D_ELEM(pointer_body_overlap, ibody, obody) = PPref.size();
						//std::cerr << " ibody= " << ibody << " obody= " << obody << " overlap= " << overlap_mask.sum() << " icc= " << PPref.size() << std::endl;
						// Extend the two vectors here!
						PPref.push_back(PPref[obody]);
						masks_bodies.push_back(overlap_mask);
						// And keep track of which ibody this entry belonged to
						pointer_body_overlap_inv.push_back(obody);

#ifdef DEBUG_OVERLAP
						It()= overlap_mask;
						fnt = "mask_ibody"+integerToString(ibody)+"_obody"+integerToString(obody)+"_overlap.spi";
						It.write(fnt);
						std::cerr << " PPref.size()= " << PPref.size() << std::endl;
#endif
					}
					else
						// if there is no overlap: just point to the original obody
						DIRECT_A2D_ELEM(pointer_body_overlap, ibody, obody) = obody;
				}
			}
		}
	}
}


void MlModel::writeBildFileBodies(FileName fn_bild)
{

	std::ofstream fh_bild;
	fh_bild.open(fn_bild.c_str(), std::ios::out);
	if (!fh_bild)
		REPORT_ERROR("HealpixSampling::writeBildFileOrientationalDistribution: cannot open " + fn_bild);

	RFLOAT xcen = -STARTINGX(Iref[0]) * pixel_size;
	RFLOAT ycen = -STARTINGY(Iref[0]) * pixel_size;
	RFLOAT zcen = -STARTINGZ(Iref[0]) * pixel_size;
	// Place a black sphere in the centre of the box
	fh_bild << ".color 0 0 0 " << std::endl;
	fh_bild << ".sphere " << xcen << " " << ycen << " " << zcen << " 3 "  << std::endl;
	for (int ibody = 0; ibody < nr_bodies; ibody++)
	{
		// Sample evenly colors from the rainbow
		RFLOAT r, g, b;
		HSL2RGB((RFLOAT)ibody/(RFLOAT)nr_bodies, 1.0, 0.5, r, g, b);
		fh_bild << ".color " << r << " " << g << " " << b << std::endl;

		// Place a sphere at the centre-of-mass
		RFLOAT x = XX(com_bodies[ibody]) * pixel_size;
		RFLOAT y = YY(com_bodies[ibody]) * pixel_size;
		RFLOAT z = ZZ(com_bodies[ibody]) * pixel_size;
		// Add the center of the box to the coordinates
		x += pixel_size + xcen;
		y += pixel_size + ycen;
		z += pixel_size + zcen;
		fh_bild << ".sphere " << x << " " << y << " " << z << " 3 "  << std::endl;
		// Add a label
		fh_bild << ".cmov " << x+5 << " " << y+5 << " " << z+5 << std::endl;
		fh_bild << "body " << ibody+1 << std::endl;
		// Add an arrow for the direction of the rotation
		RFLOAT length = 10.;
		fh_bild << ".arrow " << x << " " << y << " " << z << " "
		        << x + length*XX(rotate_direction_bodies[ibody]) * pixel_size << " "
		        << y + length*YY(rotate_direction_bodies[ibody]) * pixel_size << " "
		        << z + length*ZZ(rotate_direction_bodies[ibody]) * pixel_size << " 1 " << std::endl;
	}

	// Close and write file to disc
	fh_bild.close();

}


void MlModel::setFourierTransformMaps(bool update_tau2_spectra, int nr_threads, RFLOAT strict_lowres_exp,
		   const MultidimArray<RFLOAT> *fourier_mask)
{

	bool do_heavy(true);

	int min_ires = -1;
	if (strict_lowres_exp > 0)
	{
		min_ires = ROUND(pixel_size * ori_size / strict_lowres_exp);
//		std::cout << "MlModel::setFourierTransformMaps: strict_lowres_exp = " << strict_lowres_exp
//		          << " pixel_size = " << pixel_size << " ori_size = " << ori_size << " min_ires = " << min_ires << std::endl;;
	}

	// Note that PPref.size() can be bigger than nr_bodies in multi-body refinement, due to extra PPrefs needed for overlapping bodies
	// These only exist in PPref form, they are not needed for reconstructions, only for subtractions in getFourierTransformsAndCtfs
	for (int iclass = 0; iclass < PPref.size(); iclass++)
	{

		MultidimArray<RFLOAT> Irefp;
		if (nr_bodies > 1)
		{
			// ibody deals with overlapping bodies here, as iclass can be larger than nr_bodies when bodies overlap,
			// but there are only nr_bodies Iref; ibody is the number of the original body (max nr_bodies)
			int ibody = pointer_body_overlap_inv[iclass];
			Irefp = Iref[ibody] * masks_bodies[iclass];
			// Place each body with its center-of-mass in the center of the box
			selfTranslate(Irefp, -com_bodies[ibody], DONT_WRAP);
		}
		else
		{
			Irefp = Iref[iclass];
		}

		if(PPrefRank.size() > 1)
			do_heavy = PPrefRank[iclass];

		if (update_tau2_spectra && iclass < nr_classes * nr_bodies)
		{
			PPref[iclass].computeFourierTransformMap(Irefp, tau2_class[iclass], current_size, nr_threads, true, do_heavy, min_ires, fourier_mask, do_gpu);
		}
		else
		{
			MultidimArray<RFLOAT> dummy;
			PPref[iclass].computeFourierTransformMap(Irefp, dummy, current_size, nr_threads, true, do_heavy, min_ires, fourier_mask, do_gpu);
		}
	}

}

void MlModel::initialiseDataVersusPrior(bool fix_tau)
{

	// Calculate straightforward (i.e. don't weight by number of particles) average sigma2_noise over all optics groups
	MultidimArray<RFLOAT> avg_sigma2_noise, sum_parts;
	avg_sigma2_noise.initZeros(ori_size /2 + 1);
	sum_parts.initZeros(ori_size /2 + 1);
	RFLOAT sum = 0.;
	for (int igroup = 0; igroup < sigma2_noise.size(); igroup++)
	{
		if (sigma2_noise[igroup].sum() > 0.)
		{
			avg_sigma2_noise += sigma2_noise[igroup];
			sum += 1.;
		}
	}
	avg_sigma2_noise /= sum;

	// Get the FT of all reference structures
	// The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
	// And spectrum is squared, so ori_size*ori_size in the 3D case!
	RFLOAT normfft = (ref_dim == 3 && data_dim == 2) ? (RFLOAT)(ori_size * ori_size) : 1.;

	int nr_classes_bodies = nr_classes * nr_bodies; // also set multiple bodies!
	for (int iclass = 0; iclass < nr_classes_bodies; iclass++)
	{
		// Initialise output arrays to correct size
		tau2_class[iclass].resize(ori_size /2 + 1);

		// Get the power spectrum of the reference
		MultidimArray<RFLOAT> spectrum(ori_size /2 + 1);
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
		data_vs_prior_class[iclass].resize(ori_size /2 + 1);
		if (nr_bodies > 1)
		{
			fsc_halves_class[iclass].initZeros(ori_size /2 + 1);
		}

		long long int nr_particles = 0;
		for (long int igroup = 0; igroup < nr_particles_per_group.size(); igroup++)
		{
			nr_particles += nr_particles_per_group[igroup];
		}
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
	for (int iclass = 0; iclass < nr_classes * nr_bodies; iclass++)
	{
		int maxres = 0;
		for (int ires = 0; ires < XSIZE(data_vs_prior_class[iclass]); ires++)
		{
			if (DIRECT_A1D_ELEM(data_vs_prior_class[iclass], ires) < 1.)
				break;
			maxres = ires;
		}
		int coverwindow = maxres*2 - 1;

		estimated_resolution[iclass] = 1./getResolution(maxres);
		total_fourier_coverage[iclass] = 0.;
		RFLOAT count = 0;
		for (long int k=FIRST_XMIPP_INDEX(coverwindow); k<=LAST_XMIPP_INDEX(coverwindow); k++) \
		   for (long int i=FIRST_XMIPP_INDEX(coverwindow); i<=LAST_XMIPP_INDEX(coverwindow); i++) \
			   for (long int j=FIRST_XMIPP_INDEX(coverwindow); j<=LAST_XMIPP_INDEX(coverwindow); j++) \
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


void MlModel::reset_class(int class_idx, int to_class_idx) {
	if (to_class_idx == -1) {
		Iref[class_idx] *= 0;
		Igrad1[class_idx].initZeros();
		if (pseudo_halfsets)
			Igrad1[class_idx+nr_classes].initZeros();
		Igrad2[class_idx].initConstant(Complex(1., 1.));
		pdf_class[class_idx] = 0;
		tau2_class[class_idx] *= 0.;
		data_vs_prior_class[class_idx] *= 0.;
		pdf_class[class_idx] *= 0;
		pdf_direction[class_idx] *= 0.;
	} else {
		Iref[class_idx] = Iref[to_class_idx];
		Igrad1[class_idx] = Igrad1[to_class_idx];
		if (pseudo_halfsets)
			Igrad1[class_idx+nr_classes] = Igrad1[to_class_idx+nr_classes];
		Igrad2[class_idx] = Igrad2[to_class_idx];
		pdf_class[class_idx] = pdf_class[to_class_idx];
		tau2_class[class_idx] = tau2_class[to_class_idx];
		data_vs_prior_class[class_idx] = data_vs_prior_class[to_class_idx];
		pdf_class[class_idx] = pdf_class[to_class_idx];
		pdf_direction[class_idx] = pdf_direction[to_class_idx];

	}
}


/////////// MlWsumModel
void MlWsumModel::initialise(MlModel &_model, FileName fn_sym, bool asymmetric_padding, bool _skip_gridding, bool _pseudo_halfsets)
{
	pixel_size = _model.pixel_size;
	nr_classes = _model.nr_classes;
	nr_bodies = _model.nr_bodies;
	nr_groups = _model.nr_groups;
	nr_optics_groups = _model.nr_optics_groups;
	nr_directions = _model.nr_directions;
	ref_dim = _model.ref_dim;
	data_dim = _model.data_dim;
	ori_size = _model.ori_size;
	pdf_class = _model.pdf_class;
	class_age = _model.class_age;
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
	// Don't need scale_correction and bfactor_correction, keep wsum_signal_product and wsum_reference_power instead
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

	wsum_signal_product.resize(nr_groups);
	wsum_reference_power.resize(nr_groups);
	for (long int igroup = 0; igroup < nr_groups; igroup++)
	{
		wsum_signal_product[igroup] = 0.;
		wsum_reference_power[igroup] = 0.;
	}

	// Resize MlWsumModel-specific vectors
	BackProjector BP(ori_size, ref_dim, fn_sym, interpolator, padding_factor, r_min_nn,
					 ML_BLOB_ORDER, ML_BLOB_RADIUS, ML_BLOB_ALPHA, data_dim, _skip_gridding);
	BPref.clear();

	pseudo_halfsets = _pseudo_halfsets;
	if (_pseudo_halfsets)
		BPref.resize(2 * nr_classes * nr_bodies, BP); // also set multiple bodies
	else
		BPref.resize(nr_classes * nr_bodies, BP); // also set multiple bodies

	sumw_group.resize(nr_optics_groups);

    // For ctf_premultiplied
    MultidimArray<RFLOAT > aux;
    aux.initZeros(ori_size / 2 + 1);
    sumw_ctf2.resize(nr_optics_groups, aux);

    // For subtomogram averaging
    sumw_stMulti.resize(nr_optics_groups, aux);

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
		if (pseudo_halfsets)
			BPref[iclass + nr_classes].initZeros(current_size);

		// Assume pdf_direction is already of the right size...
		pdf_direction[iclass].initZeros();
	}

	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		pdf_class[iclass] = 0.;
		class_age[iclass] = 0.;
		if (ref_dim == 2)
			prior_offset_class[iclass].initZeros();
	}

	// Initialise sigma2_noise spectra and sumw_group
	for (int igroup = 0; igroup < nr_optics_groups; igroup++)
	{
		sigma2_noise[igroup].initZeros();
        sumw_ctf2[igroup].initZeros();
        sumw_stMulti[igroup].initZeros();
        sumw_group[igroup] = 0.;
	}

	for (int igroup = 0; igroup < nr_groups; igroup++)
	{
		wsum_signal_product[igroup] = 0.;
		wsum_reference_power[igroup] = 0.;
	}

}

//#define DEBUG_PACK
#ifdef DEBUG_PACK
#define MAX_PACK_SIZE	  100000
#else
// Approximately 1024*1024*1024/8/2 ~ 0.5 Gb
#define MAX_PACK_SIZE 67101000
#endif

void MlWsumModel::pack(MultidimArray<RFLOAT> &packed)
{
	unsigned long long packed_size = 0;
	int spectral_size = (ori_size / 2) + 1;

	// for LL & avePmax & sigma2_offset & avg_norm_correction & sigma2_rot & sigma2_tilt & sigma2_psi
	packed_size += 7 ;

	// for all optics group-related stuff
	packed_size += nr_optics_groups * 3 * spectral_size;
	packed_size += nr_optics_groups;

	// for scale correction in groups
	packed_size += 2 * nr_groups;

	// for all class-related stuff
	// data is complex: multiply by two!
	packed_size += nr_classes * nr_bodies * 2 * (unsigned long long)BPref[0].getSize();
	packed_size += nr_classes * nr_bodies * (unsigned long long)BPref[0].getSize();
	packed_size += nr_classes * nr_bodies * (unsigned long long)nr_directions;
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

	for (int igroup = 0; igroup < nr_optics_groups; igroup++)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
		{
			DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n);
		}
        sigma2_noise[igroup].clear();
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sumw_ctf2[igroup])
        {
            DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(sumw_ctf2[igroup], n);
        }
        sumw_ctf2[igroup].clear();
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sumw_stMulti[igroup])
        {
            DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(sumw_stMulti[igroup], n);
        }
		sumw_stMulti[igroup].clear();
		DIRECT_MULTIDIM_ELEM(packed, idx++) = sumw_group[igroup];
	}

	for (int igroup = 0; igroup < nr_groups; igroup++)
	{
		DIRECT_MULTIDIM_ELEM(packed, idx++) = wsum_signal_product[igroup];
		DIRECT_MULTIDIM_ELEM(packed, idx++) = wsum_reference_power[igroup];
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
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdf_direction[iclass])
		{
			DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], n);
		}
		pdf_direction[iclass].clear();
	}
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
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

	unsigned long long idx = 0;
	int spectral_size = (ori_size / 2) + 1;

	LL = DIRECT_MULTIDIM_ELEM(packed, idx++);
	ave_Pmax = DIRECT_MULTIDIM_ELEM(packed, idx++);
	sigma2_offset = DIRECT_MULTIDIM_ELEM(packed, idx++);
	avg_norm_correction = DIRECT_MULTIDIM_ELEM(packed, idx++);
	sigma2_rot = DIRECT_MULTIDIM_ELEM(packed, idx++);
	sigma2_tilt = DIRECT_MULTIDIM_ELEM(packed, idx++);
	sigma2_psi = DIRECT_MULTIDIM_ELEM(packed, idx++);

	for (int igroup = 0; igroup < nr_optics_groups; igroup++)
	{
		sigma2_noise[igroup].resize(spectral_size);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
		{
			DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
		}
        sumw_ctf2[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sumw_ctf2[igroup])
        {
            DIRECT_MULTIDIM_ELEM(sumw_ctf2[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
        }
		sumw_stMulti[igroup].resize(spectral_size);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sumw_stMulti[igroup])
		{
			DIRECT_MULTIDIM_ELEM(sumw_stMulti[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
		}
		sumw_group[igroup] = DIRECT_MULTIDIM_ELEM(packed, idx++);
	}

	for (int igroup = 0; igroup < nr_groups; igroup++)
	{
		wsum_signal_product[igroup] = DIRECT_MULTIDIM_ELEM(packed, idx++);
		wsum_reference_power[igroup] = DIRECT_MULTIDIM_ELEM(packed, idx++);
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
		pdf_direction[iclass].resize(nr_directions);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdf_direction[iclass])
		{
			DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
		}
	}
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
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
	unsigned long long nr_groups = sigma2_noise.size();
	unsigned long long nr_classes = pdf_class.size();
	unsigned long long spectral_size = (ori_size / 2) + 1;
	unsigned long long packed_size = 0;
	unsigned long long idx_start, idx_stop;

	// for LL & avePmax & sigma2_offset & avg_norm_correction & sigma2_rot & sigma2_tilt & sigma2_psi
	packed_size += 7 ;
	// for group-related spectra
	packed_size += nr_optics_groups * 3 * spectral_size; // sigma2_noise[spectral_size] and sumw_stMulti[spectral_size] and sumw_ctf2[spectral_size]
	packed_size += nr_optics_groups; // sumw_group
	// for scale correction in groups
	packed_size += 2 * nr_groups; // wsum_signal_product, wsum_reference_power
	// for all class-related stuff
	// data is complex: multiply by two!
	packed_size += BPref.size() * 2 * (unsigned long long) BPref[0].getSize(); // BPref.data
	packed_size += BPref.size() * (unsigned long long) BPref[0].getSize(); // BPref.weight
	packed_size += pdf_direction.size() * (unsigned long long) nr_directions; // pdf_directions
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
//#define DEBUG_PACK
#ifdef DEBUG_PACK
	std::cerr << " PACK: idx_start= " << idx_start << " idx_stop= " << idx_stop << " piece= " << piece << " nr_pieces= " << nr_pieces <<" packed_size= "<<packed_size<< std::endl;
	std::cerr << " nr_classes= " << nr_classes << " nr_groups= " << nr_groups << " packed_size= " << packed_size << std::endl;
	std::cerr << " MULTIDIM_SIZE(sigma2_noise[0])= " << MULTIDIM_SIZE(sigma2_noise[0]) /*<< " MULTIDIM_SIZE(wsum_signal_product_spectra[0])= " << MULTIDIM_SIZE(wsum_signal_product[0]) << " MULTIDIM_SIZE(wsum_reference_power_spectra[0])= " << MULTIDIM_SIZE(wsum_reference_power[0]) */<< std::endl;
	std::cerr << " sigma2_noise.size()= " << sigma2_noise.size() << " wsum_signal_product_spectra.size()= " << wsum_signal_product.size() << " wsum_signal_product_spectra.size()= " << wsum_signal_product.size() << std::endl;
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

	for (int igroup = 0; igroup < nr_optics_groups; igroup++)
	{

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sigma2_noise[igroup])
		{
			if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(sigma2_noise[igroup], n);
			ori_idx++;
		}
		if (idx == ori_idx && do_clear)
			sigma2_noise[igroup].clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sumw_ctf2[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(sumw_ctf2[igroup], n);
            ori_idx++;
        }
        if (idx == ori_idx && do_clear)
            sumw_ctf2[igroup].clear();

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sumw_stMulti[igroup])
		{
			if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) =DIRECT_MULTIDIM_ELEM(sumw_stMulti[igroup], n);
			ori_idx++;
		}
		if (idx == ori_idx && do_clear)
			sumw_stMulti[igroup].clear();

		if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = sumw_group[igroup];
		ori_idx++;

	}

	for (int igroup = 0; igroup < nr_groups; igroup++)
	{

		if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = wsum_signal_product[igroup];
		ori_idx++;
		if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = wsum_reference_power[igroup];
		ori_idx++;

	}

	for (int iclass = 0; iclass < BPref.size(); iclass++)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].data) {
			if (ori_idx >= idx_start && ori_idx < idx_stop)
				DIRECT_MULTIDIM_ELEM(packed, idx++) = (DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).real;
			ori_idx++;
			if (ori_idx >= idx_start && ori_idx < idx_stop)
				DIRECT_MULTIDIM_ELEM(packed, idx++) = (DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).imag;
			ori_idx++;
		}
		// Only clear after the whole array has been packed... i.e. not when we reached the pack_size halfway through
		if (idx == ori_idx && do_clear)
			BPref[iclass].data.clear();

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].weight) {
			if (ori_idx >= idx_start && ori_idx < idx_stop)
				DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(BPref[iclass].weight, n);
			ori_idx++;
		}
		if (idx == ori_idx && do_clear)
			BPref[iclass].weight.clear();
	}

	for (int iclass = 0; iclass < pdf_direction.size(); iclass++)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdf_direction[iclass])
		{
			if (ori_idx >= idx_start && ori_idx < idx_stop) DIRECT_MULTIDIM_ELEM(packed, idx++) = DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], n);
			ori_idx++;
		}
		if (idx == ori_idx && do_clear)
			pdf_direction[iclass].clear();
	}

	for (int iclass = 0; iclass < nr_classes; iclass++)
	{

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
	//std::cerr << " PACK idx= " << idx << " idx_stop-idx_start= " << idx_stop-idx_start << " idx_start= " << idx_start << " idx_stop= " << idx_stop	<< std::endl;
	if (idx != idx_stop-idx_start)
	{
		std::cerr << "idx= " << idx << "ori_idx= " << ori_idx << " idx_start= " << idx_start << " idx_stop= " << idx_stop << " packed_size= " << packed_size << std::endl;
		REPORT_ERROR("MlWsumModel::pack: idx != idx_stop-idx_start");

	}

}

void MlWsumModel::unpack(MultidimArray<RFLOAT> &packed, int piece, bool do_clear)
{


	int nr_groups = sigma2_noise.size();
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

	for (int igroup = 0; igroup < nr_optics_groups; igroup++)
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
            sumw_ctf2[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sumw_ctf2[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
                DIRECT_MULTIDIM_ELEM(sumw_ctf2[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
            ori_idx++;
        }

        if (idx == ori_idx)
            sumw_stMulti[igroup].resize(spectral_size);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sumw_stMulti[igroup])
        {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
                DIRECT_MULTIDIM_ELEM(sumw_stMulti[igroup], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
            ori_idx++;
        }

		if (ori_idx >= idx_start && ori_idx < idx_stop)
			sumw_group[igroup] = DIRECT_MULTIDIM_ELEM(packed, idx++);
		ori_idx++;
	}

	for (int igroup = 0; igroup < nr_groups; igroup++)
	{

		if (ori_idx >= idx_start && ori_idx < idx_stop)
			wsum_signal_product[igroup] = DIRECT_MULTIDIM_ELEM(packed, idx++);
		ori_idx++;
		if (ori_idx >= idx_start && ori_idx < idx_stop)
			wsum_reference_power[igroup] = DIRECT_MULTIDIM_ELEM(packed, idx++);
		ori_idx++;

	}

	for (int iclass = 0; iclass < BPref.size(); iclass++) {
		if (idx == ori_idx)
			BPref[iclass].initialiseDataAndWeight(current_size);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].data) {
			if (ori_idx >= idx_start && ori_idx < idx_stop)
				(DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).real = DIRECT_MULTIDIM_ELEM(packed, idx++);
			ori_idx++;

			if (ori_idx >= idx_start && ori_idx < idx_stop)
				(DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n)).imag = DIRECT_MULTIDIM_ELEM(packed, idx++);
			ori_idx++;
			//DIRECT_MULTIDIM_ELEM(BPref[iclass].data, n) = Complex(re, im);
		}

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(BPref[iclass].weight) {
			if (ori_idx >= idx_start && ori_idx < idx_stop)
				DIRECT_MULTIDIM_ELEM(BPref[iclass].weight, n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
			ori_idx++;
		}
	}

	for (int iclass = 0; iclass < pdf_direction.size(); iclass++)
	{
		if (idx == ori_idx)
			pdf_direction[iclass].resize(nr_directions);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pdf_direction[iclass])
		{
			if (ori_idx >= idx_start && ori_idx < idx_stop)
				DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], n) = DIRECT_MULTIDIM_ELEM(packed, idx++);
			ori_idx++;
		}

	}

	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
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
	//std::cerr << " UNPACK piece= " << piece << " idx= " << idx << " idx_stop-idx_start= " << idx_stop-idx_start << " idx_start= " << idx_start << " idx_stop= " << idx_stop	 << std::endl;
	if (idx != idx_stop-idx_start)
	{
		std::cerr << "idx= " << idx << "ori_idx= " << ori_idx << " idx_start= " << idx_start << " idx_stop= " << idx_stop << " packed_size= " << packed_size << std::endl;
		REPORT_ERROR("MlWsumModel::unpack: idx != idx_stop-idx_start");
	}


}



