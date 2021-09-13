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
#include "src/healpix_sampling.h"
//#define DEBUG_SAMPLING
//#define DEBUG_CHECKSIZES
//#define DEBUG_HELICAL_ORIENTATIONAL_SEARCH

void HealpixSampling::clear()
{
	is_3D = false;
	isRelax = false;
	fn_sym = "C1";
	fn_sym_relax = "C1";
	limit_tilt = psi_step = offset_range = offset_step = helical_offset_step = psi_step_ori = offset_range_ori = offset_step_ori = 0.;
	random_perturbation = perturbation_factor = 0.;
	// Jun19,2015 - Shaoda, Helical refinement
	helical_offset_step = -1.;
	directions_ipix.clear();
	rot_angles.clear();
	tilt_angles.clear();
	psi_angles.clear();
	translations_x.clear();
	translations_y.clear();
	translations_z.clear();
	L_repository.clear();
	R_repository.clear();
	L_repository_relax.clear();
	R_repository_relax.clear();
	pgGroup = pgOrder = 0;
	pgGroupRelaxSym = pgOrderRelaxSym = 0;

}

void HealpixSampling::initialise(
		int ref_dim,
		bool do_3d_trans,
		bool do_changepsi,
		bool do_warnpsi,
		bool do_local_searches_helical,
		bool do_helical_refine,
		RFLOAT rise_Angst,
		RFLOAT twist_deg)
{

	if (ref_dim != -1)
		is_3D = (ref_dim == 3);

	// Set the symmetry relaxation flag
	isRelax = fn_sym_relax == "" ? false : true;

	// Set flag for x,y,z-translations
	is_3d_trans = do_3d_trans;

	// By default psi_step is approximate sampling of rot,tilt in 3D; and 10 degrees in 2D
	if (psi_step < 0)
	{
		if (is_3D)
			psi_step = 360. / (6 * ROUND(std::pow(2., healpix_order)));
		else
			psi_step = 10.;
	}

	if (perturbation_factor < 0. || perturbation_factor > 1.)
		REPORT_ERROR("HealpixSampling::initialise: random perturbation factor should be between 0 and 1.");

	if (is_3D)
	{
		healpix_base.Set(healpix_order, NEST);

		// Set up symmetry
		R_repository.clear();
		L_repository.clear();
		initialiseSymMats(fn_sym, pgGroup, pgOrder, R_repository, L_repository);

		// Set up symmetry matrices for symmetry relax
		if (fn_sym_relax != "")
		{
			if (fn_sym_relax[0] != 'C' && fn_sym_relax[0] != 'c')
				REPORT_ERROR("Sorry, symmetry relaxation is currently available only for cyclic (Cn) point groups. For other symmetries, please see https://github.com/3dem/relion/issues/796.");
			R_repository_relax.clear();
			L_repository_relax.clear();
			initialiseSymMats(fn_sym_relax, pgGroupRelaxSym, pgOrderRelaxSym, R_repository_relax, L_repository_relax);
		}
	}
	else
	{
		int t_nr_psi = CEIL(360./psi_step);
		if(t_nr_psi%32!=0 && do_changepsi)
		{

			// Force-adjust psi_step to be multiples of 32 (for efficient GPU calculations)
			t_nr_psi = CEIL((float)t_nr_psi / 32.0)*32;

			if (do_warnpsi)
				std::cout << " + WARNING: Changing psi sampling rate (before oversampling) to " <<  360./(RFLOAT)t_nr_psi << " degrees, for more efficient GPU calculations" << std::endl;

		}
		psi_step = 360./(RFLOAT)t_nr_psi;
		fn_sym = "C1"; // This may not be set yet if restarting a 2D run....
	}

	// Store the not-oversampled translations, and make sure oversampled sampling is 1 pixel
	//setTranslations();
	// May06,2015 - Shaoda & Sjors, Helical translational searches
	setTranslations(-1, -1, do_local_searches_helical, do_helical_refine, -1, rise_Angst, twist_deg);

	// Store the non-oversampled projection directions
	setOrientations(-1, -1.);

	// Random perturbation and filling of the directions, psi_angles and translations vectors
	resetRandomlyPerturbedSampling();

	// SHWS 27feb2020: Set original sampling rates to allow 2D/3D classifications using coarser ones in earlier iterations
	healpix_order_ori = healpix_order;
	psi_step_ori = psi_step;
	offset_range_ori = offset_range;
	offset_step_ori = offset_step;


}

void HealpixSampling::initialiseSymMats(FileName fn_sym_, int & pgGroup_,
    		int & pgOrder_, std::vector <Matrix2D<RFLOAT> > & R_repository_,
			std::vector <Matrix2D<RFLOAT> > & L_repository_)
{
	// Set up symmetry
	SymList SL;
	SL.isSymmetryGroup(fn_sym_, pgGroup_, pgOrder_);
	SL.read_sym_file(fn_sym_);

	// Precalculate (3x3) symmetry matrices
	Matrix2D<RFLOAT>  L(4, 4), R(4, 4);
	Matrix2D<RFLOAT>  Identity(3,3);
	Identity.initIdentity();
	R_repository_.clear();
	L_repository_.clear();
	R_repository_.push_back(Identity);
	L_repository_.push_back(Identity);
	for (int isym = 0; isym < SL.SymsNo(); isym++)
	{
		SL.get_matrices(isym, L, R);
		R.resize(3, 3);
		L.resize(3, 3);
		R_repository_.push_back(R);
		L_repository_.push_back(L);
	}

}

void HealpixSampling::resetRandomlyPerturbedSampling()
{

	// Actual instance of random perturbation
	// Add to the random perturbation from the last iteration, so it keeps changing strongly...
	random_perturbation += rnd_unif(0.5*perturbation_factor, perturbation_factor);
	random_perturbation = realWRAP(random_perturbation, -perturbation_factor, perturbation_factor);

}

void HealpixSampling::read(FileName fn_in)
{

    // Open input file
    std::ifstream in(fn_in.data(), std::ios_base::in);
    if (in.fail())
        REPORT_ERROR( (std::string) "HealpixSampling::readStar: File " + fn_in + " cannot be read." );

    MetaDataTable MD;

    // Read general stuff
    MD.readStar(in, "sampling_general");
    in.close();

    if (!MD.getValue(EMDL_SAMPLING_IS_3D, is_3D) ||
    	!MD.getValue(EMDL_SAMPLING_IS_3D_TRANS, is_3d_trans) ||
		!MD.getValue(EMDL_SAMPLING_PSI_STEP, psi_step) ||
		!MD.getValue(EMDL_SAMPLING_OFFSET_RANGE, offset_range) ||
		!MD.getValue(EMDL_SAMPLING_OFFSET_STEP, offset_step) ||
		!MD.getValue(EMDL_SAMPLING_PERTURBATION_FACTOR, perturbation_factor))
		REPORT_ERROR("HealpixSampling::readStar: incorrect sampling_general table");
    // Jun19,2015 - Shaoda, Helical translational searches, backward compatibility
    if (!MD.getValue(EMDL_SAMPLING_HELICAL_OFFSET_STEP, helical_offset_step))
    	helical_offset_step = -1.;

    // SHWS 27Feb2020: backwards compatibility: older star files will not yet have original sampling parameters, just use current ones
    if (!MD.getValue(EMDL_SAMPLING_OFFSET_STEP_ORI, offset_step_ori)) offset_step_ori = offset_step;
    if (!MD.getValue(EMDL_SAMPLING_OFFSET_RANGE_ORI, offset_range_ori)) offset_range_ori = offset_range;
    if (!MD.getValue(EMDL_SAMPLING_PSI_STEP_ORI, psi_step_ori)) psi_step_ori = psi_step;

	if (is_3D)
	{
		if (!MD.getValue(EMDL_SAMPLING_HEALPIX_ORDER, healpix_order) ||
			!MD.getValue(EMDL_SAMPLING_SYMMETRY, fn_sym) ||
			!MD.getValue(EMDL_SAMPLING_LIMIT_TILT, limit_tilt) )
			REPORT_ERROR("HealpixSampling::readStar: incorrect sampling_general table for 3D sampling");

		// For 3D samplings reset psi_step to -1:
		// By default it will then be set to the healpix sampling
		// Only if the --psi_step option is given on the command line it will be set to something different!
		psi_step = -1.;

	    // SHWS 27Feb2020: backwards compatibility: older star files will not yet have original sampling parameters, just use current ones
	    if (!MD.getValue(EMDL_SAMPLING_HEALPIX_ORDER_ORI, healpix_order_ori)) healpix_order_ori = healpix_order;

	}
	else
	{
		fn_sym = "irrelevant";
		limit_tilt = 0.;
		healpix_order = 0;
	}

}

void HealpixSampling::write(FileName fn_out)
{
	MetaDataTable MD;
	std::ofstream  fh;
	FileName fn_tmp;

    fn_tmp = fn_out + "_sampling.star";
    fh.open((fn_tmp).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR( (std::string)"HealpixSampling::write: Cannot write file: " + fn_tmp);

    MD.setIsList(true);
    MD.addObject();
	MD.setName("sampling_general");
	MD.setValue(EMDL_SAMPLING_IS_3D, is_3D);
	MD.setValue(EMDL_SAMPLING_IS_3D_TRANS, is_3d_trans);
	if (is_3D)
	{
		MD.setValue(EMDL_SAMPLING_HEALPIX_ORDER, healpix_order);
		MD.setValue(EMDL_SAMPLING_SYMMETRY, fn_sym);
		MD.setValue(EMDL_SAMPLING_LIMIT_TILT, limit_tilt);
	}
	MD.setValue(EMDL_SAMPLING_PSI_STEP, psi_step);
	MD.setValue(EMDL_SAMPLING_OFFSET_RANGE, offset_range);
	MD.setValue(EMDL_SAMPLING_OFFSET_STEP, offset_step);
	// Jun19,2015 - Shaoda, Helical translational searches
	MD.setValue(EMDL_SAMPLING_HELICAL_OFFSET_STEP, helical_offset_step);
	MD.setValue(EMDL_SAMPLING_PERTURB, random_perturbation);
	MD.setValue(EMDL_SAMPLING_PERTURBATION_FACTOR, perturbation_factor);

	//27Feb2020 SHWS: write original sampling rates to allow 2D/3D classifications to use coarser ones in initial iterations
    MD.setValue(EMDL_SAMPLING_HEALPIX_ORDER_ORI, healpix_order_ori);
    MD.setValue(EMDL_SAMPLING_PSI_STEP_ORI, psi_step_ori);
    MD.setValue(EMDL_SAMPLING_OFFSET_RANGE_ORI, offset_range_ori);
    MD.setValue(EMDL_SAMPLING_OFFSET_STEP_ORI, offset_step_ori);

	MD.write(fh);

	// In the 3D case, also write a table with the sampled rot, tilt angles
	if (is_3D)
	{
		MD.clear();
		MD.setIsList(false);
		MD.setName("sampling_directions");
		for (long int idir = 0; idir < NrDirections(); idir++)
		{
			RFLOAT rot, tilt;
			getDirection(idir, rot, tilt);
			MD.addObject();
			MD.setValue(EMDL_ORIENT_ROT, rot);
			MD.setValue(EMDL_ORIENT_TILT, tilt);
		}
		MD.write(fh);
	}

	// Close the file
	fh.close();
}

void HealpixSampling::setTranslations(
		RFLOAT new_offset_step,
		RFLOAT new_offset_range,
		bool do_local_searches_helical,
		bool do_helical_refine,
		RFLOAT new_helical_offset_step,
		RFLOAT helical_rise_Angst,
		RFLOAT helical_twist_deg)
{
	// Ordinary single particles
	int maxp; // Max half nr samplings in all directions
	RFLOAT xoff, yoff, zoff, max2, old_offset_step, old_offset_range; // Offset lengths

	// Helical refinement
	int maxh; // Max half nr samplings in along helical axis
	RFLOAT h_range, old_helical_offset_step; // Translations along helical axis

	// Check old and new offsets
	old_offset_step = offset_step;  // can be < 0 ??????
	old_offset_range = offset_range;  // can be < 0 ??????
	old_helical_offset_step = helical_offset_step;  // can be < 0
	if ( (new_offset_step > 0.) && (new_offset_range >= 0.) )
	{
		offset_step = new_offset_step;
		offset_range = new_offset_range;
	}
	else
	{
		if (!(offset_step > 0.))
		{
			std::cerr << " offset_range= " << offset_range << " offset_step= " << offset_step << std::endl;
			REPORT_ERROR("HealpixSampling::setTranslations BUG %% Trying to set translations with uninitialised offset_step!");
		}
	}
	// Sometimes new offsets are set to -1, that means the old offsets remain unchanged.
	new_offset_step = offset_step;  // > 0
	new_offset_range = offset_range;  // >= 0

	// Ordinary single particles
	maxp = CEIL(offset_range / offset_step); // Perpendicular to helical axis (P1, P2)
	maxh = maxp;
	// Helical refinement
	if (do_helical_refine)
	{
		// Assume all helical parameters are valid (this should be checked before in ml_optimiser.cpp)
		helical_rise_Angst = fabs(helical_rise_Angst);
		helical_twist_deg = fabs(helical_twist_deg);

		// Search range (half) along helical axis = (-0.5 * rise, +0.5 * rise)
		h_range = (helical_rise_Angst / 2.);

		// If continue from old run or new offset is not applicable...
		if (new_helical_offset_step < 0.)
			new_helical_offset_step = old_helical_offset_step;

		// New_helical_offset_step is not OK.
		// (1) negative value
		// (2) larger than before (if there is a valid offset before)
		// (3) larger than new_offset_step
		// (4) samplings along helical axis is less than 3
		if ( (new_helical_offset_step < 0.)
				|| ( (new_helical_offset_step > old_helical_offset_step) && (old_helical_offset_step > 0.) )
				|| (new_helical_offset_step > new_offset_step)
				|| ((helical_rise_Angst / new_helical_offset_step) < 3.) )
		{
			// First try 'new_offset_step'
			new_helical_offset_step = new_offset_step;
			// Change to 'old_helical_offset_step' if the old is smaller
			if ( (old_helical_offset_step > 0.) && (new_helical_offset_step > old_helical_offset_step) )
				new_helical_offset_step = old_helical_offset_step;
			// New_helical_offset_step should be finer than 1/3 the helical rise
			if ( (helical_rise_Angst / new_helical_offset_step) < 3.)
				new_helical_offset_step = helical_rise_Angst / 3.;
		}

		maxh = CEIL(h_range / new_helical_offset_step); // Out of range samplings will be excluded next
		if (do_local_searches_helical) // Local searches along helical axis
		{
			// Local searches (2*2+1=5 samplings)
			if (maxh > 2)
				maxh = 2;
			// New helical offset step is smaller than 1/3 of the old one, samplings should be increased.
			if ( (old_helical_offset_step > 0.) && ((old_helical_offset_step / new_helical_offset_step) > 3) )
				maxh = FLOOR(old_helical_offset_step / new_helical_offset_step); // Use FLOOR here!
			// Local searches should not be wider than 1/3 of the helical rise
			if ( ((new_helical_offset_step * maxh) > (helical_rise_Angst / 6.)) )
			{
				maxh = FLOOR(helical_rise_Angst / 6. / new_helical_offset_step); // Use FLOOR here!
				if (maxh < 1) // At least we should do some searches...
					maxh = 1;
			}
		}
		// DEBUG - this should not happen
		if (maxh < 0)
			maxh = 0;
		helical_offset_step = new_helical_offset_step;
	}

	// DEBUG
	if ( (maxh < 0) || (maxp < 0) )
	{
		std::cerr << "maxh= " << maxh << " maxp= " << maxp << std::endl;
		REPORT_ERROR("HealpixSampling::setTranslations BUG %% No translations to set! ('maxh' or 'maxp' < 0)");
	}

	translations_x.clear();
	translations_y.clear();
	translations_z.clear();
	for (long int ix = -maxh; ix <= maxh; ix++)
	{
		// For helices use a different step size along helical axis X
		xoff = (do_helical_refine) ? (ix * helical_offset_step) : (ix * offset_step);
		// For helical refinement, exclude xoff outside the range of (-0.5 * rise, +0.5 * rise)
		if ( (do_helical_refine) && (ix != 0) && (fabs(xoff) > fabs(helical_rise_Angst / 2.)) )
			continue;

		for (long int iy = -maxp; iy <= maxp; iy++)
		{
			yoff = iy * offset_step;
			// For helices do not limit translations along helical axis X
			max2 = (do_helical_refine) ? (yoff * yoff) : (xoff * xoff + yoff * yoff);
			if (is_3d_trans)
			{
				for (long int iz = -maxp; iz <= maxp; iz++)
				{
					zoff = iz * offset_step;
					if ((max2 + zoff * zoff) <= (offset_range * offset_range))
					{
						translations_y.push_back(yoff);
						if (do_helical_refine) // Z axis corresponds to the helical axis in 3D subtomogram averaging !!!
						{
							translations_x.push_back(zoff);
							translations_z.push_back(xoff);
						}
						else
						{
							translations_x.push_back(xoff);
							translations_z.push_back(zoff);
						}

					}
				}
			}
			else
			{
				if (max2 < (offset_range * offset_range) + 0.001)  // +0.001 prevent precision errors in relion-3.1
				{
					translations_x.push_back(xoff);
					translations_y.push_back(yoff);
				}
			}
		}
	}
#ifdef DEBUG_SETTRANS
	std::cerr << " is_3d_trans= " << is_3d_trans << std::endl;
	for (int i = 0; i < translations_x.size(); i++)
		std::cerr << " translations_x[i]= " << translations_x[i] << std::endl;
#endif
	return;
}

/* Set only a single translation */
void HealpixSampling::addOneTranslation(
		RFLOAT offset_x,
		RFLOAT offset_y,
		RFLOAT offset_z,
		bool do_clear,
		bool do_helical_refine,
		RFLOAT rot_deg,
		RFLOAT tilt_deg,
		RFLOAT psi_deg)
{
	if (do_clear)
	{
		translations_x.clear();
		translations_y.clear();
		translations_z.clear();
	}
	if (do_helical_refine)
		transformCartesianAndHelicalCoords(offset_x, offset_y, offset_z, offset_x, offset_y, offset_z, rot_deg, tilt_deg, psi_deg, ((is_3d_trans) ? (3) : (2)), CART_TO_HELICAL_COORDS);
	translations_x.push_back(offset_x);
	translations_y.push_back(offset_y);
	if (is_3d_trans)
		translations_z.push_back(offset_z);
}

void HealpixSampling::setOrientations(int _order, RFLOAT _psi_step)
{

	// Initialise
	directions_ipix.clear();
	rot_angles.clear();
	tilt_angles.clear();
	psi_angles.clear();

	if (_order >= 0)
		healpix_order = _order;

	// Setup the HealPix object
	// For adaptive oversampling only precalculate the COARSE sampling!
	if (_order >= 0)
		healpix_base.Set(_order, NEST);

	// 3D directions
	if (is_3D)
	{
		RFLOAT rot, tilt;
		for (long int ipix = 0; ipix < healpix_base.Npix(); ipix++)
		{
			getDirectionFromHealPix(ipix, rot, tilt);

			// Push back as Matrix1D's in the vectors
			rot_angles.push_back(rot);
			tilt_angles.push_back(tilt);
			directions_ipix.push_back(ipix);


		}
//#define DEBUG_SAMPLING
#ifdef  DEBUG_SAMPLING
		writeAllOrientationsToBild("orients_all.bild", "1 0 0 ", 0.020);
#endif
		// Now remove symmetry-related pixels if not relaxing symmetry
		// TODO check size of healpix_base.max_pixrad
		if (!isRelax)
			removeSymmetryEquivalentPoints(0.5 * RAD2DEG(healpix_base.max_pixrad()));

#ifdef  DEBUG_SAMPLING
		writeAllOrientationsToBild("orients_sym.bild", "0 1 0 ", 0.021);
#endif

		// Also remove limited tilt angles
		removePointsOutsideLimitedTiltAngles();

#ifdef  DEBUG_SAMPLING
		if (ABS(limit_tilt) < 90.)
			writeAllOrientationsToBild("orients_tilt.bild", "1 1 0 ", 0.022);
#endif

	}
	else
	{
		rot_angles.push_back(0.);
		tilt_angles.push_back(0.);
		directions_ipix.push_back(-1);
	}

	// 2D in-plane angles
	// By default in 3D case: use more-or-less same psi-sampling as the 3D healpix object
	// By default in 2D case: use 5 degree
	if (_psi_step > 0.)
		psi_step = _psi_step;

	int nr_psi = CEIL(360./psi_step);
	RFLOAT psi;
	psi_step = 360./(RFLOAT)nr_psi;
	for (int ipsi = 0; ipsi < nr_psi; ipsi++)
	{
		psi = ipsi * psi_step;
		psi_angles.push_back(psi);
	}

//#define DEBUG_SAMPLING
#ifdef  DEBUG_SAMPLING
	writeAllOrientationsToBild("orients_final.bild", "1 0 0 ", 0.020);
#endif

}

/* Set only a single orientation */
void HealpixSampling::addOneOrientation(RFLOAT rot, RFLOAT tilt, RFLOAT psi, bool do_clear)
{
	if (do_clear)
	{
		directions_ipix.clear();
		rot_angles.clear();
		tilt_angles.clear();
		psi_angles.clear();
	}

	// 3D directions
	if (is_3D)
	{
		rot_angles.push_back(rot);
		tilt_angles.push_back(tilt);
		directions_ipix.push_back(-1);
	}
	else
	{
		rot_angles.push_back(0.);
		tilt_angles.push_back(0.);
		directions_ipix.push_back(-1);
	}

	// in-plane rotation
	psi_angles.push_back(psi);


}


void HealpixSampling::writeAllOrientationsToBild(FileName fn_bild, std::string rgb, RFLOAT size)
{
    std::ofstream out;
    out.open (fn_bild.c_str());
    if (!out)
        REPORT_ERROR( (std::string)"HealpixSampling::writeAllOrientationsToBild: Cannot write file: " + fn_bild);


    out << ".color 1 0 0 \n";
    out << ".arrow 0 0 0 1 0 0 0.01 \n";
    out << ".color 0 1 0 \n";
    out << ".arrow 0 0 0  0 1 0 0.01 \n";
    out << ".color 0 0 1 \n";
    out << ".arrow 0 0 0 0 0 1 0.01 \n";


    Matrix1D<RFLOAT> v(3);
	out << ".color " << rgb << std::endl;

	for (unsigned long int ipix = 0; ipix < rot_angles.size(); ipix++)
	{
		Euler_angles2direction(rot_angles[ipix], tilt_angles[ipix], v);
		out <<  ".sphere " << XX(v) << " " << YY(v) << " " << ZZ(v)  << " " <<  floatToString(size) << std::endl;
	}

	out.close();

}

void HealpixSampling::writeNonZeroPriorOrientationsToBild(FileName fn_bild, RFLOAT rot_prior, RFLOAT tilt_prior,
		std::vector<int> &pointer_dir_nonzeroprior, std::string rgb, RFLOAT size)
{
    std::ofstream out;
    out.open (fn_bild.c_str());
    if (!out)
        REPORT_ERROR( (std::string)"HealpixSampling::writeNonZeroOrientationsToBild: Cannot write file: " + fn_bild);


    out << ".color 1 0 0 \n";
    out << ".arrow 0 0 0 1 0 0 0.01 \n";
    out << ".color 0 1 0 \n";
    out << ".arrow 0 0 0  0 1 0 0.01 \n";
    out << ".color 0 0 1 \n";
    out << ".arrow 0 0 0 0 0 1 0.01 \n";

    Matrix1D<RFLOAT> v(3);

	Euler_angles2direction(rot_prior, tilt_prior, v);
	out << ".color 1 0 0 \n";
	out <<  ".sphere " << XX(v) << " " << YY(v) << " " << ZZ(v) << " " <<  floatToString(size) << std::endl;

	out << ".color " << rgb << std::endl;
	for (unsigned long int ipix = 0; ipix < pointer_dir_nonzeroprior.size(); ipix++)
	{
		long int idir = pointer_dir_nonzeroprior[ipix];
		Euler_angles2direction(rot_angles[idir], tilt_angles[idir], v);
		out <<  ".sphere " << XX(v) << " " << YY(v) << " " << ZZ(v) << " " << floatToString(size) << std::endl;
	}

	out.close();

}

RFLOAT HealpixSampling::calculateDeltaRot(Matrix1D<RFLOAT> my_direction, RFLOAT rot_prior)
{
	// Rotate the x,y-components of the direction, according to rot-prior
	Matrix1D< RFLOAT > my_rot_direction;
	Matrix2D< RFLOAT > A;
	rotation2DMatrix(rot_prior, A);
	my_rot_direction = A.inv() * my_direction;
	// Get component along the new Y-axis
	return fabs(ASIND(my_rot_direction(1)));
}

void HealpixSampling::selectOrientationsWithNonZeroPriorProbability(
		RFLOAT prior_rot, RFLOAT prior_tilt, RFLOAT prior_psi,
		RFLOAT sigma_rot, RFLOAT sigma_tilt, RFLOAT sigma_psi,
    	std::vector<int> &pointer_dir_nonzeroprior, std::vector<RFLOAT> &directions_prior,
    	std::vector<int> &pointer_psi_nonzeroprior, std::vector<RFLOAT> &psi_prior,
		bool do_bimodal_search_psi,
		RFLOAT sigma_cutoff, RFLOAT sigma_tilt_from_ninety, RFLOAT sigma_psi_from_zero)
{
	pointer_dir_nonzeroprior.clear();
	directions_prior.clear();
	// Do not check the mates again
	std::vector<bool> idir_flag(rot_angles.size(), false);

	if (is_3D)
	{
		//std::cerr<<"sigma_rot "<<sigma_rot<<" sigma_tilt "<<sigma_tilt<<std::endl;
		Matrix1D<RFLOAT> prior90_direction;
		if (sigma_tilt_from_ninety > 0.)
		{
			// pre-calculate original (0,90) direction
			Euler_angles2direction(0., 90., prior90_direction);
		}

		// Loop over all directions
		RFLOAT sumprior = 0.;
		RFLOAT sumprior_withsigmafromzero = 0.;
		// Keep track of the closest distance to prevent 0 orientations
		RFLOAT best_ang = 9999.;
		long int best_idir = -999;

		for (long int idir = 0; idir < rot_angles.size(); idir++)
		{
			// Check if this direction was met before as symmetry mate
			if (idir_flag[idir] == true)
					continue;

			bool is_nonzero_pdf = false;

			// Any prior involving BOTH rot and tilt.
			if ( (sigma_rot > 0.) && (sigma_tilt > 0.) )
			{
				// Get the direction of the prior
				Matrix1D<RFLOAT> prior_direction, my_direction, sym_direction, best_direction;
				Euler_angles2direction(prior_rot, prior_tilt, prior_direction);

				// Get the current direction in the loop
				Euler_angles2direction(rot_angles[idir], tilt_angles[idir], my_direction);
				best_direction = my_direction;

				// Loop over all symmetry operators to find the operator that brings this direction nearest to the prior if no symmetry relaxation
				if (!isRelax)
				{
					RFLOAT best_dotProduct = dotProduct(prior_direction, my_direction);
					for (int j = 0; j < R_repository.size(); j++)
					{
						sym_direction =  L_repository[j] * (my_direction.transpose() * R_repository[j]).transpose();
						RFLOAT my_dotProduct = dotProduct(prior_direction, sym_direction);
						if (my_dotProduct > best_dotProduct)
						{
							best_direction = sym_direction;
							best_dotProduct = my_dotProduct;
						}
					}
				}

				// Now that we have the best direction, find the corresponding prior probability
				RFLOAT diffang = ACOSD( dotProduct(best_direction, prior_direction) );
				if (diffang > 180.)
					diffang = ABS(diffang - 360.);
				if (do_bimodal_search_psi && (diffang > 90.))  // KThurber
					diffang = ABS(diffang - 180.);	// KThurber

				// Only consider differences within sigma_cutoff * sigma_rot
				// TODO: If sigma_rot and sigma_tilt are not the same (NOT for helices)?
				RFLOAT biggest_sigma = XMIPP_MAX(sigma_rot, sigma_tilt);
				if (diffang < sigma_cutoff * biggest_sigma)
				{
					// TODO!!! If tilt is zero then any rot will be OK!!!!!
					//std::cerr<<"Best direction index: "<<idir<<std::endl;
					pointer_dir_nonzeroprior.push_back(idir);
					RFLOAT prior = gaussian1D(diffang, biggest_sigma, 0.);
					sumprior += prior;
					if (isRelax)
					{
						idir_flag[idir] = true;
						RFLOAT my_prior = prior / R_repository_relax.size();
						directions_prior.push_back(my_prior);
						findSymmetryMate(idir, my_prior, pointer_dir_nonzeroprior, directions_prior, idir_flag);
					}
					else
						directions_prior.push_back(prior);
					is_nonzero_pdf = true;
				}


				// Keep track of the nearest direction
				if (diffang < best_ang)
				{
					best_idir = idir;
					best_ang = diffang;
				}
			}
			else if (sigma_rot > 0.)
			{
				Matrix1D<RFLOAT> my_direction, sym_direction;
				RFLOAT sym_rot, sym_tilt;

				// Get the current direction in the loop
				Euler_angles2direction(rot_angles[idir], tilt_angles[idir], my_direction);
				RFLOAT diffang = calculateDeltaRot(my_direction, prior_rot);
				RFLOAT best_diffang = diffang;
				for (int j = 0; j < R_repository.size(); j++)
				{
					sym_direction =  L_repository[j] * (my_direction.transpose() * R_repository[j]).transpose();
					diffang = calculateDeltaRot(sym_direction, prior_rot);

					if (diffang < best_diffang)
						best_diffang = diffang;
				}

				// Only consider differences within sigma_cutoff * sigma_rot
				if (best_diffang < sigma_cutoff * sigma_rot)
				{
					RFLOAT prior = gaussian1D(best_diffang, sigma_rot, 0.);
					pointer_dir_nonzeroprior.push_back(idir);
					directions_prior.push_back(prior);
					sumprior += prior;
					is_nonzero_pdf = true;
				}

				// Keep track of the nearest direction
				if (best_diffang < best_ang)
				{
					best_idir = idir;
					best_ang = diffang;
				}
			}
			else if (sigma_tilt > 0.)
			{
				Matrix1D<RFLOAT> my_direction, sym_direction;
				RFLOAT sym_rot, sym_tilt;

				// Get the current direction in the loop
				Euler_angles2direction(rot_angles[idir], tilt_angles[idir], my_direction);

				// Loop over all symmetry operators to find the operator that brings this direction nearest to the prior
				RFLOAT diffang = ABS(tilt_angles[idir] - prior_tilt);
				if (diffang > 180.)
					diffang = ABS(diffang - 360.);
				RFLOAT best_diffang = diffang;
				for (int j = 0; j < R_repository.size(); j++)
				{
					sym_direction =  L_repository[j] * (my_direction.transpose() * R_repository[j]).transpose();
					Euler_direction2angles(sym_direction, sym_rot, sym_tilt);
					diffang = ABS(sym_tilt - prior_tilt);
					if (diffang > 180.)
						diffang = ABS(diffang - 360.);
					if (diffang < best_diffang)
						best_diffang = diffang;
				}

				// Only consider differences within sigma_cutoff * sigma_tilt
				if (best_diffang < sigma_cutoff * sigma_tilt)
				{
					RFLOAT prior = gaussian1D(best_diffang, sigma_tilt, 0.);
					pointer_dir_nonzeroprior.push_back(idir);
					directions_prior.push_back(prior);
					sumprior += prior;
					is_nonzero_pdf = true;
				}

				// Keep track of the nearest direction
				if (best_diffang < best_ang)
				{
					best_idir = idir;
					best_ang = diffang;
				}
			} // end if any prior involving rot and/or tilt
			else
			{
				// If no prior on the directions: just add all of them
				pointer_dir_nonzeroprior.push_back(idir);
				directions_prior.push_back(1.);
				sumprior += 1.;
				is_nonzero_pdf = true;
			}

			// For priors on deviations from (0,90)-degree (rot,tilt) angles in multi-body refinement
			if (sigma_tilt_from_ninety > 0. && is_nonzero_pdf)
			{
				// Get the current direction in the loop (re-do, as sometimes sigma_rot and sigma_tilt are both zero!
				Matrix1D<RFLOAT> my_direction, best_direction, sym_direction;
				Euler_angles2direction(rot_angles[idir], tilt_angles[idir], my_direction);

				// Loop over all symmetry operators to find the operator that brings this direction nearest to the prior
				RFLOAT best_dotProduct = dotProduct(prior90_direction, my_direction);
				best_direction = my_direction;
				for (int j = 0; j < R_repository.size(); j++)
				{
					sym_direction =  L_repository[j] * (my_direction.transpose() * R_repository[j]).transpose();
					RFLOAT my_dotProduct = dotProduct(prior90_direction, sym_direction);
					if (my_dotProduct > best_dotProduct)
					{
						best_direction = sym_direction;
						best_dotProduct = my_dotProduct;
					}
				}

				// Now that we have the best direction, find the corresponding prior probability
				RFLOAT diffang = ACOSD( dotProduct(best_direction, prior90_direction) );
				if (diffang < -180.)
					diffang = ABS(diffang + 360.);
				else if (diffang > 180.)
					diffang = ABS(diffang - 360.);
				diffang = ABS(diffang);

				long int mypos = pointer_dir_nonzeroprior.size() - 1;
				// Check tilt angle is within 3*sigma_tilt_from_ninety
				if (diffang > sigma_cutoff * sigma_tilt_from_ninety)
				{
					pointer_dir_nonzeroprior.pop_back();
					directions_prior.pop_back();
				}
				else
				{
					RFLOAT prior = gaussian1D(diffang, sigma_tilt_from_ninety, 0.);
					directions_prior[mypos] *= prior;
					sumprior_withsigmafromzero += directions_prior[mypos];
				}
			}
			// Here add the code for relax symmetry to find the symmetry mates

		} // end for idir

		//Normalise the prior probability distribution to have sum 1 over all psi-angles
		for (long int idir = 0; idir < directions_prior.size(); idir++)
		{
			if (sigma_tilt_from_ninety > 0.)
				directions_prior[idir] /= sumprior_withsigmafromzero;
			else
				directions_prior[idir] /= sumprior;
		}

		// If there were no directions at all, just select the single nearest one:
		if (directions_prior.size() == 0)
		{
			pointer_dir_nonzeroprior.push_back(best_idir);
			//std::cerr<<"No direction has been found"<<std::endl;
			if (best_idir < 0)
				REPORT_ERROR("HealpixSampling::selectOrientationsWithNonZeroPriorProbability BUG: best_idir < 0");
			if (isRelax)
			{
				idir_flag[best_idir] = true;
				RFLOAT my_prior = 1./R_repository_relax.size();
				directions_prior.push_back(my_prior);
				findSymmetryMate(best_idir, my_prior, pointer_dir_nonzeroprior, directions_prior, idir_flag);
			}
			else
				directions_prior.push_back(1.);
		}

#ifdef  DEBUG_SAMPLING
		writeNonZeroPriorOrientationsToBild("orients_local.bild", prior_rot, prior_tilt, pointer_dir_nonzeroprior, "0 0 1", 0.023);
		std::cerr << " directions_prior.size()= " << directions_prior.size() << " pointer_dir_nonzeroprior.size()= " << pointer_dir_nonzeroprior.size() << std::endl;
		std::cerr << " sumprior= " << sumprior << std::endl;
		char c;
		std::cerr << "Written orients_local.bild for prior on angles ("<<prior_rot<<","<<prior_tilt<<") Press any key to continue.." << std::endl;
		std::cin >> c;
#endif

	}
	else
	{
		pointer_dir_nonzeroprior.push_back(0);
		directions_prior.push_back(1.);
	}

	// Psi-angles
	pointer_psi_nonzeroprior.clear();
	psi_prior.clear();

	RFLOAT sumprior = 0.;
	RFLOAT sumprior_withsigmafromzero = 0.;
	RFLOAT best_diff = 9999.;
	long int best_ipsi = -999;
	for (long int ipsi = 0; ipsi < psi_angles.size(); ipsi++)
	{
		bool is_nonzero_pdf = false;
		// Sjors 12jul2017: for small tilt-angles, rot-angle may become anything, psi-angle then follows that
		// Therefore, psi-prior may be completely wrong.... The following line would however be a very expensive fix....
		//if (sigma_psi > 0. && prior_tilt > 10.)
		if (sigma_psi > 0.)
		{
			RFLOAT diffpsi = ABS(psi_angles[ipsi] - prior_psi);
			if (diffpsi > 180.)
				diffpsi = ABS(diffpsi - 360.);
			if (do_bimodal_search_psi && (diffpsi > 90.))
				diffpsi = ABS(diffpsi - 180.);

			// Only consider differences within sigma_cutoff * sigma_psi
			if (diffpsi < sigma_cutoff * sigma_psi)
			{
				RFLOAT prior = gaussian1D(diffpsi, sigma_psi, 0.);
				pointer_psi_nonzeroprior.push_back(ipsi);
				psi_prior.push_back(prior);
				sumprior += prior;
				is_nonzero_pdf = true;

				// TMP DEBUGGING
				if (prior == 0.)
				{
					std::cerr << " psi_angles[ipsi]= " << psi_angles[ipsi] << " prior_psi= " << prior_psi << std::endl;
					std::cerr << " diffpsi= " << diffpsi << " sigma_cutoff= " << sigma_cutoff << " sigma_psi= " << sigma_psi << std::endl;
					REPORT_ERROR("prior on psi is zero!");
				}

			}

			// Keep track of the nearest sampling point
			if (diffpsi < best_diff)
			{
				best_ipsi = ipsi;
				best_diff = diffpsi;
			}
		}
		else
		{
			pointer_psi_nonzeroprior.push_back(ipsi);
			psi_prior.push_back(1.);
			sumprior += 1.;
			is_nonzero_pdf = true;
		}

		// For priors on deviations from 0 psi angles in multi-body refinement
		if (sigma_psi_from_zero > 0. && is_nonzero_pdf)
		{
			long int mypos = pointer_psi_nonzeroprior.size() - 1;
			// Check psi angle is within sigma_cutoff*sigma_psi_from_zero
			RFLOAT diff_psi = psi_angles[ipsi];
			if (diff_psi > 180.)
				diff_psi -= 360.;
			else if (diff_psi < -180.)
				diff_psi += 360.;
			diff_psi = ABS(diff_psi);
			if (diff_psi > sigma_cutoff * sigma_psi_from_zero)
			{
				pointer_psi_nonzeroprior.pop_back();
				psi_prior.pop_back();
			}
			else
			{
				RFLOAT prior = gaussian1D(diff_psi, sigma_psi_from_zero, 0.);
				psi_prior[mypos] *= prior;
				sumprior_withsigmafromzero += psi_prior[mypos];
			}
		}

	} // end for ipsi
	// Normalise the prior probability distribution to have sum 1 over all psi-angles
	for (long int ipsi = 0; ipsi < psi_prior.size(); ipsi++)
	{
		if (sigma_psi_from_zero > 0.)
			psi_prior[ipsi] /= sumprior_withsigmafromzero;
		else
			psi_prior[ipsi] /= sumprior;
	}

	// If there were no directions at all, just select the single nearest one:
	if (psi_prior.size() == 0)
	{
		if (best_ipsi < 0)
			REPORT_ERROR("HealpixSampling::selectOrientationsWithNonZeroPriorProbability BUG: best_ipsi < 0");
		pointer_psi_nonzeroprior.push_back(best_ipsi);
		psi_prior.push_back(1.);
	}

#ifdef  DEBUG_SAMPLING
	std::cerr << " psi_angles.size()= " << psi_angles.size() << " psi_step= " << psi_step << std::endl;
	std::cerr << " psi_prior.size()= " << psi_prior.size() << " pointer_psi_nonzeroprior.size()= " << pointer_psi_nonzeroprior.size() << " sumprior= " << sumprior << std::endl;
#endif
	return;
}

void HealpixSampling::findSymmetryMate(long int idir_, RFLOAT prior_,
    		std::vector<int> &pointer_dir_nonzeroprior,
			std::vector<RFLOAT> &directions_prior, std::vector<bool> &idir_flag)
{

	Matrix1D<RFLOAT> my_direction, sym_direction;
	RFLOAT angular_sampling = DEG2RAD(360. / (6 * ROUND(std::pow(2., healpix_order)))) * 2; // Calculate the search radius
	// Direction for the best-matched Healpix index
	Euler_angles2direction(rot_angles[idir_], tilt_angles[idir_], my_direction);

	// Find the best symmetry mates in the HealPix library
	for (int i = 1; i < R_repository_relax.size(); i++)
	{
		int best_direction_index;
		std::vector<int> listpix; // Array with the list of indices for the neighbors
		RFLOAT alpha;  // For Rot
		RFLOAT beta;   // For Theta

		sym_direction =  L_repository_relax[i] * (my_direction.transpose() * R_repository_relax[i]).transpose();
		Euler_direction2angles(sym_direction, alpha, beta);
		
		alpha = DEG2RAD(alpha);
		beta  = DEG2RAD(beta);
		pointing prior_direction_pointing(beta, alpha); // Object required by healpix function
     	healpix_base.query_disc(prior_direction_pointing, angular_sampling, listpix); // Search healpix for closest indices
     	best_direction_index = listpix[0];
     	// If there are more than one neighbors then select the best
		if (listpix.size() > 1)
     	{
			Matrix1D<RFLOAT> current_direction;
			Euler_angles2direction(rot_angles[best_direction_index], tilt_angles[best_direction_index], current_direction);
     		RFLOAT best_dotProduct = dotProduct(sym_direction, current_direction);
     		for (long int j = 1; j < listpix.size(); j++)
     		{
     			int current_index = listpix[j];
     			// Assuming sigma_tilt and sigma_rot are set
     			// Get the current direction
     			Euler_angles2direction(rot_angles[current_index], tilt_angles[current_index], current_direction);
     			RFLOAT my_dotProduct = dotProduct(sym_direction, current_direction);
     			if (my_dotProduct > best_dotProduct && 	idir_flag[current_index] != true)
     			{
     				best_direction_index = current_index;
     				best_dotProduct = my_dotProduct;
     			}
     		}
     	}

		// Now we have the best symmetry mate index
     	pointer_dir_nonzeroprior.push_back(best_direction_index);
		directions_prior.push_back(prior_);
		idir_flag[best_direction_index] = true;
	}
}

void HealpixSampling::selectOrientationsWithNonZeroPriorProbabilityFor3DHelicalReconstruction(
		RFLOAT prior_rot, RFLOAT prior_tilt, RFLOAT prior_psi,
		RFLOAT sigma_rot, RFLOAT sigma_tilt, RFLOAT sigma_psi,
		std::vector<int> &pointer_dir_nonzeroprior, std::vector<RFLOAT> &directions_prior,
		std::vector<int> &pointer_psi_nonzeroprior, std::vector<RFLOAT> &psi_prior,
		bool do_auto_refine_local_searches,
		RFLOAT prior_psi_flip_ratio,
		RFLOAT prior_rot_flip_ratio,  // KThurber
		RFLOAT sigma_cutoff)
{
	// Helical references are always along Z axis in 3D helical reconstructions
	// Therefore tilt priors are usually not far from 90 degrees
	// Tilt~0 problem: if tilt=0 or 180 degrees, a change in rot is equal to some changes in psi
	// If tilt priors are estimated around 0 or 180 degrees for a segment, this segment must be rubbish
	// So tilt~0 problem does not impact helical reconstruction
	// Because if user provides this rubbish segment, just let it have any orientations and it will give out rubbish
	// Just throw a warning message for this rubbish segment
	// If tilt~0 problem is ignored, I can implement 2D Gaussian priors for rot-tilt pairs
	// This can be more accurate than previous implementation
	// It also saves time when sigma_rot is much larger than sigma_tilt (during local searches in 3D auto-refine)

	RFLOAT prior_psi_flip_ratio_thres_min = 0.01;
	RFLOAT prior_rot_flip_ratio_thres_min = 0.01; 	// KThurber

	pointer_dir_nonzeroprior.clear();
	directions_prior.clear();

	if (is_3D)
	{
		// If tilt prior is less than 20 or larger than 160 degrees, print a warning message
		//if (fabs(((prior_tilt / 180.) - ROUND(prior_tilt / 180.)) * 180.) < 20.)
		//{
		//	std::cerr << " WARNING: A helical segment is found with tilt prior= " << prior_tilt
		//			<< " degrees. It will probably impact searches of orientations in 3D helical reconstruction."<< std::endl;
		//}

		// Loop over all directions
		RFLOAT sumprior = 0.;
		// Keep track of the closest distance to prevent 0 orientations
		RFLOAT best_ang = 9999.;
		long int best_idir = -999;
		for (long int idir = 0; idir < rot_angles.size(); idir++)
		{
			// Any prior involving BOTH rot and tilt.
			if ( (sigma_rot > 0.) && (sigma_tilt > 0.) )
			{
				Matrix1D<RFLOAT> prior_direction, my_direction, sym_direction, best_direction;

				// Get the direction of the prior
				Euler_angles2direction(prior_rot, prior_tilt, prior_direction);

				// Get the current direction in the loop
				Euler_angles2direction(rot_angles[idir], tilt_angles[idir], my_direction);

				// Loop over all symmetry operators to find the operator that brings this direction nearest to the prior
				RFLOAT best_dotProduct = dotProduct(prior_direction, my_direction);
				best_direction = my_direction;
				for (int j = 0; j < R_repository.size(); j++)
				{
					sym_direction =  L_repository[j] * (my_direction.transpose() * R_repository[j]).transpose();
					RFLOAT my_dotProduct = dotProduct(prior_direction, sym_direction);
					if (my_dotProduct > best_dotProduct)
					{
						best_direction = sym_direction;
						best_dotProduct = my_dotProduct;
					}
				}

				if (!do_auto_refine_local_searches)
				{
					// Assume tilt = (0, +180)
					// TODO: Check if "(tilt_angles[idir] > 0.01) && (tilt_angles[idir] < 179.99)" is needed
					//if (prior_psi_flip_ratio > prior_psi_flip_ratio_thres_min)
					if (prior_psi_flip_ratio > -1.)		// KThurber above line changed to primarily dummy if
					{
						Matrix1D<RFLOAT> my_direction2, sym_direction2, best_direction2;

						// Get the current direction in the loop
						Euler_angles2direction(rot_angles[idir], (180. - tilt_angles[idir]), my_direction2);

						// Loop over all symmetry operators to find the operator that brings this direction nearest to the prior
						RFLOAT best_dotProduct2 = dotProduct(prior_direction, my_direction2);
						best_direction2 = my_direction2;
						for (int j = 0; j < R_repository.size(); j++)
						{
							sym_direction2 =  L_repository[j] * (my_direction2.transpose() * R_repository[j]).transpose();
							RFLOAT my_dotProduct2 = dotProduct(prior_direction, sym_direction2);
							if (my_dotProduct2 > best_dotProduct2)
							{
								best_direction2 = sym_direction2;
								best_dotProduct2 = my_dotProduct2;
							}
						}

						if (best_dotProduct2 > best_dotProduct)
						{
							best_dotProduct = best_dotProduct2;
							best_direction = best_direction2;
						}
					}
				}

				// Calculate the differences from sym_rot, sym_tilt to prior_rot, prior_tilt
				RFLOAT sym_rot, sym_tilt, diff_rot, diff_tilt, diffang;
				Euler_direction2angles(best_direction, sym_rot, sym_tilt);
				diff_rot = ABS(sym_rot - prior_rot);
				if (diff_rot > 180.)
					diff_rot = ABS(diff_rot - 360.);

				// KThurber begin add
				bool is_rot_flipped = false;
				if (!do_auto_refine_local_searches)
				{
					if (diff_rot > 90.)
					{
						diff_rot = ABS(diff_rot - 180.);
						is_rot_flipped = true;
					}
				}
				// KThurber end add

				diff_tilt = ABS(sym_tilt - prior_tilt);
				if (diff_tilt > 180.)
					diff_tilt = ABS(diff_tilt - 360.);
				diffang = sqrt((diff_rot * diff_rot) + (diff_tilt * diff_tilt));
				if ( (diff_rot < sigma_cutoff * sigma_rot) && (diff_tilt < sigma_cutoff * sigma_tilt) )
				{
					// TODO!!! If tilt is zero then any rot will be OK!!!!!
					RFLOAT prior = gaussian2D(diff_rot, diff_tilt, sigma_rot, sigma_tilt, 0.);
					// KThurber begin add
					if (!do_auto_refine_local_searches)
					{
						if (is_rot_flipped)
							prior *= prior_rot_flip_ratio;
						else
							prior *= (1. - prior_rot_flip_ratio);
					}
					// KThurber end add
					pointer_dir_nonzeroprior.push_back(idir);
					directions_prior.push_back(prior);
					sumprior += prior;
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
					std::cout << "rot & tilt OK, diffang = " << diffang << std::endl;
					std::cout << " rot, tilt = " << rot_angles[idir] << ", " << tilt_angles[idir] << std::endl;
#endif
				}
				else
				{
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
					std::cout << "rot & tilt FAILED, diffang = " << diffang << std::endl;
					std::cout << " rot, tilt = " << rot_angles[idir] << ", " << tilt_angles[idir] << std::endl;
#endif
				}

				// Keep track of the nearest direction
				if (diffang < best_ang)
				{
					best_idir = idir;
					best_ang = diffang;
				}

			}
			else if (sigma_tilt > 0.)
			{
				Matrix1D<RFLOAT> my_direction, sym_direction;
				RFLOAT sym_rot, sym_tilt;

				// Get the current direction in the loop
				Euler_angles2direction(rot_angles[idir], tilt_angles[idir], my_direction);

				// Loop over all symmetry operators to find the operator that brings this direction nearest to the prior
				RFLOAT diffang = ABS(tilt_angles[idir] - prior_tilt);
				if (diffang > 180.)
					diffang = ABS(diffang - 360.);
				RFLOAT best_tilt = tilt_angles[idir];
				RFLOAT best_diffang = diffang;
				for (int j = 0; j < R_repository.size(); j++)
				{
					sym_direction =  L_repository[j] * (my_direction.transpose() * R_repository[j]).transpose();
					Euler_direction2angles(sym_direction, sym_rot, sym_tilt);
					diffang = ABS(sym_tilt - prior_tilt);
					if (diffang > 180.)
						diffang = ABS(diffang - 360.);
					if (diffang < best_diffang)
					{
						best_diffang = diffang;
						best_tilt = sym_tilt;
					}
				}

				if (!do_auto_refine_local_searches)
				{
					// Assume tilt = (0, +180)
					// TODO: Check if "(tilt_angles[idir] > 0.01) && (tilt_angles[idir] < 179.99)" is needed
					//if (prior_psi_flip_ratio > prior_psi_flip_ratio_thres_min)
					if (prior_psi_flip_ratio > -1.)		// KThurber above line changed to primarily dummy if
					{
						Matrix1D<RFLOAT> my_direction2, sym_direction2;

						// Get the current direction in the loop
						sym_tilt = 180. - tilt_angles[idir]; // Shaoda want the prior on tilt, centered around 90 degrees, so P(87) == P(93)
						Euler_angles2direction(rot_angles[idir], sym_tilt, my_direction2);

						// Loop over all symmetry operators to find the operator that brings this direction nearest to the prior
						diffang = ABS(sym_tilt - prior_tilt);
						if (diffang > 180.)
							diffang = ABS(diffang - 360.);
						if (diffang < best_diffang)
						{
							best_diffang = diffang;
							best_tilt = sym_tilt;
						}
						for (int j = 0; j < R_repository.size(); j++)
						{
							sym_direction2 =  L_repository[j] * (my_direction2.transpose() * R_repository[j]).transpose();
							Euler_direction2angles(sym_direction2, sym_rot, sym_tilt);
							diffang = ABS(sym_tilt - prior_tilt);
							if (diffang > 180.)
								diffang = ABS(diffang - 360.);
							if (diffang < best_diffang)
							{
								best_diffang = diffang;
								best_tilt = sym_tilt;
							}
						}
					}
				}

				// Only consider differences within sigma_cutoff * sigma_tilt
				if (best_diffang < sigma_cutoff * sigma_tilt)
				{
					RFLOAT prior = gaussian1D(best_diffang, sigma_tilt, 0.);
					pointer_dir_nonzeroprior.push_back(idir);
					directions_prior.push_back(prior);
					sumprior += prior;
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
					std::cout << "rot & tilt OK, diffang = " << diffang << std::endl;
					std::cout << " rot, tilt = " << rot_angles[idir] << ", " << tilt_angles[idir] << std::endl;
#endif
				}
				else
				{
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
					std::cout << "rot & tilt FAILED, diffang = " << diffang << std::endl;
					std::cout << " rot, tilt = " << rot_angles[idir] << ", " << tilt_angles[idir] << std::endl;
#endif
				}

				// Keep track of the nearest direction
				if (best_diffang < best_ang)
				{
					best_idir = idir;
					best_ang = diffang;
				}

			}
			else if (sigma_rot > 0.)
			{
				REPORT_ERROR("healpix_sampling.cpp::selectOrientationsWithNonZeroPriorProbabilityFor3DHelicalReconstruction() BUG: Rot but not tilt prior exists! It is not reasonable for 3D helical reconstruction!");
			}
			else
			{
				// If no prior on the directions: just add all of them
				pointer_dir_nonzeroprior.push_back(idir);
				directions_prior.push_back(1.);
				sumprior += 1.;
			}

		} // end for idir

		//Normalise the prior probability distribution to have sum 1 over all psi-angles
		for (long int idir = 0; idir < directions_prior.size(); idir++)
			directions_prior[idir] /= sumprior;

		// If there were no directions at all, just select the single nearest one:
		if (directions_prior.size() == 0)
		{
			if (best_idir < 0)
				REPORT_ERROR("HealpixSampling::selectOrientationsWithNonZeroPriorProbability BUG: best_idir < 0");
			pointer_dir_nonzeroprior.push_back(best_idir);
			directions_prior.push_back(1.);
		}

#ifdef  DEBUG_SAMPLING
		writeNonZeroPriorOrientationsToBild("orients_local.bild", prior_rot, prior_tilt, pointer_dir_nonzeroprior, "0 0 1", 0.023);
		std::cerr << " directions_prior.size()= " << directions_prior.size() << " pointer_dir_nonzeroprior.size()= " << pointer_dir_nonzeroprior.size() << std::endl;
		std::cerr << " sumprior= " << sumprior << std::endl;
		char c;
		std::cerr << "Written orients_local.bild for prior on angles ("<<prior_rot<<","<<prior_tilt<<") Press any key to continue.." << std::endl;
		std::cin >> c;
#endif

	}
	else
	{
		pointer_dir_nonzeroprior.push_back(0);
		directions_prior.push_back(1.);
	}

	// Psi-angles
	pointer_psi_nonzeroprior.clear();
	psi_prior.clear();

	RFLOAT sumprior = 0.;
	RFLOAT best_diff = 9999.;
	long int best_ipsi = -999;
	bool is_psi_flipped = false;
	for (long int ipsi = 0; ipsi < psi_angles.size(); ipsi++)
	{
		if (sigma_psi > 0.)
		{
			RFLOAT diffpsi = ABS(psi_angles[ipsi] - prior_psi);
			if (diffpsi > 180.)
				diffpsi = ABS(diffpsi - 360.);
			if (!do_auto_refine_local_searches)
			{
				if ( (prior_psi_flip_ratio > prior_psi_flip_ratio_thres_min) && (diffpsi > 90.))
				{
					diffpsi = ABS(diffpsi - 180.);
					is_psi_flipped = true;
				}
			}

			// Only consider differences within sigma_cutoff * sigma_psi
			if (diffpsi < sigma_cutoff * sigma_psi)
			{
				RFLOAT prior = gaussian1D(diffpsi, sigma_psi, 0.);
				if (!do_auto_refine_local_searches)
				{
					if (is_psi_flipped)
						prior *= prior_psi_flip_ratio;
					else
						prior *= (1. - prior_psi_flip_ratio);
				}
				pointer_psi_nonzeroprior.push_back(ipsi);
				psi_prior.push_back(prior);
				sumprior += prior;
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
				std::cout << "psi OK, diffang = " << diffpsi << std::endl;
				std::cout << " psi = " << psi_angles[ipsi] << std::endl;
#endif

				// TMP DEBUGGING
				if (prior == 0.)
				{
					std::cerr << " psi_angles[ipsi]= " << psi_angles[ipsi] << " prior_psi= " << prior_psi << std::endl;
					std::cerr << " diffpsi= " << diffpsi << " sigma_cutoff= " << sigma_cutoff << " sigma_psi= " << sigma_psi << std::endl;
					REPORT_ERROR("prior on psi is zero!");
				}

			}
			else
			{
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
				std::cout << "psi FAILED, diffang = " << diffpsi << std::endl;
				std::cout << " psi = " << psi_angles[ipsi] << std::endl;
#endif
			}
			// Keep track of the nearest sampling point
			if (diffpsi < best_diff)
			{
				best_ipsi = ipsi;
				best_diff = diffpsi;
			}
		}
		else
		{
			pointer_psi_nonzeroprior.push_back(ipsi);
			psi_prior.push_back(1.);
			sumprior += 1.;
		}
	}
	// Normalise the prior probability distribution to have sum 1 over all psi-angles
	for (long int ipsi = 0; ipsi < psi_prior.size(); ipsi++)
		psi_prior[ipsi] /= sumprior;

	// If there were no directions at all, just select the single nearest one:
	if (psi_prior.size() == 0)
	{
		if (best_ipsi < 0)
			REPORT_ERROR("HealpixSampling::selectOrientationsWithNonZeroPriorProbability BUG: best_ipsi < 0");
		pointer_psi_nonzeroprior.push_back(best_ipsi);
		psi_prior.push_back(1.);
	}

#ifdef  DEBUG_SAMPLING
	std::cerr << " psi_angles.size()= " << psi_angles.size() << " psi_step= " << psi_step << std::endl;
	std::cerr << " psi_prior.size()= " << psi_prior.size() << " pointer_psi_nonzeroprior.size()= " << pointer_psi_nonzeroprior.size() << " sumprior= " << sumprior << std::endl;
#endif

}

FileName HealpixSampling::symmetryGroup()
{
	return fn_sym;
}

long int HealpixSampling::getHealPixIndex(long int idir)
{
#ifdef DEBUG_CHECKSIZES
	if (idir >= directions_ipix.size())
	{
		std::cerr<< "idir= "<<idir<<" directions_ipix.size()= "<< directions_ipix.size() <<std::endl;
		REPORT_ERROR("idir >= directions_ipix.size()");
	}
#endif
	return directions_ipix[idir];
}

void HealpixSampling::checkDirection(RFLOAT &rot, RFLOAT &tilt)
{

	// The geometrical considerations about the symmetry below require that rot = [-180,180] and tilt [0,180]

	// The following was incorrect?!
	if (tilt < 0.)
	{
		tilt = -tilt;
		rot += 180.;
	}

	bool is_ok = false;
	while (!is_ok)
	{
		if (rot > 180.)
			rot -= 360.;
		else if (rot < -180.)
			rot += 360.;
		else
			is_ok = true;
	}

}

void HealpixSampling::getDirectionFromHealPix(long int ipix, RFLOAT &rot, RFLOAT &tilt)
{
    // this one always has to be double (also for SINGLE_PRECISION CALCULATIONS) for call to external library
    double zz, phi;
    healpix_base.pix2ang_z_phi(ipix, zz, phi);
    rot = RAD2DEG(phi);
    tilt = ACOSD(zz);

    // The geometrical considerations about the symmetry below require that rot = [-180,180] and tilt [0,180]
    checkDirection(rot, tilt);

}

RFLOAT HealpixSampling::getTranslationalSampling(int adaptive_oversampling)
{
	return offset_step / std::pow(2., adaptive_oversampling);
}

RFLOAT HealpixSampling::getHelicalTranslationalSampling(int adaptive_oversampling)
{
	return helical_offset_step / std::pow(2., adaptive_oversampling);
}

RFLOAT HealpixSampling::getAngularSampling(int adaptive_oversampling)
{
	if (is_3D)
	{
		int order =  healpix_order + adaptive_oversampling;
		return 360. / (6 * ROUND(std::pow(2., order)));
	}
	else
		return psi_step / std::pow(2., adaptive_oversampling);
}

long int HealpixSampling::NrDirections(int oversampling_order, const std::vector<int> *pointer_dir_nonzeroprior)
{
	long int mysize = (pointer_dir_nonzeroprior != NULL && (*pointer_dir_nonzeroprior).size() > 0) ?
			(*pointer_dir_nonzeroprior).size() : rot_angles.size();
	if (oversampling_order == 0)
		return mysize;
	else
		return ROUND(std::pow(2., oversampling_order * 2)) * mysize;
}

long int HealpixSampling::NrPsiSamplings(int oversampling_order, const std::vector<int> *pointer_psi_nonzeroprior)
{

	long int mysize = (pointer_psi_nonzeroprior != NULL && (*pointer_psi_nonzeroprior).size() > 0) ?
			(*pointer_psi_nonzeroprior).size() : psi_angles.size();
	if (oversampling_order == 0)
		return mysize;
	else
		return ROUND(std::pow(2., oversampling_order)) * mysize;
}

long int HealpixSampling::NrTranslationalSamplings(int oversampling_order)
{
	if (oversampling_order == 0)
		return translations_x.size();
	else
	{
		if (is_3d_trans)
			return ROUND(std::pow(2., oversampling_order * 3)) * translations_x.size();
		else
			return ROUND(std::pow(2., oversampling_order * 2)) * translations_x.size();
	}
}

long int HealpixSampling::NrSamplingPoints(int oversampling_order,
		const std::vector<int> *pointer_dir_nonzeroprior,
		const std::vector<int> *pointer_psi_nonzeroprior)
{
	return NrDirections(oversampling_order, pointer_dir_nonzeroprior) *
		   NrPsiSamplings(oversampling_order, pointer_psi_nonzeroprior) *
		   NrTranslationalSamplings(oversampling_order);
}

/* How often is each orientation oversampled? */
int HealpixSampling::oversamplingFactorOrientations(int oversampling_order)
{
	if (is_3D)
		return ROUND(std::pow(2., oversampling_order * 3));
	else
		return ROUND(std::pow(2., oversampling_order));
}

/* How often is each translation oversampled? */
int HealpixSampling::oversamplingFactorTranslations(int oversampling_order)
{
	if (is_3d_trans)
		return ROUND(std::pow(2., oversampling_order * 3));
	else
		return ROUND(std::pow(2., oversampling_order * 2));
}


void HealpixSampling::getDirection(long int idir, RFLOAT &rot, RFLOAT &tilt)
{
#ifdef DEBUG_CHECKSIZES
	if (idir >= rot_angles.size())
	{
		std::cerr<< "idir= "<<idir<<" rot_angles.size()= "<< rot_angles.size() <<std::endl;
		REPORT_ERROR("idir >= rot_angles.size()");
	}
#endif

	rot  = rot_angles[idir];
	tilt = tilt_angles[idir];
}

void HealpixSampling::getPsiAngle(long int ipsi, RFLOAT &psi)
{
#ifdef DEBUG_CHECKSIZES
	if (ipsi >= psi_angles.size())
	{
		std::cerr<< "ipsi= "<<ipsi<<" psi_angles.size()= "<< psi_angles.size() <<std::endl;
		REPORT_ERROR("ipsi >= psi_angles.size()");
	}
#endif
	psi = psi_angles[ipsi];
}

void HealpixSampling::getTranslationInPixel(long int itrans, RFLOAT my_pixel_size, RFLOAT &trans_x, RFLOAT &trans_y, RFLOAT &trans_z)
{
#ifdef DEBUG_CHECKSIZES
if (itrans >= translations_x.size())
{
	std::cerr<< "itrans= "<<itrans<<" translations_x.size()= "<< translations_x.size() <<std::endl;
	REPORT_ERROR("itrans >= translations_x.size()");
}
#endif
	trans_x = translations_x[itrans] / my_pixel_size;
	trans_y = translations_y[itrans] / my_pixel_size;
	if (is_3d_trans)
		trans_z = translations_z[itrans] / my_pixel_size;
}


long int HealpixSampling::getPositionSamplingPoint(int iclass, long int idir, long int ipsi, long int itrans)
{
	return iclass * rot_angles.size() * psi_angles.size() * translations_x.size()
		+ idir * psi_angles.size() * translations_x.size()
		+ ipsi * translations_x.size() + itrans;
}

long int HealpixSampling::getPositionOversampledSamplingPoint(long int ipos, int oversampling_order, int iover_rot, int iover_trans)
{
	if (oversampling_order == 0)
		return ipos;
	else
	{
		int nr_over_orient = oversamplingFactorOrientations(oversampling_order);
		int nr_over_trans = oversamplingFactorTranslations(oversampling_order);
		return ipos * nr_over_orient * nr_over_trans + nr_over_trans * iover_rot + iover_trans;
	}

}

void HealpixSampling::getTranslationsInPixel(long int itrans, int oversampling_order, RFLOAT my_pixel_size,
		std::vector<RFLOAT > &my_translations_x,
		std::vector<RFLOAT > &my_translations_y,
		std::vector<RFLOAT > &my_translations_z,
		bool do_helical_refine)
{

#ifdef DEBUG_CHECKSIZES
	if (itrans >= translations_x.size())
	{
		std::cerr<< "itrans= "<<itrans<<" translations_x.size()= "<< translations_x.size() <<std::endl;
		REPORT_ERROR("itrans >= translations_x.size()");
	}
#endif
	my_translations_x.clear();
	my_translations_y.clear();
	my_translations_z.clear();
	if (oversampling_order == 0)
	{
		my_translations_x.push_back(translations_x[itrans] / my_pixel_size);
		my_translations_y.push_back(translations_y[itrans] / my_pixel_size);
		if (is_3d_trans)
			my_translations_z.push_back(translations_z[itrans] / my_pixel_size);
	}
	else
	{
		int nr_oversamples = ROUND(std::pow(2., oversampling_order));
		// DEBUG
		if ( (nr_oversamples < 1) )
		{
			std::cerr << "oversampling_order= " << oversampling_order << " nr_oversamples= " << nr_oversamples << std::endl;
			REPORT_ERROR("HealpixSampling::getTranslations BUG %% 'nr_oversamples' should be a positive integer!");
		}

		// Helical refinement
		RFLOAT h_step = offset_step;
		if (do_helical_refine)
		{
			h_step = helical_offset_step;
		}
		if (h_step < 0.)
		{
			std::cerr << "helical_offset_step (h_step)= " << h_step << std::endl;
			REPORT_ERROR("HealpixSampling::getTranslations BUG %% 'helical_offset_step (h_step)' should be positive!");
		}

		RFLOAT over_xoff = 0., over_yoff = 0., over_zoff = 0.;
		for (int itrans_overx = 0; itrans_overx < nr_oversamples; itrans_overx++)
		{
			if ( (do_helical_refine) && (!is_3d_trans) ) // Helical reconstruction with 2D segments
				over_xoff = translations_x[itrans] - 0.5 * h_step + (0.5 + itrans_overx) * h_step / nr_oversamples;
			else
				over_xoff = translations_x[itrans] - 0.5 * offset_step + (0.5 + itrans_overx) * offset_step / nr_oversamples;
			for (int itrans_overy = 0; itrans_overy < nr_oversamples; itrans_overy++)
			{
				over_yoff = translations_y[itrans] - 0.5 * offset_step + (0.5 + itrans_overy) * offset_step / nr_oversamples;
				if (is_3d_trans)
				{
					for (int itrans_overz = 0; itrans_overz < nr_oversamples; itrans_overz++)
					{
						if (do_helical_refine) // Helical reconstruction in 3D subtomogram averaging
							over_zoff = translations_z[itrans] - 0.5 * h_step + (0.5 + itrans_overz) * h_step / nr_oversamples;
						else
							over_zoff = translations_z[itrans] - 0.5 * offset_step + (0.5 + itrans_overz) * offset_step / nr_oversamples;

						my_translations_x.push_back(over_xoff / my_pixel_size);
						my_translations_y.push_back(over_yoff / my_pixel_size);
						my_translations_z.push_back(over_zoff / my_pixel_size);
					}
				}
				else
				{
					my_translations_x.push_back(over_xoff / my_pixel_size);
					my_translations_y.push_back(over_yoff / my_pixel_size);
				}
			}
		}

		// This could rarely happen. Just in case all over-sampled xoff are excluded in helical refinement...
		// AVOID THIS BY CHOOSING AN INITIAL ANGULAR SAMPLING FINER THAN HALF OF THE TWIST !!!!!!
		if ( (do_helical_refine) && (my_translations_x.size() < 1) )
		{
			my_translations_x.push_back(translations_x[itrans] / my_pixel_size);
			my_translations_y.push_back(translations_y[itrans] / my_pixel_size);
			if (is_3d_trans)
				my_translations_z.push_back(translations_z[itrans] / my_pixel_size);
		}
	}

	if (ABS(random_perturbation) > 0.)
	{
		RFLOAT myperturb = random_perturbation * offset_step / my_pixel_size;
		// Oct31,2015 - Shaoda - TODO: Please consider this!!!
		RFLOAT myperturb_helical = random_perturbation * helical_offset_step / my_pixel_size;
		for (int iover = 0; iover < my_translations_x.size(); iover++)
		{
			// If doing helical refinement, DONT put perturbation onto translations along helical axis???
			if ( (do_helical_refine) && (!is_3d_trans) ) // Helical reconstruction with 2D segments
				my_translations_x[iover] += myperturb_helical;
			else
				my_translations_x[iover] += myperturb;
			my_translations_y[iover] += myperturb;
			if (is_3d_trans)
			{
				if (do_helical_refine)
					my_translations_z[iover] += myperturb_helical; // Helical reconstruction in 3D subtomogram averaging
				else
					my_translations_z[iover] += myperturb;
			}
		}
	}

}

void HealpixSampling::getOrientations(long int idir, long int ipsi, int oversampling_order,
		std::vector<RFLOAT > &my_rot, std::vector<RFLOAT > &my_tilt, std::vector<RFLOAT > &my_psi,
		std::vector<int> &pointer_dir_nonzeroprior, std::vector<RFLOAT> &directions_prior,
		std::vector<int> &pointer_psi_nonzeroprior, std::vector<RFLOAT> &psi_prior)
{
	my_rot.clear();
	my_tilt.clear();
	my_psi.clear();
	long int my_idir, my_ipsi;
	if (pointer_dir_nonzeroprior.size() > idir && pointer_psi_nonzeroprior.size() > ipsi)
	{
		// nonzeroprior vectors have been initialised, so use priors!
		my_idir = pointer_dir_nonzeroprior[idir];
		my_ipsi = pointer_psi_nonzeroprior[ipsi];
	}
	else
	{
		// no priors
		my_idir = idir;
		my_ipsi = ipsi;
	}

#ifdef DEBUG_CHECKSIZES
		if (my_idir >= rot_angles.size())
		{
			std::cerr<< "my_idir= "<<my_idir<<" rot_angles.size()= "<< rot_angles.size() <<std::endl;
			REPORT_ERROR("my_idir >= rot_angles.size()");
		}
		if (my_ipsi >= psi_angles.size())
		{
			std::cerr<< "my_ipsi= "<<my_ipsi<<" psi_angles.size()= "<< psi_angles.size() <<std::endl;
			REPORT_ERROR("my_ipsi >= psi_angles.size()");
		}
#endif

	if (oversampling_order == 0)
	{
		my_rot.push_back(rot_angles[my_idir]);
		my_tilt.push_back(tilt_angles[my_idir]);
		my_psi.push_back(psi_angles[my_ipsi]);
	}
	else if (!is_3D)
	{
		// for 2D sampling, only push back oversampled psi rotations
		pushbackOversampledPsiAngles(my_ipsi, oversampling_order, 0., 0., my_rot, my_tilt, my_psi);
	}
	else
	{
		// Set up oversampled grid for 3D sampling
		Healpix_Base HealPixOver(oversampling_order + healpix_order, NEST);
		int fact = HealPixOver.Nside()/healpix_base.Nside();
		int x, y, face;
		RFLOAT rot, tilt;
		// Get x, y and face for the original, coarse grid
		long int ipix = directions_ipix[my_idir];
		healpix_base.nest2xyf(ipix, x, y, face);
		// Loop over the oversampled Healpix pixels on the fine grid
		for (int j = fact * y; j < fact * (y+1); ++j)
		{
			for (int i = fact * x; i < fact * (x+1); ++i)
			{
				long int overpix = HealPixOver.xyf2nest(i, j, face);
                                // this one always has to be double (also for SINGLE_PRECISION CALCULATIONS) for call to external library
				double zz, phi;
				HealPixOver.pix2ang_z_phi(overpix, zz, phi);
				rot = RAD2DEG(phi);
				tilt = ACOSD(zz);

				// The geometrical considerations about the symmetry below require that rot = [-180,180] and tilt [0,180]
				checkDirection(rot, tilt);

				pushbackOversampledPsiAngles(my_ipsi, oversampling_order, rot, tilt, my_rot, my_tilt, my_psi);
			}
		}
	}


	// Random perturbation
	if (ABS(random_perturbation) > 0.)
	{
		RFLOAT myperturb = random_perturbation * getAngularSampling();
		for (int iover = 0; iover < my_rot.size(); iover++)
		{
			if (is_3D)
			{
				Matrix2D<RFLOAT> A(3,3), R(3,3);
				Euler_angles2matrix(my_rot[iover],
									my_tilt[iover],
									my_psi[iover],
									A);
				Euler_angles2matrix(myperturb, myperturb, myperturb, R);
				A = A * R;
				Euler_matrix2angles(A,
									my_rot[iover],
									my_tilt[iover],
									my_psi[iover]);
			}
			else
			{
				my_psi[iover] += myperturb;
			}
		}
	}



}


void HealpixSampling::pushbackOversampledPsiAngles(long int ipsi, int oversampling_order,
		RFLOAT rot, RFLOAT tilt, std::vector<RFLOAT> &oversampled_rot,
		std::vector<RFLOAT> &oversampled_tilt, std::vector<RFLOAT> &oversampled_psi)
{

	if (oversampling_order == 0)
	{
		oversampled_rot.push_back(rot);
		oversampled_tilt.push_back(tilt);
		oversampled_psi.push_back(psi_angles[ipsi]);
	}
	else
	{
		int nr_ipsi_over = ROUND(std::pow(2., oversampling_order));
		for (int ipsi_over = 0; ipsi_over < nr_ipsi_over; ipsi_over++)
		{
			RFLOAT overpsi = psi_angles[ipsi] - 0.5 * psi_step + (0.5 + ipsi_over) * psi_step / nr_ipsi_over;
			oversampled_rot.push_back(rot);
			oversampled_tilt.push_back(tilt);
                        if (!is_3D && overpsi>180.)
                            overpsi-=360.;
                        oversampled_psi.push_back(overpsi);
		}
	}

}

/* Calculate an angular distance between two sets of Euler angles */
RFLOAT HealpixSampling::calculateAngularDistance(RFLOAT rot1, RFLOAT tilt1, RFLOAT psi1,
		RFLOAT rot2, RFLOAT tilt2, RFLOAT psi2)
{

	if (is_3D)
	{
		Matrix1D<RFLOAT>  direction1(3), direction1p(3), direction2(3);
		Euler_angles2direction(rot1, tilt1, direction1);
		Euler_angles2direction(rot2, tilt2, direction2);

		// Find the symmetry operation where the Distance based on Euler axes is minimal
		RFLOAT min_axes_dist = 3600.;
		RFLOAT rot2p, tilt2p, psi2p;
		Matrix2D<RFLOAT> E1, E2;
		Matrix1D<RFLOAT> v1, v2;
		for (int j = 0; j < R_repository.size(); j++)
		{

			Euler_apply_transf(L_repository[j], R_repository[j], rot2, tilt2, psi2, rot2p, tilt2p, psi2p);

			// Distance based on Euler axes
			Euler_angles2matrix(rot1, tilt1, psi1, E1);
			Euler_angles2matrix(rot2p, tilt2p, psi2p, E2);
			RFLOAT axes_dist = 0;
			for (int i = 0; i < 3; i++)
			{
				E1.getRow(i, v1);
				E2.getRow(i, v2);
				axes_dist += ACOSD(CLIP(dotProduct(v1, v2), -1., 1.));
			}
			axes_dist /= 3.;

			if (axes_dist < min_axes_dist)
				min_axes_dist = axes_dist;

		}// for all symmetry operations j

		return min_axes_dist;
	}
	else
	{
		RFLOAT diff = ABS(psi2 - psi1);
		return realWRAP(diff, 0., 360.);
	}
}

void HealpixSampling::writeBildFileOrientationalDistribution(MultidimArray<RFLOAT> &pdf_direction,
		FileName &fn_bild, RFLOAT R, RFLOAT offset,
		const Matrix2D<RFLOAT> *Aorient, const Matrix1D<RFLOAT> *Acom,
		RFLOAT Rmax_frac, RFLOAT width_frac)
{
	if (!is_3D)
		return;

	if (XSIZE(pdf_direction) != rot_angles.size())
	{
		std::cerr << " XSIZE(pdf_direction)= " << XSIZE(pdf_direction) << " rot_angles.size()= " << rot_angles.size() << std::endl;
		REPORT_ERROR("HealpixSampling::writeBildFileOrientationalDistribution XSIZE(pdf_direction) != rot_angles.size()!");
	}


	RFLOAT pdfmax, pdfmin, pdfmean, pdfsigma;
	pdf_direction.computeStats(pdfmean, pdfsigma, pdfmin, pdfmax);

	std::ofstream fh_bild;
    fh_bild.open(fn_bild.c_str(), std::ios::out);
    if (!fh_bild)
    	REPORT_ERROR("HealpixSampling::writeBildFileOrientationalDistribution: cannot open " + fn_bild);

    // 2 * PI * R = 360 degrees, 2*radius should cover angular sampling at width_frac=1
    RFLOAT width = width_frac * PI*R*(getAngularSampling()/360.);
    Matrix1D<RFLOAT> v(3);

    for (long int iang = 0; iang < rot_angles.size(); iang++)
    {
     	RFLOAT pdf = DIRECT_A1D_ELEM(pdf_direction, iang);

     	// Don't make a cylinder for pdf==0
     	if (pdf > 0.)
     	{
			// Colour from blue to red according to deviations from sigma_pdf
			RFLOAT colscale = (pdf - pdfmean) / pdfsigma;
			colscale = XMIPP_MIN(colscale, 5.);
			colscale = XMIPP_MAX(colscale, -1.);
			colscale /= 6.;
			colscale += 1./6.; // colscale ranges from 0 (-5 sigma) to 1 (+5 sigma)

			// The length of the cylinder will depend on the pdf_direction
			RFLOAT Rp = R + Rmax_frac * R * pdf / pdfmax;

			Euler_angles2direction(rot_angles[iang], tilt_angles[iang], v);

			if (Aorient != NULL)
			{
				// In multi-body refinement, the rotations are relative to (rot,tilt)=(0,90) to prevent problems with psi-prior!!!
				Matrix2D<RFLOAT> A;
				rotation3DMatrix(90., 'Y', A, false);
				v = (*Aorient).transpose() * A * v;
			}

			Matrix1D<RFLOAT> offsetp(3);
			if (Acom != NULL)
				offsetp = *Acom;
			else
				offsetp.initZeros();


			// Don't include cylinders with zero length, as chimera will complain about that....
			if (ABS((R - Rp) * XX(v)) > 0.01 ||
					ABS((R - Rp) * YY(v)) > 0.01 ||
					ABS((R - Rp) * ZZ(v)) > 0.01)
			{
				// The width of the cylinders will be determined by the sampling:
				fh_bild << ".color " << colscale << " 0 " << 1. - colscale << std::endl;
				fh_bild << ".cylinder "
						<< R  * XX(v) + offset + XX(offsetp) << " "
						<< R  * YY(v) + offset + YY(offsetp) << " "
						<< R  * ZZ(v) + offset + ZZ(offsetp) << " "
						<< Rp * XX(v) + offset + XX(offsetp) << " "
						<< Rp * YY(v) + offset + YY(offsetp) << " "
						<< Rp * ZZ(v) + offset + ZZ(offsetp) << " "
						<< width
						<<"\n";
			}
     	}

    }

    // Close and write file to disc
    fh_bild.close();

}


///////// PRIVATE STUFF

void HealpixSampling::removePointsOutsideLimitedTiltAngles()
{

    if (ABS(limit_tilt) < 90.)
    {
    	std::vector<RFLOAT> pruned_rot_angles;
    	std::vector<RFLOAT> pruned_tilt_angles;
		std::vector<int>               pruned_directions_ipix;
		pruned_rot_angles.clear();
		pruned_tilt_angles.clear();
		pruned_directions_ipix.clear();

		for (long int i = 0; i < tilt_angles.size(); i++)
		{
			RFLOAT tilt = tilt_angles[i];
			// Let tilt angle range from -90 to 90.
			if (tilt > 90.) tilt -= 180.;

			// Keep side views || keep top views
			if ( (limit_tilt > 0. && ABS(tilt) >= ABS(limit_tilt)) || (limit_tilt < 0. && ABS(tilt) <= ABS(limit_tilt)) )
			{
				pruned_rot_angles.push_back(rot_angles[i]);
				pruned_tilt_angles.push_back(tilt_angles[i]);
				pruned_directions_ipix.push_back(directions_ipix[i]);
			}
		}
		rot_angles = pruned_rot_angles;
		tilt_angles = pruned_tilt_angles;
		directions_ipix   = pruned_directions_ipix;
    }

}


// The way symmetry is handled was copied from Xmipp.
// The original disclaimer is copied below
/***************************************************************************
 *
 * Authors:     Roberto Marabini
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

void HealpixSampling::removeSymmetryEquivalentPoints(RFLOAT max_ang)
{
    // Maximum distance
    RFLOAT cos_max_ang = cos(DEG2RAD(max_ang));
    RFLOAT my_dotProduct;
    Matrix1D<RFLOAT>  direction(3), direction1(3);
    std::vector<Matrix1D<RFLOAT> > directions_vector;

    // Calculate all vectors and fill directions_vector
    for (long int i = 0; i < rot_angles.size(); i++)
    {
    	Euler_angles2direction(rot_angles[i], tilt_angles[i], direction);
    	directions_vector.push_back(direction);
    }

    // First call to conventional remove_redundant_points
    removeSymmetryEquivalentPointsGeometric(pgGroup, pgOrder, directions_vector);

#ifdef  DEBUG_SAMPLING
    writeAllOrientationsToBild("orients_sym0.bild", "0 1 0", 0.021);
#endif

	// Only correct the seams (i.e. the borders of the asymmetric units) for small numbers of directions
    // For large numbers, the sampling is very fine and the probability distributions are probably delta functions anyway
    // Large numbers take long times to calculate...
    // Only a small fraction of the points at the border of the AU is thrown away anyway...
    if (rot_angles.size() < 4000)
    {
    	// Create no_redundant vectors
		std::vector <Matrix1D<RFLOAT> > no_redundant_directions_vector;
		std::vector <RFLOAT> no_redundant_rot_angles;
		std::vector <RFLOAT> no_redundant_tilt_angles;
		std::vector <int> no_redundant_directions_ipix;

		// Then check all points versus each other
		for (long int i = 0; i < rot_angles.size(); i++)
		{

			direction1=directions_vector[i];
			bool uniq = true;

			//for (long int k = 0; k < no_redundant_directions_vector.size(); k++)
			// i is probably closer to latest additions: loop backwards over k....
			for (long int k = no_redundant_directions_vector.size() -1; k >= 0; k--)
			{
				for (int j = 0; j < R_repository.size(); j++)
				{
					direction =  L_repository[j] *
						(no_redundant_directions_vector[k].transpose() *
						 R_repository[j]).transpose();
					//Calculate distance
					my_dotProduct = dotProduct(direction,direction1);
					if (my_dotProduct > cos_max_ang)
					{
						uniq = false;
						break;
					}
				}// for j
				if (!uniq) break;
			} // for k

			if (uniq)
			{
				no_redundant_directions_vector.push_back(directions_vector[i]);
				no_redundant_rot_angles.push_back(rot_angles[i]);
				no_redundant_tilt_angles.push_back(tilt_angles[i]);
				no_redundant_directions_ipix.push_back(directions_ipix[i]);
			}
		} // for i

		// Now overwrite the rot/tilt_angles and directions_vectors with their no_redundant counterparts
		rot_angles = no_redundant_rot_angles;
		tilt_angles = no_redundant_tilt_angles;
		directions_ipix = no_redundant_directions_ipix;
    }
}

void HealpixSampling::removeSymmetryEquivalentPointsGeometric(const int symmetry,
        int sym_order, std::vector <Matrix1D<RFLOAT> >  &directions_vector)
{
    Matrix2D<RFLOAT>  L(4, 4), R(4, 4);
    Matrix2D<RFLOAT>  aux(3, 3);
    Matrix1D<RFLOAT>  row1(3), row2(3), row(3);

    std::vector <Matrix1D<RFLOAT> > no_redundant_directions_vector;
    std::vector <RFLOAT> no_redundant_rot_angles;
    std::vector <RFLOAT> no_redundant_tilt_angles;
    std::vector <int> no_redundant_directions_ipix;

    RFLOAT my_dotProduct;
    if (symmetry == pg_CN)
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >= (-180. / sym_order) &&
                rot_angles[i] <= (180. / sym_order))
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry == pg_CI  ||
             symmetry == pg_CS )
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (tilt_angles[i] <= 90)
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_CNV )
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=    0. / sym_order &&
                rot_angles[i] <=  180. / sym_order)
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_CNH )
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >= -180. / sym_order &&
                rot_angles[i] <=  180. / sym_order &&
                tilt_angles[i] <=    90.
               )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_SN )
    {//OK
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >= -180.*2. / sym_order &&
                rot_angles[i] <=  180.*2. / sym_order &&
                tilt_angles[i] <=    90.
               )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_DN )
    {
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (sym_order == 1)
            {
            	// D1 is special!
            	if (tilt_angles[i] <=    90.)
				{
					no_redundant_rot_angles.push_back(rot_angles[i]);
					no_redundant_tilt_angles.push_back(tilt_angles[i]);
					no_redundant_directions_vector.push_back(directions_vector[i]);
					no_redundant_directions_ipix.push_back(directions_ipix[i]);
				}
            }
            else
            {

				if (rot_angles[i] >= -180. / (sym_order) + 90. &&
					rot_angles[i] <=  180. / (sym_order) + 90. &&
					tilt_angles[i] <=    90.
				   )
				{
					no_redundant_rot_angles.push_back(rot_angles[i]);
					no_redundant_tilt_angles.push_back(tilt_angles[i]);
					no_redundant_directions_vector.push_back(directions_vector[i]);
					no_redundant_directions_ipix.push_back(directions_ipix[i]);
				}
            }
        }// for i
    }
    else if (symmetry  == pg_DNV )
    {
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=   90.  &&
                rot_angles[i] <=  180. / (sym_order) + 90. &&
                tilt_angles[i] <=    90.
               )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_DNH )
    {
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=   90. &&
                rot_angles[i] <=  180. / (sym_order) + 90. &&
                tilt_angles[i] <=   90.
               )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_T )
    {//OK
        Matrix1D<RFLOAT>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3(-0.942809, 0., 0.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_2_by_3_fold_axis_3(3);
        _3_fold_axis_2_by_3_fold_axis_3 = vectorR3(0.471405, 0.272165, 0.7698);
        _3_fold_axis_2_by_3_fold_axis_3.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_3_by_3_fold_axis_1(3);
        _3_fold_axis_3_by_3_fold_axis_1 = vectorR3(0.471404, 0.816497, 0.);
        _3_fold_axis_3_by_3_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=     90. &&
                rot_angles[i] <=   150. ||
                rot_angles[i] ==     0
               )
                if (
                    dotProduct(directions_vector[i], _3_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_2_by_3_fold_axis_3) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_3_by_3_fold_axis_1) >= 0
                )
                {
                	no_redundant_rot_angles.push_back(rot_angles[i]);
                	no_redundant_tilt_angles.push_back(tilt_angles[i]);
                    no_redundant_directions_vector.push_back(directions_vector[i]);
                    no_redundant_directions_ipix.push_back(directions_ipix[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_TD )
    {//OK
        Matrix1D<RFLOAT>  _2_fold_axis_1_by_3_fold_axis_2(3);
        _2_fold_axis_1_by_3_fold_axis_2 = vectorR3(-0.942809, 0., 0.);
        _2_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_2_by_3_fold_axis_5(3);
        _3_fold_axis_2_by_3_fold_axis_5 = vectorR3(0.471405, 0.272165, 0.7698);
        _3_fold_axis_2_by_3_fold_axis_5.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_5_by_2_fold_axis_1(3);
        _3_fold_axis_5_by_2_fold_axis_1 = vectorR3(0., 0.471405, -0.666667);
        _3_fold_axis_5_by_2_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
//           if ( rot_angles[i]>=     120. &&
//                 rot_angles[i]<=   150. ||
//                 rot_angles[i]==     0
//              )
            if (
                dotProduct(directions_vector[i], _2_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_2_by_3_fold_axis_5) >= 0 &&
                dotProduct(directions_vector[i], _3_fold_axis_5_by_2_fold_axis_1) >= 0
            )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_TH )
    {//OK
        Matrix1D<RFLOAT>  _3_fold_axis_1_by_2_fold_axis_1(3);
        _3_fold_axis_1_by_2_fold_axis_1 = vectorR3(-0.816496, 0., 0.);
        _3_fold_axis_1_by_2_fold_axis_1.selfNormalize();
        Matrix1D<RFLOAT>  _2_fold_axis_1_by_2_fold_axis_2(3);
        _2_fold_axis_1_by_2_fold_axis_2 = vectorR3(0.707107, 0.408248, -0.57735);
        _2_fold_axis_1_by_2_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _2_fold_axis_2_by_3_fold_axis_1(3);
        _2_fold_axis_2_by_3_fold_axis_1 = vectorR3(-0.408248, -0.707107, 0.);
        _2_fold_axis_2_by_3_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
//           if ( rot_angles[i]>=     120. &&
//                 rot_angles[i]<=   150. ||
//                 rot_angles[i]==     0
//              )
            if (
                dotProduct(directions_vector[i], _3_fold_axis_1_by_2_fold_axis_1) >= 0 &&
                dotProduct(directions_vector[i], _2_fold_axis_1_by_2_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i], _2_fold_axis_2_by_3_fold_axis_1) >= 0
            )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_O )
    {//OK
        Matrix1D<RFLOAT>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3(0., -1., 1.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_2_by_4_fold_axis(3);
        _3_fold_axis_2_by_4_fold_axis = vectorR3(1., 1., 0.);
        _3_fold_axis_2_by_4_fold_axis.selfNormalize();
        Matrix1D<RFLOAT>  _4_fold_axis_by_3_fold_axis_1(3);
        _4_fold_axis_by_3_fold_axis_1 = vectorR3(-1., 1., 0.);
        _4_fold_axis_by_3_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if ((rot_angles[i] >=   45. &&
                 rot_angles[i] <=  135. &&
                 tilt_angles[i] <=  90.) ||
                rot_angles[i] ==  0.
               )
                if (
                    dotProduct(directions_vector[i], _3_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_2_by_4_fold_axis) >= 0 &&
                    dotProduct(directions_vector[i], _4_fold_axis_by_3_fold_axis_1) >= 0
                )
                {
                	no_redundant_rot_angles.push_back(rot_angles[i]);
                	no_redundant_tilt_angles.push_back(tilt_angles[i]);
                    no_redundant_directions_vector.push_back(directions_vector[i]);
                    no_redundant_directions_ipix.push_back(directions_ipix[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_OH )
    {//OK
        Matrix1D<RFLOAT>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3(0., -1., 1.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_2_by_4_fold_axis(3);
        _3_fold_axis_2_by_4_fold_axis = vectorR3(1., 1., 0.);
        _3_fold_axis_2_by_4_fold_axis.selfNormalize();
        Matrix1D<RFLOAT>  _4_fold_axis_by_3_fold_axis_1(3);
        _4_fold_axis_by_3_fold_axis_1 = vectorR3(-1., 1., 0.);
        _4_fold_axis_by_3_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (rot_angles[i] >=   90. &&
                rot_angles[i] <=  135. &&
                tilt_angles[i] <=  90.)
                if (
                    dotProduct(directions_vector[i], _3_fold_axis_1_by_3_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_2_by_4_fold_axis) >= 0 &&
                    dotProduct(directions_vector[i], _4_fold_axis_by_3_fold_axis_1) >= 0
                )
                {
                	no_redundant_rot_angles.push_back(rot_angles[i]);
                	no_redundant_tilt_angles.push_back(tilt_angles[i]);
                    no_redundant_directions_vector.push_back(directions_vector[i]);
                    no_redundant_directions_ipix.push_back(directions_ipix[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_I || symmetry  == pg_I2)
    {//OK
        Matrix1D<RFLOAT>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = vectorR3(-0.4999999839058737,
                                                 -0.8090170074556163,
                                                  0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = vectorR3(0.4999999839058737,
                                                -0.8090170074556163,
                                                 0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                    dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
               )
                {
					no_redundant_rot_angles.push_back(rot_angles[i]);
					no_redundant_tilt_angles.push_back(tilt_angles[i]);
                    no_redundant_directions_vector.push_back(directions_vector[i]);
                    no_redundant_directions_ipix.push_back(directions_ipix[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_I1)
    {//OK
        Matrix2D<RFLOAT>  A(3, 3);
	    Euler_angles2matrix(0, 90, 0, A);
        Matrix1D<RFLOAT>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(-0.4999999839058737,
                                                 -0.8090170074556163,
                                                  0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.4999999839058737,
                                                -0.8090170074556163,
                                                 0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();

        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                    dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
            )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I3)
    {//OK
        Matrix2D<RFLOAT>  A(3, 3);
	    Euler_angles2matrix(0,31.7174745559,0, A);
        Matrix1D<RFLOAT>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(-0.4999999839058737,
                                                 -0.8090170074556163,
                                                  0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.4999999839058737,
                                                -0.8090170074556163,
                                                 0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();

        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                    dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_by_5_fold_axis_1) >= 0
            )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I4)
    {//OK
        Matrix2D<RFLOAT>  A(3, 3);
	    Euler_angles2matrix(0,-31.7174745559,0, A);
        Matrix1D<RFLOAT>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 0., 1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(0.187592467856686,
                                        -0.303530987314591,
                                        -0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.187592467856686,
                                        0.303530987314591,
                                        -0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();

        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i],
		           _5_fold_axis_2_by_3_fold_axis)   <= 0 &&
                dotProduct(directions_vector[i],
		           _3_fold_axis_by_5_fold_axis_1)   <= 0 &&
		dotProduct(directions_vector[i],
		           _5_fold_axis_1_by_5_fold_axis_2) <= 0
            )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I5)
    {//OK
        std::cerr << "ERROR: Symmetry pg_I5 not implemented" << std::endl;
        exit(0);
    }
    else if (symmetry  == pg_IH || symmetry  == pg_I2H)
    {//OK
        Matrix1D<RFLOAT>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = vectorR3(-0.4999999839058737,
                                                 -0.8090170074556163,
                                                  0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = vectorR3(0.4999999839058737,
                                                -0.8090170074556163,
                                                 0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis =  vectorR3(1.,0.,0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                    dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_by_2_fold_axis) >= 0
               )
                {
					no_redundant_rot_angles.push_back(rot_angles[i]);
					no_redundant_tilt_angles.push_back(tilt_angles[i]);
                    no_redundant_directions_vector.push_back(directions_vector[i]);
                    no_redundant_directions_ipix.push_back(directions_ipix[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_I1H)
    {//OK
        Matrix2D<RFLOAT>  A(3, 3);
	    Euler_angles2matrix(0, 90, 0, A);
        Matrix1D<RFLOAT>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(-0.4999999839058737,
                                                 -0.8090170074556163,
                                                  0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.4999999839058737,
                                                -0.8090170074556163,
                                                 0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis =  A * vectorR3(1.,0.,0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                    dotProduct(directions_vector[i], _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                    dotProduct(directions_vector[i], _5_fold_axis_2_by_3_fold_axis) >= 0 &&
                    dotProduct(directions_vector[i], _3_fold_axis_by_2_fold_axis) >= 0
               )
                {
					no_redundant_rot_angles.push_back(rot_angles[i]);
					no_redundant_tilt_angles.push_back(tilt_angles[i]);
                    no_redundant_directions_vector.push_back(directions_vector[i]);
                    no_redundant_directions_ipix.push_back(directions_ipix[i]);
                }
        }// for i
    }
    else if (symmetry  == pg_I3H)
    {//OK
        Matrix2D<RFLOAT>  A(3, 3);
	    Euler_angles2matrix(0,31.7174745559,0, A);
        Matrix1D<RFLOAT>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 0., 1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(0.187592467856686,
                                        -0.303530987314591,
                                        -0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.187592467856686,
                                        0.303530987314591,
                                        -0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis = vectorR3(0.,1.,0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i],
		           _5_fold_axis_2_by_3_fold_axis)   >= 0 &&
                dotProduct(directions_vector[i],
		           _3_fold_axis_by_5_fold_axis_1)   >= 0 &&
		        dotProduct(directions_vector[i],
		           _5_fold_axis_1_by_5_fold_axis_2) >= 0 &&
                dotProduct(directions_vector[i],
                   _3_fold_axis_by_2_fold_axis)     >= 0
            )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I4H)
    {//OK
        Matrix2D<RFLOAT>  A(3, 3);
	Euler_angles2matrix(0,-31.7174745559,0, A);
        Matrix1D<RFLOAT>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 0., 1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<RFLOAT>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(0.187592467856686,
                                        -0.303530987314591,
                                        -0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.187592467856686,
                                        0.303530987314591,
                                        -0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<RFLOAT>  _3_fold_axis_by_2_fold_axis(3);
        _3_fold_axis_by_2_fold_axis = vectorR3(0.,1.,0.);
        _3_fold_axis_by_2_fold_axis.selfNormalize();
        for (long int i = 0; i < rot_angles.size(); i++)
        {
            if (
                dotProduct(directions_vector[i],
		           _5_fold_axis_2_by_3_fold_axis)   <= 0 &&
                dotProduct(directions_vector[i],
		           _3_fold_axis_by_5_fold_axis_1)   <= 0 &&
		        dotProduct(directions_vector[i],
		           _5_fold_axis_1_by_5_fold_axis_2) <= 0 &&
                dotProduct(directions_vector[i],
                   _3_fold_axis_by_2_fold_axis)     >= 0
            )
            {
                no_redundant_rot_angles.push_back(rot_angles[i]);
                no_redundant_tilt_angles.push_back(tilt_angles[i]);
                no_redundant_directions_vector.push_back(directions_vector[i]);
                no_redundant_directions_ipix.push_back(directions_ipix[i]);
            }
        }// for i
    }
    else if (symmetry  == pg_I5H)
    {//OK
        std::cerr << "ERROR: pg_I5H Symmetry not implemented" << std::endl;
        exit(0);
    }
    else
    {
        std::cerr << "ERROR: Symmetry " << symmetry  << "is not known" << std::endl;
        exit(0);
    }


    // Now overwrite the rot/tilt_angles and directions_vectors with their no_redundant counterparts
    rot_angles = no_redundant_rot_angles;
    tilt_angles = no_redundant_tilt_angles;
    directions_vector = no_redundant_directions_vector;
    directions_ipix = no_redundant_directions_ipix;


}


#undef DEBUG_SAMPLING
