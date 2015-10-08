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
void HealpixSampling::clear()
{
	is_3D = false;
	fn_sym = "C1";
	limit_tilt = psi_step = offset_range = offset_step = 0.;
	orientational_prior_mode = NOPRIOR;
	directions_ipix.clear();
	rot_angles.clear();
	tilt_angles.clear();
	psi_angles.clear();
	translations_x.clear();
	translations_y.clear();
	translations_z.clear();
	L_repository.clear();
	R_repository.clear();
	pgGroup = pgOrder = 0;

}

void HealpixSampling::initialise(int prior_mode, int ref_dim, bool _do_3d_trans, bool do_warnpsi)
{

	// Set the prior mode (belongs to mlmodel, but very useful inside this object)
	orientational_prior_mode = prior_mode;

	if (ref_dim != -1)
		is_3D = (ref_dim == 3);

	// Set flag for x,y,z-translations
	is_3d_trans = _do_3d_trans;

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
		SymList SL;
		SL.isSymmetryGroup(fn_sym, pgGroup, pgOrder);
		SL.read_sym_file(fn_sym);

		// Precalculate (3x3) symmetry matrices
		Matrix2D<double>  L(4, 4), R(4, 4);
		Matrix2D<double>  Identity(3,3);
		Identity.initIdentity();
		R_repository.clear();
		L_repository.clear();
		R_repository.push_back(Identity);
		L_repository.push_back(Identity);
		for (int isym = 0; isym < SL.SymsNo(); isym++)
		{
			SL.get_matrices(isym, L, R);
			R.resize(3, 3);
			L.resize(3, 3);
			R_repository.push_back(R);
			L_repository.push_back(L);
		}
	}
	else
	{
		int t_nr_psi = CEIL(360./psi_step);
		if(t_nr_psi%32!=0 && do_warnpsi)
		{
			std::cerr << std::endl << "WARNING:" << std::endl << "the chosen --psi_step (" << psi_step << " degrees) gives " << t_nr_psi << " steps (!=0 mod 32), which is not optimal when you are using GPUs" << std::endl;
			std::cerr << " May we suggest using --psi_step 11.25 (32 steps), 5.625 (64 steps), 3.75 (96 steps), etc... ? " << std::endl << "(UNLESS YOU ARE DOING LOCAL SEARCHES / AUTO-REFINE)" << std::endl << std::endl;

//			Force-adjust psi_step:
//			t_nr_psi = CEIL((float)t_nr_psi / 32.0)*32;
//			float tdiff = psi_step - 360./(double)t_nr_psi;
		}
		psi_step = 360./(double)t_nr_psi;
		fn_sym = "C1"; // This may not be set yet if restarting a 2D run....
	}


	// Store the not-oversampled translations, and make sure oversampled sampling is 1 pixel
	setTranslations();

	// Store the non-oversampled projection directions
	setOrientations();

	// Random perturbation and filling of the directions, psi_angles and translations vectors
	resetRandomlyPerturbedSampling();

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
	MD.setValue(EMDL_SAMPLING_PERTURB, random_perturbation);
	MD.setValue(EMDL_SAMPLING_PERTURBATION_FACTOR, perturbation_factor);
	MD.write(fh);

	// In the 3D case, also write a table with the sampled rot, tilt angles
	if (is_3D)
	{
		MD.clear();
		MD.setIsList(false);
		MD.setName("sampling_directions");
		for (long int idir = 0; idir < NrDirections(); idir++)
		{
			double rot, tilt;
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

void HealpixSampling::setTranslations(double _offset_step, double _offset_range)
{
	translations_x.clear();
	translations_y.clear();
	translations_z.clear();
	if (_offset_step > 0. && _offset_range >= 0.)
	{
		offset_step = _offset_step;
		offset_range = _offset_range;
	}

	int maxr = CEIL(offset_range / offset_step);
	for (long int ix = -maxr; ix <= maxr; ix++)
	{
		double xoff = ix * offset_step;
		for (long int iy = -maxr; iy <= maxr; iy++)
		{
			double yoff = iy * offset_step;

			if (is_3d_trans)
			{
				for (long int iz = -maxr; iz <= maxr; iz++)
				{
					double zoff = iz * offset_step;
					if (xoff*xoff + yoff*yoff + zoff*zoff <= offset_range * offset_range)
					{
						translations_x.push_back(xoff);
						translations_y.push_back(yoff);
						translations_z.push_back(zoff);

					}
				}
			}
			else
			{
				if (xoff*xoff + yoff*yoff <= offset_range * offset_range)
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

}
/* Set only a single translation */
void HealpixSampling::addOneTranslation(double offset_x, double offset_y, double offset_z, bool do_clear)
{
	if (do_clear)
	{
		translations_x.clear();
		translations_y.clear();
		translations_z.clear();
	}
	translations_x.push_back(offset_x);
	translations_y.push_back(offset_y);
	if (is_3d_trans)
		translations_z.push_back(offset_z);
}



void HealpixSampling::setOrientations(int _order, double _psi_step)
{

	// Initialise
	directions_ipix.clear();
	rot_angles.clear();
	tilt_angles.clear();
	psi_angles.clear();

	// Setup the HealPix object
	// For adaptive oversampling only precalculate the COARSE sampling!
	if (_order >= 0)
	{
		healpix_base.Set(_order, NEST);
		healpix_order = _order;
	}

	// 3D directions
	if (is_3D)
	{
		double rot, tilt;
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
		writeAllOrientationsToBild("orients_all.bild", "1 0 0", 0.020);
#endif
		// Now remove symmetry-related pixels
		// TODO check size of healpix_base.max_pixrad
		removeSymmetryEquivalentPoints(0.5 * RAD2DEG(healpix_base.max_pixrad()));

#ifdef  DEBUG_SAMPLING
		writeAllOrientationsToBild("orients_sym.bild", "0 1 0", 0.021);
#endif

		// Also remove limited tilt angles
		removePointsOutsideLimitedTiltAngles();

		#ifdef  DEBUG_SAMPLING
		if (ABS(limit_tilt) < 90.)
			writeAllOrientationsToBild("orients_tilt.bild", "1 1 0", 0.022);
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
	double psi;
	psi_step = 360./(double)nr_psi;
	for (int ipsi = 0; ipsi < nr_psi; ipsi++)
	{
		psi = ipsi * psi_step;
		psi_angles.push_back(psi);
	}
}

/* Set only a single orientation */
void HealpixSampling::addOneOrientation(double rot, double tilt, double psi, bool do_clear)
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


void HealpixSampling::writeAllOrientationsToBild(FileName fn_bild, std::string rgb, double size)
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


    Matrix1D<double> v(3);
	out << ".color " << rgb << std::endl;

	for (unsigned long int ipix = 0; ipix < rot_angles.size(); ipix++)
	{
		Euler_angles2direction(rot_angles[ipix], tilt_angles[ipix], v);
		out <<  ".sphere " << XX(v) << " " << YY(v) << " " << ZZ(v) <<  floatToString(size) << std::endl;
	}

	out.close();

}

void HealpixSampling::writeNonZeroPriorOrientationsToBild(FileName fn_bild, double rot_prior, double tilt_prior,
		std::vector<int> &pointer_dir_nonzeroprior, std::string rgb, double size)
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

    Matrix1D<double> v(3);

	Euler_angles2direction(rot_prior, tilt_prior, v);
	out << ".color 0 0 0 \n";
	out <<  ".sphere " << XX(v) << " " << YY(v) << " " << ZZ(v) <<  floatToString(size) << std::endl;

	out << ".color " << rgb << std::endl;
	for (unsigned long int ipix = 0; ipix < pointer_dir_nonzeroprior.size(); ipix++)
	{
		long int idir = pointer_dir_nonzeroprior[ipix];
		Euler_angles2direction(rot_angles[idir], tilt_angles[idir], v);
		out <<  ".sphere " << XX(v) << " " << YY(v) << " " << ZZ(v) <<  floatToString(size) << std::endl;
	}

	out.close();

}

void HealpixSampling::selectOrientationsWithNonZeroPriorProbability(
		double prior_rot, double prior_tilt, double prior_psi,
		double sigma_rot, double sigma_tilt, double sigma_psi,
    	std::vector<int> &pointer_dir_nonzeroprior, std::vector<double> &directions_prior,
    	std::vector<int> &pointer_psi_nonzeroprior, std::vector<double> &psi_prior,
		double sigma_cutoff)
{
	pointer_dir_nonzeroprior.clear();
	directions_prior.clear();

	if (is_3D)
	{

		// Loop over all directions
		double sumprior = 0.;
		// Keep track of the closest distance to prevent 0 orientations
		double best_ang = 9999.;
		long int best_idir = -999;
		for (long int idir = 0; idir < rot_angles.size(); idir++)
		{

			// Any prior involving rot and/or tilt.
			if (sigma_rot > 0. || sigma_tilt > 0. )
			{

				Matrix1D<double> prior_direction, my_direction, sym_direction, best_direction;
				// Get the direction of the prior
				Euler_angles2direction(prior_rot, prior_tilt, prior_direction);

				// Get the current direction in the loop
				Euler_angles2direction(rot_angles[idir], tilt_angles[idir], my_direction);

				// Loop over all symmetry operators to find the operator that brings this direction nearest to the prior
				double best_dotProduct = dotProduct(prior_direction, my_direction);
				best_direction = my_direction;
				for (int j = 0; j < R_repository.size(); j++)
				{
					sym_direction =  L_repository[j] * (my_direction.transpose() * R_repository[j]).transpose();
					double my_dotProduct = dotProduct(prior_direction, sym_direction);
					if (my_dotProduct > best_dotProduct)
					{
						best_direction = sym_direction;
						best_dotProduct = my_dotProduct;
					}
				}

				if (sigma_rot > 0. && sigma_tilt > 0.)
				{

					double diffang = ACOSD( dotProduct(best_direction, prior_direction) );
					if (diffang > 180.) diffang = ABS(diffang - 360.);

					// Only consider differences within sigma_cutoff * sigma_rot
					if (diffang < sigma_cutoff * sigma_rot)
					{
						// TODO!!! If tilt is zero then any rot will be OK!!!!!
						double prior = gaussian1D(diffang, sigma_rot, 0.);
						pointer_dir_nonzeroprior.push_back(idir);
						directions_prior.push_back(prior);
						sumprior += prior;
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
					double best_rot, best_tilt;

					Euler_direction2angles(best_direction, best_rot, best_tilt);
					double diffrot = ABS(best_rot - prior_rot);
					if (diffrot > 180.) diffrot = ABS(diffrot - 360.);

					// Only consider differences within sigma_cutoff * sigma_rot
					if (diffrot < sigma_cutoff * sigma_rot)
					{
						double prior = gaussian1D(diffrot, sigma_rot, 0.);
						pointer_dir_nonzeroprior.push_back(idir);
						directions_prior.push_back(prior);
						sumprior += prior;
					}

					// Keep track of the nearest direction
					if (diffrot < best_ang)
					{
						best_idir = idir;
						best_ang = diffrot;
					}

				}
				else if (sigma_tilt > 0.)
				{

					double best_rot, best_tilt;

					Euler_direction2angles(best_direction, best_rot, best_tilt);
					double difftilt = ABS(best_tilt - prior_tilt);
					if (difftilt > 180.) difftilt = ABS(difftilt - 360.);

					// Only consider differences within sigma_cutoff * sigma_tilt
					if (difftilt < sigma_cutoff * sigma_tilt)
					{
						double prior = gaussian1D(difftilt, sigma_tilt, 0.);
						pointer_dir_nonzeroprior.push_back(idir);
						directions_prior.push_back(prior);
						sumprior += prior;
					}

					// Keep track of the nearest direction
					if (difftilt < best_ang)
					{
						best_idir = idir;
						best_ang = difftilt;
					}

				}

			} // end if any prior involving rot and/or tilt
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

	double sumprior = 0.;
	double best_diff = 9999.;
	long int best_ipsi = -999;
	for (long int ipsi = 0; ipsi < psi_angles.size(); ipsi++)
	{
		if (sigma_psi > 0.)
		{
			double diffpsi = ABS(psi_angles[ipsi] - prior_psi);
			if (diffpsi > 180.) diffpsi = ABS(diffpsi - 360.);

			// Only consider differences within sigma_cutoff * sigma_psi
			if (diffpsi < sigma_cutoff * sigma_psi)
			{
				double prior = gaussian1D(diffpsi, sigma_psi, 0.);
				pointer_psi_nonzeroprior.push_back(ipsi);
				psi_prior.push_back(prior);
				sumprior += prior;

				// TMP DEBUGGING
				if (prior == 0.)
				{
					std::cerr << " psi_angles[ipsi]= " << psi_angles[ipsi] << " prior_psi= " << prior_psi << " orientational_prior_mode= " << orientational_prior_mode << std::endl;
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

void HealpixSampling::checkDirection(double &rot, double &tilt)
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

void HealpixSampling::getDirectionFromHealPix(long int ipix, double &rot, double &tilt)
{
	double zz, phi;
	healpix_base.pix2ang_z_phi(ipix, zz, phi);
	rot = RAD2DEG(phi);
	tilt = ACOSD(zz);

	// The geometrical considerations about the symmetry below require that rot = [-180,180] and tilt [0,180]
	checkDirection(rot, tilt);

}

double HealpixSampling::getTranslationalSampling(int adaptive_oversampling)
{
	return offset_step / std::pow(2., adaptive_oversampling);
}

double HealpixSampling::getAngularSampling(int adaptive_oversampling)
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


void HealpixSampling::getDirection(long int idir, double &rot, double &tilt)
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

void HealpixSampling::getPsiAngle(long int ipsi, double &psi)
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

void HealpixSampling::getTranslation(long int itrans, double &trans_x, double &trans_y, double &trans_z)
{
#ifdef DEBUG_CHECKSIZES
if (itrans >= translations_x.size())
{
	std::cerr<< "itrans= "<<itrans<<" translations_x.size()= "<< translations_x.size() <<std::endl;
	REPORT_ERROR("itrans >= translations_x.size()");
}
#endif
	trans_x = translations_x[itrans];
	trans_y = translations_y[itrans];
	if (is_3d_trans)
		trans_z = translations_z[itrans];
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

void HealpixSampling::getTranslations(long int itrans, int oversampling_order,
		std::vector<double > &my_translations_x,
		std::vector<double > &my_translations_y,
		std::vector<double > &my_translations_z)
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
		my_translations_x.push_back(translations_x[itrans]);
		my_translations_y.push_back(translations_y[itrans]);
		if (is_3d_trans)
			my_translations_z.push_back(translations_z[itrans]);
	}
	else
	{
		int nr_oversamples = ROUND(std::pow(2., oversampling_order));
		if (is_3d_trans)
		{
			for (int itrans_overy = 0; itrans_overy < nr_oversamples; itrans_overy++)
			{
				double over_yoff = translations_y[itrans] - 0.5 * offset_step + (0.5 + itrans_overy) * offset_step / nr_oversamples;
				for (int itrans_overx = 0; itrans_overx < nr_oversamples; itrans_overx++)
				{
					double over_xoff = translations_x[itrans] - 0.5 * offset_step + (0.5 + itrans_overx) * offset_step / nr_oversamples;
					for (int itrans_overz = 0; itrans_overz < nr_oversamples; itrans_overz++)
					{
						double over_zoff = translations_z[itrans] - 0.5 * offset_step + (0.5 + itrans_overz) * offset_step / nr_oversamples;
						my_translations_x.push_back(over_xoff);
						my_translations_y.push_back(over_yoff);
						my_translations_z.push_back(over_zoff);
					}
				}
			}
		}
		else
		{
			for (int itrans_overy = 0; itrans_overy < nr_oversamples; itrans_overy++)
			{
				double over_yoff = translations_y[itrans] - 0.5 * offset_step + (0.5 + itrans_overy) * offset_step / nr_oversamples;
				for (int itrans_overx = 0; itrans_overx < nr_oversamples; itrans_overx++)
				{
					double over_xoff = translations_x[itrans] - 0.5 * offset_step + (0.5 + itrans_overx) * offset_step / nr_oversamples;
					my_translations_x.push_back(over_xoff);
					my_translations_y.push_back(over_yoff);
				}
			}
		}
	}

	if (ABS(random_perturbation) > 0.)
	{
		double myperturb = random_perturbation * offset_step;
		for (int iover = 0; iover < my_translations_x.size(); iover++)
		{
			my_translations_x[iover] += myperturb;
			my_translations_y[iover] += myperturb;
			if (is_3d_trans)
				my_translations_z[iover] += myperturb;
		}
	}

}

void HealpixSampling::getOrientations(long int idir, long int ipsi, int oversampling_order,
		std::vector<double > &my_rot, std::vector<double > &my_tilt, std::vector<double > &my_psi,
		std::vector<int> &pointer_dir_nonzeroprior, std::vector<double> &directions_prior,
		std::vector<int> &pointer_psi_nonzeroprior, std::vector<double> &psi_prior)
{
	my_rot.clear();
	my_tilt.clear();
	my_psi.clear();
	long int my_idir, my_ipsi;
	if (orientational_prior_mode == NOPRIOR)
	{
		my_idir = idir;
		my_ipsi = ipsi;
	}
	else
	{
#ifdef DEBUG_CHECKSIZES
	if (idir >= pointer_dir_nonzeroprior.size())
	{
		std::cerr<< "idir= "<<idir<<" pointer_dir_nonzeroprior.size()= "<< pointer_dir_nonzeroprior.size() <<std::endl;
		REPORT_ERROR("idir >= pointer_dir_nonzeroprior.size()");
	}
	if (ipsi >= pointer_psi_nonzeroprior.size())
	{
		std::cerr<< "ipsi= "<<ipsi<<" pointer_psi_nonzeroprior.size()= "<< pointer_psi_nonzeroprior.size() <<std::endl;
		REPORT_ERROR("ipsi >= pointer_psi_nonzeroprior.size()");
	}
#endif
		my_idir = pointer_dir_nonzeroprior[idir];
		my_ipsi = pointer_psi_nonzeroprior[ipsi];
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
		double rot, tilt;
		// Get x, y and face for the original, coarse grid
		long int ipix = directions_ipix[my_idir];
		healpix_base.nest2xyf(ipix, x, y, face);
		// Loop over the oversampled Healpix pixels on the fine grid
		for (int j = fact * y; j < fact * (y+1); ++j)
		{
			for (int i = fact * x; i < fact * (x+1); ++i)
			{
				long int overpix = HealPixOver.xyf2nest(i, j, face);
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
		double myperturb = random_perturbation * getAngularSampling();
		for (int iover = 0; iover < my_rot.size(); iover++)
		{
			if (is_3D)
			{
				Matrix2D<double> A(3,3), R(3,3);
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
		double rot, double tilt, std::vector<double> &oversampled_rot,
		std::vector<double> &oversampled_tilt, std::vector<double> &oversampled_psi)
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
			double overpsi = psi_angles[ipsi] - 0.5 * psi_step + (0.5 + ipsi_over) * psi_step / nr_ipsi_over;
			oversampled_rot.push_back(rot);
			oversampled_tilt.push_back(tilt);
			if (!is_3D && overpsi>180.)
				overpsi-=360.;
			oversampled_psi.push_back(overpsi);
		}
	}

}

/* Calculate an angular distance between two sets of Euler angles */
double HealpixSampling::calculateAngularDistance(double rot1, double tilt1, double psi1,
		double rot2, double tilt2, double psi2)
{
	Matrix1D<double>  direction1(3), direction1p(3), direction2(3);
	Euler_angles2direction(rot1, tilt1, direction1);
	Euler_angles2direction(rot2, tilt2, direction2);

	// Find the symmetry operation where the Distance based on Euler axes is minimal
	double min_axes_dist = 3600.;
	double rot2p, tilt2p, psi2p;
	Matrix2D<double> E1, E2;
	Matrix1D<double> v1, v2;
	for (int j = 0; j < R_repository.size(); j++)
	{

        Euler_apply_transf(L_repository[j], R_repository[j], rot2, tilt2, psi2, rot2p, tilt2p, psi2p);

	    // Distance based on Euler axes
	    Euler_angles2matrix(rot1, tilt1, psi1, E1);
	    Euler_angles2matrix(rot2p, tilt2p, psi2p, E2);
	    double axes_dist = 0;
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

void HealpixSampling::writeBildFileOrientationalDistribution(MultidimArray<double> &pdf_direction,
		FileName &fn_bild, double R, double offset, double Rmax_frac, double width_frac)
{
	if (!is_3D)
		return;

	if (XSIZE(pdf_direction) != rot_angles.size())
		REPORT_ERROR("HealpixSampling::writeBildFileOrientationalDistribution XSIZE(pdf_direction) != rot_angles.size()!");


	double pdfmax, pdfmin, pdfmean, pdfsigma;
	pdf_direction.computeStats(pdfmean, pdfsigma, pdfmin, pdfmax);

	std::ofstream fh_bild;
    fh_bild.open(fn_bild.c_str(), std::ios::out);
    if (!fh_bild)
    	REPORT_ERROR("HealpixSampling::writeBildFileOrientationalDistribution: cannot open " + fn_bild);

    // 2 * PI * R = 360 degrees, 2*radius should cover angular sampling at width_frac=1
    double width = width_frac * PI*R*(getAngularSampling()/360.);
    Matrix1D<double> v(3);

    for (long int iang = 0; iang < rot_angles.size(); iang++)
    {
     	double pdf = DIRECT_A1D_ELEM(pdf_direction, iang);

     	// Don't make a cylinder for pdf==0
     	if (pdf > 0.)
     	{
			// Colour from blue to red according to deviations from sigma_pdf
			double colscale = (pdf - pdfmean) / pdfsigma;
			colscale = XMIPP_MIN(colscale, 5.);
			colscale = XMIPP_MAX(colscale, -1.);
			colscale /= 6.;
			colscale += 1./6.; // colscale ranges from 0 (-5 sigma) to 1 (+5 sigma)

			// The length of the cylinder will depend on the pdf_direction
			double Rp = R + Rmax_frac * R * pdf / pdfmax;

			Euler_angles2direction(rot_angles[iang], tilt_angles[iang], v);

			// Don't include cylinders with zero length, as chimera will complain about that....
			if (ABS((R - Rp) * XX(v)) > 0.01 ||
					ABS((R - Rp) * YY(v)) > 0.01 ||
					ABS((R - Rp) * ZZ(v)) > 0.01)
			{
				// The width of the cylinders will be determined by the sampling:
				fh_bild << ".color " << colscale << " 0 " << 1. - colscale << std::endl;
				fh_bild << ".cylinder "
						<< R  * XX(v) + offset << " "
						<< R  * YY(v) + offset << " "
						<< R  * ZZ(v) + offset << " "
						<< Rp * XX(v) + offset << " "
						<< Rp * YY(v) + offset << " "
						<< Rp * ZZ(v) + offset << " "
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
    	std::vector<double> pruned_rot_angles;
    	std::vector<double> pruned_tilt_angles;
		std::vector<int>               pruned_directions_ipix;
		pruned_rot_angles.clear();
		pruned_tilt_angles.clear();
		pruned_directions_ipix.clear();

		for (long int i = 0; i < tilt_angles.size(); i++)
		{
			double tilt = tilt_angles[i];
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

void HealpixSampling::removeSymmetryEquivalentPoints(double max_ang)
{
    // Maximum distance
    double cos_max_ang = cos(DEG2RAD(max_ang));
    double my_dotProduct;
    Matrix1D<double>  direction(3), direction1(3);
    std::vector<Matrix1D<double> > directions_vector;

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
		std::vector <Matrix1D<double> > no_redundant_directions_vector;
		std::vector <double> no_redundant_rot_angles;
		std::vector <double> no_redundant_tilt_angles;
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
        int sym_order, std::vector <Matrix1D<double> >  &directions_vector)
{
    Matrix2D<double>  L(4, 4), R(4, 4);
    Matrix2D<double>  aux(3, 3);
    Matrix1D<double>  row1(3), row2(3), row(3);

    std::vector <Matrix1D<double> > no_redundant_directions_vector;
    std::vector <double> no_redundant_rot_angles;
    std::vector <double> no_redundant_tilt_angles;
    std::vector <int> no_redundant_directions_ipix;

    double my_dotProduct;
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
        Matrix1D<double>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3(-0.942809, 0., 0.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<double>  _3_fold_axis_2_by_3_fold_axis_3(3);
        _3_fold_axis_2_by_3_fold_axis_3 = vectorR3(0.471405, 0.272165, 0.7698);
        _3_fold_axis_2_by_3_fold_axis_3.selfNormalize();
        Matrix1D<double>  _3_fold_axis_3_by_3_fold_axis_1(3);
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
        Matrix1D<double>  _2_fold_axis_1_by_3_fold_axis_2(3);
        _2_fold_axis_1_by_3_fold_axis_2 = vectorR3(-0.942809, 0., 0.);
        _2_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<double>  _3_fold_axis_2_by_3_fold_axis_5(3);
        _3_fold_axis_2_by_3_fold_axis_5 = vectorR3(0.471405, 0.272165, 0.7698);
        _3_fold_axis_2_by_3_fold_axis_5.selfNormalize();
        Matrix1D<double>  _3_fold_axis_5_by_2_fold_axis_1(3);
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
        Matrix1D<double>  _3_fold_axis_1_by_2_fold_axis_1(3);
        _3_fold_axis_1_by_2_fold_axis_1 = vectorR3(-0.816496, 0., 0.);
        _3_fold_axis_1_by_2_fold_axis_1.selfNormalize();
        Matrix1D<double>  _2_fold_axis_1_by_2_fold_axis_2(3);
        _2_fold_axis_1_by_2_fold_axis_2 = vectorR3(0.707107, 0.408248, -0.57735);
        _2_fold_axis_1_by_2_fold_axis_2.selfNormalize();
        Matrix1D<double>  _2_fold_axis_2_by_3_fold_axis_1(3);
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
        Matrix1D<double>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3(0., -1., 1.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<double>  _3_fold_axis_2_by_4_fold_axis(3);
        _3_fold_axis_2_by_4_fold_axis = vectorR3(1., 1., 0.);
        _3_fold_axis_2_by_4_fold_axis.selfNormalize();
        Matrix1D<double>  _4_fold_axis_by_3_fold_axis_1(3);
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
        Matrix1D<double>  _3_fold_axis_1_by_3_fold_axis_2(3);
        _3_fold_axis_1_by_3_fold_axis_2 = vectorR3(0., -1., 1.);
        _3_fold_axis_1_by_3_fold_axis_2.selfNormalize();
        Matrix1D<double>  _3_fold_axis_2_by_4_fold_axis(3);
        _3_fold_axis_2_by_4_fold_axis = vectorR3(1., 1., 0.);
        _3_fold_axis_2_by_4_fold_axis.selfNormalize();
        Matrix1D<double>  _4_fold_axis_by_3_fold_axis_1(3);
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
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = vectorR3(-0.4999999839058737,
                                                 -0.8090170074556163,
                                                  0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
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
        Matrix2D<double>  A(3, 3);
	    Euler_angles2matrix(0, 90, 0, A);
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(-0.4999999839058737,
                                                 -0.8090170074556163,
                                                  0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
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
        Matrix2D<double>  A(3, 3);
	    Euler_angles2matrix(0,31.7174745559,0, A);
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(-0.4999999839058737,
                                                 -0.8090170074556163,
                                                  0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
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
        Matrix2D<double>  A(3, 3);
	    Euler_angles2matrix(0,-31.7174745559,0, A);
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 0., 1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(0.187592467856686,
                                        -0.303530987314591,
                                        -0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
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
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = vectorR3(-0.4999999839058737,
                                                 -0.8090170074556163,
                                                  0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = vectorR3(0.4999999839058737,
                                                -0.8090170074556163,
                                                 0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_2_fold_axis(3);
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
        Matrix2D<double>  A(3, 3);
	    Euler_angles2matrix(0, 90, 0, A);
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 1., 0.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(-0.4999999839058737,
                                                 -0.8090170074556163,
                                                  0.3090169861701543);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.4999999839058737,
                                                -0.8090170074556163,
                                                 0.3090169861701543);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_2_fold_axis(3);
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
        Matrix2D<double>  A(3, 3);
	    Euler_angles2matrix(0,31.7174745559,0, A);
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 0., 1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(0.187592467856686,
                                        -0.303530987314591,
                                        -0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.187592467856686,
                                        0.303530987314591,
                                        -0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_2_fold_axis(3);
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
        Matrix2D<double>  A(3, 3);
	Euler_angles2matrix(0,-31.7174745559,0, A);
        Matrix1D<double>  _5_fold_axis_1_by_5_fold_axis_2(3);
        _5_fold_axis_1_by_5_fold_axis_2 = A * vectorR3(0., 0., 1.);
        _5_fold_axis_1_by_5_fold_axis_2.selfNormalize();
        Matrix1D<double>  _5_fold_axis_2_by_3_fold_axis(3);
        _5_fold_axis_2_by_3_fold_axis = A * vectorR3(0.187592467856686,
                                        -0.303530987314591,
                                        -0.491123477863004);
        _5_fold_axis_2_by_3_fold_axis.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_5_fold_axis_1(3);
        _3_fold_axis_by_5_fold_axis_1 = A * vectorR3(0.187592467856686,
                                        0.303530987314591,
                                        -0.491123477863004);
        _3_fold_axis_by_5_fold_axis_1.selfNormalize();
        Matrix1D<double>  _3_fold_axis_by_2_fold_axis(3);
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
