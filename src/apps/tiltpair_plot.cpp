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
#include <src/args.h>
#include <src/strings.h>
#include <src/funcs.h>
#include <src/memory.h>
#include <src/euler.h>
#include <src/image.h>
#include <src/projector.h>
#include <src/metadata_table.h>
#include <src/fftw.h>
#include <src/ctf.h>
#include <src/time.h>
#include <src/symmetries.h>



class tiltpair_plot_parameters
{
	public:
	FileName fn_unt, fn_til, fn_eps, fn_sym;
	MetaDataTable MDu, MDt;
	RFLOAT exp_tilt, exp_beta, dist_from_alpha, dist_from_tilt, plot_max_tilt, plot_spot_radius;
	// I/O Parser
	IOParser parser;
	SymList SL;
	std::ofstream fh_eps;


	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");
		fn_unt = parser.getOption("--u", "Input STAR file with untilted particles");
		fn_til = parser.getOption("--t", "Input STAR file with tilted particles");
		fn_eps = parser.getOption("--o", "Output EPS file ", "tiltpair.eps");
		fn_sym = parser.getOption("--sym", "Symmetry point group", "C1");
		exp_tilt = textToFloat(parser.getOption("--exp_tilt", "Choose symmetry operator that gives tilt angle closest to this value", "0."));
		exp_beta = textToFloat(parser.getOption("--exp_beta", "Choose symmetry operator that gives beta angle closest to this value", "0."));
		dist_from_alpha = textToFloat(parser.getOption("--dist_from_alpha", "Direction (alpha angle) of tilt axis from which to calculate distance", "0."));
		dist_from_tilt = textToFloat(parser.getOption("--dist_from_tilt", "Tilt angle from which to calculate distance", "0."));
		plot_max_tilt = textToFloat(parser.getOption("--max_tilt", "Maximum tilt angle to plot in the EPS file", "90."));
		plot_spot_radius = textToInteger(parser.getOption("--spot_radius", "Radius in pixels of the spots in the tiltpair plot", "3"));

		// Check for errors in the command-line option
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line, exiting...");
	}

	void initialise()
	{

		// Get the MDs for both untilted and tilted particles
		MDu.read(fn_unt);
		MDt.read(fn_til);
		if (MDu.numberOfObjects() != MDt.numberOfObjects())
			REPORT_ERROR("Tiltpair plot ERROR: untilted and tilted STAR files have unequal number of entries.");

		// Get the symmetry point group
		int pgGroup, pgOrder;
		SL.isSymmetryGroup(fn_sym, pgGroup, pgOrder);
		SL.read_sym_file(fn_sym);

		// Make postscript header
		fh_eps.open(fn_eps.c_str(), std::ios::out);
		if (!fh_eps)
			REPORT_ERROR("Tiltpair plot ERROR: Cannot open " + fn_eps + " for output");

		fh_eps << "%%!PS-Adobe-2.0\n";
		fh_eps << "%% Creator: Tilt pair analysis \n";
		fh_eps << "%% Pages: 1\n";
		fh_eps << "0 setgray\n";
		fh_eps << "0.1 setlinewidth\n";
		// Draw circles on postscript: 250pixels=plot_max_tilt
		fh_eps << "300 400 83 0 360 arc closepath stroke\n";
		fh_eps << "300 400 167 0 360 arc closepath stroke\n";
		fh_eps << "300 400 250 0 360 arc closepath stroke\n";
		fh_eps << "300 150 newpath moveto 300 650 lineto stroke\n";
		fh_eps << "50 400 newpath moveto 550 400 lineto stroke\n";
	}

	void add_to_postscript(RFLOAT tilt_angle, RFLOAT alpha, RFLOAT beta)
	{

		RFLOAT rr, th, x, y, r, g, b;

		rr = (tilt_angle / plot_max_tilt)* 250;
		x = 300. + rr * COSD(alpha);
		y = 400. + rr * SIND(alpha);
		value_to_redblue_scale(ABS(90.-beta), 0., 90., r, g, b);
		fh_eps << x << " " << y << " " << plot_spot_radius << " 0 360 arc closepath "<<r<<" "<<g<<" "<<b<<" setrgbcolor fill stroke\n";
	}

	void value_to_redblue_scale(RFLOAT val, RFLOAT minF, RFLOAT maxF, RFLOAT &r, RFLOAT &g, RFLOAT &b)
	{
		RFLOAT diff, half;
		half = (maxF - minF)/2.;
		if (val < half)
		{
			r=val/half;
			b=1.;
		}
		else
		{
			b=(maxF-val)/half;
			r=1.;
		}
		g=0.;

	}

	RFLOAT check_symmetries(RFLOAT rot1, RFLOAT tilt1, RFLOAT psi1,
			RFLOAT &rot2, RFLOAT &tilt2, RFLOAT &psi2)
	{

		int imax = SL.SymsNo() + 1;
		Matrix2D<RFLOAT>  L(4, 4), R(4, 4);  // A matrix from the list
		RFLOAT best_ang_dist = 3600;
		RFLOAT best_rot2, best_tilt2, best_psi2;
		RFLOAT tilt_angle, alpha, beta;

		for (int i = 0; i < imax; i++)
		{
			RFLOAT rot2p, tilt2p, psi2p;
			if (i == 0)
			{
				rot2p = rot2;
				tilt2p = tilt2;
				psi2p = psi2;
			}
			else
			{
				SL.get_matrices(i - 1, L, R);
				L.resize(3, 3); // Erase last row and column
				R.resize(3, 3); // as only the relative orientation is useful and not the translation
				Euler_apply_transf(L, R, rot2, tilt2, psi2, rot2p, tilt2p, psi2p);
			}

			RFLOAT ang_dist = check_tilt_pairs(rot1, tilt1, psi1, rot2p, tilt2p, psi2p);

			if (ang_dist < best_ang_dist)
			{
				best_ang_dist = ang_dist;
				best_rot2 = rot2p;
				best_tilt2 = tilt2p;
				best_psi2 = psi2p;
			}

		}

		rot2 = best_rot2;
		tilt2 = best_tilt2;
		psi2 = best_psi2;

		return best_ang_dist;
	}


	RFLOAT check_tilt_pairs(RFLOAT rot1, RFLOAT tilt1, RFLOAT psi1,
	                        RFLOAT &alpha, RFLOAT &tilt_angle, RFLOAT &beta)
	{
		// Transformation matrices
		Matrix1D<RFLOAT> axis(3);
		Matrix2D<RFLOAT> E1, E2;
		axis.resize(3);
		RFLOAT aux, sine_tilt_angle;
		RFLOAT rot2 = alpha, tilt2 = tilt_angle, psi2 = beta;

		// Calculate the transformation from one setting to the second one.
		Euler_angles2matrix(psi1, tilt1, rot1, E1);
		Euler_angles2matrix(psi2, tilt2, rot2, E2);
		E2 = E2 * E1.inv();

		// Get the tilt angle (and its sine)
		aux = ( E2(0,0) + E2(1,1) + E2(2,2) - 1. ) / 2.;
		if (ABS(aux) - 1. > XMIPP_EQUAL_ACCURACY)
			REPORT_ERROR("BUG: aux>1");
		tilt_angle = ACOSD(aux);
		sine_tilt_angle = 2. * SIND(tilt_angle);

		// Get the tilt axis direction in angles alpha and beta
		if (sine_tilt_angle > XMIPP_EQUAL_ACCURACY)
		{
			axis(0) = ( E2(2,1) - E2(1,2) ) / sine_tilt_angle;
			axis(1) = ( E2(0,2) - E2(2,0) ) / sine_tilt_angle;
			axis(2) = ( E2(1,0) - E2(0,1) ) / sine_tilt_angle;
		}
		else
		{
			axis(0) = axis(1) = 0.;
			axis(2) = 1.;
		}

		// Apply E1.inv() to the axis to get everyone in the same coordinate system again
		axis = E1.inv() * axis;

		// Convert to alpha and beta angle
		Euler_direction2angles(axis, alpha, beta);

		// Enforce positive beta: choose the other Euler angle combination to express the same direction
		if (beta < 0.)
		{
			beta = -beta;
			alpha+= 180.;
		}

		// Let alpha go from 0 to 360 degrees
		alpha = realWRAP(alpha, 0., 360.);


		// Return the value that needs to be optimized
		RFLOAT minimizer=0.;
		if (exp_beta < 999.)
			minimizer = ABS(beta - exp_beta);
		if (exp_tilt < 999.)
			minimizer += ABS(tilt_angle - exp_tilt);

		return minimizer;
	}

	void run()
	{
		int iline = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDu)
		{


			// Read input data
			RFLOAT rot1,  tilt1,  psi1;
			RFLOAT rot2,  tilt2,  psi2;
			RFLOAT rot2p, tilt2p, psi2p;
			RFLOAT best_tilt, best_alpha, best_beta;
			RFLOAT distp;

			MDu.getValue(EMDL_ORIENT_ROT, rot1);
			MDt.getValue(EMDL_ORIENT_ROT, rot2, iline);
			MDu.getValue(EMDL_ORIENT_TILT, tilt1);
			MDt.getValue(EMDL_ORIENT_TILT, tilt2, iline);
			MDu.getValue(EMDL_ORIENT_PSI, psi1);
			MDt.getValue(EMDL_ORIENT_PSI, psi2, iline);
			iline++;

			// Bring both angles to a normalized set
			rot1 = realWRAP(rot1, -180, 180);
			tilt1 = realWRAP(tilt1, -180, 180);
			psi1 = realWRAP(psi1, -180, 180);
			rot2 = realWRAP(rot2, -180, 180);
			tilt2 = realWRAP(tilt2, -180, 180);
			psi2 = realWRAP(psi2, -180, 180);

			// Apply rotations to find the minimum distance angles
			rot2p = rot2;
			tilt2p = tilt2;
			psi2p = psi2;
			distp = check_symmetries(rot1, tilt1, psi1, rot2p, tilt2p, psi2p);

			// Calculate distance to user-defined point
			RFLOAT xp, yp, x, y;
			Matrix1D<RFLOAT> aux2(4);
			xp = dist_from_tilt * COSD(dist_from_alpha);
			yp = dist_from_tilt * SIND(dist_from_alpha);
			x = tilt2p * COSD(rot2p);
			y = tilt2p * SIND(rot2p);
			aux2(3) = sqrt((xp-x)*(xp-x) + (yp-y)*(yp-y));
			aux2(0) = tilt2p;
			aux2(1) = rot2p;
			aux2(2) = psi2p;
			add_to_postscript(tilt2p, rot2p, psi2p);

		}

		// Close the EPS file to write it to disk
		fh_eps << "showpage\n";
		fh_eps.close();
	}
};


int main(int argc, char *argv[])
{
	tiltpair_plot_parameters prm;

	try
	{
		prm.read(argc, argv);
		prm.initialise();
		prm.run();
	}
	catch (RelionError XE)
	{
		std::cerr << XE;
		//prm.usage();
		return RELION_EXIT_FAILURE;
	}
	return RELION_EXIT_SUCCESS;
}

