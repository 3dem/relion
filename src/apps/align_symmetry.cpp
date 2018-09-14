/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres" and "Takanori Nakane"
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

#include <src/projector.h>
#include <src/fftw.h>
#include <src/args.h>
#include <src/euler.h>
#include <src/transformations.h>
#include <src/symmetries.h>
#include <src/time.h>

class align_symmetry
{
private:
	Matrix2D<RFLOAT> A3D;
	MultidimArray<Complex> F2D;
	MultidimArray<RFLOAT> rotated, symmetrised, dummy;
	FourierTransformer transformer;

public:

	FileName fn_in, fn_out, fn_sym;
	RFLOAT angpix, maxres, search_step;
	int nr_uniform, padding_factor, interpolator, r_min_nn, boxsize, search_range;
	bool keep_centre, only_rot;
	// I/O Parser
	IOParser parser;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("Options");
		fn_in = parser.getOption("--i", "Input map to be projected");
		fn_out = parser.getOption("--o", "Rootname for output projections", "aligned.mrc");
		fn_sym = parser.getOption("--sym", "Target point group symmetry");
		boxsize = textToInteger(parser.getOption("--box_size", "Working box size in pixels. Very small box (such that Nyquist is aroud 20 A) is usually sufficient.", "64"));
		if (boxsize % 2 != 0)
			REPORT_ERROR("The working box size (--box_size) must be an even number.");
		keep_centre = parser.checkOption("--keep_centre", "Do not re-centre the input");
		angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "-1"));
		only_rot = parser.checkOption("--only_rot", "Keep TILT and PSI fixed and search only ROT (rotation along the Z axis)");
		nr_uniform = textToInteger(parser.getOption("--nr_uniform", "Randomly search this many orientations", "400"));
		maxres = textToFloat(parser.getOption("--maxres", "Maximum resolution (in Angstrom) to consider in Fourier space (default Nyquist)", "-1"));
		search_range = textToInteger(parser.getOption("--local_search_range", "Local search range (1 + 2 * this number)", "2"));
		search_step = textToFloat(parser.getOption("--local_search_step", "Local search step (in degrees)", "2"));
		padding_factor = textToInteger(parser.getOption("--pad", "Padding factor", "2"));
		
		if (parser.checkOption("--NN", "Use nearest-neighbour instead of linear interpolation"))
			interpolator = NEAREST_NEIGHBOUR;
		else
			interpolator = TRILINEAR;

		// Hidden
		r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));

		// Check for errors in the command-line option
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	int search(MetaDataTable &MDang, Projector &projector)
	{
		init_progress_bar(MDang.numberOfObjects());
		long int best_at = 0;
		double best_diff2 = 1E99;
		RFLOAT rot, tilt, psi;

		// TODO: parallelise?
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDang)
		{
			MDang.getValue(EMDL_ORIENT_ROT, rot);
			MDang.getValue(EMDL_ORIENT_TILT, tilt);
			MDang.getValue(EMDL_ORIENT_PSI, psi);

			Euler_rotation3DMatrix(rot, tilt, psi, A3D);
			F2D.initZeros();
			projector.get2DFourierTransform(F2D, A3D, IS_NOT_INV);

			transformer.inverseFourierTransform();
			CenterFFT(rotated, false);
			symmetrised = rotated;
			symmetriseMap(symmetrised, fn_sym);

			// non-weighted real-space squared difference
			double diff2 = 0;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(rotated)	
			{                                
				diff2 += (DIRECT_MULTIDIM_ELEM(rotated, n) - DIRECT_MULTIDIM_ELEM(symmetrised, n)) * 
				         (DIRECT_MULTIDIM_ELEM(rotated, n) - DIRECT_MULTIDIM_ELEM(symmetrised, n));
			}
			if (best_diff2 > diff2)
			{
				best_diff2 = diff2;
				best_at = current_object;
			}

			if (current_object % 30 == 0) progress_bar(current_object);
#ifdef DEBUG
			std::cout << rot << " " << tilt << " " << psi << " " << diff2 << std::endl; 
#endif
		} // end search

		progress_bar(MDang.numberOfObjects());

		return best_at;
	}

	void project()
	{
		MetaDataTable MDang;
		Image<RFLOAT> vol_in, vol_work;
		int orig_size;
		RFLOAT work_angpix, r_max, rot, tilt, psi;

		std::cout << " Reading map: " << fn_in << std::endl;
		vol_in.read(fn_in);
		orig_size = XSIZE(vol_in());
		std:: cout << " The input box size: " << orig_size << std::endl;
		if (orig_size % 2 != 0)
			REPORT_ERROR("The input box size must be an even number.");
		if (orig_size < boxsize)
			REPORT_ERROR("There is no point using the working box size (--box_size) larger than the input volume.");

		if (angpix < 0.)
		{
			angpix = vol_in.samplingRateX();
			std::cout << " Using the pixel size in the input image header: " << angpix << " A/px" << std::endl;
		}

		if (!keep_centre)
		{
			selfTranslateCenterOfMassToCenter(vol_in(), DONT_WRAP, true);
			std::cout << " Re-centred to the centre of the mass" << std::endl;
		}

		vol_work = vol_in;
		resizeMap(vol_work(), boxsize);
		work_angpix = angpix * orig_size / boxsize;
		std::cout << " Downsampled to the working box size " << boxsize << " px. This corresponds to " << work_angpix << " A/px." << std::endl;
		
		if (nr_uniform > 0)
		{
			std::cout << " Generating " << nr_uniform << " projections taken randomly from a uniform angular distribution." << std::endl;
			MDang.clear();
			randomize_random_generator();
			tilt = 0;
			psi = 0;
			for (long int i = 0; i < nr_uniform; i++)
			{
				rot = rnd_unif() * 360.;

				if (!only_rot)
				{
					bool ok_tilt = false;
					while (!ok_tilt)
					{
						tilt = rnd_unif() * 180.;
						if (rnd_unif() < fabs(SIND(tilt)))
							ok_tilt = true;
					}
					psi = rnd_unif() * 360.;
				}

				MDang.addObject();
				MDang.setValue(EMDL_ORIENT_ROT, rot);
				MDang.setValue(EMDL_ORIENT_TILT, tilt);
				MDang.setValue(EMDL_ORIENT_PSI, psi);
			}
		}

		// Now that we have the size of the volume, check r_max
		if (maxres < 0.)
			r_max = boxsize;
		else
			r_max = CEIL(boxsize * work_angpix / maxres);

		// Set right size of F2D and initialize to zero
		rotated.reshape(vol_work());
		symmetrised.reshape(vol_work());
		transformer.setReal(rotated);
		transformer.getFourierAlias(F2D);

		// Set up the projector
		int data_dim = 3;
		Projector projector(boxsize, interpolator, padding_factor, r_min_nn, data_dim);
		projector.computeFourierTransformMap(vol_work(), dummy, 2* r_max);

		// Global search
		std::cout << " Searching globally ..." << std::endl;
		int best_at;
		best_at = search(MDang, projector);

		MDang.getValue(EMDL_ORIENT_ROT, rot, best_at);
		MDang.getValue(EMDL_ORIENT_TILT, tilt, best_at);
		MDang.getValue(EMDL_ORIENT_PSI, psi, best_at);
		std::cout << " The best solution is ROT = " << rot << " TILT = " << tilt << " PSI = " << psi << std::endl << std::endl;

		// Local refinement
		std::cout << " Refining locally ..." << std::endl;	
		MDang.clear();

		for (int i = -search_range; i <= search_range; i++)
		{
			for (int j = -search_range; j <= search_range; j++)
			{
				if (only_rot && j != 0) continue;

				for (int k = -search_range; k <= search_range; k++)
				{
					if (only_rot && k != 0) continue;

					MDang.addObject();
					MDang.setValue(EMDL_ORIENT_ROT, rot + i * search_step);
					MDang.setValue(EMDL_ORIENT_TILT, tilt + j * search_step);
					MDang.setValue(EMDL_ORIENT_PSI, psi + k * search_range);
				}
			}
		}
		best_at = search(MDang, projector);
		
		MDang.getValue(EMDL_ORIENT_ROT, rot, best_at);
		MDang.getValue(EMDL_ORIENT_TILT, tilt, best_at);
		MDang.getValue(EMDL_ORIENT_PSI, psi, best_at);
		std::cout << " The refined solution is ROT = " << rot << " TILT = " << tilt << " PSI = " << psi << std::endl << std::endl;

		std::cout << " Now rotating the original (full size) volume ..." << std::endl << std::endl;
		Projector full_projector(orig_size, interpolator, padding_factor, r_min_nn, data_dim);
		Image<RFLOAT> vol_out;
		FourierTransformer final_transformer;

		full_projector.computeFourierTransformMap(vol_in(), dummy, 2 * orig_size);
		Euler_rotation3DMatrix(rot, tilt, psi, A3D);
		F2D.initZeros(orig_size, orig_size, orig_size / 2 + 1);
		vol_out().reshape(vol_in());
		full_projector.get2DFourierTransform(F2D, A3D, IS_NOT_INV);

		transformer.inverseFourierTransform(F2D, vol_out());
		CenterFFT(vol_out(), false);
		vol_out.setSamplingRateInHeader(angpix);
		vol_out.write(fn_out);
		std::cout << " The aligned map has been written to " << fn_out << std::endl;

	} // end project function
};

int main(int argc, char *argv[])
{
	align_symmetry app;

	try
	{
		app.read(argc, argv);
		app.project();
	}

	catch (RelionError XE)
	{
	        //prm.usage();
        	std::cerr << XE;
	        exit(1);
	}

	return 0;
}
