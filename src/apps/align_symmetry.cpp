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
	bool keep_centre, only_rot, do_apply_sym, do_select_size, do_select_resol;
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
		fn_in = parser.getOption("--i", "Input map to be aligned (OR: input model.star, from which input map will be selected)");
		do_select_size = parser.checkOption("--select_largest_class", "Select largest class from model.star; by default map with best symmetry will be selected");
		do_select_resol = parser.checkOption("--select_highest_resol", "Select class with highest resolution from model.star; by default map with best symmetry will be selected");
		fn_out = parser.getOption("--o", "Name for output aligned map)", "aligned.mrc");
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
		do_apply_sym = parser.checkOption("--apply_sym", "Also apply the symmetry to the map");

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

	int search(MetaDataTable &MDang, Projector &projector, double &best_diff2)
	{
		init_progress_bar(MDang.numberOfObjects());
		long int best_at = 0;
		best_diff2 = 1E99;
		RFLOAT rot, tilt, psi;

		// TODO: parallelise?
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDang)
		{
			MDang.getValue(EMDL_ORIENT_ROT, rot);
			MDang.getValue(EMDL_ORIENT_TILT, tilt);
			MDang.getValue(EMDL_ORIENT_PSI, psi);

			Euler_rotation3DMatrix(rot, tilt, psi, A3D);
			F2D.initZeros();
			projector.get2DFourierTransform(F2D, A3D);

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

		// When providing a model.star file, in C1 one cannot select class based on best diff2
		if (fn_in.contains("model.star") && (fn_sym.contains("C1") || fn_sym.contains("c1")))
		{
			if (!(do_select_resol ||do_select_size ))
				REPORT_ERROR("ERROR: In symmetry C1, one has to select based on class size or resolution!");
		}

		MetaDataTable MDang;
		Image<RFLOAT> vol_in, vol_work;
		int orig_size;
		RFLOAT work_angpix, r_max, rot, tilt, psi;

		std::vector<FileName> fn_ins;
		if (fn_in.contains("model.star"))
		{
			if (do_select_resol && do_select_size)
				REPORT_ERROR("ERROR: select either on resolution or on size, not both!");
			MetaDataTable MDclasses;
			MDclasses.read(fn_in, "model_classes");
			RFLOAT largest_classsize = -1.;
			RFLOAT best_resol = 10000.;
			FileName selected_class="empty";
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDclasses)
			{
				FileName fnt;
				MDclasses.getValue(EMDL_MLMODEL_REF_IMAGE, fnt);
				RFLOAT classsize, classresol;
				MDclasses.getValue(EMDL_MLMODEL_PDF_CLASS, classsize);
				MDclasses.getValue(EMDL_MLMODEL_ESTIM_RESOL_REF, classresol);
				if (do_select_size && classsize > largest_classsize)
				{
					largest_classsize = classsize;
					selected_class = fnt;
				}
				else if (do_select_resol && classresol < best_resol)
				{
					best_resol = classresol;
					selected_class = fnt;
				}
				else if (!(do_select_size || do_select_resol)) fn_ins.push_back(fnt);
			}
			if (do_select_size || do_select_resol) fn_ins.push_back(selected_class);
		}
		else
		{
			fn_ins.push_back(fn_in);

		}

		double best_best_diff2 = 1.E99;
		int best_imap = -1;
		for (int imap = 0; imap < fn_ins.size(); imap++)
		{

			std::cout << " Reading map: " << fn_ins[imap] << std::endl;
			vol_in.read(fn_ins[imap]);


			if ( fn_sym.contains("C1") || fn_sym.contains("c1") )
			{
				// Just write out the selected map again
				vol_in.write(fn_out);
				std::cout << " The selected map has been written to " << fn_out << std::endl;
			}
			else
			{
				orig_size = XSIZE(vol_in());
				std:: cout << " The input box size: " << orig_size << std::endl;
				if (orig_size % 2 != 0)
					REPORT_ERROR("The input box size must be an even number.");
				if (orig_size < boxsize)
				{
					std::cerr << "The working box size (--box_size) is larger (" << boxsize << ") than the input box size (" << orig_size << "). The working box size is increased to the input box size." << std::endl;
					boxsize = orig_size;
				}

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
				double best_diff2;
				int best_at;
				best_at = search(MDang, projector, best_diff2);

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
				best_at = search(MDang, projector, best_diff2);

				MDang.getValue(EMDL_ORIENT_ROT, rot, best_at);
				MDang.getValue(EMDL_ORIENT_TILT, tilt, best_at);
				MDang.getValue(EMDL_ORIENT_PSI, psi, best_at);
				std::cout << " The refined solution is ROT = " << rot << " TILT = " << tilt << " PSI = " << psi << " DIFF2= " << best_diff2 << std::endl << std::endl;

				if (best_diff2 < best_best_diff2)
				{
					best_best_diff2 = best_diff2;
					std::cout << " Now rotating the original (full size) volume ..." << std::endl << std::endl;
					Projector full_projector(orig_size, interpolator, padding_factor, r_min_nn, data_dim);
					Image<RFLOAT> vol_out;
					FourierTransformer final_transformer;

					full_projector.computeFourierTransformMap(vol_in(), dummy, 2 * orig_size);
					Euler_rotation3DMatrix(rot, tilt, psi, A3D);
					F2D.initZeros(orig_size, orig_size, orig_size / 2 + 1);
					vol_out().reshape(vol_in());
					full_projector.get2DFourierTransform(F2D, A3D);

					transformer.inverseFourierTransform(F2D, vol_out());
					CenterFFT(vol_out(), false);

					if (do_apply_sym) symmetriseMap(vol_out(), fn_sym);

					vol_out.setSamplingRateInHeader(angpix);
					vol_out.write(fn_out);
					std::cout << " The aligned map has been written to " << fn_out << std::endl;
				}
			}

		} // end loop over fn_ins

		std::cout << " done!" << std::endl;
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
	        return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}
