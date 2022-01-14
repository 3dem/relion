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

#include <src/image.h>
#include <src/funcs.h>
#include <src/ctf.h>
#include <src/args.h>
#include <src/error.h>
#include <src/euler.h>
#include <src/time.h>
#include <src/jaz/single_particle/obs_model.h>

// TODO: set pixel sizes in the outputs

class stack_create_parameters
{
	public:
   	FileName fn_star, fn_root;
	MetaDataTable MD;
	// I/O Parser
	IOParser parser;
	bool do_split_per_micrograph, do_apply_trans, do_apply_trans_only, do_ignore_optics, do_one_by_one;
	ObservationModel obsModel;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");
		fn_star = parser.getOption("--i", "Input STAR file with the images (as rlnImageName) to be saved in a stack");
		fn_root = parser.getOption("--o", "Output rootname","output");
		do_split_per_micrograph = parser.checkOption("--split_per_micrograph", "Write out separate stacks for each micrograph (needs rlnMicrographName in STAR file)");
		do_apply_trans = parser.checkOption("--apply_transformation", "Apply the inplane-transformations (needs _rlnOriginX/Y and _rlnAnglePsi in STAR file) by real space interpolation");
		do_apply_trans_only = parser.checkOption("--apply_rounded_offsets_only", "Apply the rounded translations only (so-recentering without interpolation; needs _rlnOriginX/Y in STAR file)");
		do_ignore_optics = parser.checkOption("--ignore_optics", "Ignore optics groups. This allows you to read and write RELION 3.0 STAR files but does NOT allow you to convert 3.1 STAR files back to the 3.0 format.");
		do_one_by_one = parser.checkOption("--one_by_one", "Write particles one by one. This saves memory but can be slower.");

		if (do_apply_trans)
			std::cerr << "WARNING: --apply_transformation uses real space interpolation. It also invalidates CTF parameters (e.g. beam tilt & astigmatism). This can degrade the resolution. USE WITH CARE!!" << std::endl;

		// Check for errors in the command-line option
		if (parser.checkForErrors())
    			REPORT_ERROR("Errors encountered on the command line, exiting...");
	}

	void run()
	{
		if (do_ignore_optics && (do_apply_trans || do_apply_trans_only))
			REPORT_ERROR("ERROR: you cannot ignore optics and apply transformations");

		if (do_ignore_optics) MD.read(fn_star);
		else ObservationModel::loadSafely(fn_star, obsModel, MD, "particles");

		// Check for rlnImageName label
		if (!MD.containsLabel(EMDL_IMAGE_NAME))
			REPORT_ERROR("ERROR: Input STAR file does not contain the rlnImageName label. Aren't you reading RELION 3.1 STAR files with --ignore_optics?");

		if (do_split_per_micrograph && !MD.containsLabel(EMDL_MICROGRAPH_NAME))
			REPORT_ERROR("ERROR: Input STAR file does not contain the rlnMicrographName label");

		Image<RFLOAT> in;
		FileName fn_img, fn_mic;
		std::vector<FileName> fn_mics;
		std::vector<int> mics_ndims;

		// First get number of images and their size
		int ndim=0;
		bool is_first=true;
		int xdim, ydim, zdim;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			if (is_first)
			{
				MD.getValue(EMDL_IMAGE_NAME, fn_img);
				in.read(fn_img);
				xdim=XSIZE(in());
				ydim=YSIZE(in());
				zdim=ZSIZE(in());
				is_first=false;
			}

			if (do_split_per_micrograph)
			{
				MD.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
				bool have_found = false;
				for (int m = 0; m < fn_mics.size(); m++)
				{
					if (fn_mic == fn_mics[m])
					{
						have_found = true;
						mics_ndims[m]++;
						break;
					}
				}
				if (!have_found)
				{
					fn_mics.push_back(fn_mic);
					mics_ndims.push_back(1);
				}
			}
			ndim++;
		}

		// If not splitting, just fill fn_mics and mics_ndim with one entry (to re-use loop below)
		if (!do_split_per_micrograph)
		{
			fn_mics.push_back("");
			mics_ndims.push_back(ndim);
		}

		// Loop over all micrographs
		for (int m = 0; m < fn_mics.size(); m++)
		{
			ndim = mics_ndims[m];
			fn_mic = fn_mics[m];

			Image<RFLOAT> out;

			if (!do_one_by_one)
			{
				// Resize the output image
				std::cout << "Resizing the output stack to "<< ndim<<" images of size: "<<xdim<<"x"<<ydim<<"x"<<zdim << std::endl;
				RFLOAT Gb = (RFLOAT)ndim * zdim * ydim * xdim * sizeof(RFLOAT) / (1024. * 1024. * 1024.);
				std::cout << "This will require " << Gb << "Gb of memory...." << std::endl;
				std::cout << "If this runs out of memory, please try --one_by_one." << std::endl;
				out().reshape(ndim, zdim, ydim, xdim);
				// NOTE: MultidimArray::reshape takes NZYX, while Image constructor takes XYZN !!
			}

			FileName fn_out;
			if (do_split_per_micrograph)
			{
				// Remove any extensions from micrograph names....
				fn_out = fn_root + "_" + fn_mic.withoutExtension() + ".mrcs";
			}
			else
				fn_out = fn_root + ".mrcs";

			// Make all output directories if necessary
			if (fn_out.contains("/"))
			{
				mktree(fn_out.beforeLastOf("/"));
			}

			int n = 0;
			init_progress_bar(ndim);
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
			{
				FileName fn_mymic;
				if (do_split_per_micrograph)
					MD.getValue(EMDL_MICROGRAPH_NAME, fn_mymic);
				else
					fn_mymic="";

				int optics_group;
				MD.getValue(EMDL_IMAGE_OPTICS_GROUP, optics_group);
				optics_group--;
				RFLOAT angpix;
				if (do_ignore_optics) angpix = 1.0;
				else angpix = obsModel.getPixelSize(optics_group);

				if (fn_mymic == fn_mic)
				{

					MD.getValue(EMDL_IMAGE_NAME, fn_img);
					in.read(fn_img);

					if (do_apply_trans || do_apply_trans_only)
					{
						RFLOAT xoff, ori_xoff;
						RFLOAT yoff, ori_yoff;
						RFLOAT psi, ori_psi;
						MD.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, ori_xoff);
						MD.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, ori_yoff);
						MD.getValue(EMDL_ORIENT_PSI, ori_psi);
						ori_xoff /= angpix;
						ori_yoff /= angpix;

						if (do_apply_trans_only)
						{
							xoff = ROUND(ori_xoff);
							yoff = ROUND(ori_yoff);
							psi = 0.;
						}
						else
						{
							xoff = ori_xoff;
							yoff = ori_yoff;
							psi = ori_psi;
						}

						// Apply the actual transformation
						Matrix2D<RFLOAT> A;
						rotation2DMatrix(psi, A);
						MAT_ELEM(A, 0, 2) = COSD(psi) * xoff - SIND(psi) * yoff;
						MAT_ELEM(A, 1, 2) = COSD(psi) * yoff + SIND(psi) * xoff;
						selfApplyGeometry(in(), A, IS_NOT_INV, DONT_WRAP);

						MD.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, (ori_xoff - xoff)*angpix);
						MD.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, (ori_yoff - yoff)*angpix);
						MD.setValue(EMDL_ORIENT_PSI, ori_psi - psi);
					}
					FileName fn_img;
					fn_img.compose(n+1, fn_out);
					MD.setValue(EMDL_IMAGE_NAME, fn_img);

					if (!do_one_by_one)
					{
						out().setImage(n, in());
					}
					else
					{
						if (n == 0)
							in.write(fn_img, -1, false, WRITE_OVERWRITE);
						else
							in.write(fn_img, -1, true, WRITE_APPEND);
					}

					n++;
					if (n%100==0) progress_bar(n);
				}
			}
			progress_bar(ndim);

			if (!do_one_by_one)
				out.write(fn_out);
			std::cout << "Written out: " << fn_out << std::endl;
		}

		if (do_ignore_optics) MD.write(fn_root+".star");
		else obsModel.save(MD, fn_root+".star", "particles");
		std::cout << "Written out: " << fn_root << ".star" << std::endl;
		std::cout << "Done!" <<std::endl;
	}
};


int main(int argc, char *argv[])
{
	stack_create_parameters prm;

	try
	{
		prm.read(argc, argv);
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

