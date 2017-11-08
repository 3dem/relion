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

class stack_create_parameters
{
	public:
   	FileName fn_star, fn_root, fn_ext;
	MetaDataTable MD;
	// I/O Parser
	IOParser parser;
	bool do_spider, do_split_per_micrograph, do_apply_trans;

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
	    do_spider  = parser.checkOption("--spider_format", "Write out in SPIDER stack format (by default MRC stack format)");
	    do_split_per_micrograph = parser.checkOption("--split_per_micrograph", "Write out separate stacks for each micrograph (needs rlnMicrographName in STAR file)");
	    do_apply_trans = parser.checkOption("--apply_transformation", "Apply the inplane-transformations (needs _rlnOriginX/Y and _rlnAnglePsi in STAR file)");

	    fn_ext = (do_spider) ? ".spi" : ".mrcs";

      	// Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line, exiting...");
	}


	void run()
	{
		MD.read(fn_star);

		// Check for rlnImageName label
		if (!MD.containsLabel(EMDL_IMAGE_NAME))
			REPORT_ERROR("ERROR: Input STAR file does not contain the rlnImageName label");

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

			// Resize the output image
			std::cout << "Resizing the output stack to "<< ndim<<" images of size: "<<xdim<<"x"<<ydim<<"x"<<zdim << std::endl;
			RFLOAT Gb = (RFLOAT)ndim * zdim * ydim * xdim * sizeof(RFLOAT) / (1024. * 1024. * 1024.);
			std::cout << "This will require " << Gb << "Gb of memory...."<< std::endl;
			Image<RFLOAT> out(xdim, ydim, zdim, ndim);

			int n = 0;
			init_progress_bar(ndim);
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
			{
				FileName fn_mymic;
				if (do_split_per_micrograph)
					MD.getValue(EMDL_MICROGRAPH_NAME, fn_mymic);
				else
					fn_mymic="";

				if (fn_mymic == fn_mic)
				{

					MD.getValue(EMDL_IMAGE_NAME, fn_img);
					in.read(fn_img);

					if (do_apply_trans)
					{
						RFLOAT xoff = 0.;
						RFLOAT yoff = 0.;
						RFLOAT psi = 0.;
						MD.getValue(EMDL_ORIENT_ORIGIN_X, xoff);
						MD.getValue(EMDL_ORIENT_ORIGIN_Y, yoff);
						MD.getValue(EMDL_ORIENT_PSI, psi);
						// Apply the actual transformation
						Matrix2D<RFLOAT> A;
						rotation2DMatrix(psi, A);
					    MAT_ELEM(A,0, 2) = xoff;
					    MAT_ELEM(A,1, 2) = yoff;
					    selfApplyGeometry(in(), A, IS_NOT_INV, DONT_WRAP);
					}

					out().setImage(n, in());
					n++;
					if (n%100==0) progress_bar(n);

				}
			}
			progress_bar(ndim);


			FileName fn_out;
			if (do_split_per_micrograph)
			{
				// Remove any extensions from micrograph names....
				fn_out = fn_root + "_" + fn_mic.withoutExtension() + fn_ext;
			}
			else
				fn_out = fn_root + fn_ext;
			out.write(fn_out);
			std::cout << "Written out: " << fn_out << std::endl;
		}
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
        exit(1);
    }
    return 0;
}

