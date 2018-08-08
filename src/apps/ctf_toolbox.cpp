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
#include <src/args.h>
#include <src/fftw.h>
#include <src/ctf.h>
#include <src/time.h>
#include <src/symmetries.h>

class ctf_toolbox_parameters
{
	public:
   	FileName fn_in, fn_out;
   	bool do_ctf_phaseflip, do_ctf_multiply, do_1dprofile, do_2dimage, do_image_name, do_intact_ctf_first_peak;
	RFLOAT profile_angle, angpix;
   	int verb, my_size_x, my_size_y;

	// I/O Parser
	IOParser parser;

	MetaDataTable MD;
	FourierTransformer transformer;

	// Image size
	int xdim, ydim, zdim;
	long int ndim;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");
	    fn_in = parser.getOption("--i", "Input STAR file with CTF information");
	    fn_out = parser.getOption("--o", "Output rootname (for multiple images: insert this string before each image's extension)");
	    angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms (default is read from STAR file)", "-1"));
	    my_size_x = my_size_y = textToFloat(parser.getOption("--size", "Image size in pixels (default is read from STAR file)", "-1"));
	    do_image_name = parser.checkOption("--use_image_name", "Use rlnImageName  from the input STAR file? (default is use rlnMicrographName)");

	    int cst_section = parser.addSection("Options");
	    do_ctf_phaseflip = parser.checkOption("--ctf_phaseflip", "Phase-flip the image(s) in the input STAR file?");
	    do_ctf_multiply = parser.checkOption("--ctf_multiply", "Multiply the image(s) in the input STAR file with their CTF?");
	    do_2dimage  = parser.checkOption("--write_ctf_image", "Write out images with the 2D CTFs?");
	    do_1dprofile  = parser.checkOption("--write_ctf_profile", "Write out a STAR file with the 1D CTF profiles?");
	    profile_angle = textToFloat(parser.getOption("--1dprofile_angle", "Angle along which to calculate 1D CTF profiles (0=X, 90=Y)", "0"));
	    do_intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");

	    // Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

    	verb = 1;

	}



	void imageOperations(FileName fn_img, CTF &ctf, FileName my_fn_out)
	{

		Image<RFLOAT> img;
		img.read(fn_img);

		MultidimArray<Complex> Fimg;
		transformer.FourierTransform(img(), Fimg, false);

		RFLOAT xs = (RFLOAT)XSIZE(img()) * angpix;
		RFLOAT ys = (RFLOAT)YSIZE(img()) * angpix;
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(Fimg)
		{
			RFLOAT x = (RFLOAT)jp / xs;
			RFLOAT y = (RFLOAT)ip / ys;
			//bool do_abs = false, bool do_only_flip_phases = false, bool do_intact_until_first_peak = false, bool do_damping = true) const
			DIRECT_A2D_ELEM(Fimg, i, j) *= ctf.getCTF(x, y, false, do_ctf_phaseflip, do_intact_ctf_first_peak, false);
		}

		transformer.inverseFourierTransform(Fimg, img());

		// Write out the result
		// Check whether fn_out has an "@": if so REPLACE the corresponding frame in the output stack!
		long int n;
		FileName fn_tmp;
		my_fn_out.decompose(n, fn_tmp);
		n--;
		if (n >= 0) // This is a stack...
		{

			// The following assumes the images in the stack come ordered...
			if (n == 0)
				img.write(fn_tmp, n, true, WRITE_OVERWRITE); // make a new stack
			else
				img.write(fn_tmp, n, true, WRITE_APPEND);
		}
		else // individual image
		{
			img.write(my_fn_out);
		}

	}

	void writeCtf(CTF &ctf, FileName my_fn_out)
	{
		if (do_2dimage)
		{
			Image<RFLOAT> tmp;
			tmp().resize(my_size_y, my_size_x);
			ctf.getCenteredImage(tmp(), angpix, do_ctf_phaseflip, false, do_intact_ctf_first_peak, true);
			tmp.write(my_fn_out);
		}
		else if (do_1dprofile)
		{
			int mysize = XMIPP_MIN(my_size_y, my_size_x);
			MultidimArray<RFLOAT> ctf_profile;
			ctf_profile.resize(mysize);
			ctf.get1DProfile(ctf_profile, profile_angle, angpix, do_ctf_phaseflip, false, do_intact_ctf_first_peak, false);

			// Now save as a STAR file....
			MetaDataTable MDctf;
			MDctf.setName("ctf_profile");
			FOR_ALL_ELEMENTS_IN_ARRAY1D(ctf_profile)
			{
				// Only take one side of the 1D profile
				if (i >= 0)
				{
					MDctf.addObject();
					RFLOAT res = (i > 0) ? (mysize * angpix / (RFLOAT)i) : 999.;
					MDctf.setValue(EMDL_SPECTRAL_IDX, (int)i);
					MDctf.setValue(EMDL_RESOLUTION, 1./res);
					MDctf.setValue(EMDL_RESOLUTION_ANGSTROM, res);
					MDctf.setValue(EMDL_CTF_VALUE, A1D_ELEM(ctf_profile, i) );
				}
			}
			MDctf.write(my_fn_out);
		}
		else
			REPORT_ERROR("BUG: do_2dimage and do_1dprofile are both false; this should not happen!");

	}

	void run()
	{

		MD.read(fn_in);

   		if (angpix > 0.)
   		{
   			std::cout << " + Using user-provided pixel size: " << angpix << std::endl;
   		}
   		else if (MD.containsLabel(EMDL_CTF_MAGNIFICATION) && MD.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
    	{
    		RFLOAT mag, dstep;
   			MD.getValue(EMDL_CTF_MAGNIFICATION, mag);
   			MD.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
   			angpix = 10000. * dstep / mag;
   			std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << std::endl;
    	}
   		else
   			REPORT_ERROR("ERROR: provide angpix through magnification and detector pixel size in the input STAR file, or through the --angpix command line argument");

   		if (my_size_x > 0. && my_size_y > 0.)
   		{
   			std::cout << " + Using user-provided image size: " << my_size_x << std::endl;
   		}
   		else if ((do_image_name && MD.containsLabel(EMDL_IMAGE_NAME)) || MD.containsLabel(EMDL_MICROGRAPH_NAME) )
    	{
    		FileName fn_img;
    		if (do_image_name)
    			MD.getValue(EMDL_IMAGE_NAME, fn_img);
    		else
    			MD.getValue(EMDL_MICROGRAPH_NAME, fn_img);

    		Image<RFLOAT> img;
    		img.read(fn_img, false);
    		my_size_x = XSIZE(img());
    		my_size_y = YSIZE(img());
    	}
   		else
   			REPORT_ERROR("ERROR: provide a STAR file with rlnImageName or rlnMicrographName, or provide the image size through the --img_size command line argument");

  		long int i_img = 0;
		if (verb > 0)
   			init_progress_bar(MD.numberOfObjects());

  		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{

  			CTF ctf;
			ctf.read(MD, MD);

			FileName fn_img, my_fn_out;
    		if (do_image_name)
    			MD.getValue(EMDL_IMAGE_NAME, fn_img);
    		else
    			MD.getValue(EMDL_MICROGRAPH_NAME, fn_img);

    		if (do_2dimage || do_ctf_multiply || do_ctf_phaseflip)
				my_fn_out = fn_img.insertBeforeExtension("_" + fn_out);
			else
				my_fn_out = fn_img.withoutExtension() + "_" + fn_out + ".star";


			// Now do the actual work
    		if (do_ctf_multiply || do_ctf_phaseflip)
			{
				MD.getValue(EMDL_IMAGE_NAME, fn_img);
				imageOperations(fn_img, ctf, my_fn_out);
			}
			else if (do_1dprofile || do_2dimage)
			{
				writeCtf(ctf, my_fn_out);
			}

			i_img++;
			if (verb > 0)
				progress_bar(i_img);

		}

		if (verb > 0)
			progress_bar(MD.numberOfObjects());

	}

};


int main(int argc, char *argv[])
{
	ctf_toolbox_parameters prm;

	try
    {

		prm.read(argc, argv);

		prm.run();

    }
    catch (RelionError XE)
    {
        prm.usage();
        std::cerr << XE;
        exit(1);
    }
    return 0;
}



