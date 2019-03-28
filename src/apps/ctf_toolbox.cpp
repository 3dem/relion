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
	FileName fn_in, fn_out, fn_sim;
	bool do_ctf_phaseflip, do_ctf_multiply, do_1dprofile, do_2dimage, do_image_name, do_intact_ctf_first_peak;
	RFLOAT profile_angle, angpix, sim_angpix, kV, Q0, Cs, defU, defV, defAng, phase_shift;
	int verb, my_size_x, my_size_y;

	// I/O Parser
	IOParser parser;

	MetaDataTable MD;
	FourierTransformer transformer;
	ObservationModel obsModel;

	// Image size
	int xdim, ydim, zdim, sim_box, sim_box_large;
	long int ndim;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");
		fn_in = parser.getOption("--i", "Input STAR file with CTF information", "");
		fn_out = parser.getOption("--o", "Output rootname (for multiple images: insert this string before each image's extension)", "");
		angpix = textToFloat(parser.getOption("--my_angpix", "Pixel size in Angstroms (default is read from STAR file)", "-1"));
		my_size_x = my_size_y = textToFloat(parser.getOption("--size", "Image size in pixels (default is read from STAR file)", "-1"));
		do_image_name = parser.checkOption("--use_image_name", "Use rlnImageName  from the input STAR file? (default is use rlnMicrographName)");

		int cst_section = parser.addSection("Options");
		do_ctf_phaseflip = parser.checkOption("--ctf_phaseflip", "Phase-flip the image(s) in the input STAR file?");
		do_ctf_multiply = parser.checkOption("--ctf_multiply", "Multiply the image(s) in the input STAR file with their CTF?");
		do_2dimage = parser.checkOption("--write_ctf_image", "Write out images with the 2D CTFs?");
		do_1dprofile = parser.checkOption("--write_ctf_profile", "Write out a STAR file with the 1D CTF profiles?");
		profile_angle = textToFloat(parser.getOption("--1dprofile_angle", "Angle along which to calculate 1D CTF profiles (0=X, 90=Y)", "0"));
		do_intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");

		int sim_section = parser.addSection("Simulate options");
		fn_sim = parser.getOption("--simulate", "Output name for simulated CTF image","");
		sim_angpix = textToFloat(parser.getOption("--angpix", "Pixel size (A)", "1."));
		sim_box = textToInteger(parser.getOption("--box", "Box size (pix)", "256"));
		sim_box_large = textToInteger(parser.getOption("--large_box", "Firts simulate in large box, then  downscale (pix)", "-1"));
		kV = textToFloat(parser.getOption("--kV", "Voltage (kV)", "300"));
		Q0 = textToFloat(parser.getOption("--Q0", "Amplitude contrast", "0.1"));
		Cs = textToFloat(parser.getOption("--Cs", "Spherical aberration (mm)", "2.7"));
		defU = textToFloat(parser.getOption("--defU", "Defocus in U-direction (A)", "20000"));
		defV = textToFloat(parser.getOption("--defU", "Defocus in V-direction (A, default = defU)", "-1."));
		defAng = textToFloat(parser.getOption("--defAng", "Defocus angle (deg)", "0."));
		phase_shift  = textToFloat(parser.getOption("--phase_shift", "Phase shift (deg)", "0."));

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
		MultidimArray<RFLOAT> Fctf;
		transformer.FourierTransform(img(), Fimg, false);
		Fctf.resize(YSIZE(Fimg), XSIZE(Fimg));
		ctf.getFftwImage(Fctf, XSIZE(img()), YSIZE(img()), angpix, false, do_ctf_phaseflip, do_intact_ctf_first_peak, false);

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
		{
			DIRECT_MULTIDIM_ELEM(Fimg, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
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

		// CTF Simulation of a single image
		if (fn_sim != "")
		{
			if (sim_box_large < 0) sim_box_large = sim_box;
			Image<RFLOAT> Ictf(sim_box_large, sim_box_large);
			CTF ctf;
			std::cout << " + Input values: " << std::endl;
			std::cout << " +  kV= " << kV << std::endl;
			std::cout << " +  Cs= " << Cs << std::endl;
			std::cout << " +  Q0= " << Q0 << std::endl;
			std::cout << " +  defU= " << defU << std::endl;
			std::cout << " +  defV= " << defV << std::endl;
			std::cout << " +  defAng= " << defAng << std::endl;
			std::cout << " +  phase_shift = " << phase_shift << std::endl;
			std::cout << " +  angpix= " << sim_angpix<< std::endl;
			std::cout << " +  box= " << sim_box<< std::endl;
			std::cout << " +  large_box= " << sim_box_large << std::endl;
			std::cout << " + " << std::endl;
			ctf.setValues(defU, defV, defAng, kV, Cs, Q0, 0., 1., phase_shift);

			Ictf().setXmippOrigin();
			RFLOAT xs = (RFLOAT)sim_box_large * sim_angpix;
			RFLOAT ys = (RFLOAT)sim_box_large * sim_angpix;
			FOR_ALL_ELEMENTS_IN_ARRAY2D(Ictf())
			{
				RFLOAT x = (RFLOAT)j / xs;
				RFLOAT y = (RFLOAT)i / ys;

				A2D_ELEM(Ictf(), i, j) = ctf.getCTF(x, y);
			}

			resizeMap(Ictf(), sim_box);

			Ictf.write(fn_sim);
			std::cout << " + Done! written: " << fn_sim << std::endl;

		}

		else
		{
			ObservationModel::loadSafely(fn_in, obsModel, MD);


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
				ctf.readByGroup(MD, &obsModel);
				angpix = obsModel.getPixelSize(obsModel.getOpticsGroup(MD));

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
					MD.setValue(EMDL_IMAGE_NAME, my_fn_out);
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

			if (do_ctf_multiply || do_ctf_phaseflip)
			{
				obsModel.save(MD, fn_in.insertBeforeExtension("_"+fn_out));
				std::cout << " + written out new particles STAR file in: " << fn_in.insertBeforeExtension("_"+fn_out) << std::endl;
			}
		}
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



