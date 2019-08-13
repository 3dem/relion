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
#include <src/fftw.h>
#include <src/args.h>
#include <src/error.h>
#include <src/mask.h>
#include <src/time.h>
#include <src/helix.h>

class mask_create_parameters
{
	public:
	FileName fn_apply_in, fn_mask, fn_apply_out, fn_thr, fn_omask, fn_and, fn_or, fn_andnot, fn_ornot;
	RFLOAT ini_threshold, extend_ini_mask, width_soft_edge, lowpass, angpix, helical_z_percentage;
	RFLOAT inner_radius, outer_radius, center_x, center_y, center_z;
	bool do_invert, do_helix, do_denovo;
	int n_threads, box_size;
	IOParser parser;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{

		parser.setCommandLine(argc, argv);

		int create_section = parser.addSection("Mask creation options");
		fn_thr = parser.getOption("--i", "Input map to use for thresholding to generate initial binary mask","");
		fn_omask = parser.getOption("--o", "Output mask","mask.mrc");

		fn_and = parser.getOption("--and", "Pixels in the initial mask will be one if the input AND this map are above the --ini_threshold value","");
		fn_or = parser.getOption("--or", "Pixels in the initial mask will be one if the input OR this map are above the --ini_threshold value","");
		fn_andnot = parser.getOption("--and_not", "Pixels in the initial mask will be one if the input is above the --ini_threshold AND this map is below it","");
		fn_ornot = parser.getOption("--or_not", "Pixels in the initial mask will be one if the input is above the --ini_threshold OR this map is below it","");
		ini_threshold  = textToFloat(parser.getOption("--ini_threshold", "Initial threshold for binarization","0.01"));
		extend_ini_mask = textToFloat(parser.getOption("--extend_inimask", "Extend initial binary mask this number of pixels","0"));
		width_soft_edge  = textToFloat(parser.getOption("--width_soft_edge", "Width (in pixels) of the additional soft edge on the binary mask", "0"));
		do_invert = parser.checkOption("--invert", "Invert the final mask");
		do_helix = parser.checkOption("--helix", "Generate a mask for 3D helix");
		lowpass = textToFloat(parser.getOption("--lowpass", "Lowpass filter (in Angstroms) for the input map, prior to binarization (default is none)", "-1"));
		angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms) for the lowpass filter", "-1"));
		helical_z_percentage = textToFloat(parser.getOption("--z_percentage", "This box length along the center of Z axis contains good information of the helix", "0.3"));
		n_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));

		int denovo_section = parser.addSection("De novo mask creation");
		do_denovo = parser.checkOption("--denovo", "Create a mask de novo");
		box_size = textToInteger(parser.getOption("--box_size", "The box size of the mask in pixels", "-1"));
		inner_radius = textToFloat(parser.getOption("--inner_radius", "Inner radius of the masked region in pixels", "0"));
		outer_radius = textToFloat(parser.getOption("--outer_radius", "Outer radius of the mask region in pixels", "99999"));
		center_x = textToFloat(parser.getOption("--center_x", "X coordinate of the center of the mask in pixels", "0"));
		center_y = textToFloat(parser.getOption("--center_y", "Y coordinate of the center of the mask in pixels", "0"));
		center_z = textToFloat(parser.getOption("--center_z", "Z coordinate of the center of the mask in pixels", "0"));

		// Check for errors in the command-line option
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line, exiting...");

		if (fn_thr == "" && fn_apply_in == "" && !do_denovo)
			REPORT_ERROR("Either provide --i to apply a mask, OR --create_from or --denovo to create a new mask");

		if (do_denovo && box_size < 0)
			REPORT_ERROR("For de novo mask creation, please specify the box size in --box_size");
	}


	void run()
	{
		Image<RFLOAT> Iout;		

		if (do_denovo)
		{
			makeMaskFromScratch(Iout);
			if (angpix < 0)
			{
				std::cerr << "WARNING: The pixel size was not specified. 1.00 is set to the output mask." << std::endl;
				angpix = 1.0;
			}
		}
		else
		{
			makeMaskFromFile(Iout);
		}

		if (do_invert)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iout())
			{
				DIRECT_MULTIDIM_ELEM(Iout(), n) = 1. - DIRECT_MULTIDIM_ELEM(Iout(), n);
			}
		}

		// Set header and write outmap map
		Iout.setStatisticsInHeader();
		Iout.setSamplingRateInHeader(angpix, angpix);
		Iout.write(fn_omask);
		std::cout << " Done creating mask! Written out: " << fn_omask << std::endl;
	}

	void makeMaskFromScratch(Image<RFLOAT> &Iout)
	{
		Iout().reshape(box_size, box_size, box_size);

		raisedCrownMask(Iout(), inner_radius, outer_radius, width_soft_edge, center_x, center_y, center_z);
	}

	void makeMaskFromFile(Image<RFLOAT> &Iout)
	{
		Image<RFLOAT> Iin, Ip;
		std:: cout << " Creating a mask ..." << std::endl;
		Iin.read(fn_thr);

		if (angpix < 0)
		{
			angpix = Iin.samplingRateX();
			std::cerr << "WARNING: The pixel size (--angpix) was not specified." << std::endl;
			std::cerr << "         The value in the input image header (= " << angpix << ") is used instead." << std::endl;
		}

		if (lowpass > 0)
		{
			lowPassFilterMap(Iin(), lowpass, angpix);
		}

		Iin().setXmippOrigin();

		if (fn_and != "")
		{
			Ip.read(fn_and);
			Ip().setXmippOrigin();
			if (!Ip().sameShape(Iin()))
				REPORT_ERROR("ERROR: --i and --and maps are different shapes!");
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Ip())
			{
				if (DIRECT_MULTIDIM_ELEM(Ip(), n) > ini_threshold && DIRECT_MULTIDIM_ELEM(Iin(), n) > ini_threshold)
					DIRECT_MULTIDIM_ELEM(Iin(), n) = ini_threshold + 1.;
				else
					DIRECT_MULTIDIM_ELEM(Iin(), n) = ini_threshold - 1.;
			}
		}
		else if  (fn_or != "")
		{
			Ip.read(fn_or);
			Ip().setXmippOrigin();
			if (!Ip().sameShape(Iin()))
				REPORT_ERROR("ERROR: --i and --or maps are different shapes!");
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Ip())
			{
				if (DIRECT_MULTIDIM_ELEM(Ip(), n) > ini_threshold || DIRECT_MULTIDIM_ELEM(Iin(), n) > ini_threshold)
					DIRECT_MULTIDIM_ELEM(Iin(), n) = ini_threshold + 1.;
				else
					DIRECT_MULTIDIM_ELEM(Iin(), n) = ini_threshold - 1.;
			}
		}
		else if  (fn_andnot != "")
		{
			Ip.read(fn_andnot);
			Ip().setXmippOrigin();
			if (!Ip().sameShape(Iin()))
				REPORT_ERROR("ERROR: --i and --not maps are different shapes!");
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Ip())
			{
				if (DIRECT_MULTIDIM_ELEM(Iin(), n) > ini_threshold && DIRECT_MULTIDIM_ELEM(Ip(), n) < ini_threshold)
					DIRECT_MULTIDIM_ELEM(Iin(), n) = ini_threshold + 1.;
				else
					DIRECT_MULTIDIM_ELEM(Iin(), n) = ini_threshold - 1.;
			}
		}
		else if  (fn_ornot != "")
		{
			Ip.read(fn_ornot);
			Ip().setXmippOrigin();
			if (!Ip().sameShape(Iin()))
				REPORT_ERROR("ERROR: --i and --not maps are different shapes!");
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Ip())
			{
				if (DIRECT_MULTIDIM_ELEM(Iin(), n) > ini_threshold || DIRECT_MULTIDIM_ELEM(Ip(), n) < ini_threshold)
					DIRECT_MULTIDIM_ELEM(Iin(), n) = ini_threshold + 1.;
				else
					DIRECT_MULTIDIM_ELEM(Iin(), n) = ini_threshold - 1.;
			}
		}

		autoMask(Iin(), Iout(), ini_threshold, extend_ini_mask, width_soft_edge, true, n_threads); // true sets verbosity

		if (do_helix)
		{
			cutZCentralPartOfSoftMask(Iout(), helical_z_percentage, width_soft_edge);
		}
	}
};


int main(int argc, char *argv[])
{
	mask_create_parameters prm;

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


