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
	bool do_invert, do_helix;
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
	    angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms) for the lowpass filter", "1"));
	    helical_z_percentage = textToFloat(parser.getOption("--z_percentage", "This box length along the center of Z axis contains good information of the helix", "0.3"));

	    // Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line, exiting...");

    	if (fn_thr=="" && fn_apply_in == "")
	    	REPORT_ERROR("Either provide --i to apply a mask, OR --create_from to create a new mask");

	}


	void run()
	{

		Image<RFLOAT> Iin, Iout, Ip;
		std:: cout << " Creating a mask ..." << std::endl;
		Iin.read(fn_thr);

		if (lowpass > 0)
		{
			lowPassFilterMap(Iin(), lowpass, angpix);
		}

		if (fn_and != "")
		{
			Ip.read(fn_and);
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

		autoMask(Iin(), Iout(), ini_threshold, extend_ini_mask, width_soft_edge, true); // true sets verbosity

		if (do_helix)
		{
			cutZCentralPartOfSoftMask(Iout(), helical_z_percentage, width_soft_edge);
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
		Iout.setSamplingRateInHeader(Iin.samplingRateX(), Iin.samplingRateY());
		Iout.write(fn_omask);
		std::cout << " Done creating mask! Written out: " << fn_omask << std::endl;
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
        exit(1);
    }
    return 0;
}


