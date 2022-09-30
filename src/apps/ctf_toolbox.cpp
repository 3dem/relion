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
	bool do_intact_ctf_until_first_peak, do_intact_ctf_after_first_peak, do_ctf_pad, do_apply_orient;
	RFLOAT profile_angle, sim_angpix, kV, Q0, Cs, defU, defV, defAng, phase_shift;
	int verb;

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

		int general_section = parser.addSection("Pre-multiply options");
		fn_in = parser.getOption("--i", "Input STAR file with CTF information", "");
		fn_out = parser.getOption("--o", "Output rootname (for multiple images: insert this string before each image's extension)", "");
        do_apply_orient = parser.checkOption("--apply_orient", "Also apply the in-plane rotation and translation to the CTF-premultiplied images");

		int sim_section = parser.addSection("OR: simulate options");
		fn_sim = parser.getOption("--simulate", "Output name for simulated CTF image","");
		sim_angpix = textToFloat(parser.getOption("--angpix", "Pixel size (A)", "1."));
		sim_box = textToInteger(parser.getOption("--box", "Box size (pix)", "256"));
		kV = textToFloat(parser.getOption("--kV", "Voltage (kV)", "300"));
		Q0 = textToFloat(parser.getOption("--Q0", "Amplitude contrast", "0.1"));
		Cs = textToFloat(parser.getOption("--Cs", "Spherical aberration (mm)", "2.7"));
		defU = textToFloat(parser.getOption("--defU", "Defocus in U-direction (A)", "20000"));
		defV = textToFloat(parser.getOption("--defV", "Defocus in V-direction (A, default = defU)", "-1."));
		if (!checkParameter(argc, argv, "--defV"))
			defV = defU;

		defAng = textToFloat(parser.getOption("--defAng", "Defocus angle (deg)", "0."));
		phase_shift  = textToFloat(parser.getOption("--phase_shift", "Phase shift (deg)", "0."));

		int cst_section = parser.addSection("Shared options");
		do_intact_ctf_until_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");
		do_intact_ctf_after_first_peak = parser.checkOption("--ctf_intact_after_first_peak", "Leave CTFs intact after first peak");
		do_ctf_pad = parser.checkOption("--ctf_pad", "Pre-multiply with a 2x finer-sampled CTF that is then downscaled");
		// Check for errors in the command-line option
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

		verb = 1;

	}

	void run()
	{

		// CTF Simulation of a single image
		if (fn_sim != "")
		{
			if (do_ctf_pad) sim_box_large = 2 * sim_box;
			else sim_box_large = sim_box;

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
			std::cout << " +  use CTF padding? " << ((do_ctf_pad) ? "true" : "false") << std::endl;
			std::cout << " + " << std::endl;
			ctf.setValues(defU, defV, defAng, kV, Cs, Q0, 0., 1., phase_shift);

			Ictf().setXmippOrigin();
			RFLOAT xs = (RFLOAT)sim_box_large * sim_angpix;
			RFLOAT ys = (RFLOAT)sim_box_large * sim_angpix;
			FOR_ALL_ELEMENTS_IN_ARRAY2D(Ictf())
			{
				RFLOAT x = (RFLOAT)j / xs;
				RFLOAT y = (RFLOAT)i / ys;

				A2D_ELEM(Ictf(), i, j) = ctf.getCTF(x, y, false, false, do_intact_ctf_until_first_peak, true, 0.0, do_intact_ctf_after_first_peak);
			}

			resizeMap(Ictf(), sim_box);

			Ictf.write(fn_sim);
			std::cout << " + Done! written: " << fn_sim << std::endl;

		}

		else
		{

			ObservationModel::loadSafely(fn_in, obsModel, MD);
			bool do_mic_name = (obsModel.opticsMdt.getName() == "micrographs");

			long int i_img = 0;
			if (verb > 0)
				init_progress_bar(MD.numberOfObjects());

			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
			{

				CTF ctf;
				ctf.readByGroup(MD, &obsModel);
				int og = obsModel.getOpticsGroup(MD);
				RFLOAT angpix = obsModel.getPixelSize(og);

				FileName fn_img, my_fn_out;
				if (do_mic_name)
					MD.getValue(EMDL_MICROGRAPH_NAME, fn_img);
				else
					MD.getValue(EMDL_IMAGE_NAME, fn_img);

				my_fn_out = fn_img.insertBeforeExtension("_" + fn_out);

				// Now do the actual work
				MD.getValue(EMDL_IMAGE_NAME, fn_img);


				Image<RFLOAT> img;
				MultidimArray<Complex> Fimg;
				MultidimArray<RFLOAT> Fctf;

				img.read(fn_img);

				transformer.FourierTransform(img(), Fimg, false);

				Fctf.resize(YSIZE(Fimg), XSIZE(Fimg));
				ctf.getFftwImage(Fctf, XSIZE(img()), YSIZE(img()), angpix, false, false, do_intact_ctf_until_first_peak, false, do_ctf_pad, do_intact_ctf_after_first_peak);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
				{
					if (!do_intact_ctf_after_first_peak)
						DIRECT_MULTIDIM_ELEM(Fimg, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					else
					{
						// this is safe because getCTF does not return RELION_EXIT_SUCCESS.
						DIRECT_MULTIDIM_ELEM(Fimg, n) /= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
				}

				transformer.inverseFourierTransform(Fimg, img());

                if (do_apply_orient)
                {
                    RFLOAT psi;
					Matrix1D<RFLOAT> offset(2);
					Matrix2D<RFLOAT> A;
					MD.getValue(EMDL_ORIENT_PSI, psi);
					MD.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, XX(offset));
					MD.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, YY(offset));
                    offset /= angpix;
                    rotation2DMatrix(psi, A);
                    MAT_ELEM(A, 0, 2) = COSD(psi) * XX(offset) - SIND(psi) * YY(offset);
                    MAT_ELEM(A, 1, 2) = COSD(psi) * YY(offset) + SIND(psi) * XX(offset);
                    selfApplyGeometry(img(), A, IS_NOT_INV, DONT_WRAP);
                }

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
					img.write(my_fn_out);

				MD.setValue(EMDL_IMAGE_NAME, my_fn_out);
				obsModel.setCtfPremultiplied(og, true);

				i_img++;
				if (verb > 0)
					progress_bar(i_img);
			}

			if (verb > 0)
				progress_bar(MD.numberOfObjects());

			obsModel.save(MD, fn_in.insertBeforeExtension("_"+fn_out));
			std::cout << " + written out new particles STAR file in: " << fn_in.insertBeforeExtension("_"+fn_out) << std::endl;
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
		return RELION_EXIT_FAILURE;
	}
	return RELION_EXIT_SUCCESS;
}
