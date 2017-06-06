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
#include <src/ml_optimiser.h>

class particle_reposition_parameters
{
public:

	FileName fn_in, fn_opt, fn_out, fn_dat;

	RFLOAT micrograph_background;
	bool do_invert;

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

		fn_in  = parser.getOption("--i", "Input STAR file with rlnMicrographName's ");
		fn_out = parser.getOption("--o", "Output rootname, to be added to input micrograph names", "reposition");
       	fn_opt = parser.getOption("--opt", "Optimiser STAR file with the 2D classes or 3D maps to be repositioned");
       	fn_dat = parser.getOption("--data", "Data STAR file with selected particles (default is to use all particles)", "");
       	micrograph_background = textToFloat(parser.getOption("--background", "The fraction of micrograph background noise in the output micrograph", "0.1"));
       	do_invert= parser.checkOption("--invert", "Invert the contrast in the references?");

       	// Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	}

	void run()
	{

		int xdim, ydim, radius;
		MetaDataTable DFi, DFopt;
		DFi.read(fn_in);

		MlOptimiser optimiser;
		optimiser.do_preread_images = false;

		std::cerr << "Reading optimiser ..." << std::endl;
		optimiser.read(fn_opt);
		optimiser.mymodel.setFourierTransformMaps(false);
		radius = optimiser.particle_diameter / (2. * optimiser.mymodel.pixel_size);

		// Use a user-provided subset of particles instead of all of them?
		if (fn_dat != "")
		{
			MetaDataTable MDdata;
			MDdata.read(fn_dat);
			optimiser.mydata.MDimg = MDdata;
		}

		// Prepare transformer
		MultidimArray<Complex > Fref, Frefp, Fref_current_size;
		MultidimArray<RFLOAT> Mref;

		if (optimiser.mymodel.data_dim == 3)
		{
			Mref.resize(optimiser.mymodel.ori_size, optimiser.mymodel.ori_size, optimiser.mymodel.ori_size);
			Fref.resize(optimiser.mymodel.ori_size, optimiser.mymodel.ori_size, optimiser.mymodel.ori_size/2 + 1);
			Fref_current_size.resize(optimiser.mymodel.current_size, optimiser.mymodel.current_size, optimiser.mymodel.current_size/2 + 1);
		}
		else
		{
			Mref.resize(optimiser.mymodel.ori_size, optimiser.mymodel.ori_size);
			Fref.resize(optimiser.mymodel.ori_size, optimiser.mymodel.ori_size/2 + 1);
			Fref_current_size.resize(optimiser.mymodel.current_size, optimiser.mymodel.current_size/2 + 1);
		}
		Frefp.resize(Fref);

		// Loop over all micrographs
		int barstep = XMIPP_MAX(1, DFi.numberOfObjects()/ 60);
		init_progress_bar(DFi.numberOfObjects());
		long int imgno = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(DFi)
		{

			FileName fn_mic, fn_mic_out, fn_coord_out;
			DFi.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
			fn_mic_out = fn_mic.insertBeforeExtension("_" + fn_out);
			fn_coord_out = fn_mic_out.withoutExtension()+ "_coord.star";

			FourierTransformer transformer;
			MetaDataTable MDcoord;
			Image<RFLOAT> Imic_in, Imic_out;
			MultidimArray<RFLOAT> Imic_sum;
			// Read in the first micrograph
			Imic_in.read(fn_mic);
			Imic_in().setXmippOrigin();
			Imic_out().initZeros(Imic_in());
			Imic_sum.initZeros(Imic_in());
			// Get mean and stddev of the input micrograph
			RFLOAT stddev_mic, mean_mic, dummy;
			Imic_in().computeStats(mean_mic, stddev_mic, dummy, dummy);
			// Loop over all particles
			for(long int current_object2 = (optimiser.mydata.MDimg).firstObject(); \
			current_object2 != MetaDataTable::NO_MORE_OBJECTS && current_object2!= MetaDataTable::NO_OBJECTS_STORED; \
			current_object2=(optimiser.mydata.MDimg).nextObject())
			{

				FileName fn_mic2;
				optimiser.mydata.MDimg.getValue(EMDL_MICROGRAPH_NAME, fn_mic2);
				if (fn_mic2 == fn_mic)
				{
					RFLOAT rot, tilt, psi, xcoord=0., ycoord=0., zcoord=0.;
					int iclass;
					Matrix2D<RFLOAT> A;
					Matrix1D<RFLOAT> offsets(3);

					MDcoord.addObject();
		            MDcoord.setObject(optimiser.mydata.MDimg.getObject());
		            MDcoord.setValue(EMDL_MICROGRAPH_NAME,fn_mic_out);

					optimiser.mydata.MDimg.getValue(EMDL_IMAGE_COORD_X,  xcoord);
					optimiser.mydata.MDimg.getValue(EMDL_IMAGE_COORD_Y,  ycoord);
					optimiser.mydata.MDimg.getValue(EMDL_ORIENT_ROT,  rot);
					optimiser.mydata.MDimg.getValue(EMDL_ORIENT_TILT, tilt);
					optimiser.mydata.MDimg.getValue(EMDL_ORIENT_PSI,  psi);
					optimiser.mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, XX(offsets));
					optimiser.mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, YY(offsets));
					if (optimiser.mymodel.data_dim == 3)
					{
						optimiser.mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Z, ZZ(offsets));
						optimiser.mydata.MDimg.getValue(EMDL_IMAGE_COORD_Z,  zcoord);
					}
					else
					{
						ZZ(offsets) = zcoord = 0.;
					}

					optimiser.mydata.MDimg.getValue(EMDL_PARTICLE_CLASS, iclass);
					Euler_angles2matrix(rot, tilt, psi, A, false);

					// Get the 2D image (in its ori_size)
					(optimiser.mymodel.PPref[iclass-1]).get2DFourierTransform(Fref_current_size, A, IS_NOT_INV);
					windowFourierTransform(Fref_current_size, Frefp, optimiser.mymodel.ori_size);

					if (optimiser.mymodel.data_dim == 2)
						shiftImageInFourierTransform(Frefp, Fref, (RFLOAT)(optimiser.mymodel.ori_size), -XX(offsets), -YY(offsets));
					else
						shiftImageInFourierTransform(Frefp, Fref, (RFLOAT)(optimiser.mymodel.ori_size), -XX(offsets), -YY(offsets), -ZZ(offsets));

					// Take inverse transform
					transformer.inverseFourierTransform(Fref, Mref);
					CenterFFT(Mref, false);
					Mref.setXmippOrigin();

					// Reposition Mref back into the micrograph
					Imic_out().xinit = -ROUND(xcoord);
					Imic_out().yinit = -ROUND(ycoord);
					Imic_out().zinit = -ROUND(zcoord);
					Imic_sum.xinit = -ROUND(xcoord);
					Imic_sum.yinit = -ROUND(ycoord);
					Imic_sum.zinit = -ROUND(zcoord);
					FOR_ALL_ELEMENTS_IN_ARRAY3D(Mref)
					{
						long int idx = ROUND(sqrt(k*k + i*i + j*j));
						if (idx < radius)
						{
							// check the particles do not go off the side
							int kp = (k) - STARTINGZ(Imic_sum);
							int ip = (i) - STARTINGY(Imic_sum);
							int jp = (j) - STARTINGX(Imic_sum);
							if (kp >= 0 && kp < ZSIZE(Imic_sum) && ip >= 0 && ip < YSIZE(Imic_sum) && jp >= 0 && jp < XSIZE(Imic_sum) )
							{
								A3D_ELEM(Imic_out(), k, i, j) += A3D_ELEM(Mref, k, i, j);
								A3D_ELEM(Imic_sum, k, i, j) += 1.;
							}
						}
					}

				}
			} // end loop over all particles in the mydata.MDimg table


			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imic_out())
			{
				if (DIRECT_MULTIDIM_ELEM(Imic_sum, n) > 0.)
					DIRECT_MULTIDIM_ELEM(Imic_out(), n) /= DIRECT_MULTIDIM_ELEM(Imic_sum, n);
				if (do_invert)
					DIRECT_MULTIDIM_ELEM(Imic_out(), n) *= -1.;
				if (micrograph_background > 0.)
				{
					// normalize Imic_in on the fly
					DIRECT_MULTIDIM_ELEM(Imic_in(), n) -= mean_mic;
					DIRECT_MULTIDIM_ELEM(Imic_in(), n) /= stddev_mic;
					// And add a precentage to Imic_out
					DIRECT_MULTIDIM_ELEM(Imic_out(), n) *= (1. - micrograph_background);
					DIRECT_MULTIDIM_ELEM(Imic_out(), n) += micrograph_background * DIRECT_MULTIDIM_ELEM(Imic_in(), n);
				}
			}

			// Write out the new micrograph
			Imic_out.write(fn_mic_out);

			// Also write out a STAR file with the particles used
			MDcoord.write(fn_coord_out);
			MDcoord.clear();

			if (imgno%barstep==0) progress_bar(imgno);
			imgno++;

		} // end loop over input MetadataTable
		progress_bar(DFi.numberOfObjects());

		std::cout << " Done!" << std::endl;

	}// end run function



};

int main(int argc, char *argv[])
{
	time_config();
	particle_reposition_parameters prm;

	try
    {
		prm.read(argc, argv);

		prm.run();
    }

    catch (RelionError XE)
    {
        //prm.usage();
        std::cerr << XE;
        exit(1);
    }

    return 0;

}
