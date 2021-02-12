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
#include <src/jaz/single_particle/obs_model.h>
#include <stdlib.h>

class particle_reposition_parameters
{
public:

	FileName fn_in, fn_opt, fn_out, fn_dat, fn_odir;

	RFLOAT micrograph_background;
	int norm_radius;
	bool do_invert, do_ctf, do_subtract;
    ObservationModel obsModelMics;

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
		fn_out = parser.getOption("--o", "Output rootname, to be added to input micrograph names", "");
		fn_odir = parser.getOption("--odir", "Output directory (default is same as input micrographs directory", "");
		fn_opt = parser.getOption("--opt", "Optimiser STAR file with the 2D classes or 3D maps to be repositioned");
		fn_dat = parser.getOption("--data", "Data STAR file with selected particles (default is to use all particles)", "");
		micrograph_background = textToFloat(parser.getOption("--background", "The fraction of micrograph background noise in the output micrograph", "0.1"));
		do_invert= parser.checkOption("--invert", "Invert the contrast in the references?");

		do_ctf = parser.checkOption("--ctf", "Apply CTF for each particle to the references?");
		norm_radius = textToFloat(parser.getOption("--norm_radius", "Radius of the circle used for background normalisation (in pixels)", "-1"));
		do_subtract = parser.checkOption("--subtract", "Subtract repositioned micrographs from the input ones?");

		// Check for errors in the command-line option
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	void run()
	{

		if (fn_out == "" && fn_odir == "")
			REPORT_ERROR("ERROR: You need to provide either --o or --odir");

		if (fn_odir.length() > 0 && fn_odir[fn_odir.length()-1] != '/') fn_odir += "/";

		int xdim, ydim, radius;
		MetaDataTable DFi, DFopt, MDmics_out;
		ObservationModel::loadSafely(fn_in, obsModelMics, DFi, "micrographs");

		MlOptimiser optimiser;
		optimiser.do_preread_images = false;

		optimiser.read(fn_opt);
		optimiser.mymodel.setFourierTransformMaps(false);

		// Use a user-provided subset of particles instead of all of them?
		if (fn_dat != "")
		{
            std::cout <<" Reading data ..." << std::endl;
			MetaDataTable MDdata;
			MDdata.read(fn_dat);
			optimiser.mydata.MDimg = MDdata;
		}


		// Loop over all micrographs
		int barstep = XMIPP_MAX(1, DFi.numberOfObjects()/ 60);
		init_progress_bar(DFi.numberOfObjects());
		long int imgno = 0;
		FileName fn_prevdir="";
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(DFi)
		{

			FileName fn_mic, fn_mic_out;
			DFi.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
			if (fn_out != "") fn_mic_out = fn_mic.insertBeforeExtension("_" + fn_out);
			else fn_mic_out = fn_mic;
			if (fn_odir != "")
			{
				FileName fn_pre, fn_jobnr, fn_post;
				if (decomposePipelineFileName(fn_mic_out, fn_pre, fn_jobnr, fn_post))
					fn_mic_out = fn_odir + fn_post;
				else
					fn_mic_out = fn_odir + fn_mic_out;
				FileName fn_onlydir = fn_mic_out.beforeLastOf("/");
				if (fn_onlydir != fn_prevdir)
				{
					mktree(fn_onlydir);
					fn_prevdir = fn_onlydir;
				}
			}

			FourierTransformer transformer;
			MetaDataTable MDcoord;
			Image<RFLOAT> Imic_in, Imic_out;
			MultidimArray<RFLOAT> Imic_sum;
			// Read in the first micrograph
			Imic_in.read(fn_mic);
			Imic_in().setXmippOrigin();
			Imic_out().initZeros(Imic_in());
			Imic_sum.initZeros(Imic_in());
			Imic_sum.setXmippOrigin();
			// Get mean and stddev of the input micrograph
			RFLOAT stddev_mic, mean_mic, dummy;
			Imic_in().computeStats(mean_mic, stddev_mic, dummy, dummy);

			int optics_group_mic;
			DFi.getValue(EMDL_IMAGE_OPTICS_GROUP, optics_group_mic);
			RFLOAT mic_pixel_size=-1.;
			for (int i = 0; i < obsModelMics.opticsMdt.numberOfObjects(); i++)
			{
				int my_optics_group;
				obsModelMics.opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, my_optics_group);
				if (my_optics_group == optics_group_mic)
				{
					obsModelMics.opticsMdt.getValue(EMDL_MICROGRAPH_PIXEL_SIZE, mic_pixel_size);
					break;
				}
			}
			if (mic_pixel_size<0.)
				REPORT_ERROR("ERROR: could not find correct optics group in micrograph star file...");

			FileName fn_mic_pre, fn_mic_jobnr, fn_mic_post;
			decomposePipelineFileName(fn_mic, fn_mic_pre, fn_mic_jobnr, fn_mic_post);

			// Loop over all particles
			bool found_one = false;
			for (long int part_id = 0; part_id < optimiser.mydata.numberOfParticles(); part_id++)
			{
				long int ori_img_id = optimiser.mydata.particles[part_id].images[0].id;
				int optics_group = optimiser.mydata.getOpticsGroup(part_id, 0);
				RFLOAT my_pixel_size = optimiser.mydata.getImagePixelSize(part_id, 0);
				int my_image_size = optimiser.mydata.getOpticsImageSize(optics_group);

				if (do_subtract && fabs(my_pixel_size - mic_pixel_size) > 1e-6)
					REPORT_ERROR("ERROR: subtract code has only been validated with same pixel size for particles and micrographs... Sorry!");

				FileName fn_mic2;
				optimiser.mydata.MDimg.getValue(EMDL_MICROGRAPH_NAME, fn_mic2, ori_img_id);
				FileName fn_mic2_pre, fn_mic2_jobnr, fn_mic2_post;
				decomposePipelineFileName(fn_mic2, fn_mic2_pre, fn_mic2_jobnr, fn_mic2_post);

				if (fn_mic2_post == fn_mic_post)
				{

					found_one = true;

					// Prepare transformer
					MultidimArray<Complex > Fref;
					MultidimArray<RFLOAT> Mref;
					if (optimiser.mymodel.data_dim == 3)
					{
						Mref.resize(my_image_size, my_image_size, my_image_size);
						Fref.resize(my_image_size, my_image_size, my_image_size/2 + 1);
					}
					else
					{
						Mref.resize(my_image_size, my_image_size);
						Fref.resize(my_image_size, my_image_size/2 + 1);
					}


					RFLOAT rot=0., tilt=0., psi, xcoord=0., ycoord=0., zcoord=0.;
					int iclass;
					Matrix2D<RFLOAT> A;
					Matrix1D<RFLOAT> offsets(3);


					MDcoord.addObject();
					MDcoord.setObject(optimiser.mydata.MDimg.getObject(ori_img_id));
					MDcoord.setValue(EMDL_MICROGRAPH_NAME,fn_mic_out);

					optimiser.mydata.MDimg.getValue(EMDL_IMAGE_COORD_X,  xcoord, ori_img_id);
					optimiser.mydata.MDimg.getValue(EMDL_IMAGE_COORD_Y,  ycoord, ori_img_id);
					if (optimiser.mymodel.ref_dim == 3)
					{
						optimiser.mydata.MDimg.getValue(EMDL_ORIENT_ROT,  rot, ori_img_id);
						optimiser.mydata.MDimg.getValue(EMDL_ORIENT_TILT, tilt, ori_img_id);
					}
					optimiser.mydata.MDimg.getValue(EMDL_ORIENT_PSI,  psi, ori_img_id);
					optimiser.mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, XX(offsets), ori_img_id);
					optimiser.mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, YY(offsets), ori_img_id);
					if (optimiser.mymodel.data_dim == 3)
					{
						optimiser.mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, ZZ(offsets), ori_img_id);
						optimiser.mydata.MDimg.getValue(EMDL_IMAGE_COORD_Z,  zcoord, ori_img_id);
					}
					else
					{
						ZZ(offsets) = zcoord = 0.;
					}

					// Offsets in pixels
					offsets /= my_pixel_size;

					optimiser.mydata.MDimg.getValue(EMDL_PARTICLE_CLASS, iclass, ori_img_id);
					iclass--;

					Euler_angles2matrix(rot, tilt, psi, A);
					if (do_ctf)
					{
						A = optimiser.mydata.obsModel.applyAnisoMag(A, optics_group);
						A = optimiser.mydata.obsModel.applyScaleDifference(A, optics_group, optimiser.mymodel.ori_size, optimiser.mymodel.pixel_size);
					}

					// Get the 2D image (in its ori_size)
					(optimiser.mymodel.PPref[iclass]).get2DFourierTransform(Fref, A);

					if (optimiser.mymodel.data_dim == 2)
						shiftImageInFourierTransform(Fref, Fref, my_image_size, -XX(offsets), -YY(offsets));
					else
						shiftImageInFourierTransform(Fref, Fref, my_image_size, -XX(offsets), -YY(offsets), -ZZ(offsets));


					if (do_ctf)
					{
						MultidimArray<RFLOAT> Fctf;
						Fctf.resize(Fref);

						CTF ctf;
						if (optimiser.mymodel.data_dim == 3)
						{
							Image<RFLOAT> Ictf;

							FileName fn_ctf;
							optimiser.mydata.MDimg.getValue(EMDL_CTF_IMAGE, fn_ctf, ori_img_id);
							Ictf.read(fn_ctf);

							// If there is a redundant half, get rid of it
							if (XSIZE(Ictf()) == YSIZE(Ictf()))
							{
								Ictf().setXmippOrigin();
								FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf)
								{
									// Use negative kp,ip and jp indices, because the origin in the ctf_img lies half a pixel to the right of the actual center....
									DIRECT_A3D_ELEM(Fctf, k, i, j) = A3D_ELEM(Ictf(), -kp, -ip, -jp);
								}
							}
							// otherwise, just window the CTF to the current resolution
							else if (XSIZE(Ictf()) == YSIZE(Ictf()) / 2 + 1)
							{
								windowFourierTransform(Ictf(), Fctf, YSIZE(Fctf));
							}
							// if dimensions are neither cubical nor FFTW, stop
							else
							{
								REPORT_ERROR("3D CTF volume must be either cubical or adhere to FFTW format!");
							}
						}
						else
						{
							ctf.readByGroup(optimiser.mydata.MDimg, &optimiser.mydata.obsModel, ori_img_id);
							ctf.getFftwImage(Fctf, my_image_size, my_image_size, my_pixel_size,
									optimiser.ctf_phase_flipped, false, optimiser.intact_ctf_first_peak, true);
						}

						if (optimiser.mydata.obsModel.getCtfPremultiplied(optics_group))
						{
							FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
								DIRECT_MULTIDIM_ELEM(Fref, n) *= (DIRECT_MULTIDIM_ELEM(Fctf, n) * DIRECT_MULTIDIM_ELEM(Fctf, n));
						}
						else
						{
							FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
							{
								DIRECT_MULTIDIM_ELEM(Fref, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
							}
						}

						// Also do phase modulation, for beam tilt correction and other asymmetric aberrations
						optimiser.mydata.obsModel.demodulatePhase(optics_group, Fref, true); // true means do_modulate_instead
						optimiser.mydata.obsModel.divideByMtf(optics_group, Fref, true); // true means do_multiply_instead

					} // end if do_ctf

					if (optimiser.do_scale_correction)
					{
						int group_id = optimiser.mydata.getGroupId(part_id, 0);
						RFLOAT myscale = optimiser.mymodel.scale_correction[group_id];
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
						{
							DIRECT_MULTIDIM_ELEM(Fref, n) *= myscale;
						}
					}

					// Take inverse transform
					transformer.inverseFourierTransform(Fref, Mref);
					CenterFFT(Mref, false);
					Mref.setXmippOrigin();

					int mic_image_size = CEIL(my_image_size * my_pixel_size / mic_pixel_size);
					MultidimArray<RFLOAT> Mpart_mic = Mref;
					if (mic_image_size != my_image_size)
					{
						resizeMap(Mpart_mic, mic_image_size);
						Mpart_mic.setXmippOrigin();
					}

					// To keep raw micrograph and reference projections on the same scale, need to re-obtain
					// the multiplicative normalisation of the background area (outside circle) again

					RFLOAT norm_factor = 1.;
					if (norm_radius > 0)
					{
						Image<RFLOAT> Ipart;
						Ipart().resize(Mpart_mic);
						Ipart().initConstant(mean_mic); // set areas outside the micrograph to average of micrograph (just like in preprocessing)
						Imic_in().xinit = -ROUND(xcoord);
						Imic_in().yinit = -ROUND(ycoord);
						Imic_in().zinit = -ROUND(zcoord);
						FOR_ALL_ELEMENTS_IN_ARRAY3D(Mpart_mic)
						{
							// check the particles do not go off the side
							int kp = (k) - STARTINGZ(Imic_in());
							int ip = (i) - STARTINGY(Imic_in());
							int jp = (j) - STARTINGX(Imic_in());
							if (kp >= 0 && kp < ZSIZE(Imic_in()) && ip >= 0 && ip < YSIZE(Imic_in()) && jp >= 0 && jp < XSIZE(Imic_in()) )
							{
								A3D_ELEM(Ipart(), k, i, j) = A3D_ELEM(Imic_in(), k, i, j);
							}
						}

						RFLOAT psi_deg = 0., tilt_deg = 90.;
						RFLOAT part_avg, part_stdev;
						if (optimiser.do_helical_refine)
						{
							optimiser.mydata.MDimg.getValue(EMDL_ORIENT_TILT_PRIOR, tilt_deg, ori_img_id);
							optimiser.mydata.MDimg.getValue(EMDL_ORIENT_PSI_PRIOR, psi_deg, ori_img_id);
						}

						calculateBackgroundAvgStddev(Ipart, part_avg, norm_factor, norm_radius, optimiser.do_helical_refine,
								optimiser.helical_tube_outer_diameter/(2.*mic_pixel_size), tilt_deg, psi_deg);

						// Apply the per-particle norm_correction term
						if (optimiser.do_norm_correction)
						{
							RFLOAT mynorm;
							optimiser.mydata.MDimg.getValue(EMDL_IMAGE_NORM_CORRECTION, mynorm, ori_img_id);
							// TODO: check whether this is the right way around!!!
							norm_factor *= mynorm/optimiser.mymodel.avg_norm_correction;
						}
					}

					// Reposition Mpart_mic back into the micrograph
					Imic_out().xinit = -ROUND(xcoord);
					Imic_out().yinit = -ROUND(ycoord);
					Imic_out().zinit = -ROUND(zcoord);
					Imic_sum.xinit = -ROUND(xcoord);
					Imic_sum.yinit = -ROUND(ycoord);
					Imic_sum.zinit = -ROUND(zcoord);
					radius = optimiser.particle_diameter / (2. * mic_pixel_size);
					FOR_ALL_ELEMENTS_IN_ARRAY3D(Mpart_mic)
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
								A3D_ELEM(Imic_out(), k, i, j) += norm_factor * A3D_ELEM(Mpart_mic, k, i, j);
								A3D_ELEM(Imic_sum, k, i, j) += 1.;
							}
						}
					}

				}
			} // end loop over all particles in the mydata.MDimg table


			if (found_one)
			{
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imic_out())
				{
					if (DIRECT_MULTIDIM_ELEM(Imic_sum, n) > 0.)
						DIRECT_MULTIDIM_ELEM(Imic_out(), n) /= DIRECT_MULTIDIM_ELEM(Imic_sum, n);
					if (do_invert)
						DIRECT_MULTIDIM_ELEM(Imic_out(), n) *= -1.;
					if (do_subtract)
					{
						DIRECT_MULTIDIM_ELEM(Imic_out(), n) = DIRECT_MULTIDIM_ELEM(Imic_in(), n) - DIRECT_MULTIDIM_ELEM(Imic_out(), n);
					}
					else if (micrograph_background > 0.)
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

				MDmics_out.addObject();
				MDmics_out.setObject(DFi.getObject());
				MDmics_out.setValue(EMDL_MICROGRAPH_NAME, fn_mic_out);

				// Also write out a STAR file with the particles used
				FileName fn_coord_out = fn_mic_out.withoutExtension()+ "_coord.star";
				MDcoord.write(fn_coord_out);
				MDcoord.clear();
			}
			else
			{
				MDmics_out.addObject();
				MDmics_out.setObject(DFi.getObject());
			}

			if (imgno%barstep==0) progress_bar(imgno);
			imgno++;

		} // end loop over input MetadataTable
		progress_bar(DFi.numberOfObjects());


		FileName fn_star_out = fn_odir + "micrographs_reposition.star";
		if (fn_out != "") fn_star_out = fn_star_out.insertBeforeExtension("_" + fn_out);
		std::cout << "Writing out star file with the new micrographs: " << fn_star_out << std::endl;
		obsModelMics.save(MDmics_out, fn_star_out, "micrographs");

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
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}
