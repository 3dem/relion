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

// LIMITATIONS:
//  This program ignores (anisotropic) magnification and antisymmetric aberrations!

#include <src/projector.h>
#include <src/backprojector.h>
#include <src/fftw.h>
#include <src/args.h>
#include <src/ctf.h>
#include <src/strings.h>
#include <src/funcs.h>
#include <src/memory.h>
#include <src/euler.h>
#include <src/time.h>
#include <src/metadata_table.h>
#include <src/exp_model.h>
#include <src/healpix_sampling.h>
#include <src/jaz/single_particle/obs_model.h>
class project_parameters
{
public:

	FileName fn_map, fn_ang, fn_out, fn_img, fn_model, fn_sym, fn_mask, fn_ang_simulate;
	RFLOAT rot, tilt, psi, xoff, yoff, zoff, angpix, maxres, stddev_white_noise, particle_diameter, ana_prob_range, ana_prob_step, sigma_offset;
	int padding_factor;
	int r_max, r_min_nn, interpolator, nr_uniform;
	bool do_only_one, do_ctf, do_ctf2, ctf_phase_flipped, do_ctf_intact_1st_peak, do_timing, do_add_noise, do_subtract_exp, do_ignore_particle_name, do_3d_rot, write_float16;
	bool do_simulate;
	RFLOAT simulate_SNR;
	// I/O Parser
	IOParser parser;
	ObservationModel obsModel;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("Options");
		fn_map = parser.getOption("--i", "Input map to be projected");
		fn_out = parser.getOption("--o", "Rootname for output projections", "proj");
		write_float16  = parser.checkOption("--float16", "Write in half-precision 16 bit floating point numbers (MRC mode 12), instead of 32 bit (MRC mode 0).");
		do_ctf = parser.checkOption("--ctf", "Apply CTF to reference projections");
		ctf_phase_flipped = parser.checkOption("--ctf_phase_flip", "Flip phases of the CTF in the output projections");
		do_ctf_intact_1st_peak = parser.checkOption("--ctf_intact_first_peak", "Ignore CTFs until their first peak?");
		angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "-1"));
		fn_mask = parser.getOption("--mask", "Mask that will be applied to the input map prior to making projections", "");
		fn_ang = parser.getOption("--ang", "Particle STAR file with orientations and CTF for multiple projections (if None, assume single projection)", "None");
		nr_uniform = textToInteger(parser.getOption("--nr_uniform", " OR get this many random samples from a uniform angular distribution", "-1"));
		sigma_offset = textToFloat(parser.getOption("--sigma_offset", "Apply Gaussian errors (A) with this stddev to the XY-offsets", "0"));
		rot = textToFloat(parser.getOption("--rot", "First Euler angle (for a single projection)", "0"));
		tilt = textToFloat(parser.getOption("--tilt", "Second Euler angle (for a single projection)", "0"));
		psi = textToFloat(parser.getOption("--psi", "Third Euler angle (for a single projection)", "0"));
		xoff = textToFloat(parser.getOption("--xoff", "Origin X-offsets (in pixels) (for a single projection)", "0"));
		yoff = textToFloat(parser.getOption("--yoff", "Origin Y-offsets (in pixels) (for a single projection)", "0"));
		zoff = textToFloat(parser.getOption("--zoff", "Origin Z-offsets (in pixels) (for a single 3D rotation)", "0"));
		do_add_noise = parser.checkOption("--add_noise", "Add noise to the output projections (only with --ang)");
		stddev_white_noise = textToFloat(parser.getOption("--white_noise", "Standard deviation of added white Gaussian noise", "0"));
		fn_model = parser.getOption("--model_noise", "Model STAR file with power spectra for coloured Gaussian noise", "");
		do_subtract_exp = parser.checkOption("--subtract_exp", "Subtract projections from experimental images (in --ang)");
		do_ignore_particle_name = parser.checkOption("--ignore_particle_name", "Ignore the rlnParticleName column (in --ang)");
		do_only_one = (fn_ang == "None" && nr_uniform < 0);
		do_3d_rot = parser.checkOption("--3d_rot", "Perform 3D rotations instead of projection into 2D images");
		do_simulate = parser.checkOption("--simulate", "Simulate data with known ground-truth by subtracting signal and adding projection in random orientation.");
		simulate_SNR = textToFloat(parser.getOption("--adjust_simulation_SNR", "Relative SNR compared to input images for realistic simulation of data", "1."));
		fn_ang_simulate = parser.getOption("--ang_simulate", "STAR file with orientations for projections of realistic simulations (random from --ang STAR file by default)", "");

		maxres = textToFloat(parser.getOption("--maxres", "Maximum resolution (in Angstrom) to consider in Fourier space (default Nyquist)", "-1"));
		padding_factor = textToInteger(parser.getOption("--pad", "Padding factor", "2"));
		do_ctf2 = parser.checkOption("--ctf2", "Apply CTF*CTF to reference projections");
		if (parser.checkOption("--NN", "Use nearest-neighbour instead of linear interpolation"))
			interpolator = NEAREST_NEIGHBOUR;
		else
			interpolator = TRILINEAR;

		// Hidden
		r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));

		if (do_simulate)
		{
			if (!do_ctf)
			{
				std::cerr << "WARNING: with --simulate, --ctf is automatically activated." << std::endl;
				do_ctf = true;
			}
		}

		// Check for errors in the command-line option
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	void project()
	{
		int ori_size, nr_groups;
		std::vector<FileName> group_names;
		std::vector<MultidimArray<RFLOAT > > sigma2_noise;

		MetaDataTable DFo, MDang, MDang_sim;
		Matrix2D<RFLOAT> A3D;
		FileName fn_expimg;

		MultidimArray<Complex > F3D, F2D, Fexpimg;
		MultidimArray<RFLOAT> Fctf, dummy;
		Image<RFLOAT> vol, img, expimg;
		FourierTransformer transformer, transformer_expimg;

		std::cout << " Reading map: " << fn_map << std::endl;
		vol.read(fn_map);
		std::cout << " Done reading map!" << std::endl;

		if (fn_mask != "")
		{
			Image<RFLOAT> msk;
			msk.read(fn_mask);
			if (!msk().sameShape(vol()))
				REPORT_ERROR("project ERROR: mask and map have different sizes!");
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(vol())
			DIRECT_MULTIDIM_ELEM(vol(), n) *= DIRECT_MULTIDIM_ELEM(msk(), n);
		}

		if (nr_uniform > 0)
		{
			std::cout << " Generating " << nr_uniform << " projections taken randomly from a uniform angular distribution ..." << std::endl;
			MDang.clear();
			randomize_random_generator();
			for (long int i = 0; i < nr_uniform; i++)
			{
				RFLOAT rot, tilt, psi, xoff, yoff;
				rot = rnd_unif() * 360.;
				bool ok_tilt = false;
				while (!ok_tilt)
				{
					tilt = rnd_unif() * 180.;
					if (rnd_unif() < fabs(SIND(tilt)))
						ok_tilt = true;
				}
				psi = rnd_unif() * 360.;
				xoff = rnd_gaus(0., sigma_offset);
				yoff = rnd_gaus(0., sigma_offset);
				MDang.addObject();
				MDang.setValue(EMDL_ORIENT_ROT, rot);
				MDang.setValue(EMDL_ORIENT_TILT, tilt);
				MDang.setValue(EMDL_ORIENT_PSI, psi);
				MDang.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff);
				MDang.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff);
				MDang.setValue(EMDL_IMAGE_OPTICS_GROUP, 1);
			}

			std::cout << " Setting default values for optics table, though CTFs are not used in the projections ... " << std::endl;
			MetaDataTable MDopt;
			MDopt.addObject();
			MDopt.setValue(EMDL_IMAGE_OPTICS_GROUP, 1);
			std::string name = "optics1";
			MDopt.setValue(EMDL_IMAGE_OPTICS_GROUP_NAME, name);
			MDopt.setValue(EMDL_CTF_VOLTAGE, 300.);
			MDopt.setValue(EMDL_CTF_CS, 2.7);
			vol.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_X, angpix);
			MDopt.setValue(EMDL_IMAGE_PIXEL_SIZE, angpix);
			MDopt.setValue(EMDL_IMAGE_SIZE, XSIZE(vol()));
			int mydim = (do_3d_rot) ? 3 : 2;
			MDopt.setValue(EMDL_IMAGE_DIMENSIONALITY, mydim);

			obsModel = ObservationModel(MDopt);
		}
		else if (!do_only_one)
		{
			std::cout << " Reading STAR file with all angles " << fn_ang << std::endl;
			if (do_ignore_particle_name)
            {
                MDang.read(fn_ang);
            }
            else
            {
                ObservationModel::loadSafely(fn_ang, obsModel, MDang);
            }
			std::cout << " Done reading STAR file!" << std::endl;


			if (do_simulate && fn_ang_simulate != "")
			{
				std::cout << " Reading STAR file with angles for simulated images " << fn_ang << std::endl;
				MDang_sim.read(fn_ang_simulate);
				std::cout << " Done reading STAR file with angles for simulated images!" << std::endl;
				if (MDang_sim.numberOfObjects() < MDang.numberOfObjects())
				{
					REPORT_ERROR("ERROR: STAR file with angles for simulated images has fewer entries than the input STAR file with all angles.");
				}
			}
		}

		if (angpix < 0.)
		{
			if (!do_only_one)
			{
				// Get angpix from the first optics group in the obsModel
				angpix = obsModel.getPixelSize(0);
				std::cout << " + Using pixel size from the first optics group in the --ang STAR file: " << angpix << std::endl;
			}
			else
			{
				angpix = vol.samplingRateX();
				std::cerr << "WARNING: The pixel size (--angpix) was not specified." << std::endl;
				std::cerr << "         The value in the input image header (= " << angpix << ") is used instead." << std::endl;
			}
		}

		// Now that we have the size of the volume, check r_max
		if (maxres < 0.)
			r_max = XSIZE(vol());
		else
			r_max = CEIL(XSIZE(vol()) * angpix / maxres);

		// Set right size of F2D and initialize to zero
		if (do_3d_rot)
			img().resize(ZSIZE(vol()), YSIZE(vol()), XSIZE(vol()));
		else
			img().resize(YSIZE(vol()), XSIZE(vol()));
		transformer.setReal(img());
		transformer.getFourierAlias(F2D);

		// Set up the projector
		int data_dim = (do_3d_rot) ? 3 : 2;
		Projector projector((int)XSIZE(vol()), interpolator, padding_factor, r_min_nn, data_dim);
		projector.computeFourierTransformMap(vol(), dummy, 2* r_max);

		if (do_only_one)
		{
			Euler_rotation3DMatrix(rot, tilt, psi, A3D);
			F2D.initZeros();
			projector.get2DFourierTransform(F2D, A3D);
			if (ABS(xoff) > 0.001 || ABS(yoff) > 0.001 || (do_3d_rot && ABS(zoff) > 0.001) )
			{
				Matrix1D<RFLOAT> shift(2);
				XX(shift) = -xoff;
				YY(shift) = -yoff;
				if (do_3d_rot)
				{
					shift.resize(3);
					ZZ(shift) = -zoff;
					shiftImageInFourierTransform(F2D, F2D, XSIZE(vol()), XX(shift), YY(shift), ZZ(shift));
				}
				else
					shiftImageInFourierTransform(F2D, F2D, XSIZE(vol()), XX(shift), YY(shift));
			}

			// Feb 01,2017 - Shaoda, add white noise to 2D / 3D single images
			if (do_add_noise)
			{
				if ( (!(stddev_white_noise > 0.)) || (fn_model != "") )
					REPORT_ERROR("ERROR: Only add --white_noise to a single image!");
				// fftw normalization and factor sqrt(2) for two-dimensionality of complex plane
				// TODO: sqrt(2) ??? Why ???
				stddev_white_noise /= (data_dim == 3) ? (XSIZE(vol()) * XSIZE(vol())) : (XSIZE(vol()) * sqrt(2));
				// Add white noise
				FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D)
				{
					DIRECT_A3D_ELEM(F2D, k, i, j).real += rnd_gaus(0., stddev_white_noise);
					DIRECT_A3D_ELEM(F2D, k, i, j).imag += rnd_gaus(0., stddev_white_noise);
				}
			}

			transformer.inverseFourierTransform(F2D, img());
			// Shift the image back to the center...
			CenterFFT(img(), false);
			img.setSamplingRateInHeader(angpix);
			img.write(fn_out, -1, false, WRITE_OVERWRITE, write_float16 ? Float16: Float);
			std::cout<<" Done writing "<<fn_out<<std::endl;
		}
		else // not do_only_one
		{
			init_progress_bar(MDang.numberOfObjects());
			DFo.clear();
			rot = tilt = psi = xoff = yoff = zoff = 0.;

			// Can only add noise to multiple images
			// Feb 01,2017 - Shaoda, now we can add white noise to 2D / 3D single images
			if (do_add_noise)
			{
				if (fn_model != "")
				{
					std::ifstream in(fn_model.data(), std::ios_base::in);
					if (in.fail())
						REPORT_ERROR( (std::string) "MlModel::readStar: File " + fn_model + " cannot be read." );

					MetaDataTable MDlog, MDgroup, MDsigma;
					MDlog.readStar(in, "model_general");
					if (!MDlog.getValue(EMDL_MLMODEL_ORIGINAL_SIZE, ori_size) ||
					    !MDlog.getValue(EMDL_MLMODEL_NR_GROUPS, nr_groups) )
						REPORT_ERROR("MlModel::readStar: incorrect model_general table");

					MDgroup.readStar(in, "model_groups");
					group_names.resize(nr_groups, "");
					FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDgroup)
					{
						long int igroup;
						if (!MDgroup.getValue(EMDL_MLMODEL_GROUP_NO, igroup))
							REPORT_ERROR("MlModel::readStar: incorrect model_groups table");
						//Start counting of groups at 1, not at 0....
						if (!MDgroup.getValue(EMDL_MLMODEL_GROUP_NAME, group_names[igroup-1]))
							REPORT_ERROR("MlModel::readStar: incorrect model_groups table");
					}

					sigma2_noise.resize(nr_groups);
					for (int igroup = 0; igroup < nr_groups; igroup++)
					{
						sigma2_noise[igroup].resize(ori_size/2 + 1);
						MDsigma.readStar(in, "model_group_" + integerToString(igroup + 1));
						int idx;
						FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsigma)
						{
							if (!MDsigma.getValue(EMDL_SPECTRAL_IDX, idx))
								REPORT_ERROR("MlModel::readStar: incorrect table model_group_" + integerToString(igroup + 1));
							if (!MDsigma.getValue(EMDL_MLMODEL_SIGMA2_NOISE, sigma2_noise[igroup](idx)))
								REPORT_ERROR("MlModel::readStar: incorrect table model_group_" + integerToString(igroup + 1));
						}
					}
				}
				else if (stddev_white_noise > 0.)
					stddev_white_noise /= XSIZE(vol()) * sqrt(2); // fftw normalization and factor sqrt(2) for two-dimensionality of complex plane
				else
					REPORT_ERROR("ERROR: When adding noise provide either --model_noise or --white_noise");
			}

			long int imgno = 0;
			long int max_imgno = MDang.numberOfObjects() - 1;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDang)
			{
				MDang.getValue(EMDL_ORIENT_ROT, rot);
				MDang.getValue(EMDL_ORIENT_TILT, tilt);
				MDang.getValue(EMDL_ORIENT_PSI, psi);
				MDang.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff);
				MDang.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff);
				if (do_3d_rot)
					MDang.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, zoff);

				xoff /= angpix;
				yoff /= angpix;
				zoff /= angpix;

				Euler_rotation3DMatrix(rot, tilt, psi, A3D);
				F2D.initZeros();
				projector.get2DFourierTransform(F2D, A3D);

				if (ABS(xoff) > 0.001 || ABS(yoff) > 0.001 || (do_3d_rot && ABS(zoff) > 0.001) )
				{
					Matrix1D<RFLOAT> shift(2);
					XX(shift) = -xoff;
					YY(shift) = -yoff;

					if (do_3d_rot)
					{
						shift.resize(3);
						ZZ(shift) = -zoff;
						shiftImageInFourierTransform(F2D, F2D, XSIZE(vol()), XX(shift), YY(shift), ZZ(shift) );
					}
					else
						shiftImageInFourierTransform(F2D, F2D, XSIZE(vol()), XX(shift), YY(shift) );
				}

				// Apply CTF if necessary
				CTF ctf;
				if (do_ctf || do_ctf2)
				{
					if (do_3d_rot)
					{
						Image<RFLOAT> Ictf;
						FileName fn_ctf;
						MDang.getValue(EMDL_CTF_IMAGE, fn_ctf);
						Ictf.read(fn_ctf);

						// Set the CTF-image in Fctf
						Fctf.resize(F2D);

						// If there is a redundant half, get rid of it
						if (XSIZE(Ictf()) == YSIZE(Ictf()))
						{
							Ictf().setXmippOrigin();
							// Set the CTF-image in Fctf
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
						ctf.readByGroup(MDang, &obsModel); // This MDimg only contains one particle!
						Fctf.resize(F2D);
						ctf.getFftwImage(Fctf, XSIZE(vol()), XSIZE(vol()), angpix, ctf_phase_flipped, false,  do_ctf_intact_1st_peak, true);
					}
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
					{
						DIRECT_MULTIDIM_ELEM(F2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
						if (do_ctf2)
							DIRECT_MULTIDIM_ELEM(F2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
				}

				// Apply Gaussian noise
				if (do_add_noise)
				{
					if (fn_model !="")
					{
						//// 23MAY2014: for preparation of 1.3 release: removed reading a exp_model, replaced by just reading MDang
						// This does however mean that I no longer know mic_id of this image: replace by 0....
						FileName fn_group;
						if (MDang.containsLabel(EMDL_MLMODEL_GROUP_NAME))
						{
							MDang.getValue(EMDL_MLMODEL_GROUP_NAME, fn_group);
						}
						else
						{
							if (MDang.containsLabel(EMDL_MICROGRAPH_NAME))
							{
								FileName fn_orig, fn_pre, fn_jobnr;
								MDang.getValue(EMDL_MICROGRAPH_NAME, fn_orig);
								if (!decomposePipelineFileName(fn_orig, fn_pre, fn_jobnr, fn_group)) {
									fn_group = fn_orig; // Not a pipeline filename; use as is
								}
							}
							else
							{
								REPORT_ERROR("ERROR: cannot find rlnGroupName or rlnMicrographName in the input --ang file...");
							}
						}
						int my_mic_id = -1;
						for (int mic_id = 0; mic_id < group_names.size(); mic_id++)
						{
							if (fn_group == group_names[mic_id])
							{
								my_mic_id = mic_id;
								break;
							}
						}
						if (my_mic_id < 0)
							REPORT_ERROR("ERROR: cannot find " + fn_group + " in the input model file...");

						RFLOAT normcorr = 1.;
						if (MDang.containsLabel(EMDL_IMAGE_NORM_CORRECTION))
						{
							MDang.getValue(EMDL_IMAGE_NORM_CORRECTION, normcorr);
						}

						// Add coloured noise
						FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D)
						{
							int ires = ROUND( sqrt( (RFLOAT)(kp*kp + ip*ip + jp*jp) ) );
							ires = XMIPP_MIN(ires, ori_size/2); // at freqs higher than Nyquist: use last sigma2 value

							RFLOAT sigma = sqrt(DIRECT_A1D_ELEM(sigma2_noise[my_mic_id], ires));
							DIRECT_A3D_ELEM(F2D, k, i, j).real += rnd_gaus(0., sigma);
							DIRECT_A3D_ELEM(F2D, k, i, j).imag += rnd_gaus(0., sigma);
						}
					}
					else
					{
						// Add white noise
					        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D)
						{
							DIRECT_A3D_ELEM(F2D, k, i, j).real += rnd_gaus(0., stddev_white_noise);
							DIRECT_A3D_ELEM(F2D, k, i, j).imag += rnd_gaus(0., stddev_white_noise);
						}
					}
				}

				img().initZeros();
				transformer.inverseFourierTransform(F2D, img());
				// Shift the image back to the center...
				CenterFFT(img(), false);

				// Subtract the projection from the corresponding experimental image
				if (do_subtract_exp || do_simulate)
				{
					MDang.getValue(EMDL_IMAGE_NAME, fn_expimg);
					MDang.setValue(EMDL_IMAGE_ORI_NAME, fn_expimg); // Store fn_expimg in rlnOriginalParticleName
					expimg.read(fn_expimg);
					img() = expimg() - img();
				}

				// If we're simulating realistic images, then now add CTF-affected projection again
				if (do_simulate)
				{
					// Take random orientation from the input STAR file is fn_ang_simulate is empty. Otherwise, use fn_ang_simulate
					if (fn_ang_simulate == "")
					{
						long int random_imgno = -1;
						while (random_imgno < 0 || random_imgno > max_imgno)
						{
							random_imgno = rnd_unif()*max_imgno;
						}

						MDang.getValue(EMDL_ORIENT_ROT, rot, random_imgno);
						MDang.getValue(EMDL_ORIENT_TILT, tilt, random_imgno);
						MDang.getValue(EMDL_ORIENT_PSI, psi, random_imgno);
						MDang.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff, random_imgno);
						MDang.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff, random_imgno);
						if (do_3d_rot)
							MDang.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, zoff, random_imgno);

						xoff /= angpix;
						yoff /= angpix;
						zoff /= angpix;
					}
					else
					{
						MDang_sim.getValue(EMDL_ORIENT_ROT, rot, imgno);
						MDang_sim.getValue(EMDL_ORIENT_TILT, tilt, imgno);
						MDang_sim.getValue(EMDL_ORIENT_PSI, psi, imgno);
						MDang_sim.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff, imgno);
						MDang_sim.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff, imgno);
						if (do_3d_rot)
							MDang_sim.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, zoff, imgno);

						xoff /= angpix;
						yoff /= angpix;
						zoff /= angpix;
					}

					Euler_rotation3DMatrix(rot, tilt, psi, A3D);
					F2D.initZeros();
					projector.get2DFourierTransform(F2D, A3D);

					if (ABS(xoff) > 0.001 || ABS(yoff) > 0.001 || (do_3d_rot && ABS(zoff) > 0.001) )
					{
						Matrix1D<RFLOAT> shift(2);
						XX(shift) = -xoff;
						YY(shift) = -yoff;

						if (do_3d_rot)
						{
							shift.resize(3);
							ZZ(shift) = -zoff;
							shiftImageInFourierTransform(F2D, F2D, XSIZE(vol()), XX(shift), YY(shift), ZZ(shift) );
						}
						else
							shiftImageInFourierTransform(F2D, F2D, XSIZE(vol()), XX(shift), YY(shift) );
					}

					// Apply CTF
					CTF ctf;
					if (do_ctf || do_ctf2)
					{
						if (do_3d_rot)
						{
							Image<RFLOAT> Ictf;
							FileName fn_ctf;
							MDang.getValue(EMDL_CTF_IMAGE, fn_ctf);
							Ictf.read(fn_ctf);
							Ictf().setXmippOrigin();

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
							ctf.read(MDang, MDang, imgno);
							Fctf.resize(F2D);
							ctf.getFftwImage(Fctf, XSIZE(vol()), XSIZE(vol()), angpix, ctf_phase_flipped, false,  do_ctf_intact_1st_peak, true);
						}
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
						{
							DIRECT_MULTIDIM_ELEM(F2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
							if (do_ctf2)
								DIRECT_MULTIDIM_ELEM(F2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
						}
					}

					expimg().initZeros();
					transformer.inverseFourierTransform(F2D, expimg());
					// Shift the image back to the center...
					CenterFFT(expimg(), false);

					// Modify the strength of the signal
					if (fabs(simulate_SNR - 1.) > 0.000001)
					{
						expimg() *= simulate_SNR;
					}

					img() += expimg();

				}

				img.setSamplingRateInHeader(angpix);
				if (do_3d_rot)
				{
					fn_img.compose(fn_out, imgno+1,"mrc");
					img.write(fn_img, -1, false, WRITE_OVERWRITE, write_float16 ? Float16: Float);
				}
				else
				{
					// Write this particle to the stack on disc
					// First particle: write stack in overwrite mode, from then on just append to it
					fn_img.compose(imgno+1,fn_out+".mrcs");
					if (imgno == 0)
						img.write(fn_img, -1, false, WRITE_OVERWRITE, write_float16 ? Float16: Float);
					else
						img.write(fn_img, -1, false, WRITE_APPEND, write_float16 ? Float16: Float);
				}

				// Set the image name to the output STAR file
				DFo.addObject();
				DFo.setObject(MDang.getObject());
				DFo.setValue(EMDL_IMAGE_NAME,fn_img);

				if (do_simulate)
				{
					DFo.setValue(EMDL_ORIENT_ROT, rot);
					DFo.setValue(EMDL_ORIENT_TILT, tilt);
					DFo.setValue(EMDL_ORIENT_PSI, psi);
					DFo.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff * angpix);
					DFo.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff * angpix);
					if (do_3d_rot)
						DFo.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, zoff * angpix);
				}

				if (imgno%60==0) progress_bar(imgno);
				imgno++;
			}
			progress_bar(MDang.numberOfObjects());

			// Write out STAR file with all information
			fn_img = fn_out + ".star";
			obsModel.save(DFo, fn_img);
			std::cout<<" Done writing "<<imgno<<" images in "<<fn_img<<std::endl;

		} // end else do_only_one
	}// end project function
};

int main(int argc, char *argv[])
{
	time_config();
	project_parameters prm;

	try
	{
		prm.read(argc, argv);
		prm.project();
	}
	catch (RelionError XE)
	{
	        //prm.usage();
        	std::cerr << XE;
	        return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}
