/***************************************************************************
 *
 * Authors: Sjors H.W. Scheres and Jasenko Zivanov
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

#include <src/backprojector.h>
#include <src/funcs.h>
#include <src/ctf.h>
#include <src/args.h>
#include <src/error.h>
#include <src/euler.h>
#include <src/time.h>
#include <omp.h>

#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/complex_io.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/img_proc/image_op.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/new_ft.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/ctf/delocalisation_helper.h>

class reconstruct_parameters
{
	public:

		FileName fn_out, fn_sel, fn_img, fn_sym, fn_sub;

		int r_max, r_min_nn, blob_order, ref_dim, interpolator, grid_iters,
			nr_omp_threads,
			nr_helical_asu, newbox, width_mask_edge, nr_sectors;

		RFLOAT blob_radius, blob_alpha, angular_error, shift_error,
			helical_rise, helical_twist;

		bool deloc_supp, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak,
			do_fom_weighting, do_3d_rot, do_ewald;

		bool skip_gridding, debug, do_reconstruct_meas, is_positive, read_weights, div_avg;

		bool no_Wiener, writeWeights, new_Ewald_weight, Ewald_ellipsoid;

		float padding_factor, mask_diameter_ds, mask_diameter, mask_diameter_filt, flank_width;
		double padding_factor_2D;

		// I/O Parser
		IOParser parser;

		void usage()
		{
			parser.writeUsage(std::cerr);
		}

		void read(int argc, char **argv)
		{

			parser.setCommandLine(argc, argv);

			int general_section = parser.addSection("General options");
			fn_sel = parser.getOption("--i", "Input STAR file with the projection images and their orientations", "");
			fn_out = parser.getOption("--o", "Name for output reconstruction");
			fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
			padding_factor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));
			padding_factor_2D = textToDouble(parser.getOption("--pad2D", "Padding factor for 2D images", "1"));

			mask_diameter_filt = textToFloat(parser.getOption("--filter_diameter", "Diameter of filter-mask applied before division", "-1"));
			flank_width = textToFloat(parser.getOption("--filter_softness", "Width of filter-mask edge", "30"));
			nr_omp_threads = textToInteger(parser.getOption("--j", "Number of open-mp threads to use. Memory footprint is multiplied by this value.", "16"));

			int ctf_section = parser.addSection("CTF options");

			deloc_supp = parser.checkOption("--dm", "Apply delocalisation masking");
			mask_diameter_ds = textToDouble(parser.getOption("--mask_diameter_ds", "Diameter (in A) of mask for delocalisation suppression", "50"));
			intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");
			ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Images have been phase flipped");
			only_flip_phases = parser.checkOption("--only_flip_phases", "Do not correct CTF-amplitudes, only flip phases");

			read_weights = parser.checkOption("--read_weights", "Read freq. weight files");
			writeWeights = parser.checkOption("--write_weights", "Write the weights volume");
			do_ewald = parser.checkOption("--ewald", "Correct for Ewald-sphere curvature (developmental)");
			mask_diameter  = textToFloat(parser.getOption("--mask_diameter", "Diameter (in A) of mask for Ewald-sphere curvature correction", "-1."));
			width_mask_edge = textToInteger(parser.getOption("--width_mask_edge", "Width (in pixels) of the soft edge on the mask", "3"));
			is_positive = !parser.checkOption("--reverse_curvature", "Try curvature the other way around");
			newbox = textToInteger(parser.getOption("--newbox", "Box size of reconstruction after Ewald sphere correction", "-1"));
			nr_sectors = textToInteger(parser.getOption("--sectors", "Number of sectors for Ewald sphere correction", "2"));

			int helical_section = parser.addSection("Helical options");
			nr_helical_asu = textToInteger(parser.getOption("--nr_helical_asu", "Number of helical asymmetrical units", "1"));
			helical_rise = textToFloat(parser.getOption("--helical_rise", "Helical rise (in Angstroms)", "0."));
			helical_twist = textToFloat(parser.getOption("--helical_twist", "Helical twist (in degrees, + for right-handedness)", "0."));

			int expert_section = parser.addSection("Expert options");
			fn_sub = parser.getOption("--subtract","Subtract projections of this map from the images used for reconstruction", "");
			no_Wiener = parser.checkOption("--legacy", "Use gridding instead of Wiener filter");
			new_Ewald_weight = parser.checkOption("--new_Ewald_weight", "Use Ewald weight W that considers Cs as well");
			Ewald_ellipsoid = parser.checkOption("--Ewald_ellipsoid", "Allow Ewald sphere to become an ellipsoid under aniso. mag.");

			if (parser.checkOption("--NN", "Use nearest-neighbour instead of linear interpolation before gridding correction"))
			{
				interpolator = NEAREST_NEIGHBOUR;
			}
			else
			{
				interpolator = TRILINEAR;
			}

			blob_radius   = textToFloat(parser.getOption("--blob_r", "Radius of blob for gridding interpolation", "1.9"));
			blob_order    = textToInteger(parser.getOption("--blob_m", "Order of blob for gridding interpolation", "0"));
			blob_alpha    = textToFloat(parser.getOption("--blob_a", "Alpha-value of blob for gridding interpolation", "15"));
			grid_iters = textToInteger(parser.getOption("--iter", "Number of gridding-correction iterations", "10"));
			ref_dim = textToInteger(parser.getOption("--refdim", "Dimension of the reconstruction (2D or 3D)", "3"));
			angular_error = textToFloat(parser.getOption("--angular_error", "Apply random deviations with this standard deviation (in degrees) to each of the 3 Euler angles", "0."));
			shift_error = textToFloat(parser.getOption("--shift_error", "Apply random deviations with this standard deviation (in pixels) to each of the 2 translations", "0."));
			do_fom_weighting = parser.checkOption("--fom_weighting", "Weight particles according to their figure-of-merit (_rlnParticleFigureOfMerit)");
			do_3d_rot = parser.checkOption("--3d_rot", "Perform 3D rotations instead of backprojections from 2D images");
			skip_gridding = !parser.checkOption("--grid", "Perform gridding part of the reconstruction");
			div_avg = parser.checkOption("--div_avg", "Divide the per-voxel average by its weight prior to computing the preliminary FSC");

			debug = parser.checkOption("--debug", "Write out debugging data");

			// Hidden
			r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));

			// Check for errors in the command-line option
			if (parser.checkForErrors())
			{
				REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
			}
		}

		void applyCTFPandCTFQ(MultidimArray<Complex> &Fin, CTF &ctf, FourierTransformer &transformer,
		                      MultidimArray<Complex> &outP, MultidimArray<Complex> &outQ, double angpix)
		{
			//FourierTransformer transformer;
			outP.resize(Fin);
			outQ.resize(Fin);
			float angle_step = 180./nr_sectors;
			for (float angle = 0.; angle < 180.;  angle +=angle_step)
			{
				MultidimArray<Complex> CTFP(Fin), Fapp(Fin);
				MultidimArray<RFLOAT> Iapp(YSIZE(Fin), YSIZE(Fin));
				// Two passes: one for CTFP, one for CTFQ
				for (int ipass = 0; ipass < 2; ipass++)
				{
					bool is_my_positive = (ipass == 1) ? is_positive : !is_positive;

					// Get CTFP and multiply the Fapp with it
					ctf.getCTFPImage(CTFP, YSIZE(Fin), YSIZE(Fin), angpix, is_my_positive, angle);

					Fapp = Fin * CTFP; // element-wise complex multiplication!

					// inverse transform and mask out the particle....
					transformer.inverseFourierTransform(Fapp, Iapp);
					CenterFFT(Iapp, false);

					softMaskOutsideMap(Iapp, ROUND(mask_diameter/(angpix*2.)), (RFLOAT)width_mask_edge);

					// Re-box to a smaller size if necessary....
					if (newbox > 0 && newbox < YSIZE(Fin))
					{
						Iapp.setXmippOrigin();
						Iapp.window(FIRST_XMIPP_INDEX(newbox), FIRST_XMIPP_INDEX(newbox),
									LAST_XMIPP_INDEX(newbox),  LAST_XMIPP_INDEX(newbox));

					}

					// Back into Fourier-space
					CenterFFT(Iapp, true);
					transformer.FourierTransform(Iapp, Fapp, false); // false means: leave Fapp in the transformer

					// First time round: resize the output arrays
					if (ipass == 0 && fabs(angle) < XMIPP_EQUAL_ACCURACY)
					{
						outP.resize(Fapp);
						outQ.resize(Fapp);
					}

					// Now set back the right parts into outP (first pass) or outQ (second pass)
					float anglemin = angle + 90. - (0.5*angle_step);
					float anglemax = angle + 90. + (0.5*angle_step);

					// angles larger than 180
					bool is_reverse = false;
					if (anglemin >= 180.)
					{
						anglemin -= 180.;
						anglemax -= 180.;
						is_reverse = true;
					}
					MultidimArray<Complex> *myCTFPorQ, *myCTFPorQb;
					if (is_reverse)
					{
						myCTFPorQ  = (ipass == 0) ? &outQ : &outP;
						myCTFPorQb = (ipass == 0) ? &outP : &outQ;
					}
					else
					{
						myCTFPorQ  = (ipass == 0) ? &outP : &outQ;
						myCTFPorQb = (ipass == 0) ? &outQ : &outP;
					}

					// Deal with sectors with the Y-axis in the middle of the sector...
					bool do_wrap_max = false;
					if (anglemin < 180. && anglemax > 180.)
					{
						anglemax -= 180.;
						do_wrap_max = true;
					}

					// use radians instead of degrees
					anglemin = DEG2RAD(anglemin);
					anglemax = DEG2RAD(anglemax);
					FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(CTFP)
					{
						RFLOAT x = (RFLOAT)jp;
						RFLOAT y = (RFLOAT)ip;
						RFLOAT myangle = (x*x+y*y > 0) ? acos(y/sqrt(x*x+y*y)) : 0; // dot-product with Y-axis: (0,1)
						// Only take the relevant sector now...
						if (do_wrap_max)
						{
							if (myangle >= anglemin)
								DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
							else if (myangle < anglemax)
								DIRECT_A2D_ELEM(*myCTFPorQb, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
						}
						else
						{
							if (myangle >= anglemin && myangle < anglemax)
								DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
						}
					}
				}
			}
		}


		void reconstruct()
		{
			int data_dim = (do_3d_rot) ? 3 : 2;

			MultidimArray<RFLOAT> dummy;
			Image<RFLOAT> vol, sub;

			ObservationModel obsModel;
			MetaDataTable mdt0;

			ObservationModel::loadSafely(fn_sel, obsModel, mdt0);
			std::vector<double> angpix = obsModel.getPixelSizes();

			const int optGroupCount = obsModel.numberOfOpticsGroups();

			// Use pixel and box size of first opt. group for output;
			double angpixOut = angpix[0];
			int boxOut;

			// When doing Ewald-curvature correction: allow reconstructing smaller
			// box than the input images (which should have large boxes!!)
			if (do_ewald && newbox > 0)
			{
				boxOut = newbox;
			}
			else
			{
				boxOut = obsModel.getBoxSize(0);
			}

			std::vector<int> paddedSizes2D(optGroupCount);
			std::vector<int> origSizes2D(optGroupCount);

			for (int i = 0; i < optGroupCount; i++)
			{
				paddedSizes2D[i] = (int) (padding_factor_2D * obsModel.getBoxSize(i));
				origSizes2D[i] = (int) obsModel.getBoxSize(i);
			}

			// Get dimension of the images

			mdt0.getValue(EMDL_IMAGE_NAME, fn_img, 0);


			Projector subProjector(sub.data.xdim, interpolator, padding_factor, r_min_nn);

			r_max = -1;

			if (fn_sub != "")
			{
				sub.read(fn_sub);
				subProjector.computeFourierTransformMap(sub(), dummy, 2 * r_max);
			}

			std::vector<MetaDataTable> mdts = StackHelper::splitByStack(&mdt0);
			const long gc = mdts.size();

			std::vector<Image<RFLOAT>> prevRefs(2);
			std::vector<std::vector<BackProjector>> backprojectors(2);

			for (int j = 0; j < 2; j++)
			{
				backprojectors[j] = std::vector<BackProjector>(nr_omp_threads);

				for (int i = 0; i < nr_omp_threads; i++)
				{
					backprojectors[j][i] = BackProjector(
						boxOut, ref_dim, fn_sym, interpolator,
						padding_factor, r_min_nn, blob_order,
						blob_radius, blob_alpha, data_dim, skip_gridding);
				}
			}

			std::cout << "Back-projecting all images ..." << std::endl;

			time_config();
			init_progress_bar(gc/nr_omp_threads);


			#pragma omp parallel num_threads(nr_omp_threads)
			{
				int threadnum = omp_get_thread_num();

				backprojectors[0][threadnum].initZeros(2 * r_max);
				backprojectors[1][threadnum].initZeros(2 * r_max);

				RFLOAT rot, tilt, psi, fom, r_ewald_sphere;
				Matrix2D<RFLOAT> A3D;
				MultidimArray<RFLOAT> Fctf;
				Matrix1D<RFLOAT> trans(2);
				FourierTransformer transformer;

				#pragma omp for
				for (int g = 0; g < gc; g++)
				{
					std::vector<Image<RFLOAT> > obsR;

					try
					{
						obsR = StackHelper::loadStack(&mdts[g]);
					}
					catch (RelionError XE)
					{
						std::cerr << "warning: unable to load micrograph #" << (g+1) << "\n";
						continue;
					}

					const long pc = obsR.size();

					for (int p = 0; p < pc; p++)
					{
						int randSubset;
						mdts[g].getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
						randSubset = randSubset - 1;

						// Rotations
						if (ref_dim == 2)
						{
							rot = tilt = 0.;
						}
						else
						{
							mdts[g].getValue(EMDL_ORIENT_ROT, rot, p);
							mdts[g].getValue(EMDL_ORIENT_TILT, tilt, p);
						}

						psi = 0.;
						mdts[g].getValue(EMDL_ORIENT_PSI, psi, p);

						if (angular_error > 0.)
						{
							rot += rnd_gaus(0., angular_error);
							tilt += rnd_gaus(0., angular_error);
							psi += rnd_gaus(0., angular_error);
							//std::cout << rnd_gaus(0., angular_error) << std::endl;
						}

						Euler_angles2matrix(rot, tilt, psi, A3D);

						int opticsGroup = obsModel.getOpticsGroup(mdts[g], p);

						// If we are considering Ewald sphere curvature, the mag. matrix
						// has to be provided to the backprojector explicitly
						// (to avoid creating an Ewald ellipsoid)
						if (!do_ewald || Ewald_ellipsoid)
						{
							A3D = obsModel.applyAnisoMag(A3D, opticsGroup);
						}

						A3D = obsModel.applyScaleDifference(A3D, opticsGroup, boxOut, angpixOut);
						A3D /= padding_factor_2D;

						// Translations (either through phase-shifts or in real space
						trans.initZeros();
						mdts[g].getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, XX(trans), p);
						mdts[g].getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, YY(trans), p);

						XX(trans) /= angpix[opticsGroup];
						YY(trans) /= angpix[opticsGroup];

						if (shift_error > 0.)
						{
							XX(trans) += rnd_gaus(0., shift_error);
							YY(trans) += rnd_gaus(0., shift_error);
						}

						if (do_3d_rot)
						{
							trans.resize(3);
							mdts[g].getValue( EMDL_ORIENT_ORIGIN_Z, ZZ(trans), p);

							if (shift_error > 0.)
							{
								ZZ(trans) += rnd_gaus(0., shift_error);
							}
						}

						if (do_fom_weighting)
						{
							mdts[g].getValue( EMDL_PARTICLE_FOM, fom, p);
						}

						MultidimArray<Complex> Fsub, F2D, F2DP, F2DQ;

						CenterFFT(obsR[p](), true);

						const int sPad2D = paddedSizes2D[opticsGroup];

						if (padding_factor_2D > 1.0)
						{
							obsR[p] = FilterHelper::padCorner2D(obsR[p], sPad2D, sPad2D);
						}

						transformer.FourierTransform(obsR[p](), F2D);

						if (ABS(XX(trans)) > 0. || ABS(YY(trans)) > 0.)
						{
							if (do_3d_rot)
							{
								shiftImageInFourierTransform(
									F2D, F2D, sPad2D, XX(trans), YY(trans), ZZ(trans));
							}
							else
							{
								shiftImageInFourierTransform(
									F2D, F2D, sPad2D, XX(trans), YY(trans));
							}
						}

						Fctf.resize(F2D);
						Fctf.initConstant(1.);

						CTF ctf;
						ctf.readByGroup(mdts[g], &obsModel, p);

						ctf.getFftwImage(Fctf, sPad2D, sPad2D, angpix[opticsGroup],
										 ctf_phase_flipped, only_flip_phases,
										 intact_ctf_first_peak, true);

						if (deloc_supp)
						{
							DelocalisationHelper::maskOutsideBox(
								ctf, mask_diameter_ds/(2.0 * angpix[opticsGroup]),
								angpix[opticsGroup], origSizes2D[opticsGroup],
								Fctf, XX(trans), YY(trans));
						}

						obsModel.demodulatePhase(mdts[g], p, F2D);
						obsModel.divideByMtf(mdts[g], p, F2D);

						if (do_ewald)
						{
							// Ewald-sphere curvature correction
							applyCTFPandCTFQ(F2D, ctf, transformer, F2DP, F2DQ, angpix[opticsGroup]);

							// Also calculate W, store again in Fctf

							if (new_Ewald_weight)
							{
								ctf.applyWeightEwaldSphereCurvature_new(
									Fctf, sPad2D, sPad2D, angpix[opticsGroup], mask_diameter);
							}
							else
							{
								ctf.applyWeightEwaldSphereCurvature(
									Fctf, sPad2D, sPad2D, angpix[opticsGroup], mask_diameter);
							}

							// Also calculate the radius of the Ewald sphere (in pixels)
							r_ewald_sphere = boxOut * angpix[opticsGroup] / ctf.lambda;
						}

						// Subtract reference projection
						if (fn_sub != "")
						{
							obsModel.predictObservation(
								subProjector, mdts[g], p, Fsub, true, true, true);

							FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fsub)
							{
								DIRECT_MULTIDIM_ELEM(F2D, n) -= DIRECT_MULTIDIM_ELEM(Fsub, n);
							}

							// Back-project difference image
							backprojectors[randSubset][threadnum].set2DFourierTransform(
										F2D, A3D);
						}
						else
						{
							if (do_ewald)
							{
								FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
								{
									DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
								}
							}
							// "Normal" reconstruction, multiply X by CTF, and W by CTF^2
							else
							{
								FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
								{
									DIRECT_MULTIDIM_ELEM(F2D, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
									DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
								}
							}

							// Do the following after squaring the CTFs!
							if (do_fom_weighting)
							{
								FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
								{
									DIRECT_MULTIDIM_ELEM(F2D, n)  *= fom;
									DIRECT_MULTIDIM_ELEM(Fctf, n) *= fom;
								}
							}

							DIRECT_A2D_ELEM(F2D, 0, 0) = 0.0;

							if (do_ewald)
							{
								Matrix2D<RFLOAT> magMat;

								if (obsModel.hasMagMatrices && !Ewald_ellipsoid)
								{
									magMat = obsModel.getMagMatrix(opticsGroup);
								}
								else
								{
									magMat = Matrix2D<RFLOAT>(2,2);
									magMat.initIdentity();
								}

								backprojectors[randSubset][threadnum].set2DFourierTransform(
									F2DP, A3D, &Fctf, r_ewald_sphere, true, &magMat);

								backprojectors[randSubset][threadnum].set2DFourierTransform(
									F2DQ, A3D, &Fctf, r_ewald_sphere, false, &magMat);
							}
							else
							{
								backprojectors[randSubset][threadnum].set2DFourierTransform(
											F2D, A3D, &Fctf);
							}
						}

						if (threadnum == 0)
						{
							progress_bar(g);
						}
					}
				}
			}

			progress_bar(gc/nr_omp_threads);


			std::vector<BackProjector*> backprojector(2);

			for (int j = 0; j < 2; j++)
			{
				std::cerr << " + Merging volumes for half-set " << (j+1) << "...\n";

				backprojector[j] = &backprojectors[j][0];

				for (int bpi = 1; bpi < nr_omp_threads; bpi++)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(backprojector[j]->data)
					{
						DIRECT_MULTIDIM_ELEM(backprojector[j]->data, n)
								+= DIRECT_MULTIDIM_ELEM(backprojectors[j][bpi].data, n);
					}

					backprojectors[j][bpi].data.clear();

					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(backprojector[j]->weight)
					{
						DIRECT_MULTIDIM_ELEM(backprojector[j]->weight, n)
								+= DIRECT_MULTIDIM_ELEM(backprojectors[j][bpi].weight, n);
					}

					backprojectors[j][bpi].weight.clear();
				}

				std::cerr << " + Symmetrising half-set " << (j+1) << "...\n";

				backprojector[j]->symmetrise(
					nr_helical_asu, helical_twist, helical_rise/angpixOut, nr_omp_threads);
			}

			bool do_map = !no_Wiener;
			bool do_use_fsc = !no_Wiener;

			MultidimArray<RFLOAT> fsc(boxOut/2 + 1);

			if (!no_Wiener)
			{
				MultidimArray<Complex> avg0, avg1;

				backprojector[0]->getDownsampledAverage(avg0, div_avg);
				backprojector[1]->getDownsampledAverage(avg1, div_avg);
				backprojector[0]->calculateDownSampledFourierShellCorrelation(avg0, avg1, fsc);
			}

			if (debug)
			{
				std::ofstream fscNew(fn_out+"_prelim_FSC.dat");

				for (int i = 0; i < fsc.xdim; i++)
				{
					fscNew << i << " " << fsc(i) << "\n";
				}
			}

			for (int j = 0; j < 2; j++)
			{
				if (mask_diameter_filt > 0.0)
				{
					std::cout << " + Applying spherical mask of diameter " <<
							  mask_diameter_filt << " ..." << std::endl;

					const double r0 = mask_diameter_filt/2.0;
					const double r1 = r0 + flank_width;

					Image<Complex> tempC;
					Image<RFLOAT> tempR;

					BackProjector::decenterWhole(backprojector[j]->data, tempC());
					NewFFT::inverseFourierTransform(tempC(), tempR(), NewFFT::FwdOnly, false);
					tempR = FilterHelper::raisedCosEnvCorner3D(tempR, r0, r1);
					NewFFT::FourierTransform(tempR(), tempC(), NewFFT::FwdOnly);
					BackProjector::recenterWhole(tempC(), backprojector[j]->data);

					BackProjector::decenterWhole(backprojector[j]->weight, tempC());
					NewFFT::inverseFourierTransform(tempC(), tempR(), NewFFT::FwdOnly, false);
					tempR = FilterHelper::raisedCosEnvCorner3D(tempR, r0, r1);
					NewFFT::FourierTransform(tempR(), tempC(), NewFFT::FwdOnly);
					BackProjector::recenterWhole(tempC(), backprojector[j]->weight);
				}

				Image<RFLOAT> weightOut;

				std::cout << " + Starting the reconstruction ..." << std::endl;

				MultidimArray<RFLOAT> tau2;
				if (do_use_fsc) backprojector[j]->updateSSNRarrays(1., tau2, dummy, dummy, dummy, fsc, do_use_fsc, true);
				backprojector[j]->reconstruct(vol(), grid_iters, do_map, tau2,
						1., 1., -1, false, writeWeights? &weightOut : 0);

				if (writeWeights)
				{
					std::stringstream sts;
					sts << (j+1);
					std::string fnWgh = fn_out + "_half" + sts.str() + "_class001_unfil_weight.mrc";
					weightOut.write(fnWgh);
				}

				prevRefs[j] = vol;

			} // halves

			for (int j = 0; j < 2; j++)
			{
				std::stringstream sts;
				sts << (j+1);

				std::string fnFull = fn_out + "_half" + sts.str() + "_class001_unfil.mrc";

				prevRefs[j].write(fnFull);
				std::cout << " Done writing map in " << fnFull << "\n";
			}
		}
};


int main(int argc, char *argv[])
{
	reconstruct_parameters prm;

	try
	{
		prm.read(argc, argv);

		prm.reconstruct();
	}
	catch (RelionError XE)
	{
		//prm.usage();
		std::cerr << XE;
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}


