#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include "src/gpu_utils/cuda_mem_utils.h"
#include "src/complex.h"
#include <fstream>
#include <cuda_runtime.h>
#include <signal.h>

#include "src/image.h"
#include "src/autopicker.h"
#include "src/gpu_utils/cuda_autopicker.h"

#ifdef CUDA_FORCESTL
#include "src/gpu_utils/cuda_utils_stl.cuh"
#else
#include "src/gpu_utils/cuda_utils_cub.cuh"
#endif

static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

AutoPickerCuda::AutoPickerCuda(AutoPicker *basePicker, int dev_id) :
	basePckr(basePicker)
{

	/*======================================================
	                    DEVICE SETTINGS
	======================================================*/
	device_id = dev_id;
	int devCount;
	HANDLE_ERROR(cudaGetDeviceCount(&devCount));

	if(dev_id >= devCount)
	{
		std::cerr << " using device_id=" << dev_id << " (device no. " << dev_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		REPORT_ERROR("ERROR: Assigning a thread to a non-existent device (index likely too high)");
	}
	else
		HANDLE_ERROR(cudaSetDevice(dev_id));


	/*======================================================
	                    CUSTOM ALLOCATOR
	======================================================*/

#ifdef CUDA_NO_CUSTOM_ALLOCATION
	printf(" DEBUG: Custom allocator is disabled.\n");
	allocator = new CudaCustomAllocator(0, 1);
#else
	size_t allocationSize(0);

	size_t free, total;
	HANDLE_ERROR(cudaMemGetInfo( &free, &total ));

	if (basePckr->available_gpu_memory > 0)
		allocationSize = basePckr->available_gpu_memory * (1000*1000*1000);
	else
		allocationSize = (float)free * .7;

	if (allocationSize > free)
	{
		printf(" WARNING: Required memory per thread, via \"--gpu_memory_per_thread\", not available on device. (Defaulting to less)\n");
		allocationSize = (float)free * .7; //Lets leave some for other processes for now
	}

	int memAlignmentSize;
	cudaDeviceGetAttribute ( &memAlignmentSize, cudaDevAttrTextureAlignment, dev_id );

	allocator = new CudaCustomAllocator(allocationSize, memAlignmentSize);
#endif
};


void AutoPickerCuda::run()
{

	int barstep;
	if (basePckr->verb > 0)
	{
		std::cout << " Autopicking ..." << std::endl;
		init_progress_bar(basePckr->fn_micrographs.size());
		barstep = XMIPP_MAX(1,basePckr->fn_micrographs.size() / 60);
	}


	for (long int imic = 0; imic < basePckr->fn_micrographs.size(); imic++)
	{
		if (basePckr->verb > 0 && imic % barstep == 0)
			progress_bar(imic);

		autoPickOneMicrograph(basePckr->fn_micrographs[imic]);
	}

	if (basePckr->verb > 0)
		progress_bar(basePckr->fn_micrographs.size());

}


void AutoPickerCuda::autoPickOneMicrograph(FileName &fn_mic)
{
	Image<RFLOAT> Imic;
	MultidimArray<Complex > Faux, Faux2, Fmic;
	MultidimArray<RFLOAT> Maux, Mstddev, Mmean, Mdiff2, MsumX2, Mccf_best, Mpsi_best, Fctf;
	FourierTransformer transformer;
	RFLOAT sum_ref_under_circ_mask, sum_ref2_under_circ_mask;
	int my_skip_side = basePckr->autopick_skip_side + basePckr->particle_size/2;
	CTF ctf;

	int min_distance_pix = ROUND(basePckr->min_particle_distance / basePckr->angpix);

#ifdef DEBUG
	Image<RFLOAT> tt;
	tt().resize(micrograph_size, micrograph_size);
	std::cerr << " fn_mic= " << fn_mic << std::endl;
#endif
	// Read in the micrograph
	Imic.read(fn_mic);
	Imic().setXmippOrigin();

	// Let's just check the square size again....
	RFLOAT my_size, my_xsize, my_ysize;
	my_xsize = XSIZE(Imic());
	my_ysize = YSIZE(Imic());
	my_size = (my_xsize != my_ysize) ? XMIPP_MAX(my_xsize, my_ysize) : my_xsize;

	if (my_size != basePckr->micrograph_size || my_xsize != basePckr->micrograph_xsize || my_ysize != basePckr->micrograph_ysize)
	{
		Imic().printShape();
		std::cerr << " micrograph_size= " << basePckr->micrograph_size << " micrograph_xsize= " << basePckr->micrograph_xsize << " micrograph_ysize= " << basePckr->micrograph_ysize << std::endl;
		REPORT_ERROR("AutoPicker::autoPickOneMicrograph ERROR: No differently sized micrographs are allowed in one run, sorry you will have to run separately for each size...");
	}

	// Set mean to zero and stddev to 1 to prevent numerical problems with one-sweep stddev calculations....
    RFLOAT avg0, stddev0, minval0, maxval0;
	Imic().computeStats(avg0, stddev0, minval0, maxval0);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imic())
	{
		// Remove pixel values that are too far away from the mean
		if ( ABS(DIRECT_MULTIDIM_ELEM(Imic(), n) - avg0) / stddev0 > basePckr->outlier_removal_zscore)
			DIRECT_MULTIDIM_ELEM(Imic(), n) = avg0;

		DIRECT_MULTIDIM_ELEM(Imic(), n) = (DIRECT_MULTIDIM_ELEM(Imic(), n) - avg0) / stddev0;
	}

	if (basePckr->micrograph_xsize !=basePckr->micrograph_ysize)
	{
		// Window non-square micrographs to be a square with the largest side
		rewindow(Imic, basePckr->micrograph_size);

		// Fill region outside the original window with white Gaussian noise to prevent all-zeros in Mstddev
		FOR_ALL_ELEMENTS_IN_ARRAY2D(Imic())
		{
			if (i < FIRST_XMIPP_INDEX(basePckr->micrograph_ysize)
					|| i > LAST_XMIPP_INDEX(basePckr->micrograph_ysize)
					|| j < FIRST_XMIPP_INDEX(basePckr->micrograph_xsize)
					|| j > LAST_XMIPP_INDEX(basePckr->micrograph_xsize) )
				A2D_ELEM(Imic(), i, j) = rnd_gaus(0.,1.);
		}
	}

	// Read in the CTF information if needed
	if (basePckr->do_ctf)
	{
		// Search for this micrograph in the metadata table
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(basePckr->MDmic)
		{
			FileName fn_tmp;
			basePckr->MDmic.getValue(EMDL_MICROGRAPH_NAME, fn_tmp);
			if (fn_tmp==fn_mic)
			{
				ctf.read(basePckr->MDmic, basePckr->MDmic);
				Fctf.resize(basePckr->downsize_mic, basePckr->downsize_mic/2 + 1);
				ctf.getFftwImage(Fctf, basePckr->micrograph_size, basePckr->micrograph_size, basePckr->angpix, false, false, basePckr->intact_ctf_first_peak, true);
				break;
			}
		}
#ifdef DEBUG
		std::cerr << " Read CTF info from" << fn_mic.withoutExtension()<<"_ctf.star" << std::endl;
		Image<RFLOAT> Ictf;
		Ictf()=Fctf;
		Ictf.write("Mmic_ctf.spi");
#endif
	}

	Mccf_best.resize(basePckr->micrograph_size, basePckr->micrograph_size);
	Mpsi_best.resize(basePckr->micrograph_size, basePckr->micrograph_size);

	RFLOAT normfft = (RFLOAT)(basePckr->micrograph_size * basePckr->micrograph_size) / (RFLOAT)basePckr->nr_pixels_circular_mask;;
	if (basePckr->do_read_fom_maps)
	{
		FileName fn_tmp=fn_mic.withoutExtension()+"_"+basePckr->fn_out+"_stddevNoise.spi";
		Image<RFLOAT> It;
		It.read(fn_tmp);
		Mstddev = It();
	}
	else
	{
		/*
		 * Squared difference FOM:
		 * Sum ( (X-mu)/sig  - A )^2 =
		 *  = Sum((X-mu)/sig)^2 - 2 Sum (A*(X-mu)/sig) + Sum(A)^2
		 *  = (1/sig^2)*Sum(X^2) - (2*mu/sig^2)*Sum(X) + (mu^2/sig^2)*Sum(1) - (2/sig)*Sum(AX) + (2*mu/sig)*Sum(A) + Sum(A^2)
		 *
		 * However, the squared difference with an "empty" ie all-zero reference is:
		 * Sum ( (X-mu)/sig)^2
		 *
		 * The ratio of the probabilities thereby becomes:
		 * P(ref) = 1/sqrt(2pi) * exp (( (X-mu)/sig  - A )^2 / -2 )   // assuming sigma = 1!
		 * P(zero) = 1/sqrt(2pi) * exp (( (X-mu)/sig )^2 / -2 )
		 *
		 * P(ref)/P(zero) = exp(( (X-mu)/sig  - A )^2 / -2) / exp ( ( (X-mu)/sig )^2 / -2)
		 *                = exp( (- (2/sig)*Sum(AX) + (2*mu/sig)*Sum(A) + Sum(A^2)) / - 2 )
		 *
		 *                Therefore, I do not need to calculate (X-mu)/sig beforehand!!!
		 *
		 */

		// Fourier Transform (and downscale) Imic()
		CenterFFT(Imic(), true);
		transformer.FourierTransform(Imic(), Fmic);

		// Also calculate the FFT of the squared micrograph
		Maux.resize(Imic());
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Maux)
		{
			DIRECT_MULTIDIM_ELEM(Maux, n) = DIRECT_MULTIDIM_ELEM(Imic(), n) * DIRECT_MULTIDIM_ELEM(Imic(), n);
		}
		MultidimArray<Complex > Fmic2;
		transformer.FourierTransform(Maux, Fmic2);

#ifdef DEBUG
		std::cerr << " nr_pixels_circular_invmask= " << nr_pixels_circular_invmask << std::endl;
		std::cerr << " nr_pixels_circular_mask= " << nr_pixels_circular_mask << std::endl;
		windowFourierTransform(Finvmsk, Faux2, micrograph_size);
		transformer.inverseFourierTransform(Faux2, tt());
		CenterFFT(tt(), false);
		tt.write("Minvmask.spi");
		windowFourierTransform(Fmsk, Faux2, micrograph_size);
		transformer.inverseFourierTransform(Faux2, tt());
		CenterFFT(tt(), false);
		tt.write("Mmask.spi");
#endif

		// The following calculate mu and sig under the solvent area at every position in the micrograph
		basePckr->calculateStddevAndMeanUnderMask(Fmic, Fmic2, basePckr->Finvmsk,basePckr->nr_pixels_circular_invmask, Mstddev, Mmean);

		if (basePckr->do_write_fom_maps)
		{
			// TMP output
			FileName fn_tmp=fn_mic.withoutExtension()+"_"+basePckr->fn_out+"_stddevNoise.spi";
			Image<RFLOAT> It;
			It() = Mstddev;
			It.write(fn_tmp);
		}

		// From now on use downsized Fmic, as the cross-correlation with the references can be done at lower resolution
		windowFourierTransform(Fmic, Faux, basePckr->downsize_mic);
		Fmic = Faux;

	}// end if do_read_fom_maps

	// Now start looking for the peaks of all references
	// Clear the output vector with all peaks
	std::vector<Peak> peaks;
	peaks.clear();
	for (int iref = 0; iref < basePckr->Mrefs.size(); iref++)
	{
		RFLOAT expected_Pratio; // the expectedFOM for this (ctf-corrected) reference
		if (basePckr->do_read_fom_maps)
		{
			FileName fn_tmp;
			Image<RFLOAT> It;
			fn_tmp.compose(fn_mic.withoutExtension()+"_"+basePckr->fn_out+"_ref", iref,"_bestCCF.spi");
			It.read(fn_tmp);
			Mccf_best = It();
			// Retrieve expected_Pratio from the header of the image..
			It.MDMainHeader.getValue(EMDL_IMAGE_STATS_MAX, expected_Pratio);
			fn_tmp.compose(fn_mic.withoutExtension()+"_"+basePckr->fn_out+"_ref", iref,"_bestPSI.spi");
			It.read(fn_tmp);
			Mpsi_best = It();

		} //end else if do_read_fom_maps
		else
		{
			Mccf_best.initConstant(-LARGE_NUMBER);
			bool is_first_psi = true;
			for (RFLOAT psi = 0. ; psi < 360.; psi+=basePckr->psi_sampling)
			{

				// Get the Euler matrix
				Matrix2D<RFLOAT> A(3,3);
				Euler_angles2matrix(0., 0., psi, A);

				// Now get the FT of the rotated (non-ctf-corrected) template
				Faux.initZeros(basePckr->downsize_mic, basePckr->downsize_mic/2 + 1);
				basePckr->PPref[iref].get2DFourierTransform(Faux, A, IS_NOT_INV);

#ifdef DEBUG
				std::cerr << " psi= " << psi << std::endl;
				windowFourierTransform(Faux, Faux2, micrograph_size);
				transformer.inverseFourierTransform(Faux2, tt());
				CenterFFT(tt(), false);
				tt.write("Mref_rot.spi");

				windowFourierTransform(Fmic, Faux2, micrograph_size);
				transformer.inverseFourierTransform(Faux2, tt());
				CenterFFT(tt(), false);
				tt.write("Mmic.spi");

#endif

				// Apply the CTF on-the-fly (so same PPref can be used for many different micrographs)
				if (basePckr->do_ctf)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
					{
						DIRECT_MULTIDIM_ELEM(Faux, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
#ifdef DEBUG
				windowFourierTransform(Faux, Faux2, micrograph_size);
				transformer.inverseFourierTransform(Faux2, Maux);
				CenterFFT(Maux, false);
				Maux.setXmippOrigin();
				tt().resize(particle_size, particle_size);
				tt().setXmippOrigin();
				FOR_ALL_ELEMENTS_IN_ARRAY2D(tt())
				{
					A2D_ELEM(tt(), i, j) = A2D_ELEM(Maux, i, j);
				}
				tt.write("Mref_rot_ctf.spi");
#endif
				}

				if (is_first_psi)
				{
					// Calculate the expected ratio of probabilities for this CTF-corrected reference
					// and the sum_ref_under_circ_mask and sum_ref_under_circ_mask2
					// Do this also if we're not recalculating the fom maps...

					windowFourierTransform(Faux, Faux2, basePckr->micrograph_size);
					transformer.inverseFourierTransform(Faux2, Maux);
					CenterFFT(Maux, false);
					Maux.setXmippOrigin();
					// TODO: check whether I need CenterFFT(Maux, false)

					sum_ref_under_circ_mask = 0.;
					sum_ref2_under_circ_mask = 0.;
					RFLOAT suma2 = 0.;
					RFLOAT sumn = 1.;
					MultidimArray<RFLOAT> Mctfref(basePckr->particle_size, basePckr->particle_size);
					Mctfref.setXmippOrigin();
					FOR_ALL_ELEMENTS_IN_ARRAY2D(Mctfref) // only loop over smaller Mctfref, but take values from large Maux!
					{
						if (i*i + j*j < basePckr->particle_radius2)
						{
							suma2 += A2D_ELEM(Maux, i, j) * A2D_ELEM(Maux, i, j);
							suma2 += 2. * A2D_ELEM(Maux, i, j) * rnd_gaus(0., 1.);
							sum_ref_under_circ_mask += A2D_ELEM(Maux, i, j);
							sum_ref2_under_circ_mask += A2D_ELEM(Maux, i, j) * A2D_ELEM(Maux, i, j);
							sumn += 1.;
						}
#ifdef DEBUG
						A2D_ELEM(Mctfref, i, j) = A2D_ELEM(Maux, i, j);
#endif
					}
					sum_ref_under_circ_mask /= sumn;
					sum_ref2_under_circ_mask /= sumn;
					expected_Pratio = exp(suma2 / (2. * sumn));
#ifdef DEBUG
					std::cerr << " expected_Pratio["<<iref<<"]= " << expected_Pratio << std::endl;
					tt()=Mctfref;
					tt.write("Mctfref.spi");
#endif
				}

				// Now multiply template and micrograph to calculate the cross-correlation
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
				{
					DIRECT_MULTIDIM_ELEM(Faux, n) = conj(DIRECT_MULTIDIM_ELEM(Faux, n)) * DIRECT_MULTIDIM_ELEM(Fmic, n);
				}
				windowFourierTransform(Faux, Faux2, basePckr->micrograph_size);
				transformer.inverseFourierTransform(Faux2, Maux);
				CenterFFT(Maux, false);
#ifdef DEBUG
				tt()=Maux*normfft;
				tt.write("Mcc.spi");
#endif

				// Calculate ratio of prabilities P(ref)/P(zero)
				// Keep track of the best values and their corresponding iref and psi

				// So now we already had precalculated: Mdiff2 = 1/sig*Sum(X^2) - 2/sig*Sum(X) + mu^2/sig*Sum(1)
				// Still to do (per reference): - 2/sig*Sum(AX) + 2*mu/sig*Sum(A) + Sum(A^2)
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Maux)
				{
					RFLOAT diff2 = - 2. * normfft * DIRECT_MULTIDIM_ELEM(Maux, n);
					diff2 += 2. * DIRECT_MULTIDIM_ELEM(Mmean, n) * sum_ref_under_circ_mask;
					if (DIRECT_MULTIDIM_ELEM(Mstddev, n) > 1E-10)
						diff2 /= DIRECT_MULTIDIM_ELEM(Mstddev, n);
					diff2 += sum_ref2_under_circ_mask;
#ifdef DEBUG
						/*
						if (diff2 < 0. || n==28800 || n==0)
						{
							std::cerr << " n= "<<n<< "diff2= " << diff2 << " old Mdiff2=" <<DIRECT_MULTIDIM_ELEM(Mdiff2, n)
									<< " -2AX/sig " << - 2. * normfft * DIRECT_MULTIDIM_ELEM(Maux, n) / DIRECT_MULTIDIM_ELEM(Mstddev, n)
									<< " 2Amu/sig= " << 2. * DIRECT_MULTIDIM_ELEM(Mmean, n) * sum_ref_under_circ_mask[iref] / DIRECT_MULTIDIM_ELEM(Mstddev, n)
									<< " A2=" <<  sum_ref2_under_circ_mask[iref]
									<< " stddev= " <<  DIRECT_MULTIDIM_ELEM(Mstddev, n) << " avg= "<< DIRECT_MULTIDIM_ELEM(Mmean, n)
									<< std::endl;
						}
						*/
#endif
					diff2 = exp(- diff2 / 2.); // exponentiate to reflect the Gaussian error model. sigma=1 after normalization, 0.4=1/sqrt(2pi)

					// Store fraction of (1 - probability-ratio) wrt  (1 - expected Pratio)
					diff2 = (diff2 - 1.) / (expected_Pratio - 1.);
#ifdef DEBUG
					DIRECT_MULTIDIM_ELEM(Maux, n) = diff2;
#endif
					if (diff2 > DIRECT_MULTIDIM_ELEM(Mccf_best, n))
					{
						DIRECT_MULTIDIM_ELEM(Mccf_best, n) = diff2;
						DIRECT_MULTIDIM_ELEM(Mpsi_best, n) = psi;
					}
				}
#ifdef DEBUG
				std::cerr << " Maux.computeMax()= " << Maux.computeMax() << std::endl;
				tt()=Maux;
				tt.write("Mccf.spi");
			    std::cerr << " Press any key to continue... "  << std::endl;
			    char c;
			    std::cin >> c;

#endif
			    is_first_psi = false;
			} // end for psi


			if (basePckr->do_write_fom_maps)
			{
				// TMP output
				FileName fn_tmp;
				Image<RFLOAT> It;
				It() = Mccf_best;
				// Store expected_Pratio in the header of the image..
				It.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, expected_Pratio);;
				fn_tmp.compose(fn_mic.withoutExtension()+"_"+basePckr->fn_out+"_ref", iref,"_bestCCF.spi");
				It.write(fn_tmp);

				It() = Mpsi_best;
				fn_tmp.compose(fn_mic.withoutExtension()+"_"+basePckr->fn_out+"_ref", iref,"_bestPSI.spi");
				It.write(fn_tmp);
			} // end if do_write_fom_maps

		} // end if do_read_fom_maps

		// Now that we have Mccf_best and Mpsi_best, get the peaks
		std::vector<Peak> my_ref_peaks;
		Mstddev.setXmippOrigin();
		Mccf_best.setXmippOrigin();
		Mpsi_best.setXmippOrigin();
		basePckr->peakSearch(Mccf_best, Mpsi_best, Mstddev, iref, my_skip_side, my_ref_peaks);

		basePckr->prunePeakClusters(my_ref_peaks, min_distance_pix);

		// append the peaks of this reference to all the other peaks
		peaks.insert(peaks.end(), my_ref_peaks.begin(), my_ref_peaks.end());

	} // end for iref


	//Now that we have done all references, prune the list again...
	basePckr->prunePeakClusters(peaks, min_distance_pix);

	// And remove all too close neighbours
	basePckr->removeTooCloselyNeighbouringPeaks(peaks, min_distance_pix);

	// Write out a STAR file with the coordinates
	MetaDataTable MDout;
	for (int ipeak =0; ipeak < peaks.size(); ipeak++)
	{
		MDout.addObject();
		MDout.setValue(EMDL_IMAGE_COORD_X, (RFLOAT)(peaks[ipeak].x));
		MDout.setValue(EMDL_IMAGE_COORD_Y, (RFLOAT)(peaks[ipeak].y));
		MDout.setValue(EMDL_ORIENT_PSI, peaks[ipeak].psi);
		MDout.setValue(EMDL_PARTICLE_CLASS, peaks[ipeak].ref + 1); // start counting at 1
		MDout.setValue(EMDL_PARTICLE_AUTOPICK_FOM, peaks[ipeak].fom);
	}
	FileName fn_tmp = fn_mic.withoutExtension() + "_" + basePckr->fn_out + ".star";
	MDout.write(fn_tmp);

}
