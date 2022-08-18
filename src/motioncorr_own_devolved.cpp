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
#include "src/motioncorr_own_devolved.h"

void MotioncorrOwnDevolved::read(int argc, char **argv)
{
	// Define a new MpiNode
	node = new MpiNode(argc, argv);

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	movie_path = parser.getOption("--in_movie", "Path to input movie");
	micrograph_path = parser.getOption("--out_mic", "Output micrograph path");
	n_threads = textToInteger(parser.getOption("--j", "Number of threads per movie (= process)", "1"));
	max_io_threads = textToInteger(parser.getOption("--max_io_threads", "Limit the number of IO threads.", "-1"));
	grouping_for_ps = textToInteger(parser.getOption("--grouping_for_ps", "Group this number of frames and write summed power spectrum. -1 == do not write", "-1"));
	ps_size = textToInteger(parser.getOption("--ps_size", "Output size of power spectrum", "512"));
	if (ps_size % 2 != 0) REPORT_ERROR("--ps_size must be an even number.");
	first_frame_sum =  textToInteger(parser.getOption("--first_frame_sum", "First movie frame used in output sum (start at 1)", "1"));
	if (first_frame_sum < 1) first_frame_sum = 1;
	last_frame_sum =  textToInteger(parser.getOption("--last_frame_sum", "Last movie frame used in output sum (0 or negative: use all)", "-1"));
	eer_grouping = textToInteger(parser.getOption("--eer_grouping", "EER grouping", "40"));
	eer_upsampling = textToInteger(parser.getOption("--eer_upsampling", "EER upsampling (1 = 4K or 2 = 8K)", "1"));

	int doseweight_section = parser.addSection("Dose-weighting options");
	do_dose_weighting = parser.checkOption("--dose_weighting", "Use dose-weighting scheme");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms", "-1"));
	voltage = textToFloat(parser.getOption("--voltage","Voltage (in kV) for dose-weighting", "-1"));
	dose_per_frame = textToFloat(parser.getOption("--dose_per_frame", "Electron dose (in electrons/A2/frame) for dose-weighting", "1"));
	pre_exposure = textToFloat(parser.getOption("--preexposure", "Pre-exposure (in electrons/A2) for dose-weighting", "0"));

	parser.addSection("Own motion correction options");
	bin_factor =  textToFloat(parser.getOption("--bin_factor", "Binning factor (can be non-integer)", "1"));
	bfactor =  textToFloat(parser.getOption("--bfactor", "B-factor (in pix^2) that will be used inside MOTIONCOR2", "150"));
	fn_gain_reference = parser.getOption("--gainref","Location of MRC file with the gain reference to be applied","");
	gain_rotation = textToInteger(parser.getOption("--gain_rot", "Rotate the gain reference this number times 90 degrees clock-wise (in relion_display). This is same as MotionCor2's RotGain. 0, 1, 2 or 3", "0"));
	gain_flip = textToInteger(parser.getOption("--gain_flip", "Flip the gain reference. This is same as MotionCor2's FlipGain. 0, 1 (flip Y == upside down) or 2 (flip X == left to right)", "0"));
	patch_x = textToInteger(parser.getOption("--patch_x", "Patching in X-direction", "1"));
	patch_y = textToInteger(parser.getOption("--patch_y", "Patching in Y-direction", "1"));
	group = textToInteger(parser.getOption("--group_frames", "Average together this many frames before calculating the beam-induced shifts", "1"));
	write_float16  = parser.checkOption("--float16", "Write in half-precision 16 bit floating point numbers (MRC mode 12), instead of 32 bit (MRC mode 0).");
	skip_defect = parser.checkOption("--skip_defect", "Skip hot pixel detection");
	save_noDW = parser.checkOption("--save_noDW", "Save aligned but non dose weighted micrograph");
	max_iter = textToInteger(parser.getOption("--max_iter", "Maximum number of iterations for alignment. Only valid with --use_own", "5"));
	interpolate_shifts = parser.checkOption("--interpolate_shifts", "(EXPERIMENTAL) Interpolate shifts");
	ccf_downsample = textToFloat(parser.getOption("--ccf_downsample", "(EXPERT) Downsampling rate of CC map. default = 0 = automatic based on B factor", "0"));
	if (parser.checkOption("--early_binning", "Do binning before alignment to reduce memory usage. This might dampen signal near Nyquist. (ON by default)"))
		std::cerr << "Since RELION 3.1, --early_binning is on by default. Use --no_early_binning to disable it." << std::endl;

	early_binning = !parser.checkOption("--no_early_binning", "Disable --early_binning");
	if (fabs(bin_factor - 1) < 0.01)
		early_binning = false;

	if (write_float16 && grouping_for_ps <= 0)
		REPORT_ERROR("When writing in float16, you have to write power spectra for CTFFIND.");

	dose_motionstats_cutoff = textToFloat(parser.getOption("--dose_motionstats_cutoff", "Electron dose (in electrons/A2) at which to distinguish early/late global accumulated motion in output statistics", "4."));
	if (ccf_downsample > 1) REPORT_ERROR("--ccf_downsample cannot exceed 1.");

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	// Don't put any output to screen for mpi followers
	verb = (node->isLeader()) ? 1 : 0;

	// Print out MPI info
	printMpiNodesMachineNames(*node);
}

void MotioncorrOwnDevolved::run()
{
	prepareGainReference(node->isLeader());
	MPI_Barrier(MPI_COMM_WORLD); // wait for the leader to write the gain reference

	Micrograph mic(movie_path, fn_gain_reference, bin_factor, eer_upsampling, eer_grouping);
	bool result;
	result = executeOwnMotionCorrection(mic);
	if (result) saveModel(mic);

	MPI_Barrier(MPI_COMM_WORLD);

}

FileName MotioncorrOwnDevolved::getOutputFileNames(FileName fn_mic)
{
	return micrograph_path;
}
