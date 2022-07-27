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
#include "src/denoise_tomogram_runner.h"

void DenoiseTomogramRunner::read(int argc, char **argv, int rank)
{
	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "STAR file with all input tomograms");
	fn_out = parser.getOption("--o", "Directory, where all output files will be stored", "DenoiseTomo/");
    fn_denoisingwrapper_exe = parser.getOption("--wrapper_executable", "Alister Burt's wrapper script to call Denoising Programs (default is set through $RELION_DENOISING_WRAPPER_EXECUTABLE)", "");

    int cryocare_train_section = parser.addSection("CryoCARE denoising model training options");
    do_cryocare_train = parser.checkOption("--do_cryocare_train", "Use cryoCARE to generate a denoising model");
    tomograms_for_training = parser.getOption("--tomograms_for_training", "Tomograms to train denoising model on. Specify tomograms using rlnTomoName and separating tomograms with a colon, e.g.: TS_01:TS_02", "");
    number_training_subvolumes = textToInteger(parser.getOption("--number_training_subvolumes", "Number of subvolumes to be extracted per training tomogram. Number of normalisation subvolumes will be 10% of this value.", "1200"));
    subvolume_dimensions = textToInteger(parser.getOption("--subvolume_dimensions", "Dimensions of the subvolumes to be extracted for training in pixels. Do not choose a number smaller than 64.", "72"));

    int cryocare_predict_section = parser.addSection("CryoCARE denoised tomogram prediction options");
    do_cryocare_predict = parser.checkOption("--do_cryocare_predict", "Use cryoCARE to generate denoised tomograms using a given denoising model.");
    care_denoising_model = parser.getOption("--care_denoising_model", "Path to cryoCARE generated denoising model (.tar.gz)", "");
    ntiles_x = textToInteger(parser.getOption("--ntiles_x", "Number of tiles to use in denoised tomogram generation (X dimension)", ""));
    ntiles_y = textToInteger(parser.getOption("--ntiles_y", "Number of tiles to use in denoised tomogram generation (Y dimension)", ""));
    ntiles_z = textToInteger(parser.getOption("--ntiles_z", "Number of tiles to use in denoised tomogram generation (Z dimension)", ""));

    int exp_section = parser.addSection("Expert options");
    other_wrapper_args  = parser.getOption("--other_wrapper_args", "Additional command-line arguments that will be passed onto the wrapper.", "");

    // Initialise verb for non-parallel execution
	verb = 1;
	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}

void DenoiseTomogramRunner::usage()
{
	parser.writeUsage(std::cout);
}

void DenoiseTomogramRunner::initialise()
{
	// Get the denoising wrapper executable
	if (fn_denoisingwrapper_exe == "")
	{
		char *penv;
		penv = getenv("RELION_DENOISING_WRAPPER_EXECUTABLE");
		if (penv != NULL)
			fn_denoisingwrapper_exe = (std::string)penv;
	}

    int i = 0;
    if (do_cryocare_train) i++;
    if (do_cryocare_predict) i++;

    if (i != 1)
    {
        REPORT_ERROR("ERROR: you should (only) select ONE of the methods: cryoCARE:train (--do_cryocare_train), cryoCARE:predict (--do_cryocare_train)");
    }
    // Make sure fn_out ends with a slash
    if (fn_out[fn_out.length()-1] != '/')
    	fn_out += "/";

    // Check if this is a TomographyExperiment starfile, if not raise an error
    if (!tomogramSet.read(fn_in, 1))
    {
        REPORT_ERROR("ERROR: the input file is not a valid tomograms.star file");
    }

	idx_tomograms_all.clear();
	idx_tomograms.clear();
	bool warned = false;
	for (long int itomo = 0; itomo < tomogramSet.size(); itomo++)
	{
        	FileName fn_star;
        	tomogramSet.globalTable.getValue(EMDL_TOMO_TILT_SERIES_STARFILE, fn_star, itomo);
        	FileName fn_newstar = getOutputFileWithNewUniqueDate(fn_star, fn_out);
        	tomogramSet.globalTable.setValue(EMDL_TOMO_TILT_SERIES_STARFILE, fn_newstar, itomo);
		idx_tomograms.push_back(itomo);
        }
}

void DenoiseTomogramRunner::run()
{

    int barstep;
    if (verb > 0)
    {
        std::cout << " Running Denoising Wrapper ..." << std::endl;
    }

    std::vector<std::string> alltomonames;
    if (do_cryocare_train)
    {
       executeCryoCARETrainWrapper();
    }
    else if (do_cryocare_predict)
    {
       executeCryoCAREPredictWrapper();
    }

    // Abort through the pipeline_control system
    if (pipeline_control_check_abort_job())
        exit(RELION_EXIT_ABORTED);	

}


void DenoiseTomogramRunner::executeCryoCARETrainWrapper(int rank)
{

   std::string command = fn_denoisingwrapper_exe + " ";
    // Make sure the methods are the first argument to the program!
    command += " cryoCARE:train";
    command += " --tilt-series-star-file " + fn_in;
    command += " --output-directory " + fn_out;
    command += " --training-tomograms " + tomograms_for_training;
    command += " --number-training-subvolumes " + integerToString(number_training_subvolumes);
    command += " --subvolume-dimensions " + integerToString(subvolume_dimensions);

    if (other_wrapper_args.length() > 0)
        command += " " + other_wrapper_args;

    if (system(command.c_str()))
		std::cerr << "WARNING: there was an error in executing: " << command << std::endl;

}

void DenoiseTomogramRunner::executeCryoCAREPredictWrapper(int rank)
{
    std::string command = fn_denoisingwrapper_exe + " ";
    // Make sure the methods are the first argument to the program!
    command += " cryoCARE:predict";
    command += " --tilt-series-star-file " + fn_in;
    command += " --output-directory " + fn_out;
    command += " --model-name " + care_denoising_model;
if (ntiles_x > 0 && ntiles_y > 0 && ntiles_z > 0)     
    {
    command += " --n-tiles " + integerToString(ntiles_x) + " " + integerToString(ntiles_y) + " " + integerToString(ntiles_z) + " ";
    }

    if (other_wrapper_args.length() > 0)
        command += " " + other_wrapper_args;

    if (system(command.c_str()))
		std::cerr << "WARNING: there was an error in executing: " << command << std::endl;

}

