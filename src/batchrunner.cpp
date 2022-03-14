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

#include "batchrunner.h"



void BatchRunner::read(int argc, char **argv)
{
	parser.setCommandLine(argc, argv);

	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "Input text file with the jobs to execute");
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));
	is_continue = parser.checkOption("--continue", "Only execute those commands that were not done yet");
	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void BatchRunner::usage()
{
	parser.writeUsage(std::cout);
}

void BatchRunner::initialise()
{

	FileName my_fn_in = fn_in;
	std::ifstream fh_in;
	if (is_continue && exists(fn_in + "-CONTINUE")) my_fn_in = fn_in + "-CONTINUE";
	fh_in.open(my_fn_in.c_str());
	if (!fh_in) REPORT_ERROR("Cannot open " + my_fn_in);

	number_of_commands = 0;
	FileName line;
	mylines.clear();
	done_lines.clear();
	while (!fh_in.eof())
	{
		getline(fh_in, line);
		line.replaceAllSubstrings("MPI_NEWLINE","\n");
		mylines.push_back(line);
		done_lines.push_back(false);
		if (!line.contains("MPI_BARRIER"))
		{
			number_of_commands++;
		}
	}
	fh_in.close();

}

void BatchRunner::executeCommand(std::string command, int rank)
{
	// Abort through the pipeline_control system
	if (pipeline_control_check_abort_job())
		exit(RELION_EXIT_ABORTED);

	if (system(command.c_str()))
	{
		std::string errormsg = "ERROR: there was an error in executing: " + command;
		if (rank >= 0)
		{
			errormsg += " on follower: " + integerToString(rank);
		}
		REPORT_ERROR(errormsg);
	}

}

void BatchRunner::writeContinuationFile()
{
	std::ofstream  fh;
	FileName fn_out = fn_in + "-CONTINUE";
	fh.open((fn_out).c_str(), std::ios::out);
	if (!fh) REPORT_ERROR( (std::string)"ERROR: cannot write to file: " + fn_out);

	for (int iline = 0; iline < mylines.size(); iline++)
	{
		if (!done_lines[iline]) fh << mylines[iline] << std::endl;
	}
	fh.close();

}

void BatchRunner::run()
{

	int barstep, done = 0;
	if (verb > 0)
	{
		std::cout << " Executing " << number_of_commands << " commands ..." << std::endl;
		init_progress_bar(number_of_commands);
		barstep = XMIPP_MAX(1, number_of_commands / 60);
	}

    for (int iline = 0; iline < mylines.size(); iline++)
    {
    	FileName line=mylines[iline];

    	// ignore lines with MPI_BARRIER in sequential version of the code
    	std::string::size_type loc = line.find("MPI_BARRIER", 0);

        if (loc == std::string::npos)
        {
        	executeCommand(line);
        	done_lines[iline] = true;
        	writeContinuationFile();
            done++;
        }

		if (verb > 0 && done % barstep == 0) progress_bar(done);

    }

    // Remove the -CONTINUE file

	if (verb > 0)
	{
		progress_bar(number_of_commands);
		std::cout << " done! " << std::endl;
	}

}
