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

#ifndef BATCHRUNNER_H_
#define BATCHRUNNER_H_

#include "src/pipeline_control.h"
#include "src/filename.h"
#include "src/args.h"
#include "src/error.h"
#include "src/time.h"

class BatchRunner
{
public:

	// I/O Parser
	IOParser parser;

	// FileName and file handler for the list of commands to be executed
	FileName fn_in;
	std::vector<std::string> mylines;
	std::vector<bool> done_lines;

	// Is this a continuation job?
	// This is not entirely clean: the -CONTINUE file will keep having some jobs that are in fact finished, even after job completes...
	// This is hard because of how MPI code is implemented. Hopefully still useful to not to have to repeat the bulk of the work
	bool is_continue;

	// Number of commands
	int number_of_commands;

	// Verbosity
	int verb;

public:

	// Read command line arguments
	void read(int argc, char **argv);

	// Print usage instructions
	void usage();

	// Initialise some general stuff after reading
	void initialise();

	void executeCommand(std::string command, int rank = -1);

	void writeContinuationFile();

	void run();


};



#endif /* BATCHRUNNER_H_ */
