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

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include "src/pipeliner.h"
#include <src/args.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


class pipeliner_parameters
{
public:
	FileName fn_sched, fn_jobids, fn_options;
	int nr_repeat;
	long int minutes_wait, minutes_wait_before;
	std::string add_type;

	// The actual pipeline
	PipeLine pipeline;

	// I/O Parser
	IOParser parser;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{

		parser.setCommandLine(argc, argv);

		// Fill the window, but don't show it!
		int add_section = parser.addSection("Add scheduled jobs options");
		add_type = parser.getOption("--addJob", "Add a job of this type to the pipeline");
		fn_options = parser.getOption("--addJobOptions", "Options for this job");
		int run_section = parser.addSection("Run scheduled jobs options");
		fn_jobids  = parser.getOption("--RunJobs", "Run these jobs", "");
		fn_sched = parser.getOption("--schedule", "Name of the scheduler for running the scheduled jobs", "");
		nr_repeat = textToInteger(parser.getOption("--repeat", "Repeat the scheduled jobs this many times", "1"));
		minutes_wait = textToInteger(parser.getOption("--min_wait", "Wait at least this many minutes between each repeat", "0"));
		minutes_wait_before = textToInteger(parser.getOption("--min_wait_before", "Wait this many minutes before starting the running the first job", "0"));
		int expert_section = parser.addSection("Expert options");
		pipeline.name = parser.getOption("--pipeline", "Name of the pipeline", "default");

    	// Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	}

	void run()
	{

		pipeline.read(DO_LOCK);
		pipeline.write(DO_LOCK);
		if (add_type != "")
		{
			pipeline.addScheduledJob(add_type, fn_options);

		}
		else if (nr_repeat > 0)
		{
			pipeline.runScheduledJobs(fn_sched, fn_jobids, nr_repeat, minutes_wait, minutes_wait_before);
		}

	}

};

int main(int argc, char *argv[])
{

	pipeliner_parameters prm;

	try
    {

		prm.read(argc, argv);

		prm.run();

    }
    catch (RelionError XE)
    {
        std::cerr << XE;
        exit(1);
    }
    return 0;

  }
