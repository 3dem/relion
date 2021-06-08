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

#include "src/pipeliner.h"
#include <src/args.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


class pipeliner_parameters
{
public:
	FileName fn_sched, fn_jobids, fn_options, fn_alias, run_schedule, abort_schedule, add_job_star;
	FileName edit_job_in, edit_job_out, edit_label, edit_value;
	int nr_repeat;
	bool do_check_complete;
	long int minutes_wait, minutes_wait_before, seconds_wait_after, gentle_clean, harsh_clean;
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
		int check_section = parser.addSection("Check job completion options");
		do_check_complete = parser.checkOption("--check_job_completion", "Use this flag to only check whether running jobs have completed");
		int add_section = parser.addSection("Add scheduled jobs options");
		add_job_star = parser.getOption("--addJobFromStar", "Add a job with the type and options as in this job.star to the pipeline", "");
		add_type = parser.getOption("--addJob", "Add a job of this type to the pipeline", "");
		fn_options = parser.getOption("--addJobOptions", "Options for this job (either through --addJobFromStar or --addJob)", "");
		fn_alias = parser.getOption("--setJobAlias", "Set an alias to this job", "");
		int run_section = parser.addSection("Run scheduled jobs options");
		fn_jobids  = parser.getOption("--RunJobs", "Run these jobs", "");
		fn_sched = parser.getOption("--schedule", "Name of the scheduler for running the scheduled jobs", "");
		nr_repeat = textToInteger(parser.getOption("--repeat", "Run the scheduled jobs this many times", "1"));
		minutes_wait = textToInteger(parser.getOption("--min_wait", "Wait at least this many minutes between each repeat", "0"));
		minutes_wait_before = textToInteger(parser.getOption("--min_wait_before", "Wait this many minutes before starting the running the first job", "0"));
		seconds_wait_after = textToInteger(parser.getOption("--sec_wait_after", "Wait this many seconds after a process finishes (workaround for slow IO)", "10"));
		int edit_job_section = parser.addSection("Edit jobs");
		edit_job_in = parser.getOption("--editJob", "Star file of a job to be edited", "");
		edit_job_out = parser.getOption("--editJobOut", "Output star file of the edited job (default is to overwrite input)", "");
		edit_label = parser.getOption("--editOption", "The name of the joboption to be edited. This needs to be present in the input star file.", "");
		edit_value = parser.getOption("--editValue", "The value of the joboption to be set.", "");
		int expert_section = parser.addSection("Expert options");
		pipeline.name = parser.getOption("--pipeline", "Name of the pipeline", "default");
		gentle_clean = textToInteger(parser.getOption("--gentle_clean", "Gentle clean this job", "-1"));
		harsh_clean  = textToInteger(parser.getOption("--harsh_clean", "Harsh clean this job", "-1"));

		// Check for errors in the command-line option
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	void run()
	{
		if (do_check_complete)
		{
			pipeline.read(DO_LOCK);
			pipeline.write(DO_LOCK);
			pipeline.checkProcessCompletion();
		}
		else if (edit_job_in != "")
		{
			if (edit_job_out == "") edit_job_out = edit_job_in;
			if (edit_label == "") REPORT_ERROR("ERROR: no --editOption argument was found for --editJob.");

			RelionJob myjob;
			bool dummy;
			myjob.read(edit_job_in, dummy, true);

			if (myjob.joboptions.find(edit_label) != myjob.joboptions.end())
			{
				if (edit_value == "True" || edit_value == "true") edit_value = "Yes";
				if (edit_value == "False" || edit_value == "false") edit_value = "No";
				myjob.joboptions[edit_label].setString(edit_value);
			}
			else
			{
				REPORT_ERROR("ERROR: could not find the joboption " + edit_label + " in job: " + edit_job_in);
			}

			myjob.write(edit_job_out);

		}
		else if (add_job_star != "")
		{
			pipeline.read(DO_LOCK);
			pipeline.write(DO_LOCK);
			RelionJob job;
			bool is_continue;
			job.read(add_job_star, is_continue, true); // true = do_initialise
			job.is_continue = false;
			int job_num = pipeline.addScheduledJob(job, fn_options);
			if (fn_alias != "")
			{
				std::string error_message;
				if (!pipeline.setAliasJob(job_num, fn_alias, error_message))
				{
					std::cerr << "WARNING: Failed to set the job alias to " << fn_alias << ". The job name remains the default." << std::endl;
				}
			}
		}
		else if (add_type != "")
		{
			pipeline.read(DO_LOCK);
			pipeline.write(DO_LOCK);
			int job_num = pipeline.addScheduledJob(add_type, fn_options);
			if (fn_alias != "")
			{
				std::string error_message;
				if (!pipeline.setAliasJob(job_num, fn_alias, error_message))
				{
					std::cerr << "WARNING: Failed to set the job alias to " << fn_alias << ". The job name remains the default." << std::endl;
				}
			}

		}
		else if (gentle_clean > 0 || harsh_clean > 0)
		{
			pipeline.read(DO_LOCK);
			pipeline.write(DO_LOCK);
			bool found = false;
			for (int i = 0, ilim = pipeline.processList.size(); i < ilim; i++)
			{
//				std::cout << i << " " << pipeline.processList[i].name << std::endl;
				FileName fn_pre, fn_jobnr, fn_post;
				if (!decomposePipelineFileName(pipeline.processList[i].name, fn_pre, fn_jobnr, fn_post))
					continue;

				int job_nr = textToInteger(fn_jobnr.afterFirstOf("job").beforeLastOf("/"));
				if (!(job_nr == gentle_clean || job_nr == harsh_clean))
					continue;

				found = true;
//				std::cout << "Gentle clean " << pipeline.processList[i].name << std::endl;
				std::string error_message;
				if (!pipeline.cleanupJob(i, (job_nr == harsh_clean), error_message))
				{
					std::cerr << "Failed to clean!" << std::endl;
					REPORT_ERROR(error_message);
				}
				break;
			}
			if (!found)
			{
				if (gentle_clean > 0)
					std::cerr << "Could not find job to gentle clean: " << gentle_clean << std::endl;
				else
					std::cerr << "Could not find job harsh clean: " << harsh_clean << std::endl;
			}
		}
		else if (nr_repeat > 0)
		{
			pipeline.read(DO_LOCK);
			pipeline.write(DO_LOCK);
			pipeline.runScheduledJobs(fn_sched, fn_jobids, nr_repeat, minutes_wait, minutes_wait_before, seconds_wait_after);
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
		return RELION_EXIT_FAILURE;
	}
	return RELION_EXIT_SUCCESS;
}
