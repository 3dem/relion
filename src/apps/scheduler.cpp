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

#include "src/scheduler.h"
#include <src/args.h>
#include "src/strings.h"

class scheduler_parameters
{
public:
	FileName mydir, newname;
	float myconstant;
	bool do_reset, do_run;
	int verb;
	std::string add, set_var, set_mode, start_node, current_node, email, type, name, value, mode, input, input2, output, output2, boolvar;
	std::string run_pipeline;

	// The actual pipeline
	Schedule schedule;

	// I/O Parser
	IOParser parser;

	void usage()
	{
		parser.writeUsage(std::cerr);
		std::cerr << std::endl;
		std::cerr << " Different ways of using this program: " << std::endl;
		std::cerr << std::endl << " ++ Add a variable (of type float, bool or file): " << std::endl;
		std::cerr << "    --schedule Schedules/test --add variable --name iter --value 0" << std::endl;
		std::cerr << "    --schedule Schedules/test --add variable --name is_finished --value False" << std::endl;
		std::cerr << "    --schedule Schedules/test --add variable --name initial_model --value map.mrc" << std::endl;
		std::cerr << std::endl << " ++ Add an operator node (of type float, bool or file): " << std::endl;
		std::cerr << "    --schedule Schedules/test --add operator --type " << SCHEDULE_FLOAT_OPERATOR_PLUS << " --i iter --i2 iter_step --o iter" << std::endl;
		std::cerr << "    --schedule Schedules/test --add operator --type " << SCHEDULE_FLOAT_OPERATOR_PLUS << " --i iter --i2 1 --o iter" << std::endl;
		std::cerr << "    --schedule Schedules/test --add operator --type " << SCHEDULE_STRING_OPERATOR_TOUCH_FILE << " --i initial_model" << std::endl;
		std::cerr << "    --schedule Schedules/test --add operator --type " << SCHEDULE_BOOLEAN_OPERATOR_GT << " --i iter --i2 10 --o is_finished" << std::endl;
		std::cerr << "    --schedule Schedules/test --add operator --type " << SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS << " --i initial_model --o is_finished" << std::endl;
		std::cerr << std::endl << " ++ Add a job node: " << std::endl;
		std::cerr << "    --schedule Schedules/test --add job --i my_import --mode continue/new/overwrite" << std::endl;
		std::cerr << "    --schedule Schedules/test --add job --i exit" << std::endl;
		std::cerr << std::endl << " ++ Add an edge: " << std::endl;
		std::cerr << "    --schedule Schedules/test --add edge --i inputnodename --o outputnodename" << std::endl;
		std::cerr << std::endl << " ++ Add a fork: " << std::endl;
		std::cerr << "    --schedule Schedules/test --add fork --i inputnodename --bool boolvar --o outputnodename --o2 outputnodename_if_false" << std::endl;
		std::cerr << "TODO: add information about setting variables etc too!" << std::endl;
		std::cerr << std::endl;
		exit(1);
	}

	void read(int argc, char **argv)
	{

		parser.setCommandLine(argc, argv);

		// Fill the window, but don't show it!
		int gen_section = parser.addSection("General options");
		mydir = parser.getOption("--schedule", "Directory name of the schedule");
		newname = parser.getOption("--copy", "Make a copy of the schedule into this directory");
		int add_section = parser.addSection("Add elements to the schedule");
		add = parser.getOption("--add", "Specify category of element to add to the schedule (variable, operator, job, edge or fork)", "");
		type = parser.getOption("--type", "Specify type of that element to add to the schedule", "");
		input = parser.getOption("--i", "Specify input to the element ", "");
		input2 = parser.getOption("--i2", "Specify 2nd input to the element ", "");
		boolvar = parser.getOption("--bool", "Name of boolean variable (for forks)", "");
		output = parser.getOption("--o", "Specify output of the element ", "");
		output2 = parser.getOption("--o2", "Specify 2nd output of the element ", "");
		name = parser.getOption("--name", "Name of the variable or job to be added","");
		value = parser.getOption("--value", "Value of the variable to be added","");
		mode = parser.getOption("--mode", "Mode (for jobs): new, overwrite or continue","");
		int set_section = parser.addSection("Set values of variables in the schedule");
		do_reset = parser.checkOption("--reset", "Reset all variables to their original values");
		set_var = parser.getOption("--set_var", "Name of a variable to set (using also the --value argument)", "");
		set_mode = parser.getOption("--set_job_mode", "Name of a job whose mode to set (using also the --value argument)", "");
		current_node = parser.getOption("--set_current_node", "Name of a node to which to set current_node", "");
		start_node = parser.getOption("--set_start_node", "Name of a node to which to set original_start_node", "");
		email = parser.getOption("--set_email", "Email address to send messages about Scheduler status", "");
		int run_section = parser.addSection("Run the scheduler within a pipeline");
		do_run = parser.checkOption("--run", "Run the scheduler");
		verb = textToInteger(parser.getOption("--verb", "Running verbosity: 0, 1, 2 or 3)", "1"));
		run_pipeline = parser.getOption("--run_pipeline", "Name of the pipeline in which to run this schedule", "default");

		// Check for errors in the command-line option
		if (argc==1)
			usage();
		else if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	void run()
	{

		// Make sure mydir ends with a slash, and that it exists
		if (mydir[mydir.length()-1] != '/') mydir += "/";
		schedule.setName(mydir);
		schedule.do_read_only = false;

		if (do_run)
		{
			// Remove the abort signal if it exists
			FileName myabort = schedule.name + RELION_JOB_ABORT_NOW;
			if (exists(myabort)) std::remove(myabort.c_str());

			PipeLine pipeline;
			pipeline.setName(run_pipeline);

			if (exists(pipeline.name + "_pipeline.star"))
			{
				std::string lock_message = "mainGUI constructor";
				pipeline.read(DO_LOCK, lock_message);
				// With the locking system, each read needs to be followed soon with a write
				pipeline.write(DO_LOCK);
			}
			else
			{
				pipeline.write();
			}
			schedule.read(DO_LOCK); // lock for the entire duration of the run!!
			schedule.verb = verb;
			schedule.run(pipeline);
		    schedule.write(DO_LOCK);

			return; // exit now
		}

		if (!exists(mydir))
		{
			std::string command = "mkdir -p " + mydir;
			int res = system(command.c_str());

		}


		// read in schedule if it exists
		if (exists(mydir + "schedule.star"))
		{
			schedule.read(DO_LOCK);
			schedule.write(DONT_LOCK, mydir + "schedule.star.bck"); // just save another copy of the starfile ...
		}

		if (newname != "")
		{
			schedule.copy(newname);
			return;
		}
		else if (add != "")
		{
			if (add == "variable")
			{
				schedule.setVariable(name, value);
			}
			else if (add == "operator")
			{
				std::string error;
				SchedulerOperator myop = schedule.initialiseOperator(type, input, input2, output, error);
				if (error != "") REPORT_ERROR(error);
				else schedule.addOperator(myop);
			}
			else if (add == "job")
			{
				RelionJob myjob;
				bool dummy;
				myjob.read(input, dummy, true);
				schedule.addJob(myjob, input, mode);
			}
			else if (add == "edge")
			{
				schedule.addEdge(input, output);
			}
			else if (add == "fork")
			{
				schedule.addFork(input, boolvar, output, output2);
			}
		}
		else if (do_reset)
		{
			schedule.reset();
		}
		else if (set_var != "")
		{
			if (isBooleanVariable(set_var))
			{
				if (!(value == "true" || value == "True" || value == "false" || value == "False"))
					REPORT_ERROR("ERROR: invalid value for Boolean variable for --value: " + value);
				bool myval = (value == "true" || value == "True");
				schedule.setBooleanVariableValue(set_var, myval);
			}
			else if (isFloatVariable(set_var))
			{
				float floatval;
				if (!sscanf(value.c_str(), "%f", &floatval)) // is this a number?
					REPORT_ERROR("ERROR: invalid value for Float variable for --value: " + value);
				schedule.setFloatVariableValue(set_var, floatval);
			}
			else if (isStringVariable(set_var))
			{
				schedule.setStringVariableValue(set_var, value);
			}
			else
				REPORT_ERROR("ERROR: unrecognised variable whose value to set: " + set_var);
		}
		else if (set_mode != "")
		{
			if (schedule.isJob(set_mode))
			{
				if (!(value == SCHEDULE_NODE_JOB_MODE_NEW || value == SCHEDULE_NODE_JOB_MODE_CONTINUE || value == SCHEDULE_NODE_JOB_MODE_OVERWRITE))
					REPORT_ERROR("ERROR: unvalid option for job mode: " + value);
				schedule.jobs[set_mode].mode = value;
			}
			else
				REPORT_ERROR("ERROR: invalid jobname to set mode: " + set_mode);

		}
		else if (current_node != "")
		{
			schedule.current_node = current_node;
		}
		else if (start_node != "")
		{
			schedule.original_start_node = start_node;
		}
		else if (email != "")
		{
			schedule.email_address = email;
		}
		else
			REPORT_ERROR(" ERROR: nothing to do!");

		schedule.write(exists(mydir + "schedule.star")); // only LOCK if the file already existed
	}
};

int main(int argc, char *argv[])
{

	scheduler_parameters prm;

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


