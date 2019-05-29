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
	FileName mydir;
	float myconstant;
	bool do_reset, do_run;
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
		std::cerr << "    --schedule Schedules/test --add operator --type " << SCHEDULE_FLOAT_OPERATOR_PLUS_VAR << " --i iter --i2 iter_step --o iter" << std::endl;
		std::cerr << "    --schedule Schedules/test --add operator --type " << SCHEDULE_FLOAT_OPERATOR_PLUS_CONST << " --i iter --i2 1 --o iter" << std::endl;
		std::cerr << "    --schedule Schedules/test --add operator --type " << SCHEDULE_STRING_OPERATOR_TOUCH_FILE << " --i initial_model" << std::endl;
		std::cerr << "    --schedule Schedules/test --add operator --type " << SCHEDULE_BOOLEAN_OPERATOR_GT_CONST << " --i iter --i2 10 --o is_finished" << std::endl;
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
		if (!exists(mydir))
		{
			std::string command = "mkdir -p " + mydir;
			int res = system(command.c_str());

		}

		// read in schedule if it exists
		schedule.setName(mydir);
		if (exists(mydir + "schedule.star"))
		{
			schedule.read();
			schedule.write(mydir + "schedule.star.bck"); // just save another copy of the starfile ...
		}

		if (do_run)
		{
			PipeLine pipeline;
			pipeline.setName(run_pipeline);
			pipeline.read(DO_LOCK);
			pipeline.write(DO_LOCK);
			schedule.run(pipeline);
		}
		else if (add != "")
		{
			if (add == "variable")
			{
				schedule.addVariable(name, value);
			}
			else if (add == "operator")
			{
				schedule.addOperatorNode(type, input, input2, output);
			}
			else if (add == "job")
			{
				if (input == "exit")
					schedule.addExitNode();
				else
				{
					RelionJob myjob;
					bool dummy;
					myjob.read(input, dummy, true);
					schedule.addJobNode(myjob, input, mode);
				}
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
			//std::cerr << "TODO: make function calls in class!" << std::endl;
			if (isBooleanVariable(set_var))
			{
				if (!(value == "true" || value == "True" || value == "false" || value == "False"))
					REPORT_ERROR("ERROR: invalid value for Boolean variable for --value: " + value);
				bool myval = (value == "true" || value == "True");
				scheduler_bools[set_var].value = myval;
			}
			else if (isFloatVariable(set_var))
			{
				float floatval;
				if (!sscanf(value.c_str(), "%f", &floatval)) // is this a number?
					REPORT_ERROR("ERROR: invalid value for Float variable for --value: " + value);
				scheduler_floats[set_var].value = floatval;
			}
			else if (isStringVariable(set_var))
			{
				scheduler_strings[set_var].value = value;
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
				schedule.nodes[set_mode].mode = value;
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

		schedule.write();
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


