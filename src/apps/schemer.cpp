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

#include "../schemer.h"

#include <src/args.h>
#include "src/strings.h"

class schemer_parameters
{
public:
	FileName mydir, newname;
	float myconstant;
    bool do_reset, do_run, do_abort, has_ori_value;
	int verb;
    std::string add, set_var, set_mode, set_has_started, start_node, current_node, email, type, name, opname, value, ori_value, mode, input, input2, output, output2, boolvar;
	std::string run_pipeline;

	// The actual pipeline
	Scheme scheme;

	// I/O Parser
	IOParser parser;

	void usage()
	{
		parser.writeUsage(std::cerr);
		std::cerr << std::endl;
		std::cerr << " Different ways of using this program: " << std::endl;
		std::cerr << std::endl << " ++ Add a variable (of type float, bool or file): " << std::endl;
		std::cerr << "    --scheme test --add variable --name iter --value 0" << std::endl;
		std::cerr << "    --scheme test --add variable --name is_finished --value False" << std::endl;
		std::cerr << "    --scheme test --add variable --name initial_model --value map.mrc" << std::endl;
		std::cerr << std::endl << " ++ Add an operator node (of type float, bool or file): " << std::endl;
		std::cerr << "    --scheme test --add operator --name INCR_iter--type " << SCHEME_FLOAT_OPERATOR_PLUS << " --i iter --i2 iter_step --o iter" << std::endl;
		std::cerr << "    --scheme test --add operator --name INCR_iter --type " << SCHEME_FLOAT_OPERATOR_PLUS << " --i iter --i2 1 --o iter" << std::endl;
		std::cerr << "    --scheme test --add operator --name CHECK_iter --type " << SCHEME_BOOLEAN_OPERATOR_GT << " --i iter --i2 10 --o is_finished" << std::endl;
		std::cerr << "    --scheme test --add operator --name TOUCH_inimodel --type " << SCHEME_OPERATOR_TOUCH_FILE << " --i initial_model" << std::endl;
		std::cerr << "    --scheme test --add operator --name CHECK_inimodel --type " << SCHEME_BOOLEAN_OPERATOR_FILE_EXISTS << " --i initial_model --o is_finished" << std::endl;
		std::cerr << std::endl << " ++ Add a job node: " << std::endl;
		std::cerr << "    --scheme test --add job --i my_import --mode continue/new" << std::endl;
		std::cerr << "    --scheme test --add job --i exit" << std::endl;
		std::cerr << std::endl << " ++ Add an edge: " << std::endl;
		std::cerr << "    --scheme test --add edge --i inputnodename --o outputnodename" << std::endl;
		std::cerr << std::endl << " ++ Add a fork: " << std::endl;
		std::cerr << "    --scheme test --add fork --i inputnodename --bool boolvar --o outputnodename --o2 outputnodename_if_false" << std::endl;
		std::cerr << "TODO: add information about setting variables etc too!" << std::endl;
		std::cerr << std::endl;
		exit(1);
	}

	void read(int argc, char **argv)
	{

		parser.setCommandLine(argc, argv);

		// Fill the window, but don't show it!
		int gen_section = parser.addSection("General options");
		mydir = parser.getOption("--scheme", "Directory name of the scheme");
		mydir = "Schemes/" + mydir;
		newname = parser.getOption("--copy", "Make a copy of the scheme into this directory", "");
		int add_section = parser.addSection("Add elements to the scheme");
		add = parser.getOption("--add", "Specify category of element to add to the scheme (variable, operator, job, edge or fork)", "");
		type = parser.getOption("--type", "Specify type of that element to add to the scheme", "");
		input = parser.getOption("--i", "Specify input to the element ", "");
		input2 = parser.getOption("--i2", "Specify 2nd input to the element ", "");
		boolvar = parser.getOption("--bool", "Name of boolean variable (for forks)", "");
		output = parser.getOption("--o", "Specify output of the element ", "");
		output2 = parser.getOption("--o2", "Specify 2nd output of the element ", "");
		name = parser.getOption("--name", "Name of the variable, operator or job to be added","");
		value = parser.getOption("--value", "Value of the variable to be added","");
		ori_value = parser.getOption("--original_value", "Original value of the variable to be added","");
		mode = parser.getOption("--mode", "Mode (for jobs): new or continue","");
		int set_section = parser.addSection("Set values of variables in the scheme");
		do_reset = parser.checkOption("--reset", "Reset all variables to their original values");
		do_abort = parser.checkOption("--abort", "Abort a scheme that is running");
		set_var = parser.getOption("--set_var", "Name of a variable to set (using also the --value argument)", "");
		set_mode = parser.getOption("--set_job_mode", "Name of a job whose mode to set (using also the --value argument)", "");
		set_has_started = parser.getOption("--set_has_started", "Name of a job whose has_started variable to set (using also the --value argument)", "");
		current_node = parser.getOption("--set_current_node", "Name of a node to which to set current_node", "");
		int run_section = parser.addSection("Run the schemer within a pipeline");
		do_run = parser.checkOption("--run", "Run the schemer");
		verb = textToInteger(parser.getOption("--verb", "Running verbosity: 0, 1, 2 or 3)", "1"));
		run_pipeline = parser.getOption("--run_pipeline", "Name of the pipeline in which to run this scheme", "default");

		// Someone could give an empty-string ori_value....
		has_ori_value = checkParameter(argc, argv, "--original_value");

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
		scheme.setName(mydir);
		scheme.do_read_only = false;

		if (do_run)
		{
			// Remove the abort signal if it exists
			FileName myabort = scheme.name + RELION_JOB_ABORT_NOW;
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
			scheme.read(DO_LOCK); // lock for the entire duration of the run!!
			scheme.verb = verb;
			scheme.run(pipeline);
		    scheme.write(DO_LOCK);

			return; // exit now
		}

		mktree(mydir);
		if (do_abort)
		{
			scheme.read();
			scheme.abort();
			return;
		}

		// read in scheme if it exists
		if (exists(mydir + "scheme.star"))
		{
			scheme.read(DO_LOCK);
			scheme.write(DONT_LOCK, mydir + "scheme.star.bck"); // just save another copy of the starfile ...
		}

		if (newname != "")
		{
			scheme.copy(newname);
			return;
		}
		else if (add != "")
		{
			if (add == "variable")
			{
				scheme.setVariable(name, value);
				if (has_ori_value) scheme.setOriginalVariable(name, ori_value);
			}
			else if (add == "operator")
			{
				std::string error;
				SchemerOperator myop = scheme.initialiseOperator(type, input, input2, output, error);
				if (error != "") REPORT_ERROR(error);
				else scheme.addOperator(myop, name);
			}
			else if (add == "job")
			{
				RelionJob myjob;
				bool dummy;
				myjob.read(input, dummy, true);
				scheme.addJob(myjob, input, mode);
			}
			else if (add == "edge")
			{
				if (!scheme.checkUniqueInput(input))
				{
					REPORT_ERROR("ERROR: an edge/fork with this input node already exists!");
				}
				scheme.addEdge(input, output);
			}
			else if (add == "fork")
			{
				if (!scheme.checkUniqueInput(input))
				{
					REPORT_ERROR("ERROR: an edge/fork with this input node already exists!");
				}
				scheme.addFork(input, boolvar, output, output2);
			}
		}
		else if (do_reset)
		{
			scheme.reset();
		}
		else if (set_var != "")
		{
			if (isBooleanVariable(set_var))
			{
				if (!(value == "true" || value == "True" || value == "1" ||
					  value == "false" || value == "False" || value == "0"))
					REPORT_ERROR("ERROR: invalid value for Boolean variable for --value: " + value);
				bool myval = (value == "true" || value == "True" || value == "1");
				scheme.setBooleanVariableValue(set_var, myval);

				if (has_ori_value)
				{
					if (!(ori_value == "true" || ori_value == "True"  || ori_value == "1" ||
						  ori_value == "false" || ori_value == "False" || ori_value == "0"))
						REPORT_ERROR("ERROR: invalid value for Boolean variable for --original_value: " + ori_value);
					myval = (ori_value == "true" || ori_value == "True" || ori_value == "1" );
					scheme.setBooleanOriginalVariableValue(set_var, myval);
				}
			}
			else if (isFloatVariable(set_var))
			{
				float floatval;
				if (!sscanf(value.c_str(), "%f", &floatval)) // is this a number?
					REPORT_ERROR("ERROR: invalid value for Float variable for --value: " + value);
				scheme.setFloatVariableValue(set_var, floatval);

				if (has_ori_value)
				{
					if (!sscanf(ori_value.c_str(), "%f", &floatval)) // is this a number?
						REPORT_ERROR("ERROR: invalid value for Float variable for --original_value: " + ori_value);
					scheme.setFloatOriginalVariableValue(set_var, floatval);
				}
			}
			else if (isStringVariable(set_var))
			{
				scheme.setStringVariableValue(set_var, value);
				if (has_ori_value) scheme.setStringOriginalVariableValue(set_var, ori_value);
			}
			else
				REPORT_ERROR("ERROR: unrecognised variable whose value to set: " + set_var);
		}
		else if (set_mode != "")
		{
			if (scheme.isJob(set_mode))
			{
				if (!(value == SCHEME_NODE_JOB_MODE_NEW || value == SCHEME_NODE_JOB_MODE_CONTINUE))
					REPORT_ERROR("ERROR: unvalid value for setting job mode: " + value);
				scheme.jobs[set_mode].mode = value;
			}
			else
				REPORT_ERROR("ERROR: invalid jobname to set mode: " + set_mode);

		}
		else if (set_has_started != "")
		{
			if (scheme.isJob(set_has_started))
			{
                            if (value == "True" || value == "true" || value == "1")
                                scheme.jobs[set_has_started].job_has_started = true;
                            else if (value == "False" || value == "false" || value == "0")
                                scheme.jobs[set_has_started].job_has_started = false;
                            else
                                REPORT_ERROR("ERROR: unvalid value for setting has_started of a job: " + value);
			}
			else
				REPORT_ERROR("ERROR: invalid jobname to set has_started: " + set_has_started);

		}
		else if (current_node != "")
		{
			scheme.current_node = current_node;
		}
		else
			REPORT_ERROR(" ERROR: nothing to do!");

		scheme.write(exists(mydir + "scheme.star")); // only LOCK if the file already existed
	}
};

int main(int argc, char *argv[])
{

	schemer_parameters prm;

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


