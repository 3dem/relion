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
#include "src/pipeline_jobs.h"

std::vector<Node> getOutputNodesRefine(std::string outputname, int iter, int K, int dim, int nr_bodies, bool do_movies, bool do_also_rot)
{
	std::vector<Node> result;

	if (dim < 2 || dim > 3)
		REPORT_ERROR("getOutputNodesRefine ERROR: invalid dim value");

	FileName fn_out;
	if (iter < 0)
	{
		// 3D auto-refine
		fn_out = outputname;
	}
	else
	{
		// 2D or 3D classification
		fn_out.compose(outputname+"_it", iter, "", 3);
	}

	// Data and model.star files
	int node_type = (do_movies) ? NODE_MOVIE_DATA : NODE_PART_DATA;
	Node node1(fn_out + "_data.star", node_type);
	result.push_back(node1);

	if (!do_movies || do_also_rot)
	{
		if (iter > 0)
		{
			Node node2(fn_out + "_model.star", NODE_MODEL);
			result.push_back(node2);
		}

		// For 3D classification or 3D auto-refine, also use individual 3D maps as outputNodes
		if (dim == 3)
		{
			FileName fn_tmp;
			for (int iclass = 0; iclass < K; iclass++)
			{
				if (nr_bodies > 1)
					fn_tmp.compose(fn_out+"_body", iclass+1, "mrc", 3);
				else
					fn_tmp.compose(fn_out+"_class", iclass+1, "mrc", 3);

				Node node3(fn_tmp, NODE_3DREF);
				result.push_back(node3);
			}
		}

		// For auto-refine: also output the half1_class001_unfil.mrc map
		if (iter < 0)
		{
			Node node4(fn_out+"_half1_class001_unfil.mrc", NODE_HALFMAP);
			result.push_back(node4);
		}
	}
	return result;

}

// Any constructor
JobOption::JobOption(std::string _label, std::string _default_value, std::string _helptext)
{
	initialise(_label, _default_value, _helptext);
	joboption_type = JOBOPTION_ANY;
}

// FileName constructor
JobOption::JobOption(std::string _label, std::string  _default_value, std::string _pattern, std::string _directory, std::string _helptext)
{
	initialise(_label, _default_value, _helptext);
	joboption_type = JOBOPTION_FILENAME;
	pattern = _pattern;
	directory = _directory;
}

// InputNode constructor
JobOption::JobOption(std::string _label, int _nodetype, std::string _default_value, std::string _pattern, std::string _helptext)
{
	initialise(_label, _default_value, _helptext);
	joboption_type = JOBOPTION_INPUTNODE;
	pattern = _pattern;
	node_type = _nodetype;
}

// Radio constructor
JobOption::JobOption(std::string _label, int _radio_menu, int ioption,  std::string _helptext)
{
	radio_menu = _radio_menu;
	std::string defaultval;
	if (radio_menu == RADIO_SAMPLING)
		defaultval = std::string(job_sampling_options[ioption]);
	else if (radio_menu == RADIO_NODETYPE)
		defaultval = std::string(job_nodetype_options[ioption]);
	else
		REPORT_ERROR("BUG: unrecognised radio_menu type");
	initialise(_label, defaultval, _helptext);
	joboption_type = JOBOPTION_RADIO;

}

// Boolean constructor
JobOption::JobOption(std::string _label, bool _boolvalue, std::string _helptext)
{
	std::string _default_value = (_boolvalue) ? "Yes" : "No";
	initialise(_label, _default_value, _helptext);
	joboption_type = JOBOPTION_BOOLEAN;
}

// Slider constructor
JobOption::JobOption(std::string _label, float _default_value, float _min_value, float _max_value, float _step_value, std::string _helptext)
{
	initialise(_label, floatToString(_default_value), _helptext);
	joboption_type = JOBOPTION_SLIDER;
	min_value = _min_value;
	max_value = _max_value;
	step_value = _step_value;
}

void JobOption::clear()
{
	label = value = default_value = helptext = label_gui = pattern = directory = "";
	joboption_type = JOBOPTION_UNDEFINED;
	node_type = min_value = max_value = step_value = radio_menu = 0.;
}

void JobOption::initialise(std::string _label, std::string _default_value, std::string _helptext)
{
	label = label_gui = _label;
	value = default_value = _default_value;
	helptext = _helptext;
}

    // Get a string value
std::string JobOption::getString()
{
	return value;
}

// Get a numbered value
float JobOption::getNumber()
{
	if (joboption_type != JOBOPTION_SLIDER)
	{
		std::cerr << " joboption_type= " << joboption_type << " label= " << label << " value= " << value << std::endl;
		REPORT_ERROR("ERROR: this jobOption does not return a number: " + label);
	}
	else
	    return textToFloat(value);
}

// Get a boolean value
bool JobOption::getBoolean()
{
	if (joboption_type != JOBOPTION_BOOLEAN)
	{
		std::cerr << " joboption_type= " << joboption_type << " label= " << label << " value= " << value << std::endl;
		REPORT_ERROR("ERROR: this jobOption does not return a boolean: " + label);
	}
	else
		return (value == "Yes");
}

bool JobOption::readValue(std::ifstream& in)
{
    if (label != "")
    {
		// Start reading the ifstream at the top
		in.clear(); // reset eof if happened...
    	in.seekg(0, std::ios::beg);
		std::string line;
		while (getline(in, line, '\n'))
		{
			if (line.rfind(label) == 0)
			{
				// found my label
				int equalsigns = line.rfind("==");
				value = line.substr(equalsigns + 3, line.length() - equalsigns - 3);
				return true;
			}
		}
    }
    return false;
}

void JobOption::writeValue(std::ostream& out)
{
	out << label << " == " << value << std::endl;
}


bool RelionJob::read(std::string fn, bool &_is_continue, bool do_initialise)
{

	// If fn is empty, use the hidden name
	FileName myfilename = (fn=="") ? hidden_name : fn;

	std::ifstream fh;
	fh.open((myfilename+"run.job").c_str(), std::ios_base::in);
    if (fh.fail())
    	return false;
    else
    {
		std::string line;

		// Get job type from first line
		getline(fh, line, '\n');
		size_t idx = line.find("==");
		idx++;
		type = (int)textToFloat((line.substr(idx+1,line.length()-idx)).c_str());
		if (!(type >= 0 && type < NR_BROWSE_TABS))
			REPORT_ERROR("ERROR: cannot find job type in " + myfilename + "run.job");

		// Get is_continue from second line
		getline(fh, line, '\n');
		if (line.rfind("is_continue == true") == 0)
			is_continue = true;
		else
			is_continue = false;
		_is_continue = is_continue;

		if (do_initialise)
			initialise(type);

		// Read in all the stored options
		bool read_all = true;
		for (std::map<std::string,JobOption>::iterator it=joboptions.begin(); it!=joboptions.end(); ++it)
		{
			if (!(it->second).readValue(fh))
				read_all = false;
		}

		return read_all;
    }

	fh.close();

}

void RelionJob::write(std::string fn)
{
	// If fn is empty, use the hidden name
	FileName myfilename = (fn=="") ? hidden_name : fn;

	std::ofstream fh;
	fh.open((myfilename+"run.job").c_str(), std::ios::out);
    if (!fh)
    	REPORT_ERROR("ERROR: Cannot write to file: " + myfilename + "run.job");

    // Write the job type
    fh << "job_type == " << type << std::endl;

    // is_continue flag
    if (is_continue)
    	fh << "is_continue == true" << std::endl;
    else
    	fh << "is_continue == false" << std::endl;

    for (std::map<std::string,JobOption>::iterator it=joboptions.begin(); it!=joboptions.end(); ++it)
	{
		(it->second).writeValue(fh);
	}

	fh.close();


}

bool RelionJob::saveJobSubmissionScript(std::string newfilename, std::string outputname, std::vector<std::string> commands, std::string &error_message)
{

	// Open the standard job submission file
	FileName fn_qsub = joboptions["qsubscript"].getString();


	std::ofstream fo;
	std::ifstream fh;
	fh.open(fn_qsub.c_str(), std::ios_base::in);
	fo.open(newfilename.c_str(), std::ios::out);
    if (fh.fail())
    {
    	error_message = "Error reading template submission script in: " + fn_qsub;
    	return false;
    }
    else if (fo.fail())
    {
    	error_message = "Error writing to job submission script in: " + newfilename;
    	return false;
    }
    else
    {
    	int nmpi = (joboptions.find("nr_mpi") != joboptions.end()) ? joboptions["nr_mpi"].getNumber() : 1;
    	int nthr = (joboptions.find("nr_threads") != joboptions.end()) ? joboptions["nr_threads"].getNumber() : 1;
    	int ncores = nmpi * nthr;
    	int ndedi = joboptions["min_dedicated"].getNumber();
    	float fnodes = (float)ncores / (float)ndedi;
    	int nnodes = CEIL(fnodes);
    	if (fmod(fnodes, 1) > 0)
    	{
    		std:: cout << std::endl;
    		std::cout << " Warning! You're using " << nmpi << " MPI processes with " << nthr << " threads each (i.e. " << ncores << " cores), while asking for " << nnodes << " nodes with " << ndedi << " cores." << std::endl;
    		std::cout << " It is more efficient to make the number of cores (i.e. mpi*threads) a multiple of the minimum number of dedicated cores per node " << std::endl;
    	}

		fh.clear(); // reset eof if happened...
    	fh.seekg(0, std::ios::beg);
		std::string line;
		std::map<std::string, std::string> replacing;
		replacing["XXXmpinodesXXX"] = floatToString(nmpi);
		replacing["XXXthreadsXXX"] = floatToString(nthr);
		replacing["XXXcoresXXX"] = floatToString(ncores);
		replacing["XXXdedicatedXXX"] = floatToString(ndedi);
		replacing["XXXnodesXXX"] = floatToString(nnodes);
		replacing["XXXnameXXX"] = outputname;
		replacing["XXXerrfileXXX"] = outputname + "run.err";
		replacing["XXXoutfileXXX"] = outputname + "run.out";
		replacing["XXXqueueXXX"] = joboptions["queuename"].getString();
		if (joboptions.find("qsub_extra1") != joboptions.end())
			replacing["XXXextra1XXX"] = joboptions["qsub_extra1"].getString();
		if (joboptions.find("qsub_extra2") != joboptions.end())
			replacing["XXXextra2XXX"] = joboptions["qsub_extra2"].getString();

		while (getline(fh, line, '\n'))
		{

			// Replace all entries in the replacing map
			for (std::map<std::string,std::string>::iterator it=replacing.begin(); it!=replacing.end(); ++it)
			{

				std::string from = it->first;
				std::string to = it->second;

				// Replace all instances of the string on the line
				size_t start_pos = 0;
				while((start_pos = line.find(from, start_pos)) != std::string::npos)
				{
				         line.replace(start_pos, from.length(), to);
				         start_pos += to.length();
				}
			}

			if (line.find("XXXcommandXXX") == std::string::npos)
			{
				fo << line << std::endl;;
			}
			else
			{
				// Append the commands
				std::string ori_line = line;
				for (int icom = 0; icom < commands.size(); icom++)
				{
					// For multiple relion mpi commands: add multiple lines from the XXXcommandXXX template
					if ((commands[icom]).find("relion_") != std::string::npos &&
							((commands[icom]).find("_mpi`") != std::string::npos || nmpi==1) ) // if there are no MPI programs, then still use XXXcommandXXX once
					{
						std::string from = "XXXcommandXXX";
						std::string to = commands[icom];
						line.replace(line.find(from), from.length(), to);
						fo << line << std::endl;
						line = ori_line;
					}
					else
					{
						// Just add the sequential command
						fo << commands[icom] << std::endl;
					}
				}
			}

		}

		fo << std::endl;

    	fo.close();
		fh.close();
    }

	return true;

}

void RelionJob::initialisePipeline(std::string &outputname, std::string defaultname, int job_counter)
{

	outputNodes.clear();
	inputNodes.clear();

	if (outputname == "") // for continue jobs, use the same outputname
	{
		if (job_counter < 1000)
			outputname = defaultname + "/job" + integerToString(job_counter, 3) + "/";
		else
			outputname = defaultname + "/job" + integerToString(job_counter) + "/";
	}

	outputName = outputname;

}

bool RelionJob::prepareFinalCommand(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, std::string &error_message)
{
	int nr_mpi;

	// Create output directory if the outname contains a "/"
	if (do_makedir)
	{
		int last_slash = outputname.rfind("/");
		if (last_slash < outputname.size())
		{
			std::string dirs = outputname.substr(0, last_slash);
			std::string makedirs = "mkdir -p " + dirs;
			int res = system(makedirs.c_str());
		}
	}

	// Prepare full mpi commands or save jobsubmission script to disc
	if (joboptions["do_queue"].getBoolean() && do_makedir)
	{
		// Make the submission script and write it to disc
		std::string output_script = outputname + "run_submit.script";
		if (!saveJobSubmissionScript(output_script, outputname, commands, error_message))
			return false;
		final_command = joboptions["qsub"].getString() + " " + output_script + " &";
	}
	else
	{
		// If there are multiple commands, then join them all on a single line (final_command)
		// Also add mpirun in front of those commands that have relion_ and _mpi` in it (if no submission via the queue is done)
		std::string one_command;
		final_command = "";
		for (size_t icom = 0; icom < commands.size(); icom++)
		{
			// Is this a relion mpi program?
			nr_mpi = (joboptions.find("nr_mpi") != joboptions.end()) ? joboptions["nr_mpi"].getNumber() : 1;
			if (nr_mpi > 1 &&
					(commands[icom]).find("_mpi`") != std::string::npos &&
					(commands[icom]).find("relion_") != std::string::npos)
				one_command = "mpirun -n " + floatToString(nr_mpi) + " " + commands[icom] ;
			else
				one_command = commands[icom];

			// Save stdout and stderr to a .out and .err files
			// But only when a re-direct '>' is NOT already present on the command line!
			if (std::string::npos == commands[icom].find(">"))
				one_command += " >> " + outputname + "run.out 2>> " + outputname + "run.err";
			final_command += one_command;
			if (icom == commands.size() - 1)
				final_command += " & "; // end by putting composite job in the background
			else
				final_command += " && "; // execute one command after the other...
		}
	}

	char * my_warn = getenv ("RELION_ERROR_LOCAL_MPI");
	int my_nr_warn = (my_warn == NULL) ? DEFAULTWARNINGLOCALMPI : textToInteger(my_warn);

	if (nr_mpi > my_nr_warn && !joboptions["do_queue"].getBoolean())
	{
		error_message = "You're submitting a local job with " + floatToString(nr_mpi) + " parallel MPI processes. That's more than allowed by the RELION_ERROR_LOCAL_MPI environment variable.";
		return false;
	}
	else
	{
		return true;
	}
}


// Initialise
void RelionJob::initialise(int _job_type)
{
	type = _job_type;

	bool has_mpi, has_thread;
	switch (type)
	{
	case PROC_IMPORT:
	{
		has_mpi = has_thread = false;
		initialiseImportJob();
		break;
	}
	case PROC_MOTIONCORR:
	{
		has_mpi = has_thread = true;
		initialiseMotioncorrJob();
		break;
	}
	case PROC_CTFFIND:
	{
		has_mpi = true;
		has_thread = false;
		initialiseCtffindJob();
		break;
	}
	case PROC_MANUALPICK:
	{
		has_mpi = has_thread = false;
		initialiseManualpickJob();
		break;
	}
	case PROC_AUTOPICK:
	{
		has_mpi = true;
		has_thread = false;
		initialiseAutopickJob();
		break;
	}
	case PROC_EXTRACT:
	{
		has_mpi = true;
		has_thread = false;
		initialiseExtractJob();
		break;
	}
	case PROC_SORT:
	{
		has_mpi = true;
		has_thread = false;
		initialiseSortJob();
		break;
	}
	case PROC_CLASSSELECT:
	{
		has_mpi = has_thread = false;
		initialiseSelectJob();
		break;
	}
	case PROC_2DCLASS:
	{
		has_mpi = has_thread = true;
		initialiseClass2DJob();
		break;
	}
	case PROC_INIMODEL:
	{
		has_mpi = has_thread = true;
		initialiseInimodelJob();
		break;
	}
	case PROC_3DCLASS:
	{
		has_mpi = has_thread = true;
		initialiseClass3DJob();
		break;
	}
	case PROC_3DAUTO:
	{
		has_mpi = has_thread = true;
		initialiseAutorefineJob();
		break;
	}
	case PROC_MOVIEREFINE:
	{
		has_mpi = has_thread = true;
		initialiseMovierefineJob();
		break;
	}
	case PROC_POLISH:
	{
		has_mpi = has_thread = true;
		initialisePolishJob();
		break;
	}
	case PROC_MASKCREATE:
	{
		has_mpi = has_thread = false;
		initialiseMaskcreateJob();
		break;
	}
	case PROC_JOINSTAR:
	{
		has_mpi = has_thread = false;
		initialiseJoinstarJob();
		break;
	}
	case PROC_SUBTRACT:
	{
		has_mpi = has_thread = false;
		initialiseSubtractJob();
		break;
	}
	case PROC_POST:
	{
		has_mpi = has_thread = false;
		initialisePostprocessJob();
		break;
	}
	case PROC_RESMAP:
	{
		has_mpi = has_thread = true;
		initialiseLocalresJob();
		break;
	}
	default:
	{
		REPORT_ERROR("ERROR: unrecognised job-type");
	}
	}


	if (has_mpi)
		joboptions["nr_mpi"] = JobOption("Number of MPI procs:", 1, 1, 64, 1, "Number of MPI nodes to use in parallel. When set to 1, MPI will not be used.");

	if (has_thread)
		joboptions["nr_threads"] = JobOption("Number of threads:", 1, 1, 16, 1, "Number of shared-memory (POSIX) threads to use in parallel. \
		When set to 1, no multi-threading will be used. Multi-threading is often useful in 3D refinements to have more memory. 2D class averaging often proceeds more efficiently without threads.");

    joboptions["do_queue"] = JobOption("Submit to queue?", false, "If set to Yes, the job will be submit to a queue, otherwise \
the job will be executed locally. Note that only MPI jobs may be sent to a queue.");

    // Need the std::string(), as otherwise it will be overloaded and passed as a boolean....
    joboptions["queuename"] = JobOption("Queue name: ", std::string("openmpi"), "Name of the queue to which to submit the job.");

    joboptions["qsub"] = JobOption("Queue submit command:", std::string("qsub"), "Name of the command used to submit scripts to the queue, e.g. qsub or bsub.\n\n\
Note that the person who installed RELION should have made a custom script for your cluster/queue setup. Check this is the case \
(or create your own script following the RELION WIKI) if you have trouble submitting jobs.");

	// Two additional options that may be set through environment variables RELION_QSUB_EXTRA1 and RELION_QSUB_EXTRA2 (for more flexibility)
	char * extra1_text = getenv ("RELION_QSUB_EXTRA1");
	if (extra1_text != NULL)
	{
		char * extra1_default = getenv ("RELION_QSUB_EXTRA1_DEFAULT");
		char emptychar[] = "";
		if (extra1_default == NULL)
			extra1_default=emptychar;
		joboptions["qsub_extra1"] = JobOption(std::string(extra1_text), std::string(extra1_default), "Extra option to pass to the qsub template script. \
Any occurrences of XXXextra1XXX will be changed by this value.");
	}

	char * extra2_text = getenv ("RELION_QSUB_EXTRA2");
	if (extra2_text != NULL)
	{
		char * extra2_default = getenv ("RELION_QSUB_EXTRA2_DEFAULT");
		char emptychar[] = "";
		if (extra2_default == NULL)
			extra2_default=emptychar;
		joboptions["qsub_extra2"] = JobOption(std::string(extra2_text), std::string(extra2_default), "Extra option to pass to the qsub template script. \
Any occurrences of XXXextra2XXX will be changed by this value.");
	}

	// Check for environment variable RELION_QSUB_TEMPLATE
	char * default_location = getenv ("RELION_QSUB_TEMPLATE");
	if (default_location==NULL)
	{
		char mydefault[]=DEFAULTQSUBLOCATION;
		default_location=mydefault;
	}

    joboptions["qsubscript"] = JobOption("Standard submission script:", std::string(default_location), "Script Files (*.{csh,sh,bash,script})", ".",
"The template for your standard queue job submission script. \
Its default location may be changed by setting the environment variable RELION_QSUB_TEMPLATE. \
In the template script a number of variables will be replaced: \n \
XXXcommandXXX = relion command + arguments; \n \
XXXqueueXXX = The queue name; \n \
XXXmpinodesXXX = The number of MPI nodes; \n \
XXXthreadsXXX = The number of threads; \n \
XXXcoresXXX = XXXmpinodesXXX * XXXthreadsXXX; \n \
XXXdedicatedXXX = The minimum number of dedicated cores on each node; \n \
XXXnodesXXX = The number of requested nodes = CEIL(XXXcoresXXX / XXXdedicatedXXX); \n \
If these options are not enough for your standard jobs, you may define two extra variables: XXXextra1XXX and XXXextra2XXX \
Their help text is set by the environment variables RELION_QSUB_EXTRA1 and RELION_QSUB_EXTRA2 \
For example, setenv RELION_QSUB_EXTRA1 \"Max number of hours in queue\" will result in an additional (text) ein the GUI \
Any variables XXXextra1XXX in the template script will be replaced by the corresponding value.\
Likewise, default values for the extra entries can be set through environment variables RELION_QSUB_EXTRA1_DEFAULT and  RELION_QSUB_EXTRA2_DEFAULT. \
But note that (unlike all other entries in the GUI) the extra values are not remembered from one run to the other.");

	// Check for environment variable RELION_QSUB_TEMPLATE
	char * my_minimum_dedicated = getenv ("RELION_MINIMUM_DEDICATED");
	int minimum_nr_dedicated = (my_minimum_dedicated == NULL) ? DEFAULTMININIMUMDEDICATED : textToInteger(my_minimum_dedicated);
	joboptions["min_dedicated"] = JobOption("Minimum dedicated cores per node:", minimum_nr_dedicated, 1, 64, 1, "Minimum number of dedicated cores that need to be requested on each node. This is useful to force the queue to fill up entire nodes of a given size.");

	// Need the std::string(), as otherwise it will be overloaded and passed as a boolean....
	joboptions["other_args"] = JobOption("Additional arguments:", std::string(""), "In this box command-line arguments may be provided that are not generated by the GUI. \
This may be useful for testing developmental options and/or expert use of the program. \
The command 'relion_refine' will print a list of possible options.");

}


bool RelionJob::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	bool result = false;

	switch (type)
	{
	case PROC_IMPORT:
	{
		result = getCommandsImportJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_MOTIONCORR:
	{
		result = getCommandsMotioncorrJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_CTFFIND:
	{
		result = getCommandsCtffindJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_MANUALPICK:
	{
		result = getCommandsManualpickJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_AUTOPICK:
	{
		result = getCommandsAutopickJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_EXTRACT:
	{
		result = getCommandsExtractJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_SORT:
	{
		result = getCommandsSortJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_CLASSSELECT:
	{
		result = getCommandsSelectJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_2DCLASS:
	{
		result = getCommandsClass2DJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_INIMODEL:
	{
		result = getCommandsInimodelJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_3DCLASS:
	{
		result = getCommandsClass3DJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_3DAUTO:
	{
		result = getCommandsAutorefineJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_MOVIEREFINE:
	{
		result = getCommandsMovierefineJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_POLISH:
	{
		result = getCommandsPolishJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_MASKCREATE:
	{
		result = getCommandsMaskcreateJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_JOINSTAR:
	{
		result = getCommandsJoinstarJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_SUBTRACT:
	{
		result = getCommandsSubtractJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_POST:
	{
		result = getCommandsPostprocessJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	case PROC_RESMAP:
	{
		result = getCommandsLocalresJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
		break;
	}
	default:
	{
		REPORT_ERROR("ERROR: unrecognised job-type");
	}
	}

	return result;

}


void RelionJob::initialiseImportJob()
{

	hidden_name = ".gui_import";

	joboptions["fn_in"] = JobOption("Input files:", "Micrographs/*.mrcs", "Input file (*.*)", ".", "Select any file(s), possibly using Linux wildcards for 2D micrographs or 3D tomograms that you want to import into a STAR file, which will also be saved as a data Node. \n \n \
Note that for importing coordinate files, one has to give a Linux wildcard, where the *-symbol is before the coordinate-file suffix, e.g. if the micrographs are called mic1.mrc and the coordinate files mic1.box or mic1_autopick.star, one HAS to give '*.box' or '*_autopick.star', respectively.\n \n \
Also note that micrographs, movies and coordinate files all need to be in the same directory (with the same rootnames, e.g.mic1 in the example above) in order to be imported correctly. 3D masks or references can be imported from anywhere. \n \n \
Note that movie-particle STAR files cannot be imported from a previous version of RELION, as the way movies are handled has changed in RELION-2.0. \n \n \
For the import of a particle, 2D references or micrograph STAR file or of a 3D reference or mask, only a single file can be imported at a time. \n \n \
Note that due to a bug in a fltk library, you cannot import from directories that contain a substring  of the current directory, e.g. dont important from /home/betagal if your current directory is called /home/betagal_r2. In this case, just change one of the directory names.");

	joboptions["node_type"] = JobOption("Node type:", RADIO_NODETYPE, 0, "Select the type of Node this is.");

}

// Generate the correct commands
bool RelionJob::getCommandsImportJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "Import", job_counter);

	commands.push_back("echo importing...");

	std::string command;
	FileName outputstar;

	std::string fn_in = joboptions["fn_in"].getString();
	std::string node_type = joboptions["node_type"].getString();
	if (node_type == "2D micrograph movies (*.mrcs)")
	{
		outputstar = outputname+"movies.star";
		command = "relion_star_loopheader rlnMicrographMovieName > " + outputstar;;
		commands.push_back(command);
		command = "ls " + fn_in + " >> " + outputstar;
		commands.push_back(command);
		Node node(outputstar, NODE_MOVIES);
		outputNodes.push_back(node);
	}
	else if (node_type == "2D micrographs/tomograms (*.mrc)")
	{
		outputstar = outputname+"micrographs.star";
		command = "relion_star_loopheader rlnMicrographName > " + outputstar;;
		commands.push_back(command);
		command = "ls " + fn_in + " >> " + outputstar;
		commands.push_back(command);
		Node node(outputstar, NODE_MICS);
		outputNodes.push_back(node);
	}
	else if (node_type == "2D/3D particle coordinates (*.box, *_pick.star)")
	{

		// Make the same directory structure of the coordinates
		// Copy all coordinate files into the same subdirectory in the Import directory
		// But remove directory structure from pipeline if that exists
		// Dereference symbolic links if needed
        FileName fn_dir = fn_in;
		if (fn_dir.contains("/"))
			fn_dir = fn_dir.beforeLastOf("/");
		else
			fn_dir = ".";
		FileName fn_pre, fn_jobnr, fn_post;
		if (decomposePipelineSymlinkName(fn_dir, fn_pre, fn_jobnr, fn_post))
		{
			// Make the output directory
			command = "mkdir -p " + outputname + fn_post;
			commands.push_back(command);
			// Copy the coordinates there
			command = "cp " + fn_in + " " + outputname + fn_post;
			commands.push_back(command);
		}
		else
		{
			command = "cp --parents " + fn_in + " " + outputname;
			commands.push_back(command);
		}

		// Make a suffix file, which contains the actual suffix as a suffix
		// Get the coordinate-file suffix
		FileName fn_suffix = fn_in;
		FileName fn_suffix2 = fn_suffix.beforeLastOf("*");
		fn_suffix = fn_suffix.afterLastOf("*");
		fn_suffix = "coords_suffix" + fn_suffix;
		Node node(outputname + fn_suffix, NODE_MIC_COORDS);
		outputNodes.push_back(node);
		command = " echo \\\"" + fn_suffix2 + "*.mrc\\\" > " + outputname + fn_suffix;
		commands.push_back(command);

	}
	else if (node_type == "Particles STAR file (.star)" ||
			 node_type == "Movie-particles STAR file (.star)" ||
			 node_type == "Micrographs STAR file (.star)" ||
			 node_type == "2D references (.star or .mrcs)" ||
			 node_type == "3D reference (.mrc)" ||
			 node_type == "3D mask (.mrc)" ||
			 node_type == "Unfiltered half-map (unfil.mrc)")
	{
		FileName fnt = "/" + fn_in;
		fnt = fnt.afterLastOf("/");
		command = "cp " + fn_in + " " + outputname + fnt;
		commands.push_back(command);

		int mynodetype;
		if (node_type == "Particles STAR file (.star)")
			mynodetype = NODE_PART_DATA;
		else if (node_type == "Movie-particles STAR file (.star)")
			mynodetype = NODE_MOVIE_DATA;
		else if (node_type == "Micrographs STAR file (.star)")
			mynodetype = NODE_MICS;
		else if (node_type == "2D references (.star or .mrcs)")
			mynodetype = NODE_2DREFS;
		else if (node_type == "3D reference (.mrc)")
			mynodetype = NODE_3DREF;
		else if (node_type == "3D mask (.mrc)")
			mynodetype = NODE_MASK;
		else if (node_type == "Unfiltered half-map (unfil.mrc)")
			mynodetype = NODE_HALFMAP;

		Node node(outputname + fnt, mynodetype);
		outputNodes.push_back(node);

		// Also get the other half-map
		if (mynodetype == NODE_HALFMAP)
		{
			FileName fn_inb = fn_in;
			size_t pos = fn_inb.find("half1");
			if (pos != std::string::npos)
			{
				fn_inb.replace(pos, 5, "half2");

			}
			else
			{
				pos = fn_inb.find("half2");
				if (pos != std::string::npos)
				{
					fn_inb.replace(pos, 5, "half1");
				}
			}
			fnt = "/" + fn_inb;
			fnt = fnt.afterLastOf("/");
			command = "cp " + fn_inb + " " + outputname + fnt;
			commands.push_back(command);

			Node node2(outputname + fnt, mynodetype);
			outputNodes.push_back(node2);
		}

	}
	else
	{
		REPORT_ERROR("ImportJobWindow::getCommands ERROR: Unrecognized menu option for node_type= " + node_type);
	}

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);



}

void RelionJob::initialiseMotioncorrJob()
{

	hidden_name = ".gui_motioncorr";

	joboptions["input_star_mics"] = JobOption("Input movies STAR file:", NODE_MOVIES, "", "STAR files (*.star)", "A STAR file with all micrographs to run MOTIONCORR on");
	joboptions["do_save_movies"] = JobOption("Save aligned movie stacks?", true,"Save the aligned movie stacks? Say Yes if you want to perform movie-processing in RELION as well. Say No if you only want to correct motions and write out the averages.");
	joboptions["first_frame_sum"] = JobOption("First frame for corrected sum:", 1, 1, 32, 1, "First frame to use in corrected average (starts counting at 1). ");
	joboptions["last_frame_sum"] = JobOption("Last frame for corrected sum:", 0, 0, 32, 1, "Last frame to use in corrected average.");
	joboptions["angpix"] = JobOption("Pixel size (A):", 1, 0.5, 4.0, 0.1, "Provide the pixel size in Angstroms of the input movies. UNBLUR and MOTIONCOR2 use this for their bfactor and their dose-weighting.");

	// Motioncor2

	// Check for environment variable RELION_MOTIONCOR2_EXECUTABLE
	char * default_location = getenv ("RELION_MOTIONCOR2_EXECUTABLE");
	if (default_location == NULL)
	{
		char mydefault[]=DEFAULTMOTIONCOR2LOCATION;
		default_location=mydefault;
	}

	joboptions["do_motioncor2"] = JobOption("Use MOTIONCOR2?", true ,"If set to Yes, Shawn Zheng's MOTIONCOR2 will be used instead of UNBLUR.");
	joboptions["fn_motioncor2_exe"] = JobOption("MOTIONCOR2 executable:", std::string(default_location), "*.*", ".", "Location of the MOTIONCOR2 executable. You can control the default of this field by setting environment variable RELION_MOTIONCOR2_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code.");
	joboptions["fn_gain_ref"] = JobOption("Gain-reference image:", "", "*.mrc", ".", "Location of the gain-reference file to be applied to the input micrographs. Leave this empty if the movies are already gain-corrected.");
	joboptions["fn_defect"] = JobOption("Defect file:", "", "*", ".", "Location of the MOTIONCOR2-style ASCII file that describes the defect pixels on the detector (using the -DefectFile option). Leave empty if you don't have any defects, or don't want to correct for defects on your detector.");
	joboptions["patch_x"] = JobOption("Number of patches X:", std::string("1"), "Number of patches (in X and Y direction) to apply motioncor2.");
	joboptions["patch_y"] = JobOption("Number of patches Y:", std::string("1"), "Number of patches (in X and Y direction) to apply motioncor2.");
	joboptions["group_frames"] = JobOption("Group frames:", 1, 1, 5, 1, "Average together this many frames before calculating the beam-induced shifts.");
	joboptions["bin_factor"] = JobOption("Binning factor:", 1, 1, 2, 1, "Bin the micrographs this much by a windowing operation in the Fourier Tranform. Binning at this level is hard to un-do later on, but may be useful to down-scale super-resolution images. Float-values may be used for MOTIONCOR2. Do make sure though that the resulting micrograph size is even.");
	joboptions["bfactor"] = JobOption("Bfactor:", 150, 0, 1500, 50, "The B-factor (-bft) that MOTIONCOR2 will apply to the micrographs.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use: ", std::string("0"), "Provide a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':'. For example, to place two ranks on device 0 and one rank on device 1, provide '0:0:1'");
	joboptions["other_motioncor2_args"] = JobOption("Other MOTIONCOR2 arguments", std::string(""), "Additional arguments that need to be passed to MOTIONCOR2.");

	// Unblur

	joboptions["do_unblur"] = JobOption("Use UNBLUR instead?", false ,"If set to Yes, Niko Grigoerieff's UNBLUR will be used instead of MOTIONCORR. Only default settings in UNBLUR are allowed in this wrapper. Note that all options from the MOTIONCORR tab will be ignored.");

	// Check for environment variable RELION_UNBLUR_EXECUTABL
	char * default_location2 = getenv ("RELION_UNBLUR_EXECUTABLE");
	if (default_location2 == NULL)
	{
		char mydefault2[]=DEFAULTUNBLURLOCATION;
		default_location2=mydefault2;
	}
	joboptions["fn_unblur_exe"] = JobOption("UNBLUR executable:", std::string(default_location2), "*.*", ".", "Location of the UNBLUR executable. You can control the default of this field by setting environment variable RELION_UNBLUR_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code. Note that this wrapper was tested with unblur version 1.0.2.");

	// Check for environment variable RELION_SUMMOVIE_EXECUTABL
	char * default_location3 = getenv ("RELION_SUMMOVIE_EXECUTABLE");
	if (default_location3 == NULL)
	{
		char mydefault3[]=DEFAULTSUMMOVIELOCATION;
		default_location3=mydefault3;
	}
	joboptions["fn_summovie_exe"] = JobOption("SUMMOVIE executable:", std::string(default_location3), "*.*", ".", "Location of the SUMMOVIE executable. You can control the default of this field by setting environment variable RELION_SUMMOVIE_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code. Note that this wrapper was tested with summovie version 1.0.2.");

	// Dose-weight

	joboptions["do_dose_weighting"] = JobOption("Do dose-weighting?", false ,"If set to Yes, the averaged micrographs will be dose-weighted using either UNBLUR or MOTIONCOR2. Note this does not work in MOTIONCORR.");
	joboptions["voltage"] = JobOption("Voltage (kV):", 300, 80, 300, 20, "Acceleration voltage in kV.");
	joboptions["dose_per_frame"] = JobOption("Dose per frame (e/A2):", 1, 0, 5, 0.2, "Dose per movie frame (in electrons per squared Angstrom).");
	joboptions["pre_exposure"] = JobOption("Pre-exposure (e/A2):", 0, 0, 5, 0.5, "Pre-exposure dose (in electrons per squared Angstrom). This is only used by UNBLUR, not by MOTIONCOR2.");

}

bool RelionJob::getCommandsMotioncorrJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, "MotionCorr", job_counter);

	std::string command;
	if (joboptions["nr_mpi"].getNumber() > 1)
		command="`which relion_run_motioncorr_mpi`";
	else
		command="`which relion_run_motioncorr`";

	// I/O

	if (joboptions["input_star_mics"].getString() == "")
	{
		error_message = "ERROR: empty field for input STAR file...";
		return false;
	}

	command += " --i " + joboptions["input_star_mics"].getString();
	Node node(joboptions["input_star_mics"].getString(), joboptions["input_star_mics"].node_type);
	inputNodes.push_back(node);

	command += " --o " + outputname;
	outputName = outputname;
	Node node2(outputname + "corrected_micrographs.star", NODE_MICS);
	outputNodes.push_back(node2);
	if (joboptions["do_save_movies"].getBoolean())
	{
		Node node3(outputname + "corrected_micrograph_movies.star", NODE_MOVIES);
		outputNodes.push_back(node3);
	}
	Node node4(outputname + "logfile.pdf", NODE_PDF_LOGFILE);
	outputNodes.push_back(node4);

	if (joboptions["do_save_movies"].getBoolean())
		command += " --save_movies ";

	command += " --first_frame_sum " + joboptions["first_frame_sum"].getString();
	command += " --last_frame_sum " + joboptions["last_frame_sum"].getString();

	if (joboptions["do_unblur"].getBoolean())
	{
		// Unblur-specific stuff
		command += " --use_unblur";
		command += " --j " + joboptions["nr_threads"].getString();
		command += " --unblur_exe " + joboptions["fn_unblur_exe"].getString();
		command += " --summovie_exe " + joboptions["fn_summovie_exe"].getString();
		command += " --angpix " +  joboptions["angpix"].getString();
	}
	else if (joboptions["do_motioncor2"].getBoolean())
	{
		// Motioncor2-specific stuff
		command += " --use_motioncor2";
		command += " --bin_factor " + joboptions["bin_factor"].getString();
		command += " --motioncor2_exe " + joboptions["fn_motioncor2_exe"].getString();
		command += " --bfactor " + joboptions["bfactor"].getString();
		command += " --angpix " +  joboptions["angpix"].getString();
		command += " --patch_x " + joboptions["patch_x"].getString();
		command += " --patch_y " + joboptions["patch_y"].getString();
		if (joboptions["group_frames"].getNumber() > 1.)
			command += " --group_frames " + joboptions["group_frames"].getString();
		if ((joboptions["fn_gain_ref"].getString()).length() > 0)
			command += " --gainref " + joboptions["fn_gain_ref"].getString();
		if ((joboptions["fn_defect"].getString()).length() > 0)
			command += " --defect_file " + joboptions["fn_defect"].getString();


		if ((joboptions["other_motioncor2_args"].getString()).length() > 0)
			command += " --other_motioncor2_args \" " + joboptions["other_motioncor2_args"].getString() + " \"";

		// Which GPUs to use?
		command += " --gpu \"" + joboptions["gpu_ids"].getString() + "\"";

	}
	else
	{
		error_message = "ERROR: specify the use of MotionCor2 or Unblur...";
		return false;
	}

	if (joboptions["do_dose_weighting"].getBoolean())
	{
		command += " --dose_weighting";
		command += " --voltage " + joboptions["voltage"].getString();
		command += " --dose_per_frame " + joboptions["dose_per_frame"].getString();
		command += " --preexposure " + joboptions["pre_exposure"].getString();
	}

	if (is_continue)
		command += " --only_do_unfinished ";


	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);


}

void RelionJob::initialiseCtffindJob()
{

	hidden_name = ".gui_ctffind";

	char *default_location;

	joboptions["input_star_mics"] = JobOption("Input micrographs STAR file:", NODE_MICS, "", "STAR files (*.star)", "A STAR file with all micrographs to run CTFFIND or Gctf on");
	joboptions["use_noDW"] = JobOption("Use micrograph without dose-weighting?", false, "If set to Yes, the CTF estimation will be done using the micrograph without dose-weighting as in rlnMicrographNameNoDW (_noDW.mrc from MotionCor2). If set to No, the normal rlnMicrographName will be used. Note this option will not work for motion-correction with UNBLUR.");

	joboptions["cs"] = JobOption("Spherical aberration (mm):", 2, 0, 8, 0.1, "Spherical aberration of the microscope used to collect these images (in mm)");
	joboptions["kv"] = JobOption("Voltage (kV):", 300, 50, 500, 10, "Voltage the microscope was operated on (in kV)");
	joboptions["q0"] = JobOption("Amplitude contrast:", 0.1, 0, 0.3, 0.01, "Fraction of amplitude contrast. Often values around 10% work better than theoretically more accurate lower values...");
	joboptions["angpix"] = JobOption("Magnified pixel size (Angstrom):", 1.4, 0.5, 3, 0.1, "Pixel size in Angstroms. ");
	joboptions["dast"] = JobOption("Amount of astigmatism (A):", 100, 0, 2000, 100,"CTFFIND's dAst parameter, GCTFs -astm parameter");

	// CTFFIND options

	joboptions["box"] = JobOption("FFT box size (pix):", 512, 64, 1024, 8, "CTFFIND's Box parameter");
	joboptions["resmin"] = JobOption("Minimum resolution (A):", 30, 10, 200, 10, "CTFFIND's ResMin parameter");
	joboptions["resmax"] = JobOption("Maximum resolution (A):", 5, 1, 20, 1, "CTFFIND's ResMax parameter");
	joboptions["dfmin"] = JobOption("Minimum defocus value (A):", 5000, 0, 25000, 1000, "CTFFIND's dFMin parameter");
	joboptions["dfmax"] = JobOption("Maximum defocus value (A):", 50000, 20000, 100000, 1000, "CTFFIND's dFMax parameter");
	joboptions["dfstep"] = JobOption("Defocus step size (A):", 500, 200, 2000, 100,"CTFFIND's FStep parameter");

	joboptions["do_phaseshift"] = JobOption("Estimate phase shifts?", false, "If set to Yes, CTFFIND4 will estimate the phase shift, e.g. as introduced by a Volta phase-plate");
	joboptions["phase_min"] = JobOption("Phase shift (deg) - Min:", std::string("0"), "Minimum, maximum and step size (in degrees) for the search of the phase shift");
	joboptions["phase_max"] = JobOption("Phase shift (deg) - Max:", std::string("180"), "Minimum, maximum and step size (in degrees) for the search of the phase shift");
	joboptions["phase_step"] = JobOption("Phase shift (deg) - Step:", std::string("10"), "Minimum, maximum and step size (in degrees) for the search of the phase shift");

	// Check for environment variable RELION_CTFFIND_EXECUTABLE
	joboptions["use_ctffind4"] = JobOption("Use CTFFIND-4.1?", false, "If set to Yes, the wrapper will use CTFFIND4 (version 4.1) for CTF estimation. This includes thread-support, calculation of Thon rings from movie frames and phase-shift estimation for phase-plate data.");

	default_location = getenv ("RELION_CTFFIND_EXECUTABLE");
	if (default_location == NULL)
	{
		char mydefault[]=DEFAULTCTFFINDLOCATION;
		default_location=mydefault;
	}
	joboptions["fn_ctffind_exe"] = JobOption("CTFFIND-4.1 executable:", std::string(default_location), "*.exe", ".", "Location of the CTFFIND (release 4.1 or later) executable. You can control the default of this field by setting environment variable RELION_CTFFIND_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code.");
	joboptions["do_movie_thon_rings"] = JobOption("Estimate Thon rings from movies?", false, "If set to Yes, CTFFIND4 will calculate power spectra of averages of several movie frames and then average those power spectra to calculate Thon rings. This may give better rings than calculation the power spectra of averages of all movie frames, although it does come at increased costs of processing and disk access");
	joboptions["movie_rootname"] = JobOption("Movie rootname plus extension", std::string("_movie.mrcs"), "Give the movie rootname and extension for all movie files. Movies are assumed to be next to the average micrographs in the same directory.");
	joboptions["avg_movie_frames"] = JobOption("Nr of movie frames to average:", 4, 1, 20, 1,"Calculate averages over so many movie frames to calculate power spectra. Often values corresponding to an accumulated dose of ~ 4 electrons per squared Angstrom work well.");
	joboptions["ctf_win"] = JobOption("Estimate CTF on window size (pix) ", -1, -16, 4096, 16, "If a positive value is given, a squared window of this size at the center of the micrograph will be used to estimate the CTF. This may be useful to exclude parts of the micrograph that are unsuitable for CTF estimation, e.g. the labels at the edge of phtographic film. \n \n The original micrograph will be used (i.e. this option will be ignored) if a negative value is given.");


	joboptions["use_gctf"] = JobOption("Use Gctf instead?", false, "If set to Yes, Kai Zhang's Gctf program (which runs on NVIDIA GPUs) will be used instead of Niko Grigorieff's CTFFIND4.");
	// Check for environment variable RELION_CTFFIND_EXECUTABLE
	default_location = getenv ("RELION_GCTF_EXECUTABLE");
	if (default_location == NULL)
	{
		char mydefault[]=DEFAULTGCTFLOCATION;
		default_location=mydefault;
	}
	joboptions["fn_gctf_exe"] = JobOption("Gctf executable:", std::string(default_location), "*", ".", "Location of the Gctf executable. You can control the default of this field by setting environment variable RELION_GCTF_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code.");
	joboptions["do_ignore_ctffind_params"] = JobOption("Ignore 'Searches' parameters?", true, "If set to Yes, all parameters on the 'Searches' tab will be ignored, and Gctf's default parameters will be used (box.size=1024; min.resol=50; max.resol=4; min.defocus=500; max.defocus=90000; step.defocus=500; astigm=1000) \n \
\n If set to No, all parameters on the CTFFIND tab will be passed to Gctf.");
	joboptions["do_EPA"] = JobOption("Perform equi-phase averaging?", false, "If set to Yes, equi-phase averaging is used in the defocus refinement, otherwise basic rotational averaging will be performed.");
	joboptions["other_gctf_args"] = JobOption("Other Gctf options:", std::string(""), "Provide additional gctf options here.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use: ", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','. ");

}
bool RelionJob::getCommandsCtffindJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "CtfFind", job_counter);
	std::string command;

	FileName fn_outstar = outputname + "micrographs_ctf.star";
	Node node(fn_outstar, NODE_MICS);
	outputNodes.push_back(node);
	outputName = outputname;

	if (joboptions["input_star_mics"].getString() == "")
	{
		error_message = "ERROR: empty field for input STAR file...";
		return false;
	}

	Node node2(joboptions["input_star_mics"].getString(), joboptions["input_star_mics"].node_type);
	inputNodes.push_back(node2);

	if (joboptions["nr_mpi"].getNumber() > 1)
		command="`which relion_run_ctffind_mpi`";
	else
		command="`which relion_run_ctffind`";

	// Just use 10000 for magnification, so dstep==angpix
	RFLOAT magn = 10000.;

	command += " --i " + joboptions["input_star_mics"].getString();
	command += " --o " + outputname;
	command += " --CS " + joboptions["cs"].getString();
	command += " --HT " + joboptions["kv"].getString();
	command += " --AmpCnst " + joboptions["q0"].getString();
	command += " --XMAG " + floatToString(magn);
	command += " --DStep " + joboptions["angpix"].getString();
	command += " --Box " + joboptions["box"].getString();
	command += " --ResMin " + joboptions["resmin"].getString();
	command += " --ResMax " + joboptions["resmax"].getString();
	command += " --dFMin " + joboptions["dfmin"].getString();
	command += " --dFMax " + joboptions["dfmax"].getString();
	command += " --FStep " + joboptions["dfstep"].getString();
	command += " --dAst " + joboptions["dast"].getString();

	if (joboptions["use_noDW"].getBoolean())
		command += " --use_noDW ";

	if (joboptions["do_phaseshift"].getBoolean())
	{
		command += " --do_phaseshift ";
		command += " --phase_min " + joboptions["phase_min"].getString();
		command += " --phase_max " + joboptions["phase_max"].getString();
		command += " --phase_step " + joboptions["phase_step"].getString();
	}

	if (joboptions["use_gctf"].getBoolean())
	{
		command += " --use_gctf --gctf_exe " + joboptions["fn_gctf_exe"].getString();
		command += " --angpix " + joboptions["angpix"].getString();
		if (joboptions["do_ignore_ctffind_params"].getBoolean())
			command += " --ignore_ctffind_params";
		if (joboptions["do_EPA"].getBoolean())
			command += " --EPA";

		// GPU-allocation
		command += " --gpu \"" + joboptions["gpu_ids"].getString() + "\"";

		if ((joboptions["other_gctf_args"].getString()).length() > 0)
			command += " --extra_gctf_options \" " + joboptions["other_gctf_args"].getString() + " \"";

	}
	else if (joboptions["use_ctffind4"].getBoolean())
	{
		command += " --ctffind_exe " + joboptions["fn_ctffind_exe"].getString();
		command += " --ctfWin " + joboptions["ctf_win"].getString();
		command += " --is_ctffind4 ";
		if (joboptions["do_movie_thon_rings"].getBoolean())
		{
			command += " --do_movie_thon_rings --avg_movie_frames " + joboptions["avg_movie_frames"].getString();
			command += " --movie_rootname " + joboptions["movie_rootname"].getString();
		}
	}
	else
	{
		error_message = "ERROR: Please select use of CTFFIND4.1 or Gctf...";
		return false;
	}

	if (is_continue)
		command += " --only_do_unfinished ";

	// Other arguments
	command += " " + joboptions["other_args"].getString();
	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

void RelionJob::initialiseManualpickJob()
{

	hidden_name = ".gui_manualpick";

	joboptions["fn_in"] = JobOption("Input micrographs:", NODE_MICS, "", "Input micrographs (*.{star,mrc})", "Input STAR file (with or without CTF information), OR a unix-type wildcard with all micrographs in MRC format (in this case no CTFs can be used).");

	joboptions["diameter"] = JobOption("Particle diameter (A):", 100, 0, 500, 50, "The radius of the circle used around picked particles (in Angstroms). Only used for display." );
	joboptions["micscale"] = JobOption("Scale for micrographs:", 0.2, 0.1, 1, 0.05, "The micrographs will be displayed at this relative scale, i.e. a value of 0.5 means that only every second pixel will be displayed." );
	joboptions["sigma_contrast"] = JobOption("Sigma contrast:", 3, 0, 10, 0.5, "The micrographs will be displayed with the black value set to the average of all values MINUS this values times the standard deviation of all values in the micrograph, and the white value will be set \
to the average PLUS this value times the standard deviation. Use zero to set the minimum value in the micrograph to black, and the maximum value to white ");
	joboptions["white_val"] = JobOption("White value:", 0, 0, 512, 16, "Use non-zero values to set the value of the whitest pixel in the micrograph.");
	joboptions["black_val"] = JobOption("Black value:", 0, 0, 512, 16, "Use non-zero values to set the value of the blackest pixel in the micrograph.");

	joboptions["lowpass"] = JobOption("Lowpass filter (A)", 20, 10, 100, 5, "Lowpass filter that will be applied to the micrographs. Give a negative value to skip the lowpass filter.");
	joboptions["highpass"] = JobOption("Highpass filter (A)", -1, 100, 1000, 100, "Highpass filter that will be applied to the micrographs. This may be useful to get rid of background ramps due to uneven ice distributions. Give a negative value to skip the highpass filter. Useful values are often in the range of 200-400 Angstroms.");
	joboptions["angpix"] = JobOption("Pixel size (A)", -1, 0.3, 5, 0.1, "Pixel size in Angstroms. This will be used to calculate the filters and the particle diameter in pixels. If a CTF-containing STAR file is input, then the value given here will be ignored, and the pixel size will be calculated from the values in the STAR file. A negative value can then be given here.");

	joboptions["do_startend"] = JobOption("Pick start-end coordinates helices?", false, "If set to true, start and end coordinates are picked subsequently and a line will be drawn between each pair");

	joboptions["ctfscale"] = JobOption("Scale for CTF image:", 1, 0.1, 2, 0.1, "CTFFINDs CTF image (with the Thonrings) will be displayed at this relative scale, i.e. a value of 0.5 means that only every second pixel will be displayed." );

	joboptions["do_color"] = JobOption("Blue<>red color particles?", false, "If set to true, then the circles for each particles are coloured from red to blue (or the other way around) for a given metadatalabel. If this metadatalabel is not in the picked coordinates STAR file \
(basically only the rlnAutopickFigureOfMerit or rlnClassNumber) would be useful values there, then you may provide an additional STAR file (e.g. after classification/refinement below. Particles with values -999, or that are not in the additional STAR file will be coloured the default color: green");
	joboptions["color_label"] = JobOption("MetaDataLabel for color:", std::string("rlnParticleSelectZScore"), "The Metadata label of the value to plot from red<>blue. Useful examples might be: \n \
rlnParticleSelectZScore \n rlnClassNumber \n rlnAutopickFigureOfMerit \n rlnAngleTilt \n rlnLogLikeliContribution \n rlnMaxValueProbDistribution \n rlnNrOfSignificantSamples\n");
	joboptions["fn_color"] = JobOption("STAR file with color label: ", "", "STAR file (*.star)", ".", "The program will figure out which particles in this STAR file are on the current micrograph and color their circles according to the value in the corresponding column. \
Particles that are not in this STAR file, but present in the picked coordinates file will be colored green. If this field is left empty, then the color label (e.g. rlnAutopickFigureOfMerit) should be present in the coordinates STAR file.");
	joboptions["blue_value"] = JobOption("Blue value: ", 0., 0., 4., 0.1, "The value of this entry will be blue. There will be a linear scale from blue to red, according to this value and the one given below.");
	joboptions["red_value"] = JobOption("Red value: ", 2., 0., 4., 0.1, "The value of this entry will be red. There will be a linear scale from blue to red, according to this value and the one given above.");


}
bool RelionJob::getCommandsManualpickJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "ManualPick", job_counter);
	std::string command;

	command="`which relion_manualpick`";

	if (joboptions["fn_in"].getString() == "")
	{
		error_message = "ERROR: empty field for input STAR file...";
		return false;
	}

	command += " --i " + joboptions["fn_in"].getString();
	Node node(joboptions["fn_in"].getString(), joboptions["fn_in"].node_type);
	inputNodes.push_back(node);

	command += " --odir " + outputname;
	command += " --pickname manualpick";

	FileName fn_suffix = outputname + "coords_suffix_manualpick.star";
	Node node2(fn_suffix, NODE_MIC_COORDS);
	outputNodes.push_back(node2);

	// Allow saving, and always save default selection file upon launching the program
	FileName fn_outstar = outputname + "micrographs_selected.star";
	Node node3(fn_outstar, NODE_MICS);
	outputNodes.push_back(node3);
	command += " --allow_save   --fast_save --selection " + fn_outstar;

	command += " --scale " + joboptions["micscale"].getString();
	command += " --sigma_contrast " + joboptions["sigma_contrast"].getString();
	command += " --black " + joboptions["black_val"].getString();
	command += " --white " + joboptions["white_val"].getString();

	if (joboptions["lowpass"].getNumber() > 0.)
		command += " --lowpass " + joboptions["lowpass"].getString();
	if (joboptions["highpass"].getNumber() > 0.)
		command += " --highpass " + joboptions["highpass"].getString();
	if (joboptions["angpix"].getNumber() > 0.)
		command += " --angpix " + joboptions["angpix"].getString();

	command += " --ctf_scale " + joboptions["ctfscale"].getString();

	command += " --particle_diameter " + joboptions["diameter"].getString();

	if (joboptions["do_startend"].getBoolean())
	{
		command += " --pick_start_end ";
	}

	if (joboptions["do_color"].getBoolean())
	{
		command += " --color_label " + joboptions["color_label"].getString();
		command += " --blue " + joboptions["blue_value"].getString();
		command += " --red " + joboptions["red_value"].getString();
		if (joboptions["fn_color"].getString().length() > 0)
			command += " --color_star " + joboptions["fn_color"].getString();
	}

	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	// Also make the suffix file (do this after previous command was pushed back!)
	// Inside it, store the name of the micrograph STAR file, so we can display these later
	FileName fn_pre, fn_jobnr, fn_post;
	decomposePipelineSymlinkName(joboptions["fn_in"].getString(), fn_pre, fn_jobnr, fn_post);
	command = "echo " + fn_pre + fn_jobnr + fn_post + " > " + fn_suffix;
	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

void RelionJob::initialiseAutopickJob()
{

	hidden_name = ".gui_autopick";

	joboptions["fn_input_autopick"] = JobOption("Input micrographs for autopick:", NODE_MICS, "", "Input micrographs (*.{star})", "Input STAR file (preferably with CTF information) with all micrographs to pick from.");
	joboptions["fn_refs_autopick"] = JobOption("References:", NODE_2DREFS, "", "Input references (*.{star,mrcs})", "Input STAR file or MRC stack with the 2D references to be used for picking. Note that the absolute greyscale needs to be correct, so only use images created by RELION itself, e.g. by 2D class averaging or projecting a RELION reconstruction.");
	joboptions["do_gauss_ref"] = JobOption("Or use Gaussian blob?", false, "If set to Yes, a Gaussian blob will be used as a reference (you can then leave the 'References' field empty. The preferred way to autopick is by setting this to no and providing references that were generated by 2D classification from this data set in RELION. The Gaussian blob references may be useful to kickstart a new data set.");
	joboptions["gauss_max"] = JobOption("Gaussian peak value", 0.1, 0.01, 0.5, 0.01, "The peak value of the Gaussian blob. Weaker data will need lower values.");
	joboptions["angpix"] = JobOption("Pixel size in micrographs (A)", -1, 0.3, 5, 0.1, "Pixel size in Angstroms. If a CTF-containing STAR file is input, then the value given here will be ignored, and the pixel size will be calculated from the values in the STAR file. A negative value can then be given here.");
	joboptions["particle_diameter"] = JobOption("Mask diameter (A)", -1, 0, 2000, 20, "Diameter of the circular mask that will be applied around the templates in Angstroms. When set to a negative value, this value is estimated automatically from the templates themselves.");

	joboptions["lowpass"] = JobOption("Lowpass filter references (A)", 20, 10, 100, 5, "Lowpass filter that will be applied to the references before template matching. Do NOT use very high-resolution templates to search your micrographs. The signal will be too weak at high resolution anyway, and you may find Einstein from noise.... Give a negative value to skip the lowpass filter.");
	joboptions["highpass"] = JobOption("Highpass filter (A)", -1, 100, 1000, 100, "Highpass filter that will be applied to the micrographs. This may be useful to get rid of background ramps due to uneven ice distributions. Give a negative value to skip the highpass filter.  Useful values are often in the range of 200-400 Angstroms.");
	joboptions["angpix_ref"] = JobOption("Pixel size in references (A)", -1, 0.3, 5, 0.1, "Pixel size in Angstroms for the provided reference images. This will be used to calculate the filters and the particle diameter in pixels. If a negative value is given here, the pixel size in the references will be assumed to be the same as the one in the micrographs, i.e. the particles that were used to make the references were not rescaled upon extraction.");
	joboptions["psi_sampling_autopick"] = JobOption("Angular sampling (deg)", 5, 1, 30, 1, "Angular sampling in degrees for exhaustive searches of the in-plane rotations for all references.");
	joboptions["do_invert_refs"] = JobOption("References have inverted contrast?", true, "Set to Yes to indicate that the reference have inverted contrast with respect to the particles in the micrographs.");
	joboptions["do_ctf_autopick"] = JobOption("Are References CTF corrected?", true, "Set to Yes if the references were created with CTF-correction inside RELION. \n \n If set to Yes, the input micrographs can only be given as a STAR file, which should contain the CTF information for each micrograph.");
	joboptions["do_ignore_first_ctfpeak_autopick"] = JobOption("Ignore CTFs until first peak?", false,"Set this to Yes, only if this option was also used to generate the references.");

	joboptions["threshold_autopick"] = JobOption("Picking threshold:", 0.05, 0, 1., 0.01, "Use lower thresholds to pick more particles (and more junk probably)");
	joboptions["mindist_autopick"] = JobOption("Minimum inter-particle distance (A):", 100, 0, 1000, 20, "Particles closer together than this distance will be consider to be a single cluster. From each cluster, only one particle will be picked. \
\n\nThis option takes no effect for picking helical segments. The inter-box distance is calculated with the number of asymmetrical units and the helical rise on 'Helix' tab.");
	joboptions["maxstddevnoise_autopick"] = JobOption("Maximum stddev noise:", 1.1, 0.9, 1.5, 0.02, "This is useful to prevent picking in carbon areas, or areas with big contamination features. Peaks in areas where the background standard deviation in the normalized micrographs is higher than this value will be ignored. Useful values are probably in the range 1.0 to 1.2. Set to -1 to switch off the feature to eliminate peaks due to high background standard deviations.");
	joboptions["do_write_fom_maps"] = JobOption("Write FOM maps?", false, "If set to Yes, intermediate probability maps will be written out, which (upon reading them back in) will speed up tremendously the optimization of the threshold and inter-particle distance parameters. However, with this option, one cannot run in parallel, as disc I/O is very heavy with this option set.");
	joboptions["do_read_fom_maps"] = JobOption("Read FOM maps?", false, "If written out previously, read the FOM maps back in and re-run the picking to quickly find the optimal threshold and inter-particle distance parameters");

	joboptions["shrink"] = JobOption("Shrink factor:", 1, 0, 1, 0.1, "This is useful to speed up the calculations, and to make them less memory-intensive. The micrographs will be downscaled (shrunk) to calculate the cross-correlations, and peak searching will be done in the downscaled FOM maps. When set to 0, the micrographs will de downscaled to the lowpass filter of the references, a value between 0 and 1 will downscale the micrographs by that factor. Note that the results will not be exactly the same when you shrink micrographs!");
	joboptions["use_gpu"] = JobOption("Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':'. For example: 0:1:0:1:0:1");
	joboptions["do_pick_helical_segments"] = JobOption("Pick 2D helical segments?", false, "Set to Yes if you want to pick 2D helical segments.");

	joboptions["helical_tube_outer_diameter"] = JobOption("Tube diameter (A): ", 200, 100, 1000, 10, "Outer diameter (in Angstroms) of helical tubes. \
This value should be slightly larger than the actual width of the tubes.");
	joboptions["helical_nr_asu"] = JobOption("Number of asymmetrical units:", 1, 1, 100, 1, "Number of helical asymmetrical units in each segment box. This integer should not be less than 1. The inter-box distance (pixels) = helical rise (Angstroms) * number of asymmetrical units / pixel size (Angstroms). \
The optimal inter-box distance might also depend on the box size, the helical rise and the flexibility of the structure. In general, an inter-box distance of ~10% * the box size seems appropriate.");
	joboptions["helical_rise"] = JobOption("Helical rise (A):", -1, 0, 100, 0.01, "Helical rise in Angstroms. (Please click '?' next to the option above for details about how the inter-box distance is calculated.)");
	joboptions["helical_tube_kappa_max"] = JobOption("Maximum curvature (kappa): ", 0.1, 0.05, 0.5, 0.01, "Maximum curvature allowed for picking helical tubes. \
Kappa = 0.3 means that the curvature of the picked helical tubes should not be larger than 30% the curvature of a circle (diameter = particle mask diameter). \
Kappa ~ 0.05 is recommended for long and straight tubes (e.g. TMV, VipA/VipB and AChR tubes) while 0.20 ~ 0.40 seems suitable for flexible ones (e.g. ParM and MAVS-CARD filaments).");
	joboptions["helical_tube_length_min"] = JobOption("Minimum length (A): ", -1, 100, 1000, 10, "Minimum length (in Angstroms) of helical tubes for auto-picking. \
Helical tubes with shorter lengths will not be picked. Note that a long helical tube seen by human eye might be treated as short broken pieces due to low FOM values or high picking threshold.");

}

bool RelionJob::getCommandsAutopickJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "AutoPick", job_counter);
	std::string command;
	if (joboptions["nr_mpi"].getNumber() > 1)
		command="`which relion_autopick_mpi`";
	else
		command="`which relion_autopick`";

	// Input
	if (joboptions["fn_input_autopick"].getString() == "")
	{
		error_message = "ERROR: empty field for input STAR file...";
		return false;
	}
	command += " --i " + joboptions["fn_input_autopick"].getString();
	Node node(joboptions["fn_input_autopick"].getString(), joboptions["fn_input_autopick"].node_type);
	inputNodes.push_back(node);

	if (joboptions["do_gauss_ref"].getBoolean())
	{
		command += " --ref gauss ";
		command += " --gauss_max " + joboptions["gauss_max"].getString();
	}
	else
	{
		if (joboptions["fn_refs_autopick"].getString() == "")
		{
			error_message ="ERROR: empty field for references...";
			return false;
		}
		command += " --ref " + joboptions["fn_refs_autopick"].getString();
		Node node2(joboptions["fn_refs_autopick"].getString(), joboptions["fn_refs_autopick"].node_type);
		inputNodes.push_back(node2);
	}

	// Output
	Node node3(outputname + "coords_suffix_autopick.star", NODE_MIC_COORDS);
	outputNodes.push_back(node3);

	command += " --odir " + outputname;
	command += " --pickname autopick";

	if (joboptions["do_invert_refs"].getBoolean())
		command += " --invert ";

	if (joboptions["do_ctf_autopick"].getBoolean())
	{
		command += " --ctf ";
		if (joboptions["do_ignore_first_ctfpeak_autopick"].getBoolean())
			command += " --ctf_intact_first_peak ";
	}
	command += " --ang " + joboptions["psi_sampling_autopick"].getString();

	command += " --shrink " + joboptions["shrink"].getString();
	if (joboptions["lowpass"].getNumber() > 0.)
		command += " --lowpass " + joboptions["lowpass"].getString();
	if (joboptions["highpass"].getNumber() > 0.)
		command += " --highpass " + joboptions["highpass"].getString();
	if (joboptions["angpix"].getNumber() > 0.)
		command += " --angpix " + joboptions["angpix"].getString();
	if (joboptions["angpix_ref"].getNumber() > 0.)
		command += " --angpix_ref " + joboptions["angpix_ref"].getString();
	if (joboptions["particle_diameter"].getNumber() > 0.)
		command += " --particle_diameter " + joboptions["particle_diameter"].getString();

	if (joboptions["do_write_fom_maps"].getBoolean())
		command += " --write_fom_maps ";

	if (joboptions["do_read_fom_maps"].getBoolean())
		command += " --read_fom_maps ";

	command += " --threshold " + joboptions["threshold_autopick"].getString();
	if (joboptions["do_pick_helical_segments"].getBoolean())
		command += " --min_distance " + floatToString(joboptions["helical_nr_asu"].getNumber() * joboptions["helical_rise"].getNumber());
	else
		command += " --min_distance " + joboptions["mindist_autopick"].getString();
	command += " --max_stddev_noise " + joboptions["maxstddevnoise_autopick"].getString();

	// Helix
	if (joboptions["do_pick_helical_segments"].getBoolean())
	{
		command += " --helix";
		command += " --helical_tube_outer_diameter " + joboptions["helical_tube_outer_diameter"].getString();
		command += " --helical_tube_kappa_max " + joboptions["helical_tube_kappa_max"].getString();
		command += " --helical_tube_length_min " + joboptions["helical_tube_length_min"].getString();
	}

	if (is_continue && !(joboptions["do_read_fom_maps"].getBoolean() || joboptions["do_write_fom_maps"].getBoolean()))
		command += " --only_do_unfinished ";

	// GPU-stuff
	if (joboptions["use_gpu"].getBoolean())
	{
		// for the moment always use --shrink 0 with GPUs ...
		command += " --gpu \"" + joboptions["gpu_ids"].getString() + "\"";
	}

	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	// Also touch the suffix file. Do this after the first command had completed
	// Instead of the symlink from the alias, use the original jobnr filename
	FileName fn_pre, fn_jobnr, fn_post;
	decomposePipelineSymlinkName(joboptions["fn_input_autopick"].getString(), fn_pre, fn_jobnr, fn_post);
	command = "echo " + fn_pre + fn_jobnr + fn_post + " > " +  outputname + "coords_suffix_autopick.star";
	commands.push_back(command.c_str());

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);


}

void RelionJob::initialiseExtractJob()
{

	hidden_name = ".gui_extract";

    joboptions["star_mics"]= JobOption("micrograph STAR file: ", NODE_MICS, "", "Input STAR file (*.{star})", "Filename of the STAR file that contains all micrographs from which to extract particles.");
	joboptions["coords_suffix"] = JobOption("Input coordinates: ", NODE_MIC_COORDS, "", "Input coords_suffix file ({coords_suffix}*)", "Filename of the coords_suffix file with the directory structure and the suffix of all coordinate files.");
	joboptions["do_reextract"] = JobOption("OR re-extract refined particles? ", false, "If set to Yes, the input Coordinates above will be ignored. Instead, one uses a _data.star file from a previous 2D or 3D refinement to re-extract the particles in that refinement, possibly re-centered with their refined origin offsets. This is particularly useful when going from binned to unbinned particles.");
	joboptions["fndata_reextract"] = JobOption("Refined particles STAR file: ", NODE_PART_DATA, "", "Input STAR file (*.{star})", "Filename of the STAR file with the refined particle coordinates, e.g. from a previous 2D or 3D classification or auto-refine run.");
	joboptions["do_recenter"] = JobOption("Re-center refined coordinates? ", true, "If set to Yes, the input coordinates will be re-centered according to the refined origin offsets in the provided _data.star file .");
	joboptions["do_set_angpix"] = JobOption("Manually set pixel size? ", false, "If set to Yes, the rlnMagnification and rlnDetectorPixelSize will be set in the resulting STAR file. Only use this option if CTF information is NOT coming from the input coordinate STAR file(s). For example, because you decided not to estimate the CTF for your micrographs.");
	joboptions["angpix"] = JobOption("Pixel size (A)", 1, 0.3, 5, 0.1, "Provide the pixel size in Angstroms in the micrograph (so before any re-scaling).  If you provide input CTF parameters, then leave this value to the default of -1.");

	joboptions["extract_size"] = JobOption("Particle box size (pix):", 128, 64, 512, 8, "Size of the extracted particles (in pixels). This should be an even number!");
	joboptions["do_invert"] = JobOption("Invert contrast?", true, "If set to Yes, the contrast in the particles will be inverted.");

	joboptions["do_norm"] = JobOption("Normalize particles?", true, "If set to Yes, particles will be normalized in the way RELION prefers it.");
	joboptions["bg_diameter"] = JobOption("Diameter background circle (pix): ", -1, -1, 600, 10, "Particles will be normalized to a mean value of zero and a standard-deviation of one for all pixels in the background area.\
The background area is defined as all pixels outside a circle with this given diameter in pixels (before rescaling). When specifying a negative value, a default value of 75% of the Particle box size will be used.");
	joboptions["white_dust"] = JobOption("Stddev for white dust removal: ", -1, -1, 10, 0.1, "Remove very white pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
	joboptions["black_dust"] = JobOption("Stddev for black dust removal: ", -1, -1, 10, 0.1, "Remove very black pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
	joboptions["do_rescale"] = JobOption("Rescale particles?", false, "If set to Yes, particles will be re-scaled. Note that the particle diameter below will be in the down-scaled images.");
	joboptions["rescale"] = JobOption("Re-scaled size (pixels): ", 128, 64, 512, 8, "The re-scaled value needs to be an even number");

	joboptions["do_extract_helix"] = JobOption("Extract helical segments?", false, "Set to Yes if you want to extract helical segments. RELION (.star), EMAN2 (.box) and XIMDISP (.coords) formats of tube or segment coordinates are supported.");
	joboptions["helical_tube_outer_diameter"] = JobOption("Tube diameter (A): ", 200, 100, 1000, 10, "Outer diameter (in Angstroms) of helical tubes. \
This value should be slightly larger than the actual width of helical tubes.");
	joboptions["helical_bimodal_angular_priors"] = JobOption("Use bimodal angular priors?", true, "Normally it should be set to Yes and bimodal angular priors will be applied in the following classification and refinement jobs. \
Set to No if the 3D helix looks the same when rotated upside down.");
	joboptions["do_extract_helical_tubes"] = JobOption("Coordinates are start-end only?", true, "Set to Yes if you want to extract helical segments from manually picked tube coordinates (starting and end points of helical tubes in RELION, EMAN or XIMDISP format). \
Set to No if segment coordinates (RELION auto-picked results or EMAN / XIMDISP segments) are provided.");
	joboptions["do_cut_into_segments"] = JobOption("Cut helical tubes into segments?", true, "Set to Yes if you want to extract multiple helical segments with a fixed inter-box distance. \
If it is set to No, only one box at the center of each helical tube will be extracted.");
	joboptions["helical_nr_asu"] = JobOption("Number of asymmetrical units:", 1, 1, 100, 1, "Number of helical asymmetrical units in each segment box. This integer should not be less than 1. The inter-box distance (pixels) = helical rise (Angstroms) * number of asymmetrical units / pixel size (Angstroms). \
The optimal inter-box distance might also depend on the box size, the helical rise and the flexibility of the structure. In general, an inter-box distance of ~10% * the box size seems appropriate.");
	joboptions["helical_rise"] = JobOption("Helical rise (A):", 1, 0, 100, 0.01, "Helical rise in Angstroms. (Please click '?' next to the option above for details about how the inter-box distance is calculated.)");

}
bool RelionJob::getCommandsExtractJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "Extract", job_counter);
	std::string command;
	if (joboptions["nr_mpi"].getNumber() > 1)
		command="`which relion_preprocess_mpi`";
	else
		command="`which relion_preprocess`";

	// Input
	if (joboptions["star_mics"].getString() == "")
	{
		error_message = "ERROR: empty field for input STAR file...";
		return false;
	}
	command += " --i " + joboptions["star_mics"].getString();
	Node node(joboptions["star_mics"].getString(), joboptions["star_mics"].node_type);
	inputNodes.push_back(node);

	if (joboptions["do_reextract"].getBoolean())
	{
		if (joboptions["fndata_reextract"].getString() == "")
		{
			error_message = "ERROR: empty field for refined particles STAR file...";
			return false;
		}
		command += " --reextract_data_star " + joboptions["fndata_reextract"].getString();
		Node node2(joboptions["fndata_reextract"].getString(), joboptions["fndata_reextract"].node_type);
		inputNodes.push_back(node2);
		if (joboptions["do_recenter"].getBoolean())
		{
			command += " --recenter";
		}
	}
	else
	{
		FileName mysuffix = joboptions["coords_suffix"].getString();
		command += " --coord_dir " + mysuffix.beforeLastOf("/") + "/";
		command += " --coord_suffix " + (mysuffix.afterLastOf("/")).without("coords_suffix");
		Node node2(joboptions["coords_suffix"].getString(), joboptions["coords_suffix"].node_type);
		inputNodes.push_back(node2);
	}

	// Output
	FileName fn_ostar = outputname + "particles.star";
	Node node3(fn_ostar, NODE_PART_DATA);
	outputNodes.push_back(node3);
	command += " --part_star " + fn_ostar;

	command += " --part_dir " + outputname;
	command += " --extract";
	command += " --extract_size " + joboptions["extract_size"].getString();

	// Operate stuff
	// Get an integer number for the bg_radius
	RFLOAT bg_radius = (joboptions["bg_diameter"].getNumber() < 0.) ? 0.75 * joboptions["extract_size"].getNumber() : joboptions["bg_diameter"].getNumber();
	bg_radius /= 2.; // Go from diameter to radius
	if (joboptions["do_rescale"].getBoolean())
	{
		command += " --scale " + joboptions["rescale"].getString();
		bg_radius *= joboptions["rescale"].getNumber() / joboptions["extract_size"].getNumber();
	}
	if (joboptions["do_norm"].getBoolean())
	{
		// Get an integer number for the bg_radius
		bg_radius = (int)bg_radius;
		command += " --norm --bg_radius " + floatToString(bg_radius);
		command += " --white_dust " + joboptions["white_dust"].getString();
		command += " --black_dust " + joboptions["black_dust"].getString();
	}
	if (joboptions["do_invert"].getBoolean())
		command += " --invert_contrast ";

	if (joboptions["do_set_angpix"].getBoolean())
	{
		command += " --set_angpix " + joboptions["angpix"].getString();
	}

	// Helix
	if (joboptions["do_extract_helix"].getBoolean())
	{
		command += " --helix";
		command += " --helical_outer_diameter " + joboptions["helical_tube_outer_diameter"].getString();
		if (joboptions["helical_bimodal_angular_priors"].getBoolean())
			command += " --helical_bimodal_angular_priors";
		if (joboptions["do_extract_helical_tubes"].getBoolean())
		{
			command += " --helical_tubes";
			if (joboptions["do_cut_into_segments"].getBoolean())
			{
				command += " --helical_cut_into_segments";
				command += " --helical_nr_asu " + joboptions["helical_nr_asu"].getString();
				command += " --helical_rise " + joboptions["helical_rise"].getString();
			}
			else
				command += " --helical_nr_asu 1 --helical_rise 1";
		}
	}

	if (is_continue)
		command += " --only_extract_unfinished ";


	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();
	commands.push_back(command);

	if (joboptions["do_reextract"].getBoolean() || (joboptions["do_extract_helix"].getBoolean() && joboptions["do_extract_helical_tubes"].getBoolean()) )
	{
		// Also touch the suffix file. Do this after the first command had completed
		command = "echo " + joboptions["star_mics"].getString() + " > " +  outputname + "coords_suffix_extract.star";
		commands.push_back(command.c_str());

		Node node(outputname + "coords_suffix_extract.star", NODE_MIC_COORDS);
		outputNodes.push_back(node);

	}

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

void RelionJob::initialiseSortJob()
{
	hidden_name = ".gui_sort";

	joboptions["input_star"] = JobOption("Input particles to be sorted:", NODE_PART_DATA, "", "Input particles(*.{star})", "This STAR file should contain in-plane rotations, in-plane translations and a class number that were obtained by alignment (class2D/class3D or auto3D) OR auto-picking. A column called rlnParticleSelectZScore will be added to this same STAR file with the sorting result. This column can then be used in the display programs to sort the particles on.");
	joboptions["is_autopick"] = JobOption("Are these from an Extract job?", false, "Set the particles are from an Extract job (as opposed to a 2D/3D Classification or 3D auto-refine job), then you can only perform sorting if you used coordinates from Auto-picking in the particle extraction. In that case, provide the 2D references from the auto-pick run below.");
	joboptions["autopick_refs"] = JobOption("Autopicking references:", NODE_2DREFS, "", "References(*.{star})", "This STAR file should contain the 2D references that were used for the auto-picking.");

	joboptions["angpix_ref"] = JobOption("Pixel size in references (A)", -1, 0.3, 5, 0.1, "Pixel size in Angstroms for the provided reference images. This will be used to calculate the filters and the particle diameter in pixels. If a negative value is given here, the pixel size in the references will be assumed to be the same as the one in the micrographs, i.e. the particles that were used to make the references were not rescaled upon extraction.");
	joboptions["do_ctf"] = JobOption("Are References CTF corrected?", true, "Set to Yes if the references were created with CTF-correction inside RELION. \n ");
	joboptions["do_ignore_first_ctfpeak"] = JobOption("Ignore CTFs until first peak?", false,"Set this to Yes, only if this option was also used to generate the references.");


}

bool RelionJob::getCommandsSortJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "Sort", job_counter);
	std::string command;

	if (joboptions["nr_mpi"].getNumber() > 1)
		command="`which relion_particle_sort_mpi`";
	else
		command="`which relion_particle_sort`";

	if (joboptions["input_star"].getString() == "")
	{
		error_message = "ERROR: empty field for input STAR file...";
		return false;
	}
	command += " --i " + joboptions["input_star"].getString();
	Node node(joboptions["input_star"].getString(), joboptions["input_star"].node_type);
	inputNodes.push_back(node);

	if (joboptions["angpix_ref"].getNumber() > 0.)
		command += " --angpix_ref " + joboptions["angpix_ref"].getString();

	// Determine the --ref automatically, from the particle input filename
	FileName fn_ref, fn_in = joboptions["input_star"].getString();
	int node_type;
	if (fn_in.contains("_data.star") && (fn_in.contains("Class2D/") || fn_in.contains("Class3D/")) )
	{
		fn_ref = fn_in.without("_data.star") + "_model.star";
		node_type= NODE_MODEL;
	}
	else if (fn_in.contains("_data.star") && fn_in.contains("Refine3D/"))
	{
		if (fn_in.contains("_it0") || fn_in.contains("_it1"))
			fn_ref = fn_in.without("_data.star") + "_half1_model.star";
		else
			fn_ref = fn_in.without("_data.star") + "_model.star";
		node_type= NODE_MODEL;
	}
	else if (fn_in.contains("Extract/"))
	{
		if (joboptions["is_autopick"].getBoolean())
		{

			if (joboptions["autopick_refs"].getString() == "")
			{
				error_message = "ERROR: empty field for autopicking references. This is compulsory for Extract jobs...";
				return false;
			}

			fn_ref = joboptions["autopick_refs"].getString();
			node_type= NODE_2DREFS;
		}
		else
		{
			error_message = "ERROR: these particles are from an Extract job. Without auto-picking references you cannot run sorting!";
			return false;
		}
	}
	command += " --ref " + fn_ref;
	Node node2(fn_ref, node_type);
	inputNodes.push_back(node2);

	command += " --o " + outputname + "particles_sort.star";
	Node node3(outputname + "particles_sort.star", NODE_PART_DATA);
	outputNodes.push_back(node3);

	if (joboptions["do_ctf"].getBoolean())
	{
		command += " --ctf ";
		if (joboptions["do_ignore_first_ctfpeak"].getBoolean())
			command += " --ctf_intact_first_peak ";
	}

	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

void RelionJob::initialiseSelectJob()
{

	hidden_name = ".gui_select";

	joboptions["fn_model"] = JobOption("Select classes from model.star:", NODE_MODEL, "", "STAR files (*.star)", "A _model.star file from a previous 2D or 3D classification run to select classes from.");
	joboptions["fn_mic"] = JobOption("OR select from micrographs.star:", NODE_MICS, "", "STAR files (*.star)", "A micrographs.star file to select micrographs from.");
	joboptions["fn_data"] = JobOption("OR select from particles.star:", NODE_PART_DATA, "", "STAR files (*.star)", "A particles.star file to select individual particles from.");
	joboptions["fn_coords"] = JobOption("OR select from picked coords:", NODE_MIC_COORDS, "", "STAR files (coords_suffix*.star)", "A coordinate suffix .star file to select micrographs while inspecting coordinates (and/or CTFs).");

	joboptions["do_recenter"] = JobOption("Re-center the class averages?", true, "This option is only used when selecting particles from 2D classes. The selected class averages will all re-centered on their center-of-mass. This is useful when you plane to use these class averages as templates for auto-picking.");
	joboptions["do_regroup"] = JobOption("Regroup the particles?", false, "If set to Yes, then the program will regroup the selected particles in 'more-or-less' the number of groups indicated below. For re-grouping from individual particle _data.star files, a _model.star file with the same prefix should exist, i.e. the particle star file should be generated by relion_refine");
	joboptions["nr_groups"] = JobOption("Approximate nr of groups: ", 1, 50, 20, 1, "It is normal that the actual number of groups may deviate a little from this number. ");

}

bool RelionJob::getCommandsSelectJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "Select", job_counter);
	std::string command;
	command="`which relion_display`";

	// I/O
	if (joboptions["fn_model"].getString() != "")
	{
		command += " --gui --i " + joboptions["fn_model"].getString();
		Node node(joboptions["fn_model"].getString(), joboptions["fn_model"].node_type);
		inputNodes.push_back(node);

		FileName fn_parts = outputname+"particles.star";
		command += " --allow_save --fn_parts " + fn_parts;
		Node node2(fn_parts, NODE_PART_DATA);
		outputNodes.push_back(node2);

		// Only save the 2D class averages for 2D jobs
		FileName fnt = joboptions["fn_model"].getString();
		if (fnt.contains("Class2D/"))
		{
			FileName fn_imgs = outputname+"class_averages.star";
			command += " --fn_imgs " + fn_imgs;
			Node node3(fn_imgs, NODE_2DREFS);
			outputNodes.push_back(node3);

			if (joboptions["do_recenter"].getBoolean())
			{
				command += " --recenter ";
			}
		}
	}
	else if (joboptions["fn_mic"].getString() != "")
	{
		command += " --gui --i " + joboptions["fn_mic"].getString();
		Node node(joboptions["fn_mic"].getString(), joboptions["fn_mic"].node_type);
		inputNodes.push_back(node);

		FileName fn_mics = outputname+"micrographs.star";
		command += " --allow_save --fn_imgs " + fn_mics;
		Node node2(fn_mics, NODE_MICS);
		outputNodes.push_back(node2);
	}
	else if (joboptions["fn_data"].getString() != "")
	{
		command += " --gui --i " + joboptions["fn_data"].getString();
		Node node(joboptions["fn_data"].getString(), joboptions["fn_data"].node_type);
		inputNodes.push_back(node);

		FileName fn_parts = outputname+"particles.star";
		command += " --allow_save --fn_imgs " + fn_parts;
		Node node2(fn_parts, NODE_PART_DATA);
		outputNodes.push_back(node2);
	}
	else if  (joboptions["fn_coords"].getString() != "")
	{

    	RelionJob manualpickjob;

    	FileName fn_job = ".gui_manualpick";
		bool iscont=false;
		if (exists(fn_job+"run.job"))
		{
			manualpickjob.read(fn_job.c_str(), iscont, true); // true means do initialise
		}
		else
		{
			error_message = "You need to save 'Manual picking' job settings (using the Jobs menu) before you can display coordinate files.";
			return false;
		}

		// Get the name of the micrograph STAR file from reading the suffix file
	    FileName fn_suffix = joboptions["fn_coords"].getString();
		FileName fn_star;
	    if (is_continue)
	    {
	    	fn_star = outputname + "micrographs_selected.star";
	    }
	    else
	    {
	    	std::ifstream in(fn_suffix.data(), std::ios_base::in);
	    	in >> fn_star ;
	    	in.close();
	    }
	    FileName fn_dirs = fn_suffix.beforeLastOf("/")+"/";
		fn_suffix = fn_suffix.afterLastOf("/").without("coords_suffix_");
		fn_suffix = fn_suffix.withoutExtension();

		// Launch the manualpicker...
		command="`which relion_manualpick` --i " + fn_star;
		Node node4(joboptions["fn_coords"].getString(), joboptions["fn_coords"].node_type);
		inputNodes.push_back(node4);

		command += " --odir " + fn_dirs;
		command += " --pickname " + fn_suffix;

		// The output selection
		FileName fn_outstar = outputname + "micrographs_selected.star";
		Node node3(fn_outstar, NODE_MICS);
		outputNodes.push_back(node3);
		command += " --allow_save  --selection " + fn_outstar;

		// All the stuff from the saved manualpickjob
		command += " --scale " + manualpickjob.joboptions["micscale"].getString();
		command += " --sigma_contrast " + manualpickjob.joboptions["sigma_contrast"].getString();
		command += " --black " + manualpickjob.joboptions["black_val"].getString();
		command += " --white " + manualpickjob.joboptions["white_val"].getString();

		if (manualpickjob.joboptions["lowpass"].getNumber() > 0.)
			command += " --lowpass " + manualpickjob.joboptions["lowpass"].getString();
		if (manualpickjob.joboptions["highpass"].getNumber() > 0.)
			command += " --highpass " + manualpickjob.joboptions["highpass"].getString();
		if (manualpickjob.joboptions["angpix"].getNumber() > 0.)
			command += " --angpix " + manualpickjob.joboptions["angpix"].getString();

		command += " --ctf_scale " + manualpickjob.joboptions["ctfscale"].getString();

		command += " --particle_diameter " + manualpickjob.joboptions["diameter"].getString();


		if (manualpickjob.joboptions["do_color"].getBoolean())
		{
			command += " --color_label " + manualpickjob.joboptions["color_label"].getString();
			command += " --blue " + manualpickjob.joboptions["blue_value"].getString();
			command += " --red " + manualpickjob.joboptions["red_value"].getString();
			if (manualpickjob.joboptions["fn_color"].getString().length() > 0)
				command += " --color_star " + manualpickjob.joboptions["fn_color"].getString();
		}

		// Other arguments for extraction
		command += " " + manualpickjob.joboptions["other_args"].getString() + " &";
	}
	else
	{
		// Nothing was selected...
		error_message = "Please select an input file.";
		return false;
	}

	// Re-grouping
	if (joboptions["do_regroup"].getBoolean() && joboptions["fn_coords"].getString() == "")
	{
		if (joboptions["fn_model"].getString() == "")
		{
			error_message = "Re-grouping only works for model.star files...";
			return false;
    	}
        command += " --regroup " + joboptions["nr_groups"].getString();
	}

	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

void RelionJob::initialiseClass2DJob()
{

	hidden_name = ".gui_class2d";

	joboptions["fn_img"] = JobOption("Input images STAR file:", NODE_PART_DATA, "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");
	joboptions["fn_cont"] = JobOption("Continue from here: ", std::string(""), "STAR Files (*_optimiser.star)", "CURRENT_ODIR",  "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");

	joboptions["do_ctf_correction"] = JobOption("Do CTF-correction?", true, "If set to Yes, CTFs will be corrected inside the MAP refinement. \
The resulting algorithm intrinsically implements the optimal linear, or Wiener filter. \
Note that CTF parameters for all images need to be given in the input STAR file. \
The command 'relion_refine --print_metadata_labels' will print a list of all possible metadata labels for that STAR file. \
See the RELION Wiki for more details.\n\n Also make sure that the correct pixel size (in Angstrom) is given above!)");
	joboptions["ctf_phase_flipped"] = JobOption("Have data been phase-flipped?", false, "Set this to Yes if the images have been \
ctf-phase corrected during the pre-processing steps. \
Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION, \
as this can be done inside the internal CTF-correction. \
However, if the phases have been flipped, you should tell the program about it by setting this option to Yes.");
	joboptions["ctf_intact_first_peak"] = JobOption("Ignore CTFs until first peak?", false, "If set to Yes, then CTF-amplitude correction will \
only be performed from the first peak of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution. \
Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results. \
Therefore, this option is not generally recommended: try increasing amplitude contrast (in your input STAR file) first!");

	joboptions["nr_classes"] = JobOption("Number of classes:", 1, 1, 50, 1, "The number of classes (K) for a multi-reference refinement. \
These classes will be made in an unsupervised manner from a single reference by division of the data into random subsets during the first iteration.");
	joboptions["tau_fudge"] = JobOption("Regularisation parameter T:", 2 , 0.1, 10, 0.1, "Bayes law strictly determines the relative weight between \
the contribution of the experimental data and the prior. However, in practice one may need to adjust this weight to put slightly more weight on \
the experimental data to allow optimal results. Values greater than 1 for this regularisation parameter (T in the JMB2011 paper) put more \
weight on the experimental data. Values around 2-4 have been observed to be useful for 3D refinements, values of 1-2 for 2D refinements. \
Too small values yield too-low resolution structures; too high values result in over-estimated resolutions, mostly notable by the apparition of high-frequency noise in the references.");
	joboptions["nr_iter"] = JobOption("Number of iterations:", 25, 1, 50, 1, "Number of iterations to be performed. \
Note that the current implementation of 2D class averaging and 3D classification does NOT comprise a convergence criterium. \
Therefore, the calculations will need to be stopped by the user if further iterations do not yield improvements in resolution or classes. \n\n \
Also note that upon restarting, the iteration number continues to be increased, starting from the final iteration in the previous run. \
The number given here is the TOTAL number of iterations. For example, if 10 iterations have been performed previously and one restarts to perform \
an additional 5 iterations (for example with a finer angular sampling), then the number given here should be 10+5=15.");

	joboptions["do_subsets"] = JobOption("Use subsets for initial updates?", false, "If set to True, multiple maximization updates (as many as defined by the 'Number of subset updates') will be performed during the first iteration(s): each time after the number of particles in a subset has been processed. \
By using subsets with much fewer particles than the entire data set, the initial updates will be much faster, while the very low resolution class averages will not be notably worse than with the entire data set. \n\n \
This will greatly speed up 2D classifications with very many (hundreds of thousands or more) particles. A useful subset size is probably in the order of ten thousand particles. If the data set only comprises (tens of) thousands of particles, this option may be less useful.");
	joboptions["subset_size"] = JobOption("Initial subset size:", 10000, 1000, 50000, 1000, "Number of individual particles after which one will perform a maximization update in the first iteration(s). \
A useful subset size is probably in the order of ten thousand particles.");
	joboptions["max_subsets"] = JobOption("Number of subset updates:", 3, 1, 10, 1, "This option is only used when a positive number is given for the 'Initial subset size'. In that case, in the first iteration, maximization updates are performed over a smaller subset of the particles to speed up calculations.\
Useful values are probably in the range of 2-5 subset updates. Using more might speed up further, but with the risk of affecting the results. If the number of subsets times the subset size is larger than the number of particles in the data set, then more than 1 iteration will be split into subsets.");
	joboptions["particle_diameter"] = JobOption("Mask diameter (A):", 200, 0, 1000, 10, "The experimental images will be masked with a soft \
circular mask with this diameter. Make sure this radius is not set too small because that may mask away part of the signal! \
If set to a value larger than the image size no masking will be performed.\n\n\
The same diameter will also be used for a spherical mask of the reference structures if no user-provided mask is specified.");
	joboptions["do_zero_mask"] = JobOption("Mask individual particles with zeros?", true, "If set to Yes, then in the individual particles, \
the area outside a circle with the radius of the particle will be set to zeros prior to taking the Fourier transform. \
This will remove noise and therefore increase sensitivity in the alignment and classification. However, it will also introduce correlations \
between the Fourier components that are not modelled. When set to No, then the solvent area is filled with random noise, which prevents introducing correlations.\
High-resolution refinements (e.g. ribosomes or other large complexes in 3D auto-refine) tend to work better when filling the solvent area with random noise (i.e. setting this option to No), refinements of smaller complexes and most classifications go better when using zeros (i.e. setting this option to Yes).");
	joboptions["highres_limit"] = JobOption("Limit resolution E-step to (A): ", -1, -1, 20, 1, "If set to a positive number, then the expectation step (i.e. the alignment) will be done only including the Fourier components up to this resolution (in Angstroms). \
This is useful to prevent overfitting, as the classification runs in RELION are not to be guaranteed to be 100% overfitting-free (unlike the 3D auto-refine with its gold-standard FSC). In particular for very difficult data sets, e.g. of very small or featureless particles, this has been shown to give much better class averages. \
In such cases, values in the range of 7-12 Angstroms have proven useful.");

	joboptions["dont_skip_align"] = JobOption("Perform image alignment?", true, "If set to No, then rather than \
performing both alignment and classification, only classification will be performed. This allows the use of very focused masks.\
This requires that the optimal orientations of all particles are already stored in the input STAR file. ");
	joboptions["psi_sampling"] = JobOption("In-plane angular sampling:", 6., 0.5, 20, 0.5, "The sampling rate for the in-plane rotation angle (psi) in degrees. \
Using fine values will slow down the program. Recommended value for most 2D refinements: 5 degrees.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");
	joboptions["offset_range"] = JobOption("Offset search range (pix):", 5, 0, 30, 1, "Probabilities will be calculated only for translations \
in a circle with this radius (in pixels). The center of this circle changes at every iteration and is placed at the optimal translation \
for each image in the previous iteration.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");
	joboptions["offset_step"] = JobOption("Offset search step (pix):", 1, 0.1, 5, 0.1, "Translations will be sampled with this step-size (in pixels).\
Translational sampling is also done using the adaptive approach. \
Therefore, if adaptive=1, the translations will first be evaluated on a 2x coarser grid.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");

	joboptions["do_helix"] = JobOption("Classify 2D helical segments?", false, "Set to Yes if you want to classify 2D helical segments. Note that the helical segments should come with priors of psi angles");
	joboptions["helical_tube_outer_diameter"] = JobOption("Tube diameter (A): ", 200, 100, 1000, 10, "Outer diameter (in Angstroms) of helical tubes. \
This value should be slightly larger than the actual width of the tubes. You may want to copy the value from previous particle extraction job. \
If negative value is provided, this option is disabled and ordinary circular masks will be applied. Sometimes '--dont_check_norm' option is useful to prevent errors in normalisation of helical segments.");
	joboptions["do_bimodal_psi"] = JobOption("Do bimodal angular searches?", true, "Do bimodal search for psi angles? \
Set to Yes if you want to classify 2D helical segments with priors of psi angles. The priors should be bimodal due to unknown polarities of the segments. \
Set to No if the 3D helix looks the same when rotated upside down. If it is set to No, ordinary angular searches will be performed.\n\nThis option will be invalid if you choose not to perform image alignment on 'Sampling' tab.");
	joboptions["range_psi"] = JobOption("Angular search range - psi (deg):", 6, 3, 30, 1, "Local angular searches will be performed \
within +/- the given amount (in degrees) from the psi priors estimated through helical segment picking. \
A range of 15 degrees is the same as sigma = 5 degrees. Note that the ranges of angular searches should be much larger than the sampling.\
\n\nThis option will be invalid if you choose not to perform image alignment on 'Sampling' tab.");

	joboptions["nr_pool"] = JobOption("Number of pooled particles:", 3, 1, 16, 1, "Particles are processed in individual batches by MPI slaves. During each batch, a stack of particle images is only opened and closed once to improve disk access times. \
All particle images of a single batch are read into memory together. The size of these batches is at least one particle per thread used. The nr_pooled_particles parameter controls how many particles are read together for each thread. If it is set to 3 and one uses 8 threads, batches of 3x8=24 particles will be read together. \
This may improve performance on systems where disk access, and particularly metadata handling of disk access, is a problem. It has a modest cost of increased RAM usage.");
	joboptions["do_parallel_discio"] = JobOption("Use parallel disc I/O?", true, "If set to Yes, all MPI slaves will read images from disc. \
Otherwise, only the master will read images and send them through the network to the slaves. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many slaves reading in parallel.");
	joboptions["do_preread_images"] = JobOption("Pre-read all particles into RAM?", false, "If set to Yes, all particle images will be read into computer memory, which will greatly speed up calculations on systems with slow disk access. However, one should of course be careful with the amount of RAM available. \
Because particles are read in float-precision, it will take ( N * box_size * box_size * 4 / (1024 * 1024 * 1024) ) Giga-bytes to read N particles into RAM. For 100 thousand 200x200 images, that becomes 15Gb, or 60 Gb for the same number of 400x400 particles. \
Remember that running a single MPI slave on each node that runs as many threads as available cores will have access to all available RAM. \n \n If parallel disc I/O is set to No, then only the master reads all particles into RAM and sends those particles through the network to the MPI slaves during the refinement iterations.");
	joboptions["scratch_dir"] = JobOption("Copy particles to scratch directory:", std::string(""), "If a directory is provided here, then the job will create a sub-directory in it called relion_volatile. If that relion_volatile directory already exists, it will be wiped. Then, the program will copy all input particles into a large stack inside the relion_volatile subdirectory. \
Provided this directory is on a fast local drive (e.g. an SSD drive), processing in all the iterations will be faster. If the job finishes correctly, the relion_volatile directory will be wiped. If the job crashes, you may want to remove it yourself.");
	joboptions["do_combine_thru_disc"] = JobOption("Combine iterations through disc?", false, "If set to Yes, at the end of every iteration all MPI slaves will write out a large file with their accumulated results. The MPI master will read in all these files, combine them all, and write out a new file with the combined results. \
All MPI salves will then read in the combined results. This reduces heavy load on the network, but increases load on the disc I/O. \
This will affect the time it takes between the progress-bar in the expectation step reaching its end (the mouse gets to the cheese) and the start of the ensuing maximisation step. It will depend on your system setup which is most efficient.");

	joboptions["use_gpu"] = JobOption("Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','. For example: '0,0:1,1:0,0:1,1'");

}

bool RelionJob::getCommandsClass2DJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "Class2D", job_counter);
	std::string command;

	if (joboptions["nr_mpi"].getNumber() > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";

    FileName fn_run = "run";
	if (is_continue)
    {
		if (joboptions["fn_cont"].getString() == "")
		{
			error_message = "ERROR: empty field for continuation STAR file...";
			return false;
		}
		int pos_it = joboptions["fn_cont"].getString().rfind("_it");
		int pos_op = joboptions["fn_cont"].getString().rfind("_optimiser");
		if (pos_it < 0 || pos_op < 0)
		{
			error_message = "Warning: invalid optimiser.star filename provided for continuation run!";
			return false;
		}
		int it = (int)textToFloat((joboptions["fn_cont"].getString().substr(pos_it+3, 6)).c_str());
		fn_run += "_ct" + floatToString(it);
		command += " --continue " + joboptions["fn_cont"].getString();
    }

    command += " --o " + outputname + fn_run;
	outputNodes = getOutputNodesRefine(outputname + fn_run, (int)joboptions["nr_iter"].getNumber(), (int)joboptions["nr_classes"].getNumber(), 2, 1);

	if (!is_continue)
	{
		if (joboptions["fn_img"].getString() == "")
		{
			error_message = "ERROR: empty field for input STAR file...";
			return false;
		}
		command += " --i " + joboptions["fn_img"].getString();
		Node node(joboptions["fn_img"].getString(), joboptions["fn_img"].node_type);
		inputNodes.push_back(node);
	}

	// Always do compute stuff
	if (!joboptions["do_combine_thru_disc"].getBoolean())
		command += " --dont_combine_weights_via_disc";
	if (!joboptions["do_parallel_discio"].getBoolean())
		command += " --no_parallel_disc_io";
	if (joboptions["do_preread_images"].getBoolean())
		command += " --preread_images " ;
	else if (joboptions["scratch_dir"].getString() != "")
            command += " --scratch_dir " +  joboptions["scratch_dir"].getString();
        command += " --pool " + joboptions["nr_pool"].getString();

	// CTF stuff
	if (!is_continue)
	{

		if (joboptions["do_ctf_correction"].getBoolean())
		{
			command += " --ctf ";
			if (joboptions["ctf_phase_flipped"].getBoolean())
				command += " --ctf_phase_flipped ";
			if (joboptions["ctf_intact_first_peak"].getBoolean())
				command += " --ctf_intact_first_peak ";
		}
	}

	// Optimisation
	command += " --iter " + joboptions["nr_iter"].getString();
	if (joboptions["do_subsets"].getBoolean())
	{
		command += " --write_subsets 1 --subset_size " + joboptions["subset_size"].getString();
		command += " --max_subsets " + joboptions["max_subsets"].getString();
	}

	command += " --tau2_fudge " + joboptions["tau_fudge"].getString();
    command += " --particle_diameter " + joboptions["particle_diameter"].getString();
	if (!is_continue)
	{
		command += " --K " + joboptions["nr_classes"].getString();
		// Always flatten the solvent
		command += " --flatten_solvent ";
		if (joboptions["do_zero_mask"].getBoolean())
			command += " --zero_mask ";
		if (joboptions["highres_limit"].getNumber() > 0)
			command += " --strict_highres_exp " + joboptions["highres_limit"].getString();
	}

	// Sampling
	int iover = 1;
	command += " --oversampling " + floatToString((float)iover);

	if (!joboptions["dont_skip_align"].getBoolean())
	{
		command += " --skip_align ";
	}
	else
	{
		// The sampling given in the GUI will be the oversampled one!
		command += " --psi_step " + floatToString(joboptions["psi_sampling"].getNumber() * pow(2., iover));
		// Offset range
		command += " --offset_range " + joboptions["offset_range"].getString();
		// The sampling given in the GUI will be the oversampled one!
		command += " --offset_step " + floatToString(joboptions["offset_step"].getNumber() * pow(2., iover));
	}

	// Helix
	if (joboptions["do_helix"].getBoolean())
	{
		command += " --helical_outer_diameter " + joboptions["helical_tube_outer_diameter"].getString();

		if (joboptions["dont_skip_align"].getBoolean())
		{
			if (joboptions["do_bimodal_psi"].getBoolean())
				command += " --bimodal_psi";

			RFLOAT val = joboptions["range_psi"].getNumber();
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_psi " + floatToString(val / 3.);
		}
	}

	// Always do norm and scale correction
	if (!is_continue)
		command += " --norm --scale ";

	// Running stuff
	command += " --j " + joboptions["nr_threads"].getString();

	// GPU-stuff
	if (joboptions["use_gpu"].getBoolean())
	{
		command += " --gpu \"" + joboptions["gpu_ids"].getString() +"\"";
	}

	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);
	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

// Constructor for initial model job
void RelionJob::initialiseInimodelJob()
{

	hidden_name = ".gui_inimodel";

	joboptions["fn_img"] = JobOption("Input images STAR file:", NODE_PART_DATA, "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \
In SGD, it is very important that there are particles from enough different orientations. One only needs a few thousand to 10k particles. When selecting good 2D classes in the Subset Selection jobtype, use the option to select a maximum number of particles from each class to generate more even angular distributions for SGD.\
\n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");
	joboptions["fn_cont"] = JobOption("Continue from here: ", std::string(""), "STAR Files (*_optimiser.star)", "CURRENT_ODIR", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");

	joboptions["nr_iter"] = JobOption("Number of iterations:", 1, 1, 10, 1, "Number of iterations to be performed. Often 1 or 2 iterations with approximately ten thousand particles, or 5-10 iterations with several thousand particles is enough.");
	joboptions["sgd_subset_size"] = JobOption("SGD subset size:", 200, 100, 1000, 100, "How many particles will be processed for each SGD step. Often 200 seems to work well.");
	joboptions["sgd_write_subsets"] = JobOption("Write-out frequency subsets:", 10, -1, 25, 1, "Every how many subsets do you want to write the model to disk. Negative value means only write out model after entire iteration. ");
	joboptions["sgd_highres_limit"] = JobOption("Limit resolution SGD to (A): ", 20, -1, 40, 1, "If set to a positive number, then the SGD will be done only including the Fourier components up to this resolution (in Angstroms). \
This is essential in SGD, as there is very little regularisation, i.e. overfitting will start to happen very quickly. \
Values in the range of 15-30 Angstroms have proven useful.");
	joboptions["sgd_sigma2fudge_halflife"] = JobOption("SGD increased noise variance half-life:", -1, -100, 10000, 100, "When set to a positive value, the initial estimates of the noise variance will internally be multiplied by 8, and then be gradually reduced, \
having 50% after this many particles have been processed. By default, this option is switched off by setting this value to a negative number. \
In some difficult cases, switching this option on helps. In such cases, values around 1000 have found to be useful. Change the factor of eight with the additional argument --sgd_sigma2fudge_ini");

	//joboptions["nr_classes"] = JobOption("Number of classes:", 1, 1, 50, 1, "The number of classes (K) for a multi-reference refinement. \
These classes will be made in an unsupervised manner from a single reference by division of the data into random subsets during the first iteration.");
	joboptions["sym_name"] = JobOption("Symmetry:", std::string("C1"), "Initial SGD runs are often performed in C1. If a particle is confirmed to have symmetry, the SGD can also be repeated with the corresponding \
point group symmetry. That has the advantage that the symetry axes in the reference will be aligned correctly.");
	joboptions["particle_diameter"] = JobOption("Mask diameter (A):", 200, 0, 1000, 10, "The experimental images will be masked with a soft \
circular mask with this diameter. Make sure this radius is not set too small because that may mask away part of the signal! \
If set to a value larger than the image size no masking will be performed.\n\n\
The same diameter will also be used for a spherical mask of the reference structures if no user-provided mask is specified.");
	//joboptions["do_zero_mask"] = JobOption("Mask individual particles with zeros?", true, "If set to Yes, then in the individual particles, \
the area outside a circle with the radius of the particle will be set to zeros prior to taking the Fourier transform. \
This will remove noise and therefore increase sensitivity in the alignment and classification. However, it will also introduce correlations \
between the Fourier components that are not modelled. When set to No, then the solvent area is filled with random noise, which prevents introducing correlations.\
High-resolution refinements (e.g. ribosomes or other large complexes in 3D auto-refine) tend to work better when filling the solvent area with random noise (i.e. setting this option to No), refinements of smaller complexes and most classifications go better when using zeros (i.e. setting this option to Yes).");

	joboptions["do_ctf_correction"] = JobOption("Do CTF-correction?", true, "If set to Yes, CTFs will be corrected inside the MAP refinement. \
The resulting algorithm intrinsically implements the optimal linear, or Wiener filter. \
Note that CTF parameters for all images need to be given in the input STAR file. \
The command 'relion_refine --print_metadata_labels' will print a list of all possible metadata labels for that STAR file. \
See the RELION Wiki for more details.\n\n Also make sure that the correct pixel size (in Angstrom) is given above!)");
	joboptions["ctf_phase_flipped"] = JobOption("Have data been phase-flipped?", false, "Set this to Yes if the images have been \
ctf-phase corrected during the pre-processing steps. \
Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION, \
as this can be done inside the internal CTF-correction. \
However, if the phases have been flipped, you should tell the program about it by setting this option to Yes.");
	joboptions["ctf_intact_first_peak"] = JobOption("Ignore CTFs until first peak?", false, "If set to Yes, then CTF-amplitude correction will \
only be performed from the first peak of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution. \
Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results. \
Therefore, this option is not generally recommended: try increasing amplitude contrast (in your input STAR file) first!");


	joboptions["sampling"] = JobOption("Angular sampling interval:", RADIO_SAMPLING, 1, "There are only a few discrete \
angular samplings possible because we use the HealPix library to generate the sampling of the first two Euler angles on the sphere. \
The samplings are approximate numbers and vary slightly over the sphere.\n\n For initial model generation at low resolutions, coarser angular samplings can often be used than in normal 3D classifications/refinements, e.g. 15 degrees. ");
	joboptions["offset_range"] = JobOption("Offset search range (pix):", 6, 0, 30, 1, "Probabilities will be calculated only for translations \
in a circle with this radius (in pixels). The center of this circle changes at every iteration and is placed at the optimal translation \
for each image in the previous iteration.\n\n");
	joboptions["offset_step"] = JobOption("Offset search step (pix):", 2, 0.1, 5, 0.1, "Translations will be sampled with this step-size (in pixels).\
Translational sampling is also done using the adaptive approach. \
Therefore, if adaptive=1, the translations will first be evaluated on a 2x coarser grid.\n\n ");

	joboptions["do_parallel_discio"] = JobOption("Use parallel disc I/O?", true, "If set to Yes, all MPI slaves will read their own images from disc. \
Otherwise, only the master will read images and send them through the network to the slaves. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many slaves reading in parallel.");
	joboptions["nr_pool"] = JobOption("Number of pooled particles:", 3, 1, 16, 1, "Particles are processed in individual batches by MPI slaves. During each batch, a stack of particle images is only opened and closed once to improve disk access times. \
All particle images of a single batch are read into memory together. The size of these batches is at least one particle per thread used. The nr_pooled_particles parameter controls how many particles are read together for each thread. If it is set to 3 and one uses 8 threads, batches of 3x8=24 particles will be read together. \
This may improve performance on systems where disk access, and particularly metadata handling of disk access, is a problem. It has a modest cost of increased RAM usage.");
	joboptions["do_preread_images"] = JobOption("Pre-read all particles into RAM?", false, "If set to Yes, all particle images will be read into computer memory, which will greatly speed up calculations on systems with slow disk access. However, one should of course be careful with the amount of RAM available. \
Because particles are read in float-precision, it will take ( N * box_size * box_size * 4 / (1024 * 1024 * 1024) ) Giga-bytes to read N particles into RAM. For 100 thousand 200x200 images, that becomes 15Gb, or 60 Gb for the same number of 400x400 particles. \
Remember that running a single MPI slave on each node that runs as many threads as available cores will have access to all available RAM. \n \n If parallel disc I/O is set to No, then only the master reads all particles into RAM and sends those particles through the network to the MPI slaves during the refinement iterations.");
	joboptions["scratch_dir"] = JobOption("Copy particles to scratch directory:", std::string(""), "If a directory is provided here, then the job will create a sub-directory in it called relion_volatile. If that relion_volatile directory already exists, it will be wiped. Then, the program will copy all input particles into a large stack inside the relion_volatile subdirectory. \
Provided this directory is on a fast local drive (e.g. an SSD drive), processing in all the iterations will be faster. If the job finishes correctly, the relion_volatile directory will be wiped. If the job crashes, you may want to remove it yourself.");
	joboptions["do_combine_thru_disc"] = JobOption("Combine iterations through disc?", false, "If set to Yes, at the end of every iteration all MPI slaves will write out a large file with their accumulated results. The MPI master will read in all these files, combine them all, and write out a new file with the combined results. \
All MPI salves will then read in the combined results. This reduces heavy load on the network, but increases load on the disc I/O. \
This will affect the time it takes between the progress-bar in the expectation step reaching its end (the mouse gets to the cheese) and the start of the ensuing maximisation step. It will depend on your system setup which is most efficient.");

	joboptions["use_gpu"] = JobOption("Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','. For example: '0,0:1,1:0,0:1,1'");

}

bool RelionJob::getCommandsInimodelJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();

	initialisePipeline(outputname, "InitialModel", job_counter);

	std::string command;
	if (joboptions["nr_mpi"].getNumber() > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";

    FileName fn_run = "run";
	if (is_continue)
    {
		if (joboptions["fn_cont"].getString() == "")
		{
			error_message = "ERROR: empty field for continuation STAR file...";
			return false;
		}
		int pos_it = joboptions["fn_cont"].getString().rfind("_it");
		int pos_op = joboptions["fn_cont"].getString().rfind("_optimiser");
		if (pos_it < 0 || pos_op < 0)
			std::cerr << "Warning: invalid optimiser.star filename provided for continuation run: " << joboptions["fn_cont"].getString() << std::endl;
		int it = (int)textToFloat((joboptions["fn_cont"].getString().substr(pos_it+3, 6)).c_str());
		fn_run += "_ct" + floatToString(it);
		command += " --continue " + joboptions["fn_cont"].getString();
    }

    command += " --o " + outputname + fn_run;
	outputNodes = getOutputNodesRefine(outputname + fn_run, (int)joboptions["nr_iter"].getNumber(), 1, 3, 1);

	command += " --sgd ";
	command += " --subset_size " + joboptions["sgd_subset_size"].getString();
	command += " --strict_highres_sgd " + joboptions["sgd_highres_limit"].getString();
	command += " --write_subsets " + joboptions["sgd_write_subsets"].getString();

	if (!is_continue)
	{
		if (joboptions["fn_img"].getString() == "")
		{
			error_message = "ERROR: empty field for input STAR file...";
			return false;
		}
		command += " --denovo_3dref --i " + joboptions["fn_img"].getString();
		Node node(joboptions["fn_img"].getString(), joboptions["fn_img"].node_type);
		inputNodes.push_back(node);

		// CTF stuff
		if (joboptions["do_ctf_correction"].getBoolean())
		{
			command += " --ctf";
			if (joboptions["ctf_phase_flipped"].getBoolean())
				command += " --ctf_phase_flipped";
			if (joboptions["ctf_intact_first_peak"].getBoolean())
				command += " --ctf_intact_first_peak";
		}

		//command += " --K " + joboptions["nr_classes"].getString();
		command += " --sym " + joboptions["sym_name"].getString();

		//if (joboptions["do_zero_mask"].getBoolean())
			command += " --zero_mask";
	}

	// Always do compute stuff
	if (!joboptions["do_combine_thru_disc"].getBoolean())
		command += " --dont_combine_weights_via_disc";
	if (!joboptions["do_parallel_discio"].getBoolean())
		command += " --no_parallel_disc_io";
	if (joboptions["do_preread_images"].getBoolean())
            command += " --preread_images " ;
	else if (joboptions["scratch_dir"].getString() != "")
            command += " --scratch_dir " +  joboptions["scratch_dir"].getString();
        command += " --pool " + joboptions["nr_pool"].getString();

	// Optimisation
	command += " --iter " + joboptions["nr_iter"].getString();
    command += " --particle_diameter " + joboptions["particle_diameter"].getString();

    // Sampling
    int iover = 1;
	command += " --oversampling " + floatToString((float)iover);
	for (int i = 0; i < 10; i++)
	{
		if (strcmp((joboptions["sampling"].getString()).c_str(), job_sampling_options[i]) == 0)
		{
			// The sampling given in the GUI will be the oversampled one!
			command += " --healpix_order " + floatToString((float)i + 1 - iover);
			break;
		}
	}
	// Offset range
	command += " --offset_range " + joboptions["offset_range"].getString();
	// The sampling given in the GUI will be the oversampled one!
	command += " --offset_step " + floatToString(joboptions["offset_step"].getNumber() * pow(2., iover));

	// Running stuff
	command += " --j " + joboptions["nr_threads"].getString();

	// GPU-stuff
	if (joboptions["use_gpu"].getBoolean())
	{
		command += " --gpu \"" + joboptions["gpu_ids"].getString() +"\"";
	}

	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);


	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

void RelionJob::initialiseClass3DJob()
{

	hidden_name = ".gui_class3d";

	joboptions["fn_img"] = JobOption("Input images STAR file:", NODE_PART_DATA, "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");
	joboptions["fn_cont"] = JobOption("Continue from here: ", std::string(""), "STAR Files (*_optimiser.star)", "CURRENT_ODIR", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");
	joboptions["fn_ref"] = JobOption("Reference map:", NODE_3DREF, "", "Image Files (*.{spi,vol,mrc})", "A 3D map in MRC/Spider format. \
	Make sure this map has the same dimensions and the same pixel size as your input images.");
	joboptions["fn_mask"] = JobOption("Reference mask (optional):", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "\
If no mask is provided, a soft spherical mask based on the particle diameter will be used.\n\
\n\
Otherwise, provide a Spider/mrc map containing a (soft) mask with the same \
dimensions as the reference(s), and values between 0 and 1, with 1 being 100% protein and 0 being 100% solvent. \
The reconstructed reference map will be multiplied by this mask.\n\
\n\
In some cases, for example for non-empty icosahedral viruses, it is also useful to use a second mask. For all white (value 1) pixels in this second mask \
the corresponding pixels in the reconstructed map are set to the average value of these pixels. \
Thereby, for example, the higher density inside the virion may be set to a constant. \
Note that this second mask should have one-values inside the virion and zero-values in the capsid and the solvent areas. \
To use a second mask, use the additional option --solvent_mask2, which may given in the Additional arguments line (in the Running tab).");

	joboptions["ref_correct_greyscale"] = JobOption("Ref. map is on absolute greyscale?", false, "Probabilities are calculated based on a Gaussian noise model, \
which contains a squared difference term between the reference and the experimental image. This has a consequence that the \
reference needs to be on the same absolute intensity grey-scale as the experimental images. \
RELION and XMIPP reconstruct maps at their absolute intensity grey-scale. \
Other packages may perform internal normalisations of the reference density, which will result in incorrect grey-scales. \
Therefore: if the map was reconstructed in RELION or in XMIPP, set this option to Yes, otherwise set it to No. \
If set to No, RELION will use a (grey-scale invariant) cross-correlation criterion in the first iteration, \
and prior to the second iteration the map will be filtered again using the initial low-pass filter. \
This procedure is relatively quick and typically does not negatively affect the outcome of the subsequent MAP refinement. \
Therefore, if in doubt it is recommended to set this option to No.");
	joboptions["ini_high"] = JobOption("Initial low-pass filter (A):", 60, 0, 200, 5, "It is recommended to strongly low-pass filter your initial reference map. \
If it has not yet been low-pass filtered, it may be done internally using this option. \
If set to 0, no low-pass filter will be applied to the initial reference(s).");
	joboptions["sym_name"] = JobOption("Symmetry:", std::string("C1"), "If the molecule is asymmetric, \
set Symmetry group to C1. Note their are multiple possibilities for icosahedral symmetry: \n \
* I1: No-Crowther 222 (standard in Heymann, Chagoyen & Belnap, JSB, 151 (2005) 196–207) \n \
* I2: Crowther 222 \n \
* I3: 52-setting (as used in SPIDER?)\n \
* I4: A different 52 setting \n \
The command 'relion_refine --sym D2 --print_symmetry_ops' prints a list of all symmetry operators for symmetry group D2. \
RELION uses XMIPP's libraries for symmetry operations. \
Therefore, look at the XMIPP Wiki for more details:  http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/WebHome?topic=Symmetry");

	joboptions["do_ctf_correction"] = JobOption("Do CTF-correction?", true, "If set to Yes, CTFs will be corrected inside the MAP refinement. \
The resulting algorithm intrinsically implements the optimal linear, or Wiener filter. \
Note that CTF parameters for all images need to be given in the input STAR file. \
The command 'relion_refine --print_metadata_labels' will print a list of all possible metadata labels for that STAR file. \
See the RELION Wiki for more details.\n\n Also make sure that the correct pixel size (in Angstrom) is given above!)");
	joboptions["ctf_corrected_ref"] = JobOption("Has reference been CTF-corrected?", false, "Set this option to Yes if the reference map \
represents density that is unaffected by CTF phases and amplitudes, e.g. it was created using CTF correction (Wiener filtering) inside RELION or from a PDB. \n\n\
If set to No, then in the first iteration, the Fourier transforms of the reference projections are not multiplied by the CTFs.");
	joboptions["ctf_phase_flipped"] = JobOption("Have data been phase-flipped?", false, "Set this to Yes if the images have been \
ctf-phase corrected during the pre-processing steps. \
Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION, \
as this can be done inside the internal CTF-correction. \
However, if the phases have been flipped, you should tell the program about it by setting this option to Yes.");
	joboptions["ctf_intact_first_peak"] = JobOption("Ignore CTFs until first peak?", false, "If set to Yes, then CTF-amplitude correction will \
only be performed from the first peak of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution. \
Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results. \
Therefore, this option is not generally recommended: try increasing amplitude contrast (in your input STAR file) first!");

	joboptions["nr_classes"] = JobOption("Number of classes:", 1, 1, 50, 1, "The number of classes (K) for a multi-reference refinement. \
These classes will be made in an unsupervised manner from a single reference by division of the data into random subsets during the first iteration.");
	joboptions["tau_fudge"] = JobOption("Regularisation parameter T:", 4 , 0.1, 10, 0.1, "Bayes law strictly determines the relative weight between \
the contribution of the experimental data and the prior. However, in practice one may need to adjust this weight to put slightly more weight on \
the experimental data to allow optimal results. Values greater than 1 for this regularisation parameter (T in the JMB2011 paper) put more \
weight on the experimental data. Values around 2-4 have been observed to be useful for 3D refinements, values of 1-2 for 2D refinements. \
Too small values yield too-low resolution structures; too high values result in over-estimated resolutions, mostly notable by the apparition of high-frequency noise in the references.");
	joboptions["nr_iter"] = JobOption("Number of iterations:", 25, 1, 50, 1, "Number of iterations to be performed. \
Note that the current implementation of 2D class averaging and 3D classification does NOT comprise a convergence criterium. \
Therefore, the calculations will need to be stopped by the user if further iterations do not yield improvements in resolution or classes. \n\n \
Also note that upon restarting, the iteration number continues to be increased, starting from the final iteration in the previous run. \
The number given here is the TOTAL number of iterations. For example, if 10 iterations have been performed previously and one restarts to perform \
an additional 5 iterations (for example with a finer angular sampling), then the number given here should be 10+5=15.");
	joboptions["do_subsets"] = JobOption("Use subsets for initial updates?", false, "If set to True, multiple maximization updates (as many as defined by the 'Number of subset updates') will be performed during the first iteration(s): each time after the number of particles in a subset has been processed. \
By using subsets with much fewer particles than the entire data set, the initial updates will be much faster, while the very low resolution class averages will not be notably worse than with the entire data set. \n\n \
This will greatly speed up classifications with very many (hundreds of thousands or more) particles. A useful subset size is probably in the order of ten thousand particles. If the data set only comprises (tens of) thousands of particles, this option may be less useful.");
	joboptions["subset_size"] = JobOption("Initial subset size:", 10000, 1000, 50000, 1000, "Number of individual particles after which one will perform a maximization update in the first iteration(s). \
A useful subset size is probably in the order of ten thousand particles.");
	joboptions["max_subsets"] = JobOption("Number of subset updates:", 3, 1, 10, 1, "This option is only used when a positive number is given for the 'Initial subset size'. In that case, in the first iteration, maximization updates are performed over a smaller subset of the particles to speed up calculations.\
Useful values are probably in the range of 2-5 subset updates. Using more might speed up further, but with the risk of affecting the results. If the number of subsets times the subset size is larger than the number of particles in the data set, then more than 1 iteration will be split into subsets.");
	joboptions["particle_diameter"] = JobOption("Mask diameter (A):", 200, 0, 1000, 10, "The experimental images will be masked with a soft \
circular mask with this diameter. Make sure this radius is not set too small because that may mask away part of the signal! \
If set to a value larger than the image size no masking will be performed.\n\n\
The same diameter will also be used for a spherical mask of the reference structures if no user-provided mask is specified.");
	joboptions["do_zero_mask"] = JobOption("Mask individual particles with zeros?", true, "If set to Yes, then in the individual particles, \
the area outside a circle with the radius of the particle will be set to zeros prior to taking the Fourier transform. \
This will remove noise and therefore increase sensitivity in the alignment and classification. However, it will also introduce correlations \
between the Fourier components that are not modelled. When set to No, then the solvent area is filled with random noise, which prevents introducing correlations.\
High-resolution refinements (e.g. ribosomes or other large complexes in 3D auto-refine) tend to work better when filling the solvent area with random noise (i.e. setting this option to No), refinements of smaller complexes and most classifications go better when using zeros (i.e. setting this option to Yes).");
	joboptions["highres_limit"] = JobOption("Limit resolution E-step to (A): ", -1, -1, 20, 1, "If set to a positive number, then the expectation step (i.e. the alignment) will be done only including the Fourier components up to this resolution (in Angstroms). \
This is useful to prevent overfitting, as the classification runs in RELION are not to be guaranteed to be 100% overfitting-free (unlike the 3D auto-refine with its gold-standard FSC). In particular for very difficult data sets, e.g. of very small or featureless particles, this has been shown to give much better class averages. \
In such cases, values in the range of 7-12 Angstroms have proven useful.");

	joboptions["dont_skip_align"] = JobOption("Perform image alignment?", true, "If set to No, then rather than \
performing both alignment and classification, only classification will be performed. This allows the use of very focused masks.\
This requires that the optimal orientations of all particles are already stored in the input STAR file. ");
	joboptions["sampling"] = JobOption("Angular sampling interval:", RADIO_SAMPLING, 2, "There are only a few discrete \
angular samplings possible because we use the HealPix library to generate the sampling of the first two Euler angles on the sphere. \
The samplings are approximate numbers and vary slightly over the sphere.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");
	joboptions["offset_range"] = JobOption("Offset search range (pix):", 5, 0, 30, 1, "Probabilities will be calculated only for translations \
in a circle with this radius (in pixels). The center of this circle changes at every iteration and is placed at the optimal translation \
for each image in the previous iteration.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");
	joboptions["offset_step"] = JobOption("Offset search step (pix):", 1, 0.1, 5, 0.1, "Translations will be sampled with this step-size (in pixels).\
Translational sampling is also done using the adaptive approach. \
Therefore, if adaptive=1, the translations will first be evaluated on a 2x coarser grid.\n\n \
If auto-sampling is used, this will be the value for the first iteration(s) only, and the sampling rate will be increased automatically after that.");
	joboptions["do_local_ang_searches"] = JobOption("Perform local angular searches?", false, "If set to Yes, then rather than \
performing exhaustive angular searches, local searches within the range given below will be performed. \
A prior Gaussian distribution centered at the optimal orientation in the previous iteration and \
with a stddev of 1/3 of the range given below will be enforced.");
	joboptions["sigma_angles"] = JobOption("Local angular search range:", 5., 0, 15, 0.1, "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.");

	joboptions["do_helix"] = JobOption("Do helical reconstruction?", false, "If set to Yes, then perform 3D helical reconstruction.");
	joboptions["helical_tube_inner_diameter"] = JobOption("Tube diameter - inner (A):", std::string("-1"),"Inner and outer diameter (in Angstroms) of the reconstructed helix spanning across Z axis. \
Set the inner diameter to negative value if the helix is not hollow in the center. The outer diameter should be slightly larger than the actual width of helical tubes because it also decides the shape of 2D \
particle mask for each segment. If the psi priors of the extracted segments are not accurate enough due to high noise level or flexibility of the structure, then set the outer diameter to a large value.");
	joboptions["helical_tube_outer_diameter"] = JobOption("Tube diameter - outer (A):", std::string("-1"),"Inner and outer diameter (in Angstroms) of the reconstructed helix spanning across Z axis. \
Set the inner diameter to negative value if the helix is not hollow in the center. The outer diameter should be slightly larger than the actual width of helical tubes because it also decides the shape of 2D \
particle mask for each segment. If the psi priors of the extracted segments are not accurate enough due to high noise level or flexibility of the structure, then set the outer diameter to a large value.");
	joboptions["range_tilt"] = JobOption("Angular search range - tilt (deg):", std::string("15"), "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.\n\nThese ranges will only be applied to the \
tilt and psi angles in the first few iterations (global searches for orientations) in 3D helical reconstruction. \
Values of 9 or 15 degrees are commonly used. Higher values are recommended for more flexible structures and more memory and computation time will be used. \
A range of 15 degrees means sigma = 5 degrees.\n\nThese options will be invalid if you choose to perform local angular searches or not to perform image alignment on 'Sampling' tab.");
	joboptions["range_psi"] = JobOption("Angular search range - psi (deg):", std::string("10"), "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.\n\nThese ranges will only be applied to the \
tilt and psi angles in the first few iterations (global searches for orientations) in 3D helical reconstruction. \
Values of 9 or 15 degrees are commonly used. Higher values are recommended for more flexible structures and more memory and computation time will be used. \
A range of 15 degrees means sigma = 5 degrees.\n\nThese options will be invalid if you choose to perform local angular searches or not to perform image alignment on 'Sampling' tab.");
	joboptions["do_apply_helical_symmetry"] = JobOption("Apply helical symmetry?", true, "If set to Yes, helical symmetry will be applied in every iteration. Set to No if you have just started a project, helical symmetry is unknown or not yet estimated.");
	joboptions["helical_nr_asu"] = JobOption("Number of asymmetrical units:", 1, 1, 100, 1, "Number of helical asymmetrical units in each segment box. If the inter-box distance (set in segment picking step) \
is 100 Angstroms and the estimated helical rise is ~20 Angstroms, then set this value to 100 / 20 = 5 (nearest integer). This integer should not be less than 1. The correct value is essential in measuring the \
signal to noise ratio in helical reconstruction.");
	joboptions["helical_twist_initial"] =  JobOption("Initial helical twist (deg):", std::string("0"),"Initial helical symmetry. Set helical twist (in degrees) to positive value if it is a right-handed helix. \
Helical rise is a positive value in Angstroms. If local searches of helical symmetry are planned, initial values of helical twist and rise should be within their respective ranges.");
	joboptions["helical_rise_initial"] = JobOption("Initial helical rise (A):", std::string("0"), "Initial helical symmetry. Set helical twist (in degrees) to positive value if it is a right-handed helix. \
Helical rise is a positive value in Angstroms. If local searches of helical symmetry are planned, initial values of helical twist and rise should be within their respective ranges.");
	joboptions["helical_z_percentage"] = JobOption("Central Z length (%):", 30., 5., 80., 1., "Reconstructed helix suffers from inaccuracies of orientation searches. \
The central part of the box contains more reliable information compared to the top and bottom parts along Z axis, where Fourier artefacts are also present if the \
number of helical asymmetrical units is larger than 1. Therefore, information from the central part of the box is used for searching and imposing \
helical symmetry in real space. Set this value (%) to the central part length along Z axis divided by the box size. Values around 30% are commonly used.");
	joboptions["do_local_search_helical_symmetry"] = JobOption("Do local searches of symmetry?", false, "If set to Yes, then perform local searches of helical twist and rise within given ranges.");
	joboptions["helical_twist_min"] = JobOption("Helical twist search (deg) - Min:", std::string("0"), "Minimum, maximum and initial step for helical twist search. Set helical twist (in degrees) \
to positive value if it is a right-handed helix. Generally it is not necessary for the user to provide an initial step (less than 1 degree, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.");
	joboptions["helical_twist_max"] = JobOption("Helical twist search (deg) - Max:", std::string("0"), "Minimum, maximum and initial step for helical twist search. Set helical twist (in degrees) \
to positive value if it is a right-handed helix. Generally it is not necessary for the user to provide an initial step (less than 1 degree, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.");
	joboptions["helical_twist_inistep"] = JobOption("Helical twist search (deg) - Step:", std::string("0"), "Minimum, maximum and initial step for helical twist search. Set helical twist (in degrees) \
to positive value if it is a right-handed helix. Generally it is not necessary for the user to provide an initial step (less than 1 degree, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.");
	joboptions["helical_rise_min"] = JobOption("Helical rise search (A) - Min:", std::string("0"), "Minimum, maximum and initial step for helical rise search. Helical rise is a positive value in Angstroms. \
Generally it is not necessary for the user to provide an initial step (less than 1% the initial helical rise, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.");
	joboptions["helical_rise_max"] = JobOption("Helical rise search (A) - Max:", std::string("0"), "Minimum, maximum and initial step for helical rise search. Helical rise is a positive value in Angstroms. \
Generally it is not necessary for the user to provide an initial step (less than 1% the initial helical rise, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.");
	joboptions["helical_rise_inistep"] = JobOption("Helical rise search (A) - Step:", std::string("0"), "Minimum, maximum and initial step for helical rise search. Helical rise is a positive value in Angstroms. \
Generally it is not necessary for the user to provide an initial step (less than 1% the initial helical rise, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.");
	joboptions["helical_range_distance"] = JobOption("Range factor of local averaging:", -1., 1., 5., 0.1, "Local averaging of orientations and translations will be performed within a range of +/- this value * the box size. Polarities are also set to be the same for segments coming from the same tube during local refinement. \
Values of ~ 2.0 are recommended for flexible structures such as MAVS-CARD filaments, ParM, MamK, etc. This option might not improve the reconstructions of helices formed from curled 2D lattices (TMV and VipA/VipB). Set to negative to disable this option.");

	joboptions["do_parallel_discio"] = JobOption("Use parallel disc I/O?", true, "If set to Yes, all MPI slaves will read their own images from disc. \
Otherwise, only the master will read images and send them through the network to the slaves. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many slaves reading in parallel.");
	joboptions["nr_pool"] = JobOption("Number of pooled particles:", 3, 1, 16, 1, "Particles are processed in individual batches by MPI slaves. During each batch, a stack of particle images is only opened and closed once to improve disk access times. \
All particle images of a single batch are read into memory together. The size of these batches is at least one particle per thread used. The nr_pooled_particles parameter controls how many particles are read together for each thread. If it is set to 3 and one uses 8 threads, batches of 3x8=24 particles will be read together. \
This may improve performance on systems where disk access, and particularly metadata handling of disk access, is a problem. It has a modest cost of increased RAM usage.");
	joboptions["do_preread_images"] = JobOption("Pre-read all particles into RAM?", false, "If set to Yes, all particle images will be read into computer memory, which will greatly speed up calculations on systems with slow disk access. However, one should of course be careful with the amount of RAM available. \
Because particles are read in float-precision, it will take ( N * box_size * box_size * 4 / (1024 * 1024 * 1024) ) Giga-bytes to read N particles into RAM. For 100 thousand 200x200 images, that becomes 15Gb, or 60 Gb for the same number of 400x400 particles. \
Remember that running a single MPI slave on each node that runs as many threads as available cores will have access to all available RAM. \n \n If parallel disc I/O is set to No, then only the master reads all particles into RAM and sends those particles through the network to the MPI slaves during the refinement iterations.");
	joboptions["scratch_dir"] = JobOption("Copy particles to scratch directory:", std::string(""), "If a directory is provided here, then the job will create a sub-directory in it called relion_volatile. If that relion_volatile directory already exists, it will be wiped. Then, the program will copy all input particles into a large stack inside the relion_volatile subdirectory. \
Provided this directory is on a fast local drive (e.g. an SSD drive), processing in all the iterations will be faster. If the job finishes correctly, the relion_volatile directory will be wiped. If the job crashes, you may want to remove it yourself.");
	joboptions["do_combine_thru_disc"] = JobOption("Combine iterations through disc?", false, "If set to Yes, at the end of every iteration all MPI slaves will write out a large file with their accumulated results. The MPI master will read in all these files, combine them all, and write out a new file with the combined results. \
All MPI salves will then read in the combined results. This reduces heavy load on the network, but increases load on the disc I/O. \
This will affect the time it takes between the progress-bar in the expectation step reaching its end (the mouse gets to the cheese) and the start of the ensuing maximisation step. It will depend on your system setup which is most efficient.");

	joboptions["use_gpu"] = JobOption("Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','.  For example: '0,0:1,1:0,0:1,1'");

}

bool RelionJob::getCommandsClass3DJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "Class3D", job_counter);
	std::string command;

	if (joboptions["nr_mpi"].getNumber() > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";

    FileName fn_run = "run";
	if (is_continue)
    {
		if (joboptions["fn_cont"].getString() == "")
		{
			error_message = "ERROR: empty field for continuation STAR file...";
			return false;
		}
		int pos_it = joboptions["fn_cont"].getString().rfind("_it");
		int pos_op = joboptions["fn_cont"].getString().rfind("_optimiser");
		if (pos_it < 0 || pos_op < 0)
			std::cerr << "Warning: invalid optimiser.star filename provided for continuation run: " << joboptions["fn_cont"].getString() << std::endl;
		int it = (int)textToFloat((joboptions["fn_cont"].getString().substr(pos_it+3, 6)).c_str());
		fn_run += "_ct" + floatToString(it);;
		command += " --continue " + joboptions["fn_cont"].getString();
    }

    command += " --o " + outputname + fn_run;
	outputNodes = getOutputNodesRefine(outputname + fn_run, (int)joboptions["nr_iter"].getNumber(), (int)joboptions["nr_classes"].getNumber(), 3, 1);

	if (!is_continue)
	{
		if (joboptions["fn_img"].getString() == "")
		{
			error_message = "ERROR: empty field for input STAR file...";
			return false;
		}
		command += " --i " + joboptions["fn_img"].getString();
		Node node(joboptions["fn_img"].getString(), joboptions["fn_img"].node_type);
		inputNodes.push_back(node);

		if (joboptions["fn_ref"].getString() == "")
		{
			error_message = "ERROR: empty field for reference. Type None for de-novo subtomogram averaging, provide reference for single-particle analysis.";
			return false;
		}

		if (joboptions["fn_ref"].getString() != "None")
		{
			command += " --ref " + joboptions["fn_ref"].getString();
			Node node(joboptions["fn_ref"].getString(), joboptions["fn_ref"].node_type);
			inputNodes.push_back(node);

			if (!joboptions["ref_correct_greyscale"].getBoolean()) // dont do firstiter_cc when giving None
				command += " --firstiter_cc";
		}

		if (joboptions["ini_high"].getNumber() > 0.)
			command += " --ini_high " + joboptions["ini_high"].getString();
	}

	// Always do compute stuff
	if (!joboptions["do_combine_thru_disc"].getBoolean())
		command += " --dont_combine_weights_via_disc";
	if (!joboptions["do_parallel_discio"].getBoolean())
		command += " --no_parallel_disc_io";
	if (joboptions["do_preread_images"].getBoolean())
            command += " --preread_images " ;
	else if (joboptions["scratch_dir"].getString() != "")
            command += " --scratch_dir " +  joboptions["scratch_dir"].getString();
        command += " --pool " + joboptions["nr_pool"].getString();

	// CTF stuff
	if (!is_continue)
	{

		if (joboptions["do_ctf_correction"].getBoolean())
		{
			command += " --ctf";
			if (joboptions["ctf_corrected_ref"].getBoolean())
				command += " --ctf_corrected_ref";
			if (joboptions["ctf_phase_flipped"].getBoolean())
				command += " --ctf_phase_flipped";
			if (joboptions["ctf_intact_first_peak"].getBoolean())
				command += " --ctf_intact_first_peak";
		}
	}

	// Optimisation
	command += " --iter " + joboptions["nr_iter"].getString();
	if (joboptions["do_subsets"].getBoolean())
	{
		command += " --write_subsets 1 --subset_size " + joboptions["subset_size"].getString();
		command += " --max_subsets " + joboptions["max_subsets"].getString();
	}
	command += " --tau2_fudge " + joboptions["tau_fudge"].getString();
    command += " --particle_diameter " + joboptions["particle_diameter"].getString();
	if (!is_continue)
	{
		command += " --K " + joboptions["nr_classes"].getString();
		// Always flatten the solvent
		command += " --flatten_solvent";
		if (joboptions["do_zero_mask"].getBoolean())
			command += " --zero_mask";
		if (joboptions["highres_limit"].getNumber() > 0)
			command += " --strict_highres_exp " + joboptions["highres_limit"].getString();
	}
	if (joboptions["fn_mask"].getString().length() > 0)
	{
		command += " --solvent_mask " + joboptions["fn_mask"].getString();
		Node node(joboptions["fn_mask"].getString(), joboptions["fn_mask"].node_type);
		inputNodes.push_back(node);
	}

	// Sampling
	if (!joboptions["dont_skip_align"].getBoolean())
	{
		command += " --skip_align ";
	}
	else
	{
		int iover = 1;
		command += " --oversampling " + floatToString((float)iover);
		for (int i = 0; i < 10; i++)
		{
			if (strcmp((joboptions["sampling"].getString()).c_str(), job_sampling_options[i]) == 0)
			{
				// The sampling given in the GUI will be the oversampled one!
				command += " --healpix_order " + floatToString((float)i + 1 - iover);
				break;
			}
		}
		// Manually input local angular searches
		if (joboptions["do_local_ang_searches"].getBoolean())
			command += " --sigma_ang " + floatToString(joboptions["sigma_angles"].getNumber() / 3.);

		// Offset range
		command += " --offset_range " + joboptions["offset_range"].getString();
		// The sampling given in the GUI will be the oversampled one!
		command += " --offset_step " +  floatToString(joboptions["offset_step"].getNumber() * pow(2., iover));
	}

	// Provide symmetry, and always do norm and scale correction
	if (!is_continue)
	{
		command += " --sym " + joboptions["sym_name"].getString();
		command += " --norm --scale ";
	}

	if ( (!is_continue) && (joboptions["do_helix"].getBoolean()) )
	{
		command += " --helix";
		if (textToFloat(joboptions["helical_tube_inner_diameter"].getString()) > 0.)
			command += " --helical_inner_diameter " + joboptions["helical_tube_inner_diameter"].getString();
		command += " --helical_outer_diameter " + joboptions["helical_tube_outer_diameter"].getString();
		if (joboptions["do_apply_helical_symmetry"].getBoolean())
		{
			command += " --helical_nr_asu " + joboptions["helical_nr_asu"].getString();
			command += " --helical_twist_initial " + joboptions["helical_twist_initial"].getString();
			command += " --helical_rise_initial " + joboptions["helical_rise_initial"].getString();
			command += " --helical_z_percentage " + floatToString(joboptions["helical_z_percentage"].getNumber() / 100.);
			if (joboptions["do_local_search_helical_symmetry"].getBoolean())
			{
				command += " --helical_symmetry_search";
				command += " --helical_twist_min " + joboptions["helical_twist_min"].getString();
				command += " --helical_twist_max " + joboptions["helical_twist_max"].getString();
				if (textToFloat(joboptions["helical_twist_inistep"].getString()) > 0.)
					command += " --helical_twist_inistep " + joboptions["helical_twist_inistep"].getString();
				command += " --helical_rise_min " + joboptions["helical_rise_min"].getString();
				command += " --helical_rise_max " + joboptions["helical_rise_max"].getString();
				if (textToFloat(joboptions["helical_rise_inistep"].getString()) > 0.)
					command += " --helical_rise_inistep " + joboptions["helical_rise_inistep"].getString();
			}
		}
		else
			command += " --ignore_helical_symmetry";
		if ( (joboptions["dont_skip_align"].getBoolean()) && (!joboptions["do_local_ang_searches"].getBoolean()) )
		{
			RFLOAT val = textToFloat(joboptions["range_tilt"].getString());
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_tilt " + floatToString(val / 3.);
			val = textToFloat(joboptions["range_psi"].getString());
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_psi " + floatToString(val / 3.);
			if (joboptions["helical_range_distance"].getNumber() > 0.)
				command += " --helical_sigma_distance " + floatToString(joboptions["helical_range_distance"].getNumber() / 3.);
		}
	}

	// Running stuff
	command += " --j " + joboptions["nr_threads"].getString();

	// GPU-stuff
	if (joboptions["use_gpu"].getBoolean())
	{
		command += " --gpu \"" + joboptions["gpu_ids"].getString() + "\"";
	}

	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

void RelionJob::initialiseAutorefineJob()
{
	type = PROC_3DAUTO;

	hidden_name = ".gui_auto3d";

	joboptions["fn_img"] = JobOption("Input images STAR file:", NODE_PART_DATA, "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");
	joboptions["fn_cont"] = JobOption("Continue from here: ", std::string(""), "STAR Files (*_optimiser.star)", "CURRENT_ODIR", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");
	joboptions["fn_ref"] = JobOption("Reference map:", NODE_3DREF, "", "Image Files (*.{spi,vol,mrc})", "A 3D map in MRC/Spider format. \
	Make sure this map has the same dimensions and the same pixel size as your input images.");
	joboptions["fn_mask"] = JobOption("Reference mask (optional):", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "\
If no mask is provided, a soft spherical mask based on the particle diameter will be used.\n\
\n\
Otherwise, provide a Spider/mrc map containing a (soft) mask with the same \
dimensions as the reference(s), and values between 0 and 1, with 1 being 100% protein and 0 being 100% solvent. \
The reconstructed reference map will be multiplied by this mask.\n\
\n\
In some cases, for example for non-empty icosahedral viruses, it is also useful to use a second mask. For all white (value 1) pixels in this second mask \
the corresponding pixels in the reconstructed map are set to the average value of these pixels. \
Thereby, for example, the higher density inside the virion may be set to a constant. \
Note that this second mask should have one-values inside the virion and zero-values in the capsid and the solvent areas. \
To use a second mask, use the additional option --solvent_mask2, which may given in the Additional arguments line (in the Running tab).");

	joboptions["ref_correct_greyscale"] = JobOption("Ref. map is on absolute greyscale?", false, "Probabilities are calculated based on a Gaussian noise model, \
which contains a squared difference term between the reference and the experimental image. This has a consequence that the \
reference needs to be on the same absolute intensity grey-scale as the experimental images. \
RELION and XMIPP reconstruct maps at their absolute intensity grey-scale. \
Other packages may perform internal normalisations of the reference density, which will result in incorrect grey-scales. \
Therefore: if the map was reconstructed in RELION or in XMIPP, set this option to Yes, otherwise set it to No. \
If set to No, RELION will use a (grey-scale invariant) cross-correlation criterion in the first iteration, \
and prior to the second iteration the map will be filtered again using the initial low-pass filter. \
This procedure is relatively quick and typically does not negatively affect the outcome of the subsequent MAP refinement. \
Therefore, if in doubt it is recommended to set this option to No.");
	joboptions["ini_high"] = JobOption("Initial low-pass filter (A):", 60, 0, 200, 5, "It is recommended to strongly low-pass filter your initial reference map. \
If it has not yet been low-pass filtered, it may be done internally using this option. \
If set to 0, no low-pass filter will be applied to the initial reference(s).");
	joboptions["sym_name"] = JobOption("Symmetry:", std::string("C1"), "If the molecule is asymmetric, \
set Symmetry group to C1. Note their are multiple possibilities for icosahedral symmetry: \n \
* I1: No-Crowther 222 (standard in Heymann, Chagoyen & Belnap, JSB, 151 (2005) 196–207) \n \
* I2: Crowther 222 \n \
* I3: 52-setting (as used in SPIDER?)\n \
* I4: A different 52 setting \n \
The command 'relion_refine --sym D2 --print_symmetry_ops' prints a list of all symmetry operators for symmetry group D2. \
RELION uses XMIPP's libraries for symmetry operations. \
Therefore, look at the XMIPP Wiki for more details:  http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/WebHome?topic=Symmetry");

	joboptions["do_ctf_correction"] = JobOption("Do CTF-correction?", true, "If set to Yes, CTFs will be applied to the projections of the map. This requires that CTF information is present in the input STAR file.");
	joboptions["ctf_corrected_ref"] = JobOption("Has reference been CTF-corrected?", false, "Set this option to Yes if the reference map \
represents density that is unaffected by CTF phases and amplitudes, e.g. it was created using CTF correction (Wiener filtering) inside RELION or from a PDB. \n\n\
If set to No, then in the first iteration, the Fourier transforms of the reference projections are not multiplied by the CTFs.");
	joboptions["ctf_phase_flipped"] = JobOption("Have data been phase-flipped?", false, "Set this to Yes if the images have been \
ctf-phase corrected during the pre-processing steps. \
Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION, \
as this can be done inside the internal CTF-correction. \
However, if the phases have been flipped, you should tell the program about it by setting this option to Yes.");
	joboptions["ctf_intact_first_peak"] = JobOption("Ignore CTFs until first peak?", false, "If set to Yes, then CTF-amplitude correction will \
only be performed from the first peak of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution. \
Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results. \
Therefore, this option is not generally recommended: try increasing amplitude contrast (in your input STAR file) first!");

	joboptions["particle_diameter"] = JobOption("Mask diameter (A):", 200, 0, 1000, 10, "The experimental images will be masked with a soft \
circular mask with this diameter. Make sure this radius is not set too small because that may mask away part of the signal! \
If set to a value larger than the image size no masking will be performed.\n\n\
The same diameter will also be used for a spherical mask of the reference structures if no user-provided mask is specified.");
	joboptions["do_zero_mask"] = JobOption("Mask individual particles with zeros?", true, "If set to Yes, then in the individual particles, \
the area outside a circle with the radius of the particle will be set to zeros prior to taking the Fourier transform. \
This will remove noise and therefore increase sensitivity in the alignment and classification. However, it will also introduce correlations \
between the Fourier components that are not modelled. When set to No, then the solvent area is filled with random noise, which prevents introducing correlations.\
High-resolution refinements (e.g. ribosomes or other large complexes in 3D auto-refine) tend to work better when filling the solvent area with random noise (i.e. setting this option to No), refinements of smaller complexes and most classifications go better when using zeros (i.e. setting this option to Yes).");
	joboptions["do_solvent_fsc"] = JobOption("Use solvent-flattened FSCs?", false, "If set to Yes, then instead of using unmasked maps to calculate the gold-standard FSCs during refinement, \
masked half-maps are used and a post-processing-like correction of the FSC curves (with phase-randomisation) is performed every iteration. This only works when a reference mask is provided on the I/O tab. \
This may yield higher-resolution maps, especially when the mask contains only a relatively small volume inside the box.");

	joboptions["sampling"] = JobOption("Initial angular sampling:", RADIO_SAMPLING, 2, "There are only a few discrete \
angular samplings possible because we use the HealPix library to generate the sampling of the first two Euler angles on the sphere. \
The samplings are approximate numbers and vary slightly over the sphere.\n\n \
Note that this will only be the value for the first few iteration(s): the sampling rate will be increased automatically after that.");
	joboptions["offset_range"] = JobOption("Initial offset range (pix):", 5, 0, 30, 1, "Probabilities will be calculated only for translations \
in a circle with this radius (in pixels). The center of this circle changes at every iteration and is placed at the optimal translation \
for each image in the previous iteration.\n\n \
Note that this will only be the value for the first few iteration(s): the sampling rate will be increased automatically after that.");
	joboptions["offset_step"] = JobOption("Initial offset step (pix):", 1, 0.1, 5, 0.1, "Translations will be sampled with this step-size (in pixels).\
Translational sampling is also done using the adaptive approach. \
Therefore, if adaptive=1, the translations will first be evaluated on a 2x coarser grid.\n\n \
Note that this will only be the value for the first few iteration(s): the sampling rate will be increased automatically after that.");
	joboptions["auto_local_sampling"] = JobOption("Local searches from auto-sampling:", RADIO_SAMPLING, 4, "In the automated procedure to \
increase the angular samplings, local angular searches of -6/+6 times the sampling rate will be used from this angular sampling rate onwards. For most \
lower-symmetric particles a value of 1.8 degrees will be sufficient. Perhaps icosahedral symmetries may benefit from a smaller value such as 0.9 degrees.");


	joboptions["do_helix"] = JobOption("Do helical reconstruction?", false, "If set to Yes, then perform 3D helical reconstruction.");
	joboptions["helical_tube_inner_diameter"] = JobOption("Tube diameter - inner (A):", std::string("-1"),"Inner and outer diameter (in Angstroms) of the reconstructed helix spanning across Z axis. \
Set the inner diameter to negative value if the helix is not hollow in the center. The outer diameter should be slightly larger than the actual width of helical tubes because it also decides the shape of 2D \
particle mask for each segment. If the psi priors of the extracted segments are not accurate enough due to high noise level or flexibility of the structure, then set the outer diameter to a large value.");
	joboptions["helical_tube_outer_diameter"] = JobOption("Tube diameter - outer (A):", std::string("-1"),"Inner and outer diameter (in Angstroms) of the reconstructed helix spanning across Z axis. \
Set the inner diameter to negative value if the helix is not hollow in the center. The outer diameter should be slightly larger than the actual width of helical tubes because it also decides the shape of 2D \
particle mask for each segment. If the psi priors of the extracted segments are not accurate enough due to high noise level or flexibility of the structure, then set the outer diameter to a large value.");
	joboptions["range_tilt"] = JobOption("Angular search range - tilt (deg):", std::string("15"), "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.\n\nThese ranges will only be applied to the \
tilt and psi angles in the first few iterations (global searches for orientations) in 3D helical reconstruction. \
Values of 9 or 15 degrees are commonly used. Higher values are recommended for more flexible structures and more memory and computation time will be used. \
A range of 15 degrees means sigma = 5 degrees.\n\nThese options will be invalid if you choose to perform local angular searches or not to perform image alignment on 'Sampling' tab.");
	joboptions["range_psi"] = JobOption("Angular search range - psi (deg):", std::string("10"), "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.\n\nThese ranges will only be applied to the \
tilt and psi angles in the first few iterations (global searches for orientations) in 3D helical reconstruction. \
Values of 9 or 15 degrees are commonly used. Higher values are recommended for more flexible structures and more memory and computation time will be used. \
A range of 15 degrees means sigma = 5 degrees.\n\nThese options will be invalid if you choose to perform local angular searches or not to perform image alignment on 'Sampling' tab.");
	joboptions["do_apply_helical_symmetry"] = JobOption("Apply helical symmetry?", true, "If set to Yes, helical symmetry will be applied in every iteration. Set to No if you have just started a project, helical symmetry is unknown or not yet estimated.");
	joboptions["helical_nr_asu"] = JobOption("Number of asymmetrical units:", 1, 1, 100, 1, "Number of helical asymmetrical units in each segment box. If the inter-box distance (set in segment picking step) \
is 100 Angstroms and the estimated helical rise is ~20 Angstroms, then set this value to 100 / 20 = 5 (nearest integer). This integer should not be less than 1. The correct value is essential in measuring the \
signal to noise ratio in helical reconstruction.");
	joboptions["helical_twist_initial"] =  JobOption("Initial helical twist (deg):", std::string("0"),"Initial helical symmetry. Set helical twist (in degrees) to positive value if it is a right-handed helix. \
Helical rise is a positive value in Angstroms. If local searches of helical symmetry are planned, initial values of helical twist and rise should be within their respective ranges.");
	joboptions["helical_rise_initial"] = JobOption("Initial helical rise (A):", std::string("0"), "Initial helical symmetry. Set helical twist (in degrees) to positive value if it is a right-handed helix. \
Helical rise is a positive value in Angstroms. If local searches of helical symmetry are planned, initial values of helical twist and rise should be within their respective ranges.");
	joboptions["helical_z_percentage"] = JobOption("Central Z length (%):", 30., 5., 80., 1., "Reconstructed helix suffers from inaccuracies of orientation searches. \
The central part of the box contains more reliable information compared to the top and bottom parts along Z axis, where Fourier artefacts are also present if the \
number of helical asymmetrical units is larger than 1. Therefore, information from the central part of the box is used for searching and imposing \
helical symmetry in real space. Set this value (%) to the central part length along Z axis divided by the box size. Values around 30% are commonly used.");
	joboptions["do_local_search_helical_symmetry"] = JobOption("Do local searches of symmetry?", false, "If set to Yes, then perform local searches of helical twist and rise within given ranges.");
	joboptions["helical_twist_min"] = JobOption("Helical twist search (deg) - Min:", std::string("0"), "Minimum, maximum and initial step for helical twist search. Set helical twist (in degrees) \
to positive value if it is a right-handed helix. Generally it is not necessary for the user to provide an initial step (less than 1 degree, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.");
	joboptions["helical_twist_max"] = JobOption("Helical twist search (deg) - Max:", std::string("0"), "Minimum, maximum and initial step for helical twist search. Set helical twist (in degrees) \
to positive value if it is a right-handed helix. Generally it is not necessary for the user to provide an initial step (less than 1 degree, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.");
	joboptions["helical_twist_inistep"] = JobOption("Helical twist search (deg) - Step:", std::string("0"), "Minimum, maximum and initial step for helical twist search. Set helical twist (in degrees) \
to positive value if it is a right-handed helix. Generally it is not necessary for the user to provide an initial step (less than 1 degree, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.");
	joboptions["helical_rise_min"] = JobOption("Helical rise search (A) - Min:", std::string("0"), "Minimum, maximum and initial step for helical rise search. Helical rise is a positive value in Angstroms. \
Generally it is not necessary for the user to provide an initial step (less than 1% the initial helical rise, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.");
	joboptions["helical_rise_max"] = JobOption("Helical rise search (A) - Max:", std::string("0"), "Minimum, maximum and initial step for helical rise search. Helical rise is a positive value in Angstroms. \
Generally it is not necessary for the user to provide an initial step (less than 1% the initial helical rise, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.");
	joboptions["helical_rise_inistep"] = JobOption("Helical rise search (A) - Step:", std::string("0"), "Minimum, maximum and initial step for helical rise search. Helical rise is a positive value in Angstroms. \
Generally it is not necessary for the user to provide an initial step (less than 1% the initial helical rise, 5~1000 samplings as default). But it needs to be set manually if the default value \
does not guarantee convergence. The program cannot find a reasonable symmetry if the true helical parameters fall out of the given ranges. Note that the final reconstruction can still converge if wrong helical and point group symmetry are provided.");
	joboptions["helical_range_distance"] = JobOption("Range factor of local averaging:", -1., 1., 5., 0.1, "Local averaging of orientations and translations will be performed within a range of +/- this value * the box size. Polarities are also set to be the same for segments coming from the same tube during local refinement. \
Values of ~ 2.0 are recommended for flexible structures such as MAVS-CARD filaments, ParM, MamK, etc. This option might not improve the reconstructions of helices formed from curled 2D lattices (TMV and VipA/VipB). Set to negative to disable this option.");


	joboptions["do_parallel_discio"] = JobOption("Use parallel disc I/O?", true, "If set to Yes, all MPI slaves will read their own images from disc. \
Otherwise, only the master will read images and send them through the network to the slaves. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many slaves reading in parallel.");
	joboptions["nr_pool"] = JobOption("Number of pooled particles:", 3, 1, 16, 1, "Particles are processed in individual batches by MPI slaves. During each batch, a stack of particle images is only opened and closed once to improve disk access times. \
All particle images of a single batch are read into memory together. The size of these batches is at least one particle per thread used. The nr_pooled_particles parameter controls how many particles are read together for each thread. If it is set to 3 and one uses 8 threads, batches of 3x8=24 particles will be read together. \
This may improve performance on systems where disk access, and particularly metadata handling of disk access, is a problem. It has a modest cost of increased RAM usage.");
	joboptions["do_preread_images"] = JobOption("Pre-read all particles into RAM?", false, "If set to Yes, all particle images will be read into computer memory, which will greatly speed up calculations on systems with slow disk access. However, one should of course be careful with the amount of RAM available. \
Because particles are read in float-precision, it will take ( N * box_size * box_size * 8 / (1024 * 1024 * 1024) ) Giga-bytes to read N particles into RAM. For 100 thousand 200x200 images, that becomes 15Gb, or 60 Gb for the same number of 400x400 particles. \
Remember that running a single MPI slave on each node that runs as many threads as available cores will have access to all available RAM. \n \n If parallel disc I/O is set to No, then only the master reads all particles into RAM and sends those particles through the network to the MPI slaves during the refinement iterations.");
	joboptions["scratch_dir"] = JobOption("Copy particles to scratch directory:", std::string(""), "If a directory is provided here, then the job will create a sub-directory in it called relion_volatile. If that relion_volatile directory already exists, it will be wiped. Then, the program will copy all input particles into a large stack inside the relion_volatile subdirectory. \
Provided this directory is on a fast local drive (e.g. an SSD drive), processing in all the iterations will be faster. If the job finishes correctly, the relion_volatile directory will be wiped. If the job crashes, you may want to remove it yourself.");
	joboptions["do_combine_thru_disc"] = JobOption("Combine iterations through disc?", false, "If set to Yes, at the end of every iteration all MPI slaves will write out a large file with their accumulated results. The MPI master will read in all these files, combine them all, and write out a new file with the combined results. \
All MPI salves will then read in the combined results. This reduces heavy load on the network, but increases load on the disc I/O. \
This will affect the time it takes between the progress-bar in the expectation step reaching its end (the mouse gets to the cheese) and the start of the ensuing maximisation step. It will depend on your system setup which is most efficient.");
	joboptions["use_gpu"] = JobOption("Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','.  For example: '0,0:1,1:0,0:1,1'");


}

bool RelionJob::getCommandsAutorefineJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "Refine3D", job_counter);
	std::string command;

	if (joboptions["nr_mpi"].getNumber() > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";

    FileName fn_run = "run";
	if (is_continue)
    {
		if (joboptions["fn_cont"].getString() == "")
		{
			error_message = "ERROR: empty field for continuation STAR file...";
			return false;
		}
		int pos_it = joboptions["fn_cont"].getString().rfind("_it");
		int pos_op = joboptions["fn_cont"].getString().rfind("_optimiser");
		if (pos_it < 0 || pos_op < 0)
			std::cerr << "Warning: invalid optimiser.star filename provided for continuation run: " << joboptions["fn_cont"].getString() << std::endl;
		int it = (int)textToFloat((joboptions["fn_cont"].getString().substr(pos_it+3, 6)).c_str());
		fn_run += "_ct" + floatToString(it);
		command += " --continue " + joboptions["fn_cont"].getString();
    }

    command += " --o " + outputname + fn_run;
	// TODO: add bodies!! (probably in next version)
	outputNodes = getOutputNodesRefine(outputname + fn_run, -1, 1, 3, 1, false, false); // false false means dont do movies

	if (!is_continue)
	{
		command += " --auto_refine --split_random_halves --i " + joboptions["fn_img"].getString();
		if (joboptions["fn_img"].getString() == "")
		{
			error_message = "ERROR: empty field for input STAR file...";
			return false;
		}
		Node node(joboptions["fn_img"].getString(), joboptions["fn_img"].node_type);
		inputNodes.push_back(node);
		if (joboptions["fn_ref"].getString() == "")
		{
			error_message = "ERROR: empty field for input reference...";
			return false;
		}
		if (joboptions["fn_ref"].getString() != "None")
		{
			command += " --ref " + joboptions["fn_ref"].getString();
			Node node(joboptions["fn_ref"].getString(), joboptions["fn_ref"].node_type);
			inputNodes.push_back(node);

			if (!joboptions["ref_correct_greyscale"].getBoolean())
				command += " --firstiter_cc";
		}
		if (joboptions["ini_high"].getNumber() > 0.)
			command += " --ini_high " + joboptions["ini_high"].getString();

	}

	// Always do compute stuff
	if (!joboptions["do_combine_thru_disc"].getBoolean())
		command += " --dont_combine_weights_via_disc";
	if (!joboptions["do_parallel_discio"].getBoolean())
		command += " --no_parallel_disc_io";
	if (joboptions["do_preread_images"].getBoolean())
		command += " --preread_images " ;
	else if (joboptions["scratch_dir"].getString() != "")
                command += " --scratch_dir " +  joboptions["scratch_dir"].getString();
	command += " --pool " + joboptions["nr_pool"].getString();

	// CTF stuff
	if (!is_continue)
	{

		if (joboptions["do_ctf_correction"].getBoolean())
		{
			command += " --ctf";
			if (joboptions["ctf_corrected_ref"].getBoolean())
				command += " --ctf_corrected_ref";
			if (joboptions["ctf_phase_flipped"].getBoolean())
				command += " --ctf_phase_flipped";
			if (joboptions["ctf_intact_first_peak"].getBoolean())
				command += " --ctf_intact_first_peak";
		}
	}

	// Optimisation
    command += " --particle_diameter " + joboptions["particle_diameter"].getString();
	if (!is_continue)
	{
		// Always flatten the solvent
		command += " --flatten_solvent";
		if (joboptions["do_zero_mask"].getBoolean())
			command += " --zero_mask";
	}
	if (joboptions["fn_mask"].getString().length() > 0)
	{
		command += " --solvent_mask " + joboptions["fn_mask"].getString();

		if (joboptions["do_solvent_fsc"].getBoolean())
			command += " --solvent_correct_fsc ";

		Node node(joboptions["fn_mask"].getString(), joboptions["fn_mask"].node_type);
		inputNodes.push_back(node);
	}

	if (!is_continue)
	{
		// Sampling
		int iover = 1;
		command += " --oversampling " + floatToString((float)iover);
		for (int i = 0; i < 10; i++)
		{
			if (strcmp((joboptions["sampling"].getString()).c_str(), job_sampling_options[i]) == 0)
			{
				// The sampling given in the GUI will be the oversampled one!
				command += " --healpix_order " + floatToString((float)i + 1 - iover);
				break;
			}
		}
		// Minimum sampling rate to perform local searches (may be changed upon continuation
		for (int i = 0; i < 10; i++)
		{
			if (strcmp((joboptions["auto_local_sampling"].getString()).c_str(), job_sampling_options[i]) == 0)
			{
				command += " --auto_local_healpix_order " + floatToString((float)i + 1 - iover);
				break;
			}
		}

		// Offset range
		command += " --offset_range " + joboptions["offset_range"].getString();
		// The sampling given in the GUI will be the oversampled one!
		command += " --offset_step " + floatToString(joboptions["offset_step"].getNumber() * pow(2., iover));

		command += " --sym " + joboptions["sym_name"].getString();
		// Always join low-res data, as some D&I point group refinements may fall into different hands!
		command += " --low_resol_join_halves 40";
		command += " --norm --scale ";

		// Helix
		if (joboptions["do_helix"].getBoolean())
		{
			command += " --helix";
			if (textToFloat(joboptions["helical_tube_inner_diameter"].getString()) > 0.)
				command += " --helical_inner_diameter " + joboptions["helical_tube_inner_diameter"].getString();
			command += " --helical_outer_diameter " + joboptions["helical_tube_outer_diameter"].getString();
			if (joboptions["do_apply_helical_symmetry"].getBoolean())
			{
				command += " --helical_nr_asu " + joboptions["helical_nr_asu"].getString();
				command += " --helical_twist_initial " + joboptions["helical_twist_initial"].getString();
				command += " --helical_rise_initial " + joboptions["helical_rise_initial"].getString();
				command += " --helical_z_percentage " + floatToString(joboptions["helical_z_percentage"].getNumber() / 100.);
				if (joboptions["do_local_search_helical_symmetry"].getBoolean())
				{
					command += " --helical_symmetry_search";
					command += " --helical_twist_min " + joboptions["helical_twist_min"].getString();
					command += " --helical_twist_max " + joboptions["helical_twist_max"].getString();
					if (textToFloat(joboptions["helical_twist_inistep"].getString()) > 0.)
						command += " --helical_twist_inistep " + joboptions["helical_twist_inistep"].getString();
					command += " --helical_rise_min " + joboptions["helical_rise_min"].getString();
					command += " --helical_rise_max " + joboptions["helical_rise_max"].getString();
					if (textToFloat(joboptions["helical_rise_inistep"].getString()) > 0.)
						command += " --helical_rise_inistep " + joboptions["helical_rise_inistep"].getString();
				}
			}
			else
				command += " --ignore_helical_symmetry";
			RFLOAT val = textToFloat(joboptions["range_tilt"].getString());
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_tilt " + floatToString(val / 3.);
			val = textToFloat(joboptions["range_psi"].getString());
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_psi " + floatToString(val / 3.);
			if (joboptions["helical_range_distance"].getNumber() > 0.)
				command += " --helical_sigma_distance " + floatToString(joboptions["helical_range_distance"].getNumber() / 3.);
		}
	}

	// Running stuff
	command += " --j " + joboptions["nr_threads"].getString();

	// GPU-stuff
	if (joboptions["use_gpu"].getBoolean())
	{
		command += " --gpu \"" + joboptions["gpu_ids"].getString() + "\"";
	}

	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);


}

void RelionJob::initialiseMovierefineJob()
{

	hidden_name = ".gui_movierefine";

	joboptions["fn_movie_star"] = JobOption("Input movies STAR file:", NODE_MOVIES, "", "STAR Files (*.{star})", "Select the STAR file with the input movie micrographs, either from an Import or a Motion correction job.");
    joboptions["movie_rootname"] = JobOption("Rootname of movies files:", std::string("movie"), "rootname to relate each movie to the single-frame averaged micropgraph. With a rootname of 'movie', the movie for mic001.mrc should be called mic001_movie.mrcs. If you've run the MOTIONCORR wrapper in RELION, then the correct rootname is 'movie'.");

	joboptions["fn_cont"] = JobOption("Continue 3D auto-refine from: ", std::string(""), "STAR Files (*_optimiser.star)", "Refine3D/", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run. \n \
Besides restarting jobs that were somehow stopped before convergence, also use the continue-option after the last iteration to do movie processing.");

	joboptions["join_nr_mics"] = JobOption("Process micrographs in batches of: ", 50, 10, 200, 10, "All movie-particles will be extracted from each micrograph separately, but the resulting STAR files will be joined together into batches of the size specified here. The movie-refinement will be performed on these batches. \
Very large batches cost more RAM, but the parallelisation in smaller batches is poorer. If a negative number is given, all micrographs are processed in one batch. Note that movie-refinement that INCLUDE rotational searches cannot be performed in batches!");

	joboptions["first_movie_frame"] = JobOption("First movie frame to extract: ", 1, 1, 20, 1, "Extract from this movie frame onwards. The first frame is number 1.");
	joboptions["last_movie_frame"] = JobOption("Last movie frame to extract: ", 0, 0, 64, 1, "Extract until this movie frame. Zero means: extract all frames in the movie. You may want to specify the last frame number though, as it will be useful to detect movies which accidentally have fewer frames.");
	joboptions["avg_movie_frames"] = JobOption("Average every so many frames: ", 1, 1, 8, 1, "Average every so many movie frames together upon the extraction. For example, 32-frame movies may be reduced to 16-frame movie-particles when provding a value of 2 here. This will reduce computational costs in movie-refinement and polishing, but too large values will affect the results. Default is a value of 1, so no averaging");
	joboptions["max_mpi_nodes"] = JobOption("Maximum number of MPI nodes: ", 8, 2, 24, 1, "The number of MPI nodes used by the relion_preprocess program will be limited to this value, regardless of the number of MPI nodes requested on the Running tab (which is also used for the refinement step). This is useful to protect the file system from too heavy disk I/O.");
	joboptions["extract_size"] = JobOption("Particle box size (pix):", 128, 64, 512, 8, "Size of the extracted particles (in pixels). Use the same as for the refined particles!");
	joboptions["do_invert"] = JobOption("Invert contrast?", true, "If set to Yes, the contrast in the particles will be inverted. Use the same as for the refined particles!");

	joboptions["do_rescale"] = JobOption("Rescale particles?", false, "If set to Yes, particles will be re-scaled. Note that the particle diameter below will be in the down-scaled images.");
	joboptions["rescale"] = JobOption("Re-scaled size (pixels): ", 128, 64, 512, 8, "The re-scaled value needs to be an even number. Use the same as for the refined particles! ");

	joboptions["do_norm"] = JobOption("Normalize movie-particles?", true, "If set to Yes, particles will be normalized in the way RELION prefers it. Use the same values as for the refined particles!");
	joboptions["bg_diameter"] = JobOption("Diameter background circle (pix): ", -1, -1, 600, 10, " Use the same as for the refined particles! Particles will be normalized to a mean value of zero and a standard-deviation of one for all pixels in the background area.\
The background area is defined as all pixels outside a circle with this given diameter in pixels (before rescaling). When specifying a negative value, a default value of 75% of the Particle box size will be used.");
	joboptions["white_dust"] = JobOption("Stddev for white dust removal: ", -1, -1, 10, 0.1, " Use the same as for the refined particles! Remove very white pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
	joboptions["black_dust"] = JobOption("Stddev for black dust removal: ", -1, -1, 10, 0.1, " Use the same as for the refined particles! Remove very black pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");

	joboptions["movie_runavg_window"] = JobOption("Running average window:", 5, 1, 15, 1, "The individual movie frames will be averaged using a running \
average window with the specified width. Use an odd number. The optimal value will depend on the SNR in the individual movie frames. For ribosomes, we used a value of 5, where \
each movie frame integrated approximately 1 electron per squared Angstrom.");
	joboptions["movie_sigma_offset"] = JobOption("Stddev on the translations (pix):", 1., 0.5, 10, 0.5, "A Gaussian prior with the specified standard deviation \
will be centered at the rotations determined for the corresponding particle where all movie-frames were averaged. For ribosomes, we used a value of 2 pixels");
	joboptions["do_alsorot_movies"] = JobOption("Also include rotational searches?", false, "If set to Yes, then running averages of the individual frames of recorded movies will also be aligned rotationally. \n \
If one wants to perform particle polishing, then rotational alignments of the movie frames is NOT necessary and will only take more computing time.");
	joboptions["movie_sigma_angles"] = JobOption("Stddev on the rotations (deg):", 1., 0.5, 10, 0.5, "A Gaussian prior with the specified standard deviation \
will be centered at the rotations determined for the corresponding particle where all movie-frames were averaged. For ribosomes, we used a value of 1 degree");


}

bool RelionJob::getCommandsMovierefineJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "MovieRefine", job_counter);
	std::string command;

	// A. First get the extract command
	if (joboptions["nr_mpi"].getNumber() > 1)
		command="`which relion_preprocess_mpi`";
	else
		command="`which relion_preprocess`";

	// Input
	if (joboptions["fn_movie_star"].getString() == "")
	{
		error_message = "ERROR: empty field for input STAR file...";
		return false;
	}
	command += " --i " + joboptions["fn_movie_star"].getString();
	Node node(joboptions["fn_movie_star"].getString(), joboptions["fn_movie_star"].node_type);
	inputNodes.push_back(node);

	// Get the data.star to be used for re-extraction from the optimiser name
	if (joboptions["fn_cont"].getString() == "")
	{
		error_message = "ERROR: empty field for continuation STAR file...";
		return false;
	}
	MetaDataTable MDopt;
	MDopt.read(joboptions["fn_cont"].getString(), "optimiser_general");
	FileName fn_data;
	MDopt.getValue(EMDL_OPTIMISER_DATA_STARFILE, fn_data);
	command += " --reextract_data_star " + fn_data;
	Node node2(joboptions["fn_cont"].getString(), NODE_OPTIMISER);
	inputNodes.push_back(node2);

	// Output files of the extraction: to be used in the second step: movie-refinement
	command += " --part_dir " + outputname;

	FileName fn_ostar = outputname + "particles_" + joboptions["movie_rootname"].getString() + ".star";
	FileName fn_olist = outputname + "micrographs_" + joboptions["movie_rootname"].getString() + "_list.star";
	command += " --list_star " + fn_olist;
	command += " --join_nr_mics " + joboptions["join_nr_mics"].getString();
	if (joboptions["join_nr_mics"].getNumber() <= 0)
	{
		// Only write out a STAR file with all movie-particles if we're NOT processing on a per-micrograph basis
		command += " --part_star " + fn_ostar;
	}

	// Extraction parameters
	command += " --extract --extract_movies";
	command += " --extract_size " + joboptions["extract_size"].getString();
	command += " --movie_name " + joboptions["movie_rootname"].getString();
	command += " --first_movie_frame " + joboptions["first_movie_frame"].getString();
	command += " --last_movie_frame " + joboptions["last_movie_frame"].getString();
	command += " --avg_movie_frames " + joboptions["avg_movie_frames"].getString();
	// Limit MPI nodes
	if (joboptions["nr_mpi"].getNumber() > ROUND(joboptions["max_mpi_nodes"].getNumber()))
		command += " --max_mpi_nodes " + joboptions["max_mpi_nodes"].getString();

	// Operate stuff
	// Get an integer number for the bg_radius
	RFLOAT bg_radius = (joboptions["bg_diameter"].getNumber() < 0.) ? 0.75 * joboptions["extract_size"].getNumber() : joboptions["bg_diameter"].getNumber();
	bg_radius /= 2.; // Go from diameter to radius
	if (joboptions["do_rescale"].getBoolean())
	{
		command += " --scale " + joboptions["rescale"].getString();
		bg_radius *= joboptions["rescale"].getNumber() / joboptions["extract_size"].getNumber();
	}
	if (joboptions["do_norm"].getBoolean())
	{
		// Get an integer number for the bg_radius
		bg_radius = (int)bg_radius;
		command += " --norm --bg_radius " + floatToString(bg_radius);
		command += " --white_dust " + joboptions["white_dust"].getString();
		command += " --black_dust " + joboptions["black_dust"].getString();
	}
	if (joboptions["do_invert"].getBoolean())
		command += " --invert_contrast ";

	if (is_continue)
		command += " --only_extract_unfinished ";

	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();
	commands.push_back(command);

	// B. Then get the actual movie-refinement command
	if (joboptions["nr_mpi"].getNumber() > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";

	// Make specific run-names for different variables of ravg, sigma_offset and sigma_angles,
	// This way you don't have to re-extract all movie particles each time you try a different parameter setting
	std::string runname = "run_ravg" + joboptions["movie_runavg_window"].getString() + "_off" + joboptions["movie_sigma_offset"].getString();
	if (joboptions["do_alsorot_movies"].getBoolean())
		runname += "_ang" + joboptions["movie_sigma_angles"].getString();

	command += " --o " + outputname + runname;
	outputNodes = getOutputNodesRefine(outputname + runname, -1, 1, 3, 1, true, joboptions["do_alsorot_movies"].getBoolean() );

	command += " --continue " + joboptions["fn_cont"].getString();

	if (joboptions["join_nr_mics"].getNumber() > 0)
	{
		if (joboptions["do_alsorot_movies"].getBoolean())
		{
			error_message = "You cannot process micrographs in batches and perform rotational searches!";
			return false;
		}

		command += " --process_movies_in_batches --realign_movie_frames " + fn_olist;

		if (is_continue)
			command += " --only_do_unfinished_movies ";
	}
	else
	{
		command += " --realign_movie_frames " + fn_ostar;
	}

	command += " --movie_frames_running_avg " + joboptions["movie_runavg_window"].getString();
	command += " --movie_name " + joboptions["movie_rootname"].getString();
	command += " --sigma_off " + joboptions["movie_sigma_offset"].getString();

	if (joboptions["do_alsorot_movies"].getBoolean())
	{
		command += " --sigma_ang " + joboptions["movie_sigma_angles"].getString();
	}
	else
	{
		command += " --skip_rotate --skip_maximize ";
	}

	// Running stuff
	command += " --j " + joboptions["nr_threads"].getString();

	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);


}

void RelionJob::initialisePolishJob()
{

	hidden_name = ".gui_polish";

	joboptions["fn_in"] = JobOption("Input STAR file with aligned movies:", NODE_MOVIE_DATA, "", "STAR files (*_data.star)",  "Provide the data.star file that was output by the movie-processing option in the auto-refine job.");
	joboptions["fn_mask"] = JobOption("Mask for the reconstructions", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "A continuous mask with values between 0 (solvent) and 1 (protein). You may provide the same map that was obtained in the post-processing of the corresponding auto-refine jobs before the movie processing.");

	joboptions["do_fit_movement"] = JobOption("Linear fit particle movements?", true, "If set to Yes, then the program will fit linear tracks (in X-Y and in time) through \
the estimated movie tracks in the input STAR file. For small particles (e.g. < 1MDa) this will lead to more robust beam-induced movement modelling. Because particles that are close to each other on a \
micrograph often move in similar directions, the estimated tracks from neighbouring particles may also be included in the fitting of each particle. Again, in particular for smaller particles \
this may improve the robustness of the fits.");
	joboptions["sigma_nb"] = JobOption("Stddev on particle distance (pix)", 100, 0, 1000, 50, "This value determines how much neighbouring particles contribute to the fit of the movements of each particle. \
This value is the standard deviation of a Gaussian on the inter-particle distance. Larger values mean that particles that are further away still contribute more. Particles beyond 3 standard deviations are excluded \
from the fit. Very large values will lead to all fitted tracks pointing in the same direction. A value of zero means that each particle is fitted independently.");

	joboptions["do_bfactor_weighting"] = JobOption("Perform B-factor weighting?", true, "If set to Yes, then the program will estimate a resolution-dependent weight for each movie frames by calculating independent half-reconstructions for each movie frame separately. \
Gold-standard FSCs between these are then converted into relative Guinier plots, through which straight lines are fitted. Linear fits are often suitable beyond 20A resolution. Small particles may not yield such high resolutions for the individual-frame reconstructions. \
Therefore, in some cases it may be better to skip this step. It is however recommended to always try and perform B-factor weighting, and to inspect the output bfactors.star and guinier.star files, as adequate weighting may significantly improve resolution in the final map.");
	joboptions["perframe_highres"] = JobOption("Highres-limit per-frame maps (A)", 6, 1, 25, 1, "To estimate the resolution and dose dependency of the radiation damage, the program \
will calculate reconstructions from all first movie frames, second movie frames, etc. These per-frame reconstructions will have lower resolution than the reconstruction from all-frames. \
To speed up the calculations (and reduce memory requirements), the per-frame reconstructions may be limited in resolution using this parameter. One can inspect the output STAR files of the per-frame reconstructions \
to check afterwards that this value was not chosen lower than the actual resolution of these reconstructions");
	joboptions["perframe_bfac_lowres"] = JobOption("Lowres-limit B-factor estimation (A)", 20 , 1, 40, 1, "This value describes the lowest resolution that is included in the B-factor estimation of the per-frame reconstructions. \
Because the power spectrum of per-frame reconstructions is compared to the power spectrum of the reconstruction from all frames, a much lower value than the 10A described in the Rosenthal and Henderson (2003) paper in JMB can be used. Probably a value around 20A is still OK.");
	joboptions["average_frame_bfactor"] = JobOption("Average frames B-factor estimation", 1 , 1, 7, 1, "B-factors for each movie frame will be estimated from reconstructions of all particles for that movie frame. Single-frame reconstructions sometimes give not enough signal to estimate reliable B-factors. \
This option allows one to calculate the B-factors from running averages of movie frames. The value specified should be an odd number. Calculating B-factors from multiple movie frames improves the SNR in the reconstructions, but limits the estimation of sudden changes in B-factors throughout the movie, for example in the first few frames when beam-induced movement is very rapid. \
Therefore, one should not use higher values than strictly necessary.");
	joboptions["sym_name"] = JobOption("Symmetry:", std::string("C1"), "If the molecule is asymmetric, \
set Symmetry group to C1. Note their are multiple possibilities for icosahedral symmetry: \n \
* I1: No-Crowther 222 (standard in Heymann, Chagoyen & Belnap, JSB, 151 (2005) 196–207) \n \
* I2: Crowther 222 \n \
* I3: 52-setting (as used in SPIDER?)\n \
* I4: A different 52 setting \n \
The command 'relion_refine --sym D2 --print_symmetry_ops' prints a list of all symmetry operators for symmetry group D2. \
RELION uses XMIPP's libraries for symmetry operations. \
Therefore, look at the XMIPP Wiki for more details:  http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/WebHome?topic=Symmetry");

	joboptions["bg_diameter"] = JobOption("Diameter background circle (pix): ", -1, -1, 600, 10, "Particles will be normalized to a mean value of zero and a standard-deviation of one for all pixels in the background area.\
The background area is defined as all pixels outside a circle with this given diameter in pixels (before rescaling). When specifying a negative value, a default value of 75% of the Particle box size will be used.");
	joboptions["white_dust"] = JobOption("Stddev for white dust removal: ", -1, -1, 10, 0.1, "Remove very white pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
	joboptions["black_dust"] = JobOption("Stddev for black dust removal: ", -1, -1, 10, 0.1, "Remove very black pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");

	joboptions["do_helix"] = JobOption("Do helical reconstruction?", false, "If set to Yes, then perform 3D helical reconstruction.");
	joboptions["helical_nr_asu"] = JobOption("Number of asymmetrical units:", 1, 1, 100, 1, "Number of helical asymmetrical units in each segment box. If the interbox distance (set in segment picking step) \
is 100 Angstroms and the estimated helical rise is ~20 Angstroms, then set this value to 100 / 20 = 5 (nearest integer). This integer should not be less than 1. The correct value is essential in measuring the \
signal to noise ratio in helical reconstruction.\n\nPlease copy this value from previous 3D refinement.");
	joboptions["helical_twist"] = JobOption("Helical twist (deg):", std::string("0"), "Helical symmetry. Set helical twist (in degrees) to positive value if it is a right-handed helix. \
Helical rise is a positive value in Angstroms.\n\nPlease copy the refined helical symmetry from previous 3D refinement.");
	joboptions["helical_rise"] = JobOption("Helical rise (A):", std::string("0"), "Helical symmetry. Set helical twist (in degrees) to positive value if it is a right-handed helix. \
Helical rise is a positive value in Angstroms.\n\nPlease copy the refined helical symmetry from previous 3D refinement.");

}

bool RelionJob::getCommandsPolishJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "Polish", job_counter);
	std::string command;

	if (joboptions["nr_mpi"].getNumber() > 1)
		command="`which relion_particle_polish_mpi`";
	else
		command="`which relion_particle_polish`";

	// General
	if (joboptions["fn_in"].getString() == "")
	{
		error_message = "ERROR: empty field for input STAR file...";
		return false;
	}
	command += " --i " + joboptions["fn_in"].getString();
	Node node(joboptions["fn_in"].getString(), joboptions["fn_in"].node_type);
	inputNodes.push_back(node);

	if (joboptions["fn_mask"].getString() == "")
	{
		error_message = "ERROR: empty field for input mask...";
		return false;
	}
	command += " --mask " + joboptions["fn_mask"].getString();
	Node node2(joboptions["fn_mask"].getString(), joboptions["fn_mask"].node_type);
	inputNodes.push_back(node2);

	command += " --o " + outputname;
	Node node3(outputname + "shiny.star", NODE_PART_DATA);
	outputNodes.push_back(node3);

	Node node4(outputname + "logfile.pdf", NODE_PDF_LOGFILE);
	outputNodes.push_back(node4);

	Node node5(outputname + "shiny_post.mrc", NODE_FINALMAP);
	outputNodes.push_back(node5);

	// If this is not a continue job, then re-start from scratch....
	if (is_continue)
		command += " --only_do_unfinished ";

	// Beam-induced movement fitting options
	if (joboptions["do_fit_movement"].getBoolean())
		command += " --sigma_nb " + joboptions["sigma_nb"].getString();
	else
		command += " --no_fit ";

	// Damage
	if (joboptions["do_bfactor_weighting"].getBoolean())
	{
		command += " --perframe_highres " + joboptions["perframe_highres"].getString();
		command += " --autob_lowres " + joboptions["perframe_bfac_lowres"].getString();
		command += " --bfactor_running_avg " + joboptions["average_frame_bfactor"].getString();
	}
	else
	{
		command += " --skip_bfactor_weighting ";
	}

	// Symmetry group
	command += " --sym " + joboptions["sym_name"].getString();

	// Normalisation
	command += " --bg_radius " + floatToString(FLOOR(joboptions["bg_diameter"].getNumber()/2.));
	command += " --white_dust " + joboptions["white_dust"].getString();
	command += " --black_dust " + joboptions["black_dust"].getString();

	// Helix
	if (joboptions["do_helix"].getBoolean())
	{
		command += " --helix --helical_nr_asu " + joboptions["helical_nr_asu"].getString();
		command += " --helical_twist " + joboptions["helical_twist"].getString();
		command += " --helical_rise " + joboptions["helical_rise"].getString();
	}

	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);


}

void RelionJob::initialiseMaskcreateJob()
{
	hidden_name = ".gui_maskcreate";

	joboptions["fn_in"] = JobOption("Input 3D map:", NODE_3DREF, "", "MRC map files (*.mrc)", "Provide an input MRC map from which to start binarizing the map.");

	joboptions["lowpass_filter"] = JobOption("Lowpass filter map (A)", -1, 10, 100, 5, "Lowpass filter that will be applied to the input map, prior to binarization. To calculate solvent masks, a lowpass filter of 15-20A may work well.");
	joboptions["angpix"] = JobOption("Pixel size (A)", 1, 0.3, 5, 0.1, "Provide the pixel size in Angstroms of the input map (to calculate the low-pass filter)");

	joboptions["inimask_threshold"] = JobOption("Initial binarisation threshold:", 0.02, 0., 0.5, 0.01, "This threshold is used to make an initial binary mask from the average of the two unfiltered half-reconstructions. \
If you don't know what value to use, display one of the unfiltered half-maps in a 3D surface rendering viewer and find the lowest threshold that gives no noise peaks outside the reconstruction.");
	joboptions["extend_inimask"] = JobOption("Extend binary map this many pixels:", 3, 0, 20, 1, "The initial binary mask is extended this number of pixels in all directions." );
	joboptions["width_mask_edge"] = JobOption("Add a soft-edge of this many pixels:", 3, 0, 20, 1, "The extended binary mask is further extended with a raised-cosine soft edge of the specified width." );

	joboptions["do_helix"] = JobOption("Mask a 3D helix?", false, "Generate a mask for 3D helix which spans across Z axis of the box.");
	joboptions["helical_z_percentage"] = JobOption("Central Z length (%):", 30., 5., 80., 1., "Reconstructed helix suffers from inaccuracies of orientation searches. \
The central part of the box contains more reliable information compared to the top and bottom parts along Z axis. Set this value (%) to the central part length along Z axis divided by the box size. Values around 30% are commonly used but you may want to try different lengths.");

}

bool RelionJob::getCommandsMaskcreateJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "MaskCreate", job_counter);
	std::string command;

	command="`which relion_mask_create`";

	// I/O
	if (joboptions["fn_in"].getString() == "")
	{
		error_message = "ERROR: empty field for input STAR file...";
		return false;
	}
	command += " --i " + joboptions["fn_in"].getString();
	Node node(joboptions["fn_in"].getString(), joboptions["fn_in"].node_type);
	inputNodes.push_back(node);

	command += " --o " + outputname + "mask.mrc";
	Node node2(outputname + "mask.mrc", NODE_MASK);
	outputNodes.push_back(node2);

	if (joboptions["lowpass_filter"].getNumber() > 0)
	{
		command += " --lowpass " + joboptions["lowpass_filter"].getString();
		command += " --angpix " + joboptions["angpix"].getString();
	}
	command += " --ini_threshold " + joboptions["inimask_threshold"].getString();
	command += " --extend_inimask " + joboptions["extend_inimask"].getString();
	command += " --width_soft_edge " + joboptions["width_mask_edge"].getString();

	if (joboptions["do_helix"].getBoolean())
	{
		command += " --helix --z_percentage " + floatToString(joboptions["helical_z_percentage"].getNumber() / 100.);
	}

	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);


}

void RelionJob::initialiseJoinstarJob()
{

	hidden_name = ".gui_joinstar";

	joboptions["do_part"] = JobOption("Combine particle STAR files?", false, "");
	joboptions["fn_part1"] = JobOption("Particle STAR file 1: ", NODE_PART_DATA, "", "particle STAR file (*.star)", "The first of the particle STAR files to be combined.");
	joboptions["fn_part2"] = JobOption("Particle STAR file 2: ", NODE_PART_DATA, "", "particle STAR file (*.star)", "The second of the particle STAR files to be combined.");
	joboptions["fn_part3"] = JobOption("Particle STAR file 3: ", NODE_PART_DATA, "", "particle STAR file (*.star)", "The third of the particle STAR files to be combined. Leave empty if there are only two files to be combined.");
	joboptions["fn_part4"] = JobOption("Particle STAR file 4: ", NODE_PART_DATA, "", "particle STAR file (*.star)", "The fourth of the particle STAR files to be combined. Leave empty if there are only two or three files to be combined.");

	joboptions["do_mic"] = JobOption("Combine micrograph STAR files?", false, "");
	joboptions["fn_mic1"] = JobOption("Micrograph STAR file 1: ", NODE_MICS, "", "micrograph STAR file (*.star)", "The first of the micrograph STAR files to be combined.");
	joboptions["fn_mic2"] = JobOption("Micrograph STAR file 2: ", NODE_MICS, "", "micrograph STAR file (*.star)", "The second of the micrograph STAR files to be combined.");
	joboptions["fn_mic3"] = JobOption("Micrograph STAR file 3: ", NODE_MICS, "", "micrograph STAR file (*.star)", "The third of the micrograph STAR files to be combined. Leave empty if there are only two files to be combined.");
	joboptions["fn_mic4"] = JobOption("Micrograph STAR file 4: ", NODE_MICS, "", "micrograph STAR file (*.star)", "The fourth of the micrograph STAR files to be combined. Leave empty if there are only two or three files to be combined.");

	joboptions["do_mov"] = JobOption("Combine movie STAR files?", false, "");
	joboptions["fn_mov1"] = JobOption("Movie STAR file 1: ", NODE_MOVIES, "", "movie STAR file (*.star)", "The first of the micrograph movie STAR files to be combined.");
	joboptions["fn_mov2"] = JobOption("Movie STAR file 2: ", NODE_MOVIES, "", "movie STAR file (*.star)", "The second of the micrograph movie STAR files to be combined.");
	joboptions["fn_mov3"] = JobOption("Movie STAR file 3: ", NODE_MOVIES, "", "movie STAR file (*.star)", "The third of the micrograph movie STAR files to be combined. Leave empty if there are only two files to be combined.");
	joboptions["fn_mov4"] = JobOption("Movie STAR file 4: ", NODE_MOVIES, "", "movie STAR file (*.star)", "The fourth of the micrograph movie STAR files to be combined. Leave empty if there are only two or three files to be combined.");

}

bool RelionJob::getCommandsJoinstarJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "JoinStar", job_counter);
	std::string command;
	command="`which relion_star_combine`";

	int ii = 0;
	if (joboptions["do_part"].getBoolean())
		ii++;
	if (joboptions["do_mic"].getBoolean())
		ii++;
	if (joboptions["do_mov"].getBoolean())
		ii++;

	if (ii == 0)
	{
		error_message = "You've selected no type of files for joining. Select a single type!";
		return false;
	}

	if (ii > 1)
	{
		error_message = "You've selected more than one type of files for joining. Only select a single type!";
		return false;
	}

	// I/O
	if (joboptions["do_part"].getBoolean())
	{
		if (joboptions["fn_part1"].getString() == "" || joboptions["fn_part2"].getString() == "")
		{
			error_message = "ERROR: empty field for first or second input STAR file...";
			return false;
		}
		command += " --i \" " + joboptions["fn_part1"].getString();
		Node node(joboptions["fn_part1"].getString(), joboptions["fn_part1"].node_type);
		inputNodes.push_back(node);
		command += " " + joboptions["fn_part2"].getString();
		Node node2(joboptions["fn_part2"].getString(), joboptions["fn_part2"].node_type);
		inputNodes.push_back(node2);
		if (joboptions["fn_part3"].getString() != "")
		{
			command += " " + joboptions["fn_part3"].getString();
			Node node3(joboptions["fn_part3"].getString(), joboptions["fn_part3"].node_type);
			inputNodes.push_back(node3);
		}
		if (joboptions["fn_part4"].getString() != "")
		{
			command += " " + joboptions["fn_part4"].getString();
			Node node4(joboptions["fn_part4"].getString(), joboptions["fn_part4"].node_type);
			inputNodes.push_back(node4);
		}
		command += " \" ";

		// Check for duplicates
		command += " --check_duplicates rlnImageName ";
		command += " --o " + outputname + "join_particles.star";
		Node node5(outputname + "join_particles.star", joboptions["fn_part1"].node_type);
		outputNodes.push_back(node5);

	}
	else if (joboptions["do_mic"].getBoolean())
	{
		if (joboptions["fn_mic1"].getString() == "" || joboptions["fn_mic2"].getString() == "")
		{
			error_message = "ERROR: empty field for first or second input STAR file...";
			return false;
		}
		command += " --i \" " + joboptions["fn_mic1"].getString();
		Node node(joboptions["fn_mic1"].getString(), joboptions["fn_mic1"].node_type);
		inputNodes.push_back(node);
		command += " " + joboptions["fn_mic2"].getString();
		Node node2(joboptions["fn_mic2"].getString(), joboptions["fn_mic2"].node_type);
		inputNodes.push_back(node2);
		if (joboptions["fn_mic3"].getString() != "")
		{
			command += " " + joboptions["fn_mic3"].getString();
			Node node3(joboptions["fn_mic3"].getString(), joboptions["fn_mic3"].node_type);
			inputNodes.push_back(node3);
		}
		if (joboptions["fn_mic4"].getString() != "")
		{
			command += " " + joboptions["fn_mic4"].getString();
			Node node4(joboptions["fn_mic4"].getString(), joboptions["fn_mic4"].node_type);
			inputNodes.push_back(node4);
		}
		command += " \" ";

		// Check for duplicates
		command += " --check_duplicates rlnMicrographName ";
		command += " --o " + outputname + "join_mics.star";
		Node node5(outputname + "join_mics.star", joboptions["fn_mic1"].node_type);
		outputNodes.push_back(node5);

	}
	else if (joboptions["do_mov"].getBoolean())
	{
		if (joboptions["fn_mov1"].getString() == "" || joboptions["fn_mov2"].getString() == "")
		{
			error_message = "ERROR: empty field for first or second input STAR file...";
			return false;
		}
		command += " --i \" " + joboptions["fn_mov1"].getString();
		Node node(joboptions["fn_mov1"].getString(), joboptions["fn_mov1"].node_type);
		inputNodes.push_back(node);
		command += " " + joboptions["fn_mov2"].getString();
		Node node2(joboptions["fn_mov2"].getString(), joboptions["fn_mov2"].node_type);
		inputNodes.push_back(node2);
		if (joboptions["fn_mov3"].getString() != "")
		{
			command += " " + joboptions["fn_mov3"].getString();
			Node node3(joboptions["fn_mov3"].getString(), joboptions["fn_mov3"].node_type);
			inputNodes.push_back(node3);
		}
		if (joboptions["fn_mov4"].getString() != "")
		{
			command += " " + joboptions["fn_mov4"].getString();
			Node node4(joboptions["fn_mov4"].getString(), joboptions["fn_mov4"].node_type);
			inputNodes.push_back(node4);
		}
		command += " \" ";

		// Check for duplicates
		command += " --check_duplicates rlnMicrographMovieName ";
		command += " --o " + outputname + "join_movies.star";
		Node node5(outputname + "join_movies.star", joboptions["fn_mov1"].node_type);
		outputNodes.push_back(node5);

	}

	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

void RelionJob::initialiseSubtractJob()
{

	hidden_name = ".gui_subtract";

	joboptions["fn_data"] = JobOption("Input particles:", NODE_PART_DATA, "", "STAR files (*.star)", "The input STAR file with the metadata of all particles.");
	joboptions["do_subtract"] = JobOption("Subtract partial signal?", true, "If set to Yes, a copy of the entire set of particle images in this STAR file will be made that contains the subtracted particle images.");
	joboptions["fn_in"] = JobOption("Map to be projected:", NODE_3DREF, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide the map that will be used to calculate projections, which will be subtracted from the experimental particles. Make sure this map was calculated by RELION from the same particles as above, and preferably with those orientations, as it is crucial that the absolute greyscale is the same as in the experimental particles.");
    joboptions["fn_mask"] = JobOption("Mask to apply to this map:", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a soft mask where the protein density you wish to subtract from the experimental particles is white (1) and the rest of the protein and the solvent is black (0).");
    joboptions["do_fliplabel"] = JobOption("OR revert to original particles?", false, "If set to Yes, no signal subtraction is performed. Instead, the labels of rlnImageName and rlnImageOriginalName are flipped in the input STAR file. This will make the STAR file point back to the original, non-subtracted images.");

	joboptions["do_ctf_correction"] = JobOption("Do apply CTFs?", true, "If set to Yes, CTFs will be applied to the projections to subtract from the experimental particles");
	joboptions["ctf_phase_flipped"] = JobOption("Have data been phase-flipped?", false, "Set this to Yes if the experimental particles have been \
ctf-phase corrected during the pre-processing steps. \
Note that CTF-phase flipping is NOT a necessary pre-processing step for MAP-refinement in RELION, \
as this can be done inside the internal CTF-correction. \
However, if the phases have been flipped, you should tell the program about it by setting this option to Yes.");
	joboptions["ctf_intact_first_peak"] = JobOption("Ignore CTFs until first peak?", false, "Set this to Yes when you have used the same option in the generation of the input map to be projected.");

}

bool RelionJob::getCommandsSubtractJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, "Subtract", job_counter);
	std::string command;

	if (joboptions["do_subtract"].getBoolean())
	{
		command="`which relion_project`";


		// I/O
		if (joboptions["fn_in"].getString() == "")
		{
			error_message = "ERROR: empty field for input map...";
			return false;
		}
		command += " --subtract_exp --i " + joboptions["fn_in"].getString();
		Node node(joboptions["fn_in"].getString(), joboptions["fn_in"].node_type);
		inputNodes.push_back(node);
		if (joboptions["fn_mask"].getString() == "")
		{
			error_message = "ERROR: empty field for input mask...";
			return false;
		}
		command += " --mask " + joboptions["fn_mask"].getString();
		Node node2(joboptions["fn_mask"].getString(), joboptions["fn_mask"].node_type);
		inputNodes.push_back(node2);
		if (joboptions["fn_data"].getString() == "")
		{
			error_message = "ERROR: empty field for input STAR file...";
			return false;
		}
		command += " --ang " + joboptions["fn_data"].getString();
		Node node3(joboptions["fn_data"].getString(), joboptions["fn_data"].node_type);
		inputNodes.push_back(node3);

		command += " --o " + outputname + "subtracted";
		Node node4(outputname + "subtracted.star", NODE_PART_DATA);
		outputNodes.push_back(node4);

		if (joboptions["do_ctf_correction"].getBoolean())
		{
			command += " --ctf --angpix -1";
			if (joboptions["ctf_phase_flipped"].getBoolean())
				command += " --ctf_phase_flip";
			if (joboptions["ctf_intact_first_peak"].getBoolean())
				command += " --ctf_intact_first_peak";
		}

		// Other arguments
		command += " " + joboptions["other_args"].getString();
	}
	else if (joboptions["do_fliplabel"].getBoolean())
	{
		Node node(joboptions["fn_data"].getString(), joboptions["fn_data"].node_type);
		inputNodes.push_back(node);

		Node node2(outputname + "original.star", NODE_PART_DATA);
		outputNodes.push_back(node2);

		command = "awk '{if  ($1==\"_rlnImageName\") {$1=\"_rlnImageOriginalName\"} else if ($1==\"_rlnImageOriginalName\") {$1=\"_rlnImageName\"}; print }' < ";
		command += joboptions["fn_data"].getString() + " > " + outputname + "original.star";
	}
	else
	{
		error_message = "Choose either to subtract particle signal, or to revert to original particles";
		return false;
	}

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}

void RelionJob::initialisePostprocessJob()
{

	hidden_name = ".gui_post";

	joboptions["fn_in"] = JobOption("One of the 2 unfiltered half-maps:", NODE_HALFMAP, "", "MRC map files (*half1_class001_unfil.mrc)",  "Provide one of the two unfiltered half-reconstructions that were output upon convergence of a 3D auto-refine run.");
	joboptions["fn_mask"] = JobOption("Solvent mask:", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a soft mask where the protein is white (1) and the solvent is black (0). Often, the softer the mask the higher resolution estimates you will get. A soft edge of 5-10 pixels is often a good edge width.");
	joboptions["angpix"] = JobOption("Calibrated pixel size (A)", 1, 0.3, 5, 0.1, "Provide the final, calibrated pixel size in Angstroms. This value may be different from the pixel-size used thus far, e.g. when you have recalibrated the pixel size using the fit to a PDB model. The X-axis of the output FSC plot will use this calibrated value.");

	joboptions["fn_mtf"] = JobOption("MTF of the detector (STAR file)", "", "STAR Files (*.star)", ".", "The MTF of the detector is used in the (later) post-processing and particle polishing stages of refinement.  \
If you know the MTF of your detector, provide it here. Curves for some well-known detectors may be downloaded from the RELION Wiki. Also see there for the exact format \
\n If you do not know the MTF of your detector and do not want to measure it, then by leaving this entry empty, you include the MTF of your detector in your overall estimated B-factor upon sharpening the map.\
Although that is probably slightly less accurate, the overall quality of your map will probably not suffer very much.");
	joboptions["do_auto_bfac"] = JobOption("Estimate B-factor automatically?", true, "If set to Yes, then the program will use the automated procedure described by Rosenthal and Henderson (2003, JMB) to estimate an overall B-factor for your map, and sharpen it accordingly. \
Note that your map must extend well beyond the lowest resolution included in the procedure below, which should not be set to resolutions much lower than 10 Angstroms. ");
	joboptions["autob_lowres"] = JobOption("Lowest resolution for auto-B fit (A):", 10, 8, 15, 0.5, "This is the lowest frequency (in Angstroms) that will be included in the linear fit of the Guinier plot as described in Rosenthal and Henderson (2003, JMB). Dont use values much lower or higher than 10 Angstroms. If your map does not extend beyond 10 Angstroms, then instead of the automated procedure use your own B-factor.");
	joboptions["do_adhoc_bfac"] = JobOption("Use your own B-factor?", false, "Instead of using the automated B-factor estimation, provide your own value. Use negative values for sharpening the map. \
This option is useful if your map does not extend beyond the 10A needed for the automated procedure, or when the automated procedure does not give a suitable value (e.g. in more disordered parts of the map).");
	joboptions["adhoc_bfac"] = JobOption("User-provided B-factor:", -1000, -2000, 0, -50, "Use negative values for sharpening. Be careful: if you over-sharpen your map, you may end up interpreting noise for signal!");

	joboptions["do_skip_fsc_weighting"] = JobOption("Skip FSC-weighting?", false, "If set to No (the default), then the output map will be low-pass filtered according to the mask-corrected, gold-standard FSC-curve. \
Sometimes, it is also useful to provide an ad-hoc low-pass filter (option below), as due to local resolution variations some parts of the map may be better and other parts may be worse than the overall resolution as measured by the FSC. \
In such cases, set this option to Yes and provide an ad-hoc filter as described below.");
	joboptions["low_pass"] = JobOption("Ad-hoc low-pass filter (A):",5,1,40,1,"This option allows one to low-pass filter the map at a user-provided frequency (in Angstroms). When using a resolution that is higher than the gold-standard FSC-reported resolution, take care not to interpret noise in the map for signal...");

}

bool RelionJob::getCommandsPostprocessJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "PostProcess", job_counter);
	std::string command;

	command="`which relion_postprocess`";

	// Input mask
	if (joboptions["fn_mask"].getString() == "")
	{
		error_message = "ERROR: empty field for input mask...";
		return false;
	}
	command += " --mask " + joboptions["fn_mask"].getString();
	Node node3(joboptions["fn_mask"].getString(), joboptions["fn_mask"].node_type);
	inputNodes.push_back(node3);

	// Get the input rootname from the half-map name
	// run1_half1_class001_unfil.mrc -> run1
	if (joboptions["fn_in"].getString() == "")
	{
		error_message = "ERROR: empty field for input half-map...";
		return false;
	}
	Node node(joboptions["fn_in"].getString(), joboptions["fn_in"].node_type);
	inputNodes.push_back(node);
	int pos_half = joboptions["fn_in"].getString().rfind("_half");
	if (pos_half < joboptions["fn_in"].getString().size())
	{
		command += " --i " + joboptions["fn_in"].getString().substr(0, pos_half);
	}
	else
	{
		error_message = "PostJobWindow::getCommands ERROR: cannot find _half substring in input filename!";
		return false;
	}

	// The output name contains a directory: use it for output
	command += " --o " + outputname + "postprocess";
	command += "  --angpix " + joboptions["angpix"].getString();
	Node node1(outputname+"postprocess.mrc", NODE_FINALMAP);
	outputNodes.push_back(node1);
	Node node2(outputname+"postprocess_masked.mrc", NODE_FINALMAP);
	outputNodes.push_back(node2);

	Node node2b(outputname+"logfile.pdf", NODE_PDF_LOGFILE);
	outputNodes.push_back(node2b);


	// Sharpening
	if (joboptions["fn_mtf"].getString().length() > 0)
	{
		command += " --mtf " + joboptions["fn_mtf"].getString();
	}
	if (joboptions["do_auto_bfac"].getBoolean())
	{
		command += " --auto_bfac ";
		command += " --autob_lowres " + joboptions["autob_lowres"].getString();
	}
	if (joboptions["do_adhoc_bfac"].getBoolean())
	{
		command += " --adhoc_bfac " + joboptions["adhoc_bfac"].getString();
	}

	// Filtering
	if (joboptions["do_skip_fsc_weighting"].getBoolean())
	{
		command += " --skip_fsc_weighting ";
		command += " --low_pass "  + joboptions["low_pass"].getString();
	}

	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);
	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

void RelionJob::initialiseLocalresJob()
{

	hidden_name = ".gui_localres";

	joboptions["fn_in"] = JobOption("One of the 2 unfiltered half-maps:", NODE_HALFMAP, "", "MRC map files (*_unfil.mrc)",  "Provide one of the two unfiltered half-reconstructions that were output upon convergence of a 3D auto-refine run.");
	joboptions["angpix"] = JobOption("Calibrated pixel size (A)", 1, 0.3, 5, 0.1, "Provide the final, calibrated pixel size in Angstroms. This value may be different from the pixel-size used thus far, e.g. when you have recalibrated the pixel size using the fit to a PDB model. The X-axis of the output FSC plot will use this calibrated value.");

	// Check for environment variable RELION_RESMAP_TEMPLATE
	char * default_location = getenv ("RELION_RESMAP_EXECUTABLE");
	if (default_location == NULL)
	{
		char mydefault[] = DEFAULTRESMAPLOCATION;
		default_location = mydefault;
	}

	joboptions["do_resmap_locres"] = JobOption("Use ResMap?", true, "If set to Yes, then ResMap will be used for local resolution estimation.");
	joboptions["fn_resmap"] = JobOption("ResMap executable:", std::string(default_location), "ResMap*", ".", "Location of the ResMap executable. You can control the default of this field by setting environment variable RELION_RESMAP_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code. \n \n Note that the ResMap wrapper cannot use MPI.");
	joboptions["fn_mask"] = JobOption("User-provided solvent mask:", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a mask with values between 0 and 1 around all domains of the complex.");
	joboptions["pval"] = JobOption("P-value:", 0.05, 0., 1., 0.01, "This value is typically left at 0.05. If you change it, report the modified value in your paper!");
	joboptions["minres"] = JobOption("Highest resolution (A): ", 0., 0., 10., 0.1, "ResMaps minRes parameter. By default (0), the program will start at just above 2x the pixel size");
	joboptions["maxres"] = JobOption("Lowest resolution (A): ", 0., 0., 10., 0.1, "ResMaps maxRes parameter. By default (0), the program will stop at 4x the pixel size");
	joboptions["stepres"] = JobOption("Resolution step size (A)", 1., 0.1, 3, 0.1, "ResMaps stepSize parameter." );

	joboptions["do_relion_locres"] = JobOption("Use Relion?", false, "If set to Yes, then relion_postprocess will be used for local-rtesolution estimation. This program basically performs a series of post-processing operations with a small soft, spherical mask that is moved over the entire map, while using phase-randomisation to estimate the convolution effects of that mask. \
\n \n The output relion_locres.mrc map can be used to color the surface of a map in UCSF Chimera according to its local resolution. The output relion_locres_filtered.mrc is a composite map that is locally filtered to the estimated resolution. \
This is a developmental feature in need of further testing, but initial results indicate it may be useful. \n \n Note that only this program can use MPI, the ResMap wrapper cannot use MPI.");

	//joboptions["locres_sampling"] = JobOption("Sampling rate (A):", 25, 5, 50, 5, "The local-resolution map will be calculated every so many Angstroms, by placing soft spherical masks on a cubic grid with this spacing. Very fine samplings (e.g. < 15A?) may take a long time to compute and give spurious estimates!");
	//joboptions["randomize_at"] = JobOption("Frequency for phase-randomisation (A): ", 10., 5, 20., 1, "From this frequency onwards, the phases for the mask-corrected FSC-calculation will be randomized. Make sure this is a lower resolution (i.e. a higher number) than the local resolutions you are after in your map.");
	joboptions["adhoc_bfac"] = JobOption("User-provided B-factor:", -100, -500, 0, -25, "Probably, the overall B-factor as was estimated in the postprocess is a useful value for here. Use negative values for sharpening. Be careful: if you over-sharpen your map, you may end up interpreting noise for signal!");
	joboptions["fn_mtf"] = JobOption("MTF of the detector (STAR file)", "", "STAR Files (*.star)", ".", "The MTF of the detector is used to complement the user-provided B-factor in the sharpening. If you don't have this curve, you can leave this field empty.");

}

bool RelionJob::getCommandsLocalresJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
	initialisePipeline(outputname, "LocalRes", job_counter);
	std::string command;

	if (joboptions["do_resmap_locres"].getBoolean() == joboptions["do_relion_locres"].getBoolean())
	{
		error_message = "ERROR: choose either ResMap or Relion for local resolution estimation";
		return false;
	}

	if (joboptions["fn_in"].getString() == "")
	{
		error_message = "ERROR: empty field for input half-map...";
		return false;
	}
	// Get the two half-reconstruction names from the single one
	std::string fn_half1, fn_half2;
	int pos_half = joboptions["fn_in"].getString().rfind("_half");
	if (pos_half < joboptions["fn_in"].getString().size())
	{
		fn_half1 = joboptions["fn_in"].getString().substr(0, pos_half) + "_half1_class001_unfil.mrc";
		fn_half2 = joboptions["fn_in"].getString().substr(0, pos_half) + "_half2_class001_unfil.mrc";
	}
	else
	{
		error_message = "ERROR: cannot find _half substring in input filenames...";
		return false;
	}

	Node node(joboptions["fn_in"].getString(), joboptions["fn_in"].node_type);
	inputNodes.push_back(node);

	if (joboptions["do_resmap_locres"].getBoolean())
	{

		// ResMap wrapper
		if (joboptions["fn_resmap"].getString().length() == 0)
		{
			error_message = "ERROR: please provide an executable for the ResMap program.";
			return false;
		}

		if (joboptions["fn_mask"].getString() == "")
		{
			error_message = "ERROR: Please provide an input mask for ResMap local-resolution estimation.";
			return false;
		}

		if (joboptions["do_queue"].getBoolean())
		{
			error_message = "ERROR: You cannot submit a ResMap job to the queue, as it needs user interaction.";
			return false;
		}

		if (joboptions["nr_mpi"].getNumber() > 1)
		{
			error_message = "You cannot use more than 1 MPI processor for the ResMap wrapper...";
			return false;
		}
		// Make symbolic links to the half-maps in the output directory
		commands.push_back("ln -s ../../" + fn_half1 + " " + outputname + "half1.mrc");
		commands.push_back("ln -s ../../" + fn_half2 + " " + outputname + "half2.mrc");

		Node node2(joboptions["fn_mask"].getString(), joboptions["fn_mask"].node_type);
		inputNodes.push_back(node2);

		Node node3(outputname + "half1_resmap.mrc", NODE_RESMAP);
		outputNodes.push_back(node3);

		command = joboptions["fn_resmap"].getString();
		command += " --maskVol=" + joboptions["fn_mask"].getString();
		command += " --noguiSplit " + outputname + "half1.mrc " +  outputname + "half2.mrc";
		command += " --vxSize=" + joboptions["angpix"].getString();
		command += " --pVal=" + joboptions["pval"].getString();
		command += " --minRes=" + joboptions["minres"].getString();
		command += " --maxRes=" + joboptions["maxres"].getString();
		command += " --stepRes=" + joboptions["stepres"].getString();

	}
	else if (joboptions["do_relion_locres"].getBoolean())
	{
		// Relion postprocessing

		if (joboptions["nr_mpi"].getNumber() > 1)
			command="`which relion_postprocess_mpi`";
		else
			command="`which relion_postprocess`";

		command += " --locres --i " + joboptions["fn_in"].getString().substr(0, pos_half);
		command += " --o " + outputname + "relion";
		command += " --angpix " + joboptions["angpix"].getString();
		//command += " --locres_sampling " + joboptions["locres_sampling"].getString();
		//command += " --locres_randomize_at " + joboptions["randomize_at"].getString();
		command += " --adhoc_bfac " + joboptions["adhoc_bfac"].getString();
		if (joboptions["fn_mtf"].getString() != "")
			command += " --mtf " + joboptions["fn_mtf"].getString();

		Node node1(outputname+"relion_locres_filtered.mrc", NODE_FINALMAP);
		outputNodes.push_back(node1);
		Node node2(outputname+"relion_locres.mrc", NODE_RESMAP);
		outputNodes.push_back(node2);

	}

	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();
	commands.push_back(command);


	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

