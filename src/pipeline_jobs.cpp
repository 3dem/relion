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

std::vector<Node> getOutputNodesRefine(std::string outputname, int iter, int K, int dim, int nr_bodies, std::string jobtype)
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
	if (nr_bodies > 1)
	{
		FileName fn_tmp;
		for (int ibody = 0; ibody < nr_bodies; ibody++)
		{
			fn_tmp.compose(fn_out+"_half1_body", ibody+1, "", 3);
			fn_tmp += "_unfil.mrc";
			Node node4(fn_tmp, LABEL_MULTIBODY_HALFMAP);
			result.push_back(node4);
		}
	}
	else // normal refinements/classifications
	{
		if (jobtype == "Refine3D")
		{
			Node node1(fn_out + "_data.star", LABEL_REFINE3D_PARTS);
            Node node2(fn_out + "_optimiser.star", LABEL_REFINE3D_OPT);
    		result.push_back(node1);
            result.push_back(node2);
		}
		if (jobtype == "Class3D")
		{
			Node node1(fn_out + "_data.star", LABEL_CLASS3D_PARTS);
            Node node2(fn_out + "_optimiser.star", LABEL_CLASS3D_OPT);
    		result.push_back(node1);
            result.push_back(node2);
		}
		if (jobtype == "Class2D")
		{
			Node node1(fn_out + "_data.star", LABEL_CLASS2D_PARTS);
            Node node2(fn_out + "_optimiser.star", LABEL_CLASS2D_OPT);
    		result.push_back(node1);
            result.push_back(node2);
		}


                // For auto-refine: also output the run_half1_class001_unfil.mrc map
		if (iter < 0)
		{
			if (jobtype == "Refine3D")
			{
				Node node4(fn_out+"_half1_class001_unfil.mrc", LABEL_REFINE3D_HALFMAP);
				result.push_back(node4);
			}
			if (jobtype == "MultiBody")
			{
				Node node4(fn_out+"_half1_class001_unfil.mrc", LABEL_MULTIBODY_HALFMAP);
				result.push_back(node4);
			}
		}

		// For 3D classification or 3D auto-refine, also use individual 3D maps as outputNodes
		if (dim == 3)
		{
			FileName fn_tmp;
			for (int iclass = 0; iclass < K; iclass++)
			{
				fn_tmp.compose(fn_out+"_class", iclass+1, "mrc", 3);
				if (jobtype == "Refine3D")
				{
					Node node3(fn_tmp, LABEL_REFINE3D_MAP);
					result.push_back(node3);
				}
				if (jobtype == "Class3D")
				{
					Node node3(fn_tmp, LABEL_CLASS3D_MAP);
					result.push_back(node3);
				}
			}
		}
	}

	return result;
}

// Any constructor
JobOption::JobOption(std::string _label, std::string _default_value, std::string _helptext)
{
	clear();
	initialise(_label, _default_value, _helptext);
	joboption_type = JOBOPTION_ANY;
}

// FileName constructor
JobOption::JobOption(std::string _label, std::string  _default_value, std::string _pattern, std::string _directory, std::string _helptext)
{
	clear();
	initialise(_label, _default_value, _helptext);
	joboption_type = JOBOPTION_FILENAME;
	pattern = _pattern;
	directory = _directory;
}

// InputNode constructor
JobOption::JobOption(std::string _label, int _nodetype, std::string _default_value, std::string _pattern, std::string _helptext)
{
	clear();
	initialise(_label, _default_value, _helptext);
	joboption_type = JOBOPTION_INPUTNODE;
	pattern = _pattern;
	node_type = get_node_label(_nodetype);
}

// Radio constructor
JobOption::JobOption(std::string _label, std::vector<std::string> _radio_options, int ioption,  std::string _helptext)
{
	clear();
	std::string defaultval;
	radio_options = _radio_options;

	defaultval = radio_options[ioption];

	initialise(_label, defaultval, _helptext);
	joboption_type = JOBOPTION_RADIO;
}

// Boolean constructor
JobOption::JobOption(std::string _label, bool _boolvalue, std::string _helptext)
{
	clear();
	std::string _default_value = (_boolvalue) ? "Yes" : "No";
	initialise(_label, _default_value, _helptext);
	joboption_type = JOBOPTION_BOOLEAN;
}

// Slider constructor
JobOption::JobOption(std::string _label, float _default_value, float _min_value, float _max_value, float _step_value, std::string _helptext)
{
	clear();
	initialise(_label, floatToString(_default_value), _helptext);
	joboption_type = JOBOPTION_SLIDER;
	min_value = _min_value;
	max_value = _max_value;
	step_value = _step_value;
}

void JobOption::writeToMetaDataTable(MetaDataTable& MD) const
{

	MD.addObject();
	MD.setValue(EMDL_JOBOPTION_VARIABLE, variable);
	MD.setValue(EMDL_JOBOPTION_VALUE, value);

	return;
}


void JobOption::clear()
{
	label = value = default_value = helptext = label_gui = pattern = directory = "undefined";
	joboption_type = JOBOPTION_UNDEFINED;
	radio_options = job_undefined_options;
	node_type = min_value = max_value = step_value = 0.;
}

void JobOption::initialise(std::string _label, std::string _default_value, std::string _helptext)
{
	label = label_gui = _label;
	value = default_value = _default_value;
	helptext = _helptext;
}

bool JobOption::isSchedulerVariable()
{
	std::string pat="$$";
	return (value.find(pat) != std::string::npos);
}

// Get a string value
std::string JobOption::getString()
{
	return value;
}

// Set a string value
void JobOption::setString(std::string set_to)
{
	value = set_to;
}

int JobOption::getHealPixOrder(std::string s)
{
	for (int i = 0; i < 9; i++)
	{
		if (s == job_sampling_options[i])
			return i + 1;
	}

	return -1;
}

std::string JobOption::getCtfFitString(std::string s)
{
	if (s == job_ctffit_options[0]) return "f";
	else if (s == job_ctffit_options[1]) return "m";
	else if (s == job_ctffit_options[2]) return "p";
	else return "";
}

// Get a numbered value
float JobOption::getNumber(std::string &errmsg)
{
	errmsg = "";
	if (value.substr(0, 2) == "$$")
	{
		return 0;
	}
	else
	{
		float retval;
		int ok = sscanf(value.c_str(), "%f", &retval);

		if (ok)
		{
			return retval;
		}
		else
		{
			errmsg = "Error in textToFloat of " + value;
			return 0.;
		}
	}
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
		std::string search_for = label;
		if (label == "Estimate beamtilt?") // 3.0 compatibility
			search_for = "Perform beamtilt estimation?";
		else if (label == "Perform MTF correction?")
		{
			std::cerr << "A legacy job option \"Perform MTF correction?\" is ignored. If an MTF file name is supplied, MTF correction will be applied." << std::endl;
			return false;
		}

		// Start reading the ifstream at the top
		in.clear(); // reset eof if happened...
		in.seekg(0, std::ios::beg);
		std::string line;
		while (getline(in, line, '\n'))
		{
			if (line.rfind(search_for) == 0)
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

bool RelionJob::containsLabel(std::string _label, std::string &option)
{
	for (std::map<std::string,JobOption>::iterator it=joboptions.begin(); it!=joboptions.end(); ++it)
	{
		if ((it->second).label == _label)
		{
			option = it->first;
			return true;
		}
	}
	return false;
}

void RelionJob::setOption(std::string setOptionLine)
{
	std::size_t equalsigns = setOptionLine.find("==");
	if (equalsigns == std::string::npos)
		REPORT_ERROR(" ERROR: no '==' entry on JobOptionLine: " + setOptionLine);

	std::string label, value, option;
	label = setOptionLine.substr(0, equalsigns - 1);
	value = setOptionLine.substr(equalsigns + 3, setOptionLine.length() - equalsigns - 3);

	if (joboptions.find(label) != joboptions.end())
	{
		joboptions[label].setString(value);
	}
	else if (containsLabel(label, option))
	{
		joboptions[option].setString(value);
	}
	else
	{
		REPORT_ERROR(" ERROR: Job does not contain label: " + label);
	}
}

bool RelionJob::read(std::string fn, bool &_is_continue, bool do_initialise)
{
	// If fn is empty, use the hidden name
	FileName myfilename = (fn=="") ? hidden_name : fn;
	bool have_read = false;

	// For backwards compatibility
	if (!exists(myfilename + "job.star") && exists(myfilename + "run.job"))
	{
		std::ifstream fh;
		fh.open((myfilename+"run.job").c_str(), std::ios_base::in);
		if (fh.fail())
		{
			REPORT_ERROR("ERROR reading file: " + myfilename + "run.job");
		}
		else
		{
			std::string line;

			// Get job type from first line
			getline(fh, line, '\n');
			size_t idx = line.find("==");
			idx++;

			type = (int)textToFloat((line.substr(idx+1,line.length()-idx)).c_str());

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
			have_read = true;
		}

		fh.close();
	}

	if (!have_read)
	{
		// Read from STAR
		MetaDataTable MDhead;
		MetaDataTable MDvals;

		FileName fn_star = myfilename;
		if (fn_star.getExtension() != "star" || !exists(fn_star)) // full name was given
		{
			fn_star += "job.star"; // "Refine3D/job123" OR ".gui_auto3d"
			if (!exists(fn_star))
				return false;
		}

		MDhead.read(fn_star, "job");
		if (MDhead.containsLabel(EMDL_JOB_TYPE_LABEL))
		{
			MDhead.getValue(EMDL_JOB_TYPE_LABEL, label);
			type = get_proc_type(label);
		}
		else
		{
			// backwards compatibility with 3.0
			MDhead.getValue(EMDL_JOB_TYPE, type);
		}

		MDhead.getValue(EMDL_JOB_IS_CONTINUE, is_continue);
		_is_continue = is_continue;
		MDhead.getValue(EMDL_JOB_IS_TOMO, is_tomo);
		if (do_initialise)
			initialise(type);

		MDvals.read(fn_star, "joboptions_values");
		std::string label, value;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDvals)
		{
			MDvals.getValue(EMDL_JOBOPTION_VARIABLE, label);
			MDvals.getValue(EMDL_JOBOPTION_VALUE, value);
			if (joboptions.find(label) == joboptions.end())
			{
				std::cerr << "WARNING: cannot find " << label << " in the defined joboptions. Ignoring it ..." <<std::endl;
			}
			else
			{
				joboptions[label].value = value;
			}
		}
		have_read = true;
	}

	if (have_read)
	{
		// Just check that went OK
		if (type != PROC_IMPORT &&
		    type != PROC_MOTIONCORR &&
		    type != PROC_CTFFIND &&
		    type != PROC_MANUALPICK &&
		    type != PROC_AUTOPICK &&
		    type != PROC_EXTRACT &&
		    type != PROC_CLASSSELECT &&
		    type != PROC_2DCLASS &&
		    type != PROC_3DCLASS &&
		    type != PROC_3DAUTO &&
		    type != PROC_MULTIBODY &&
		    type != PROC_MASKCREATE &&
		    type != PROC_JOINSTAR &&
		    type != PROC_SUBTRACT &&
		    type != PROC_POST &&
		    type != PROC_RESMAP &&
		    type != PROC_INIMODEL &&
		    type != PROC_MOTIONREFINE &&
		    type != PROC_CTFREFINE &&
			type != PROC_TOMO_IMPORT &&
			type != PROC_TOMO_SUBTOMO &&
			type != PROC_TOMO_CTFREFINE &&
			type != PROC_TOMO_ALIGN &&
			type != PROC_TOMO_RECONSTRUCT &&
		    type != PROC_EXTERNAL)
			return false;

		return true;
	}
	else
	{
		return false;
	}
}

void RelionJob::write(std::string fn)
{
	// If fn is empty, use the hidden name
	FileName myfilename = (fn=="") ? hidden_name : fn;

	FileName fn_star = myfilename;
	if (fn_star.getExtension() != "star") // full name was given
	{
		fn_star += "job.star"; // "Refine3D/job123" OR ".gui_auto3d"
	}

	std::ofstream fh;
	fh.open((fn_star).c_str(), std::ios::out);
	if (!fh)
		REPORT_ERROR("ERROR: Cannot write to file: " + fn_star);

	MetaDataTable MDhead;
	MetaDataTable MDvals, MDopts;

	MDhead.setIsList(true);
	MDhead.addObject();
	// as of 3.1-beta do not write integer into the STAR files anymore....
	// MDhead.setValue(EMDL_JOB_TYPE, type);
	MDhead.setValue(EMDL_JOB_TYPE_LABEL, label);
	MDhead.setValue(EMDL_JOB_IS_CONTINUE, is_continue);
	MDhead.setValue(EMDL_JOB_IS_TOMO, is_tomo);
	// TODO: add name for output directories!!! make a std:;map between type and name for all options!
	//MDhead.setValue(EMDL_JOB_TYPE_NAME, type);
	MDhead.setName("job");
	MDhead.write(fh);

	// Now make a table with all the values
	for (std::map<std::string,JobOption>::iterator it=joboptions.begin(); it!=joboptions.end(); ++it)
	{
		(it->second).writeToMetaDataTable(MDvals);
	}
	MDvals.setName("joboptions_values");
	MDvals.write(fh);

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
		int nmpi = (joboptions.find("nr_mpi") != joboptions.end()) ? joboptions["nr_mpi"].getNumber(error_message) : 1;
		if (error_message != "") return false;

		int nthr = (joboptions.find("nr_threads") != joboptions.end()) ? joboptions["nr_threads"].getNumber(error_message) : 1;
		if (error_message != "") return false;

		int ncores = nmpi * nthr;
		int ndedi = joboptions["min_dedicated"].getNumber(error_message);
		if (error_message != "") return false;

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
		char *extra_count_text = getenv("RELION_QSUB_EXTRA_COUNT");
		const char extra_count_val = (extra_count_text ? atoi(extra_count_text) : 2);
		for (int i=1; i<=extra_count_val; i++)
		{
			std::stringstream out;
			out<<i;
			const std::string i_str=out.str();
			if (joboptions.find(std::string("qsub_extra")+i_str) != joboptions.end())
			{
				replacing[std::string("XXXextra")+i_str+"XXX"] = joboptions[std::string("qsub_extra")+i_str].getString();
			}
		}

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

void RelionJob::initialisePipeline(std::string &outputname, int job_counter)
{
	outputNodes.clear();
	inputNodes.clear();

	FileName dirname = proc_type2dirname.at(type);
	// TODO: insert "relion." prefix to dirname when using the ccpem-pipeliner...

	if (outputname == "") // for continue jobs, use the same outputname
	{
		if (job_counter < 1000)
			outputname = dirname + "/job" + integerToString(job_counter, 3) + "/";
		else
			outputname = dirname + "/job" + integerToString(job_counter) + "/";
	}

	// This is the default label, deeper levels can be added for specific jobs
	label = get_proc_label(type);
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
			mktree(outputname.substr(0, last_slash));
		}
	}

	// Add the --pipeline_control argument to all relion_ programs
	for (int icom = 0; icom < commands.size(); icom++)
	{
		if ((commands[icom]).find("relion_") != std::string::npos)
			commands[icom] += " --pipeline_control " + outputname;

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
			nr_mpi = (joboptions.find("nr_mpi") != joboptions.end()) ? joboptions["nr_mpi"].getNumber(error_message) : 1;
			if (error_message != "") return false;

			if (nr_mpi > 1 &&
					(commands[icom]).find("_mpi`") != std::string::npos &&
					(commands[icom]).find("relion_") != std::string::npos)
			{

				const char *default_mpirun = getenv("RELION_MPIRUN");
				if (default_mpirun == NULL)
				{
					default_mpirun = DEFAULTMPIRUN;
				}
				one_command = std::string(default_mpirun) + " -n " + floatToString(nr_mpi) + " " + commands[icom] ;
			}
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

	char * my_warn = getenv("RELION_ERROR_LOCAL_MPI");
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
	if (type == PROC_IMPORT)
	{
		has_mpi = has_thread = false;
		initialiseImportJob();
	}
	else if (type == PROC_MOTIONCORR)
	{
		has_mpi = has_thread = true;
		initialiseMotioncorrJob();
	}
	else if (type == PROC_CTFFIND)
	{
		has_mpi = true;
		has_thread = false;
		initialiseCtffindJob();
	}
	else if (type == PROC_MANUALPICK)
	{
		has_mpi = has_thread = false;
		initialiseManualpickJob();
	}
	else if (type == PROC_AUTOPICK)
	{
		has_mpi = true;
		has_thread = false;
		initialiseAutopickJob();
	}
	else if (type == PROC_EXTRACT)
	{
		has_mpi = true;
		has_thread = false;
		initialiseExtractJob();
	}
	else if (type == PROC_CLASSSELECT)
	{
		has_mpi = has_thread = false;
		initialiseSelectJob();
	}
	else if (type == PROC_2DCLASS)
	{
		has_mpi = has_thread = true;
		initialiseClass2DJob();
	}
	else if (type == PROC_INIMODEL)
	{
		has_mpi = has_thread = true;
		initialiseInimodelJob();
	}
	else if (type == PROC_3DCLASS)
	{
		has_mpi = has_thread = true;
		initialiseClass3DJob();
	}
	else if (type == PROC_3DAUTO)
	{
		has_mpi = has_thread = true;
		initialiseAutorefineJob();
	}
	else if (type == PROC_MULTIBODY)
	{
		has_mpi = has_thread = true;
		initialiseMultiBodyJob();
	}
	else if (type == PROC_MASKCREATE)
	{
		has_mpi = false;
		has_thread = true;
		initialiseMaskcreateJob();
	}
	else if (type == PROC_JOINSTAR)
	{
		has_mpi = has_thread = false;
		initialiseJoinstarJob();
	}
	else if (type == PROC_SUBTRACT)
	{
		has_mpi = true;
		has_thread = false;
		initialiseSubtractJob();
	}
	else if (type == PROC_POST)
	{
		has_mpi = has_thread = false;
		initialisePostprocessJob();
	}
	else if (type == PROC_RESMAP)
	{
            has_mpi = true;
            has_thread = false;
		initialiseLocalresJob();
	}
	else if (type == PROC_MOTIONREFINE)
	{
		has_mpi = has_thread = true;
		initialiseMotionrefineJob();
	}
	else if (type == PROC_CTFREFINE)
	{
		has_mpi = has_thread = true;
		initialiseCtfrefineJob();
	}
	else if (type == PROC_TOMO_IMPORT)
	{
		has_mpi = has_thread = false;
		initialiseTomoImportJob();
	}
	else if (type == PROC_TOMO_SUBTOMO)
	{
		has_mpi = has_thread = true;
		initialiseTomoSubtomoJob();
	}
	else if (type == PROC_TOMO_CTFREFINE)
	{
		has_mpi = has_thread = true;
		initialiseTomoCtfRefineJob();
	}
	else if (type == PROC_TOMO_ALIGN)
	{
		has_mpi = has_thread = true;
		initialiseTomoAlignJob();
	}
	else if (type == PROC_TOMO_RECONSTRUCT)
	{
		has_mpi = has_thread = true;
		initialiseTomoReconPartJob();
	}
	else if (type == PROC_EXTERNAL)
	{
		has_mpi = false;
		has_thread = true;
		initialiseExternalJob();
	}
	else
		REPORT_ERROR("ERROR: unrecognised job-type");

	// Check for environment variable RELION_MPI_MAX and RELION_QSUB_NRMPI
	const char *mpi_max_input = getenv("RELION_MPI_MAX");
	int mpi_max = (mpi_max_input == NULL) ? DEFAULTMPIMAX : textToInteger(mpi_max_input);
	char * qsub_nrmpi_text = getenv("RELION_QSUB_NRMPI");
	const char qsub_nrmpi_val = (qsub_nrmpi_text ? atoi(qsub_nrmpi_text) : DEFAULTNRMPI);
	if (has_mpi)
	{
		joboptions["nr_mpi"] = JobOption("Number of MPI procs:", qsub_nrmpi_val , 1, mpi_max, 1, "Number of MPI nodes to use in parallel. When set to 1, MPI will not be used. The maximum can be set through the environment variable RELION_MPI_MAX.");
	}

	const char *thread_max_input = getenv("RELION_THREAD_MAX");
	int thread_max = (thread_max_input == NULL) ? DEFAULTTHREADMAX : textToInteger(thread_max_input);
	char * qsub_nrthr_text = getenv("RELION_QSUB_NRTHREADS");
	const char qsub_nrthreads_val = (qsub_nrthr_text ? atoi(qsub_nrthr_text) : DEFAULTNRTHREADS);
	if (has_thread)
	{
		joboptions["nr_threads"] = JobOption("Number of threads:", qsub_nrthreads_val, 1, thread_max, 1, "Number of shared-memory (POSIX) threads to use in parallel. \
When set to 1, no multi-threading will be used. The maximum can be set through the environment variable RELION_THREAD_MAX.");
	}

	const char * use_queue_input = getenv("RELION_QUEUE_USE");
	bool use_queue = (use_queue_input == NULL) ? DEFAULTQUEUEUSE : textToBool(use_queue_input);
	joboptions["do_queue"] = JobOption("Submit to queue?", use_queue, "If set to Yes, the job will be submit to a queue, otherwise \
the job will be executed locally. Note that only MPI jobs may be sent to a queue. The default can be set through the environment variable RELION_QUEUE_USE.");

	// Check for environment variable RELION_QUEUE_NAME
	const char * default_queue = getenv("RELION_QUEUE_NAME");
	if (default_queue==NULL)
	{
		default_queue = DEFAULTQUEUENAME;
	}

	// Need the std::string(), as otherwise it will be overloaded and passed as a boolean....
	joboptions["queuename"] = JobOption("Queue name: ", std::string(default_queue), "Name of the queue to which to submit the job. The default name can be set through the environment variable RELION_QUEUE_NAME.");

	// Check for environment variable RELION_QSUB_COMMAND
	const char * default_command = getenv("RELION_QSUB_COMMAND");
	if (default_command==NULL)
	{
		default_command = DEFAULTQSUBCOMMAND;
	}

	joboptions["qsub"] = JobOption("Queue submit command:", std::string(default_command), "Name of the command used to submit scripts to the queue, e.g. qsub or bsub.\n\n\
Note that the person who installed RELION should have made a custom script for your cluster/queue setup. Check this is the case \
(or create your own script following the RELION Wiki) if you have trouble submitting jobs. The default command can be set through the environment variable RELION_QSUB_COMMAND.");

	// additional options that may be set through environment variables RELION_QSUB_EXTRAi and RELION_QSUB_EXTRAi (for more flexibility)
	char * extra_count_text = getenv("RELION_QSUB_EXTRA_COUNT");
	const char extra_count_val = (extra_count_text ? atoi(extra_count_text) : 2);
	for (int i=1; i<=extra_count_val; i++)
	{
		std::stringstream out;
		out<<i;
		const std::string i_str=out.str();
		char * extra_text = getenv((std::string("RELION_QSUB_EXTRA")+i_str).c_str());
		if (extra_text != NULL)
		{
			const std::string query_default=std::string("RELION_QSUB_EXTRA")+i_str+"_DEFAULT";
			char *extra_default = getenv(query_default.c_str());
			char emptychar[] = "";
			if (extra_default == NULL)
			{
				extra_default=emptychar;
			}
			const std::string query_help=std::string("RELION_QSUB_EXTRA")+i_str+"_HELP";
			char *extra_help = getenv(query_help.c_str());
			std::string txt;
			if (extra_help == NULL)
			{
				txt = std::string("Extra option to pass to the qsub template script. Any occurrences of XXXextra")+i_str+"XXX will be changed by this value.";
			}
			else
			{
				txt=std::string(extra_help);
			}
			joboptions[std::string("qsub_extra")+i_str] = JobOption(std::string(extra_text), std::string(extra_default), txt.c_str());
                }
	}

	// Check for environment variable RELION_QSUB_TEMPLATE
	char * default_location = getenv("RELION_QSUB_TEMPLATE");
	char default_qsub[] = DEFAULTQSUBLOCATION;
	if (default_location == NULL)
	{
		default_location = default_qsub;
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
If these options are not enough for your standard jobs, you may define a user-specified number of extra variables: XXXextra1XXX, XXXextra2XXX, etc. \
The number of extra variables is controlled through the environment variable RELION_QSUB_EXTRA_COUNT. \
Their help text is set by the environment variables RELION_QSUB_EXTRA1, RELION_QSUB_EXTRA2, etc \
For example, setenv RELION_QSUB_EXTRA_COUNT 1, together with setenv RELION_QSUB_EXTRA1 \"Max number of hours in queue\" will result in an additional (text) ein the GUI \
Any variables XXXextra1XXX in the template script will be replaced by the corresponding value.\
Likewise, default values for the extra entries can be set through environment variables RELION_QSUB_EXTRA1_DEFAULT, RELION_QSUB_EXTRA2_DEFAULT, etc. \
But note that (unlike all other entries in the GUI) the extra values are not remembered from one run to the other.");

	// Check for environment variable RELION_QSUB_TEMPLATE
	char * my_minimum_dedicated = getenv ("RELION_MINIMUM_DEDICATED");
	int minimum_nr_dedicated = (my_minimum_dedicated == NULL) ? DEFAULTMININIMUMDEDICATED : textToInteger(my_minimum_dedicated);
	joboptions["min_dedicated"] = JobOption("Minimum dedicated cores per node:", minimum_nr_dedicated, 1, 64, 1, "Minimum number of dedicated cores that need to be requested on each node. This is useful to force the queue to fill up entire nodes of a given size. The default can be set through the environment variable RELION_MINIMUM_DEDICATED.");

	// Need the std::string(), as otherwise it will be overloaded and passed as a boolean....
	joboptions["other_args"] = JobOption("Additional arguments:", std::string(""), "In this box command-line arguments may be provided that are not generated by the GUI. \
This may be useful for testing developmental options and/or expert use of the program. \
To print a list of possible options, run the corresponding program from the command line without any arguments.");

	// Set the variable name in all joboptions

	std::map<std::string, JobOption>::iterator it;
	for ( it = joboptions.begin(); it != joboptions.end(); it++ )
	{
	    (it->second).variable = it->first;
	}
}

bool RelionJob::getCommands(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	bool result = false;

	if (type == PROC_IMPORT)
	{
		result = getCommandsImportJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_MOTIONCORR)
	{
		result = getCommandsMotioncorrJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_CTFFIND)
	{
		result = getCommandsCtffindJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_MANUALPICK)
	{
		result = getCommandsManualpickJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_AUTOPICK)
	{
		result = getCommandsAutopickJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_EXTRACT)
	{
		result = getCommandsExtractJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_CLASSSELECT)
	{
		result = getCommandsSelectJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_2DCLASS)
	{
		result = getCommandsClass2DJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_INIMODEL)
	{
		result = getCommandsInimodelJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_3DCLASS)
	{
		result = getCommandsClass3DJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_3DAUTO)
	{
		result = getCommandsAutorefineJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_MULTIBODY)
	{
		result = getCommandsMultiBodyJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_MASKCREATE)
	{
		result = getCommandsMaskcreateJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_JOINSTAR)
	{
		result = getCommandsJoinstarJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_SUBTRACT)
	{
		result = getCommandsSubtractJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_POST)
	{
		result = getCommandsPostprocessJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_RESMAP)
	{
		result = getCommandsLocalresJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_MOTIONREFINE)
	{
		result = getCommandsMotionrefineJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_CTFREFINE)
	{
		result = getCommandsCtfrefineJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_TOMO_IMPORT)
	{
		result = getCommandsTomoImportJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_TOMO_SUBTOMO)
	{
		result = getCommandsTomoSubtomoJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_TOMO_CTFREFINE)
	{
		result = getCommandsTomoCtfRefineJob(outputname, commands, final_command, do_makedir, job_counter,
											 error_message);
	}
	else if (type == PROC_TOMO_ALIGN)
	{
		result = getCommandsTomoAlignJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else if (type == PROC_TOMO_RECONSTRUCT)
	{
		result = getCommandsTomoReconPartJob(outputname, commands, final_command, do_makedir, job_counter,
											 error_message);
	}
	else if (type == PROC_EXTERNAL)
	{
		result = getCommandsExternalJob(outputname, commands, final_command, do_makedir, job_counter, error_message);
	}
	else
	{
		REPORT_ERROR("ERROR: unrecognised job-type: type = " + integerToString(type));
	}

	return result;
}

void RelionJob::initialiseImportJob()
{
	hidden_name = ".gui_import";

	joboptions["do_raw"] = JobOption("Import raw movies/micrographs?", true, "Set this to Yes if you plan to import raw movies or micrographs");
	joboptions["fn_in_raw"] = JobOption("Raw input files:", "Micrographs/*.tif", (std::string)"Movie or Image (*.{mrc,mrcs,tif,tiff})", ".", "Provide a Linux wildcard that selects all raw movies or micrographs to be imported. The path must be a relative path from the project directory. To import files outside the project directory, first make a symbolic link by an absolute path and then specify the link by a relative path. See the FAQ page on RELION wiki (https://www3.mrc-lmb.cam.ac.uk/relion/index.php/FAQs#What_is_the_right_way_to_import_files_outside_the_project_directory.3F) for details.");
	joboptions["is_multiframe"] = JobOption("Are these multi-frame movies?", true, "Set to Yes for multi-frame movies, set to No for single-frame micrographs.");

	joboptions["optics_group_name"] = JobOption("Optics group name:", (std::string)"opticsGroup1", "Name of this optics group. Each group of movies/micrographs with different optics characteristics for CTF refinement should have a unique name.");
	joboptions["fn_mtf"] = JobOption("MTF of the detector:", "", "STAR Files (*.star)", ".", "As of release-3.1, the MTF of the detector is used in the refinement stages of refinement.  \
If you know the MTF of your detector, provide it here. Curves for some well-known detectors may be downloaded from the RELION Wiki. Also see there for the exact format \
\n If you do not know the MTF of your detector and do not want to measure it, then by leaving this entry empty, you include the MTF of your detector in your overall estimated B-factor upon sharpening the map.\
Although that is probably slightly less accurate, the overall quality of your map will probably not suffer very much. \n \n Note that when combining data from different detectors, the differences between their MTFs can no longer be absorbed in a single B-factor, and providing the MTF here is important!");

	joboptions["angpix"] = JobOption("Pixel size (Angstrom):", 1.4, 0.5, 3, 0.1, "Pixel size in Angstroms. ");
	joboptions["kV"] = JobOption("Voltage (kV):", 300, 50, 500, 10, "Voltage the microscope was operated on (in kV)");
	joboptions["Cs"] = JobOption("Spherical aberration (mm):", 2.7, 0, 8, 0.1, "Spherical aberration of the microscope used to collect these images (in mm). Typical values are 2.7 (FEI Titan & Talos, most JEOL CRYO-ARM), 2.0 (FEI Polara), 1.4 (some JEOL CRYO-ARM) and 0.01 (microscopes with a Cs corrector).");
	joboptions["Q0"] = JobOption("Amplitude contrast:", 0.1, 0, 0.3, 0.01, "Fraction of amplitude contrast. Often values around 10% work better than theoretically more accurate lower values...");
	joboptions["beamtilt_x"] = JobOption("Beamtilt in X (mrad):", 0.0, -1.0, 1.0, 0.1, "Known beamtilt in the X-direction (in mrad). Set to zero if unknown.");
	joboptions["beamtilt_y"] = JobOption("Beamtilt in Y (mrad):", 0.0, -1.0, 1.0, 0.1, "Known beamtilt in the Y-direction (in mrad). Set to zero if unknown.");


	joboptions["do_other"] = JobOption("Import other node types?", false, "Set this to Yes  if you plan to import anything else than movies or micrographs");

	joboptions["fn_in_other"] = JobOption("Input file:", "ref.mrc", "Input file (*.*)", ".", "Select any file(s) to import. \n \n \
Note that for importing coordinate files, one has to give a Linux wildcard, where the *-symbol is before the coordinate-file suffix, e.g. if the micrographs are called mic1.mrc and the coordinate files mic1.box or mic1_autopick.star, one HAS to give '*.box' or '*_autopick.star', respectively.\n \n \
Also note that micrographs, movies and coordinate files all need to be in the same directory (with the same rootnames, e.g.mic1 in the example above) in order to be imported correctly. 3D masks or references can be imported from anywhere. \n \n \
Note that movie-particle STAR files cannot be imported from a previous version of RELION, as the way movies are handled has changed in RELION-2.0. \n \n \
For the import of a particle, 2D references or micrograph STAR file or of a 3D reference or mask, only a single file can be imported at a time. \n \n \
Note that due to a bug in a fltk library, you cannot import from directories that contain a substring  of the current directory, e.g. dont important from /home/betagal if your current directory is called /home/betagal_r2. In this case, just change one of the directory names.");

	joboptions["node_type"] = JobOption("Node type:", job_nodetype_options, 0, "Select the type of Node this is.");
	joboptions["optics_group_particles"] = JobOption("Rename optics group for particles:", (std::string)"", "Only for the import of a particles STAR file with a single, or no, optics groups defined: rename the optics group for the imported particles to this string.");
}

// Generate the correct commands
bool RelionJob::getCommandsImportJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);

	std::string command;
	FileName fn_out, fn_in;
	command = "relion_import ";

	bool do_raw = joboptions["do_raw"].getBoolean();
	bool do_other = joboptions["do_other"].getBoolean();

	if (do_raw && do_other)
	{
		error_message = "ERROR: you cannot import BOTH raw movies/micrographs AND other node types at the same time...";
		return false;
	}
	if ((!do_raw) && (!do_other))
	{
		error_message = "ERROR: nothing to do... ";
		return false;
	}

	if (do_raw)
	{
		label += ".movies";

		fn_in = joboptions["fn_in_raw"].getString();

		if (fn_in.rfind("../") != std::string::npos) // Forbid at any place
		{
			error_message = "ERROR: don't import files outside the project directory.\nPlease make a symbolic link by an absolute path before importing.";
			return false;
		}

		if (fn_in.rfind("/", 0) == 0) // Forbid only at the beginning
		{
			error_message = "ERROR: please import files by a relative path.\nIf you want to import files outside the project directory, make a symbolic link by an absolute path and\nimport the symbolic link by a relative path.";
			return false;
		}

		if (joboptions["is_multiframe"].getBoolean())
		{
			fn_out = "movies.star";
			Node node(outputname + fn_out, LABEL_IMPORT_MOVIES);
			outputNodes.push_back(node);
			command += " --do_movies ";
		}
		else
		{
			fn_out = "micrographs.star";
			Node node(outputname + fn_out, LABEL_IMPORT_MICS);
			outputNodes.push_back(node);
			command += " --do_micrographs ";
		}

		FileName optics_group = joboptions["optics_group_name"].getString();
		if (optics_group == "")
		{
			error_message = "ERROR: please specify an optics group name.";
			return false;
		}

		if (!optics_group.validateCharactersStrict(true)) // true means: do_allow_double_dollar (for scheduler)
		{
			error_message = "ERROR: an optics group name may contain only numbers, alphabets and hyphen(-).";
			return false;
		}

		command += " --optics_group_name \"" + optics_group + "\"";
		if (joboptions["fn_mtf"].getString().length() > 0)
		{
			command += " --optics_group_mtf " + joboptions["fn_mtf"].getString();
		}
		command += " --angpix " + joboptions["angpix"].getString();
		command += " --kV " + joboptions["kV"].getString();
		command += " --Cs " + joboptions["Cs"].getString();
		command += " --Q0 " + joboptions["Q0"].getString();
		command += " --beamtilt_x " + joboptions["beamtilt_x"].getString();
		command += " --beamtilt_y " + joboptions["beamtilt_y"].getString();

	}
	else if (do_other)
	{
		label += ".other";

		fn_in = joboptions["fn_in_other"].getString();
		std::string node_type = joboptions["node_type"].getString();
		if (node_type == "Particle coordinates (*.box, *_pick.star)")
		{
			// Make a suffix file, which contains the actual suffix as a suffix
			// Get the coordinate-file suffix
			fn_out = "coords_suffix" + fn_in.afterLastOf("*");
			Node node(outputname + fn_out, LABEL_IMPORT_COORDS);
			outputNodes.push_back(node);
			command += " --do_coordinates ";
		}
		else
		{
			fn_out = "/" + fn_in;
			fn_out = fn_out.afterLastOf("/");

			std::string mynodetype;
			if (node_type == "Particles STAR file (.star)")
				mynodetype = LABEL_IMPORT_PARTS;
			else if (node_type == "Multiple (2D or 3D) references (.star or .mrcs)")
				mynodetype = LABEL_IMPORT_2DIMG;
			else if (node_type == "3D reference (.mrc)")
				mynodetype = LABEL_IMPORT_MAP;
			else if (node_type == "3D mask (.mrc)")
				mynodetype = LABEL_IMPORT_MASK;
			else if (node_type == "Micrographs STAR file (.star)")
				mynodetype = LABEL_IMPORT_MICS;
			else if (node_type == "Unfiltered half-map (unfil.mrc)")
				mynodetype = LABEL_IMPORT_HALFMAP;
			else
			{
				error_message = "Unrecognized menu option for node_type = " + node_type;
				return false;
			}

			Node node(outputname + fn_out, mynodetype);
			outputNodes.push_back(node);

			// Also get the other half-map
			if (mynodetype == LABEL_HALFMAP_CPIPE)
			{
				FileName fn_inb = "/" + fn_in;
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
				fn_inb = fn_inb.afterLastOf("/");
				Node node2(outputname + fn_inb, mynodetype);
				outputNodes.push_back(node2);
				command += " --do_halfmaps";
			}
			else if (mynodetype == LABEL_PARTS_CPIPE)
			{
				command += " --do_particles";
				FileName optics_group = joboptions["optics_group_particles"].getString();
				if (optics_group != "")
				{
					if (!optics_group.validateCharactersStrict())
					{
						error_message = "ERROR: an optics group name may contain only numbers, alphabets and hyphen(-).";
						return false;
					}
					command += " --particles_optics_group_name \"" + optics_group + "\"";
				}
			}
			else
			{
				command += " --do_other";
			}
		}
	}

	// Now finish the command call to relion_import program, which does the actual copying
	command += " --i \"" + fn_in + "\"";
	command += " --odir " + outputname;
	command += " --ofile " + fn_out;

	if (is_continue)
		command += " --continue ";

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}

void RelionJob::initialiseMotioncorrJob()
{
	hidden_name = ".gui_motioncorr";

	joboptions["input_star_mics"] = JobOption("Input movies STAR file:", NODE_MOVIES_CPIPE, "", "STAR files (*.star)", "A STAR file with all micrographs to run MOTIONCORR on");
	joboptions["first_frame_sum"] = JobOption("First frame for corrected sum:", 1, 1, 32, 1, "First frame to use in corrected average (starts counting at 1). ");
	joboptions["last_frame_sum"] = JobOption("Last frame for corrected sum:", -1, 0, 32, 1, "Last frame to use in corrected average. Values equal to or smaller than 0 mean 'use all frames'.");
	joboptions["eer_grouping"] = JobOption("EER fractionation:", 32, 1, 100, 1, "The number of hardware frames to group into one fraction. This option is relevant only for Falcon4 movies in the EER format. Note that all 'frames' in the GUI (e.g. first and last frame for corrected sum, dose per frame) refer to fractions, not raw detector frames. See https://www3.mrc-lmb.cam.ac.uk/relion/index.php/Image_compression#Falcon4_EER for detailed guidance on EER processing.");
	joboptions["do_float16"] = JobOption("Write output in float16?", true ,"If set to Yes, RelionCor2 will write output images in float16 MRC format. This will save a factor of two in disk space compared to the default of writing in float32. Note that RELION and CCPEM will read float16 images, but other programs may not (yet) do so. For example, Gctf will not work with float16 images. Also note that this option does not work with UCSF MotionCor2. For CTF estimation, use CTFFIND-4.1 with pre-calculated power spectra (activate the 'Save sum of power spectra' option).");

	// Motioncor2
	char *default_location = getenv ("RELION_MOTIONCOR2_EXECUTABLE");
	char default_motioncor2[] = DEFAULTMOTIONCOR2LOCATION;
	if (default_location == NULL)
	{
		default_location = default_motioncor2;
	}

	// Common arguments RELION and UCSF implementation
	joboptions["bfactor"] = JobOption("Bfactor:", 150, 0, 1500, 50, "The B-factor that will be applied to the micrographs.");
	joboptions["patch_x"] = JobOption("Number of patches X:", std::string("1"), "Number of patches (in X and Y direction) to apply motioncor2.");
	joboptions["patch_y"] = JobOption("Number of patches Y:", std::string("1"), "Number of patches (in X and Y direction) to apply motioncor2.");
	joboptions["group_frames"] = JobOption("Group frames:", 1, 1, 5, 1, "Average together this many frames before calculating the beam-induced shifts.");
	joboptions["bin_factor"] = JobOption("Binning factor:", 1, 1, 2, 1, "Bin the micrographs this much by a windowing operation in the Fourier Tranform. Binning at this level is hard to un-do later on, but may be useful to down-scale super-resolution images. Float-values may be used. Do make sure though that the resulting micrograph size is even.");
	joboptions["fn_gain_ref"] = JobOption("Gain-reference image:", "", "*.mrc", ".", "Location of the gain-reference file to be applied to the input micrographs. Leave this empty if the movies are already gain-corrected.");
	joboptions["gain_rot"] = JobOption("Gain rotation:", job_gain_rotation_options, 0, "Rotate the gain reference by this number times 90 degrees clockwise in relion_display. This is the same as -RotGain in MotionCor2. Note that MotionCor2 uses a different convention for rotation so it says 'counter-clockwise'. Valid values are 0, 1, 2 and 3.");
	joboptions["gain_flip"] = JobOption("Gain flip:", job_gain_flip_options, 0, "Flip the gain reference after rotation. This is the same as -FlipGain in MotionCor2. 0 means do nothing, 1 means flip Y (upside down) and 2 means flip X (left to right).");

	// UCSF-wrapper
	joboptions["do_own_motioncor"] = JobOption("Use RELION's own implementation?", true ,"If set to Yes, use RELION's own implementation of a MotionCor2-like algorithm by Takanori Nakane. Otherwise, wrap to the UCSF implementation. Note that Takanori's program only runs on CPUs but uses multiple threads, while the UCSF-implementation needs a GPU but uses only one CPU thread. Takanori's implementation is most efficient when the number of frames is divisible by the number of threads (e.g. 12 or 18 threads per MPI process for 36 frames). On some machines, setting the OMP_PROC_BIND environmental variable to TRUE accelerates the program.\n\
When running on 4k x 4k movies and using 6 to 12 threads, the speeds should be similar. Note that Takanori's program uses the same model as the UCSF program and gives results that are almost identical.\n\
Whichever program you use, 'Motion Refinement' is highly recommended to get the most of your dataset.");
	joboptions["fn_motioncor2_exe"] = JobOption("MOTIONCOR2 executable:", std::string(default_location), "*.*", ".", "Location of the MOTIONCOR2 executable. You can control the default of this field by setting environment variable RELION_MOTIONCOR2_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code.");
	joboptions["fn_defect"] = JobOption("Defect file:", "", "*", ".", "Location of a UCSF MotionCor2-style defect text file or a defect map that describe the defect pixels on the detector. Each line of a defect text file should contain four numbers specifying x, y, width and height of a defect region. A defect map is an image (MRC or TIFF), where 0 means good and 1 means bad pixels. The coordinate system is the same as the input movie before application of binning, rotation and/or flipping.\nNote that the format of the defect text is DIFFERENT from the defect text produced by SerialEM! One can convert a SerialEM-style defect file into a defect map using IMOD utilities e.g. \"clip defect -D defect.txt -f tif movie.mrc defect_map.tif\". See explanations in the SerialEM manual.\n\nLeave empty if you don't have any defects, or don't want to correct for defects on your detector.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string("0"), "Provide a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':'. For example, to place one rank on device 0 and one rank on device 1, provide '0:1'.\n\
Note that multiple MotionCor2 processes should not share a GPU; otherwise, it can lead to crash or broken outputs (e.g. black images) .");
	joboptions["other_motioncor2_args"] = JobOption("Other MOTIONCOR2 arguments", std::string(""), "Additional arguments that need to be passed to MOTIONCOR2.");

	// Dose-weight
	joboptions["do_dose_weighting"] = JobOption("Do dose-weighting?", true ,"If set to Yes, the averaged micrographs will be dose-weighted.");
	joboptions["do_save_noDW"] = JobOption("Save non-dose weighted as well?", false, "Aligned but non-dose weighted images are sometimes useful in CTF estimation, although there is no difference in most cases. Whichever the choice, CTF refinement job is always done on dose-weighted particles.");
	joboptions["dose_per_frame"] = JobOption("Dose per frame (e/A2):", 1, 0, 5, 0.2, "Dose per movie frame (in electrons per squared Angstrom).");
	joboptions["pre_exposure"] = JobOption("Pre-exposure (e/A2):", 0, 0, 5, 0.5, "Pre-exposure dose (in electrons per squared Angstrom).");

	joboptions["do_save_ps"] = JobOption("Save sum of power spectra?", true, "Sum of non-dose weighted power spectra provides better signal for CTF estimation. The power spectra can be used by CTFFIND4 but not by GCTF. This option is not available for UCSF MotionCor2. You must use this option when writing in float16.");
	joboptions["group_for_ps"] = JobOption("Sum power spectra every e/A2:", 4, 0, 10, 0.5, "McMullan et al (Ultramicroscopy, 2015) sugggest summing power spectra every 4.0 e/A2 gives optimal Thon rings");
}

bool RelionJob::getCommandsMotioncorrJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);

	std::string command;
	if (joboptions["nr_mpi"].getNumber(error_message) > 1)
		command="`which relion_run_motioncorr_mpi`";
	else
		command="`which relion_run_motioncorr`";
	if (error_message != "") return false;

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
	Node node2(outputname + "corrected_micrographs.star", LABEL_MOCORR_MICS);
	outputNodes.push_back(node2);
	Node node4(outputname + "logfile.pdf", LABEL_MOCORR_LOG);
	outputNodes.push_back(node4);

	command += " --first_frame_sum " + joboptions["first_frame_sum"].getString();
	command += " --last_frame_sum " + joboptions["last_frame_sum"].getString();

	if (joboptions["do_own_motioncor"].getBoolean())
	{
		label += ".own";

		command += " --use_own ";
		command += " --j " + joboptions["nr_threads"].getString();
		if (joboptions["do_float16"].getBoolean())
		{
			if (!joboptions["do_save_ps"].getBoolean())
			{
				error_message = "When writing to float16, you have to write power spectra for CTFFIND-4.1.";
				return false;
			}

			command +=" --float16";
		}
	}
	else
	{
		label += ".motioncor2";

		command += " --use_motioncor2 ";
		command += " --motioncor2_exe " + joboptions["fn_motioncor2_exe"].getString();

		if (joboptions["do_float16"].getBoolean())
		{
			error_message = "ERROR: MotionCor2 cannot write float16 files.";
			return false;
		}

		if ((joboptions["other_motioncor2_args"].getString()).length() > 0)
			command += " --other_motioncor2_args \" " + joboptions["other_motioncor2_args"].getString() + " \"";

		// Which GPUs to use?
		command += " --gpu \"" + joboptions["gpu_ids"].getString() + "\"";
	}

	if ((joboptions["fn_defect"].getString()).length() > 0)
		command += " --defect_file " + joboptions["fn_defect"].getString();

	command += " --bin_factor " + joboptions["bin_factor"].getString();
	command += " --bfactor " + joboptions["bfactor"].getString();
	command += " --dose_per_frame " + joboptions["dose_per_frame"].getString();
	command += " --preexposure " + joboptions["pre_exposure"].getString();
	command += " --patch_x " + joboptions["patch_x"].getString();
	command += " --patch_y " + joboptions["patch_y"].getString();
	command += " --eer_grouping " + joboptions["eer_grouping"].getString();

	if (joboptions["group_frames"].getNumber(error_message) > 1.)
		command += " --group_frames " + joboptions["group_frames"].getString();
	if (error_message != "") return false;

	if ((joboptions["fn_gain_ref"].getString()).length() > 0)
	{

		int gain_rot = -1, gain_flip = -1;
		for (int i = 0; i <= 3; i++)
		{
			if (strcmp((joboptions["gain_rot"].getString()).c_str(), job_gain_rotation_options[i].c_str()) == 0)
			{
				gain_rot = i;
				break;
			}
		}

		for (int i = 0; i <= 2; i++)
		{
			if (strcmp((joboptions["gain_flip"].getString()).c_str(), job_gain_flip_options[i].c_str()) == 0)
			{
				gain_flip = i;
				break;
			}
		}

		if (gain_rot == -1 || gain_flip == -1)
			REPORT_ERROR("Illegal gain_rot and/or gain_flip.");

		command += " --gainref " + joboptions["fn_gain_ref"].getString();
		command += " --gain_rot " + integerToString(gain_rot);
		command += " --gain_flip " + integerToString(gain_flip);
	}

	if (joboptions["do_dose_weighting"].getBoolean())
	{
		command += " --dose_weighting ";
		if (joboptions["do_save_noDW"].getBoolean())
		{
			command += " --save_noDW ";
		}
	}

	if (joboptions["do_save_ps"].getBoolean())
	{
		if (!joboptions["do_own_motioncor"].getBoolean())
		{
			error_message = "'Save sum of power spectra' is not available with UCSF MotionCor2.";
			return false;
		}

		float dose_for_ps = joboptions["group_for_ps"].getNumber(error_message);
		if (error_message != "") return false;
		float dose_rate = joboptions["dose_per_frame"].getNumber(error_message);
		if (error_message != "") return false;
		if (dose_rate <= 0)
		{
			error_message = "Please specify the dose rate to calculate the grouping for power spectra.";
			return false;
		}
		if (dose_for_ps <= 0)
		{
			error_message = "Invalid dose for the grouping for power spectra.";
			return false;
		}

		int grouping_for_ps = ROUND(dose_for_ps / dose_rate);
		if (grouping_for_ps == 0)
			grouping_for_ps = 1;

		command += " --grouping_for_ps " + integerToString(grouping_for_ps) + " ";
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

	joboptions["input_star_mics"] = JobOption("Input micrographs STAR file:", NODE_MICS_CPIPE, "", "STAR files (*.star)", "A STAR file with all micrographs to run CTFFIND or Gctf on");
	joboptions["use_noDW"] = JobOption("Use micrograph without dose-weighting?", false, "If set to Yes, the CTF estimation will be done using the micrograph without dose-weighting as in rlnMicrographNameNoDW (_noDW.mrc from MotionCor2). If set to No, the normal rlnMicrographName will be used.");

	joboptions["do_phaseshift"] = JobOption("Estimate phase shifts?", false, "If set to Yes, CTFFIND4 will estimate the phase shift, e.g. as introduced by a Volta phase-plate");
	joboptions["phase_min"] = JobOption("Phase shift (deg) - Min:", std::string("0"), "Minimum, maximum and step size (in degrees) for the search of the phase shift");
	joboptions["phase_max"] = JobOption("Phase shift (deg) - Max:", std::string("180"), "Minimum, maximum and step size (in degrees) for the search of the phase shift");
	joboptions["phase_step"] = JobOption("Phase shift (deg) - Step:", std::string("10"), "Minimum, maximum and step size (in degrees) for the search of the phase shift");

	joboptions["dast"] = JobOption("Amount of astigmatism (A):", 100, 0, 2000, 100,"CTFFIND's dAst parameter, GCTFs -astm parameter");

	// CTFFIND options

	// Check for environment variable RELION_CTFFIND_EXECUTABLE
	joboptions["use_ctffind4"] = JobOption("Use CTFFIND-4.1?", false, "If set to Yes, the wrapper will use CTFFIND4 (version 4.1) for CTF estimation. This includes thread-support, calculation of Thon rings from movie frames and phase-shift estimation for phase-plate data.");
	joboptions["use_given_ps"] = JobOption("Use power spectra from MotionCorr job?", true, "If set to Yes, the CTF estimation will be done using power spectra calculated during motion correction. You must use this option if you used float16 in motion correction.");
	default_location = getenv ("RELION_CTFFIND_EXECUTABLE");
	char default_ctffind[] = DEFAULTCTFFINDLOCATION;
	if (default_location == NULL)
	{
		default_location = default_ctffind;
	}
	joboptions["fn_ctffind_exe"] = JobOption("CTFFIND-4.1 executable:", std::string(default_location), "*", ".", "Location of the CTFFIND (release 4.1 or later) executable. You can control the default of this field by setting environment variable RELION_CTFFIND_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code.");
	joboptions["slow_search"] = JobOption("Use exhaustive search?", false, "If set to Yes, CTFFIND4 will use slower but more exhaustive search. This option is recommended for CTFFIND version 4.1.8 and earlier, but probably not necessary for 4.1.10 and later. It is also worth trying this option when astigmatism and/or phase shifts are difficult to fit.");

	joboptions["box"] = JobOption("FFT box size (pix):", 512, 64, 1024, 8, "CTFFIND's Box parameter");
	joboptions["resmin"] = JobOption("Minimum resolution (A):", 30, 10, 200, 10, "CTFFIND's ResMin parameter");
	joboptions["resmax"] = JobOption("Maximum resolution (A):", 5, 1, 20, 1, "CTFFIND's ResMax parameter");
	joboptions["dfmin"] = JobOption("Minimum defocus value (A):", 5000, 0, 25000, 1000, "CTFFIND's dFMin parameter");
	joboptions["dfmax"] = JobOption("Maximum defocus value (A):", 50000, 20000, 100000, 1000, "CTFFIND's dFMax parameter");
	joboptions["dfstep"] = JobOption("Defocus step size (A):", 500, 200, 2000, 100,"CTFFIND's FStep parameter");

	joboptions["ctf_win"] = JobOption("Estimate CTF on window size (pix) ", -1, -16, 4096, 16, "If a positive value is given, a squared window of this size at the center of the micrograph will be used to estimate the CTF. This may be useful to exclude parts of the micrograph that are unsuitable for CTF estimation, e.g. the labels at the edge of phtographic film. \n \n The original micrograph will be used (i.e. this option will be ignored) if a negative value is given.");

	joboptions["use_gctf"] = JobOption("Use Gctf instead?", false, "If set to Yes, Kai Zhang's Gctf program (which runs on NVIDIA GPUs) will be used instead of Niko Grigorieff's CTFFIND4.");
	default_location = getenv("RELION_GCTF_EXECUTABLE");
	char default_gctf[] = DEFAULTGCTFLOCATION;
	if (default_location == NULL)
	{
		default_location = default_gctf;
	}
	joboptions["fn_gctf_exe"] = JobOption("Gctf executable:", std::string(default_location), "*", ".", "Location of the Gctf executable. You can control the default of this field by setting environment variable RELION_GCTF_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code.");
	joboptions["do_ignore_ctffind_params"] = JobOption("Ignore 'Searches' parameters?", true, "If set to Yes, all parameters EXCEPT for phase shift search and its ranges on the 'Searches' tab will be ignored, and Gctf's default parameters will be used (box.size=1024; min.resol=50; max.resol=4; min.defocus=500; max.defocus=90000; step.defocus=500; astigm=1000) \n \
\nIf set to No, all parameters on the CTFFIND tab will be passed to Gctf.");
	joboptions["do_EPA"] = JobOption("Perform equi-phase averaging?", false, "If set to Yes, equi-phase averaging is used in the defocus refinement, otherwise basic rotational averaging will be performed.");
	joboptions["other_gctf_args"] = JobOption("Other Gctf options:", std::string(""), "Provide additional gctf options here.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','. ");
}

bool RelionJob::getCommandsCtffindJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;

	FileName fn_outstar = outputname + "micrographs_ctf.star";
	Node node(fn_outstar, LABEL_CTFFIND_MICS);
	outputNodes.push_back(node);
	outputName = outputname;

	// PDF with histograms of the eigenvalues
	Node node3(outputname + "logfile.pdf", LABEL_CTFFIND_LOG);
	outputNodes.push_back(node3);

	if (joboptions["input_star_mics"].getString() == "")
	{
		error_message = "ERROR: empty field for input STAR file...";
		return false;
	}

	Node node2(joboptions["input_star_mics"].getString(), joboptions["input_star_mics"].node_type);
	inputNodes.push_back(node2);

	if (joboptions["nr_mpi"].getNumber(error_message) > 1)
		command="`which relion_run_ctffind_mpi`";
	else
		command="`which relion_run_ctffind`";
	if (error_message != "") return false;

	command += " --i " + joboptions["input_star_mics"].getString();
	command += " --o " + outputname;
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
		label += ".gctf";

		command += " --use_gctf --gctf_exe " + joboptions["fn_gctf_exe"].getString();
		if (joboptions["do_ignore_ctffind_params"].getBoolean())
			command += " --ignore_ctffind_params";
		if (joboptions["do_EPA"].getBoolean())
			command += " --EPA";

		// GPU-allocation
		command += " --gpu \"" + joboptions["gpu_ids"].getString() + "\"";

		if (joboptions["other_gctf_args"].getString().find("--phase_shift_H") != std::string::npos ||
	            joboptions["other_gctf_args"].getString().find("--phase_shift_L") != std::string::npos ||
		    joboptions["other_gctf_args"].getString().find("--phase_shift_S") != std::string::npos)
		{
			error_message = "Please don't specify --phase_shift_L, H, S in 'Other Gctf options'. Use 'Estimate phase shifts' and 'Phase shift - Min, Max, Step' instead.";
			return false;
		}

		if ((joboptions["other_gctf_args"].getString()).length() > 0)
			command += " --extra_gctf_options \" " + joboptions["other_gctf_args"].getString() + " \"";

	}
	else if (joboptions["use_ctffind4"].getBoolean())
	{
		label += ".ctffind4";

		command += " --ctffind_exe " + joboptions["fn_ctffind_exe"].getString();
		command += " --ctfWin " + joboptions["ctf_win"].getString();
		command += " --is_ctffind4 ";
		if (!joboptions["slow_search"].getBoolean())
		{
			command += " --fast_search ";
		}
		if (joboptions["use_given_ps"].getBoolean())
		{
			command += " --use_given_ps ";
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

	joboptions["fn_in"] = JobOption("Input micrographs:", NODE_MICS_CPIPE, "", "Input micrographs (*.{star,mrc})", "Input STAR file (with or without CTF information), OR a unix-type wildcard with all micrographs in MRC format (in this case no CTFs can be used).");

	joboptions["diameter"] = JobOption("Particle diameter (A):", 100, 0, 500, 50, "The diameter of the circle used around picked particles (in Angstroms). Only used for display." );
	joboptions["micscale"] = JobOption("Scale for micrographs:", 0.2, 0.1, 1, 0.05, "The micrographs will be displayed at this relative scale, i.e. a value of 0.5 means that only every second pixel will be displayed." );
	joboptions["sigma_contrast"] = JobOption("Sigma contrast:", 3, 0, 10, 0.5, "The micrographs will be displayed with the black value set to the average of all values MINUS this values times the standard deviation of all values in the micrograph, and the white value will be set \
to the average PLUS this value times the standard deviation. Use zero to set the minimum value in the micrograph to black, and the maximum value to white ");
	joboptions["white_val"] = JobOption("White value:", 0, 0, 512, 16, "Use non-zero values to set the value of the whitest pixel in the micrograph.");
	joboptions["black_val"] = JobOption("Black value:", 0, 0, 512, 16, "Use non-zero values to set the value of the blackest pixel in the micrograph.");

	joboptions["lowpass"] = JobOption("Lowpass filter (A)", 20, 10, 100, 5, "Lowpass filter that will be applied to the micrographs. Give a negative value to skip the lowpass filter.");
	joboptions["highpass"] = JobOption("Highpass filter (A)", -1, 100, 1000, 100, "Highpass filter that will be applied to the micrographs. This may be useful to get rid of background ramps due to uneven ice distributions. Give a negative value to skip the highpass filter. Useful values are often in the range of 200-400 Angstroms.");
	joboptions["angpix"] = JobOption("Pixel size (A)", -1, 0.3, 5, 0.1, "Pixel size in Angstroms. This will be used to calculate the filters and the particle diameter in pixels. If a CTF-containing STAR file is input, then the value given here will be ignored, and the pixel size will be calculated from the values in the STAR file. A negative value can then be given here.");
	joboptions["do_topaz_denoise"] = JobOption("OR: use Topaz denoising?", false, "If set to true, Topaz denoising will be performed instead of lowpass filtering.");
	char *default_location = getenv ("RELION_TOPAZ_EXECUTABLE");
	char default_topaz[] = DEFAULTTOPAZLOCATION;
	if (default_location == NULL)
	{
		default_location = default_topaz;
	}
	joboptions["fn_topaz_exec"] = JobOption("Topaz executable", std::string(default_location), "The location of the Topaz executable. If you need to activate conda environment, please make a wrapper shell script to do so and specify it. You can control the default of this field by setting environment variable RELION_TOPAZ_EXECUTABLE.");

	joboptions["do_startend"] = JobOption("Pick start-end coordinates helices?", false, "If set to true, start and end coordinates are picked subsequently and a line will be drawn between each pair");

	joboptions["do_fom_threshold"] = JobOption("Use autopick FOM threshold?", false, "If set to Yes, only particles with rlnAutopickFigureOfMerit values below the threshold below will be extracted.");
	joboptions["minimum_pick_fom"] = JobOption("Minimum autopick FOM: ", 0, -5, 10, 0.1, "The minimum value for the rlnAutopickFigureOfMerit for particles to be extracted.");

	joboptions["do_color"] = JobOption("Blue<>red color particles?", false, "If set to true, then the circles for each particles are coloured from red to blue (or the other way around) for a given metadatalabel. If this metadatalabel is not in the picked coordinates STAR file \
(basically only the rlnAutopickFigureOfMerit or rlnClassNumber) would be useful values there, then you may provide an additional STAR file (e.g. after classification/refinement below. Particles with values -999, or that are not in the additional STAR file will be coloured the default color: green");
	joboptions["color_label"] = JobOption("MetaDataLabel for color:", std::string("rlnAutopickFigureOfMerit"), "The Metadata label of the value to plot from red<>blue. Useful examples might be: \n \
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
	initialisePipeline(outputname, job_counter);
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

	// Allow saving, and always save default selection file upon launching the program
	FileName fn_outstar = outputname + "micrographs_selected.star";
	Node node3(fn_outstar, LABEL_MANPICK_MICS);
	outputNodes.push_back(node3);
	command += " --allow_save   --fast_save --selection " + fn_outstar;

	command += " --scale " + joboptions["micscale"].getString();
	command += " --sigma_contrast " + joboptions["sigma_contrast"].getString();
	command += " --black " + joboptions["black_val"].getString();
	command += " --white " + joboptions["white_val"].getString();

	if (joboptions["do_topaz_denoise"].getBoolean())
	{
		command += " --topaz_denoise";
		command += " --topaz_exe " + joboptions["fn_topaz_exec"].getString();
	}
	else
	{
		if (joboptions["lowpass"].getNumber(error_message) > 0.)
			command += " --lowpass " + joboptions["lowpass"].getString();
		if (error_message != "") return false;

		if (joboptions["highpass"].getNumber(error_message) > 0.)
			command += " --highpass " + joboptions["highpass"].getString();
		if (error_message != "") return false;

		if (joboptions["angpix"].getNumber(error_message) > 0.)
			command += " --angpix " + joboptions["angpix"].getString();
		if (error_message != "") return false;
	}

	if (joboptions["do_fom_threshold"].getBoolean())
	{
		command += " --minimum_pick_fom " + joboptions["minimum_pick_fom"].getString();
	}

	command += " --particle_diameter " + joboptions["diameter"].getString();

	if (joboptions["do_startend"].getBoolean())
	{
		label += ".helical";

		command += " --pick_start_end ";

		// new version: no longer save coords_suffix nodetype, but 2-column list of micrographs and coordinate files
		Node node2(outputname + "manualpick.star", LABEL_MANPICK_COORDS_HELIX);
		outputNodes.push_back(node2);
	}

	else
	{
		// new version: no longer save coords_suffix nodetype, but 2-column list of micrographs and coordinate files
		Node node2(outputname + "manualpick.star", LABEL_MANPICK_COORDS);
		outputNodes.push_back(node2);
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

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}

void RelionJob::initialiseAutopickJob()
{
	hidden_name = ".gui_autopick";

	joboptions["fn_input_autopick"] = JobOption("Input micrographs for autopick:", NODE_MICS_CPIPE, "", "Input micrographs (*.{star})", "Input STAR file (preferably with CTF information) with all micrographs to pick from.");
	joboptions["angpix"] = JobOption("Pixel size in micrographs (A)", -1, 0.3, 5, 0.1, "Pixel size in Angstroms. If a CTF-containing STAR file is input, then the value given here will be ignored, and the pixel size will be calculated from the values in the STAR file. A negative value can then be given here.");
	joboptions["continue_manual"] = JobOption("OR: continue manually?", false, "If set to Yes, an Autopick job can be continued as a manualpick job, so that incorrect picks can be corrected interactively.");

	joboptions["do_log"] = JobOption("OR: use Laplacian-of-Gaussian?", false, "If set to Yes, a Laplacian-of-Gaussian blob detection will be used (you can then leave the 'References' field empty. The preferred way to autopick is by setting this to no and providing references that were generated by 2D classification from this data set in RELION. The Laplacian-of-Gaussian method may be useful to kickstart a new data set. Please note that some options in the autopick tab are ignored in this method. See help messages of each option for details.");
	joboptions["log_diam_min"] = JobOption("Min. diameter for LoG filter (A)", 200, 50, 500, 10, "The smallest allowed diameter for the blob-detection algorithm. This should correspond to the smallest size of your particles in Angstroms.");
	joboptions["log_diam_max"] = JobOption("Max. diameter for LoG filter (A)", 250, 50, 500, 10, "The largest allowed diameter for the blob-detection algorithm. This should correspond to the largest size of your particles in Angstroms.");
	joboptions["log_invert"] = JobOption("Are the particles white?", false, "Set this option to No if the particles are black, and to Yes if the particles are white.");
	joboptions["log_maxres"] = JobOption("Maximum resolution to consider (A)", 20, 10, 100, 5, "The Laplacian-of-Gaussian filter will be applied to downscaled micrographs with the corresponding size. Give a negative value to skip downscaling.");
	joboptions["log_adjust_thr"] = JobOption("Adjust default threshold (stddev):", 0, -1., 1., 0.05, "Use this to pick more (negative number -> lower threshold) or less (positive number -> higher threshold) particles compared to the default setting. The threshold is moved this many standard deviations away from the average.");
	joboptions["log_upper_thr"] = JobOption("Upper threshold (stddev):", 999., 0., 10., 0.5, "Use this to discard picks with LoG thresholds that are this many standard deviations above the average, e.g. to avoid high contrast contamination like ice and ethane droplets. Good values depend on the contrast of micrographs and need to be interactively explored; for low contrast micrographs, values of ~ 1.5 may be reasonable, but the same value will be too low for high-contrast micrographs.");

	joboptions["do_topaz"] = JobOption("OR: use Topaz?", false, "If set to Yes, topaz will be used for autopicking. Run 2 separate jobs from the Topaz tab: one for training the model and for the actual picking.");
	char *default_location = getenv ("RELION_TOPAZ_EXECUTABLE");
	char default_topaz[] = DEFAULTTOPAZLOCATION;
	if (default_location == NULL)
	{
		default_location = default_topaz;
	}
	joboptions["fn_topaz_exec"] = JobOption("Topaz executable", std::string(default_location), "The location of the Topaz executable. If you need to activate conda environment, please make a wrapper shell script to do so and specify it. You can control the default of this field by setting environment variable RELION_TOPAZ_EXECUTABLE.");
	joboptions["do_topaz_train"] = JobOption("Perform topaz training?", false, "Set this option to Yes if you want to train a topaz model.");
	joboptions["topaz_train_picks"] = JobOption("Input picked coordinates for training:", NODE_COORDS_CPIPE, "", "Input micrographs (*.{star})", "Input STAR file (preferably with CTF information) with all micrographs to pick from.");
	joboptions["do_topaz_train_parts"] = JobOption("OR train on a set of particles? ", false, "If set to Yes, the input Coordinates above will be ignored. Instead, one uses a _data.star file from a previous 2D or 3D refinement or selection to use those particle positions for training.");
	joboptions["topaz_train_parts"] = JobOption("Particles STAR file for training: ", NODE_PARTS_CPIPE, "", "Input STAR file (*.{star})", "Filename of the STAR file with the particle coordinates to be used for training, e.g. from a previous 2D or 3D classification or selection.");
	joboptions["do_topaz_pick"] = JobOption("Perform topaz picking?", false, "Set this option to Yes if you want to use a topaz model for autopicking.");
	joboptions["topaz_particle_diameter"] = JobOption("Particle diameter (A) ", -1, 0, 2000, 20, "Diameter of the particle (to be used to infer topaz downscale factor and particle radius)");
	joboptions["topaz_nr_particles"] = JobOption("Nr of particles per micrograph: ", -1, 0, 2000, 20, "Expected average number of particles per micrograph");
	joboptions["topaz_model"] = JobOption("Trained topaz model: ", "", "SAV Files (*.sav)", ".", "Trained topaz model for topaz-based picking. Use on job for training and a next job for picking. Leave this empty to use the default (general) model.");
	joboptions["topaz_other_args"]= JobOption("Additional topaz arguments:", std::string(""), "These additional arguments will be passed onto all topaz programs.");

	joboptions["do_refs"] = JobOption("Use reference-based template-matching?", false, "If set to Yes, 2D or 3D references, as defined on the References tab will be used for autopicking.");
	joboptions["fn_refs_autopick"] = JobOption("2D references:", NODE_2DIMGS_CPIPE, "", "Input references (*.{star,mrcs})", "Input STAR file or MRC stack with the 2D references to be used for picking. Note that the absolute greyscale needs to be correct, so only use images created by RELION itself, e.g. by 2D class averaging or projecting a RELION reconstruction.");
	joboptions["do_ref3d"]= JobOption("OR: provide a 3D reference?", false, "Set this option to Yes if you want to provide a 3D map, which will be projected into multiple directions to generate 2D references.");
	joboptions["fn_ref3d_autopick"] = JobOption("3D reference:", NODE_MAP_CPIPE, "", "Input reference (*.{mrc})", "Input MRC file with the 3D reference maps, from which 2D references will be made by projection. Note that the absolute greyscale needs to be correct, so only use maps created by RELION itself from this data set.");
	joboptions["ref3d_symmetry"] = JobOption("Symmetry:", std::string("C1"), "Symmetry point group of the 3D reference. Only projections in the asymmetric part of the sphere will be generated.");
	joboptions["ref3d_sampling"] = JobOption("3D angular sampling:", job_sampling_options, 0, "There are only a few discrete \
angular samplings possible because we use the HealPix library to generate the sampling of the first two Euler angles on the sphere. \
The samplings are approximate numbers and vary slightly over the sphere.\n\n For autopicking, 30 degrees is usually fine enough, but for highly symmetrical objects one may need to go finer to adequately sample the asymmetric part of the sphere.");

	joboptions["lowpass"] = JobOption("Lowpass filter references (A)", 20, 10, 100, 5, "Lowpass filter that will be applied to the references before template matching. Do NOT use very high-resolution templates to search your micrographs. The signal will be too weak at high resolution anyway, and you may find Einstein from noise.... Give a negative value to skip the lowpass filter.");
	joboptions["highpass"] = JobOption("Highpass filter (A)", -1, 100, 1000, 100, "Highpass filter that will be applied to the micrographs. This may be useful to get rid of background ramps due to uneven ice distributions. Give a negative value to skip the highpass filter.  Useful values are often in the range of 200-400 Angstroms.");
	joboptions["angpix_ref"] = JobOption("Pixel size in references (A)", -1, 0.3, 5, 0.1, "Pixel size in Angstroms for the provided reference images. This will be used to calculate the filters and the particle diameter in pixels. If a negative value is given here, the pixel size in the references will be assumed to be the same as the one in the micrographs, i.e. the particles that were used to make the references were not rescaled upon extraction.");
	joboptions["psi_sampling_autopick"] = JobOption("In-plane angular sampling (deg)", 5, 1, 30, 1, "Angular sampling in degrees for exhaustive searches of the in-plane rotations for all references.");

	joboptions["do_invert_refs"] = JobOption("References have inverted contrast?", true, "Set to Yes to indicate that the reference have inverted contrast with respect to the particles in the micrographs.");
	joboptions["do_ctf_autopick"] = JobOption("Are References CTF corrected?", true, "Set to Yes if the references were created with CTF-correction inside RELION. \n \n If set to Yes, the input micrographs can only be given as a STAR file, which should contain the CTF information for each micrograph.");
	joboptions["do_ignore_first_ctfpeak_autopick"] = JobOption("Ignore CTFs until first peak?", false,"Set this to Yes, only if this option was also used to generate the references.");

	joboptions["threshold_autopick"] = JobOption("Picking threshold:", 0.05, 0, 1., 0.01, "Use lower thresholds to pick more particles (and more junk probably).\
\n\nThis option is ignored in the Laplacian-of-Gaussian picker. Please use 'Adjust default threshold' in the 'Laplacian' tab instead.");
	joboptions["mindist_autopick"] = JobOption("Minimum inter-particle distance (A):", 100, 0, 1000, 20, "Particles closer together than this distance will be consider to be a single cluster. From each cluster, only one particle will be picked. \
\n\nThis option takes no effect for picking helical segments. The inter-box distance is calculated with the number of asymmetrical units and the helical rise on 'Helix' tab. This option is also ignored in the Laplacian-of-Gaussian picker. The inter-box distance is calculated from particle diameters.");
	joboptions["maxstddevnoise_autopick"] = JobOption("Maximum stddev noise:", 1.1, 0.9, 1.5, 0.02, "This is useful to prevent picking in carbon areas, or areas with big contamination features. Peaks in areas where the background standard deviation in the normalized micrographs is higher than this value will be ignored. Useful values are probably in the range 1.0 to 1.2. Set to -1 to switch off the feature to eliminate peaks due to high background standard deviations.\
\n\nThis option is ignored in the Laplacian-of-Gaussian picker.");
	joboptions["minavgnoise_autopick"] = JobOption("Minimum avg noise:", -999., -2, 0.5, 0.05, "This is useful to prevent picking in carbon areas, or areas with big contamination features. Peaks in areas where the background standard deviation in the normalized micrographs is higher than this value will be ignored. Useful values are probably in the range -0.5 to 0. Set to -999 to switch off the feature to eliminate peaks due to low average background densities.\
\n\nThis option is ignored in the Laplacian-of-Gaussian picker.");
	joboptions["do_write_fom_maps"] = JobOption("Write FOM maps?", false, "If set to Yes, intermediate probability maps will be written out, which (upon reading them back in) will speed up tremendously the optimization of the threshold and inter-particle distance parameters. However, with this option, one cannot run in parallel, as disc I/O is very heavy with this option set.");
	joboptions["do_read_fom_maps"] = JobOption("Read FOM maps?", false, "If written out previously, read the FOM maps back in and re-run the picking to quickly find the optimal threshold and inter-particle distance parameters");

	joboptions["shrink"] = JobOption("Shrink factor:", 0, 0, 1, 0.1, "This is useful to speed up the calculations, and to make them less memory-intensive. The micrographs will be downscaled (shrunk) to calculate the cross-correlations, and peak searching will be done in the downscaled FOM maps. When set to 0, the micrographs will de downscaled to the lowpass filter of the references, a value between 0 and 1 will downscale the micrographs by that factor. Note that the results will not be exactly the same when you shrink micrographs!\
\n\nIn the Laplacian-of-Gaussian picker, this option is ignored and the shrink factor always becomes 0.");
	joboptions["use_gpu"] = JobOption("Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration. The Laplacian-of-Gaussian picker does not support GPU.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':'. For example: 0:1:0:1:0:1");

	joboptions["do_pick_helical_segments"] = JobOption("Pick 2D helical segments?", false, "Set to Yes if you want to pick 2D helical segments.");
	joboptions["do_amyloid"] = JobOption("Pick amyloid segments?", false, "Set to Yes if you want to use the algorithm that was developed specifically for picking amyloids.");

	joboptions["helical_tube_outer_diameter"] = JobOption("Tube diameter (A): ", 200, 100, 1000, 10, "Outer diameter (in Angstroms) of helical tubes. \
This value should be slightly larger than the actual width of the tubes.");
	joboptions["helical_nr_asu"] = JobOption("Number of unique asymmetrical units:", 1, 1, 100, 1, "Number of unique helical asymmetrical units in each segment box. This integer should not be less than 1. The inter-box distance (pixels) = helical rise (Angstroms) * number of asymmetrical units / pixel size (Angstroms). \
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
	initialisePipeline(outputname, job_counter);

	std::string command;
	if (is_continue && joboptions["continue_manual"].getBoolean())
	{

		label += ".continuemanual";

		command="`which relion_manualpick`";

		command += " --i " + joboptions["fn_input_autopick"].getString();
		command += " --odir " + outputname;
		command += " --pickname autopick";

		Node node(joboptions["fn_input_autopick"].getString(), joboptions["fn_input_autopick"].node_type);
		inputNodes.push_back(node);

		// Output new version: no longer save coords_suffix nodetype, but 2-column list of micrographs and coordinate files
		Node node2(outputname + "autopick.star", LABEL_AUTOPICK_COORDS);
		outputNodes.push_back(node2);

		// The output micrographs selection
		FileName fn_outstar = outputname + "micrographs_selected.star";
		Node node3(fn_outstar, LABEL_AUTOPICK_MICS);
		outputNodes.push_back(node3);
		command += " --allow_save  --selection " + fn_outstar;

		// A manualpicker jobwindow for display of micrographs....
		FileName fn_job = ".gui_manualpick";
		if (exists(fn_job+"job.star") || exists(fn_job+"run.job"))
		{
			RelionJob manualpickjob;
			bool iscont = false;
			manualpickjob.read(fn_job.c_str(), iscont, true); // true means do initialise

			command += " --scale " + manualpickjob.joboptions["micscale"].getString();
			command += " --sigma_contrast " + manualpickjob.joboptions["sigma_contrast"].getString();
			command += " --black " + manualpickjob.joboptions["black_val"].getString();
			command += " --white " + manualpickjob.joboptions["white_val"].getString();

			if (manualpickjob.joboptions["do_startend"].getBoolean())
			{
				command += " --pick_start_end ";
			}
			if (manualpickjob.joboptions["do_topaz_denoise"].getBoolean())
			{
				command += " --topaz_denoise --topaz_exe " + manualpickjob.joboptions["fn_topaz_exec"].getString();
			}
			else
			{
				std::string error_message = "";
				float mylowpass = manualpickjob.joboptions["lowpass"].getNumber(error_message);
				if (mylowpass > 0.)
					command += " --lowpass " + manualpickjob.joboptions["lowpass"].getString();

				float myhighpass = manualpickjob.joboptions["highpass"].getNumber(error_message);
				if (myhighpass > 0.)
					command += " --highpass " + manualpickjob.joboptions["highpass"].getString();

				float myangpix = manualpickjob.joboptions["angpix"].getNumber(error_message);
				if (myangpix > 0.)
					command += " --angpix " + manualpickjob.joboptions["angpix"].getString();
			}

			command += " --particle_diameter " + manualpickjob.joboptions["diameter"].getString();
			if (manualpickjob.joboptions["do_fom_threshold"].getBoolean())
			{
				command += " --minimum_pick_fom " + manualpickjob.joboptions["minimum_pick_fom"].getString();
			}

			if (manualpickjob.joboptions["do_color"].getBoolean())
			{
				command += " --color_label " + manualpickjob.joboptions["color_label"].getString();
				command += " --blue " + manualpickjob.joboptions["blue_value"].getString();
				command += " --red " + manualpickjob.joboptions["red_value"].getString();
				if (manualpickjob.joboptions["fn_color"].getString().length() > 0)
					command += " --color_star " + manualpickjob.joboptions["fn_color"].getString();
			}

		}
		else
		{
			// Just use some defaults if no .gui_manualpickjob.star exists
			command += " --scale 0.25";
			command += " --sigma_contrast 3";
			command += " --lowpass 20";
			command += " --particle_diameter 100";
		}

	}
	else
	{
		// Run autopicking
		if (joboptions["nr_mpi"].getNumber(error_message) > 1)
			command="`which relion_autopick_mpi`";
		else
			command="`which relion_autopick`";
		if (error_message != "") return false;

		// Input
		int icheck = 0;
		if (joboptions["do_log"].getBoolean()) icheck++;
		if (joboptions["do_topaz"].getBoolean()) icheck++;
		if (joboptions["do_refs"].getBoolean()) icheck++;

		if ( icheck != 1)
		{
			error_message = "ERROR: On the I/O tab specify (only) one of three methods: template-matching, LoG or topaz ...";
			return false;
		}

		if (joboptions["fn_input_autopick"].getString() == "" )
		{
			error_message = "ERROR: empty field for input STAR file...";
			return false;
		}

		command += " --i " + joboptions["fn_input_autopick"].getString();
		Node node(joboptions["fn_input_autopick"].getString(), joboptions["fn_input_autopick"].node_type);
		inputNodes.push_back(node);

		if (!(joboptions["do_topaz"].getBoolean() && joboptions["do_topaz_train"].getBoolean()))
		{
			// Output new version: no longer save coords_suffix nodetype, but 2-column list of micrographs and coordinate files
			Node node3(outputname + "autopick.star", LABEL_AUTOPICK_COORDS);
			outputNodes.push_back(node3);

			// PDF with histograms of the eigenvalues
			Node node3b(outputname + "logfile.pdf", LABEL_AUTOPICK_LOG);
			outputNodes.push_back(node3b);
		}

		command += " --odir " + outputname;
		command += " --pickname autopick";

		if (joboptions["do_topaz"].getBoolean())
		{

			label += ".topaz";

			icheck = 0;
			if (joboptions["do_topaz_train"].getBoolean()) icheck++;
			if (joboptions["do_topaz_pick"].getBoolean()) icheck++;
			if ( icheck != 1)
			{
				error_message = "ERROR: On the Topaz tab specify (only) one of two methods: training or picking...";
				return false;
			}

			command += " --topaz_exe " + joboptions["fn_topaz_exec"].getString();

			if (joboptions["topaz_particle_diameter"].getNumber(error_message) > 0.)
				command += " --particle_diameter " + joboptions["topaz_particle_diameter"].getString();
			if (error_message != "") return false;

			if (joboptions["do_topaz_train"].getBoolean())
			{

				label += ".train";

				if (!joboptions["use_gpu"].getBoolean())
				{
					error_message ="ERROR: For Topaz training, specify which GPUs to use on the autopicking tab; for Topaz picking GPU usage is optional";
					return false;
				}

				command += " --topaz_train";

				if (joboptions["topaz_nr_particles"].getNumber(error_message) > 0.)
					command += " --topaz_nr_particles " + joboptions["topaz_nr_particles"].getString();
				if (error_message != "") return false;

				if (joboptions["do_topaz_train_parts"].getBoolean())
				{
					command += " --topaz_train_parts " + joboptions["topaz_train_parts"].getString();
					// Output new version: no longer save coords_suffix nodetype, but 2-column list of micrographs and coordinate files
					Node nodet(outputname + "input_training_coords.star", LABEL_COORDS_CPIPE);
					outputNodes.push_back(nodet);

				}
				else
				{
					command += " --topaz_train_picks " + joboptions["topaz_train_picks"].getString();
				}

			}
			else if (joboptions["do_topaz_pick"].getBoolean())
			{
				label += ".pick";

				command += " --topaz_extract";
				if (joboptions["topaz_model"].getString() != "")
					command += " --topaz_model " + joboptions["topaz_model"].getString();
			}

			if ((joboptions["topaz_other_args"].getString()).length() > 0)
				command += " --topaz_args \" " + joboptions["topaz_other_args"].getString() + " \"";

			// GPU-stuff
			if (joboptions["use_gpu"].getBoolean())
			{
				command += " --gpu \"" + joboptions["gpu_ids"].getString() + "\"";
			}

		}
		else if (joboptions["do_log"].getBoolean())
		{
			if (joboptions["use_gpu"].getBoolean())
			{
				error_message ="ERROR: The Laplacian-of-Gaussian picker does not support GPU.";
				return false;
			}

			label += ".log";

			command += " --LoG ";
			command += " --LoG_diam_min " + joboptions["log_diam_min"].getString();
			command += " --LoG_diam_max " + joboptions["log_diam_max"].getString();
			command += " --shrink 0 --lowpass " + joboptions["log_maxres"].getString();
			command += " --LoG_adjust_threshold " + joboptions["log_adjust_thr"].getString();
			if (joboptions["log_upper_thr"].getNumber(error_message) < 999.)
				command += " --LoG_upper_threshold " + joboptions["log_upper_thr"].getString();
			if (error_message != "") return false;

			if (joboptions["log_invert"].getBoolean())
				command += " --Log_invert ";
		}
		else if (joboptions["do_refs"].getBoolean())
		{
			if (joboptions["do_ref3d"].getBoolean())
			{

				if (joboptions["fn_ref3d_autopick"].getString() == "")
				{
					error_message ="ERROR: empty field for 3D reference...";
					return false;
				}

				label += ".ref3d";

				command += " --ref " + joboptions["fn_ref3d_autopick"].getString();
				Node node2(joboptions["fn_ref3d_autopick"].getString(), LABEL_MAP_CPIPE);
				inputNodes.push_back(node2);
				command += " --sym " + joboptions["ref3d_symmetry"].getString();

				// Sampling
				int ref3d_sampling = JobOption::getHealPixOrder(joboptions["ref3d_sampling"].getString());
				if (ref3d_sampling <= 0)
				{
					error_message = "Wrong choice for ref3d_sampling";
					return false;
				}

				command += " --healpix_order " + integerToString(ref3d_sampling);
			}
			else
			{
				if (joboptions["fn_refs_autopick"].getString() == "")
				{
					error_message ="ERROR: empty field for references...";
					return false;
				}

				label += ".ref2d";

				command += " --ref " + joboptions["fn_refs_autopick"].getString();
				Node node2(joboptions["fn_refs_autopick"].getString(), LABEL_2DIMGS_CPIPE);
				inputNodes.push_back(node2);
			}

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
			if (joboptions["lowpass"].getNumber(error_message) > 0.)
				command += " --lowpass " + joboptions["lowpass"].getString();
			if (error_message != "") return false;

			if (joboptions["highpass"].getNumber(error_message) > 0.)
				command += " --highpass " + joboptions["highpass"].getString();
			if (error_message != "") return false;

			if (joboptions["angpix"].getNumber(error_message) > 0.)
				command += " --angpix " + joboptions["angpix"].getString();
			if (error_message != "") return false;

			if (joboptions["angpix_ref"].getNumber(error_message) > 0.)
				command += " --angpix_ref " + joboptions["angpix_ref"].getString();
			if (error_message != "") return false;

			command += " --threshold " + joboptions["threshold_autopick"].getString();
			if (joboptions["do_pick_helical_segments"].getBoolean())
				command += " --min_distance " + floatToString(joboptions["helical_nr_asu"].getNumber(error_message) * joboptions["helical_rise"].getNumber(error_message));
			else
				command += " --min_distance " + joboptions["mindist_autopick"].getString();
			if (error_message != "") return false;

			command += " --max_stddev_noise " + joboptions["maxstddevnoise_autopick"].getString();
			if (joboptions["minavgnoise_autopick"].getNumber(error_message) > -900.)
				command += " --min_avg_noise " + joboptions["minavgnoise_autopick"].getString();
			if (error_message != "") return false;

			// Helix
			if (joboptions["do_pick_helical_segments"].getBoolean())
			{
				command += " --helix";
				if (joboptions["do_amyloid"].getBoolean())
					command += " --amyloid";
				command += " --helical_tube_outer_diameter " + joboptions["helical_tube_outer_diameter"].getString();
				command += " --helical_tube_kappa_max " + joboptions["helical_tube_kappa_max"].getString();
				command += " --helical_tube_length_min " + joboptions["helical_tube_length_min"].getString();
			}

			// GPU-stuff
			if (joboptions["use_gpu"].getBoolean())
			{
				// for the moment always use --shrink 0 with GPUs ...
				command += " --gpu \"" + joboptions["gpu_ids"].getString() + "\"";
			}

		}

		if (joboptions["do_refs"].getBoolean() || joboptions["do_log"].getBoolean())
		{

			// Although mainly for debugging, LoG-picking does have write/read_fom_maps...
			if (joboptions["do_write_fom_maps"].getBoolean())
				command += " --write_fom_maps ";

			if (joboptions["do_read_fom_maps"].getBoolean())
				command += " --read_fom_maps ";

			if (is_continue && !(joboptions["do_read_fom_maps"].getBoolean() || joboptions["do_write_fom_maps"].getBoolean()))
				command += " --only_do_unfinished ";
		}
		else if (joboptions["do_topaz"].getBoolean())
		{
			if (is_continue)
				command += " --only_do_unfinished ";
		}
	}

	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}

void RelionJob::initialiseExtractJob()
{
	hidden_name = ".gui_extract";

    joboptions["star_mics"]= JobOption("micrograph STAR file: ", NODE_MICS_CPIPE, "", "Input STAR file (*.{star})", "Filename of the STAR file that contains all micrographs from which to extract particles.");
    // TO DOL set helical option for this
    joboptions["coords_suffix"] = JobOption("Input coordinates: ", NODE_COORDS_CPIPE, "", "Input coordinates list file (*.star)", "Starfile with a 2-column list of micrograph names and corresponding coordinate filenames (in .star, .box or as 2 or 3-column free text format)");
	joboptions["do_reextract"] = JobOption("OR re-extract refined particles? ", false, "If set to Yes, the input Coordinates above will be ignored. Instead, one uses a _data.star file from a previous 2D or 3D refinement to re-extract the particles in that refinement, possibly re-centered with their refined origin offsets. This is particularly useful when going from binned to unbinned particles.");
	joboptions["fndata_reextract"] = JobOption("Refined particles STAR file: ", NODE_PARTS_CPIPE, "", "Input STAR file (*.{star})", "Filename of the STAR file with the refined particle coordinates, e.g. from a previous 2D or 3D classification or auto-refine run.");
	joboptions["do_reset_offsets"] = JobOption("Reset the refined offsets to zero? ", false, "If set to Yes, the input origin offsets will be reset to zero. This may be useful after 2D classification of helical segments, where one does not want neighbouring segments to be translated on top of each other for a subsequent 3D refinement or classification.");
	joboptions["do_recenter"] = JobOption("OR: re-center refined coordinates? ", false, "If set to Yes, the input coordinates will be re-centered according to the refined origin offsets in the provided _data.star file. The unit is pixel, not angstrom. The origin is at the center of the box, not at the corner.");
	joboptions["recenter_x"] = JobOption("Re-center on X-coordinate (in pix): ", std::string("0"), "Re-extract particles centered on this X-coordinate (in pixels in the reference)");
	joboptions["recenter_y"] = JobOption("Re-center on Y-coordinate (in pix): ", std::string("0"), "Re-extract particles centered on this Y-coordinate (in pixels in the reference)");
	joboptions["recenter_z"] = JobOption("Re-center on Z-coordinate (in pix): ", std::string("0"), "Re-extract particles centered on this Z-coordinate (in pixels in the reference)");
	joboptions["extract_size"] = JobOption("Particle box size (pix):", 128, 64, 512, 8, "Size of the extracted particles (in pixels). This should be an even number!");
	joboptions["do_invert"] = JobOption("Invert contrast?", true, "If set to Yes, the contrast in the particles will be inverted.");
	joboptions["do_float16"] = JobOption("Write output in float16?", true ,"If set to Yes, this program will write output images in float16 MRC format. This will save a factor of two in disk space compared to the default of writing in float32. Note that RELION and CCPEM will read float16 images, but other programs may not (yet) do so.");

	joboptions["do_norm"] = JobOption("Normalize particles?", true, "If set to Yes, particles will be normalized in the way RELION prefers it.");
	joboptions["bg_diameter"] = JobOption("Diameter background circle (pix): ", -1, -1, 600, 10, "Particles will be normalized to a mean value of zero and a standard-deviation of one for all pixels in the background area.\
The background area is defined as all pixels outside a circle with this given diameter in pixels (before rescaling). When specifying a negative value, a default value of 75% of the Particle box size will be used.");
	joboptions["white_dust"] = JobOption("Stddev for white dust removal: ", -1, -1, 10, 0.1, "Remove very white pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
	joboptions["black_dust"] = JobOption("Stddev for black dust removal: ", -1, -1, 10, 0.1, "Remove very black pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
	joboptions["do_rescale"] = JobOption("Rescale particles?", false, "If set to Yes, particles will be re-scaled. Note that the particle diameter below will be in the down-scaled images.");
	joboptions["rescale"] = JobOption("Re-scaled size (pixels): ", 128, 64, 512, 8, "The re-scaled value needs to be an even number");
	joboptions["do_fom_threshold"] = JobOption("Use autopick FOM threshold?", false, "If set to Yes, only particles with rlnAutopickFigureOfMerit values below the threshold below will be extracted.");
	joboptions["minimum_pick_fom"] = JobOption("Minimum autopick FOM: ", 0, -5, 10, 0.1, "The minimum value for the rlnAutopickFigureOfMerit for particles to be extracted.");

	joboptions["do_extract_helix"] = JobOption("Extract helical segments?", false, "Set to Yes if you want to extract helical segments. RELION (.star), EMAN2 (.box) and XIMDISP (.coords) formats of tube or segment coordinates are supported.");
	joboptions["helical_tube_outer_diameter"] = JobOption("Tube diameter (A): ", 200, 100, 1000, 10, "Outer diameter (in Angstroms) of helical tubes. \
This value should be slightly larger than the actual width of helical tubes.");
	joboptions["helical_bimodal_angular_priors"] = JobOption("Use bimodal angular priors?", true, "Normally it should be set to Yes and bimodal angular priors will be applied in the following classification and refinement jobs. \
Set to No if the 3D helix looks the same when rotated upside down.");
	joboptions["do_extract_helical_tubes"] = JobOption("Coordinates are start-end only?", true, "Set to Yes if you want to extract helical segments from manually picked tube coordinates (starting and end points of helical tubes in RELION, EMAN or XIMDISP format). \
Set to No if segment coordinates (RELION auto-picked results or EMAN / XIMDISP segments) are provided.");
	joboptions["do_cut_into_segments"] = JobOption("Cut helical tubes into segments?", true, "Set to Yes if you want to extract multiple helical segments with a fixed inter-box distance. \
If it is set to No, only one box at the center of each helical tube will be extracted.");
	joboptions["helical_nr_asu"] = JobOption("Number of unique asymmetrical units:", 1, 1, 100, 1, "Number of unique helical asymmetrical units in each segment box. This integer should not be less than 1. The inter-box distance (pixels) = helical rise (Angstroms) * number of asymmetrical units / pixel size (Angstroms). \
The optimal inter-box distance might also depend on the box size, the helical rise and the flexibility of the structure. In general, an inter-box distance of ~10% * the box size seems appropriate.");
	joboptions["helical_rise"] = JobOption("Helical rise (A):", 1, 0, 100, 0.01, "Helical rise in Angstroms. (Please click '?' next to the option above for details about how the inter-box distance is calculated.)");

}

bool RelionJob::getCommandsExtractJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;
	if (joboptions["nr_mpi"].getNumber(error_message) > 1)
		command="`which relion_preprocess_mpi`";
	else
		command="`which relion_preprocess`";
	if (error_message != "") return false;

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

		if (joboptions["do_reset_offsets"].getBoolean() && joboptions["do_recenter"].getBoolean())
		{
			error_message = "ERROR: you cannot both reset refined offsets and recenter on refined coordinates, choose one...";
			return false;
		}

		label += ".reextract";

		command += " --reextract_data_star " + joboptions["fndata_reextract"].getString();
		Node node2(joboptions["fndata_reextract"].getString(), joboptions["fndata_reextract"].node_type);
		inputNodes.push_back(node2);
		if (joboptions["do_reset_offsets"].getBoolean())
		{
			command += " --reset_offsets";
		}
		else if (joboptions["do_recenter"].getBoolean())
		{
			command += " --recenter";
			command += " --recenter_x " + joboptions["recenter_x"].getString();
			command += " --recenter_y " + joboptions["recenter_y"].getString();
			command += " --recenter_z " + joboptions["recenter_z"].getString();
		}
	}
	else
	{
		FileName mylist = joboptions["coords_suffix"].getString();
		if (mylist == "")
		{
			error_message = "ERROR: empty field for coordinate STAR file...";
			return false;
		}
		// Attempt at backwards compatibility
		if (mylist.contains("coords_suffix"))
		{
			command += " --coord_dir " + mylist.beforeLastOf("/") + "/";
			command += " --coord_suffix " + (mylist.afterLastOf("/")).without("coords_suffix");
		}
		else
		{
			command += " --coord_list " + mylist;
		}
		Node node2(mylist, joboptions["coords_suffix"].node_type);
		inputNodes.push_back(node2);
	}

	// Output
	FileName fn_ostar = outputname + "particles.star";

	command += " --part_star " + fn_ostar;

	if (joboptions["do_reextract"].getBoolean())
	{
		FileName fn_pickstar = outputname + "extractpick.star";
		Node node(fn_pickstar, LABEL_EXTRACT_COORDS_REEX);
		outputNodes.push_back(node);
		command += " --pick_star " + fn_pickstar;
	}

	if (joboptions["do_extract_helix"].getBoolean() && joboptions["do_extract_helical_tubes"].getBoolean())
	{
		FileName fn_pickstar = outputname + "extractpick.star";
		Node node(fn_pickstar, LABEL_EXTRACT_COORDS_HELIX);
		outputNodes.push_back(node);
		command += " --pick_star " + fn_pickstar;
	}


	command += " --part_dir " + outputname;
	command += " --extract";
	command += " --extract_size " + joboptions["extract_size"].getString();

	if (joboptions["do_fom_threshold"].getBoolean())
	{
		command += " --minimum_pick_fom " + joboptions["minimum_pick_fom"].getString();
	}

	if (joboptions["do_float16"].getBoolean())
	{
		command += " --float16 ";
	}

	// Operate stuff
	// Get an integer number for the bg_radius
	RFLOAT bg_radius = (joboptions["bg_diameter"].getNumber(error_message) < 0.) ? 0.75 * joboptions["extract_size"].getNumber(error_message) : joboptions["bg_diameter"].getNumber(error_message);
	if (error_message != "") return false;

	bg_radius /= 2.; // Go from diameter to radius
	if (joboptions["do_rescale"].getBoolean())
	{
		command += " --scale " + joboptions["rescale"].getString();
		bg_radius *= joboptions["rescale"].getNumber(error_message);
		if (error_message != "") return false;

		bg_radius /= joboptions["extract_size"].getNumber(error_message);
		if (error_message != "") return false;
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

	// Helix
	if (joboptions["do_extract_helix"].getBoolean())
	{
		Node node3(fn_ostar, LABEL_EXTRACT_PARTS_HELIX);
		outputNodes.push_back(node3);

		label += ".helical";

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
	else
	{
		Node node3(fn_ostar, LABEL_EXTRACT_PARTS);
		outputNodes.push_back(node3);
	}


	if (is_continue)
		command += " --only_do_unfinished ";


	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();
	commands.push_back(command);


	if (joboptions["do_reextract"].getBoolean())
	{
		Node node(outputname + "reextract.star", LABEL_EXTRACT_COORDS_REEX);
		outputNodes.push_back(node);
	}

	if (joboptions["do_extract_helix"].getBoolean() && joboptions["do_extract_helical_tubes"].getBoolean())
	{
		Node node(outputname + "helix_segments.star", LABEL_EXTRACT_COORDS_HELIX);
		outputNodes.push_back(node);
	}


	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}

void RelionJob::initialiseSelectJob()
{
	hidden_name = ".gui_select";

	joboptions["fn_model"] = JobOption("Select classes from job:", NODE_OPTIMISER_CPIPE, "", "STAR files (*_optimiser.star)", "A _optimiser.star (or for backwards compatibility also a _model.star) file from a previous 2D or 3D classification run to select classes from.");
	joboptions["fn_mic"] = JobOption("OR select from micrographs.star:", NODE_MICS_CPIPE, "", "STAR files (*.star)", "A micrographs.star file to select micrographs from.");
	joboptions["fn_data"] = JobOption("OR select from particles.star:", NODE_PARTS_CPIPE, "", "STAR files (*.star)", "A particles.star file to select individual particles from.");

	joboptions["do_class_ranker"] = JobOption("Automatically select 2D classes?", false, "If set to True, the class_ranker program will be used to make an automated class selection, based on the parameters below. This option only works when selecting classes from a relion_refine job (input optimiser.star on the I.O tab)");
	joboptions["rank_threshold"] = JobOption("Minimum threshold for auto-selection: ", 0.5, 0, 1, 0.05, "Only classes with a pre dicted threshold above this value will be selected.");
	joboptions["select_nr_parts"] = JobOption("Select at least this many particles: ", -1, -1, 10000, 500, "Even if they have scores below the minimum threshold, select at least this many particles with the best scores.");
	joboptions["select_nr_classes"] = JobOption("OR: select at least this many classes: ", -1, -1, 24, 1, "Even if they have scores below the minimum threshold, select at least this many classes with the best scores.");

    char *default_location;
    default_location = getenv("RELION_PYTHON_EXECUTABLE");
    char default_python[] = DEFAULTPYTHONLOCATION;
    if (default_location == NULL)
    {
        default_location = default_python;
    }
	joboptions["python_exe"] = JobOption("Python executable: ", std::string(default_location), "This version of python should include torch and numpy. We have found that the one from topaz (which is also used for auto-picking) works well. At the LMB, it is here: /public/EM/anaconda3/envs/topaz/bin/python");

	joboptions["do_recenter"] = JobOption("Re-center the class averages?", true, "This option is only used when selecting particles from 2D classes. The selected class averages will all re-centered on their center-of-mass. This is useful when you plane to use these class averages as templates for auto-picking.");
	joboptions["do_regroup"] = JobOption("Regroup the particles?", false, "If set to Yes, then the program will regroup the selected particles in 'more-or-less' the number of groups indicated below. For re-grouping from individual particle _data.star files, a _model.star file with the same prefix should exist, i.e. the particle star file should be generated by relion_refine");
	joboptions["nr_groups"] = JobOption("Approximate nr of groups: ", 1, 50, 20, 1, "It is normal that the actual number of groups may deviate a little from this number. ");

	joboptions["do_select_values"] = JobOption("Select based on metadata values?", false, "If set to Yes, the job will be non-interactive and the selected star file will be based only on the value of the corresponding metadata label. Note that this option is only valid for micrographs or particles STAR files.");
	joboptions["select_label"] = JobOption("Metadata label for subset selection:", (std::string)"rlnCtfMaxResolution", "This column from the input STAR file will be used for the subset selection.");
	joboptions["select_minval"] = JobOption("Minimum metadata value:",  (std::string)"-9999.", "Only lines in the input STAR file with the corresponding metadata value larger than or equal to this value will be included in the subset.");
	joboptions["select_maxval"] = JobOption("Maximum metadata value:",  (std::string)"9999.", "Only lines in the input STAR file with the corresponding metadata value smaller than or equal to this value will be included in the subset.");

	joboptions["do_discard"] = JobOption("OR: select on image statistics?", false, "If set to Yes, the job will be non-interactive and all images in the input star file that have average and/or stddev pixel values that are more than the specified sigma-values away from the ensemble mean will be discarded.");
	joboptions["discard_label"] = JobOption("Metadata label for images:", (std::string)"rlnImageName", "Specify which column from the input STAR contains the names of the images to be used to calculate the average and stddev values.");
	joboptions["discard_sigma"] = JobOption("Sigma-value for discarding images:", 4, 1, 10, 0.1, "Images with average and/or stddev values that are more than this many times the ensemble stddev away from the ensemble mean will be discarded.");

	joboptions["do_split"] = JobOption("OR: split into subsets?", false, "If set to Yes, the job will be non-interactive and the star file will be split into subsets as defined below.");
	joboptions["do_random"] = JobOption("Randomise order before making subsets?:", false, "If set to Yes, the input STAR file order will be randomised. If set to No, the original order in the input STAR file will be maintained.");
	joboptions["split_size"] = JobOption("Subset size: ", 100, 100, 10000, 100, "The number of lines in each of the output subsets. When this is -1, items are divided into a number of subsets specified in the next option.");
	joboptions["nr_split"] = JobOption("OR: number of subsets: ", -1, 1, 50, 1, "Give a positive integer to specify into how many equal-sized subsets the data will be divided. When the subset size is also specified, only this number of subsets, each with the specified size, will be written, possibly missing some items. When this is -1, all items are used, generating as many subsets as necessary.");

	joboptions["do_remove_duplicates"] = JobOption("OR: remove duplicates?", false, "If set to Yes, duplicated particles that are within a given distance are removed leaving only one. Duplicated particles are sometimes generated when particles drift into the same position during alignment. They inflate and invalidate gold-standard FSC calculation.");
	joboptions["duplicate_threshold"] = JobOption("Minimum inter-particle distance (A)", 30, 0, 1000, 1, "Particles within this distance are removed leaving only one.");
	joboptions["image_angpix"] = JobOption("Pixel size before extraction (A)", -1, -1, 10, 0.01, "The pixel size of particles (relevant to rlnOriginX/Y) is read from the STAR file. When the pixel size of the original micrograph used for auto-picking and extraction (relevant to rlnCoordinateX/Y) is different, specify it here. In other words, this is the pixel size after binning during motion correction, but before down-sampling during extraction.");
}

bool RelionJob::getCommandsSelectJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;

	if (joboptions["fn_model"].getString() == "" &&
			joboptions["fn_mic"].getString() == "" && joboptions["fn_data"].getString() == "")
	{
		// Nothing was selected...
		error_message = "Please select an input file.";
		return false;
	}

	int c = 0;
	if (joboptions["do_select_values"].getBoolean()) c++;
	if (joboptions["do_discard"].getBoolean()) c++;
	if (joboptions["do_split"].getBoolean()) c++;
	if (joboptions["do_remove_duplicates"].getBoolean()) c++;
	if (c > 1)
	{
		error_message = "You cannot do many tasks simultaneously...";
		return false;
	}

	if (joboptions["do_remove_duplicates"].getBoolean())
	{

		label += ".removeduplicates";

		// Remove duplicates
		command="`which relion_star_handler`";

		if (joboptions["fn_mic"].getString() != "" || joboptions["fn_model"].getString() != "")
		{
			error_message = "ERROR: Duplicate removal is only possible for particle STAR files...";
			return false;
		}

		if (joboptions["fn_data"].getString() == "")
		{
			error_message = "ERROR: Duplicate removal needs a particle STAR file...";
			return false;
		}

		Node node(joboptions["fn_data"].getString(), joboptions["fn_data"].node_type);
		inputNodes.push_back(node);
		command += " --i " + joboptions["fn_data"].getString();

		FileName fn_out = outputname+"particles.star";
		Node node2(fn_out, LABEL_SELECT_PARTS);
		outputNodes.push_back(node2);
		command += " --o " + fn_out;

		command += " --remove_duplicates " + joboptions["duplicate_threshold"].getString();
		if (joboptions["image_angpix"].getNumber(error_message) > 0)
			command += " --image_angpix " + joboptions["image_angpix"].getString();
		if (error_message != "") return false;

	}
	else if (joboptions["do_select_values"].getBoolean() || joboptions["do_discard"].getBoolean() || joboptions["do_split"].getBoolean())
	{
		// Value-based selection
		command="`which relion_star_handler`";

		if (joboptions["fn_model"].getString() != "")
		{
			error_message = "ERROR: Value-selection or subset splitting is only possible for micrograph or particle STAR files...";
			return false;
		}

		FileName fn_out;
		if (joboptions["fn_mic"].getString() != "")
		{
			Node node(joboptions["fn_mic"].getString(), joboptions["fn_mic"].node_type);
			inputNodes.push_back(node);
			command += " --i " + joboptions["fn_mic"].getString();
			fn_out = outputname+"micrographs.star";

		}
		else if (joboptions["fn_data"].getString() != "")
		{
			Node node(joboptions["fn_data"].getString(), joboptions["fn_data"].node_type);
			inputNodes.push_back(node);
			command += " --i " + joboptions["fn_data"].getString();
			fn_out = outputname+"particles.star";
		}
		command += " --o " + fn_out;

		if (joboptions["do_select_values"].getBoolean() || joboptions["do_discard"].getBoolean())
		{

			if (joboptions["fn_mic"].getString() != "")
			{
				Node node2(fn_out, LABEL_SELECT_MICS);
				outputNodes.push_back(node2);
			}
			else if (joboptions["fn_data"].getString() != "")
			{
				Node node2(fn_out, LABEL_SELECT_PARTS);
				outputNodes.push_back(node2);
			}

			if (joboptions["do_select_values"].getBoolean())
			{
				label += ".onvalue";

				command += " --select " + joboptions["select_label"].getString();
				command += " --minval " + joboptions["select_minval"].getString();
				command += " --maxval " + joboptions["select_maxval"].getString();
			}
			else if (joboptions["do_discard"].getBoolean())
			{
				label += ".discard";

				command += " --discard_on_stats ";
				command += " --discard_label " + joboptions["discard_label"].getString();
				command += " --discard_sigma " + joboptions["discard_sigma"].getString();
			}

		}
		else if (joboptions["do_split"].getBoolean())
		{

			label += ".split";

			int nr_split=0;
			command += " --split ";
			if (joboptions["do_random"].getBoolean())
			{
				command += " --random_order ";
			}

			if (joboptions["nr_split"].getNumber(error_message) <= 0 && joboptions["split_size"].getNumber(error_message) <= 0
					&& !joboptions["nr_split"].isSchedulerVariable() && !joboptions["split_size"].isSchedulerVariable())
			{
				error_message = "ERROR: When splitting the input STAR file into subsets, set nr_split and/or split_size to a positive value";
				return false;
			}

			if (joboptions["nr_split"].getNumber(error_message) > 0 && !joboptions["nr_split"].isSchedulerVariable())
			{
				if (error_message != "") return false;

				nr_split = joboptions["nr_split"].getNumber(error_message);
				command += " --nr_split " + joboptions["nr_split"].getString();
			}
			if (joboptions["split_size"].getNumber(error_message) > 0 && !joboptions["split_size"].isSchedulerVariable())
			{
				if (error_message != "") return false;

				command += " --size_split " + joboptions["split_size"].getString();
			}

			// As of relion-3.1, star_handler will write out a star file with the output nodes, which will be read by the pipeliner
		}
	}
	else
	{

		// Automated 2D class selection through the class_ranker
		if (joboptions["do_class_ranker"].getBoolean())
		{

			label += ".class2dauto";

			if (joboptions["fn_model"].getString() == "")
			{
				error_message = "ERROR: When using automatically selecting 2D classes, one needs to provide an optimiser.star file";
				return false;
			}

			if (joboptions["do_regroup"].getBoolean() || joboptions["do_recenter"].getBoolean())
			{
				error_message = "ERROR: regrouping and recentering have not been implemented in class_ranker.";
				return false;
			}

			command = "`which relion_class_ranker`";

			// input
			command += " --opt " + joboptions["fn_model"].getString();
			Node node(joboptions["fn_model"].getString(), joboptions["fn_model"].node_type);
			inputNodes.push_back(node);

			//output
			command += " --o " + outputname + " --fn_sel_parts particles.star --fn_sel_classavgs class_averages.star";

			command += " --python " + joboptions["python_exe"].getString();

			if (joboptions["select_nr_parts"].getNumber(error_message) > 0)
			{
				command += " --select_min_nr_particles " + joboptions["select_nr_parts"].getString();
			}
			else if (joboptions["select_nr_classes"].getNumber(error_message) > 0)
			{
				command += " --select_min_nr_classes " + joboptions["select_nr_classes"].getString();
			}

			FileName fn_parts = outputname+"particles.star";
			Node node2(fn_parts, LABEL_SELECT_PARTS);
			outputNodes.push_back(node2);

			FileName fn_imgs = outputname+"class_averages.star";
			Node node3(fn_imgs, LABEL_SELECT_CLAVS);
			outputNodes.push_back(node3);

			// Also save optimiser.star, which could be used for next manual selection (but ordered for examples on the new scores)
			command += " --fn_root rank";

			// Only save the 2D class averages for 2D jobs
			FileName fn_opt = outputname+"rank_optimiser.star";
			Node node4(fn_opt, LABEL_SELECT_OPT);
			outputNodes.push_back(node4);

			// perform the actual prediction and selection
			command += " --do_granularity_features ";
			command += " --auto_select ";
			command += " --min_score " + joboptions["rank_threshold"].getString();
		}
		else
		{

			// Interactive selection
			label += ".interactive";

			command="`which relion_display`";

			// I/O
			if (joboptions["fn_model"].getString() != "")
			{

				command += " --gui --i " + joboptions["fn_model"].getString();
				Node node(joboptions["fn_model"].getString(), joboptions["fn_model"].node_type);
				inputNodes.push_back(node);

				FileName fn_parts = outputname+"particles.star";
				command += " --allow_save --fn_parts " + fn_parts;
				Node node2(fn_parts, LABEL_SELECT_PARTS);
				outputNodes.push_back(node2);

				// Only save the 2D class averages for 2D jobs
				FileName fnt = joboptions["fn_model"].getString();
				if (fnt.contains("Class2D/"))
				{
					FileName fn_imgs = outputname+"class_averages.star";
					command += " --fn_imgs " + fn_imgs;
					Node node3(fn_imgs, LABEL_SELECT_CLAVS);
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
				Node node2(fn_mics, LABEL_SELECT_MICS);
				outputNodes.push_back(node2);
			}
			else if (joboptions["fn_data"].getString() != "")
			{
				command += " --gui --i " + joboptions["fn_data"].getString();
				Node node(joboptions["fn_data"].getString(), joboptions["fn_data"].node_type);
				inputNodes.push_back(node);

				FileName fn_parts = outputname+"particles.star";
				command += " --allow_save --fn_imgs " + fn_parts;
				Node node2(fn_parts, LABEL_SELECT_PARTS);
				outputNodes.push_back(node2);
			}
		}
	}

	// Re-grouping
	if (joboptions["do_regroup"].getBoolean())
	{
		if (joboptions["fn_model"].getString() == "")
		{
			error_message = "Re-grouping only works for model.star/optimiser.star files...";
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

	joboptions["fn_img"] = JobOption("Input images STAR file:", NODE_PARTS_CPIPE, "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
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


	joboptions["do_em"] = JobOption("Use EM algorithm?", false, "If set to Yes, the slower expectation-maximization algorithm will be used. This was the default option in releases prior to 4.0-beta. If set to No, then one needs to use the (faster) VDAM (variable metric gradient descent with adaptive moments) algorithm below. will be used.");
	joboptions["nr_iter_em"] = JobOption("Number of EM iterations:", 25, 1, 50, 1, "Number of EM iterations to be performed. \
Note that the current implementation of 2D class averaging and 3D classification does NOT comprise a convergence criterium. \
Therefore, the calculations will need to be stopped by the user if further iterations do not yield improvements in resolution or classes. \n\n \
Also note that upon restarting, the iteration number continues to be increased, starting from the final iteration in the previous run. \
The number given here is the TOTAL number of iterations. For example, if 10 iterations have been performed previously and one restarts to perform \
an additional 5 iterations (for example with a finer angular sampling), then the number given here should be 10+5=15.");


	joboptions["do_grad"] = JobOption("Use VDAM algorithm?", true, "If set to Yes, the faster VDAM algorithm will be used. This algorithm was introduced with relion-4.0. If set to No, then the slower EM algorithm needs to be used.");
	joboptions["nr_iter_grad"] = JobOption("Number of VDAM mini-batches:", 200, 50, 500, 10, "Number of mini-batches to be processed using the VDAM algorithm. Using 200 has given good results for many data sets. Using 100 will run faster, at the expense of some quality in the results.");

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
	joboptions["do_center"] = JobOption("Center class averages?", true, "If set to Yes, every iteration the class average images will be centered on their center-of-mass. This will only work for positive signals, so the particles should be white.");

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
	joboptions["allow_coarser"] = JobOption("Allow coarser sampling?", false, "If set to Yes, the program will use coarser angular and translational samplings if the estimated accuracies of the assignments is still low in the earlier iterations. This may speed up the calculations.");

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
	joboptions["do_restrict_xoff"] = JobOption("Restrict helical offsets to rise:", true, "Set to Yes if you want to restrict the translational offsets along the helices to the rise of the helix given below. Set to No to allow free (conventional) translational offsets.");
	joboptions["helical_rise"] = JobOption("Helical rise (A):", 4.75, -1, 100, 1, "The helical rise (in Angstroms). Translational offsets along the helical axis will be limited from -rise/2 to +rise/2, with a flat prior.");


	joboptions["nr_pool"] = JobOption("Number of pooled particles:", 3, 1, 16, 1, "Particles are processed in individual batches by MPI followers. During each batch, a stack of particle images is only opened and closed once to improve disk access times. \
All particle images of a single batch are read into memory together. The size of these batches is at least one particle per thread used. The nr_pooled_particles parameter controls how many particles are read together for each thread. If it is set to 3 and one uses 8 threads, batches of 3x8=24 particles will be read together. \
This may improve performance on systems where disk access, and particularly metadata handling of disk access, is a problem. It has a modest cost of increased RAM usage.");
	joboptions["do_parallel_discio"] = JobOption("Use parallel disc I/O?", true, "If set to Yes, all MPI followers will read images from disc. \
Otherwise, only the leader will read images and send them through the network to the followers. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many followers reading in parallel. If your datasets contain particles with different box sizes, you have to say Yes.");
	joboptions["do_preread_images"] = JobOption("Pre-read all particles into RAM?", false, "If set to Yes, all particle images will be read into computer memory, which will greatly speed up calculations on systems with slow disk access. However, one should of course be careful with the amount of RAM available. \
Because particles are read in float-precision, it will take ( N * box_size * box_size * 4 / (1024 * 1024 * 1024) ) Giga-bytes to read N particles into RAM. For 100 thousand 200x200 images, that becomes 15Gb, or 60 Gb for the same number of 400x400 particles. \
Remember that running a single MPI follower on each node that runs as many threads as available cores will have access to all available RAM. \n \n If parallel disc I/O is set to No, then only the leader reads all particles into RAM and sends those particles through the network to the MPI followers during the refinement iterations.");
	const char *default_scratch = getenv("RELION_SCRATCH_DIR");
	if (default_scratch == NULL)
	{
		default_scratch = DEFAULTSCRATCHDIR;
	}
	joboptions["scratch_dir"] = JobOption("Copy particles to scratch directory:", std::string(default_scratch), "If a directory is provided here, then the job will create a sub-directory in it called relion_volatile. If that relion_volatile directory already exists, it will be wiped. Then, the program will copy all input particles into a large stack inside the relion_volatile subdirectory. \
Provided this directory is on a fast local drive (e.g. an SSD drive), processing in all the iterations will be faster. If the job finishes correctly, the relion_volatile directory will be wiped. If the job crashes, you may want to remove it yourself.");
	joboptions["do_combine_thru_disc"] = JobOption("Combine iterations through disc?", false, "If set to Yes, at the end of every iteration all MPI followers will write out a large file with their accumulated results. The MPI leader will read in all these files, combine them all, and write out a new file with the combined results. \
All MPI salves will then read in the combined results. This reduces heavy load on the network, but increases load on the disc I/O. \
This will affect the time it takes between the progress-bar in the expectation step reaching its end (the mouse gets to the cheese) and the start of the ensuing maximisation step. It will depend on your system setup which is most efficient.");

	joboptions["use_gpu"] = JobOption("Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','. For example: '0,0:1,1:0,0:1,1'");
}

bool RelionJob::getCommandsClass2DJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;

	if (joboptions["nr_mpi"].getNumber(error_message) > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";
	if (error_message != "") return false;

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
		// SHWS 10dec2020: switch off using run_ctXX output for continue jobs, as this will affect Schedulers
		//int it = (int)textToFloat((joboptions["fn_cont"].getString().substr(pos_it+3, 6)).c_str());
		//fn_run += "_ct" + floatToString(it);
		command += " --continue " + joboptions["fn_cont"].getString();
	}

	command += " --o " + outputname + fn_run;

	int my_classes = (int)joboptions["nr_classes"].getNumber(error_message);
	if (error_message != "") return false;

	// Optimisation
	int my_iter;
	if (joboptions["do_em"].getBoolean())
        {

            if (joboptions["do_grad"].getBoolean())
            {
                error_message = "You cannot specify to use both the EM and the VDAM algorithm!";
                return false;
            }

            command += " --iter " + joboptions["nr_iter_em"].getString();

            my_iter = (int)joboptions["nr_iter_em"].getNumber(error_message);
            if (error_message != "") return false;
        }
        else if (joboptions["do_grad"].getBoolean())
	{
            if (joboptions["nr_mpi"].getNumber(error_message) > 1)
            {
                error_message = "Gradient refinement (running the VDAM algorithm) is not supported together with MPI.";
                return false;
	    }

            command += " --grad --class_inactivity_threshold 0.1 --grad_write_iter 10";
            command += " --iter " + joboptions["nr_iter_grad"].getString();

            my_iter = (int)joboptions["nr_iter_grad"].getNumber(error_message);
            if (error_message != "") return false;
	}
        else
        {
            error_message = "You need to specify to use either the EM or the VDAM algorithm";
            return false;
        }

	outputNodes = getOutputNodesRefine(outputname + fn_run, my_iter, my_classes, 2, 1, "Class2D");

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
	// Takanori observed bad 2D classifications with pad1, so use pad2 always. Memory isnt a problem here anyway.
	command += " --pad 2 ";

	// CTF stuff
	if (!is_continue)
	{
		if (joboptions["do_ctf_correction"].getBoolean())
		{
			command += " --ctf ";
			if (joboptions["ctf_intact_first_peak"].getBoolean())
				command += " --ctf_intact_first_peak ";
		}
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
		if (joboptions["highres_limit"].getNumber(error_message) > 0)
			command += " --strict_highres_exp " + joboptions["highres_limit"].getString();
		if (error_message != "") return false;

	}

	if (joboptions["do_center"].getBoolean())
	{
		command += " --center_classes ";
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
		command += " --psi_step " + floatToString(joboptions["psi_sampling"].getNumber(error_message) * pow(2., iover));
		if (error_message != "") return false;

		// Offset range
		command += " --offset_range " + joboptions["offset_range"].getString();
		// The sampling given in the GUI will be the oversampled one!
		command += " --offset_step " + floatToString(joboptions["offset_step"].getNumber(error_message) * pow(2., iover));
		if (error_message != "") return false;

		if (joboptions["allow_coarser"].getBoolean())
		{
			command += " --allow_coarser_sampling";
		}

	}

	// Helix
	if (joboptions["do_helix"].getBoolean())
	{
		label += ".helical";

		command += " --helical_outer_diameter " + joboptions["helical_tube_outer_diameter"].getString();

		if (joboptions["dont_skip_align"].getBoolean())
		{
			if (joboptions["do_bimodal_psi"].getBoolean())
				command += " --bimodal_psi";

			RFLOAT val = joboptions["range_psi"].getNumber(error_message);
			if (error_message != "") return false;

			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_psi " + floatToString(val / 3.);

			if (joboptions["do_restrict_xoff"].getBoolean())
			{
				command += " --helix --helical_rise_initial " + joboptions["helical_rise"].getString();
			}
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

	joboptions["fn_img"] = JobOption("Input images STAR file:", NODE_PARTS_CPIPE, "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \
In Gradient optimisation, it is very important that there are particles from enough different orientations. One only needs a few thousand to 10k particles. When selecting good 2D classes in the Subset Selection jobtype, use the option to select a maximum number of particles from each class to generate more even angular distributions for SGD.\
\n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");
	joboptions["fn_cont"] = JobOption("Continue from here: ", std::string(""), "STAR Files (*_optimiser.star)", "CURRENT_ODIR", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");

	joboptions["nr_iter"] = JobOption("Number of VDAM mini-batches:", 200, 50, 500, 10, "How many iterations (i.e. mini-batches) to perform with the VDAM algorithm?");
	joboptions["tau_fudge"] = JobOption("Regularisation parameter T:", 4 , 0.1, 10, 0.1, "Bayes law strictly determines the relative weight between \
the contribution of the experimental data and the prior. However, in practice one may need to adjust this weight to put slightly more weight on \
the experimental data to allow optimal results. Values greater than 1 for this regularisation parameter (T in the JMB2011 paper) put more \
weight on the experimental data. Values around 2-4 have been observed to be useful for 3D initial model calculations");

	joboptions["nr_classes"] = JobOption("Number of classes:", 1, 1, 50, 1, "The number of classes (K) for a multi-reference ab initio SGD refinement. \
These classes will be made in an unsupervised manner, starting from a single reference in the initial iterations of the SGD, and the references will become increasingly dissimilar during the inbetween iterations.");
	joboptions["sym_name"] = JobOption("Symmetry:", std::string("C1"), "The initial model is always generated in C1 and then aligned to and symmetrized with the specified point group. If the automatic alignment fails, please manually rotate run_itNNN_class001.mrc (NNN is the number of iterations) so that it conforms the symmetry convention.");
	joboptions["do_run_C1"] = JobOption("Run in C1 and apply symmetry later? ", true, "If set to Yes, the gradient-driven optimisation is run in C1 and the symmetry orientation is searched and applied later. If set to No, the entire optimisation is run in the symmetry point group indicated above.");
	joboptions["particle_diameter"] = JobOption("Mask diameter (A):", 200, 0, 1000, 10, "The experimental images will be masked with a soft \
circular mask with this diameter. Make sure this radius is not set too small because that may mask away part of the signal! \
If set to a value larger than the image size no masking will be performed.\n\n\
The same diameter will also be used for a spherical mask of the reference structures if no user-provided mask is specified.");
	joboptions["do_solvent"] = JobOption("Flatten and enforce non-negative solvent?", true, "If set to Yes, the job will apply a spherical mask and enforce all values in the reference to be non-negative.");

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
	joboptions["ctf_intact_first_peak"] = JobOption("Ignore CTFs until first peak?", false, "If set to Yes, then CTF-amplitude correction will \
only be performed from the first peak of each CTF onward. This can be useful if the CTF model is inadequate at the lowest resolution. \
Still, in general using higher amplitude contrast on the CTFs (e.g. 10-20%) often yields better results. \
Therefore, this option is not generally recommended: try increasing amplitude contrast (in your input STAR file) first!");

	joboptions["do_parallel_discio"] = JobOption("Use parallel disc I/O?", true, "If set to Yes, all MPI followers will read their own images from disc. \
Otherwise, only the leader will read images and send them through the network to the followers. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many followers reading in parallel. If your datasets contain particles with different box sizes, you have to say Yes.");
	joboptions["nr_pool"] = JobOption("Number of pooled particles:", 3, 1, 16, 1, "Particles are processed in individual batches by MPI followers. During each batch, a stack of particle images is only opened and closed once to improve disk access times. \
All particle images of a single batch are read into memory together. The size of these batches is at least one particle per thread used. The nr_pooled_particles parameter controls how many particles are read together for each thread. If it is set to 3 and one uses 8 threads, batches of 3x8=24 particles will be read together. \
This may improve performance on systems where disk access, and particularly metadata handling of disk access, is a problem. It has a modest cost of increased RAM usage.");
	joboptions["do_preread_images"] = JobOption("Pre-read all particles into RAM?", false, "If set to Yes, all particle images will be read into computer memory, which will greatly speed up calculations on systems with slow disk access. However, one should of course be careful with the amount of RAM available. \
Because particles are read in float-precision, it will take ( N * box_size * box_size * 4 / (1024 * 1024 * 1024) ) Giga-bytes to read N particles into RAM. For 100 thousand 200x200 images, that becomes 15Gb, or 60 Gb for the same number of 400x400 particles. \
Remember that running a single MPI follower on each node that runs as many threads as available cores will have access to all available RAM. \n \n If parallel disc I/O is set to No, then only the leader reads all particles into RAM and sends those particles through the network to the MPI followers during the refinement iterations.");
	const char *default_scratch = getenv("RELION_SCRATCH_DIR");
	if (default_scratch == NULL)
	{
		default_scratch = DEFAULTSCRATCHDIR;
	}
	joboptions["scratch_dir"] = JobOption("Copy particles to scratch directory:", std::string(default_scratch), "If a directory is provided here, then the job will create a sub-directory in it called relion_volatile. If that relion_volatile directory already exists, it will be wiped. Then, the program will copy all input particles into a large stack inside the relion_volatile subdirectory. \
Provided this directory is on a fast local drive (e.g. an SSD drive), processing in all the iterations will be faster. If the job finishes correctly, the relion_volatile directory will be wiped. If the job crashes, you may want to remove it yourself.");
	joboptions["do_combine_thru_disc"] = JobOption("Combine iterations through disc?", false, "If set to Yes, at the end of every iteration all MPI followers will write out a large file with their accumulated results. The MPI leader will read in all these files, combine them all, and write out a new file with the combined results. \
All MPI salves will then read in the combined results. This reduces heavy load on the network, but increases load on the disc I/O. \
This will affect the time it takes between the progress-bar in the expectation step reaching its end (the mouse gets to the cheese) and the start of the ensuing maximisation step. It will depend on your system setup which is most efficient.");

	joboptions["use_gpu"] = JobOption("Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','. For example: '0,0:1,1:0,0:1,1'");
}

bool RelionJob::getCommandsInimodelJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();

	initialisePipeline(outputname, job_counter);

	if (joboptions["nr_mpi"].getNumber(error_message) > 1)
	{
		error_message = "Gradient refinement is not supported together with MPI.";
		return false;
	}
	if (error_message != "") return false;

	std::string command;
	command="`which relion_refine`";

	FileName fn_sym = joboptions["sym_name"].getString();

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
		// SHWS 10dec2020: switch off using run_ctXX output for continue jobs, as this will affect Schedulers
		//int it = (int)textToFloat((joboptions["fn_cont"].getString().substr(pos_it+3, 6)).c_str());
		//fn_run += "_ct" + floatToString(it);
		command += " --continue " + joboptions["fn_cont"].getString();
	}

	command += " --o " + outputname + fn_run;
        command += " --iter " + joboptions["nr_iter"].getString();

	int total_nr_iter = joboptions["nr_iter"].getNumber(error_message);
	if (error_message != "") return false;
    int nr_classes = joboptions["nr_classes"].getNumber(error_message);
	if (error_message != "") return false;


	if (!is_continue)
	{
		command += " --grad --denovo_3dref ";

		if (joboptions["fn_img"].getString() == "")
		{
			error_message = "ERROR: empty field for input STAR file...";
			return false;
		}
		command += " --i " + joboptions["fn_img"].getString();
		Node node(joboptions["fn_img"].getString(), joboptions["fn_img"].node_type);
		inputNodes.push_back(node);

		// CTF stuff
		if (joboptions["do_ctf_correction"].getBoolean())
		{
			command += " --ctf";
			if (joboptions["ctf_intact_first_peak"].getBoolean())
				command += " --ctf_intact_first_peak";
		}

		command += " --K " + joboptions["nr_classes"].getString();
		if (joboptions["do_run_C1"].getBoolean())
		{
			command += " --sym C1 ";
		}
		else
		{
			command += " --sym " + fn_sym;
		}

		if (joboptions["do_solvent"].getBoolean())
			command += " --flatten_solvent ";
		command += " --zero_mask ";
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
	command += " --pad 1 ";

	// Optimisation
	command += " --particle_diameter " + joboptions["particle_diameter"].getString();
	command += " --oversampling 1  --healpix_order 1  --offset_range 6  --offset_step 2 --auto_sampling ";
	command += " --tau2_fudge " + joboptions["tau_fudge"].getString();

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

	// Quickly remove RELION_JOB_EXIT_SUCCESS
	std::string command0 = "rm -f " + outputname + RELION_JOB_EXIT_SUCCESS;
	commands.push_back(command0);


	FileName fn_model;
	fn_model.compose(outputname + fn_run + "_it", total_nr_iter,"",3);
	fn_model+="_model.star";

	// Align with symmetry axes and apply symmetry
	std::string command2 = "`which relion_align_symmetry`";
	command2 += " --i " + fn_model;
	command2 += " --o " + outputname + "initial_model.mrc";

	if ( joboptions["do_run_C1"].getBoolean() && !(fn_sym.contains("C1") || fn_sym.contains("c1")) )
	{
		command2 += " --sym " + joboptions["sym_name"].getString();
	}
	else
	{
		command2 += " --sym C1 ";
	}
	command2 += " --apply_sym --select_largest_class ";
	commands.push_back(command2);

	// And re-introduce RELION_JOB_EXIT_SUCCESS
	std::string commandF = "touch " + outputname + RELION_JOB_EXIT_SUCCESS;
	commands.push_back(commandF);

	// Output nodes
	Node node2(outputname + "initial_model.mrc", LABEL_INIMOD_MAP);
    outputNodes.push_back(node2);

    // If doing more than 1 class, make them all available (one of them will be the same as initial_model.mrc)
    if (nr_classes > 1)
    {
        for (int iclass = 0; iclass < nr_classes; iclass++)
        {
			FileName fn_tmp;
			fn_tmp.compose(outputname + fn_run + "_it", total_nr_iter, "", 3);
			fn_tmp.compose(fn_tmp + "_class", iclass+1, "mrc", 3);
			Node node3(fn_tmp, LABEL_INIMOD_MAP);
			outputNodes.push_back(node3);
        }
    }

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}

void RelionJob::initialiseClass3DJob()
{
	hidden_name = ".gui_class3d";

	joboptions["fn_img"] = JobOption("Input images STAR file:", NODE_PARTS_CPIPE, "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");
	joboptions["fn_cont"] = JobOption("Continue from here: ", std::string(""), "STAR Files (*_optimiser.star)", "CURRENT_ODIR", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");
	joboptions["fn_ref"] = JobOption("Reference map:", NODE_MAP_CPIPE, "", "Image Files (*.{spi,vol,mrc})", "A 3D map in MRC/Spider format. \
	Make sure this map has the same dimensions and the same pixel size as your input images.");
	joboptions["fn_mask"] = JobOption("Reference mask (optional):", NODE_MASK_CPIPE, "", "Image Files (*.{spi,vol,msk,mrc})", "\
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
	joboptions["do_fast_subsets"] = JobOption("Use fast subsets (for large data sets)?", false, "If set to Yes, the first 5 iterations will be done with random subsets of only K*1500 particles (K being the number of classes); the next 5 with K*4500 particles, the next 5 with 30% of the data set; and the final ones with all data. This was inspired by a cisTEM implementation by Niko Grigorieff et al.");

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
	joboptions["sampling"] = JobOption("Angular sampling interval:", job_sampling_options, 2, "There are only a few discrete \
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
	joboptions["allow_coarser"] = JobOption("Allow coarser sampling?", false, "If set to Yes, the program will use coarser angular and translational samplings if the estimated accuracies of the assignments is still low in the earlier iterations. This may speed up the calculations.");
	joboptions["relax_sym"] = JobOption("Relax symmetry:", std::string(""), "With this option, poses related to the standard local angular search range by the given point group will also be explored. For example, if you have a pseudo-symmetric dimer A-A', refinement or classification in C1 with symmetry relaxation by C2 might be able to improve distinction between A and A'. Note that the reference must be more-or-less aligned to the convention of (pseudo-)symmetry operators. For details, see Ilca et al 2019 and Abrishami et al 2020 cited in the About dialog.");

	joboptions["do_helix"] = JobOption("Do helical reconstruction?", false, "If set to Yes, then perform 3D helical reconstruction.");
	joboptions["helical_tube_inner_diameter"] = JobOption("Tube diameter - inner (A):", std::string("-1"),"Inner and outer diameter (in Angstroms) of the reconstructed helix spanning across Z axis. \
Set the inner diameter to negative value if the helix is not hollow in the center. The outer diameter should be slightly larger than the actual width of helical tubes because it also decides the shape of 2D \
particle mask for each segment. If the psi priors of the extracted segments are not accurate enough due to high noise level or flexibility of the structure, then set the outer diameter to a large value.");
	joboptions["helical_tube_outer_diameter"] = JobOption("Tube diameter - outer (A):", std::string("-1"),"Inner and outer diameter (in Angstroms) of the reconstructed helix spanning across Z axis. \
Set the inner diameter to negative value if the helix is not hollow in the center. The outer diameter should be slightly larger than the actual width of helical tubes because it also decides the shape of 2D \
particle mask for each segment. If the psi priors of the extracted segments are not accurate enough due to high noise level or flexibility of the structure, then set the outer diameter to a large value.");
	joboptions["range_rot"] = JobOption("Angular search range - rot (deg):", std::string("-1"), "Local angular searches will be performed \
within +/- of the given amount (in degrees) from the optimal orientation in the previous iteration. The default negative value means that no local searches will be performed. \
A Gaussian prior will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.\n\nThese ranges will only be applied to the \
rot, tilt and psi angles in the first few iterations (global searches for orientations) in 3D helical reconstruction. \
Values of 9 or 15 degrees are commonly used. Higher values are recommended for more flexible structures and more memory and computation time will be used. \
A range of 15 degrees means sigma = 5 degrees.\n\nThese options will be invalid if you choose to perform local angular searches or not to perform image alignment on 'Sampling' tab.");
	joboptions["range_tilt"] = JobOption("Angular search range - tilt (deg):", std::string("15"), "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.\n\nThese ranges will only be applied to the \
rot, tilt and psi angles in the first few iterations (global searches for orientations) in 3D helical reconstruction. \
Values of 9 or 15 degrees are commonly used. Higher values are recommended for more flexible structures and more memory and computation time will be used. \
A range of 15 degrees means sigma = 5 degrees.\n\nThese options will be invalid if you choose to perform local angular searches or not to perform image alignment on 'Sampling' tab.");
	joboptions["range_psi"] = JobOption("Angular search range - psi (deg):", std::string("10"), "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.\n\nThese ranges will only be applied to the \
rot, tilt and psi angles in the first few iterations (global searches for orientations) in 3D helical reconstruction. \
Values of 9 or 15 degrees are commonly used. Higher values are recommended for more flexible structures and more memory and computation time will be used. \
A range of 15 degrees means sigma = 5 degrees.\n\nThese options will be invalid if you choose to perform local angular searches or not to perform image alignment on 'Sampling' tab.");
	joboptions["do_apply_helical_symmetry"] = JobOption("Apply helical symmetry?", true, "If set to Yes, helical symmetry will be applied in every iteration. Set to No if you have just started a project, helical symmetry is unknown or not yet estimated.");
	joboptions["helical_nr_asu"] = JobOption("Number of unique asymmetrical units:", 1, 1, 100, 1, "Number of unique helical asymmetrical units in each segment box. If the inter-box distance (set in segment picking step) \
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
	joboptions["keep_tilt_prior_fixed"] = JobOption("Keep tilt-prior fixed:", true, "If set to yes, the tilt prior will not change during the optimisation. If set to No, at each iteration the tilt prior will move to the optimal tilt value for that segment from the previous iteration.");

	joboptions["do_parallel_discio"] = JobOption("Use parallel disc I/O?", true, "If set to Yes, all MPI followers will read their own images from disc. \
Otherwise, only the leader will read images and send them through the network to the followers. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many followers reading in parallel. If your datasets contain particles with different box sizes, you have to say Yes.");
	joboptions["nr_pool"] = JobOption("Number of pooled particles:", 3, 1, 16, 1, "Particles are processed in individual batches by MPI followers. During each batch, a stack of particle images is only opened and closed once to improve disk access times. \
All particle images of a single batch are read into memory together. The size of these batches is at least one particle per thread used. The nr_pooled_particles parameter controls how many particles are read together for each thread. If it is set to 3 and one uses 8 threads, batches of 3x8=24 particles will be read together. \
This may improve performance on systems where disk access, and particularly metadata handling of disk access, is a problem. It has a modest cost of increased RAM usage.");
	joboptions["do_pad1"] = JobOption("Skip padding?", false, "If set to Yes, the calculations will not use padding in Fourier space for better interpolation in the references. Otherwise, references are padded 2x before Fourier transforms are calculated. Skipping padding (i.e. use --pad 1) gives nearly as good results as using --pad 2, but some artifacts may appear in the corners from signal that is folded back.");
	joboptions["do_preread_images"] = JobOption("Pre-read all particles into RAM?", false, "If set to Yes, all particle images will be read into computer memory, which will greatly speed up calculations on systems with slow disk access. However, one should of course be careful with the amount of RAM available. \
Because particles are read in float-precision, it will take ( N * box_size * box_size * 4 / (1024 * 1024 * 1024) ) Giga-bytes to read N particles into RAM. For 100 thousand 200x200 images, that becomes 15Gb, or 60 Gb for the same number of 400x400 particles. \
Remember that running a single MPI follower on each node that runs as many threads as available cores will have access to all available RAM. \n \n If parallel disc I/O is set to No, then only the leader reads all particles into RAM and sends those particles through the network to the MPI followers during the refinement iterations.");
	const char *default_scratch = getenv("RELION_SCRATCH_DIR");
	if (default_scratch == NULL)
	{
		default_scratch = DEFAULTSCRATCHDIR;
	}
	joboptions["scratch_dir"] = JobOption("Copy particles to scratch directory:", std::string(default_scratch), "If a directory is provided here, then the job will create a sub-directory in it called relion_volatile. If that relion_volatile directory already exists, it will be wiped. Then, the program will copy all input particles into a large stack inside the relion_volatile subdirectory. \
Provided this directory is on a fast local drive (e.g. an SSD drive), processing in all the iterations will be faster. If the job finishes correctly, the relion_volatile directory will be wiped. If the job crashes, you may want to remove it yourself.");
	joboptions["do_combine_thru_disc"] = JobOption("Combine iterations through disc?", false, "If set to Yes, at the end of every iteration all MPI followers will write out a large file with their accumulated results. The MPI leader will read in all these files, combine them all, and write out a new file with the combined results. \
All MPI salves will then read in the combined results. This reduces heavy load on the network, but increases load on the disc I/O. \
This will affect the time it takes between the progress-bar in the expectation step reaching its end (the mouse gets to the cheese) and the start of the ensuing maximisation step. It will depend on your system setup which is most efficient.");

	joboptions["use_gpu"] = JobOption("Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','.  For example: '0,0:1,1:0,0:1,1'");
}

bool RelionJob::getCommandsClass3DJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;

	if (joboptions["nr_mpi"].getNumber(error_message) > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";
	if (error_message != "") return false;

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
		// SHWS 10dec2020: switch off using run_ctXX output for continue jobs, as this will affect Schedulers
		//int it = (int)textToFloat((joboptions["fn_cont"].getString().substr(pos_it+3, 6)).c_str());
		//fn_run += "_ct" + floatToString(it);;
		command += " --continue " + joboptions["fn_cont"].getString();
	}

	command += " --o " + outputname + fn_run;

	int my_iter = (int)joboptions["nr_iter"].getNumber(error_message);
	if (error_message != "") return false;

	int my_classes = (int)joboptions["nr_classes"].getNumber(error_message);
	if (error_message != "") return false;

	outputNodes = getOutputNodesRefine(outputname + fn_run, my_iter, my_classes, 3, 1, "Class3D");

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

		if (joboptions["ini_high"].getNumber(error_message) > 0.)
			command += " --ini_high " + joboptions["ini_high"].getString();
		if (error_message != "") return false;

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
	if (joboptions["do_pad1"].getBoolean())
		command += " --pad 1 ";
	else
		command += " --pad 2 ";

	// CTF stuff
	if (!is_continue)
	{

		if (joboptions["do_ctf_correction"].getBoolean())
		{
			command += " --ctf";
			if (joboptions["ctf_intact_first_peak"].getBoolean())
				command += " --ctf_intact_first_peak";
		}
	}

	// Optimisation
	command += " --iter " + joboptions["nr_iter"].getString();
	command += " --tau2_fudge " + joboptions["tau_fudge"].getString();
        command += " --particle_diameter " + joboptions["particle_diameter"].getString();
	if (!is_continue)
	{
		if (joboptions["do_fast_subsets"].getBoolean())
			command += " --fast_subsets ";

		command += " --K " + joboptions["nr_classes"].getString();
		// Always flatten the solvent
		command += " --flatten_solvent";
		if (joboptions["do_zero_mask"].getBoolean())
			command += " --zero_mask";
		if (joboptions["highres_limit"].getNumber(error_message) > 0)
			command += " --strict_highres_exp " + joboptions["highres_limit"].getString();
		if (error_message != "") return false;

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
		int sampling = JobOption::getHealPixOrder(joboptions["sampling"].getString());
		if (sampling <= 0)
		{
			error_message = "Wrong choice for sampling";
			return false;
		}
		// The sampling given in the GUI will be the oversampled one!
		command += " --healpix_order " + integerToString(sampling - iover);

		// Manually input local angular searches
		if (joboptions["do_local_ang_searches"].getBoolean())
		{
			command += " --sigma_ang " + floatToString(joboptions["sigma_angles"].getNumber(error_message) / 3.);
			if (joboptions["relax_sym"].getString().length() > 0)
				command += " --relax_sym " + joboptions["relax_sym"].getString();

			if (error_message != "") return false;
		}

		// Offset range
		command += " --offset_range " + joboptions["offset_range"].getString();
		// The sampling given in the GUI will be the oversampled one!
		command += " --offset_step " +  floatToString(joboptions["offset_step"].getNumber(error_message) * pow(2., iover));
		if (error_message != "") return false;

		if (joboptions["allow_coarser"].getBoolean())
		{
			command += " --allow_coarser_sampling";
		}

	}

	// Provide symmetry, and always do norm and scale correction
	if (!is_continue)
	{
		command += " --sym " + joboptions["sym_name"].getString();
		command += " --norm --scale ";
	}

	if ( (!is_continue) && (joboptions["do_helix"].getBoolean()) )
	{
		label += ".helical";

		command += " --helix";

		float inner_diam = joboptions["helical_tube_inner_diameter"].getNumber(error_message);
		if (error_message != "") return false;
		if (inner_diam > 0.)
			command += " --helical_inner_diameter " + joboptions["helical_tube_inner_diameter"].getString();

		command += " --helical_outer_diameter " + joboptions["helical_tube_outer_diameter"].getString();
		if (joboptions["do_apply_helical_symmetry"].getBoolean())
		{
			command += " --helical_nr_asu " + joboptions["helical_nr_asu"].getString();
			command += " --helical_twist_initial " + joboptions["helical_twist_initial"].getString();
			command += " --helical_rise_initial " + joboptions["helical_rise_initial"].getString();

			float myz = joboptions["helical_z_percentage"].getNumber(error_message) / 100.;
			if (error_message != "") return false;
			command += " --helical_z_percentage " + floatToString(myz);

			if (joboptions["do_local_search_helical_symmetry"].getBoolean())
			{
				command += " --helical_symmetry_search";
				command += " --helical_twist_min " + joboptions["helical_twist_min"].getString();
				command += " --helical_twist_max " + joboptions["helical_twist_max"].getString();

				float twist_inistep = joboptions["helical_twist_inistep"].getNumber(error_message);
				if (error_message != "") return false;
				if (twist_inistep > 0.)
					command += " --helical_twist_inistep " + joboptions["helical_twist_inistep"].getString();

				command += " --helical_rise_min " + joboptions["helical_rise_min"].getString();
				command += " --helical_rise_max " + joboptions["helical_rise_max"].getString();

				float rise_inistep = joboptions["helical_rise_inistep"].getNumber(error_message);
				if (error_message != "") return false;
				if (rise_inistep > 0.)
					command += " --helical_rise_inistep " + joboptions["helical_rise_inistep"].getString();
			}
		}
		else
			command += " --ignore_helical_symmetry";
		if (joboptions["keep_tilt_prior_fixed"].getBoolean())
			command += " --helical_keep_tilt_prior_fixed";
		if ( (joboptions["dont_skip_align"].getBoolean()) && (!joboptions["do_local_ang_searches"].getBoolean()) )
		{
			float val = joboptions["range_tilt"].getNumber(error_message);
			if (error_message != "") return false;
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_tilt " + floatToString(val / 3.);

			val = joboptions["range_psi"].getNumber(error_message);
			if (error_message != "") return false;
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_psi " + floatToString(val / 3.);

			val = joboptions["range_rot"].getNumber(error_message);
			if (error_message != "") return false;
			val = (val < 0.) ? (0.) : (val);
			val = (val > 90.) ? (90.) : (val);
			command += " --sigma_rot " + floatToString(val / 3.);

			val = joboptions["helical_range_distance"].getNumber(error_message);
			if (error_message != "") return false;
			if (val > 0.)
				command += " --helical_sigma_distance " + floatToString(val / 3.);
		}
	}

	// Running stuff
	command += " --j " + joboptions["nr_threads"].getString();

	// GPU-stuff
	if (joboptions["use_gpu"].getBoolean())
	{
		if (!joboptions["dont_skip_align"].getBoolean())
		{
			error_message = "ERROR: you cannot use GPUs when skipping image alignments.";
			return false;
		}
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

	if (is_tomo)
	{
		joboptions["in_optimisation"] = JobOption("Input optimisation set: ", OUTNODE_TOMO_OPTIMISATION, "", "Optimisation set STAR file (*.star)", "Input tomo optimisation set. Input images STAR file, reference halfmaps and reference mask files will be extracted. If input files are specified below, then they will override the components in this optimisation set.");
	}
	joboptions["fn_img"] = JobOption("Input images STAR file:", NODE_PARTS_CPIPE, "", "STAR files (*.star) \t Image stacks (not recommended, read help!) (*.{spi,mrcs})", "A STAR file with all images (and their metadata). \n \n Alternatively, you may give a Spider/MRC stack of 2D images, but in that case NO metadata can be included and thus NO CTF correction can be performed, \
nor will it be possible to perform noise spectra estimation or intensity scale corrections in image groups. Therefore, running RELION with an input stack will in general provide sub-optimal results and is therefore not recommended!! Use the Preprocessing procedure to get the input STAR file in a semi-automated manner. Read the RELION wiki for more information.");
	joboptions["fn_cont"] = JobOption("Continue from here: ", std::string(""), "STAR Files (*_it*_optimiser.star)", "CURRENT_ODIR", "Select the *_optimiser.star file for the iteration \
from which you want to continue a previous run. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");
	joboptions["fn_ref"] = JobOption("Reference map:", NODE_MAP_CPIPE, "", "Image Files (*.{spi,vol,mrc})", "A 3D map in MRC/Spider format. \
	Make sure this map has the same dimensions and the same pixel size as your input images.");
	joboptions["fn_mask"] = JobOption("Reference mask (optional):", NODE_MASK_CPIPE, "", "Image Files (*.{spi,vol,msk,mrc})", "\
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

	joboptions["sampling"] = JobOption("Initial angular sampling:", job_sampling_options, 2, "There are only a few discrete \
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
	joboptions["auto_local_sampling"] = JobOption("Local searches from auto-sampling:", job_sampling_options, 4, "In the automated procedure to \
increase the angular samplings, local angular searches of -6/+6 times the sampling rate will be used from this angular sampling rate onwards. For most \
lower-symmetric particles a value of 1.8 degrees will be sufficient. Perhaps icosahedral symmetries may benefit from a smaller value such as 0.9 degrees.");
	joboptions["relax_sym"] = JobOption("Relax symmetry:", std::string(""), "With this option, poses related to the standard local angular search range by the given point group will also be explored. For example, if you have a pseudo-symmetric dimer A-A', refinement or classification in C1 with symmetry relaxation by C2 might be able to improve distinction between A and A'. Note that the reference must be more-or-less aligned to the convention of (pseudo-)symmetry operators. For details, see Ilca et al 2019 and Abrishami et al 2020 cited in the About dialog.");
	joboptions["auto_faster"] = JobOption("Use finer angular sampling faster?", false, "If set to Yes, then let auto-refinement proceed faster with finer angular samplings. Two additional command-line options will be passed to the refine program: \n \n \
--auto_ignore_angles lets angular sampling go down despite changes still happening in the angles \n \n \
--auto_resol_angles lets angular sampling go down if the current resolution already requires that sampling at the edge of the particle.  \n\n \
This option will make the computation faster, but hasn't been tested for many cases for potential loss in reconstruction quality upon convergence.");

	joboptions["do_helix"] = JobOption("Do helical reconstruction?", false, "If set to Yes, then perform 3D helical reconstruction.");
	joboptions["helical_tube_inner_diameter"] = JobOption("Tube diameter - inner (A):", std::string("-1"),"Inner and outer diameter (in Angstroms) of the reconstructed helix spanning across Z axis. \
Set the inner diameter to negative value if the helix is not hollow in the center. The outer diameter should be slightly larger than the actual width of helical tubes because it also decides the shape of 2D \
particle mask for each segment. If the psi priors of the extracted segments are not accurate enough due to high noise level or flexibility of the structure, then set the outer diameter to a large value.");
	joboptions["helical_tube_outer_diameter"] = JobOption("Tube diameter - outer (A):", std::string("-1"),"Inner and outer diameter (in Angstroms) of the reconstructed helix spanning across Z axis. \
Set the inner diameter to negative value if the helix is not hollow in the center. The outer diameter should be slightly larger than the actual width of helical tubes because it also decides the shape of 2D \
particle mask for each segment. If the psi priors of the extracted segments are not accurate enough due to high noise level or flexibility of the structure, then set the outer diameter to a large value.");
	joboptions["range_rot"] = JobOption("Angular search range - rot (deg):", std::string("-1"), "Local angular searches will be performed \
within +/- of the given amount (in degrees) from the optimal orientation in the previous iteration. The default negative value means that no local searches will be performed. \
A Gaussian prior will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.\n\nThese ranges will only be applied to the \
rot, tilt and psi angles in the first few iterations (global searches for orientations) in 3D helical reconstruction. \
Values of 9 or 15 degrees are commonly used. Higher values are recommended for more flexible structures and more memory and computation time will be used. \
A range of 15 degrees means sigma = 5 degrees.\n\nThese options will be invalid if you choose to perform local angular searches or not to perform image alignment on 'Sampling' tab.");
	joboptions["range_tilt"] = JobOption("Angular search range - tilt (deg):", std::string("15"), "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.\n\nThese ranges will only be applied to the \
rot, tilt and psi angles in the first few iterations (global searches for orientations) in 3D helical reconstruction. \
Values of 9 or 15 degrees are commonly used. Higher values are recommended for more flexible structures and more memory and computation time will be used. \
A range of 15 degrees means sigma = 5 degrees.\n\nThese options will be invalid if you choose to perform local angular searches or not to perform image alignment on 'Sampling' tab.");
	joboptions["range_psi"] = JobOption("Angular search range - psi (deg):", std::string("10"), "Local angular searches will be performed \
within +/- the given amount (in degrees) from the optimal orientation in the previous iteration. \
A Gaussian prior (also see previous option) will be applied, so that orientations closer to the optimal orientation \
in the previous iteration will get higher weights than those further away.\n\nThese ranges will only be applied to the \
rot, tilt and psi angles in the first few iterations (global searches for orientations) in 3D helical reconstruction. \
Values of 9 or 15 degrees are commonly used. Higher values are recommended for more flexible structures and more memory and computation time will be used. \
A range of 15 degrees means sigma = 5 degrees.\n\nThese options will be invalid if you choose to perform local angular searches or not to perform image alignment on 'Sampling' tab.");
	joboptions["do_apply_helical_symmetry"] = JobOption("Apply helical symmetry?", true, "If set to Yes, helical symmetry will be applied in every iteration. Set to No if you have just started a project, helical symmetry is unknown or not yet estimated.");
	joboptions["helical_nr_asu"] = JobOption("Number of unique asymmetrical units:", 1, 1, 100, 1, "Number of unique helical asymmetrical units in each segment box. If the inter-box distance (set in segment picking step) \
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
	joboptions["keep_tilt_prior_fixed"] = JobOption("Keep tilt-prior fixed:", true, "If set to yes, the tilt prior will not change during the optimisation. If set to No, at each iteration the tilt prior will move to the optimal tilt value for that segment from the previous iteration.");

	joboptions["do_parallel_discio"] = JobOption("Use parallel disc I/O?", true, "If set to Yes, all MPI followers will read their own images from disc. \
Otherwise, only the leader will read images and send them through the network to the followers. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many followers reading in parallel. If your datasets contain particles with different box sizes, you have to say Yes.");
	joboptions["nr_pool"] = JobOption("Number of pooled particles:", 3, 1, 16, 1, "Particles are processed in individual batches by MPI followers. During each batch, a stack of particle images is only opened and closed once to improve disk access times. \
All particle images of a single batch are read into memory together. The size of these batches is at least one particle per thread used. The nr_pooled_particles parameter controls how many particles are read together for each thread. If it is set to 3 and one uses 8 threads, batches of 3x8=24 particles will be read together. \
This may improve performance on systems where disk access, and particularly metadata handling of disk access, is a problem. It has a modest cost of increased RAM usage.");
	joboptions["do_pad1"] = JobOption("Skip padding?", false, "If set to Yes, the calculations will not use padding in Fourier space for better interpolation in the references. Otherwise, references are padded 2x before Fourier transforms are calculated. Skipping padding (i.e. use --pad 1) gives nearly as good results as using --pad 2, but some artifacts may appear in the corners from signal that is folded back.");
	joboptions["do_preread_images"] = JobOption("Pre-read all particles into RAM?", false, "If set to Yes, all particle images will be read into computer memory, which will greatly speed up calculations on systems with slow disk access. However, one should of course be careful with the amount of RAM available. \
Because particles are read in float-precision, it will take ( N * box_size * box_size * 8 / (1024 * 1024 * 1024) ) Giga-bytes to read N particles into RAM. For 100 thousand 200x200 images, that becomes 15Gb, or 60 Gb for the same number of 400x400 particles. \
Remember that running a single MPI follower on each node that runs as many threads as available cores will have access to all available RAM. \n \n If parallel disc I/O is set to No, then only the leader reads all particles into RAM and sends those particles through the network to the MPI followers during the refinement iterations.");
	const char *default_scratch = getenv("RELION_SCRATCH_DIR");
	if (default_scratch == NULL)
	{
		default_scratch = DEFAULTSCRATCHDIR;
	}
	joboptions["scratch_dir"] = JobOption("Copy particles to scratch directory:", std::string(default_scratch), "If a directory is provided here, then the job will create a sub-directory in it called relion_volatile. If that relion_volatile directory already exists, it will be wiped. Then, the program will copy all input particles into a large stack inside the relion_volatile subdirectory. \
Provided this directory is on a fast local drive (e.g. an SSD drive), processing in all the iterations will be faster. If the job finishes correctly, the relion_volatile directory will be wiped. If the job crashes, you may want to remove it yourself.");
	joboptions["do_combine_thru_disc"] = JobOption("Combine iterations through disc?", false, "If set to Yes, at the end of every iteration all MPI followers will write out a large file with their accumulated results. The MPI leader will read in all these files, combine them all, and write out a new file with the combined results. \
All MPI salves will then read in the combined results. This reduces heavy load on the network, but increases load on the disc I/O. \
This will affect the time it takes between the progress-bar in the expectation step reaching its end (the mouse gets to the cheese) and the start of the ensuing maximisation step. It will depend on your system setup which is most efficient.");
	joboptions["use_gpu"] = JobOption("Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','.  For example: '0,0:1,1:0,0:1,1'");
}

bool RelionJob::getCommandsAutorefineJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;

	if (joboptions["nr_mpi"].getNumber(error_message) > 1)
		command="`which relion_refine_mpi`";
	else
		command="`which relion_refine`";
	if (error_message != "") return false;

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
			error_message = "Invalid optimiser.star filename provided for auto-refine continuation run: " + joboptions["fn_cont"].getString();
                        return false;
                }

		// SHWS 10dec2020: switch off using run_ctXX output for continue jobs, as this will affect Schedulers
		//int it = (int)textToFloat((joboptions["fn_cont"].getString().substr(pos_it+3, 6)).c_str());
		//fn_run += "_ct" + floatToString(it);
		command += " --continue " + joboptions["fn_cont"].getString();

		// If is_continue we still need tomo optimisation set to create output optimisation set
		if (is_tomo && joboptions["in_optimisation"].getString() != "")
		{
			FileName fn_OS = joboptions["in_optimisation"].getString();
			Node node(fn_OS, joboptions["in_optimisation"].node_type);
			inputNodes.push_back(node);
			command += " --ios " + fn_OS;
		}
	}

	command += " --o " + outputname + fn_run;
	// TODO: add bodies!! (probably in next version)
	outputNodes = getOutputNodesRefine(outputname + fn_run, -1, 1, 3, 1, "Refine3D");

	if (is_tomo) label += ".tomo";

	if (!is_continue)
	{
		command += " --auto_refine --split_random_halves";

		// If tomo optimiser set is passed, fn_img and fn_ref can be empty
		if (is_tomo && joboptions["in_optimisation"].getString() != "")
		{
			// Optimiser set should contain particles, halfmap and refmask or they should be set especifically
			// If Optimiset set is passed without halfmaps or refmask, they cannot be set as "None" in the GUI.
			FileName fn_OS = joboptions["in_optimisation"].getString();
			Node node(fn_OS, joboptions["in_optimisation"].node_type);
			inputNodes.push_back(node);
			command += " --ios " + fn_OS;

			Node node1( outputname + fn_run + "_optimisation_set.star", LABEL_TOMO_OPTIMISATION);
			outputNodes.push_back(node1);

			if (joboptions["fn_mask"].getString() == "" && joboptions["do_solvent_fsc"].getBoolean())
				command += " --solvent_correct_fsc ";
			if (joboptions["fn_ref"].getString() == "" && !joboptions["ref_correct_greyscale"].getBoolean())
				command += " --firstiter_cc";
		}
		else if (joboptions["fn_img"].getString() == "")
		{
			error_message = "ERROR: empty field for input STAR file...";
			return false;
		}
		else if (joboptions["fn_ref"].getString() == "")
		{
			error_message = "ERROR: empty field for input reference...";
			return false;
		}

		if (joboptions["fn_img"].getString() != "")
		{
			command += " --i " + joboptions["fn_img"].getString();
			Node node(joboptions["fn_img"].getString(), joboptions["fn_img"].node_type);
			inputNodes.push_back(node);
		}

		FileName fn_ref = joboptions["fn_ref"].getString();
		if (fn_ref != "" && fn_ref != "None")
		{
			command += " --ref " + fn_ref;
			Node node(fn_ref, joboptions["fn_ref"].node_type);
			inputNodes.push_back(node);

			if (!joboptions["ref_correct_greyscale"].getBoolean())
				command += " --firstiter_cc";
		}
		if (joboptions["ini_high"].getNumber(error_message) > 0.)
		{
			if (error_message != "") return false;
			command += " --ini_high " + joboptions["ini_high"].getString();
		}

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
	if (joboptions["do_pad1"].getBoolean())
		command += " --pad 1 ";
	else
		command += " --pad 2 ";
	if (joboptions["auto_faster"].getBoolean())
	{
		command += " --auto_ignore_angles --auto_resol_angles";
	}

	// CTF stuff
	if (!is_continue)
	{
		if (joboptions["do_ctf_correction"].getBoolean())
		{
			command += " --ctf";
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

		int sampling = JobOption::getHealPixOrder(joboptions["sampling"].getString());
		if (sampling <= 0)
		{
			error_message = "Wrong choice for sampling";
			return false;
		}
		// The sampling given in the GUI will be the oversampled one!
		command += " --healpix_order " + integerToString(sampling - iover);

		// Minimum sampling rate to perform local searches (may be changed upon continuation
		int auto_local_sampling = JobOption::getHealPixOrder(joboptions["auto_local_sampling"].getString());
		if (auto_local_sampling <= 0)
		{
			error_message = "Wrong choice for auto_local_sampling";
			return false;
		}
		// The sampling given in the GUI will be the oversampled one!
		command += " --auto_local_healpix_order " + integerToString(auto_local_sampling - iover);

		// Offset range
		command += " --offset_range " + joboptions["offset_range"].getString();
		// The sampling given in the GUI will be the oversampled one!
		command += " --offset_step " + floatToString(joboptions["offset_step"].getNumber(error_message) * pow(2., iover));
		if (error_message != "") return false;

		command += " --sym " + joboptions["sym_name"].getString();
		// Always join low-res data, as some D&I point group refinements may fall into different hands!
		command += " --low_resol_join_halves 40";
		command += " --norm --scale ";

		// Helix
		if (joboptions["do_helix"].getBoolean())
		{
			label += ".helical";

			command += " --helix";

			float inner_diam = joboptions["helical_tube_inner_diameter"].getNumber(error_message);
			if (error_message != "") return false;
			if (inner_diam > 0.)
				command += " --helical_inner_diameter " + joboptions["helical_tube_inner_diameter"].getString();

			command += " --helical_outer_diameter " + joboptions["helical_tube_outer_diameter"].getString();
			if (joboptions["do_apply_helical_symmetry"].getBoolean())
			{
				command += " --helical_nr_asu " + joboptions["helical_nr_asu"].getString();
				command += " --helical_twist_initial " + joboptions["helical_twist_initial"].getString();
				command += " --helical_rise_initial " + joboptions["helical_rise_initial"].getString();

				float myz = joboptions["helical_z_percentage"].getNumber(error_message) / 100.;
				if (error_message != "") return false;
				command += " --helical_z_percentage " + floatToString(myz);

				if (joboptions["do_local_search_helical_symmetry"].getBoolean())
				{
					command += " --helical_symmetry_search";
					command += " --helical_twist_min " + joboptions["helical_twist_min"].getString();
					command += " --helical_twist_max " + joboptions["helical_twist_max"].getString();

					float twist_inistep = joboptions["helical_twist_inistep"].getNumber(error_message);
					if (error_message != "") return false;
					if (twist_inistep > 0.)
						command += " --helical_twist_inistep " + joboptions["helical_twist_inistep"].getString();

					command += " --helical_rise_min " + joboptions["helical_rise_min"].getString();
					command += " --helical_rise_max " + joboptions["helical_rise_max"].getString();

					float rise_inistep = joboptions["helical_rise_inistep"].getNumber(error_message);
					if (error_message != "") return false;
					if (rise_inistep > 0.)
						command += " --helical_rise_inistep " + joboptions["helical_rise_inistep"].getString();
				}
			}
			else
				command += " --ignore_helical_symmetry";

			float val;
			if (sampling != auto_local_sampling)
			{
				val = joboptions["range_tilt"].getNumber(error_message);
				if (error_message != "") return false;
				val = (val < 0.) ? (0.) : (val);
				val = (val > 90.) ? (90.) : (val);
				command += " --sigma_tilt " + floatToString(val / 3.);

				val = joboptions["range_psi"].getNumber(error_message);
				if (error_message != "") return false;
				val = (val < 0.) ? (0.) : (val);
				val = (val > 90.) ? (90.) : (val);
				command += " --sigma_psi " + floatToString(val / 3.);

				val = joboptions["range_rot"].getNumber(error_message);
				if (error_message != "") return false;
				val = (val < 0.) ? (0.) : (val);
				val = (val > 90.) ? (90.) : (val);
				command += " --sigma_rot " + floatToString(val / 3.);
			}

			val = joboptions["helical_range_distance"].getNumber(error_message);
			if (error_message != "") return false;
			if (val > 0.)
				command += " --helical_sigma_distance " + floatToString(val / 3.);

			if (joboptions["keep_tilt_prior_fixed"].getBoolean())
				command += " --helical_keep_tilt_prior_fixed";
		}
	}

	if (joboptions["relax_sym"].getString().length() > 0)
		command += " --relax_sym " + joboptions["relax_sym"].getString();

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

void RelionJob::initialiseMultiBodyJob()
{
	type = PROC_MULTIBODY;

	hidden_name = ".gui_multibody";

	joboptions["fn_in"] = JobOption("Consensus refinement optimiser.star: ", std::string(""), "STAR Files (run_it*_optimiser.star)", "Refine3D/.", "Select the *_optimiser.star file for the iteration of the consensus refinement \
from which you want to start multi-body refinement.");

	joboptions["fn_cont"] = JobOption("Continue from here: ", std::string(""), "STAR Files (*_optimiser.star)", "CURRENT_ODIR", "Select the *_optimiser.star file for the iteration \
from which you want to continue this multi-body refinement. \
Note that the Output rootname of the continued run and the rootname of the previous run cannot be the same. \
If they are the same, the program will automatically add a '_ctX' to the output rootname, \
with X being the iteration from which one continues the previous run.");

	joboptions["fn_bodies"] = JobOption("Body STAR file:", std::string(""), "STAR Files (*.{star})", ".", " Provide the STAR file with all information about the bodies to be used in multi-body refinement. \
An example for a three-body refinement would look like this: \n\
\n\
data_\n\
loop_\n\
_rlnBodyMaskName\n\
_rlnBodyRotateRelativeTo\n\
_rlnBodySigmaAngles\n\
_rlnBodySigmaOffset\n\
large_body_mask.mrc 2 10 2\n\
small_body_mask.mrc 1 10 2\n\
head_body_mask.mrc 2 10 2\n\
\n\
Where each data line represents a different body, and: \n \
 - rlnBodyMaskName contains the name of a soft-edged mask with values in [0,1] that define the body; \n\
 - rlnBodyRotateRelativeTo defines relative to which other body this body rotates (first body is number 1); \n\
 - rlnBodySigmaAngles and _rlnBodySigmaOffset are the standard deviations (widths) of Gaussian priors on the consensus rotations and translations; \n\
\n \
Optionally, there can be a fifth column with _rlnBodyReferenceName. Entries can be 'None' (without the ''s) or the name of a MRC map with an initial reference for that body. In case the entry is None, the reference will be taken from the density in the consensus refinement.\n \n\
Also note that larger bodies should be above smaller bodies in the STAR file. For more information, see the multi-body paper.");

	joboptions["do_subtracted_bodies"] = JobOption("Reconstruct subtracted bodies?", true, "If set to Yes, then the reconstruction of each of the bodies will use the subtracted images. This may give \
useful insights about how well the subtraction worked. If set to No, the original particles are used for reconstruction (while the subtracted ones are still used for alignment). This will result in fuzzy densities for bodies outside the one used for refinement.");

	joboptions["sampling"] = JobOption("Initial angular sampling:", job_sampling_options, 4, "There are only a few discrete \
angular samplings possible because we use the HealPix library to generate the sampling of the first two Euler angles on the sphere. \
The samplings are approximate numbers and vary slightly over the sphere.\n\n \
Note that this will only be the value for the first few iteration(s): the sampling rate will be increased automatically after that.");
	joboptions["offset_range"] = JobOption("Initial offset range (pix):", 3, 0, 30, 1, "Probabilities will be calculated only for translations \
in a circle with this radius (in pixels). The center of this circle changes at every iteration and is placed at the optimal translation \
for each image in the previous iteration.\n\n \
Note that this will only be the value for the first few iteration(s): the sampling rate will be increased automatically after that.");
	joboptions["offset_step"] = JobOption("Initial offset step (pix):", 0.75, 0.1, 5, 0.1, "Translations will be sampled with this step-size (in pixels).\
Translational sampling is also done using the adaptive approach. \
Therefore, if adaptive=1, the translations will first be evaluated on a 2x coarser grid.\n\n \
Note that this will only be the value for the first few iteration(s): the sampling rate will be increased automatically after that.");


	joboptions["do_analyse"] = JobOption("Run flexibility analysis?", true, "If set to Yes, after the multi-body refinement has completed, a PCA analysis will be run on the orientations all all bodies in the data set. This can be set to No initially, and then the job can be continued afterwards to only perform this analysis.");
	joboptions["nr_movies"] = JobOption("Number of eigenvector movies:", 3, 0, 16, 1, "Series of ten output maps will be generated along this many eigenvectors. These maps can be opened as a 'Volume Series' in UCSF Chimera, and then displayed as a movie. They represent the principal motions in the particles.");
	joboptions["do_select"] = JobOption("Select particles based on eigenvalues?", false, "If set to Yes, a particles.star file is written out with all particles that have the below indicated eigenvalue in the selected range.");
	joboptions["select_eigenval"] = JobOption("Select on eigenvalue:", 1, 1, 20, 1, "This is the number of the eigenvalue to be used in the particle subset selection (start counting at 1).");
	joboptions["eigenval_min"] = JobOption("Minimum eigenvalue:", -999., -50, 50, 1, "This is the minimum value for the selected eigenvalue; only particles with the selected eigenvalue larger than this value will be included in the output particles.star file");
	joboptions["eigenval_max"] = JobOption("Maximum eigenvalue:", 999., -50, 50, 1, "This is the maximum value for the selected eigenvalue; only particles with the selected eigenvalue less than this value will be included in the output particles.star file");

	joboptions["do_parallel_discio"] = JobOption("Use parallel disc I/O?", true, "If set to Yes, all MPI followers will read their own images from disc. \
Otherwise, only the leader will read images and send them through the network to the followers. Parallel file systems like gluster of fhgfs are good at parallel disc I/O. NFS may break with many followers reading in parallel. If your datasets contain particles with different box sizes, you have to say Yes.");
	joboptions["nr_pool"] = JobOption("Number of pooled particles:", 3, 1, 16, 1, "Particles are processed in individual batches by MPI followers. During each batch, a stack of particle images is only opened and closed once to improve disk access times. \
All particle images of a single batch are read into memory together. The size of these batches is at least one particle per thread used. The nr_pooled_particles parameter controls how many particles are read together for each thread. If it is set to 3 and one uses 8 threads, batches of 3x8=24 particles will be read together. \
This may improve performance on systems where disk access, and particularly metadata handling of disk access, is a problem. It has a modest cost of increased RAM usage.");
	joboptions["do_pad1"] = JobOption("Skip padding?", false, "If set to Yes, the calculations will not use padding in Fourier space for better interpolation in the references. Otherwise, references are padded 2x before Fourier transforms are calculated. Skipping padding (i.e. use --pad 1) gives nearly as good results as using --pad 2, but some artifacts may appear in the corners from signal that is folded back.");
	joboptions["do_preread_images"] = JobOption("Pre-read all particles into RAM?", false, "If set to Yes, all particle images will be read into computer memory, which will greatly speed up calculations on systems with slow disk access. However, one should of course be careful with the amount of RAM available. \
Because particles are read in float-precision, it will take ( N * box_size * box_size * 8 / (1024 * 1024 * 1024) ) Giga-bytes to read N particles into RAM. For 100 thousand 200x200 images, that becomes 15Gb, or 60 Gb for the same number of 400x400 particles. \
Remember that running a single MPI follower on each node that runs as many threads as available cores will have access to all available RAM. \n \n If parallel disc I/O is set to No, then only the leader reads all particles into RAM and sends those particles through the network to the MPI followers during the refinement iterations.");
	const char *default_scratch = getenv("RELION_SCRATCH_DIR");
	if (default_scratch == NULL)
	{
		default_scratch = DEFAULTSCRATCHDIR;
	}
	joboptions["scratch_dir"] = JobOption("Copy particles to scratch directory:", std::string(default_scratch), "If a directory is provided here, then the job will create a sub-directory in it called relion_volatile. If that relion_volatile directory already exists, it will be wiped. Then, the program will copy all input particles into a large stack inside the relion_volatile subdirectory. \
Provided this directory is on a fast local drive (e.g. an SSD drive), processing in all the iterations will be faster. If the job finishes correctly, the relion_volatile directory will be wiped. If the job crashes, you may want to remove it yourself.");
	joboptions["do_combine_thru_disc"] = JobOption("Combine iterations through disc?", false, "If set to Yes, at the end of every iteration all MPI followers will write out a large file with their accumulated results. The MPI leader will read in all these files, combine them all, and write out a new file with the combined results. \
All MPI salves will then read in the combined results. This reduces heavy load on the network, but increases load on the disc I/O. \
This will affect the time it takes between the progress-bar in the expectation step reaching its end (the mouse gets to the cheese) and the start of the ensuing maximisation step. It will depend on your system setup which is most efficient.");
	joboptions["use_gpu"] = JobOption("Use GPU acceleration?", false, "If set to Yes, the job will try to use GPU acceleration.");
	joboptions["gpu_ids"] = JobOption("Which GPUs to use:", std::string(""), "This argument is not necessary. If left empty, the job itself will try to allocate available GPU resources. You can override the default allocation by providing a list of which GPUs (0,1,2,3, etc) to use. MPI-processes are separated by ':', threads by ','.  For example: '0,0:1,1:0,0:1,1'");
}

bool RelionJob::getCommandsMultiBodyJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;

	if (!exists(joboptions["fn_bodies"].getString()))
	{
		error_message = "ERROR: you have to specify an existing body STAR file.";
		return false;
	}

	if (is_continue && joboptions["fn_cont"].getString() == "" && !joboptions["do_analyse"].getBoolean())
	{
		error_message = "ERROR: either specify a optimiser file to continue multibody refinement from; OR run flexibility analysis...";
		return false;
	}

	FileName fn_run = "";
	if (!is_continue || (is_continue && joboptions["fn_cont"].getString() != ""))
	{

		if (joboptions["nr_mpi"].getNumber(error_message) > 1)
			command="`which relion_refine_mpi`";
		else
			command="`which relion_refine`";
		if (error_message != "") return false;

		MetaDataTable MD;
		MD.read(joboptions["fn_bodies"].getString());
		int nr_bodies = MD.numberOfObjects();

		if (is_continue)
		{
			int pos_it = joboptions["fn_cont"].getString().rfind("_it");
			int pos_op = joboptions["fn_cont"].getString().rfind("_optimiser");
			if (pos_it < 0 || pos_op < 0)
				std::cerr << "Warning: invalid optimiser.star filename provided for continuation run: " << joboptions["fn_cont"].getString() << std::endl;
			int it = (int)textToFloat((joboptions["fn_cont"].getString().substr(pos_it+3, 6)).c_str());
			fn_run = "run_ct" + floatToString(it);
			command += " --continue " + joboptions["fn_cont"].getString();
			command += " --o " + outputname + fn_run;
			outputNodes = getOutputNodesRefine(outputname + fn_run, -1, 1, 3, nr_bodies, "MultiBody");

		}
		else
		{
			fn_run = "run";
			command += " --continue " + joboptions["fn_in"].getString();
			command += " --o " + outputname + fn_run;
			outputNodes = getOutputNodesRefine(outputname + "run", -1, 1, 3, nr_bodies, "MultiBody");
			command += " --solvent_correct_fsc --multibody_masks " + joboptions["fn_bodies"].getString();

			Node node(joboptions["fn_in"].getString(), LABEL_REFINE3D_OPT);
			inputNodes.push_back(node);

			// Sampling
			int iover = 1;
			command += " --oversampling " + floatToString((float)iover);
			int sampling = JobOption::getHealPixOrder(joboptions["sampling"].getString());
			if (sampling <= 0)
			{
				error_message = "Wrong choice for sampling";
				return false;
			}
			// The sampling given in the GUI will be the oversampled one!
			command += " --healpix_order " + integerToString(sampling - iover);
			// Always perform local searches!
			command += " --auto_local_healpix_order " + integerToString(sampling - iover);

			// Offset range
			command += " --offset_range " + joboptions["offset_range"].getString();
			// The sampling given in the GUI will be the oversampled one!
			command += " --offset_step " + floatToString(joboptions["offset_step"].getNumber(error_message) * pow(2., iover));
			if (error_message != "") return false;
		}

		if (joboptions["do_subtracted_bodies"].getBoolean())
			command += " --reconstruct_subtracted_bodies ";

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
		if (joboptions["do_pad1"].getBoolean())
			command += " --pad 1 ";
		else
			command += " --pad 2 ";

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
	} // end if (!is_continue || (is_continue && joboptions["fn_cont"].getString() != ""))

	if (joboptions["do_analyse"].getBoolean())
	{
		command = "`which relion_flex_analyse`";

		// If we had performed relion_refine command, then fn_run would be set now
		// Otherwise, we have to search for _model.star files that do NOT have a _it??? specifier
		if (fn_run == "")
		{
			FileName fn_wildcard = outputname + "run*_model.star";
			std::vector<FileName> fns_model;
			std::vector<FileName> fns_ok;
			fn_wildcard.globFiles(fns_model);
			for (int i = 0; i < fns_model.size(); i++)
			{
				if (!fns_model[i].contains("_it"))
					fns_ok.push_back(fns_model[i]);
			}
			if (fns_ok.size() == 0)
			{
				error_message = "ERROR: cannot find appropriate model.star file in the output directory";
				return false;
			}
			if (fns_ok.size() > 1)
			{
				error_message = "ERROR: there are more than one model.star files (without '_it' specifiers) in the output directory. Move all but one out of the way.";
				return false;
			}
			fn_run = fns_ok[0].beforeFirstOf("_model.star");
		}
		else
			fn_run = outputname + fn_run;

		// General I/O
		command += " --PCA_orient ";
		command += " --model " + fn_run + "_model.star";
		command += " --data " + fn_run + "_data.star";
		command += " --bodies " + joboptions["fn_bodies"].getString();
		command += " --o " + outputname + "analyse";

		// Eigenvector movie maps
		if (joboptions["nr_movies"].getNumber(error_message) > 0)
		{
			command += " --do_maps ";
			command += " --k " + joboptions["nr_movies"].getString();
		}
		if (error_message != "") return false;

		// Selection
		if (joboptions["do_select"].getBoolean())
		{
			float minval = joboptions["eigenval_min"].getNumber(error_message);
			if (error_message != "") return false;

			float maxval = joboptions["eigenval_max"].getNumber(error_message);
			if (error_message != "") return false;

			if ( minval >= maxval)
			{
				error_message = "ERROR: the maximum eigenvalue should be larger than the minimum one!";
				return false;
			}

			command += " --select_eigenvalue " + joboptions["select_eigenval"].getString();
			command += " --select_eigenvalue_min " + joboptions["eigenval_min"].getString();
			command += " --select_eigenvalue_max " + joboptions["eigenval_max"].getString();

			// Add output node: selected particles star file
			FileName fnt = outputname + "analyse_eval"+integerToString(joboptions["select_eigenval"].getNumber(error_message),3)+"_select";
			if (error_message != "") return false;

			int min = ROUND(joboptions["eigenval_min"].getNumber(error_message));
			if (error_message != "") return false;

			int max = ROUND(joboptions["eigenval_max"].getNumber(error_message));
			if (error_message != "") return false;

			if (min > -99998)
				fnt += "_min"+integerToString(min);
			if (max < 99998)
				fnt += "_max"+integerToString(max);
			fnt += ".star";
			Node node2(fnt, LABEL_MULTIBODY_SEL_PARTS);
			outputNodes.push_back(node2);

		}

		// PDF with histograms of the eigenvalues
		Node node3(outputname + "analyse_logfile.pdf", LABEL_MULTIBODY_FLEXLOG);
		outputNodes.push_back(node3);

		commands.push_back(command);
	}

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}

void RelionJob::initialiseMaskcreateJob()
{
	hidden_name = ".gui_maskcreate";

	joboptions["fn_in"] = JobOption("Input 3D map:", NODE_MAP_CPIPE, "", "MRC map files (*.mrc)", "Provide an input MRC map from which to start binarizing the map.");

	joboptions["lowpass_filter"] = JobOption("Lowpass filter map (A)", 15, 10, 100, 5, "Lowpass filter that will be applied to the input map, prior to binarization. To calculate solvent masks, a lowpass filter of 15-20A may work well.");
	joboptions["angpix"] = JobOption("Pixel size (A)", -1, 0.3, 5, 0.1, "Provide the pixel size of the input map in Angstroms to calculate the low-pass filter. This value is also used in the output image header.");

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
	initialisePipeline(outputname, job_counter);
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
	Node node2(outputname + "mask.mrc", LABEL_MASK3D_MASK);
	outputNodes.push_back(node2);

	if (joboptions["lowpass_filter"].getNumber(error_message) > 0)
	{
		command += " --lowpass " + joboptions["lowpass_filter"].getString();
	}
	if (error_message != "") return false;

	if (joboptions["angpix"].getNumber(error_message) > 0)
	{
		command += " --angpix " + joboptions["angpix"].getString();
	}
	if (error_message != "") return false;

	command += " --ini_threshold " + joboptions["inimask_threshold"].getString();
	command += " --extend_inimask " + joboptions["extend_inimask"].getString();
	command += " --width_soft_edge " + joboptions["width_mask_edge"].getString();

	if (joboptions["do_helix"].getBoolean())
	{
		command += " --helix --z_percentage " + floatToString(joboptions["helical_z_percentage"].getNumber(error_message) / 100.);
		if (error_message != "") return false;
	}

	// Running stuff
	command += " --j " + joboptions["nr_threads"].getString();

	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}

void RelionJob::initialiseJoinstarJob()
{
	hidden_name = ".gui_joinstar";

	joboptions["do_part"] = JobOption("Combine particle STAR files?", false, "");
	joboptions["fn_part1"] = JobOption("Particle STAR file 1: ", NODE_PARTS_CPIPE, "", "particle STAR file (*.star)", "The first of the particle STAR files to be combined.");
	joboptions["fn_part2"] = JobOption("Particle STAR file 2: ", NODE_PARTS_CPIPE, "", "particle STAR file (*.star)", "The second of the particle STAR files to be combined.");
	joboptions["fn_part3"] = JobOption("Particle STAR file 3: ", NODE_PARTS_CPIPE, "", "particle STAR file (*.star)", "The third of the particle STAR files to be combined. Leave empty if there are only two files to be combined.");
	joboptions["fn_part4"] = JobOption("Particle STAR file 4: ", NODE_PARTS_CPIPE, "", "particle STAR file (*.star)", "The fourth of the particle STAR files to be combined. Leave empty if there are only two or three files to be combined.");

	joboptions["do_mic"] = JobOption("Combine micrograph STAR files?", false, "");
	joboptions["fn_mic1"] = JobOption("Micrograph STAR file 1: ", NODE_MICS_CPIPE, "", "micrograph STAR file (*.star)", "The first of the micrograph STAR files to be combined.");
	joboptions["fn_mic2"] = JobOption("Micrograph STAR file 2: ", NODE_MICS_CPIPE, "", "micrograph STAR file (*.star)", "The second of the micrograph STAR files to be combined.");
	joboptions["fn_mic3"] = JobOption("Micrograph STAR file 3: ", NODE_MICS_CPIPE, "", "micrograph STAR file (*.star)", "The third of the micrograph STAR files to be combined. Leave empty if there are only two files to be combined.");
	joboptions["fn_mic4"] = JobOption("Micrograph STAR file 4: ", NODE_MICS_CPIPE, "", "micrograph STAR file (*.star)", "The fourth of the micrograph STAR files to be combined. Leave empty if there are only two or three files to be combined.");

	joboptions["do_mov"] = JobOption("Combine movie STAR files?", false, "");
	joboptions["fn_mov1"] = JobOption("Movie STAR file 1: ", NODE_MOVIES_CPIPE, "", "movie STAR file (*.star)", "The first of the micrograph movie STAR files to be combined.");
	joboptions["fn_mov2"] = JobOption("Movie STAR file 2: ", NODE_MOVIES_CPIPE, "", "movie STAR file (*.star)", "The second of the micrograph movie STAR files to be combined.");
	joboptions["fn_mov3"] = JobOption("Movie STAR file 3: ", NODE_MOVIES_CPIPE, "", "movie STAR file (*.star)", "The third of the micrograph movie STAR files to be combined. Leave empty if there are only two files to be combined.");
	joboptions["fn_mov4"] = JobOption("Movie STAR file 4: ", NODE_MOVIES_CPIPE, "", "movie STAR file (*.star)", "The fourth of the micrograph movie STAR files to be combined. Leave empty if there are only two or three files to be combined.");
}

bool RelionJob::getCommandsJoinstarJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;
	command="`which relion_star_handler`";

	int ii = 0;
	if (joboptions["do_part"].getBoolean())
	{
		ii++;
		label += ".particles";
	}
	if (joboptions["do_mic"].getBoolean())
	{
		ii++;
		label += ".micrographs";
	}
	if (joboptions["do_mov"].getBoolean())
	{
		ii++;
		label += ".movies";
	}

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
		command += " --combine --i \" " + joboptions["fn_part1"].getString();
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
		command += " --combine --i \" " + joboptions["fn_mic1"].getString();
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
		command += " --combine --i \" " + joboptions["fn_mov1"].getString();
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

	joboptions["fn_opt"] = JobOption("Input optimiser.star: ", NODE_OPTIMISER_CPIPE, "", "STAR Files (*_optimiser.star)", "Select the *_optimiser.star file for the iteration of the 3D refinement/classification \
which you want to use for subtraction. It will use the maps from this run for the subtraction, and of no particles input STAR file is given below, it will use all of the particles from this run.");
	joboptions["fn_mask"] = JobOption("Mask of the signal to keep:", NODE_MASK_CPIPE, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a soft mask where the protein density you wish to subtract from the experimental particles is black (0) and the density you wish to keep is white (1).");
	joboptions["do_data"] = JobOption("Use different particles?", false, "If set to Yes, subtraction will be performed on the particles in the STAR file below, instead of on all the particles of the 3D refinement/classification from the optimiser.star file.");
	joboptions["fn_data"] = JobOption("Input particle star file:", NODE_PARTS_CPIPE, "", "particle STAR file (*.star)", "The particle STAR files with particles that will be used in the subtraction. Leave this field empty if all particles from the input refinement/classification run are to be used.");
	joboptions["do_float16"] = JobOption("Write output in float16?", true ,"If set to Yes, this program will write output images in float16 MRC format. This will save a factor of two in disk space compared to the default of writing in float32. Note that RELION and CCPEM will read float16 images, but other programs may not (yet) do so.");

	joboptions["do_fliplabel"] = JobOption("OR revert to original particles?", false, "If set to Yes, no signal subtraction is performed. Instead, the labels of rlnImageName and rlnImageOriginalName are flipped in the input STAR file given in the field below. This will make the STAR file point back to the original, non-subtracted images.");
	joboptions["fn_fliplabel"] = JobOption("revert this particle star file:", NODE_PARTS_CPIPE, "", "particle STAR file (*.star)", "The particle STAR files with particles that will be used for label reversion.");

	joboptions["do_center_mask"] = JobOption("Do center subtracted images on mask?", true, "If set to Yes, the subtracted particles will be centered on projections of the center-of-mass of the input mask.");
	joboptions["do_center_xyz"] = JobOption("Do center on my coordinates?", false, "If set to Yes, the subtracted particles will be centered on projections of the x,y,z coordinates below. The unit is pixel, not angstrom. The origin is at the center of the box, not at the corner.");
	joboptions["center_x"] = JobOption("Center coordinate (pix) - X:", std::string("0"), "X-coordinate of the 3D center (in pixels).");
	joboptions["center_y"] = JobOption("Center coordinate (pix) - Y:", std::string("0"), "Y-coordinate of the 3D center (in pixels).");
	joboptions["center_z"] = JobOption("Center coordinate (pix) - Z:", std::string("0"), "Z-coordinate of the 3D center (in pixels).");

	joboptions["new_box"] = JobOption("New box size:", -1, 64, 512, 32, "Provide a non-negative value to re-window the subtracted particles in a smaller box size." );
}

bool RelionJob::getCommandsSubtractJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;

	if (joboptions["do_fliplabel"].getBoolean())
	{
		if (joboptions["nr_mpi"].getNumber(error_message) > 1)
		{
			error_message = "You cannot use MPI parallelization to revert particle labels.";
			return false;
		}

		Node node(joboptions["fn_fliplabel"].getString(), joboptions["fn_fliplabel"].node_type);
		inputNodes.push_back(node);

		Node node2(outputname + "original.star", LABEL_SUBTRACT_REVERTED);
		outputNodes.push_back(node2);

		label += ".revert";

		command = "`which relion_particle_subtract`";
		command += " --revert " + joboptions["fn_fliplabel"].getString() + " --o " + outputname;
	}
	else
	{
		if (joboptions["nr_mpi"].getNumber(error_message) > 1)
			command="`which relion_particle_subtract_mpi`";
		else
			command="`which relion_particle_subtract`";
		if (error_message != "") return false;

		// I/O
		if (joboptions["fn_opt"].getString() == "")
		{
			error_message = "ERROR: empty field for input optimiser.star...";
			return false;
		}
		command += " --i " + joboptions["fn_opt"].getString();
		Node node(joboptions["fn_opt"].getString(), LABEL_OPTIMISER_CPIPE);
		inputNodes.push_back(node);

		if (joboptions["fn_mask"].getString() != "")
		{
			command += " --mask " + joboptions["fn_mask"].getString();
			Node node2(joboptions["fn_mask"].getString(), joboptions["fn_mask"].node_type);
			inputNodes.push_back(node2);
		}
		if (joboptions["do_data"].getBoolean())
		{
			if (joboptions["fn_data"].getString() == "")
			{
				error_message = "ERROR: empty field for the input particle STAR file...";
				return false;
			}
			command += " --data " + joboptions["fn_data"].getString();
			Node node3(joboptions["fn_data"].getString(), joboptions["fn_data"].node_type);
			inputNodes.push_back(node3);
		}

		command += " --o " + outputname;
		Node node4(outputname + "particles_subtracted.star", LABEL_SUBTRACT_SUBTRACTED);
		outputNodes.push_back(node4);

		if (joboptions["do_center_mask"].getBoolean())
		{
			command += " --recenter_on_mask";
		}
		else if (joboptions["do_center_xyz"].getBoolean())
		{
			command += " --center_x " + joboptions["center_x"].getString();
			command += " --center_y " + joboptions["center_y"].getString();
			command += " --center_z " + joboptions["center_z"].getString();
		}

		if (joboptions["do_float16"].getBoolean())
		{
			command += " --float16 ";
		}

		if (joboptions["new_box"].getNumber(error_message) > 0)
		{
			command += " --new_box " + joboptions["new_box"].getString();
		}
		if (error_message != "") return false;

	}

	// Other arguments
	command += " " + joboptions["other_args"].getString();

	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}

void RelionJob::initialisePostprocessJob()
{
	hidden_name = ".gui_post";

	if (is_tomo)
	{
		joboptions["in_optimisation"] = JobOption("Input optimisation set: ", OUTNODE_TOMO_OPTIMISATION, "", "Optimisation set STAR file (*.star)", "Input tomo optimisation set. Half map files will be extracted. If half maps are specified below, then they will override the components in this optimisation set.");
	}
	joboptions["fn_in"] = JobOption("One of the 2 unfiltered half-maps:", NODE_HALFMAP_CPIPE, "", "MRC map files (*half1*.mrc)",  "Provide one of the two unfiltered half-reconstructions that were output upon convergence of a 3D auto-refine run.");
	joboptions["fn_mask"] = JobOption("Solvent mask:", NODE_MASK_CPIPE, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a soft mask where the protein is white (1) and the solvent is black (0). Often, the softer the mask the higher resolution estimates you will get. A soft edge of 5-10 pixels is often a good edge width.");
	joboptions["angpix"] = JobOption("Calibrated pixel size (A)", -1, 0.3, 5, 0.1, "Provide the final, calibrated pixel size in Angstroms. This value may be different from the pixel-size used thus far, e.g. when you have recalibrated the pixel size using the fit to a PDB model. The X-axis of the output FSC plot will use this calibrated value.");

	joboptions["do_auto_bfac"] = JobOption("Estimate B-factor automatically?", true, "If set to Yes, then the program will use the automated procedure described by Rosenthal and Henderson (2003, JMB) to estimate an overall B-factor for your map, and sharpen it accordingly. \
Note that your map must extend well beyond the lowest resolution included in the procedure below, which should not be set to resolutions much lower than 10 Angstroms. ");
	joboptions["autob_lowres"] = JobOption("Lowest resolution for auto-B fit (A):", 10, 8, 15, 0.5, "This is the lowest frequency (in Angstroms) that will be included in the linear fit of the Guinier plot as described in Rosenthal and Henderson (2003, JMB). Dont use values much lower or higher than 10 Angstroms. If your map does not extend beyond 10 Angstroms, then instead of the automated procedure use your own B-factor.");
	joboptions["do_adhoc_bfac"] = JobOption("Use your own B-factor?", false, "Instead of using the automated B-factor estimation, provide your own value. Use negative values for sharpening the map. \
This option is useful if your map does not extend beyond the 10A needed for the automated procedure, or when the automated procedure does not give a suitable value (e.g. in more disordered parts of the map).");
	joboptions["adhoc_bfac"] = JobOption("User-provided B-factor:", -1000, -2000, 0, -50, "Use negative values for sharpening. Be careful: if you over-sharpen your map, you may end up interpreting noise for signal!");

	joboptions["fn_mtf"] = JobOption("MTF of the detector (STAR file)", "", "STAR Files (*.star)", ".", "If you know the MTF of your detector, provide it here. Curves for some well-known detectors may be downloaded from the RELION Wiki. Also see there for the exact format \
\n If you do not know the MTF of your detector and do not want to measure it, then by leaving this entry empty, you include the MTF of your detector in your overall estimated B-factor upon sharpening the map.\
Although that is probably slightly less accurate, the overall quality of your map will probably not suffer very much.");
	joboptions["mtf_angpix"] = JobOption("Original detector pixel size:", 1.0, 0.3, 2.0, 0.1, "This is the original pixel size (in Angstroms) in the raw (non-super-resolution!) micrographs.");

	joboptions["do_skip_fsc_weighting"] = JobOption("Skip FSC-weighting?", false, "If set to No (the default), then the output map will be low-pass filtered according to the mask-corrected, gold-standard FSC-curve. \
Sometimes, it is also useful to provide an ad-hoc low-pass filter (option below), as due to local resolution variations some parts of the map may be better and other parts may be worse than the overall resolution as measured by the FSC. \
In such cases, set this option to Yes and provide an ad-hoc filter as described below.");
	joboptions["low_pass"] = JobOption("Ad-hoc low-pass filter (A):",5,1,40,1,"This option allows one to low-pass filter the map at a user-provided frequency (in Angstroms). When using a resolution that is higher than the gold-standard FSC-reported resolution, take care not to interpret noise in the map for signal...");
}

bool RelionJob::getCommandsPostprocessJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
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

	// Input half map (one of them)
	FileName fn_half1 = joboptions["fn_in"].getString();
	FileName fn_half2;

	if (is_tomo && joboptions["in_optimisation"].getString() != "")
	{
		FileName fn_OS = joboptions["in_optimisation"].getString();
		Node node(fn_OS, joboptions["in_optimisation"].node_type);
		inputNodes.push_back(node);
		command += " --ios " + fn_OS;

		Node node1(outputname + "postprocess_optimisation_set.star", LABEL_TOMO_OPTIMISATION);
		outputNodes.push_back(node1);
	}
	else if (fn_half1 == "")
	{
		error_message = "ERROR: empty field for input half-map...";
		return false;
	}

	if (fn_half1 != "")
	{
		if (!fn_half1.getTheOtherHalf(fn_half2))
		{
			error_message = "ERROR: cannot find 'half' substring in the input filename...";
			return false;
		}

		Node node(fn_half1, joboptions["fn_in"].node_type);
		inputNodes.push_back(node);
		command += " --i " + fn_half1;
	}

	// The output name contains a directory: use it for output
	command += " --o " + outputname + "postprocess";
	command += "  --angpix " + joboptions["angpix"].getString();
	Node node1(outputname+"postprocess.mrc", LABEL_POST_MAP);
	outputNodes.push_back(node1);
	Node node2(outputname+"postprocess_masked.mrc", LABEL_POST_MASKED);
	outputNodes.push_back(node2);

	Node node2b(outputname+"logfile.pdf", LABEL_POST_LOG);
	outputNodes.push_back(node2b);

	Node node2c(outputname+"postprocess.star", LABEL_POST);
	outputNodes.push_back(node2c);

	// Sharpening
	if (joboptions["fn_mtf"].getString().length() > 0)
	{
		command += " --mtf " + joboptions["fn_mtf"].getString();
		command += " --mtf_angpix " + joboptions["mtf_angpix"].getString();
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

	joboptions["fn_in"] = JobOption("One of the 2 unfiltered half-maps:", NODE_HALFMAP_CPIPE, "", "MRC map files (*half1*.mrc)",  "Provide one of the two unfiltered half-reconstructions that were output upon convergence of a 3D auto-refine run.");
	joboptions["angpix"] = JobOption("Calibrated pixel size (A)", 1, 0.3, 5, 0.1, "Provide the final, calibrated pixel size in Angstroms. This value may be different from the pixel-size used thus far, e.g. when you have recalibrated the pixel size using the fit to a PDB model. The X-axis of the output FSC plot will use this calibrated value.");

	// Check for environment variable RELION_RESMAP_TEMPLATE
	char *default_location = getenv("RELION_RESMAP_EXECUTABLE");
	char default_resmap[] = DEFAULTRESMAPLOCATION;
	if (default_location == NULL)
	{
		default_location = default_resmap;
	}

	joboptions["do_resmap_locres"] = JobOption("Use ResMap?", true, "If set to Yes, then ResMap will be used for local resolution estimation.");
	joboptions["fn_resmap"] = JobOption("ResMap executable:", std::string(default_location), "ResMap*", ".", "Location of the ResMap executable. You can control the default of this field by setting environment variable RELION_RESMAP_EXECUTABLE, or by editing the first few lines in src/gui_jobwindow.h and recompile the code. \n \n Note that the ResMap wrapper cannot use MPI.");
	joboptions["fn_mask"] = JobOption("User-provided solvent mask:", NODE_MASK_CPIPE, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a mask with values between 0 and 1 around all domains of the complex. ResMap uses this mask for local resolution calculation. RELION does NOT use this mask for calculation, but makes a histogram of local resolution within this mask.");
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
	initialisePipeline(outputname, job_counter);
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
	FileName fn_half1 = joboptions["fn_in"].getString();
	FileName fn_half2;
	if (!fn_half1.getTheOtherHalf(fn_half2))
	{
		error_message = "ERROR: cannot find 'half' substring in the input filename...";
		return false;
	}

	Node node(joboptions["fn_in"].getString(), joboptions["fn_in"].node_type);
	inputNodes.push_back(node);

	if (joboptions["do_resmap_locres"].getBoolean())
	{

		label += ".resmap";

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

		if (joboptions["nr_mpi"].getNumber(error_message) > 1)
		{
			error_message = "You cannot use more than 1 MPI processor for the ResMap wrapper...";
			return false;
		}
		if (error_message != "") return false;

		// Make symbolic links to the half-maps in the output directory
		commands.push_back("ln -s ../../" + fn_half1 + " " + outputname + "half1.mrc");
		commands.push_back("ln -s ../../" + fn_half2 + " " + outputname + "half2.mrc");

		Node node2(joboptions["fn_mask"].getString(), joboptions["fn_mask"].node_type);
		inputNodes.push_back(node2);

		Node node3(outputname + "half1_resmap.mrc", LABEL_LOCRES_RESMAP);
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
		label += ".own";

		if (joboptions["nr_mpi"].getNumber(error_message) > 1)
			command="`which relion_postprocess_mpi`";
		else
			command="`which relion_postprocess`";
		if (error_message != "") return false;

		command += " --locres --i " + joboptions["fn_in"].getString();
		command += " --o " + outputname + "relion";
		command += " --angpix " + joboptions["angpix"].getString();
		//command += " --locres_sampling " + joboptions["locres_sampling"].getString();
		//command += " --locres_randomize_at " + joboptions["randomize_at"].getString();
		command += " --adhoc_bfac " + joboptions["adhoc_bfac"].getString();
		if (joboptions["fn_mtf"].getString() != "")
			command += " --mtf " + joboptions["fn_mtf"].getString();

		if (joboptions["fn_mask"].getString() != "")
		{
			command += " --mask " + joboptions["fn_mask"].getString();
			Node node0(outputname+"histogram.pdf", LABEL_LOCRES_LOG);
			outputNodes.push_back(node0);
		}

		Node node1(outputname+"relion_locres_filtered.mrc", LABEL_LOCRES_FILTMAP);
		outputNodes.push_back(node1);
		Node node2(outputname+"relion_locres.mrc", LABEL_LOCRES_RESMAP);
		outputNodes.push_back(node2);
	}

	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();
	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}

void RelionJob::initialiseMotionrefineJob()
{
	hidden_name = ".gui_bayespolish";

	// I/O
	joboptions["fn_mic"] = JobOption("Micrographs (from MotionCorr):", NODE_MICS_CPIPE,  "", "STAR files (*.star)", "The input STAR file with the micrograph (and their movie metadata) from a MotionCorr job.");
	joboptions["fn_data"] = JobOption("Particles (from Refine3D or CtfRefine):", NODE_PARTS_CPIPE,  "", "STAR files (*.star)", "The input STAR file with the metadata of all particles.");
	joboptions["fn_post"] = JobOption("Postprocess STAR file:", OUTNODE_POST,  "", "STAR files (postprocess.star)", "The STAR file generated by a PostProcess job. \
The mask used for this postprocessing will be applied to the unfiltered half-maps and should encompass the entire complex. The resulting FSC curve will be used for weighting the different frequencies.");
	joboptions["do_float16"] = JobOption("Write output in float16?", true ,"If set to Yes, this program will write output images in float16 MRC format. This will save a factor of two in disk space compared to the default of writing in float32. Note that RELION and CCPEM will read float16 images, but other programs may not (yet) do so.");

	// Frame range
	joboptions["first_frame"] = JobOption("First movie frame: ", 1., 1., 10., 1, "First movie frame to take into account in motion fit and combination step");
	joboptions["last_frame"] = JobOption("Last movie frame: ", -1., 5., 50., 1, "Last movie frame to take into account in motion fit and combination step. Values equal to or smaller than 0 mean 'use all frames'.");

	joboptions["extract_size"] = JobOption("Extraction size (pix in unbinned movie):", -1, 64, 1024, 8, "Size of the extracted particles in the unbinned original movie(in pixels). This should be an even number.");
	joboptions["rescale"] = JobOption("Re-scaled size (pixels): ", -1, 64, 1024, 8, "The re-scaled value needs to be an even number.");

	// Parameter optimisation
	joboptions["do_param_optim"] = JobOption("Train optimal parameters?", false, "If set to Yes, then relion_motion_refine will estimate optimal parameter values for the three sigma values above on a subset of the data (determined by the minimum number of particles to be used below).");
	joboptions["eval_frac"] = JobOption("Fraction of Fourier pixels for testing: ", 0.5, 0, 1., 0.01, "This fraction of Fourier pixels (at higher resolution) will be used for evaluation of the parameters (test set), whereas the rest (at lower resolution) will be used for parameter estimation itself (work set).");
	joboptions["optim_min_part"] = JobOption("Use this many particles: ", 10000, 5000, 50000, 1000, "Use at least this many particles for the meta-parameter optimisation. The more particles the more expensive in time and computer memory the calculation becomes, but the better the results may get.");

	// motion_fit
	joboptions["do_polish"] = JobOption("Perform particle polishing?", true, "If set to Yes, then relion_motion_refine will be run to estimate per-particle motion-tracks using the parameters below, and polished particles will be generated.");
	joboptions["opt_params"] = JobOption("Optimised parameter file:", OUTNODE_POLISH_PARAMS,  "", "TXT files (*.txt)", "The output TXT file from a previous Bayesian polishing job in which the optimal parameters were determined.");
	joboptions["do_own_params"] = JobOption("OR use your own parameters?", false, "If set to Yes, then the field for the optimised parameter file will be ignored and the parameters specified below will be used instead.");
	joboptions["sigma_vel"] = JobOption("Sigma for velocity (A/dose): ", 0.2, 1., 10., 0.1, "Standard deviation for the velocity regularisation. Smaller values requires the tracks to be shorter.");
	joboptions["sigma_div"] = JobOption("Sigma for divergence (A): ", 5000, 0, 10000, 10000, "Standard deviation for the divergence of tracks across the micrograph. Smaller values requires the tracks to be spatially more uniform in a micrograph.");
	joboptions["sigma_acc"] = JobOption("Sigma for acceleration (A/dose): ", 2, -1, 7, 0.1, "Standard deviation for the acceleration regularisation. Smaller values requires the tracks to be straighter.");

	//combine_frames
	joboptions["minres"] = JobOption("Minimum resolution for B-factor fit (A): ", 20, 8, 40, 1, "The minimum spatial frequency (in Angstrom) used in the B-factor fit.");
	joboptions["maxres"] = JobOption("Maximum resolution for B-factor fit (A): ", -1, -1, 15, 1, "The maximum spatial frequency (in Angstrom) used in the B-factor fit. If a negative value is given, the maximum is determined from the input FSC curve.");
}

bool RelionJob::getCommandsMotionrefineJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;

	if (joboptions["nr_mpi"].getNumber(error_message) > 1)
		command="`which relion_motion_refine_mpi`";
	else
		command="`which relion_motion_refine`";
	if (error_message != "") return false;

	if (joboptions["fn_data"].getString() == "")
	{
		error_message = "ERROR: empty field for input particle STAR file...";
		return false;
	}
	if (joboptions["fn_mic"].getString() == "")
	{
		error_message = "ERROR: empty field for input micrograph STAR file...";
		return false;
	}
	if (joboptions["fn_post"].getString() == "")
	{
		error_message = "ERROR: empty field for input PostProcess STAR file...";
		return false;
	}

	if (joboptions["do_param_optim"].getBoolean() && joboptions["do_polish"].getBoolean())
	{
		error_message = "ERROR: Choose either parameter training or polishing, not both.";
		return false;
	}

	if (!joboptions["do_param_optim"].getBoolean() && !joboptions["do_polish"].getBoolean())
	{
		error_message = "ERROR: nothing to do, choose either parameter training or polishing.";
		return false;
	}

	if ((joboptions["eval_frac"].getNumber(error_message) <= 0.1 || joboptions["eval_frac"].getNumber(error_message) > 0.9 )
			&& !joboptions["eval_frac"].isSchedulerVariable() )
	{
		error_message = "ERROR: the fraction of Fourier pixels used for evaluation should be between 0.1 and 0.9.";
		return false;
	}
	if (error_message != "") return false;

	Node node(joboptions["fn_data"].getString(), joboptions["fn_data"].node_type);
	inputNodes.push_back(node);

	Node node2(joboptions["fn_post"].getString(), joboptions["fn_post"].node_type);
	inputNodes.push_back(node);

	command += " --i " + joboptions["fn_data"].getString();
	command += " --f " + joboptions["fn_post"].getString();
	command += " --corr_mic " + joboptions["fn_mic"].getString();
	command += " --first_frame " + joboptions["first_frame"].getString();
	command += " --last_frame " + joboptions["last_frame"].getString();
	command += " --o " + outputname;

	if (joboptions["do_float16"].getBoolean())
	{
		command += " --float16 ";
	}

	if (joboptions["do_param_optim"].getBoolean())
	{

		label += ".train";

		// Estimate meta-parameters
		RFLOAT align_frac = 1.0 - joboptions["eval_frac"].getNumber(error_message);
		if (error_message != "") return false;
		command += " --min_p " + joboptions["optim_min_part"].getString();
		command += " --eval_frac " + joboptions["eval_frac"].getString();
		command += " --align_frac " + floatToString(align_frac);

		if (joboptions["sigma_acc"].getNumber(error_message) < 0)
		{
			command += " --params2 ";
		}
		else
		{
			command += " --params3 ";
		}
		if (error_message != "") return false;

		Node node5(outputname+"opt_params_all_groups.txt", LABEL_POLISH_PARAMS);
		outputNodes.push_back(node5);
	}
	else if (joboptions["do_polish"].getBoolean())
	{
		if (joboptions["do_own_params"].getBoolean())
		{
			// User-specified Parameters
			command += " --s_vel " + joboptions["sigma_vel"].getString();
			command += " --s_div " + joboptions["sigma_div"].getString();
			command += " --s_acc " + joboptions["sigma_acc"].getString();
		}
		else
		{
			if (joboptions["opt_params"].getString() == "")
			{
				error_message = "ERROR: Please specify an optimised parameter file OR choose 'use own paramaeters' and set three sigma values.";
				return false;
			}
			command += " --params_file " + joboptions["opt_params"].getString();
		}

		command += " --combine_frames";
		command += " --bfac_minfreq " + joboptions["minres"].getString();
		command += " --bfac_maxfreq " + joboptions["maxres"].getString();

		const int window = ROUND(joboptions["extract_size"].getNumber(error_message));
		if (error_message != "") return false;

		const int scale = ROUND(joboptions["rescale"].getNumber(error_message));
		if (error_message != "") return false;

		if (window * scale <= 0)
		{
			error_message = "ERROR: Please specify both the extraction box size and the downsampled size, or leave both the default (-1)";
			return false;
		}

		if (window > 0 && scale > 0)
		{
			if (window % 2 != 0)
			{
				error_message = "ERROR: The extraction box size must be an even number";
				return false;
			}
			command += " --window " + joboptions["extract_size"].getString();

			if (scale % 2 != 0)
			{
				error_message = "ERROR: The downsampled box size must be an even number.";
				return false;
			}

			if (scale > window)
			{
				error_message = "ERROR: The downsampled box size cannot be larger than the extraction size.";
				return false;
			}
			command += " --scale " + joboptions["rescale"].getString();
		}

		Node node6(outputname+"logfile.pdf", LABEL_POLISH_LOG);
		outputNodes.push_back(node6);

		Node node7(outputname+"shiny.star", LABEL_POLISH_PARTS);
		outputNodes.push_back(node7);
	}

	// If this is a continue job, then only process unfinished micrographs
	if (is_continue)
		command += " --only_do_unfinished ";

	// Running stuff
	command += " --j " + joboptions["nr_threads"].getString();

	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();
	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}

void RelionJob::initialiseCtfrefineJob()
{
	hidden_name = ".gui_ctfrefine";

	// I/O
	joboptions["fn_data"] = JobOption("Particles (from Refine3D):", NODE_PARTS_CPIPE,  "", "STAR files (*.star)", "The input STAR file with the metadata of all particles.");
	joboptions["fn_post"] = JobOption("Postprocess STAR file:", OUTNODE_POST,  "", "STAR files (postprocess.star)", "The STAR file generated by a PostProcess job. \
The mask used for this postprocessing will be applied to the unfiltered half-maps and should encompass the entire complex. The resulting FSC curve will be used for weighting the different frequencies. \n \n Note that for helices it is common practice to use a mask only encompassing the central 30% or so of the box. \
This gives higher resolution estimates, as it disregards ill-defined regions near the box edges. However, for ctf_refine it is better to use a mask encompassing (almost) the entire box, as otherwise there may not be enough signal.");

	joboptions["minres"] = JobOption("Minimum resolution for fits (A): ", 30, 8, 40, 1, "The minimum spatial frequency (in Angstrom) used in the beamtilt fit.");

	// Defocus fit
	joboptions["do_ctf"] = JobOption("Perform CTF parameter fitting?", true, "If set to Yes, then relion_ctf_refine will be used to estimate the selected parameters below.");
	joboptions["do_defocus"] = JobOption("Fit defocus?", job_ctffit_options, 0, "If set to per-particle or per-micrograph, then relion_ctf_refine will estimate defocus values.");
	joboptions["do_astig"] = JobOption("Fit astigmatism?", job_ctffit_options, 0, "If set to per-particle or per-micrograph, then relion_ctf_refine will estimate astigmatism.");
	joboptions["do_bfactor"] = JobOption("Fit B-factor?", job_ctffit_options, 0, "If set to per-particle or per-micrograph, then relion_ctf_refine will estimate B-factors that describe the signal falloff.");
	joboptions["do_phase"] = JobOption("Fit phase-shift?", job_ctffit_options, 0, "If set to per-particle or per-micrograph, then relion_ctf_refine will estimate (VPP?) phase shift values.");

	// aberrations
	joboptions["do_aniso_mag"] = JobOption("Estimate (anisotropic) magnification?", false, "If set to Yes, then relion_ctf_refine will also estimate the (anisotropic) magnification per optics group. \
This option cannot be done simultaneously with higher-order aberration estimation. It's probably best to estimate the one that is most off first, and the other one second. It might be worth repeating the estimation if both are off.");

	joboptions["do_tilt"] = JobOption("Estimate beamtilt?", false, "If set to Yes, then relion_ctf_refine will also estimate the beamtilt per optics group. This option is only recommended for data sets that extend beyond 4.5 Angstrom resolution.");
	joboptions["do_trefoil"] = JobOption("Also estimate trefoil?", false, "If set to Yes, then relion_ctf_refine will also estimate the trefoil (3-fold astigmatism) per optics group. This option is only recommended for data sets that extend beyond 3.5 Angstrom resolution.");

	joboptions["do_4thorder"] = JobOption("Estimate 4th order aberrations?", false, "If set to Yes, then relion_ctf_refine will also estimate the Cs and the tetrafoil (4-fold astigmatism) per optics group. This option is only recommended for data sets that extend beyond 3 Angstrom resolution.");
}



bool RelionJob::getCommandsCtfrefineJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;


	if (joboptions["nr_mpi"].getNumber(error_message) > 1)
		command="`which relion_ctf_refine_mpi`";
	else
		command="`which relion_ctf_refine`";
	if (error_message != "") return false;

	if (joboptions["fn_data"].getString() == "")
	{
		error_message = "ERROR: empty field for input particle STAR file...";
		return false;
	}
	if (joboptions["fn_post"].getString() == "")
	{
		error_message = "ERROR: empty field for input PostProcess STAR file...";
		return false;
	}

	if (!joboptions["do_aniso_mag"].getBoolean() &&
	    !joboptions["do_ctf"].getBoolean() &&
	    !joboptions["do_tilt"].getBoolean() &&
	    !joboptions["do_4thorder"].getBoolean())
	{
		error_message = "ERROR: you haven't selected to fit anything...";
		return false;
	}

	if (!joboptions["do_aniso_mag"].getBoolean() && joboptions["do_ctf"].getBoolean() &&
	    joboptions["do_defocus"].getString() == job_ctffit_options[0] &&
	    joboptions["do_astig"].getString() == job_ctffit_options[0] &&
	    joboptions["do_bfactor"].getString() == job_ctffit_options[0] &&
	    joboptions["do_phase"].getString() == job_ctffit_options[0])
	{
		error_message = "ERROR: you did not select any CTF parameter to fit. Either switch off CTF parameter fitting, or select one to fit.";
		return false;
	}

	Node node(joboptions["fn_data"].getString(), joboptions["fn_data"].node_type);
	inputNodes.push_back(node);

	Node node2(joboptions["fn_post"].getString(), joboptions["fn_post"].node_type);
	inputNodes.push_back(node);

	Node node6(outputname+"logfile.pdf", LABEL_CTFREFINE_LOG);
	outputNodes.push_back(node6);

	command += " --i " + joboptions["fn_data"].getString();
	command += " --f " + joboptions["fn_post"].getString();
	command += " --o " + outputname;

	// Always either do anisotropic magnification, or CTF,tilt-odd,even
	if (joboptions["do_aniso_mag"].getBoolean())
	{
		label += ".anisomag";

		command += " --fit_aniso";
		command += " --kmin_mag " + joboptions["minres"].getString();

		Node node5(outputname+"particles_ctf_refine.star", LABEL_CTFREFINE_ANISOPARTS);
		outputNodes.push_back(node5);

	}
	else
	{
		Node node5(outputname+"particles_ctf_refine.star", LABEL_CTFREFINE_REFINEPARTS);
		outputNodes.push_back(node5);

		if (joboptions["do_ctf"].getBoolean())
		{
			command += " --fit_defocus --kmin_defocus " + joboptions["minres"].getString();
			std::string fit_options = "";

			fit_options += JobOption::getCtfFitString(joboptions["do_phase"].getString());
			fit_options += JobOption::getCtfFitString(joboptions["do_defocus"].getString());
			fit_options += JobOption::getCtfFitString(joboptions["do_astig"].getString());
			fit_options += "f"; // always have Cs refinement switched off
			fit_options += JobOption::getCtfFitString(joboptions["do_bfactor"].getString());

			if (fit_options.size() != 5)
			{
				error_message = "Wrong CTF fitting options";
				return false;
			}

			command += " --fit_mode " + fit_options;
		}

		// do not allow anisotropic magnification to be done simultaneously with higher-order aberrations
		if (joboptions["do_tilt"].getBoolean())
		{
			command += " --fit_beamtilt";
			command += " --kmin_tilt " + joboptions["minres"].getString();

			if (joboptions["do_trefoil"].getBoolean())
			{
				command += " --odd_aberr_max_n 3";
			}
		}

		if (joboptions["do_4thorder"].getBoolean())
		{
			command += " --fit_aberr";
		}
	}

	// If this is a continue job, then only process unfinished micrographs
	if (is_continue)
	{
		command += " --only_do_unfinished ";
	}

	// Running stuff
	command += " --j " + joboptions["nr_threads"].getString();

	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();
	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}


void RelionJob::initialiseExternalJob()
{
	hidden_name = ".gui_external";

	// I/O
	joboptions["fn_exe"] = JobOption("External executable:", "", "", ".", "Location of the script that will launch the external program. This script should write all its output in the directory specified with --o. Also, it should write in that same directory a file called RELION_JOB_EXIT_SUCCESS upon successful exit, and RELION_JOB_EXIT_FAILURE upon failure.");

	// Optional input nodes
	joboptions["in_mov"] = JobOption("Input movies: ", NODE_MOVIES_CPIPE, "", "movie STAR file (*.star)", "Input movies. This will be passed with a --in_movies argument to the executable.");
	joboptions["in_mic"] = JobOption("Input micrographs: ", NODE_MICS_CPIPE, "", "micrographs STAR file (*.star)", "Input micrographs. This will be passed with a --in_mics argument to the executable.");
	joboptions["in_part"] = JobOption("Input particles: ", NODE_PARTS_CPIPE, "", "particles STAR file (*.star)", "Input particles. This will be passed with a --in_parts argument to the executable.");
	joboptions["in_coords"] = JobOption("Input coordinates: ", NODE_COORDS_CPIPE, "", "STAR files (coords_suffix*.star)", "Input coordinates. This will be passed with a --in_coords argument to the executable.");
	joboptions["in_3dref"] = JobOption("Input 3D reference: ", NODE_MAP_CPIPE, "", "MRC files (*.mrc)", "Input 3D reference map. This will be passed with a --in_3dref argument to the executable.");
	joboptions["in_mask"] = JobOption("Input 3D mask: ", NODE_MASK_CPIPE, "", "MRC files (*.mrc)", "Input 3D mask. This will be passed with a --in_mask argument to the executable.");

	// Optional parameters
	joboptions["param1_label"] = JobOption("Param1 - label:", std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param1_value"] = JobOption("Param1 - value:" , std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param2_label"] = JobOption("Param2 - label:", std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param2_value"] = JobOption("Param2 - value:" , std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param3_label"] = JobOption("Param3 - label:", std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param3_value"] = JobOption("Param3 - value:" , std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param4_label"] = JobOption("Param4 - label:", std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param4_value"] = JobOption("Param4 - value:" , std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param5_label"] = JobOption("Param5 - label:", std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param5_value"] = JobOption("Param5 - value:" , std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param6_label"] = JobOption("Param6 - label:", std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param6_value"] = JobOption("Param6 - value:" , std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param7_label"] = JobOption("Param7 - label:", std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param7_value"] = JobOption("Param7 - value:" , std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param8_label"] = JobOption("Param8 - label:", std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param8_value"] = JobOption("Param8 - value:" , std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param9_label"] = JobOption("Param9 - label:", std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param9_value"] = JobOption("Param9 - value:" , std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param10_label"] = JobOption("Param10 - label:", std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
	joboptions["param10_value"] = JobOption("Param10 - value:" , std::string(""), "Define label and value for optional parameters to the script. These will be passed as an argument --label value");
}

bool RelionJob::getCommandsExternalJob(std::string &outputname, std::vector<std::string> &commands,
		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;

	if (joboptions["fn_exe"].getString() == "")
	{
		error_message = "ERROR: empty field for the external executable script...";
		return false;
	}

	command=joboptions["fn_exe"].getString();
	command += " --o " + outputname;

	// Optional input nodes
	if (joboptions["in_mov"].getString() != "")
	{
		Node node(joboptions["in_mov"].getString(), joboptions["in_mov"].node_type);
		inputNodes.push_back(node);
		command += " --in_movies " + joboptions["in_mov"].getString();
	}
	if (joboptions["in_mic"].getString() != "")
	{
		Node node(joboptions["in_mic"].getString(), joboptions["in_mic"].node_type);
		inputNodes.push_back(node);
		command += " --in_mics " + joboptions["in_mic"].getString();
	}
	if (joboptions["in_part"].getString() != "")
	{
		Node node(joboptions["in_part"].getString(), joboptions["in_part"].node_type);
		inputNodes.push_back(node);
		command += " --in_parts " + joboptions["in_part"].getString();
	}
	if (joboptions["in_coords"].getString() != "")
	{
		Node node(joboptions["in_coords"].getString(), joboptions["in_coords"].node_type);
		inputNodes.push_back(node);
		command += " --in_coords " + joboptions["in_coords"].getString();
	}
	if (joboptions["in_3dref"].getString() != "")
	{
		Node node(joboptions["in_3dref"].getString(), joboptions["in_3dref"].node_type);
		inputNodes.push_back(node);
		command += " --in_3dref " + joboptions["in_3dref"].getString();
	}
	if (joboptions["in_mask"].getString() != "")
	{
		Node node(joboptions["in_mask"].getString(), joboptions["in_mask"].node_type);
		inputNodes.push_back(node);
		command += " --in_mask " + joboptions["in_mask"].getString();
	}

	// Optional arguments
	if (joboptions["param1_label"].getString() != "")
	{
		command += " --" + joboptions["param1_label"].getString() + " " + joboptions["param1_value"].getString();
	}
	if (joboptions["param2_label"].getString() != "")
	{
		command += " --" + joboptions["param2_label"].getString() + " " + joboptions["param2_value"].getString();
	}
	if (joboptions["param3_label"].getString() != "")
	{
		command += " --" + joboptions["param3_label"].getString() + " " + joboptions["param3_value"].getString();
	}
	if (joboptions["param4_label"].getString() != "")
	{
		command += " --" + joboptions["param4_label"].getString() + " " + joboptions["param4_value"].getString();
	}
	if (joboptions["param5_label"].getString() != "")
	{
		command += " --" + joboptions["param5_label"].getString() + " " + joboptions["param5_value"].getString();
	}
	if (joboptions["param6_label"].getString() != "")
	{
		command += " --" + joboptions["param6_label"].getString() + " " + joboptions["param6_value"].getString();
	}
	if (joboptions["param7_label"].getString() != "")
	{
		command += " --" + joboptions["param7_label"].getString() + " " + joboptions["param7_value"].getString();
	}
	if (joboptions["param8_label"].getString() != "")
	{
		command += " --" + joboptions["param8_label"].getString() + " " + joboptions["param8_value"].getString();
	}
	if (joboptions["param9_label"].getString() != "")
	{
		command += " --" + joboptions["param9_label"].getString() + " " + joboptions["param9_value"].getString();
	}
	if (joboptions["param10_label"].getString() != "")
	{
		command += " --" + joboptions["param10_label"].getString() + " " + joboptions["param10_value"].getString();
	}

	// Running stuff
	command += " --j " + joboptions["nr_threads"].getString();

	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();
	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}

void RelionJob::addTomoInputOptions(bool has_tomograms, bool has_particles,
		bool has_trajectories, bool has_manifolds, bool has_halfmaps, bool has_postprocess)
{
    // Optional input nodes
     joboptions["in_optimisation"] = JobOption("Input optimisation set: ", OUTNODE_TOMO_OPTIMISATION, "", "Optimisation set STAR file (*.star)", "Input optimisation set. This will be passed with a --i argument to the executable. If any inidividual components of the optimisation set are specified below, then they will override the components in this optimisation set.");
     if (has_particles) joboptions["in_particles"] = JobOption("Input particle set: ", OUTNODE_TOMO_PARTS, "", "Particle STAR file (*.star)", "Input particle set. This will be passed with a --p argument to the executable. If specified, this will override the entry in the input optimisation set. If left empty, the entry from the optimisation set will be used.");
     if (has_tomograms) joboptions["in_tomograms"] = JobOption("Input tomogram set: ", OUTNODE_TOMO_TOMOGRAMS, "", "Tomogram set STAR file (*.star)", "Input tomogram set. This will be passed with a --m argument to the executable. If specified, this will override the entry in the input optimisation set. If left empty, the entry from the optimisation set will be used.");
     if (has_trajectories) joboptions["in_trajectories"] = JobOption("Input trajectory set: ", OUTNODE_TOMO_TRAJECTORIES, "", "Trajectory set STAR file (*.star)", "Input trajectory set. This will be passed with a --mot argument to the executable. If specified, this will override the entry in the input optimisation set. If left empty, the entry from the optimisation set will be used.");
     if (has_manifolds) joboptions["in_manifolds"] = JobOption("Input manifold set: ", OUTNODE_TOMO_MANIFOLDS, "", "Manifold set STAR file (*.star)", "Input manifold set. This will be passed with a --man argument to the executable. If specified, this will override the entry in the input optimisation set. If left empty, the entry from the optimisation set will be used.");
	 if (has_halfmaps) joboptions["in_halfmaps"] = JobOption("One of the 2 reference half-maps:", OUTNODE_TOMO_HALFMAP, "", "MRC map files (*half1*.mrc)", "Provide one of the two reference half-reconstructions. Both maps will be passed with a --ref1 and --ref2 arguments to the executable. If specified, this will override the entry in the input optimisation set. If left empty, the entry from the optimisation set will be used.");
	 if (has_postprocess)
	 {
	 	joboptions["in_refmask"] = JobOption("Reference mask: ", NODE_MASK_CPIPE, "", "Image Files (*.mrc)", "Input reference mask. This will be passed with a --mask argument to the executable. If specified, this will override the entry in the input optimisation set. If left empty, the entry from the optimisation set will be used.");
	 	joboptions["in_post"] = JobOption("Input postprocess STAR: ", OUTNODE_TOMO_POST, "", "Postprocess STAR file (postprocess.star)", "Input STAR file from a relion_postprocess job. This will be passed with a --fsc argument to the executable. If specified, this will override the entry in the input optimisation set. If left empty, the entry from the optimisation set will be used.");
	 }
}

std::string RelionJob::getTomoInputCommmand(std::string &command, int has_tomograms, int has_particles,
		int has_trajectories, int has_manifolds, bool has_halfmaps, int has_postprocess)
{
	std::string error_message = "";

	// if no optimisation set is given, check all other necessary files are present
	if (joboptions["in_optimisation"].getString() == "")
	{
		if (has_tomograms == HAS_COMPULSORY && joboptions["in_tomograms"].getString() == "")
		{
			error_message = "ERROR: no optimisation set is specified, yet also no tomogram set is specified";
			return error_message;
		}
		if (has_particles == HAS_COMPULSORY && joboptions["in_particles"].getString() == "")
		{
			error_message = "ERROR: no optimisation set is specified, yet also no particle set is specified";
			return error_message;
		}
		if (has_trajectories == HAS_COMPULSORY && joboptions["in_trajectories"].getString() == "")
		{
			error_message = "ERROR: no optimisation set is specified, yet also no trajectory set is specified";
			return error_message;
		}
		if (has_manifolds == HAS_COMPULSORY && joboptions["in_manifolds"].getString() == "")
		{
			error_message = "ERROR: no optimisation set is specified, yet also no manifold set is specified";
			return error_message;
		}
		if (has_halfmaps == HAS_COMPULSORY && joboptions["in_halfmaps"].getString() == "")
		{
			error_message = "ERROR: no optimisation set is specified, yet also no reference half map file is specified";
			return error_message;
		}
		if (has_postprocess == HAS_COMPULSORY)
		{
			if (joboptions["in_refmask"].getString() == "")
			{
				error_message = "ERROR: no optimisation set is specified, yet also no reference mask file is specified";
				return error_message;
			}
			if (joboptions["in_post"].getString() == "")
			{
				error_message = "ERROR: no optimisation set is specified, yet also no postprocess star file is specified";
				return error_message;
			}
		}
	}

	if (joboptions["in_optimisation"].getString() != "")
	{
		Node node(joboptions["in_optimisation"].getString(), joboptions["in_optimisation"].node_type);
		inputNodes.push_back(node);
    	command += " --i " + joboptions["in_optimisation"].getString();
	}
	if (has_tomograms != HAS_NOT && joboptions["in_tomograms"].getString() != "")
	{
		Node node(joboptions["in_tomograms"].getString(), joboptions["in_tomograms"].node_type);
		inputNodes.push_back(node);
    	command += " --t " + joboptions["in_tomograms"].getString();
	}
	if (has_particles != HAS_NOT && joboptions["in_particles"].getString() != "")
	{
		Node node(joboptions["in_particles"].getString(), joboptions["in_particles"].node_type);
		inputNodes.push_back(node);
    	command += " --p " + joboptions["in_particles"].getString();
	}
	if (has_trajectories != HAS_NOT && joboptions["in_trajectories"].getString() != "")
	{
		Node node(joboptions["in_trajectories"].getString(), joboptions["in_trajectories"].node_type);
		inputNodes.push_back(node);
    	command += " --mot " + joboptions["in_trajectories"].getString();
	}
	if (has_manifolds != HAS_NOT && joboptions["in_manifolds"].getString() != "")
	{
		Node node(joboptions["in_manifolds"].getString(), joboptions["in_manifolds"].node_type);
		inputNodes.push_back(node);
    	command += " --man " + joboptions["in_manifolds"].getString();
	}
	if (has_halfmaps != HAS_NOT && joboptions["in_halfmaps"].getString() != "")
	{
		// Input half map (one of them)
		FileName fn_half1 = joboptions["in_halfmaps"].getString();
		FileName fn_half2;
		if (!fn_half1.getTheOtherHalf(fn_half2))
		{
			error_message = "ERROR: cannot find 'half' substring in the halfmap filename...";
			return error_message;
		}
		Node node(fn_half1, joboptions["in_halfmaps"].node_type);
		inputNodes.push_back(node);
		command += " --ref1 " + fn_half1;
		command += " --ref2 " + fn_half2;
	}
	if (has_postprocess != HAS_NOT && joboptions["in_refmask"].getString() != "")
	{
		Node node(joboptions["in_refmask"].getString(), joboptions["in_refmask"].node_type);
		inputNodes.push_back(node);
		command += " --mask " + joboptions["in_refmask"].getString();
	}
	if (has_postprocess != HAS_NOT && joboptions["in_post"].getString() != "")
	{
		Node node(joboptions["in_post"].getString(), joboptions["in_post"].node_type);
		inputNodes.push_back(node);
		command += " --fsc " + joboptions["in_post"].getString();
	}

	return error_message;
}

std::string RelionJob::setTomoOutputCommand(std::string &command, std::string optimisationSet,	std::string tomograms,
								std::string particles, std::string trajectories, std::string manifolds,
								std::string halfmap1, std::string postprocess, std::string refmask,
								std::string optimisationSetOut)
{
	std::string error_message = "";

	// Create output optimisation set
	command = "`which relion_tomo_make_optimisation_set`";
	command += " --o " + optimisationSetOut;

	if (optimisationSet != "") command += " --i " + optimisationSet;
	if (tomograms != "") command += " --t " + tomograms;
	if (particles != "") command += " --p " + particles;
	if (trajectories != "") command += " --mot " + trajectories;
	if (manifolds != "") command += " --man " + manifolds;
	if (halfmap1 != "")
	{
		FileName fn_half1 = halfmap1;
		FileName halfmap2;
		if (!fn_half1.getTheOtherHalf(halfmap2))
		{
			error_message = "ERROR: cannot find 'half' substring in the input filename...";
			return error_message;
		}
		command += " --ref1 " + halfmap1;
		command += " --ref2 " + halfmap2;
	}
	if (postprocess != "") command += " --fsc " + postprocess;
	if (refmask != "") command += " --mask " + refmask;

	Node node1(optimisationSetOut, LABEL_TOMO_OPTIMISATION);
	outputNodes.push_back(node1);

	return error_message;
}


void RelionJob::initialiseTomoImportJob()
{
        hidden_name = ".gui_tomo_import";

       	joboptions["do_tomo"] = JobOption("Import tomograms?", true, "Set this to Yes for importing tomogram directories from IMOD.");
        joboptions["io_tomos"] = JobOption("Append to tomograms set: ", OUTNODE_TOMO_TOMOGRAMS, "", "Tomogram set STAR file (*.star)", "The imported tomograms will be output into this tomogram set. If any tomograms were already in this tomogram set, then the newly imported ones will be added to those.");
        joboptions["tomo_star"] = JobOption("STAR file with tomograms description: ", "", "Input file (*.star)", ".", "Provide a STAR file with the basic following information to import tomogsrams: \n\n"
                  " - rlnTomoImportImodDir: path to the IMOD directory.\n"
                  " - rlnTomoImportCtfFindFile or rlnTomoImportCtfPlotterFile: path to the initial CTF estimate from either CTFFind or CtfPlotter, respectively.\n"
                  " - rlnTomoTiltSeriesName: path to the actual tilt series file. Note if the filename ends with .st, this needs to be specificed with an .st:mrc ending to tell RELION  to interpret it as an mrc file.\n\n"
                  "The following additional columns may also be present:\n\n"
                  " - rlnTomoName: The tomogram name. If not specified rlnTomoTiltSeriesName will be used instead.\n"
                  " - rlnTomoImportFractionalDose: the electron dose corresponding to one tilt image. If omitted, the value of the --fd argument will be used.\n"
                  " - rlnTomoImportOrderList: path to a two-column csv text file specifying the chronological order in which the images were acquired. The first number counts up from 1, while the second describes the sequence of tilt angles.\n"
                  " - rlnOpticsGroupName: an arbitrary name for an optics group. This allows the set of tilt series to be separated into subsets that share the same optical aberrations. This is useful if the data have been collected in multiple sessions that might exhibit different aberrations. If omitted, all tilt series will be assigned to the same default optics group.\n"
                  " - rlnTomoImportOffset<X/Y/Z>: an arbitrary offset to the 3D coordinate system. This is useful if particles have already been picked in tomograms that have been cropped after reconstruction by IMOD. If the IMOD-internal SHIFT command has been used to apply offsets, then this will be handled internally and does not need to be specified here. If omitted, then the values of the --off<x/y/z> command line arguments will be used instead (which default to 0).\n"
                  " - rlnTomoImportCulledFile: output file name for a new tilt series with the excluded frames missing. This is only needed if tilt images have been excluded using IMOD’s EXCLUDE, EXCLUDELIST or EXCLUDELIST2 commands. In that case, this becomes a mandatory parameter.");
    	joboptions["angpix"] = JobOption("Pixel size (Angstrom):", (std::string)"", "Pixel size in Angstroms. If this values varies among the input tomograms, then specify it using its own column (rlnTomoTiltSeriesPixelSize) in the input tomogram description STAR file.");
    	joboptions["kV"] = JobOption("Voltage (kV):", (std::string)"", "Voltage the microscope was operated on (in kV; default=300). If this values varies among the input tomograms, then specify it using its own column (rlnVoltage) in the input tomogram description STAR file.");
    	joboptions["Cs"] = JobOption("Spherical aberration (mm):", (std::string)"", "Spherical aberration of the microscope used to collect these images (in mm; default=2.7). Typical values are 2.7 (FEI Titan & Talos, most JEOL CRYO-ARM), 2.0 (FEI Polara), 1.4 (some JEOL CRYO-ARM) and 0.01 (microscopes with a Cs corrector). If this values varies among the input tomograms, then specify it using its own column (rlnSphericalAberration) in the input tomogram description STAR file.");
    	joboptions["Q0"] = JobOption("Amplitude contrast:", (std::string)"", "Fraction of amplitude contrast (default=0.1). Often values around 10% work better than theoretically more accurate lower values.  If this values varies among the input tomograms, then specify it using its own column (rlnAmplitudeContrast) in the input tomogram description STAR file.");
    	joboptions["dose"] = JobOption("Frame dose:", (std::string)"", "Electron dose (in e/A^2) per frame (image) in the tilt series.  If this values varies among the input tomograms, then specify it using its own column (rlnTomoImportFractionalDose) in the input tomogram description STAR file.");
    	joboptions["order_list"] = JobOption("Ordered list:", (std::string)"", "", ".", "A 2-column, comma-separated file with the frame-order list of the tilt series, where the first column is the frame (image) number (starting at 1) and the second column is the tilt angle (in degrees). If this values varies among the input tomograms, then specify it using its own column (rlnTomoImportOrderList) in the input tomogram description STAR file.");
    	joboptions["do_flipYZ"] = JobOption("Flip YZ?", true, "Set this to Yes if you want to interchange the Y and Z coordinates.  If this values varies among the input tomograms, then append opposite values to tomogram set using another Import tomo job.");
    	joboptions["do_flipZ"] = JobOption("Flip Z?", true, "Set this to Yes if you want to change the sign of the Z coordinates.  If this values varies among the input tomograms, then append opposite values to tomogram set using another Import tomo job.");
    	joboptions["hand"] = JobOption("Tilt handedness:", (std::string)"", "Set this to indicate the handedness of the tilt geometry (default=-1). The value of this parameter is either +1 or -1, and it describes whether the focus increases or decreases as a function of Z distance. It has to be determined experimentally. In our experiments, it has always been -1. Y If this values varies among the input tomograms, then append opposite values to tomogram set using another Import tomo job.");
       	joboptions["do_coords"] = JobOption("Import coordinates?", false, "Set this to Yes for importing particle coordinates.");
        joboptions["part_star"] = JobOption("STAR file with coordinates: ", "", "Input file (*.star)", ".", "Provide a STAR file with the following information to input particles: \n \n TODO TODO TODO ");
        joboptions["part_tomos"] = JobOption("Tomograms set: ", OUTNODE_TOMO_TOMOGRAMS, "", "Tomogram set STAR file (*.star)", "The tomograms set from which these particles were picked.");
        joboptions["do_coords_flipZ"] = JobOption("Flip Z coordinates?", false, "Set this to Yes if you want to flip particles Z coordinate. Use it in case imported tomograms Z axis are flipped compared to tomograms used for picking.");
    	joboptions["do_other"] = JobOption("Import other node types?", false, "Set this to Yes  if you plan to import anything else than movies or micrographs");
    	joboptions["fn_in_other"] = JobOption("Input file:", "ref.mrc", "Input file (*.*)", ".", "Select any file(s) to import. \n \n \
    Note that for importing coordinate files, one has to give a Linux wildcard, where the *-symbol is before the coordinate-file suffix, e.g. if the micrographs are called mic1.mrc and the coordinate files mic1.box or mic1_autopick.star, one HAS to give '*.box' or '*_autopick.star', respectively.\n \n \
    Also note that micrographs, movies and coordinate files all need to be in the same directory (with the same rootnames, e.g.mic1 in the example above) in order to be imported correctly. 3D masks or references can be imported from anywhere. \n \n \
    Note that movie-particle STAR files cannot be imported from a previous version of RELION, as the way movies are handled has changed in RELION-2.0. \n \n \
    For the import of a particle, 2D references or micrograph STAR file or of a 3D reference or mask, only a single file can be imported at a time. \n \n \
    Note that due to a bug in a fltk library, you cannot import from directories that contain a substring  of the current directory, e.g. dont important from /home/betagal if your current directory is called /home/betagal_r2. In this case, just change one of the directory names.");

    	joboptions["node_type"] = JobOption("Node type:", job_nodetype_options_tomo, 0, "Select the type of Node this is.");
    	joboptions["optics_group_particles"] = JobOption("Rename optics group for particles:", (std::string)"", "Only for the import of a particles STAR file with a single, or no, optics groups defined: rename the optics group for the imported particles to this string.");



}

bool RelionJob::getCommandsTomoImportJob(std::string &outputname, std::vector<std::string> &commands,
		 std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
    initialisePipeline(outputname, job_counter);
    std::string command;

	// Some code here was copied from the SPA import job...
	bool do_tomo = joboptions["do_tomo"].getBoolean();
	bool do_coords = joboptions["do_coords"].getBoolean();
	bool do_other = joboptions["do_other"].getBoolean();

	int i = 0;
	if (do_tomo) i++;
	if (do_coords) i++;
	if (do_other) i++;

	if (i != 1)
	{
		error_message = "ERROR: you can only select to import tomograms, import particles, OR import other nodes.";
	return false;
	}

	if (do_tomo)
	{

		if (joboptions["tomo_star"].getString() == "")
		{
			error_message = "ERROR: you need to provide an input STAR file with information about the tomograms to be imported";
			return false;
		}

		// TODO: insert call to relion_tomo_import_tomograms here
		command = "relion_tomo_import_tomograms ";

		command += " --i " + joboptions["tomo_star"].getString();
		command += " --o " + outputname+"tomograms.star";
                if (joboptions["io_tomos"].getString() != "") command += " --t " + joboptions["io_tomos"].getString();

		Node node(outputname+"tomograms.star", LABEL_TOMO_TOMOGRAMS);
		outputNodes.push_back(node);

		if (joboptions["angpix"].getString() != "") command += " --angpix " + joboptions["angpix"].getString();
		if (joboptions["kV"].getString() != "") command += " --voltage " + joboptions["kV"].getString();
		if (joboptions["Cs"].getString() != "") command += " --Cs " + joboptions["Cs"].getString();
		if (joboptions["Q0"].getString() != "") command += " --Q0 " + joboptions["Q0"].getString();
		if (joboptions["dose"].getString() != "") command += " --fd " + joboptions["dose"].getString();
		if (joboptions["order_list"].getString() != "") command += " --ol " + joboptions["order_list"].getString();
		if (joboptions["do_flipYZ"].getBoolean()) command += " --flipYZ ";
		if (joboptions["do_flipZ"].getBoolean()) command += " --flipZ ";
		if (joboptions["hand"].getString() != "") command += " --hand " + joboptions["hand"].getString();


	}
	else if (do_coords)
	{

		if (joboptions["part_star"].getString() == "")
		{
			error_message = "ERROR: you need to provide an input STAR file with information about the tomograms to be imported.";
			return false;
		}

		if (joboptions["part_tomos"].getString() == "")
		{
			error_message = "ERROR: you need to provide an input tomograms set with information about the tomograms from which they particles originate.";
			return false;
		}

		command = "relion_tomo_import_particles ";

		command += " --i " + joboptions["part_star"].getString();
		command += " --o " + outputname;
		command += " --t " + joboptions["part_tomos"].getString();

		if (joboptions["do_coords_flipZ"].getBoolean())
		{
			command += " --flipZ";
		}

		Node node(outputname+"particles.star", LABEL_TOMO_PARTS);
		outputNodes.push_back(node);
		Node node2(outputname+"optimisation_set.star", LABEL_TOMO_OPTIMISATION);
		outputNodes.push_back(node2);
	}
	else if (do_other)
	{
		FileName fn_out, fn_in;
		command = "relion_import ";

		fn_in = joboptions["fn_in_other"].getString();
		std::string node_type = joboptions["node_type"].getString();

		fn_out = "/" + fn_in;
		fn_out = fn_out.afterLastOf("/");

		std::string mynodetype;
		if (node_type == "Particles STAR file (.star)")
			mynodetype = LABEL_TOMO_PARTS;
		else if (node_type == "Set of tomograms STAR file (.star)")
			mynodetype = LABEL_TOMO_TOMOGRAMS;
		else if (node_type == "Multiple (2D or 3D) references (.star or .mrcs)")
			mynodetype = LABEL_2DIMGS_CPIPE;
		else if (node_type == "3D reference (.mrc)")
			mynodetype = LABEL_MAP_CPIPE;
		else if (node_type == "3D mask (.mrc)")
			mynodetype = LABEL_MASK_CPIPE;
		else if (node_type == "Unfiltered half-map (unfil.mrc)")
			mynodetype = LABEL_TOMO_HALFMAP;
		else
		{
			error_message = "Unrecognized menu option for node_type = " + node_type;
			return false;
		}

		Node node(outputname + fn_out, mynodetype);
		outputNodes.push_back(node);

		// Also get the other half-map
		if (mynodetype == LABEL_TOMO_HALFMAP)
		{
			FileName fn_inb = "/" + fn_in;
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
			fn_inb = fn_inb.afterLastOf("/");
			Node node2(outputname + fn_inb, mynodetype);
			outputNodes.push_back(node2);
			command += " --do_halfmaps";
		}
		else if (mynodetype == LABEL_TOMO_PARTS)
		{
			command += " --do_particles";
			FileName optics_group = joboptions["optics_group_particles"].getString();
			if (optics_group != "")
			{
				if (!optics_group.validateCharactersStrict())
				{
					error_message = "ERROR: an optics group name may contain only numbers, alphabets and hyphen(-).";
					return false;
				}
				command += " --particles_optics_group_name \"" + optics_group + "\"";
			}
		}
		else
		{
			command += " --do_other";
		}

		// Now finish the command call to relion_import program, which does the actual copying
		command += " --i \"" + fn_in + "\"";
		command += " --odir " + outputname;
		command += " --ofile " + fn_out;

	}


    // Other arguments for extraction
    command += " " + joboptions["other_args"].getString();
    commands.push_back(command);

    return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

void RelionJob::initialiseTomoSubtomoJob()
{

	hidden_name = ".gui_tomo_subtomo";

	addTomoInputOptions(true, true, true, false, false, false);

	joboptions["box_size"] = JobOption("Box size (pix):", 128, 32, 512, 16, "The initial box size of the reconstruction. A sufficiently large box size allows more of the high-frequency signal to be captured that has been delocalised by the CTF.");
	joboptions["crop_size"] = JobOption("Cropped box size (pix):", -1, -1, 512, 16, "If set to a positive value, after construction, the resulting pseudo subtomograms are cropped to this size. A smaller box size allows the (generally expensive) refinement using relion_refine to proceed more rapidly.");
	joboptions["binning"] = JobOption("Binning factor:", 1, 1, 16, 1, "The tilt series images will be binned by this (real-valued) factor and then reconstructed in the specified box size above. Note that thereby the reconstructed region becomes larger when specifying binning factors larger than one.");

	joboptions["do_cone_weight"] = JobOption("Use cone weight?", false, "If set to Yes, then downweight a cone in Fourier space along the Z axis (as defined by the coordinate system of the particle). This is useful for particles embedded in a membrane, as it can prevent the alignment from being driven by the membrane signal (the signal of a planar membrane is localised within one line in 3D Fourier space). Note that the coordinate system of a particle is given by both the subtomogram orientation (if defined) and the particle orientation (see particle set). This allows the user to first obtain a membrane-driven alignment, and to then specifically suppress the signal in that direction.");
	joboptions["cone_angle"] = JobOption("Cone angle:", 10, 1, 50, 1, "The (full) opening angle of the cone to be suppressed, given in degrees. This angle should include both the uncertainty about the membrane orientation and its variation across the region represented in the subtomogram.");

	joboptions["do_float16"] = JobOption("Write output in float16?", true ,"If set to Yes, this program will write output images in float16 MRC format. This will save a factor of two in disk space compared to the default of writing in float32. Note that RELION and CCPEM will read float16 images, but other programs may not (yet) do so.");

}



bool RelionJob::getCommandsTomoSubtomoJob(std::string &outputname, std::vector<std::string> &commands,
				std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{
	commands.clear();
	initialisePipeline(outputname, job_counter);
	std::string command;

	if (joboptions["nr_mpi"].getNumber(error_message) > 1)
		command="`which relion_tomo_subtomo_mpi`";
	else
		command="`which relion_tomo_subtomo`";
	if (error_message != "") return false;

	// I/O
	error_message = getTomoInputCommmand(command, HAS_COMPULSORY, HAS_COMPULSORY, HAS_OPTIONAL, HAS_NOT, HAS_NOT,
										 HAS_NOT);
	if (error_message != "") return false;

	command += " --theme classic --o " + outputname;

	Node node1(outputname+"optimisation_set.star", LABEL_TOMO_OPTIMISATION);
	outputNodes.push_back(node1);
	Node node2(outputname+"particles.star", LABEL_TOMO_PARTS);
	outputNodes.push_back(node2);

	// Job-specific stuff goes here

	command += " --b " + joboptions["box_size"].getString();

	int crop_size = joboptions["crop_size"].getNumber(error_message);
    if (error_message != "") return false;
    if (crop_size > 0.) command += " --crop " + joboptions["crop_size"].getString();

    command += " --bin " + joboptions["binning"].getString();

    if (joboptions["do_cone_weight"].getBoolean())
    {
    	command += " --cone_weight --cone_angle " + joboptions["cone_angle"].getString();
    }

	if (joboptions["do_float16"].getBoolean())
	{
		command += " --float16 ";
	}

	if (is_continue)
	{
		command += " --only_do_unfinished ";
	}

	// Running stuff
	command += " --j " + joboptions["nr_threads"].getString();

	// Other arguments for extraction
	command += " " + joboptions["other_args"].getString();
	commands.push_back(command);

	return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);
}


void RelionJob::initialiseTomoCtfRefineJob()
{
        hidden_name = ".gui_tomo_refine_ctf";

	addTomoInputOptions(true, true, true, false, true, true);

    	joboptions["box_size"] = JobOption("Box size for estimation (pix):", 128, 32, 512, 16, "Box size to be used for the estimation. Note that this can be larger than the box size of the reference map. A sufficiently large box size allows more of the high-frequency signal to be captured that has been delocalised by the CTF.");
    	joboptions["do_defocus"] = JobOption("Refine defocus?", true, "If set to Yes, then estimate the defoci of the individual tilt images.");
    	joboptions["focus_range"] = JobOption("Defocus search range (A):", 3000, 0, 10000, 500, "Defocus search range (in A). This search range will be, by default, sampled in 100 steps. Use the additional argument --ds to change the number of sampling points.");
     	joboptions["do_reg_def"] = JobOption("Do defocus regularisation?", false, "Apply defocus regularisation. " \
		"High-tilt images do not offer enough signal to recover the defocus value precisely. The regularisation " \
		"forces the estimated defoci to assume similar values within a given tilt series, which prevents those " \
		"high-tilt images from overfitting.");
    	joboptions["lambda"] = JobOption("Defocus regularisation lambda:", 0.1, 0, 1, 0.05, "Defocus regularisation scale. ");
    	joboptions["do_scale"] = JobOption("Refine contrast scale?", true, "If set to Yes, then estimate the signal " \
    	"scale or ice thickness.");
    	joboptions["do_frame_scale"] = JobOption("Refine scale per frame?", true, "If set to Yes, then estimate the " \
		"signal-scale parameter independently for each tilt. If not specified, the ice thickness, beam luminance and " \
		"surface normal are estimated instead. Those three parameters then imply the signal intensity for each frame. "\
		"Due to the smaller number of parameters, the ice thickness model is more robust to noise. By default, the "\
		"ice thickness and surface normal will be estimated per tilt-series, and the beam luminance globally.");
    	joboptions["do_tomo_scale"] = JobOption("Refine scale per tomogram?", false, "If set to Yes, then estimate "\
    	"the beam luminance separately for each tilt series. This is not recommended.");

    	joboptions["do_even_aberr"] = JobOption("Refine even aberrations?", true, "If set to Yes, then estimates the even higher-order aberrations.");
    	joboptions["nr_even_aberr"] = JobOption("Order of even aberrations:", 4, 4, 8, 2, "The maximum order for the even aberrations to be estimated.");
    	joboptions["do_odd_aberr"] = JobOption("Refine odd aberrations?", true, "If set to Yes, then estimates the odd higher-order aberrations.");
    	joboptions["nr_odd_aberr"] = JobOption("Order of odd aberrations:", 3, 3, 7, 2, "The maximum order for the odd aberrations to be estimated.");

}

bool RelionJob::getCommandsTomoCtfRefineJob(std::string &outputname, std::vector<std::string> &commands,
				std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
    initialisePipeline(outputname, job_counter);
    std::string command;

    if (joboptions["nr_mpi"].getNumber(error_message) > 1)
            command="`which relion_tomo_refine_ctf_mpi`";
    else
            command="`which relion_tomo_refine_ctf`";
    if (error_message != "") return false;

    // I/O
    error_message = getTomoInputCommmand(command, HAS_COMPULSORY, HAS_COMPULSORY, HAS_OPTIONAL, HAS_NOT, HAS_COMPULSORY,
										 HAS_OPTIONAL);
	if (error_message != "") return false;

	command += " --theme classic --o " + outputname;

	Node node1(outputname+"optimisation_set.star", LABEL_TOMO_OPTIMISATION);
	outputNodes.push_back(node1);
	Node node2(outputname+"tomograms.star", LABEL_TOMO_TOMOGRAMS);
	outputNodes.push_back(node2);
	Node node3(outputname + "logfile.pdf", LABEL_TOMO_CTFREFINE_LOG);
	outputNodes.push_back(node3);

	// Job-specific stuff goes here

	command += " --b " + joboptions["box_size"].getString();

	if (joboptions["do_defocus"].getBoolean())
	{
		command += " --do_defocus";
		command += " --d0 -" + joboptions["focus_range"].getString();
		command += " --d1 " + joboptions["focus_range"].getString();

		if (joboptions["do_reg_def"].getBoolean())
			command += " --do_reg_defocus --lambda " + joboptions["lambda"].getString();
	}

	if (joboptions["do_scale"].getBoolean())
	{
		command += " --do_scale";
		if (joboptions["do_frame_scale"].getBoolean() && joboptions["do_tomo_scale"].getBoolean())
		{
			error_message = "ERROR: per-tomogram scale estimation and per-frame scale estimation are mutually exclusive";
			return false;
		}
		if (joboptions["do_frame_scale"].getBoolean()) command += " --per_frame_scale";
		if (joboptions["do_tomo_scale"].getBoolean()) command += " --per_tomo_scale";

	}

	if (joboptions["do_even_aberr"].getBoolean())
	{
		command += " --do_even_aberrations --ne " + joboptions["nr_even_aberr"].getString();
	}

	if (joboptions["do_odd_aberr"].getBoolean())
	{
		command += " --do_odd_aberrations --no " + joboptions["nr_odd_aberr"].getString();
	}

	if (is_continue)
	{
		command += " --only_do_unfinished ";
	}

	// Running stuff
    command += " --j " + joboptions["nr_threads"].getString();

    // Other arguments for extraction
    command += " " + joboptions["other_args"].getString();
    commands.push_back(command);

    return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

void RelionJob::initialiseTomoAlignJob()
{
	hidden_name = ".gui_tomo_align";

	addTomoInputOptions(true, true, true, false, true, true);

	joboptions["box_size"] = JobOption("Box size for estimation (pix):", 128, 32, 512, 16, "Box size to be used for the estimation. Note that this can be larger than the box size of the reference map. A sufficiently large box size allows more of the high-frequency signal to be captured that has been delocalised by the CTF.");
	joboptions["max_error"] = JobOption("Max position error (pix):", 5, 0, 64, 1, "maximal assumed error in the initial 2D particle-positions (distances between the projected 3D positions and their true positions in the images), given in pixels.");

	joboptions["do_shift_align"] = JobOption("Align by shift only?", false, "If set to Yes, tilt series projection shifts are refined based on cross-correlation. Useful for very badly aligned frames. No iterative optimisation.");
	joboptions["shift_align_type"] = JobOption("Alignment model: ", job_tomo_align_shiftonly_options, 0, "If set to \"Only particles\", it estimates rigid shift by aligning only the particles instead of by predicting entire micrographs. In this case, only misalignments smaller than half the box size of the particle can be corrected.");

	joboptions["do_motion"] = JobOption("Fit per-particle motion?", false, "If set to Yes, then the subtomogram version of Bayesian polishing will be used to fit per-particle (3D) motion tracks, besides the rigid part of the motion in the tilt series.");
	joboptions["sigma_vel"] = JobOption("Sigma for velocity (A/dose): ", 0.2, 1., 10., 0.1, "The expected amount of motion (i.e. the std. deviation of particle positions in Angstroms after 1 electron per A^2 of radiation)");
	joboptions["sigma_div"] = JobOption("Sigma for divergence (A): ", 5000, 0, 10000, 10000, "The expected spatial smoothness of the particle trajectories in A (a greater value means spatially smoother motion");
	joboptions["do_sq_exp_ker"] = JobOption("Use Gaussian decay?", false, "If set to Yes, then assume that the correlation of the velocities of two particles decays as a Gaussian over their distance, instead of as an exponential. This will produce spatially smoother motion and result in a shorter program runtime.");

	joboptions["do_deform"] = JobOption("Estimate 2D deformations?", false, "If set to Yes, then the subtomogram version of Bayesian polishing will be used to fit per-particle (3D) motion tracks, besides the rigid part of the motion in the tilt series.");
	joboptions["def_w"] = JobOption("Horizontal sampling points: ", 3, 0, 10, 1, "Number of horizontal sampling points for the deformation grid.");
	joboptions["def_h"] = JobOption("Vertical sampling points: ", 3, 0, 10, 1, "Number of vertical sampling points for the deformation grid.");
	joboptions["def_model"] = JobOption("Deformation Model:", job_tomo_align_def_model, 1, "Type of model to use (linear, spline or Fourier).");
	joboptions["lambda"] = JobOption("Deformation regularisation lambda:", 0., 0, 1, 0.05, "Deformation regularisation scale.");
	joboptions["do_frame_def"] = JobOption("Refine deformations per frame?", false, "If set to Yes, it models deformations per tilt frame instead of per tilt series.");

}

bool RelionJob::getCommandsTomoAlignJob(std::string &outputname, std::vector<std::string> &commands,
				std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
    initialisePipeline(outputname, job_counter);
    std::string command;

	if (joboptions["nr_mpi"].getNumber(error_message) > 1)
		command="`which relion_tomo_align_mpi`";
	else
		command="`which relion_tomo_align`";
	if (error_message != "") return false;

	// I/O
    error_message = getTomoInputCommmand(command, HAS_COMPULSORY, HAS_COMPULSORY, HAS_OPTIONAL, HAS_NOT, HAS_COMPULSORY,
					HAS_OPTIONAL);
    if (error_message != "") return false;

	command += " --theme classic --o " + outputname;

	Node node1(outputname+"optimisation_set.star", LABEL_TOMO_OPTIMISATION);
	outputNodes.push_back(node1);
	Node node2(outputname+"tomograms.star", LABEL_TOMO_TOMOGRAMS);
	outputNodes.push_back(node2);
	Node node3(outputname+"particles.star", LABEL_TOMO_PARTS);
	outputNodes.push_back(node3);
	if (joboptions["do_motion"].getBoolean())
	{
		Node node4(outputname+"motion.star", LABEL_TOMO_TRAJECTORIES);
		outputNodes.push_back(node4);
	}
	Node node5(outputname + "logfile.pdf", LABEL_TOMO_FRAMEALIGN_LOG);
	outputNodes.push_back(node5);

	// Job-specific stuff goes here
	command += " --b " + joboptions["box_size"].getString();
	command += " --r " + joboptions["max_error"].getString();

	bool do_shift_align = joboptions["do_shift_align"].getBoolean();
	bool do_motion = joboptions["do_motion"].getBoolean();

	int i = 0;
	if (do_shift_align) i++;
	if (do_motion) i++;

	if (i > 1)
	{
		error_message = "ERROR: Per-particle motion and shift only corrections cannot be applied simultaneously.";
		return false;
	}

	if (do_shift_align)
	{
		command += " --shift_only ";
		if (joboptions["shift_align_type"].getString() == "Only particles")
		command += " --shift_only_by_particles ";
	}

    if (do_motion)
	{
    	command += " --motion ";
    	command += " --s_vel " + joboptions["sigma_vel"].getString();
    	command += " --s_div " + joboptions["sigma_div"].getString();
    	if (joboptions["do_sq_exp_ker"].getBoolean())
    	{
    		command += " --sq_exp_ker ";
    	}
	}
	if (joboptions["do_deform"].getBoolean())
	{
		command += " --deformation ";
		command += " --def_w " + joboptions["def_w"].getString();
		command += " --def_h " + joboptions["def_h"].getString();
		command += " --def_model " + joboptions["def_model"].getString();
		command += " --def_reg " + joboptions["lambda"].getString();

		if (joboptions["do_frame_def"].getBoolean())
		{
			command += " --per_frame_deformation ";
		}
	}

	if (is_continue)
	{
		command += " --only_do_unfinished ";
	}

	// Running stuff
    command += " --j " + joboptions["nr_threads"].getString();

    // Other arguments for extraction
    command += " " + joboptions["other_args"].getString();
    commands.push_back(command);

    return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

void RelionJob::initialiseTomoReconPartJob()
{
	hidden_name = ".gui_tomo_reconstruct_particle";

	addTomoInputOptions(true, true, true, false, false, false);

   	joboptions["do_from2d"] = JobOption("Average from 2D tilt series?", true, "If set to Yes, then relion_tomo_reconstruct_particle is used, with the options below, to calculate the new average from the original 2D tilt series images. This yields the best results. If set to No, then relion_reconstruct is used to calculate the average of the 3D subtomogram images in the particle set on the I/O tab. This is quicker, but gives worse results.");
	joboptions["box_size"] = JobOption("Box size (pix):", 128, 32, 512, 16, "Box size of the reconstruction. Note that this is independent of the box size that has been used to refine the particle. This allows the user to construct a 3D map of arbitrary size to gain an overview of the structure surrounding the particle. A sufficiently large box size also allows more of the high-frequency signal to be captured that has been delocalised by the CTF.");
	joboptions["crop_size"] = JobOption("Cropped box size (pix):", -1, -1, 512, 16, "If set to a positive value, the program will output an additional set of maps that have been cropped to this size. This is useful if a map is desired that is smaller than the box size required to retrieve the CTF-delocalised signal.");
	joboptions["binning"] = JobOption("Binning factor:", 1, 1, 16, 1, "The tilt series images will be binned by this (real-valued) factor and then reconstructed in the specified box size above. Note that thereby the reconstructed region becomes larger when specifying binning factors larger than one.");
	joboptions["snr"] = JobOption("Wiener SNR constant:", 0, 0, 0.0001, 0.00001, "If set to a positive value, apply a Wiener filter with this signal-to-noise ratio. If omitted, the reconstruction will use a heuristic to prevent divisions by excessively small numbers. Please note that using a low (even though realistic) SNR might wash out the higher frequencies, which could make the map unsuitable to be used for further refinement.");
	joboptions["fn_mask"] = JobOption("FSC Solvent mask:", NODE_MASK_CPIPE, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a soft mask to automatically estimate the postprocess FSC. It will also create an optimisation set file to be used in other tomo protocols.");
	joboptions["sym_name"] = JobOption("Symmetry:", std::string("C1"), "If the molecule is asymmetric, \
set Symmetry group to C1. Note their are multiple possibilities for icosahedral symmetry: \n \
* I1: No-Crowther 222 (standard in Heymann, Chagoyen & Belnap, JSB, 151 (2005) 196–207) \n \
* I2: Crowther 222 \n \
* I3: 52-setting (as used in SPIDER?)\n \
* I4: A different 52 setting \n \
The command 'relion_refine --sym D2 --print_symmetry_ops' prints a list of all symmetry operators for symmetry group D2. \
RELION uses XMIPP's libraries for symmetry operations. \
Therefore, look at the XMIPP Wiki for more details:  http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/WebHome?topic=Symmetry");

}

bool RelionJob::getCommandsTomoReconPartJob(std::string &outputname, std::vector<std::string> &commands,
				std::string &final_command, bool do_makedir, int job_counter, std::string &error_message)
{

	commands.clear();
    initialisePipeline(outputname, job_counter);
    std::string command, command2;

    if (joboptions["do_from2d"].getBoolean())
    {
		if (joboptions["nr_mpi"].getNumber(error_message) > 1)
			command="`which relion_tomo_reconstruct_particle_mpi`";
		else
			command="`which relion_tomo_reconstruct_particle`";
		if (error_message != "") return false;

		// I/O
		error_message = getTomoInputCommmand(command, HAS_COMPULSORY, HAS_COMPULSORY, HAS_OPTIONAL, HAS_NOT, HAS_NOT,
						HAS_NOT);
		if (error_message != "") return false;

		command += " --theme classic --o " + outputname;

		Node node1(outputname+"merged.mrc", LABEL_TOMO_MAP);
		outputNodes.push_back(node1);
		Node node2(outputname+"half1.mrc", LABEL_TOMO_HALFMAP);
		outputNodes.push_back(node2);
		Node node3(outputname+"optimisation_set.star", LABEL_TOMO_OPTIMISATION);
		outputNodes.push_back(node3);

		// Job-specific stuff goes here
		command += " --b " + joboptions["box_size"].getString();

		int crop_size = joboptions["crop_size"].getNumber(error_message);
		if (error_message != "") return false;
		if (crop_size > 0.) command += " --crop " + joboptions["crop_size"].getString();

		command += " --bin " + joboptions["binning"].getString();

		float SNR = joboptions["snr"].getNumber(error_message);
		if (error_message != "") return false;
		if (SNR > 0.) command += " --SNR " + joboptions["snr"].getString();

		// Running stuff
		command += " --j " + joboptions["nr_threads"].getString();
		command += " --j_out " + joboptions["nr_threads"].getString();
		command += " --j_in 1 ";

		if (is_continue)
		{
			command += " --only_do_unfinished ";
		}

		// Estimate FSC
		if (joboptions["fn_mask"].getString() != "")
		{
			command2 = "`which relion_tomo_make_reference`";
			command2 += " --rec "+ outputname;
			command2 += " --o "+ outputname;
			command2 += " --mask "+ joboptions["fn_mask"].getString();
			error_message = getTomoInputCommmand(command2, HAS_COMPULSORY, HAS_COMPULSORY, HAS_OPTIONAL, HAS_NOT,
												 HAS_NOT,
												 HAS_NOT);
			Node node4(outputname+"PostProcess/logfile.pdf", LABEL_TOMO_POST_LOG);
			outputNodes.push_back(node4);
			Node node5(outputname+"PostProcess/postprocess.star", LABEL_TOMO_POST);
			outputNodes.push_back(node5);
		}
	}
    else
    {
		if (joboptions["in_particles"].getString() == "")
		{
			error_message = "ERROR: when not reconstructing from the 2D tilt series images, you need to provide a particle set on the I/O tab.";
			return false;
		}

		if (joboptions["nr_mpi"].getNumber(error_message) > 1)
			command="`which relion_reconstruct_mpi`";
		else
			command="`which relion_reconstruct`";
		if (error_message != "") return false;

		Node node(joboptions["in_particles"].getString(), joboptions["in_particles"].node_type);
		inputNodes.push_back(node);
    	command += " --i " + joboptions["in_particles"].getString();

		Node node1(outputname+"reconstruct.mrc", LABEL_TOMO_MAP);
		outputNodes.push_back(node1);
    	command += " --o " + outputname + "reconstruct.mrc";
    	command += " --ctf ";
    }

	command += " --sym " + joboptions["sym_name"].getString();

	// Other arguments for extraction
    command += " " + joboptions["other_args"].getString();
    commands.push_back(command);

	if (command2 != "")
	{
		commands.push_back(command2);
	}

    return prepareFinalCommand(outputname, commands, final_command, do_makedir, error_message);

}

