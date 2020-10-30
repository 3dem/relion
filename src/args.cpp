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
/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "src/args.h"
#include "src/gcc_version.h"
#include "src/matrix1d.h"
#include <algorithm>

// Get parameters from the command line ====================================
std::string getParameter(int argc, char **argv, const std::string param, const std::string option)
{
	int i = 0;
	int i_found = -1;

	while (i < argc)
	{
//		std::cout << i << " " << i_found << " " << argv[i] << " looking for " << param << std::endl;
		if (strcmp(param.c_str(), argv[i]) == 0)
		{
			if (i_found != -1)
			{
				std::cerr << "WARNING: Command-line option " << param << 
					     " was specified more than once. The value specified the last is used." << std::endl;
			}
			i_found = i;
		}
		i++;
	}

	if (i_found > 0 && i_found < argc - 1)
	{
		return argv[i_found + 1];
	}
	else
	{
		if (option != "NULL")
		{
			return (std::string)option;
		}
		else
		{
			std::string auxstr;
			auxstr = (std::string)"Argument " + param + " not found or invalid argument";
			REPORT_ERROR(auxstr);
		}
	}
}

// Checks if a boolean parameter was included the command line =============
bool checkParameter(int argc, char **argv, std::string param)
{
	int i = 0;

	while ((i < argc) && (strcmp(param.c_str(), argv[i]) != 0))
		i++;

	if (i < argc)
		return(true);
	else
		return(false);
}

IOParser::IOParser()
{
	clear();

}

IOParser::IOParser(const IOParser &in)
{
	copy(in);
}

IOParser& IOParser::operator= (const IOParser &in)
{
	copy(in);
	return (*this);
}

IOParser::~IOParser()
{
	clear();
}


void IOParser::copy(const IOParser &in)
{
	options = in.options;
	usages = in.usages;
	optionals = in.optionals;
	defaultvalues = in.defaultvalues;
	argc = in.argc;
	argv = in.argv;
	error_messages = in.error_messages;
	warning_messages = in.warning_messages;
	current_section = in.current_section;
	section_names = in.section_names;
	section_numbers = in.section_numbers;
}

void IOParser::clear()
{
	argc = 0;
	argv = NULL;
	options.clear();
	usages.clear();
	optionals.clear();
	defaultvalues.clear();
	error_messages.clear();
	warning_messages.clear();
	section_names.clear();
	section_numbers.clear();
	current_section = 0;
}

void IOParser::setCommandLine(int _argc, char** _argv)
{
	argc = _argc;
	argv = _argv;

	// Print version of software and exit
	if ( checkParameter(argc, argv, "--version"))
	{
		PRINT_VERSION_INFO();
		exit(0);
	}
	// Dirty hack to get pipeline control for all programs...
	if (checkParameter(argc, argv, "--pipeline_control"))
		pipeline_control_outputname = getParameter(argc, argv, "--pipeline_control");
	else
		pipeline_control_outputname = "";
}

void IOParser::addOption(std::string option, std::string usage, std::string defaultvalue, bool hidden)
{
	if (hidden)
	{
		hiddenOptions.push_back(option);
	}
	else
	{
		if (section_names.size() == 0)
			REPORT_ERROR("IOParser::addOption: ERROR First add a section to the parser, then the options!");
		options.push_back(option);
		usages.push_back(usage);
		section_numbers.push_back(current_section);
		if (defaultvalue == "NULL")
		{
			optionals.push_back(false);
			defaultvalues.push_back(" ");
		}
		else
		{
			optionals.push_back(true);
			defaultvalues.push_back((std::string)defaultvalue);
		}
	}
}

int IOParser::addSection(std::string name)
{
	current_section = section_names.size();
	section_names.push_back(name);
	return current_section;
}

/** Set the current section to this number */
void IOParser::setSection(int number)
{
	current_section = number;
}

bool IOParser::optionExists(std::string option)
{
	for (int ii = 0; ii < options.size(); ii++)
		if (strcmp((options[ii]).c_str(), option.c_str()) == 0)
			return true;

	for (int ii = 0; ii < hiddenOptions.size(); ii++)
		if (strcmp((hiddenOptions[ii]).c_str(), option.c_str()) == 0)
			return true;

	return false;
}

std::string IOParser::getOption(std::string option, std::string usage, std::string defaultvalue, bool hidden)
{

	// If this option did not exist yet, add it to the list
	if (!optionExists(option))
		addOption(option, usage, defaultvalue, hidden);

	int i = 0;
	int i_found = -1;

	while (i < argc)
	{
		if (strcmp(option.c_str(), argv[i]) == 0)
		{
			if (i_found != -1)
			{
				std::cerr << "WARNING: Command-line option " << option << 
					     " was specified more than once. The value specified the last is used." << std::endl;
			}
			i_found = i;
		}
		i++;
	}

	if (i_found > 0 && i_found < argc - 1)
	{
		return argv[i_found + 1];
	}
	else
	{
		if (defaultvalue != "NULL")
		{
			return (std::string)defaultvalue;
		}
		else
		{
			std::string auxstr;
			auxstr = (std::string)"ERROR: Argument " + option + " not found or invalid argument";
			error_messages.push_back(auxstr);
			return "";
		}
	}
}

// Checks if a boolean parameter was included the command line =============
bool IOParser::checkOption(std::string option, std::string usage, std::string defaultvalue, bool hidden)
{
	// If this option did not exist yet, add it to the list
	if (!optionExists(option))
		addOption(option, usage, defaultvalue, hidden);

	return checkParameter(argc, argv, option);
}

void IOParser::writeCommandLine(std::ostream &out)
{
	for (int i = 1; i < argc; i++)
		out << argv[i] << " ";
	out << std::endl;

}

bool IOParser::checkForErrors(int verb)
{

	if(checkParameter(argc, argv, "--version"))
	{
		std::cout << "RELION version " << g_RELION_VERSION << std::endl;
		exit(0);
	}
	if(argc==1 || (argc==2 && checkParameter(argc, argv, "--continue")) || checkParameter(argc, argv, "--help") || checkParameter(argc, argv, "-h"))
	{
		writeUsage(std::cout);
	 	exit(0);
	}

	// First check the command line for unknown arguments
	checkForUnknownArguments();

	// First print warning messages
	if (warning_messages.size() > 0)
	{
		if (verb > 0)
		{
			std::cerr << "The following warnings were encountered upon command-line parsing: " << std::endl;
			for (unsigned int i = 0; i < warning_messages.size(); ++i)
				std::cerr << warning_messages[i] << std::endl;
		}
	}

	// Then check for error messages
	if (error_messages.size() > 0)
	{
		if (verb > 0)
		{
			std::cerr << "The following errors were encountered upon command-line parsing: " << std::endl;
			for (unsigned int i = 0; i < error_messages.size(); ++i)
				std::cerr << error_messages[i] << std::endl;
		}
		return true;
	}
	else
	{
		return false;
	}

}

void IOParser::checkForUnknownArguments()
{
	for (int i = 1; i < argc; i++)
	{
		// Valid options should start with "--"
		bool is_ok = true;
		if (strncmp("--", argv[i], 2) == 0)
		{
			if (!optionExists((std::string)argv[i]) && !(strncmp("--pipeline_control", argv[i], 18) == 0) )
			{
				is_ok = false;
			}
		}
		// If argv[i] starts with one "-": check it is a number and argv[i-1] is a valid option
		// or whether this is perhaps
		else if (strncmp("--", argv[i], 1) == 0)
		{
			float testval;
			// test whether this is a number
			int is_a_number = sscanf(argv[i], "%f", &testval);
			if (is_a_number)
			{
			// check whether  argv[i-1] is a valid option
			if (!optionExists(argv[i-1]))
				is_ok = false;
		 	}
			else
				is_ok = false;
		}

		if (!is_ok)
		{
			std::string auxstr;
			auxstr = (std::string)"WARNING: Option " + argv[i] + "\tis not a valid RELION argument";

			warning_messages.push_back(auxstr);
		}
	}
}

void IOParser::writeUsageOneLine(int i, std::ostream &out)
{
	std::string aux = "  ";
	aux += options[i];

	if (optionals[i])
	{
		aux += " (";
		aux += defaultvalues[i];
		aux += ")";
	}

	out << std::setw(35) << aux;
	out << " : ";
	out << usages[i];
	out << std::endl;

}

void IOParser::writeUsageOneSection(int section, std::ostream &out)
{
	// First write all compulsory options
	//out << "+++ Compulsory:" << std::endl;
	for (int i = 0; i < options.size(); i++)
	{
		if (!optionals[i] && section_numbers[i] == section)
			writeUsageOneLine(i, out);
	}

	// Then write optional ones
	//out << "+++ Optional (defaults between parentheses):" << std::endl;
	for (int i = 0; i < options.size(); i++)
	{
		if (optionals[i] && section_numbers[i] == section)
			writeUsageOneLine(i, out);
	}
}

void IOParser::writeUsage(std::ostream &out)
{
	out << "+++ RELION: command line arguments (with defaults for optional ones between parantheses) +++"<<std::endl;

	for (int section = 0; section < section_names.size(); section++)
	{
		out << "====== " << section_names[section] << " ===== " << std::endl;
		writeUsageOneSection(section, out);
	}
	out << std::setw(35) << "--version";
	out << " : Print RELION version and exit" << std::endl;
}

void IOParser::reportError(const std::string& message)
{
	error_messages.push_back(message);
}

void untangleDeviceIDs(std::string &tangled, std::vector < std::vector < std::string > > &untangled)
{
	// Handle GPU (device) assignments for each rank, if speficied
	size_t pos = 0;
	std::string delim = ":";
	std::vector < std::string > allRankIDs;
	std::string thisRankIDs, thisThreadID;
	while ((pos = tangled.find(delim)) != std::string::npos)
	{
		thisRankIDs = tangled.substr(0, pos);
//		std::cout << "in loop " << thisRankIDs << std::endl;
		tangled.erase(0, pos + delim.length());
		allRankIDs.push_back(thisRankIDs);
	}
	allRankIDs.push_back(tangled);

	untangled.resize(allRankIDs.size());
	//Now handle the thread assignements in each rank
	for (int i = 0; i < allRankIDs.size(); i++)
	{
		pos=0;
		delim = ",";
//		std::cout  << "in 2nd loop "<< allRankIDs[i] << std::endl;
		while ((pos = allRankIDs[i].find(delim)) != std::string::npos)
		{
			thisThreadID = allRankIDs[i].substr(0, pos);
//			std::cout << "in 3rd loop " << thisThreadID << std::endl;
			allRankIDs[i].erase(0, pos + delim.length());
			untangled[i].push_back(thisThreadID);
		}
		untangled[i].push_back(allRankIDs[i]);
	}
#ifdef DEBUG
	std::cout << "untangled.size() == " << untangled.size() << std::endl;
	for (int irank = 0; irank < untangled.size(); irank++)
	{
		std::cout << "untangled[" << irank << "]: ";
		for (int ithread = 0; ithread < untangled[irank].size(); ithread++)
			std::cout << untangled[irank][ithread] << " ";
		std::cout << std::endl;
	}
#endif
}

