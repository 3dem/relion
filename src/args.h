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

#ifndef ARGS_H
#define ARGS_H

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
//ROB
#include <string.h>

#include "src/funcs.h"
#include "src/matrix1d.h"
template <typename T> class Matrix1D;


/** @defgroup Arguments Functions for parsing the command line
 *
 * These functions help you to manage the command line parameters
 */

/** Get parameters from the command line.
 * @ingroup CommandLineFunctions
 *
 * This function assumes that the command line is structured in such a way that
 * for each parameter a block like "-param <param_value>" is defined. The label
 * "param" can be substituted by any other one you like. If the parameter is
 * optional then this function allows you to define a default value. If no
 * default value is supplied and the parameter is not specified in the command
 * line, then an exception is thrown. You may change the default exception.
 *
 * You may also indicate that in case of error no exception is raised and force
 * the program to abort (use the exit variable).
 *
 * @code
 * m_param = textToFloat(getParameter(argc, argv, "-m"));
 *
 * // Get compulsory parameter "-m"
 * m_param = textToFloat(getParameter(argc, argv, "-m","2.65"));
 *
 * // Optional parameter, if no parameter is given it takes 2.65 by default
 * m_param = textToFloat(getParameter(argc, argv, "-m", NULL, 6001, "-m parameter not \
 *     found. I'm going out", TRUE);
 *
 * // Compulsory parameter, if not found give an special error message and exit
 * // the program
 *
 * @endcode
 */
std::string getParameter(int argc,
                char** argv,
                std::string param,
                std::string option = "NULL");

/** Get boolean parameters from the command line.
 * @ingroup CommandLineFunctions
 *
 * This function assumes that the command line is structured in such a way that
 * for each parameter a block like "-param" is defined. The label "param" can be
 * substituted by any other one you like. It might be used to look for a boolean
 * parameter, for instance:
 *
 *     -verbose means that verbose functionality is set (TRUE)
 *
 * @code
 * verbose = checkParameter(argc, argv, "-verbose"));
 *
 * // checks if "-verbose" was supplied in the command line. If -verbose was
 * // supplied the function returns TRUE (1), otherwise returns FALSE (0)
 * @endcode
 */
bool checkParameter(int argc, char** argv, std::string param);


class IOParser
{
private:

	std::vector<std::string> options;
	std::vector<std::string> hiddenOptions;
	std::vector<std::string> usages;
	std::vector<bool>        optionals;
	std::vector<std::string> defaultvalues;
	std::vector<int>         section_numbers;
	std::vector<std::string> section_names;
	std::vector<std::string> error_messages;
	std::vector<std::string> warning_messages;

	int current_section;

	// The original command line
	int argc;
	char** argv;

public:

    /** Constructor */
    IOParser();

    /** Copy constructor */
    IOParser(const IOParser &in);

    /**Assignment operator */
    IOParser& operator= (const IOParser &in);

    /** Destructor */
    ~IOParser();

    /** Copy everything from input to this */
    void copy(const IOParser &in);

    /** Clear object */
    void clear();

    /** Store pointer to command line */
    void setCommandLine(int _argc, char** _argv);

    /** Check whether option exists in the stored options */
    bool optionExists(std::string option);

    /** Add a section to the parser, and set the current section to the newly created one, returns number of current section */
    int addSection(std::string name);

    /** Set the current section to this number */
    void setSection(int number);

    /** Get the current section to this number */
    int getSection()
    {
    	return current_section;
    }

    /** Add an option to the object list */
    void addOption(std::string option, std::string usage, std::string defaultvalue = "NULL", bool hidden = false);

    /** Get the value from the command line, and adds option to the list if it did not yet exist */
    std::string getOption(std::string option, std::string usage, std::string defaultvalue = "NULL", bool hidden = false);

    /** Returns true if option was given and false if not, and adds option to the list if it did not yet exist */
	bool checkOption(std::string option, std::string usage, std::string defaultvalue = "false", bool hidden = false);

    /** Checks the whole command line and reports an error if it contains an undefined option */
    bool commandLineContainsUndefinedOption();

    /** Write the stored command line to outstream */
    void writeCommandLine(std::ostream &outstream);

    /** Returns true is there were any error messages (and prints them if verb>0 */
    bool checkForErrors(int verb = 1);

    /** Check the whole command line for invalid arguments, if found add to the error messages */
    void checkForUnknownArguments();

    /** Write one line of the usage to outstream */
    void writeUsageOneLine(int i, std::ostream &out);

    /** Write one section of the usage to outstream */
   void writeUsageOneSection(int section, std::ostream &out);

    /** Write the usage for all options to outstream */
    void writeUsage(std::ostream &outstream);

	/** Add an error message to be displayed after parsing */
	void reportError(const std::string& message);
};

/*
 * Takes a string with device-indices which is
 *  	: delimited for ranks
 *  	, delimited for threads within each rank
 * and outputs a rank-major array which supplies
 * a mapping as input for distribution of ranks
 * and threads over the availiable/specfied GPUs.
 */
void untangleDeviceIDs(std::string &tangled, std::vector < std::vector < std::string > > &untangled);

#endif
