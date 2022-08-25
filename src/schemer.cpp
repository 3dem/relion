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
#include "schemer.h"

// one global timestamp...
static time_t annotated_time;
static time_t exit_time;
bool has_annotated_time = false;
bool has_exit_time = false;


// Global variables, but only with reach within this file!
std::map<std::string, SchemerBooleanVariable> schemer_global_bools;
std::map<std::string, SchemerFloatVariable> schemer_global_floats;
std::map<std::string, SchemerStringVariable> schemer_global_strings;
std::map<std::string, SchemerOperator> schemer_global_operators;

bool isBooleanVariable(std::string _name)
{
	return (schemer_global_bools.find(_name) != schemer_global_bools.end());
}

bool isFloatVariable(std::string _name)
{
	return (schemer_global_floats.find(_name) != schemer_global_floats.end());
}

bool isStringVariable(std::string _name)
{
	return (schemer_global_strings.find(_name) != schemer_global_strings.end());
}

bool isSchemeOperator(std::string _name)
{
	return (schemer_global_operators.find(_name) != schemer_global_operators.end());
}

SchemerOperator::SchemerOperator(std::string _type, std::string _input1, std::string _input2, std::string _output)
{
	std::string myerror = initialise(_type, _input1, _input2, _output);
	if (myerror != "") REPORT_ERROR(myerror);
}

std::string SchemerOperator::initialise(std::string _type, std::string _input1, std::string _input2, std::string _output)
{
	type = _type;

	// Check output
	if ((type == SCHEME_BOOLEAN_OPERATOR_GT ||
		 type == SCHEME_BOOLEAN_OPERATOR_LT ||
		 type == SCHEME_BOOLEAN_OPERATOR_EQ ||
		 type == SCHEME_BOOLEAN_OPERATOR_GE ||
		 type == SCHEME_BOOLEAN_OPERATOR_LE ||
		 type == SCHEME_BOOLEAN_OPERATOR_SET ||
		 type == SCHEME_BOOLEAN_OPERATOR_AND ||
		 type == SCHEME_BOOLEAN_OPERATOR_OR ||
		 type == SCHEME_BOOLEAN_OPERATOR_FILE_EXISTS ||
		 type == SCHEME_BOOLEAN_OPERATOR_READ_STAR
		) && !isBooleanVariable(_output))
		return "ERROR: boolean operator does not have valid boolean output: " + _output;
	if ((type == SCHEME_FLOAT_OPERATOR_SET ||
		 type == SCHEME_FLOAT_OPERATOR_PLUS ||
		 type == SCHEME_FLOAT_OPERATOR_MINUS ||
		 type == SCHEME_FLOAT_OPERATOR_MULT ||
		 type == SCHEME_FLOAT_OPERATOR_DIVIDE ||
		 type == SCHEME_FLOAT_OPERATOR_ROUND ||
		 type == SCHEME_FLOAT_OPERATOR_COUNT_IMAGES ||
		 type == SCHEME_FLOAT_OPERATOR_COUNT_WORDS ||
		 type == SCHEME_FLOAT_OPERATOR_READ_STAR ||
		 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MAX ||
		 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MIN ||
		 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_AVG ||
		 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX
		 ) && !isFloatVariable(_output))
		return "ERROR: float operator does not have valid float output: " + _output;
	if ((type == SCHEME_STRING_OPERATOR_READ_STAR ||
		 type == SCHEME_STRING_OPERATOR_SET ||
		 type == SCHEME_STRING_OPERATOR_JOIN ||
		 type == SCHEME_STRING_OPERATOR_BEFORE_FIRST ||
		 type == SCHEME_STRING_OPERATOR_AFTER_FIRST ||
		 type == SCHEME_STRING_OPERATOR_BEFORE_LAST ||
		 type == SCHEME_STRING_OPERATOR_AFTER_LAST ||
		 type == SCHEME_STRING_OPERATOR_GLOB ||
		 type == SCHEME_STRING_OPERATOR_NTH_WORD
		 )	&& ! isStringVariable(_output))
		return "ERROR: string operator does not have valid string output: " + _output;

	// Check input1
	if ((type == SCHEME_BOOLEAN_OPERATOR_SET ||
		 type == SCHEME_BOOLEAN_OPERATOR_AND ||
		 type == SCHEME_BOOLEAN_OPERATOR_OR )  && !isBooleanVariable(_input1))
		return "ERROR: operator does not have valid boolean input1: " + _input1;
	if ((type == SCHEME_FLOAT_OPERATOR_SET ||
		 type == SCHEME_FLOAT_OPERATOR_PLUS ||
		 type == SCHEME_FLOAT_OPERATOR_MINUS ||
		 type == SCHEME_FLOAT_OPERATOR_MULT ||
		 type == SCHEME_FLOAT_OPERATOR_DIVIDE ||
		 type == SCHEME_FLOAT_OPERATOR_ROUND ||
		 type == SCHEME_EXIT_MAXTIME
		 ) && !isFloatVariable(_input1))
		return "ERROR: operator does not have valid float input1: " + _input1;
	if ((type == SCHEME_BOOLEAN_OPERATOR_READ_STAR ||
		 type == SCHEME_BOOLEAN_OPERATOR_FILE_EXISTS ||
		 type == SCHEME_FLOAT_OPERATOR_COUNT_IMAGES ||
		 type == SCHEME_FLOAT_OPERATOR_COUNT_WORDS ||
		 type == SCHEME_FLOAT_OPERATOR_READ_STAR ||
		 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MAX ||
		 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MIN ||
		 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_AVG ||
		 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX ||
		 type == SCHEME_STRING_OPERATOR_READ_STAR ||
		 type == SCHEME_STRING_OPERATOR_SET ||
		 type == SCHEME_STRING_OPERATOR_JOIN ||
		 type == SCHEME_STRING_OPERATOR_BEFORE_FIRST ||
		 type == SCHEME_STRING_OPERATOR_AFTER_FIRST ||
		 type == SCHEME_STRING_OPERATOR_BEFORE_LAST ||
		 type == SCHEME_STRING_OPERATOR_AFTER_LAST ||
		 type == SCHEME_STRING_OPERATOR_GLOB ||
		 type == SCHEME_STRING_OPERATOR_NTH_WORD ||
		 type == SCHEME_OPERATOR_TOUCH_FILE ||
		 type == SCHEME_OPERATOR_COPY_FILE ||
		 type == SCHEME_OPERATOR_MOVE_FILE ||
		 type == SCHEME_OPERATOR_DELETE_FILE
		 ) && ! isStringVariable(_input1))
		return "ERROR: operator does not have valid string input1: " + _input1;

	// Check input2
	if ((type == SCHEME_BOOLEAN_OPERATOR_AND ||
		 type == SCHEME_BOOLEAN_OPERATOR_OR
		 ) && !isBooleanVariable(_input2))
		return "ERROR: operator does not have valid boolean input2: " + _input2;
	if ((type == SCHEME_BOOLEAN_OPERATOR_GT ||
		 type == SCHEME_BOOLEAN_OPERATOR_LT ||
		 type == SCHEME_BOOLEAN_OPERATOR_EQ ||
		 type == SCHEME_BOOLEAN_OPERATOR_GE ||
		 type == SCHEME_BOOLEAN_OPERATOR_LE ||
		 type == SCHEME_FLOAT_OPERATOR_PLUS ||
		 type == SCHEME_FLOAT_OPERATOR_MINUS ||
		 type == SCHEME_FLOAT_OPERATOR_MULT ||
		 type == SCHEME_FLOAT_OPERATOR_DIVIDE ||
		 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX ||
		 type == SCHEME_STRING_OPERATOR_NTH_WORD
		 ) && !(isFloatVariable(_input2) || isNumber(_input2)))
		return "ERROR: operator does not have valid number (float variable or text) input2: " + _input2;
	if ((type == SCHEME_OPERATOR_COPY_FILE ||
		 type == SCHEME_OPERATOR_MOVE_FILE ||
		 type == SCHEME_STRING_OPERATOR_BEFORE_FIRST ||
		 type == SCHEME_STRING_OPERATOR_AFTER_FIRST ||
		 type == SCHEME_STRING_OPERATOR_BEFORE_LAST ||
		 type == SCHEME_STRING_OPERATOR_AFTER_LAST ||
		 type == SCHEME_STRING_OPERATOR_JOIN
	 ) && ! isStringVariable(_input2))
	return "ERROR: operator does not have valid string input2: " + _input2;

	input1 = (_input1 == "") ? "undefined" : _input1;
	input2 = (_input2 == "") ? "undefined" : _input2;
	output = (_output == "") ? "undefined" : _output;

	return "";
}

// Separate comma-separated labels for table, input and output
void  SchemerOperator::readFromStarFile() const
{
	MetaDataTable MD;
	std::string mystring, mystarfile, mytable;
	EMDLabel mylabel;

	// The localtion is always in input1
	mystring = schemer_global_strings[input1].value;
	std::vector< std::string > splits;
	int nr_splits = splitString(mystring, ",", splits);
	if (splits.size() < 3) REPORT_ERROR("Need at least three comma-separated values for starfilename, tablename and labelname");
	mystarfile = splits[0];
	mytable = splits[1];
	mylabel = EMDL::str2Label(splits[2]);

	// Read the correct table from the STAR file
	MD.read(mystarfile, mytable);

	int ival;
	long idxmin, idxmax;
	RFLOAT mymin=99.e99;
	RFLOAT mymax = -mymin;
	RFLOAT mysum = 0.;
	RFLOAT myval;
	long ii = 0;
	MultidimArray<RFLOAT> for_sorting(MD.numberOfObjects());

	long idx = (isFloatVariable(input2)) ? ROUND(schemer_global_floats[input2].value) : 0;
	if (type == SCHEME_BOOLEAN_OPERATOR_READ_STAR ||
		type == SCHEME_FLOAT_OPERATOR_READ_STAR ||
		type == SCHEME_STRING_OPERATOR_READ_STAR)
	{

		if (EMDL::isDouble(mylabel))
		{
			RFLOAT fval;
			MD.getValue(mylabel, fval, idx);
			schemer_global_floats[output].value = fval;
		}
		else if (EMDL::isInt(mylabel))
		{
			int ival;
			MD.getValue(mylabel, ival, idx);
			schemer_global_floats[output].value = ROUND(ival);
		}
		else if (EMDL::isString(mylabel))
		{
			std::string val;
			MD.getValue(mylabel, val, idx);
			schemer_global_strings[output].value = val;
		}
		else if (EMDL::isBool(mylabel))
		{
			bool val;
			MD.getValue(mylabel, val, idx);
			schemer_global_bools[output].value = val;
		}
	}
	else if (type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MAX ||
			type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MIN ||
			type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_AVG ||
			type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX)
	{
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			if (EMDL::isDouble(mylabel))
			{
				RFLOAT fval;
				MD.getValue(mylabel, fval);
				myval = fval;
			}
			else if (EMDL::isInt(mylabel))
			{
				int ival;
				MD.getValue(mylabel, ival);
				myval = ival;
			}
			else
				REPORT_ERROR("ERROR: metadata label " + EMDL::label2Str(mylabel) + " is not of a number type!");

			DIRECT_MULTIDIM_ELEM(for_sorting, ii) = myval;

			if (myval < mymin)
			{
				mymin = myval;
				idxmin = ii;
			}
			if (myval > mymax)
			{
				mymax = myval;
				idxmax = ii;
			}
			mymax = XMIPP_MAX(myval, mymax);
			mysum += myval;
			ii++;
		}

		if (ii > 0)
			mysum /= (RFLOAT)ii;

		if (type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MAX) schemer_global_floats[output].value = mymax;
		else if (type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MIN) schemer_global_floats[output].value = mymin;
		else if (type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_AVG) schemer_global_floats[output].value = mysum;
		else if (type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX)
		{
			if (idx == 0) REPORT_ERROR("Give a positive or a negative value for input2 in sorted_idx: 1 is largest, -1 is smallest value");
			if (idx < 0) idx = ii + idx; // smallest value is numberOfObjects - 1
			else idx--; // now start counting at 0
			MultidimArray<long> sorted_idx;
			for_sorting.sorted_index(sorted_idx);
			schemer_global_floats[output].value = DIRECT_MULTIDIM_ELEM(sorted_idx, idx);
		}
	}
}

bool SchemerOperator::performOperation() const
{
	RFLOAT val2 = (isFloatVariable(input2)) ? schemer_global_floats[input2].value : 0;

	if (type == SCHEME_BOOLEAN_OPERATOR_SET)
	{
		schemer_global_bools[output].value = schemer_global_bools[input1].value;
	}
	else if (type == SCHEME_BOOLEAN_OPERATOR_AND)
	{
		schemer_global_bools[output].value = (schemer_global_bools[input1].value && schemer_global_bools[input2].value);
	}
	else if (type == SCHEME_BOOLEAN_OPERATOR_OR)
	{
		schemer_global_bools[output].value = (schemer_global_bools[input1].value || schemer_global_bools[input2].value);
	}
	else if (type == SCHEME_BOOLEAN_OPERATOR_NOT)
	{
		schemer_global_bools[output].value = (!(schemer_global_bools[input1].value));
	}
	else if (type == SCHEME_BOOLEAN_OPERATOR_GT)
	{
		schemer_global_bools[output].value = (schemer_global_floats[input1].value > val2);
	}
	else if (type == SCHEME_BOOLEAN_OPERATOR_LT)
	{
		schemer_global_bools[output].value = (schemer_global_floats[input1].value < val2);
	}
	else if (type == SCHEME_BOOLEAN_OPERATOR_EQ)
	{
		schemer_global_bools[output].value = (fabs(schemer_global_floats[input1].value - val2) < 1E-8);
	}
	else if (type == SCHEME_BOOLEAN_OPERATOR_GE)
	{
		schemer_global_bools[output].value = (schemer_global_floats[input1].value >= val2);
	}
	else if (type == SCHEME_BOOLEAN_OPERATOR_LE)
	{
		schemer_global_bools[output].value = (schemer_global_floats[input1].value <= val2);
	}
	else if (type == SCHEME_BOOLEAN_OPERATOR_FILE_EXISTS)
	{
		schemer_global_bools[output].value = (exists(schemer_global_strings[input1].value));
	}
	else if (type == SCHEME_BOOLEAN_OPERATOR_READ_STAR ||
			 type == SCHEME_FLOAT_OPERATOR_READ_STAR ||
			 type == SCHEME_STRING_OPERATOR_READ_STAR ||
			 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MAX ||
			 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MIN ||
			 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_AVG ||
			 type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX
			 )
	{
		readFromStarFile();
	}
	else if (type == SCHEME_FLOAT_OPERATOR_SET)
	{
		schemer_global_floats[output].value = schemer_global_floats[input1].value;
	}
	else if (type == SCHEME_FLOAT_OPERATOR_PLUS)
	{
		schemer_global_floats[output].value = schemer_global_floats[input1].value + val2;
	}
	else if (type ==SCHEME_FLOAT_OPERATOR_MINUS)
	{
		schemer_global_floats[output].value = schemer_global_floats[input1].value - val2;
	}
	else if (type == SCHEME_FLOAT_OPERATOR_MULT)
	{
		schemer_global_floats[output].value = schemer_global_floats[input1].value * val2;
	}
	else if (type == SCHEME_FLOAT_OPERATOR_DIVIDE)
	{
		schemer_global_floats[output].value = schemer_global_floats[input1].value / val2;
	}
	else if (type == SCHEME_FLOAT_OPERATOR_ROUND)
	{
		schemer_global_floats[output].value = ROUND(schemer_global_floats[input1].value);
	}
	else if (type == SCHEME_FLOAT_OPERATOR_COUNT_IMAGES)
	{
		if (!exists(schemer_global_strings[input1].value))
		{
			schemer_global_floats[output].value = 0;
		}
		else
		{
			ObservationModel obsmodel;
			MetaDataTable MDimg;
			std::string mytablename = (isStringVariable(input2)) ? schemer_global_strings[input2].value : "particles";
			ObservationModel::loadSafely(schemer_global_strings[input1].value, obsmodel, MDimg, mytablename);
			schemer_global_floats[output].value = MDimg.numberOfObjects();
		}
	}
	else if (type == SCHEME_FLOAT_OPERATOR_COUNT_WORDS)
	{
		// return zero for an undefined string
		if (schemer_global_strings[input1].value == "undefined")
		{
			schemer_global_floats[output].value = 0;
		}
		else
		{
			std::vector< std::string > splits;
			int nr_splits = splitString(schemer_global_strings[input1].value, ",", splits);
			schemer_global_floats[output].value = splits.size();
		}
	}
	else if (type == SCHEME_STRING_OPERATOR_SET)
	{
		schemer_global_strings[output].value = schemer_global_strings[input1].value;
	}
	else if (type == SCHEME_STRING_OPERATOR_JOIN)
	{
		schemer_global_strings[output].value = schemer_global_strings[input1].value + schemer_global_strings[input2].value;
	}
	else if (type == SCHEME_STRING_OPERATOR_BEFORE_FIRST)
	{
		schemer_global_strings[output].value = schemer_global_strings[input1].value.beforeFirstOf(schemer_global_strings[input2].value);
	}
	else if (type == SCHEME_STRING_OPERATOR_AFTER_FIRST)
	{
		schemer_global_strings[output].value = schemer_global_strings[input1].value.afterFirstOf(schemer_global_strings[input2].value);
	}
	else if (type == SCHEME_STRING_OPERATOR_BEFORE_LAST)
	{
		schemer_global_strings[output].value = schemer_global_strings[input1].value.beforeLastOf(schemer_global_strings[input2].value);
	}
	else if (type == SCHEME_STRING_OPERATOR_AFTER_LAST)
	{
		schemer_global_strings[output].value = schemer_global_strings[input1].value.afterLastOf(schemer_global_strings[input2].value);
	}
	else if (type == SCHEME_STRING_OPERATOR_GLOB)
	{
		FileName input=schemer_global_strings[input1].value;
		std::vector<FileName> files;
		input.globFiles(files);
		if (files.size() == 0)
		{
			schemer_global_strings[output].value = "undefined";
		}
		else
		{
			schemer_global_strings[output].value = files[0];
			for (int i = 1; i < files.size(); i++)
			{
				schemer_global_strings[output].value += "," + files[i];
			}
		}
	}
	else if (type == SCHEME_STRING_OPERATOR_NTH_WORD)
	{

		std::vector< std::string > splits;
		int nr_splits = splitString(schemer_global_strings[input1].value, ",", splits);
		int mypos = ROUND(val2);
		// for negative Ns, count from the back
		if (mypos < 0) mypos = splits.size() - mypos + 1;
		// Started counting at 1, but first element of vector is zero!
		mypos--;
		if (mypos >= splits.size() || mypos < 0)
		{
			schemer_global_strings[output].value = "undefined";
		}
		else
		{
			schemer_global_strings[output].value = splits[mypos];
		}
	}
	else if (type == SCHEME_OPERATOR_TOUCH_FILE)
	{
		std::cout << " + Touching: " << schemer_global_strings[input1].value << std::endl;
		touch(schemer_global_strings[input1].value);
	}
	else if (type == SCHEME_OPERATOR_COPY_FILE || type == SCHEME_OPERATOR_MOVE_FILE)
	{
		std::string mycommand;
		if (type == SCHEME_OPERATOR_COPY_FILE)
		{
			std::cout << " + Copying: " << schemer_global_strings[input1].value << " to " << schemer_global_strings[input2].value << std::endl;
			mycommand = "cp ";
		}
		else
		{
			std::cout << " + Moving: " << schemer_global_strings[input1].value << " to " << schemer_global_strings[input2].value << std::endl;
			mycommand = "mv ";
		}
		// Make output directory if it doesn't exist
		if (schemer_global_strings[input2].value.contains("/"))
		{
			FileName mydirs = schemer_global_strings[input2].value.beforeLastOf("/");
			mktree(mydirs);
		}
		// Execute the command
		mycommand += schemer_global_strings[input1].value + " " + schemer_global_strings[input2].value;
		int res = system(mycommand.c_str());
	}
	else if (type == SCHEME_OPERATOR_DELETE_FILE)
	{
		std::string mycommand = "rm -f " + schemer_global_strings[input1].value;
		std::cout << " + Deleting: " << schemer_global_strings[input1].value << std::endl;
		int res = system(mycommand.c_str());
	}
	else if (type == SCHEME_WAIT_OPERATOR_SINCE_LAST_TIME)
	{
		if (has_annotated_time)
		{
			time_t current_time = time(NULL);
			RFLOAT elapsed = current_time - annotated_time;
			// Also set output to elapsed seconds if it is a floatVariable
			if (isFloatVariable(output)) schemer_global_floats[output].value = elapsed;
			RFLOAT wait_seconds =  schemer_global_floats[input1].value - elapsed;
			if (wait_seconds > 0)
			{
				std::cout << " + Waiting for " << wait_seconds << " seconds ..." << std::endl;
				RFLOAT waited = 0.;
				while (waited < wait_seconds)
				{
					sleep(10);
					waited += 10.;
					// Abort mechanism
					if (pipeline_control_check_abort_job())
					{
						std::cout << " + Interrupted waiting due to a abort signal." << std::endl;
						break;
					}
				}

				std::cout << " + Finished waiting." << std::endl;
			}
			else
			{
				std::cout << " + Not waiting, as more than " << schemer_global_floats[input1].value << " seconds have passed." << std::endl;
			}
		}
		annotated_time = time(NULL);
		has_annotated_time =true;
	}
	else if (type == SCHEME_EMAIL_OPERATOR)
	{

		time_t my_time = time(NULL);
		std::string mymessage = std::string(ctime(&my_time)) + "\n";
		mymessage += "input1: " + input1 + " = ";
		if (isStringVariable(input1)) mymessage += schemer_global_strings[input1].value + "\n";
		else if (isBooleanVariable(input1)) mymessage += (schemer_global_bools[input1].value) ? "True \n" : "False \n";
		else if (isFloatVariable(input1)) mymessage += floatToString(schemer_global_floats[input1].value) + "\n";
		if (isBooleanVariable(input2) || isFloatVariable(input2) || isStringVariable(input2))
		{
			mymessage += "input2: " + input2 + " = ";
			if (isStringVariable(input2)) mymessage += schemer_global_strings[input2].value + "\n";
			else if (isBooleanVariable(input2)) mymessage += (schemer_global_bools[input2].value) ? "True \n" : "False \n";
			else if (isFloatVariable(input2)) mymessage += floatToString(schemer_global_floats[input2].value) + "\n";
		}
		std::cout << " + Sending e-mail." << std::endl;

		schemerSendEmail(mymessage);
	}
	else if (type == SCHEME_EXIT_MAXTIME)
	{
		// The first time round, set the exit time, else check whether it has passed
		if (!has_exit_time)
		{
			time_t my_time = time(NULL);
			// input1 maximum time is in hours, add to my_time in seconds
			exit_time = my_time + schemer_global_floats[input1].value*3600;
			has_exit_time = true;
			tm *gmt_time = gmtime(&exit_time);
			std::cout << " + Setting exit time at: " << asctime(gmt_time);
		}
		else
		{
			time_t my_time = time(NULL);
			if (my_time >= exit_time)
			{
				tm *gmt_time = gmtime(&my_time);
				std::cout << " + It is now: " << asctime(gmt_time);
				std::cout << " + The scheme has reached its exit time. Exiting ..." << std::endl;
				return false; // to exit the scheme
			}
		}

	}
	else if (type == SCHEME_EXIT_OPERATOR)
	{
		std::cout << " + The scheme has reached an exit point. Exiting ..." << std::endl;
		return false; // to exit the scheme
	}
	else
		REPORT_ERROR("ERROR: unrecognised Operator type:" + type);

	return true;
}

std::string SchemerOperator::getName()
{
	if (type == SCHEME_BOOLEAN_OPERATOR_GT) return output + "=" + input1 + "_GT_" + input2;
	if (type == SCHEME_BOOLEAN_OPERATOR_LT) return output + "=" + input1 + "_LT_" + input2;
	if (type == SCHEME_BOOLEAN_OPERATOR_EQ) return output + "=" + input1 + "_EQ_" + input2;
	if (type == SCHEME_BOOLEAN_OPERATOR_GE) return output + "=" + input1 + "_GE_" + input2;
	if (type == SCHEME_BOOLEAN_OPERATOR_LE) return output + "=" + input1 + "_LE_" + input2;
	if (type == SCHEME_BOOLEAN_OPERATOR_AND) return output + "=" + input1 + "_AND_" + input2;
	if (type == SCHEME_BOOLEAN_OPERATOR_OR) return output + "=" + input1 + "_OR_" + input2;
	if (type == SCHEME_BOOLEAN_OPERATOR_NOT) return output + "=" + "NOT_" + input1;
	if (type == SCHEME_BOOLEAN_OPERATOR_FILE_EXISTS) return output + "=" + "EXISTS_" + input1;
	if (type == SCHEME_BOOLEAN_OPERATOR_READ_STAR) return output + "=" + "STAR_" + input1 + "_" + input2;
	if (type == SCHEME_FLOAT_OPERATOR_SET) return output + "=" + "SET_" + input1;
	if (type == SCHEME_FLOAT_OPERATOR_PLUS) return output + "=" + input1 + "_PLUS_" + input2;
	if (type == SCHEME_FLOAT_OPERATOR_MINUS) return output + "=" + input1 + "_MINUS_" + input2;
	if (type == SCHEME_FLOAT_OPERATOR_MULT) return output + "=" + input1 + "_MULT_" + input2;
	if (type == SCHEME_FLOAT_OPERATOR_DIVIDE) return output + "=" + input1 + "_DIV_" + input2;
	if (type == SCHEME_FLOAT_OPERATOR_ROUND) return output + "=" + "ROUND_" + input1;
	if (type == SCHEME_FLOAT_OPERATOR_COUNT_IMAGES) return output + "=" + "COUNT_IMGS_" + input1 + "_" + input2;
	if (type == SCHEME_FLOAT_OPERATOR_COUNT_WORDS) return output + "=" + "COUNT_WORDS_" + input1;
	if (type == SCHEME_FLOAT_OPERATOR_READ_STAR) return output + "=" + "STAR_" + input1 + "_" + input2;
	if (type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MAX) return output + "=" + "STAR_MAX_" + input1;
	if (type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MIN) return output + "=" + "STAR_MIN_" + input1;
	if (type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_AVG) return output + "=" + "STAR_AVG_" + input1;
	if (type == SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX) return output + "=" + "STAR_SORT_IDX_" + input1 + "_" + input2;
	if (type == SCHEME_STRING_OPERATOR_JOIN) return output + "=" + "JOIN_" + input1 + "_" + input2;
	if (type == SCHEME_STRING_OPERATOR_BEFORE_FIRST) return output + "=" + "BEF_FIRST_" + input1 + "_" + input2;
	if (type == SCHEME_STRING_OPERATOR_AFTER_FIRST) return output + "=" + "AFT_FIRST_" + input1 + "_" + input2;
	if (type == SCHEME_STRING_OPERATOR_BEFORE_LAST) return output + "=" + "BEF_LAST_" + input1 + "_" + input2;
	if (type == SCHEME_STRING_OPERATOR_AFTER_LAST) return output + "=" + "AFT_LAST_" + input1 + "_" + input2;
	if (type == SCHEME_STRING_OPERATOR_READ_STAR) return output + "=" + "STAR_" + input1 + "_" + input2;
	if (type == SCHEME_STRING_OPERATOR_GLOB) return output + "=" + "GLOB_" + input1;
	if (type == SCHEME_STRING_OPERATOR_NTH_WORD) return output + "=" + "NTH_WORD_" + input1 + "_" + input2;
	if (type == SCHEME_OPERATOR_TOUCH_FILE) return "TOUCH_" + input1;
	if (type == SCHEME_OPERATOR_COPY_FILE) return "COPY_" + input1 + "_TO_" + input2;
	if (type == SCHEME_OPERATOR_MOVE_FILE) return "MOVE_" + input1 + "_TO_" + input2;
	if (type == SCHEME_OPERATOR_DELETE_FILE) return "DELETE_" + input1;
	if (type == SCHEME_WAIT_OPERATOR_SINCE_LAST_TIME) return "WAIT_" + input1;
	if (type == SCHEME_EMAIL_OPERATOR) return "EMAIL_" + input1 + "_" + input2;
	if (type == SCHEME_EXIT_OPERATOR) return "EXIT";
	else
		REPORT_ERROR("ERROR: unrecognised Operator type:" + type);
}

std::string SchemerEdge::getOutputNode() const
{
	if (is_fork) return (schemer_global_bools[myBooleanVariable].value) ?  outputNodeTrue : outputNode;
	else return outputNode;
}

void Scheme::clear()
{
	verb = 1;
	current_node = "undefined";
	name = "undefined";
	do_read_only = false;
	schemer_global_bools.clear();
	schemer_global_floats.clear();
	schemer_global_strings.clear();
	jobs.clear();
	edges.clear();
	schemer_global_operators.clear();
}

std::string Scheme::findJobByCurrentName(std::string _name)
{
	std::map<std::string, SchemerJob>::iterator it;

	for ( it = jobs.begin(); it != jobs.end(); it++ )
	{
	    if (it->second.current_name == _name)
	    	return it->second.current_name;
	}
	REPORT_ERROR("ERROR: cannot find job: " + _name);
	return "";
}

void Scheme::read(bool do_lock, FileName fn)
{
#ifdef DEBUG_LOCK
	std::cerr << "entering read lock_message=" << lock_message << std::endl;
#endif

	if (fn == "") fn = name + "scheme.star";

	FileName name_wo_dir = fn.beforeLastOf("/");
	FileName dir_lock=".relion_lock_scheme_" + name_wo_dir.afterLastOf("/"), fn_lock=dir_lock + "/lock_scheme";;
	if (do_lock && !do_read_only)
	{
		int iwait =0;
		int status = mkdir(dir_lock.c_str(), S_IRWXU);

#ifdef DEBUG_LOCK
		std::cerr <<  " A status= " << status << std::endl;
#endif
		while (status != 0)
		{
			if (errno == EACCES) // interestingly, not EACCESS!
				REPORT_ERROR("ERROR: Scheme::read cannot create a lock directory " + dir_lock + ". You don't have write permission to this project. If you want to look at other's project directory (but run nothing there), please start RELION with --readonly.");

			// If the lock exists: wait 3 seconds and try again
			// First time round, print a warning message
			if (iwait == 0)
			{
				std::cout << " WARNING: trying to read scheme.star, but directory " << dir_lock << " exists (which protects against simultaneous writing)" << std::endl;
			}
			sleep(3);
			status =  mkdir(dir_lock.c_str(), S_IRWXU);
#ifdef DEBUG_LOCK
			std::cerr <<  " B status= " << status << std::endl;
#endif

			iwait++;
			if (iwait > 40)
			{

				REPORT_ERROR("ERROR: Scheme::read has waited for 2 minutes for lock directory to disappear. Make sure this schemer is not running, and then manually remove the file: " + fn_lock);
			}

		}
		// Generate the lock file
		std::ofstream  fh;
		fh.open(fn_lock.c_str(), std::ios::out);
		if (!fh)
			REPORT_ERROR( (std::string)"ERROR: Cannot open file: " + fn_lock);
		std::string lock_message = "lock mechanism from Schemer";
		fh << lock_message << std::endl;
		fh.close();
	}

	// Clear current model
	clear();

	// Open input file
	std::ifstream in(fn.data(), std::ios_base::in);
	if (in.fail())
		REPORT_ERROR( (std::string) "Scheme::read: File " + fn + " cannot be read." );

	// For reading: do the nodes before the general table, in order to set current_node
	MetaDataTable MD;
	MD.readStar(in, "scheme_general");
	MD.getValue(EMDL_SCHEME_GENERAL_NAME, name);
	MD.getValue(EMDL_SCHEME_GENERAL_CURRENT_NODE, current_node);
	MD.clear();

	MD.readStar(in, "scheme_floats");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string myname;
		RFLOAT value, original_value;
		MD.getValue(EMDL_SCHEME_VAR_FLOAT_NAME, myname);
		MD.getValue(EMDL_SCHEME_VAR_FLOAT_VALUE, value);
		MD.getValue(EMDL_SCHEME_VAR_FLOAT_ORI_VALUE, original_value);
		SchemerFloatVariable myval(value, original_value);
		schemer_global_floats[myname] = myval;
	}
	MD.clear();

	MD.readStar(in, "scheme_bools");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string myname;
		bool value, original_value;
		MD.getValue(EMDL_SCHEME_VAR_BOOL_NAME, myname);
		MD.getValue(EMDL_SCHEME_VAR_BOOL_VALUE, value);
		MD.getValue(EMDL_SCHEME_VAR_BOOL_ORI_VALUE, original_value);
		SchemerBooleanVariable myval(value, original_value);
		schemer_global_bools[myname] = myval;
	}
	MD.clear();

	MD.readStar(in, "scheme_strings");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string myname;
		FileName value, original_value;
		MD.getValue(EMDL_SCHEME_VAR_STRING_NAME, myname);
		MD.getValue(EMDL_SCHEME_VAR_STRING_VALUE, value);
		MD.getValue(EMDL_SCHEME_VAR_STRING_ORI_VALUE, original_value);
		SchemerStringVariable myval(value, original_value);
		schemer_global_strings[myname] = myval;
	}
	MD.clear();

	MD.readStar(in, "scheme_operators");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string myname, type, input1, input2, output;
		RFLOAT constant;

		MD.getValue(EMDL_SCHEME_OPERATOR_NAME, myname);
		MD.getValue(EMDL_SCHEME_OPERATOR_TYPE, type);
		MD.getValue(EMDL_SCHEME_OPERATOR_INPUT1, input1);
		MD.getValue(EMDL_SCHEME_OPERATOR_INPUT2, input2);
		MD.getValue(EMDL_SCHEME_OPERATOR_OUTPUT, output);
		SchemerOperator myval(type, input1, input2, output);
		schemer_global_operators[myname] = myval;
	}
	MD.clear();

	MD.readStar(in, "scheme_jobs");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string myname, ori_name, mode, type;
		bool has_started;

		MD.getValue(EMDL_SCHEME_JOB_NAME, myname);
		MD.getValue(EMDL_SCHEME_JOB_ORI_NAME, ori_name);
		MD.getValue(EMDL_SCHEME_JOB_MODE, mode);
		MD.getValue(EMDL_SCHEME_JOB_HAS_STARTED, has_started);

		SchemerJob myval(myname, mode, has_started);
		jobs[ori_name] = myval;
	}
	MD.clear();


	MD.readStar(in, "scheme_edges");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		int number;
		std::string inputname, outputname, outputname_true, bool_name;
		bool is_fork;

		MD.getValue(EMDL_SCHEME_EDGE_INPUT, inputname);
		MD.getValue(EMDL_SCHEME_EDGE_OUTPUT, outputname);
		MD.getValue(EMDL_SCHEME_EDGE_IS_FORK, is_fork);
		MD.getValue(EMDL_SCHEME_EDGE_OUTPUT_TRUE, outputname_true);
		MD.getValue(EMDL_SCHEME_EDGE_BOOLEAN, bool_name);
		SchemerEdge myval(inputname, outputname, is_fork, bool_name, outputname_true);
		edges.push_back(myval);
	}
	MD.clear();

	// Never let current_node be undefined...
	if (current_node == "undefined" && edges.size() > 0)
	{
		current_node = edges[0].inputNode;
	}

	// Close file handler
	in.close();

}

bool Scheme::isWriteLocked()
{
	FileName name_wo_dir = name;
	name_wo_dir = name_wo_dir.beforeLastOf("/");
	FileName dir_lock=".relion_lock_scheme_" + name_wo_dir.afterLastOf("/"), fn_lock=dir_lock + "/lock_scheme";;
	return exists(dir_lock);
}

void Scheme::write(bool do_lock, FileName fn)
{
	if (do_read_only)
		return;

	FileName name_wo_dir = name;
	name_wo_dir = name_wo_dir.beforeLastOf("/");
	FileName dir_lock=".relion_lock_scheme_" + name_wo_dir.afterLastOf("/"), fn_lock=dir_lock + "/lock_scheme";;
	if (do_lock)
	{

#ifdef DEBUG_LOCK
		if (exists(fn_lock))
		{
			std::cerr << "writing pipeline: " << fn_lock << " exists as expected" << std::endl;
		}
#endif

		int iwait =0;
		while( !exists(fn_lock) )
		{
			// If the lock exists: wait 3 seconds and try again
			// First time round, print a warning message
			if (iwait == 0)
			{
				std::cerr << " WARNING: was expecting a file called "+fn_lock+ " but it isn't there. Will wait for 1 minute to see whether it appears" << std::endl;
			}
			sleep(3);
			iwait++;
			if (iwait > 40)
			{
				REPORT_ERROR("ERROR: PipeLine::write has waited for 2 minutes for lock file to appear, but it doesn't....");
			}
		}
	}

	if (fn == "") fn = name + "scheme.star";

	// B. Write STAR file with the entire scheme
	std::ofstream  fh;
	fh.open((fn).c_str(), std::ios::out);
	if (!fh)
		REPORT_ERROR( (std::string)"Scheme::write: Cannot write file: " + fn);

	MetaDataTable MDgeneral;
	MDgeneral.setName("scheme_general");
	MDgeneral.setIsList(true);
	MDgeneral.addObject();
	MDgeneral.setValue(EMDL_SCHEME_GENERAL_NAME, name);
	MDgeneral.setValue(EMDL_SCHEME_GENERAL_CURRENT_NODE, current_node);
	MDgeneral.write(fh);

	if (schemer_global_floats.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("scheme_floats");
		std::map<std::string, SchemerFloatVariable>::iterator it;
		for ( it = schemer_global_floats.begin(); it != schemer_global_floats.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEME_VAR_FLOAT_NAME, it->first);
			MD.setValue(EMDL_SCHEME_VAR_FLOAT_VALUE, it->second.value);
			MD.setValue(EMDL_SCHEME_VAR_FLOAT_ORI_VALUE, it->second.original_value);
		}
		MD.write(fh);
	}

	if (schemer_global_bools.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("scheme_bools");
		std::map<std::string, SchemerBooleanVariable>::iterator it;
		for ( it = schemer_global_bools.begin(); it != schemer_global_bools.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEME_VAR_BOOL_NAME, it->first);
			MD.setValue(EMDL_SCHEME_VAR_BOOL_VALUE, it->second.value);
			MD.setValue(EMDL_SCHEME_VAR_BOOL_ORI_VALUE, it->second.original_value);
		}
		MD.write(fh);
	}

	if (schemer_global_strings.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("scheme_strings");
		std::map<std::string, SchemerStringVariable>::iterator it;
		for ( it = schemer_global_strings.begin(); it != schemer_global_strings.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEME_VAR_STRING_NAME, it->first);
			MD.setValue(EMDL_SCHEME_VAR_STRING_VALUE, it->second.value);
			MD.setValue(EMDL_SCHEME_VAR_STRING_ORI_VALUE, it->second.original_value);
		}
		MD.write(fh);
	}

	if (schemer_global_operators.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("scheme_operators");
		std::map<std::string, SchemerOperator>::iterator it;
		for ( it = schemer_global_operators.begin(); it != schemer_global_operators.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEME_OPERATOR_NAME, it->first);
			MD.setValue(EMDL_SCHEME_OPERATOR_TYPE, it->second.type);
			MD.setValue(EMDL_SCHEME_OPERATOR_OUTPUT, it->second.output );
			MD.setValue(EMDL_SCHEME_OPERATOR_INPUT1, it->second.input1 );
			MD.setValue(EMDL_SCHEME_OPERATOR_INPUT2, it->second.input2 );
		}
		MD.write(fh);
	}

	if (jobs.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("scheme_jobs");
		std::map<std::string, SchemerJob>::iterator it;
		for ( it = jobs.begin(); it != jobs.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEME_JOB_ORI_NAME, it->first);
			MD.setValue(EMDL_SCHEME_JOB_NAME, it->second.current_name);
			MD.setValue(EMDL_SCHEME_JOB_MODE, it->second.mode);
			MD.setValue(EMDL_SCHEME_JOB_HAS_STARTED, it->second.job_has_started);
		}
		MD.write(fh);
	}

	if (edges.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("scheme_edges");
		for (int i = 0; i < edges.size(); i++)
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEME_EDGE_INPUT, edges[i].inputNode);
			MD.setValue(EMDL_SCHEME_EDGE_OUTPUT, edges[i].outputNode);
			MD.setValue(EMDL_SCHEME_EDGE_IS_FORK, edges[i].is_fork);
			MD.setValue(EMDL_SCHEME_EDGE_OUTPUT_TRUE, edges[i].outputNodeTrue);
			MD.setValue(EMDL_SCHEME_EDGE_BOOLEAN, edges[i].myBooleanVariable);
		}
		MD.write(fh);
	}

	// Close the file handler
	fh.close();

	if (do_lock)
	{

#ifdef DEBUG_LOCK
		std::cerr << " write scheme: now deleting " << fn_lock << std::endl;
#endif

		if (!exists(fn_lock))
			REPORT_ERROR("ERROR: Scheme::write was expecting a file called "+fn_lock+ " but it is no longer there.");
		unlock();
	}

	// Touch a file to indicate to the GUI that the pipeline has just changed
	//FileName mychanged = name+SCHEME_HAS_CHANGED;
	//touch(mychanged);
}

// Reset all variables to their original value
void Scheme::reset()
{
	{
		std::map<std::string, SchemerFloatVariable>::iterator it;
		for ( it = schemer_global_floats.begin(); it != schemer_global_floats.end(); it++ )
			it->second.value = it->second.original_value;
	}

	{
		std::map<std::string, SchemerBooleanVariable>::iterator it;
		for ( it = schemer_global_bools.begin(); it != schemer_global_bools.end(); it++ )
	        it->second.value = it->second.original_value;
	}

	{
		std::map<std::string, SchemerStringVariable>::iterator it;
		for ( it = schemer_global_strings.begin(); it != schemer_global_strings.end(); it++ )
	        it->second.value = it->second.original_value;
	}

	{
		std::map<std::string, SchemerJob>::iterator it;
		for ( it = jobs.begin(); it != jobs.end(); it++ )
		{
			it->second.current_name = it->first;
			it->second.job_has_started = false;
		}
	}

	if (edges.size() > 0)
		current_node = edges[0].inputNode;
	else
		current_node = "undefined";
}

bool Scheme::isNode(std::string _name)
{
	// is this either an operator or a job?
	return (jobs.find(_name) != jobs.end() || schemer_global_operators.find(_name) != schemer_global_operators.end());
}

bool Scheme::isJob(std::string _name)
{
	return (jobs.find(_name) != jobs.end());
}

bool Scheme::isOperator(std::string _name)
{
	return isSchemeOperator(_name);
}

void Scheme::setVariable(std::string name, FileName value)
{
	float floatval;
	if (value != "" && sscanf(value.c_str(), "%f", &floatval)) // is this a number?
	{
		if (isFloatVariable(name)) setFloatVariableValue(name, floatval);
		else addFloatVariable(name, floatval);
	}
	else if (value == "true" || value == "True" || value == "false" || value == "False") // or a boolean?
	{
		bool myval = (value == "true" || value == "True");
		if (isBooleanVariable(name)) setBooleanVariableValue(name, myval);
		else addBooleanVariable(name, myval);
	}
	else
	{
		if (isStringVariable(name)) setStringVariableValue(name, value);
		else addStringVariable(name, value);
	}
}

void Scheme::setOriginalVariable(std::string name, FileName value)
{
	float floatval;
	if (value != "" && sscanf(value.c_str(), "%f", &floatval)) // is this a number?
	{
		if (isFloatVariable(name)) setFloatOriginalVariableValue(name, floatval);
		else addFloatVariable(name, floatval);
	}
	else if (value == "true" || value == "True" || value == "false" || value == "False") // or a boolean?
	{
		bool myval = (value == "true" || value == "True");
		if (isBooleanVariable(name)) setBooleanOriginalVariableValue(name, myval);
		else addBooleanVariable(name, myval);
	}
	else
	{
		if (isStringVariable(name)) setStringOriginalVariableValue(name, value);
		else addStringVariable(name, value);
	}
}

void Scheme::addFloatVariable(std::string name, RFLOAT value)
{
	if (isFloatVariable(name))
		REPORT_ERROR("ERROR: trying to add a float variable with a name that already exists: " + name);

	SchemerFloatVariable myvar(value, value);
	schemer_global_floats[name] = myvar;
}

void Scheme::addBooleanVariable(std::string name, bool value)
{
	if (isBooleanVariable(name))
		REPORT_ERROR("ERROR: trying to add a boolean variable with a name that already exists: " + name);

	SchemerBooleanVariable myvar(value, value);
	schemer_global_bools[name] = myvar;
}

void Scheme::addStringVariable(std::string name, FileName value)
{
	if (isStringVariable(name))
		REPORT_ERROR("ERROR: trying to add a string variable with a name that already exists: " + name);

	SchemerStringVariable myvar(value, value);
	schemer_global_strings[name] = myvar;
}

float Scheme::getFloatVariableValue(std::string name)
{
	if (!isFloatVariable(name))
		REPORT_ERROR("ERROR: cannot find float variable with name:" + name);
	return schemer_global_floats[name].value;
}

float Scheme::getFloatOriginalVariableValue(std::string name)
{
	if (!isFloatVariable(name))
		REPORT_ERROR("ERROR: cannot find float variable with name:" + name);
	return schemer_global_floats[name].original_value;
}

void Scheme::setFloatVariableValue(std::string name, RFLOAT val)
{
	if (!isFloatVariable(name))
		REPORT_ERROR("ERROR: cannot find float variable with name:" + name);
	schemer_global_floats[name].value = val;
	schemer_global_floats[name].value = val;
}

void Scheme::setFloatOriginalVariableValue(std::string name, RFLOAT val)
{
	if (!isFloatVariable(name))
		REPORT_ERROR("ERROR: cannot find float variable with name:" + name);
	schemer_global_floats[name].original_value = val;
}

bool Scheme::getBooleanVariableValue(std::string name)
{
	if (!isBooleanVariable(name))
		REPORT_ERROR("ERROR: cannot find boolean variable with name:" + name);
	return schemer_global_bools[name].value;
}

bool Scheme::getBooleanOriginalVariableValue(std::string name)
{
	if (!isBooleanVariable(name))
		REPORT_ERROR("ERROR: cannot find boolean variable with name:" + name);
	return schemer_global_bools[name].original_value;
}

void Scheme::setBooleanVariableValue(std::string name, bool val)
{
	if (!isBooleanVariable(name))
		REPORT_ERROR("ERROR: cannot find boolean variable with name:" + name);
	schemer_global_bools[name].value = val;
}

void Scheme::setBooleanOriginalVariableValue(std::string name, bool val)
{
	if (!isBooleanVariable(name))
		REPORT_ERROR("ERROR: cannot find boolean variable with name:" + name);
	schemer_global_bools[name].original_value = val;
}

std::string Scheme::getStringVariableValue(std::string name)
{
	if (!isStringVariable(name))
		REPORT_ERROR("ERROR: cannot find string variable with name:" + name);
	return schemer_global_strings[name].value;
}

std::string Scheme::getStringOriginalVariableValue(std::string name)
{
	if (!isStringVariable(name))
		REPORT_ERROR("ERROR: cannot find string variable with name:" + name);
	return schemer_global_strings[name].original_value;
}

void Scheme::setStringVariableValue(std::string name, std::string val)
{
	if (!isStringVariable(name))
		REPORT_ERROR("ERROR: cannot find string variable with name:" + name);
	schemer_global_strings[name].value = val;
}

void Scheme::setStringOriginalVariableValue(std::string name, std::string val)
{
	if (!isStringVariable(name))
		REPORT_ERROR("ERROR: cannot find string variable with name:" + name);
	schemer_global_strings[name].original_value = val;
}

std::string Scheme::getVariableValueAsString(std::string name)
{
	if (isStringVariable(name)) return schemer_global_strings[name].value;
	else if (isBooleanVariable(name)) return (schemer_global_bools[name].value) ? "True" : "False";
	else if (isFloatVariable(name)) return floatToString(schemer_global_floats[name].value);
	else REPORT_ERROR("Scheme::getVariableValueAsString: no variable named " + name);
}

void Scheme::setOperatorParameters(std::string name, std::string _type, std::string _input1, std::string _input2, std::string _output)
{
	if (!isOperator(name))
		REPORT_ERROR("ERROR: cannot find operator with name:" + name);
	// Just make a new one, so all the check are done automatically...
	SchemerOperator myop(_type, _input1, _input2, _output);
	schemer_global_operators[name] = myop;
}

void Scheme::getOperatorParameters(std::string name, std::string &_type, std::string &_input1, std::string &_input2, std::string &_output)
{
	if (!isOperator(name))
		REPORT_ERROR("ERROR: cannot find operator with name:" + name);
	_type = schemer_global_operators[name].type;
	_input1 = schemer_global_operators[name].input1;
	_input2 = schemer_global_operators[name].input2;
	_output = schemer_global_operators[name].output;
}

std::map<std::string, SchemerFloatVariable> Scheme::getCurrentFloatVariables()
{
	return schemer_global_floats;
}

std::map<std::string, SchemerBooleanVariable> Scheme::getCurrentBooleanVariables()
{
	return schemer_global_bools;
}

std::map<std::string, SchemerStringVariable> Scheme::getCurrentStringVariables()
{
	return schemer_global_strings;
}

std::map<std::string, SchemerOperator> Scheme::getCurrentOperators()
{
	return schemer_global_operators;
}

SchemerOperator Scheme::initialiseOperator(std::string type, std::string input_name, std::string input2_name,
		std::string output_name, std::string &error_message)
{
	SchemerOperator myop;
	error_message = myop.initialise(type, input_name, input2_name, output_name);
	return myop;
}

void Scheme::addOperator(SchemerOperator &myop, std::string &myname)
{
	if (name == "") myname = myop.getName();
	schemer_global_operators[myname] = myop;
}

void Scheme::addJob(RelionJob &myjob, std::string jobname, std::string mode)
{

	//remove spaces from jobname
	for (int i = 0; i < jobname.length(); i++)
	{
		if (jobname[i] == ' ') jobname[i] = '_';
	}

	// Check whether the jobname is unique
	if (isNode(jobname))
		REPORT_ERROR("ERROR: trying to add a JobNode that already exists: " + jobname);
	// Now add this job to the local Scheme
	std::string error_message;
	std::vector<std::string> commands;
	std::string final_command;
	std::string output_name = name + jobname + '/';

	// Save a copy of the job in the Schemes directory
	mktree(name + jobname);

	myjob.write(output_name);

	if (!myjob.getCommands(output_name, commands, final_command, false, 1, error_message))
		REPORT_ERROR("ERROR in getting commands for schemed job: " + error_message);

	SchemerJob mynode(jobname, mode, false);
	jobs[jobname] = mynode;
}

void Scheme::removeVariable(std::string name)
{
	// Remove any operators with this variable in it
	removeOperatorsWithThisInputOrOutput(name);

	if (isBooleanVariable(name))
	{
		schemer_global_bools.erase(name);
		// Also remove forks with this boolean variable
		removeEdgesWithThisInputOutputOrBoolean(name);
	}
	else if (isFloatVariable(name)) schemer_global_floats.erase(name);
	else if (isStringVariable(name)) schemer_global_strings.erase(name);
	else REPORT_ERROR("ERROR: cannot find variable to erase: " + name);
}

void Scheme::removeEdgesWithThisInputOutputOrBoolean(std::string name)
{
	std::vector<SchemerEdge> new_edges;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i].inputNode != name && edges[i].outputNode != name &&
			edges[i].outputNodeTrue != name && edges[i].myBooleanVariable != name)
		{
			new_edges.push_back(edges[i]);
		}
	}
	edges= new_edges;
}

void Scheme::removeOperator(std::string name)
{
	if (isOperator(name)) schemer_global_operators.erase(name);
	else REPORT_ERROR("ERROR: cannot find operator to erase: " + name);

	// Also remove any edges that input/output with this operator
	removeEdgesWithThisInputOutputOrBoolean(name);
}

void Scheme::removeOperatorsWithThisInputOrOutput(std::string name)
{
	std::map<std::string, SchemerOperator> new_operators;

	std::map<std::string, SchemerOperator>::iterator it;
	for ( it = schemer_global_operators.begin(); it != schemer_global_operators.end(); it++ )
		if (!(it->second.input1 == name || it->second.input2 == name ||  it->second.output == name))
		{
			// keep
			new_operators[it->first] = it->second;
		}
		else
		{
			// also remove edges with this operator in it
			removeEdgesWithThisInputOutputOrBoolean(it->first);
		}

	schemer_global_operators = new_operators;
}

void Scheme::removeJob(std::string name)
{
	if (isJob(name)) jobs.erase(name);
	else REPORT_ERROR("ERROR: cannot find job to erase: " + name);

	// Also remove any edges that input/output with this job
	removeEdgesWithThisInputOutputOrBoolean(name);
}

void Scheme::removeEdge(int idx)
{
	edges.erase(edges.begin()+idx);
}

void Scheme::copy(FileName newname)
{
	// Make sure newname ends with a slash
	if (newname[newname.length()-1] != '/') newname += "/";

	// Make the output directory,
	mktree(newname);

	// Replace original name in all stringVariables
	for (std::map<std::string, SchemerStringVariable>::iterator it = schemer_global_strings.begin(); it != schemer_global_strings.end(); it++ )
	{
		(it->second.value).replaceAllSubstrings(name, newname);
		(it->second.original_value).replaceAllSubstrings(name, newname);
	}

	// Replace all names in the pipeliner jobs
	for (std::map<std::string, SchemerJob>::iterator it = jobs.begin(); it != jobs.end(); it++ )
	{
		RelionJob myjob;
		bool dummy;
		if (!myjob.read(name + it->first + '/', dummy, true))
			REPORT_ERROR("There was an error reading job: " + it->first);

		for (std::map<std::string,JobOption>::iterator it2 = myjob.joboptions.begin(); it2 != myjob.joboptions.end(); ++it2)
		{
			FileName mystring = (it2->second).value;
			if (mystring.contains(name))
			{
				mystring.replaceAllSubstrings(name, newname);
				(it2->second).value = mystring;
			}
		}

		// Write the new job in the new directory
		std::string mydir = newname + it->first + '/';
		mktree(mydir);
		myjob.write(mydir);
	}

	// Change the name itself
	setName(newname);

	// And write the new scheme.star and scheule_pipeline.star files
	write();
}

void schemerSendEmail(std::string message, std::string subject)
{
	if (isStringVariable("email"))
	{
		std::string command = "echo \"" + message + "\" | mail -s \"" + subject + "\" -r RELION " + schemer_global_strings["email"].value;
		int res = system(command.c_str());
	}
}

bool Scheme::checkUniqueInput(std::string inputnode_name)
{
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i].inputNode == inputnode_name) return false;
	}

	return true;
}

void Scheme::addEdge(std::string inputnode_name, std::string outputnode_name)
{
	SchemerEdge myval(inputnode_name, outputnode_name);
	edges.push_back(myval);
}

void Scheme::addFork(std::string inputnode_name, std::string mybool_name, std::string outputnode_name_if_true, std::string outputnode_name)
{
	SchemerEdge myval(inputnode_name, outputnode_name, true, mybool_name, outputnode_name_if_true);
	edges.push_back(myval);
}

bool Scheme::isValid()
{
	// TODO: check if duplicate edges, forks or schemer_global_operators exist....

	// Check Schemer ends with an exit

	return false; // to be implemented
}

std::string Scheme::getNextNode()
{
	std::string result= "undefined";
	if (current_node == "undefined")
	{
		REPORT_ERROR("ERROR: the current_node is not defined...");
	}
	else
	{
		for (int i = 0; i < edges.size(); i++)
		{
			if (edges[i].inputNode == current_node)
			{
				result = edges[i].getOutputNode();
			}
		}
	}
	return result;
}

std::string Scheme::getPreviousNode()
{
	std::string result= "undefined";
	if (current_node == "undefined")
		REPORT_ERROR("ERROR: cannot return previous node, as the current node is undefined or equal to the original start node...");

	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i].getOutputNode() == current_node)
		{
			result = edges[i].inputNode;
			return result;
		}
	}
	return result;
}

bool Scheme::gotoNextNode()
{
	current_node = getNextNode();
	if (verb > 1) std::cout << " ++ Setting current node to: " << current_node << std::endl;

	// Write out current status, but maintain lock on the directory!
	write();

	return (current_node == "undefined") ? false : true;
}

bool Scheme::gotoNextJob()
{
	// This loops through the next Nodes until encountering a JOB
	while (gotoNextNode())
	{
		if (pipeline_control_check_abort_job())
		{
			write(DO_LOCK);
			exit(RELION_EXIT_ABORTED);
		}

		if (isOperator(current_node))
		{
			if (!executeOperator(current_node)) return false;
		}
		else // this is a job, get its current_name and options
		{
			return true;
		}
	}

	return false;
}


bool Scheme::changeStringForJobnames(FileName &mystring, FileName current_node)
{

	// Check for any strings containing the 'name' of this Scheme; if so, replace by current_name of the corresponding job
	FileName original_string = mystring;
	if (mystring.contains(name))
	{

		// Remove leading directory and tailing slash to get the process current_name in the pipeline_schemer
		FileName my_ori_name = (mystring.afterFirstOf(name)).beforeLastOf("/");
		// find that process in the nodes, and get its current current_name
		std::string my_current_name =jobs[my_ori_name].current_name;
		std::string to_replace = name  + my_ori_name + '/';
		mystring.replaceAllSubstrings(to_replace, my_current_name);
		std::cout << " + " << current_node << ": " << original_string  << " -> " << mystring << std::endl;

		return true;

	}
	// Also check for any other Scheme that might run in parallel
	else if (mystring.contains("Schemes/"))
	{

		// Remove leading directory and tailing slash to get the process current_name in the pipeline_schemer
		FileName my_scheme_name = (mystring.afterFirstOf("Schemes/")).beforeFirstOf("/");
		FileName my_scheme_star = "Schemes/"+my_scheme_name+"/scheme.star";
		if (exists(my_scheme_star)) {
            // Remove leading directory and tailing slash to get the process current_name in the pipeline_schemer
            FileName my_ori_name = (mystring.afterFirstOf(my_scheme_name + "/")).beforeLastOf("/");

            // Read only the jobs table from the other scheme;
            // otherwise global variables like schemer_global_floats, strings etc get overwritten!
            // But be careful, as STAR file might just be written out be the other Schemer
            // Therefore try at least 3 times
            MetaDataTable MDjobs;

            int itry = 0;
            bool have_jobs = false;
            while (!have_jobs)
            {

                MDjobs.read(my_scheme_star, "scheme_jobs");
                have_jobs = (MDjobs.numberOfObjects() > 0);

                if (!have_jobs)
                {
                    // Wait one second before trying again
                    sleep(1);

                    itry++;
                    if (itry > 5) REPORT_ERROR("ERROR: trying to read scheme_jobs table from " + my_scheme_star +
                                               " but cannot find any jobs in it....");
                }
            }

            FileName job_current_name, job_ori_name, my_current_name = "undefined";
            FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDjobs)
            {
                MDjobs.getValue(EMDL_SCHEME_JOB_ORI_NAME, job_ori_name);
                MDjobs.getValue(EMDL_SCHEME_JOB_NAME, job_current_name);
                if (job_ori_name == my_ori_name)
                {
                    my_current_name = job_current_name;
                }
            }

            if (my_current_name[my_current_name.length()-1] != '/')
                    my_current_name += "/";

            if (my_current_name == "undefined")
            {
                    REPORT_ERROR("ERROR: cannot find job " + my_ori_name + " in scheme: " + my_scheme_star);
            }
            std::string to_replace = "Schemes/" + my_scheme_name + "/" + my_ori_name + "/";
            mystring.replaceAllSubstrings(to_replace, my_current_name);

            std::cout << " ++ in " << current_node << ": " << original_string  << " -> " << mystring << std::endl;

            return true;
		}
		else
		{
			REPORT_ERROR("ERROR: cannot find Scheme STAR file: " + my_scheme_star);
		}

	}

	return false;

}

// Modify an operator to set variables from the Schemer
bool Scheme::executeOperator(FileName current_node)
{
	SchemerOperator my_op = schemer_global_operators[current_node];

	// Use temporary variables to store the current values, then perform the operation and delete the temporary ones again
	// This is so that the original values in input1,2 or output are retained
	if (isStringVariable(my_op.input1))
	{
		FileName mystring = schemer_global_strings[my_op.input1].value;
		if (changeStringForJobnames(mystring, current_node))
		{
			addStringVariable("xxx_tmp_input1", mystring);
			my_op.input1 = "xxx_tmp_input1";
		}
	}
	if (isStringVariable(my_op.input2))
	{
		FileName mystring = schemer_global_strings[my_op.input2].value;
		if (changeStringForJobnames(mystring, current_node))
		{
			addStringVariable("xxx_tmp_input2", mystring);
			my_op.input1 = "xxx_tmp_input2";
		}
	}
	if (isStringVariable(my_op.output))
	{
		FileName mystring = schemer_global_strings[my_op.output].value;
		if (changeStringForJobnames(mystring, current_node))
		{
			addStringVariable("xxx_tmp_output", mystring);
			my_op.output = "xxx_tmp_output";
		}
	}

	bool my_success = my_op.performOperation();

	if (isStringVariable("xxx_tmp_input1")) removeVariable("xxx_tmp_input1");
	if (isStringVariable("xxx_tmp_input2")) removeVariable("xxx_tmp_input2");
	if (isStringVariable("xxx_tmp_output"))
	{
		// Set output in the variable from the original operator, and then remove tmp_output
		std::string myoutputvariable = schemer_global_operators[current_node].output;
		schemer_global_strings[myoutputvariable].value = schemer_global_strings["xxx_tmp_output"].value;
		removeVariable("xxx_tmp_output");
	}

	if (verb > 0 && schemer_global_operators[current_node].output != "undefined")
	{
		std::cout << " + " << current_node << ": " << schemer_global_operators[current_node].output << " -> " <<getVariableValueAsString(schemer_global_operators[current_node].output) << std::endl;
	}

	return my_success;

}


// Modify a job to set variables from the Schemer
RelionJob Scheme::prepareJob(FileName current_node)
{

	RelionJob myjob;
	bool dummy, is_continue, needs_to_add_job = false;

	// Always re-read from the original job.star from the Scheme directory,
	// because $$ variables and input jobnames may have been replaced when continuing existing jobs
	if (!myjob.read(name + current_node + '/', dummy, true))
        REPORT_ERROR("ERROR: there was a problem reading job " + name + current_node);

    // Check whether there are any joboption values with a jobname from one of the processes in this Schemer
	// And replace these by their corresponding 'current_name'
	for (std::map<std::string,JobOption>::iterator it=myjob.joboptions.begin(); it!=myjob.joboptions.end(); ++it)
	{
		FileName mystring = (it->second).value;

		// Check for options with a value containing $$, and replace with the current value of the corresponding Variable
		if (mystring.contains("$$"))
		{
			std::vector< std::string > myvars;
			bool has_found = false;
			while (mystring.contains("$$"))
			{
				has_found = true;
				FileName before = mystring.beforeFirstOf("$$");
				FileName after = mystring.afterFirstOf("$$");
				std::vector< std::string > splits;
				int nr_splits = splitString(after, " ", splits);
				if (splits.size() == 0)
				{
					REPORT_ERROR(" ERROR: cannot find anything after $$ sign in string: " + mystring);
				}
				std::string mypat = splits[0];
				myvars.push_back(mypat);

				// Found an option that needs replacement! Now find which variable to insert
				std::string my_value;
				if (isBooleanVariable(mypat))
				{
					if (myjob.joboptions[it->first].joboption_type != JOBOPTION_BOOLEAN)
					{
						REPORT_ERROR(" ERROR: trying to set a BooleanVariable: " + mypat + " into a non-boolean option: " + it->first);
					}

					my_value = (schemer_global_bools[mypat].value) ? "Yes" : "No";
				}
				else if (isFloatVariable(mypat))
				{
					if (myjob.joboptions[it->first].joboption_type == JOBOPTION_BOOLEAN)
					{
						REPORT_ERROR(" ERROR: trying to set FloatVariable: " + mypat + " into a boolean option: " + it->first);
					}

					my_value = floatToString(schemer_global_floats[mypat].value);
				}
				else if (isStringVariable(splits[0]))
				{
					if (myjob.joboptions[it->first].joboption_type == JOBOPTION_BOOLEAN)
					{
						REPORT_ERROR(" ERROR: trying to set StringVariable: " + mypat + " into a boolean option: " + it->first);
					}
					if (myjob.joboptions[it->first].joboption_type == JOBOPTION_SLIDER)
					{
						REPORT_ERROR(" ERROR: trying to set StringVariable: " + mypat + " into a slider option: " + it->first);
					}

					my_value = schemer_global_strings[mypat].value;
				}
				else
				{
					REPORT_ERROR(" ERROR: variable in job is not part of this Scheme: " + mypat);
				}

				mystring = before + my_value;
				for (int i = 1; i < splits.size(); i++)
					mystring += " " + splits[i];
			}

			if (has_found)
			{
				myjob.joboptions[it->first].value = mystring;
				std::string myvarsstr = "";
				for (int i = 0; i < myvars.size(); i++)
					myvarsstr+= myvars[i]+ " ";
				if (verb > 2) std::cout << " +++ Setting joboption " << it->first << " to " << mystring << " based on variable(s): " << myvarsstr<< std::endl;
			}

		} //end if mystring contains $$

		// Also check for jobnames
		if (changeStringForJobnames(mystring, current_node)) myjob.joboptions[it->first].value = mystring;

	}

	return myjob;
}



void Scheme::run(PipeLine &pipeline)
{
	time_config();

	pipeline_control_delete_exit_files();


	if (current_node == "undefined")
	{
		if (edges.size() > 0)
			current_node = edges[0].inputNode;
		else
		REPORT_ERROR("No edges defined yet...");
	}

    if (verb > 0)
    {
    	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++"  << std::endl;
    	std::cout << " + Starting execution of scheme at: " << current_node << std::endl;
    }

    // If we start from operators instead of jobs, then execute the operator and proceed until the next Job
	bool has_more_jobs = true;
    if (isOperator(current_node))
    {
		if (!executeOperator(current_node)) REPORT_ERROR("ERROR: something went wrong with execution of the initial operator...");
		has_more_jobs = gotoNextJob();
    }

	// go through all nodes
	bool is_ok = true;
	while (has_more_jobs)
	{
		// Abort mechanism
		if (pipeline_control_check_abort_job())
		{
			write(DO_LOCK);
			exit(RELION_EXIT_ABORTED);
		}

		RelionJob myjob = prepareJob(current_node);

		bool dummy, is_continue = false;
		int current_job;
		if (!jobs[current_node].job_has_started || jobs[current_node].mode == SCHEME_NODE_JOB_MODE_NEW)
		{

			// Now add this job to the pipeline we will actually be running in
			current_job = pipeline.addScheduledJob(myjob, "", false); // false means dont write hidden guifile

			// Set the current_name of the current node now
			jobs[current_node].current_name = pipeline.processList[current_job].name;
			if (verb > 0) std::cout << " + Creating new Job: " << jobs[current_node].current_name << " from Node: " << current_node << std::endl;

		}
		else if (jobs[current_node].mode == SCHEME_NODE_JOB_MODE_CONTINUE)
		{

			is_continue = (jobs[current_node].mode == SCHEME_NODE_JOB_MODE_CONTINUE);
			current_job = pipeline.findProcessByName(jobs[current_node].current_name);

			if (current_job < 0)
				REPORT_ERROR("ERROR: RunScheme cannot find process with name: " + jobs[current_node].current_name);

		}
		else
		{
			REPORT_ERROR("ERROR: unrecognised mode for running a new process: " + jobs[current_node].mode);
		}

		// Check whether the input nodes are there, before executing the job
		for (long int inode = 0; inode < pipeline.processList[current_job].inputNodeList.size(); inode++)
		{
			long int mynode = pipeline.processList[current_job].inputNodeList[inode];
			int itry = 0;
			while (!exists(pipeline.nodeList[mynode].name))
			{
				std::cerr << " + Warning " << pipeline.nodeList[mynode].name << " does not exist. Waiting 10 seconds ... " << std::endl;
				sleep(10);

				// Abort mechanism
				if (pipeline_control_check_abort_job())
				{
					write(DO_LOCK);
					exit(RELION_EXIT_ABORTED);
				}
				else if (itry > 3)
				{
					std::cout << " + Gave up on waiting for " << pipeline.nodeList[mynode].name << ". Aborting ... " << std::endl;
					write(DO_LOCK);
					exit(RELION_EXIT_ABORTED);
				}

				itry++;
			}
		}

		// Now actually run the Schemed job
		std::string error_message;
		if (verb > 0)
		{
			time_t my_time = time(NULL);
			std::cout << " + Executing Job: " << jobs[current_node].current_name << " at " << ctime(&my_time);
		}
		jobs[current_node].job_has_started = true;

		// last false: don't write hidden GUI files, so defaults in Schemes don't mess up defaults for environment variables in later execution of RELION GUI
		if (!pipeline.runJob(myjob, current_job, false, is_continue, true, error_message, false))
			REPORT_ERROR(error_message);

		// Write out current status, but maintain lock on the directory!
		write();

		// Wait for job to finish
		bool is_failure = false;
		bool is_aborted = false;
		pipeline.waitForJobToFinish(current_job, is_failure, is_aborted);


		std::string message = "";
		if (is_failure) message = " + Stopping scheme due to job " + jobs[current_node].current_name + " failing with an error ...";
		else if (is_aborted) message = " + Stopping scheme due to user abort of job " + jobs[current_node].current_name + " ...";
		if (message != "")
		{
			schemerSendEmail(message, "Scheme: " + name);
			std::cout << message << std::endl;
			is_ok = false;
			break;
		}

		has_more_jobs = gotoNextJob();
	} // end while has_more_jobs

	if (is_ok) schemerSendEmail("Finished successfully!", "Scheme: " + name);

	if (verb > 0)
	{
		if (exists(name + RELION_JOB_ABORT_NOW))
			std::cout << " + Found an ABORT signal... " << std::endl;
		std::cout << " + Schemer " << name << " stops now... " << std::endl;
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++"  << std::endl << std::endl;
	}
}

void Scheme::unlock()
{
	FileName name_wo_dir = name;
	name_wo_dir = name_wo_dir.beforeLastOf("/");
	FileName dir_lock = ".relion_lock_scheme_" + name_wo_dir.afterLastOf("/");
	FileName fn_lock = dir_lock + "/lock_scheme";;

	if (exists(fn_lock))
	{
		if (std::remove(fn_lock.c_str()))
			REPORT_ERROR("ERROR: in removing lock file "+fn_lock);
		if (rmdir(dir_lock.c_str()))
			REPORT_ERROR("ERROR: in removing lock directory "+dir_lock);
	}
}

void Scheme::abort()
{
	std::cout << " + Aborting scheme while at: " << current_node << std::endl;
	// Only abort the job if current_name is no longer the original name from the Schemer!
	if (isJob(current_node) && current_node != jobs[current_node].current_name)
	{
		touch(jobs[current_node].current_name + RELION_JOB_ABORT_NOW);
		std::cout << " ++ Touched file: " << jobs[current_node].current_name << RELION_JOB_ABORT_NOW << std::endl;
	}
	touch(name + RELION_JOB_ABORT_NOW);
	std::cout << " ++ Touched file: " << name << RELION_JOB_ABORT_NOW << std::endl;
	std::cout << " ++ Now wait for the job and the schemer to abort ..." << std::endl;
}
