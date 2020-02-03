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

// one global timestamp...
static time_t annotated_time;
bool has_annotated_time = false;

// Global variables, but only with reach within this file!
std::map<std::string, SchedulerBooleanVariable> scheduler_global_bools;
std::map<std::string, SchedulerFloatVariable> scheduler_global_floats;
std::map<std::string, SchedulerStringVariable> scheduler_global_strings;
std::map<std::string, SchedulerOperator> scheduler_global_operators;

bool isBooleanVariable(std::string _name)
{
	return (scheduler_global_bools.find(_name) != scheduler_global_bools.end());
}

bool isFloatVariable(std::string _name)
{
	return (scheduler_global_floats.find(_name) != scheduler_global_floats.end());
}

bool isStringVariable(std::string _name)
{
	return (scheduler_global_strings.find(_name) != scheduler_global_strings.end());
}

bool isScheduleOperator(std::string _name)
{
	return (scheduler_global_operators.find(_name) != scheduler_global_operators.end());
}

SchedulerOperator::SchedulerOperator(std::string _type, std::string _input1, std::string _input2, std::string _output)
{
	std::string myerror = initialise(_type, _input1, _input2, _output);
	if (myerror != "") REPORT_ERROR(myerror);
}

std::string SchedulerOperator::initialise(std::string _type, std::string _input1, std::string _input2, std::string _output)
{
	type = _type;

	// Check output
	if ((type == SCHEDULE_BOOLEAN_OPERATOR_GT ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_LT ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_EQ ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_GE ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_LE ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_AND ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_OR ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_READ_STAR
		) && !isBooleanVariable(_output))
		return "ERROR: boolean operator does not have valid boolean output: " + _output;
	if ((type == SCHEDULE_FLOAT_OPERATOR_SET ||
		 type == SCHEDULE_FLOAT_OPERATOR_PLUS ||
		 type == SCHEDULE_FLOAT_OPERATOR_MINUS ||
		 type == SCHEDULE_FLOAT_OPERATOR_MULT ||
		 type == SCHEDULE_FLOAT_OPERATOR_DIVIDE ||
		 type == SCHEDULE_FLOAT_OPERATOR_ROUND ||
		 type == SCHEDULE_FLOAT_OPERATOR_COUNT_IMAGES ||
		 type == SCHEDULE_FLOAT_OPERATOR_COUNT_WORDS ||
		 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR ||
		 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MAX ||
		 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MIN ||
		 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_AVG ||
		 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX
		 ) && !isFloatVariable(_output))
		return "ERROR: float operator does not have valid float output: " + _output;
	if ((type == SCHEDULE_STRING_OPERATOR_READ_STAR ||
		 type == SCHEDULE_STRING_OPERATOR_JOIN ||
		 type == SCHEDULE_STRING_OPERATOR_BEFORE_FIRST ||
		 type == SCHEDULE_STRING_OPERATOR_AFTER_FIRST ||
		 type == SCHEDULE_STRING_OPERATOR_BEFORE_LAST ||
		 type == SCHEDULE_STRING_OPERATOR_AFTER_LAST ||
		 type == SCHEDULE_STRING_OPERATOR_GLOB ||
		 type == SCHEDULE_STRING_OPERATOR_NTH_WORD
		 )	&& ! isStringVariable(_output))
		return "ERROR: string operator does not have valid string output: " + _output;

	// Check input1
	if ((type == SCHEDULE_BOOLEAN_OPERATOR_AND ||
		type == SCHEDULE_BOOLEAN_OPERATOR_OR )  && !isBooleanVariable(_input1))
		return "ERROR: boolean operator does not have valid boolean input1: " + _input1;
	if ((type == SCHEDULE_FLOAT_OPERATOR_SET ||
		 type == SCHEDULE_FLOAT_OPERATOR_PLUS ||
		 type == SCHEDULE_FLOAT_OPERATOR_MINUS ||
		 type == SCHEDULE_FLOAT_OPERATOR_MULT ||
		 type == SCHEDULE_FLOAT_OPERATOR_DIVIDE ||
		 type == SCHEDULE_FLOAT_OPERATOR_ROUND
		 ) && !isFloatVariable(_input1))
		return "ERROR: float operator does not have valid float input1: " + _input1;
	if ((type == SCHEDULE_BOOLEAN_OPERATOR_READ_STAR ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS ||
		 type == SCHEDULE_FLOAT_OPERATOR_COUNT_IMAGES ||
		 type == SCHEDULE_FLOAT_OPERATOR_COUNT_WORDS ||
		 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR ||
		 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MAX ||
		 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MIN ||
		 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_AVG ||
		 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX ||
		 type == SCHEDULE_STRING_OPERATOR_READ_STAR ||
		 type == SCHEDULE_STRING_OPERATOR_JOIN ||
		 type == SCHEDULE_STRING_OPERATOR_BEFORE_FIRST ||
		 type == SCHEDULE_STRING_OPERATOR_AFTER_FIRST ||
		 type == SCHEDULE_STRING_OPERATOR_BEFORE_LAST ||
		 type == SCHEDULE_STRING_OPERATOR_AFTER_LAST ||
		 type == SCHEDULE_STRING_OPERATOR_GLOB ||
		 type == SCHEDULE_STRING_OPERATOR_NTH_WORD ||
		 type == SCHEDULE_OPERATOR_TOUCH_FILE ||
		 type == SCHEDULE_OPERATOR_COPY_FILE ||
		 type == SCHEDULE_OPERATOR_MOVE_FILE ||
		 type == SCHEDULE_OPERATOR_DELETE_FILE
		 ) && ! isStringVariable(_input1))
		return "ERROR: operator does not have valid string input1: " + _input1;

	// Check input2
	if ((type == SCHEDULE_BOOLEAN_OPERATOR_AND ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_OR
		 ) && !isBooleanVariable(_input2))
		return "ERROR: boolean operator does not have valid boolean input2: " + _input2;
	if ((type == SCHEDULE_BOOLEAN_OPERATOR_GT ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_LT ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_EQ ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_GE ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_LE ||
		 type == SCHEDULE_FLOAT_OPERATOR_PLUS ||
		 type == SCHEDULE_FLOAT_OPERATOR_MINUS ||
		 type == SCHEDULE_FLOAT_OPERATOR_MULT ||
		 type == SCHEDULE_FLOAT_OPERATOR_DIVIDE ||
		 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX ||
		 type == SCHEDULE_STRING_OPERATOR_NTH_WORD
		 ) && !(isFloatVariable(_input2) || isNumber(_input2)))
		return "ERROR: operator does not have valid number (float variable or text) input2: " + _input2;
	if ((type == SCHEDULE_OPERATOR_COPY_FILE ||
		 type == SCHEDULE_OPERATOR_MOVE_FILE ||
		 type == SCHEDULE_STRING_OPERATOR_BEFORE_FIRST ||
		 type == SCHEDULE_STRING_OPERATOR_AFTER_FIRST ||
		 type == SCHEDULE_STRING_OPERATOR_BEFORE_LAST ||
		 type == SCHEDULE_STRING_OPERATOR_AFTER_LAST ||
		 type == SCHEDULE_STRING_OPERATOR_JOIN
	 ) && ! isStringVariable(_input2))
	return "ERROR: operator does not have valid string input2: " + _input2;

	input1 = (_input1 == "") ? "undefined" : _input1;
	input2 = (_input2 == "") ? "undefined" : _input2;
	output = (_output == "") ? "undefined" : _output;

	return "";
}

// Separate comma-separated labels for table, input and output
void  SchedulerOperator::readFromStarFile() const
{
	MetaDataTable MD;
	std::string mystring, mystarfile, mytable;
	EMDLabel mylabel;

	// The localtion is always in input1
	mystring = scheduler_global_strings[input1].value;
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

	long idx = (isFloatVariable(input2)) ? ROUND(scheduler_global_floats[input2].value) : 0;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_READ_STAR ||
		type == SCHEDULE_FLOAT_OPERATOR_READ_STAR ||
		type == SCHEDULE_STRING_OPERATOR_READ_STAR)
	{

		if (EMDL::isDouble(mylabel))
		{
			RFLOAT fval;
			MD.getValue(mylabel, fval, idx);
			scheduler_global_floats[output].value = fval;
		}
		else if (EMDL::isInt(mylabel))
		{
			int ival;
			MD.getValue(mylabel, ival, idx);
			scheduler_global_floats[output].value = ROUND(ival);
		}
		else if (EMDL::isString(mylabel))
		{
			std::string val;
			MD.getValue(mylabel, val, idx);
			scheduler_global_strings[output].value = val;
		}
		else if (EMDL::isBool(mylabel))
		{
			bool val;
			MD.getValue(mylabel, val, idx);
			scheduler_global_bools[output].value = val;
		}
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MAX ||
			type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MIN ||
			type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_AVG ||
			type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX)
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

		if (type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MAX) scheduler_global_floats[output].value = mymax;
		else if (type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MIN) scheduler_global_floats[output].value = mymin;
		else if (type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_AVG) scheduler_global_floats[output].value = mysum;
		else if (type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX)
		{
			if (idx == 0) REPORT_ERROR("Give a positive or a negative value for input2 in sorted_idx: 1 is largest, -1 is smallest value");
			if (idx < 0) idx = ii + idx; // smallest value is numberOfObjects - 1
			else idx--; // now start counting at 0
			MultidimArray<long> sorted_idx;
			for_sorting.sorted_index(sorted_idx);
			scheduler_global_floats[output].value = DIRECT_MULTIDIM_ELEM(sorted_idx, idx);
		}
	}
}

bool SchedulerOperator::performOperation() const
{
	RFLOAT val2 = (isFloatVariable(input2)) ? scheduler_global_floats[input2].value : 0;

	if (type == SCHEDULE_BOOLEAN_OPERATOR_AND)
	{
		scheduler_global_bools[output].value = (scheduler_global_bools[input1].value && scheduler_global_bools[input2].value);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_OR)
	{
		scheduler_global_bools[output].value = (scheduler_global_bools[input1].value || scheduler_global_bools[input2].value);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_NOT)
	{
		scheduler_global_bools[output].value = (!(scheduler_global_bools[input1].value));
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_GT)
	{
		scheduler_global_bools[output].value = (scheduler_global_floats[input1].value > val2);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_LT)
	{
		scheduler_global_bools[output].value = (scheduler_global_floats[input1].value < val2);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_EQ)
	{
		scheduler_global_bools[output].value = (fabs(scheduler_global_floats[input1].value - val2) < 1E-8);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_GE)
	{
		scheduler_global_bools[output].value = (scheduler_global_floats[input1].value >= val2);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_LE)
	{
		scheduler_global_bools[output].value = (scheduler_global_floats[input1].value <= val2);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS)
	{
		scheduler_global_bools[output].value = (exists(scheduler_global_strings[input1].value));
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_READ_STAR ||
			 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR ||
			 type == SCHEDULE_STRING_OPERATOR_READ_STAR ||
			 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MAX ||
			 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MIN ||
			 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_AVG ||
			 type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX
			 )
	{
		readFromStarFile();
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_SET)
	{
		scheduler_global_floats[output].value = scheduler_global_floats[input1].value;
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_PLUS)
	{
		scheduler_global_floats[output].value = scheduler_global_floats[input1].value + val2;
	}
	else if (type ==SCHEDULE_FLOAT_OPERATOR_MINUS)
	{
		scheduler_global_floats[output].value = scheduler_global_floats[input1].value - val2;
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_MULT)
	{
		scheduler_global_floats[output].value = scheduler_global_floats[input1].value * val2;
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE)
	{
		scheduler_global_floats[output].value = scheduler_global_floats[input1].value / val2;
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_ROUND)
	{
		scheduler_global_floats[output].value = ROUND(scheduler_global_floats[input1].value);
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_COUNT_IMAGES)
	{
		ObservationModel obsmodel;
		MetaDataTable MDimg;
		std::string mytablename = (isStringVariable(input2)) ? scheduler_global_strings[input2].value : "particles";
		ObservationModel::loadSafely(scheduler_global_strings[input1].value, obsmodel, MDimg, mytablename);
		scheduler_global_floats[output].value = MDimg.numberOfObjects();
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_COUNT_WORDS)
	{
		// return zero for an undefined string
		if (scheduler_global_strings[input1].value == "undefined")
		{
			scheduler_global_floats[output].value = 0;
		}
		else
		{
			std::vector< std::string > splits;
			int nr_splits = splitString(scheduler_global_strings[input1].value, ",", splits);
			scheduler_global_floats[output].value = splits.size();
		}
	}
	else if (type == SCHEDULE_STRING_OPERATOR_JOIN)
	{
		scheduler_global_strings[output].value = scheduler_global_strings[input1].value + scheduler_global_strings[input2].value;
	}
	else if (type == SCHEDULE_STRING_OPERATOR_BEFORE_FIRST)
	{
		scheduler_global_strings[output].value = scheduler_global_strings[input1].value.beforeFirstOf(scheduler_global_strings[input2].value);
	}
	else if (type == SCHEDULE_STRING_OPERATOR_AFTER_FIRST)
	{
		scheduler_global_strings[output].value = scheduler_global_strings[input1].value.afterFirstOf(scheduler_global_strings[input2].value);
	}
	else if (type == SCHEDULE_STRING_OPERATOR_BEFORE_LAST)
	{
		scheduler_global_strings[output].value = scheduler_global_strings[input1].value.beforeLastOf(scheduler_global_strings[input2].value);
	}
	else if (type == SCHEDULE_STRING_OPERATOR_AFTER_LAST)
	{
		scheduler_global_strings[output].value = scheduler_global_strings[input1].value.afterLastOf(scheduler_global_strings[input2].value);
	}
	else if (type == SCHEDULE_STRING_OPERATOR_GLOB)
	{
		FileName input=scheduler_global_strings[input1].value;
		std::vector<FileName> files;
		input.globFiles(files);
		if (files.size() == 0)
		{
			scheduler_global_strings[output].value = "undefined";
		}
		else
		{
			scheduler_global_strings[output].value = files[0];
			for (int i = 1; i < files.size(); i++)
			{
				scheduler_global_strings[output].value += "," + files[i];
			}
		}
	}
	else if (type == SCHEDULE_STRING_OPERATOR_NTH_WORD)
	{

		std::vector< std::string > splits;
		int nr_splits = splitString(scheduler_global_strings[input1].value, ",", splits);
		int mypos = ROUND(val2);
		// for negative Ns, count from the back
		if (mypos < 0) mypos = splits.size() - mypos + 1;
		// Started counting at 1, but first element of vector is zero!
		mypos--;
		if (mypos >= splits.size() || mypos < 0)
		{
			scheduler_global_strings[output].value = "undefined";
		}
		else
		{
			scheduler_global_strings[output].value = splits[mypos];
		}
	}
	else if (type == SCHEDULE_OPERATOR_TOUCH_FILE)
	{
		std::cout << " + Touching: " << scheduler_global_strings[input1].value << std::endl;
		touch(scheduler_global_strings[input1].value);
	}
	else if (type == SCHEDULE_OPERATOR_COPY_FILE || type == SCHEDULE_OPERATOR_MOVE_FILE)
	{
		std::string mycommand;
		if (type == SCHEDULE_OPERATOR_COPY_FILE)
		{
			std::cout << " + Copying: " << scheduler_global_strings[input1].value << " to " << scheduler_global_strings[input2].value << std::endl;
			mycommand = "cp ";
		}
		else
		{
			std::cout << " + Moving: " << scheduler_global_strings[input1].value << " to " << scheduler_global_strings[input2].value << std::endl;
			mycommand = "mv ";
		}
		// Make output directory if it doesn't exist
		if (scheduler_global_strings[input2].value.contains("/"))
		{
			FileName mydirs = scheduler_global_strings[input2].value.beforeLastOf("/");
			std::string mycommand = "mkdir -p " + mydirs;
 			int res = system(mycommand.c_str());
		}
		// Execute the command
		mycommand += scheduler_global_strings[input1].value + " " + scheduler_global_strings[input2].value;
		int res = system(mycommand.c_str());
	}
	else if (type == SCHEDULE_OPERATOR_DELETE_FILE)
	{
		std::string mycommand = "rm -f " + scheduler_global_strings[input1].value;
		std::cout << " + Deleting: " << scheduler_global_strings[input1].value << std::endl;
		int res = system(mycommand.c_str());
	}
	else if (type == SCHEDULE_WAIT_OPERATOR_SINCE_LAST_TIME)
	{
		if (has_annotated_time)
		{
			time_t current_time = time(NULL);
			RFLOAT elapsed = current_time - annotated_time;
			// Also set output to elapsed seconds if it is a floatVariable
			if (isFloatVariable(output)) scheduler_global_floats[output].value = elapsed;
			RFLOAT wait_seconds =  scheduler_global_floats[input1].value - elapsed;
			if (wait_seconds > 0)
			{
				std::cout << " + Waiting for " << wait_seconds << " seconds ..." << std::endl;
				sleep(wait_seconds);
				std::cout << " + Finished waiting." << std::endl;
			}
			else
			{
				std::cout << " + Not waiting, as more than " << scheduler_global_floats[input1].value << " seconds have passed." << std::endl;
			}
		}
		annotated_time = time(NULL);
		has_annotated_time =true;
	}
	else if (type == SCHEDULE_EMAIL_OPERATOR)
	{

		time_t my_time = time(NULL);
		std::string mymessage = std::string(ctime(&my_time)) + "\n";
		mymessage += "input1: " + input1 + " = ";
		if (isStringVariable(input1)) mymessage += scheduler_global_strings[input1].value + "\n";
		else if (isBooleanVariable(input1)) mymessage += (scheduler_global_bools[input1].value) ? "True \n" : "False \n";
		else if (isFloatVariable(input1)) mymessage += floatToString(scheduler_global_floats[input1].value) + "\n";
		if (isBooleanVariable(input2) || isFloatVariable(input2) || isStringVariable(input2))
		{
			mymessage += "input2: " + input2 + " = ";
			if (isStringVariable(input2)) mymessage += scheduler_global_strings[input2].value + "\n";
			else if (isBooleanVariable(input2)) mymessage += (scheduler_global_bools[input2].value) ? "True \n" : "False \n";
			else if (isFloatVariable(input2)) mymessage += floatToString(scheduler_global_floats[input2].value) + "\n";
		}
		std::cout << " + Sending e-mail." << std::endl;

		schedulerSendEmail(mymessage);
	}
	else if (type == SCHEDULE_EXIT_OPERATOR)
	{
		std::cout << " + The schedule has reached an exit point ..." << std::endl;
		return false; // to exit the schedule
	}
	else
		REPORT_ERROR("ERROR: unrecognised Operator type:" + type);

	return true;
}

std::string SchedulerOperator::getName()
{
	if (type == SCHEDULE_BOOLEAN_OPERATOR_GT) return output + "=" + input1 + "_GT_" + input2;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_LT) return output + "=" + input1 + "_LT_" + input2;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_EQ) return output + "=" + input1 + "_EQ_" + input2;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_GE) return output + "=" + input1 + "_GE_" + input2;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_LE) return output + "=" + input1 + "_LE_" + input2;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_AND) return output + "=" + input1 + "_AND_" + input2;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_OR) return output + "=" + input1 + "_OR_" + input2;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_NOT) return output + "=" + "NOT_" + input1;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS) return output + "=" + "EXISTS_" + input1;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_READ_STAR) return output + "=" + "STAR_" + input1 + "_" + input2;
	if (type == SCHEDULE_FLOAT_OPERATOR_SET) return output + "=" + "SET_" + input1;
	if (type == SCHEDULE_FLOAT_OPERATOR_PLUS) return output + "=" + input1 + "_PLUS_" + input2;
	if (type == SCHEDULE_FLOAT_OPERATOR_MINUS) return output + "=" + input1 + "_MINUS_" + input2;
	if (type == SCHEDULE_FLOAT_OPERATOR_MULT) return output + "=" + input1 + "_MULT_" + input2;
	if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE) return output + "=" + input1 + "_DIV_" + input2;
	if (type == SCHEDULE_FLOAT_OPERATOR_ROUND) return output + "=" + "ROUND_" + input1;
	if (type == SCHEDULE_FLOAT_OPERATOR_COUNT_IMAGES) return output + "=" + "COUNT_IMGS_" + input1 + "_" + input2;
	if (type == SCHEDULE_FLOAT_OPERATOR_COUNT_WORDS) return output + "=" + "COUNT_WORDS_" + input1;
	if (type == SCHEDULE_FLOAT_OPERATOR_READ_STAR) return output + "=" + "STAR_" + input1 + "_" + input2;
	if (type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MAX) return output + "=" + "STAR_MAX_" + input1;
	if (type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MIN) return output + "=" + "STAR_MIN_" + input1;
	if (type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_AVG) return output + "=" + "STAR_AVG_" + input1;
	if (type == SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX) return output + "=" + "STAR_SORT_IDX_" + input1 + "_" + input2;
	if (type == SCHEDULE_STRING_OPERATOR_JOIN) return output + "=" + "JOIN_" + input1 + "_" + input2;
	if (type == SCHEDULE_STRING_OPERATOR_BEFORE_FIRST) return output + "=" + "BEF_FIRST_" + input1 + "_" + input2;
	if (type == SCHEDULE_STRING_OPERATOR_AFTER_FIRST) return output + "=" + "AFT_FIRST_" + input1 + "_" + input2;
	if (type == SCHEDULE_STRING_OPERATOR_BEFORE_LAST) return output + "=" + "BEF_LAST_" + input1 + "_" + input2;
	if (type == SCHEDULE_STRING_OPERATOR_AFTER_LAST) return output + "=" + "AFT_LAST_" + input1 + "_" + input2;
	if (type == SCHEDULE_STRING_OPERATOR_READ_STAR) return output + "=" + "STAR_" + input1 + "_" + input2;
	if (type == SCHEDULE_STRING_OPERATOR_GLOB) return output + "=" + "GLOB_" + input1;
	if (type == SCHEDULE_STRING_OPERATOR_NTH_WORD) return output + "=" + "NTH_WORD_" + input1 + "_" + input2;
	if (type == SCHEDULE_OPERATOR_TOUCH_FILE) return "TOUCH_" + input1;
	if (type == SCHEDULE_OPERATOR_COPY_FILE) return "COPY_" + input1 + "_TO_" + input2;
	if (type == SCHEDULE_OPERATOR_MOVE_FILE) return "MOVE_" + input1 + "_TO_" + input2;
	if (type == SCHEDULE_OPERATOR_DELETE_FILE) return "DELETE_" + input1;
	if (type == SCHEDULE_WAIT_OPERATOR_SINCE_LAST_TIME) return "WAIT_" + input1;
	if (type == SCHEDULE_EMAIL_OPERATOR) return "EMAIL_" + input1 + "_" + input2;
	if (type == SCHEDULE_EXIT_OPERATOR) return "EXIT";
	else
		REPORT_ERROR("ERROR: unrecognised Operator type:" + type);
}

std::string SchedulerEdge::getOutputNode() const
{
	if (is_fork) return (scheduler_global_bools[myBooleanVariable].value) ?  outputNodeTrue : outputNode;
	else return outputNode;
}

void Schedule::clear()
{
	verb = 1;
	current_node = "undefined";
	name = "undefined";
	do_read_only = false;
	scheduler_global_bools.clear();
	scheduler_global_floats.clear();
	scheduler_global_strings.clear();
	jobs.clear();
	edges.clear();
	scheduler_global_operators.clear();
}

std::string Schedule::findJobByCurrentName(std::string _name)
{
	std::map<std::string, SchedulerJob>::iterator it;

	for ( it = jobs.begin(); it != jobs.end(); it++ )
	{
	    if (it->second.current_name == _name)
	    	return it->second.current_name;
	}
	REPORT_ERROR("ERROR: cannot find job: " + _name);
	return "";
}

void Schedule::read(bool do_lock, FileName fn)
{
#ifdef DEBUG_LOCK
	std::cerr << "entering read lock_message=" << lock_message << std::endl;
#endif
	FileName name_wo_dir = name;
	name_wo_dir = name_wo_dir.beforeLastOf("/");
	FileName dir_lock=".relion_lock_schedule_" + name_wo_dir.afterLastOf("/"), fn_lock=dir_lock + "/lock_schedule";;
	if (do_lock && !do_read_only)
	{
		int iwait =0;
		int status = mkdir(dir_lock.c_str(), S_IRWXU);

#ifdef DEBUG_LOCK
		std::cerr <<  " A status= " << status << std::endl;
#endif
		while (!status == 0)
		{
			if (errno == EACCES) // interestingly, not EACCESS!
				REPORT_ERROR("ERROR: Schedule::read cannot create a lock directory " + dir_lock + ". You don't have write permission to this project. If you want to look at other's project directory (but run nothing there), please start RELION with --readonly.");

			// If the lock exists: wait 3 seconds and try again
			// First time round, print a warning message
			if (iwait == 0)
			{
				std::cout << " WARNING: trying to read schedule.star, but directory " << dir_lock << " exists (which protects against simultaneous writing)" << std::endl;
			}
			sleep(3);
			status =  mkdir(dir_lock.c_str(), S_IRWXU);
#ifdef DEBUG_LOCK
			std::cerr <<  " B status= " << status << std::endl;
#endif

			iwait++;
			if (iwait > 40)
			{

				REPORT_ERROR("ERROR: Schedule::read has waited for 2 minutes for lock directory to disappear. Make sure this scheduler is not running, and then manually remove the file: " + fn_lock);
			}

		}
		// Generate the lock file
		std::ofstream  fh;
		fh.open(fn_lock.c_str(), std::ios::out);
		if (!fh)
			REPORT_ERROR( (std::string)"ERROR: Cannot open file: " + fn_lock);
		std::string lock_message = "lock mechanism from Scheduler";
		fh << lock_message << std::endl;
		fh.close();
	}

	if (fn == "") fn = name + "schedule.star";

	// Clear current model
	clear();

	// Open input file
	std::ifstream in(fn.data(), std::ios_base::in);
	if (in.fail())
		REPORT_ERROR( (std::string) "Schedule::read: File " + fn + " cannot be read." );

	// For reading: do the nodes before the general table, in order to set current_node
	MetaDataTable MD;
	MD.readStar(in, "schedule_general");
	MD.getValue(EMDL_SCHEDULE_GENERAL_NAME, name);
	MD.getValue(EMDL_SCHEDULE_GENERAL_CURRENT_NODE, current_node);
	MD.clear();

	MD.readStar(in, "schedule_floats");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string myname;
		RFLOAT value, original_value;
		MD.getValue(EMDL_SCHEDULE_VAR_FLOAT_NAME, myname);
		MD.getValue(EMDL_SCHEDULE_VAR_FLOAT_VALUE, value);
		MD.getValue(EMDL_SCHEDULE_VAR_FLOAT_ORI_VALUE, original_value);
		SchedulerFloatVariable myval(value, original_value);
		scheduler_global_floats[myname] = myval;
	}
	MD.clear();

	MD.readStar(in, "schedule_bools");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string myname;
		bool value, original_value;
		MD.getValue(EMDL_SCHEDULE_VAR_BOOL_NAME, myname);
		MD.getValue(EMDL_SCHEDULE_VAR_BOOL_VALUE, value);
		MD.getValue(EMDL_SCHEDULE_VAR_BOOL_ORI_VALUE, original_value);
		SchedulerBooleanVariable myval(value, original_value);
		scheduler_global_bools[myname] = myval;
	}
	MD.clear();

	MD.readStar(in, "schedule_strings");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string myname;
		FileName value, original_value;
		MD.getValue(EMDL_SCHEDULE_VAR_STRING_NAME, myname);
		MD.getValue(EMDL_SCHEDULE_VAR_STRING_VALUE, value);
		MD.getValue(EMDL_SCHEDULE_VAR_STRING_ORI_VALUE, original_value);
		SchedulerStringVariable myval(value, original_value);
		scheduler_global_strings[myname] = myval;
	}
	MD.clear();

	MD.readStar(in, "schedule_operators");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string myname, type, input1, input2, output;
		RFLOAT constant;

		MD.getValue(EMDL_SCHEDULE_OPERATOR_NAME, myname);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_TYPE, type);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_INPUT1, input1);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_INPUT2, input2);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_OUTPUT, output);
		SchedulerOperator myval(type, input1, input2, output);
		scheduler_global_operators[myname] = myval;
	}
	MD.clear();

	MD.readStar(in, "schedule_jobs");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string myname, ori_name, mode, type;
		bool has_started;

		MD.getValue(EMDL_SCHEDULE_JOB_NAME, myname);
		MD.getValue(EMDL_SCHEDULE_JOB_ORI_NAME, ori_name);
		MD.getValue(EMDL_SCHEDULE_JOB_MODE, mode);
		MD.getValue(EMDL_SCHEDULE_JOB_HAS_STARTED, has_started);

		SchedulerJob myval(myname, mode, has_started);
		jobs[ori_name] = myval;
	}
	MD.clear();


	MD.readStar(in, "schedule_edges");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		int number;
		std::string inputname, outputname, outputname_true, bool_name;
		bool is_fork;

		MD.getValue(EMDL_SCHEDULE_EDGE_INPUT, inputname);
		MD.getValue(EMDL_SCHEDULE_EDGE_OUTPUT, outputname);
		MD.getValue(EMDL_SCHEDULE_EDGE_IS_FORK, is_fork);
		MD.getValue(EMDL_SCHEDULE_EDGE_OUTPUT_TRUE, outputname_true);
		MD.getValue(EMDL_SCHEDULE_EDGE_BOOLEAN, bool_name);
		SchedulerEdge myval(inputname, outputname, is_fork, bool_name, outputname_true);
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

	// Also read in the schedule_pipeline (no need to lock now?)
	schedule_pipeline.read();

}

bool Schedule::isWriteLocked()
{
	FileName name_wo_dir = name;
	name_wo_dir = name_wo_dir.beforeLastOf("/");
	FileName dir_lock=".relion_lock_schedule_" + name_wo_dir.afterLastOf("/"), fn_lock=dir_lock + "/lock_schedule";;
	return exists(dir_lock);
}

void Schedule::write(bool do_lock, FileName fn)
{
	if (do_read_only)
		return;

	FileName name_wo_dir = name;
	name_wo_dir = name_wo_dir.beforeLastOf("/");
	FileName dir_lock=".relion_lock_schedule_" + name_wo_dir.afterLastOf("/"), fn_lock=dir_lock + "/lock_schedule";;
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

	if (fn == "") fn = name + "schedule.star";

	// B. Write STAR file with the entire schedule
	std::ofstream  fh;
	fh.open((fn).c_str(), std::ios::out);
	if (!fh)
		REPORT_ERROR( (std::string)"Schedule::write: Cannot write file: " + fn);

	MetaDataTable MDgeneral;
	MDgeneral.setName("schedule_general");
	MDgeneral.setIsList(true);
	MDgeneral.addObject();
	MDgeneral.setValue(EMDL_SCHEDULE_GENERAL_NAME, name);
	MDgeneral.setValue(EMDL_SCHEDULE_GENERAL_CURRENT_NODE, current_node);
	MDgeneral.write(fh);

	if (scheduler_global_floats.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_floats");
		std::map<std::string, SchedulerFloatVariable>::iterator it;
		for ( it = scheduler_global_floats.begin(); it != scheduler_global_floats.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_VAR_FLOAT_NAME, it->first);
			MD.setValue(EMDL_SCHEDULE_VAR_FLOAT_VALUE, it->second.value);
			MD.setValue(EMDL_SCHEDULE_VAR_FLOAT_ORI_VALUE, it->second.original_value);
		}
		MD.write(fh);
	}

	if (scheduler_global_bools.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_bools");
		std::map<std::string, SchedulerBooleanVariable>::iterator it;
		for ( it = scheduler_global_bools.begin(); it != scheduler_global_bools.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_VAR_BOOL_NAME, it->first);
			MD.setValue(EMDL_SCHEDULE_VAR_BOOL_VALUE, it->second.value);
			MD.setValue(EMDL_SCHEDULE_VAR_BOOL_ORI_VALUE, it->second.original_value);
		}
		MD.write(fh);
	}

	if (scheduler_global_strings.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_strings");
		std::map<std::string, SchedulerStringVariable>::iterator it;
		for ( it = scheduler_global_strings.begin(); it != scheduler_global_strings.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_VAR_STRING_NAME, it->first);
			MD.setValue(EMDL_SCHEDULE_VAR_STRING_VALUE, it->second.value);
			MD.setValue(EMDL_SCHEDULE_VAR_STRING_ORI_VALUE, it->second.original_value);
		}
		MD.write(fh);
	}

	if (scheduler_global_operators.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_operators");
		std::map<std::string, SchedulerOperator>::iterator it;
		for ( it = scheduler_global_operators.begin(); it != scheduler_global_operators.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_OPERATOR_NAME, it->first);
			MD.setValue(EMDL_SCHEDULE_OPERATOR_TYPE, it->second.type);
			MD.setValue(EMDL_SCHEDULE_OPERATOR_OUTPUT, it->second.output );
			MD.setValue(EMDL_SCHEDULE_OPERATOR_INPUT1, it->second.input1 );
			MD.setValue(EMDL_SCHEDULE_OPERATOR_INPUT2, it->second.input2 );
		}
		MD.write(fh);
	}

	if (jobs.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_jobs");
		std::map<std::string, SchedulerJob>::iterator it;
		for ( it = jobs.begin(); it != jobs.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_JOB_ORI_NAME, it->first);
			MD.setValue(EMDL_SCHEDULE_JOB_NAME, it->second.current_name);
			MD.setValue(EMDL_SCHEDULE_JOB_MODE, it->second.mode);
			MD.setValue(EMDL_SCHEDULE_JOB_HAS_STARTED, it->second.job_has_started);
		}
		MD.write(fh);
	}

	if (edges.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_edges");
		for (int i = 0; i < edges.size(); i++)
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_EDGE_INPUT, edges[i].inputNode);
			MD.setValue(EMDL_SCHEDULE_EDGE_OUTPUT, edges[i].outputNode);
			MD.setValue(EMDL_SCHEDULE_EDGE_IS_FORK, edges[i].is_fork);
			MD.setValue(EMDL_SCHEDULE_EDGE_OUTPUT_TRUE, edges[i].outputNodeTrue);
			MD.setValue(EMDL_SCHEDULE_EDGE_BOOLEAN, edges[i].myBooleanVariable);
		}
		MD.write(fh);
	}

	// Close the file handler
	fh.close();

	// Also write out the schedule_pipeline (no need to lock now?)
	schedule_pipeline.write();

	if (do_lock)
	{

#ifdef DEBUG_LOCK
		std::cerr << " write schedule: now deleting " << fn_lock << std::endl;
#endif

		if (!exists(fn_lock))
			REPORT_ERROR("ERROR: Schedule::write was expecting a file called "+fn_lock+ " but it is no longer there.");
		unlock();
	}

	// Touch a file to indicate to the GUI that the pipeline has just changed
	FileName mychanged = name+SCHEDULE_HAS_CHANGED;
	touch(mychanged);
}

// Reset all variables to their original value
void Schedule::reset()
{
	{
		std::map<std::string, SchedulerFloatVariable>::iterator it;
		for ( it = scheduler_global_floats.begin(); it != scheduler_global_floats.end(); it++ )
			it->second.value = it->second.original_value;
	}

	{
		std::map<std::string, SchedulerBooleanVariable>::iterator it;
		for ( it = scheduler_global_bools.begin(); it != scheduler_global_bools.end(); it++ )
	        it->second.value = it->second.original_value;
	}

	{
		std::map<std::string, SchedulerStringVariable>::iterator it;
		for ( it = scheduler_global_strings.begin(); it != scheduler_global_strings.end(); it++ )
	        it->second.value = it->second.original_value;
	}

	{
		std::map<std::string, SchedulerJob>::iterator it;
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

bool Schedule::isNode(std::string _name)
{
	// is this either an operator or a job?
	return (jobs.find(_name) != jobs.end() || scheduler_global_operators.find(_name) != scheduler_global_operators.end());
}

bool Schedule::isJob(std::string _name)
{
	return (jobs.find(_name) != jobs.end());
}

bool Schedule::isOperator(std::string _name)
{
	return isScheduleOperator(_name);
}

void Schedule::setVariable(std::string name, FileName value)
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

void Schedule::setOriginalVariable(std::string name, FileName value)
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

void Schedule::addFloatVariable(std::string name, RFLOAT value)
{
	if (isFloatVariable(name))
		REPORT_ERROR("ERROR: trying to add a float variable with a name that already exists: " + name);

	SchedulerFloatVariable myvar(value, value);
	scheduler_global_floats[name] = myvar;
}

void Schedule::addBooleanVariable(std::string name, bool value)
{
	if (isBooleanVariable(name))
		REPORT_ERROR("ERROR: trying to add a boolean variable with a name that already exists: " + name);

	SchedulerBooleanVariable myvar(value, value);
	scheduler_global_bools[name] = myvar;
}

void Schedule::addStringVariable(std::string name, FileName value)
{
	if (isStringVariable(name))
		REPORT_ERROR("ERROR: trying to add a string variable with a name that already exists: " + name);

	SchedulerStringVariable myvar(value, value);
	scheduler_global_strings[name] = myvar;
}

float Schedule::getFloatVariableValue(std::string name)
{
	if (!isFloatVariable(name))
		REPORT_ERROR("ERROR: cannot find float variable with name:" + name);
	return scheduler_global_floats[name].value;
}

float Schedule::getFloatOriginalVariableValue(std::string name)
{
	if (!isFloatVariable(name))
		REPORT_ERROR("ERROR: cannot find float variable with name:" + name);
	return scheduler_global_floats[name].original_value;
}

void Schedule::setFloatVariableValue(std::string name, RFLOAT val)
{
	if (!isFloatVariable(name))
		REPORT_ERROR("ERROR: cannot find float variable with name:" + name);
	scheduler_global_floats[name].value = val;
	scheduler_global_floats[name].value = val;
}

void Schedule::setFloatOriginalVariableValue(std::string name, RFLOAT val)
{
	if (!isFloatVariable(name))
		REPORT_ERROR("ERROR: cannot find float variable with name:" + name);
	scheduler_global_floats[name].original_value = val;
}

bool Schedule::getBooleanVariableValue(std::string name)
{
	if (!isBooleanVariable(name))
		REPORT_ERROR("ERROR: cannot find boolean variable with name:" + name);
	return scheduler_global_bools[name].value;
}

bool Schedule::getBooleanOriginalVariableValue(std::string name)
{
	if (!isBooleanVariable(name))
		REPORT_ERROR("ERROR: cannot find boolean variable with name:" + name);
	return scheduler_global_bools[name].original_value;
}

void Schedule::setBooleanVariableValue(std::string name, bool val)
{
	if (!isBooleanVariable(name))
		REPORT_ERROR("ERROR: cannot find boolean variable with name:" + name);
	scheduler_global_bools[name].value = val;
}

void Schedule::setBooleanOriginalVariableValue(std::string name, bool val)
{
	if (!isBooleanVariable(name))
		REPORT_ERROR("ERROR: cannot find boolean variable with name:" + name);
	scheduler_global_bools[name].original_value = val;
}

std::string Schedule::getStringVariableValue(std::string name)
{
	if (!isStringVariable(name))
		REPORT_ERROR("ERROR: cannot find string variable with name:" + name);
	return scheduler_global_strings[name].value;
}

std::string Schedule::getStringOriginalVariableValue(std::string name)
{
	if (!isStringVariable(name))
		REPORT_ERROR("ERROR: cannot find string variable with name:" + name);
	return scheduler_global_strings[name].original_value;
}

void Schedule::setStringVariableValue(std::string name, std::string val)
{
	if (!isStringVariable(name))
		REPORT_ERROR("ERROR: cannot find string variable with name:" + name);
	scheduler_global_strings[name].value = val;
}

void Schedule::setStringOriginalVariableValue(std::string name, std::string val)
{
	if (!isStringVariable(name))
		REPORT_ERROR("ERROR: cannot find string variable with name:" + name);
	scheduler_global_strings[name].original_value = val;
}

std::string Schedule::getVariableValueAsString(std::string name)
{
	if (isStringVariable(name)) return scheduler_global_strings[name].value;
	else if (isBooleanVariable(name)) return (scheduler_global_bools[name].value) ? "True" : "False";
	else if (isFloatVariable(name)) return floatToString(scheduler_global_floats[name].value);
	else REPORT_ERROR("Schedule::getVariableValueAsString: no variable named " + name);
}

void Schedule::setOperatorParameters(std::string name, std::string _type, std::string _input1, std::string _input2, std::string _output)
{
	if (!isOperator(name))
		REPORT_ERROR("ERROR: cannot find operator with name:" + name);
	// Just make a new one, so all the check are done automatically...
	SchedulerOperator myop(_type, _input1, _input2, _output);
	scheduler_global_operators[name] = myop;
}

void Schedule::getOperatorParameters(std::string name, std::string &_type, std::string &_input1, std::string &_input2, std::string &_output)
{
	if (!isOperator(name))
		REPORT_ERROR("ERROR: cannot find operator with name:" + name);
	_type = scheduler_global_operators[name].type;
	_input1 = scheduler_global_operators[name].input1;
	_input2 = scheduler_global_operators[name].input2;
	_output = scheduler_global_operators[name].output;
}

std::map<std::string, SchedulerFloatVariable> Schedule::getCurrentFloatVariables()
{
	return scheduler_global_floats;
}

std::map<std::string, SchedulerBooleanVariable> Schedule::getCurrentBooleanVariables()
{
	return scheduler_global_bools;
}

std::map<std::string, SchedulerStringVariable> Schedule::getCurrentStringVariables()
{
	return scheduler_global_strings;
}

std::map<std::string, SchedulerOperator> Schedule::getCurrentOperators()
{
	return scheduler_global_operators;
}

SchedulerOperator Schedule::initialiseOperator(std::string type, std::string input_name, std::string input2_name,
		std::string output_name, std::string &error_message)
{
	SchedulerOperator myop;
	error_message = myop.initialise(type, input_name, input2_name, output_name);
	return myop;
}

void Schedule::addOperator(SchedulerOperator &myop)
{
	std::string myname = myop.getName();
	scheduler_global_operators[myname] = myop;
}

void Schedule::addJob(RelionJob &myjob, std::string jobname, std::string mode)
{

	//remove spaces from jobname
	for (int i = 0; i < jobname.length(); i++)
	{
		if (jobname[i] == ' ') jobname[i] = '_';
	}

	// Check whether the jobname is unique
	if (isNode(jobname))
		REPORT_ERROR("ERROR: trying to add a JobNode that already exists: " + jobname);
	// Now add this job to the local schedule_pipeline
	std::string error_message;
	std::vector<std::string> commands;
	std::string final_command;
	std::string output_name = name + jobname + '/';

	// Save a copy of the job in the Schedules directory
	std::string command = "mkdir -p " + name + jobname;
	int res = system(command.c_str());

	myjob.write(output_name);

	if (!myjob.getCommands(output_name, commands, final_command, false, schedule_pipeline.job_counter, error_message))
		REPORT_ERROR("ERROR in getting commands for scheduled job: " + error_message);

	int current_job = schedule_pipeline.addJob(myjob, PROC_SCHEDULED, false, false); // 1st false is do_overwrite, 2nd false is do_write_minipipeline

	if (current_job < 0)
		REPORT_ERROR("ERROR: current job should not be negative now ...");

	SchedulerJob mynode(jobname, mode, false);
	jobs[jobname] = mynode;
}

void Schedule::removeVariable(std::string name)
{
	// Remove any operators with this variable in it
	removeOperatorsWithThisInputOrOutput(name);

	if (isBooleanVariable(name))
	{
		scheduler_global_bools.erase(name);
		// Also remove forks with this boolean variable
		removeEdgesWithThisInputOutputOrBoolean(name);
	}
	else if (isFloatVariable(name)) scheduler_global_floats.erase(name);
	else if (isStringVariable(name)) scheduler_global_strings.erase(name);
	else REPORT_ERROR("ERROR: cannot find variable to erase: " + name);
}

void Schedule::removeEdgesWithThisInputOutputOrBoolean(std::string name)
{
	std::vector<SchedulerEdge> new_edges;
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

void Schedule::removeOperator(std::string name)
{
	if (isOperator(name)) scheduler_global_operators.erase(name);
	else REPORT_ERROR("ERROR: cannot find operator to erase: " + name);

	// Also remove any edges that input/output with this operator
	removeEdgesWithThisInputOutputOrBoolean(name);
}

void Schedule::removeOperatorsWithThisInputOrOutput(std::string name)
{
	std::map<std::string, SchedulerOperator> new_operators;

	std::map<std::string, SchedulerOperator>::iterator it;
	for ( it = scheduler_global_operators.begin(); it != scheduler_global_operators.end(); it++ )
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

	scheduler_global_operators = new_operators;
}

void Schedule::removeJob(std::string name)
{
	if (isJob(name)) jobs.erase(name);
	else REPORT_ERROR("ERROR: cannot find job to erase: " + name);

	// Also remove any edges that input/output with this job
	removeEdgesWithThisInputOutputOrBoolean(name);
}

void Schedule::removeEdge(int idx)
{
	edges.erase(edges.begin()+idx);
}

void Schedule::copy(FileName newname)
{
	// Make sure newname ends with a slash
	if (newname[newname.length()-1] != '/') newname += "/";

	// Make the output directory,
	std::string command = "mkdir -p " + newname;
	int res = system(command.c_str());

	// Replace original name in all stringVariables
	for (std::map<std::string, SchedulerStringVariable>::iterator it = scheduler_global_strings.begin(); it != scheduler_global_strings.end(); it++ )
	{
		(it->second.value).replaceAllSubstrings(name, newname);
		(it->second.original_value).replaceAllSubstrings(name, newname);
	}

	// Also replace all names of Nodes and Processes in the Pipeliner
	for (int i = 0; i < schedule_pipeline.nodeList.size(); i++)
	{
		FileName myname = schedule_pipeline.nodeList[i].name;
		if (myname.contains(name)) myname.replaceAllSubstrings(name, newname);
		schedule_pipeline.nodeList[i].name = myname;
	}

	for (int i = 0; i < schedule_pipeline.processList.size(); i++)
	{
		FileName myname = schedule_pipeline.processList[i].name;
		if (myname.contains(name)) myname.replaceAllSubstrings(name, newname);
		schedule_pipeline.processList[i].name = myname;
	}

	// Replace all names in the pipeliner jobs
	for (std::map<std::string, SchedulerJob>::iterator it = jobs.begin(); it != jobs.end(); it++ )
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
		std::string command = "mkdir -p " + mydir;
		int res = system(command.c_str());
		myjob.write(mydir);
	}

	// Change the name itself
	setName(newname);

	// And write the new schedule.star and scheule_pipeline.star files
	write();
}

void schedulerSendEmail(std::string message, std::string subject)
{
	if (isStringVariable("email"))
	{
		std::string command = "echo \"" + message + "\" | mail -s \"" + subject + "\" -r RELION " + scheduler_global_strings["email"].value;
		int res = system(command.c_str());
	}
}

void Schedule::addEdge(std::string inputnode_name, std::string outputnode_name)
{
	SchedulerEdge myval(inputnode_name, outputnode_name);
	edges.push_back(myval);
}

void Schedule::addFork(std::string inputnode_name, std::string mybool_name, std::string outputnode_name_if_true, std::string outputnode_name)
{
	SchedulerEdge myval(inputnode_name, outputnode_name, true, mybool_name, outputnode_name_if_true);
	edges.push_back(myval);
}

bool Schedule::isValid()
{
	// TODO: check if duplicate edges, forks or scheduler_global_operators exist....

	// Check Scheduler ends with an exit

	return false; // to be implemented
}

std::string Schedule::getNextNode()
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

std::string Schedule::getPreviousNode()
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

bool Schedule::gotoNextNode()
{
	current_node = getNextNode();
	if (verb > 1) std::cout << " ++ Setting current node to: " << current_node << std::endl;

	// Write out current status, but maintain lock on the directory!
	write();

	return (current_node == "undefined") ? false : true;
}

bool Schedule::gotoNextJob()
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
			SchedulerOperator my_op = scheduler_global_operators[current_node];

			// Now change any original job names in input/output for their corresponding current_names
			std::map<std::string, SchedulerJob>::iterator it;
			for ( it = jobs.begin(); it != jobs.end(); it++ )
			{
				if (isStringVariable(my_op.input1) && scheduler_global_strings[my_op.input1].value.contains(it->first))
				{
					// Make a new temporary stringVariable, act on it with my_op and then delete again?
					FileName newval = scheduler_global_strings[my_op.input1].value;
					newval.replaceAllSubstrings(name+it->first+"/", it->second.current_name);
					addStringVariable("xxx_tmp_input1", newval);
					my_op.input1 = "xxx_tmp_input1";
				}
				if (isStringVariable(my_op.input2) && scheduler_global_strings[my_op.input2].value.contains(it->first))
				{
					FileName newval = scheduler_global_strings[my_op.input2].value;
					newval.replaceAllSubstrings(name+it->first+"/", it->second.current_name);
					addStringVariable("xxx_tmp_input2", newval);
					my_op.input2 = "xxx_tmp_input2";
				}
				if (isStringVariable(my_op.output) && scheduler_global_strings[my_op.output].value.contains(it->first))
				{
					FileName newval = scheduler_global_strings[my_op.output].value;
					newval.replaceAllSubstrings(name+it->first+"/", it->second.current_name);
					addStringVariable("xxx_tmp_output", newval);
					my_op.output = "xxx_tmp_input1";
				}
			}

			bool op_success = my_op.performOperation();
			if (isStringVariable("xxx_tmp_input1")) removeVariable("xxx_tmp_input1");
			if (isStringVariable("xxx_tmp_input2")) removeVariable("xxx_tmp_input2");
			if (isStringVariable("xxx_tmp_output"))
			{
				// Set output in the variable from the original operator, and then remove tmp_output
				std::string myname = scheduler_global_operators[current_node].output;
				scheduler_global_strings[myname].value = scheduler_global_strings["xxx_tmp_output"].value;
				removeVariable("xxx_tmp_output");
			}

			if (verb > 0 && op_success)
			{
				if (scheduler_global_operators[current_node].output != "undefined")
				{
					std::cout << " + " << current_node << " => " <<getVariableValueAsString(scheduler_global_operators[current_node].output) << std::endl;
				}
			}

			if (!op_success) return false;
		}
		else // this is a job, get its current_name and options
		{
			return true;
		}
	}

	return false;
}

// Modify a job to set variables from the Scheduler
void Schedule::setVariablesInJob(RelionJob &job, FileName original_job_name, bool &needs_a_restart)
{
	needs_a_restart = false;

	RelionJob ori_job;
	bool dummy;
	ori_job.read(name + original_job_name + '/', dummy, true);

	// Check where this job gets its input from: change names from local scheduler ones to the current pipeline
	int ori_process = schedule_pipeline.findProcessByName(name + original_job_name + '/');
	// Loop over all input nodes to this job
	for (int inode = 0; inode <  schedule_pipeline.processList[ori_process].inputNodeList.size(); inode++)
	{
		int mynode = schedule_pipeline.processList[ori_process].inputNodeList[inode];
		// find from which pipeline_scheduler job this jobs gets its input nodes
		int output_from_process =  schedule_pipeline.nodeList[mynode].outputFromProcess;

		if (output_from_process < 0)
		{
			// This was not a process, just continue
			break;
		}

		// Get the original current_name in the pipeline_scheduler of this job
		FileName my_ori_name = schedule_pipeline.processList[output_from_process].name;
		// Remove leading directory and tailing slash to get the process current_name in the pipeline_scheduler
		FileName my_process_name = (my_ori_name.afterFirstOf(name)).beforeLastOf("/");
		// find that process in the nodes, and get its current current_name
		std::string my_current_name =jobs[my_process_name].current_name;

		// Change all instances of the my_ori_name to my_current_name, take from ori_job and set into this job
		for (std::map<std::string,JobOption>::iterator it=ori_job.joboptions.begin(); it!=ori_job.joboptions.end(); ++it)
		{
			FileName mystring = (it->second).value;
			if (mystring.contains(my_ori_name))
			{
				mystring.replaceAllSubstrings(my_ori_name, my_current_name);
				FileName myval = job.joboptions[it->first].value;
				myval = myval.beforeLastOf("/") + "/";
				// If any of the input nodes are not the same as my_current_name (from continuation jobs) or my_ori_name (from new jobs)
				if (myval != my_current_name && myval != my_ori_name)
				{
					std::cerr << "restart needed!" << std::endl;
					needs_a_restart = true;
				}
				job.joboptions[it->first].value = mystring;
			}
		}
	}

	// Check whether there are any options with a value containing $$, which is the sign for inserting Scheduler variables
	for (std::map<std::string,JobOption>::iterator it=ori_job.joboptions.begin(); it!=ori_job.joboptions.end(); ++it)
	{
		FileName mystring = (it->second).value;
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
				REPORT_ERROR(" ERROR: cannot find anything after $$ sign in string: " + mystring);
			std::string mypat = splits[0];
			myvars.push_back(mypat);

			// Found an option that needs replacement! Now find which variable to insert
			std::string my_value;
			if (isBooleanVariable(mypat))
			{
				if (job.joboptions[it->first].joboption_type != JOBOPTION_BOOLEAN)
					REPORT_ERROR(" ERROR: trying to set a BooleanVariable: " + mypat + " into a non-boolean option: " + it->first);

				my_value = (scheduler_global_bools[mypat].value) ? "Yes" : "No";
			}
			else if (isFloatVariable(mypat))
			{
				if (job.joboptions[it->first].joboption_type == JOBOPTION_BOOLEAN)
					REPORT_ERROR(" ERROR: trying to set FloatVariable: " + mypat + " into a boolean option: " + it->first);

				my_value = floatToString(scheduler_global_floats[mypat].value);
			}
			else if (isStringVariable(splits[0]))
			{
				if (job.joboptions[it->first].joboption_type == JOBOPTION_BOOLEAN)
					REPORT_ERROR(" ERROR: trying to set StringVariable: " + mypat + " into a boolean option: " + it->first);
				if (job.joboptions[it->first].joboption_type == JOBOPTION_SLIDER)
					REPORT_ERROR(" ERROR: trying to set StringVariable: " + mypat + " into a slider option: " + it->first);

				my_value = scheduler_global_strings[mypat].value;
			}
			else
				REPORT_ERROR(" ERROR: variable in job is not part of this Schedule: " + mypat);

			mystring = before + my_value;
			for (int i = 1; i < splits.size(); i++)
				mystring += " " + splits[i];
		}

		if (has_found)
		{
			job.joboptions[it->first].value = mystring;
			std::string myvarsstr = "";
			for (int i = 0; i < myvars.size(); i++)
				myvarsstr+= myvars[i]+ " ";
			if (verb > 2) std::cout << " +++ Setting joboption " << it->first << " to " << mystring << " based on variable(s): " << myvarsstr<< std::endl;
		}
	}
}

void Schedule::run(PipeLine &pipeline)
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
    	std::cout << " + Starting execution of schedule at: " << current_node << std::endl;
    }

    // If we start from operators instead of jobs, then execute the operator and proceed until the next Job
	bool has_more_jobs = true;
    if (isOperator(current_node))
    {
		bool op_success = scheduler_global_operators[current_node].performOperation();
		if (op_success)
		{
			if (scheduler_global_operators[current_node].output != "undefined" && verb > 0)
				std::cout << " + " << current_node << " => " <<getVariableValueAsString(scheduler_global_operators[current_node].output) << std::endl;
		}
		else
		{
			REPORT_ERROR("ERROR: something went wrong with execution of the initial operator...");
		}

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

		RelionJob myjob;
		bool is_continue, do_overwrite_current, dummy;
		int current_job;
		if (!jobs[current_node].job_has_started || jobs[current_node].mode == SCHEDULE_NODE_JOB_MODE_NEW)
		{
			// Read the job from inside the Schedule schedule_pipeline to a new job
			bool dummy;
			myjob.read(name + current_node + '/', dummy, true); // true means initialise the job

			// This function replaces variable calls starting with a '$$' from the original_job into the current_job
			// It will also take care of dealing with inheritance of the correct inputNode names
			setVariablesInJob(myjob, current_node, dummy);

			// Now add this job to the pipeline we will actually be running in
			current_job = pipeline.addScheduledJob(myjob);
			is_continue = false;
			do_overwrite_current = false;

			// Set the current_name of the current node now
			jobs[current_node].current_name = pipeline.processList[current_job].name;
				if (verb > 0) std::cout << " + Creating new Job: " << jobs[current_node].current_name << " from Node: " << current_node << std::endl;
		}
		else if (jobs[current_node].mode == SCHEDULE_NODE_JOB_MODE_CONTINUE || jobs[current_node].mode == SCHEDULE_NODE_JOB_MODE_OVERWRITE)
		{
			is_continue = (jobs[current_node].mode == SCHEDULE_NODE_JOB_MODE_CONTINUE);
			do_overwrite_current = (jobs[current_node].mode == SCHEDULE_NODE_JOB_MODE_OVERWRITE);

			current_job = pipeline.findProcessByName(jobs[current_node].current_name);
			if (current_job < 0)
				REPORT_ERROR("ERROR: RunSchedule cannot find process with name: " + jobs[current_node].current_name);

			// Read the job from the pipeline we are running in
			if (!myjob.read(pipeline.processList[current_job].name, dummy, true)) // true means also initialise the job
				REPORT_ERROR("There was an error reading job: " + pipeline.processList[current_job].name);

			// This function replaces variable calls starting with a '$$' from the original_job into the current_job
			// It will also take care of dealing with inheritance of the correct inputNode names
			bool needs_a_restart = false;
			setVariablesInJob(myjob, current_node, needs_a_restart);

			if (needs_a_restart)
			{
				// Now add this job to the pipeline we will actually be running in
				current_job = pipeline.addScheduledJob(myjob);
				is_continue = false;
				do_overwrite_current = false;

				// Set the current_name of the current node now
				jobs[current_node].current_name = pipeline.processList[current_job].name;
				if (verb > 0) std::cout << " + Creating new Job: " << jobs[current_node].current_name << " from node " << current_node << std::endl;
			}
		}
		else
			REPORT_ERROR("ERROR: unrecognised mode for running a new process: " + jobs[current_node].mode);

		// Check whether the input nodes are there, before executing the job
		for (long int inode = 0; inode < pipeline.processList[current_job].inputNodeList.size(); inode++)
		{
			long int mynode = pipeline.processList[current_job].inputNodeList[inode];
			while (!exists(pipeline.nodeList[mynode].name))
			{
				std::cerr << " + -- Warning " << pipeline.nodeList[mynode].name << " does not exist. Waiting 10 seconds ... " << std::endl;
				sleep(10);

				// Abort mechanism
				if (pipeline_control_check_abort_job())
				{
					write(DO_LOCK);
					exit(RELION_EXIT_ABORTED);
				}

			}
		}

		// Now actually run the Scheduled job
		std::string error_message;
		if (verb > 0)
		{
			time_t my_time = time(NULL);
			std::cout << " + Executing Job: " << jobs[current_node].current_name << " at " << ctime(&my_time);
		}
		jobs[current_node].job_has_started = true;

		if (!pipeline.runJob(myjob, current_job, false, is_continue, true, do_overwrite_current, error_message))
			REPORT_ERROR(error_message);

		// Write out current status, but maintain lock on the directory!
		write();

		// Wait for job to finish
		bool is_failure = false;
		bool is_aborted = false;
		pipeline.waitForJobToFinish(current_job, is_failure, is_aborted);


		std::string message = "";
		if (is_failure) message = " + Stopping schedule due to job " + jobs[current_node].current_name + " failing with an error ...";
		else if (is_aborted) message = " + Stopping schedule due to user abort of job " + jobs[current_node].current_name + " ...";
		if (message != "")
		{
			schedulerSendEmail(message, "Schedule: " + name);
			std::cout << message << std::endl;
			is_ok = false;
			break;
		}

		has_more_jobs = gotoNextJob();
	} // end while has_more_jobs

	if (is_ok) schedulerSendEmail("Finished successfully!", "Schedule: " + name);

	if (verb > 0)
	{
		if (exists(name + RELION_JOB_ABORT_NOW))
			std::cout << " + Found an ABORT signal... " << std::endl;
		std::cout << " + Scheduler " << name << " stops now... " << std::endl;
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++"  << std::endl << std::endl;
	}
}

void Schedule::unlock()
{
	FileName name_wo_dir = name;
	name_wo_dir = name_wo_dir.beforeLastOf("/");
	FileName dir_lock = ".relion_lock_schedule_" + name_wo_dir.afterLastOf("/");
	FileName fn_lock = dir_lock + "/lock_schedule";;

	if (exists(fn_lock))
	{
		if (std::remove(fn_lock.c_str()))
			REPORT_ERROR("ERROR: in removing lock file "+fn_lock);
		if (rmdir(dir_lock.c_str()))
			REPORT_ERROR("ERROR: in removing lock directory "+dir_lock);
	}
}

void Schedule::abort()
{
	std::cout << " Aborting schedule while at: " << current_node << std::endl;
	if (isJob(current_node))
	{
		touch(jobs[current_node].current_name + RELION_JOB_ABORT_NOW);
		std::cerr << " Touched file: " << jobs[current_node].current_name << RELION_JOB_ABORT_NOW << std::endl;
	}
	touch(name + RELION_JOB_ABORT_NOW);
	std::cerr << " Touched file: " << name << RELION_JOB_ABORT_NOW << std::endl;
	std::cerr << " Now wait for the job and the scheduler to abort ..." << std::endl;
}
