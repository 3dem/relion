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
TimeStamp global_timestamp;
bool has_annotated_time = false;

std::map<std::string, SchedulerBooleanVariable> scheduler_bools;
std::map<std::string, SchedulerFloatVariable> scheduler_floats;
std::map<std::string, SchedulerStringVariable> scheduler_strings;
std::map<std::string, SchedulerOperator> operators;

bool isBooleanVariable(std::string _name)
{
	return (scheduler_bools.find(_name) != scheduler_bools.end());
}

bool isFloatVariable(std::string _name)
{
	return (scheduler_floats.find(_name) != scheduler_floats.end());
}

bool isStringVariable(std::string _name)
{
	return (scheduler_strings.find(_name) != scheduler_strings.end());
}

bool isOperator(std::string _name)
{
	return (operators.find(_name) != operators.end());
}

SchedulerOperator::SchedulerOperator(std::string _type, std::string _input1, std::string _input2, std::string _output)
{

	type = _type;

	// Check output
	if ((type == SCHEDULE_BOOLEAN_OPERATOR_GT_CONST || type == SCHEDULE_BOOLEAN_OPERATOR_GT_VAR ||
		type == SCHEDULE_BOOLEAN_OPERATOR_LT_CONST || type == SCHEDULE_BOOLEAN_OPERATOR_LT_VAR ||
		type == SCHEDULE_BOOLEAN_OPERATOR_EQ_CONST || type == SCHEDULE_BOOLEAN_OPERATOR_EQ_VAR ||
		type == SCHEDULE_BOOLEAN_OPERATOR_AND || type == SCHEDULE_BOOLEAN_OPERATOR_OR ||
		type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS) && !isBooleanVariable(_output))
		REPORT_ERROR("ERROR: boolean operator does not have valid boolean output: " + _output);
	if ((type == SCHEDULE_FLOAT_OPERATOR_PLUS_VAR || type == SCHEDULE_FLOAT_OPERATOR_PLUS_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_MINUS_VAR || type == SCHEDULE_FLOAT_OPERATOR_MINUS_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_MULT_VAR || type == SCHEDULE_FLOAT_OPERATOR_MULT_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_VAR || type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST_INV) && !isFloatVariable(_output))
		REPORT_ERROR("ERROR: float operator does not have valid float output: " + _output);
	if ((type == SCHEDULE_STRING_OPERATOR_TOUCH_FILE ||
		 type == SCHEDULE_STRING_OPERATOR_COPY_FILE ||
		 type == SCHEDULE_STRING_OPERATOR_MOVE_FILE ||
		 type == SCHEDULE_STRING_OPERATOR_DELETE_FILE)	&& ! isStringVariable(_output))
		REPORT_ERROR("ERROR: string operator does not have valid string output: " + _output);

	// Check input1
	if ((type == SCHEDULE_BOOLEAN_OPERATOR_AND ||
		type == SCHEDULE_BOOLEAN_OPERATOR_OR )  && !isBooleanVariable(_input1))
		REPORT_ERROR("ERROR: boolean operator does not have valid boolean input1: " + _input1);
	if ((type == SCHEDULE_FLOAT_OPERATOR_PLUS_VAR || type == SCHEDULE_FLOAT_OPERATOR_PLUS_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_MINUS_VAR || type == SCHEDULE_FLOAT_OPERATOR_MINUS_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_MULT_VAR || type == SCHEDULE_FLOAT_OPERATOR_MULT_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_VAR || type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST_INV ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_GT_CONST ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_GT_VAR ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_LT_CONST ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_LT_VAR) && !isFloatVariable(_input1))
		REPORT_ERROR("ERROR: float operator does not have valid float input1: " + _input1);
	if ((type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS ||
		 type == SCHEDULE_STRING_OPERATOR_COPY_FILE ||
		 type == SCHEDULE_STRING_OPERATOR_MOVE_FILE) && ! isStringVariable(_input1))
		REPORT_ERROR("ERROR: operator does not have valid string input1: " + _input1);

	// Check input2
	if ((type == SCHEDULE_BOOLEAN_OPERATOR_AND ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_OR) && !isBooleanVariable(_input2))
		REPORT_ERROR("ERROR: boolean operator does not have valid boolean input2: " + _input2);
	if ((type == SCHEDULE_BOOLEAN_OPERATOR_GT_VAR ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_LT_VAR ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_EQ_VAR ||
		 type == SCHEDULE_FLOAT_OPERATOR_PLUS_VAR ||
		 type == SCHEDULE_FLOAT_OPERATOR_MINUS_VAR ||
		 type == SCHEDULE_FLOAT_OPERATOR_MULT_VAR ||
		 type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_VAR) && !isFloatVariable(_input2))
		REPORT_ERROR("ERROR: operator does not have valid float input2: " + _input2);
	float floatval;
	if ((type == SCHEDULE_BOOLEAN_OPERATOR_GT_CONST ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_LT_CONST ||
		 type == SCHEDULE_BOOLEAN_OPERATOR_EQ_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_PLUS_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_MINUS_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_MULT_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST ||
		 type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST_INV) && !sscanf(_input2.c_str(), "%f", &floatval))
		REPORT_ERROR("ERROR: operator does not have valid number input2: " + _input2);

	input1 = _input1;
	input2 = _input2;
	output = _output;
}

void SchedulerOperator::performOperation() const
{
	if (type == SCHEDULE_BOOLEAN_OPERATOR_AND)
	{
		scheduler_bools[output].value = (scheduler_bools[input1].value && scheduler_bools[input2].value);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_OR)
	{
		scheduler_bools[output].value = (scheduler_bools[input1].value || scheduler_bools[input2].value);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_NOT)
	{
		scheduler_bools[output].value = (!(scheduler_bools[input1].value));
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_VAR)
	{
		scheduler_bools[output].value = (scheduler_floats[input1].value > scheduler_floats[input2].value);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_LT_VAR)
	{
		scheduler_bools[output].value = (scheduler_floats[input1].value < scheduler_floats[input2].value);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_EQ_VAR)
	{
		scheduler_bools[output].value = (fabs(scheduler_floats[input1].value - scheduler_floats[input2].value) < 1E-8);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_CONST)
	{
		scheduler_bools[output].value = (scheduler_floats[input1].value > textToFloat(input2));
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_LT_CONST)
	{
		scheduler_bools[output].value = (scheduler_floats[input1].value < textToFloat(input2));
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_EQ_CONST)
	{
		scheduler_bools[output].value = (fabs(scheduler_floats[input1].value - textToFloat(input2)) < 1E-8);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS)
	{
		scheduler_bools[output].value = (exists(scheduler_strings[input1].value));
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_VAR)
	{
		scheduler_floats[output].value = scheduler_floats[input1].value + scheduler_floats[input2].value;
	}
	else if (type ==SCHEDULE_FLOAT_OPERATOR_MINUS_VAR)
	{
		scheduler_floats[output].value = scheduler_floats[input1].value - scheduler_floats[input2].value;
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_MULT_VAR)
	{
		scheduler_floats[output].value = scheduler_floats[input1].value * scheduler_floats[input2].value;
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_VAR)
	{
		scheduler_floats[output].value = scheduler_floats[input1].value / scheduler_floats[input2].value;
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_CONST)
	{
		scheduler_floats[output].value = scheduler_floats[input1].value + textToFloat(input2);
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_MINUS_CONST)
	{
		scheduler_floats[output].value = scheduler_floats[input1].value - textToFloat(input2);
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_MULT_CONST)
	{
		scheduler_floats[output].value = scheduler_floats[input1].value * textToFloat(input2);
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST)
	{
		scheduler_floats[output].value = scheduler_floats[input1].value / textToFloat(input2);
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST_INV)
	{
		scheduler_floats[output].value = textToFloat(input2) / scheduler_floats[input1].value;
	}
	else if (type == SCHEDULE_STRING_OPERATOR_TOUCH_FILE)
	{
		std::cout << " Touching file " << scheduler_strings[output].value << std::endl;
		touch(scheduler_strings[output].value);
	}
	else if (type == SCHEDULE_STRING_OPERATOR_COPY_FILE)
	{
		std::cout << " Copying file " << scheduler_strings[input1].value << " to " << scheduler_strings[output].value << std::endl;
		copy(scheduler_strings[input1].value, scheduler_strings[output].value);
	}
	else if (type == SCHEDULE_STRING_OPERATOR_MOVE_FILE)
	{
		std::cout << " Moving file " << scheduler_strings[input1].value << " to " << scheduler_strings[output].value << std::endl;
		move(scheduler_strings[input1].value, scheduler_strings[output].value);
	}
	else if (type == SCHEDULE_STRING_OPERATOR_DELETE_FILE)
	{
		std::cout << " Deleting file " << scheduler_strings[output].value << std::endl;
		delete(scheduler_strings[output].value.c_str());
	}
	else if (type == SCHEDULE_WAIT_OPERATOR_SINCE_LAST_TIME)
	{
		if (has_annotated_time)
		{
			RFLOAT elapsed = elapsed_time(global_timestamp);
			RFLOAT wait_seconds =  textToFloat(input1)-elapsed;
			std::cout << " Waiting for " << wait_seconds << " seconds ..." << std::endl;
			sleep(wait_seconds);
			std::cout << " Finished waiting!" << std::endl;
			annotate_time(&global_timestamp);
			has_annotated_time =true;
		}
	}
	else
		REPORT_ERROR("ERROR: unrecognised Operator type:" + type);

}


std::string SchedulerOperator::getName()
{
	if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_CONST || type == SCHEDULE_BOOLEAN_OPERATOR_GT_VAR) return input1 + "_GT_" + input2;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_LT_CONST || type == SCHEDULE_BOOLEAN_OPERATOR_LT_VAR) return input1 + "_LT_" + input2;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_EQ_CONST || type == SCHEDULE_BOOLEAN_OPERATOR_EQ_VAR) return input1 + "_EQ_" + input2;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_AND) return input1 + "_AND_" + input2;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_OR) return input1 + "_OR_" + input2;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_NOT) return "NOT_" + input1;
	if (type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS) return "EXISTS_" + input1;
	if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_VAR || type == SCHEDULE_FLOAT_OPERATOR_PLUS_CONST) return input1 + "_PLUS_" + input2;
	if (type == SCHEDULE_FLOAT_OPERATOR_MINUS_VAR || type == SCHEDULE_FLOAT_OPERATOR_MINUS_CONST) return input1 + "_MINUS_" + input2;
	if (type == SCHEDULE_FLOAT_OPERATOR_MULT_VAR || type == SCHEDULE_FLOAT_OPERATOR_MULT_CONST) return input1 + "_MULT_" + input2;
	if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_VAR || type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST) return input1 + "_DIV_" + input2;
	if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST_INV) return input2 + "_DIV_" + input1;
	if (type == SCHEDULE_STRING_OPERATOR_TOUCH_FILE) return "TOUCH_" + input1;
	if (type == SCHEDULE_STRING_OPERATOR_COPY_FILE) return "COPY_" + input1 + "_" + input2;
	if (type == SCHEDULE_STRING_OPERATOR_MOVE_FILE) return "MOVE_" + input1 + "_" + input2;
	if (type == SCHEDULE_STRING_OPERATOR_DELETE_FILE) return "DELETE_" + input1;
	else
		REPORT_ERROR("ERROR: unrecognised Operator type:" + type);
}

std::string SchedulerEdge::getOutputNode() const
{
	if (is_fork) return (scheduler_bools[myBooleanVariable].value) ?  outputNode : outputNodeFalse;
	else return outputNode;
}

void Schedule::clear()
{
	current_node = "undefined";
	original_start_node = "undefined";
	name = "undefined";
	email_address = "undefined";
	scheduler_bools.clear();
	scheduler_floats.clear();
	scheduler_strings.clear();
	nodes.clear();
	edges.clear();
	operators.clear();
}

std::string Schedule::findNodeByCurrentName(std::string _name)
{
	std::map<std::string, SchedulerNode>::iterator it;

	for ( it = nodes.begin(); it != nodes.end(); it++ )
	{
	    if (it->second.current_name == _name)
	    	return it->second.current_name;
	}
	REPORT_ERROR("ERROR: cannot find node: " + _name);
	return "";
}

bool SchedulerNode::performOperation()
{
	if (type == SCHEDULE_NODE_TYPE_OPERATOR)
	{
		operators[myOperator].performOperation();
		return true;
	}
	else if (type == SCHEDULE_NODE_TYPE_EXIT)
	{
		std::cout << " Reached an exit of the Schedule ..." << std::endl;
		return true;
	}
    else if (type == SCHEDULE_NODE_TYPE_JOB)
    {
		return false;
	}
    else
    	REPORT_ERROR("ERROR: unrecognised node type: " + type);

}

void Schedule::read(FileName fn)
{
	if (fn == "") fn = name + "schedule.star";

	// Clear current model
	clear();

	// Open input file
	std::ifstream in(fn.data(), std::ios_base::in);
	if (in.fail())
		REPORT_ERROR( (std::string) "Schedule::read: File " + fn + " cannot be read." );

	// For reading: do the nodes before the general table, in order to set current_node and original_start_node
	MetaDataTable MD;
	MD.readStar(in, "schedule_general");
	std::string current_node_name, original_start_node_name;
	MD.getValue(EMDL_SCHEDULE_GENERAL_NAME, name);
	MD.getValue(EMDL_SCHEDULE_GENERAL_CURRENT_NODE, current_node);
	MD.getValue(EMDL_SCHEDULE_GENERAL_ORIGINAL_START_NODE, original_start_node);
	MD.getValue(EMDL_SCHEDULE_GENERAL_EMAIL, email_address);

	MD.clear();

	MD.readStar(in, "schedule_floats");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string name;
		RFLOAT value, original_value;
		MD.getValue(EMDL_SCHEDULE_VAR_FLOAT_NAME, name);
		MD.getValue(EMDL_SCHEDULE_VAR_FLOAT_VALUE, value);
		MD.getValue(EMDL_SCHEDULE_VAR_FLOAT_ORI_VALUE, original_value);
		SchedulerFloatVariable myval(value, original_value);
		scheduler_floats[name] = myval;
	}
	MD.clear();

	MD.readStar(in, "schedule_bools");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string name;
		bool value, original_value;
		MD.getValue(EMDL_SCHEDULE_VAR_BOOL_NAME, name);
		MD.getValue(EMDL_SCHEDULE_VAR_BOOL_VALUE, value);
		MD.getValue(EMDL_SCHEDULE_VAR_BOOL_ORI_VALUE, original_value);
		SchedulerBooleanVariable myval(value, original_value);
		scheduler_bools[name] = myval;
	}
	MD.clear();

	MD.readStar(in, "schedule_strings");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string name;
		FileName value, original_value;
		MD.getValue(EMDL_SCHEDULE_VAR_STRING_NAME, name);
		MD.getValue(EMDL_SCHEDULE_VAR_STRING_VALUE, value);
		MD.getValue(EMDL_SCHEDULE_VAR_STRING_ORI_VALUE, original_value);
		SchedulerStringVariable myval(value, original_value);
		scheduler_strings[name] = myval;
	}
	MD.clear();

	MD.readStar(in, "schedule_operators");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string name, type, input1, input2, output;
		RFLOAT constant;

		MD.getValue(EMDL_SCHEDULE_OPERATOR_NAME, name);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_TYPE, type);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_INPUT1, input1);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_INPUT2, input2);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_OUTPUT, output);
		SchedulerOperator myval(type, input1, input2, output);
		operators[name] = myval;
	}
	MD.clear();

	MD.readStar(in, "schedule_nodes");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string name, ori_name, mode, type;
		bool has_started;

		MD.getValue(EMDL_SCHEDULE_NODE_NAME, name);
		MD.getValue(EMDL_SCHEDULE_NODE_ORI_NAME, ori_name);
		MD.getValue(EMDL_SCHEDULE_NODE_JOB_MODE, mode);
		MD.getValue(EMDL_SCHEDULE_NODE_JOB_HAS_STARTED, has_started);
		MD.getValue(EMDL_SCHEDULE_NODE_TYPE, type);

		SchedulerNode myval(name, type, mode, has_started);
		nodes[ori_name] = myval;
	}
	MD.clear();


	MD.readStar(in, "schedule_edges");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		int number;
		std::string inputname, outputname, outputname_false, bool_name;
		bool is_fork;

		MD.getValue(EMDL_SCHEDULE_EDGE_NUMBER, number);
		MD.getValue(EMDL_SCHEDULE_EDGE_INPUT, inputname);
		MD.getValue(EMDL_SCHEDULE_EDGE_OUTPUT, outputname);
		MD.getValue(EMDL_SCHEDULE_EDGE_IS_FORK, is_fork);
		MD.getValue(EMDL_SCHEDULE_EDGE_OUTPUT_FALSE, outputname_false);
		MD.getValue(EMDL_SCHEDULE_EDGE_BOOLEAN, bool_name);
		SchedulerEdge myval(inputname, outputname, is_fork, bool_name, outputname_false);
		edges.push_back(myval);
	}
	MD.clear();

	// Close file handler
	in.close();

	// Also read in the schedule_pipeline (no need to lock now?)
	schedule_pipeline.read();

}


void Schedule::write(FileName fn)
{

	if (fn == "") fn = name + "schedule.star";

	// B. Write STAR file with the entire schedule
	std::ofstream  fh;
	fh.open((fn).c_str(), std::ios::out);
	if (!fh)
		REPORT_ERROR( (std::string)"Schedule::write: Cannot write file: " + fn);

	// For undefined values
	std::string str_aux="undefined";
	RFLOAT float_aux=0.;
	bool bool_aux=false;

	MetaDataTable MDgeneral;
	MDgeneral.setName("schedule_general");
	MDgeneral.setIsList(true);
	MDgeneral.addObject();
	MDgeneral.setValue(EMDL_SCHEDULE_GENERAL_NAME, name);
	MDgeneral.setValue(EMDL_SCHEDULE_GENERAL_CURRENT_NODE, current_node);
	MDgeneral.setValue(EMDL_SCHEDULE_GENERAL_ORIGINAL_START_NODE, original_start_node);
	MDgeneral.setValue(EMDL_SCHEDULE_GENERAL_EMAIL, email_address);
	MDgeneral.write(fh);

	if (scheduler_floats.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_floats");
		std::map<std::string, SchedulerFloatVariable>::iterator it;
		for ( it = scheduler_floats.begin(); it != scheduler_floats.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_VAR_FLOAT_NAME, it->first);
			MD.setValue(EMDL_SCHEDULE_VAR_FLOAT_VALUE, it->second.value);
			MD.setValue(EMDL_SCHEDULE_VAR_FLOAT_ORI_VALUE, it->second.original_value);
		}
		MD.write(fh);
	}

	if (scheduler_bools.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_bools");
		std::map<std::string, SchedulerBooleanVariable>::iterator it;
		for ( it = scheduler_bools.begin(); it != scheduler_bools.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_VAR_BOOL_NAME, it->first);
			MD.setValue(EMDL_SCHEDULE_VAR_BOOL_VALUE, it->second.value);
			MD.setValue(EMDL_SCHEDULE_VAR_BOOL_ORI_VALUE, it->second.original_value);
		}
		MD.write(fh);
	}

	if (scheduler_strings.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_strings");
		std::map<std::string, SchedulerStringVariable>::iterator it;
		for ( it = scheduler_strings.begin(); it != scheduler_strings.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_VAR_STRING_NAME, it->first);
			MD.setValue(EMDL_SCHEDULE_VAR_STRING_VALUE, it->second.value);
			MD.setValue(EMDL_SCHEDULE_VAR_STRING_ORI_VALUE, it->second.original_value);
		}
		MD.write(fh);
	}

	if (operators.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_operators");
		std::map<std::string, SchedulerOperator>::iterator it;
		for ( it = operators.begin(); it != operators.end(); it++ )
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

	if (nodes.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_nodes");
		std::map<std::string, SchedulerNode>::iterator it;
		for ( it = nodes.begin(); it != nodes.end(); it++ )
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_NODE_ORI_NAME, it->first);
			MD.setValue(EMDL_SCHEDULE_NODE_TYPE, it->second.type);
			MD.setValue(EMDL_SCHEDULE_NODE_NAME, it->second.current_name);
			MD.setValue(EMDL_SCHEDULE_NODE_JOB_MODE, it->second.mode);
			MD.setValue(EMDL_SCHEDULE_NODE_JOB_HAS_STARTED, it->second.job_has_started);
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
			MD.setValue(EMDL_SCHEDULE_EDGE_NUMBER, i+1); // start counting at 1!
			MD.setValue(EMDL_SCHEDULE_EDGE_INPUT, edges[i].inputNode);
			MD.setValue(EMDL_SCHEDULE_EDGE_OUTPUT, edges[i].outputNode);
			MD.setValue(EMDL_SCHEDULE_EDGE_IS_FORK, edges[i].is_fork);

			if (edges[i].is_fork)
			{
				MD.setValue(EMDL_SCHEDULE_EDGE_OUTPUT_FALSE, edges[i].outputNodeFalse);
				MD.setValue(EMDL_SCHEDULE_EDGE_BOOLEAN, edges[i].myBooleanVariable);
			}
			else
			{
				MD.setValue(EMDL_SCHEDULE_EDGE_OUTPUT_FALSE, str_aux);
				MD.setValue(EMDL_SCHEDULE_EDGE_BOOLEAN, str_aux);
			}
		}
		MD.write(fh);
	}

	// Close the file handler
	fh.close();

	// Also write out the schedule_pipeline (no need to lock now?)
	schedule_pipeline.write();

}

// Reset all variables to their original value
void Schedule::reset()
{

	{
		std::map<std::string, SchedulerFloatVariable>::iterator it;
		for ( it = scheduler_floats.begin(); it != scheduler_floats.end(); it++ )
			it->second.value = it->second.original_value;
	}

	{
		std::map<std::string, SchedulerBooleanVariable>::iterator it;
		for ( it = scheduler_bools.begin(); it != scheduler_bools.end(); it++ )
	        it->second.value = it->second.original_value;
	}

	{
		std::map<std::string, SchedulerStringVariable>::iterator it;
		for ( it = scheduler_strings.begin(); it != scheduler_strings.end(); it++ )
	        it->second.value = it->second.original_value;
	}

	{
		std::map<std::string, SchedulerNode>::iterator it;
		for ( it = nodes.begin(); it != nodes.end(); it++ )
		{
			it->second.current_name = it->first;
			it->second.job_has_started = false;
		}
    }

    current_node = "undefined";
}

bool Schedule::gotoNextNode()
{
    if (current_node == "undefined")
    {
        if (original_start_node == "undefined")
        	REPORT_ERROR("ERROR: the starting node was not defined...");
    	current_node = original_start_node; // go to first node in the list
        std::cout << " Setting current_node to original_start_node: " << original_start_node << std::endl;
    	return true;
    }

    for (int i = 0; i < edges.size(); i++)
    {
        if (edges[i].inputNode == current_node)
        {
            current_node = edges[i].getOutputNode();
            return (current_node == "undefined") ? false : true;
        }
    }

    return false;
}

bool Schedule::gotoNextJob(FileName &job_name, FileName &original_name, std::string &mode, bool &has_started)
{

    // This loops through the next Nodes until encountering a JOB
    while (gotoNextNode())
    {
    	// If this node is an operator, perform its operation, else get the Job
		if (nodes[current_node].performOperation())
		{
			if (nodes[current_node].type == SCHEDULE_NODE_TYPE_EXIT)
			{
				sendEmail("Finished successfully!");
				return false;
			}
		}
		else // this is a job, get its current_name and options
		{
            std::cout << " Now preparing Job: " << nodes[current_node].current_name << std::endl;
			job_name = nodes[current_node].current_name;
            original_name = current_node;
            mode = nodes[current_node].mode;
            has_started = nodes[current_node].job_has_started;
            return true;
        }
    }

    return false;
}

bool Schedule::isNode(std::string _name)
{
	return (nodes.find(_name) != nodes.end());
}


bool Schedule::isJob(std::string _name)
{
	if (isNode(_name))
		return (nodes[_name].type == SCHEDULE_NODE_TYPE_JOB);
	return false;
}

void Schedule::addVariable(std::string name, FileName value)
{
	if (value == "true" || value == "True" || value == "false" || value == "False") // or a boolean?
	{
		bool myval = (value == "true" || value == "True");
		addBooleanVariable(name, myval);
	}
	else
	{
		float floatval;
		if (sscanf(value.c_str(), "%f", &floatval)) // is this a number?
		{
			std::cerr << name << " = " << floatval << std::endl;
			addFloatVariable(name, floatval);
		}
		else
		{
			addStringVariable(name, value);
		}
	}
}

void Schedule::addFloatVariable(std::string name, RFLOAT value)
{
	if (isFloatVariable(name))
		REPORT_ERROR("ERROR: trying to add a float variable with a name that already exists: " + name);

	SchedulerFloatVariable myvar(value, value);
	scheduler_floats[name] = myvar;
}

void Schedule::addBooleanVariable(std::string name, bool value)
{
	if (isBooleanVariable(name))
		REPORT_ERROR("ERROR: trying to add a boolean variable with a name that already exists: " + name);

	SchedulerBooleanVariable myvar(value, value);
	scheduler_bools[name] = myvar;
}

void Schedule::addStringVariable(std::string name, FileName value)
{
	if (isStringVariable(name))
		REPORT_ERROR("ERROR: trying to add a string variable with a name that already exists: " + name);

	SchedulerStringVariable myvar(value, value);
	scheduler_strings[name] = myvar;
}

void Schedule::addOperatorNode(std::string type, std::string input_name, std::string input2_name, std::string output_name)
{
	SchedulerOperator myop(type, input_name, input2_name, output_name);
	std::string myname = myop.getName();
	operators[myname] = myop;
	SchedulerNode mynode(myname, SCHEDULE_NODE_TYPE_OPERATOR, type);
}

void Schedule::addJobNode(RelionJob &myjob, std::string jobname, std::string mode)
{

	// Check whether the jobname is unique
	if (isNode(jobname))
		REPORT_ERROR("ERROR: trying to add a JobNode that already exists...");

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

	SchedulerNode mynode(jobname, SCHEDULE_NODE_TYPE_JOB, mode, false);
	nodes[jobname] = mynode;
}

void Schedule::addExitNode()
{
	SchedulerNode mynode;
	nodes["exit"] = mynode;
}

void Schedule::sendEmail(std::string message)
{
	if (email_address != "undefined")
	{
		std::string mysubject = "Schedule: " + name;
		std::string command = "echo \"" + message + "\" | mail -s \"" + mysubject + "\" -r RELION " + email_address;
		int res = system(command.c_str());
	}
}

void Schedule::addEdge(std::string inputnode_name, std::string outputnode_name)
{
	SchedulerEdge myval(inputnode_name, outputnode_name);
	edges.push_back(myval);
}

void Schedule::addFork(std::string inputnode_name, std::string mybool_name, std::string outputnode_name, std::string outputnode_name_if_false)
{
	SchedulerEdge myval(inputnode_name, outputnode_name, true, mybool_name, outputnode_name_if_false);
	edges.push_back(myval);
}

bool Schedule::isValid()
{
	// TODO: check if duplicate edges, forks or operators exist....

	// Check original_start node was set

	// Check Scheduler ends with an exit

}

RelionJob Schedule::copyNewJobFromSchedulePipeline(FileName original_job_name)
{


	RelionJob job;
	bool dummy;
	job.read(name + original_job_name + '/', dummy, true); // true means initialise the job

	// Check where this job gets its input from: change names from local scheduler ones to the current pipeline
	int local_process = schedule_pipeline.findProcessByName(name + original_job_name + '/');
	// Loop over all input nodes to this job
	for (int inode = 0; inode <  schedule_pipeline.processList[local_process].inputNodeList.size(); inode++)
	{
		int mynode = schedule_pipeline.processList[local_process].inputNodeList[inode];
		// find from which pipeline_scheduler job this jobs gets its input nodes
		int output_from_process =  schedule_pipeline.nodeList[mynode].outputFromProcess;
		if (output_from_process < 0)
			REPORT_ERROR("ERROR: cannot find outputProcess of node: " +  schedule_pipeline.nodeList[mynode].name);
		// Get the original current_name in the pipeline_scheduler of this job
		FileName my_ori_name = schedule_pipeline.processList[output_from_process].name;
		// Remove leading directory and tailing slash to get the process current_name in the pipeline_scheduler
		FileName my_process_name = (my_ori_name.afterFirstOf(name)).beforeLastOf("/");
		// find that process in the nodes, and get its current current_name
		std::string my_current_name =nodes[my_process_name].current_name;
		// Change all instances of the my_ori_name to my_current_name in the joboptions of this job
		for (std::map<std::string,JobOption>::iterator it=job.joboptions.begin(); it!=job.joboptions.end(); ++it)
		{
			FileName mystring = (it->second).value;
			mystring.replaceAllSubstrings(my_ori_name, my_current_name);
			(it->second).value = mystring;
		}
	}

	return job;
}

// Modify a job to set variables from the Scheduler
void Schedule::setVariablesInJob(RelionJob &job, FileName original_job_name)
{

	RelionJob ori_job;
	bool dummy;
	ori_job.read(name + original_job_name + '/', dummy, true);

	// Check whether there are any options with a value containing &&, which is the sign for inserting Scheduler variables
	for (std::map<std::string,JobOption>::iterator it=ori_job.joboptions.begin(); it!=ori_job.joboptions.end(); ++it)
	{
		FileName mystring = (it->second).value;
		if (mystring.contains("$$"))
		{
			mystring = mystring.afterFirstOf("$$");

			// Found an option that needs replacement! Now find which variable to insert
			std::string my_value;
			if (isBooleanVariable(mystring))
			{
				if (job.joboptions[it->first].joboption_type != JOBOPTION_BOOLEAN)
					REPORT_ERROR(" ERROR: trying to set a BooleanVariable: " + mystring + " into a non-boolean option: " + it->first);

				my_value = (scheduler_bools[mystring].value) ? "Yes" : "No";
			}
			else if (isFloatVariable(mystring))
			{
				if (job.joboptions[it->first].joboption_type == JOBOPTION_BOOLEAN)
					REPORT_ERROR(" ERROR: trying to set FloatVariable: " + mystring + " into a boolean option: " + it->first);

				my_value = floatToString(scheduler_floats[mystring].value);
			}
			else if (isStringVariable(mystring))
			{
				if (job.joboptions[it->first].joboption_type == JOBOPTION_BOOLEAN)
					REPORT_ERROR(" ERROR: trying to set StringVariable: " + mystring + " into a boolean option: " + it->first);
				if (job.joboptions[it->first].joboption_type == JOBOPTION_SLIDER)
					REPORT_ERROR(" ERROR: trying to set StringVariable: " + mystring + " into a slider option: " + it->first);

				my_value = scheduler_strings[mystring].value;
			}
			else
				REPORT_ERROR(" ERROR: variable in job is not part of this Schedule: " + mystring);

			job.joboptions[it->first].value = my_value;
			std::cout << " Setting joboption " << it->first << " to " << my_value << " based on variable: " << mystring << std::endl;
		}
	}

}

void Schedule::run(PipeLine &pipeline)
{
    // go through all nodes
    FileName job_name, original_job_name, mode;
    bool job_has_started;
    while (gotoNextJob(job_name, original_job_name, mode, job_has_started))
    {
        RelionJob myjob;
        bool is_continue, do_overwrite_current, dummy;
    	int current_job;
    	if (!job_has_started || mode == SCHEDULE_NODE_JOB_MODE_NEW)
    	{

    		// Copy the job inside the schedule_pipeline to a new job inside the pipeline we are actually running in
    		// This function also takes care of fixing the names of the inputNodes
    		myjob = copyNewJobFromSchedulePipeline(original_job_name);
        	// Now add this job to the pipeline we will actually be running in
        	current_job = pipeline.addScheduledJob(myjob);
        	is_continue = false;
        	do_overwrite_current = false;

        	// Set the current_name of the current node now
        	nodes[current_node].current_name = pipeline.processList[current_job].name;
    	}
    	else if (mode == SCHEDULE_NODE_JOB_MODE_CONTINUE || mode == SCHEDULE_NODE_JOB_MODE_OVERWRITE)
    	{
    		is_continue = (mode == SCHEDULE_NODE_JOB_MODE_CONTINUE);
    		do_overwrite_current = (mode == SCHEDULE_NODE_JOB_MODE_OVERWRITE);

    		current_job = pipeline.findProcessByName(job_name);
    		if (current_job < 0)
				REPORT_ERROR("ERROR: RunSchedule cannot find process with name: " + job_name);

    		// Read the job from the pipeline we are running in
    		if (!myjob.read(pipeline.processList[current_job].name, dummy, true)) // true means also initialise the job
					REPORT_ERROR("There was an error reading job: " + pipeline.processList[current_job].name);

    	}
    	else
    		REPORT_ERROR("ERROR: unrecognised mode for running a new process: " + mode);


    	// This function replaces variable calls starting with a '&&' from the original_job into the current_job
    	setVariablesInJob(myjob, original_job_name);

    	// Check whether the input nodes are there, before executing the job
		for (long int inode = 0; inode < pipeline.processList[current_job].inputNodeList.size(); inode++)
		{
			long int mynode = pipeline.processList[current_job].inputNodeList[inode];
			while (!exists(pipeline.nodeList[mynode].name))
			{
				std::cerr << " + -- Warning " << pipeline.nodeList[mynode].name << " does not exist. Waiting 10 seconds ... " << std::endl;
				sleep(10);
			}
		}

		// Now actually run the Scheduled job
		std::string error_message;
		if (!pipeline.runJob(myjob, current_job, false, is_continue, true, do_overwrite_current, error_message))
			REPORT_ERROR(error_message);

		// Wait for job to finish
		bool is_failure = false;
		bool is_aborted = false;
		pipeline.waitForJobToFinish(current_job, is_failure, is_aborted);

		std::string message = "";
		if (is_failure) message = " Stopping schedule due to job " + job_name + " failing with an error ...";
		else if (is_aborted) message = " Stopping schedule due to user abort of job " + job_name + " ...";
		if (message != "")
		{
			sendEmail(message);
			std::cerr << message << std::endl;
			break;
		}
		else
		{
			// job finished successfully, write out the updated scheduler file
			nodes[current_node].job_has_started = true;
			write();
		}
    } // end while gotoNextJob

    std::cout << " Scheduler " << name << " has finished now... " << std::endl;
}

void Schedule::abort()
{
    // TODO: make an abort mechanism with a touch file, specific for the Scheduler!!

	// TODO: also make a mechanism that the same scheduler cannot be run twice simultaneously!




    /*
	// Abort by aborting the current job
	FileName job_name, original_job_name, mode;
    bool has_done = false, has_started = false;

    // Go to the next job in the Schedule (starting from current_job)
    // Place an abort in each of the Jobs directories
    gotoNextJob(job_name, original_job_name, mode, has_started);
    int current_job = pipeline.findProcessByName(job_name);
    if (pipeline.processList[current_job].status == PROC_RUNNING)
    {
        has_done = true;
    	touch(processList[current_job].current_name + RELION_JOB_ABORT_NOW);
        std::cout << " Marking job " << processList[current_job].current_name << " for abortion; schedule should abort too ... " << std::endl;
    }


    while ()
    {
        int current_job = findProcessByName(job_name);
        if (processList[current_job].status == PROC_RUNNING)
        {
            has_done = true;
        	touch(processList[current_job].current_name + RELION_JOB_ABORT_NOW);
            std::cout << " Marking job " << processList[current_job].current_name << " for abortion; schedule should abort too ... " << std::endl;
        }
    }

    if (!has_done)
    	std::cout << " This schedule seems to have finished already ..." << std::endl;
	*/

}

