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

SchedulerBooleanVariable* Schedule::findBooleanVariable(std::string _name)
{
	if (_name == "NULL")
		return NULL;
	for (int i = 0; i < bools.size(); i++)
		if (bools[i].name == _name) return &bools[i];
	REPORT_ERROR("ERROR: cannot find boolean variable: " + _name);
	return NULL;
}

SchedulerFloatVariable* Schedule::findFloatVariable(std::string _name)
{
	for (int i = 0; i < floats.size(); i++)
		if (floats[i].name == _name) return &floats[i];
	REPORT_ERROR("ERROR: cannot find float variable: " + _name);
	return NULL;
}

SchedulerStringVariable* Schedule::findStringVariable(std::string _name)
{
	for (int i = 0; i < strings.size(); i++)
		if (strings[i].name == _name) return &strings[i];
	REPORT_ERROR("ERROR: cannot find file variable: " + _name);
	return NULL;
}

SchedulerBooleanOperator* Schedule::findBooleanOperator(std::string _name)
{
	for (int i = 0; i < bool_ops.size(); i++)
		if (bool_ops[i].name == _name) return &bool_ops[i];
	REPORT_ERROR("ERROR: cannot find boolean operator: " + _name);
	return NULL;
}

SchedulerFloatOperator* Schedule::findFloatOperator(std::string _name)
{
	for (int i = 0; i < float_ops.size(); i++)
		if (float_ops[i].name == _name) return &float_ops[i];
	REPORT_ERROR("ERROR: cannot find float operator: " + _name);
	return NULL;
}

SchedulerStringOperator* Schedule::findStringOperator(std::string _name)
{
	for (int i = 0; i < string_ops.size(); i++)
		if (string_ops[i].name == _name) return &string_ops[i];
	REPORT_ERROR("ERROR: cannot find file operator: " + _name);
	return NULL;
}

SchedulerNode* Schedule::findNode(std::string _name)
{
	for (int i = 0; i < nodes.size(); i++)
	{
		if (nodes[i].name == _name) return &nodes[i];
	}
	REPORT_ERROR("ERROR: cannot find node: " + _name);
	return NULL;
}

SchedulerNode* Schedule::findNodeByOriginalName(std::string _name)
{
	for (int i = 0; i < nodes.size(); i++)
	{
		if (nodes[i].original_name == _name) return &nodes[i];
	}
	REPORT_ERROR("ERROR: cannot find node: " + _name);
	return NULL;
}

bool SchedulerNode::performOperation()
{
	if (type == SCHEDULE_NODE_TYPE_BOOL_OPERATOR)
	{
		(myBooleanOperator)->performOperation();
		return true;
	}
	else if (type == SCHEDULE_NODE_TYPE_FLOAT_OPERATOR)
	{
		(myFloatOperator)->performOperation();
		return true;
	}
	else if (type == SCHEDULE_NODE_TYPE_STRING_OPERATOR)
	{
		(myStringOperator)->performOperation();
		return true;
	}
	else if (type == SCHEDULE_NODE_TYPE_TIMER_WAIT)
	{
		if (!has_annotated_time)
		{
			return true;
		}
		else
		{
			RFLOAT elapsed = elapsed_time(global_timestamp);
			std::cout << " Now waiting for another " << wait_seconds-elapsed << " seconds ..." << std::endl;
			sleep(wait_seconds-elapsed);
			std::cout << " Finished waiting!" << std::endl;
			annotate_time(&global_timestamp);
			has_annotated_time =true;
			return true;
		}
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
	MD.getValue(EMDL_SCHEDULE_GENERAL_CURRENT_NODE, current_node_name);
	MD.getValue(EMDL_SCHEDULE_GENERAL_ORIGINAL_START_NODE, original_start_node_name);
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
		SchedulerFloatVariable myval(name, value, original_value);
		floats.push_back(myval);
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
		SchedulerBooleanVariable myval(name, value, original_value);
		bools.push_back(myval);
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
		SchedulerStringVariable myval(name, value, original_value);
		strings.push_back(myval);
	}
	MD.clear();

	MD.readStar(in, "schedule_boolean_operators");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string name, type, input1_name, input2_name, float1_name, float2_name, file_name, output_name;
		RFLOAT constant;

		MD.getValue(EMDL_SCHEDULE_OPERATOR_BOOL_NAME, name);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_BOOL_TYPE, type);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_BOOL_INPUT1, input1_name);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_BOOL_INPUT1, input2_name);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_BOOL_FLOAT1, float1_name);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_BOOL_FLOAT2, float2_name);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_BOOL_CONST, constant);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_BOOL_FILE, file_name);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_BOOL_OUTPUT, output_name);
		if (type == SCHEDULE_BOOLEAN_OPERATOR_AND || type == SCHEDULE_BOOLEAN_OPERATOR_OR)
		{
			SchedulerBooleanOperator myval(type, findBooleanVariable(input1_name),
					findBooleanVariable(input2_name), findBooleanVariable(output_name));
			bool_ops.push_back(myval);
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_NOT)
		{
			SchedulerBooleanOperator myval(type, findBooleanVariable(input1_name), findBooleanVariable(output_name));
			bool_ops.push_back(myval);
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_VAR || type == SCHEDULE_BOOLEAN_OPERATOR_LT_VAR || type == SCHEDULE_BOOLEAN_OPERATOR_EQ_VAR)
		{
			SchedulerBooleanOperator myval(type, findFloatVariable(float1_name), findFloatVariable(float2_name), findBooleanVariable(output_name));
			bool_ops.push_back(myval);
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_CONST || type == SCHEDULE_BOOLEAN_OPERATOR_LT_CONST || type == SCHEDULE_BOOLEAN_OPERATOR_EQ_CONST)
		{
			SchedulerBooleanOperator myval(type, findFloatVariable(float1_name), constant, findBooleanVariable(output_name));
			bool_ops.push_back(myval);
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS)
		{
			SchedulerBooleanOperator myval(type, findStringVariable(file_name), findBooleanVariable(output_name));
			bool_ops.push_back(myval);
		}
		else
			REPORT_ERROR("ERROR: unrecognised boolean operator type:" + type);
	}
	MD.clear();

	MD.readStar(in, "schedule_float_operators");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string name, type, input1_name, input2_name, output_name;
		RFLOAT myconstant;

		MD.getValue(EMDL_SCHEDULE_OPERATOR_FLOAT_NAME, name);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_FLOAT_TYPE, type);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_FLOAT_INPUT1, input1_name);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_FLOAT_INPUT1, input2_name);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_FLOAT_CONSTANT, myconstant);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_FLOAT_OUTPUT, output_name);

		if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_VAR ||
			type == SCHEDULE_FLOAT_OPERATOR_MINUS_VAR ||
			type == SCHEDULE_FLOAT_OPERATOR_MULT_VAR ||
			type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_VAR)
		{
			SchedulerFloatOperator myval(type, findFloatVariable(input1_name),
					findFloatVariable(input2_name), findFloatVariable(output_name));
			float_ops.push_back(myval);
		}
		else if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_CONST ||
				type == SCHEDULE_FLOAT_OPERATOR_MINUS_CONST ||
				type == SCHEDULE_FLOAT_OPERATOR_MULT_CONST ||
				type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST ||
				type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST_INV)
		{
			SchedulerFloatOperator myval(type, findFloatVariable(input1_name),
					myconstant, findFloatVariable(output_name));
			float_ops.push_back(myval);
		}
		else
			REPORT_ERROR("ERROR: unrecognised float operator type:" + type);
	}
	MD.clear();

	MD.readStar(in, "schedule_string_operators");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string name, type, input_name, output_name;
		RFLOAT myconstant;

		MD.getValue(EMDL_SCHEDULE_OPERATOR_STRING_NAME, name);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_STRING_TYPE, type);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_STRING_INPUT, input_name);
		MD.getValue(EMDL_SCHEDULE_OPERATOR_STRING_OUTPUT, output_name);

		if (type == SCHEDULE_STRING_OPERATOR_COPY_FILE ||
			type == SCHEDULE_STRING_OPERATOR_MOVE_FILE)
		{
			SchedulerStringOperator myval(type, findStringVariable(input_name), findStringVariable(output_name));
			string_ops.push_back(myval);
		}
		else if (type == SCHEDULE_STRING_OPERATOR_TOUCH_FILE ||
				 type == SCHEDULE_STRING_OPERATOR_DELETE_FILE)
		{
			SchedulerStringOperator myval(type, findStringVariable(input_name));
			string_ops.push_back(myval);
		}
		else
			REPORT_ERROR("ERROR: unrecognised string operator type:" + type);
	}
	MD.clear();

	MD.readStar(in, "schedule_nodes");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		std::string name, ori_name, mode, type;
		RFLOAT wait_seconds;
		bool has_started;

		MD.getValue(EMDL_SCHEDULE_NODE_NAME, name);
		MD.getValue(EMDL_SCHEDULE_NODE_ORI_NAME, ori_name);
		MD.getValue(EMDL_SCHEDULE_NODE_JOB_MODE, mode);
		MD.getValue(EMDL_SCHEDULE_NODE_JOB_HAS_STARTED, has_started);
		MD.getValue(EMDL_SCHEDULE_NODE_TYPE, type);
		MD.getValue(EMDL_SCHEDULE_NODE_WAIT_TIME, wait_seconds);

		if (type == SCHEDULE_NODE_TYPE_JOB)
		{
			SchedulerNode myval(name, ori_name, mode, has_started);
			nodes.push_back(myval);
		}
		else if (type == SCHEDULE_NODE_TYPE_BOOL_OPERATOR)
		{
			SchedulerNode myval(findBooleanOperator(name));
			nodes.push_back(myval);
		}
		else if (type == SCHEDULE_NODE_TYPE_FLOAT_OPERATOR)
		{
			SchedulerNode myval(findFloatOperator(name));
			nodes.push_back(myval);
		}
		else if (type == SCHEDULE_NODE_TYPE_STRING_OPERATOR)
		{
			SchedulerNode myval(findStringOperator(name));
			nodes.push_back(myval);
		}
		else if (type == SCHEDULE_NODE_TYPE_TIMER_WAIT)
		{
			SchedulerNode myval(wait_seconds);
			nodes.push_back(myval);
		}
		else if (type == SCHEDULE_NODE_TYPE_EXIT)
		{
			SchedulerNode myval;
			nodes.push_back(myval);
		}
	}
	MD.clear();

	// Now that we have the nodes, also set current_node and original_start_node
	if (current_node_name == "undefined") current_node = NULL;
	else setCurrentNode(current_node_name);
	if (original_start_node_name == "undefined") original_start_node = NULL;
	else setOriginalStartNode(original_start_node_name);

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
		if (is_fork)
		{
			SchedulerEdge myval(findNode(inputname), findBooleanVariable(bool_name), findNode(outputname), findNode(outputname_false));
			edges.push_back(myval);
		}
		else
		{
			SchedulerEdge myval(findNode(inputname), findNode(outputname));
			edges.push_back(myval);
		}
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
	if (current_node != NULL)
		MDgeneral.setValue(EMDL_SCHEDULE_GENERAL_CURRENT_NODE, current_node->name);
	else
		MDgeneral.setValue(EMDL_SCHEDULE_GENERAL_CURRENT_NODE, str_aux);
	if (original_start_node != NULL)
		MDgeneral.setValue(EMDL_SCHEDULE_GENERAL_ORIGINAL_START_NODE, original_start_node->name);
	else
		MDgeneral.setValue(EMDL_SCHEDULE_GENERAL_ORIGINAL_START_NODE, str_aux);
	MDgeneral.setValue(EMDL_SCHEDULE_GENERAL_EMAIL, email_address);
	MDgeneral.write(fh);

	if (floats.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_floats");
		for (int i = 0; i < floats.size(); i++)
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_VAR_FLOAT_NAME, floats[i].name);
			MD.setValue(EMDL_SCHEDULE_VAR_FLOAT_VALUE, (floats[i]).value);
			MD.setValue(EMDL_SCHEDULE_VAR_FLOAT_ORI_VALUE, (floats[i]).original_value);
		}
		MD.write(fh);
	}

	if (bools.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_bools");
		for (int i = 0; i < bools.size(); i++)
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_VAR_BOOL_NAME, bools[i].name);
			MD.setValue(EMDL_SCHEDULE_VAR_BOOL_VALUE, bools[i].value);
			MD.setValue(EMDL_SCHEDULE_VAR_BOOL_ORI_VALUE, bools[i].original_value);
		}
		MD.write(fh);
	}
	if (strings.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_strings");
		for (int i = 0; i < strings.size(); i++)
		{
			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_VAR_STRING_NAME, strings[i].name);
			MD.setValue(EMDL_SCHEDULE_VAR_STRING_VALUE, strings[i].value);
			MD.setValue(EMDL_SCHEDULE_VAR_STRING_ORI_VALUE, strings[i].original_value);
		}
		MD.write(fh);
	}


	if (bool_ops.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_boolean_operators");
		for (int i = 0; i < bool_ops.size(); i++)
		{
			std::string type = bool_ops[i].type;

			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_NAME, bool_ops[i].name);
			MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_TYPE, bool_ops[i].type);
			MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_OUTPUT, bool_ops[i].output->name);

			if (type == SCHEDULE_BOOLEAN_OPERATOR_AND ||
				type == SCHEDULE_BOOLEAN_OPERATOR_OR ||
				type == SCHEDULE_BOOLEAN_OPERATOR_NOT)
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_INPUT1, bool_ops[i].input1->name);
			}
			else
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_INPUT1, bool_aux);
			}
			if (type == SCHEDULE_BOOLEAN_OPERATOR_AND ||
				type == SCHEDULE_BOOLEAN_OPERATOR_OR )
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_INPUT2, bool_ops[i].input2->name);
			}
			else
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_INPUT2, bool_aux);
			}

			MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_FLOAT1, bool_ops[i].float1->name);

			if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_CONST ||
				type == SCHEDULE_BOOLEAN_OPERATOR_LT_CONST ||
				type == SCHEDULE_BOOLEAN_OPERATOR_EQ_CONST ||
				type == SCHEDULE_BOOLEAN_OPERATOR_GT_VAR ||
				type == SCHEDULE_BOOLEAN_OPERATOR_LT_VAR ||
				type == SCHEDULE_BOOLEAN_OPERATOR_EQ_VAR)
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_FLOAT1, bool_ops[i].float1->name);
			}
			else
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_FLOAT1, str_aux);
			}

			if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_VAR ||
				type == SCHEDULE_BOOLEAN_OPERATOR_LT_VAR ||
				type == SCHEDULE_BOOLEAN_OPERATOR_EQ_VAR)
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_FLOAT2, bool_ops[i].float2->name);
			}
			else
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_FLOAT2, str_aux);
			}

			if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_CONST ||
				type == SCHEDULE_BOOLEAN_OPERATOR_LT_CONST ||
				type == SCHEDULE_BOOLEAN_OPERATOR_EQ_CONST)
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_CONST, bool_ops[i].myconstant);
			}
			else
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_CONST, float_aux);
			}

			if (type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS)
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_FILE, bool_ops[i].file->name);
			}
			else
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_BOOL_FILE, str_aux);
			}
		}
		MD.write(fh);
	}

	if (float_ops.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_float_operators");
		for (int i = 0; i < float_ops.size(); i++)
		{
			MD.addObject();
			std::string type = float_ops[i].type;

			MD.setValue(EMDL_SCHEDULE_OPERATOR_FLOAT_NAME, float_ops[i].name);
			MD.setValue(EMDL_SCHEDULE_OPERATOR_FLOAT_TYPE, float_ops[i].type);
			MD.setValue(EMDL_SCHEDULE_OPERATOR_FLOAT_OUTPUT, float_ops[i].output->name);

			MD.setValue(EMDL_SCHEDULE_OPERATOR_FLOAT_INPUT1, float_ops[i].input1->name);

			if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_VAR ||
				type == SCHEDULE_FLOAT_OPERATOR_MINUS_VAR ||
				type == SCHEDULE_FLOAT_OPERATOR_MULT_VAR ||
				type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_VAR)
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_FLOAT_INPUT2, float_ops[i].input2->name);
			}
			else
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_FLOAT_INPUT2, str_aux);
			}

			if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_CONST ||
				type == SCHEDULE_FLOAT_OPERATOR_MINUS_CONST ||
				type == SCHEDULE_FLOAT_OPERATOR_MULT_CONST ||
				type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST ||
				type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST_INV)
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_FLOAT_CONSTANT, float_ops[i].myconstant);
			}
			else
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_FLOAT_CONSTANT, float_aux);
			}
		}
		MD.write(fh);
	}

	if (string_ops.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_string_operators");
		for (int i = 0; i < string_ops.size(); i++)
		{
			MD.addObject();

			std::string type = string_ops[i].type;

			MD.setValue(EMDL_SCHEDULE_OPERATOR_STRING_NAME, string_ops[i].name);
			MD.setValue(EMDL_SCHEDULE_OPERATOR_STRING_TYPE, string_ops[i].type);
			MD.setValue(EMDL_SCHEDULE_OPERATOR_STRING_INPUT, string_ops[i].input->name);

			if (type == SCHEDULE_STRING_OPERATOR_COPY_FILE ||
				type == SCHEDULE_STRING_OPERATOR_MOVE_FILE)
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_STRING_OUTPUT, string_ops[i].output->name);
			}
			else
			{
				MD.setValue(EMDL_SCHEDULE_OPERATOR_STRING_OUTPUT, str_aux);
			}
		}
		MD.write(fh);
	}

	if (nodes.size() > 0)
	{
		MetaDataTable MD;
		MD.setName("schedule_nodes");
		for (int i = 0; i < nodes.size(); i++)
		{
			std::string type = nodes[i].type;

			MD.addObject();
			MD.setValue(EMDL_SCHEDULE_NODE_NAME, nodes[i].name);
			MD.setValue(EMDL_SCHEDULE_NODE_TYPE, nodes[i].type);


			if (type == SCHEDULE_NODE_TYPE_EXIT)
			{
				MD.setValue(EMDL_SCHEDULE_NODE_ORI_NAME, nodes[i].original_name);
			}
			else if (type == SCHEDULE_NODE_TYPE_JOB)
			{
				MD.setValue(EMDL_SCHEDULE_NODE_ORI_NAME, nodes[i].original_name);
				MD.setValue(EMDL_SCHEDULE_NODE_JOB_MODE, nodes[i].mode);
				MD.setValue(EMDL_SCHEDULE_NODE_JOB_HAS_STARTED, nodes[i].job_has_started);
			}
			else
			{
				MD.setValue(EMDL_SCHEDULE_NODE_ORI_NAME, str_aux);
				MD.setValue(EMDL_SCHEDULE_NODE_JOB_MODE, str_aux);
				MD.setValue(EMDL_SCHEDULE_NODE_JOB_HAS_STARTED, bool_aux);
			}
			if (type == SCHEDULE_NODE_TYPE_TIMER_WAIT)
			{
				MD.setValue(EMDL_SCHEDULE_NODE_WAIT_TIME, nodes[i].wait_seconds);
			}
			else
			{
				MD.setValue(EMDL_SCHEDULE_NODE_WAIT_TIME, float_aux);
			}

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
			MD.setValue(EMDL_SCHEDULE_EDGE_INPUT, edges[i].inputNode->name);
			MD.setValue(EMDL_SCHEDULE_EDGE_OUTPUT, edges[i].outputNode->name);
			MD.setValue(EMDL_SCHEDULE_EDGE_IS_FORK, edges[i].is_fork);

			if (edges[i].is_fork)
			{
				MD.setValue(EMDL_SCHEDULE_EDGE_OUTPUT_FALSE, edges[i].outputNodeFalse->name);
				MD.setValue(EMDL_SCHEDULE_EDGE_BOOLEAN, edges[i].myBooleanVariable->name);
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

    for (int i = 0; i < floats.size(); i++)
        floats[i].value = floats[i].original_value;

    for (int i = 0; i < bools.size(); i++)
        bools[i].value = bools[i].original_value;

    for (int i = 0; i < strings.size(); i++)
        strings[i].value = strings[i].original_value;

    for (int i = 0; i < nodes.size(); i++)
    {
    	nodes[i].name = nodes[i].original_name;
    	nodes[i].job_has_started = false;
    }

    current_node = NULL;
}

void Schedule::setCurrentNode(std::string _name)
{
    current_node = findNode(_name);
}
void Schedule::setOriginalStartNode(std::string _name)
{
    original_start_node = findNode(_name);
}

bool Schedule::gotoNextNode()
{
    if (current_node == NULL)
    {
        if (original_start_node == NULL)
        {
        	std::cerr << " WARNING: start node was not defined, assuming it is the first node in the list ..." << std::endl;
        	original_start_node = &nodes[0];
        }
    	current_node = original_start_node; // go to first node in the list
        std::cout << " Setting current_node to original_start_node: " << current_node->name << std::endl;
    	return true;
    }

    for (int i = 0; i < edges.size(); i++)
    {
        if (edges[i].inputNode == current_node)
        {
            current_node = edges[i].getOutputNode();
            return (current_node == NULL) ? false : true;
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
		if (current_node->performOperation())
		{
			if (current_node->type == SCHEDULE_NODE_TYPE_EXIT)
			{
				sendEmail("Finished successfully!");
				return false;
			}
		}
		else // this is a job, get its name and options
		{
            std::cout << " Now preparing Job: " << current_node->name << std::endl;
			job_name = current_node->name;
            original_name = current_node->original_name;
            mode = current_node->mode;
            has_started = current_node->job_has_started;
            return true;
        }
    }

    return false;
}

bool Schedule::isBooleanVariable(std::string _name)
{
	for (int i = 0; i < bools.size(); i++)
		if (bools[i].name == _name) return true;
	return false;
}

bool Schedule::isFloatVariable(std::string _name)
{
	for (int i = 0; i < floats.size(); i++)
		if (floats[i].name == _name) return true;
	return false;
}

bool Schedule::isStringVariable(std::string _name)
{
	for (int i = 0; i < strings.size(); i++)
		if (strings[i].name == _name) return true;
	return false;
}

bool Schedule::isNode(std::string _name)
{
	for (int i = 0; i < nodes.size(); i++)
		if (nodes[i].name == _name) return true;
	return false;
}

bool Schedule::isJob(std::string _name)
{
	for (int i = 0; i < nodes.size(); i++)
		if (nodes[i].name == _name) return (nodes[i].type == SCHEDULE_NODE_TYPE_JOB);
	return false;
}


void Schedule::addFloatVariable(std::string name, RFLOAT value)
{
	if (isFloatVariable(name))
		REPORT_ERROR("ERROR: trying to add a float variable with a name that already exists: " + name);

	SchedulerFloatVariable myvar(name, value, value);
	floats.push_back(myvar);
}

void Schedule::addBooleanVariable(std::string name, bool value)
{
	if (isBooleanVariable(name))
		REPORT_ERROR("ERROR: trying to add a boolean variable with a name that already exists: " + name);

	SchedulerBooleanVariable myvar(name, value, value);
	bools.push_back(myvar);
}

void Schedule::addStringVariable(std::string name, FileName value)
{
	if (isStringVariable(name))
		REPORT_ERROR("ERROR: trying to add a string variable with a name that already exists: " + name);

	SchedulerStringVariable myvar(name, value, value);
	strings.push_back(myvar);
}

void Schedule::addOperatorNode(std::string type, std::string input_name, std::string input2_name, std::string output_name)
{
	if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_CONST ||
		type == SCHEDULE_FLOAT_OPERATOR_MINUS_CONST ||
		type == SCHEDULE_FLOAT_OPERATOR_MULT_CONST ||
		type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST ||
		type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST_INV)
	{
		SchedulerFloatOperator myop(type, findFloatVariable(input_name), textToFloat(input2_name), findFloatVariable(output_name));
		float_ops.push_back(myop);
		SchedulerNode mynode(&myop);
		nodes.push_back(mynode);
	}
	else if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_VAR ||
			 type == SCHEDULE_FLOAT_OPERATOR_MINUS_VAR ||
			 type == SCHEDULE_FLOAT_OPERATOR_MULT_VAR ||
			 type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_VAR)
	{
		SchedulerFloatOperator myop(type, findFloatVariable(input_name), findFloatVariable(input2_name), findFloatVariable(output_name));
		float_ops.push_back(myop);
		SchedulerNode mynode(&myop);
		nodes.push_back(mynode);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_NOT)
	{
		SchedulerBooleanOperator myop(type, findBooleanVariable(input_name), findBooleanVariable(output_name));
		bool_ops.push_back(myop);
		SchedulerNode mynode(&myop);
		nodes.push_back(mynode);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_AND || type == SCHEDULE_BOOLEAN_OPERATOR_OR)
	{
		SchedulerBooleanOperator myop(type, findBooleanVariable(input_name), findBooleanVariable(input2_name), findBooleanVariable(output_name));
		bool_ops.push_back(myop);
		SchedulerNode mynode(&myop);
		nodes.push_back(mynode);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_VAR ||
			 type == SCHEDULE_BOOLEAN_OPERATOR_LT_VAR ||
			 type == SCHEDULE_BOOLEAN_OPERATOR_EQ_VAR)
	{
		SchedulerBooleanOperator myop(type, findFloatVariable(input_name), findFloatVariable(input2_name), findBooleanVariable(output_name));
		bool_ops.push_back(myop);
		SchedulerNode mynode(&myop);
		nodes.push_back(mynode);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_CONST ||
		     type == SCHEDULE_BOOLEAN_OPERATOR_LT_CONST ||
		     type == SCHEDULE_BOOLEAN_OPERATOR_EQ_CONST)
	{
		SchedulerBooleanOperator myop(type, findFloatVariable(input_name), textToFloat(input2_name), findBooleanVariable(output_name));
		bool_ops.push_back(myop);
		SchedulerNode mynode(&myop);
		nodes.push_back(mynode);
	}
	else if (type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS)
	{
		SchedulerBooleanOperator myop(type, findStringVariable(input_name), findBooleanVariable(output_name));
		bool_ops.push_back(myop);
		SchedulerNode mynode(&myop);
		nodes.push_back(mynode);
	}
	else if (type == SCHEDULE_STRING_OPERATOR_COPY_FILE || type == SCHEDULE_STRING_OPERATOR_MOVE_FILE )
	{
		SchedulerStringOperator myop(type, findStringVariable(input_name), findStringVariable(output_name));
		string_ops.push_back(myop);
		SchedulerNode mynode(&myop);
		nodes.push_back(mynode);
	}
	else if (type == SCHEDULE_STRING_OPERATOR_TOUCH_FILE || type == SCHEDULE_STRING_OPERATOR_DELETE_FILE )
	{
		SchedulerStringOperator myop(type, findStringVariable(input_name));
		string_ops.push_back(myop);
		SchedulerNode mynode(&myop);
		nodes.push_back(mynode);
	}
	else
		REPORT_ERROR("ERROR: unrecognised operator type: " + type);

}

void Schedule::addJobNode(FileName jobname, std::string mode)
{
	// Check whether the jobname is unique
	if (isNode(jobname))
		REPORT_ERROR("ERROR: trying to add a JobNode that already exists...");

	// Copy the input job.star file into the Schedule directory as jobname + "job.star"
	copy(jobname + "job.star", name + jobname + "job.star");

	// And add this file to the internal schedule_pipeline
	RelionJob job;
	bool dummy;
	job.read(jobname, dummy, true);

	// Now add this job to the local schedule_pipeline
	std::string error_message;
	std::vector<std::string> commands;
	std::string final_command;
	std::string output_name = name + jobname + '/';

	if (!job.getCommands(output_name, commands, final_command, false, schedule_pipeline.job_counter, error_message))
		REPORT_ERROR("ERROR in getting commands for scheduled job: " + error_message);

	int current_job = schedule_pipeline.addJob(job, PROC_SCHEDULED, false, false); // 1st false is do_overwrite, 2nd false is do_write_minipipeline

	if (current_job < 0)
		REPORT_ERROR("ERROR: current job should not be negative now ...");

	SchedulerNode mynode(jobname, jobname, mode);
	nodes.push_back(mynode);
}

void Schedule::addExitNode()
{
	SchedulerNode mynode;
	nodes.push_back(mynode);
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
	SchedulerEdge myval(findNode(inputnode_name), findNode(outputnode_name));
	edges.push_back(myval);
}

void Schedule::addFork(std::string inputnode_name, std::string mybool_name, std::string outputnode_name, std::string outputnode_name_if_false)
{
	SchedulerEdge myval(findNode(inputnode_name), findBooleanVariable(mybool_name), findNode(outputnode_name), findNode(outputnode_name_if_false));
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
	job.read(name + original_job_name, dummy, true); // true means initialise the job

	// Check where this job gets its input from: change names from local scheduler ones to the current pipeline
	int local_process = schedule_pipeline.findProcessByName(name + original_job_name + '/');
	// Loop over all input nodes to this job
	for (int inode = 0; inode <  schedule_pipeline.processList[local_process].inputNodeList.size(); inode++)
	{
		int mynode = schedule_pipeline.processList[local_process].inputNodeList[inode];
		// find from which pipeline_scheduler job this jobs gets its input nodes
		std::cerr << " schedule_pipeline.nodeList[mynode].name= " << schedule_pipeline.nodeList[mynode].name << std::endl;
		int output_from_process =  schedule_pipeline.nodeList[mynode].outputFromProcess;
		if (output_from_process < 0)
			REPORT_ERROR("ERROR: cannot find outputProcess of node: " +  schedule_pipeline.nodeList[mynode].name);
		// Get the original name in the pipeline_scheduler of this job
		FileName my_ori_name = schedule_pipeline.processList[output_from_process].name;
		// Remove leading directory and tailing slash to get the process name in the pipeline_scheduler
		FileName my_process_name = (my_ori_name.afterFirstOf(name)).beforeLastOf("/");
		// find that process in the nodes, and get its current name
		std::string my_current_name =findNodeByOriginalName(my_process_name)->name;
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
	ori_job.read(name + original_job_name, dummy, true);

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

				my_value = (findBooleanVariable(mystring)->value) ? "Yes" : "No";
			}
			else if (isFloatVariable(mystring))
			{
				if (job.joboptions[it->first].joboption_type == JOBOPTION_BOOLEAN)
					REPORT_ERROR(" ERROR: trying to set FloatVariable: " + mystring + " into a boolean option: " + it->first);

				my_value = floatToString(findFloatVariable(mystring)->value);
			}
			else if (isStringVariable(mystring))
			{
				if (job.joboptions[it->first].joboption_type == JOBOPTION_BOOLEAN)
					REPORT_ERROR(" ERROR: trying to set StringVariable: " + mystring + " into a boolean option: " + it->first);
				if (job.joboptions[it->first].joboption_type == JOBOPTION_SLIDER)
					REPORT_ERROR(" ERROR: trying to set StringVariable: " + mystring + " into a slider option: " + it->first);

				my_value = findStringVariable(mystring)->value;
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

        	// Set the name of the current node now
        	current_node->name = pipeline.processList[current_job].name;

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
			current_node->job_has_started = true;
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
    	touch(processList[current_job].name + RELION_JOB_ABORT_NOW);
        std::cout << " Marking job " << processList[current_job].name << " for abortion; schedule should abort too ... " << std::endl;
    }


    while ()
    {
        int current_job = findProcessByName(job_name);
        if (processList[current_job].status == PROC_RUNNING)
        {
            has_done = true;
        	touch(processList[current_job].name + RELION_JOB_ABORT_NOW);
            std::cout << " Marking job " << processList[current_job].name << " for abortion; schedule should abort too ... " << std::endl;
        }
    }

    if (!has_done)
    	std::cout << " This schedule seems to have finished already ..." << std::endl;
	*/

}

