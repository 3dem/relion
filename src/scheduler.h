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

#ifndef SCHEDULER_H_
#define SCHEDULER_H_

#include "src/time.h"
#include "src/pipeliner.h"

class SchedulerFloatVariable
{
	public:
	std::string name;
	RFLOAT value, original_value;

	SchedulerFloatVariable(std::string _myname, RFLOAT _my_value, RFLOAT _my_ori_value)
	{
		name = _myname;
		value = _my_value;
		original_value = _my_ori_value;
	}
};

class SchedulerBooleanVariable
{
	public:
	std::string name;
	bool value, original_value;

	SchedulerBooleanVariable(std::string _myname, bool _my_value, bool _my_ori_value)
	{
		name = _myname;
		value = _my_value;
		original_value = _my_ori_value;
	}
};

class SchedulerStringVariable
{
	public:
	std::string name;
	FileName value, original_value;

	SchedulerStringVariable(std::string _myname, FileName _my_value, FileName _my_ori_value)
	{
		name = _myname;
		value = _my_value;
		original_value = _my_ori_value;
	}
};

#define SCHEDULE_BOOLEAN_OPERATOR_AND "bool_op_and"
#define SCHEDULE_BOOLEAN_OPERATOR_OR  "bool_op_or"
#define SCHEDULE_BOOLEAN_OPERATOR_NOT "bool_op_not"
#define SCHEDULE_BOOLEAN_OPERATOR_GT_VAR "bool_op_gt_var"
#define SCHEDULE_BOOLEAN_OPERATOR_LT_VAR "bool_op_lt_var"
#define SCHEDULE_BOOLEAN_OPERATOR_EQ_VAR "bool_op_eq_var"
#define SCHEDULE_BOOLEAN_OPERATOR_GT_CONST "bool_op_gt_const"
#define SCHEDULE_BOOLEAN_OPERATOR_LT_CONST "bool_op_lt_const"
#define SCHEDULE_BOOLEAN_OPERATOR_EQ_CONST "bool_op_eq_const"
#define SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS "bool_op_file_exists"

// A class that performs boolean logic
class SchedulerBooleanOperator
{
	public:
	std::string type, name;
	SchedulerBooleanVariable *input1, *input2;
	SchedulerFloatVariable *float1, *float2;
	RFLOAT myconstant;
	SchedulerStringVariable *file;
	SchedulerBooleanVariable *output;

	public:

	void performOperation() const
	{
		if (type == SCHEDULE_BOOLEAN_OPERATOR_AND)
		{
			output->value = (input1->value && input2->value);
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_OR)
		{
			output->value = (input1->value || input2->value);
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_NOT)
		{
			output->value = (!(input1->value));
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_VAR)
		{
			output->value = (float1->value > float2->value);
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_LT_VAR)
		{
			output->value = (float1->value < float2->value);
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_EQ_VAR)
		{
			output->value = (fabs(float1->value - float2->value) < 1E-8);
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_CONST)
		{
			output->value = (float1->value > myconstant);
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_LT_CONST)
		{
			output->value = (float1->value < myconstant);
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_EQ_CONST)
		{
			output->value = (fabs(float1->value - myconstant) < 1E-8);
		}
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS)
		{
			output->value = (exists(file->value));
		}
		else
			REPORT_ERROR("ERROR: undefined boolean operator type");

		if (output->value)
			std::cout << " Setting Boolean variable " << output->name << " to: True " << std::endl;
		else
			std::cout << " Setting Boolean variable " << output->name << " to: False " << std::endl;

	}

	SchedulerBooleanOperator(std::string _type, SchedulerBooleanVariable *_input1,
			SchedulerBooleanVariable*_input2, SchedulerBooleanVariable *_output)
	{

		type = _type;
		input1 = _input1;
		input2 = _input2;
		float1 = NULL;
		float2 = NULL;
		myconstant = 0.;
		file = NULL;
		output = _output;

		if (type == SCHEDULE_BOOLEAN_OPERATOR_AND)
			name = input1->name + "_AND_" + input2->name;
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_OR)
			name = input1->name + "_AND_" + input2->name;
		else
			REPORT_ERROR("BUG: incorrect boolean operator type ....");
	}

	SchedulerBooleanOperator(std::string _type, SchedulerBooleanVariable *_input1, SchedulerBooleanVariable *_output)
	{
		type = _type;
		input1 = _input1;
		input2 = NULL;
		float1 = NULL;
		float2 = NULL;
		myconstant = 0.;
		file = NULL;
		output = _output;

		if (type == SCHEDULE_BOOLEAN_OPERATOR_NOT)
			name = "NOT_" + input1->name;
		else
			REPORT_ERROR("BUG: incorrect boolean operator type ....");
	}

	SchedulerBooleanOperator(std::string _type, SchedulerFloatVariable *_float1, SchedulerFloatVariable *_float2, SchedulerBooleanVariable *_output)
	{
		type = _type;
		input1 = NULL;
		input2 = NULL;
		float1 = _float1;
		float2 = _float2;
		myconstant = 0.;
		file = NULL;
		output = _output;

		if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_VAR)
			name = float1->name + "_GT_" + float2->name;
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_LT_VAR)
			name = float1->name + "_LT_" + float2->name;
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_EQ_VAR)
			name = float1->name + "_EQ_" + float2->name;
		else
			REPORT_ERROR("BUG: incorrect boolean operator type ....");

	}

	SchedulerBooleanOperator(std::string _type, SchedulerFloatVariable *_float1, RFLOAT _constant, SchedulerBooleanVariable *_output)
	{
		type = _type;
		input1 = NULL;
		input2 = NULL;
		float1 = _float1;
		float2 = NULL;
		myconstant = _constant;
		file = NULL;
		output = _output;

		if (type == SCHEDULE_BOOLEAN_OPERATOR_GT_CONST)
			name = float1->name + "_GT_" + floatToString(myconstant);
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_LT_CONST)
			name = float1->name + "_LT_" + floatToString(myconstant);
		else if (type == SCHEDULE_BOOLEAN_OPERATOR_EQ_CONST)
			name = float1->name + "_EQ_" + floatToString(myconstant);
		else
			REPORT_ERROR("BUG: incorrect boolean operator type ....");

	}

	SchedulerBooleanOperator(std::string _type, SchedulerStringVariable *_string, SchedulerBooleanVariable *_output)
	{
		type = _type;
		input1 = NULL;
		input2 = NULL;
		float1 = NULL;
		float2 = NULL;
		myconstant = 0.;
		file = _string;
		output = _output;

		if (type == SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS)
			name = "EXISTS_" + file->name;
		else
			REPORT_ERROR("BUG: incorrect boolean operator type ....");

	}

};

#define SCHEDULE_FLOAT_OPERATOR_PLUS_VAR "float_op_plus_float"
#define SCHEDULE_FLOAT_OPERATOR_MINUS_VAR "float_op_minus_float"
#define SCHEDULE_FLOAT_OPERATOR_MULT_VAR "float_op_mult_float"
#define SCHEDULE_FLOAT_OPERATOR_DIVIDE_VAR "float_op_divide_float"
#define SCHEDULE_FLOAT_OPERATOR_PLUS_CONST "float_op_plus_const"
#define SCHEDULE_FLOAT_OPERATOR_MINUS_CONST "float_op_minus_const"
#define SCHEDULE_FLOAT_OPERATOR_MULT_CONST "float_op_mult_const"
#define SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST "float_op_div_by_const"
#define SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST_INV "float_op_div_const_by"

// A class that performs basic operations on float variables
class SchedulerFloatOperator
{
	public:
	SchedulerFloatVariable *input1, *input2;
	SchedulerFloatVariable *output;
	RFLOAT myconstant;
	std::string name, type;

	public:

	void performOperation() const
	{
		if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_VAR)
		{
			output->value = input1->value + input2->value;
		}
		else if (type ==SCHEDULE_FLOAT_OPERATOR_MINUS_VAR)
		{
			output->value = input1->value - input2->value;
		}
		else if (type == SCHEDULE_FLOAT_OPERATOR_MULT_VAR)
		{
			output->value = input1->value * input2->value;
		}
		else if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_VAR)
		{
			output->value = input1->value / input2->value;
		}
		else if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_CONST)
		{
			output->value = input1->value + myconstant;
		}
		else if (type == SCHEDULE_FLOAT_OPERATOR_MINUS_CONST)
		{
			output->value = input1->value - myconstant;
		}
		else if (type == SCHEDULE_FLOAT_OPERATOR_MULT_CONST)
		{
			output->value = input1->value * myconstant;
		}
		else if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST)
		{
			output->value = input1->value / myconstant;
		}
		else if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST_INV)
		{
			output->value = myconstant / input1->value;
		}
		else
		{
			REPORT_ERROR("BUG: unrecognised float operator type... ");
		}

		std::cout << " Setting Float variable " << output->name << " to: " << output->value << std::endl;
	}

	SchedulerFloatOperator(std::string _type, SchedulerFloatVariable *_input1,
			SchedulerFloatVariable *_input2,
			SchedulerFloatVariable *_output)
	{
		type = _type;
		input1 = _input1;
		input2 = _input2;
		output = _output;
		myconstant = 0.;

		if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_VAR)
			name = input1->name + "_PLUS_" + input2->name;
		else if (type == SCHEDULE_FLOAT_OPERATOR_MINUS_VAR)
			name = input1->name + "_MINUS_" + input2->name;
		else if (type == SCHEDULE_FLOAT_OPERATOR_MULT_VAR)
			name = input1->name + "_MULT_" + input2->name;
		else if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_VAR)
			name = input1->name + "_DIV_" + input2->name;
		else
			REPORT_ERROR("BUG: incorrect float operator type for two variables ...");
	}

	SchedulerFloatOperator(std::string _type, SchedulerFloatVariable *_input1,
			RFLOAT _constant, SchedulerFloatVariable *_output)
	{
		type = _type;
		myconstant = _constant;
		input1 = _input1;
		input2 = NULL;
		output = _output;

		if (type == SCHEDULE_FLOAT_OPERATOR_PLUS_CONST)
			name = input1->name + "_PLUS_" + floatToString(myconstant);
		else if (type == SCHEDULE_FLOAT_OPERATOR_MINUS_CONST)
			name = input1->name + "_MINUS_" + floatToString(myconstant);
		else if (type == SCHEDULE_FLOAT_OPERATOR_MULT_CONST)
			name = input1->name + "_MULT_" + floatToString(myconstant);
		else if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST)
			name = input1->name + "_DIV_" + floatToString(myconstant);
		else if (type == SCHEDULE_FLOAT_OPERATOR_DIVIDE_CONST_INV)
			name = floatToString(myconstant) + "_DIV_" + input1->name;
		else
			REPORT_ERROR("BUG: incorrect float operator type for constant ...");
	}

};

#define SCHEDULE_STRING_OPERATOR_TOUCH_FILE "string_op_touch_file"
#define SCHEDULE_STRING_OPERATOR_COPY_FILE "string_op_copy_file"
#define SCHEDULE_STRING_OPERATOR_MOVE_FILE "string_op_move_file"
#define SCHEDULE_STRING_OPERATOR_DELETE_FILE "string_op_delete_file"

// A class that performs basic operations on string variables
class SchedulerStringOperator
{
	public:
	SchedulerStringVariable *input;;
	SchedulerStringVariable *output;
	std::string name, type;

	public:

	void performOperation() const
	{
		if (type == SCHEDULE_STRING_OPERATOR_TOUCH_FILE)
		{
			std::cout << " Touching file " << input->value << std::endl;
			touch(input->value);
		}
		if (type == SCHEDULE_STRING_OPERATOR_COPY_FILE)
		{
			std::cout << " Copying file " << input->value << " to " << output->value << std::endl;
			copy(input->value, output->value);
		}
		if (type == SCHEDULE_STRING_OPERATOR_MOVE_FILE)
		{
			std::cout << " Moving file from " << input->value << " to " << output->value << std::endl;
			move(input->value, output->value);
		}
		if (type == SCHEDULE_STRING_OPERATOR_DELETE_FILE)
		{
			std::cout << " Deleting file " << input->value << std::endl;
			delete((input->value).c_str());
		}
		else
		{
			REPORT_ERROR("BUG: unrecognised string operator type... ");
		}
	}

	SchedulerStringOperator(std::string _type, SchedulerStringVariable *_input,
			SchedulerStringVariable *_output)
	{
		type = _type;
		input = _input;
		output = _output;

		if (type == SCHEDULE_STRING_OPERATOR_COPY_FILE)
			name = "COPY_" + input->name + "_" + output->name;
		else if (type == SCHEDULE_STRING_OPERATOR_MOVE_FILE)
			name = "MOVE_" + input->name + "_" + output->name;
		else
			REPORT_ERROR("BUG: incorrect string operator type for two string variables ...");
	}

	SchedulerStringOperator(std::string _type, SchedulerStringVariable *_input)
	{
		type = _type;
		input = _input;

		if (type == SCHEDULE_STRING_OPERATOR_TOUCH_FILE)
			name = "TOUCH_" + input->name;
		else if (type == SCHEDULE_STRING_OPERATOR_DELETE_FILE)
			name = "DELETE_" + input->name;
		else
			REPORT_ERROR("BUG: incorrect float operator type for one string variable ...");

	}
};


#define SCHEDULE_NODE_TYPE_JOB "job"
#define SCHEDULE_NODE_TYPE_BOOL_OPERATOR "bool_op"
#define SCHEDULE_NODE_TYPE_FLOAT_OPERATOR "float_op"
#define SCHEDULE_NODE_TYPE_STRING_OPERATOR "string_op"
#define SCHEDULE_NODE_TYPE_TIMER_WAIT "wait"
#define SCHEDULE_NODE_TYPE_EXIT "exit"

#define SCHEDULE_NODE_JOB_MODE_NEW "new"
#define SCHEDULE_NODE_JOB_MODE_CONTINUE "continue"
#define SCHEDULE_NODE_JOB_MODE_OVERWRITE "overwrite"

class SchedulerNode
{
	public:
	std::string type, name, original_name, mode;
	bool job_has_started;
	SchedulerBooleanOperator *myBooleanOperator;
	SchedulerFloatOperator *myFloatOperator;
	SchedulerStringOperator *myStringOperator;
	RFLOAT wait_seconds;

	public:

	SchedulerNode(std::string _job_name, std::string _ori_name, std::string _mode, bool _has_started = false)
	{
		type = SCHEDULE_NODE_TYPE_JOB;
		name = _job_name;
		original_name = _ori_name;
		mode = _mode;
		myBooleanOperator = NULL;
		myFloatOperator = NULL;
		myStringOperator = NULL;
		wait_seconds = 0.;
		job_has_started = _has_started;
	}

	SchedulerNode(SchedulerBooleanOperator* _mybooloperator)
	{
		type = SCHEDULE_NODE_TYPE_BOOL_OPERATOR;
		name = original_name = _mybooloperator->name;
		myBooleanOperator = _mybooloperator;
		myFloatOperator = NULL;
		myStringOperator = NULL;
		wait_seconds = 0.;
		job_has_started = false;
	}

	SchedulerNode(SchedulerFloatOperator* _myfloatoperator)
	{
		type = SCHEDULE_NODE_TYPE_FLOAT_OPERATOR;
		name = original_name = _myfloatoperator->name;
		myFloatOperator = _myfloatoperator;
		myBooleanOperator = NULL;
		myStringOperator = NULL;
		wait_seconds = 0.;
		job_has_started = false;
	}

	SchedulerNode(SchedulerStringOperator* _mystringoperator)
	{
		type = SCHEDULE_NODE_TYPE_STRING_OPERATOR;
		name = original_name = _mystringoperator->name;
		mode = "undefined";
		myStringOperator = _mystringoperator;
		myBooleanOperator = NULL;
		myFloatOperator = NULL;
		wait_seconds = 0.;
		job_has_started = false;
	}

	SchedulerNode(RFLOAT _wait_seconds)
	{
		type = SCHEDULE_NODE_TYPE_TIMER_WAIT;
		name = original_name = "WAIT_" + floatToString( _wait_seconds);
		wait_seconds =  _wait_seconds;
		mode = "undefined";
		myFloatOperator = NULL;
		myBooleanOperator = NULL;
		myStringOperator = NULL;
		job_has_started = false;
	}

	SchedulerNode()
	{
		type = SCHEDULE_NODE_TYPE_EXIT;
		name = original_name = "exit";
		wait_seconds =  0.;
		mode = "undefined";
		myFloatOperator = NULL;
		myBooleanOperator = NULL;
		myStringOperator = NULL;
		job_has_started = false;
	}

	// Perform operation and return TRUE if not a JOB; just return FALSE if a JOB
	bool performOperation();
};


// A class that defines the edges between a graph that defines execution order, where the nodes are individual JOB instances
// An edge can also be a fork, where the output is controlled through a boolean variable
class SchedulerEdge
{
	public:
	SchedulerNode *inputNode, *outputNode, *outputNodeFalse;
	SchedulerBooleanVariable *myBooleanVariable;
	bool is_fork;

	SchedulerNode* getOutputNode() const
	{
		if (is_fork) return (myBooleanVariable->value) ?  outputNode : outputNodeFalse;
		else return outputNode;
	}

	SchedulerEdge(SchedulerNode *_input, SchedulerNode *_output)
	{
		inputNode = _input;
		outputNode= _output;
		is_fork = false;
		outputNodeFalse = NULL;
		myBooleanVariable = NULL;
	}

	SchedulerEdge(SchedulerNode *_input, SchedulerBooleanVariable *_mybool, SchedulerNode *_output = NULL,
			SchedulerNode* _output_if_false = NULL)
	{
		inputNode = _input;
		outputNode= _output;
		is_fork = true;
		outputNodeFalse = _output_if_false;
		myBooleanVariable = _mybool;
	}

};

class Schedule
{
	public:

	SchedulerNode *current_node, *original_start_node;
	std::string name, email_address;

	std::vector<SchedulerBooleanVariable> bools;
	std::vector<SchedulerFloatVariable> floats;
	std::vector<SchedulerStringVariable> strings;

	std::vector<SchedulerBooleanOperator> bool_ops;
	std::vector<SchedulerFloatOperator> float_ops;
	std::vector<SchedulerStringOperator> string_ops;

	std::vector<SchedulerNode> nodes;
	std::vector<SchedulerEdge> edges;

	PipeLine schedule_pipeline;

public:

	Schedule()
	{
		clear();
	}

	void clear()
	{
		current_node = NULL;
		original_start_node = NULL;
		name = "undefined";
		email_address = "undefined";
		bools.clear();
		floats.clear();
		strings.clear();
		nodes.clear();
		edges.clear();
		bool_ops.clear();
		float_ops.clear();
		string_ops.clear();
	}

	void setName(std::string _name)
	{
		name = _name;
		schedule_pipeline.setName(_name + "schedule");
	}

	void read(FileName fn = "");

	void write(FileName fn = "");

    void reset();

    void setCurrentNode(std::string _name);
    void setOriginalStartNode(std::string _name);

    bool gotoNextNode();

    bool gotoNextJob(FileName &job_name, FileName &original_name, std::string &mode, bool &has_started);

    bool isBooleanVariable(std::string name);
	bool isFloatVariable(std::string name);
	bool isStringVariable(std::string name);
	bool isNode(std::string name);
	bool isJob(std::string name);

	SchedulerBooleanVariable* findBooleanVariable(std::string name);
	SchedulerFloatVariable* findFloatVariable(std::string name);
	SchedulerStringVariable* findStringVariable(std::string name);
	SchedulerBooleanOperator* findBooleanOperator(std::string name);
	SchedulerFloatOperator* findFloatOperator(std::string name);
	SchedulerStringOperator* findStringOperator(std::string name);
	SchedulerNode* findNode(std::string name);
	SchedulerNode* findNodeByOriginalName(std::string name);

    // Build a new Schedule
    // Add variables
    void addFloatVariable(std::string name, RFLOAT value);
    void addBooleanVariable(std::string name, bool value);
    void addStringVariable(std::string name, FileName value);

    // Add operators (of any kind), also adds its corresponding node
    void addOperatorNode(std::string type, std::string input_name, std::string input2_name, std::string output_name);

    // Add a new job, also adds its corresponding node
    void addJobNode(RelionJob &myjob, std::string jobname, std::string mode);

    void addExitNode();

    void sendEmail(std::string message);

    // Add edges and forks in between the nodes
    void addEdge(std::string inputnode_name, std::string outputnode_name);
    void addFork(std::string inputnode_name, std::string mybool_name, std::string outputnode_name, std::string outputnode_name_if_false );

    // Test integrity of the Schedule. Warn for unused variables, nodes, etc.
    bool isValid();

    // This function fixes the dependency of newly generated jobs, as determined by the pipeline_schedule
    RelionJob copyNewJobFromSchedulePipeline(FileName original_job_name);

    // Modify a job to set variables from the Scheduler
    void setVariablesInJob(RelionJob &job, FileName original_job_name);

    // Run the Schedule
    void run(PipeLine &pipeline);

    // Abort a running schedule
    void abort();

};


#endif /* SCHEDULER_H_ */
