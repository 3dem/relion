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
#include "src/jaz/obs_model.h"
#define SCHEDULE_HAS_CHANGED ".schedule_has_changed";

class SchedulerFloatVariable
{
	public:
	RFLOAT value, original_value;

	SchedulerFloatVariable() {};

	SchedulerFloatVariable(RFLOAT _value, RFLOAT _original_value)
	{
		value = _value;
		original_value = _original_value;
	}
};

class SchedulerBooleanVariable
{
	public:
	bool value, original_value;

	SchedulerBooleanVariable() {};

	SchedulerBooleanVariable(bool _value, bool _original_value)
	{
		value = _value;
		original_value = _original_value;
	}
};

class SchedulerStringVariable
{
	public:
	FileName value, original_value;

	SchedulerStringVariable() {};

	SchedulerStringVariable(FileName _value, FileName _original_value)
	{
		value = (_value == "") ? "undefined" : _value;
		original_value = (_original_value == "") ? "undefined" : _original_value;
	}
};

bool isBooleanVariable(std::string name);
bool isFloatVariable(std::string name);
bool isStringVariable(std::string name);
bool isScheduleOperator(std::string name);

#define SCHEDULE_BOOLEAN_OPERATOR_AND "bool_op_and"
#define SCHEDULE_BOOLEAN_OPERATOR_OR  "bool_op_or"
#define SCHEDULE_BOOLEAN_OPERATOR_NOT "bool_op_not"
#define SCHEDULE_BOOLEAN_OPERATOR_GT "bool_op_gt"
#define SCHEDULE_BOOLEAN_OPERATOR_LT "bool_op_lt"
#define SCHEDULE_BOOLEAN_OPERATOR_EQ "bool_op_eq"
#define SCHEDULE_BOOLEAN_OPERATOR_FILE_EXISTS "bool_op_file_exists"
#define SCHEDULE_BOOLEAN_OPERATOR_READ_STAR "bool_op_read_star"
#define SCHEDULE_FLOAT_OPERATOR_PLUS "float_op_plus"
#define SCHEDULE_FLOAT_OPERATOR_MINUS "float_op_minus"
#define SCHEDULE_FLOAT_OPERATOR_MULT "float_op_mult"
#define SCHEDULE_FLOAT_OPERATOR_DIVIDE "float_op_divide"
#define SCHEDULE_FLOAT_OPERATOR_INVDIV "float_op_invdiv"
#define SCHEDULE_FLOAT_OPERATOR_COUNT_IMAGES "float_op_count_images"
#define SCHEDULE_FLOAT_OPERATOR_COUNT_WORDS "float_op_count_words"
#define SCHEDULE_FLOAT_OPERATOR_READ_STAR "float_op_star"
#define SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MAX "float_op_star_table_max"
#define SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MIN "float_op_star_table_min"
#define SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_AVG "float_op_star_table_avg"
#define SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MAX_IDX "float_op_star_table_max_idx"
#define SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_MIN_IDX "float_op_star_table_min_idx"
#define SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_IDX "float_op_star_table_idx"
#define SCHEDULE_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX "float_op_star_table_sort_idx"
#define SCHEDULE_STRING_OPERATOR_JOIN "string_op_join"
#define SCHEDULE_STRING_OPERATOR_BEFORE_FIRST "string_op_before_first"
#define SCHEDULE_STRING_OPERATOR_AFTER_FIRST "string_op_after_first"
#define SCHEDULE_STRING_OPERATOR_BEFORE_LAST "string_op_before_last"
#define SCHEDULE_STRING_OPERATOR_AFTER_LAST "string_op_after_last"
#define SCHEDULE_STRING_OPERATOR_TOUCH_FILE "string_op_touch_file"
#define SCHEDULE_STRING_OPERATOR_COPY_FILE "string_op_copy_file"
#define SCHEDULE_STRING_OPERATOR_MOVE_FILE "string_op_move_file"
#define SCHEDULE_STRING_OPERATOR_DELETE_FILE "string_op_delete_file"
#define SCHEDULE_STRING_OPERATOR_READ_STAR "string_op_read_star"
#define SCHEDULE_STRING_OPERATOR_GLOB "string_op_glob"
#define SCHEDULE_STRING_OPERATOR_NTH_WORD "string_op_nth_word"
#define SCHEDULE_WAIT_OPERATOR_SINCE_LAST_TIME "wait_since_last_time"
#define SCHEDULE_EXIT_OPERATOR "exit"

// A class that performs operators on variables
class SchedulerOperator
{
	public:
	std::string type, input1, input2, output;

	public:


	SchedulerOperator() {};

	SchedulerOperator(std::string _type, std::string _input1="undefined", std::string _input2="undefined", std::string _output="undefined");

	std::string initialise(std::string _type, std::string _input1="undefined", std::string _input2="undefined", std::string _output="undefined");

	// Generate a meaningful current_name for the operator
	std::string getName();

	// Read a specific value from a STAR file
	void readFromStarFile() const;

	bool performOperation() const;

};

#define SCHEDULE_NODE_JOB_MODE_NEW "new"
#define SCHEDULE_NODE_JOB_MODE_CONTINUE "continue"
#define SCHEDULE_NODE_JOB_MODE_OVERWRITE "overwrite"

class SchedulerJob
{
	public:
	std::string current_name, mode;
	bool job_has_started;

	public:

	SchedulerJob() {};

	SchedulerJob(std::string _name, std::string _mode, bool _has_started = false)
	{
		current_name = _name;
		mode = _mode;
		job_has_started = _has_started;
	}

	// Perform operation and return TRUE if not a JOB; just return FALSE if a JOB
	bool performOperation();
};


// A class that defines the edges between a graph that defines execution order, where the nodes are individual JOB instances
// An edge can also be a fork, where the output is controlled through a boolean variable
class SchedulerEdge
{
	public:
	std::string inputNode, outputNode, outputNodeTrue;
	std::string myBooleanVariable;
	bool is_fork;

	std::string getOutputNode() const;

	SchedulerEdge(std::string _input, std::string _output, bool _is_fork, std::string _mybool, std::string _output_if_true)
	{
		inputNode = _input;
		outputNode= _output;
		is_fork = _is_fork;
		outputNodeTrue = _output_if_true;
		myBooleanVariable = _mybool;
	}

	SchedulerEdge(std::string _input, std::string _output)
	{
		inputNode = _input;
		outputNode= _output;
		is_fork = false;
		outputNodeTrue = "undefined";
		myBooleanVariable = "undefined";
	}

};


class Schedule
{

public:

	std::string name, current_node, original_start_node, email_address;
	bool do_read_only;
	int verb;

	std::map<std::string, SchedulerJob> jobs;
	std::vector<SchedulerEdge> edges;

	PipeLine schedule_pipeline;

public:

	Schedule()
	{
		clear();
	}

	void clear();

	void setName(std::string _name)
	{
		name = _name;
		schedule_pipeline.setName(_name + "schedule");
	}


	void read(bool do_lock = false, FileName fn = "");

	bool isWriteLocked();
	void write(bool do_lock = false, FileName fn = "");

    void reset();

    void setCurrentNode(std::string _name);
    void setOriginalStartNode(std::string _name);

	bool isNode(std::string name);
	bool isJob(std::string name);
	bool isOperator(std::string name);

	std::string findJobByCurrentName(std::string name);

    // Get/set Variables and Operators(scheduler_floats is only visible in this file!)
    float getFloatVariableValue(std::string name);
    float getFloatOriginalVariableValue(std::string name);
    void setFloatVariableValue(std::string name, RFLOAT val);
    void setFloatOriginalVariableValue(std::string name, RFLOAT val);

    bool getBooleanVariableValue(std::string name);
    bool getBooleanOriginalVariableValue(std::string name);
    void setBooleanVariableValue(std::string name, bool val);
    void setBooleanOriginalVariableValue(std::string name, bool val);

    std::string getStringVariableValue(std::string name);
    std::string getStringOriginalVariableValue(std::string name);
    void setStringVariableValue(std::string name, std::string val);
    void setStringOriginalVariableValue(std::string name, std::string val);

    std::string getVariableValueAsString(std::string name);

    std::string getOperatorName(std::string type, std::string input1, std::string input2, std::string output)
    {
    	SchedulerOperator op(type, input1, input2, output);
    	return op.getName();
    }
    void setOperatorParameters(std::string name, std::string type, std::string input1, std::string input2, std::string output);
    void getOperatorParameters(std::string name, std::string &type, std::string &input1, std::string &input2, std::string &output);

    // Get vectors with current Variables / Operators
    std::map<std::string, SchedulerFloatVariable> getCurrentFloatVariables();
    std::map<std::string, SchedulerBooleanVariable> getCurrentBooleanVariables();
    std::map<std::string, SchedulerStringVariable> getCurrentStringVariables();
    std::map<std::string, SchedulerOperator> getCurrentOperators();

    // Get/set operators

    // Add variables
    void setVariable(std::string name, FileName value); // (Add new one if exists, otherwise set value)
	void addFloatVariable(std::string name, RFLOAT value);
    void addBooleanVariable(std::string name, bool value);
    void addStringVariable(std::string name, FileName value);

    // Add operators (of any kind), also adds its corresponding node
    SchedulerOperator initialiseOperator(std::string type, std::string input_name, std::string input2_name,
    		std::string output_name, std::string &error_message);
    void addOperator(SchedulerOperator &op);

    // Add a new job, also adds its corresponding node
    void addJob(RelionJob &myjob, std::string jobname, std::string mode);

    void addExitNode();

    // Remove variables/operators/jobs
    void removeVariable(std::string name);
    void removeEdgesWithThisInputOutputOrBoolean(std::string name);
    void removeOperator(std::string name);
    void removeOperatorsWithThisInputOrOutput(std::string name);
    void removeJob(std::string name);
    void removeEdge(int idx);


    // Rename this scheduler into a new directory
    void copy(FileName newname);

    // Send an email
    void sendEmail(std::string message);

    // Add edges and forks in between the nodes
    void addEdge(std::string inputnode_name, std::string outputnode_name);
    void addFork(std::string inputnode_name, std::string mybool_name, std::string outputnode_name, std::string outputnode_name_if_false );

    // Test integrity of the Schedule. Warn for unused variables, nodes, etc.
    bool isValid();

    std::string getNextNode();
    std::string getPreviousNode();

    bool gotoNextNode();
    bool gotoNextJob();


    // Modify a job to set variables and input nodes from the Scheduler
    void setVariablesInJob(RelionJob &job, FileName original_job_name, bool &needs_a_restart);

    // Run the Schedule
    void run(PipeLine &pipeline);

    // Abort a running schedule
    void abort();

};


#endif /* SCHEDULER_H_ */
