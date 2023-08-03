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
#include "src/pipeline_jobs.h"
#include <src/jaz/single_particle/obs_model.h>
//#define SCHEME_HAS_CHANGED ".scheme_has_changed";

class SchemerFloatVariable
{
	public:
	RFLOAT value, original_value;

	SchemerFloatVariable() {};

	SchemerFloatVariable(RFLOAT _value, RFLOAT _original_value)
	{
		value = _value;
		original_value = _original_value;
	}
};

class SchemerBooleanVariable
{
	public:
	bool value, original_value;

	SchemerBooleanVariable() {};

	SchemerBooleanVariable(bool _value, bool _original_value)
	{
		value = _value;
		original_value = _original_value;
	}
};

class SchemerStringVariable
{
	public:
	FileName value, original_value;

	SchemerStringVariable() {};

	SchemerStringVariable(FileName _value, FileName _original_value)
	{
		value = _value;
		original_value = _original_value;
	}
};

bool isBooleanVariable(std::string name);
bool isFloatVariable(std::string name);
bool isStringVariable(std::string name);
bool isSchemeOperator(std::string name);

#define SCHEME_BOOLEAN_OPERATOR_SET "bool=set"
#define SCHEME_BOOLEAN_OPERATOR_AND "bool=and"
#define SCHEME_BOOLEAN_OPERATOR_OR  "bool=or"
#define SCHEME_BOOLEAN_OPERATOR_NOT "bool=not"
#define SCHEME_BOOLEAN_OPERATOR_GT "bool=gt"
#define SCHEME_BOOLEAN_OPERATOR_LT "bool=lt"
#define SCHEME_BOOLEAN_OPERATOR_GE "bool=ge"
#define SCHEME_BOOLEAN_OPERATOR_LE "bool=le"
#define SCHEME_BOOLEAN_OPERATOR_EQ "bool=eq"
#define SCHEME_BOOLEAN_OPERATOR_FILE_EXISTS "bool=file_exists"
#define SCHEME_BOOLEAN_OPERATOR_READ_STAR "bool=read_star"
#define SCHEME_FLOAT_OPERATOR_SET "float=set"
#define SCHEME_FLOAT_OPERATOR_PLUS "float=plus"
#define SCHEME_FLOAT_OPERATOR_MINUS "float=minus"
#define SCHEME_FLOAT_OPERATOR_MULT "float=mult"
#define SCHEME_FLOAT_OPERATOR_DIVIDE "float=divide"
#define SCHEME_FLOAT_OPERATOR_ROUND "float=round"
#define SCHEME_FLOAT_OPERATOR_COUNT_IMAGES "float=count_images"
#define SCHEME_FLOAT_OPERATOR_COUNT_WORDS "float=count_words"
#define SCHEME_FLOAT_OPERATOR_READ_STAR "float=read_star"
#define SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MAX "float=star_table_max"
#define SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MIN "float=star_table_min"
#define SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_AVG "float=star_table_avg"
#define SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MAX_IDX "float=star_table_max_idx"
#define SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_MIN_IDX "float=star_table_min_idx"
#define SCHEME_FLOAT_OPERATOR_READ_STAR_TABLE_SORT_IDX "float=star_table_sort_idx"
#define SCHEME_STRING_OPERATOR_SET "string=set"
#define SCHEME_STRING_OPERATOR_JOIN "string=join"
#define SCHEME_STRING_OPERATOR_BEFORE_FIRST "string=before_first"
#define SCHEME_STRING_OPERATOR_AFTER_FIRST "string=after_first"
#define SCHEME_STRING_OPERATOR_BEFORE_LAST "string=before_last"
#define SCHEME_STRING_OPERATOR_AFTER_LAST "string=after_last"
#define SCHEME_STRING_OPERATOR_READ_STAR "string=read_star"
#define SCHEME_STRING_OPERATOR_GLOB "string=glob"
#define SCHEME_STRING_OPERATOR_NTH_WORD "string=nth_word"
#define SCHEME_OPERATOR_TOUCH_FILE "touch_file"
#define SCHEME_OPERATOR_COPY_FILE "copy_file"
#define SCHEME_OPERATOR_MOVE_FILE "move_file"
#define SCHEME_OPERATOR_DELETE_FILE "delete_file"
#define SCHEME_WAIT_OPERATOR_SINCE_LAST_TIME "wait"
#define SCHEME_EMAIL_OPERATOR "email"
#define SCHEME_EXIT_MAXTIME "exit_maxtime"
#define SCHEME_EXIT_OPERATOR "exit"

// A class that performs operators on variables
class SchemerOperator
{
	public:
	std::string type, input1, input2, output;

	public:


	SchemerOperator() {};

	SchemerOperator(std::string _type, std::string _input1="undefined", std::string _input2="undefined", std::string _output="undefined");

	std::string initialise(std::string _type, std::string _input1="undefined", std::string _input2="undefined", std::string _output="undefined");

	// Generate a meaningful current_name for the operator
	std::string getName();

	// Read a specific value from a STAR file
	void readFromStarFile() const;

	bool performOperation() const;

};

#define SCHEME_NODE_JOB_MODE_NEW "new"
#define SCHEME_NODE_JOB_MODE_CONTINUE "continue"

class SchemerJob
{
	public:
	std::string current_name, mode;
	bool job_has_started;

	public:

	SchemerJob() {};

	SchemerJob(std::string _name, std::string _mode, bool _has_started = false)
	{
		current_name = _name;
		mode = _mode;
		job_has_started = _has_started;
	}

	// Perform operation and return TRUE if not a JOB; just return FALSE if a JOB
	bool performOperation();


};

// Send an email
void schemerSendEmail(std::string message, std::string subject = "Schemer");


// A class that defines the edges between a graph that defines execution order, where the nodes are individual JOB instances
// An edge can also be a fork, where the output is controlled through a boolean variable
class SchemerEdge
{
	public:
	std::string inputNode, outputNode, outputNodeTrue;
	std::string myBooleanVariable;
	bool is_fork;

	std::string getOutputNode() const;

	SchemerEdge(std::string _input, std::string _output, bool _is_fork, std::string _mybool, std::string _output_if_true)
	{
		inputNode = _input;
		outputNode= _output;
		is_fork = _is_fork;
		outputNodeTrue = _output_if_true;
		myBooleanVariable = _mybool;
	}

	SchemerEdge(std::string _input, std::string _output)
	{
		inputNode = _input;
		outputNode= _output;
		is_fork = false;
		outputNodeTrue = "undefined";
		myBooleanVariable = "undefined";
	}

};

class Scheme
{

public:

	std::string name, current_node;
	bool do_read_only;
	int verb;

	std::map<std::string, SchemerJob> jobs;
	std::vector<SchemerEdge> edges;

public:

	Scheme()
	{
		clear();
	}

	void clear();

	void setName(std::string _name)
	{
		if (name[name.length()-1] != '/') name = name + '/';
        name = _name;
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

    // Get/set Variables and Operators(schemer_floats is only visible in this file!)
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
    	SchemerOperator op(type, input1, input2, output);
    	return op.getName();
    }
    void setOperatorParameters(std::string name, std::string type, std::string input1, std::string input2, std::string output);
    void getOperatorParameters(std::string name, std::string &type, std::string &input1, std::string &input2, std::string &output);

    // Get vectors with current Variables / Operators
    std::map<std::string, SchemerFloatVariable> getCurrentFloatVariables();
    std::map<std::string, SchemerBooleanVariable> getCurrentBooleanVariables();
    std::map<std::string, SchemerStringVariable> getCurrentStringVariables();
    std::map<std::string, SchemerOperator> getCurrentOperators();

    // Get/set operators

    // Add variables
    void setVariable(std::string name, FileName value); // (Add new one if exists, otherwise set value)
    void setOriginalVariable(std::string name, FileName value); // (Add new one if exists, otherwise set original_value)
	void addFloatVariable(std::string name, RFLOAT value);
    void addBooleanVariable(std::string name, bool value);
    void addStringVariable(std::string name, FileName value);

    // Add operators (of any kind), also adds its corresponding node
    SchemerOperator initialiseOperator(std::string type, std::string input_name, std::string input2_name,
    		std::string output_name, std::string &error_message);
    void addOperator(SchemerOperator &op, std::string &myname);

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


    // Rename this schemer into a new directory
    void copy(FileName newname);

    // Add edges and forks in between the nodes
    bool checkUniqueInput(std::string inputnode_name);
    void addEdge(std::string inputnode_name, std::string outputnode_name);
    void addFork(std::string inputnode_name, std::string mybool_name, std::string outputnode_name, std::string outputnode_name_if_false );

    // Test integrity of the Scheme. Warn for unused variables, nodes, etc.
    bool isValid();

    std::string getNextNode();
    std::string getPreviousNode();

    bool gotoNextNode();
    bool gotoNextJob();

    // Checks a string for any jobnames that are part of this or any other Schemes and changes for their corresponding current_name
    // Returns true if my string was changed
    bool changeStringForJobnames(FileName &mystring, FileName current_node);

    // Execute an operator from the Schemer and return true if successful
    bool executeOperator(FileName current_node);

    // This read in the original job.star, checks for pipeliner node dependencies and $$ variables and adds the new job to the pipeliner if neccesary
    RelionJob prepareJob(FileName current_node);

    // Run the Scheme
    void run(PipeLine &pipeline);

    // Remove the lock file for read/write protection
    void unlock();

    // Abort a running scheme
    void abort();

};


#endif /* SCHEDULER_H_ */
