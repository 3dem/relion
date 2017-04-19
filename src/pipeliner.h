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

#ifndef PIPELINER_H_
#define PIPELINER_H_
#include <iostream>
#include <sstream>
#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include "src/metadata_table.h"
#include "src/pipeline_jobs.h"
#define DEFAULTPDFVIEWER "evince"

class Process
{

	public:
	std::string name;
	std::string alias;
	int type;
	int status;
	std::vector<long int> inputNodeList;  // List of Nodes of input to this process
	std::vector<long int> outputNodeList; // List of Nodes of output from this process

	// Constructor
	Process(std::string _name, int _type, int _status, std::string _alias="None")
	{
		name = _name;
		type = _type;
		status = _status;
		alias = _alias;
	}

	// Destructor
	~Process()
	{
		inputNodeList.clear();
		outputNodeList.clear();
	}

};

#define DO_LOCK true
#define DONT_LOCK false
// Forward definition
class PipeLineFlowChart;

#define PIPELINE_HAS_CHANGED ".pipeline_has_changed"
class PipeLine
{

	public:
	int job_counter;
	bool do_read_only;

	std::string name;
	std::vector<Node> nodeList; //list of all Nodes in the pipeline
	std::vector<Process> processList; //list of all Processes in the pipeline

	PipeLine()
	{
		name = "default";
		job_counter = 1;
	}

	~PipeLine()
	{
		clear();
	}

	void clear()
	{
		nodeList.clear();
		processList.clear();
		job_counter = 1;
	}

	void setName(std::string _name)
	{
		name = _name;
	}

	// Add a new input Edge to the list
	// Check whether Node with that name already exists in the Node list, and if so update that one
	// The input_for_process will be added to the inputForProcessList of this Node
	//
	void addNewInputEdge(Node &_Node, long int input_for_process);

	// Add a new output Edge to the list
	// Check whether Node with that name already exists in the Node list, and if so update that one
	// The output_from_process will be added to the outputFromProcessList of this Node
	//
	void addNewOutputEdge(long int output_from_process, Node &_Node);

	// Check whether Node already exists in the nodeList. If not add and return pointer to new node, otherwise return pointer to existing node
	// Also touch entry in .Nodes directory, use touch_if_not_exist for scheduled jobs
	long int addNode(Node &_Node, bool touch_if_not_exist = false);

	// Add a new Process to the list (no checks are performed)
	long int addNewProcess(Process &_Process, bool do_overwrite = false);

	// Find nodes or process (by name or alias)
	long int findNodeByName(std::string name);
	long int findProcessByName(std::string name);
	long int findProcessByAlias(std::string name);

	// Touch each individual Node name in the temporary Nodes directory
	// Return true if Node output file exists and temporary file is written, false otherwise
	bool touchTemporaryNodeFile(Node &node, bool touch_even_if_not_exist=false);
	void touchTemporaryNodeFiles(Process &process);

	// And delete these temporary files
	void deleteTemporaryNodeFile(Node &node);
	void deleteTemporaryNodeFiles(Process &process);

	// Re-make entries of all NodeNames in the hidden .Nodes directory (for file browsing for InputNode I/O)
	void remakeNodeDirectory();

	// Check for process completion by checking for the presence of all outputNode filenames
	// Returns true if any of the running processes has completed, false otherwise
	bool checkProcessCompletion();

	// Get the command line arguments for thisjob
	bool getCommandLineJob(RelionJob &thisjob, int current_job, bool is_main_continue, bool is_scheduled, bool do_makedir, std::vector<std::string> &commands,
			std::string &final_command, std::string &error_message);

	// Adds _job to the pipeline and return the id of the newprocess
	long int addJob(RelionJob &_job, int as_status, bool do_overwrite);

	// Runs a job and adds it to the pipeline
	bool runJob(RelionJob &_job, int &current_job, bool only_schedule, bool is_main_continue, bool is_scheduled, std::string &error_message);

	// Runs a series of scheduled jobs, possibly in a loop
	void runScheduledJobs(FileName fn_sched, FileName fn_jobids, int nr_repeat, long int minutes_wait);

	// If I'm deleting this_job from the pipeline, which Nodes and which Processes need to be deleted?
	void deleteJobGetNodesAndProcesses(int this_job, bool do_recursive, std::vector<bool> &deleteNodes, std::vector<bool> &deleteProcesses);

	// Given the lists of which Nodes and Processes to delete, now do the actual deleting
	void deleteNodesAndProcesses(std::vector<bool> &deleteNodes, std::vector<bool> &deleteProcesses);

	// Changes the status of this_job to finished in the pipeline, returns false is job hadn't started yet
	bool markAsFinishedJob(int this_job, std::string &error_message);

	// Set the alias for a job, return true for success, false otherwise
	bool setAliasJob(int this_job, std::string alias, std::string &error_message);

	// Make the flowchart for this job
	bool makeFlowChart(long int current_job, bool do_display_pdf, std::string &error_message);

	// Undelete a JOb from the pipeline
	void undeleteJob(FileName fn_undel);

	// Clean up intermediate files from this_job
	bool cleanupJob(int this_job, bool do_harsh, std::string &error_message);

	// Clean upintermediate files from all jobs in the pipeline
	bool cleanupAllJobs(bool do_harsh, std::string &error_message);

	void replaceFilesForImportExportOfScheduledJobs(FileName fn_in_dir, FileName fn_out_dir,
			std::vector<std::string> &find_pattern, std::vector<std::string> &replace_pattern);

	// Export all scheduled jobs
	bool exportAllScheduledJobs(std::string mydir, std::string &error_message);

	// Import previously exported jobs
	void importJobs(FileName fn_export);

	// Import a job into the pipeline
	// Return true if at least one process is imported correctly
	bool importPipeline(std::string _name);

	// Write out the pipeline to a STAR file
	void write(bool do_lock = false, FileName fn_del="", std::vector<bool> deleteNode = std::vector<bool>(), std::vector<bool> deleteProcess = std::vector<bool>());

	// Read in the pipeline from a STAR file
	void read(bool do_lock = false);

};



class PipeLineFlowChart
{
public:

	// Use short process names, or original, full ones
	bool do_short_names;

	// Also make upwardsFlowCharts for all branches?
	bool do_branches;

	// All the processes for which a upwardFlowChart will be made
	std::vector<long int> todo_list;

public:

	PipeLineFlowChart()
	{
		do_branches= true;
		do_short_names = false;
	}

	// Write how many particles or classes or whatever the node is that represents a downward arrow
	std::string getDownwardsArrowLabel(PipeLine &pipeline, long int lower_process, long int new_process);

	// The process will be added to the top
	// The function returns the parent process from which the upper_node came
	// It will return a negative value if there was no parent process
	long int addProcessToUpwardsFlowChart(std::ofstream &fh, PipeLine &pipeline, long int lower_process,
			long int new_process, std::vector<long int> &branched_procs);

	void makeOneUpwardsFlowChart(std::ofstream &fh, PipeLine &pipeline, long int from_node,
			std::vector<long int> &all_branches, bool is_main_flow);

	void makeAllUpwardsFlowCharts(FileName &fn_out, PipeLine &pipeline, long int from_process);

	// Open and close a new flowchart picture
	void openTikZPicture(std::ofstream &fh, bool is_main_flow);

	void closeTikZPicture(std::ofstream &fh, bool is_main_flow);

	void adaptNamesForTikZ(FileName &name);

	// Open and close a new output file
	void openFlowChartFile(FileName &fn_out, std::ofstream &fh);

	void closeFlowChartFile(std::ofstream &fh);


};

#endif /* PIPELINER_H_ */
