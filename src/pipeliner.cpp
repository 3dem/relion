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

#include "src/pipeliner.h"
#include <unistd.h>

//#define DEBUG


long int PipeLine::addNode(Node &_Node, bool touch_if_not_exist)
{

	if (_Node.name=="")
		REPORT_ERROR("PipeLine::addNode ERROR: Adding an empty nodename. Did you fill in all Node names correctly?");

	// Check if _Node has an aliased name, and if so, revert back to the original name!
	FileName fn_node = _Node.name;
	for (size_t i = 0; i < processList.size(); i++)
	{
		if (fn_node.contains(processList[i].alias))
		{
			// Replace the alias by the name
			fn_node = processList[i].name + fn_node.without(processList[i].alias);
			_Node.name = fn_node;
			break;
		}
	}

	bool is_found = false;
	long int i;
	for (i=0; i < nodeList.size(); i++)
	{
		if (nodeList[i].name == _Node.name)
		{
			is_found = true;
			break;
		}
	}
	if (!is_found)
	{
		nodeList.push_back(_Node);
	}

	touchTemporaryNodeFile(nodeList[i], touch_if_not_exist);

	return i;

}


void PipeLine::addNewInputEdge(Node &_Node, long int myProcess)
{

	// 1. Check whether Node with that name already exists in the Node list
	long int oldsize = nodeList.size();
	long int  myNode = addNode(_Node);
	long int newsize = nodeList.size();

	// 2. Set the input_for_process in the inputForProcessList of this Node
	// But only if it doesn't exist yet
	bool found = false;
	for (size_t i = 0; i < (nodeList[myNode]).inputForProcessList.size(); i++)
	{
		if ((nodeList[myNode]).inputForProcessList[i] == myProcess)
		{
			found=true;
			break;
		}
	}
	if (!found)
	{
		(nodeList[myNode]).inputForProcessList.push_back(myProcess);
		(processList[myProcess]).inputNodeList.push_back(myNode);
	}

	if (newsize > oldsize)
	{
		// This is a previously unobserved Node, being used as input in a new Process. Check whether it came from an old Process
		for (long int i=0; i < processList.size(); i++)
		{
			// TODO: check this name comparison!!!
			FileName nodename = (nodeList[myNode]).name;
			if (nodename.contains(processList[i].name))
			{
#ifdef DEBUG
				std::cerr << " nodename.find(processList[i].name) = " << nodename.find(processList[i].name) << " processList[i].name= " << processList[i].name << std::endl;
				std::cerr << " found! " <<  nodename << " as output from " << processList[i].name << std::endl;
#endif
				(processList[i].outputNodeList).push_back(myNode);
				nodeList[myNode].outputFromProcess = i;
				break;
			}
		}
	}
}

void PipeLine::addNewOutputEdge(long int myProcess, Node &_Node)
{

	long int old_size = nodeList.size();

	// 1. Check whether Node with that name already exists in the Node list
	// Touch .Nodes entries even if they don't exist for scheduled jobs
	bool touch_if_not_exist = (processList[myProcess].status == PROC_SCHEDULED);

	// 2. Set the output_from_process of this Node
	_Node.outputFromProcess = myProcess;
	long int myNode = addNode(_Node, touch_if_not_exist);

	// 3. Only for new Nodes, add this Node to the outputNodeList of myProcess
	if (myNode == old_size)
		processList[myProcess].outputNodeList.push_back(myNode);

}

long int PipeLine::addNewProcess(Process &_Process, bool do_overwrite)
{
	// Check whether Process with the same name already exists in the processList
	bool is_found = false;
	long int i;
	for (i=0; i < processList.size(); i++)
	{
		if (processList[i].name == _Process.name)
		{
			is_found = true;
			processList[i].status = _Process.status;
			break;
		}
	}
	if (!is_found)
	{
		processList.push_back(_Process);
		job_counter++;
	}
	else if (!do_overwrite)
	{
		REPORT_ERROR("PipeLine::addNewProcess: ERROR: trying to add existing Process to the pipeline, while overwriting is not allowed.");
	}
	return i;
}

long int PipeLine::findNodeByName(std::string name)
{
	for(long int ipos = 0 ; ipos < nodeList.size() ; ipos++)
	{
		if(nodeList[ipos].name == name)
			return ipos;
	}
	return -1;
}


long int PipeLine::findProcessByName(std::string name)
{
	for(long int ipos = 0 ; ipos < processList.size() ; ipos++)
	{
		if(processList[ipos].name == name)
			return ipos;
	}
	return -1;
}

long int PipeLine::findProcessByAlias(std::string name)
{
	if (name == "None")
		return -1;

	for(long int ipos = 0 ; ipos < processList.size() ; ipos++)
	{
		if(processList[ipos].alias == name)
			return ipos;
	}
	return -1;
}

bool PipeLine::touchTemporaryNodeFile(Node &node, bool touch_even_if_not_exist)
{
	if (do_read_only)
		return false;

	FileName fn_dir = ".Nodes/";
	FileName fnt;

	// Check whether there is an alias for the corresponding process
	FileName fn_alias = (node.outputFromProcess < 0) ? "None" : processList[node.outputFromProcess].alias;

	if (fn_alias != "None")
	{
		// Make sure fn_alias ends with a slash
		if (fn_alias[fn_alias.length()-1] != '/')
			fn_alias += "/";
		FileName fn_pre, fn_jobnr, fn_post;
		if (decomposePipelineFileName(node.name, fn_pre, fn_jobnr, fn_post))
			fnt = fn_alias + fn_post;
		else
			REPORT_ERROR("PipeLine::touchTemporaryNodeFile ERROR: invalid node name: " + node.name);
	}
	else
		fnt = node.name;

	if (exists(node.name) || touch_even_if_not_exist)
	{
		// Make subdirectory for each type of node
		FileName fn_type = integerToString(node.type) + "/";
		FileName mydir = fn_dir + fn_type + fnt.substr(0, fnt.rfind("/") + 1);
		FileName mynode = fn_dir + fn_type + fnt;
		std::string command;
		if (!exists(mydir))
		{
			command = "mkdir -p " + mydir;
			int res = system(command.c_str());
		}
		touch(mynode);
		return true;
	}
	else
		return false;
}

void PipeLine::touchTemporaryNodeFiles(Process &process)
{

	if (do_read_only)
		return;

	bool touch_if_not_exist = (process.status == PROC_SCHEDULED);

	for (int j = 0; j < process.outputNodeList.size(); j++)
	{
		long int mynode = process.outputNodeList[j];
		touchTemporaryNodeFile(nodeList[mynode], touch_if_not_exist);
	}

}

void PipeLine::deleteTemporaryNodeFile(Node &node)
{
	if (do_read_only)
		return;

	FileName fn_dir = ".Nodes/";
	FileName fnt;

	// Check whether there is an alias for the corresponding process
	FileName fn_alias = (node.outputFromProcess < 0) ? "None" : processList[node.outputFromProcess].alias;

	if (fn_alias != "None")
	{
		// Make sure fn_alias ends with a slash
		if (fn_alias[fn_alias.length()-1] != '/')
			fn_alias += "/";
		FileName fn_pre, fn_jobnr, fn_post;
		if (decomposePipelineFileName(node.name, fn_pre, fn_jobnr, fn_post))
			fnt = fn_alias + fn_post;
		else
			REPORT_ERROR("PipeLine::deleteTemporaryNodeFile ERROR: invalid node name: " + node.name);
	}
	else
		fnt = node.name;

	FileName fn_type = integerToString(node.type) + "/";
	FileName fn = fn_dir + fn_type + fnt;
	int res = remove(fn.c_str());

	// Also remove the directory if it is empty
	// TODO: check what happens if the directory is not empty yet....
	fn = fn.beforeLastOf("/");
	res = remove(fn.c_str());

}

void PipeLine::deleteTemporaryNodeFiles(Process &process)
{
	if (do_read_only)
		return;

	for (int j = 0; j < process.outputNodeList.size(); j++)
	{
		long int mynode = process.outputNodeList[j];
		deleteTemporaryNodeFile(nodeList[mynode]);
	}

}

void PipeLine::remakeNodeDirectory()
{
	if (do_read_only)
		return;

	// Clear existing directory
	FileName fn_dir = ".Nodes/";
	std::string command = " rm -rf " + fn_dir;
	int res = system(command.c_str());

	for (long int i = 0; i < nodeList.size(); i++)
	{
		int myproc = nodeList[i].outputFromProcess;
		bool touch_if_not_exist = (myproc < 0) ? false : (processList[myproc].status == PROC_SCHEDULED);
		touchTemporaryNodeFile(nodeList[i], touch_if_not_exist);
	}
	command = "chmod 777 -R " + fn_dir;
	res = system(command.c_str());
}


bool PipeLine::checkProcessCompletion()
{
	if (do_read_only)
		return false;

	std::vector<long int> finished_success_processes;
	std::vector<long int> finished_failure_processes;
	std::vector<long int> finished_aborted_processes;

	for (long int i=0; i < processList.size(); i++)
	{
		// Only check running processes for file existence
		if (processList[i].status == PROC_RUNNING)
		{

			if (exists(processList[i].name + RELION_JOB_EXIT_SUCCESS))
				finished_success_processes.push_back(i);
			else if (exists(processList[i].name + RELION_JOB_EXIT_FAILURE))
				finished_failure_processes.push_back(i);
			else if (exists(processList[i].name + RELION_JOB_EXIT_ABORTED))
				finished_aborted_processes.push_back(i);

		}
	}

	// Only do read/write cycle in case a process was finished, otherwise the GUI slows down too much
	if (finished_success_processes.size() > 0 ||
		finished_failure_processes.size() > 0 ||
		finished_aborted_processes.size() > 0 ||
		exists(PIPELINE_HAS_CHANGED) )
	{
		// Read in the latest version of the pipeline, just in case anyone else made a change meanwhile...
		std::string lock_message = "";
		if (finished_success_processes.size() > 0)
		{
			lock_message += "checkProcessCompletion: the following jobs have successfully finished: ";
			for (long int i = 0; i < finished_success_processes.size(); i++)
				lock_message += " " + processList[finished_success_processes[i]].name;
			lock_message += "\n";
		}
		if (finished_failure_processes.size() > 0)
		{
			lock_message += "checkProcessCompletion: the following jobs have failed with an error: ";
			for (long int i = 0; i < finished_failure_processes.size(); i++)
				lock_message += " " + processList[finished_failure_processes[i]].name;
			lock_message += "\n";
		}
		if (finished_aborted_processes.size() > 0)
		{
			lock_message += "checkProcessCompletion: the following jobs have been aborted: ";
			for (long int i = 0; i < finished_aborted_processes.size(); i++)
				lock_message += " " + processList[finished_aborted_processes[i]].name;
			lock_message += "\n";
		}

		read(DO_LOCK, lock_message);

		// Set the new status of all the finished processes
		for (int i=0; i < finished_success_processes.size(); i++)
		{
			int myproc = finished_success_processes[i];
			processList[myproc].status = PROC_FINISHED_SUCCESS;

			// Also see whether there was an output nodes starfile
			getOutputNodesFromStarFile(myproc);

			// Also touch the outputNodes in the .Nodes directory
			for (long int j = 0; j < processList[myproc].outputNodeList.size(); j++)
			{
				int myNode = (processList[myproc]).outputNodeList[j];
				if (myNode < 0 || myNode >= nodeList.size())
					REPORT_ERROR("pipeline checkProcessCompletion ERROR: " + integerToString(j) + "th output node of " + processList[myproc].name + " is invalid: " + integerToString(myNode));
				if (!touchTemporaryNodeFile(nodeList[myNode]))
					REPORT_ERROR("ERROR: output node " + nodeList[myNode].name + " does not exist, while job " + processList[myproc].name +" has finished successfully...");
			}

		}
		for (int i=0; i < finished_failure_processes.size(); i++)
		{
			processList[finished_failure_processes[i]].status = PROC_FINISHED_FAILURE;
		}
		for (int i=0; i < finished_aborted_processes.size(); i++)
		{
			processList[finished_aborted_processes[i]].status = PROC_FINISHED_ABORTED;
		}

		// Always couple read/write with DO_LOCK
		// This is to make sure two different windows do not get out-of-sync
		write(DO_LOCK);
		if (exists(PIPELINE_HAS_CHANGED))
			std::remove(PIPELINE_HAS_CHANGED);
		return true;
	}

	return false;

}

bool PipeLine::getCommandLineJob(RelionJob &thisjob, int current_job, bool is_main_continue,
                                 bool is_scheduled, bool do_makedir, bool do_overwrite_current,
                                 std::vector<std::string> &commands, std::string &final_command, std::string &error_message)
{

	if (do_overwrite_current) is_main_continue = false;

	// Except for continuation or scheduled jobs, all jobs get a new unique directory
	std::string my_outputname;
	if ((is_main_continue || is_scheduled || do_overwrite_current) && current_job < processList.size())
	{
		if (current_job < 0)
			REPORT_ERROR("BUG: current_job < 0");
		FileName fn_settings = processList[current_job].name;
		my_outputname = fn_settings.beforeLastOf("/") + "/";
	}
	else
		my_outputname = "";

	// Set is_continue flag inside the job
	thisjob.is_continue = is_main_continue;

	bool result = thisjob.getCommands(my_outputname, commands, final_command, do_makedir, job_counter, error_message);

	if (result && commands.size()==0)
	{
		error_message = " PipeLine::getCommandLineJob: Nothing to do...";
		return false;
	}

	return result;
}

// Adds thisjob to the pipeline and returns the id of the newprocess
long int PipeLine::addJob(RelionJob &thisjob, int as_status, bool do_overwrite, bool do_write_minipipeline)
{

	// Also write a mini-pipeline in the output directory
	PipeLine mini_pipeline;
	mini_pipeline.setName(thisjob.outputName+"job");

	// Add Process to the processList of the pipeline
	Process process(thisjob.outputName, thisjob.type, as_status);
	long int myProcess = addNewProcess(process, do_overwrite);
	mini_pipeline.addNewProcess(process);

	// Add all input nodes
	for (int i=0; i < thisjob.inputNodes.size(); i++)
	{
		addNewInputEdge(thisjob.inputNodes[i], myProcess);
		mini_pipeline.addNewInputEdge(thisjob.inputNodes[i], 0);
	}
	// Add all output nodes
	for (int i=0; i < thisjob.outputNodes.size(); i++)
	{
		addNewOutputEdge(myProcess, thisjob.outputNodes[i]);
		mini_pipeline.addNewOutputEdge(0, thisjob.outputNodes[i]);
	}

	if (do_write_minipipeline)
	{
		// Write the mini-pipeline to an updated STAR file
		mini_pipeline.write();
	}
	// Writing of the overall pipeline is done in the function calling addToPipeLine

	return myProcess;

}

bool PipeLine::runJob(RelionJob &_job, int &current_job, bool only_schedule, bool is_main_continue,
                      bool is_scheduled, bool do_overwrite_current, std::string &error_message)
{

	std::vector<std::string> commands;
	std::string final_command;

	// Remove run.out and run.err when overwriting a job
	if (do_overwrite_current) is_main_continue = false;

	// true means makedir
	if (!getCommandLineJob(_job, current_job, is_main_continue, is_scheduled, true, do_overwrite_current, commands, final_command, error_message))
	{
		return false;
	}

	// Remove run.out and run.err when overwriting a job
	if (do_overwrite_current)
	{
		// Completely empty the output directory, NOTE that  _job.outputName+ is not defined until AFTER calling getCommandLineJob!!!
		std::string command = " rm -rf " + _job.outputName + "*";
		int res = system(command.c_str());

		// Above deletes run_submit.script too, so we have to call this again ...
		if (!getCommandLineJob(_job, current_job, is_main_continue, is_scheduled, true, do_overwrite_current, commands, final_command, error_message))
			return false;
	}

	// Read in the latest version of the pipeline, just in case anyone else made a change meanwhile...
	std::string lock_message = "runJob: " + _job.outputName;
	read(DO_LOCK, lock_message);

	// Save temporary hidden file with this jobs settings as default for a new job
	_job.write("");

	// Also save a copy of the GUI settings with the current output name
	_job.write(_job.outputName);

	// Make sure none of the exit or abort files from the pipeline_control system are here from before
	std::remove((_job.outputName+RELION_JOB_ABORT_NOW).c_str());
	std::remove((_job.outputName+RELION_JOB_EXIT_ABORTED).c_str());
	std::remove((_job.outputName+RELION_JOB_EXIT_SUCCESS).c_str());
	std::remove((_job.outputName+RELION_JOB_EXIT_FAILURE).c_str());


	/*
	// If this is a continuation job, check whether output files exist and move away!
	// This is to ensure that the continuation job goes OK and will show up as 'running' in the GUI
	bool do_move_output_nodes_to_old = false;
	if (!only_schedule && is_main_continue)
	{
		do_move_output_nodes_to_old = !(processList[current_job].type == PROC_2DCLASS ||
		                                processList[current_job].type == PROC_3DCLASS ||
		                                processList[current_job].type == PROC_INIMODEL ||
		                                processList[current_job].type == PROC_3DAUTO ||
		                                processList[current_job].type == PROC_MULTIBODY ||
		                                processList[current_job].type == PROC_MANUALPICK ||
		                                processList[current_job].type == PROC_CLASSSELECT);

		// For continuation of relion_refine jobs, remove the original output nodes from the list
		if (processList[current_job].type == PROC_2DCLASS ||
		    processList[current_job].type == PROC_3DCLASS ||
		    processList[current_job].type == PROC_3DAUTO ||
		    processList[current_job].type == PROC_MULTIBODY ||
		    processList[current_job].type == PROC_INIMODEL)
		{

			std::vector<bool> deleteNodes, deleteProcesses;
			deleteNodes.resize(nodeList.size(), false);
			deleteProcesses.resize(processList.size(), false);

			for (long int inode = 0; inode < (processList[current_job]).outputNodeList.size(); inode++)
			{
				long int mynode = (processList[current_job]).outputNodeList[inode];
				if(!exists(nodeList[mynode].name))
					deleteNodes[mynode] = true;
			}

			FileName fn_del = "tmp";
			write(DO_LOCK, fn_del, deleteNodes, deleteProcesses);
			std::remove("tmpdeleted_pipeline.star");

			// Read the updated pipeline back in again
			lock_message += " part 2";
			read(DO_LOCK, lock_message);

		}
	} // end if !only_schedule && is_main_continue

	// If a job is executed with a non-continue scheduled job, then also move away any existing output node files
	if (current_job >= 0 && (is_scheduled && !is_main_continue) || do_overwrite_current)
		do_move_output_nodes_to_old = true;

	// Move away existing output nodes
	if (do_move_output_nodes_to_old)
	{

		for (int i = 0; i < processList[current_job].outputNodeList.size(); i++)
		{
			int j = processList[current_job].outputNodeList[i];
			std::string fn_node = nodeList[j].name;
			if (exists(fn_node))
			{
				std::string path2 =  fn_node + ".old";
				rename(fn_node.c_str(), path2.c_str());
			}
		}
	}
	*/

	// For continuation of relion_refine jobs, remove the original output nodes from the list
	if (!only_schedule && is_main_continue)
	{
		if (processList[current_job].type == PROC_2DCLASS ||
		    processList[current_job].type == PROC_3DCLASS ||
		    processList[current_job].type == PROC_3DAUTO ||
		    processList[current_job].type == PROC_MULTIBODY ||
		    processList[current_job].type == PROC_INIMODEL)
		{

			std::vector<bool> deleteNodes, deleteProcesses;
			deleteNodes.resize(nodeList.size(), false);
			deleteProcesses.resize(processList.size(), false);

			for (long int inode = 0; inode < (processList[current_job]).outputNodeList.size(); inode++)
			{
				long int mynode = (processList[current_job]).outputNodeList[inode];
				if(!exists(nodeList[mynode].name))
					deleteNodes[mynode] = true;
			}

			FileName fn_del = "tmp";
			write(DO_LOCK, fn_del, deleteNodes, deleteProcesses);
			std::remove("tmpdeleted_pipeline.star");

			// Read the updated pipeline back in again
			lock_message += " part 2";
			read(DO_LOCK, lock_message);

		}
	} // end if !only_schedule && is_main_continue

	// Now save the job (and its status) to the PipeLine
	int mynewstatus;
	if (only_schedule)
		mynewstatus = PROC_SCHEDULED;
	else
		mynewstatus = PROC_RUNNING;
	bool allow_overwrite = is_main_continue || is_scheduled; // continuation and scheduled jobs always allow overwriting into the existing directory

	// Add the job to the pipeline, and set current_job to the new one
	current_job = addJob(_job, mynewstatus, allow_overwrite || do_overwrite_current);

	// Write out the new pipeline
	write(DO_LOCK);

	// Now actually execute the Job
	if (!only_schedule)
	{
		//std::cout << "Executing: " << final_command << std::endl;
		int res = system(final_command.c_str());

		// Also print the final_command to the note for future reference
		FileName fn_note = processList[current_job].name + "note.txt";
		std::ofstream ofs;
		ofs.open (fn_note.c_str(), std::ofstream::out | std::ofstream::app);

		// current date/time based on current system
		time_t now = time(0);
		ofs << std::endl << " ++++ Executing new job on " << ctime(&now);
		ofs <<  " ++++ with the following command(s): " << std::endl;
		for (size_t icom = 0; icom < commands.size(); icom++)
			ofs << commands[icom] << std::endl;
		ofs <<  " ++++ " << std::endl;
		ofs.close();
	}

	// Copy pipeline star file as backup to the output directory
	FileName fn_pipe = name + "_pipeline.star";
	if (exists(fn_pipe))
		copy(fn_pipe, processList[current_job].name + fn_pipe);

	return true;
}

// Adds a scheduled job to the pipeline from the command line
int PipeLine::addScheduledJob(std::string typestring, std::string fn_options)
{

	int type;
	if (typestring == PROC_IMPORT_NAME)
		type = PROC_IMPORT;
	else if (typestring == PROC_MOTIONCORR_NAME)
		type = PROC_MOTIONCORR;
	else if (typestring == PROC_CTFFIND_NAME)
		type = PROC_CTFFIND;
	else if (typestring == PROC_MANUALPICK_NAME)
		type = PROC_MANUALPICK;
	else if (typestring == PROC_AUTOPICK_NAME)
		type = PROC_AUTOPICK;
	else if (typestring == PROC_EXTRACT_NAME)
		type = PROC_EXTRACT;
	else if (typestring == PROC_CLASSSELECT_NAME)
		type = PROC_CLASSSELECT;
	else if (typestring == PROC_2DCLASS_NAME)
		type = PROC_2DCLASS;
	else if (typestring == PROC_3DCLASS_NAME)
		type = PROC_3DCLASS;
	else if (typestring == PROC_3DAUTO_NAME)
		type = PROC_3DAUTO;
	else if (typestring == PROC_MASKCREATE_NAME)
		typestring = PROC_MASKCREATE;
	else if (typestring == PROC_JOINSTAR_NAME)
		type = PROC_JOINSTAR;
	else if (typestring == PROC_SUBTRACT_NAME)
		type = PROC_SUBTRACT;
	else if (typestring == PROC_POST_NAME)
		type = PROC_POST;
	else if (typestring == PROC_RESMAP_NAME)
		type = PROC_RESMAP;
	else if (typestring == PROC_INIMODEL_NAME)
		type = PROC_INIMODEL;
	else
		REPORT_ERROR("ERROR: unrecognised string for job type: " + typestring);

	return addScheduledJob(type, fn_options);

}

// Adds a scheduled job to the pipeline from the command line
int PipeLine::addScheduledJob(int job_type, std::string fn_options)
{
	RelionJob job;
	job.initialise(job_type);
	std::vector<std::string> options;
	splitString(fn_options, ";", options);
	for (long int iopt = 0; iopt < options.size(); iopt++)
		job.setOption(options[iopt]);

	// Always add Pre-processing jobs as continuation ones (for convenient on-the-fly processing)
	if (job_type == PROC_MOTIONCORR || job_type == PROC_CTFFIND || job_type == PROC_AUTOPICK  || job_type == PROC_EXTRACT)
		job.is_continue = true;

	std::string error_message;
	int current_job = processList.size();
	if (!runJob(job, current_job, true, job.is_continue, false, false, error_message)) // true is only_schedule, false means !is_scheduled, 2nd false means dont overwrite current
		REPORT_ERROR(error_message.c_str());

	return current_job;
}

// Adds a scheduled job to the pipeline from the command line
int PipeLine::addScheduledJob(RelionJob &job)
{
	std::string error_message;
	int current_job = processList.size();
	if (!runJob(job, current_job, true, job.is_continue, false, false, error_message)) // true is only_schedule, false means !is_scheduled, 2nd false means dont overwrite current
		REPORT_ERROR(error_message.c_str());

	return current_job;
}

void PipeLine::waitForJobToFinish(int current_job, bool &is_failure, bool &is_aborted)
{

	while (true)
    {
        sleep(10);
        checkProcessCompletion();
        if (processList[current_job].status == PROC_FINISHED_SUCCESS ||
            processList[current_job].status == PROC_FINISHED_ABORTED ||
            processList[current_job].status == PROC_FINISHED_FAILURE)
        {
            // Prepare a string for a more informative .lock file
            std::string lock_message = " pipeliner noticed that " + processList[current_job].name + " finished and is trying to update the pipeline";

            // Read in existing pipeline, in case some other window had changed something else
            read(DO_LOCK, lock_message);

            if (processList[current_job].status == PROC_FINISHED_FAILURE)
                is_failure = true;
            else if (processList[current_job].status == PROC_FINISHED_ABORTED)
                is_aborted = true;

            // Write out the modified pipeline with the new status of current_job
            write(DO_LOCK);
            break;

        } // endif something has happened
    } // while true, waiting for job to finish


}

void PipeLine::runScheduledJobs(FileName fn_sched, FileName fn_jobids, int nr_repeat,
		long int minutes_wait, long int minutes_wait_before, long int seconds_wait_after, bool do_overwrite_current)
{

	std::vector<FileName> my_scheduled_processes;
	std::vector<std::string> jobids;
	int njobs = splitString(fn_jobids, " ", jobids);
	if (njobs == 0)
	{
		REPORT_ERROR("PipeLine::runScheduledJobs: Nothing to do...");
	}
	else
	{
		for (int i = 0; i < njobs; i++)
			my_scheduled_processes.push_back(jobids[i]);
	}

	FileName fn_log = "pipeline_" + fn_sched + ".log";
	std::ofstream  fh;
	fh.open((fn_log).c_str(), std::ios::app);

	if (nr_repeat > 1)
	{
		std::cout << " PIPELINER: writing out information in logfile " << fn_log << std::endl;
	}

	// Touch the fn_check file
	FileName fn_check = "RUNNING_PIPELINER_" + fn_sched;
	bool fn_check_exists = false;
	if (nr_repeat > 1)
	{
		touch(fn_check);
		fn_check_exists = true;
	}

	fh << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	fh << " Starting a new scheduler execution called " << fn_sched << std::endl;
	fh << " The scheduled jobs are: " << std::endl;
	for (long int i = 0; i < my_scheduled_processes.size(); i++)
		fh << " - " << my_scheduled_processes[i] << std::endl;
	if (nr_repeat > 1)
	{
		if (minutes_wait_before > 0)
			fh << " Will wait " << minutes_wait_before << " minutes before running the first job." << std::endl;
		fh << " Will execute the scheduled jobs " << nr_repeat << " times." << std::endl;
		if (nr_repeat > 1)
			fh << " Will wait until at least " << minutes_wait << " minutes have passed between each repeat." << std::endl;
		fh << " Will be checking for existence of file " << fn_check << "; if it no longer exists, the scheduler will stop." << std::endl;
	}
	fh << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

	// Wait this many minutes before starting the repeat cycle...
	if (minutes_wait_before > 0)
		sleep(minutes_wait_before * 60);

	bool is_failure = false;
	bool is_aborted = false;

	int repeat = 0;
	time_t now = time(0);
	for (repeat = 0 ; repeat < nr_repeat; repeat++)
	{
		if (nr_repeat > 1)
		{
			fh << " + " << ctime(&now) << " -- Starting the " << repeat+1 << "th repeat" << std::endl;
		}

		// Get starting time of the repeat cycle
		timeval time_start, time_end;
		gettimeofday(&time_start, NULL);

		for (long int i = 0; i < my_scheduled_processes.size(); i++)
		{
			int current_job = findProcessByName(my_scheduled_processes[i]);
			if (current_job < 0)
			{
				// Also try finding it by alias
				current_job = findProcessByAlias(my_scheduled_processes[i]);
				if (current_job < 0)
					REPORT_ERROR("ERROR: cannot find process with name: " + my_scheduled_processes[i]);
			}
			RelionJob myjob;
			bool is_continue;
			if (!myjob.read(processList[current_job].name, is_continue, true)) // true means also initialise the job
				REPORT_ERROR("There was an error reading job: " + processList[current_job].name);

			// Check whether the input nodes are there, before executing the job
			for (long int inode = 0; inode < processList[current_job].inputNodeList.size(); inode++)
			{
				long int mynode = processList[current_job].inputNodeList[inode];
				while (!exists(nodeList[mynode].name))
				{
					fh << " + -- Warning " << nodeList[mynode].name << " does not exist. Waiting 60 seconds ... " << std::endl;
					sleep(60);
				}
			}
			now = time(0);
			fh << " + " << ctime(&now) << " ---- Executing " << processList[current_job].name  << std::endl;
			std::string error_message;

			if (!runJob(myjob, current_job, false, is_continue, true, do_overwrite_current, error_message)) // true means is_scheduled; false=dont overwrite current
				REPORT_ERROR(error_message);

			// Now wait until that job is done!
			while (true)
			{
				if (nr_repeat > 1 && !exists(fn_check))
				{
					fn_check_exists = false;
					break;
				}

				sleep(seconds_wait_after);
				checkProcessCompletion();
				if (processList[current_job].status == PROC_FINISHED_SUCCESS ||
					processList[current_job].status == PROC_FINISHED_ABORTED ||
					processList[current_job].status == PROC_FINISHED_FAILURE)
				{
					// Prepare a string for a more informative .lock file
					std::string lock_message = " Scheduler " + fn_sched + " noticed that " + processList[current_job].name +
							" finished and is trying to update the pipeline";

					// Read in existing pipeline, in case some other window had changed something else
					read(DO_LOCK, lock_message);

					if (processList[current_job].status == PROC_FINISHED_SUCCESS)
					{
						// Will we go on to do another repeat?
						if (repeat + 1 != nr_repeat)
						{
							int mytype = processList[current_job].type;
							// The following jobtypes have functionality to only do the unfinished part of the job
							if (mytype == PROC_MOTIONCORR || mytype == PROC_CTFFIND || mytype == PROC_AUTOPICK || mytype == PROC_EXTRACT
									|| mytype == PROC_CLASSSELECT )
							{
								myjob.is_continue = true;
								// Write the job again, now with the updated is_continue status
								myjob.write(processList[current_job].name);
							}
							processList[current_job].status = PROC_SCHEDULED;
						}
						else
						{
							processList[current_job].status = PROC_FINISHED_SUCCESS;
						}
					}
					else if (processList[current_job].status == PROC_FINISHED_FAILURE)
					{
						is_failure = true;
					}
					else if (processList[current_job].status == PROC_FINISHED_ABORTED)
					{
						is_aborted = true;
					}

					// Write out the modified pipeline with the new status of current_job
					write(DO_LOCK);
					break;
				}
			}

			// break out of scheduled processes loop
			if (is_failure || is_aborted)
				break;

			if (nr_repeat > 1 && !fn_check_exists)
				break;

		} //end loop my_scheduled_processes


		// break out of repeat loop
		if (nr_repeat > 1 && !fn_check_exists)
			break;

		if (is_failure || is_aborted)
			break;

		// Wait at least until 'minutes_wait' minutes have passed from the beginning of the repeat cycle
		gettimeofday(&time_end, NULL);
		long int passed_minutes = (time_end.tv_sec - time_start.tv_sec)/60;
		long int still_wait = minutes_wait - passed_minutes;
		if (still_wait > 0 && repeat+1 != nr_repeat)
		{
			fh << " + -- Waiting " << still_wait << " minutes until next repeat .."<< std::endl;
			sleep(still_wait * 60);
		}

	}

	if (is_failure)
	{
		fh << " + Stopping pipeliner due to a job that failed with an error ..." << std::endl;
	}
	else if (is_aborted)
	{
		fh << " + Stopping pipeliner due to user abort .. " << std::endl;
	}
	else if (repeat == nr_repeat)
	{
		if (nr_repeat > 1)
		{
			fh << " + performed all requested repeats in scheduler " << fn_sched << ". Stopping pipeliner now ..." << std::endl;
		}
		else
		{
			fh << " + All jobs have finished, so stopping pipeliner now ..." << std::endl;
		}

		// Read in existing pipeline, in case some other window had changed it
		std::string lock_message = " Scheduler " + fn_sched + " has finished and is trying to update the pipeline";
		read(DO_LOCK, lock_message);

		// After breaking out of repeat, set status of the jobs to finished
		for (long int i = 0; i < my_scheduled_processes.size(); i++)
		{
			int current_job = findProcessByName(my_scheduled_processes[i]);
			processList[current_job].status = PROC_FINISHED_SUCCESS;
		}

		// Write the pipeline to an updated STAR file
		write(DO_LOCK);

		// Remove the temporary file
		std::remove(fn_check.c_str());
	}
	else if (!fn_check_exists && nr_repeat > 1)
	{
		fh << " + File " << fn_check << " was removed. Stopping now .." << std::endl;
		std::cout << " PIPELINER: the " << fn_check << " file was removed. Stopping now ..." << std::endl;
	}

	fh << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	fh.close();

}

void PipeLine::deleteJobGetNodesAndProcesses(int this_job, bool do_recursive, std::vector<bool> &deleteNodes, std::vector<bool> &deleteProcesses)
{

	deleteProcesses.resize(processList.size(), false);
	deleteNodes.resize(nodeList.size(), false);

	std::vector<long int> to_delete_processes;
	to_delete_processes.push_back(this_job);

	bool is_done = false;
	size_t istart = 0;
	while (!is_done)
	{
		size_t imax = to_delete_processes.size();
		for (long int i = istart; i < imax; i++)
		{
			// re-set istart for next recursive round
			istart = imax;
			long int idel = to_delete_processes[i];
			deleteProcesses[idel] = true;
			is_done = true;
			for (size_t inode = 0; inode < (processList[idel]).outputNodeList.size(); inode++)
			{
				long int mynode = (processList[idel]).outputNodeList[inode];
				deleteNodes[mynode] = true;

				if (do_recursive)
				{
					// Check whether this node is being used as input for another process, and if so, delete those as well
					for (size_t ii = 0; ii < (nodeList[mynode]).inputForProcessList.size(); ii++)
					{
						long int iproc = (nodeList[mynode]).inputForProcessList[ii];
						// See if this process is not already in the list to be deleted
						bool already_in = false;
						for (size_t j = 0; j < to_delete_processes.size(); j++)
						{
							if (to_delete_processes[j] == iproc)
								already_in = true;
						}
						if (!already_in)
						{
							to_delete_processes.push_back(iproc);
							is_done = false;
						}
					}
				}
			}
		}
	}

}

void PipeLine::deleteNodesAndProcesses(std::vector<bool> &deleteNodes, std::vector<bool> &deleteProcesses)
{

	// Read in existing pipeline, in case some other window had changed it
	std::string lock_message = "deleteNodesAndProcesses";
	read(DO_LOCK, lock_message);

	// Write new pipeline without the deleted processes and nodes to disc and read in again
	FileName fn_del;
	for (int i = 0; i < deleteProcesses.size(); i++)
		if (deleteProcesses[i])
		{
			fn_del = processList[i].name;
			break;
		}
	write(DO_LOCK, fn_del, deleteNodes, deleteProcesses);

	// Delete the output directories for all selected processes from the hard disk
	// Do this after pipeline.write to get the deleted_pipeline.star still in the correct directory
	for (int i = 0; i < deleteProcesses.size(); i++)
	{
		if (deleteProcesses[i])
		{
			FileName alldirs = processList[i].name;
			alldirs = alldirs.beforeLastOf("/");
			// Move entire output directory (with subdirectory structure) to the Trash folder
			FileName firstdirs = alldirs.beforeLastOf("/");
			FileName fn_tree="Trash/" + firstdirs;
			int res = mktree(fn_tree);
			std::string command = "mv -f " + alldirs + " " + "Trash/" + firstdirs+"/.";
			res = system(command.c_str());
			// Also remove the symlink if it exists
			FileName fn_alias = (processList[i]).alias;
			if (fn_alias != "None")
			{
				int res = unlink((fn_alias.beforeLastOf("/")).c_str());
			}

			deleteTemporaryNodeFiles(processList[i]);
		}
	}

	// Read new pipeline back in again
	lock_message += " part 2";
	read(DO_LOCK, lock_message);
	write(DO_LOCK);

}

void PipeLine::getOutputNodesFromStarFile(int this_job)
{

	// See if a STAR file called RELION_OUTPUT_NODES.star exists, and if so, read in which output nodes were created
	FileName outnodes = processList[this_job].name + "RELION_OUTPUT_NODES.star";
	if (exists(outnodes))
	{

		MetaDataTable MDnodes;
		MDnodes.read(outnodes, "output_nodes");

		FileName nodename;
		int nodetype;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDnodes)
		{
			MDnodes.getValue(EMDL_PIPELINE_NODE_NAME, nodename);
			MDnodes.getValue(EMDL_PIPELINE_NODE_TYPE, nodetype);

			// if this node does not exist yet, then add it to the pipeline
			if (findNodeByName(nodename) < 0 )
			{
				Node node(nodename, nodetype);
				addNewOutputEdge(this_job, node);
			}
		}

	}

}

bool PipeLine::markAsFinishedJob(int this_job, std::string &error_message)
{

	// Read in existing pipeline, in case some other window had changed it
	std::string lock_message = "markAsFinishedJob";
	read(DO_LOCK, lock_message);

	processList[this_job].status = PROC_FINISHED_SUCCESS;

	// For relion_refine jobs, add last iteration optimiser.star, data.star, model.star and class???.mrc to the pipeline
	if (processList[this_job].type == PROC_2DCLASS ||
		processList[this_job].type == PROC_3DCLASS ||
		processList[this_job].type == PROC_3DAUTO  ||
		processList[this_job].type == PROC_INIMODEL )
	{
		// Get the last iteration optimiser file
		FileName fn_opt;
		FileName fn_root1 = (processList[this_job].alias != "None") ? processList[this_job].alias : processList[this_job].name;

		std::vector<FileName> fn_opts;
		fn_opt = fn_root1 + "run_it*optimiser.star";
		fn_opt.globFiles(fn_opts);
		// It could also be a continuation
		fn_opt = fn_root1 + "run_ct?_it???_optimiser.star";
		fn_opt.globFiles(fn_opts, false); // false means: don't clear fn_opts vector
		// It could also be a continuation
		fn_opt = fn_root1 + "run_ct??_it???_optimiser.star";
		fn_opt.globFiles(fn_opts, false); // false means: don't clear fn_opts vector
		if (fn_opts.size() > 0)
		{

			fn_opt = fn_opts[fn_opts.size()-1]; // the last one

			// Also get data.star
			FileName fn_data = fn_opt.without("_optimiser.star") + "_data.star";
			Node node2(fn_data, NODE_PART_DATA);
			addNewOutputEdge(this_job, node2);

			FileName fn_root = fn_opt.without("_optimiser.star");
			if (processList[this_job].type == PROC_3DAUTO)
				fn_root += "_half1";

			FileName fn_model = fn_root + "_model.star";
			Node node3(fn_model, NODE_MODEL);
			addNewOutputEdge(this_job, node3);


			FileName fn_map = fn_root + "_class???.mrc";
			std::vector<FileName> fn_maps;
			fn_map.globFiles(fn_maps);
			for (int i = 0; i < fn_maps.size(); i++)
			{
				Node node4(fn_maps[i], NODE_3DREF);
				addNewOutputEdge(this_job, node4);
			}
		}
		else
		{
			error_message = "You are trying to mark a relion_refine job as finished that hasn't even started. \n This will be ignored. Perhaps you wanted to delete it instead?";
			processList[this_job].status = PROC_RUNNING;
			write(DO_LOCK);
			return false;
		}
	}

	// Also see whether there is an output nodes star file...
	getOutputNodesFromStarFile(this_job);

	// Remove any of the expected output nodes from the pipeline if the corresponding file doesn't already exist
	std::vector<bool> deleteNodes, deleteProcesses;
	deleteNodes.resize(nodeList.size(), false);
	deleteProcesses.resize(processList.size(), false);

	for (long int inode = 0; inode < (processList[this_job]).outputNodeList.size(); inode++)
	{
		long int mynode = (processList[this_job]).outputNodeList[inode];
		if(!exists(nodeList[mynode].name))
			deleteNodes[mynode] = true;
	}
	FileName fn_del = "tmp";
	write(DO_LOCK, fn_del, deleteNodes, deleteProcesses);
	std::remove("tmpdeleted_pipeline.star");

	// Read the updated pipeline back in and write it out again
	// With the locking mechanism, each pipeline.read(bool, DO_LOCK) needs to be followed soon by a pipeline.write(DO_LOCK)!
	lock_message += " part 2";
	read(DO_LOCK, lock_message);
	write(DO_LOCK);

	return true;

}

// Set the alias for a job
bool PipeLine::setAliasJob(int this_job, std::string alias, std::string &error_message)
{

	FileName fn_pre, fn_jobnr, fn_post;
	if (!decomposePipelineFileName(processList[this_job].name, fn_pre, fn_jobnr, fn_post))
		REPORT_ERROR("PipeLine::setAlias ERROR: invalid pipeline process name: " + processList[this_job].name);

	if (alias.length() == 0)
	{
		alias = "None";
	}
	else if (alias.length() < 2)
	{
		 error_message = "Alias cannot be less than 2 characters, please provide another one";
		 return false;
	}
	else if (alias.length() > 2 && alias[0]=='j' && alias[1]=='o' && alias[2]=='b')
	{
		 error_message = "Alias cannot start with 'job', please provide another one";
		 return false;
	}
	else if (alias.find("*") != std::string::npos || alias.find("?") != std::string::npos || alias.find("(") != std::string::npos || alias.find(")") != std::string::npos ||
	         alias.find("/") != std::string::npos || alias.find("\"") != std::string::npos || alias.find("\\") != std::string::npos || alias.find("|") != std::string::npos ||
		 alias.find("#") != std::string::npos || alias.find("<") != std::string::npos || alias.find(">") != std::string::npos || alias.find("&") != std::string::npos ||
		 alias.find("%") != std::string::npos || alias.find("{") != std::string::npos || alias.find("}") != std::string::npos || alias.find("$") != std::string::npos)
	{
		error_message = "Alias cannot contain following symbols: *, ?, (, ), /, \", \\, |, #, <, >, &, %, {, }, $";
		return false;
	}
	else
	{

		//remove spaces from any potential alias
		for (int i = 0; i < alias.length(); i++)
		{
			if (alias[i] == ' ')
				alias[i] = '_';
		}

		// Make sure the alias ends with a slash
		if (alias[alias.length()-1] != '/')
			alias += "/";

		// Check uniqueness of the alias
		bool is_unique = true;
		for (size_t i = 0; i < processList.size(); i++)
		{
			if ( processList[i].alias == fn_pre + alias && alias != "None")
			{
				is_unique = false;
				break;
			}
		}
		if (!is_unique || alias.length() < 1)
		{
			 error_message ="Alias is not unique, please provide another one";
			 return false;
		}
	}


	// Read in existing pipeline, in case some other window had changed it
	std::string lock_message = "setAliasJob";
	read(DO_LOCK, lock_message);

	// Remove the original .Nodes entry
	deleteTemporaryNodeFiles(processList[this_job]);

	FileName fn_old_alias = processList[this_job].alias;
	if (fn_old_alias != "None")
		fn_old_alias = fn_old_alias.beforeLastOf("/");

	// No alias if the alias contains a unique jobnr string because of continuation of relion_refine jobs
	if (alias == "None" )
	{
		processList[this_job].alias = "None";
	}
	else
	{
		// If this was already an alias: remove the old symbolic link
		if (fn_old_alias != "None")
			int res2 = unlink(fn_old_alias.c_str());

		// Set the alias in the pipeline
		processList[this_job].alias = fn_pre + alias;

		//Make the new symbolic link
		FileName path1 = "../" + processList[this_job].name;
		FileName path2 = processList[this_job].alias;
		int res = symlink(path1.c_str(), path2.beforeLastOf("/").c_str());

	}

	// Remake the new .Nodes entry
	touchTemporaryNodeFiles(processList[this_job]);

	// Write new pipeline to disc
	write(DO_LOCK);

	return true;

}

bool PipeLine::makeFlowChart(long int current_job, bool do_display_pdf, std::string &error_message)
{
	if (current_job < 0)
	{
		error_message = " You can only make flowcharts for existing jobs ... ";
		return false;
	}

	const char * default_pdf_viewer = getenv ("RELION_PDFVIEWER_EXECUTABLE");
	char mydefault[]=DEFAULTPDFVIEWER;
	if (default_pdf_viewer == NULL)
	{
		default_pdf_viewer=mydefault;
	}
	std::string myviewer(default_pdf_viewer);

	PipeLineFlowChart flowchart;
	FileName fn_dir = processList[current_job].name;
	FileName fn_out = "flowchart.tex";
	flowchart.makeAllUpwardsFlowCharts(fn_out, *this, current_job);
	std::string command = "latex flowchart.tex > flowchart.log && dvipdf flowchart >>flowchart.log && mv flowchart* " + fn_dir;
	std:: cout << " Executing: " << command << std::endl;
	int res = std::system(command.c_str());
	if (do_display_pdf)
	{
		command = myviewer + " " + fn_dir + "flowchart.pdf &";
		res = std::system(command.c_str());
	}

	// Read in existing pipeline, in case some other window had changed it
	std::string lock_message = "makeFlowChart";
	read(DO_LOCK, lock_message);

	// Add the PDF file as a logfile to the outputnodes of this job, so it can be visualised from the Display button
	Node node(fn_dir+"flowchart.pdf", NODE_PDF_LOGFILE);
	addNewOutputEdge(current_job, node);

	write(DO_LOCK);

	return true;
}
// Undelete a Job from the pipeline, move back from Trash and insert back into the graph
void PipeLine::undeleteJob(FileName fn_undel)
{

	// Read in existing pipeline, in case some other window had changed it
	std::string lock_message = "undeleteJob";
	read(DO_LOCK, lock_message);

	importPipeline(fn_undel.beforeLastOf("_pipeline.star"));

	// Copy all processes in the STAR file back into the ProjectDirectory
	MetaDataTable MDproc;
	MDproc.read(fn_undel, "pipeline_processes");
	std::cout <<"  Undeleting from Trash ... " << std::endl;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDproc)
	{
		FileName fn_proc;
		MDproc.getValue(EMDL_PIPELINE_PROCESS_NAME, fn_proc);

		// Copy the job back from the Trash folder
		FileName fn_dest = fn_proc.beforeLastOf("/"); //gets rid of ending "/"
		FileName fn_dir_dest = fn_dest.beforeLastOf("/"); // Now only get the job-type directory
		if (!exists(fn_dir_dest))
		{
			mktree(fn_dir_dest);
		}
		std::string command = "mv Trash/" + fn_dest + " " + fn_dest;
		std::cout << command << std::endl;
		int res = system(command.c_str());

		// Also re-make all entries in the .Nodes directory
		long int myproc = findProcessByName(fn_proc);
		touchTemporaryNodeFiles(processList[myproc]);
	}
	std::cout << " Done undeleting! " << std::endl;

	// Write the new pipeline to disk and reread it back in again
	write(DO_LOCK);

}

bool PipeLine::cleanupJob(int this_job, bool do_harsh, std::string &error_message)
{

	if (this_job < 0 || processList[this_job].status != PROC_FINISHED_SUCCESS)
	{
		error_message =" You can only clean up finished jobs ... ";
		return false;
	}

	// These job types do not have cleanup:
	if (processList[this_job].type == PROC_IMPORT ||
		processList[this_job].type == PROC_MANUALPICK ||
		processList[this_job].type == PROC_CLASSSELECT ||
		processList[this_job].type == PROC_MASKCREATE ||
		processList[this_job].type == PROC_JOINSTAR ||
		processList[this_job].type == PROC_RESMAP)
		return true;

	// Find any subdirectories
	std::vector<FileName> fns_subdir;
	FileName fn_curr_dir = "";
	int idir = -1, istop = 0;
	bool is_first = true;
	// Recursively find all subdirectories
	while (idir < istop)
	{
		FileName fn_curr_dir = (is_first) ? processList[this_job].name : processList[this_job].name + fns_subdir[idir];
		DIR *dir = opendir(fn_curr_dir.c_str());
		struct dirent *entry = readdir(dir);
		while (entry != NULL)
		{
			// Only want directories, and not '.' or '..'
			if (entry->d_type == DT_DIR && (entry->d_name[0] != '.'))
			{
				FileName fnt = (is_first) ? entry->d_name : fns_subdir[idir] + entry->d_name;
				fns_subdir.push_back(fnt + "/");
			}
			entry = readdir(dir);
		}
		closedir(dir);
		istop = fns_subdir.size();
		idir++;
		is_first = false;
	}

	std::vector<FileName> fns_del;
	FileName fn_pattern;

	// In all jobs cleanup the .old files (from continuation runs)
	//fn_pattern = processList[this_job].name + "*.old";
	//fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del

	////////// Now see which jobs needs cleaning up
	if (processList[this_job].type == PROC_MOTIONCORR)
	{

		for (int idir = 0; idir < fns_subdir.size(); idir++)
		{
			if (do_harsh)
			{
				//remove entire movies directory
				fns_del.push_back(processList[this_job].name + fns_subdir[idir]);
			}
			else
			{
				fn_pattern = processList[this_job].name + fns_subdir[idir] + "*.com";
				fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
				fn_pattern = processList[this_job].name + fns_subdir[idir] + "*.err";
				fn_pattern.globFiles(fns_del, false);
				fn_pattern = processList[this_job].name + fns_subdir[idir] + "*.out";
				fn_pattern.globFiles(fns_del, false);
				fn_pattern = processList[this_job].name + fns_subdir[idir] + "*.log";
				fn_pattern.globFiles(fns_del, false);
			}
		}

	} // end if motioncorr
	else if (processList[this_job].type == PROC_CTFFIND)
	{

		fn_pattern = processList[this_job].name + "gctf*.out";
		fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
		fn_pattern = processList[this_job].name + "gctf*.err";
		fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
		for (int idir = 0; idir < fns_subdir.size(); idir++)
		{
			//remove entire Micrographs directory structure
			fns_del.push_back(processList[this_job].name + fns_subdir[idir]);
		}

	} // end if ctffind
	else if (processList[this_job].type == PROC_AUTOPICK)
	{

		for (int idir = 0; idir < fns_subdir.size(); idir++)
		{
			// remove the Spider files with the FOM maps
			fn_pattern = processList[this_job].name + fns_subdir[idir] + "*.spi";
			fn_pattern.globFiles(fns_del, false);
		}

	} // end if autopick
	else if (processList[this_job].type == PROC_EXTRACT)
	{

		for (int idir = 0; idir < fns_subdir.size(); idir++)
		{
			if (do_harsh)
			{
				//remove entire directory (STAR files and particle stacks!
				fns_del.push_back(processList[this_job].name + fns_subdir[idir]);
			}
			else
			{
				// only remove the STAR files with the metadata (this will only give moderate file savings)
				fn_pattern = processList[this_job].name + fns_subdir[idir] + "*_extract.star";
				fn_pattern.globFiles(fns_del, false);
			}
		}

	} // end if extract
	else if (processList[this_job].type == PROC_2DCLASS ||
			 processList[this_job].type == PROC_3DCLASS ||
			 processList[this_job].type == PROC_3DAUTO ||
			 processList[this_job].type == PROC_INIMODEL ||
			 processList[this_job].type == PROC_MULTIBODY)
	{

		// First find the _data.star from each iteration
		std::vector<FileName> fns_iter;
		fn_pattern = processList[this_job].name + "run_it[0-9][0-9][0-9]_data.star";
		fn_pattern.globFiles(fns_iter);
		fn_pattern = processList[this_job].name + "run_ct[0-9]_it[0-9][0-9][0-9]_data.star";
		fn_pattern.globFiles(fns_iter, false);
		fn_pattern = processList[this_job].name + "run_ct[0-9][0-9]_it[0-9][0-9][0-9]_data.star";
		fn_pattern.globFiles(fns_iter, false);
		fn_pattern = processList[this_job].name + "run_ct[0-9][0-9][0-9]_it[0-9][0-9][0-9]_data.star";
		fn_pattern.globFiles(fns_iter, false);
		// Keep everything for the last iteration, such thatone could for example still do a multibody refinement after gentle cleaning
		for (int ifile = 0; ifile < (signed int)(fns_iter.size())-1; ifile++)
		{
			FileName fn_file = (fns_iter[ifile]).without("_data.star");
			// Find the iterations to keep: i.e. those that are part of the pipeline
			bool is_in_pipeline = false;
			for (long int inode = 0; inode < nodeList.size(); inode++)
			{
				FileName fn_node = nodeList[inode].name;
				if (fn_node.contains(fn_file))
				{
					is_in_pipeline = true;
					break;
				}
			}
			// Delete all files from this iteration
			if (!is_in_pipeline)
			{
				fn_pattern = fn_file + "*";
				fn_pattern.globFiles(fns_del, false);
			}

			// Also clean up maps for PCA movies when doing harsh cleaning
			if (do_harsh && processList[this_job].type == PROC_MULTIBODY)
			{
				fn_pattern = processList[this_job].name + "analyse_component???_bin???.mrc";
				fn_pattern.globFiles(fns_del, false);
			}

		} //end loop over ifile (i.e. the _data.star files from all iterations)

	} // end if refine job
	else if (processList[this_job].type == PROC_CTFREFINE)
	{

		for (int idir = 0; idir < fns_subdir.size(); idir++)
		{
			// remove the temporary output files
			fn_pattern = processList[this_job].name + fns_subdir[idir] + "*_wAcc.mrc";
			fn_pattern.globFiles(fns_del, false);
			fn_pattern = processList[this_job].name + fns_subdir[idir] + "*_xyAcc_real.mrc";
			fn_pattern.globFiles(fns_del, false);
			fn_pattern = processList[this_job].name + fns_subdir[idir] + "*_xyAcc_imag.mrc";
			fn_pattern.globFiles(fns_del, false);
		}

	} // end if ctf_refine
	else if (processList[this_job].type == PROC_MOTIONREFINE)
	{

		for (int idir = 0; idir < fns_subdir.size(); idir++)
		{
			// remove the temporary output files
			fn_pattern = processList[this_job].name + fns_subdir[idir] + "*_FCC_cc.mrc";
			fn_pattern.globFiles(fns_del, false);
			fn_pattern = processList[this_job].name + fns_subdir[idir] + "*_FCC_w0.mrc";
			fn_pattern.globFiles(fns_del, false);
			fn_pattern = processList[this_job].name + fns_subdir[idir] + "*_FCC_w1.mrc";
			fn_pattern.globFiles(fns_del, false);

			if (do_harsh)
			{
				fn_pattern = processList[this_job].name + fns_subdir[idir] + "*_shiny.mrcs";
				fn_pattern.globFiles(fns_del, false);
				fn_pattern = processList[this_job].name + fns_subdir[idir] + "*_shiny.star";
				fn_pattern.globFiles(fns_del, false);
			}
		}

	} // end if motion_refine
	else if (processList[this_job].type == PROC_SUBTRACT)
	{

		if (do_harsh)
		{
			fn_pattern = processList[this_job].name + "subtracted.*";
			fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
		}

	} // end if subtract
	else if (processList[this_job].type == PROC_POST)
	{

		fn_pattern = processList[this_job].name + "*masked.mrc";
		fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del

	} // end if postprocess


	// Now actually move all the files
	FileName fn_old_dir = "";
	for (long int idel = 0; idel < fns_del.size(); idel++)
	{
		FileName fn_dest = "Trash/" + fns_del[idel];
		FileName fn_dir = fn_dest.beforeLastOf("/");
		if (fn_dir != fn_old_dir && ! exists(fn_dir))
			int res = mktree(fn_dir);
		// by removing entire directories, it could be the file is gone already
		if (exists(fns_del[idel]))
		{
			std::string command = "mv -f " + fns_del[idel] + " "+ fn_dir;
			int res = system(command.c_str());
		}
	} // end loop over all files to be deleted

	return true;

}

// Clean upintermediate files from all jobs in the pipeline
bool PipeLine::cleanupAllJobs(bool do_harsh, std::string &error_message)
{

	for (int myjob = 0; myjob < processList.size(); myjob++)
	{
		if (processList[myjob].status == PROC_FINISHED_SUCCESS)
		{
			if (do_harsh && exists(processList[myjob].name + "NO_HARSH_CLEAN"))
				continue;
			if (!cleanupJob(myjob, do_harsh, error_message))
				return false;
		}
	}

	return true;

}

void PipeLine::replaceFilesForImportExportOfScheduledJobs(FileName fn_in_dir, FileName fn_out_dir, std::vector<std::string> &find_pattern, std::vector<std::string> &replace_pattern)
{
	int res;
	std::string command;
	std::vector<std::string> myfiles;
	myfiles.push_back("run.job");
	myfiles.push_back("note.txt");
	myfiles.push_back("job_pipeline.star");

	// Copy the run.job, the note.txt and the job_pipeline.star
	// Replace all instances of all find_pattern's with the replace_pattern's
	for (int ifile = 0; ifile < myfiles.size(); ifile++)
	{
		for (int ipatt = 0; ipatt < find_pattern.size(); ipatt++)
		{
			FileName outfile = fn_out_dir + myfiles[ifile];
			FileName tmpfile = fn_out_dir + "tmp";
			FileName infile = (ipatt == 0) ? fn_in_dir + myfiles[ifile] : tmpfile;
			// Create directory first time round
			if (ipatt == 0)
			{
				FileName dirs = outfile.beforeLastOf("/");
				command =  "mkdir -p " + dirs;
				res = system(command.c_str());
			}
			command =  "sed 's|" + find_pattern[ipatt] + "|" + replace_pattern[ipatt] + "|g' < " + infile + " > " + outfile;
			//std::cerr << " Executing: " << command<<std::endl;
			res = system(command.c_str());
			if (ipatt+1 < find_pattern.size())
			{
				std::rename(outfile.c_str(), tmpfile.c_str());
				//std::cerr << " Excuting: mv " << outfile<<" "<<tmpfile<<std::endl;
			}
		}
	}

}

bool PipeLine::exportAllScheduledJobs(std::string mydir, std::string &error_message)
{

	// Make sure the directory name ends with a slash
	mydir += "/";
	std::string command = "mkdir -p ExportJobs/" + mydir;
	int res = system(command.c_str());

	MetaDataTable MDexported;

	// Loop through all the Scheduled jobs and export them one-by-one
	int iexp =0;
	std::vector<std::string> find_pattern, replace_pattern;
	for (long int i = 0; i < processList.size(); i++)
	{
		if (processList[i].status == PROC_SCHEDULED)
		{
			iexp++;
			if (processList[i].alias != "None")
			{
				error_message = "ERROR: aliases are not allowed on Scheduled jobs that are to be exported! Make sure all scheduled jobs are made with unaliases names.";
				return false;
			}

			// A general name for the exported job:
			FileName expname = processList[i].name;
			expname = expname.beforeFirstOf("/") + "/exp"+integerToString(iexp, 3)+"/";
			find_pattern.push_back(processList[i].name);
			replace_pattern.push_back(expname);

			MDexported.addObject();
			MDexported.setValue(EMDL_PIPELINE_PROCESS_NAME, expname);

			// Copy the run.job, the note.txt and the job_pipeline.star and replace patterns
			replaceFilesForImportExportOfScheduledJobs(processList[i].name, "ExportJobs/" + mydir + expname, find_pattern, replace_pattern);
		}
	}

	MDexported.write("ExportJobs/" + mydir + "exported.star");
	return true;

}

void PipeLine::importJobs(FileName fn_export)
{

	FileName fn_export_dir = fn_export.beforeLastOf("/")+"/";

	//FileName fn_dir_export = fn_export.beforeLastOf("/")+"/";
	MetaDataTable MDexported;
	MDexported.read(fn_export);

	// Read in existing pipeline, in case some other window had changed it
	std::string lock_message = "importJobs";
	read(DO_LOCK, lock_message);

	std::vector<std::string> find_pattern, replace_pattern;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDexported)
	{
		FileName expname;
		MDexported.getValue(EMDL_PIPELINE_PROCESS_NAME, expname);
		find_pattern.push_back(expname);
		// Make a new name for this job
		FileName newname = expname.beforeFirstOf("/")+"/job"+integerToString(job_counter, 3)+"/";
		//std::cerr << " expname= " << expname << " newname= " << newname << std::endl;
		replace_pattern.push_back(newname);
		replaceFilesForImportExportOfScheduledJobs(fn_export_dir + expname, newname, find_pattern, replace_pattern);
		// Import the job into the pipeline
	    importPipeline(newname+"job");
	    job_counter++;
	}

	// Write the new pipeline to disk
	write(DO_LOCK);

}

// Import a job into the pipeline
bool PipeLine::importPipeline(std::string _name)
{

	if (_name == name)
	{
		std::cerr << " importPipeline WARNING: ignoring request to import myself! "<<std::endl;
		return false;
	}

	PipeLine mini_pipeline;
	mini_pipeline.setName(_name);
	mini_pipeline.read();

	// TODO: vectors that map imported process and nodes numbers to the new ones!!
	std::vector<long int> imported_process_nr;
	std::vector<long int> imported_node_nr;
	long int ori_nr_proc = processList.size();
	long int ori_nr_node = nodeList.size();

	bool imported = false;
	for (int iproc = 0; iproc < mini_pipeline.processList.size(); iproc++)
	{
		// Check that the new processes all have unique names
		if (findProcessByName(mini_pipeline.processList[iproc].name) >= 0)
			REPORT_ERROR("importPipeline ERROR: cannot import pipeline with non-unique job name: " + mini_pipeline.processList[iproc].name);

		if (findProcessByAlias(mini_pipeline.processList[iproc].alias) >= 0)
		{
			std::cerr << "importPipeline WARNING: resetting non-unique imported alias: " << mini_pipeline.processList[iproc].alias << std::endl;
			mini_pipeline.processList[iproc].alias = "None";
		}

		imported_process_nr.push_back(processList.size());
		processList.push_back(mini_pipeline.processList[iproc]);
		imported = true;

	}

	if (imported)
	{
		for (int inode = 0; inode < mini_pipeline.nodeList.size(); inode++)
		{
			// Only push_back nodes that weren't present yet
			int mynode = findNodeByName(mini_pipeline.nodeList[inode].name);
			if (mynode < 0)
			{
				imported_node_nr.push_back(nodeList.size());
				nodeList.push_back(mini_pipeline.nodeList[inode]);
			}
			else
			{
				//std::cerr << "mynode= "<<mynode << " name=" << nodeList[mynode].name << std::endl;
				imported_node_nr.push_back(mynode);
			}
		}

		// Now fix the numbers of the lists in the imported processes and nodes
		for (int iproc = ori_nr_proc; iproc < processList.size(); iproc++)
		{
			for (int inode = 0; inode < processList[iproc].inputNodeList.size(); inode++)
			{
				int ori_node = processList[iproc].inputNodeList[inode];
				//std::cerr << " ori_node=" << ori_node << " name= "<<nodeList[ori_node].name << std::endl;
				//std::cerr << " imported_node_nr[ori_node]= " << imported_node_nr[ori_node] << " name= "<<nodeList[imported_node_nr[ori_node]].name<<std::endl;
				processList[iproc].inputNodeList[inode] = imported_node_nr[ori_node];
			}
			for (int inode = 0; inode < processList[iproc].outputNodeList.size(); inode++)
			{
				int ori_node = processList[iproc].outputNodeList[inode];
				processList[iproc].outputNodeList[inode] = imported_node_nr[ori_node];
			}
		}
		for (int inode = ori_nr_node; inode < nodeList.size(); inode++)
		{
			for (int iproc = 0; iproc < nodeList[inode].inputForProcessList.size(); iproc++)
			{
				int ori_proc = nodeList[inode].inputForProcessList[iproc];
				nodeList[inode].inputForProcessList[iproc] = imported_process_nr[ori_proc];
			}
			int ori_proc2 = nodeList[inode].outputFromProcess;
			nodeList[inode].outputFromProcess = imported_process_nr[ori_proc2];
		}

	}

	return imported;

}

// Read pipeline from STAR file
void PipeLine::read(bool do_lock, std::string lock_message)
{

#ifdef DEBUG_LOCK
	std::cerr << "entering read lock_message=" << lock_message << std::endl;
#endif
	FileName name_wo_dir = name;
	FileName dir_lock=".relion_lock", fn_lock=".relion_lock/lock_" + name_wo_dir.afterLastOf("/") + "_pipeline.star";;
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
				REPORT_ERROR("ERROR: PipeLine::read cannot create a lock directory " + dir_lock + ". You don't have write permission to this project. If you want to look at other's project directory (but run nothing there), please start RELION with --readonly.");

			// If the lock exists: wait 3 seconds and try again
			// First time round, print a warning message
			if (iwait == 0)
			{
				std::cout << " WARNING: trying to read pipeline.star, but directory " << dir_lock << " exists (which protects against simultaneous writing by multiple instances of the GUI)" << std::endl;
			}
			sleep(3);
			status =  mkdir(dir_lock.c_str(), S_IRWXU);
#ifdef DEBUG_LOCK
			std::cerr <<  " B status= " << status << std::endl;
#endif

			iwait++;
			if (iwait > 40)
			{

				REPORT_ERROR("ERROR: PipeLine::read has waited for 2 minutes for lock directory to disappear. You may want to manually remove the file: " + fn_lock);
			}

		}
		// Generate the lock file
		std::ofstream  fh;
	    fh.open(fn_lock.c_str(), std::ios::out);
	    if (!fh)
	        REPORT_ERROR( (std::string)"ERROR: Cannot open file: " + fn_lock);
	    fh << lock_message << std::endl;
	    fh.close();

	}


	// Start from scratch
	clear();

	FileName fn = name + "_pipeline.star";
	std::ifstream in(fn.c_str(), std::ios_base::in);

	if (in.fail())
		REPORT_ERROR( (std::string) "PipeLine::read: File " + fn + " cannot be read." );

	MetaDataTable MDgen, MDnode, MDproc, MDedge1, MDedge2;

	// This if allows for older version of the pipeline without the jobcounter
	// TODO: remove after alpha-testing
	if (MDgen.readStar(in, "pipeline_general"))
	{
		MDgen.getValue(EMDL_PIPELINE_JOB_COUNTER, job_counter);
	}

	MDnode.readStar(in, "pipeline_nodes");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDnode)
	{
		std::string name;
		int type;
		if (!MDnode.getValue(EMDL_PIPELINE_NODE_NAME, name) ||
			!MDnode.getValue(EMDL_PIPELINE_NODE_TYPE, type)	)
			REPORT_ERROR("PipeLine::read: cannot find name or type in pipeline_nodes table");

		Node newNode(name, type);
		nodeList.push_back(newNode);
	}

	MDproc.readStar(in, "pipeline_processes");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDproc)
	{
		std::string name, alias;
		int type, status;
		if (!MDproc.getValue(EMDL_PIPELINE_PROCESS_NAME, name) ||
			!MDproc.getValue(EMDL_PIPELINE_PROCESS_ALIAS, alias) ||
			!MDproc.getValue(EMDL_PIPELINE_PROCESS_TYPE, type) ||
			!MDproc.getValue(EMDL_PIPELINE_PROCESS_STATUS, status)	)
			REPORT_ERROR("PipeLine::read: cannot find name or type in pipeline_processes table");

		Process newProcess(name, type, status, alias);
		processList.push_back(newProcess);

		// Make a symbolic link to the alias if it isn't there...
		if (alias != "None")
		{
			// Also make a symbolic link for the output directory!
			// Make sure it doesn't end in a slash
			FileName fn_alias = alias;
			if (fn_alias[fn_alias.length()-1] == '/')
				fn_alias = fn_alias.beforeLastOf("/");

			// Only make the alias if it doesn't exist yet, otherwise you end up with recursive ones.
			if (!exists(fn_alias))
			{
				//std::string command = " ln -s ../" + name + " " + fn_alias;
				//int res= system(command.c_str());
				std::string path1 = "../" + name;
				int res = symlink(path1.c_str(), fn_alias.c_str());
			}
		}
	}

	// Read in all input (Node->Process) edges
	MDedge1.readStar(in, "pipeline_input_edges");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDedge1)
	{
		std::string fromnodename, procname;
		if (!MDedge1.getValue(EMDL_PIPELINE_EDGE_PROCESS, procname) ||
			!MDedge1.getValue(EMDL_PIPELINE_EDGE_FROM, fromnodename) )
			REPORT_ERROR("PipeLine::read: cannot find procname or fromnodename in pipeline_edges table");

		// Now fill in all To and FromEdgeLists of all Nodes
		long int myProcess = findProcessByName(procname);
		bool found_both = true;
		if (myProcess < 0 || myProcess >= processList.size())
		{
			std::cerr << "PipeLine WARNING: cannot find child process with name: " << procname << std::endl;
			found_both = false;
			//REPORT_ERROR("PipeLine::read ERROR: cannot find to-process with name: " + procname);
		}
		long int fromNode = findNodeByName(fromnodename);
		if (fromNode < 0 || fromNode >= nodeList.size())
		{
			std::cerr << "PipeLine WARNING: cannot find parent node with name: " << fromnodename << std::endl;
			found_both = false;
			//REPORT_ERROR("PipeLine::read ERROR: cannot find from-node with name: " + fromnodename);
		}
		if (found_both)
		{
			processList[myProcess].inputNodeList.push_back(fromNode);
			nodeList[fromNode].inputForProcessList.push_back(myProcess);
		}
	}

	// Read in all output (Process->Node) edges
	MDedge2.readStar(in, "pipeline_output_edges");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDedge2)
	{
		std::string tonodename, procname;
		if (!MDedge2.getValue(EMDL_PIPELINE_EDGE_TO, tonodename) ||
			!MDedge2.getValue(EMDL_PIPELINE_EDGE_PROCESS, procname) )
			REPORT_ERROR("PipeLine::read: cannot find procname or tonodename in pipeline_edges table");

		// Now fill in all To and FromEdgeLists of all Nodes
		long int myProcess = findProcessByName(procname);
		bool found_both = true;
		if (myProcess < 0 || myProcess >= processList.size())
		{
			std::cerr << "PipeLine WARNING: cannot find parent process with name: " << procname << std::endl;
			found_both = false;
			//REPORT_ERROR("PipeLine::read ERROR: cannot find from-process with name: " + procname);
		}
		long int toNode = findNodeByName(tonodename);
		if (toNode < 0 || toNode  >= nodeList.size())
		{
			std::cerr << "PipeLine WARNING: cannot find child node with name: " << tonodename << std::endl;
			found_both = false;
			//REPORT_ERROR("PipeLine::read ERROR: cannot find to-node with name: " + tonodename);
		}
		if (found_both)
		{
			processList[myProcess].outputNodeList.push_back(toNode);
			nodeList[toNode].outputFromProcess = myProcess;
		}
	}

	// Close file handler
	in.close();

}

void PipeLine::write(bool do_lock, FileName fn_del, std::vector<bool> deleteNode, std::vector<bool> deleteProcess)
{
	if (do_read_only)
		return;

	FileName name_wo_dir = name;
	FileName dir_lock=".relion_lock", fn_lock=".relion_lock/lock_" + name_wo_dir.afterLastOf("/") + "_pipeline.star";;
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
				REPORT_ERROR("ERROR: PipeLine::read has waited for 2 minutes for lock file to appear, but it doesn't. This should not happen. Is something wrong with the disk access?");
			}
		}
	}

	std::ofstream  fh, fh_del;
    FileName fn = name + "_pipeline.star";
    fh.open(fn.c_str(), std::ios::out);
	if (fh.fail())
		REPORT_ERROR("ERROR: cannot write to pipeline file: " + fn);

    if (fn_del != "")
    {
        FileName fnt = fn_del + "deleted_pipeline.star";
        fh_del.open(fnt.c_str(), std::ios::out);
        if (deleteNode.size() != nodeList.size())
        	REPORT_ERROR("PipeLine::write BUG: not enough entries in deleteNode vector!");
        if (deleteProcess.size() != processList.size())
        	REPORT_ERROR("PipeLine::write BUG: not enough entries in deleteProcess vector!");
    }

    MetaDataTable MDgen, MDnode, MDproc, MDedge1, MDedge2;
    MetaDataTable MDgen_del, MDnode_del, MDproc_del, MDedge1_del, MDedge2_del;

#ifdef DEBUG
    std::cerr << " writing pipeline as " << fn << std::endl;
#endif

    MDgen.setName("pipeline_general");
    MDgen.setIsList(true);
    MDgen.addObject();
    MDgen.setValue(EMDL_PIPELINE_JOB_COUNTER, job_counter);
    MDgen.write(fh);

    if (fn_del != "")
    {
		MDgen_del.setName("pipeline_general");
		MDgen_del.setIsList(true);
		MDgen_del.addObject();
		MDgen_del.setValue(EMDL_PIPELINE_JOB_COUNTER, job_counter);
		MDgen_del.write(fh_del);
    }

    MDproc.setName("pipeline_processes");
    MDproc_del.setName("pipeline_processes");
    for(long int i=0 ; i < processList.size() ; i++)
    {
    	if (fn_del == "" || !deleteProcess[i])
    	{
			MDproc.addObject();
			MDproc.setValue(EMDL_PIPELINE_PROCESS_NAME, processList[i].name);
			MDproc.setValue(EMDL_PIPELINE_PROCESS_ALIAS, processList[i].alias);
			MDproc.setValue(EMDL_PIPELINE_PROCESS_TYPE, processList[i].type);
			MDproc.setValue(EMDL_PIPELINE_PROCESS_STATUS, processList[i].status);
    	}
    	else
    	{
			MDproc_del.addObject();
			MDproc_del.setValue(EMDL_PIPELINE_PROCESS_NAME, processList[i].name);
			MDproc_del.setValue(EMDL_PIPELINE_PROCESS_ALIAS, processList[i].alias);
			MDproc_del.setValue(EMDL_PIPELINE_PROCESS_TYPE, processList[i].type);
			MDproc_del.setValue(EMDL_PIPELINE_PROCESS_STATUS, processList[i].status);
    	}

    }
#ifdef DEBUG
    MDproc.write(std::cerr);
#endif
    MDproc.write(fh);
    if (fn_del != "")
    	MDproc_del.write(fh_del);

    MDnode.setName("pipeline_nodes");
    MDnode_del.setName("pipeline_nodes");
    for(long int i=0 ; i < nodeList.size() ; i++)
    {
    	if (fn_del == "" || !deleteNode[i])
    	{
			MDnode.addObject();
			MDnode.setValue(EMDL_PIPELINE_NODE_NAME, nodeList[i].name);
			MDnode.setValue(EMDL_PIPELINE_NODE_TYPE, nodeList[i].type);
    	}
    	else
    	{
			MDnode_del.addObject();
			MDnode_del.setValue(EMDL_PIPELINE_NODE_NAME, nodeList[i].name);
			MDnode_del.setValue(EMDL_PIPELINE_NODE_TYPE, nodeList[i].type);

    	}
    }
#ifdef DEBUG
    MDnode.write(std::cerr);
#endif
    MDnode.write(fh);
    if (fn_del != "")
    	MDnode_del.write(fh_del);

    // Also write all (Node->Process) edges to a single table
    MDedge1.setName("pipeline_input_edges");
    MDedge1_del.setName("pipeline_input_edges");
    for(long int i=0 ; i < processList.size() ; i++)
    {
		for (long int j=0; j < ((processList[i]).inputNodeList).size(); j++)
		{
			long int inputNode = ((processList[i]).inputNodeList)[j];
	    	if (fn_del == "" || (!deleteProcess[i] && !deleteNode[inputNode]) )
	    	{
				MDedge1.addObject();
				MDedge1.setValue(EMDL_PIPELINE_EDGE_FROM, nodeList[inputNode].name);
				MDedge1.setValue(EMDL_PIPELINE_EDGE_PROCESS, processList[i].name);
	    	}
	    	else
	    	{
				MDedge1_del.addObject();
				MDedge1_del.setValue(EMDL_PIPELINE_EDGE_FROM, nodeList[inputNode].name);
				MDedge1_del.setValue(EMDL_PIPELINE_EDGE_PROCESS, processList[i].name);
	    	}
    	}
    }
#ifdef DEBUG
    MDedge1.write(std::cerr);
#endif
    MDedge1.write(fh);
    if (fn_del != "")
        MDedge1_del.write(fh_del);

    // Also write all (Process->Node) edges to a single table
    MDedge2.setName("pipeline_output_edges");
    MDedge2_del.setName("pipeline_output_edges");
    for(long int i=0 ; i < processList.size() ; i++)
    {
    	for (long int j=0; j < ((processList[i]).outputNodeList).size(); j++)
    	{
    		long int outputNode = ((processList[i]).outputNodeList)[j];
	    	if (fn_del == "" || (!deleteProcess[i] && !deleteNode[outputNode]) )
	    	{
				MDedge2.addObject();
				MDedge2.setValue(EMDL_PIPELINE_EDGE_PROCESS,  processList[i].name);
				MDedge2.setValue(EMDL_PIPELINE_EDGE_TO, nodeList[outputNode].name);
	    	}
	    	else
	    	{
				MDedge2_del.addObject();
				MDedge2_del.setValue(EMDL_PIPELINE_EDGE_PROCESS,  processList[i].name);
				MDedge2_del.setValue(EMDL_PIPELINE_EDGE_TO, nodeList[outputNode].name);
	    	}
    	}
    }
    MDedge2.write(fh);
    if (fn_del != "")
        MDedge2_del.write(fh_del);

#ifdef DEBUG
    MDedge2.write(std::cerr);
#endif

    fh.close();
    if (fn_del != "")
    	fh_del.close();

	if (do_lock)
	{

#ifdef DEBUG_LOCK
		std::cerr << " write pipeline: now deleting " << fn_lock << std::endl;
#endif

		if (!exists(fn_lock))
			REPORT_ERROR("ERROR: PipeLine::write was expecting a file called "+fn_lock+ " but it is no longer there.");
		if (std::remove(fn_lock.c_str()))
			REPORT_ERROR("ERROR: PipeLine::write reported error in removing file "+fn_lock);
		if (rmdir(dir_lock.c_str()))
			REPORT_ERROR("ERROR: PipeLine::write reported error in removing directory "+dir_lock);
	}

	// Touch a file to indicate to the GUI that the pipeline has just changed
	touch(PIPELINE_HAS_CHANGED);

}

std::string PipeLineFlowChart::getDownwardsArrowLabel(PipeLine &pipeline, long int lower_process, long int upper_process)
{
	// What is the type of the node between upper_process and lower_process?
	bool is_found = false;
	long int mynode = -1;
	for (long int i = 0; i < pipeline.processList[lower_process].inputNodeList.size(); i++)
	{
		long int inode= pipeline.processList[lower_process].inputNodeList[i];
		// Find this one in the outputNodeList of the upper_process
		if (pipeline.nodeList[inode].outputFromProcess == upper_process)
		{
			is_found = true;
			mynode = inode;
			break;
		}
	}

	if (!is_found)
		REPORT_ERROR("PipeLineFlowChart::getDownwardsArrowLabel ERROR: cannot find node connecting " + pipeline.processList[upper_process].name + " -> " + pipeline.processList[lower_process].name);

	std::string mylabel = "";
	MetaDataTable MD;
	long int nr_obj;

	switch (pipeline.nodeList[mynode].type)
	{
		case NODE_MOVIES:
		{
			nr_obj = MD.read(pipeline.nodeList[mynode].name, "", NULL, "", true); // true means: only count nr entries;
			mylabel = integerToString(nr_obj) + " movies";
			break;
		}
		case NODE_MICS:
		{
			nr_obj = MD.read(pipeline.nodeList[mynode].name, "", NULL, "", true); // true means: only count nr entries;
			mylabel = integerToString(nr_obj) + " micrographs";
			break;
		}
		case NODE_PART_DATA:
		{
			nr_obj = MD.read(pipeline.nodeList[mynode].name, "", NULL, "", true); // true means: only count nr entries;
			mylabel = integerToString(nr_obj) + " particles";
			break;
		}
		case NODE_2DREFS:
		{
			mylabel = "2Drefs";
			break;
		}
		case NODE_3DREF:
		{
			mylabel = "3D ref";
			break;
		}
		case NODE_MASK:
		{
			mylabel = "mask";
			break;
		}
		case NODE_MODEL:
		{
			nr_obj = MD.read(pipeline.nodeList[mynode].name, "model_classes", NULL, "", true); // true means: only count nr entries;
			mylabel = integerToString(nr_obj) + " classes";
			break;
		}
		case NODE_OPTIMISER:
		{
			mylabel = "continue";
			break;
		}
		case NODE_HALFMAP:
		{
			mylabel = "half-map";
			break;
		}
		case NODE_FINALMAP:
		{
			mylabel = "final map";
			break;
		}
		case NODE_RESMAP:
		{
			mylabel = "local-res map";
			break;
		}
		default:
		{
			mylabel = "";
			break;
		}
	}

	return mylabel;

}
void PipeLineFlowChart::adaptNamesForTikZ(FileName &name)
{
	name.replaceAllSubstrings((std::string)"_", (std::string)"\\_");
	name.replaceAllSubstrings((std::string)".", (std::string)"-");
	name.replaceAllSubstrings((std::string)",", (std::string)"-");
	name.replaceAllSubstrings((std::string)"^", (std::string)"\\textasciicircum ");
	name.replaceAllSubstrings((std::string)"~", (std::string)"\\textasciitilde ");
}

long int PipeLineFlowChart::addProcessToUpwardsFlowChart(std::ofstream &fh, PipeLine &pipeline,
		long int lower_process, long int new_process, std::vector<long int> &branched_procs)
{


	branched_procs.clear();
	FileName procname;
	if (pipeline.processList[new_process].alias != "None")
		procname = pipeline.processList[new_process].alias;
	else
    	procname = pipeline.processList[new_process].name;

	if (do_short_names)
		procname = procname.beforeFirstOf("/");
	else
	{
		FileName longname = (procname.afterFirstOf("/")).beforeLastOf("/");
		adaptNamesForTikZ(longname);
		procname = procname.beforeFirstOf("/") + "\\\\" + longname;
	}

	FileName new_nodename= pipeline.processList[new_process].name;
	adaptNamesForTikZ(new_nodename);
	FileName lower_nodename;

	// First put the box of the process
	// If this is the lowest process, don't use "above-of" statement, and don't draw an arrow
	if (lower_process < 0)
	{
			fh << "\\node [block] (" << new_nodename << ") {" << procname << "};" << std::endl;
	}
	else
	{
		lower_nodename = pipeline.processList[lower_process].name;
		adaptNamesForTikZ(lower_nodename);

		fh << "\\node [block, above of="<< lower_nodename <<"] (" << new_nodename << ") {" << procname << "};" << std::endl;
		std::string mylabel = getDownwardsArrowLabel(pipeline, lower_process, new_process);
		// Make an arrow from the box to the node it came from
		fh << "\\path [line] ("<< new_nodename <<") -- node[right] {" << mylabel << "} ("<< lower_nodename <<");" << std::endl;
	}

	// See if there are any branchings side-wards, e.g. masks, 2D/3D references, coords, model, optimiser, etc
	long int result = -1;
	if (pipeline.processList[new_process].inputNodeList.size() == 0)
	{
		// Reached the top of the tree!
		return -1;
	}
	if (pipeline.processList[new_process].inputNodeList.size() > 1)
	{

		std::string rightname, leftname;
		for (int inode = 0; inode < pipeline.processList[new_process].inputNodeList.size(); inode++)
		{
			bool is_left = false;
			bool is_right = false;
			bool is_upper_left = false;
			bool is_upper_right = false;
			std::string right_label="", left_label="";

			long int inputnode = pipeline.processList[new_process].inputNodeList[inode];
			int mynodetype = pipeline.nodeList[inputnode].type;

	        if (pipeline.processList[new_process].type == PROC_AUTOPICK)
	        {
	        	is_right = (mynodetype == NODE_2DREFS);
	        	right_label="2D refs";
	        }
	        else if (pipeline.processList[new_process].type == PROC_EXTRACT)
	        {

	        	// If the coordinates come from NODE_MIC_COORDS, then straight up is the CTF info
	        	// If the coordinates come from NODE_PART_DATA, then that should be straight up
	        	// therefore first check whether this node has NODE_PART_DATA input
	        	bool has_part_data = false;
	        	for (int inode2 = 0; inode2 < pipeline.processList[new_process].inputNodeList.size(); inode2++)
	        	{
	    			long int inputnode2 = pipeline.processList[new_process].inputNodeList[inode2];
	    			if (pipeline.nodeList[inputnode2].type == NODE_PART_DATA)
	    			{
	    				has_part_data = true;
	    				break;
	    			}
	        	}
	        	if (has_part_data)
	        	{
	        		is_right = (mynodetype == NODE_MICS);
	        		right_label = "mics";
	        	}
	        	else
	        	{
					is_right = (mynodetype == NODE_MIC_COORDS);
					right_label = "coords";
	        	}
	        }
	        else if (pipeline.processList[new_process].type == PROC_3DCLASS)
	        {
	        	is_right = (mynodetype == NODE_3DREF);
	        	right_label = "3D ref";
	        	is_left = (mynodetype == NODE_MASK);
	        	left_label = "mask";
	        }
	        else if (pipeline.processList[new_process].type == PROC_3DAUTO)
	        {
	        	is_right = (mynodetype == NODE_3DREF);
	        	right_label = "3D ref";
	        	is_left = (mynodetype == NODE_MASK);
	        	left_label = "mask";
	        }
	        else if (pipeline.processList[new_process].type == PROC_JOINSTAR)
	        {
	        	// For joinstar: there will be no parent process that returns a postive value!
	        	// Thereby, joinstar will always end in the 2-4 input processes, each of for which a new flowchart will be made on a new tikZpicture
	        	if (mynodetype == NODE_MOVIES)
	        		right_label = left_label = "mics";
	        	else if (mynodetype == NODE_PART_DATA)
	        		right_label = left_label = "parts";
	        	is_right = (inode == 0);
	        	is_left = (inode == 1);
	        	is_upper_right = (inode == 2);
	        	is_upper_left = (inode == 3);
	        }
	        else if (pipeline.processList[new_process].type == PROC_SUBTRACT)
	        {
	        	is_right = (mynodetype == NODE_3DREF);
	        	right_label = "3D ref";
	        	is_left = (mynodetype == NODE_MASK);
	        	left_label = "mask";
	        }
	        else if (pipeline.processList[new_process].type == PROC_POST)
	        {
	        	is_left = (mynodetype == NODE_MASK);
	        	left_label = "mask";
	        }
	        else if (pipeline.processList[new_process].type == PROC_RESMAP)
	        {
	        	is_left = (mynodetype == NODE_MASK);
	        	left_label = "mask";
	        }

			if (is_right || is_left || is_upper_right || is_upper_left)
			{
				FileName hyperrefname;
				FileName parent_nodename, newprocname;
				long int parent_process = pipeline.nodeList[inputnode].outputFromProcess;
				if (parent_process < 0)
				{
					std::cout << " WARNING: cannot get parent of node: " << pipeline.nodeList[inputnode].name << std::endl;
					parent_nodename = (is_right || is_upper_right) ? new_nodename + "_rigth" : new_nodename + "_left";
					newprocname = "unknown";
				}
				else
				{
					// Keep track of all the side-wards branches
					branched_procs.push_back(parent_process);
					if (pipeline.processList[parent_process].alias != "None")
						newprocname = pipeline.processList[parent_process].alias;
					else
						newprocname = pipeline.processList[parent_process].name;
					if (do_short_names)
						newprocname = newprocname.beforeFirstOf("/");
					else
					{
						FileName longname2 = (newprocname.afterFirstOf("/")).beforeLastOf("/");
						adaptNamesForTikZ(longname2);
						hyperrefname = "sec:" + newprocname.beforeFirstOf("/") + "/" + longname2;
						if (pipeline.processList[parent_process].type==PROC_IMPORT)
							newprocname = newprocname.beforeFirstOf("/") + "\\\\" + longname2;
						else
							newprocname = " \\hyperlink{" + hyperrefname + "}{" + newprocname.beforeFirstOf("/") + "}\\\\" + longname2;
					}

					parent_nodename = pipeline.processList[parent_process].name;
					adaptNamesForTikZ(parent_nodename);
					std::string labelpos;
					if (is_right)
						rightname = parent_nodename;
					if (is_left)
						leftname = parent_nodename;
				}

				if (is_right || is_left)
				{
					std::string pos = (is_right) ? "right" : "left";
					fh << "\\node [block2, "<< pos<<" of="<< new_nodename<<"] (" << parent_nodename << ") {" << newprocname << "};" << std::endl;
				}
				else if (is_upper_right || is_upper_left)
				{
					std::string abovename = (is_upper_right) ? rightname : leftname;
					fh << "\\node [block2b, above of="<< abovename <<"] (" << parent_nodename << ") {" << newprocname << "};" << std::endl;
				}

				// Make an arrow from the box to the process it came from
				std::string arrowlabel = (is_right || is_upper_right) ? right_label : left_label;
				fh << "\\path [line] ("<< parent_nodename <<") -- node[above] {" << arrowlabel << "} ("<< new_nodename <<");" << std::endl;
			}
			else
				result = pipeline.nodeList[inputnode].outputFromProcess;
		}

		return result;
	}
	else
	{
		// Only a single input node: return the process that one came from
		long int inputnode = pipeline.processList[new_process].inputNodeList[0];
		return pipeline.nodeList[inputnode].outputFromProcess;
	}


}

void PipeLineFlowChart::makeOneUpwardsFlowChart(std::ofstream &fh, PipeLine &pipeline, long int from_process,
		std::vector<long int> &all_branches, bool is_main_flow)
{

	openTikZPicture(fh, is_main_flow);
	long int prev_process = -1;
	long int current_process = from_process;
	bool do_stop = false;
	int counter = 0;
	while (!do_stop)
	{

		std::vector<long int> branched_procs;
		long int next_process = addProcessToUpwardsFlowChart(fh, pipeline, prev_process, current_process, branched_procs);

		if (counter > 10)
		{
			closeTikZPicture(fh, false);
			counter = 0;
			next_process = current_process;
			current_process = prev_process;
			prev_process = -1;
			openTikZPicture(fh, false);
		}

		if (next_process < 0)
		{
			do_stop = true;
		}
		else
		{
			prev_process = current_process;
			current_process = next_process;
		}

		// See if there are any new branches, and if so add them to the all_branches vector
		if (do_branches)
		{
			for (int i = 0; i <  branched_procs.size(); i++)
			{
				long int mybranch = branched_procs[i];
				bool already_exists= false;
				for (int j = 0; j < all_branches.size(); j++)
				{
					if (all_branches[j] == mybranch)
					{
						already_exists = true;
						break;
					}
				}
				if (!already_exists && pipeline.processList[mybranch].type != PROC_IMPORT)
				{
					all_branches.push_back(mybranch);
				}
			}
		}

		counter++;

	}
	closeTikZPicture(fh, is_main_flow);

}
void PipeLineFlowChart::makeAllUpwardsFlowCharts(FileName &fn_out, PipeLine &pipeline, long int from_process)
{
	std::ofstream fh;
	openFlowChartFile(fn_out, fh);

	// At the beginning of the flowchart file, first make an overview flowchart with short names
	do_short_names = true;
	do_branches = false;
	FileName myorititle = (pipeline.processList[from_process].alias != "None") ?
			pipeline.processList[from_process].alias : pipeline.processList[from_process].name;
	myorititle=myorititle.beforeLastOf("/");
	adaptNamesForTikZ(myorititle);
	fh << "\\section*{Overview flowchart for " << myorititle << "}" << std::endl;
	std::vector<long int> dummy;
	makeOneUpwardsFlowChart(fh, pipeline, from_process, dummy, true);

	// Then, make fully branched flowcharts below
	do_short_names = false;
	do_branches = true;
	std::vector<long int> all_branches;
	int i = 0;
	all_branches.push_back(from_process);
	while (i < all_branches.size())
	{
		FileName mytitle = (pipeline.processList[all_branches[i]].alias != "None") ? pipeline.processList[all_branches[i]].alias : pipeline.processList[all_branches[i]].name;
		mytitle=mytitle.beforeLastOf("/");
		adaptNamesForTikZ(mytitle);
		if (i == 0)
		{
			std::cout << " Making main branched flowchart ... " <<std::endl;
			fh << "\\section*{Branched flowchart for " << mytitle << "}" << std::endl;
		}
		else
		{
			std::cout << " Making flowchart for branch: " << integerToString(i) << " ... " << std::endl;
			std::string hypertarget = "sec:" + mytitle;
			fh << "\\subsection*{Flowchart for branch " << integerToString(i)<< ": "<< mytitle << "\\hypertarget{"<<hypertarget<<"}{}}" << std::endl;
		}

		makeOneUpwardsFlowChart(fh, pipeline, all_branches[i], all_branches, (i==0) );

		i++;
	}


	closeFlowChartFile(fh);

}

void PipeLineFlowChart::openTikZPicture(std::ofstream &fh, bool is_main_flow)
{
	if (is_main_flow)
	{
		fh << "% For large flowcharts: try reducing the fraction on the next line." << std::endl;
		fh << "\\resizebox{!}{0.75\\textheight}{" << std::endl;
	}
	fh << "\\begin{tikzpicture}[scale=1, auto]" << std::endl;
	// Override the long-name styles with the shorter ones
	if (do_short_names)
	{
		fh << "\\tikzstyle{block} = [rectangle, draw, fill=white,text width=2.5cm, node distance = 1.6cm, text centered, rounded corners, minimum height=0.8cm]" << std::endl;
		fh << "\\tikzstyle{block2} = [rectangle, draw, fill=white,text width=2.5cm, node distance = 4cm, text centered, rounded corners, minimum height=0.8cm]" << std::endl;
		fh << "\\tikzstyle{block2b} = [rectangle, draw, fill=white,text width=2.5cm, node distance = 1.6cm, text centered, rounded corners, minimum height=0.8cm]" << std::endl;
	}
}

void PipeLineFlowChart::closeTikZPicture(std::ofstream &fh, bool is_main_flow)
{
	fh << "\\end{tikzpicture}" << std::endl;
	if (is_main_flow)
	{
		fh << "% For large flowcharts: close resizebox here..." << std::endl;
		fh << "}" << std::endl; // closes resizebox
	}
}

void PipeLineFlowChart::openFlowChartFile(FileName &fn_out, std::ofstream &fh)
{

	fh.open((fn_out).c_str(), std::ios::out);
	if (!fh)
		REPORT_ERROR( (std::string)"PipeLineFlowChart ERROR: Cannot write to file: " + fn_out);

	// Set up the LaTex header
	fh << "\\documentclass{article}" << std::endl;
	fh << "\\usepackage{tikz,hyperref}" << std::endl;
	fh << "\\usetikzlibrary{shapes,arrows}" << std::endl;
	fh << "\\begin{document}" << std::endl;
	// These are the styles for the long names!
	fh << "\\tikzstyle{block} = [rectangle, draw, fill=white,text width=3.5cm, node distance = 1.8cm, text centered, rounded corners]" << std::endl;
	fh << "\\tikzstyle{block2} = [rectangle, draw, fill=blue!20,text width=3.5cm, node distance = 5cm, text centered, rounded corners]" << std::endl;
	fh << "\\tikzstyle{block2b} = [rectangle, draw, fill=blue!20,text width=3.5cm, node distance = 1.8cm, text centered, rounded corners]" << std::endl;


	fh << "\\tikzstyle{line} = [draw, very thick, color=black!50, -latex']" << std::endl << std::endl;
}

void PipeLineFlowChart::closeFlowChartFile(std::ofstream &fh)
{
	fh << "\\end{document}" << std::endl;
	fh.close();
}
