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
	bool touch_if_not_exist = (processList[myProcess].status == PROC_SCHEDULED_CONT || processList[myProcess].status == PROC_SCHEDULED_NEW);

	// 2. Set the output_from_process of this Node
	_Node.outputFromProcess = myProcess;
	long int myNode = addNode(_Node, touch_if_not_exist);

	// 3. Only for new Nodes, add this Node to the outputNodeList of myProcess
	if (myNode == old_size)
		processList[myProcess].outputNodeList.push_back(myNode);

}

long int PipeLine::addNewProcess(Process &_Process, bool do_overwrite)
{
	// Check whether Process  with the same name already exists in the processList
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

	bool touch_if_not_exist = (process.status == PROC_SCHEDULED_CONT || process.status == PROC_SCHEDULED_NEW);

	for (int j = 0; j < process.outputNodeList.size(); j++)
	{
		long int mynode = process.outputNodeList[j];
		touchTemporaryNodeFile(nodeList[mynode], touch_if_not_exist);
	}

}

void PipeLine::deleteTemporaryNodeFile(Node &node)
{
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

	for (int j = 0; j < process.outputNodeList.size(); j++)
	{
		long int mynode = process.outputNodeList[j];
		deleteTemporaryNodeFile(nodeList[mynode]);
	}

}

void PipeLine::remakeNodeDirectory()
{
	// Clear existing directory
	FileName fn_dir = ".Nodes/";
	std::string command = " rm -rf " + fn_dir;
	int res = system(command.c_str());

	for (long int i = 0; i < nodeList.size(); i++)
	{
		int myproc = nodeList[i].outputFromProcess;
		bool touch_if_not_exist = (myproc < 0) ? false : (processList[myproc].status == PROC_SCHEDULED_CONT ||
								   	   	   	   	   	   	  processList[myproc].status == PROC_SCHEDULED_NEW);
		touchTemporaryNodeFile(nodeList[i], touch_if_not_exist);
	}
	command = "chmod 777 -R " + fn_dir;
	res = system(command.c_str());
}


void PipeLine::checkProcessCompletion()
{

	for (long int i=0; i < processList.size(); i++)
	{
		// Only check running processes for file existence
		if (processList[i].status == PROC_RUNNING)
		{
			bool all_exist = true;
			for (long int j = 0; j < processList[i].outputNodeList.size(); j++)
			{
				int myNode = (processList[i]).outputNodeList[j];
				if (myNode < 0 || myNode >= nodeList.size())
					REPORT_ERROR("pipeline checkProcessCompletion ERROR: " + integerToString(j) + "th output node of " + processList[i].name + " is invalid: " + integerToString(myNode));
				if (!touchTemporaryNodeFile(nodeList[myNode]))
				{
					all_exist = false;
					break;
				}
			}
			if (all_exist)
			{
				processList[i].status = PROC_FINISHED;
			}
		}
	}

	// Don't write out the updated pipeline here.
	// Always make sure to read in the existing pipeline inside gui_mainwindow.cpp before writing a new one to disk
	// This is to make sure two different windows do not get out-of-sync

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
void PipeLine::read(bool do_lock)
{

	FileName fn_lock = ".lock_" + name + "_pipeline.star";
	if (do_lock)
	{
		int iwait =0;
		while( exists(fn_lock) )
		{
			// If the lock exists: wait 3 seconds and try again
			// First time round, print a warning message
			if (iwait == 0)
			{
				std::cout << "WARNING: trying to read pipeline.star, but " << fn_lock << " exists!" << std::endl;
				std::cout << " This is a protection against simultaneous writing to the pipeline by multiple instances of the GUI." << std::endl;
				std::cout << " You can override this by manually deleting the " << fn_lock << " file." << std::endl;
				std::cout << " This instance of the GUI will try for 1 minute to see whether the lock disappears." << std::endl;
			}
			sleep(3);
			iwait++;
			if (iwait > 20)
				REPORT_ERROR("ERROR: PipeLine::read has waited for 1 minute for lock file to disappear. You may want to manually remove the file: " + fn_lock);
		}
		// Generate the lock file
		touch(fn_lock);
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

	FileName fn_lock = ".lock_" + name + "_pipeline.star";
	if (do_lock)
	{
		if (!exists(fn_lock))
			REPORT_ERROR("ERROR: PipeLine::write was expecting a file called "+fn_lock+ " but it isn't there.");
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
		if (!exists(fn_lock))
			REPORT_ERROR("ERROR: PipeLine::write was expecting a file called "+fn_lock+ " but it is no longer there.");
		std::remove(fn_lock.c_str());
	}

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
    case NODE_MOVIE_DATA:
    {
    	mylabel = "particle movie-frames";
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
	name.replaceAllSubstrings((std::string)"_", (std::string)"-");
	name.replaceAllSubstrings((std::string)".", (std::string)"-");
	name.replaceAllSubstrings((std::string)",", (std::string)"-");
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

	        switch (pipeline.processList[new_process].type)
	        {
	        case PROC_AUTOPICK:
	        {
	        	is_right = (mynodetype == NODE_2DREFS);
	        	right_label="2D refs";
	            break;
	        }
	        case PROC_EXTRACT:
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
	        	break;
	        }
	        case PROC_SORT:
	        {
	        	is_right = (mynodetype == NODE_MODEL || mynodetype == NODE_2DREFS);
	        	right_label = "refs";
	        	break;
	        }
	        case PROC_3DCLASS:
	        {
	        	is_right = (mynodetype == NODE_3DREF);
	        	right_label = "3D ref";
	        	is_left = (mynodetype == NODE_MASK);
	        	left_label = "mask";
	        	break;
	        }
	        case PROC_3DAUTO:
	        {
	        	is_right = (mynodetype == NODE_3DREF);
	        	right_label = "3D ref";
	        	is_left = (mynodetype == NODE_MASK);
	        	left_label = "mask";
	        	break;
	        }
	        case PROC_MOVIEREFINE:
	        {
	        	is_right = (mynodetype == NODE_MOVIES);
	        	right_label = "movies";
	        	break;
	        }
	        case PROC_POLISH:
	        {
	        	is_left = (mynodetype == NODE_MASK);
	        	left_label = "mask";
	        	break;
	        }
	        case PROC_JOINSTAR:
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
	        	break;
	        }
	        case PROC_SUBTRACT:
	        {
	        	is_right = (mynodetype == NODE_3DREF);
	        	right_label = "3D ref";
	        	is_left = (mynodetype == NODE_MASK);
	        	left_label = "mask";
	        	break;
	        }
	        case PROC_POST:
	        {
	        	is_left = (mynodetype == NODE_MASK);
	        	left_label = "mask";
	        	break;
	        }
	        case PROC_RESMAP:
	        {
	        	is_left = (mynodetype == NODE_MASK);
	        	left_label = "mask";
	        	break;
	        }
	    	default:
	    	{
	    		break;
	    	}
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
	std::vector<long int> all_branches;



	std::string myorititle;
	all_branches.push_back(from_process);
	int i = 0;
	while (i < all_branches.size())
	{
		FileName mytitle = (pipeline.processList[all_branches[i]].alias != "None") ? pipeline.processList[all_branches[i]].alias : pipeline.processList[all_branches[i]].name;
		mytitle=mytitle.beforeLastOf("/");
		adaptNamesForTikZ(mytitle);
		if (i == 0)
		{
			std::cout << " Making first flowchart ... " <<std::endl;
			fh << "\\section*{Branched flowchart for " << mytitle << "}" << std::endl;
			myorititle = mytitle;
		}
		else
		{
			std::cout << " Making flowchart for branch: " << integerToString(i) << " ... " << std::endl;
			std::string hypertarget = "sec:" + mytitle;
			fh << "\\subsection*{Flowchart for branch " << integerToString(i)<< ": "<< mytitle << "\\hypertarget{"<<hypertarget<<"}{}}" << std::endl;
			//fh << "\\hypertarget{" << hypertarget<<"}{.}"<<std::endl;
			//fh << "\\newline" << std::endl;
		}

		makeOneUpwardsFlowChart(fh, pipeline, all_branches[i], all_branches, (i==0) );

		i++;
	}

	// At the end of the flowchart file, also make one with short names
	do_short_names = true;
	do_branches = false;
	fh << "\\section{Overview flowchart for " << myorititle << "}" << std::endl;
	makeOneUpwardsFlowChart(fh, pipeline,all_branches[0], all_branches, true);

	closeFlowChartFile(fh);

}

void PipeLineFlowChart::openTikZPicture(std::ofstream &fh, bool is_main_flow)
{
	if (is_main_flow)
	{
		fh << "% For large flowcharts: try removing percent sign on next line, and on line below." << std::endl;
		fh << "%\\resizebox{!}{0.95\\textheight}{" << std::endl;
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
		fh << "% For large flowcharts: also remove percent sign on next line." << std::endl;
		fh << "%}" << std::endl; // closes resizebox
	}
}

void PipeLineFlowChart::openFlowChartFile(FileName &fn_out, std::ofstream &fh)
{

    fh.open((fn_out).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR( (std::string)"PipeLineFlowChart ERROR: Cannot write to file: " + fn_out);

    // Set up the LaTex header
    fh << "\\documentclass{article}" << std::endl;
    fh << "\\usepackage{tikz,hyperref, sectsty,xcolor}" << std::endl;
    fh << "\\usetikzlibrary{shapes,arrows}" << std::endl;
    fh << "\\begin{document}" << std::endl;
    fh << "\\subsectionfont{\\color{blue!50}}" << std::endl;
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
