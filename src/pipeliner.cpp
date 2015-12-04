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

//#define DEBUG

long int PipeLine::addNode(Node &_Node)
{

	if (_Node.name=="")
		REPORT_ERROR("PipeLine::addNode ERROR: Adding an empty nodename. Did you fill in all Node named correctly?");
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
	return i;

}


void PipeLine::addNewInputEdge(Node &_Node, long int myProcess)
{

	// 1. Check whether Node with that name already exists in the Node list
	long int oldsize = nodeList.size();
	long int  myNode = addNode(_Node);
	long int newsize = nodeList.size();

	// 2. Set the input_for_process in the inputForProcessList of this Node
	(nodeList[myNode]).inputForProcessList.push_back(myProcess);
	(processList[myProcess]).inputNodeList.push_back(myNode);

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

	// 1. Check whether Node with that name already exists in the Node list
	long int myNode = addNode(_Node);

	// 1. Set the output_from_process in the inputForProcessList of this Node
	nodeList[myNode].outputFromProcess = myProcess;
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
			break;
		}
	}
	if (!is_found)
	{
		processList.push_back(_Process);
	}
	else if (!do_overwrite)
	{
		REPORT_ERROR("PipeLine::addNewProcess: ERROR: trying to add existing Process to the pipeline, while overwriting is not allowed.");
	}
	return i;
}

void PipeLine::deleteProcess(int ipos, bool recursive)
{

	std::vector<bool> deleteProcesses, deleteNodes;
	deleteProcesses.resize(processList.size(), false);
	deleteNodes.resize(nodeList.size(), false);

	std::vector<long int> to_delete_processes;
	to_delete_processes.push_back(ipos);

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
			for (size_t inode = 0; inode < (processList[idel]).outputNodeList.size(); inode++)
			{
				long int mynode = (processList[ipos]).outputNodeList[inode];
				deleteNodes[mynode] = true;
				is_done = true;
				if (recursive)
				{
					// Check whether this node is being used as input for another process, and if so, delete those as well
					for (size_t ii = 0; ii < (nodeList[inode]).inputForProcessList.size(); ii++)
					{
						long int iproc = (nodeList[inode]).inputForProcessList[ii];
						to_delete_processes.push_back(iproc);
						is_done = false;
					}
				}
			}
		}
	}

	// Write new pipeline to disc and read in again
	write(deleteNodes, deleteProcesses);
	read();

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

bool PipeLine::touchTemporaryNodeFile(Node &node, bool touch_even_if_not_exist)
{
	FileName fn_dir = ".Nodes/";
	FileName fnt = node.name;
	if (exists(fnt) || touch_even_if_not_exist)
	{
		// Make subdirectory for each type of node
		FileName fn_type = integerToString(node.type) + "/";
		std::string command = "mkdir -p " + fn_dir + fn_type + fnt.substr(0, fnt.rfind("/") + 1);
		int res = system(command.c_str());
		command = "touch " + fn_dir + fn_type + fnt;
		res = system(command.c_str());
		return true;
	}
	else
		return false;
}

void PipeLine::makeNodeDirectory()
{
	// Clear existing directory
	FileName fn_dir = ".Nodes/";
	std::string command = " rm -rf " + fn_dir;
	int res = system(command.c_str());

	for (long int i = 0; i < nodeList.size(); i++)
	{
		touchTemporaryNodeFile(nodeList[i]);
	}

}

void PipeLine::checkProcessCompletion()
{
	for (long int i=0; i < processList.size(); i++)
	{
		// Only check running processes for file existance
		if (processList[i].status == PROC_RUNNING)
		{
			bool all_exist = true;
			for (long int j = 0; j < processList[i].outputNodeList.size(); j++)
			{
				int myNode = (processList[i]).outputNodeList[j];
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
		else if (processList[i].status == PROC_SCHEDULED)
		{
			// Also touch the output nodes of the scheduled jobs (even though they don't exist yet. That way can link future jobs together!)
			for (long int j = 0; j < processList[i].outputNodeList.size(); j++)
			{
				int myNode = (processList[i]).outputNodeList[j];
				touchTemporaryNodeFile(nodeList[myNode], true); // true means touch even if output node does not exist yet
			}
		}
	}

}

void PipeLine::read()
{

	// Start from scratch
	clear();

	FileName fn = name + "_pipeline.star";
	std::ifstream in(fn.c_str(), std::ios_base::in);
    if (in.fail())
        REPORT_ERROR( (std::string) "PipeLine::read: File " + fn + " cannot be read." );

    MetaDataTable MDnode, MDproc, MDedge1, MDedge2;
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
		std::string name;
		int type, status;
		if (!MDproc.getValue(EMDL_PIPELINE_PROCESS_NAME, name) ||
			!MDproc.getValue(EMDL_PIPELINE_PROCESS_TYPE, type) ||
			!MDproc.getValue(EMDL_PIPELINE_PROCESS_STATUS, status)	)
			REPORT_ERROR("PipeLine::read: cannot find name or type in pipeline_processes table");

		Process newProcess(name, type, status);
		processList.push_back(newProcess);
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
		if (myProcess < 0)
			REPORT_ERROR("PipeLine::read ERROR: cannot find to-process with name: " + procname);
		long int fromNode = findNodeByName(fromnodename);
		if (fromNode < 0)
			REPORT_ERROR("PipeLine::read ERROR: cannot find from-node with name: " + fromnodename);
		processList[myProcess].inputNodeList.push_back(fromNode);
		nodeList[fromNode].inputForProcessList.push_back(myProcess);
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
		if (myProcess < 0)
			REPORT_ERROR("PipeLine::read ERROR: cannot find from-process with name: " + procname);
		long int toNode = findNodeByName(tonodename);
		if (toNode < 0)
			REPORT_ERROR("PipeLine::read ERROR: cannot find to-node with name: " + tonodename);

		processList[myProcess].outputNodeList.push_back(toNode);
		nodeList[toNode].outputFromProcess = myProcess;
	}

	// Close file handler
	in.close();

	// Re-make the .Node directory
	makeNodeDirectory();

}

void PipeLine::write(std::vector<bool> &deleteNode, std::vector<bool> &deleteProcess)
{
    std::ofstream  fh;
    FileName fn = name + "_pipeline.star";
    fh.open(fn.c_str(), std::ios::out);

    bool do_delete = (deleteProcess.size() == processList.size());
    if (do_delete && (deleteNode.size() != nodeList.size()) )
    	REPORT_ERROR("PipeLine::write BUG: not enough entries in deleteNode vector!");

    MetaDataTable MDnode, MDproc, MDedge1, MDedge2;
#ifdef DEBUG
    std::cerr << " writing pipeline as " << fn << std::endl;
#endif

    MDproc.setName("pipeline_processes");
    for(long int i=0 ; i < processList.size() ; i++)
    {
    	if (!do_delete || !deleteProcess[i])
    	{
			MDproc.addObject();
			MDproc.setValue(EMDL_PIPELINE_PROCESS_NAME, processList[i].name);
			MDproc.setValue(EMDL_PIPELINE_PROCESS_TYPE, processList[i].type);
			MDproc.setValue(EMDL_PIPELINE_PROCESS_STATUS, processList[i].status);
    	}
    }
#ifdef DEBUG
    MDproc.write(std::cerr);
#endif
    MDproc.write(fh);

    MDnode.setName("pipeline_nodes");
    for(long int i=0 ; i < nodeList.size() ; i++)
    {
    	if (!do_delete || !deleteNode[i])
    	{
			MDnode.addObject();
			MDnode.setValue(EMDL_PIPELINE_NODE_NAME, nodeList[i].name);
			MDnode.setValue(EMDL_PIPELINE_NODE_TYPE, nodeList[i].type);
    	}
    }
#ifdef DEBUG
    MDnode.write(std::cerr);
#endif
    MDnode.write(fh);

    // Also write all (Node->Process) edges to a single table
    MDedge1.setName("pipeline_input_edges");
    for(long int i=0 ; i < processList.size() ; i++)
    {
		for (long int j=0; j < ((processList[i]).inputNodeList).size(); j++)
		{
			long int inputNode = ((processList[i]).inputNodeList)[j];
	    	if (!do_delete || (!deleteProcess[i] && !deleteNode[inputNode]) )
	    	{
				MDedge1.addObject();
				MDedge1.setValue(EMDL_PIPELINE_EDGE_FROM, nodeList[inputNode].name);
				MDedge1.setValue(EMDL_PIPELINE_EDGE_PROCESS, processList[i].name);
	    	}
    	}
    }
#ifdef DEBUG
    MDedge1.write(std::cerr);
#endif
    MDedge1.write(fh);

    // Also write all (Process->Node) edges to a single table
    MDedge2.setName("pipeline_output_edges");
    for(long int i=0 ; i < processList.size() ; i++)
    {
    	for (long int j=0; j < ((processList[i]).outputNodeList).size(); j++)
    	{
    		long int outputNode = ((processList[i]).outputNodeList)[j];
	    	if (!do_delete || (!deleteProcess[i] && !deleteNode[outputNode]) )
	    	{
				MDedge2.addObject();
				MDedge2.setValue(EMDL_PIPELINE_EDGE_PROCESS,  processList[i].name);
				MDedge2.setValue(EMDL_PIPELINE_EDGE_TO, nodeList[outputNode].name);
	    	}
    	}
    }
    MDedge2.write(fh);
#ifdef DEBUG
    MDedge2.write(std::cerr);
#endif

    fh.close();

}


