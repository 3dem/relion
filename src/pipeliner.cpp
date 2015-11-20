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

#define DEBUG

std::string Node::printType()
{
    switch ( type )
    {
		case NODE_3DREF: 		return "3D-reference"; break;
		case NODE_2DREF:		return "2D-reference"; break;
		case NODE_HALFMAP:		return "Half-map"; break;
		case NODE_2DMASK:		return "2D-mask"; break;
		case NODE_3DMASK:		return "3D-mask"; break;
		case NODE_MOVIE:		return "movie"; break;
		case NODE_2DMIC:		return "micrograph"; break;
		case NODE_TOMO:			return "tomogram"; break;
		case NODE_MIC_COORDS:	return "Micrograph coords"; break;
		case NODE_MIC_CTFS:		return "Micrograph CTF-info"; break;
		case NODE_PART_DATA:	return "Metadata List of particles"; break;
		case NODE_OPTIMISER:	return "Optimiser"; break;
		default: return "unrecognized";
    }

}


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
			std::string nodename = (nodeList[myNode]).name;
			if (nodename.find(processList[i].name))
			{
#ifdef DEBUG
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

void PipeLine::eraseNode(int ipos)
{
	REPORT_ERROR("PipeLine::eraseNode needs to be re-implemented with vector indices instead of pointers!");

	/*
	Node* eraseNode = &nodeList[ipos];
	nodeList.erase(nodeList.begin()+ipos);
	// Erase all references to this node in the inputNodeList and outputNodeList of all Processes
	for (int i = 0; i < processList.size(); i++)
	{
		for (int j = 0; j < (processList[i]).inputNodeList.size(); j++)
			if (((processList[i]).inputNodeList)[j] == eraseNode)
				(processList[i]).inputNodeList.erase((processList[i]).inputNodeList.begin()+j);
		for (int j = 0; j < (processList[i]).outputNodeList.size(); j++)
			if (((processList[i]).outputNodeList)[j] == eraseNode)
				(processList[i]).outputNodeList.erase((processList[i]).outputNodeList.begin()+j);
	}
	*/
}

void PipeLine::eraseProcess(int ipos)
{
	REPORT_ERROR("PipeLine::eraseProcess needs to be re-implemented with vector indices instead of pointers!");

	/*
	Process* eraseProcess = &processList[ipos];
	processList.erase(processList.begin()+ipos);
	// Erase all references to this process in the inputNodeList and outputNodeList of all Processes
	for (int i = 0; i < nodeList.size(); i++)
	{
		for (int j = 0; j < (nodeList[i]).inputForProcessList.size(); j++)
			if (((nodeList[i]).inputForProcessList)[j] == eraseProcess)
				(nodeList[i]).inputForProcessList.erase((nodeList[i]).inputForProcessList.begin()+j);
		if ((nodeList[i]).outputFromProcess == eraseProcess)
			nodeList[i].outputFromProcess = NULL;
	}
	*/

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

void PipeLine::makeNodeDirectory()
{
	// Clear existing directory
	FileName fn_dir = ".Nodes/";
	std::string command = " rm -rf " + fn_dir;
	std::cerr << command << std::endl;
	int res = system(command.c_str());

	for (long int i = 0; i < nodeList.size(); i++)
	{
		FileName fnt = nodeList[i].name;
		if (exists(fnt))
		{
			// Make subdirectory for each type of node
			FileName fn_type = integerToString(nodeList[i].type) + "/";
			command = "mkdir -p " + fn_dir + fn_type + fnt.substr(0, fnt.rfind("/") + 1);
			std::cerr << command << std::endl;
			res = system(command.c_str());
			command = "touch " + fn_dir + fn_type + fnt;
			std::cerr << command << std::endl;
			res = system(command.c_str());
		}
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
				FileName fnt = nodeList[myNode].name;
				if (!exists(fnt))
				{
					all_exist = false;
					break;
				}
			}
			if (all_exist)
			{
#ifdef DEBUG
				std::cerr << " Updating status of job " << processList[i].name << " to finished " << std::endl;
#endif
				processList[i].status = PROC_FINISHED;
			}
		}
	}

}

void PipeLine::read()
{

	FileName fn = name + "_pipeline.star";
	std::ifstream in(fn.c_str(), std::ios_base::in);
    if (in.fail())
        REPORT_ERROR( (std::string) "PipeLine::read: File " + fn + " cannot be read." );

    MetaDataTable MDnode, MDproc, MDedge1, MDedge2;
    if (!MDnode.readStar(in, "pipeline_nodes"))
    	REPORT_ERROR("PipeLine::read: cannot find table pipeline_nodes in " + name);

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

    if (!MDproc.readStar(in, "pipeline_processes"))
    	REPORT_ERROR("PipeLine::read: cannot find table pipeline_processes in " + name);

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
    if (!MDedge1.readStar(in, "pipeline_input_edges"))
    	REPORT_ERROR("PipeLine::read: cannot find table pipeline_input_edges in " + name);

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
    if (!MDedge2.readStar(in, "pipeline_output_edges"))
    	REPORT_ERROR("PipeLine::read: cannot find table pipeline_output_edges in " + name);

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

}

void PipeLine::write()
{
    std::ofstream  fh;
    FileName fn = name + "_pipeline.star";
    fh.open(fn.c_str(), std::ios::out);

    MetaDataTable MDnode, MDproc, MDedge1, MDedge2;
#ifdef DEBUG
    std::cerr << " writing pipeline as " << fn << std::endl;
#endif

    MDproc.setName("pipeline_processes");
    for(long int i=0 ; i < processList.size() ; i++)
    {
    	MDproc.addObject();
    	MDproc.setValue(EMDL_PIPELINE_PROCESS_NAME, processList[i].name);
    	MDproc.setValue(EMDL_PIPELINE_PROCESS_TYPE, processList[i].type);
    	MDproc.setValue(EMDL_PIPELINE_PROCESS_STATUS, processList[i].status);
    }
#ifdef DEBUG
    MDproc.write(std::cerr);
#endif
    MDproc.write(fh);

    MDnode.setName("pipeline_nodes");
    for(long int i=0 ; i < nodeList.size() ; i++)
    {
    	MDnode.addObject();
    	MDnode.setValue(EMDL_PIPELINE_NODE_NAME, nodeList[i].name);
    	MDnode.setValue(EMDL_PIPELINE_NODE_TYPE, nodeList[i].type);
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
    		MDedge1.addObject();
    		MDedge1.setValue(EMDL_PIPELINE_EDGE_FROM, nodeList[inputNode].name);
    		MDedge1.setValue(EMDL_PIPELINE_EDGE_PROCESS, processList[i].name);
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
    		MDedge2.addObject();
    		MDedge2.setValue(EMDL_PIPELINE_EDGE_PROCESS,  processList[i].name);
    		MDedge2.setValue(EMDL_PIPELINE_EDGE_TO, nodeList[outputNode].name);
    	}
    }
    MDedge2.write(fh);
#ifdef DEBUG
    MDedge2.write(std::cerr);
#endif

    fh.close();

}


