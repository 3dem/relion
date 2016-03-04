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

	// 1. Check whether Node with that name already exists in the Node list
	long int old_size = nodeList.size();
	long int myNode = addNode(_Node);

	// 2. Set the output_from_process of this Node
	nodeList[myNode].outputFromProcess = myProcess;

	// 3. Only for new Nodes, add this Node to the outputNodeList of myProcess
	if (myNode == old_size)
		processList[myProcess].outputNodeList.push_back(myNode);

	// Touch .Nodes file, even if it doesn't exist yet for scheduled jobs
	bool touch_if_not_exist = (processList[myProcess].status == PROC_SCHEDULED_CONT || processList[myProcess].status == PROC_SCHEDULED_NEW);
	touchTemporaryNodeFile(nodeList[myNode], touch_if_not_exist);

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

void PipeLine::deleteProcess(int ipos, bool recursive)
{

	FileName fn_del = processList[ipos].name;
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
	write(fn_del, deleteNodes, deleteProcesses);
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
	//if (node.outputFromProcess < 0)
	//	REPORT_ERROR("Pipeline ERROR: node " + node.name + " does not seem to come from any process in the pipeline..." );
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
		std::string command = "mkdir -p " + fn_dir + fn_type + fnt.substr(0, fnt.rfind("/") + 1);
		int res = system(command.c_str());
		command = "touch " + fn_dir + fn_type + fnt;
		res = system(command.c_str());
		command = "chmod 777 " + fn_dir + fn_type + fnt;
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
		int myproc = nodeList[i].outputFromProcess;
		//if (myproc < 0)
		//	REPORT_ERROR("PipeLine::makeNodeDirectory ERROR: cannot get from which process node " + nodeList[i].name + " is coming.");
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
void PipeLine::read(bool only_read_if_file_exists)
{

	// Start from scratch
	clear();

	FileName fn = name + "_pipeline.star";
	std::ifstream in(fn.c_str(), std::ios_base::in);

	if (in.fail())
	{
		if (only_read_if_file_exists)
			return;
		else
			REPORT_ERROR( (std::string) "PipeLine::read: File " + fn + " cannot be read." );
	}

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
				std::string command = " ln -s ../" + name + " " + fn_alias;
				int res= system(command.c_str());

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

	// Re-make the .Node directory
	makeNodeDirectory();

}

void PipeLine::write(FileName fn_del, std::vector<bool> deleteNode, std::vector<bool> deleteProcess)
{
    std::ofstream  fh, fh_del;
    FileName fn = name + "_pipeline.star";
    fh.open(fn.c_str(), std::ios::out);

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
}


