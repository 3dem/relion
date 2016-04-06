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
#include "src/metadata_table.h"
#include <iostream>
#include <sstream>


//forward declaration
class Process;

/*
 * The Node class represents data and metadata that are either input to or output from Processes
 * Nodes are connected to each by Edges:
 * - the fromEdgeList are connections with Nodes earlier (higher up) in the pipeline
 * - the toEdgeList are connections with Nodes later (lower down) in the pipeline
 *
 * Nodes could be of the following types:
 */

#define NODE_MOVIES			0 // 2D micrograph movie(s), e.g. Falcon001_movie.mrcs or micrograph_movies.star
#define NODE_MICS			1 // 2D micrograph(s), possibly with CTF information as well, e.g. Falcon001.mrc or micrographs.star
#define NODE_MIC_COORDS		2 // Suffix for particle coordinates in micrographs (e.g. autopick.star or .box)
#define NODE_PART_DATA		3 // A metadata (STAR) file with particles (e.g. particles.star or run1_data.star)
#define NODE_MOVIE_DATA		4 // A metadata (STAR) file with particle movie-frames (e.g. particles_movie.star or run1_ct27_data.star)
#define NODE_2DREFS       	5 // A STAR file with one or multiple 2D references, e.g. autopick_references.star
#define NODE_3DREF       	6 // A single 3D-reference, e.g. map.mrc
#define NODE_MASK			7 // 3D mask, e.g. mask.mrc or masks.star
#define NODE_MODEL		    8 // A model STAR-file for class selection
#define NODE_OPTIMISER		9 // An optimiser STAR-file for job continuation
#define NODE_HALFMAP		10// Unfiltered half-maps from 3D auto-refine, e.g. run1_half?_class001_unfil.mrc
#define NODE_FINALMAP		11// Sharpened final map from post-processing (cannot be used as input)
#define NODE_RESMAP			12// Resmap with local resolution (cannot be used as input)
#define NODE_PDF_LOGFILE    13//PDF logfile

class Node
{
	public:
	std::string name; // what's my name?
	int type; // which type of node am I?
	std::vector<long int> inputForProcessList; 	  //list of processes that use this Node as input
	long int outputFromProcess;   //Which process made this Node

	// Constructor
	Node(std::string _name, int _type)
	{
		name = _name;
		type = _type;
		outputFromProcess = -1;
	}

	// Destructor
	// Do not delete the adjacent nodes here... They will be deleted by graph destructor
	~Node()
	{
		inputForProcessList.clear();
	}

};

/*
 * The Process class represents tasks/jobs
 * A Process converts input Nodes into output Nodes, thereby generating an Edge between them (so each Edge has an associated Process)
 * A Process connects one or more input Nodes to one or more output Nodes
 * An Edge is added from each input Node to each Output node
 *
 * Processes could be of the following types:
 */


// This order defines the order of the process browser in the GUI!
// TODO:#define PROC_IMPORT      	1 // Import any node into the pipeline (by reading them from external files)
#define PROC_IMPORT         0 // Import any file as a Node of a given type
#define PROC_MOTIONCORR 	1 // Import any file as a Node of a given type
#define PROC_CTFFIND	    2 // Estimate CTF parameters from micrographs for either entire micrographs and/or particles
#define PROC_MANUALPICK		3 // Manually pick particle coordinates from micrographs
#define PROC_AUTOPICK		4 // Automatically pick particle coordinates from micrographs, their CTF and 2D references
#define PROC_EXTRACT		5 // Window particles, normalize, downsize etc from micrographs (also combine CTF into metadata file)
#define PROC_SORT           6 // Sort particles based on their Z-scores
#define PROC_CLASSSELECT    7 // Read in model.star file, and let user interactively select classes through the display (later: auto-selection as well)
#define PROC_2DCLASS		8 // 2D classification (from input particles)
#define PROC_3DCLASS		9 // 3D classification (from input 2D/3D particles, an input 3D-reference, and possibly a 3D mask)
#define PROC_3DAUTO	        10 // 3D auto-refine (from input particles, an input 3Dreference, and possibly a 3D mask)
#define PROC_POLISH			11// Particle-polishing (from movie-particles)
#define PROC_MASKCREATE     12// Process to create masks from input maps
#define PROC_JOINSTAR       13// Process to create masks from input maps
#define PROC_SUBTRACT       14// Process to subtract projections of parts of the reference from experimental images
#define PROC_POST			15// Post-processing (from unfiltered half-maps and a possibly a 3D mask)
#define PROC_RESMAP			16// Local resolution estimation (from unfiltered half-maps and a 3D mask)
#define PROC_MOVIEREFINE    17// Movie-particle extraction and refinement combined
#define NR_BROWSE_TABS      18

// Status a Process may have
#define PROC_RUNNING   0
#define PROC_SCHEDULED_NEW 1
#define PROC_FINISHED  2
#define PROC_SCHEDULED_CONT 3


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

class PipeLine
{

	public:
	int job_counter;

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

	// Delete a process and its output nodes (and all input edges) from the pipeline
	void deleteProcess(int ipos, bool recursive = false);

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

	// Check for process completion by cheking for the presence of all outputNode filenames
	void checkProcessCompletion();

	// Import a job into the pipeline
	// Return true if at least one process is imported correctly
	bool importPipeline(std::string _name);

	// Write out the pipeline to a STAR file
	void write(FileName fn_del="", std::vector<bool> deleteNode = std::vector<bool>(), std::vector<bool> deleteProcess = std::vector<bool>());

	// Read in the pipeline from a STAR file
	void read(bool only_read_if_file_exists=false);

	// Make LaTeX and TikZ-based flowcharts
	void makeUpwardsFlowChart(long int from_process);

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
			std::vector<long int> &all_branches);

	void makeAllUpwardsFlowCharts(FileName &fn_out, PipeLine &pipeline, long int from_process);

	// Open and close a new flowchart picture
	void openTikZPicture(std::ofstream &fh);

	void closeTikZPicture(std::ofstream &fh);

	void adaptNamesForTikZ(FileName &name);

	// Open and close a new output file
	void openFlowChartFile(FileName &fn_out, std::ofstream &fh);

	void closeFlowChartFile(std::ofstream &fh);


};


#endif /* PIPELINER_H_ */
