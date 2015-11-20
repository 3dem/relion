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

#define NODE_3DREF       	0 // 3D-reference(s), e.g. map.mrc or refs.star
#define NODE_2DREF          1 // 2D reference(s), e.g. 1@refs.star, ref.mrc or refs.star
#define NODE_HALFMAP		2 // Unfiltered half-maps from 3D auto-refine, e.g. run1_half?_class001_unfil.mrc
#define NODE_2DMASK			3 // 2D mask(s), e.g. mask.mrc or masks.star
#define NODE_3DMASK			4 // 3D mask(s), e.g. mask.mrc or masks.star
#define NODE_MOVIE			6 // 2D micrograph movie(s), e.g. Falcon001_movie.mrcs or micrograph_movies.star
#define NODE_2DMIC			5 // 2D micrograph(s), e.g. Falcon001.mrc or micrographs.star
#define NODE_TOMO			7 // 3D tomogram(s), e.g. tomo001.mrc or tomograms.star
#define NODE_MIC_COORDS		8 // (A list of) coordinates for micrographs (*_autopick.star) // TODO: make a STAR file with all autopick.star names?
#define NODE_MIC_CTFS		9 // (A list of) CTF parameters for micrographs (micrographs_ctf.star)
#define NODE_PART_DATA		10 // A metadata (STAR) file with particles (e.g. particles.star or run1_data.star)
#define NODE_MOVIE_DATA		11 // A metadata (STAR) file with particle movie-frames (e.g. particles_movie.star or run1_ct27_data.star)
#define NODE_OPTIMISER		12 // A optimiser STAR-file for job continuation
#define NODE_FINALMAP		13 // Sharpened final map from post-processing (cannot be used as input)
#define NODE_RESMAP			14 // Resmap with local resolution (cannot be used as input)

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
		outputFromProcess = NULL;
	}

	// Destructor
	// Do not delete the adjacent nodes here... They will be deleted by graph destructor
	~Node()
	{
		inputForProcessList.clear();
	}

	std::string printType();

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
#define PROC_GENERAL        1 // old general tab; TODO: remove!
#define PROC_MANUALPICK		2 // Manually pick particle coordinates from micrographs
// TODO #define PROC_MOTIONCORR	    2 // Motion-correct movies into aligned movies and average micrographs
#define PROC_CTFFIND	    3 // Estimate CTF parameters from micrographs for either entire micrographs and/or particles
#define PROC_AUTOPICK		4 // Automatically pick particle coordinates from micrographs, their CTF and 2D references
#define PROC_SORT           5 // Sort particles based on their Z-scores
#define PROC_EXTRACT		6 // Window particles, normalize, downsize etc from micrographs (also combine CTF into metadata file)
#define PROC_2DCLASS		7 // 2D classification (from input particles)
#define PROC_3DCLASS		8 // 3D classification (from input 2D/3D particles, an input 3D-reference, and possibly a 3D mask)
// TODO #define PROC_CLASSSELECT    9
// TODO #define PROC_MASKCREATE     10
// TODO #define PROC_UNBIN          11
#define PROC_3DAUTO	        9  // 3D auto-refine (from input particles, an input 3Dreference, and possibly a 3D mask)
#define PROC_POLISH			10 // Particle-polishing (from movie-particles)
#define PROC_POST			11 // Post-processing (from unfiltered half-maps and a possibly a 3D mask)
#define PROC_RESMAP			12 // Local resolution estimation (from unfiltered half-maps and a 3D mask)
#define PROC_PUBLISH		13 // TODO: remove! Local resolution estimation (from unfiltered half-maps and a 3D mask)
#define NR_BROWSE_TABS PROC_PUBLISH

// Status a Process may have
#define PROC_RUNNING   0
#define PROC_SCHEDULED 1
#define PROC_FINISHED  2
#define PROC_CANCELLED 3


class Process
{

	public:
	std::string name;
	int type;
	int status;
	std::vector<long int> inputNodeList;  // List of Nodes of input to this process
	std::vector<long int> outputNodeList; // List of Nodes of output from this process

	// Constructor
	Process(std::string _name, int _type, int _status)
	{
		name = _name;
		type = _type;
		status = _status;
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
	std::string name;
	std::vector<Node> nodeList; //list of all Nodes in the pipeline
	std::vector<Process> processList; //list of all Processes in the pipeline

	PipeLine()
	{
		name = "default";
	}

	~PipeLine()
	{
		nodeList.clear();
		processList.clear();
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
	long int addNode(Node &_Node);

	// Add a new Process to the list (no checks are performed)
	long int addNewProcess(Process &_Process, bool do_overwrite = false);

	// Erase nodes and processes from the pipeline
	// At the ipos'th position in the nodeList or processList
	void eraseNode(int ipos);
	void eraseProcess(int ipos);

	// Find nodes or process
	long int findNodeByName(std::string name);
	long int findProcessByName(std::string name);

	// Make empty entries of all NodeNames in a hidden directory (useful for file browsing for InputNode I/O)
	void makeNodeDirectory();

	// Check for process completion by cheking for the presence of all outputNode filenames
	void checkProcessCompletion();

	// Write out the pipeline to a STAR file
	void write();

	// Read in the pipeline from a STAR file
	void read();

};



#endif /* PIPELINER_H_ */
