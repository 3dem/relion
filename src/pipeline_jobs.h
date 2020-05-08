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


#ifndef SRC_PIPELINE_JOBS_H_
#define SRC_PIPELINE_JOBS_H_

#define JOBOPTION_UNDEFINED 0
#define JOBOPTION_ANY 1
#define JOBOPTION_FILENAME 2
#define JOBOPTION_INPUTNODE 3
#define JOBOPTION_RADIO 4
#define JOBOPTION_BOOLEAN 5
#define JOBOPTION_SLIDER 6
#define JOBOPTION_ONLYTEXT 7
#include "src/macros.h"
#include "src/metadata_table.h"
#include "src/filename.h"
#include <string>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#define YSTEP 20

#define TOGGLE_DEACTIVATE 0
#define TOGGLE_REACTIVATE 1
#define TOGGLE_ALWAYS_DEACTIVATE 2
#define TOGGLE_LEAVE_ACTIVE 3


#define HAS_MPI true
#define HAS_NOT_MPI false
#define HAS_THREAD true
#define HAS_NOT_THREAD false

#define RADIO_SAMPLING 0
#define RADIO_NODETYPE 1
#define RADIO_GAIN_ROTATION 2
#define RADIO_GAIN_FLIP 3

// Our own defaults at LMB are the hard-coded ones
#define DEFAULTQSUBLOCATION "/public/EM/RELION/relion/bin/relion_qsub.csh"
#define DEFAULTCTFFINDLOCATION "/public/EM/ctffind/ctffind.exe"
#define DEFAULTMOTIONCOR2LOCATION "/public/EM/MOTIONCOR2/MotionCor2"
#define DEFAULTGCTFLOCATION "/public/EM/Gctf/bin/Gctf"
#define DEFAULTRESMAPLOCATION "/public/EM/ResMap/ResMap-1.1.4-linux64"
#define DEFAULTQSUBCOMMAND "qsub"
#define DEFAULTQUEUENAME "openmpi"
#define DEFAULTMININIMUMDEDICATED 1
#define DEFAULTWARNINGLOCALMPI 32
#define DEFAULTALLOWCHANGEMINDEDICATED true
#define DEFAULTQUEUEUSE false
#define DEFAULTNRMPI 1
#define DEFAULTMPIMAX 64
#define DEFAULTNRTHREADS 1
#define DEFAULTTHREADMAX 16
#define DEFAULTMPIRUN "mpirun"
#define DEFAULTSCRATCHDIR ""

static const std::vector<std::string> job_undefined_options{
	"undefined"
};

static const std::vector<std::string> job_boolean_options{
	"Yes",
	"No"
};

static const std::vector<std::string> job_sampling_options{
	"30 degrees",
	"15 degrees",
	"7.5 degrees",
	"3.7 degrees",
	"1.8 degrees",
	"0.9 degrees",
	"0.5 degrees",
	"0.2 degrees",
	"0.1 degrees"
};

static const std::vector<std::string> job_nodetype_options{
	"Particle coordinates (*.box, *_pick.star)",
	"Particles STAR file (.star)",
	"Movie-particles STAR file (.star)",
	"2D references (.star or .mrcs)",
	"Micrographs STAR file (.star)",
	"3D reference (.mrc)",
	"3D mask (.mrc)",
	"Unfiltered half-map (unfil.mrc)"
};

static const std::vector<std::string> job_gain_rotation_options{
	"No rotation (0)",
	"90 degrees (1)",
	"180 degrees (2)",
	"270 degrees (3)"
};

static const std::vector<std::string> job_gain_flip_options{
	"No flipping (0)",
	"Flip upside down (1)",
	"Flip left to right (2)"
};

static const std::vector<std::string> job_ctffit_options{
	"No",
	"Per-micrograph",
	"Per-particle"
};

static std::string getStringFitOption(std::string option)
{
	if (option == job_ctffit_options[0]) return "f";
	if (option == job_ctffit_options[1]) return "m";
	if (option == job_ctffit_options[2]) return "p";
	REPORT_ERROR("ERROR: unrecognised fit_mode for ctf_refine");
}

// To have a line on the GUI to change the minimum number of dedicated in a job
static bool do_allow_change_minimum_dedicated;

// Optional output file for any jobtype that explicitly defines the output nodes
#define RELION_OUTPUT_NODES "RELION_OUTPUT_NODES.star"

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
//#define NODE_MOVIE_DATA		4 // A metadata (STAR) file with particle movie-frames (e.g. particles_movie.star or run1_ct27_data.star)
#define NODE_2DREFS       	5 // A STAR file with one or multiple 2D references, e.g. autopick_references.star
#define NODE_3DREF       	6 // A single 3D-reference, e.g. map.mrc
#define NODE_MASK			7 // 3D mask, e.g. mask.mrc or masks.star
#define NODE_MODEL		    8 // A model STAR-file for class selection
#define NODE_OPTIMISER		9 // An optimiser STAR-file for job continuation
#define NODE_HALFMAP		10// Unfiltered half-maps from 3D auto-refine, e.g. run1_half?_class001_unfil.mrc
#define NODE_FINALMAP		11// Sharpened final map from post-processing (cannot be used as input)
#define NODE_RESMAP			12// Resmap with local resolution (cannot be used as input)
#define NODE_PDF_LOGFILE    13// PDF logfile
#define NODE_POST           14// Postprocess STAR file (with FSC curve, unfil half-maps, masks etc in it: used by Jasenko's programs
#define NODE_POLISH_PARAMS  15// Txt file with optimal parameters for Bayesian polishing

// All the directory names of the different types of jobs defined inside the pipeline
#define PROC_IMPORT_NAME        "Import"       // Import any file as a Node of a given type
#define PROC_MOTIONCORR_NAME 	"MotionCorr"   // Import any file as a Node of a given type
#define PROC_CTFFIND_NAME	    "CtfFind"  	   // Estimate CTF parameters from micrographs for either entire micrographs and/or particles
#define PROC_MANUALPICK_NAME    "ManualPick"   // Manually pick particle coordinates from micrographs
#define PROC_AUTOPICK_NAME		"AutoPick"     // Automatically pick particle coordinates from micrographs, their CTF and 2D references
#define PROC_EXTRACT_NAME		"Extract"      // Window particles, normalize, downsize etc from micrographs (also combine CTF into metadata file)
#define PROC_CLASSSELECT_NAME   "Select" 	   // Read in model.star file, and let user interactively select classes through the display (later: auto-selection as well)
#define PROC_2DCLASS_NAME 		"Class2D"      // 2D classification (from input particles)
#define PROC_3DCLASS_NAME		"Class3D"      // 3D classification (from input 2D/3D particles, an input 3D-reference, and possibly a 3D mask)
#define PROC_3DAUTO_NAME        "Refine3D"     // 3D auto-refine (from input particles, an input 3Dreference, and possibly a 3D mask)
//#define PROC_POLISH_NAME	    "Polish"       // Particle-polishing (from movie-particles)
#define PROC_MASKCREATE_NAME    "MaskCreate"   // Process to create masks from input maps
#define PROC_JOINSTAR_NAME      "JoinStar"     // Process to create masks from input maps
#define PROC_SUBTRACT_NAME      "Subtract"     // Process to subtract projections of parts of the reference from experimental images
#define PROC_POST_NAME			"PostProcess"  // Post-processing (from unfiltered half-maps and a possibly a 3D mask)
#define PROC_RESMAP_NAME  	    "LocalRes"     // Local resolution estimation (from unfiltered half-maps and a 3D mask)
//#define PROC_MOVIEREFINE_NAME   "MovieRefine"  // Movie-particle extraction and refinement combined
#define PROC_INIMODEL_NAME		"InitialModel" // De-novo generation of 3D initial model (using SGD)
#define PROC_MULTIBODY_NAME		"MultiBody"    // Multi-body refinement
#define PROC_MOTIONREFINE_NAME  "Polish"       // Jasenko's motion fitting program for Bayesian polishing (to replace MovieRefine?)
#define PROC_CTFREFINE_NAME     "CtfRefine"    // Jasenko's program for defocus and beamtilt optimisation
#define PROC_EXTERNAL_NAME      "External"     // For running non-relion programs

#define PROC_IMPORT         0 // Import any file as a Node of a given type
#define PROC_MOTIONCORR 	1 // Import any file as a Node of a given type
#define PROC_CTFFIND	    2 // Estimate CTF parameters from micrographs for either entire micrographs and/or particles
#define PROC_MANUALPICK 	3 // Manually pick particle coordinates from micrographs
#define PROC_AUTOPICK		4 // Automatically pick particle coordinates from micrographs, their CTF and 2D references
#define PROC_EXTRACT		5 // Window particles, normalize, downsize etc from micrographs (also combine CTF into metadata file)
//#define PROC_SORT         6 // Sort particles based on their Z-scores
#define PROC_CLASSSELECT    7 // Read in model.star file, and let user interactively select classes through the display (later: auto-selection as well)
#define PROC_2DCLASS		8 // 2D classification (from input particles)
#define PROC_3DCLASS		9 // 3D classification (from input 2D/3D particles, an input 3D-reference, and possibly a 3D mask)
#define PROC_3DAUTO         10// 3D auto-refine (from input particles, an input 3Dreference, and possibly a 3D mask)
//#define PROC_POLISH  		11// Particle-polishing (from movie-particles)
#define PROC_MASKCREATE     12// Process to create masks from input maps
#define PROC_JOINSTAR       13// Process to create masks from input maps
#define PROC_SUBTRACT       14// Process to subtract projections of parts of the reference from experimental images
#define PROC_POST			15// Post-processing (from unfiltered half-maps and a possibly a 3D mask)
#define PROC_RESMAP 		16// Local resolution estimation (from unfiltered half-maps and a 3D mask)
//#define PROC_MOVIEREFINE    17// Movie-particle extraction and refinement combined
#define PROC_INIMODEL		18// De-novo generation of 3D initial model (using SGD)
#define PROC_MULTIBODY      19// Multi-body refinement
#define PROC_MOTIONREFINE   20// Jasenko's motion_refine
#define PROC_CTFREFINE      21// Jasenko's ctf_refine
#define PROC_EXTERNAL       99// External scripts
#define NR_BROWSE_TABS      20

// Status a Process may have
#define PROC_RUNNING          0 // (hopefully) running
#define PROC_SCHEDULED        1 // scheduled for future execution
#define PROC_FINISHED_SUCCESS 2 // successfully finished
#define PROC_FINISHED_FAILURE 3 // reported an error
#define PROC_FINISHED_ABORTED 4 // aborted by the user

struct gui_layout
{
    /// Name for the tab
    std::string tabname;
    /// y-position
    int ypos;
    ///
    RFLOAT w;
};

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

// Helper function to get the outputnames of refine jobs
std::vector<Node> getOutputNodesRefine(std::string outputname, int iter, int K, int dim, int nr_bodies=1);

// One class to store any type of Option for a GUI entry
class JobOption
{
public:

	std::string label;
	std::string label_gui;
	int joboption_type;
	std::string variable;
	std::string value;
	std::string default_value;
	std::string helptext;
	float min_value;
	float max_value;
	float step_value;
	int node_type;
	std::string pattern;
	std::string directory;
	std::vector<std::string> radio_options;

public:

	// Any constructor
	JobOption(std::string _label, std::string _default_value, std::string _helptext);

	// FileName constructor
	JobOption(std::string _label, std::string  _default_value, std::string _pattern, std::string _directory, std::string _helptext);

	// InputNode constructor
	JobOption(std::string _label, int _nodetype, std::string _default_value, std::string _pattern, std::string _helptext);

	// Radio constructor
	JobOption(std::string _label, std::vector<std::string> radio_options, int ioption,  std::string _helptext);

	// Boolean constructor
	JobOption(std::string _label, bool _boolvalue, std::string _helptext);

	// Slider constructor
	JobOption(std::string _label, float _default_value, float _min_value, float _max_value, float _step_value, std::string _helptext);

	// Write to a STAR file
	void writeToMetaDataTable(MetaDataTable& MD) const;

	// Empty constructor
	JobOption() { clear(); }

	// Empty destructor
	~JobOption() { clear(); }

	void clear();

	// Set values of label, value, default_value and helptext (common for all types)
	void initialise(std::string _label, std::string _default_value, std::string _helptext);

	// Contains $$ for SchedulerVariable
	bool isSchedulerVariable();

	// Get a string value
	std::string getString();

	// Set a string value
	void setString(std::string set_to);

	// Get a string value
	Node getNode();

	// Get a numbered value
	float getNumber(std::string &errmsg);

	// Get a boolean value
	bool getBoolean();

	// Read value from an ifstream. Return false if cannot find it
	bool readValue(std::ifstream& in);

	// Write value to an ostream
	void writeValue(std::ostream& out);
};

class RelionJob
{

public:

	// The name of this job
	std::string outputName;

	// The alias to this job
	std::string alias;

	// Name of the hidden file
	std::string hidden_name;

	// Which job type is this?
	int type;

	// Is this a continuation job?
	bool is_continue;

	// List of Nodes of input to this process
	std::vector<Node> inputNodes;

	// List of Nodes of output from this process
	std::vector<Node> outputNodes;

	// All the options to this job
	std::map<std::string, JobOption > joboptions;


public:
	// Constructor
	RelionJob() { clear(); };

	// Empty Destructor
	~RelionJob() { clear(); };

	// Clear everything
	void clear()
	{
		outputName = alias = "";
		type = -1;
		inputNodes.clear();
		outputNodes.clear();
		joboptions.clear();
		is_continue = false;
	}

	// Returns true if the option is present in joboptions
	bool containsLabel(std::string label, std::string &option);

	// Set this option in the job
	void setOption(std::string setOptionLine);

	// write/read settings to disc
	// fn is a directory name (e.g. Refine3D/job123/) or a STAR file
	bool read(std::string fn, bool &_is_continue, bool do_initialise = false); // return false if unsuccessful
	void write(std::string fn);

	// Write the job submission script
	bool saveJobSubmissionScript(std::string newfilename, std::string outputname, std::vector<std::string> commands, std::string &error_message);

	// Initialise pipeline stuff for each job, return outputname
	void initialisePipeline(std::string &outputname, std::string defaultname, int job_counter);

	// Prepare the final (job submission or combined (mpi) command of possibly multiple lines)
	// Returns true to go ahead, and false to cancel
	bool prepareFinalCommand(std::string &outputname, std::vector<std::string> &commands, std::string &final_command,
			bool do_makedir, std::string &warning_message);

	// Initialise the generic RelionJob
	void initialise(int job_type);

	// Generic getCommands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands,
	 		std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	// Now all the specific job types are defined
	void initialiseImportJob();
	bool getCommandsImportJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseMotioncorrJob();
	bool getCommandsMotioncorrJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseCtffindJob();
	bool getCommandsCtffindJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseManualpickJob();
	bool getCommandsManualpickJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseAutopickJob();
	bool getCommandsAutopickJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseExtractJob();
	bool getCommandsExtractJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseSelectJob();
	bool getCommandsSelectJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseClass2DJob();
	bool getCommandsClass2DJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseInimodelJob();
	bool getCommandsInimodelJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseClass3DJob();
	bool getCommandsClass3DJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseAutorefineJob();
	bool getCommandsAutorefineJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseMultiBodyJob();
	bool getCommandsMultiBodyJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseMaskcreateJob();
	bool getCommandsMaskcreateJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseJoinstarJob();
	bool getCommandsJoinstarJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseSubtractJob();
	bool getCommandsSubtractJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialisePostprocessJob();
	bool getCommandsPostprocessJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseLocalresJob();
	bool getCommandsLocalresJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseMotionrefineJob();
	bool getCommandsMotionrefineJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseCtfrefineJob();
	bool getCommandsCtfrefineJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseExternalJob();
	bool getCommandsExternalJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);
};

#endif /* SRC_PIPELINE_JOBS_H_ */
