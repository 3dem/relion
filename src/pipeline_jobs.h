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

#define RADIO_SAMPLING 0
#define RADIO_NODETYPE 1
#define RADIO_GAIN_ROTATION 2
#define RADIO_GAIN_FLIP 3

// Our own defaults at LMB are the hard-coded ones
#define DEFAULTQSUBLOCATION "/public/EM/RELION/relion/bin/relion_qsub.csh"
#define DEFAULTCTFFINDLOCATION "/public/EM/ctffind/ctffind.exe"
#define DEFAULTMOTIONCOR2LOCATION "/public/EM/MOTIONCOR2/MotionCor2"
#define DEFAULTGCTFLOCATION "/public/EM/Gctf/bin/Gctf"
#define DEFAULTTOPAZLOCATION "/public/EM/RELION/topaz"
#define DEFAULTRESMAPLOCATION "/public/EM/ResMap/ResMap-1.1.4-linux64"
#define DEFAULTPYTHONLOCATION "python"
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
// Modify the loop in JobOption::getHealPixOrder when you add more choices!

static const std::vector<std::string> job_nodetype_options{
	"Particle coordinates (*.box, *_pick.star)",
	"Particles STAR file (.star)",
	"Multiple (2D or 3D) references (.star or .mrcs)",
	"Micrographs STAR file (.star)",
	"3D reference (.mrc)",
	"3D mask (.mrc)",
	"Unfiltered half-map (unfil.mrc)"
};

static const std::vector<std::string> job_nodetype_options_tomo{
	"Set of tomograms STAR file (.star)",
	"Particles STAR file (.star)",
	"Multiple (2D or 3D) references (.star or .mrcs)",
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

static const std::vector<std::string> job_tomo_align_shiftonly_options{
		"Entire micrographs",
		"Only particles"
};

static const std::vector<std::string> job_tomo_align_def_model{
"linear",
"spline",
"Fourier"
};

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
#define NODE_REFS       	5 // A STAR file with one or multiple references, e.g. autopick_references.star
#define NODE_3DREF       	6 // A single 3D-reference, e.g. map.mrc
#define NODE_MASK			7 // 3D mask, e.g. mask.mrc or masks.star
// SHWS 28nov2019: NODE_MODEL should disappear now, only here for backwards compatibility
#define NODE_MODEL		    8 // A model STAR-file for class selection
#define NODE_OPTIMISER		9 // An optimiser STAR-file for job continuation
#define NODE_HALFMAP		10// Unfiltered half-maps from 3D auto-refine, e.g. run1_half?_class001_unfil.mrc
// SHWS 28nov2019: NODE_FINALMAP should disappear now, only here for backwards compatibility
#define NODE_FINALMAP		11// Sharpened final map from post-processing (cannot be used as input)
#define NODE_RESMAP			12// Resmap with local resolution (cannot be used as input)
#define NODE_PDF_LOGFILE    13// PDF logfile
#define NODE_POST           14// Postprocess STAR file (with FSC curve, unfil half-maps, masks etc in it: used by Jasenko's programs
#define NODE_POLISH_PARAMS  15// Txt file with optimal parameters for Bayesian polishing

#define NODE_TOMO_OPTIMISATION 50 // Jasenko's combined optimisation set for subtomogram averaging
#define NODE_TOMO_TOMOGRAMS    51 // Jasenko's set of tomograms
#define NODE_TOMO_TRAJECTORIES 52 // Jasenko's definition of subtomogram motion trajectories
#define NODE_TOMO_MANIFOLDS    53 // Jasenko's definition of 3D shapes (for picking of subtomograms)

// These labels were temporarily in use during alpha-testing of relion-4.0. To be removed once relion-4 is stable
// Now replaced with the CCPEM-pipeliner compatble node names

#define NODE_MOVIES_LABEL   	     "relion.MovieStar"
#define NODE_MICS_LABEL			     "relion.MicrographStar"
#define NODE_MIC_COORDS_LABEL	     "relion.CoordinateStar"
#define NODE_PART_DATA_LABEL	     "relion.ParticleStar"
#define NODE_REFS_LABEL              "relion.ReferenceStar"
#define NODE_3DREF_LABEL       	     "relion.DensityMap"
#define NODE_MASK_LABEL			     "relion.Mask"
#define NODE_OPTIMISER_LABEL	     "relion.OptimiserStar"
#define NODE_HALFMAP_LABEL		     "relion.HalfMap"
#define NODE_RESMAP_LABEL		     "relion.LocalResolutionMap"
#define NODE_PDF_LOGFILE_LABEL       "relion.PdfLogfile"
#define NODE_POST_LABEL              "relion.PostprocessStar"
#define NODE_POLISH_PARAMS_LABEL     "relion.PolishParams"
#define NODE_TOMO_OPTIMISATION_LABEL "relion.TomoOptimisationSet"
#define NODE_TOMO_TOMOGRAMS_LABEL    "relion.TomoTomogramSet"
#define NODE_TOMO_TRAJECTORIES_LABEL "relion.TomoTrajectorySet"
#define NODE_TOMO_MANIFOLDS_LABEL    "relion.TomoManifoldSet"


// nodes compatible with the CCPEM pipeliner
// General types that are used as input nodes in the pipeliner
#define NODE_MOVIES_CPIPE				1
#define NODE_MICS_CPIPE					2
#define NODE_2DIMGS_CPIPE				3
#define NODE_MAP_CPIPE					4
#define NODE_PARTS_CPIPE				5
#define NODE_COORDS_CPIPE				6
#define NODE_COORDS_HELIX_CPIPE			7
#define NODE_PARTS_HELIX_CPIPE			8
#define NODE_OPTIMISER_CPIPE			9
#define NODE_MASK_CPIPE					10
#define NODE_HALFMAP_CPIPE				11
#define NODE_RESMAP_CPIPE				12
#define NODE_LOGFILE_CPIPE				13

// Job-specific output nodes
// Import
#define OUTNODE_IMPORT_MOVIES				21
#define OUTNODE_IMPORT_MICS 				22
#define OUTNODE_IMPORT_COORDS				23
#define OUTNODE_IMPORT_PARTS				24
#define OUTNODE_IMPORT_2DIMG				25
#define OUTNODE_IMPORT_MAP					26
#define OUTNODE_IMPORT_MASK					27
#define OUTNODE_IMPORT_HALFMAP				28

// MotionCorr
#define OUTNODE_MOCORR_MICS					31
#define OUTNODE_MOCORR_LOG					32

// CtfFind
#define OUTNODE_CTFFIND_MICS				41
#define OUTNODE_CTFFIND_LOG 				42

// ManualPick
#define OUTNODE_MANPICK_MICS				51
#define OUTNODE_MANPICK_COORDS				52
#define OUTNODE_MANPICK_COORDS_HELIX		53

// AutoPick
#define OUTNODE_AUTOPICK_COORDS				61
#define OUTNODE_AUTOPICK_LOG				62
#define OUTNODE_AUTOPICK_TOPAZMODEL			63
#define OUTNODE_AUTOPICK_MICS				64

// Extract
#define OUTNODE_EXTRACT_PARTS				71
#define OUTNODE_EXTRACT_PARTS_HELIX			72
#define OUTNODE_EXTRACT_COORDS_HELIX		73
#define OUTNODE_EXTRACT_PARTS_REEX			74
#define OUTNODE_EXTRACT_COORDS_REEX			75

// Class2D

#define OUTNODE_CLASS2D_PARTS				81
#define OUTNODE_CLASS2D_OPT					82
#define OUTNODE_CLASS2D_PARTS_HELIX			83

// Select
#define OUTNODE_SELECT_MICS				 	91
#define OUTNODE_SELECT_PARTS				92
#define OUTNODE_SELECT_OPT					93
#define OUTNODE_SELECT_CLAVS				94

// Initial model
#define OUTNODE_INIMOD_MAP					101

// Class3D
#define OUTNODE_CLASS3D_OPT					111
#define OUTNODE_CLASS3D_MAP					112
#define OUTNODE_CLASS3D_PARTS				113
#define OUTNODE_CLASS3D_PARTS_HELIX			114

// Refine3D
#define OUTNODE_REFINE3D_HALFMAP			121
#define OUTNODE_REFINE3D_OPT				122
#define OUTNODE_REFINE3D_MAP				123
#define OUTNODE_REFINE3D_PARTS				124
#define OUTNODE_REFINE3D_PARTS_HELIX		125

// MultiBody
#define OUTNODE_MULTIBODY_HALFMAP			131
#define OUTNODE_MULTIBODY_PARTS				132
#define OUTNODE_MULTIBODY_OPT				133
#define OUTNODE_MULTIBODY_FLEXLOG			134
#define OUTNODE_MULTIBODY_SEL_PARTS			135

// MaskCreate
#define OUTNODE_MASK3D_MASK					141

//JoinStar
// input nodetypes are used as output node types
// so none needed

// Subtract
#define OUTNODE_SUBTRACT_SUBTRACTED			151
#define OUTNODE_SUBTRACT_REVERTED			152

//Local Ref
#define OUTNODE_LOCRES_OWN				161
#define OUTNODE_LOCRES_RESMAP			162
#define OUTNODE_LOCRES_FILTMAP			163
#define OUTNODE_LOCRES_LOG				164

// CtfRefine
#define OUTNODE_CTFREFINE_REFINEPARTS	171
#define OUTNODE_CTFREFINE_LOG			172
#define OUTNODE_CTFREFINE_ANISOPARTS	173

// Polish
#define OUTNODE_POLISH_PARTS			181
#define OUTNODE_POLISH_LOG				182
#define OUTNODE_POLISH_PARAMS			183

// PostProcess
#define OUTNODE_POST					191
#define OUTNODE_POST_MAP				192
#define OUTNODE_POST_MASKED				193
#define OUTNODE_POST_LOG				194

// Tomo
#define OUTNODE_TOMO_OPTIMISATION		201
#define OUTNODE_TOMO_TOMOGRAMS   		202
#define OUTNODE_TOMO_TRAJECTORIES 		203
#define OUTNODE_TOMO_MANIFOLDS    		204
#define OUTNODE_TOMO_PARTS				205
#define OUTNODE_TOMO_MAP				206
#define OUTNODE_TOMO_HALFMAP			207
#define OUTNODE_TOMO_POST				208
#define OUTNODE_TOMO_POST_LOG			209
#define OUTNODE_TOMO_FRAMEALIGN_LOG	    210
#define OUTNODE_TOMO_CTFREFINE_LOG		211


#define LABEL_MOVIES_CPIPE             "MicrographMoviesData.star.relion"
#define LABEL_MICS_CPIPE               "MicrographsData.star.relion"
#define LABEL_2DIMGS_CPIPE             "ImagesData.star.relion"
#define LABEL_MAP_CPIPE                "DensityMap.mrc"
#define LABEL_PARTS_CPIPE              "ParticlesData.star.relion"
#define LABEL_COORDS_CPIPE             "MicrographsCoords.star.relion"
#define LABEL_COORDS_HELIX_CPIPE       "MicrographsCoords.star.relion.helixstartend"
#define LABEL_PARTS_HELIX_CPIPE        "ParticlesData.star.relion.helicalsegments"
#define LABEL_OPTIMISER_CPIPE          "ProcessData.star.relion.optimiser"
#define LABEL_MASK_CPIPE               "Mask3D.mrc"
#define LABEL_HALFMAP_CPIPE	           "DensityMap.mrc.halfmap"
#define LABEL_RESMAP_CPIPE             "Image3D.mrc.localresmap"
#define LABEL_LOGFILE_CPIPE            "LogFile.pdf.relion"
#define LABEL_IMPORT_MOVIES            "MicrographMoviesData.star.relion"
#define LABEL_IMPORT_MICS              "MicrographsData.star.relion"
#define LABEL_IMPORT_COORDS            "MicrographsCoords.star.relion"
#define LABEL_IMPORT_PARTS             "ParticlesData.star.relion"
#define LABEL_IMPORT_2DIMG             "ImagesData.star.relion"
#define LABEL_IMPORT_MAP               "DensityMap.mrc"
#define LABEL_IMPORT_MASK              "Mask3D.mrc"
#define LABEL_IMPORT_HALFMAP           "DensityMap.mrc.halfmap"
#define LABEL_MOCORR_MICS              "MicrographsData.star.relion.motioncorr"
#define LABEL_MOCORR_LOG               "LogFile.pdf.relion.motioncorr"
#define LABEL_CTFFIND_MICS             "MicrographsData.star.relion.ctf"
#define LABEL_CTFFIND_LOG              "LogFile.pdf.relion.ctffind"
#define LABEL_MANPICK_MICS             "MicrographsData.star.relion"
#define LABEL_MANPICK_COORDS           "MicrographsCoords.star.relion.manualpick"
#define LABEL_MANPICK_COORDS_HELIX     "MicrographsCoords.star.relion.manualpick.helixstartend"
#define LABEL_AUTOPICK_COORDS          "MicrographsCoords.star.relion.autopick"
#define LABEL_AUTOPICK_LOG             "LogFile.pdf.relion.autopick"
#define LABEL_AUTOPICK_TOPAZMODEL      "ProcessData.sav.topaz.model"	// to be added?
#define LABEL_AUTOPICK_MICS            "MicrographsData.star.relion"
#define LABEL_EXTRACT_PARTS            "ParticlesData.star.relion"
#define LABEL_EXTRACT_PARTS_HELIX      "ParticlesData.star.relion.helicalsegments"
#define LABEL_EXTRACT_COORDS_HELIX     "MicrographsCoords.star.relion.helixstartend"
#define LABEL_EXTRACT_PARTS_REEX       "ParticlesData.star.relion.reextract"
#define LABEL_EXTRACT_COORDS_REEX      "MicrographsCoords.star.relion.reextract"
#define LABEL_CLASS2D_PARTS            "ParticlesData.star.relion.class2d"
#define LABEL_CLASS2D_OPT              "ProcessData.star.relion.optimiser.class2d"
#define LABEL_CLASS2D_PARTS_HELIX      "ParticlesData.star.relion.class2d.helicalsegments"
#define LABEL_SELECT_MICS              "MicrographsData.star.relion"
#define LABEL_SELECT_MOVS              "MicrographMoviesData.star.relion"
#define LABEL_SELECT_PARTS             "ParticlesData.star.relion"
#define LABEL_SELECT_OPT               "ProcessData.star.relion.optimiser.autoselect"
#define LABEL_SELECT_CLAVS             "ImagesData.star.relion.classaverages"
#define LABEL_INIMOD_MAP               "DensityMap.mrc.relion.initialmodel"
#define LABEL_CLASS3D_OPT              "ProcessData.star.relion.optimiser.class3d"
#define LABEL_CLASS3D_MAP              "DensityMap.mrc.relion.class3d"
#define LABEL_CLASS3D_PARTS            "ParticlesData.star.relion.class3d"
#define LABEL_CLASS3D_PARTS_HELIX      "ParticlesData.star.relion.class3d.helicalsegments"
#define LABEL_REFINE3D_HALFMAP         "DensityMap.mrc.relion.halfmap.refine3d"
#define LABEL_REFINE3D_OPT             "ProcessData.star.relion.optimiser.refine3d"
#define LABEL_REFINE3D_MAP             "DensityMap.mrc.relion.refine3d"
#define LABEL_REFINE3D_PARTS           "ParticlesData.star.relion.refine3d"
#define LABEL_REFINE3D_PARTS_HELIX     "ParticlesData.star.relion.refine3d.helicalsegements"
#define LABEL_MULTIBODY_HALFMAP        "DensityMap.mrc.relion.halfmap.multibody"
#define LABEL_MULTIBODY_PARTS          "ParticlesData.star.relion.multibody"
#define LABEL_MULTIBODY_OPT            "ProcessData.star.relion.optimiser.multibody"
#define LABEL_MULTIBODY_FLEXLOG        "LogFile.pdf.relion.flexanalysis"
#define LABEL_MULTIBODY_SEL_PARTS      "ParticlesData.star.relion.flexanalysis.eigenselected"
#define LABEL_MASK3D_MASK              "Mask3D.mrc.relion"
#define LABEL_SUBTRACT_SUBTRACTED      "ParticlesData.star.relion.subtracted"
#define LABEL_SUBTRACT_REVERTED        "ParticlesData.star.relion"
#define LABEL_LOCRES_OWN               "Image3D.mrc.relion.localresmap"
#define LABEL_LOCRES_RESMAP            "Image3D.mrc.resmap.localresmap"
#define LABEL_LOCRES_FILTMAP           "DensityMap.mrc.relion.localresfiltered"
#define LABEL_LOCRES_LOG               "LogFile.pdf.relion.localres"
#define LABEL_CTFREFINE_REFINEPARTS    "ParticlesData.star.relion.ctfrefine"
#define LABEL_CTFREFINE_LOG            "LogFile.pdf.relion.ctfrefine"
#define LABEL_CTFREFINE_ANISOPARTS     "ParticlesData.star.relion.anisomagrefine"
#define LABEL_POLISH_PARTS             "ParticlesData.star.relion.polished"
#define LABEL_POLISH_LOG               "LogFile.pdf.relion.polish"
#define LABEL_POLISH_PARAMS            "ProcessData.txt.relion.polish.params"
#define LABEL_POST                     "ProcessData.star.relion.postprocess"
#define LABEL_POST_MAP                 "DensityMap.mrc.relion.postprocess"
#define LABEL_POST_MASKED              "DensityMap.mrc.relion.postprocess.masked"
#define LABEL_POST_LOG                 "LogFile.pdf.relion.postprocess"
#define LABEL_TOMO_OPTIMISATION        "ProcessData.star.relion.tomo.optimisation_set"
#define LABEL_TOMO_TOMOGRAMS           "ProcessData.star.relion.tomo.relion.tomogram_set"
#define LABEL_TOMO_TRAJECTORIES        "ProcessData.star.relion.tomo.relion.trajectory_set"
#define LABEL_TOMO_MANIFOLDS           "ProcessData.star.relion.tomo.manifoldset"
#define LABEL_TOMO_PARTS               "Particles.star.relion.tomo"
#define LABEL_TOMO_MAP                 "DensityMap.mrc.relion.tomo.subvolume"
#define LABEL_TOMO_HALFMAP             "DensityMap.mrc.relion.tomo.halfmap"
#define LABEL_TOMO_POST                "ProcessData.star.relion.tomo.postprocess"
#define LABEL_TOMO_POST_LOG            "LogFile.pdf.relion.tomo.postprocess"
#define LABEL_TOMO_FRAMEALIGN_LOG      "LogFile.pdf.relion.tomo.framealign"
#define LABEL_TOMO_CTFREFINE_LOG       "LogFile.pdf.relion.tomo.ctfrefine"


static std::map<int, std::string> node_type2pipeliner_label = {{NODE_MOVIES_CPIPE, LABEL_MOVIES_CPIPE},
	{NODE_MICS_CPIPE, LABEL_MICS_CPIPE},
	{NODE_2DIMGS_CPIPE, LABEL_2DIMGS_CPIPE},
	{NODE_MAP_CPIPE, LABEL_MAP_CPIPE},
	{NODE_PARTS_CPIPE, LABEL_PARTS_CPIPE},
	{NODE_COORDS_CPIPE, LABEL_COORDS_CPIPE},
	{NODE_COORDS_HELIX_CPIPE, LABEL_COORDS_HELIX_CPIPE},
	{NODE_PARTS_HELIX_CPIPE, LABEL_PARTS_HELIX_CPIPE},
	{NODE_OPTIMISER_CPIPE, LABEL_OPTIMISER_CPIPE},
	{NODE_MASK_CPIPE, LABEL_MASK_CPIPE},
	{NODE_HALFMAP_CPIPE	, LABEL_HALFMAP_CPIPE},
	{NODE_RESMAP_CPIPE, LABEL_RESMAP_CPIPE},
	{NODE_LOGFILE_CPIPE, LABEL_LOGFILE_CPIPE},
	{OUTNODE_IMPORT_MOVIES, LABEL_IMPORT_MOVIES},
	{OUTNODE_IMPORT_MICS, LABEL_IMPORT_MICS},
	{OUTNODE_IMPORT_COORDS, LABEL_IMPORT_COORDS},
	{OUTNODE_IMPORT_PARTS, LABEL_IMPORT_PARTS},
	{OUTNODE_IMPORT_2DIMG, LABEL_IMPORT_2DIMG},
	{OUTNODE_IMPORT_MAP, LABEL_IMPORT_MAP},
	{OUTNODE_IMPORT_MASK, LABEL_IMPORT_MASK},
	{OUTNODE_IMPORT_HALFMAP, LABEL_IMPORT_HALFMAP},
	{OUTNODE_MOCORR_MICS, LABEL_MOCORR_MICS},
	{OUTNODE_MOCORR_LOG, LABEL_MOCORR_LOG},
	{OUTNODE_CTFFIND_MICS, LABEL_CTFFIND_MICS},
	{OUTNODE_CTFFIND_LOG, LABEL_CTFFIND_LOG},
	{OUTNODE_MANPICK_MICS, LABEL_MANPICK_MICS},
	{OUTNODE_MANPICK_COORDS, LABEL_MANPICK_COORDS},
	{OUTNODE_MANPICK_COORDS_HELIX, LABEL_MANPICK_COORDS_HELIX},
	{OUTNODE_AUTOPICK_COORDS, LABEL_AUTOPICK_COORDS},
	{OUTNODE_AUTOPICK_LOG, LABEL_AUTOPICK_LOG},
	{OUTNODE_AUTOPICK_TOPAZMODEL, LABEL_AUTOPICK_TOPAZMODEL},	// to be added?
	{OUTNODE_AUTOPICK_MICS, LABEL_AUTOPICK_MICS},
	{OUTNODE_EXTRACT_PARTS, LABEL_EXTRACT_PARTS},
	{OUTNODE_EXTRACT_PARTS_HELIX, LABEL_EXTRACT_PARTS_HELIX},
	{OUTNODE_EXTRACT_COORDS_HELIX, LABEL_EXTRACT_COORDS_HELIX},
	{OUTNODE_EXTRACT_PARTS_REEX, LABEL_EXTRACT_PARTS_REEX},
	{OUTNODE_EXTRACT_COORDS_REEX, LABEL_EXTRACT_COORDS_REEX},
	{OUTNODE_CLASS2D_PARTS, LABEL_CLASS2D_PARTS},
	{OUTNODE_CLASS2D_OPT, LABEL_CLASS2D_OPT},
	{OUTNODE_CLASS2D_PARTS_HELIX, LABEL_CLASS2D_PARTS_HELIX},
	{OUTNODE_SELECT_MICS, LABEL_SELECT_MICS},
	{OUTNODE_SELECT_PARTS, LABEL_SELECT_PARTS},
	{OUTNODE_SELECT_OPT, LABEL_SELECT_OPT},
	{OUTNODE_SELECT_CLAVS, LABEL_SELECT_CLAVS},
	{OUTNODE_INIMOD_MAP, LABEL_INIMOD_MAP},
	{OUTNODE_CLASS3D_OPT, LABEL_CLASS3D_OPT},
	{OUTNODE_CLASS3D_MAP, LABEL_CLASS3D_MAP},
	{OUTNODE_CLASS3D_PARTS, LABEL_CLASS3D_PARTS},
	{OUTNODE_CLASS3D_PARTS_HELIX, LABEL_CLASS3D_PARTS_HELIX},
	{OUTNODE_REFINE3D_HALFMAP, LABEL_REFINE3D_HALFMAP},
	{OUTNODE_REFINE3D_OPT, LABEL_REFINE3D_OPT},
	{OUTNODE_REFINE3D_MAP, LABEL_REFINE3D_MAP},
	{OUTNODE_REFINE3D_PARTS, LABEL_REFINE3D_PARTS},
	{OUTNODE_REFINE3D_PARTS_HELIX, LABEL_REFINE3D_PARTS_HELIX},
	{OUTNODE_MULTIBODY_HALFMAP, LABEL_MULTIBODY_HALFMAP},
	{OUTNODE_MULTIBODY_PARTS, LABEL_MULTIBODY_PARTS},
	{OUTNODE_MULTIBODY_OPT, LABEL_MULTIBODY_OPT},
	{OUTNODE_MULTIBODY_FLEXLOG, LABEL_MULTIBODY_FLEXLOG},
	{OUTNODE_MULTIBODY_SEL_PARTS, LABEL_MULTIBODY_SEL_PARTS},
	{OUTNODE_MASK3D_MASK, LABEL_MASK3D_MASK},
	{OUTNODE_SUBTRACT_SUBTRACTED, LABEL_SUBTRACT_SUBTRACTED},
	{OUTNODE_SUBTRACT_REVERTED, LABEL_SUBTRACT_REVERTED},
	{OUTNODE_LOCRES_OWN, LABEL_LOCRES_OWN},
	{OUTNODE_LOCRES_RESMAP, LABEL_LOCRES_RESMAP},
	{OUTNODE_LOCRES_FILTMAP, LABEL_LOCRES_FILTMAP},
	{OUTNODE_LOCRES_LOG, LABEL_LOCRES_LOG},
	{OUTNODE_CTFREFINE_REFINEPARTS, LABEL_CTFREFINE_REFINEPARTS},
	{OUTNODE_CTFREFINE_LOG, LABEL_CTFREFINE_LOG},
	{OUTNODE_CTFREFINE_ANISOPARTS, LABEL_CTFREFINE_ANISOPARTS},
	{OUTNODE_POLISH_PARTS, LABEL_POLISH_PARTS},
	{OUTNODE_POLISH_LOG, LABEL_POLISH_LOG},
	{OUTNODE_POLISH_PARAMS, LABEL_POLISH_PARAMS},
	{OUTNODE_POST, LABEL_POST},
	{OUTNODE_POST_MAP, LABEL_POST_MAP},
	{OUTNODE_POST_MASKED, LABEL_POST_MASKED},
	{OUTNODE_POST_LOG, LABEL_POST_LOG},
	{OUTNODE_TOMO_OPTIMISATION, LABEL_TOMO_OPTIMISATION},
	{OUTNODE_TOMO_TOMOGRAMS, LABEL_TOMO_TOMOGRAMS},
	{OUTNODE_TOMO_TRAJECTORIES, LABEL_TOMO_TRAJECTORIES},
	{OUTNODE_TOMO_MANIFOLDS, LABEL_TOMO_MANIFOLDS},
	{OUTNODE_TOMO_PARTS, LABEL_TOMO_PARTS},
	{OUTNODE_TOMO_MAP, LABEL_TOMO_MAP},
	{OUTNODE_TOMO_HALFMAP, LABEL_TOMO_HALFMAP},
	{OUTNODE_TOMO_POST, LABEL_TOMO_POST},
	{OUTNODE_TOMO_POST_LOG, LABEL_TOMO_POST_LOG},
    {OUTNODE_TOMO_FRAMEALIGN_LOG, LABEL_TOMO_FRAMEALIGN_LOG},
    {OUTNODE_TOMO_CTFREFINE_LOG, LABEL_TOMO_CTFREFINE_LOG}};


static std::map<std::string, int> pipeliner_label2type = {{LABEL_MOVIES_CPIPE, NODE_MOVIES_CPIPE},
	{LABEL_MICS_CPIPE, NODE_MICS_CPIPE},
	{LABEL_2DIMGS_CPIPE, NODE_2DIMGS_CPIPE},
	{LABEL_MAP_CPIPE, NODE_MAP_CPIPE},
	{LABEL_PARTS_CPIPE, NODE_PARTS_CPIPE},
	{LABEL_COORDS_CPIPE, NODE_COORDS_CPIPE},
	{LABEL_COORDS_HELIX_CPIPE, NODE_COORDS_HELIX_CPIPE},
	{LABEL_PARTS_HELIX_CPIPE, NODE_PARTS_HELIX_CPIPE},
	{LABEL_OPTIMISER_CPIPE, NODE_OPTIMISER_CPIPE},
	{LABEL_MASK_CPIPE, NODE_MASK_CPIPE},
	{LABEL_HALFMAP_CPIPE	, NODE_HALFMAP_CPIPE},
	{LABEL_RESMAP_CPIPE, NODE_RESMAP_CPIPE},
	{LABEL_LOGFILE_CPIPE, NODE_LOGFILE_CPIPE},
	{LABEL_IMPORT_MOVIES, OUTNODE_IMPORT_MOVIES},
	{LABEL_IMPORT_MICS, OUTNODE_IMPORT_MICS},
	{LABEL_IMPORT_COORDS, OUTNODE_IMPORT_COORDS},
	{LABEL_IMPORT_PARTS, OUTNODE_IMPORT_PARTS},
	{LABEL_IMPORT_2DIMG, OUTNODE_IMPORT_2DIMG},
	{LABEL_IMPORT_MAP, OUTNODE_IMPORT_MAP},
	{LABEL_IMPORT_MASK, OUTNODE_IMPORT_MASK},
	{LABEL_IMPORT_HALFMAP, OUTNODE_IMPORT_HALFMAP},
	{LABEL_MOCORR_MICS, OUTNODE_MOCORR_MICS},
	{LABEL_MOCORR_LOG, OUTNODE_MOCORR_LOG},
	{LABEL_CTFFIND_MICS, OUTNODE_CTFFIND_MICS},
	{LABEL_CTFFIND_LOG, OUTNODE_CTFFIND_LOG},
	{LABEL_MANPICK_MICS, OUTNODE_MANPICK_MICS},
	{LABEL_MANPICK_COORDS, OUTNODE_MANPICK_COORDS},
	{LABEL_MANPICK_COORDS_HELIX, OUTNODE_MANPICK_COORDS_HELIX},
	{LABEL_AUTOPICK_COORDS, OUTNODE_AUTOPICK_COORDS},
	{LABEL_AUTOPICK_LOG, OUTNODE_AUTOPICK_LOG},
	{LABEL_AUTOPICK_TOPAZMODEL, OUTNODE_AUTOPICK_TOPAZMODEL},	// to be added?
	{LABEL_AUTOPICK_MICS, OUTNODE_AUTOPICK_MICS},
	{LABEL_EXTRACT_PARTS, OUTNODE_EXTRACT_PARTS},
	{LABEL_EXTRACT_PARTS_HELIX, OUTNODE_EXTRACT_PARTS_HELIX},
	{LABEL_EXTRACT_COORDS_HELIX, OUTNODE_EXTRACT_COORDS_HELIX},
	{LABEL_EXTRACT_PARTS_REEX, OUTNODE_EXTRACT_PARTS_REEX},
	{LABEL_EXTRACT_COORDS_REEX, OUTNODE_EXTRACT_COORDS_REEX},
	{LABEL_CLASS2D_PARTS, OUTNODE_CLASS2D_PARTS},
	{LABEL_CLASS2D_OPT, OUTNODE_CLASS2D_OPT},
	{LABEL_CLASS2D_PARTS_HELIX, OUTNODE_CLASS2D_PARTS_HELIX},
	{LABEL_SELECT_MICS, OUTNODE_SELECT_MICS},
	{LABEL_SELECT_PARTS, OUTNODE_SELECT_PARTS},
	{LABEL_SELECT_OPT, OUTNODE_SELECT_OPT},
	{LABEL_SELECT_CLAVS, OUTNODE_SELECT_CLAVS},
	{LABEL_INIMOD_MAP, OUTNODE_INIMOD_MAP},
	{LABEL_CLASS3D_OPT, OUTNODE_CLASS3D_OPT},
	{LABEL_CLASS3D_MAP, OUTNODE_CLASS3D_MAP},
	{LABEL_CLASS3D_PARTS, OUTNODE_CLASS3D_PARTS},
	{LABEL_CLASS3D_PARTS_HELIX, OUTNODE_CLASS3D_PARTS_HELIX},
	{LABEL_REFINE3D_HALFMAP, OUTNODE_REFINE3D_HALFMAP},
	{LABEL_REFINE3D_OPT, OUTNODE_REFINE3D_OPT},
	{LABEL_REFINE3D_MAP, OUTNODE_REFINE3D_MAP},
	{LABEL_REFINE3D_PARTS, OUTNODE_REFINE3D_PARTS},
	{LABEL_REFINE3D_PARTS_HELIX, OUTNODE_REFINE3D_PARTS_HELIX},
	{LABEL_MULTIBODY_HALFMAP, OUTNODE_MULTIBODY_HALFMAP},
	{LABEL_MULTIBODY_PARTS, OUTNODE_MULTIBODY_PARTS},
	{LABEL_MULTIBODY_OPT, OUTNODE_MULTIBODY_OPT},
	{LABEL_MULTIBODY_FLEXLOG, OUTNODE_MULTIBODY_FLEXLOG},
	{LABEL_MULTIBODY_SEL_PARTS, OUTNODE_MULTIBODY_SEL_PARTS},
	{LABEL_MASK3D_MASK, OUTNODE_MASK3D_MASK},
	{LABEL_SUBTRACT_SUBTRACTED, OUTNODE_SUBTRACT_SUBTRACTED},
	{LABEL_SUBTRACT_REVERTED, OUTNODE_SUBTRACT_REVERTED},
	{LABEL_LOCRES_OWN, OUTNODE_LOCRES_OWN},
	{LABEL_LOCRES_RESMAP, OUTNODE_LOCRES_RESMAP},
	{LABEL_LOCRES_FILTMAP, OUTNODE_LOCRES_FILTMAP},
	{LABEL_LOCRES_LOG, OUTNODE_LOCRES_LOG},
	{LABEL_CTFREFINE_REFINEPARTS, OUTNODE_CTFREFINE_REFINEPARTS},
	{LABEL_CTFREFINE_LOG, OUTNODE_CTFREFINE_LOG},
	{LABEL_CTFREFINE_ANISOPARTS, OUTNODE_CTFREFINE_ANISOPARTS},
	{LABEL_POLISH_PARTS, OUTNODE_POLISH_PARTS},
	{LABEL_POLISH_LOG, OUTNODE_POLISH_LOG},
	{LABEL_POLISH_PARAMS, OUTNODE_POLISH_PARAMS},
	{LABEL_POST, OUTNODE_POST},
	{LABEL_POST_MAP, OUTNODE_POST_MAP},
	{LABEL_POST_MASKED, OUTNODE_POST_MASKED},
	{LABEL_POST_LOG, OUTNODE_POST_LOG},
	{LABEL_TOMO_OPTIMISATION, OUTNODE_TOMO_OPTIMISATION},
	{LABEL_TOMO_TOMOGRAMS, OUTNODE_TOMO_TOMOGRAMS},
	{LABEL_TOMO_TRAJECTORIES, OUTNODE_TOMO_TRAJECTORIES},
	{LABEL_TOMO_MANIFOLDS, OUTNODE_TOMO_MANIFOLDS},
	{LABEL_TOMO_PARTS, OUTNODE_TOMO_PARTS},
	{LABEL_TOMO_MAP, OUTNODE_TOMO_MAP},
	{LABEL_TOMO_HALFMAP, OUTNODE_TOMO_HALFMAP},
	{LABEL_TOMO_POST, OUTNODE_TOMO_POST},
	{LABEL_TOMO_POST_LOG, OUTNODE_TOMO_POST_LOG},
	{LABEL_TOMO_FRAMEALIGN_LOG, OUTNODE_TOMO_FRAMEALIGN_LOG},
    {LABEL_TOMO_CTFREFINE_LOG, OUTNODE_TOMO_CTFREFINE_LOG} };

//// Conversion dict for CCPEM-pipeliner compatibility
static std::map<std::string, int> node_label2type = {{NODE_MOVIES_LABEL, NODE_MOVIES},
		{NODE_MICS_LABEL, NODE_MICS_CPIPE},
		{NODE_MIC_COORDS_LABEL, NODE_COORDS_CPIPE},
		{NODE_PART_DATA_LABEL, NODE_PARTS_CPIPE},
		{NODE_REFS_LABEL, NODE_2DIMGS_CPIPE},
		{NODE_3DREF_LABEL, NODE_MAP_CPIPE},
		{NODE_MASK_LABEL, NODE_MASK_CPIPE},
		{NODE_OPTIMISER_LABEL, NODE_OPTIMISER_CPIPE},
		{NODE_HALFMAP_LABEL, NODE_HALFMAP_CPIPE},
		{NODE_RESMAP_LABEL, NODE_RESMAP_CPIPE},
		{NODE_PDF_LOGFILE_LABEL, NODE_LOGFILE_CPIPE},
		{NODE_POST_LABEL, OUTNODE_POST},
		{NODE_POLISH_PARAMS_LABEL, OUTNODE_POLISH_PARAMS},
		{NODE_TOMO_OPTIMISATION_LABEL, OUTNODE_TOMO_OPTIMISATION},
		{NODE_TOMO_TOMOGRAMS_LABEL,    OUTNODE_TOMO_TOMOGRAMS},
		{NODE_TOMO_TRAJECTORIES_LABEL, OUTNODE_TOMO_TRAJECTORIES},
		{NODE_TOMO_MANIFOLDS_LABEL,    OUTNODE_TOMO_MANIFOLDS} };

static std::string get_node_label(int type)
{
	if (node_type2pipeliner_label.find(type) != node_type2pipeliner_label.end())
	{
		return node_type2pipeliner_label.at(type);;
	}
	else
	{
		std::cerr << " WARNING: unrecognised node type: " << type << std::endl;
	}
	return "";
}

static int get_node_type(std::string label)
{

	if (pipeliner_label2type.find(label) != pipeliner_label2type.end())
	{
		return pipeliner_label2type.at(label);;
	}
	else if (node_label2type.find(label) != node_label2type.end())
	{
		return node_label2type.at(label);
	}
	else
	{
		std::cerr << " WARNING: unrecognised node label: " << label << std::endl;
	}
	return -1;

}


// All the directory names of the different types of jobs defined inside the pipeline
#define PROC_IMPORT_DIRNAME           "Import"       // Import any file as a Node of a given type
#define PROC_MOTIONCORR_DIRNAME 	  "MotionCorr"   // Import any file as a Node of a given type
#define PROC_CTFFIND_DIRNAME	      "CtfFind"  	   // Estimate CTF parameters from micrographs for either entire micrographs and/or particles
#define PROC_MANUALPICK_DIRNAME       "ManualPick"   // Manually pick particle coordinates from micrographs
#define PROC_AUTOPICK_DIRNAME		  "AutoPick"     // Automatically pick particle coordinates from micrographs, their CTF and 2D references
#define PROC_EXTRACT_DIRNAME		  "Extract"      // Window particles, normalize, downsize etc from micrographs (also combine CTF into metadata file)
#define PROC_CLASSSELECT_DIRNAME      "Select" 	   // Read in model.star file, and let user interactively select classes through the display (later: auto-selection as well)
#define PROC_2DCLASS_DIRNAME 		  "Class2D"      // 2D classification (from input particles)
#define PROC_3DCLASS_DIRNAME		  "Class3D"      // 3D classification (from input 2D/3D particles, an input 3D-reference, and possibly a 3D mask)
#define PROC_3DAUTO_DIRNAME           "Refine3D"     // 3D auto-refine (from input particles, an input 3Dreference, and possibly a 3D mask)
#define PROC_MASKCREATE_DIRNAME       "MaskCreate"   // Process to create masks from input maps
#define PROC_JOINSTAR_DIRNAME         "JoinStar"     // Process to create masks from input maps
#define PROC_SUBTRACT_DIRNAME         "Subtract"     // Process to subtract projections of parts of the reference from experimental images
#define PROC_POST_DIRNAME			  "PostProcess"  // Post-processing (from unfiltered half-maps and a possibly a 3D mask)
#define PROC_RESMAP_DIRNAME  	      "LocalRes"     // Local resolution estimation (from unfiltered half-maps and a 3D mask)
#define PROC_INIMODEL_DIRNAME		  "InitialModel" // De-novo generation of 3D initial model (using SGD)
#define PROC_MULTIBODY_DIRNAME	      "MultiBody"    // Multi-body refinement
#define PROC_MOTIONREFINE_DIRNAME     "Polish"       // Jasenko's motion fitting program for Bayesian polishing (to replace MovieRefine?)
#define PROC_CTFREFINE_DIRNAME        "CtfRefine"    // Jasenko's program for defocus and beamtilt optimisation
#define PROC_TOMO_IMPORT_DIRNAME      "ImportTomo"              // Import for tomography GUI
#define PROC_TOMO_SUBTOMO_DIRNAME     "PseudoSubtomo"           // Creation of pseudo-subtomograms from tilt series images
#define PROC_TOMO_CTFREFINE_DIRNAME   "CtfRefineTomo"           // CTF refinement (defocus & aberrations) for tomography
#define PROC_TOMO_ALIGN_DIRNAME       "FrameAlignTomo"          // Frame alignment and particle polishing for subtomography
#define PROC_TOMO_RECONSTRUCT_DIRNAME "ReconstructParticleTomo" // Calculation of particle average from the individual tilt series images
#define PROC_EXTERNAL_DIRNAME         "External"     // For running non-relion programs

// All the directory names of the different types of jobs defined inside the pipeline
#define PROC_IMPORT_LABELNEW           "relion.import"       // Import any file as a Node of a given type
#define PROC_MOTIONCORR_LABELNEW 	   "relion.motioncorr"   // Import any file as a Node of a given type
#define PROC_CTFFIND_LABELNEW	       "relion.ctffind"  	   // Estimate CTF parameters from micrographs for either entire micrographs and/or particles
#define PROC_MANUALPICK_LABELNEW       "relion.manualpick"   // Manually pick particle coordinates from micrographs
#define PROC_AUTOPICK_LABELNEW		   "relion.autopick"     // Automatically pick particle coordinates from micrographs, their CTF and 2D references
#define PROC_EXTRACT_LABELNEW	       "relion.extract"      // Window particles, normalize, downsize etc from micrographs (also combine CTF into metadata file)
#define PROC_CLASSSELECT_LABELNEW      "relion.select" 	   // Read in model.star file, and let user interactively select classes through the display (later: auto-selection as well)
#define PROC_2DCLASS_LABELNEW 		   "relion.class2d"      // 2D classification (from input particles)
#define PROC_3DCLASS_LABELNEW		   "relion.class3d"      // 3D classification (from input 2D/3D particles, an input 3D-reference, and possibly a 3D mask)
#define PROC_3DAUTO_LABELNEW           "relion.refine3d"     // 3D auto-refine (from input particles, an input 3Dreference, and possibly a 3D mask)
#define PROC_MASKCREATE_LABELNEW       "relion.maskcreate"   // Process to create masks from input maps
#define PROC_JOINSTAR_LABELNEW         "relion.joinstar"     // Process to create masks from input maps
#define PROC_SUBTRACT_LABELNEW         "relion.subtract"     // Process to subtract projections of parts of the reference from experimental images
#define PROC_POST_LABELNEW			   "relion.postprocess"  // Post-processing (from unfiltered half-maps and a possibly a 3D mask)
#define PROC_RESMAP_LABELNEW  	       "relion.localres"     // Local resolution estimation (from unfiltered half-maps and a 3D mask)
#define PROC_INIMODEL_LABELNEW		   "relion.initialmodel" // De-novo generation of 3D initial model (using SGD)
#define PROC_MULTIBODY_LABELNEW	       "relion.multibody"    // Multi-body refinement
#define PROC_MOTIONREFINE_LABELNEW     "relion.polish"       // Jasenko's motion fitting program for Bayesian polishing (to replace MovieRefine?)
#define PROC_CTFREFINE_LABELNEW        "relion.ctfrefine"    // Jasenko's program for defocus and beamtilt optimisation
#define PROC_TOMO_IMPORT_LABELNEW      "relion.importtomo"              // Import for tomography GUI
#define PROC_TOMO_SUBTOMO_LABELNEW     "relion.pseudosubtomo"           // Creation of pseudo-subtomograms from tilt series images
#define PROC_TOMO_CTFREFINE_LABELNEW   "relion.ctfrefinetomo"           // CTF refinement (defocus & aberrations) for tomography
#define PROC_TOMO_ALIGN_LABELNEW       "relion.framealigntomo"          // Frame alignment and particle polishing for subtomography
#define PROC_TOMO_RECONSTRUCT_LABELNEW "relion.reconstructparticletomo" // Calculation of particle average from the individual tilt series images
#define PROC_EXTERNAL_LABELNEW         "relion.external"     // For running non-relion programs


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
#define PROC_TOMO_IMPORT        50// Import for tomography GUI
#define PROC_TOMO_SUBTOMO   51// Creation of pseudo-subtomograms from tilt series images
#define PROC_TOMO_CTFREFINE     52// CTF refinement (defocus & aberrations for tomography)
#define PROC_TOMO_ALIGN        53// Frame alignment and particle polishing for subtomography
#define PROC_TOMO_RECONSTRUCT       54// Calculation of particle average from the individual tilt series images
#define PROC_EXTERNAL       99// External scripts


static std::map<int, std::string> proc_type2dirname = {{PROC_IMPORT, PROC_IMPORT_DIRNAME},
		{PROC_MOTIONCORR, PROC_MOTIONCORR_DIRNAME},
		{PROC_CTFFIND, PROC_CTFFIND_DIRNAME},
		{PROC_MANUALPICK, PROC_MANUALPICK_DIRNAME},
		{PROC_AUTOPICK, PROC_AUTOPICK_DIRNAME},
		{PROC_EXTRACT, PROC_EXTRACT_DIRNAME},
		{PROC_CLASSSELECT, PROC_CLASSSELECT_DIRNAME},
		{PROC_2DCLASS, PROC_2DCLASS_DIRNAME},
		{PROC_3DCLASS, PROC_3DCLASS_DIRNAME},
		{PROC_3DAUTO, PROC_3DAUTO_DIRNAME},
		{PROC_MASKCREATE,        PROC_MASKCREATE_DIRNAME},
		{PROC_JOINSTAR,        PROC_JOINSTAR_DIRNAME},
		{PROC_SUBTRACT,        PROC_SUBTRACT_DIRNAME},
		{PROC_POST,             PROC_POST_DIRNAME},
		{PROC_RESMAP,           PROC_RESMAP_DIRNAME},
		{PROC_INIMODEL,         PROC_INIMODEL_DIRNAME},
		{PROC_MULTIBODY,        PROC_MULTIBODY_DIRNAME},
		{PROC_MOTIONREFINE,     PROC_MOTIONREFINE_DIRNAME},
		{PROC_CTFREFINE,        PROC_CTFREFINE_DIRNAME},
		{PROC_TOMO_IMPORT,      PROC_TOMO_IMPORT_DIRNAME},
		{PROC_TOMO_SUBTOMO,     PROC_TOMO_SUBTOMO_DIRNAME},
		{PROC_TOMO_CTFREFINE,   PROC_TOMO_CTFREFINE_DIRNAME},
		{PROC_TOMO_ALIGN,       PROC_TOMO_ALIGN_DIRNAME},
		{PROC_TOMO_RECONSTRUCT, PROC_TOMO_RECONSTRUCT_DIRNAME},
		{PROC_EXTERNAL,         PROC_EXTERNAL_DIRNAME}};

static std::map<int, std::string> proc_type2labelnew = {{PROC_IMPORT, PROC_IMPORT_LABELNEW},
		{PROC_MOTIONCORR, PROC_MOTIONCORR_LABELNEW},
		{PROC_CTFFIND, PROC_CTFFIND_LABELNEW},
		{PROC_MANUALPICK, PROC_MANUALPICK_LABELNEW},
		{PROC_AUTOPICK, PROC_AUTOPICK_LABELNEW},
		{PROC_EXTRACT, PROC_EXTRACT_LABELNEW},
		{PROC_CLASSSELECT, PROC_CLASSSELECT_LABELNEW},
		{PROC_2DCLASS, PROC_2DCLASS_LABELNEW},
		{PROC_3DCLASS, PROC_3DCLASS_LABELNEW},
		{PROC_3DAUTO, PROC_3DAUTO_LABELNEW},
		{PROC_MASKCREATE,        PROC_MASKCREATE_LABELNEW},
		{PROC_JOINSTAR,        PROC_JOINSTAR_LABELNEW},
		{PROC_SUBTRACT,        PROC_SUBTRACT_LABELNEW},
		{PROC_POST,             PROC_POST_LABELNEW},
		{PROC_RESMAP,           PROC_RESMAP_LABELNEW},
		{PROC_INIMODEL,         PROC_INIMODEL_LABELNEW},
		{PROC_MULTIBODY,        PROC_MULTIBODY_LABELNEW},
		{PROC_MOTIONREFINE,     PROC_MOTIONREFINE_LABELNEW},
		{PROC_CTFREFINE,        PROC_CTFREFINE_LABELNEW},
		{PROC_TOMO_IMPORT,      PROC_TOMO_IMPORT_LABELNEW},
		{PROC_TOMO_SUBTOMO,     PROC_TOMO_SUBTOMO_LABELNEW},
		{PROC_TOMO_CTFREFINE,   PROC_TOMO_CTFREFINE_LABELNEW},
		{PROC_TOMO_ALIGN,       PROC_TOMO_ALIGN_LABELNEW},
		{PROC_TOMO_RECONSTRUCT, PROC_TOMO_RECONSTRUCT_LABELNEW},
		{PROC_EXTERNAL,         PROC_EXTERNAL_LABELNEW}};

static std::map<std::string, int> proc_dirname2type = {
		{PROC_IMPORT_DIRNAME,           PROC_IMPORT},
		{PROC_MOTIONCORR_DIRNAME,       PROC_MOTIONCORR},
		{PROC_CTFFIND_DIRNAME,          PROC_CTFFIND},
		{PROC_MANUALPICK_DIRNAME,       PROC_MANUALPICK},
		{PROC_AUTOPICK_DIRNAME,         PROC_AUTOPICK},
		{PROC_EXTRACT_DIRNAME,          PROC_EXTRACT},
		{PROC_CLASSSELECT_DIRNAME,      PROC_CLASSSELECT},
		{PROC_2DCLASS_DIRNAME,          PROC_2DCLASS},
		{PROC_3DCLASS_DIRNAME,          PROC_3DCLASS},
		{PROC_3DAUTO_DIRNAME,           PROC_3DAUTO},
		{PROC_MASKCREATE_DIRNAME,       PROC_MASKCREATE},
		{PROC_JOINSTAR_DIRNAME,         PROC_JOINSTAR},
		{PROC_SUBTRACT_DIRNAME,         PROC_SUBTRACT},
		{PROC_POST_DIRNAME,             PROC_POST},
		{PROC_RESMAP_DIRNAME,           PROC_RESMAP},
		{PROC_INIMODEL_DIRNAME,         PROC_INIMODEL},
		{PROC_MULTIBODY_DIRNAME,        PROC_MULTIBODY},
		{PROC_MOTIONREFINE_DIRNAME,     PROC_MOTIONREFINE},
		{PROC_CTFREFINE_DIRNAME,        PROC_CTFREFINE},
		{PROC_TOMO_IMPORT_DIRNAME,      PROC_TOMO_IMPORT},
		{PROC_TOMO_SUBTOMO_DIRNAME,     PROC_TOMO_SUBTOMO},
		{PROC_TOMO_CTFREFINE_DIRNAME,   PROC_TOMO_CTFREFINE},
		{PROC_TOMO_ALIGN_DIRNAME,       PROC_TOMO_ALIGN},
		{PROC_TOMO_RECONSTRUCT_DIRNAME, PROC_TOMO_RECONSTRUCT},
		{PROC_EXTERNAL_DIRNAME,         PROC_EXTERNAL}};

static std::map<std::string, int> proc_labelnew2type = {
		{PROC_IMPORT_LABELNEW,           PROC_IMPORT},
		{PROC_MOTIONCORR_LABELNEW,       PROC_MOTIONCORR},
		{PROC_CTFFIND_LABELNEW,          PROC_CTFFIND},
		{PROC_MANUALPICK_LABELNEW,       PROC_MANUALPICK},
		{PROC_AUTOPICK_LABELNEW,         PROC_AUTOPICK},
		{PROC_EXTRACT_LABELNEW,          PROC_EXTRACT},
		{PROC_CLASSSELECT_LABELNEW,      PROC_CLASSSELECT},
		{PROC_2DCLASS_LABELNEW,          PROC_2DCLASS},
		{PROC_3DCLASS_LABELNEW,          PROC_3DCLASS},
		{PROC_3DAUTO_LABELNEW,           PROC_3DAUTO},
		{PROC_MASKCREATE_LABELNEW,       PROC_MASKCREATE},
		{PROC_JOINSTAR_LABELNEW,         PROC_JOINSTAR},
		{PROC_SUBTRACT_LABELNEW,         PROC_SUBTRACT},
		{PROC_POST_LABELNEW,             PROC_POST},
		{PROC_RESMAP_LABELNEW,           PROC_RESMAP},
		{PROC_INIMODEL_LABELNEW,         PROC_INIMODEL},
		{PROC_MULTIBODY_LABELNEW,        PROC_MULTIBODY},
		{PROC_MOTIONREFINE_LABELNEW,     PROC_MOTIONREFINE},
		{PROC_CTFREFINE_LABELNEW,        PROC_CTFREFINE},
		{PROC_TOMO_IMPORT_LABELNEW,      PROC_TOMO_IMPORT},
		{PROC_TOMO_SUBTOMO_LABELNEW,     PROC_TOMO_SUBTOMO},
		{PROC_TOMO_CTFREFINE_LABELNEW,   PROC_TOMO_CTFREFINE},
		{PROC_TOMO_ALIGN_LABELNEW,       PROC_TOMO_ALIGN},
		{PROC_TOMO_RECONSTRUCT_LABELNEW, PROC_TOMO_RECONSTRUCT},
		{PROC_EXTERNAL_LABELNEW,         PROC_EXTERNAL}};


static std::string get_proc_label(int type)
{
	if (proc_type2labelnew.find(type) != proc_type2labelnew.end())
	{
		return proc_type2labelnew.at(type);;
	}
	else
	{
		std::cerr << " WARNING: unrecognised process type: " << type << std::endl;
	}
	return "";
}

static int get_proc_type(std::string label)
{

	FileName mylabel = label;
	bool keep_trying = true;
	while (keep_trying)
	{
		if (proc_labelnew2type.find(mylabel) != proc_labelnew2type.end())
		{
			return proc_labelnew2type.at(mylabel);;
		}
		else
		{
			// remove deepest layers of the process label, from ccpem-pipeliner
			if (mylabel.contains("."))
			{
				mylabel = mylabel.beforeLastOf(".");
			}
			else
			{
				keep_trying = false;
			}
		}
	}

	// Backward compatibility with alpha-versions of relion-4.0
	if (proc_dirname2type.find(label) != proc_dirname2type.end())
	{
		return proc_dirname2type.at(label);
	}
	else
	{
		std::cerr << " WARNING: unrecognised process label: " << label << std::endl;
	}
	return -1;

}


// Status a Process may have
#define PROC_RUNNING          0 // (hopefully) running
#define PROC_SCHEDULED        1 // scheduled for future execution
#define PROC_FINISHED_SUCCESS 2 // successfully finished
#define PROC_FINISHED_FAILURE 3 // reported an error
#define PROC_FINISHED_ABORTED 4 // aborted by the user

static std::map<std::string, int> procstatus_label2type = {
		{"Running", PROC_RUNNING},
		{"Scheduled", PROC_SCHEDULED},
		{"Succeeded", PROC_FINISHED_SUCCESS},
		{"Failed", PROC_FINISHED_FAILURE},
		{"Aborted", PROC_FINISHED_ABORTED}};

static std::map<int, std::string> procstatus_type2label = {
		{PROC_RUNNING, "Running", },
		{PROC_SCHEDULED, "Scheduled", },
		{PROC_FINISHED_SUCCESS, "Succeeded"},
		{PROC_FINISHED_FAILURE, "Failed"},
		{PROC_FINISHED_ABORTED, "Aborted"}};


#define HAS_NOT 0
#define HAS_OPTIONAL 1
#define HAS_COMPULSORY 2

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
	std::string type; // which type of node am I
	std::vector<long int> inputForProcessList; 	  //list of processes that use this Node as input
	long int outputFromProcess;   //Which process made this Node

	// Constructor
	Node(std::string _name, std::string _type)
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
	// Get HealPix order from string. Returns -1 when failed.
	static int getHealPixOrder(std::string s);

	// Get a f/p/m character for CTF fitting. Returns "" when failed.
	static std::string getCtfFitString(std::string option);

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
	std::string node_type;
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

	// Whats is my type-label? (can have deeper levels in ccpem-pipeliner)
	FileName label;

	// Is this a continuation job?
	bool is_continue;

	// Activate custom stuff for tomo
	bool is_tomo;

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
		is_tomo = false;
	}

	// Returns true if the option is present in joboptions
	bool containsLabel(std::string label, std::string &option);

	// Set this option in the job
	void setOption(std::string setOptionLine);

	// Activate specific options for tomo
	void setTomo(bool _is_tomo)
	{
		is_tomo = _is_tomo;
	}

	// write/read settings to disc
	// fn is a directory name (e.g. Refine3D/job123/) or a STAR file
	bool read(std::string fn, bool &_is_continue, bool do_initialise = false); // return false if unsuccessful
	void write(std::string fn);

	// Write the job submission script
	bool saveJobSubmissionScript(std::string newfilename, std::string outputname, std::vector<std::string> commands, std::string &error_message);

	// Initialise pipeline stuff for each job, return outputname
	void initialisePipeline(std::string &outputname, int job_counter);

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


	// relion-3.2: add subtomogram averaging programs by Jasenko
	void addTomoInputOptions(bool has_tomograms, bool has_particles,
							 bool has_trajectories, bool has_manifolds, bool has_halfmaps, bool has_postprocess);
	std::string getTomoInputCommmand(std::string &command, int has_tomograms, int has_particles,
									 int has_trajectories, int has_manifolds, bool has_halfmaps, int has_postprocess);

	std::string setTomoOutputCommand(std::string &command, std::string optimisationSet, std::string tomograms,
									 std::string particles, std::string trajectories, std::string manifolds,
									 std::string halfmap1, std::string postprocess, std::string refmask,
									 std::string optimisationSetOut);

	void initialiseTomoImportJob();
	bool getCommandsTomoImportJob(std::string &outputname, std::vector<std::string> &commands,
								  std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseTomoSubtomoJob();
	bool getCommandsTomoSubtomoJob(std::string &outputname, std::vector<std::string> &commands,
								   std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseTomoCtfRefineJob();
	bool getCommandsTomoCtfRefineJob(std::string &outputname, std::vector<std::string> &commands,
									 std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseTomoAlignJob();
	bool getCommandsTomoAlignJob(std::string &outputname, std::vector<std::string> &commands,
								 std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

	void initialiseTomoReconPartJob();
	bool getCommandsTomoReconPartJob(std::string &outputname, std::vector<std::string> &commands,
									 std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

};

#endif /* SRC_PIPELINE_JOBS_H_ */
