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
#define DEFAULTARETOMOLOCATION "/public/EM/AreTomo/AreTomo2/AreTomo2"
#define DEFAULTBATCHTOMOLOCATION "/public/EM/imod/IMOD/bin/batchruntomo"
//#define DEFAULTIMODWRAPPERLOCATION "relion_tomo_align_tilt_series"
//#define DEFAULTDENOISINGWRAPPERLOCATION "relion_tomo_denoise"
#define DEFAULTMOTIONCOR2LOCATION "/public/EM/MOTIONCOR2/MotionCor2"
#define DEFAULTGCTFLOCATION "/public/EM/Gctf/bin/Gctf"
#define DEFAULTRESMAPLOCATION "/public/EM/ResMap/ResMap-1.1.4-linux64"
#define DEFAULTQSUBCOMMAND "sbatch"
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
    "Set of tiltseries STAR file (.star)",
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

static const std::vector<std::string> job_tomo_pick_mode{
"particles",
"spheres",
"surfaces",
"filaments"
};

static const std::vector<std::string> job_modelangelo_alphabet_options{
"amino",
"DNA",
"RNA"
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

// General nodes used for input into relion jobs have a _CPIPE suffix

#define LABEL_MOVIES_CPIPE             "MicrographMovieGroupMetadata.star.relion"
#define LABEL_MICS_CPIPE               "MicrographGroupMetadata.star.relion"
#define LABEL_2DIMGS_CPIPE             "Image2DGroupMetadata.star.relion"
#define LABEL_MAP_CPIPE                "DensityMap.mrc"
#define LABEL_PARTS_CPIPE              "ParticleGroupMetadata.star.relion"
#define LABEL_COORDS_CPIPE             "MicrographCoordsGroup.star.relion"
#define LABEL_COORDS_HELIX_CPIPE       "MicrographCoordsGroup.star.relion.helixstartend"
#define LABEL_PARTS_HELIX_CPIPE        "ParticleGroupMetadata.star.relion.helicalsegments"
#define LABEL_OPTIMISER_CPIPE          "OptimiserData.star.relion"
#define LABEL_MASK_CPIPE               "Mask3D.mrc"
#define LABEL_HALFMAP_CPIPE	           "DensityMap.mrc.halfmap"
#define LABEL_RESMAP_CPIPE             "Image3D.mrc.localresmap"
#define LABEL_LOGFILE_CPIPE            "LogFile.pdf.relion"
#define LABEL_SEQUENCE_CPIPE           "Sequence.fasta"
#define LABEL_SEQUENCEALIGNMENT_CPIPE  "SequenceAlignment.hmm"
#define LABEL_ATOMCOORDS_CPIPE         "AtomCoords.cif"
#define LABEL_TOMO_OPTSET_CPIPE        "TomoOptimisationSet.star.relion"
#define LABEL_TOMOGRAMS_CPIPE          "TomogramGroupMetadata.star.relion"
#define LABEL_TRAJECTORIES_CPIPE       "TomoTrajectoryData.star.relion"
#define LABEL_MANIFOLDS_CPIPE          "TomoManifoldData.star.relion"
#define LABEL_POSTPROCESS_CPIPE        "ProcessData.star.relion.postprocess"

// More specific output nodes are below

#define LABEL_IMPORT_MOVIES            "MicrographMovieGroupMetadata.star.relion"
#define LABEL_IMPORT_MICS              "MicrographGroupMetadata.star.relion"
#define LABEL_IMPORT_COORDS            "MicrographCoordsGroup.star.relion"
#define LABEL_IMPORT_PARTS             "ParticleGroupMetadata.star.relion"
#define LABEL_IMPORT_2DIMG             "Image2DGroupMetadata.star.relion"
#define LABEL_IMPORT_MAP               "DensityMap.mrc"
#define LABEL_IMPORT_MASK              "Mask3D.mrc"
#define LABEL_IMPORT_HALFMAP           "DensityMap.mrc.halfmap"
#define LABEL_MOCORR_MICS              "MicrographGroupMetadata.star.relion.motioncorr"
#define LABEL_MOCORR_LOG               "LogFile.pdf.relion.motioncorr"
#define LABEL_CTFFIND_MICS             "MicrographGroupMetadata.star.relion.ctf"
#define LABEL_CTFFIND_LOG              "LogFile.pdf.relion.ctffind"
#define LABEL_CTFFIND_POWER_SPECTRA    "Image2DGroupMetadata.star.relion.ctffind.power_spectra"
#define LABEL_MANPICK_MICS             "MicrographGroupMetadata.star.relion"
#define LABEL_MANPICK_COORDS           "MicrographCoordsGroup.star.relion.manualpick"
#define LABEL_MANPICK_COORDS_HELIX     "MicrographCoordsGroup.star.relion.manualpick.helixstartend"
#define LABEL_AUTOPICK_COORDS          "MicrographCoordsGroup.star.relion.autopick"
#define LABEL_AUTOPICK_LOG             "LogFile.pdf.relion.autopick"
#define LABEL_AUTOPICK_TOPAZMODEL      "ParamsData.sav.topaz.model" // to be added?
#define LABEL_AUTOPICK_MICS            "MicrographGroupMetadata.star.relion"
#define LABEL_EXTRACT_PARTS            "ParticleGroupMetadata.star.relion"
#define LABEL_EXTRACT_PARTS_HELIX      "ParticleGroupMetadata.star.relion.helicalsegments"
#define LABEL_EXTRACT_COORDS_HELIX     "MicrographCoordsGroup.star.relion.helixstartend"
#define LABEL_EXTRACT_PARTS_REEX       "ParticleGroupMetadata.star.relion.reextract"
#define LABEL_EXTRACT_COORDS_REEX      "MicrographCoordsGroup.star.relion.reextract"
#define LABEL_CLASS2D_PARTS            "ParticleGroupMetadata.star.relion.class2d"
#define LABEL_CLASS2D_OPT              "OptimiserData.star.relion.class2d"
#define LABEL_CLASS2D_PARTS_HELIX      "ParticleGroupMetadata.star.relion.class2d.helicalsegments"
#define LABEL_SELECT_MICS              "MicrographGroupMetadata.star.relion"
#define LABEL_SELECT_MOVS              "MicrographMovieGroupMetadata.star.relion"
#define LABEL_SELECT_PARTS             "ParticleGroupMetadata.star.relion"
#define LABEL_SELECT_OPT               "OptimiserData.star.relion.select"
#define LABEL_SELECT_CLAVS             "Image2DGroupMetadata.star.relion.classaverages"
#define LABEL_SELECT_LOG               "LogFile.pdf.relion.select"
#define LABEL_INIMOD_MAP               "DensityMap.mrc.relion.initialmodel"
#define LABEL_INIMOD_OPTSET            "TomoOptimisationSet.star.relion.initialmodel"
#define LABEL_CLASS3D_OPT              "OptimiserData.star.relion.class3d"
#define LABEL_CLASS3D_MAP              "DensityMap.mrc.relion.class3d"
#define LABEL_CLASS3D_PARTS            "ParticleGroupMetadata.star.relion.class3d"
#define LABEL_CLASS3D_PARTS_HELIX      "ParticleGroupMetadata.star.relion.class3d.helicalsegments"
#define LABEL_CLASS3D_OPTSET           "TomoOptimisationSet.star.relion.class3d"
#define LABEL_REFINE3D_HALFMAP         "DensityMap.mrc.relion.halfmap.refine3d"
#define LABEL_REFINE3D_OPT             "OptimiserData.star.relion.refine3d"
#define LABEL_REFINE3D_MAP             "DensityMap.mrc.relion.refine3d"
#define LABEL_REFINE3D_PARTS           "ParticleGroupMetadata.star.relion.refine3d"
#define LABEL_REFINE3D_PARTS_HELIX     "ParticleGroupMetadata.star.relion.refine3d.helicalsegements"
#define LABEL_REFINE3D_OPTSET          "TomoOptimisationSet.star.relion.refine3d"
#define LABEL_MULTIBODY_HALFMAP        "DensityMap.mrc.relion.halfmap.multibody"
#define LABEL_MULTIBODY_PARTS          "ParticleGroupMetadata.star.relion.multibody"
#define LABEL_MULTIBODY_OPT            "OptimiserData.star.relion.multibody"
#define LABEL_MULTIBODY_FLEXLOG        "LogFile.pdf.relion.flexanalysis"
#define LABEL_MULTIBODY_SEL_PARTS      "ParticleGroupMetadata.star.relion.flexanalysis.eigenselected"
#define LABEL_MULTIBODY_OPTSET         "TomoOptimisationSet.star.relion.multibody"
#define LABEL_DYNAMIGHT_HALFMAP        "DensityMap.mrc.relion.halfmap.dynamight"
#define LABEL_MASK3D_MASK              "Mask3D.mrc.relion"
#define LABEL_SUBTRACT_SUBTRACTED      "ParticleGroupMetadata.star.relion.subtracted"
#define LABEL_SUBTRACT_REVERTED        "ParticleGroupMetadata.star.relion"
#define LABEL_LOCRES_OWN               "Image3D.mrc.relion.localresmap"
#define LABEL_LOCRES_RESMAP            "Image3D.mrc.resmap.localresmap"
#define LABEL_LOCRES_FILTMAP           "DensityMap.mrc.relion.localresfiltered"
#define LABEL_LOCRES_LOG               "LogFile.pdf.relion.localres"
#define LABEL_CTFREFINE_REFINEPARTS    "ParticleGroupMetadata.star.relion.ctfrefine"
#define LABEL_CTFREFINE_LOG            "LogFile.pdf.relion.ctfrefine"
#define LABEL_CTFREFINE_ANISOPARTS     "ParticleGroupMetadata.star.relion.anisomagrefine"
#define LABEL_POLISH_PARTS             "ParticleGroupMetadata.star.relion.polished"
#define LABEL_POLISH_LOG               "LogFile.pdf.relion.polish"
#define LABEL_POLISH_PARAMS            "ParamsData.txt.relion.polish"
#define LABEL_POST_POST                "ProcessData.star.relion.postprocess"
#define LABEL_POST_MAP                 "DensityMap.mrc.relion.postprocess"
#define LABEL_POST_MASKED              "DensityMap.mrc.relion.postprocess.masked"
#define LABEL_POST_LOG                 "LogFile.pdf.relion.postprocess"

// Tomography-specific jobs
#define LABEL_IMPORT_TOMOGRAMS         "TomogramGroupMetadata.star.relion.tomo.import"
#define LABEL_IMPORT_TOMO_COORDS       "ParticleGroupMetadata.star.relion.tomo.import"
#define LABEL_MOCORR_TOMOGRAMS         "TomogramGroupMetadata.star.relion.tomo.motioncorr"
#define LABEL_CTFFIND_TOMOGRAMS        "TomogramGroupMetadata.star.relion.tomo.ctffind"
#define LABEL_TILTALIGN_TOMOGRAMS      "TomogramGroupMetadata.star.relion.tomo.aligntiltseries"
#define LABEL_TILTALIGN_LOG            "LogFile.pdf.relion.tomo.aligntiltseries"
#define LABEL_RECONSTRUCT_TOMOGRAMS    "TomogramGroupMetadata.star.relion.tomo.reconstruct"
#define LABEL_DENOISE_TOMOGRAMS        "TomogramGroupMetadata.star.relion.tomo.denoise"
#define LABEL_TOMOPICK_PARTS_PARTS     "ParticleGroupMetadata.star.relion.tomo.manualpick.particles"
#define LABEL_TOMOPICK_PARTS_SPHERE    "ParticleGroupMetadata.star.relion.tomo.manualpick.spheres"
#define LABEL_TOMOPICK_PARTS_FILAMENT  "ParticleGroupMetadata.star.relion.tomo.manualpick.filaments"
#define LABEL_TOMOPICK_PARTS_SURFACE   "ParticleGroupMetadata.star.relion.tomo.manualpick.surfaces"
#define LABEL_TOMOPICK_OPTSET          "TomoOptimisationSet.star.relion.tomo.manualpick"
#define LABEL_EXCLUDE_TOMOGRAMS        "TomogramGroupMetadata.star.relion.tomo.excludeimages"
#define LABEL_SUBTOMO_PARTS            "ParticleGroupMetadata.star.relion.tomo.extract"
#define LABEL_SUBTOMO_OPTSET           "TomoOptimisationSet.star.relion.tomo.extract"
#define LABEL_CTFREFINE_TOMOGRAMS      "TomogramGroupMetadata.star.relion.tomo.ctfrefine"
#define LABEL_CTFREFINE_OPTSET         "TomoOptimisationSet.star.relion.tomo.ctfrefine"
#define LABEL_CTFREFINE_TOMO_LOG       "LogFile.pdf.relion.tomo.ctfrefine"
#define LABEL_FRAMEALIGN_TOMOGRAMS     "TomogramGroupMetadata.star.relion.tomo.polish"
#define LABEL_FRAMEALIGN_OPTSET        "TomoOptimisationSet.star.relion.tomo.polish"
#define LABEL_FRAMEALIGN_LOG           "LogFile.pdf.relion.tomo.polish"
#define LABEL_FRAMEALIGN_PARTS         "ParticleGroupMetadata.star.relion.tomo.polish"
#define LABEL_FRAMEALIGN_TRAJS         "TomoTrajectoryData.star.relion.polish"
#define LABEL_RECONSPART_OPTSET        "TomoOptimisationSet.star.relion.tomo.reconstruct"
#define LABEL_RECONSPART_HALFMAP       "DensityMap.mrc.relion.tomo.halfmap.reconstruct"
#define LABEL_RECONSPART_MAP           "DensityMap.mrc.relion.tomo.map.reconstruct"


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
#define PROC_DYNAMIGHT_DIRNAME        "DynaMight"    // Johannes' DynaMight for modelling continuous heterogeneity
#define PROC_MODELANGELO_DIRNAME      "ModelAngelo"  // Kiarash's ModelAngelo for automated model building
#define PROC_TOMO_IMPORT_DIRNAME      "Import"              // Import for tomography GUI
#define PROC_TOMO_SUBTOMO_DIRNAME     "Extract"           // Creation of pseudo-subtomograms from tilt series images
#define PROC_TOMO_CTFREFINE_DIRNAME   "CtfRefine"           // CTF refinement (defocus & aberrations) for tomography
#define PROC_TOMO_EXCLUDE_TILT_IMAGES_DIRNAME   "ExcludeTiltImages"  // Exclusion of bad tilt-images from tilt-series
#define PROC_TOMO_ALIGN_DIRNAME       "Polish"  // Frame alignment and particle polishing for subtomography
#define PROC_TOMO_RECONSTRUCT_DIRNAME "Reconstruct" // Calculation of particle average from the individual tilt series images
#define PROC_TOMO_DENOISE_DIRNAME "Denoise" // Denoise tomograms
#define PROC_TOMO_PICK_DIRNAME "Picks" // Pick particles in tomograms
#define PROC_EXTERNAL_DIRNAME         "External"     // For running non-relion programs
#define PROC_TOMO_ALIGN_TILTSERIES_DIRNAME "AlignTiltSeries" // Tilt series alignment for tomogram reconstruction
#define PROC_TOMO_RECONSTRUCT_TOMOGRAM_DIRNAME "Tomograms" // Reconstruction of tomograms for particle picking

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
#define PROC_DYNAMIGHT_LABELNEW        "dynamight"           // Johannes' DynaMight for modelling continuous heterogeneity
#define PROC_MODELANGELO_LABELNEW      "modelangelo"         // Kiarash's ModelAngelo for automated model building
#define PROC_TOMO_IMPORT_LABELNEW      "relion.importtomo"              // Import for tomography GUI
#define PROC_TOMO_EXCLUDE_TILT_IMAGES_LABELNEW "relion.excludetilts"    // Exclude bad tilt-images
#define PROC_TOMO_SUBTOMO_LABELNEW     "relion.pseudosubtomo"           // Creation of pseudo-subtomograms from tilt series images
#define PROC_TOMO_CTFREFINE_LABELNEW   "relion.ctfrefinetomo"           // CTF refinement (defocus & aberrations) for tomography
#define PROC_TOMO_ALIGN_LABELNEW       "relion.framealigntomo"          // Frame alignment and particle polishing for subtomography
#define PROC_TOMO_RECONSTRUCT_LABELNEW "relion.reconstructparticletomo" // Calculation of particle average from the individual tilt series images
#define PROC_TOMO_DENOISE_LABELNEW     "relion.denoisetomo"  // Denoise tomograms
#define PROC_TOMO_PICK_LABELNEW        "relion.picktomo"     // Pick tomograms
#define PROC_EXTERNAL_LABELNEW         "relion.external"     // For running non-relion programs
#define PROC_TOMO_ALIGN_TILTSERIES_LABELNEW "relion.aligntiltseries" // Tilt series alignment for tomogram reconstruction
#define PROC_TOMO_RECONSTRUCT_TOMOGRAM_LABELNEW "relion.reconstructtomograms" // Reconstruction of tomograms for particle picking


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
#define PROC_DYNAMIGHT      22// wrapper to Johannes' DynaMight
#define PROC_MODELANGELO    23// wrapper to Kiasrash's ModelAngelo
#define PROC_TOMO_IMPORT    50// Import for tomography GUI
#define PROC_TOMO_SUBTOMO   51// Creation of pseudo-subtomograms from tilt series images
#define PROC_TOMO_CTFREFINE     52// CTF refinement (defocus & aberrations for tomography)
#define PROC_TOMO_ALIGN        53// Frame alignment and particle polishing for subtomography
#define PROC_TOMO_RECONSTRUCT       54// Calculation of particle average from the individual tilt series images
#define PROC_TOMO_ALIGN_TILTSERIES 55// Tilt series alignment for tomogram reconstruction
#define PROC_TOMO_RECONSTRUCT_TOMOGRAM 56 // Reconstruction of tomograms for particle picking
#define PROC_TOMO_EXCLUDE_TILT_IMAGES 57 // Exclude bad tilt-images from tilt-series
#define PROC_TOMO_DENOISE_TOMOGRAM 58 // Denoise tomograms
#define PROC_TOMO_PICK_TOMOGRAM 59 // Denoise tomograms
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
        {PROC_DYNAMIGHT,        PROC_DYNAMIGHT_DIRNAME},
        {PROC_MODELANGELO,      PROC_MODELANGELO_DIRNAME},
		{PROC_TOMO_IMPORT,      PROC_TOMO_IMPORT_DIRNAME},
		{PROC_TOMO_SUBTOMO,     PROC_TOMO_SUBTOMO_DIRNAME},
		{PROC_TOMO_CTFREFINE,   PROC_TOMO_CTFREFINE_DIRNAME},
		{PROC_TOMO_ALIGN,       PROC_TOMO_ALIGN_DIRNAME},
		{PROC_TOMO_RECONSTRUCT, PROC_TOMO_RECONSTRUCT_DIRNAME},
        {PROC_TOMO_ALIGN_TILTSERIES,     PROC_TOMO_ALIGN_TILTSERIES_DIRNAME},
        {PROC_TOMO_RECONSTRUCT_TOMOGRAM, PROC_TOMO_RECONSTRUCT_TOMOGRAM_DIRNAME},
 	    {PROC_TOMO_DENOISE_TOMOGRAM, PROC_TOMO_DENOISE_DIRNAME},
        {PROC_TOMO_PICK_TOMOGRAM, PROC_TOMO_PICK_DIRNAME},
        {PROC_TOMO_EXCLUDE_TILT_IMAGES, PROC_TOMO_EXCLUDE_TILT_IMAGES_DIRNAME},
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
        {PROC_DYNAMIGHT,        PROC_DYNAMIGHT_LABELNEW},
        {PROC_MODELANGELO,      PROC_MODELANGELO_LABELNEW},
		{PROC_TOMO_IMPORT,      PROC_TOMO_IMPORT_LABELNEW},
		{PROC_TOMO_SUBTOMO,     PROC_TOMO_SUBTOMO_LABELNEW},
		{PROC_TOMO_CTFREFINE,   PROC_TOMO_CTFREFINE_LABELNEW},
		{PROC_TOMO_ALIGN,       PROC_TOMO_ALIGN_LABELNEW},
		{PROC_TOMO_RECONSTRUCT, PROC_TOMO_RECONSTRUCT_LABELNEW},
        {PROC_TOMO_ALIGN_TILTSERIES,     PROC_TOMO_ALIGN_TILTSERIES_LABELNEW},
        {PROC_TOMO_RECONSTRUCT_TOMOGRAM, PROC_TOMO_RECONSTRUCT_TOMOGRAM_LABELNEW},
 	    {PROC_TOMO_DENOISE_TOMOGRAM, PROC_TOMO_DENOISE_LABELNEW},
 	    {PROC_TOMO_PICK_TOMOGRAM, PROC_TOMO_PICK_LABELNEW},
	    {PROC_TOMO_EXCLUDE_TILT_IMAGES, PROC_TOMO_EXCLUDE_TILT_IMAGES_LABELNEW},
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
        {PROC_DYNAMIGHT_DIRNAME,        PROC_DYNAMIGHT},
        {PROC_MODELANGELO_DIRNAME,      PROC_MODELANGELO},
		{PROC_TOMO_IMPORT_DIRNAME,      PROC_TOMO_IMPORT},
		{PROC_TOMO_SUBTOMO_DIRNAME,     PROC_TOMO_SUBTOMO},
		{PROC_TOMO_CTFREFINE_DIRNAME,   PROC_TOMO_CTFREFINE},
		{PROC_TOMO_ALIGN_DIRNAME,       PROC_TOMO_ALIGN},
		{PROC_TOMO_RECONSTRUCT_DIRNAME, PROC_TOMO_RECONSTRUCT},
        {PROC_TOMO_ALIGN_TILTSERIES_DIRNAME,     PROC_TOMO_ALIGN_TILTSERIES},
        {PROC_TOMO_RECONSTRUCT_TOMOGRAM_DIRNAME, PROC_TOMO_RECONSTRUCT_TOMOGRAM},
        {PROC_TOMO_DENOISE_DIRNAME, PROC_TOMO_DENOISE_TOMOGRAM},
        {PROC_TOMO_PICK_DIRNAME, PROC_TOMO_PICK_TOMOGRAM},
	    {PROC_TOMO_EXCLUDE_TILT_IMAGES_DIRNAME, PROC_TOMO_EXCLUDE_TILT_IMAGES},
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
        {PROC_DYNAMIGHT_LABELNEW,        PROC_DYNAMIGHT},
        {PROC_MODELANGELO_LABELNEW,      PROC_MODELANGELO},
		{PROC_TOMO_IMPORT_LABELNEW,      PROC_TOMO_IMPORT},
		{PROC_TOMO_SUBTOMO_LABELNEW,     PROC_TOMO_SUBTOMO},
		{PROC_TOMO_CTFREFINE_LABELNEW,   PROC_TOMO_CTFREFINE},
		{PROC_TOMO_ALIGN_LABELNEW,       PROC_TOMO_ALIGN},
		{PROC_TOMO_RECONSTRUCT_LABELNEW, PROC_TOMO_RECONSTRUCT},
        {PROC_TOMO_ALIGN_TILTSERIES_LABELNEW,     PROC_TOMO_ALIGN_TILTSERIES},
        {PROC_TOMO_RECONSTRUCT_TOMOGRAM_LABELNEW, PROC_TOMO_RECONSTRUCT_TOMOGRAM},
 	    {PROC_TOMO_DENOISE_LABELNEW, PROC_TOMO_DENOISE_TOMOGRAM},
 	    {PROC_TOMO_PICK_LABELNEW, PROC_TOMO_PICK_TOMOGRAM},
	    {PROC_TOMO_EXCLUDE_TILT_IMAGES_LABELNEW, PROC_TOMO_EXCLUDE_TILT_IMAGES},
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
    int type_depth; // how far deep does the .Nodes directory go
	std::vector<long int> inputForProcessList; 	  //list of processes that use this Node as input
	long int outputFromProcess;   //Which process made this Node

	// Constructor
	Node(std::string _name, std::string _type, int _type_depth = 1)
	{
		name = _name;
		type = _type;
        type_depth = _type_depth;
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
std::vector<Node> getOutputNodesRefine(std::string outputname, std::string jobtype, int iter, int K, int dim, int nr_bodies=1, bool _is_tomo = false);

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
    int node_type_depth;
	std::string pattern;
	std::string directory;
	std::vector<std::string> radio_options;

public:

	// Any constructor
	JobOption(std::string _label, std::string _default_value, std::string _helptext);

	// FileName constructor
	JobOption(std::string _label, std::string  _default_value, std::string _pattern, std::string _directory, std::string _helptext);

	// InputNode constructor
	JobOption(std::string _label, std::string _nodetype, int _node_type_depth, std::string _default_value, std::string _pattern, std::string _helptext);

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
			bool do_makedir, std::string &warning_message, bool do_dash_for_python = false);

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

  	void initialiseDynaMightJob();
	bool getCommandsDynaMightJob(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

  	void initialiseModelAngeloJob();
	bool getCommandsModelAngeloJob(std::string &outputname, std::vector<std::string> &commands,
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
							 bool has_trajectories, bool has_manifolds=false);

	std::string getTomoInputCommmand(bool is_for_refine, std::string &command, int has_tomograms, int has_particles,
									 int has_trajectories, int has_manifolds=HAS_NOT);

	void initialiseTomoImportJob();
	bool getCommandsTomoImportJob(std::string &outputname, std::vector<std::string> &commands,
								  std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

    void initialiseTomoAlignTiltSeriesJob();
    bool getCommandsTomoAlignTiltSeriesJob(std::string &outputname, std::vector<std::string> &commands,
                                  std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

    void initialiseTomoReconstructTomogramsJob();
    bool getCommandsTomoReconstructTomogramsJob(std::string &outputname, std::vector<std::string> &commands,
                                  std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

    void initialiseTomoDenoiseTomogramsJob();
    bool getCommandsTomoDenoiseTomogramsJob(std::string &outputname, std::vector<std::string> &commands,
                                  std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

    void initialiseTomoPickTomogramsJob();
    bool getCommandsTomoPickTomogramsJob(std::string &outputname, std::vector<std::string> &commands,
                                            std::string &final_command, bool do_makedir, int job_counter, std::string &error_message);

    void initialiseTomoExcludeTiltImagesJob();
    bool getCommandsTomoExcludeTiltImagesJob(std::string &outputname, std::vector<std::string> &commands,
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
