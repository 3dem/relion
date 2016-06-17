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

#ifndef GUI_JOBWINDOW_H_
#define GUI_JOBWINDOW_H_

#define HAS_MPI true
#define HAS_NOT_MPI false
#define HAS_THREAD true
#define HAS_NOT_THREAD false
#define HAS_PARALLEL_DISCIO true
#define HAS_NOT_PARALLEL_DISCIO false

#include <ctime>
#include "src/gui_entries.h"
#include "src/metadata_table.h"
#include "src/pipeliner.h"

// Our own defaults at LMB are the hard-coded ones
#define DEFAULTQSUBLOCATION "/public/EM/RELION/relion/bin/qsub.csh"
#define DEFAULTCTFFINDLOCATION "\"/public/EM/ctffind/ctffind.exe  --omp-num-threads 1 --old-school-input\""
#define DEFAULTMOTIONCORRLOCATION "/public/EM/MOTIONCORR/bin/motioncorr"
#define DEFAULTUNBLURLOCATION "/public/EM/UNBLUR/unblur.exe"
#define DEFAULTSUMMOVIELOCATION "/public/EM/SUMMOVIE/summovie.exe"
#define DEFAULTGCTFLOCATION "/public/EM/Gctf/bin/Gctf"
#define DEFAULTRESMAPLOCATION "/public/EM/ResMap/ResMap-1.1.4-linux64"
#define DEFAULTMININIMUMDEDICATED 1
#define DEFAULTWARNINGLOCALMPI 32
#define DEFAULTALLOWCHANGEMINDEDICATED true

static Fl_Menu_Item sampling_options[] = {
		      {"30 degrees"},
		      {"15 degrees"},
		      {"7.5 degrees"},
		      {"3.7 degrees"},
		      {"1.8 degrees"},
		      {"0.9 degrees"},
		      {"0.5 degrees"},
		      {"0.2 degrees"},
		      {"0.1 degrees"},
		      {0} // this should be the last entry
};

static Fl_Menu_Item node_type_options[] = {
		      {"2D micrograph movies (*.mrcs)"},
	          {"2D micrographs/tomograms (*.mrc)"},
	          {"2D/3D particle coordinates (*.box, *_pick.star)"},
	          {"Particles STAR file (.star)"},
	          {"Movie-particles STAR file (.star)"},
	          {"2D references (.star or .mrcs)"},
	          {"Micrographs STAR file (.star)"},
		      {"3D reference (.mrc)"},
		      {"3D mask (.mrc)"},
		      {"Unfiltered half-map (unfil.mrc)"},
		      {0} // this should be the last entry
};

static int minimum_nr_dedicated;
static bool do_allow_change_minimum_dedicated;


// Output filenames of relion_refine (use iter=-1 for auto-refine)
std::vector<Node> getOutputNodesRefine(std::string outputname, int iter, int K, int dim, int nr_bodies = 1, bool do_movies = false, bool do_also_rot = false);

class RelionJobWindow : public Fl_Box
{
protected:

	// Which process type is this?
	int type;
	bool is_continue;

public:
	// Position of current cursor to place new elements
	int start_y, current_y;

	// Menu
	Fl_Menu_Bar *menubar;

	// Choice for new/continue
	Fl_Choice *choice_continue;
    Fl_Menu_ *menu_continue;

	// Tabs
    Fl_Tabs *tabs;
	Fl_Group *tab1, *tab2, *tab3, *tab4, *tab5, *tab6, *tab7, *runtab;

	// Running
	Fl_Group *queue_group;
	SliderEntry nr_mpi;
	SliderEntry nr_threads;
    BooleanEntry do_queue;
	AnyEntry queuename;
	AnyEntry qsub;
	SliderEntry min_dedicated;
	FileNameEntry qsubscript;
	bool have_extra1, have_extra2;
	AnyEntry qsub_extra1;
	AnyEntry qsub_extra2;
	AnyEntry other_args;

	// For the run tab
	bool has_mpi;
	bool has_thread;

	// General pipeline I/O (the same for all jobtypes)
	std::vector<Node> pipelineInputNodes;
	std::vector<Node> pipelineOutputNodes;
	std::string pipelineOutputName;


public:
	// Constructor with x, y, w, h and a title
	RelionJobWindow(int nr_tabs, bool _has_mpi, bool _has_thread,
			int x = WCOL0, int y = 2, int w = GUIWIDTH - WCOL0 - 10, int h = GUIHEIGHT_OLD-65, const char* title = "");

    // Destructor
    ~RelionJobWindow() {};

    void resetHeight();

    void setupRunTab();

    // by default nothing happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// write/read settings to disc
	void openWriteFile(std::string fn, std::ofstream &fh);
	bool openReadFile(std::string fn, std::ifstream &fh);
	void closeWriteFile(std::ofstream& fh, std::string fn);
	void closeReadFile(std::ifstream& fh);

	// Write the job submission script
	void saveJobSubmissionScript(std::string newfilename, std::string outputname, std::vector<std::string> commands);

	// Initialise pipeiline stuff for each job, return outputname
	void initialisePipeline(std::string &outputname, std::string defaultname, int job_counter);

	// Prepare the final (job submission or combined (mpi) command of possibly multiple lines)
	// Returns true to go ahead, and false to cancel
	bool prepareFinalCommand(std::string &outputname, std::vector<std::string> &commands, std::string &final_command, bool do_makedir = true);

private:

    static void cb_menu_continue(Fl_Widget*, void*);
    inline void cb_menu_continue_i();

};

/*
class XXXXJobWindow : public RelionJobWindow
{
public:

	// I/O

public:

	// Constructor
	XXXXJobWindow();

	// Destructor
	~XXXXJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	void getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command);

};
*/

class ImportJobWindow : public RelionJobWindow
{
public:

	// I/O
	FileNameEntry fn_in;
	RadioEntry node_type;

public:

	// Constructor
	ImportJobWindow();

	// Destructor
	~ImportJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command, bool do_makedir, int job_counter);

};

class MotioncorrJobWindow : public RelionJobWindow
{
public:

	InputNodeEntry input_star_mics;
	FileNameEntry fn_motioncorr_exe;
	FileNameEntry fn_summovie_exe;
	SliderEntry bin_factor;
	SliderEntry first_frame_ali;
	SliderEntry last_frame_ali;
	SliderEntry first_frame_sum;
	SliderEntry last_frame_sum;
	SliderEntry bfactor;
	AnyEntry other_motioncorr_args;
	BooleanEntry do_save_movies;

	BooleanEntry do_unblur;
	FileNameEntry fn_unblur_exe;
	SliderEntry angpix;

    AnyEntry gpu_ids;

    Fl_Group *unblur_group;

public:

	// Constructor
	MotioncorrJobWindow();

	// Destructor
	~MotioncorrJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	bool getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter);

};

class CtffindJobWindow : public RelionJobWindow
{
public:

	InputNodeEntry input_star_mics;
	FileNameEntry fn_ctffind_exe, fn_gctf_exe;
	SliderEntry ctf_win;
	SliderEntry cs, kv, q0, angpix, dast;
	SliderEntry box, resmin, resmax, dfmin, dfmax, dfstep;
	BooleanEntry is_ctffind4, do_movie_thon_rings, do_phaseshift;
	SliderEntry avg_movie_frames;
	AnyEntry phase_min, phase_max, phase_step;
	AnyEntry movie_rootname;
	BooleanEntry use_gctf, do_ignore_ctffind_params, do_EPA;
    AnyEntry gpu_ids;

	AnyEntry other_gctf_args;

	Fl_Group *gctf_group, *ctffind4_group, *movie_group, *phaseshift_group;

public:

	// Constructor
	CtffindJobWindow();

	// Destructor
	~CtffindJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	bool getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter);

};

class AutopickJobWindow : public RelionJobWindow
{
public:

	InputNodeEntry fn_input_autopick;
	InputNodeEntry fn_refs_autopick;
	BooleanEntry do_invert_refs;
	BooleanEntry do_ctf_autopick;
	BooleanEntry do_ignore_first_ctfpeak_autopick;
	SliderEntry lowpass;
	SliderEntry highpass;
	SliderEntry angpix;
	SliderEntry angpix_ref;
	SliderEntry particle_diameter;
	SliderEntry psi_sampling_autopick;
	SliderEntry shrink;
	BooleanEntry do_write_fom_maps, do_read_fom_maps;
	SliderEntry threshold_autopick;
	SliderEntry mindist_autopick;
	SliderEntry maxstddevnoise_autopick;
	BooleanEntry do_pick_helical_segments;
	SliderEntry helical_tube_kappa_max;
	SliderEntry helical_tube_outer_diameter;
	SliderEntry helical_tube_length_min;
	SliderEntry helical_nr_asu;
	SliderEntry helical_rise;
    BooleanEntry use_gpu;
    AnyEntry gpu_ids;

	Fl_Group *autopick_ctf_group, *autopick_helix_group, *gpu_group;

public:

	// Constructor
	AutopickJobWindow();

	// Destructor
	~AutopickJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	bool getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter);

};


class ManualpickJobWindow : public RelionJobWindow
{
public:

	InputNodeEntry fn_in;
	SliderEntry lowpass;
	SliderEntry highpass;
	SliderEntry angpix;
	SliderEntry diameter;
	SliderEntry micscale;
	SliderEntry ctfscale;
	SliderEntry sigma_contrast;
	SliderEntry white_val;
	SliderEntry black_val;
	BooleanEntry do_color;
	AnyEntry color_label;
	FileNameEntry fn_color;
	SliderEntry blue_value;
	SliderEntry red_value;

	Fl_Group *color_group;

public:

	// Constructor
	ManualpickJobWindow();

	// Destructor
	~ManualpickJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	bool getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter);

};


class ExtractJobWindow : public RelionJobWindow
{
public:

	// I/O
	InputNodeEntry star_mics;
	InputNodeEntry coords_suffix;
	BooleanEntry do_reextract;
	InputNodeEntry fndata_reextract;
	BooleanEntry do_recenter;

	// extract
	SliderEntry extract_size;
	BooleanEntry do_set_angpix;
	SliderEntry angpix;
	BooleanEntry do_rescale;
	SliderEntry rescale;
	BooleanEntry do_norm;
	SliderEntry bg_diameter;
	SliderEntry white_dust;
	SliderEntry black_dust;
	BooleanEntry do_invert;

	// Helix
	BooleanEntry do_extract_helix;
	BooleanEntry do_extract_helical_tubes;
	BooleanEntry do_cut_into_segments;
	SliderEntry helical_nr_asu;
	SliderEntry helical_rise;
	SliderEntry helical_tube_outer_diameter;
	BooleanEntry helical_bimodal_angular_priors;

	Fl_Group *reextract_group, *rescale_group, *set_angpix_group, *norm_group, *movie_extract_group, *helix_group, *helical_tubes_group;

public:

	// Constructor
	ExtractJobWindow();

	// Destructor
	~ExtractJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter);

};

class SortJobWindow : public RelionJobWindow
{
public:

	// I/O
	InputNodeEntry input_star;
	BooleanEntry is_autopick;
	InputNodeEntry autopick_refs;
	BooleanEntry do_ctf;
	BooleanEntry do_ignore_first_ctfpeak;

	Fl_Group *ctf_group, *autopick_group;

public:

	// Constructor
	SortJobWindow();

	// Destructor
	~SortJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter);

};

class Class2DJobWindow : public RelionJobWindow
{
public:

	// I/O
	FileNameEntry fn_cont;
	InputNodeEntry fn_img;

	// CTF
	BooleanEntry do_ctf_correction;
	BooleanEntry ctf_phase_flipped;
	BooleanEntry ctf_intact_first_peak;

	// Optimisation
	SliderEntry nr_iter;
	SliderEntry tau_fudge;
	SliderEntry particle_diameter;
	BooleanEntry do_zero_mask;
	SliderEntry highres_limit;

	// Sampling
	SliderEntry nr_classes;
	BooleanEntry dont_skip_align;
	SliderEntry psi_sampling;
	SliderEntry offset_range;
	SliderEntry offset_step;

	// Helix
	BooleanEntry do_helix;
	BooleanEntry do_bimodal_psi;
	SliderEntry range_psi;
	SliderEntry helical_tube_outer_diameter;

	// Compute
	BooleanEntry do_combine_thru_disc;
	BooleanEntry do_parallel_discio;
	SliderEntry nr_pool;
	BooleanEntry do_preread_images;
    BooleanEntry use_gpu;
    AnyEntry gpu_ids;
    AnyEntry scratch_dir;


	Fl_Group *ctf_group, *dont_skip_align_group, *helix_group, *gpu_group;

public:

	// Constructor
	Class2DJobWindow();

	// Destructor
	~Class2DJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter);

};


class Class3DJobWindow : public RelionJobWindow
{
public:

	// I/O
	FileNameEntry fn_cont;
	InputNodeEntry fn_img;
	InputNodeEntry fn_ref;
	InputNodeEntry fn_mask;
	SliderEntry nr_classes;

	// Reference
	BooleanEntry ref_correct_greyscale;
	SliderEntry ini_high;
	AnyEntry sym_name;

	// CTF
	BooleanEntry do_ctf_correction;
	BooleanEntry ctf_corrected_ref;
	BooleanEntry ctf_phase_flipped;
	BooleanEntry ctf_intact_first_peak;

	// Optimisation
	SliderEntry nr_iter;
	SliderEntry tau_fudge;
	SliderEntry particle_diameter;
	BooleanEntry do_zero_mask;
	SliderEntry highres_limit;

	// Sampling
	BooleanEntry dont_skip_align;
	RadioEntry sampling;
	SliderEntry offset_range;
	SliderEntry offset_step;
	BooleanEntry do_local_ang_searches;
	SliderEntry sigma_angles;

	// Helix
	BooleanEntry do_helix;
	AnyEntry helical_tube_inner_diameter;
	AnyEntry helical_tube_outer_diameter;
	SliderEntry helical_nr_asu;
	AnyEntry helical_twist_initial;
	AnyEntry helical_rise_initial;
	BooleanEntry do_local_search_helical_symmetry;
	AnyEntry helical_twist_min;
	AnyEntry helical_twist_max;
	AnyEntry helical_twist_inistep;
	AnyEntry helical_rise_min;
	AnyEntry helical_rise_max;
	AnyEntry helical_rise_inistep;
	SliderEntry helical_z_percentage;
	AnyEntry range_tilt;
	AnyEntry range_psi;
	SliderEntry helical_range_distance;

	// Compute
	BooleanEntry do_combine_thru_disc;
	BooleanEntry do_parallel_discio;
	SliderEntry nr_pool;
	BooleanEntry do_preread_images;
    BooleanEntry use_gpu;
    AnyEntry gpu_ids;
    AnyEntry scratch_dir;

	Fl_Group *ctf_group, *dont_skip_align_group, *localsearch_group, *helix_group, *helix_symmetry_search_group, *gpu_group;

public:

	// Constructor
	Class3DJobWindow();

	// Destructor
	~Class3DJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter);

};


class Auto3DJobWindow : public RelionJobWindow
{
public:

	// I/O
	FileNameEntry fn_cont;
	InputNodeEntry fn_img;
	InputNodeEntry fn_ref;
	InputNodeEntry fn_mask;

	// Reference
	BooleanEntry ref_correct_greyscale;
	SliderEntry ini_high;
	AnyEntry sym_name;

	// CTF
	BooleanEntry do_ctf_correction;
	BooleanEntry ctf_corrected_ref;
	BooleanEntry ctf_phase_flipped;
	BooleanEntry ctf_intact_first_peak;

	// Optimisation
	SliderEntry particle_diameter;
	BooleanEntry do_zero_mask;
	BooleanEntry do_solvent_fsc;

	// Sampling
	textOnlyEntry autosample_text;
	RadioEntry sampling;
	SliderEntry offset_range;
	SliderEntry offset_step;
	RadioEntry auto_local_sampling;

	// Helix
	BooleanEntry do_helix;
	AnyEntry helical_tube_inner_diameter;
	AnyEntry helical_tube_outer_diameter;
	SliderEntry helical_nr_asu;
	AnyEntry helical_twist_initial;
	AnyEntry helical_rise_initial;
	BooleanEntry do_local_search_helical_symmetry;
	AnyEntry helical_twist_min;
	AnyEntry helical_twist_max;
	AnyEntry helical_twist_inistep;
	AnyEntry helical_rise_min;
	AnyEntry helical_rise_max;
	AnyEntry helical_rise_inistep;
	SliderEntry helical_z_percentage;
	AnyEntry range_tilt;
	AnyEntry range_psi;
	SliderEntry helical_range_distance;

	// Compute
	BooleanEntry do_combine_thru_disc;
	BooleanEntry do_parallel_discio;
	SliderEntry nr_pool;
	BooleanEntry do_preread_images;
    BooleanEntry use_gpu;
    AnyEntry gpu_ids;
    AnyEntry scratch_dir;

	Fl_Group *ctf_group, *helix_group, *helix_symmetry_search_group, *gpu_group;

public:

	// Constructor
	Auto3DJobWindow();

	// Destructor
	~Auto3DJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter);

};

class MovieRefineJobWindow : public RelionJobWindow
{
public:

	// I/O
	FileNameEntry fn_cont;
	InputNodeEntry fn_movie_star;
	AnyEntry movie_rootname;

	// Extract movie-particles
	SliderEntry first_movie_frame;
	SliderEntry last_movie_frame;
	SliderEntry avg_movie_frames;
	SliderEntry max_mpi_nodes;
	SliderEntry extract_size;
	BooleanEntry do_set_angpix;
	SliderEntry angpix;
	BooleanEntry do_rescale;
	SliderEntry rescale;
	BooleanEntry do_norm;
	SliderEntry bg_diameter;
	SliderEntry white_dust;
	SliderEntry black_dust;
	BooleanEntry do_invert;

	// Movies
	SliderEntry join_nr_mics;
	SliderEntry movie_runavg_window;
	SliderEntry movie_sigma_offset;
	BooleanEntry do_alsorot_movies;
	SliderEntry movie_sigma_angles;

	Fl_Group *norm_group, *rescale_group, *alsorot_movie_group;

public:

	// Constructor
	MovieRefineJobWindow();

	// Destructor
	~MovieRefineJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, bool do_makedir, int job_counter);

};

class PolishJobWindow : public RelionJobWindow
{

public:
	// I/O
	InputNodeEntry fn_in;
	InputNodeEntry fn_mask;

	// Movements
	BooleanEntry do_fit_movement;
	SliderEntry sigma_nb;

	// Damage-model
	BooleanEntry do_bfactor_weighting;
	SliderEntry perframe_highres;
	SliderEntry perframe_bfac_lowres;
	SliderEntry average_frame_bfactor;
	AnyEntry sym_name;

	// Normalisation
	SliderEntry bg_diameter;
	SliderEntry white_dust;
	SliderEntry black_dust;

	// Helix
	BooleanEntry do_helix;
	SliderEntry helical_nr_asu;
	AnyEntry helical_twist;
	AnyEntry helical_rise;

	Fl_Group *fit_group, *weight_group, *helix_group;

public:
	// Constructor
	PolishJobWindow();

	// Destructor
	~PolishJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command, bool do_makedir, int job_counter);

};

class ClassSelectJobWindow : public RelionJobWindow
{
public:

	// I/O
	InputNodeEntry fn_model;
	InputNodeEntry fn_data;
	InputNodeEntry fn_mic;
	InputNodeEntry fn_coords;

	// Allow clearing of input fields for a new job
	FileName ori_fn_model, ori_fn_data, ori_fn_mic, ori_fn_coords;

	BooleanEntry do_recenter;
	BooleanEntry do_regroup;
	SliderEntry nr_groups;

	Fl_Group *regroup_group;

public:

	// Constructor
	ClassSelectJobWindow();

	// Destructor
	~ClassSelectJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command, bool do_makedir, int job_counter);

};

class MaskCreateJobWindow : public RelionJobWindow
{
public:

	// I/O
	InputNodeEntry fn_in;
	SliderEntry lowpass_filter;
	SliderEntry angpix;
	SliderEntry inimask_threshold;
	SliderEntry extend_inimask;
	SliderEntry width_mask_edge;
	BooleanEntry do_helix;
	SliderEntry helical_z_percentage;

	Fl_Group *helix_group;

public:

	// Constructor
	MaskCreateJobWindow();

	// Destructor
	~MaskCreateJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command, bool do_makedir, int job_counter);

};

class JoinStarJobWindow : public RelionJobWindow
{
public:

	// I/O
	BooleanEntry do_part;
	InputNodeEntry fn_part1;
	InputNodeEntry fn_part2;
	InputNodeEntry fn_part3;
	InputNodeEntry fn_part4;

	BooleanEntry do_mic;
	InputNodeEntry fn_mic1;
	InputNodeEntry fn_mic2;
	InputNodeEntry fn_mic3;
	InputNodeEntry fn_mic4;

	BooleanEntry do_mov;
	InputNodeEntry fn_mov1;
	InputNodeEntry fn_mov2;
	InputNodeEntry fn_mov3;
	InputNodeEntry fn_mov4;


	// Allow clearing of input fields for a new job
	FileName ori_fn_part1, ori_fn_part2, ori_fn_part3, ori_fn_part4;
	FileName ori_fn_mic1, ori_fn_mic2, ori_fn_mic3, ori_fn_mic4;
	FileName ori_fn_mov1, ori_fn_mov2, ori_fn_mov3, ori_fn_mov4;

	Fl_Group *part_group, *mic_group, *mov_group;

public:

	// Constructor
	JoinStarJobWindow();

	// Destructor
	~JoinStarJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command, bool do_makedir, int job_counter);

};

class SubtractJobWindow : public RelionJobWindow
{
public:

	// I/O
	InputNodeEntry fn_in;
	BooleanEntry do_subtract;
	InputNodeEntry fn_data;
	InputNodeEntry fn_mask;
	BooleanEntry do_fliplabel;

	// CTF
	BooleanEntry do_ctf_correction;
	BooleanEntry ctf_phase_flipped;
	BooleanEntry ctf_intact_first_peak;

	Fl_Group *ctf_group, *subtract_group;

public:

	// Constructor
	SubtractJobWindow();

	// Destructor
	~SubtractJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command, bool do_makedir, int job_counter);

};



class PostJobWindow : public RelionJobWindow
{
public:

	// I/O
	InputNodeEntry fn_in;
	InputNodeEntry fn_mask;
	SliderEntry angpix;

	// Sharpening
	BooleanEntry do_auto_bfac;
	SliderEntry autob_lowres;
	BooleanEntry do_adhoc_bfac;
	SliderEntry adhoc_bfac;
	FileNameEntry fn_mtf;

	// Filtering
	BooleanEntry do_skip_fsc_weighting;
	SliderEntry low_pass;

	Fl_Group *autobfac_group, *adhocbfac_group, *skipweight_group;

public:

	// Constructor
	PostJobWindow();

	// Destructor
	~PostJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command, bool do_makedir, int job_counter);

};

class ResmapJobWindow : public RelionJobWindow
{
public:

	// I/O
	InputNodeEntry fn_in;
	InputNodeEntry fn_mask;

	// Params
	FileNameEntry fn_resmap;
	SliderEntry angpix;
	SliderEntry pval;
	SliderEntry minres;
	SliderEntry maxres;
	SliderEntry stepres;


public:

	// Constructor
	ResmapJobWindow();

	// Destructor
	~ResmapJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	bool getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command, bool do_makedir, int job_counter);

};


#endif /* GUI_JOBWINDOW_H_ */
