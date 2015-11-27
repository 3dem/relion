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

#define MENUHEIGHT 30
#define TABHEIGHT 25
#define HAS_MPI true
#define HAS_NOT_MPI false
#define HAS_THREAD true
#define HAS_NOT_THREAD false
#define HAS_PARALLEL_DISCIO true
#define HAS_NOT_PARALLEL_DISCIO false
#define HAS_RUN true
#define HAS_NOT_RUN false

#include <ctime>
#include "src/gui_entries.h"
#include "src/pipeliner.h"

// Our own defaults at LMB are the hard-coded ones
#define DEFAULTQSUBLOCATION "/public/EM/RELION/relion/bin/qsub.csh"
#define DEFAULTCTFFINDLOCATION "\"/public/EM/ctffind/ctffind.exe  --omp-num-threads 1 --old-school-input\""
#define DEFAULTGCTFLOCATION "/public/EM/Gctf/bin/Gctf"
#define DEFAULTRESMAPLOCATION "/public/EM/ResMap/ResMap-1.1.4-linux64"
#define DEFAULTMININIMUMDEDICATED 1
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

static int minimum_nr_dedicated;
static bool do_allow_change_minimum_dedicated;


void getContinueOutname(std::string &outputname, FileNameEntry &fn_cont);
// Output filenames of relion_refine (use iter=-1 for auto-refine)
std::vector<Node> getOutputNodesRefine(std::string outputname, int iter, int K, int dim, int nr_bodies = 1);

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
	Fl_Group *tab1, *tab2, *tab3, *tab4, *tab5, *tab6, *runtab;

	// Running
	Fl_Group *queue_group;
	SliderEntry nr_mpi;
	SliderEntry nr_threads;
	SliderEntry ram_per_thread;
    BooleanEntry do_queue;
	AnyEntry queuename;
	AnyEntry qsub;
	SliderEntry min_dedicated;
	FileNameEntry qsubscript;
	bool have_extra1, have_extra2;
	AnyEntry qsub_extra1;
	AnyEntry qsub_extra2;
	AnyEntry other_args;

    // Run buttons
    Fl_Button *run_button, *print_CL_button, *cite_button;

	// For the run tab
	bool has_mpi;
	bool has_thread;

	// General pipeline I/O (the same for all jobtypes)
	std::vector<Node> pipelineInputNodes;
	std::vector<Node> pipelineOutputNodes;
	std::string pipelineOutputName;


public:
	// Constructor with x, y, w, h and a title
	RelionJobWindow(int nr_tabs, bool _has_mpi, bool _has_thread, bool _has_run = true,
			int x = WCOL0, int y = MENUHEIGHT+10, int w = GUIWIDTH - WCOL0 - 10, int h = GUIHEIGHT_OLD-70, const char* title = "");

    // Destructor
    ~RelionJobWindow() {};

    void resetHeight();

    void setupRunTab();

    // by default nothing happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// write/read settings to disc
	void openWriteFile(std::string fn, std::ofstream &fh);
	bool openReadFile(std::string fn, std::ifstream &fh);
	void closeWriteFile(std::ofstream& fh);
	void closeReadFile(std::ifstream& fh);

	// Write the job submission script
	void saveJobSubmissionScript(std::string newfilename, std::string outputname, std::vector<std::string> commands);

	// See if outputname contains DATENTIMENRUN, and if so, replace with current date and time
	void changeDateNTimeInOutputname(std::string &outputname);

	// Prepare the final (job submission or combined (mpi) command of possibly multiple lines)
	void prepareFinalCommand(std::string &outputname, std::vector<std::string> &commands, std::string &final_command);

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

class GeneralJobWindow : public RelionJobWindow
{
public:

	// I/O
	SliderEntry angpix, particle_diameter;

public:

	// Constructor
	GeneralJobWindow();

	// Destructor
	~GeneralJobWindow(){};

	// write/read settings to disc
	void write(std::string fn);
	void read(std::string fn, bool &_is_continue);

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	void getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command);

};

class CtffindJobWindow : public RelionJobWindow
{
public:

	FileNameEntry mic_names;
	AnyEntry output_star_ctf_mics;
	FileNameEntry fn_ctffind_exe;
	SliderEntry ctf_win;
	SliderEntry cs, kv, q0, angpix, dstep, dast;
	SliderEntry box, resmin, resmax, dfmin, dfmax, dfstep;
	BooleanEntry use_gctf, do_ignore_ctffind_params, do_EPA;
	FileNameEntry fn_gctf_exe;

	Fl_Group *gctf_group;

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

	void getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, RFLOAT angpix);

};

class AutopickJobWindow : public RelionJobWindow
{
public:

	FileNameEntry fn_input_autopick;
	FileNameEntry fn_refs_autopick;
	AnyEntry autopick_rootname;
	BooleanEntry do_invert_refs;
	BooleanEntry do_ctf_autopick;
	BooleanEntry do_ignore_first_ctfpeak_autopick;
	SliderEntry lowpass_autopick;
	SliderEntry psi_sampling_autopick;
	BooleanEntry do_write_fom_maps, do_read_fom_maps;
	SliderEntry threshold_autopick;
	SliderEntry mindist_autopick;
	SliderEntry maxstddevnoise_autopick;

	Fl_Group *autopick_ctf_group;

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

	void getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, RFLOAT angpix, RFLOAT particle_diameter);

};


class ManualpickJobWindow : public RelionJobWindow
{
public:

	FileNameEntry fn_in;
	FileNameEntry fn_out;
	AnyEntry manualpick_rootname;
	SliderEntry lowpass;
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

	void getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, RFLOAT angpix, RFLOAT particle_diameter);

};


class ExtractJobWindow : public RelionJobWindow
{
public:

	// I/O
	FileNameEntry star_mics;
	AnyEntry pick_suffix;
	AnyEntry extract_rootname;

	// extract
	BooleanEntry do_extract;
	SliderEntry extract_size;
	BooleanEntry do_rescale;
	SliderEntry rescale;
	BooleanEntry do_norm;
	SliderEntry white_dust;
	SliderEntry black_dust;
	BooleanEntry do_invert;

	// movies
	BooleanEntry do_movie_extract;
	AnyEntry movie_rootname;
	SliderEntry first_movie_frame;
	SliderEntry last_movie_frame;

	Fl_Group *extract_group, *rescale_group, *norm_group, *movie_extract_group;

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
	void getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, RFLOAT angpix, RFLOAT particle_diameter);

};

class SortJobWindow : public RelionJobWindow
{
public:

	// I/O
	InputNodeEntry input_star;
	AnyEntry fn_out;
	InputNodeEntry fn_refs;
	BooleanEntry do_ctf;
	BooleanEntry do_ignore_first_ctfpeak;

	Fl_Group *ctf_group;

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
	void getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, RFLOAT angpix, RFLOAT particle_diameter);

};



class Class2DJobWindow : public RelionJobWindow
{
public:

	// I/O
	AnyEntry fn_out;
	FileNameEntry fn_cont;
	InputNodeEntry fn_img;
	BooleanEntry do_parallel_discio;

	// CTF
	BooleanEntry do_ctf_correction;
	BooleanEntry ctf_phase_flipped;
	BooleanEntry ctf_intact_first_peak;

	// Optimisation
	SliderEntry nr_iter;
	SliderEntry tau_fudge;
	BooleanEntry do_zero_mask;
	SliderEntry highres_limit;

	// Sampling
	SliderEntry nr_classes;
	BooleanEntry dont_skip_align;
	SliderEntry psi_sampling;
	SliderEntry offset_range;
	SliderEntry offset_step;

	Fl_Group *ctf_group, *dont_skip_align_group;

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
	void getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, RFLOAT angpix, RFLOAT particle_diameter);

};


class Class3DJobWindow : public RelionJobWindow
{
public:

	// I/O
	AnyEntry fn_out;
	FileNameEntry fn_cont;
	InputNodeEntry fn_img;
	InputNodeEntry fn_ref;
	InputNodeEntry fn_mask;
	SliderEntry nr_classes;
	BooleanEntry do_parallel_discio;

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
	BooleanEntry do_zero_mask;
	SliderEntry highres_limit;

	// Sampling
	BooleanEntry dont_skip_align;
	RadioEntry sampling;
	SliderEntry offset_range;
	SliderEntry offset_step;
	BooleanEntry do_local_ang_searches;
	SliderEntry sigma_angles;

	Fl_Group *ctf_group, *dont_skip_align_group, *localsearch_group;

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
	void getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, RFLOAT angpix, RFLOAT particle_diameter);

};


class Auto3DJobWindow : public RelionJobWindow
{
public:

	// I/O
	AnyEntry fn_out;
	FileNameEntry fn_cont;
	InputNodeEntry fn_img;
	InputNodeEntry fn_ref;
	InputNodeEntry fn_mask;
	BooleanEntry do_parallel_discio;

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
	BooleanEntry do_zero_mask;

	// Sampling
	textOnlyEntry autosample_text;
	RadioEntry sampling;
	SliderEntry offset_range;
	SliderEntry offset_step;
	RadioEntry auto_local_sampling;

	// Movies
	BooleanEntry do_movies;
	FileNameEntry fn_movie_star;
	SliderEntry movie_runavg_window;
	SliderEntry movie_sigma_offset;
	BooleanEntry do_alsorot_movies;
	SliderEntry movie_sigma_angles;

	Fl_Group *ctf_group, *movie_group, *alsorot_movie_group;

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
	void getCommands(std::string &outputname, std::vector<std::string> &commands,
			std::string &final_command, RFLOAT angpix, RFLOAT particle_diameter);

};

class PostJobWindow : public RelionJobWindow
{
public:

	// I/O
	InputNodeEntry fn_in;
	AnyEntry fn_out;

	// TODO: take away auto-masking from postprocessing and do this outside the postprocessing in a new jobtype

	// Masking
	BooleanEntry do_automask;
	SliderEntry inimask_threshold;
	SliderEntry extend_inimask;
	SliderEntry width_mask_edge;
	BooleanEntry do_usermask;
	FileNameEntry fn_mask;

	// Sharpening
	BooleanEntry do_auto_bfac;
	SliderEntry autob_lowres;
	BooleanEntry do_adhoc_bfac;
	SliderEntry adhoc_bfac;
	FileNameEntry fn_mtf;

	// Filtering
	BooleanEntry do_skip_fsc_weighting;
	SliderEntry low_pass;

	Fl_Group *automask_group, *usermask_group, *autobfac_group, *adhocbfac_group, *skipweight_group;

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
	void getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command,
			RFLOAT angpix);

};

class PolishJobWindow : public RelionJobWindow
{

public:
	// I/O
	InputNodeEntry fn_in;
	InputNodeEntry fn_mask;
	AnyEntry fn_out;

	// Movements
	SliderEntry movie_runavg_window;
	BooleanEntry do_fit_movement;
	SliderEntry sigma_nb;

	// Damage-model
	BooleanEntry do_bfactor_weighting;
	SliderEntry perframe_highres;
	SliderEntry perframe_bfac_lowres;
	AnyEntry sym_name;

	Fl_Group *fit_group, *weight_group;


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
	void getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command,
			RFLOAT angpix, RFLOAT particle_diameter, RFLOAT black_dust, RFLOAT white_dust);

};

class ResmapJobWindow : public RelionJobWindow
{
public:

	// I/O
	InputNodeEntry fn_in;
	InputNodeEntry fn_mask;

	// Params
	FileNameEntry fn_resmap;
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
	void getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command,
			RFLOAT angpix);

};


class PublishJobWindow : public RelionJobWindow
{
public:

	// I/O
	textOnlyEntry cite_text;
	textOnlyEntry cite_external_text;

public:

	// Constructor
	PublishJobWindow();

	// Destructor
	~PublishJobWindow(){};

	// what happens if you change continue old run radiobutton
	void toggle_new_continue(bool is_continue);

	// Generate the correct commands
	void getCommands(std::string &outputname, std::vector<std::string> &commands, std::string &final_command);

};


#endif /* GUI_JOBWINDOW_H_ */
