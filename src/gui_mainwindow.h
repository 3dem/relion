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

#ifndef GUI_MAINWINDOW_H_
#define GUI_MAINWINDOW_H_
#include "src/gui_jobwindow.h"
#include "src/gui_entries.h"
#include "src/pipeliner.h"
#include <signal.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/time.h>
#include <algorithm>
#include <iostream>
#include <vector>
#define DO_WRITE true
#define DONT_WRITE false
#define DO_READ true
#define DONT_READ false
#define DO_TOGGLE_CONT true
#define DONT_TOGGLE_CONT false
#define DO_GET_CL true
#define DONT_GET_CL false
#define DO_MKDIR true
#define DONT_MKDIR false
// font size of browser windows on the main GUI
#define RLN_FONTSIZE 14

// Maximum number of jobs in the job-browsers in the pipeline-part of the GUI
#define MAX_JOBS_BROWSER 50

// This class organises the main winfow of the relion GUI
static Fl_Hold_Browser *browser;
static Fl_Group *browse_grp[NR_BROWSE_TABS];
static int browse_jobtype[NR_BROWSE_TABS]; // this allow non-consecutive numbering of jobtypes in the job browser
static Fl_Choice *display_io_node;
static Fl_Select_Browser *finished_job_browser, *running_job_browser, *scheduled_job_browser, *input_job_browser, *output_job_browser;
static Fl_Box *image_box;
static Fl_JPEG_Image *jpeg_image;
// For keeping track of which process to use in the process browser on the GUI
static std::vector<long int> running_processes, finished_processes, scheduled_processes, input_processes, output_processes, io_nodes;
static bool is_main_continue;
static ImportJobWindow *job_import;
static MotioncorrJobWindow *job_motioncorr;
static CtffindJobWindow *job_ctffind;
static ManualpickJobWindow *job_manualpick;
static AutopickJobWindow *job_autopick;
static ExtractJobWindow *job_extract;
static SortJobWindow *job_sort;
static Class2DJobWindow *job_class2d;
static Class3DJobWindow *job_class3d;
static Auto3DJobWindow *job_auto3d;
static MovieRefineJobWindow *job_movierefine;
static ClassSelectJobWindow *job_classselect;
static MaskCreateJobWindow *job_maskcreate;
static JoinStarJobWindow *job_joinstar;
static SubtractJobWindow *job_subtract;
static PostJobWindow *job_post;
static PolishJobWindow *job_polish;
static ResmapJobWindow *job_resmap;
// Run button
static Fl_Button *run_button;
static Fl_Button *print_CL_button;
static Fl_Button *schedule_button;
static Fl_Input *alias_current_job;
// Stdout and stderr display
static Fl_Text_Display *disp_stdout;
static Fl_Text_Display *disp_stderr;
static Fl_Text_Buffer *textbuff_stdout;
static Fl_Text_Buffer *textbuff_stderr;

static FileName fn_settings;
// Initial screen
static bool show_initial_screen;

// A manualpicker jobwindow for display of micrographs....
static ManualpickJobWindow global_manualpickjob;

// Store all the history
static PipeLine pipeline;
// Which is the current job being displayed?
static int current_job;
static FileName global_outputname;

// Order jobs in finished window alphabetically?
static bool do_order_alphabetically;

class NoteEditorWindow : public Fl_Window
{

public:

	FileName fn_note;
	Fl_Text_Editor *editor;
	Fl_Text_Buffer *textbuff_note;
	NoteEditorWindow(int w, int h, const char* t, FileName _fn_note);

	~NoteEditorWindow() {};

private:

    static void cb_save(Fl_Widget*, void*);
    inline void cb_save_i();

    static void cb_cancel(Fl_Widget*, void*);
    inline void cb_cancel_i();

};




class RelionMainWindow : public Fl_Window
{

public:

	// For Tabs
	Fl_Menu_Bar *menubar, *menubar2;
	Fl_Tabs *tabs;
	Fl_Group *tab0, *tab1, *tab2, *tab3, *tab4, *tab5;

    // For job submission
    std::string final_command;
    std::vector<std::string> commands;

    // Constructor with w x h size of the window and a title
	RelionMainWindow(int w, int h, const char* title, FileName fn_pipe);

    // Destructor
    ~RelionMainWindow(){};

    // Communicate with the different jobtype objects
    void jobCommunicate(bool do_write, bool do_read, bool do_toggle_continue, bool do_commandline, bool do_makedir, int this_job = 0);

    // Add a process to the PipeLine, return the number of the process
    long int addToPipeLine(int as_status, bool do_overwrite = false, int this_job = 0);

    // Update the content of the finished, running and scheduled job lists
    void fillRunningJobLists();

    // Update the content of the input and output job lists for the current job
    void fillToAndFromJobLists();

    // Update all job lists (running, scheduled, finished, as well as to/from)
    void updateJobLists();

    // When a job is selected from the job browsers at the bottom: set current_job there, load that one in the current window
    // and update all job lists at the bottom
    void loadJobFromPipeline();

    // Run scheduled jobs from the pipeliner
    void runScheduledJobs(int nr_repeat, long int minutes_wait);

private:


    // Vertical distance from the top
    int start_y;

    // Current height
    int current_y;


    /** Call-back functions for the Run button
     *  The method of using two functions of static void and inline void was copied from:
     *  http://www3.telus.net/public/robark/
     */

    static void cb_select_browsegroup(Fl_Widget*, void*);
    inline void cb_select_browsegroup_i();

    static void cb_select_finished_job(Fl_Widget*, void*);
    inline void cb_select_finished_job_i();

    static void cb_select_running_job(Fl_Widget*, void*);
    inline void cb_select_running_job_i();

    static void cb_select_scheduled_job(Fl_Widget*, void*);
    inline void cb_select_scheduled_job_i();

    static void cb_select_input_job(Fl_Widget*, void*);
    inline void cb_select_input_job_i();

    static void cb_select_output_job(Fl_Widget*, void*);
    inline void cb_select_output_job_i();

    static void cb_display_io_node(Fl_Widget*, void*);
    inline void cb_display_io_node_i();

    static void cb_display(Fl_Widget*, void*);
    inline void cb_display_i();

    inline void cb_toggle_continue_i();

    static void cb_run(Fl_Widget*, void*);
    static void cb_schedule(Fl_Widget*, void*);
    inline void cb_run_i(bool only_schedule = false, bool do_open_edit = true);

    static void cb_delete(Fl_Widget*, void*);
    inline void cb_delete_i(bool do_ask = true, bool do_recursive = true);

    static void cb_gently_clean_all_jobs(Fl_Widget*, void*);
    static void cb_harshly_clean_all_jobs(Fl_Widget*, void*);
    inline void cb_clean_all_jobs_i(bool do_harsh);

    static void cb_gentle_cleanup(Fl_Widget*, void*);
    static void cb_harsh_cleanup(Fl_Widget*, void*);
    inline void cb_cleanup_i(int myjob = -1, bool do_verb = true, bool do_harsh = false);

    static void cb_set_alias(Fl_Widget*, void*);
    inline void cb_set_alias_i(std::string newalias = "");

    static void cb_mark_as_finished(Fl_Widget*, void*);
    inline void cb_mark_as_finished_i();

    static void cb_edit_project_note(Fl_Widget*, void*);
    static void cb_edit_note(Fl_Widget*, void*);
    inline void cb_edit_note_i(bool is_project_note = false);

    inline void cb_fill_stdout_i();

    static void cb_print_cl(Fl_Widget*, void*);
    inline void cb_print_cl_i();

    static void cb_save(Fl_Widget*, void*);
    inline void cb_save_i();

    static void cb_load(Fl_Widget*, void*);
    inline void cb_load_i();

    static void cb_import(Fl_Widget*, void*);
    static void cb_undelete_job(Fl_Widget*, void*);
    inline void cb_import_i(bool is_undelete);

    static void cb_order_jobs_alphabetically(Fl_Widget*, void*);
    static void cb_order_jobs_chronologically(Fl_Widget*, void*);

    static void cb_empty_trash(Fl_Widget*, void*);
    inline void cb_empty_trash_i();

    static void cb_print_notes(Fl_Widget*, void*);
    inline void cb_print_notes_i();

    static void cb_reread_pipeline(Fl_Widget*, void*);
    inline void cb_reread_pipeline_i();

    static void cb_reactivate_runbutton(Fl_Widget*, void*);
    inline void cb_reactivate_runbutton_i();

    static void cb_show_initial_screen(Fl_Widget*, void*);
    inline void cb_show_initial_screen_i();

    static void cb_start_pipeliner(Fl_Widget*, void*);
    inline void cb_start_pipeliner_i();

    static void cb_stop_pipeliner(Fl_Widget*, void*);
    inline void cb_stop_pipeliner_i();

    static void cb_about(Fl_Widget*, void*);
    inline void cb_about_i();

public:
    static void cb_quit(Fl_Widget*, void*);
private:
    inline void cb_quit_i();
};

#endif /* GUI_MAINWINDOW_H_ */
