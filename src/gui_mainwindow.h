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

// Maximum number of jobs in the job-browsers in the pipeline-part of the GUI
#define MAX_JOBS_BROWSER 50

// This class organises the main winfow of the relion GUI
static Fl_Choice *add_new_job, *display_io_node;
static Fl_Text_Display *text_current_job;
static Fl_Text_Buffer *textbuff_current_job;
static Fl_Select_Browser *finished_job_browser, *running_job_browser, *scheduled_job_browser, *input_job_browser, *output_job_browser;
// For keeping track of which process to use in the process browser on the GUI
static std::vector<long int> running_processes, finished_processes, scheduled_processes, input_processes, output_processes, io_nodes;
static Fl_Group        *browse_grp[NR_BROWSE_TABS];
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
static ClassSelectJobWindow *job_classselect;
static MaskCreateJobWindow *job_maskcreate;
static SubtractJobWindow *job_subtract;
static PostJobWindow *job_post;
static PolishJobWindow *job_polish;
static ResmapJobWindow *job_resmap;
static PublishJobWindow *job_publish;
// Run button
static Fl_Button *run_button;
static FileName fn_settings;

// A manualpicker jobwindow for display of micrographs....
static ManualpickJobWindow global_manualpickjob;

// Store all the history
static PipeLine pipeline;
// Which is the current job being displayed?
static int current_job;
FileName global_outputname;

static Fl_Menu_Item new_job_options[] = {
		{"Import"},
		{"Motion correction"},
		{"CTF estimation"},
		{"Manual picking"},
		{"Auto-picking"},
		{"Particle extraction"},
		{"Particle sorting"},
		{"2D classification"},
		{"3D classification"},
		{"3D auto-refine"},
		{"Particle polishing"},
		{"Class selection"},
		{"Mask creation"},
		{"Image subtraction"},
		{"Post-processing"},
		{"Local-resolution"},
		{0} // this should be the last entry
};

class RelionMainWindow : public Fl_Window
{

public:

	// For Tabs
	Fl_Menu_Bar *menubar, *menubar2;
	Fl_Tabs *tabs;
	Fl_Group *tab0, *tab1, *tab2, *tab3, *tab4, *tab5;

    // Run button
    Fl_Button *print_CL_button, *cite_button;
    Fl_Button *schedule_button;

    // For job submission
    std::string final_command;
    std::vector<std::string> commands;

    // Constructor with w x h size of the window and a title
	RelionMainWindow(int w, int h, const char* title, FileName fn_pipe);

    // Destructor
    ~RelionMainWindow(){};

    // Communicate with the different jobtype objects
    void jobCommunicate(bool do_write, bool do_read, bool do_toggle_continue, bool do_commandline, bool do_makedir, int this_job = 0);

    // Add a process to the PipeLine
    void addToPipeLine(int as_status, bool do_overwrite = false, int this_job = 0);

    // Update the content of the finished, running and scheduled job lists
    void fillRunningJobLists();

    // Update the content of the input and output job lists for the current job
    void fillToAndFromJobLists();

    // Update all job lists (running, scheduled, finished, as well as to/from)
    void updateJobLists();

    // When a job is selected from the job browsers at the bottom: set current_job there, load that one in the current window
    // and update all job lists at the bottom
    void loadJobFromPipeline();


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
    inline void cb_run_i();

    static void cb_schedule(Fl_Widget*, void*);
    inline void cb_schedule_i();

    static void cb_run_scheduled(Fl_Widget*, void*);
    inline void cb_run_scheduled_i();

    static void cb_delete(Fl_Widget*, void*);
    inline void cb_delete_i(bool do_ask = true, bool do_recursive = true);

    static void cb_cleanup(Fl_Widget*, void*);
    inline void cb_cleanup_i();

    static void cb_mark_as_finished(Fl_Widget*, void*);
    inline void cb_mark_as_finished_i();

    static void cb_print_cl(Fl_Widget*, void*);
    inline void cb_print_cl_i();

    static void cb_menubar_save(Fl_Widget*, void*);
    inline void cb_menubar_save_i();

    static void cb_menubar_reactivate_runbutton(Fl_Widget*, void*);
    inline void cb_menubar_reactivate_runbutton_i();

    static void cb_menubar_about(Fl_Widget*, void*);
    inline void cb_menubar_about_i();

    static void cb_menubar_quit(Fl_Widget*, void*);
    inline void cb_menubar_quit_i();
};

#endif /* GUI_MAINWINDOW_H_ */
