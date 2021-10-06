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

#ifndef SRC_GUI_MAINWINDOW_H_
#define SRC_GUI_MAINWINDOW_H_
#define Complex tmpComplex
#include <FL/Fl_Scroll.H>
#include "src/gui_jobwindow.h"
#undef Complex
#include "src/pipeliner.h"


#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <sys/time.h>
#include <algorithm>
#include <iostream>
#include <vector>
// Sizing
#define JOBCOLWIDTH (250)
#define XJOBCOL1 (10)
#define XJOBCOL2 (JOBCOLWIDTH + 25)
#define XJOBCOL3 (2*JOBCOLWIDTH + 40)
#define JOBHEIGHT (170)
#define JOBHALFHEIGHT ( (JOBHEIGHT) / (2) )
#define STDOUT_Y (60)
#define STDERR_Y (170)

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
#define RLN_FONTSIZE 13

// Maximum number of jobs in the job-browsers in the pipeline-part of the GUI
#define MAX_JOBS_BROWSER 50


// This class organises the main window of the relion GUI
static Fl_Hold_Browser *browser;
static Fl_Group *browse_grp[MAX_JOBS_BROWSER];
static Fl_Group *background_grp;
static Fl_Group *pipeliner_jobs_grp;
static Fl_Group *pipeliner_grp;
static Fl_Group *expand_stdout_grp;
static Fl_Choice *display_io_node;
static Fl_Select_Browser *finished_job_browser, *running_job_browser, *scheduled_job_browser, *input_job_browser, *output_job_browser;
static Fl_Box *image_box;
static Fl_Pixmap *xpm_image;
// For keeping track of which process to use in the process browser on the GUI
static std::vector<long int> running_processes, finished_processes, scheduled_processes, input_processes, output_processes, io_nodes;
static bool is_main_continue;
static bool do_overwrite_continue;

static JobWindow *gui_jobwindows[MAX_JOBS_BROWSER];

// Run button
// Sjors 16feb2018: somehow suddenly this run_button needs to be a non-static: otherwise it doesn't change to 'continue now' and doesnt grey out...
static Fl_Button *run_button;
static Fl_Button *print_CL_button;
static Fl_Button *schedule_button;
static Fl_Button *expand_stdout_button;
static Fl_Input *alias_current_job;

static Fl_Text_Buffer *textbuff_stdout;
static Fl_Text_Buffer *textbuff_stderr;

static void Gui_Timer_CB(void *userdata);

// Read-only GUI?
static bool maingui_do_read_only;
// Show expand stdout view
extern bool show_expand_stdout;
// Use ccpem-pipeliner?
static bool use_ccpem_pipeliner;

// The pipeline this GUI is acting on
static PipeLine pipeline;

// Which is the current job being displayed?
static int current_job;
static FileName global_outputname;

// Order jobs in finished window alphabetically?
static bool do_order_alphabetically;

// The last time something changed
static time_t time_last_change;

class SchedulerWindow : public Fl_Window
{
public:

       FileName pipeline_name; // Name of this pipeline (e.g. default)
       std::vector<Fl_Check_Button*> check_buttons;
       Fl_Input *repeat, *wait_before, *wait, *schedule_name, *wait_after;
       std::vector<FileName> my_jobs; // Which jobs to execute

       SchedulerWindow(int w, int h, const char* title): Fl_Window(w, h, title){}

       ~SchedulerWindow() {};

       int fill(FileName _pipeline_name, std::vector<FileName> _scheduled_jobs);

private:

       static void cb_execute(Fl_Widget*, void*);
       inline void cb_execute_i();

       static void cb_cancel(Fl_Widget*, void*);
       inline void cb_cancel_i();


};

// Stdout and stderr display
class StdOutDisplay : public Fl_Text_Display
{
public:
	std::string fn_file;
	StdOutDisplay(int X, int Y, int W, int H, const char *l = 0) : Fl_Text_Display(X, Y, W, H, l){};
	~StdOutDisplay() {};
	int handle(int ev);
};

static StdOutDisplay *disp_stdout, *disp_expand_stdout;
static StdOutDisplay *disp_stderr, *disp_expand_stderr;

class NoteEditorWindow : public Fl_Window
{

public:

	FileName fn_note;
	Fl_Text_Editor *editor;
	Fl_Text_Buffer *textbuff_note;
	bool allow_save;
	NoteEditorWindow(int w, int h, const char* t, FileName _fn_note, bool _allow_save = true);

	~NoteEditorWindow() {};

private:

    static void cb_save(Fl_Widget*, void*);
    inline void cb_save_i();

    static void cb_cancel(Fl_Widget*, void*);
    inline void cb_cancel_i();

};

class GuiMainWindow : public Fl_Window
{

public:

	// For Tabs
	Fl_Menu_Bar *menubar, *menubar2;

	// For clicking in stdout/err windows
	StdOutDisplay *stdoutbox, *stderrbox;

	// Update GUI every how many seconds
	int update_every_sec;

	// Exit GUI after how many seconds idle?
	float exit_after_sec;

	// Numbr of browse tabs
	int nr_browse_tabs;

	// For job submission
    std::string final_command;
    std::vector<std::string> commands;

    // Constructor with w x h size of the window and a title
	GuiMainWindow(int w, int h, const char* title, FileName fn_pipe,
			int _update_every_sec, int _exit_after_sec, bool _do_read_only = false,
			bool _do_tomo = false, bool _use_ccpem_pipeliner = false, bool _do_projdir = false);

    // Destructor
    ~GuiMainWindow(){ clear(); };

    // Clear stuff
    void clear();

    // How will jobs be displayed in the GUI job running, finished, in, out & scheduled job lists
    std::string getJobNameForDisplay(Process &job);

    // Update the content of the finished, running and scheduled job lists
    void fillRunningJobLists();

    // Update the content of the input and output job lists for the current job
    void fillToAndFromJobLists();

    // Need public access for auto-updating the GUI
    void fillStdOutAndErr();

    // Touch the TimeStamp of the last change
    void tickTimeLastChanged();

    // Update all job lists (running, scheduled, finished, as well as to/from)
    void updateJobLists();

    // When a job is selected from the job browsers at the bottom: set current_job there, load that one in the current window
    // and update all job lists at the bottom
    void loadJobFromPipeline(int this_job);

    // Run scheduled jobs from the pipeliner
    void runScheduledJobs(FileName fn_sched, FileName fn_jobids, int nr_repeat, long int minutes_wait);

private:


    // Vertical distance from the top
    int start_y;

    // Current height
    int current_y;


    /** Call-back functions
     *  The method of using two functions of static void and inline void was copied from:
     *  http://www3.telus.net/public/robark/
     */

    static void cb_select_browsegroup(Fl_Widget*, void*);
    inline void cb_select_browsegroup_i(bool is_initial = false);

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

    static void cb_abort(Fl_Widget*, void*);
    inline void cb_abort_i(std::string newalias = "");

    static void cb_mark_as_finished(Fl_Widget*, void*);
    static void cb_mark_as_failed(Fl_Widget*, void*);
    inline void cb_mark_as_finished_i(bool is_failed = false);

    static void cb_edit_project_note(Fl_Widget*, void*);
    static void cb_edit_note(Fl_Widget*, void*);
    inline void cb_edit_note_i(bool is_project_note = false);

    static void cb_print_cl(Fl_Widget*, void*);
    inline void cb_print_cl_i();

    static void cb_save(Fl_Widget*, void*);
    inline void cb_save_i();

    static void cb_load(Fl_Widget*, void*);
    inline void cb_load_i();

    static void cb_undelete_job(Fl_Widget*, void*);
    inline void cb_undelete_job_i();

    static void cb_order_jobs_alphabetically(Fl_Widget*, void*);
    static void cb_order_jobs_chronologically(Fl_Widget*, void*);

    static void cb_empty_trash(Fl_Widget*, void*);
    inline void cb_empty_trash_i();

    static void cb_print_notes(Fl_Widget*, void*);
    inline void cb_print_notes_i();

    static void cb_remake_nodesdir(Fl_Widget*, void*);
    inline void cb_remake_nodesdir_i();

    static void cb_reread_pipeline(Fl_Widget*, void*);
    inline void cb_reread_pipeline_i();

    static void cb_reactivate_runbutton(Fl_Widget*, void*);
    inline void cb_reactivate_runbutton_i();

    static void cb_toggle_overwrite_continue(Fl_Widget*, void*);
    inline void cb_toggle_overwrite_continue_i();

    static void cb_show_initial_screen(Fl_Widget*, void*);
    inline void cb_show_initial_screen_i();

    static void cb_start_pipeliner(Fl_Widget*, void*);
    inline void cb_start_pipeliner_i();

    static void cb_stop_pipeliner(Fl_Widget*, void*);
    inline void cb_stop_pipeliner_i();

    static void cb_toggle_expand_stdout(Fl_Widget*, void*);
    inline void cb_toggle_expand_stdout_i();

    static void cb_about(Fl_Widget*, void*);
    inline void cb_about_i();

public:
    static void cb_quit(Fl_Widget*, void*);
private:
    inline void cb_quit_i();
};




#endif /* SRC_NEWGUI_MAINWINDOW_CPP_ */
