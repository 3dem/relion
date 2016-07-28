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
#include "src/gui_mainwindow.h"

#define DEBUG

// The StdOutDisplay allows looking at the entire stdout or stderr file
int StdOutDisplay::handle(int ev)
{

	if (ev==FL_PUSH && Fl::event_clicks())
	{
		// double-click
		if (Fl::event_clicks())
		{
			FileName fn = current_browse_directory + fn_file;
			std::string command;
			if (exists(fn))
			{
				if (fn_file == "run.out")
				{
					std::string command = "awk -F\"\r\" '{if (NF>1) {print $NF} else {print}}' < " + fn + " > .gui_tmpstd";
					int res = system(command.c_str());
					NoteEditorWindow* w = new NoteEditorWindow(800, 400, fn.c_str(), ".gui_tmpstd", false); // false means dont_allow_save
					w->show();
					return 1;
				}
				else
				{
					NoteEditorWindow* w = new NoteEditorWindow(800, 400, fn.c_str(), fn, false); // false means dont_allow_save
					w->show();
					return 1;
				}

			}
		} // end if double click
	} // end if FL_PUSH

	return 0;
}

int SchedulerWindow::fill(FileName _pipeline_name, std::vector<FileName> _scheduled_jobs, std::vector<long int> _scheduled_job_ids)
{
	//color(GUI_BACKGROUND_COLOR);
    int current_y = 2, max_y = 2;
    int ystep = 35;

    int xcol = w()-120;

    // Scroll bars
    Fl_Scroll scroll(0, current_y, w(), h());
    scroll.type(Fl_Scroll::VERTICAL);

    my_jobs.clear();
    pipeline_name = _pipeline_name;
    for (int ijob = 0; ijob < _scheduled_jobs.size(); ijob++)
    {
    	my_jobs.push_back(_scheduled_job_ids[ijob]);
    	int xcoor = (ijob < 1+_scheduled_jobs.size()/2) ? 20 : w()-170;
    	if (ijob == 1+_scheduled_jobs.size()/2)
    		current_y = 2;
    	Fl_Check_Button *mycheck = new Fl_Check_Button(xcoor, current_y, ystep-8, ystep-8, _scheduled_jobs[ijob].c_str());
    	mycheck->labelsize(ENTRY_FONTSIZE);
    	check_buttons.push_back(mycheck);
		mycheck->value(1);
		current_y += ystep;
		if (current_y > max_y)
			max_y = current_y;
    }
    current_y = max_y;
    schedule_name = new Fl_Input(xcol, current_y, 100, ystep-8, "Provide a name for this schedule: ");
    current_y += ystep;
    repeat = new Fl_Input(xcol, current_y, 100, ystep-8, "Run the jobs how many times?");
    current_y += ystep;
    wait = new Fl_Input(xcol, current_y, 100, ystep-8, "Wait at least in between (in minutes)?");
    current_y += ystep;

    // Set the input value
    schedule_name->value("schedule1");
	schedule_name->color(GUI_INPUT_COLOR);
	schedule_name->textsize(ENTRY_FONTSIZE);
	schedule_name->labelsize(ENTRY_FONTSIZE);
    repeat->value("1");
	repeat->color(GUI_INPUT_COLOR);
	repeat->textsize(ENTRY_FONTSIZE);
	repeat->labelsize(ENTRY_FONTSIZE);
    wait->value("15");
	wait->color(GUI_INPUT_COLOR);
	wait->textsize(ENTRY_FONTSIZE);
	wait->labelsize(ENTRY_FONTSIZE);

	// Button to execute
	Fl_Button *execute_button = new Fl_Button(w()-200, current_y, 80, 30, "Execute");
	execute_button->color(GUI_RUNBUTTON_COLOR);
	execute_button->labelsize(12);
	execute_button->callback( cb_execute, this);

	// Button to cancel
	Fl_Button *cancel_button = new Fl_Button(w()-100, current_y, 80, 30, "Cancel");
	cancel_button->color(GUI_RUNBUTTON_COLOR);
	cancel_button->labelsize(12);
	cancel_button->callback( cb_cancel, this);

	resizable(*this);
	show();

	return Fl::run();

}

void SchedulerWindow::cb_cancel(Fl_Widget*, void* v)
{
    SchedulerWindow* T=(SchedulerWindow*)v;
    T->hide();
}

void SchedulerWindow::cb_execute(Fl_Widget*, void* v)
{
    SchedulerWindow* T=(SchedulerWindow*)v;
    T->cb_execute_i();
    T->hide();
}

void SchedulerWindow::cb_execute_i()
{
	FileName fn_sched(schedule_name->value());
	FileName fn_check = "RUNNING_PIPELINER_" + pipeline_name + "_" + fn_sched;
	if (exists(fn_check))
	{
		std::string msg =  "ERROR: a file called " + fn_check + " already exists. \n This implies another set of scheduled jobs with this name is already running. \n Cancelling job execution...";
		fl_message(msg.c_str());
	}
	else
	{
		// Make a string with all job-ids to process
		std::string jobids="\"";
		for (int ijob = 0; ijob < my_jobs.size(); ijob++)
		{
			if (check_buttons[ijob]->value())
				jobids += integerToString(my_jobs[ijob]) + " ";
		}
		jobids += "\"";
		std::cerr << " jobids= " << jobids << std::endl;

		std::string myrepeat(repeat->value());
		std::string mywait(wait->value());

		std::string command = "relion_pipeliner --pipeline " + pipeline_name;
		command += " --schedule " + fn_sched;
		command += " --repeat " + myrepeat;
		command += " --min_wait " + mywait;
		command += " --jobids " + jobids;
		// Run this in the background, so control returns to the window
		command += " &";
		int res = system(command.c_str());
		std::cout << " Launching: " << command << std::endl;
		std::cout << " Stop execution of this set of scheduled jobs by deleting file: " << fn_check << std::endl;

	}

}


NoteEditorWindow::NoteEditorWindow(int w, int h, const char* title, FileName _fn_note, bool _allow_save):Fl_Window(w,h,title)
{
	allow_save = _allow_save;
	editor = new Fl_Text_Editor(0, 0, w, h-50);
    editor->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS,10);
	textbuff_note = new Fl_Text_Buffer;
	editor->buffer(textbuff_note);
	textbuff_note->transcoding_warning_action=NULL;
	fn_note = _fn_note;
	if (exists(fn_note))
		int err = textbuff_note->loadfile(fn_note.c_str());
	else
		textbuff_note->text("Describe what this job or project is about here...");
	editor->insert_position(editor->buffer()->length());
	editor->show_insert_position();

	if (allow_save)
	{
		// Button to save and exit
		Fl_Button *save_button = new Fl_Button(w-200, h-40, 80, 30, "Save");
		save_button->color(GUI_RUNBUTTON_COLOR);
		save_button->labelsize(12);
		save_button->callback( cb_save, this);
	}

	// Button to exit
	Fl_Button *cancel_button = new Fl_Button(w-100, h-40, 80, 30, "Cancel");
	cancel_button->color(GUI_RUNBUTTON_COLOR);
	cancel_button->labelsize(12);
	cancel_button->callback( cb_cancel, this);


}

void NoteEditorWindow::cb_cancel(Fl_Widget*, void* v)
{
    NoteEditorWindow* T=(NoteEditorWindow*)v;
    T->hide();
}

void NoteEditorWindow::cb_save(Fl_Widget*, void* v)
{
    NoteEditorWindow* T=(NoteEditorWindow*)v;
    T->cb_save_i();
    T->hide();
}

void NoteEditorWindow::cb_save_i()
{
	int err = textbuff_note->savefile(fn_note.c_str());
}



RelionMainWindow::RelionMainWindow(int w, int h, const char* title, FileName fn_pipe):Fl_Window(w,h,title)
{

	show_initial_screen = true;
	do_order_alphabetically = false;

	FileName fn_lock=".gui_projectdir";
	if (!exists(fn_lock))
	{
		std::cout << " Only run the relion GUI from your ProjectDirectory. Do you want to start a new project here [y/n]? ";
		char c;
		std::cin >> c;
		if (c == 'y' || c == 'Y')
			touch(".gui_projectdir");
		else
		{
			std::cout << " Exiting ... " << std::endl;
			exit(0);
		}
	}

	// First setup the old part of the GUI
	h = GUIHEIGHT_OLD;

	// TODO: control file location and use better figure
	background_grp = new Fl_Group(WCOL0-10, 0 ,w-WCOL0, h-55);

    // First look for image in the binary install directory,
    // then in the source tree.
    FileName fn_bg = std::string(INSTALL_LIBRARY_DIR) + std::string("gui_background.xpm");

    if(!exists(fn_bg))
    {
        fn_bg = std::string(SOURCE_DIR) + std::string("gui_background.xpm");
    }

	if (exists(fn_bg))
	{
		// Initial screen picture with some explanation on how to use the GUI
		image_box = new Fl_Box(WCOL0-10, 0 ,w-WCOL0, h-55); // widget that will contain image
		xpm_image = new Fl_XPM_Image(fn_bg.c_str());
		image_box->image(xpm_image); // attach xpm image to box
		forgot_button = new Fl_Button(450, 143, 10, 32, "?");
		forgot_button->color(GUI_BUTTON_COLOR);
		forgot_button->labelsize(12);
		forgot_button->callback( cb_forgot, this);
	 }
	background_grp->end();

	// Read in the pipeline STAR file if it exists
	pipeline.name = fn_pipe;
	if (exists(fn_pipe + "_pipeline.star"))
	{
		pipeline.read(DO_LOCK);
		// With the locking system, each read needs to be followed soon with a write
		pipeline.write(DO_LOCK);
	}
	else
	{
		pipeline.write();
	}

	// Check which jobs have finished
	pipeline.checkProcessCompletion();

    color(GUI_BACKGROUND_COLOR);
    menubar = new Fl_Menu_Bar(-3, 0, WCOL0-7, MENUHEIGHT);
    menubar->add("File/Re-read pipeline",  FL_ALT+'r', cb_reread_pipeline, this);
    menubar->add("File/Edit project note",  FL_ALT+'e', cb_edit_project_note, this);
    menubar->add("File/Print all notes",  FL_ALT+'p', cb_print_notes, this);
    menubar->add("File/Remake .Nodes\\/",  FL_ALT+'n', cb_remake_nodesdir, this);
    menubar->add("File/Display",  FL_ALT+'d', cb_display, this);
    menubar->add("File/_Show initial screen",  FL_ALT+'z', cb_show_initial_screen, this);
    menubar->add("File/_Empty trash",  FL_ALT+'t', cb_empty_trash, this);
    menubar->add("File/About", 0, cb_about, this);
    menubar->add("File/Quit", FL_ALT+'q', cb_quit, this);

    menubar->add("Jobs/Save job settings",  FL_ALT+'s', cb_save, this);
    menubar->add("Jobs/_Load job settings",  FL_ALT+'l', cb_load, this);
    menubar->add("Jobs/Order alphabetically",  FL_ALT+'a', cb_order_jobs_alphabetically, this);
    menubar->add("Jobs/_Order chronologically",  FL_ALT+'c', cb_order_jobs_chronologically, this);
    menubar->add("Jobs/_Undelete job(s)",  FL_ALT+'u', cb_undelete_job, this);
    menubar->add("Jobs/Export scheduled job(s)",  FL_ALT+'x', cb_export_jobs, this);
    menubar->add("Jobs/_Import scheduled job(s)",  FL_ALT+'i', cb_import_jobs, this);
    menubar->add("Jobs/Gently clean all jobs",  FL_ALT+'g', cb_gently_clean_all_jobs, this);
    menubar->add("Jobs/Harshly clean all jobs",  FL_ALT+'h', cb_harshly_clean_all_jobs, this);

    menubar->add("Autorun/Run scheduled jobs", 0, cb_start_pipeliner, this);
    menubar->add("Autorun/Stop running scheduled jobs", 0, cb_stop_pipeliner, this);

    current_y = MENUHEIGHT + 10;

    // Add run buttons on the menubar as well
	print_CL_button = new Fl_Button(GUIWIDTH - 330, h-90, 100, 32, "Print command");
	print_CL_button->color(GUI_RUNBUTTON_COLOR);
	print_CL_button->labelsize(12);
	print_CL_button->callback( cb_print_cl, this);

	schedule_button = new Fl_Button(GUIWIDTH - 220 , h-90, 100, 32, "Schedule");
	schedule_button->color(GUI_RUNBUTTON_COLOR);
	schedule_button->labelfont(FL_ITALIC);
	schedule_button->labelsize(14);
	schedule_button->callback( cb_schedule, this);

	run_button = new Fl_Button(GUIWIDTH - 110 , h-90, 100, 32, "Run now");
	run_button->color(GUI_RUNBUTTON_COLOR);
	run_button->labelfont(FL_ITALIC);
	run_button->labelsize(14);
	run_button->callback( cb_run, this);

	// Fill browser in the right order
	browser = new Fl_Hold_Browser(10,MENUHEIGHT+5,WCOL0-20,h-MENUHEIGHT-60);
    current_job = -1;

    browse_grp[0] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[0] = PROC_IMPORT;
    browser->add("Import");
    job_import = new ImportJobWindow();
    browse_grp[0]->end();

    browse_grp[1] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[1] = PROC_MOTIONCORR;
    browser->add("Motion correction");
	job_motioncorr = new MotioncorrJobWindow();
    browse_grp[1]->end();

    browse_grp[2] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[2] = PROC_CTFFIND;
    browser->add("CTF estimation");
    job_ctffind = new CtffindJobWindow();
    browse_grp[2]->end();

    browse_grp[3] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[3] = PROC_MANUALPICK;
	browser->add("Manual picking");
	job_manualpick = new ManualpickJobWindow();
    browse_grp[3]->end();

    browse_grp[4] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[4] = PROC_AUTOPICK;
	browser->add("Auto-picking");
	job_autopick = new AutopickJobWindow();
    browse_grp[4]->end();

    browse_grp[5] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[5] = PROC_EXTRACT;
	browser->add("Particle extraction");
	job_extract = new ExtractJobWindow();
    browse_grp[5]->end();

    browse_grp[6] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[6] = PROC_SORT;
	browser->add("Particle sorting");
	job_sort = new SortJobWindow();
    browse_grp[6]->end();

    browse_grp[7] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[7] = PROC_CLASSSELECT;
	browser->add("Subset selection");
	job_classselect = new ClassSelectJobWindow();
    browse_grp[7]->end();

    browse_grp[8] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[8] = PROC_2DCLASS;
	browser->add("2D classification");
	job_class2d = new Class2DJobWindow();
    browse_grp[8]->end();

    browse_grp[9] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[9] = PROC_3DCLASS;
	browser->add("3D classification");
	job_class3d = new Class3DJobWindow();
    browse_grp[9]->end();

    browse_grp[10] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[10] = PROC_3DAUTO;
	browser->add("3D auto-refine");
	job_auto3d = new Auto3DJobWindow();
    browse_grp[10]->end();

    browse_grp[11] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[11] = PROC_MOVIEREFINE;
	browser->add("Movie refinement");
	job_movierefine = new MovieRefineJobWindow();
    browse_grp[11]->end();

    browse_grp[12] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[12] = PROC_POLISH;
	browser->add("Particle polishing");
	job_polish = new PolishJobWindow();
    browse_grp[12]->end();

    browse_grp[13] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[13] = PROC_MASKCREATE;
	browser->add("Mask creation");
	job_maskcreate = new MaskCreateJobWindow();
    browse_grp[13]->end();

    browse_grp[14] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[14] = PROC_JOINSTAR;
	browser->add("Join star files");
	job_joinstar = new JoinStarJobWindow();
    browse_grp[14]->end();

    browse_grp[15] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[15] = PROC_SUBTRACT;
	browser->add("Particle subtraction");
	job_subtract = new SubtractJobWindow();
    browse_grp[15]->end();

    browse_grp[16] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[16] = PROC_POST;
	browser->add("Post-processing");
	job_post = new PostJobWindow();
    browse_grp[16]->end();

    browse_grp[17] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browse_jobtype[17] = PROC_RESMAP;
	browser->add("Local resolution");
	job_resmap = new ResmapJobWindow();
    browse_grp[17]->end();

    browser->callback(cb_select_browsegroup);
    browser->textsize(RLN_FONTSIZE);
    browser->end();
    browser->select(1); // just start from the beginning

    // Pipeline part of the GUI

    menubar2 = new Fl_Menu_Bar(XJOBCOL1, GUIHEIGHT_EXT_START, 100, MENUHEIGHT);
    menubar2->color(GUI_BUTTON_COLOR);
    menubar2->add("Job actions/Edit Note", 0, cb_edit_note, this);
    menubar2->add("Job actions/Alias", 0, cb_set_alias, this);
    menubar2->add("Job actions/Mark as finished", 0, cb_mark_as_finished, this);
    menubar2->add("Job actions/Make flowchart", 0, cb_make_flowchart, this);
    menubar2->add("Job actions/Gentle clean", 0, cb_gentle_cleanup, this);
    menubar2->add("Job actions/Harsh clean", 0, cb_harsh_cleanup, this);
    menubar2->add("Job actions/Delete", 0, cb_delete, this);

    // Fl_input with the alias of the new job (or the name of an existing one)
    alias_current_job = new Fl_Input(XJOBCOL2-50 , GUIHEIGHT_EXT_START+3, JOBCOLWIDTH, MENUHEIGHT-6, "Current job:");

    // Left-hand side browsers for input/output nodes and processes
	display_io_node  = new Fl_Choice(XJOBCOL3, GUIHEIGHT_EXT_START+3, 250, MENUHEIGHT-6);
    display_io_node->label("Display:");
    display_io_node->color(GUI_BUTTON_COLOR);
	display_io_node->callback(cb_display_io_node, this);

	// Add browsers for finished, running and scheduled jobs
    Fl_Text_Buffer *textbuff1 = new Fl_Text_Buffer();
    Fl_Text_Buffer *textbuff2 = new Fl_Text_Buffer();
    Fl_Text_Buffer *textbuff3 = new Fl_Text_Buffer();
    Fl_Text_Buffer *textbuff4 = new Fl_Text_Buffer();
    Fl_Text_Buffer *textbuff5 = new Fl_Text_Buffer();
    textbuff1->text("Finished jobs");
    textbuff2->text("Running jobs");
    textbuff3->text("Scheduled jobs");
    textbuff4->text("Input to this job");
    textbuff5->text("Output from this job");
    Fl_Text_Display* textdisp1 = new Fl_Text_Display(XJOBCOL1, GUIHEIGHT_EXT_START2, JOBCOLWIDTH, 25);
	Fl_Text_Display* textdisp2 = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_EXT_START2, JOBCOLWIDTH, 25);
	Fl_Text_Display* textdisp3 = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_EXT_START2 + JOBHALFHEIGHT + 25, JOBCOLWIDTH, 25);
	Fl_Text_Display* textdisp4 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_EXT_START2, JOBCOLWIDTH, 25);
	Fl_Text_Display* textdisp5 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_EXT_START2 + JOBHALFHEIGHT + 25, JOBCOLWIDTH, 25);
	textdisp1->buffer(textbuff1);
	textdisp2->buffer(textbuff2);
	textdisp3->buffer(textbuff3);
	textdisp4->buffer(textbuff4);
	textdisp5->buffer(textbuff5);
	textdisp1->color(GUI_BACKGROUND_COLOR);
	textdisp2->color(GUI_BACKGROUND_COLOR);
	textdisp3->color(GUI_BACKGROUND_COLOR);
	textdisp4->color(GUI_BACKGROUND_COLOR);
	textdisp5->color(GUI_BACKGROUND_COLOR);

    finished_job_browser  = new Fl_Select_Browser(XJOBCOL1, GUIHEIGHT_EXT_START2 + 25, JOBCOLWIDTH, JOBHEIGHT+25);
    running_job_browser   = new Fl_Select_Browser(XJOBCOL2, GUIHEIGHT_EXT_START2 + 25, JOBCOLWIDTH, JOBHALFHEIGHT);
    scheduled_job_browser = new Fl_Select_Browser(XJOBCOL2, GUIHEIGHT_EXT_START2 + 25 + JOBHALFHEIGHT + 25, JOBCOLWIDTH, JOBHALFHEIGHT);
    input_job_browser    = new Fl_Select_Browser(XJOBCOL3,  GUIHEIGHT_EXT_START2 + 25, JOBCOLWIDTH, JOBHALFHEIGHT);
    output_job_browser   = new Fl_Select_Browser(XJOBCOL3,  GUIHEIGHT_EXT_START2 + 25 + JOBHALFHEIGHT + 25, JOBCOLWIDTH, JOBHALFHEIGHT);

    // Fill the actual browsers
    fillRunningJobLists();

    // Set the callbacks
    finished_job_browser->callback(cb_select_finished_job);
    running_job_browser->callback(cb_select_running_job);
    scheduled_job_browser->callback(cb_select_scheduled_job);
    input_job_browser->callback(cb_select_input_job);
    output_job_browser->callback(cb_select_output_job);
    finished_job_browser->textsize(RLN_FONTSIZE);
    running_job_browser->textsize(RLN_FONTSIZE);
    scheduled_job_browser->textsize(RLN_FONTSIZE);
    input_job_browser->textsize(RLN_FONTSIZE);
    output_job_browser->textsize(RLN_FONTSIZE);

    finished_job_browser->end();
    running_job_browser->end();
    scheduled_job_browser->end();
    input_job_browser->end();
    output_job_browser->end();

    // Display stdout and stderr of jobs
    textbuff_stdout = new Fl_Text_Buffer();
    textbuff_stderr = new Fl_Text_Buffer();
    // Disable warning message about UTF-8 transcoding
	disp_stdout = new StdOutDisplay(XJOBCOL1, GUIHEIGHT_EXT_START2 + JOBHEIGHT + STDOUT_Y, w-20, 110);
    disp_stderr = new StdOutDisplay(XJOBCOL1, GUIHEIGHT_EXT_START2 + JOBHEIGHT + STDERR_Y, w-20, 60);
    disp_stdout->fn_file = "run.out";
    disp_stderr->fn_file = "run.err";
    textbuff_stdout->text("stdout will go here; double-click this window to open stdout in a separate window");
    textbuff_stderr->text("stderr will go here; double-click this window to open stderr in a separate window");
    disp_stdout->buffer(textbuff_stdout);
    disp_stderr->buffer(textbuff_stderr);
    disp_stderr->textcolor(FL_RED);
    disp_stdout->textsize(RLN_FONTSIZE);
    disp_stderr->textsize(RLN_FONTSIZE);
    disp_stdout->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS,0);
    disp_stderr->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS,0);
    disp_stdout->scrollbar_width(0);
    disp_stderr->scrollbar_width(0);

    // Set and activate current selection from side-browser
    cb_select_browsegroup_i(); // make default active
    is_main_continue = false; // default is a new run

    Fl::add_timeout(2, Timer_CB, (void*)this);

}

static void Timer_CB(void *userdata) {
	RelionMainWindow *o = (RelionMainWindow*)userdata;
    // Update the stdout and stderr windows if we're currently pointing at a running job
    if (current_job >= 0 && pipeline.processList[current_job].status == PROC_RUNNING)
    	o->fillStdOutAndErr();
    // Always check for job completion
    o->updateJobLists();
    // Refresh every 2 seconds
    Fl::repeat_timeout(2, Timer_CB, userdata);
}

// Update the content of the finished, running and scheduled job lists
void RelionMainWindow::fillRunningJobLists()
{
	// Go back to the same positions in the vertical scroll bars of the job lists after updating...
	int mypos_running = running_job_browser->position();
	int mypos_scheduled = scheduled_job_browser->position();
	int mypos_finished = finished_job_browser->position();

    // Clear whatever was in there
	finished_job_browser->clear();
	finished_processes.clear();
	running_job_browser->clear();
	running_processes.clear();
	scheduled_job_browser->clear();
	scheduled_processes.clear();

	// Fill the finished Jobs browsers
	if (do_order_alphabetically)
	{
		// Only re-order the finished jobs!
        std::vector<std::pair<std::string,long int> > vp;
		for (long int i = pipeline.processList.size() -1; i >= 0; i--)
		{
			if (pipeline.processList[i].alias != "None")
				vp.push_back(std::make_pair(pipeline.processList[i].alias, i));
			else
				vp.push_back(std::make_pair(pipeline.processList[i].name, i));
		}
        // Sort on the first elements of the pairs
        std::sort(vp.begin(), vp.end());

		for (long int ip = 0; ip < vp.size(); ip++)
		{
			long int i = vp[ip].second;
			if (pipeline.processList[i].status == PROC_FINISHED)
			{
				finished_processes.push_back(i);
				finished_job_browser->add(vp[ip].first.c_str());
			}
		}
	}
	else
	{
		// For finished jobs search backwards, so that last jobs are at the top
		for (long int i = pipeline.processList.size() -1; i >= 0; i--)
		{
			if (pipeline.processList[i].status == PROC_FINISHED)
			{
				finished_processes.push_back(i);
				if (pipeline.processList[i].alias != "None")
					finished_job_browser->add(pipeline.processList[i].alias.c_str());
				else
					finished_job_browser->add(pipeline.processList[i].name.c_str());
			}
		}
	}

	// For running and scheduled jobs search forwards, so that last jobs are at the bottom
	for (long int i = 0; i < pipeline.processList.size(); i++)
	{
		if (pipeline.processList[i].status == PROC_RUNNING)
		{
			running_processes.push_back(i);
			if (pipeline.processList[i].alias != "None")
				running_job_browser->add(pipeline.processList[i].alias.c_str());
			else
				running_job_browser->add(pipeline.processList[i].name.c_str());
		}
		else if (pipeline.processList[i].status == PROC_SCHEDULED_NEW ||
				 pipeline.processList[i].status == PROC_SCHEDULED_CONT)
		{
			scheduled_processes.push_back(i);
			if (pipeline.processList[i].alias != "None")
				scheduled_job_browser->add(pipeline.processList[i].alias.c_str());
			else
				scheduled_job_browser->add(pipeline.processList[i].name.c_str());
		}
	}

	running_job_browser->position(mypos_running);
	scheduled_job_browser->position(mypos_scheduled);
	finished_job_browser->position(mypos_finished);
}

void RelionMainWindow::fillToAndFromJobLists()
{
	display_io_node->clear();
	input_job_browser->clear();
	output_job_browser->clear();
	io_nodes.clear();
	input_processes.clear();
	output_processes.clear();

	if (current_job >= 0)
	{
		// Where do the input nodes come from?
		for (long int inode = 0; inode < (pipeline.processList[current_job]).inputNodeList.size(); inode++)
		{
			long int mynode = (pipeline.processList[current_job]).inputNodeList[inode];

			if (pipeline.nodeList[mynode].type != NODE_MOVIES) // no display for movie rootname
			{
				FileName fnt = pipeline.nodeList[mynode].name;
				if (exists(fnt))
				{
					fnt = "in: " + fnt.afterLastOf("/");
					display_io_node->add(fnt.c_str());
					io_nodes.push_back(mynode);
				}
			}

			long int myproc = (pipeline.nodeList[mynode]).outputFromProcess;
			if (myproc >= 0)
			{
				// Check if this process was already there
				bool already_there = false;
				for (long int i = 0; i < input_processes.size(); i++)
				{
					if (myproc == input_processes[i])
					{
						already_there=true;
						break;
					}
				}
				if (!already_there)
				{
					input_processes.push_back(myproc);
					if (pipeline.processList[myproc].alias != "None")
						input_job_browser->add(pipeline.processList[myproc].alias.c_str());
					else
						input_job_browser->add(pipeline.processList[myproc].name.c_str());
				}
			}
		}
		// Where do the output nodes lead to?
		for (long int inode = 0; inode < (pipeline.processList[current_job]).outputNodeList.size(); inode++)
		{
			long int mynode = (pipeline.processList[current_job]).outputNodeList[inode];
			FileName fnt = pipeline.nodeList[mynode].name;
			if (exists(fnt))
			{
				fnt = "out: " + fnt.afterLastOf("/");
				display_io_node->add(fnt.c_str());
				io_nodes.push_back(mynode);
			}

			long int nr_outputs = (pipeline.nodeList[mynode]).inputForProcessList.size();
			for (long int iproc = 0; iproc < nr_outputs; iproc++)
			{
				long int myproc =  (pipeline.nodeList[mynode]).inputForProcessList[iproc];
				// Check if this process was already there
				bool already_there = false;
				for (long int i = 0; i < output_processes.size(); i++)
				{
					if (myproc == output_processes[i])
					{
						already_there=true;
						break;
					}
				}
				if (!already_there)
				{
					output_processes.push_back(myproc);
					if (pipeline.processList[myproc].alias != "None")
						output_job_browser->add(pipeline.processList[myproc].alias.c_str());
					else
						output_job_browser->add(pipeline.processList[myproc].name.c_str());
				}
			}
		}
	}
}

void RelionMainWindow::updateJobLists()
{
	pipeline.checkProcessCompletion();
	fillRunningJobLists();
	fillToAndFromJobLists();
}

void RelionMainWindow::loadJobFromPipeline()
{

	int this_job = current_job;
	int itype = pipeline.processList[current_job].type;
	fn_settings = pipeline.processList[current_job].name;

	// The following line allows certain browse buttons to only open the current directory (using CURRENT_ODIR)
	current_browse_directory = pipeline.processList[current_job].name;

	// What type of job is this?
	for ( int t=0; t<NR_BROWSE_TABS; t++ )
	{
		if ( browse_jobtype[t] == itype )
			browser->value(t+1);
	}

	// the new job browser has reset current_job to -1....
	current_job = this_job;

	// Re-read the settings for this job
	cb_select_browsegroup_i(); // change to the corresponding jobwindow
	jobCommunicate(DONT_WRITE, DO_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DONT_MKDIR);

	// Update all job lists in the main GUI
	updateJobLists();

	// If a finished or running job was loaded from the pipeline: set this to be a continuation job
	// If a scheduled job was loaded, only set is_main_continue to true when it is PROC_SCHEDULED_NEW
    if (pipeline.processList[current_job].status != PROC_SCHEDULED_NEW)
    	is_main_continue = true;
    else
    	is_main_continue = false;

    cb_toggle_continue_i();

    if (pipeline.processList[current_job].alias != "None")
    	alias_current_job->value(pipeline.processList[current_job].alias.c_str());
    else
    	alias_current_job->value(pipeline.processList[current_job].name.c_str());

	fillStdOutAndErr();
}

long int RelionMainWindow::addToPipeLine(int as_status, bool do_overwrite, int this_job)
{
	int itype = (this_job > 0) ? this_job : (browse_jobtype[browser->value() - 1]); // browser starts counting at 1 ...
	std::vector<Node> inputnodes;
	std::vector<Node> outputnodes;
	std::string oname;

	switch (itype)
	{
	case PROC_IMPORT:
	{
		inputnodes = job_import->pipelineInputNodes;
		outputnodes= job_import->pipelineOutputNodes;
		oname = job_import->pipelineOutputName;
		break;
	}
	case PROC_MOTIONCORR:
	{
		inputnodes = job_motioncorr->pipelineInputNodes;
		outputnodes= job_motioncorr->pipelineOutputNodes;
		oname = job_motioncorr->pipelineOutputName;
		break;
	}
	case PROC_MANUALPICK:
	{
		inputnodes = job_manualpick->pipelineInputNodes;
		outputnodes= job_manualpick->pipelineOutputNodes;
		oname = job_manualpick->pipelineOutputName;
		break;
	}
	case PROC_CTFFIND:
	{
		inputnodes = job_ctffind->pipelineInputNodes;
		outputnodes= job_ctffind->pipelineOutputNodes;
		oname = job_ctffind->pipelineOutputName;
		break;
	}
	case PROC_AUTOPICK:
	{
		inputnodes = job_autopick->pipelineInputNodes;
		outputnodes= job_autopick->pipelineOutputNodes;
		oname = job_autopick->pipelineOutputName;
		break;
	}
	case PROC_EXTRACT:
	{
		inputnodes = job_extract->pipelineInputNodes;
		outputnodes= job_extract->pipelineOutputNodes;
		oname = job_extract->pipelineOutputName;
		break;
	}
	case PROC_SORT:
	{
		inputnodes = job_sort->pipelineInputNodes;
		outputnodes= job_sort->pipelineOutputNodes;
		oname = job_sort->pipelineOutputName;
		break;
	}
	case PROC_CLASSSELECT:
	{
		inputnodes = job_classselect->pipelineInputNodes;
		outputnodes= job_classselect->pipelineOutputNodes;
		oname = job_classselect->pipelineOutputName;
		break;
	}
	case PROC_2DCLASS:
	{

		inputnodes = job_class2d->pipelineInputNodes;
		outputnodes= job_class2d->pipelineOutputNodes;
		oname = job_class2d->pipelineOutputName;
		break;
	}
	case PROC_3DCLASS:
	{

		inputnodes = job_class3d->pipelineInputNodes;
		outputnodes= job_class3d->pipelineOutputNodes;
		oname = job_class3d->pipelineOutputName;
		break;
	}
	case PROC_3DAUTO:
	{

		inputnodes = job_auto3d->pipelineInputNodes;
		outputnodes= job_auto3d->pipelineOutputNodes;
		oname = job_auto3d->pipelineOutputName;
		break;
	}
	case PROC_MOVIEREFINE:
	{
		inputnodes = job_movierefine->pipelineInputNodes;
		outputnodes= job_movierefine->pipelineOutputNodes;
		oname = job_movierefine->pipelineOutputName;
		break;
	}
	case PROC_POLISH:
	{
		inputnodes = job_polish->pipelineInputNodes;
		outputnodes= job_polish->pipelineOutputNodes;
		oname = job_polish->pipelineOutputName;
		break;
	}
	case PROC_MASKCREATE:
	{
		inputnodes = job_maskcreate->pipelineInputNodes;
		outputnodes= job_maskcreate->pipelineOutputNodes;
		oname = job_maskcreate->pipelineOutputName;
		break;
	}
	case PROC_JOINSTAR:
	{
		inputnodes = job_joinstar->pipelineInputNodes;
		outputnodes= job_joinstar->pipelineOutputNodes;
		oname = job_joinstar->pipelineOutputName;
		break;
	}
	case PROC_SUBTRACT:
	{
		inputnodes = job_subtract->pipelineInputNodes;
		outputnodes= job_subtract->pipelineOutputNodes;
		oname = job_subtract->pipelineOutputName;
		break;
	}
	case PROC_POST:
	{
		inputnodes = job_post->pipelineInputNodes;
		outputnodes= job_post->pipelineOutputNodes;
		oname = job_post->pipelineOutputName;
		break;
	}
	case PROC_RESMAP:
	{

		inputnodes = job_resmap->pipelineInputNodes;
		outputnodes= job_resmap->pipelineOutputNodes;
		oname = job_resmap->pipelineOutputName;
		break;
	}
	default:
	{
		REPORT_ERROR("ERROR: unrecognised job-type to add to the pipeline");
	}
	}

	// Also write a mini-pipeline in the output directory
	PipeLine mini_pipeline;
	mini_pipeline.setName(oname+"job");

	// Add Process to the processList of the pipeline
	Process process(oname, itype, as_status);
	long int myProcess = pipeline.addNewProcess(process, do_overwrite);
	mini_pipeline.addNewProcess(process);

	// Add all input nodes
	for (int i=0; i < inputnodes.size(); i++)
	{
		pipeline.addNewInputEdge(inputnodes[i], myProcess);
		mini_pipeline.addNewInputEdge(inputnodes[i], 0);
	}
	// Add all output nodes
	for (int i=0; i < outputnodes.size(); i++)
	{
		pipeline.addNewOutputEdge(myProcess, outputnodes[i]);
		mini_pipeline.addNewOutputEdge(0, outputnodes[i]);
	}

	// Write the mini-pipeline to an updated STAR file
	mini_pipeline.write();
	// Writing of the overall pipeline is done in the function calling addToPipeLine

	return myProcess;
}

bool RelionMainWindow::jobCommunicate(bool do_write, bool do_read, bool do_toggle_continue, bool do_commandline, bool do_makedir, int this_job)
{
	int itype = (this_job > 0) ? this_job : (browse_jobtype[browser->value() - 1]); // browser starts counting at 1 ....
	show_initial_screen = false;

	bool is_scheduled=false;
	if (current_job >= 0)
		is_scheduled= (pipeline.processList[current_job].status == PROC_SCHEDULED_CONT || pipeline.processList[current_job].status == PROC_SCHEDULED_NEW);

	if (do_commandline)
	{
		// Except for continuation or scheduled jobs, all jobs get a new unique directory
		if (is_main_continue || is_scheduled)
			global_outputname = fn_settings.beforeLastOf("/") + "/";
		else
			global_outputname = "";
	}

	bool result = true;

	switch (itype)
	{
	case PROC_IMPORT:
	{
		if (do_write)
			job_import->write(fn_settings);
		if (do_read)
			job_import->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_import->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_import->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_MOTIONCORR:
	{
		if (do_write)
			job_motioncorr->write(fn_settings);
		if (do_read)
			job_motioncorr->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_motioncorr->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_motioncorr->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_MANUALPICK:
	{
		if (do_write)
			job_manualpick->write(fn_settings);
		if (do_read)
			job_manualpick->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_manualpick->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_manualpick->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_CTFFIND:
	{
		if (do_write)
			job_ctffind->write(fn_settings);
		if (do_read)
			job_ctffind->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_ctffind->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_ctffind->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_AUTOPICK:
	{
		if (do_write)
			job_autopick->write(fn_settings);
		if (do_read)
			job_autopick->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_autopick->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_autopick->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_EXTRACT:
	{
		if (do_write)
			job_extract->write(fn_settings);
		if (do_read)
			job_extract->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_extract->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_extract->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_SORT:
	{
		if (do_write)
			job_sort->write(fn_settings);
		if (do_read)
			job_sort->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_sort->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_sort->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_CLASSSELECT:
	{
		if (do_write)
			job_classselect->write(fn_settings);
		if (do_read)
			job_classselect->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_classselect->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_classselect->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_2DCLASS:
	{
		if (do_write)
			job_class2d->write(fn_settings);
		if (do_read)
			job_class2d->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_class2d->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_class2d->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_3DCLASS:
	{
		if (do_write)
			job_class3d->write(fn_settings);
		if (do_read)
			job_class3d->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_class3d->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_class3d->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_3DAUTO:
	{
		if (do_write)
			job_auto3d->write(fn_settings);
		if (do_read)
			job_auto3d->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_auto3d->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_auto3d->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_MOVIEREFINE:
	{
		if (do_write)
			job_movierefine->write(fn_settings);
		if (do_read)
			job_movierefine->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_movierefine->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_movierefine->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_POLISH:
	{
		if (do_write)
			job_polish->write(fn_settings);
		if (do_read)
			job_polish->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_polish->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_polish->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_MASKCREATE:
	{
		if (do_write)
			job_maskcreate->write(fn_settings);
		if (do_read)
			job_maskcreate->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_maskcreate->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_maskcreate->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_JOINSTAR:
	{
		if (do_write)
			job_joinstar->write(fn_settings);
		if (do_read)
			job_joinstar->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_joinstar->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_joinstar->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_SUBTRACT:
	{
		if (do_write)
			job_subtract->write(fn_settings);
		if (do_read)
			job_subtract->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_subtract->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_subtract->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_POST:
	{
		if (do_write)
			job_post->write(fn_settings);
		if (do_read)
			job_post->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_post->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_post->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	case PROC_RESMAP:
	{
		if (do_write)
			job_resmap->write(fn_settings);
		if (do_read)
			job_resmap->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_resmap->toggle_new_continue(is_main_continue);
		if (do_commandline)
			result = job_resmap->getCommands(global_outputname, commands, final_command, do_makedir, pipeline.job_counter);
		break;
	}
	} // end switch

	// set the continue button correct upon reading of old settings
	if (do_read)
	{
		// Make the choice active
		cb_toggle_continue_i();
	}

	// getCommands may be cancelled
	return result;
}

void RelionMainWindow::runScheduledJobs(FileName fn_sched, FileName fn_jobids, int nr_repeat, long int minutes_wait)
{

	std::vector<long int> my_scheduled_processes;
	std::vector<std::string> jobids;
	int njobs = splitString(fn_jobids, " ", jobids);
	if (njobs == 0)
		my_scheduled_processes = scheduled_processes;
	else
		for (int i = 0; i < njobs; i++)
			my_scheduled_processes.push_back(textToInteger(jobids[i]));

	FileName fn_log = "pipeline_" + pipeline.name + "_" + fn_sched + ".log";
	std::ofstream  fh;
	fh.open((fn_log).c_str(), std::ios::app);

	std::cout << " PIPELINER: writing out information in logfile " << fn_log << std::endl;

	// Touch the fn_check file
	FileName fn_check = "RUNNING_PIPELINER_"+pipeline.name + "_" + fn_sched;
	touch(fn_check);
	bool fn_check_exists = true;

	if (my_scheduled_processes.size() < 1)
		REPORT_ERROR("relion_pipeliner ERROR: there are no scheduled jobs. Exiting...");

	fh << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	fh << " Starting a new scheduler execution called " << fn_sched << std::endl;
	fh << " The scheduled jobs are: " << std::endl;
	for (long int i = 0; i < scheduled_processes.size(); i++)
		fh << " - " << pipeline.processList[scheduled_processes[i]].name << std::endl;
	fh << " Will execute the scheduled jobs " << nr_repeat << " times." << std::endl;
	if (nr_repeat > 1)
		fh << " Will wait until at least " << minutes_wait << " minutes have passed between each repeat." << std::endl;
	fh << " Will be checking for existence of file " << fn_check << "; if it no longer exists, the scheduler will stop." << std::endl;
	fh << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

	int repeat = 0;
	for (repeat = 0 ; repeat < nr_repeat; repeat++)
	{
		fh << " + Starting the " << repeat+1 << "th repeat .." << std::endl;

		// Get starting time of the repeat cycle
		timeval time_start, time_end;
		gettimeofday(&time_start, NULL);

		for (long int i = 0; i < my_scheduled_processes.size(); i++)
		{
			current_job = my_scheduled_processes[i];
			loadJobFromPipeline();
			fh << " + -- Executing " << pipeline.processList[current_job].name << " .. " << std::endl;
			cb_run_i(false, false); //dont only schedule and dont open the editor window

			// Now wait until that job is done!
			while (true)
			{
				if (!exists(fn_check))
				{
					fn_check_exists = false;
					break;
				}

				sleep(10);
				pipeline.checkProcessCompletion();
				if (pipeline.processList[current_job].status == PROC_FINISHED)
				{
					// Read in existing pipeline, in case some other window had changed something else
					pipeline.read(DO_LOCK);
					// Re-set set status of current_job to finished
					pipeline.processList[current_job].status = PROC_FINISHED;
					// Write out the modified pipeline with the status of current_job to finished
					pipeline.write(DO_LOCK);
					break;
				}
			}

			if (!fn_check_exists)
				break;

			if (repeat + 1 != nr_repeat)
			{

				// Read in existing pipeline, in case some other window had changed it
				pipeline.read(DO_LOCK);

				// Set the current job back into the job list of the repeating cycle
				// Do we want to run this as NEW or CONTINUED NEXT TIME?
				int mytype = pipeline.processList[current_job].type;
				// The following jobtypes have functionality to only do the unfinished part of the job
				if (mytype == PROC_MOTIONCORR || mytype == PROC_CTFFIND || mytype == PROC_AUTOPICK || mytype == PROC_EXTRACT || mytype == PROC_MOVIEREFINE)
					pipeline.processList[current_job].status = PROC_SCHEDULED_CONT;
				else
					pipeline.processList[current_job].status = PROC_SCHEDULED_NEW;

				// Write the pipeline to an updated STAR file
				pipeline.write(DO_LOCK);
			}
		}

		if (!fn_check_exists)
			break;

		// Wait at least until 'minutes_wait' minutes have passed from the beginning of the repeat cycle
		gettimeofday(&time_end, NULL);
		long int passed_minutes = (time_end.tv_sec - time_start.tv_sec)/60;
		long int still_wait = minutes_wait - passed_minutes;
		if (still_wait > 0 && repeat+1 != nr_repeat)
		{
			fh << " + -- Waiting " << still_wait << " minutes until next repeat .."<< std::endl;
			sleep(still_wait * 60);
		}

	}

	if (repeat == nr_repeat)
	{
		fh << " + performed all requested repeats in scheduler " << fn_sched << ". Stopping now ..." << std::endl;
		std::cout << " PIPELINER: performed all requested repeats, stopping now ..." << std::endl;
		std::cout << " PIPELINER: you may want to re-read the pipeline (from the File menu) to update the job lists." << std::endl;

		// Read in existing pipeline, in case some other window had changed it
		pipeline.read(DO_LOCK);

		// After breaking out of repeat, set status of the jobs to finished
		for (long int i = 0; i < my_scheduled_processes.size(); i++)
		{
			pipeline.processList[my_scheduled_processes[i]].status = PROC_FINISHED;
		}

		// Write the pipeline to an updated STAR file
		pipeline.write(DO_LOCK);

		// Remove the temporary file
		std::remove(fn_check.c_str());
		exit(0);
	}
	else if (!fn_check_exists)
	{
		fh << " + File " << fn_check << " was removed. Stopping now .." << std::endl;
		std::cout << " PIPELINER: the " << fn_check << " file was removed. Stopping now ..." << std::endl;
		std::cout << " PIPELINER: you may want to re-read the pipeline (from the File menu) to update the job lists." << std::endl;
	}
	else
	{
		REPORT_ERROR("PIPELINER BUG: This shouldn't happen, either fn_check should not exist or we should reach end of repeat cycles...");
	}

	fh << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	fh.close();

	cb_quit_i();
}

void RelionMainWindow::cb_select_browsegroup(Fl_Widget* o, void* v)
{
	RelionMainWindow* T=(RelionMainWindow*)v;
	// When clicking the job browser on the left: reset current_job to -1 (i.e. a new job, not yet in the pipeline)
	current_job = -1;
	T->cb_select_browsegroup_i();
	run_button->activate();

}

void RelionMainWindow::cb_select_browsegroup_i()
{

	// Hide the initial screen
	if (show_initial_screen)
		background_grp->show();
	else
		background_grp->hide();

	// Show the 'selected' group, hide the others
    for ( int t=0; t<NR_BROWSE_TABS; t++ )
    {
    	// During the initial screen: show a nice picture with some explanations
    	if ( t == (browser->value() - 1) && !show_initial_screen ) // browser starts counting at 1...
        {
        	browse_grp[t]->show();
        }
        else
        {
        	browse_grp[t]->hide();
        }
    }

    // Toggle the new tab according to the continue toggle button
    jobCommunicate(DONT_WRITE, DONT_READ, DO_TOGGLE_CONT, DONT_GET_CL, DONT_MKDIR);

	// Update all job lists in the main GUI
	updateJobLists();

	// If an new job was loaded from the job-browser: set the GUI to be a new job
    is_main_continue = false;
    cb_toggle_continue_i(); // make default active

    alias_current_job->value("Give_alias_here");
	//current_job = -1;
    textbuff_stdout->text("stdout will go here; double-click this window to open stdout in a separate window");
    textbuff_stderr->text("stderr will go here; double-click this window to open stderr in a separate window");

}


void RelionMainWindow::cb_select_finished_job(Fl_Widget* o, void* v)
{
	RelionMainWindow* T=(RelionMainWindow*)v;
	T->cb_select_finished_job_i();
	run_button->activate();
}

void RelionMainWindow::cb_select_finished_job_i()
{
	// Show the 'selected' group, hide the others
    int idx = finished_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
    {
    	current_job = finished_processes[idx];
		loadJobFromPipeline();
    }
}

void RelionMainWindow::cb_select_running_job(Fl_Widget* o, void* v)
{
	RelionMainWindow* T=(RelionMainWindow*)v;
	T->cb_select_running_job_i();
	run_button->activate();
}

void RelionMainWindow::cb_select_running_job_i()
{
	// Show the 'selected' group, hide the others
    int idx = running_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
    {
    	current_job = running_processes[idx];
    	loadJobFromPipeline();
    }
}

void RelionMainWindow::cb_select_scheduled_job(Fl_Widget* o, void* v)
{
	RelionMainWindow* T=(RelionMainWindow*)v;
	T->cb_select_scheduled_job_i();
	run_button->activate();
}

void RelionMainWindow::cb_select_scheduled_job_i()
{
	// Show the 'selected' group, hide the others
    int idx = scheduled_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
    {
    	current_job = scheduled_processes[idx];
		loadJobFromPipeline();
    }
}

void RelionMainWindow::cb_select_input_job(Fl_Widget* o, void* v)
{
	RelionMainWindow* T=(RelionMainWindow*)v;
	T->cb_select_input_job_i();
	run_button->activate();
}

void RelionMainWindow::cb_select_input_job_i()
{

	// Show the 'selected' group, hide the others
    int idx = input_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
    {
    	current_job = input_processes[idx];
		loadJobFromPipeline();
    }

}

void RelionMainWindow::cb_select_output_job(Fl_Widget* o, void* v)
{
	RelionMainWindow* T=(RelionMainWindow*)v;
	T->cb_select_output_job_i();
	run_button->activate();
}

void RelionMainWindow::cb_select_output_job_i()
{

	// Show the 'selected' group, hide the others
    int idx = output_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
    {
    	current_job = output_processes[idx];
		loadJobFromPipeline();
    }

}

void RelionMainWindow::cb_display_io_node(Fl_Widget* o, void* v)
{
	RelionMainWindow* T=(RelionMainWindow*)v;
	T->cb_display_io_node_i();
	run_button->activate();
}

void RelionMainWindow::cb_display_io_node_i()
{

	// Run relion_display on the output node
	int idx = display_io_node->value();
	long int mynode = io_nodes[idx];
	std::string command;

	if (pipeline.nodeList[mynode].type == NODE_MIC_COORDS)
	{
		FileName fn_job = ".gui_manualpickrun.job";
		bool iscont=false;
		if (exists(fn_job))
		{
			global_manualpickjob.read(fn_job.beforeLastOf("run.job").c_str(), iscont);
		}
		else
		{
			fl_message("ERROR: Save a Manual picking job parameter file (using the Save jobs settings option from the Jobs menu) before displaying coordinate files. ");
			return;
		}

		// Get the name of the micrograph STAR file from reading the suffix file
	    FileName fn_suffix = pipeline.nodeList[mynode].name;
	    if (fn_suffix.getExtension() == "star")
	    {
			std::ifstream in(fn_suffix.data(), std::ios_base::in);
			FileName fn_star;
			in >> fn_star ;
			in.close();
			if (fn_star != "")
			{
				FileName fn_dirs = fn_suffix.beforeLastOf("/")+"/";
				fn_suffix = fn_suffix.afterLastOf("/").without("coords_suffix_");
				fn_suffix = fn_suffix.withoutExtension();
				// Launch the manualpicker...
				command="`which relion_manualpick` --i " + fn_star;
				command += " --odir " + fn_dirs;
				command += " --pickname " + fn_suffix;
				command += " --scale " + floatToString(global_manualpickjob.micscale.getValue());
				command += " --sigma_contrast " + floatToString(global_manualpickjob.sigma_contrast.getValue());
				command += " --black " + floatToString(global_manualpickjob.black_val.getValue());
				command += " --white " + floatToString(global_manualpickjob.white_val.getValue());

				if (global_manualpickjob.lowpass.getValue() > 0.)
					command += " --lowpass " + floatToString(global_manualpickjob.lowpass.getValue());
				if (global_manualpickjob.highpass.getValue() > 0.)
					command += " --highpass " + floatToString(global_manualpickjob.highpass.getValue());
				if (global_manualpickjob.angpix.getValue() > 0.)
					command += " --angpix " + floatToString(global_manualpickjob.angpix.getValue());

				command += " --ctf_scale " + floatToString(global_manualpickjob.ctfscale.getValue());

				command += " --particle_diameter " + floatToString(global_manualpickjob.diameter.getValue());

				if (global_manualpickjob.do_color.getValue())
				{
					command += " --color_label " + global_manualpickjob.color_label.getValue();
					command += " --blue " + floatToString(global_manualpickjob.blue_value.getValue());
					command += " --red " + floatToString(global_manualpickjob.red_value.getValue());
					if (global_manualpickjob.fn_color.getValue().length() > 0)
						command += " --color_star " + global_manualpickjob.fn_color.getValue();
				}

				// Other arguments for extraction
				command += " " + global_manualpickjob.other_args.getValue() + " &";
			}
			else
			{
				fl_message("Only coordinates in .star format, generated in the pipeline, can be displayed here.");
			}
	    }
	    else
	    {
	    	fl_message("Only coordinates in .star format, generated in the pipeline, can be displayed here.");
	    }
	}
	else if (pipeline.nodeList[mynode].type == NODE_PDF_LOGFILE)
	{
		const char * default_pdf_viewer = getenv ("RELION_PDFVIEWER_EXECUTABLE");
		if (default_pdf_viewer == NULL)
		{
			char mydefault[]=DEFAULTPDFVIEWER;
			default_pdf_viewer=mydefault;
		}
		std::string myviewer(default_pdf_viewer);
		command = myviewer + " " + pipeline.nodeList[mynode].name + "&";
	}
	else
	{
		command = "relion_display --gui --i " + pipeline.nodeList[mynode].name + " &";
	}
	//std::cerr << " command= " << command << std::endl;
	int res= system(command.c_str());

}

void RelionMainWindow::cb_display(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_display_i();
}


void RelionMainWindow::cb_display_i()
{
        std::string command = " relion_display --gui &" ;
        int res = system(command.c_str());
}

void RelionMainWindow::cb_toggle_continue_i()
{

	if (is_main_continue)
	{
		run_button->label("Continue now");
		run_button->color(GUI_BUTTON_COLOR);
		run_button->labelfont(FL_ITALIC);
		run_button->labelsize(13);
		alias_current_job->deactivate();
	}
	else
	{
		run_button->label("Run now!");
		run_button->color(GUI_RUNBUTTON_COLOR);
		run_button->labelfont(FL_ITALIC);
		run_button->labelsize(16);
		alias_current_job->activate();
	}

	jobCommunicate(DONT_WRITE, DONT_READ, DO_TOGGLE_CONT, DONT_GET_CL, DONT_MKDIR);

}

void RelionMainWindow::fillStdOutAndErr()
{

	FileName fn_out = (current_job >= 0) ? pipeline.processList[current_job].name + "run.out" : "";
	FileName fn_err = (current_job >= 0) ? pipeline.processList[current_job].name + "run.err" : "";
	if (exists(fn_out))
	{
		// Remove annoying carriage returns
		std::string command = "awk -F\"\r\" '{if (NF>1) {print $NF} else {print}}' < " + fn_out + " | tail -6 > .gui_tmpout";
		int res = system(command.c_str());
		std::ifstream in(".gui_tmpout", std::ios_base::in);
		if (in.fail())
			REPORT_ERROR( (std::string) "MetaDataTable::read: File .gui_tmpout does not exists" );
		int err = textbuff_stdout->loadfile(".gui_tmpout");
		in.close();
		// Scroll to the bottom
		disp_stdout->insert_position(textbuff_stdout->length()-1);
		disp_stdout->show_insert_position();
	}
	else
		textbuff_stdout->text("stdout will go here; double-click this window to open stdout in a separate window");

	if (exists(fn_err))
	{
		std::string command = "tail -3 " + fn_err + " > .gui_tmperr";
		int res = system(command.c_str());
		std::ifstream in(".gui_tmperr", std::ios_base::in);
		if (in.fail())
			REPORT_ERROR( (std::string) "MetaDataTable::read: File .gui_tmperr does not exists" );
		int err = textbuff_stderr->loadfile(".gui_tmperr");
		in.close();
		// Scroll to the bottom
		disp_stderr->insert_position(textbuff_stderr->length()-1);
		disp_stderr->show_insert_position();
	}
	else
		textbuff_stderr->text("stderr will go here; double-click this window to open stderr in a separate window");

}

void RelionMainWindow::cb_print_cl(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_print_cl_i();
}

void RelionMainWindow::cb_print_cl_i()
{
	// Update all job lists in the main GUI
	updateJobLists();

	jobCommunicate(DONT_WRITE, DONT_READ, DONT_TOGGLE_CONT, DO_GET_CL, DONT_MKDIR);
    std::cout << " *** The command is:" << std::endl;
    for (int icom = 0; icom < commands.size(); icom++)
    	std::cout << commands[icom] << std::endl;

}

// Run button call-back functions
void RelionMainWindow::cb_forgot(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_forgot_i(); // 1st true means only_schedule, do not run, 2nd true means open the note editor window
}

void RelionMainWindow::cb_forgot_i()
{
	fl_message("Really?! Perhaps you should spend fewer nights at the microscope and try to sleep a bit more...");
}

// Run button call-back functions
void RelionMainWindow::cb_run(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
	// Deactivate Run button to prevent the user from accidentally submitting many jobs
	run_button->deactivate();
	// Run the job
	T->cb_run_i(false, true); // false means dont only_schedule, true means open the note editor window
}

// Run button call-back functions
void RelionMainWindow::cb_schedule(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_run_i(true, true); // 1st true means only_schedule, do not run, 2nd true means open the note editor window
}

void RelionMainWindow::cb_run_i(bool only_schedule, bool do_open_edit)
{

	// Make sure that whenever we write out a pipeline, we first read in the existing version on disk
	// This prevents sync errors with multiple windows acting on the same pipeline

	// Read in the latest version of the pipeline, just in case anyone else made a change meanwhile...
	pipeline.read(DO_LOCK); // true means: only_read if_file_exists

	// Get the command line arguments from the currently active jobwindow,
	// If the job gets cancelled, just exit this function
	if (!jobCommunicate(DONT_WRITE, DONT_READ, DONT_TOGGLE_CONT, DO_GET_CL, DO_MKDIR))
	{
		std::cout << " Cancelling job" << std::endl;
		return;
	}

	// Save temporary hidden file with this jobs settings as default for a new job
	fn_settings = "";
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DONT_MKDIR);

	// Also save a copy of the GUI settings with the current output name
	fn_settings = global_outputname;
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DONT_MKDIR);

	if (commands.size()==0)
	{
		std::cout << " Nothing to do..."<< std::endl;
		return;
	}

	// If this is a continuation job, check whether output files exist and move away!
	// This is to ensure that the continuation job goes OK and will show up as 'running' in the GUI
	if (!only_schedule && is_main_continue)
	{
		bool skip_this = (pipeline.processList[current_job].type == PROC_2DCLASS ||
				pipeline.processList[current_job].type == PROC_3DCLASS ||
				pipeline.processList[current_job].type == PROC_3DAUTO ||
				pipeline.processList[current_job].type == PROC_MANUALPICK ||
				pipeline.processList[current_job].type == PROC_CLASSSELECT ||
				pipeline.processList[current_job].type == PROC_MOVIEREFINE ||
				pipeline.processList[current_job].type == PROC_POLISH);
		if (!skip_this )
		{
			for (int i = 0; i < pipeline.processList[current_job].outputNodeList.size(); i++)
			{
				int j = pipeline.processList[current_job].outputNodeList[i];
				std::string fn_node = pipeline.nodeList[j].name;
				if (exists(fn_node))
				{
					std::string path2 =  fn_node + ".old";
					rename(fn_node.c_str(), path2.c_str());
				}
			}
		}

		// For continuation of relion_refine jobs, remove the original output nodes from the list
		if (pipeline.processList[current_job].type == PROC_2DCLASS ||
				pipeline.processList[current_job].type == PROC_3DCLASS ||
				pipeline.processList[current_job].type == PROC_3DAUTO)
		{

			std::vector<bool> deleteNodes, deleteProcesses;
			deleteNodes.resize(pipeline.nodeList.size(), false);
			deleteProcesses.resize(pipeline.processList.size(), false);

			for (long int inode = 0; inode < (pipeline.processList[current_job]).outputNodeList.size(); inode++)
			{
				long int mynode = (pipeline.processList[current_job]).outputNodeList[inode];
				if(!exists(pipeline.nodeList[mynode].name))
					deleteNodes[mynode] = true;
			}

			FileName fn_del = "tmp";
			pipeline.write(DO_LOCK, fn_del, deleteNodes, deleteProcesses);
			std::remove("tmpdeleted_pipeline.star");

			// Read the updated pipeline back in again
			pipeline.read(DO_LOCK);

		}

	}

	// Now save the job (and its status) to the PipeLine
	int mynewstatus;
	if (only_schedule)
		mynewstatus = (is_main_continue) ? PROC_SCHEDULED_CONT : PROC_SCHEDULED_NEW;
	else
		mynewstatus = PROC_RUNNING;
	bool allow_overwrite = is_main_continue; // continuation jobs always allow overwriting into the existing directory
	if (current_job >= 0 && !allow_overwrite) // and so do scheduled jobs
		allow_overwrite = (pipeline.processList[current_job].status == PROC_SCHEDULED_NEW ||
				 pipeline.processList[current_job].status == PROC_SCHEDULED_CONT) ;
	current_job = addToPipeLine(mynewstatus, allow_overwrite);

	// Update all job lists in the main GUI
	updateJobLists();

	// Write out the new pipeline
	pipeline.write(DO_LOCK);

	if (!only_schedule)
	{
		//std::cout << "Executing: " << final_command << std::endl;
		int res = system(final_command.c_str());

		// Also print the final_command to the note for future reference
		FileName fn_note = pipeline.processList[current_job].name + "note.txt";
		std::ofstream ofs;
		ofs.open (fn_note.c_str(), std::ofstream::out | std::ofstream::app);

		// current date/time based on current system
		time_t now = time(0);
		ofs << std::endl << " ++++ Executing new job on " << ctime(&now);
		ofs <<  " ++++ with the following command(s): " << std::endl;
		for (size_t icom = 0; icom < commands.size(); icom++)
			ofs << commands[icom] << std::endl;
		ofs <<  " ++++ " << std::endl;
		ofs.close();
	}

	// Open the edit note window
	if (do_open_edit)
	{
		// Open the note editor window
		cb_edit_note_i();

		// Also set alias from the alias_current_job input
		if (!is_main_continue)
		{
			std::string alias= (std::string)alias_current_job->value();
			if (alias != "Give_alias_here" && alias != pipeline.processList[current_job].name)
				cb_set_alias_i(alias);
		}
	}

	// Copy pipeline star file as backup to the output directory
	FileName fn_pipe = pipeline.name + "_pipeline.star";
	if (exists(fn_pipe))
		copy(fn_pipe, pipeline.processList[current_job].name + fn_pipe);

}


// Run button call-back functions
void RelionMainWindow::cb_delete(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_delete_i();
}

void RelionMainWindow::cb_delete_i(bool do_ask, bool do_recursive)
{

	if (current_job < 0)
	{
		std::cout << " You can only delete existing jobs ... " << std::endl;
		return;
	}

	std::vector<bool> deleteProcesses, deleteNodes;
	deleteProcesses.resize(pipeline.processList.size(), false);
	deleteNodes.resize(pipeline.nodeList.size(), false);
	FileName fn_del = pipeline.processList[current_job].name;

	std::vector<long int> to_delete_processes;
	to_delete_processes.push_back(current_job);

	bool is_done = false;
	size_t istart = 0;
	while (!is_done)
	{
		size_t imax = to_delete_processes.size();
		for (long int i = istart; i < imax; i++)
		{
			// re-set istart for next recursive round
			istart = imax;
			long int idel = to_delete_processes[i];
			deleteProcesses[idel] = true;
                        is_done = true;
                        for (size_t inode = 0; inode < (pipeline.processList[idel]).outputNodeList.size(); inode++)
			{
				long int mynode = (pipeline.processList[idel]).outputNodeList[inode];
				deleteNodes[mynode] = true;

				if (do_recursive)
				{
					// Check whether this node is being used as input for another process, and if so, delete those as well
					for (size_t ii = 0; ii < (pipeline.nodeList[mynode]).inputForProcessList.size(); ii++)
					{
						long int iproc = (pipeline.nodeList[mynode]).inputForProcessList[ii];
						// See if this process is not already in the list to be deleted
						bool already_in = false;
						for (size_t j = 0; j < to_delete_processes.size(); j++)
						{
							if (to_delete_processes[j] == iproc)
								already_in = true;
						}
						if (!already_in)
						{
							to_delete_processes.push_back(iproc);
							is_done = false;
						}
					}
				}
			}
		}
	}

	// Before we do anything: confirm this is really what the user wants to do....
	int proceed;
	if (do_ask)
	{
		std::string ask;
		ask = "Are you sure you want to move the following processes to Trash? \n";
		for (size_t i = 0; i < deleteProcesses.size(); i++)
		{
			if (deleteProcesses[i])
			{
				std::string name = (pipeline.processList[i].alias == "None") ? pipeline.processList[i].name : pipeline.processList[i].alias;
				ask += " - " + name + "\n";
			}
		}
		proceed =  fl_choice("%s", "Don't move", "Move", NULL, ask.c_str());
	}
	else
	{
		proceed = 1;
	}

	if (proceed)
	{

		// Read in existing pipeline, in case some other window had changed it
		pipeline.read(DO_LOCK);

		// Write new pipeline without the deleted processes and nodes to disc and read in again
		pipeline.write(DO_LOCK, fn_del, deleteNodes, deleteProcesses);

		// Delete the output directories for all selected processes from the hard disk
		// Do this after pipeline.write to get the deleted_pipeline.star still in the correct directory
		for (int i = 0; i < deleteProcesses.size(); i++)
		{
			if (deleteProcesses[i])
			{
				FileName alldirs = pipeline.processList[i].name;
				alldirs = alldirs.beforeLastOf("/");
				// Move entire output directory (with subdirectory structure) to the Trash folder
				FileName firstdirs = alldirs.beforeLastOf("/");
				FileName fn_tree="Trash/" + firstdirs;
				int res = mktree(fn_tree);
				std::string command = "mv -f " + alldirs + " " + "Trash/" + firstdirs+"/.";
				res = system(command.c_str());
				// Also remove the symlink if it exists
				FileName fn_alias = (pipeline.processList[i]).alias;
				if (fn_alias != "None")
				{
					int res = unlink((fn_alias.beforeLastOf("/")).c_str());
				}

				pipeline.deleteTemporaryNodeFiles(pipeline.processList[i]);

			}
		}

		// Read new pipeline back in again
		pipeline.read(DO_LOCK);

		// Update all job lists in the main GUI
		updateJobLists();

		pipeline.write(DO_LOCK);
	}

}

// Run button call-back functions
void RelionMainWindow::cb_gently_clean_all_jobs(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_clean_all_jobs_i(false);
}

// Run button call-back functions
void RelionMainWindow::cb_harshly_clean_all_jobs(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_clean_all_jobs_i(true);
}

void RelionMainWindow::cb_clean_all_jobs_i(bool do_harsh)
{


	int proceed = 1;
	std::string ask;
	if (do_harsh)
	{
		ask = "Are you sure you want to harshly clean up intermediate files from the entire pipeline? \n\n\
Harsh cleaning will remove micrographs, movies and particle stacks from all MotionCorr, Extract, MovieRefine, \n\
Polish and Subtract directories. This means you will NOT be able to use those images in subsequent runs anymore, \n\
although you could always recreate the data by continuing the job (possibly at considerable computing costs).\n \n \
You can protect specific jobs from harsh cleaning by creating a file called \"NO_HARSH_CLEAN\" inside their directory,\n\
e.g. by using \"touch Polish/job045/NO_HARSH_CLEAN\". Below is a list of currently protected jobs (if any):\n \n";

		for (int myjob = 0; myjob < pipeline.processList.size(); myjob++)
		{
			if (pipeline.processList[myjob].status == PROC_FINISHED &&
					(pipeline.processList[myjob].type == PROC_MOTIONCORR ||
					pipeline.processList[myjob].type == PROC_EXTRACT ||
					pipeline.processList[myjob].type == PROC_MOVIEREFINE ||
					pipeline.processList[myjob].type == PROC_POLISH ||
					pipeline.processList[myjob].type == PROC_SUBTRACT) )
			{
				if (exists(pipeline.processList[myjob].name + "NO_HARSH_CLEAN"))
					ask += pipeline.processList[myjob].name + " \n";
			}
		}
	}
	else
	{
		ask = "Are you sure you want to gently clean up intermediate files from the entire pipeline?";
	}

	proceed = fl_choice("%s", "Don't clean up", "Clean up", NULL, ask.c_str());
	if (proceed)
	{
		std::string how = (do_harsh) ? "Harshly" : "Gently";
		std::cout << how << " cleaning all finished jobs ..." << std::endl;
		for (int myjob = 0; myjob < pipeline.processList.size(); myjob++)
		{
			if (pipeline.processList[myjob].status == PROC_FINISHED)
			{
				if (do_harsh && exists(pipeline.processList[myjob].name + "NO_HARSH_CLEAN"))
					continue;
				cb_cleanup_i(myjob, false, do_harsh); // false means no verb
			}
		}
		fl_message("Done cleaning! Don't forget the files are all still in the Trash folder. Use the \"Empty Trash\" option from the File menu to permanently delete them.");
	}
}

// Run button call-back functions
void RelionMainWindow::cb_gentle_cleanup(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_cleanup_i(-1, true, false);
}

void RelionMainWindow::cb_harsh_cleanup(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_cleanup_i(-1, true, true);
}

void RelionMainWindow::cb_cleanup_i(int myjob, bool do_verb, bool do_harsh)
{
	// Allow cleaning the currently selected job from the GUI
	if (myjob < 0)
		myjob = current_job;

	if (myjob < 0 || pipeline.processList[myjob].status != PROC_FINISHED)
	{
		fl_message(" You can only clean up finished jobs ... ");
		return;
	}

	// These job types do not have cleanup:
	if (pipeline.processList[myjob].type == PROC_IMPORT ||
		pipeline.processList[myjob].type == PROC_MANUALPICK ||
		pipeline.processList[myjob].type == PROC_SORT ||
		pipeline.processList[myjob].type == PROC_CLASSSELECT ||
		pipeline.processList[myjob].type == PROC_MASKCREATE ||
		pipeline.processList[myjob].type == PROC_JOINSTAR ||
		pipeline.processList[myjob].type == PROC_RESMAP)
		return;

	int proceed = 1;
	if (do_verb)
	{
		std::string ask;
		ask = "Are you sure you want to clean up intermediate files from " + pipeline.processList[myjob].name + "?";
		proceed = fl_choice("%s", "Don't clean up", "Clean up", NULL, ask.c_str());
	}

	if (proceed)
	{

		// Find any subdirectories
		std::vector<FileName> fns_subdir;
		FileName fn_curr_dir = "";
		int idir = -1, istop = 0;
		bool is_first = true;
		// Recursively find all subdirectories
		while (idir < istop)
		{
			FileName fn_curr_dir = (is_first) ? pipeline.processList[myjob].name : pipeline.processList[myjob].name + fns_subdir[idir];
			DIR *dir = opendir(fn_curr_dir.c_str());
			struct dirent *entry = readdir(dir);
			while (entry != NULL)
			{
				// Only want directories, and not '.' or '..'
				if (entry->d_type == DT_DIR && (entry->d_name[0] != '.'))
				{
					FileName fnt = (is_first) ? entry->d_name : fns_subdir[idir] + entry->d_name;
					fns_subdir.push_back(fnt + "/");
				}
				entry = readdir(dir);
			}
			closedir(dir);
			istop = fns_subdir.size();
			idir++;
			is_first = false;
		}

		std::vector<FileName> fns_del;
		FileName fn_pattern;

		// In all jobs cleanup the .old files (from continuation runs)
		fn_pattern = pipeline.processList[myjob].name + "*.old";
		fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del

		////////// Now see which jobs needs cleaning up
		if (pipeline.processList[myjob].type == PROC_MOTIONCORR)
		{

			for (int idir = 0; idir < fns_subdir.size(); idir++)
			{
				if (do_harsh)
				{
					//remove entire directory
					fns_del.push_back(pipeline.processList[myjob].name + fns_subdir[idir]);
				}
				else
				{
					fn_pattern = pipeline.processList[myjob].name + fns_subdir[idir] + "*.com";
					fn_pattern.globFiles(fns_del, false);
					fn_pattern = pipeline.processList[myjob].name + fns_subdir[idir] + "*.err";
					fn_pattern.globFiles(fns_del, false);
					fn_pattern = pipeline.processList[myjob].name + fns_subdir[idir] + "*.out";
					fn_pattern.globFiles(fns_del, false);
					fn_pattern = pipeline.processList[myjob].name + fns_subdir[idir] + "*.log";
					fn_pattern.globFiles(fns_del, false);
				}
			}

		} // end if motioncorr
		else if (pipeline.processList[myjob].type == PROC_CTFFIND)
		{

			fn_pattern = pipeline.processList[myjob].name + "gctf*.out";
			fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
			fn_pattern = pipeline.processList[myjob].name + "gctf*.err";
			fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
			for (int idir = 0; idir < fns_subdir.size(); idir++)
			{
				//remove entire Micrographs directory structure
				fns_del.push_back(pipeline.processList[myjob].name + fns_subdir[idir]);
			}

		} // end if ctffind
		else if (pipeline.processList[myjob].type == PROC_AUTOPICK)
		{

			for (int idir = 0; idir < fns_subdir.size(); idir++)
			{
				// remove the Spider files with the FOM maps
				fn_pattern = pipeline.processList[myjob].name + fns_subdir[idir] + "*.spi";
				fn_pattern.globFiles(fns_del, false);
			}

		} // end if autopick
		else if (pipeline.processList[myjob].type == PROC_EXTRACT)
		{

			for (int idir = 0; idir < fns_subdir.size(); idir++)
			{
				if (do_harsh)
				{
					//remove entire directory (STAR files and particle stacks!
					fns_del.push_back(pipeline.processList[myjob].name + fns_subdir[idir]);
				}
				else
				{
					// only remove the STAR files with the metadata (this will only give moderate file savings)
					fn_pattern = pipeline.processList[myjob].name + fns_subdir[idir] + "*_extract.star";
					fn_pattern.globFiles(fns_del, false);
				}
			}

		} // end if extract
		else if (pipeline.processList[myjob].type == PROC_2DCLASS ||
				 pipeline.processList[myjob].type == PROC_3DCLASS ||
				 pipeline.processList[myjob].type == PROC_3DAUTO)
	    {

			// First find the _data.star from each iteration
			std::vector<FileName> fns_iter;
			fn_pattern = pipeline.processList[myjob].name + "*_it[0-9][0-9][0-9]_data.star";
			fn_pattern.globFiles(fns_iter);
			for (int ifile = 0; ifile < fns_iter.size(); ifile++)
			{
				FileName fn_file = (fns_iter[ifile]).without("_data.star");
				// Find the iterations to keep: i.e. those that are part of the pipeline
				bool is_in_pipeline = false;
				for (long int inode = 0; inode < pipeline.nodeList.size(); inode++)
				{
					FileName fn_node = pipeline.nodeList[inode].name;
					if (fn_node.contains(fn_file))
					{
						is_in_pipeline = true;
						break;
					}
				}
				// Delete all files from this iteration
				if (!is_in_pipeline)
				{
					fn_pattern = fn_file + "*";
					fn_pattern.globFiles(fns_del, false);
				}

			} //end loop over ifile (i.e. the _data.star files from all iterations)

		} // end if refine job
		else if (pipeline.processList[myjob].type == PROC_MOVIEREFINE)
		{

			fn_pattern = pipeline.processList[myjob].name + "batch*mics_nr[0-9][0-9][0-9].star";
			fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
			fn_pattern = pipeline.processList[myjob].name + "run_it[0-9][0-9][0-9]*";
			fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
			for (int idir = 0; idir < fns_subdir.size(); idir++)
			{
				if (do_harsh)
				{
					//remove entire Micrographs directory (STAR files and particle stacks!)
					fns_del.push_back(pipeline.processList[myjob].name + fns_subdir[idir]);
				}
				else
				{
					// only remove the STAR files with the metadata (this will only give moderate file savings)
					fn_pattern = pipeline.processList[myjob].name + fns_subdir[idir] + "*_extract.star";
					fn_pattern.globFiles(fns_del, false);
				}
			}

		} // end if movierefine
		else if (pipeline.processList[myjob].type == PROC_POLISH)
		{

			fn_pattern = pipeline.processList[myjob].name + "*.mrc";
			fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
			fn_pattern = pipeline.processList[myjob].name + "shiny*xml";
			fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
			for (int idir = 0; idir < fns_subdir.size(); idir++)
			{
				if (do_harsh)
				{
					//remove entire Micrographs directory (STAR files and particle stacks!)
					fns_del.push_back(pipeline.processList[myjob].name + fns_subdir[idir]);
				}
				else
				{
					// only remove the STAR files with the metadata (this will only give moderate file savings)
					fn_pattern = pipeline.processList[myjob].name + fns_subdir[idir] + "*.star";
					fn_pattern.globFiles(fns_del, false);
				}
			}

		} // end if polish
		else if (pipeline.processList[myjob].type == PROC_SUBTRACT)
		{

			if (do_harsh)
			{
				fn_pattern = pipeline.processList[myjob].name + "subtracted.*";
				fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del
			}

		} // end if subtract
		else if (pipeline.processList[myjob].type == PROC_POST)
		{

			fn_pattern = pipeline.processList[myjob].name + "*masked.mrc";
			fn_pattern.globFiles(fns_del, false); // false means do not clear fns_del

		} // end if postprocess


		// Now actually move all the files
		FileName fn_old_dir = "";
		for (long int idel = 0; idel < fns_del.size(); idel++)
		{
			FileName fn_dest = "Trash/" + fns_del[idel];
			FileName fn_dir = fn_dest.beforeLastOf("/");
			if (fn_dir != fn_old_dir && ! exists(fn_dir))
				int res = mktree(fn_dir);
			// by removing entire directories, it could be the file is gone already
			if (exists(fns_del[idel]))
			{
				std::string command = "mv -f " + fns_del[idel] + " "+ fn_dir;
				int res = system(command.c_str());
			}
		} // end loop over all files to be deleted

	} // end if proceed

}


// Run button call-back functions
void RelionMainWindow::cb_set_alias(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_set_alias_i();
}

void RelionMainWindow::cb_set_alias_i(std::string alias)
{

	FileName fn_pre, fn_jobnr, fn_post, fn_dummy, default_ask;
	if (!decomposePipelineFileName(pipeline.processList[current_job].name, fn_pre, fn_jobnr, fn_post))
		REPORT_ERROR("RelionMainWindow::cb_set_alias_i ERROR: invalid pipeline process name: " + pipeline.processList[current_job].name);

	// If alias already exists: remove that symlink
	FileName fn_alias = pipeline.processList[current_job].alias;
	FileName fn_old_alias="";
	if (fn_alias != "None")
	{
		fn_old_alias = fn_alias.beforeLastOf("/");
		default_ask = fn_alias.without(fn_pre);
		if (default_ask[default_ask.length()-1] == '/')
			default_ask = default_ask.beforeLastOf("/");
	}
	else
		default_ask = fn_jobnr.beforeLastOf("/");

	bool is_done = false;
	while (!is_done)
	{
		// If the alias already contains a uniquedate string it may be a continuation of a relion_refine job
		// (where alias_current_job contains a different uniqdate than the outputname of the job)
		if (alias == "" || decomposePipelineFileName(alias, fn_dummy, fn_dummy, fn_dummy) ) // if an alias is provided, just check it is unique, otherwise ask
		{
			const char * palias;
			palias =  fl_input("Rename to: ", default_ask.c_str());
			// TODO: check for cancel
			if (palias == NULL)
				return;
			std::string al2(palias);
			alias = al2;
		}

		if (alias.length() < 2)
		{
			 fl_message("Alias cannot be less than 2 characters, please provide another one");
			 alias="";
		}
		else if (alias.length() > 2 && alias[0]=='j' && alias[1]=='o' && alias[2]=='b')
		{
			 fl_message("Alias cannot start with 'job', please provide another one");
			 alias="";
		}
		else
		{

			// Read in existing pipeline, in case some other window had changed it
			pipeline.read(DO_LOCK);

			//remove spaces from any potential alias
			for (int i = 0; i < alias.length(); i++)
			{
				if (alias[i] == ' ')
					alias[i] = '_';
			}

			// Make sure the alias ends with a slash
			if (alias[alias.length()-1] != '/')
				alias += "/";

			// Check uniqueness of the alias
			bool is_unique = true;
			for (size_t i = 0; i < pipeline.processList.size(); i++)
			{
				if ( pipeline.processList[i].alias == fn_pre + alias && alias != "None")
				{
					is_unique = false;
					break;
				}
			}
			if (!is_unique || alias.length() < 1)
			{
				 fl_message("Alias is not unique, please provide another one");
				 alias="";
			}
			else
				is_done = true;
		}
	}

	// Remove the original .Nodes entry
	pipeline.deleteTemporaryNodeFiles(pipeline.processList[current_job]);

	// No alias if the alias contains a unique jobnr string because of continuation of relion_refine jobs
	// (where alias_current_job contains a different uniq jobnr than the outputname of the job)
	if (alias == "None" )
	{
		pipeline.processList[current_job].alias = "None";
	}
	else
	{
		// If this was already an alias: remove the old symbolic link
		if (fn_old_alias != "")
		{
			int res2 = unlink(fn_old_alias.c_str());
		}

		// Set the alias in the pipeline
		pipeline.processList[current_job].alias = fn_pre + alias;

		//Make the new symbolic link
		FileName path1 = "../" + pipeline.processList[current_job].name;
		FileName path2 = pipeline.processList[current_job].alias;
		int res = symlink(path1.c_str(), path2.beforeLastOf("/").c_str());

	}

	// Remake the new .Nodes entry
	pipeline.touchTemporaryNodeFiles(pipeline.processList[current_job]);

	// Update the name in the lists
	updateJobLists();

	// Write new pipeline to disc
	pipeline.write(DO_LOCK);

}



// Run button call-back functions
void RelionMainWindow::cb_mark_as_finished(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_mark_as_finished_i();
}

void RelionMainWindow::cb_mark_as_finished_i()
{

	if (current_job < 0)
	{
		fl_message("You can only mark existing jobs as finished!");
		return;
	}

	// Read in existing pipeline, in case some other window had changed it
	pipeline.read(DO_LOCK);

	pipeline.processList[current_job].status = PROC_FINISHED;

	// For relion_refine jobs, add last iteration optimiser.star, data.star, model.star and class???.mrc to the pipeline
	if (pipeline.processList[current_job].type == PROC_2DCLASS ||
		pipeline.processList[current_job].type == PROC_3DCLASS ||
		pipeline.processList[current_job].type == PROC_3DAUTO)
	{
		// Get the last iteration optimiser file
		FileName fn_opt;
		FileName fn_root1 = (pipeline.processList[current_job].alias != "None") ? pipeline.processList[current_job].alias : pipeline.processList[current_job].name;

		std::vector<FileName> fn_opts;
		fn_opt = fn_root1 + "run_it*optimiser.star";
		fn_opt.globFiles(fn_opts);
		// It could also be a continuation
		fn_opt = fn_root1 + "run_ct?_it???_optimiser.star";
		fn_opt.globFiles(fn_opts, false); // false means: don't clear fn_opts vector
		// It could also be a continuation
		fn_opt = fn_root1 + "run_ct??_it???_optimiser.star";
		fn_opt.globFiles(fn_opts, false); // false means: don't clear fn_opts vector
		if (fn_opts.size() > 0)
		{

			fn_opt = fn_opts[fn_opts.size()-1]; // the last one

			// Also get data.star
			FileName fn_data = fn_opt.without("_optimiser.star") + "_data.star";
			Node node2(fn_data, NODE_PART_DATA);
			pipeline.addNewOutputEdge(current_job, node2);

			FileName fn_root = fn_opt.without("_optimiser.star");
			if (pipeline.processList[current_job].type == PROC_3DAUTO)
				fn_root += "_half1";

			FileName fn_model = fn_root + "_model.star";
			Node node3(fn_model, NODE_MODEL);
			pipeline.addNewOutputEdge(current_job, node3);


			FileName fn_map = fn_root + "_class???.mrc";
			std::vector<FileName> fn_maps;
			fn_map.globFiles(fn_maps);
			for (int i = 0; i < fn_maps.size(); i++)
			{
				Node node4(fn_maps[i], NODE_3DREF);
				pipeline.addNewOutputEdge(current_job, node4);
			}
		}
		else
		{
			fl_message(" You are trying to mark a relion_refine job as finished that hasn't even started. \n This will be ignored. Perhaps you wanted to delete it instead?");
			pipeline.processList[current_job].status = PROC_RUNNING;
		}
	}

	// Remove any of the expected output nodes from the pipeline if the corresponding file doesn't already exist
	std::vector<bool> deleteNodes, deleteProcesses;
	deleteNodes.resize(pipeline.nodeList.size(), false);
	deleteProcesses.resize(pipeline.processList.size(), false);

	for (long int inode = 0; inode < (pipeline.processList[current_job]).outputNodeList.size(); inode++)
	{
		long int mynode = (pipeline.processList[current_job]).outputNodeList[inode];
		if(!exists(pipeline.nodeList[mynode].name))
			deleteNodes[mynode] = true;
	}
	FileName fn_del = "tmp";
	pipeline.write(DO_LOCK, fn_del, deleteNodes, deleteProcesses);
	std::remove("tmpdeleted_pipeline.star");

	// Read the updated pipeline back in again
	pipeline.read(DO_LOCK);

	// Update all job lists in the main GUI
	updateJobLists();

	// With the locking mechanism, each pipeline.read(bool, DO_LOCK) needs to be followed soon by a pipeline.write(DO_LOCK)!
	pipeline.write(DO_LOCK);

}

// Run button call-back functions
void RelionMainWindow::cb_make_flowchart(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_make_flowchart_i();
}

void RelionMainWindow::cb_make_flowchart_i()
{
	if (current_job < 0)
	{
		fl_message(" You can only make flowcharts for existing jobs ... ");
		return;
	}

	const char * default_pdf_viewer = getenv ("RELION_CTFFIND_EXECUTABLE");
	if (default_pdf_viewer == NULL)
	{
		char mydefault[]=DEFAULTPDFVIEWER;
		default_pdf_viewer=mydefault;
	}
	std::string myviewer(default_pdf_viewer);


	PipeLineFlowChart flowchart;
	FileName fn_dir = pipeline.processList[current_job].name;
	FileName fn_out = "flowchart.tex";
	flowchart.makeAllUpwardsFlowCharts(fn_out, pipeline, current_job);
	std::string command = "latex flowchart.tex > flowchart.log && dvipdf flowchart >>flowchart.log && mv flowchart* " + fn_dir;
	std:: cout << " Executing: " << command << std::endl;
	int res = std::system(command.c_str());
	command = myviewer + " " + fn_dir + "flowchart.pdf &";
	res = std::system(command.c_str());

	// Read in existing pipeline, in case some other window had changed it
	pipeline.read(DO_LOCK);

	// Add the PDF file as a logfile to the outputnodes of this job, so it can be visualised from the Display button
	Node node(fn_dir+"flowchart.pdf", NODE_PDF_LOGFILE);
	pipeline.addNewOutputEdge(current_job, node);
	updateJobLists();

	pipeline.write(DO_LOCK);

	return;
}

void RelionMainWindow::cb_edit_note(Fl_Widget*, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_edit_note_i();

}

void RelionMainWindow::cb_edit_project_note(Fl_Widget*, void* v)
{
  	RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_edit_note_i(true); // true means is_project_note
}

void RelionMainWindow::cb_edit_note_i(bool is_project_note)
{

	FileName fn_note;
	std::string title;
	if (is_project_note)
	{
		fn_note = "project_note.txt";
		title = "Overall project notes";
	}
	else
	{
		if (current_job < 0)
		{
			fl_message(" You can only edit the note for existing jobs ... ");
			return;
		}
		fn_note = pipeline.processList[current_job].name + "note.txt";
		title = (pipeline.processList[current_job].alias == "None") ? pipeline.processList[current_job].name : pipeline.processList[current_job].alias;
	}
	NoteEditorWindow* w = new NoteEditorWindow(660, 400, title.c_str(), fn_note);
	w->show();

}


// Save button call-back function
void RelionMainWindow::cb_save(Fl_Widget* o, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_save_i();
}

void RelionMainWindow::cb_save_i()
{
	// For scheduled jobs, also allow saving the .job file in the output directory
	if (current_job >= 0 && (pipeline.processList[current_job].status == PROC_SCHEDULED_NEW ||
			                pipeline.processList[current_job].status == PROC_SCHEDULED_CONT))
	{
		fn_settings = pipeline.processList[current_job].name;
		jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);
	}

	fn_settings = "";
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

}

// Load button call-back function
void RelionMainWindow::cb_load(Fl_Widget* o, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_load_i();
}

void RelionMainWindow::cb_load_i()
{

	fn_settings = "";
	jobCommunicate(DONT_WRITE, DO_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DONT_MKDIR);

}

// Load button call-back function
void RelionMainWindow::cb_undelete_job(Fl_Widget* o, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_undelete_job_i();
}

void RelionMainWindow::cb_undelete_job_i()
{

	std::string fn_dir = "./Trash/.";
	std::string fn_filter = "Pipeline STAR files (job_pipeline.star)";
	Fl_File_Chooser chooser(fn_dir.c_str(),  fn_filter.c_str(), Fl_File_Chooser::SINGLE, "Choose pipeline STAR file to import");
	chooser.show();
	// Block until user picks something.
	while(chooser.shown())
		{ Fl::wait(); }

	// User hit cancel?
	if ( chooser.value() == NULL )
		return;

	char relname[FL_PATH_MAX];
    fl_filename_relative(relname,sizeof(relname),chooser.value());
	FileName fn_pipe(relname);


	// Read in existing pipeline, in case some other window had changed it
	pipeline.read(DO_LOCK);

    pipeline.importPipeline(fn_pipe.beforeLastOf("_pipeline.star"));

	// Copy all processes in the STAR file back into the ProjectDirectory
	MetaDataTable MDproc;
	MDproc.read(fn_pipe, "pipeline_processes");
	std::cout <<"  Undeleting from Trash ... " << std::endl;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDproc)
	{
		FileName fn_proc;
		MDproc.getValue(EMDL_PIPELINE_PROCESS_NAME, fn_proc);

		// Copy the job back from the Trash folder
		FileName fn_dest = fn_proc.beforeLastOf("/"); //gets rid of ending "/"
		FileName fn_dir_dest = fn_dest.beforeLastOf("/"); // Now only get the job-type directory
		if (!exists(fn_dir_dest))
		{
			mktree(fn_dir_dest);
		}
		std::string command = "mv Trash/" + fn_dest + " " + fn_dest;
		std::cout << command << std::endl;
		int res = system(command.c_str());

		// Also re-make all entries in the .Nodes directory
		long int myproc = pipeline.findProcessByName(fn_proc);
		pipeline.touchTemporaryNodeFiles(pipeline.processList[myproc]);
	}
	std::cout << " Done undeleting! " << std::endl;

	// Write the new pipeline to disk and reread it back in again
	pipeline.write(DO_LOCK);

}


void replaceFilesForImportExportOfScheduledJobs(FileName fn_in_dir, FileName fn_out_dir, std::vector<std::string> &find_pattern, std::vector<std::string> &replace_pattern)
{
	int res;
	std::string command;
	std::vector<std::string> myfiles;
	myfiles.push_back("run.job");
	myfiles.push_back("note.txt");
	myfiles.push_back("job_pipeline.star");

	// Copy the run.job, the note.txt and the job_pipeline.star
	// Replace all instances of all find_pattern's with the replace_pattern's
	for (int ifile = 0; ifile < myfiles.size(); ifile++)
	{
		for (int ipatt = 0; ipatt < find_pattern.size(); ipatt++)
		{
			FileName outfile = fn_out_dir + myfiles[ifile];
			FileName tmpfile = fn_out_dir + "tmp";
			FileName infile = (ipatt == 0) ? fn_in_dir + myfiles[ifile] : tmpfile;
			// Create directory first time round
			if (ipatt == 0)
			{
				FileName dirs = outfile.beforeLastOf("/");
				command =  "mkdir -p " + dirs;
				res = system(command.c_str());
			}
			command =  "sed 's|" + find_pattern[ipatt] + "|" + replace_pattern[ipatt] + "|g' < " + infile + " > " + outfile;
			//std::cerr << " Executing: " << command<<std::endl;
			res = system(command.c_str());
			if (ipatt+1 < find_pattern.size())
			{
				std::rename(outfile.c_str(), tmpfile.c_str());
				//std::cerr << " Excuting: mv " << outfile<<" "<<tmpfile<<std::endl;
			}
		}
	}

}

void RelionMainWindow::cb_export_jobs(Fl_Widget* o, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_export_jobs_i();
}

void RelionMainWindow::cb_export_jobs_i()
{
	// Get the name of this block of exported jobs and make the corresponding directory
	const char * answer;
	std::string default_answer="export1";
	answer =  fl_input("Name of the exported block of jobs? ", default_answer.c_str());
	std::string mydir(answer);
	mydir += "/";
	std::string command = "mkdir -p ExportJobs/" + mydir;
	int res = system(command.c_str());

	MetaDataTable MDexported;

	// Loop through all the Scheduled jobs and export them one-by-one
	int iexp =0;
	std::vector<std::string> find_pattern, replace_pattern;
	for (long int i = 0; i < pipeline.processList.size(); i++)
	{
		if (pipeline.processList[i].status == PROC_SCHEDULED_NEW ||
				pipeline.processList[i].status == PROC_SCHEDULED_CONT)
		{
			iexp++;
			if (pipeline.processList[i].alias != "None")
			{
				fl_message("ERROR: aliases are not allowed on Scheduled jobs that are to be exported! Make sure all scheduled jobs are made with unaliases names.");
				return;
			}

			// A general name for the exported job:
			FileName expname = pipeline.processList[i].name;
			expname = expname.beforeFirstOf("/") + "/exp"+integerToString(iexp, 3)+"/";
			//std::cerr << " pipeline.processList[i].name= " << pipeline.processList[i].name << " expname= " << expname << std::endl;
			find_pattern.push_back(pipeline.processList[i].name);
			replace_pattern.push_back(expname);

			MDexported.addObject();
			MDexported.setValue(EMDL_PIPELINE_PROCESS_NAME, expname);

			// Copy the run.job, the note.txt and the job_pipeline.star and replace patterns
			replaceFilesForImportExportOfScheduledJobs(pipeline.processList[i].name, "ExportJobs/" + mydir + expname, find_pattern, replace_pattern);
		}
	}

	MDexported.write("ExportJobs/" + mydir + "exported.star");

}

void RelionMainWindow::cb_import_jobs(Fl_Widget* o, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_import_jobs_i();
}


void RelionMainWindow::cb_import_jobs_i()
{

	// Get the directory with the Exported jobs
	std::string fn_dir = ".";
	std::string fn_filter = "Export STAR file (exported.star)";
	Fl_File_Chooser chooser(fn_dir.c_str(),  fn_filter.c_str(), Fl_File_Chooser::SINGLE, "Choose pipeline STAR file to import");
	chooser.show();
	// Block until user picks something.
	while(chooser.shown())
		{ Fl::wait(); }

	// User hit cancel?
	if ( chooser.value() == NULL )
		return;

	FileName fn_export(chooser.value());
	FileName fn_export_dir = fn_export.beforeLastOf("/")+"/";

	//FileName fn_dir_export = fn_export.beforeLastOf("/")+"/";
	MetaDataTable MDexported;
	MDexported.read(fn_export);

	// Read in existing pipeline, in case some other window had changed it
	pipeline.read(DO_LOCK);

	std::vector<std::string> find_pattern, replace_pattern;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDexported)
	{
		FileName expname;
		MDexported.getValue(EMDL_PIPELINE_PROCESS_NAME, expname);
		find_pattern.push_back(expname);
		// Make a new name for this job
		FileName newname = expname.beforeFirstOf("/")+"/job"+integerToString(pipeline.job_counter, 3)+"/";
		//std::cerr << " expname= " << expname << " newname= " << newname << std::endl;
		replace_pattern.push_back(newname);
		replaceFilesForImportExportOfScheduledJobs(fn_export_dir + expname, newname, find_pattern, replace_pattern);
		// Import the job into the pipeline
	    pipeline.importPipeline(newname+"job");
	    pipeline.job_counter++;
	}

	// Write the new pipeline to disk
	fillRunningJobLists();
	pipeline.write(DO_LOCK);

}

// Re-order running and finished job lists
void RelionMainWindow::cb_order_jobs_alphabetically(Fl_Widget* o, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    do_order_alphabetically = true;
    T->fillRunningJobLists();
}

// Re-order running and finished job lists
void RelionMainWindow::cb_order_jobs_chronologically(Fl_Widget* o, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    do_order_alphabetically = false;
    T->fillRunningJobLists();
}

// Empty-trash button call-back function
void RelionMainWindow::cb_empty_trash(Fl_Widget* o, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_empty_trash_i();
}

void RelionMainWindow::cb_empty_trash_i()
{
	std::string ask = "Are you sure you want to remove the entire Trash folder?";
	int proceed =  fl_choice("%s", "Don't empty trash", "Empty Trash", NULL, ask.c_str());
	if (proceed)
	{
		std::string command = "rm -rf Trash";
		std::cout << " Executing: " << command << std::endl;
		int res = system(command.c_str());
	}
}

void RelionMainWindow::cb_print_notes(Fl_Widget*, void* v)
{
  	 RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_print_notes_i();
}

void RelionMainWindow::cb_print_notes_i()
{
	std::ofstream  fh;
	FileName fn_tmp = pipeline.name + "_all_notes.txt";
	fh.open((fn_tmp).c_str(), std::ios::out);

	for (size_t i = 0; i < pipeline.processList.size(); i++)
	{
		FileName fn_note = pipeline.processList[i].name+"note.txt";
		fh << " ################################################################ " << std::endl;
		fh << " # Job= " << pipeline.processList[i].name;
		if (pipeline.processList[i].alias != "None")
			fh <<" alias: " << pipeline.processList[i].alias;
		fh	<< std::endl;
		if (exists(fn_note))
		{
			std::ifstream in(fn_note.data(), std::ios_base::in);
			std::string line;
			if (in.fail())
				REPORT_ERROR( (std::string) "ERROR: cannot read file " + fn_note);
    	    in.seekg(0);
    	    while (getline(in, line, '\n'))
    	    {
    	    	fh << line << std::endl;
    	    }
			in.close();
		}
	}
	fh.close();

	fl_message("Done writing all notes into file: %s" , fn_tmp.c_str());

}

void RelionMainWindow::cb_remake_nodesdir(Fl_Widget*, void* v)
{
  	 RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_remake_nodesdir_i();
}

void RelionMainWindow::cb_remake_nodesdir_i()
{
	pipeline.remakeNodeDirectory();
}

void RelionMainWindow::cb_reread_pipeline(Fl_Widget*, void* v)
{
  	 RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_reread_pipeline_i();
}

void RelionMainWindow::cb_reread_pipeline_i()
{
	pipeline.read(DO_LOCK);
	// With the locking system, each read needs to be followed soon with a write
	pipeline.write(DO_LOCK);
}


void RelionMainWindow::cb_reactivate_runbutton(Fl_Widget* o, void* v)
{

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_reactivate_runbutton_i();
}

void RelionMainWindow::cb_reactivate_runbutton_i()
{
	run_button->activate();
}

void RelionMainWindow::cb_show_initial_screen(Fl_Widget* o, void* v)
{

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_show_initial_screen_i();
}

void RelionMainWindow::cb_show_initial_screen_i()
{
	show_initial_screen = true;
	cb_select_browsegroup_i();
}

void RelionMainWindow::cb_start_pipeliner(Fl_Widget* o, void* v)
{

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_start_pipeliner_i();
}

void RelionMainWindow::cb_start_pipeliner_i()
{

	std::vector<FileName> job_names;
	std::vector<long int> job_ids;

	for (long int ii =0; ii < scheduled_processes.size(); ii++)
	{
		long int id = scheduled_processes[ii];
		job_ids.push_back(id);
		job_names.push_back(pipeline.processList[id].name);
	}
	std::vector<long int> my_scheduled_processes = scheduled_processes;
	SchedulerWindow* w = new SchedulerWindow(400, 300, "Select which jobs to execute");
	w->fill(pipeline.name, job_names, job_ids);

}

void RelionMainWindow::cb_stop_pipeliner(Fl_Widget* o, void* v)
{

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_stop_pipeliner_i();
}

void RelionMainWindow::cb_stop_pipeliner_i()
{
	std::string fn_filter = "Pipeline scheduled file (RUNNING_PIPELINER_" + pipeline.name + "_*)";
	Fl_File_Chooser chooser(".",  fn_filter.c_str(), Fl_File_Chooser::SINGLE, "Choose which scheduler to stop");
	chooser.show();
	// Block until user picks something.
	while(chooser.shown())
		{ Fl::wait(); }

	// User hit cancel?
	if ( chooser.value() == NULL )
		return;

	FileName fn_del(chooser.value());
	std::cout <<" Deleting file : " << fn_del<< std::endl;
	std::remove(fn_del.c_str());
}

void RelionMainWindow::cb_about(Fl_Widget* o, void* v)
{

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_about_i();
}

void RelionMainWindow::cb_about_i()
{
	ShowHelpText *help = new ShowHelpText("\
RELION is written by Sjors Scheres at the MRC Laboratory of Molecular Biology (scheres@mrc-lmb.cam.ac.uk).\n \
\n\
Note that RELION is completely free, open-source software. You can redistribute it and/or modify it for your own purposes, but please do make sure \
the contribution of Sjors Scheres is acknowledged appropriately. In order to maintain an overview of existing versions, he would also appreciate being \
notified of any redistribution of (modified versions of) the code. \n \n \n \
If RELION is useful in your work, please cite us. Relevant papers are:\n \n \
 * General Bayesian approach (and first mention of RELION): \n \
     Scheres (2012) J. Mol. Biol. (PMID: 22100448)	 \n \n\
 * RELION implementation details and the 3D auto-refine procedure: \n \
     Scheres (2012) J. Struct. Biol. (PMID: 23000701)	 \n \n\
 * Gold-standard FSC and the relevance of the 0.143 criterion: \n \
     Scheres & Chen (2012) Nat. Meth. (PMID: 22842542)	 \n \n\
 * Movie-processing procedure: \n \
     Bai et al. (2013) eLife (PMID: 23427024 )	 \n \n\
 * Correction of mask effects on the FSC curve by randomised phases: \n \
     Chen et al. (2013) Ultramicroscopy (PMID: 23872039)	 \n \n\
 * Particle-polishing: \n \
     Scheres (2014) eLife (PMID: 25122622)	 \n \n\
 * Auto-picking : \n \
     Scheres (2014) J. Struct. Biol. (PMID: 25486611) \n \n \
 * Sub-tomogram averaging : \n \
     Bharat et al. (2015) Structure (PMID: 26256537) \n \n \
Please also cite the following EXTERNAL programs: \n \n \
* CTFFIND for CTF-estimation: \n \
    Mindell & Grigorieff (2003) J. Mol. Biol. (PMID: 12781660) \n \n\
* MOTIONCORR for beam-induced motion correction: \n \
    Li et al (2013) Nat. Methods (PMID: 23644547) \n \n\
* ResMap for local-resolution estimation:  \n\
    Kucukelbir et al. (2014) Nat. Meth. (PMID: 24213166) \n \n\
* CTFFIND4 for CTF-estimation: \n \
    Rohou & Grigorieff (2015) J. Struct. Biol. (PMID: 26278980) \n \n\
* Gctf for CTF-estimation: \n \
    Zhang (2016) J. Struct. Biol. (PMID: 2659270) \n \n\
* Postscript plots are made using CPlot2D from  www.amzsaki.com\n \
");
}


void RelionMainWindow::cb_quit(Fl_Widget* o, void* v)
{
	RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_quit_i();
}

void RelionMainWindow::cb_quit_i()
{
	exit(0);
}
