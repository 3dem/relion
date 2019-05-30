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
#include "src/gui_background.xpm"

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

// The StdOutDisplay allows looking at the entire stdout or stderr file
int StdOutDisplay::handle(int ev)
{

	if (ev==FL_PUSH && Fl::event_clicks())
	{
		// double-click
		if (Fl::event_clicks())
		{
			if (current_job < 0) return 0;
			current_browse_directory = pipeline.processList[current_job].name;
			FileName fn = current_browse_directory + fn_file;
			std::string command;
			if (exists(fn))
			{
				if (fn_file == "run.out")
				{
					if (maingui_do_read_only)
					{
						NoteEditorWindow* w = new NoteEditorWindow(800, 400, fn.c_str(), fn.c_str(), false); // false means dont_allow_save
						w->show();
						return 1;
					}
					else
					{
						std::string command = "awk -F\"\r\" '{if (NF>1) {print $NF} else {print}}' < " + fn + " > .gui_tmpstd";
						int res = system(command.c_str());
						NoteEditorWindow* w = new NoteEditorWindow(800, 400, fn.c_str(), ".gui_tmpstd", false); //false means dont_allow_save, as its temp file anyway
						w->show();
						return 1;
					}
				}
				else
				{
					NoteEditorWindow* w = new NoteEditorWindow(800, 400, fn.c_str(), fn, true); // true means allow_save, this is useful to remove past errors
					w->show();
					return 1;
				}

			}
		} // end if double click
	} // end if FL_PUSH

	return 0;
}

int SchedulerWindow::fill(FileName _pipeline_name, std::vector<FileName> _scheduled_jobs)
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
    	my_jobs.push_back(_scheduled_jobs[ijob]);
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
    wait_before = new Fl_Input(xcol, current_y, 100, ystep-8, "Wait this many minutes before starting?");
    current_y += ystep;
    repeat = new Fl_Input(xcol, current_y, 100, ystep-8, "Run the jobs how many times?");
    current_y += ystep;
    wait = new Fl_Input(xcol, current_y, 100, ystep-8, "Wait at least in between (in minutes)?");
    current_y += ystep;
    wait_after = new Fl_Input(xcol, current_y, 100, ystep-8, "Wait at least after each job (in seconds)?");
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
    wait_before->value("0");
	wait_before->color(GUI_INPUT_COLOR);
	wait_before->textsize(ENTRY_FONTSIZE);
	wait_before->labelsize(ENTRY_FONTSIZE);
    wait_after->value("10");
	wait_after->color(GUI_INPUT_COLOR);
	wait_after->textsize(ENTRY_FONTSIZE);
	wait_after->labelsize(ENTRY_FONTSIZE);

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
		fl_message("%s",msg.c_str());
	}
	else
	{
		// Make a string with all job-ids to process
		std::string jobids="\"";
		for (int ijob = 0; ijob < my_jobs.size(); ijob++)
		{
			if (check_buttons[ijob]->value())
				jobids += my_jobs[ijob] + " ";
		}
		jobids += "\"";

		std::string myrepeat(repeat->value());
		std::string mywait(wait->value());
		std::string mywait_before(wait_before->value());
		std::string mywait_after(wait_after->value());

		std::string command = "relion_pipeliner --pipeline " + pipeline_name;
		command += " --schedule " + fn_sched;
		command += " --repeat " + myrepeat;
		command += " --min_wait " + mywait;
		command += " --min_wait_before " + mywait_before;
		command += " --sec_wait_after " + mywait_after;
		command += " --RunJobs " + jobids;
		// Run this in the background, so control returns to the window
		command += " &";
		int res = system(command.c_str());
		std::cout << " Launching: " << command << std::endl;
		std::cout << " Stop execution of this set of scheduled jobs by deleting file: " << fn_check << std::endl;

	}

}
/*
int SchedulerAddVariableOperatorWindow::fill(bool is_variable, bool is_add)
{
	std::cerr << "in fill" << std::endl;

	//color(GUI_BACKGROUND_COLOR);
    int current_y = 2, max_y = 2;
    int ystep = 35;

    int xcol1 = w()-300;
    int xcol2 = w()-120;

    if (is_variable)
    {
    	scheduler_add_variable_name = new Fl_Input(xcol1, current_y, 100, ystep-8, "Name: ");
    	scheduler_add_variable_value = new Fl_Input(xcol2, current_y, 100, ystep-8, "Value:");
		current_y += ystep;
    }
    else
    {
    	std::cerr << "todo: set/add scheduler_operators" << std::endl;
    }

	// Button to execute
	std::string mybuttontext = (is_add) ? "Add" : "Set";
    Fl_Button *execute_button = new Fl_Button(w()-100, current_y, 80, 30, mybuttontext.c_str());
	execute_button->color(GUI_RUNBUTTON_COLOR);
	execute_button->labelsize(12);
	execute_button->callback( cb_add, this);

	// Button to cancel
	Fl_Button *cancel_button = new Fl_Button(w()-200, current_y, 80, 30, "Cancel");
	cancel_button->color(GUI_RUNBUTTON_COLOR);
	cancel_button->labelsize(12);
	cancel_button->callback( cb_cancel, this);

	resizable(*this);
	show();

	return Fl::run();

}
void SchedulerAddVariableOperatorWindow::cb_cancel(Fl_Widget*, void* v)
{
    SchedulerAddVariableOperatorWindow* T=(SchedulerAddVariableOperatorWindow*)v;
    scheduler_add_variable_name->value("");
    scheduler_add_variable_value->value("");
    T->hide();
}

void SchedulerAddVariableOperatorWindow::cb_add(Fl_Widget*, void* v)
{
    SchedulerAddVariableOperatorWindow* T=(SchedulerAddVariableOperatorWindow*)v;
    T->cb_add_i();
    T->hide();
}

void SchedulerAddVariableOperatorWindow::cb_add_i()
{
	std::string myname = scheduler_add_variable_name->value();
	std::string myvalue = scheduler_add_variable_value->value();

	// Check wether a variable with this name already exists
	if (!(myvalue == "true" || myvalue == "True" || myvalue == "false" || myvalue == "False"))
	{
		if (schedule.isBooleanVariable(myname))
		{
			fl_message("ERROR: a boolean variable with this name already exists");
			return;
		}
	}
	else
	{
		float floatval;
		if (sscanf(myvalue.c_str(), "%f", &floatval)) // is this a number?
		{
			if (schedule.isFloatVariable(myname))
			{
				fl_message("ERROR: a float variable with this name already exists");
				return;
			}
		}
		else
		{
			if (schedule.isStringVariable(myname))
			{
				fl_message("ERROR: a string variable with this name already exists");
				return;
			}
		}
	}

	// Now set the global variables that will be used back in the main GUI

}
*/

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
	resizable(*this);

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


GuiMainWindow::GuiMainWindow(int w, int h, const char* title, FileName fn_pipe, FileName fn_sched, int _update_every_sec, int _exit_after_sec, bool _do_read_only):Fl_Window(w,h,title)
{

	// Set initial Timer
	tickTimeLastChanged();

	// Setup read_only
	maingui_do_read_only = _do_read_only;
	pipeline.do_read_only = _do_read_only;
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

	// Initial screen picture with some explanation on how to use the GUI
        //image_box = new Fl_Box(WCOL0-8, 0 ,w-WCOL0, h-35); // widget that will contain image
	image_box = new Fl_Box(WCOL0-8, 50 ,w-WCOL0, h-120); // widget that will contain image
	xpm_image = new Fl_Pixmap(gui_background);
	image_box->image(xpm_image); // attach xpm image to box

	background_grp->end();

	// read in schedule if it exists, otherwise just initialise schedule with its name
	if (fn_sched != "")
	{
		show_scheduler = true;
		schedule.setName("Schedules/"+fn_sched+"/");
		pipeline.name = "Schedules/"+fn_sched+"/schedule";
		if (exists(schedule.name+"schedule.star"))
		{
			schedule.read();
			pipeline.name = "Schedules/"+fn_sched+"/schedule";
		}
		else
		{
			std::string command = "mkdir -p Schedules/" + fn_sched;
			int res = system(command.c_str());
			pipeline.write(); // empty write
		}

	}
	else
	{
		// Read in the pipeline STAR file if it exists
		pipeline.name = fn_pipe;
	}
	if (exists(fn_pipe + "_pipeline.star"))
	{
		std::string lock_message = "mainGUI constructor";
		pipeline.read(DO_LOCK, lock_message);
		// With the locking system, each read needs to be followed soon with a write
		pipeline.write(DO_LOCK);
	}
	else
	{
		pipeline.write();
	}

 	color(GUI_BACKGROUND_COLOR);
    menubar = new Fl_Menu_Bar(-3, 0, WCOL0-7, MENUHEIGHT);
	menubar->add("File/Re-read pipeline",  FL_ALT+'r', cb_reread_pipeline, this);
	menubar->add("File/Edit project note",  FL_ALT+'e', cb_edit_project_note, this);
	if (!maingui_do_read_only)
		menubar->add("File/Print all notes",  FL_ALT+'p', cb_print_notes, this);
	if (!maingui_do_read_only)
		menubar->add("File/Remake .Nodes\\/",  FL_ALT+'n', cb_remake_nodesdir, this);
	menubar->add("File/Display",  FL_ALT+'d', cb_display, this);
	menubar->add("File/_Overwrite continue",  FL_ALT+'o', cb_toggle_overwrite_continue, this);
	menubar->add("File/_Show initial screen",  FL_ALT+'z', cb_show_initial_screen, this);
	if (!maingui_do_read_only)
		menubar->add("File/_Empty trash",  FL_ALT+'t', cb_empty_trash, this);
    menubar->add("File/About", 0, cb_about, this);
    menubar->add("File/Quit", FL_ALT+'q', cb_quit, this);
	if (!maingui_do_read_only)
	{	menubar->add("Jobs/Save job settings",  FL_ALT+'s', cb_save, this);
		menubar->add("Jobs/_Load job settings",  FL_ALT+'l', cb_load, this);
	}
	menubar->add("Jobs/Order alphabetically",  FL_ALT+'a', cb_order_jobs_alphabetically, this);
	menubar->add("Jobs/_Order chronologically",  FL_ALT+'c', cb_order_jobs_chronologically, this);
	if (!maingui_do_read_only)
	{	menubar->add("Jobs/_Undelete job(s)",  FL_ALT+'u', cb_undelete_job, this);
		menubar->add("Jobs/Export scheduled job(s)",  FL_ALT+'x', cb_export_jobs, this);
		menubar->add("Jobs/_Import scheduled job(s)",  FL_ALT+'i', cb_import_jobs, this);
		menubar->add("Jobs/Gently clean all jobs",  FL_ALT+'g', cb_gently_clean_all_jobs, this);
		menubar->add("Jobs/Harshly clean all jobs",  FL_ALT+'h', cb_harshly_clean_all_jobs, this);
		menubar->add("Autorun/Run scheduled jobs", 0, cb_start_pipeliner, this);
		menubar->add("Autorun/Stop running scheduled jobs", 0, cb_stop_pipeliner, this);
	}
    current_y = MENUHEIGHT + 10;

    // Add run buttons on the menubar as well
	print_CL_button = new Fl_Button(GUIWIDTH - 330, h-90, 100, 32, "Print command");
	print_CL_button->color(GUI_RUNBUTTON_COLOR);
	print_CL_button->labelsize(12);
	print_CL_button->callback( cb_print_cl, this);


	// Fill browser in the right order
	browser = new Fl_Hold_Browser(10,MENUHEIGHT+5,WCOL0-20,h-MENUHEIGHT-60);
    browser->textsize(RLN_FONTSIZE-1);
    current_job = -1;

    int i = 0;
    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browser->add("Import");
    gui_jobwindows[i] = new JobWindow();
    gui_jobwindows[i]->initialise(PROC_IMPORT);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browser->add("Motion correction");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_MOTIONCORR);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
    browser->add("CTF estimation");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_CTFFIND);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("Manual picking");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_MANUALPICK);
	browse_grp[i]->end();
	i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("Auto-picking");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_AUTOPICK);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("Particle extraction");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_EXTRACT);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("Subset selection");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_CLASSSELECT);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("2D classification");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_2DCLASS);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("3D initial model");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_INIMODEL);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("3D classification");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_3DCLASS);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("3D auto-refine");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_3DAUTO);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("3D multi-body");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_MULTIBODY);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("CTF refinement");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_CTFREFINE);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("Bayesian polishing");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_MOTIONREFINE);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("Mask creation");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_MASKCREATE);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("Join star files");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_JOINSTAR);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("Particle subtraction");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_SUBTRACT);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("Post-processing");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_POST);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("Local resolution");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_RESMAP);
    browse_grp[i]->end();
    i++;

    browse_grp[i] = new Fl_Group(WCOL0, 2, 550, 615-MENUHEIGHT);
	browser->add("External");
	gui_jobwindows[i] = new JobWindow();
	gui_jobwindows[i]->initialise(PROC_EXTERNAL);
    browse_grp[i]->end();

    browser->callback(cb_select_browsegroup);
    browser->end();
    browser->select(1); // just start from the beginning


    // 0) Both pipeliner and scheduler GUI

	Fl_Text_Buffer *textbuff3 = new Fl_Text_Buffer();
	Fl_Text_Buffer *textbuff4 = new Fl_Text_Buffer();
	Fl_Text_Buffer *textbuff5 = new Fl_Text_Buffer();
	textbuff3->text("Scheduled jobs");
	textbuff4->text("Input to this job");
	textbuff5->text("Output from this job");
	Fl_Text_Display* textdisp3, *textdisp4, *textdisp5;
	if (show_scheduler)
	{
		textdisp3 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_EXT_START, JOBCOLWIDTH, 25);
		textdisp4 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_EXT_START + JOBHALFHEIGHT +25, JOBCOLWIDTH, 25);
		textdisp5 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_EXT_START + 1.5*JOBHALFHEIGHT + 50, JOBCOLWIDTH, 25);
		scheduled_job_browser = new Fl_Select_Browser(XJOBCOL3, GUIHEIGHT_EXT_START + 25 , JOBCOLWIDTH, JOBHALFHEIGHT);
		input_job_browser     = new Fl_Select_Browser(XJOBCOL3, GUIHEIGHT_EXT_START + JOBHALFHEIGHT + 50, JOBCOLWIDTH, 0.5*JOBHALFHEIGHT);
		output_job_browser    = new Fl_Select_Browser(XJOBCOL3, GUIHEIGHT_EXT_START + 1.5*JOBHALFHEIGHT + 75, JOBCOLWIDTH, 0.5*JOBHALFHEIGHT);

	}
	else
	{
		textdisp3 = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_EXT_START2 + JOBHALFHEIGHT + 25, JOBCOLWIDTH, 25);
		textdisp4 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_EXT_START2, JOBCOLWIDTH, 25);
		textdisp5 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_EXT_START2 + JOBHALFHEIGHT + 25, JOBCOLWIDTH, 25);
		scheduled_job_browser = new Fl_Select_Browser(XJOBCOL2, GUIHEIGHT_EXT_START2 + 25 + JOBHALFHEIGHT + 25, JOBCOLWIDTH, JOBHALFHEIGHT);
		input_job_browser    = new Fl_Select_Browser(XJOBCOL3,  GUIHEIGHT_EXT_START2 + 25, JOBCOLWIDTH, JOBHALFHEIGHT);
		output_job_browser   = new Fl_Select_Browser(XJOBCOL3,  GUIHEIGHT_EXT_START2 + 25 + JOBHALFHEIGHT + 25, JOBCOLWIDTH, JOBHALFHEIGHT);
	}

	textdisp3->buffer(textbuff3);
	textdisp4->buffer(textbuff4);
	textdisp5->buffer(textbuff5);
	textdisp3->color(GUI_BACKGROUND_COLOR);
	textdisp4->color(GUI_BACKGROUND_COLOR);
	textdisp5->color(GUI_BACKGROUND_COLOR);

	scheduled_job_browser->callback(cb_select_scheduled_job);
	input_job_browser->callback(cb_select_input_job);
	output_job_browser->callback(cb_select_output_job);
	scheduled_job_browser->textsize(RLN_FONTSIZE-1);
	input_job_browser->textsize(RLN_FONTSIZE-1);
	output_job_browser->textsize(RLN_FONTSIZE-1);

	scheduled_job_browser->end();
	input_job_browser->end();
	output_job_browser->end();

    // A) Pipeliner part of the GUI
	pipeliner_grp = new Fl_Group(0, 0, 2*w, 2*h);
	pipeliner_grp->begin();

	run_button = new Fl_Button(GUIWIDTH - 110 , h-90, 100, 32, "Run!");
	run_button->color(GUI_RUNBUTTON_COLOR);
	run_button->labelfont(FL_ITALIC);
	run_button->labelsize(14);
	run_button->callback( cb_run, this);
	if (maingui_do_read_only)
		run_button->deactivate();

	schedule_button = new Fl_Button(GUIWIDTH - 220 , h-90, 100, 32, "Schedule");
	schedule_button->color(GUI_RUNBUTTON_COLOR);
	schedule_button->labelfont(FL_ITALIC);
	schedule_button->labelsize(14);
	schedule_button->callback( cb_schedule, this);
	if (maingui_do_read_only)
		schedule_button->deactivate();


	menubar2 = new Fl_Menu_Bar(XJOBCOL1, GUIHEIGHT_EXT_START, 100, MENUHEIGHT);
	menubar2->color(GUI_BUTTON_COLOR);
	menubar2->add("Job actions/Edit Note", 0, cb_edit_note, this);
	if (!maingui_do_read_only)
	{
		menubar2->add("Job actions/Alias", 0, cb_set_alias, this);
		menubar2->add("Job actions/Abort running", 0, cb_abort, this);
		menubar2->add("Job actions/Mark as finished", 0, cb_mark_as_finished, this);
		menubar2->add("Job actions/Make flowchart", 0, cb_make_flowchart, this);
		menubar2->add("Job actions/Gentle clean", 0, cb_gentle_cleanup, this);
		menubar2->add("Job actions/Harsh clean", 0, cb_harsh_cleanup, this);
		menubar2->add("Job actions/Delete", 0, cb_delete, this);
	}

	// Fl_input with the alias of the new job (or the name of an existing one)
	alias_current_job = new Fl_Input(XJOBCOL2-50 , GUIHEIGHT_EXT_START+3, JOBCOLWIDTH, MENUHEIGHT-6, "Current job:");

	// Left-hand side browsers for input/output nodes and processes
	display_io_node  = new Fl_Choice(XJOBCOL3, GUIHEIGHT_EXT_START+3, 250, MENUHEIGHT-6);
	display_io_node->label("Display:");
	display_io_node->color(GUI_BUTTON_COLOR);
	display_io_node->callback(cb_display_io_node, this);

	// Add browsers for finished and running jobs
	Fl_Text_Buffer *textbuff1 = new Fl_Text_Buffer();
	textbuff1->text("Finished jobs");
	Fl_Text_Display* textdisp1 = new Fl_Text_Display(XJOBCOL1, GUIHEIGHT_EXT_START2, JOBCOLWIDTH, 25);
	textdisp1->buffer(textbuff1);
	textdisp1->color(GUI_BACKGROUND_COLOR);
	finished_job_browser  = new Fl_Select_Browser(XJOBCOL1, GUIHEIGHT_EXT_START2 + 25, JOBCOLWIDTH, JOBHEIGHT+25);
	finished_job_browser->callback(cb_select_finished_job);
	finished_job_browser->textsize(RLN_FONTSIZE-1);
	finished_job_browser->end();

	Fl_Text_Buffer *textbuff2 = new Fl_Text_Buffer();
	textbuff2->text("Running jobs");
	Fl_Text_Display* textdisp2 = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_EXT_START2, JOBCOLWIDTH, 25);
	textdisp2->buffer(textbuff2);
	textdisp2->color(GUI_BACKGROUND_COLOR);
	running_job_browser   = new Fl_Select_Browser(XJOBCOL2, GUIHEIGHT_EXT_START2 + 25, JOBCOLWIDTH, JOBHALFHEIGHT);
	running_job_browser->callback(cb_select_running_job);
	running_job_browser->textsize(RLN_FONTSIZE-1);
	running_job_browser->end();

	// Display stdout and stderr of jobs
	textbuff_stdout = new Fl_Text_Buffer();
	textbuff_stderr = new Fl_Text_Buffer();
	// Disable warning message about UTF-8 transcoding
	textbuff_stdout->transcoding_warning_action=NULL;
	textbuff_stderr->transcoding_warning_action=NULL;
	disp_stdout = new StdOutDisplay(XJOBCOL1, GUIHEIGHT_EXT_START2 + JOBHEIGHT + STDOUT_Y-5, w-20, 105);
	disp_stderr = new StdOutDisplay(XJOBCOL1, GUIHEIGHT_EXT_START2 + JOBHEIGHT + STDERR_Y-5, w-20, 50);
	disp_stdout->fn_file = "run.out";
	disp_stderr->fn_file = "run.err";
	textbuff_stdout->text("stdout will go here; double-click this window to open stdout in a separate window");
	textbuff_stderr->text("stderr will go here; double-click this window to open stderr in a separate window");
	disp_stdout->buffer(textbuff_stdout);
	disp_stderr->buffer(textbuff_stderr);
	disp_stderr->textcolor(FL_RED);
	disp_stdout->textsize(RLN_FONTSIZE-1);
	disp_stderr->textsize(RLN_FONTSIZE-1);
	disp_stdout->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS,0);
	disp_stderr->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS,0);
	disp_stdout->scrollbar_width(0);
	disp_stderr->scrollbar_width(0);

	pipeliner_grp->end();

	// B) Scheduler part of the GUI
	scheduler_grp = new Fl_Group(0, 0, 4*w, 4*h);
	scheduler_grp->begin();


	scheduler_job_name = new Fl_Input(GUIWIDTH - 500, h-90, 150, 32, "Name:");
	scheduler_job_name->color(GUI_INPUT_COLOR);

	// Select one of three modes for adding a new job
	job_mode  = new Fl_Choice(GUIWIDTH - 220 , h-90, 100, 32);
	job_mode->label("");
	job_mode->color(GUI_RUNBUTTON_COLOR);
	//job_mode->callback(cb_display_io_node, this);
	job_mode->menu(job_mode_options);
	// TODO: fill options for this choice!

	add_job_button = new Fl_Button(GUIWIDTH - 110 , h-90, 100, 32, "Add job");
	add_job_button->color(GUI_RUNBUTTON_COLOR);
	add_job_button->labelfont(FL_ITALIC);
	add_job_button->labelsize(14);
	add_job_button->callback( cb_scheduler_add_job, this);

	// Add browsers for variable and nodes
	Fl_Text_Buffer *textbuffvar = new Fl_Text_Buffer();
	textbuffvar->text("Variables");
	Fl_Text_Display* textdispvar = new Fl_Text_Display(XJOBCOL1, GUIHEIGHT_EXT_START, JOBCOLWIDTH, 25);
	textdispvar->buffer(textbuffvar);
	textdispvar->color(GUI_BACKGROUND_COLOR);
	scheduler_set_variable_name = new Fl_Input(XJOBCOL1, GUIHEIGHT_EXT_START+25, JOBCOLWIDTH/2, 25);
    scheduler_set_variable_name->color(GUI_INPUT_COLOR);
    scheduler_set_variable_value = new Fl_Input(XJOBCOL1+JOBCOLWIDTH/2, GUIHEIGHT_EXT_START+25, JOBCOLWIDTH/2, 25);
    scheduler_set_variable_value->color(GUI_INPUT_COLOR);
	delete_scheduler_variable_button = new Fl_Button(XJOBCOL1+JOBCOLWIDTH-105, GUIHEIGHT_EXT_START, 50, 25);
	delete_scheduler_variable_button->color(GUI_RUNBUTTON_COLOR);
	delete_scheduler_variable_button->labelfont(FL_ITALIC);
	delete_scheduler_variable_button->labelsize(RLN_FONTSIZE);
	delete_scheduler_variable_button->label("Del");
	delete_scheduler_variable_button->callback(cb_delete_scheduler_variable, this);
	set_scheduler_variable_button = new Fl_Button(XJOBCOL1+JOBCOLWIDTH-50, GUIHEIGHT_EXT_START, 50, 25);
	set_scheduler_variable_button->color(GUI_RUNBUTTON_COLOR);
	set_scheduler_variable_button->labelfont(FL_ITALIC);
	set_scheduler_variable_button->labelsize(RLN_FONTSIZE);
	set_scheduler_variable_button->label("Set");
	set_scheduler_variable_button->callback(cb_set_scheduler_variable, this);
	scheduler_variable_browser  = new Fl_Hold_Browser(XJOBCOL1, GUIHEIGHT_EXT_START + 50, JOBCOLWIDTH, JOBHEIGHT-16-25);
	scheduler_variable_browser->callback(cb_select_scheduler_variable);
	scheduler_variable_browser->textsize(RLN_FONTSIZE-1);
	scheduler_variable_browser->end();

	Fl_Text_Buffer *textbuffnode = new Fl_Text_Buffer();
	textbuffnode->text("Operators");
	Fl_Text_Display* textdispnode = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_EXT_START, JOBCOLWIDTH, 25);
	textdispnode->buffer(textbuffnode);
	textdispnode->color(GUI_BACKGROUND_COLOR);
	delete_scheduler_operator_button = new Fl_Button(XJOBCOL2+JOBCOLWIDTH-105, GUIHEIGHT_EXT_START, 50, 25);
	delete_scheduler_operator_button->color(GUI_RUNBUTTON_COLOR);
	delete_scheduler_operator_button->labelfont(FL_ITALIC);
	delete_scheduler_operator_button->labelsize(RLN_FONTSIZE);
	delete_scheduler_operator_button->label("Del");
	delete_scheduler_operator_button->callback( cb_delete_scheduler_operator, this);
	set_scheduler_operator_button = new Fl_Button(XJOBCOL2+JOBCOLWIDTH-50, GUIHEIGHT_EXT_START, 50, 25);
	set_scheduler_operator_button->color(GUI_RUNBUTTON_COLOR);
	set_scheduler_operator_button->labelfont(FL_ITALIC);
	set_scheduler_operator_button->labelsize(RLN_FONTSIZE);
	set_scheduler_operator_button->label("Set");
	set_scheduler_operator_button->callback( cb_set_scheduler_operator, this);
	scheduler_operator_browser  = new Fl_Hold_Browser(XJOBCOL2, GUIHEIGHT_EXT_START + 25, JOBCOLWIDTH, JOBHEIGHT-16);
	scheduler_operator_browser->callback(cb_select_scheduler_node);
	scheduler_operator_browser->textsize(RLN_FONTSIZE-1);
	scheduler_operator_browser->end();

	scheduler_grp->end();




	if (show_scheduler)
	{
		pipeliner_grp->hide();
		scheduler_grp->show();
		fillSchedulerNodesAndVariables();
	}
	else
	{
		scheduler_grp->hide();
		pipeliner_grp->show();
	}

	// Fill the actual browsers
	fillRunningJobLists();

	// Mechanism to update stdout and stderr continuously and also update the JobLists
	// Also exit the GUI if it has been idle for too long
	update_every_sec = _update_every_sec;
	exit_after_sec = (float)_exit_after_sec;
	if (update_every_sec > 0)
		Fl::add_timeout(update_every_sec, Gui_Timer_CB, (void*)this);

    cb_show_initial_screen_i();

    // Set and activate current selection from side-browser
	cb_select_browsegroup_i(true); // make default active; true is used to show_initial_screen
	is_main_continue = false; // default is a new run

}

static void Gui_Timer_CB(void *userdata)
{
	GuiMainWindow *o = (GuiMainWindow*)userdata;

	time_t now;
	time (&now);

	double dif = difftime (now, time_last_change);
    // If the GUI has been idle for too long, then exit
	if (dif > o->exit_after_sec)
    {
		std::cout << " The relion GUI has been idle for more than " << o->exit_after_sec << " seconds, exiting now... " << std::endl;
		exit(0);
    }

	// Update the stdout and stderr windows if we're currently pointing at a running job
	if (current_job >= 0 && pipeline.processList[current_job].status == PROC_RUNNING)
    	o->fillStdOutAndErr();

    // Check for job completion if the pipeline has been changed
	if (exists(PIPELINE_HAS_CHANGED))
		o->updateJobLists();

    // Refresh every so many seconds
    Fl::repeat_timeout(o->update_every_sec, Gui_Timer_CB, userdata);
}

void GuiMainWindow::clear()
{
	if (menubar != NULL)
	{
		delete menubar;
		menubar = NULL;
	}
	if (menubar2 != NULL)
	{
		delete menubar2;
		menubar2 = NULL;
	}
}

std::string GuiMainWindow::getJobNameForDisplay(Process &job)
{
	FileName result;
	FileName fn_pre, fn_jobnr, fn_post;

	if (show_scheduler)
	{
		result = job.name;
		result = result.afterFirstOf(schedule.name);
		return result.beforeLastOf("/");
	}
	else if (!decomposePipelineFileName(job.name, fn_pre, fn_jobnr, fn_post))
	{
		result = job.name;
	}
	else
	{
		std::string numberonly = (fn_jobnr.afterFirstOf("b")).beforeFirstOf("/");
		if (job.alias != "None")
			result = numberonly + ": " + job.alias;
		else
			result = numberonly + ": " + job.name;
	}

	return result;
}


// Update the content of the finished, running and scheduled job lists
void GuiMainWindow::fillRunningJobLists()
{
	// Go back to the same positions in the vertical scroll bars of the job lists after updating...
	int mypos_running = running_job_browser->position();
	int mypos_scheduled = scheduled_job_browser->position();
	int mypos_finished = finished_job_browser->position();
	int myhpos_running = running_job_browser->hposition();
	int myhpos_scheduled = scheduled_job_browser->hposition();
	int myhpos_finished = finished_job_browser->hposition();

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
			if (pipeline.processList[i].status == PROC_FINISHED_SUCCESS ||
				pipeline.processList[i].status == PROC_FINISHED_FAILURE ||
				pipeline.processList[i].status == PROC_FINISHED_ABORTED)
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
			if (pipeline.processList[i].status == PROC_FINISHED_SUCCESS ||
				pipeline.processList[i].status == PROC_FINISHED_FAILURE ||
				pipeline.processList[i].status == PROC_FINISHED_ABORTED)
			{
				finished_processes.push_back(i);
				finished_job_browser->add((getJobNameForDisplay(pipeline.processList[i])).c_str());
			}
		}
	}

	// For running and scheduled jobs search forwards, so that last jobs are at the bottom
	for (long int i = 0; i < pipeline.processList.size(); i++)
	{
		if (pipeline.processList[i].status == PROC_RUNNING)
		{
			running_processes.push_back(i);
			running_job_browser->add((getJobNameForDisplay(pipeline.processList[i])).c_str());
		}
		else if (pipeline.processList[i].status == PROC_SCHEDULED)
		{
			scheduled_processes.push_back(i);
			scheduled_job_browser->add((getJobNameForDisplay(pipeline.processList[i])).c_str());
		}
	}

	running_job_browser->position(mypos_running);
	scheduled_job_browser->position(mypos_scheduled);
	finished_job_browser->position(mypos_finished);
	running_job_browser->hposition(myhpos_running);
	scheduled_job_browser->hposition(myhpos_scheduled);
	finished_job_browser->hposition(myhpos_finished);
}

void GuiMainWindow::fillToAndFromJobLists()
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
					input_job_browser->add((getJobNameForDisplay(pipeline.processList[myproc])).c_str());
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
					output_job_browser->add((getJobNameForDisplay(pipeline.processList[myproc])).c_str());
				}
			}
		}
	}
}

// Update the content of the finished, running and scheduled job lists
void GuiMainWindow::fillSchedulerNodesAndVariables()
{
	// Go back to the same positions in the vertical scroll bars of the job lists after updating...
	int mypos_scheduler_variable = scheduler_variable_browser->position();
	int mypos_scheduler_node = scheduler_operator_browser->position();

    // Clear whatever was in there
	scheduler_variable_browser->clear();
	scheduler_operator_browser->clear();
	scheduler_floats_var_browser_idx.clear();
	scheduler_bools_var_browser_idx.clear();
	scheduler_strings_var_browser_idx.clear();
	scheduler_operator_browser_idx.clear();

	// Fill variables browser
	{
		std::map<std::string, SchedulerFloatVariable> scheduler_floats = schedule.getCurrentFloatVariables();
		std::map<std::string, SchedulerFloatVariable>::iterator it;
		int i = 0;
		for ( it = scheduler_floats.begin(); it != scheduler_floats.end(); it++ )
		{
			std::string mylabel = it->first +
					" = " + floatToString(it->second.value) +
					" (" + floatToString(it->second.original_value) + ")";
			scheduler_variable_browser->add(mylabel.c_str());
			scheduler_floats_var_browser_idx.push_back(i);
			scheduler_bools_var_browser_idx.push_back(-1);
			scheduler_strings_var_browser_idx.push_back(-1);
			i++;
		}
	}

	{
		std::map<std::string, SchedulerBooleanVariable> scheduler_bools = schedule.getCurrentBooleanVariables();
		std::map<std::string, SchedulerBooleanVariable>::iterator it;
		int i = 0;
		for ( it = scheduler_bools.begin(); it != scheduler_bools.end(); it++ )
		{
			std::string myval = (it->second.value) ? "True" : "False";
			std::string myorival = (it->second.original_value) ? "True" : "False";

			std::string mylabel = it->first +
					" = " + myval +
					" (" + myorival + ")";
			scheduler_floats_var_browser_idx.push_back(-1);
			scheduler_variable_browser->add(mylabel.c_str());
			scheduler_bools_var_browser_idx.push_back(i);
			scheduler_strings_var_browser_idx.push_back(-1);
			i++;
		}
	}

	{
		std::map<std::string, SchedulerStringVariable> scheduler_strings = schedule.getCurrentStringVariables();
		std::map<std::string, SchedulerStringVariable>::iterator it;
		int i = 0;
		for ( it = scheduler_strings.begin(); it != scheduler_strings.end(); it++ )
		{
			std::string mylabel = it->first +
					" = " + it->second.value +
					" (" + it->second.original_value + ")";
			scheduler_variable_browser->add(mylabel.c_str());
			scheduler_floats_var_browser_idx.push_back(-1);
			scheduler_bools_var_browser_idx.push_back(-i);
			scheduler_strings_var_browser_idx.push_back(i);
			i++;
		}
	}

	// Fill operator browser
	{
		std::map<std::string, SchedulerOperator> scheduler_operators = schedule.getCurrentOperators();
		std::map<std::string, SchedulerOperator>::iterator it;
		int i = 0;
		for ( it = scheduler_operators.begin(); it != scheduler_operators.end(); it++ )
		{
			std::string mylabel = it->first +
					" -> " + it->second.output;
			scheduler_operator_browser->add(mylabel.c_str());
			scheduler_operator_browser_idx.push_back(i);
			i++;
		}
	}

	scheduler_variable_browser->position(mypos_scheduler_variable);
	scheduler_operator_browser->position(mypos_scheduler_node);

}


void GuiMainWindow::fillStdOutAndErr()
{

	FileName fn_out = "";
	FileName fn_err = "";
	FileName fn_outtail, fn_errtail;
	if (current_job >= 0)
	{
		fn_out = pipeline.processList[current_job].name + "run.out";
		fn_err = pipeline.processList[current_job].name + "run.err";
		fn_outtail = pipeline.processList[current_job].name + ".run.out.tail";
		fn_errtail = pipeline.processList[current_job].name + ".run.err.tail";
	}

	if (exists(fn_out))
	{
		if (maingui_do_read_only)
		{
			int err = textbuff_stdout->loadfile(fn_out.c_str());
		}
		else
		{
			// Remove annoying carriage returns
			std::string command = "tail -n 6 < " + fn_out + " | awk -F\"\r\" '{if (NF>1) {print $NF} else {print}}' > " + fn_outtail;
			int res = system(command.c_str());
			std::ifstream in(fn_outtail.c_str(), std::ios_base::in);
			if (in.fail())
				REPORT_ERROR( (std::string) "MetaDataTable::read: File " + fn_outtail + " does not exists" );
			int err = textbuff_stdout->loadfile(fn_outtail.c_str());
			in.close();
		}
		// Scroll to the bottom
		disp_stdout->insert_position(textbuff_stdout->length()-1);
		disp_stdout->show_insert_position();
	}
	else
		textbuff_stdout->text("stdout will go here; double-click this window to open stdout in a separate window");

	if (exists(fn_err))
	{
		if (maingui_do_read_only)
		{
			int err = textbuff_stderr->loadfile(fn_err.c_str());
		}
		else
		{
			std::string command = "tail -3 " + fn_err + " > " + fn_errtail;
			int res = system(command.c_str());
			std::ifstream in(fn_errtail.c_str(), std::ios_base::in);
			if (in.fail())
				REPORT_ERROR( (std::string) "MetaDataTable::read: File " + fn_errtail + " does not exists" );
			int err = textbuff_stderr->loadfile(fn_errtail.c_str());
			in.close();
		}
		// Scroll to the bottom
		disp_stderr->insert_position(textbuff_stderr->length()-1);
		disp_stderr->show_insert_position();
	}
	else
		textbuff_stderr->text("stderr will go here; double-click this window to open stderr in a separate window");

}

void GuiMainWindow::tickTimeLastChanged()
{
	time(&time_last_change);
}

void GuiMainWindow::updateJobLists()
{
	pipeline.checkProcessCompletion();
	fillRunningJobLists();
	fillToAndFromJobLists();
}


void GuiMainWindow::loadJobFromPipeline(int this_job)
{
	// Set the "static int" to which job we're currently pointing
	current_job = this_job;
	int itype = pipeline.processList[current_job].type;
	// The following line allows certain browse buttons to only open the current directory (using CURRENT_ODIR)
	current_browse_directory = pipeline.processList[current_job].name;

	// What type of job is this?
	for ( int t=0; t<NR_BROWSE_TABS; t++ )
	{
		if ( gui_jobwindows[t]->myjob.type == itype )
			browser->value(t+1);
	}

	// change GUI to the corresponding jobwindow
	cb_select_browsegroup_i();

	// Re-read the settings for this job and update the values inside the GUI
	int iwin = (browser->value() - 1);
	gui_jobwindows[iwin]->myjob.read(pipeline.processList[current_job].name, is_main_continue);
	gui_jobwindows[iwin]->updateMyGui();

	// If a finished or running job was loaded from the pipeline: set this to be a continuation job
	// If a scheduled job was loaded, only set is_main_continue to true when it is PROC_SCHEDULED
    //if (pipeline.processList[current_job].status == PROC_SCHEDULED && !gui_jobwindows[iwin]->myjob.is_continue)
    //	is_main_continue = false;
    //else
    //	is_main_continue = true;

    // Any job loaded from the pipeline will initially be set as a continuation job
    // but for show-scheduler, no job should be a continuation
	if (show_scheduler)
	{
		is_main_continue = false;
		do_overwrite_continue = true;
	}
	else
	{
		is_main_continue = true;
	}
    cb_toggle_continue_i();

    // Set the alias in the window
    alias_current_job->value((getJobNameForDisplay(pipeline.processList[current_job])).c_str());

	// Update all job lists in the main GUI
	updateJobLists();

	// File the out and err windows
	fillStdOutAndErr();

}

void GuiMainWindow::cb_select_browsegroup(Fl_Widget* o, void* v)
{

	GuiMainWindow* T=(GuiMainWindow*)v;

	// When clicking the job browser on the left: reset current_job to -1 (i.e. a new job, not yet in the pipeline)
	current_job = -1;
	T->cb_select_browsegroup_i();
	run_button->activate();

}

void GuiMainWindow::cb_select_browsegroup_i(bool show_initial_screen)
{

	// Update timer
	tickTimeLastChanged();

	// Hide the initial screen
	if (show_initial_screen)
		background_grp->show();
	else
		background_grp->hide();

	int iwin = (browser->value() - 1);
	if (iwin < 0 || iwin >= NR_BROWSE_TABS) return;
	// Show the 'selected' group, hide the others
	for ( int t=0; t<NR_BROWSE_TABS; t++ )
    {
    	// During the initial screen: show a nice picture with some explanations
    	if ( t == iwin && !show_initial_screen) // browser starts counting at 1...
        {
    		browse_grp[t]->show();
        }
        else
        {
        	browse_grp[t]->hide();
        }
    }

	// Update all job lists in the main GUI
	updateJobLists();

    is_main_continue = false;
    do_overwrite_continue = false;

	// If the GUI got changed, put that change into the joboption now
    gui_jobwindows[iwin]->updateMyJob();

	// toggle the continue status of this job
	cb_toggle_continue_i();

	alias_current_job->value("Give_alias_here");

	scheduler_job_name->value("");
	scheduler_job_name->activate();

	textbuff_stdout->text("stdout will go here; double-click this window to open stdout in a separate window");
	textbuff_stderr->text("stderr will go here; double-click this window to open stderr in a separate window");

}

void GuiMainWindow::cb_select_finished_job(Fl_Widget* o, void* v)
{
	GuiMainWindow* T=(GuiMainWindow*)v;
	T->cb_select_finished_job_i();
	run_button->activate();
}

void GuiMainWindow::cb_select_finished_job_i()
{
	// Update timer
	tickTimeLastChanged();

	// Show the 'selected' group, hide the others
    int idx = finished_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
		loadJobFromPipeline(finished_processes[idx]);
}

void GuiMainWindow::cb_select_running_job(Fl_Widget* o, void* v)
{
	GuiMainWindow* T=(GuiMainWindow*)v;
	T->cb_select_running_job_i();
	run_button->activate();
}

void GuiMainWindow::cb_select_running_job_i()
{
	// Update timer
	tickTimeLastChanged();

	// Show the 'selected' group, hide the others
    int idx = running_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
    	loadJobFromPipeline(running_processes[idx]);

}

void GuiMainWindow::cb_select_scheduled_job(Fl_Widget* o, void* v)
{
	GuiMainWindow* T=(GuiMainWindow*)v;
	T->cb_select_scheduled_job_i();
	run_button->activate();
}

void GuiMainWindow::cb_select_scheduled_job_i()
{
	// Update timer
	tickTimeLastChanged();

	// Show the 'selected' group, hide the others
    int idx = scheduled_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
		loadJobFromPipeline(scheduled_processes[idx]);

    if (show_scheduler)
    {
		FileName jobname = getJobNameForDisplay(pipeline.processList[current_job]);
    	scheduler_job_name->value(jobname.c_str());
		scheduler_job_name->deactivate();
		bool found = false;
		for (int i = 0; i < 3; i++)
		{
			if (schedule.jobs[jobname].mode == job_mode_options[i].label())
			{
				found = true;
				job_mode->picked(&job_mode_options[i]);
				job_mode_options[job_mode->value()].label();
			}
		}
		if (!found) REPORT_ERROR("ERROR: unrecognised job_mode ...");
    }
}

void GuiMainWindow::cb_select_input_job(Fl_Widget* o, void* v)
{
	GuiMainWindow* T=(GuiMainWindow*)v;
	T->cb_select_input_job_i();
	run_button->activate();
}

void GuiMainWindow::cb_select_input_job_i()
{
	// Update timer
	tickTimeLastChanged();

	// Show the 'selected' group, hide the others
    int idx = input_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
		loadJobFromPipeline(input_processes[idx]);

}

void GuiMainWindow::cb_select_output_job(Fl_Widget* o, void* v)
{
	GuiMainWindow* T=(GuiMainWindow*)v;
	T->cb_select_output_job_i();
	run_button->activate();
}

void GuiMainWindow::cb_select_output_job_i()
{
	// Update timer
	tickTimeLastChanged();

	// Show the 'selected' group, hide the others
    int idx = output_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
		loadJobFromPipeline(output_processes[idx]);

}

void GuiMainWindow::cb_display_io_node(Fl_Widget* o, void* v)
{
	GuiMainWindow* T=(GuiMainWindow*)v;
	T->cb_display_io_node_i();
	run_button->activate();
}

void GuiMainWindow::cb_display_io_node_i()
{

	// Run relion_display on the output node
	int idx = display_io_node->value();
	long int mynode = io_nodes[idx];
	std::string command;

	if (pipeline.nodeList[mynode].type == NODE_MIC_COORDS)
	{
		// A manualpicker jobwindow for display of micrographs....
		RelionJob manualpickjob;
		FileName fn_job = ".gui_manualpickrun.job";
		bool iscont=false;
		if (exists(fn_job))
		{
			manualpickjob.read(fn_job.beforeLastOf("run.job").c_str(), iscont, true); // true means do initialise
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
				command += " --scale " + manualpickjob.joboptions["micscale"].getString();
				command += " --sigma_contrast " + manualpickjob.joboptions["sigma_contrast"].getString();
				command += " --black " + manualpickjob.joboptions["black_val"].getString();
				command += " --white " + manualpickjob.joboptions["white_val"].getString();

				if (manualpickjob.joboptions["lowpass"].getNumber() > 0.)
					command += " --lowpass " + manualpickjob.joboptions["lowpass"].getString();
				if (manualpickjob.joboptions["highpass"].getNumber() > 0.)
					command += " --highpass " + manualpickjob.joboptions["highpass"].getString();
				if (manualpickjob.joboptions["angpix"].getNumber() > 0.)
					command += " --angpix " + manualpickjob.joboptions["angpix"].getString();

				command += " --ctf_scale " + manualpickjob.joboptions["ctfscale"].getString();

				command += " --particle_diameter " + manualpickjob.joboptions["diameter"].getString();

				if (manualpickjob.joboptions["do_color"].getBoolean())
				{
					command += " --color_label " + manualpickjob.joboptions["color_label"].getString();
					command += " --blue " + manualpickjob.joboptions["blue_value"].getString();
					command += " --red " + manualpickjob.joboptions["red_value"].getString();
					if (manualpickjob.joboptions["fn_color"].getString().length() > 0)
						command += " --color_star " + manualpickjob.joboptions["fn_color"].getString();
				}

				// Other arguments for extraction
				command += " " + manualpickjob.joboptions["other_args"].getString() + " &";
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
		char mydefault[]=DEFAULTPDFVIEWER;
		if (default_pdf_viewer == NULL)
		{
			default_pdf_viewer=mydefault;
		}
		std::string myviewer(default_pdf_viewer);
		command = myviewer + " " + pipeline.nodeList[mynode].name + "&";
	}
	else if (pipeline.nodeList[mynode].type != NODE_POST)
	{
		command = "relion_display --gui --i " + pipeline.nodeList[mynode].name + " &";
	}
	//std::cerr << " command= " << command << std::endl;
	int res= system(command.c_str());

}

void GuiMainWindow::cb_set_scheduler_variable(Fl_Widget* o, void*v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_set_scheduler_variable_i();
}

void GuiMainWindow::cb_set_scheduler_variable_i()
{
	schedule.setVariable(scheduler_set_variable_name->value(), scheduler_set_variable_value->value());
	fillSchedulerNodesAndVariables();
	schedule.write();

}

void GuiMainWindow::cb_delete_scheduler_variable(Fl_Widget* o, void*v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_delete_scheduler_variable_i();
}

void GuiMainWindow::cb_delete_scheduler_variable_i()
{
	std::string ask = "Are you sure you want to delete this variable?";
	int proceed =  fl_choice("%s", "Cancel", "Delete!", NULL, ask.c_str());
	if (!proceed)
	{
		do_overwrite_continue = false;
		return;
	}

	schedule.removeVariable(scheduler_set_variable_name->value());
	fillSchedulerNodesAndVariables();
	schedule.write();
}

void GuiMainWindow::cb_set_scheduler_operator(Fl_Widget* o, void*v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_set_scheduler_operator_i();
}

void GuiMainWindow::cb_set_scheduler_operator_i()
{

}

void GuiMainWindow::cb_delete_scheduler_operator(Fl_Widget* o, void*v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_delete_scheduler_operator_i();
}

void GuiMainWindow::cb_delete_scheduler_operator_i()
{
	std::string ask = "Are you sure you want to delete this operator?";
	int proceed =  fl_choice("%s", "Cancel", "Delete!", NULL, ask.c_str());
	if (!proceed)
	{
		do_overwrite_continue = false;
		return;
	}

	schedule.removeOperator(scheduler_set_variable_name->value());
	fillSchedulerNodesAndVariables();
	schedule.write();
}

void GuiMainWindow::cb_select_scheduler_variable(Fl_Widget* o, void*v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_select_scheduler_variable_i();
}

void GuiMainWindow::cb_select_scheduler_variable_i()
{
	// Get position of the browser:
	int idx = scheduler_variable_browser->value();
	FileName mytext = scheduler_variable_browser->text(idx);
	FileName myname = mytext.beforeFirstOf(" = ");
	FileName myval = mytext.afterFirstOf(" = ");
	myval = myval.beforeFirstOf(" (");
	std::cerr << scheduler_variable_browser->text(idx) << std::endl;
	std::cerr << scheduler_variable_browser->selected(idx) << std::endl;
	scheduler_set_variable_name->value(myname.c_str());
	scheduler_set_variable_value->value(myval.c_str());

}

void GuiMainWindow::cb_select_scheduler_node(Fl_Widget *o, void* v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_select_scheduler_node_i();

}

void GuiMainWindow::cb_select_scheduler_node_i()
{

}


void GuiMainWindow::cb_display(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_display_i();
}


void GuiMainWindow::cb_display_i()
{
        std::string command = " relion_display --gui &" ;
        int res = system(command.c_str());
}

void GuiMainWindow::cb_toggle_continue_i()
{

	if (is_main_continue || do_overwrite_continue)
	{
		if (do_overwrite_continue)
		{
			run_button->label("Overwrite!");
			add_job_button->label("Save");
		}
		else
		{
			run_button->label("Continue!");
		}
		run_button->color(GUI_BUTTON_COLOR);
		run_button->labelfont(FL_ITALIC);
		run_button->labelsize(13);
		alias_current_job->deactivate();
	}
	else
	{
		run_button->label("Run!");
		add_job_button->label("Add job");
		run_button->color(GUI_RUNBUTTON_COLOR);
		run_button->labelfont(FL_ITALIC);
		run_button->labelsize(16);
		alias_current_job->activate();
	}

	int my_window = (browser->value() - 1);
	gui_jobwindows[my_window]->toggle_new_continue(is_main_continue && !do_overwrite_continue);

}


void GuiMainWindow::cb_print_cl(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_print_cl_i();
}

void GuiMainWindow::cb_print_cl_i()
{

	int iwin = browser->value() - 1;
	// And update the job inside it
	gui_jobwindows[iwin]->updateMyJob();

	std::string error_message;
	if (!pipeline.getCommandLineJob(gui_jobwindows[iwin]->myjob, current_job, is_main_continue, false,
			DONT_MKDIR, do_overwrite_continue, commands, final_command, error_message))
	{
		fl_message("%s",error_message.c_str());
	}
	else
	{
		std::cout << " *** The command is:" << std::endl;
		for (int icom = 0; icom < commands.size(); icom++)
			std::cout << commands[icom] << std::endl;
	}

}

// Run button call-back functions
void GuiMainWindow::cb_run(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
	// Deactivate Run button to prevent the user from accidentally submitting many jobs
	run_button->deactivate();
	// Run the job
	T->cb_run_i(false, false); // 1st false means dont only_schedule, 2nd false means dont open the note editor window
}

// Run button call-back functions
void GuiMainWindow::cb_schedule(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_run_i(true, false); // 1st true means only_schedule, do not run, 2nd false means dont open the note editor window
}

void GuiMainWindow::cb_run_i(bool only_schedule, bool do_open_edit)
{

	if (do_overwrite_continue)
	{
		std::string ask = "Are you sure you want to overwrite this job?";
		int proceed =  fl_choice("%s", "Cancel", "Overwrite!", NULL, ask.c_str());
		if (!proceed)
		{
			do_overwrite_continue = false;
			return;
		}
	}

	// Get which jobtype the GUI is on now
	int iwin = browser->value() - 1;
	// And update the job inside it
	gui_jobwindows[iwin]->updateMyJob();

	// Update timer
	tickTimeLastChanged();

	std::string error_message;
	if (!pipeline.runJob(gui_jobwindows[iwin]->myjob, current_job, only_schedule, is_main_continue, false, do_overwrite_continue, error_message))
	{
		fl_message("%s",error_message.c_str());
		// Allow the user to fix the error and submit this job again
		run_button->activate();
		return;
	}

	// Update all job lists in the main GUI
	updateJobLists();

	// Open the edit note window
	if (do_open_edit)
	{
		// Open the note editor window
		cb_edit_note_i();
	}

	// Also set alias from the alias_current_job input
	if (!(is_main_continue || do_overwrite_continue))
	{
		std::string alias= (std::string)alias_current_job->value();
		if (alias != "Give_alias_here" && alias != pipeline.processList[current_job].name)
			cb_set_alias_i(alias);
	}

	do_overwrite_continue = false;

	// Select this job now
	loadJobFromPipeline(current_job);

}

void GuiMainWindow::cb_scheduler_add_job(Fl_Widget* o, void* v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_scheduler_add_job_i();
}

void GuiMainWindow::cb_scheduler_add_job_i()
{

	// Get which jobtype the GUI is on now
	int iwin = browser->value() - 1;
	// And update the job inside it
	gui_jobwindows[iwin]->updateMyJob();

	std::string mode = job_mode_options[job_mode->value()].label();
	std::string jobname = scheduler_job_name->value();

	if (do_overwrite_continue)
	{
		// Write the possibly updated job settings
		gui_jobwindows[iwin]->myjob.write(pipeline.processList[current_job].name);

		// Also write the possibly updated job_mode
		std::string mode = job_mode_options[job_mode->value()].label();
		schedule.jobs[jobname].mode = mode;
		schedule.write();

		std::cerr << "TODO: write has_started variable as well, add button to edit this!" << std::endl;

	}
	else
	{
		// Add job to the schedule
		// Get the mode, and the jobname
		if (jobname == "")
		{
			fl_message("%s","You need to provide a Name for this job in the scheduler.");
			return;
		}

		schedule.addJob(gui_jobwindows[iwin]->myjob, jobname, mode);

		scheduler_job_name->value("");
		schedule.write();
		updateJobLists();
	}

}

// Run button call-back functions
void GuiMainWindow::cb_delete(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_delete_i();
}

void GuiMainWindow::cb_delete_i(bool do_ask, bool do_recursive)
{

	if (current_job < 0)
	{
		std::cout << " You can only delete existing jobs ... " << std::endl;
		return;
	}

	std::vector<bool> deleteProcesses, deleteNodes;
	pipeline.deleteJobGetNodesAndProcesses(current_job, do_recursive, deleteNodes, deleteProcesses);

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
		proceed =  fl_choice("%s", "Cancel", "Move", NULL, ask.c_str());
	}
	else
	{
		proceed = 1;
	}

	if (proceed)
	{

		pipeline.deleteNodesAndProcesses(deleteNodes, deleteProcesses);

		// Reset current_job
		current_job = -1;
		fillStdOutAndErr();

		// Update all job lists in the main GUI
		updateJobLists();

	}

}

// Run button call-back functions
void GuiMainWindow::cb_gently_clean_all_jobs(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_clean_all_jobs_i(false);
}

// Run button call-back functions
void GuiMainWindow::cb_harshly_clean_all_jobs(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_clean_all_jobs_i(true);
}

void GuiMainWindow::cb_clean_all_jobs_i(bool do_harsh)
{


	int proceed = 1;
	std::string ask;
	if (do_harsh)
	{
		ask = "Are you sure you want to harshly clean up intermediate files from the entire pipeline? \n\n\
Harsh cleaning will remove micrographs, movies and particle stacks from all MotionCorr, Extract, \n\
Polish and Subtract directories. This means you will NOT be able to use those images in subsequent runs anymore, \n\
although you could always recreate the data by continuing the job (possibly at considerable computing costs).\n \n \
You can protect specific jobs from harsh cleaning by creating a file called \"NO_HARSH_CLEAN\" inside their directory,\n\
e.g. by using \"touch Polish/job045/NO_HARSH_CLEAN\". Below is a list of currently protected jobs (if any):\n \n";

		for (int myjob = 0; myjob < pipeline.processList.size(); myjob++)
		{
			if (pipeline.processList[myjob].status == PROC_FINISHED_SUCCESS &&
					(pipeline.processList[myjob].type == PROC_MOTIONCORR ||
					pipeline.processList[myjob].type == PROC_EXTRACT ||
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

	proceed = fl_choice("%s", "Cancel", "Clean up", NULL, ask.c_str());
	if (proceed)
	{
		std::string how = (do_harsh) ? "Harshly" : "Gently";
		std::cout << how << " cleaning all finished jobs ..." << std::endl;

		std::string error_message;
		if (!pipeline.cleanupAllJobs(do_harsh, error_message))
			fl_message("%s",error_message.c_str());

		fl_message("Done cleaning! Don't forget the files are all still in the Trash folder. Use the \"Empty Trash\" option from the File menu to permanently delete them.");
	}
}

// Run button call-back functions
void GuiMainWindow::cb_gentle_cleanup(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_cleanup_i(-1, true, false);
}

void GuiMainWindow::cb_harsh_cleanup(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_cleanup_i(-1, true, true);
}

void GuiMainWindow::cb_cleanup_i(int myjob, bool do_verb, bool do_harsh)
{
	// Allow cleaning the currently selected job from the GUI
	if (myjob < 0)
		myjob = current_job;

	int proceed = 1;
	if (do_verb)
	{
		std::string ask;
		ask = "Are you sure you want to clean up intermediate files from " + pipeline.processList[current_job].name + "?";
		proceed = fl_choice("%s", "Cancel", "Clean up", NULL, ask.c_str());
	}

	if (proceed)
	{
		std::string error_message;
		if (!pipeline.cleanupJob(myjob, do_harsh, error_message))
			fl_message("%s",error_message.c_str());
	}

}


// Run button call-back functions
void GuiMainWindow::cb_set_alias(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_set_alias_i();
}

void GuiMainWindow::cb_set_alias_i(std::string alias)
{

	FileName fn_pre, fn_jobnr, fn_post, fn_dummy, default_ask;
	if (!decomposePipelineFileName(pipeline.processList[current_job].name, fn_pre, fn_jobnr, fn_post))
		REPORT_ERROR("GuiMainWindow::cb_set_alias_i ERROR: invalid pipeline process name: " + pipeline.processList[current_job].name);

	// Start the asking window with the current alias
	std::string error_message;
	FileName fn_alias = pipeline.processList[current_job].alias;
	if (fn_alias != "None")
	{
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
			if (palias == NULL)
				return;
			std::string al2(palias);
			alias = al2;
		}

		if (pipeline.setAliasJob(current_job, alias, error_message))
			is_done = true;
		else
		{
			alias = "";
			fl_message("%s",error_message.c_str());
		}
	}

}

// Run button call-back functions
void GuiMainWindow::cb_abort(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_abort_i();
}

void GuiMainWindow::cb_abort_i(std::string alias)
{

	if (pipeline.processList[current_job].status != PROC_RUNNING)
	{
		std::string error_message = "You can only abort running jobs ... ";
		fl_message("%s",error_message.c_str());
	}
	else
	{
		std::string ask = "Are you sure you want to abort job: " + pipeline.processList[current_job].name + " ?";
		int proceed =  fl_choice("%s", "Cancel", "Abort!", NULL, ask.c_str());
		if (proceed)
		{
			touch(pipeline.processList[current_job].name + RELION_JOB_ABORT_NOW);
		}
	}
}



// Run button call-back functions
void GuiMainWindow::cb_mark_as_finished(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_mark_as_finished_i();
}

void GuiMainWindow::cb_mark_as_finished_i()
{

	if (current_job < 0)
	{
		fl_message("You can only mark existing jobs as finished!");
		return;
	}

	std::string error_message;
	if (!pipeline.markAsFinishedJob(current_job, error_message))
		fl_message("%s",error_message.c_str());
	else
		updateJobLists();

}

// Run button call-back functions
void GuiMainWindow::cb_make_flowchart(Fl_Widget* o, void* v) {

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_make_flowchart_i();
}

void GuiMainWindow::cb_make_flowchart_i()
{

	std::string error_message;
	if (!pipeline.makeFlowChart(current_job, true, error_message))
		fl_message("%s",error_message.c_str());
	else
		updateJobLists();

}

void GuiMainWindow::cb_edit_note(Fl_Widget*, void* v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_edit_note_i();

}

void GuiMainWindow::cb_edit_project_note(Fl_Widget*, void* v)
{
  	GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_edit_note_i(true); // true means is_project_note
}

void GuiMainWindow::cb_edit_note_i(bool is_project_note)
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
	NoteEditorWindow* w = new NoteEditorWindow(660, 400, title.c_str(), fn_note, !maingui_do_read_only);
	w->show();

}


// Save button call-back function
void GuiMainWindow::cb_save(Fl_Widget* o, void* v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_save_i();
}

void GuiMainWindow::cb_save_i()
{
	// Get which job we're dealing with, and update it from the GUI
	int iwin = browser->value() - 1;
	gui_jobwindows[iwin]->updateMyJob();

	// For scheduled jobs, also allow saving the .job file in the output directory
	if (current_job >= 0 && (pipeline.processList[current_job].status == PROC_SCHEDULED))
	{
		gui_jobwindows[iwin]->myjob.write(pipeline.processList[current_job].name);
	}
	// Write the hidden file
	gui_jobwindows[iwin]->myjob.write("");

}

// Load button call-back function
void GuiMainWindow::cb_load(Fl_Widget* o, void* v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_load_i();
}

void GuiMainWindow::cb_load_i()
{
	int iwin = browser->value() - 1;
	gui_jobwindows[iwin]->myjob.read("", is_main_continue);
	alias_current_job->value("Give_alias_here");
	gui_jobwindows[iwin]->updateMyGui();

	// Make the current continue-setting active
	cb_toggle_continue_i();
}

// Load button call-back function
void GuiMainWindow::cb_undelete_job(Fl_Widget* o, void* v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_undelete_job_i();
}

void GuiMainWindow::cb_undelete_job_i()
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

	pipeline.undeleteJob(fn_pipe);

}

void GuiMainWindow::cb_export_jobs(Fl_Widget* o, void* v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_export_jobs_i();
}

void GuiMainWindow::cb_export_jobs_i()
{

	// Get the name of this block of exported jobs and make the corresponding directory
	const char * answer;
	std::string default_answer="export1";
	answer =  fl_input("Name of the exported block of jobs? ", default_answer.c_str());
	std::string mydir(answer);

	std::string error_message;
	if (!pipeline.exportAllScheduledJobs(mydir, error_message))
		fl_message("%s",error_message.c_str());

}

void GuiMainWindow::cb_import_jobs(Fl_Widget* o, void* v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_import_jobs_i();
}


void GuiMainWindow::cb_import_jobs_i()
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

	pipeline.importJobs(fn_export);

	// refresh the joblists
	updateJobLists();
}

// Re-order running and finished job lists
void GuiMainWindow::cb_order_jobs_alphabetically(Fl_Widget* o, void* v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    do_order_alphabetically = true;
    T->fillRunningJobLists();
}

// Re-order running and finished job lists
void GuiMainWindow::cb_order_jobs_chronologically(Fl_Widget* o, void* v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    do_order_alphabetically = false;
    T->fillRunningJobLists();
}

// Empty-trash button call-back function
void GuiMainWindow::cb_empty_trash(Fl_Widget* o, void* v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_empty_trash_i();
}

void GuiMainWindow::cb_empty_trash_i()
{
	std::string ask = "Are you sure you want to remove the entire Trash folder?";
	int proceed =  fl_choice("%s", "Cancel", "Empty Trash", NULL, ask.c_str());
	if (proceed)
	{
		std::string command = "rm -rf Trash";
		std::cout << " Executing: " << command << std::endl;
		int res = system(command.c_str());
	}
}

void GuiMainWindow::cb_print_notes(Fl_Widget*, void* v)
{
  	 GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_print_notes_i();
}

void GuiMainWindow::cb_print_notes_i()
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

void GuiMainWindow::cb_remake_nodesdir(Fl_Widget*, void* v)
{
  	 GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_remake_nodesdir_i();
}

void GuiMainWindow::cb_remake_nodesdir_i()
{
	pipeline.remakeNodeDirectory();
}

void GuiMainWindow::cb_reread_pipeline(Fl_Widget*, void* v)
{
  	 GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_reread_pipeline_i();
}

void GuiMainWindow::cb_reread_pipeline_i()
{
	std::string lock_message = " mainGUI reread_pipeline_i";
	pipeline.read(DO_LOCK, lock_message);
	// With the locking system, each read needs to be followed soon with a write
	pipeline.write(DO_LOCK);
}


void GuiMainWindow::cb_reactivate_runbutton(Fl_Widget* o, void* v)
{

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_reactivate_runbutton_i();
}

void GuiMainWindow::cb_reactivate_runbutton_i()
{
	run_button->activate();
}

void GuiMainWindow::cb_toggle_overwrite_continue(Fl_Widget* o, void* v)
{

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_toggle_overwrite_continue_i();
}

void GuiMainWindow::cb_toggle_overwrite_continue_i()
{
    do_overwrite_continue = !do_overwrite_continue;

	cb_toggle_continue_i();
}

void GuiMainWindow::cb_show_initial_screen(Fl_Widget* o, void* v)
{

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_show_initial_screen_i();
}

void GuiMainWindow::cb_show_initial_screen_i()
{
    run_button->deactivate();

	cb_select_browsegroup_i(true);
}

void GuiMainWindow::cb_toggle_pipeliner_scheduler(Fl_Widget* o, void* v)
{
    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_toggle_pipeliner_scheduler_i();
}

void GuiMainWindow::cb_toggle_pipeliner_scheduler_i()
{
	if (show_scheduler)
	{
		pipeliner_grp->hide();
		scheduler_grp->show();
	}
	else
	{
		scheduler_grp->hide();
		pipeliner_grp->show();
	}
}

void GuiMainWindow::cb_start_pipeliner(Fl_Widget* o, void* v)
{

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_start_pipeliner_i();
}

void GuiMainWindow::cb_start_pipeliner_i()
{

	std::vector<FileName> job_names;

	for (long int ii =0; ii < scheduled_processes.size(); ii++)
	{
		long int id = scheduled_processes[ii];
		job_names.push_back(pipeline.processList[id].name);
	}
	SchedulerWindow* w = new SchedulerWindow(400, 300, "Select which jobs to execute");
	w->fill(pipeline.name, job_names);

}

void GuiMainWindow::cb_stop_pipeliner(Fl_Widget* o, void* v)
{

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_stop_pipeliner_i();
}

void GuiMainWindow::cb_stop_pipeliner_i()
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

void GuiMainWindow::cb_about(Fl_Widget* o, void* v)
{

    GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_about_i();
}

void GuiMainWindow::cb_about_i()
{

#define HELPTEXT ("RELION " RELION_SHORT_VERSION "\n \n \
RELION is is developed in the groups of\n\n \
Sjors H.W. Scheres at the MRC Laboratory of Molecular Biology\n \n \
- Sjors H.W. Scheres\n \
- Shaoda He\n \
- Takanori Nakane\n \
- Jasenko Zivanov\n \
- Liyi Dong\n \
\n \n \
and Erik Lindahl at Stockholm University\n \n \
- Erik Lindahl\n \
- Björn O. Forsberg\n \
- Dari Kimanius\n \
\n\
Note that RELION is completely free, open-source software. You can redistribute it and/or modify it for your own purposes, but please do make sure \
the contribution of the developers are acknowledged appropriately. In order to maintain an overview of existing versions, a notification regarding  \
any redistribution of (modified versions of) the code is appreciated (contact Sjors directly). \n \n \n \
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
 * Auto-picking : \n \
     Scheres (2014) J. Struct. Biol. (PMID: 25486611) \n \n \
 * Sub-tomogram averaging : \n \
     Bharat et al. (2015) Structure (PMID: 26256537) \n \n \
 * v.2.0 GPU capability and autopicking acceleration : \n \
     Kimanius et al. (2016) eLife (PMID: 27845625) \n \n \
 * Helical reconstruction : \n \
     He & Scheres (2017) J, Struct. Biol. (PMID: 28193500) \n \n \
Please also cite the following EXTERNAL programs: \n \n \
* MOTIONCOR2 for beam-induced motion correction: \n \
    Zheng et al (2017) Nat. Meth. (PMID: 28250466) \n \n\
* UNBLUR for beam-induced motion correction: \n \
    Grabnt & Grigorieff eLife (PMID: 26023829) \n \n\
* CTFFIND4 for CTF-estimation: \n \
    Rohou & Grigorieff (2015) J. Struct. Biol. (PMID: 26278980) \n \n\
* Gctf for CTF-estimation: \n \
    Zhang (2016) J. Struct. Biol. (PMID: 2659270) \n \n\
* Stochastic Gradient Descent initial model generation:  \n\
    Punjani et al. (2017) Nat. Meth. (PMID: 28165473) \n \n\
* ResMap for local-resolution estimation:  \n\
    Kucukelbir et al. (2014) Nat. Meth. (PMID: 24213166) \n \n\
* Postscript plots are made using CPlot2D from  www.amzsaki.com\n ")

	ShowHelpText *help = new ShowHelpText(HELPTEXT);
}



void GuiMainWindow::cb_quit(Fl_Widget* o, void* v)
{
	GuiMainWindow* T=(GuiMainWindow*)v;
    T->cb_quit_i();
}

void GuiMainWindow::cb_quit_i()
{
	exit(0);
}

