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


NoteEditorWindow::NoteEditorWindow(int w, int h, const char* title, FileName _fn_note):Fl_Window(w,h,title)
{
	editor = new Fl_Text_Editor(0, 0, w, h-50);
    editor->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS,10);
	textbuff_note = new Fl_Text_Buffer;
	editor->buffer(textbuff_note);
	fn_note = _fn_note;
	if (exists(fn_note))
		int errno = textbuff_note->loadfile(fn_note.c_str());
	else
		textbuff_note->text("Describe what this job is about here...");

	// Button to exit
	Fl_Button *cancel_button = new Fl_Button(w-200, h-40, 80, 30, "Cancel");
	cancel_button->color(GUI_RUNBUTTON_COLOR);
	cancel_button->labelsize(12);
	cancel_button->callback( cb_cancel, this);

	// Button to save and exit
	Fl_Button *save_button = new Fl_Button(w-100, h-40, 80, 30, "Save");
	save_button->color(GUI_RUNBUTTON_COLOR);
	save_button->labelsize(12);
	save_button->callback( cb_save, this);


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
	int errno = textbuff_note->savefile(fn_note.c_str());
}



RelionMainWindow::RelionMainWindow(int w, int h, const char* title, FileName fn_pipe):Fl_Window(w,h,title)
{

	show_initial_screen = true;

	FileName fn_lock=".gui_projectdir";
	if (!exists(fn_lock))
	{
		std::cout << " Only run the relion GUI from your ProjectDirectory. Do you want to start a new project here [y/n]? ";
		char c;
		std::cin >> c;
		if (c == 'y' || c == 'Y')
		{
			std::string command = " touch .gui_projectdir ";
			int res= system(command.c_str());
		}
		else
		{
			std::cout << " Exiting ... " << std::endl;
			exit(0);
		}
	}

	// First setup the old part of the GUI
	h = GUIHEIGHT_OLD;

	// Initial screen picture with some density and some explanation
	//fl_register_images(); // initialize image lib
	//image_box = new Fl_Box(0,100,w,h); // widget that will contain image
	// TODO: control file location and use better figure
	//jpeg_image = new Fl_JPEG_Image("/lmb/home/scheres/bg.jpg"); // load jpeg image into ram
	//image_box->image(jpeg_image); // attach jpg image to box

	// Read in the pipeline STAR file if it exists
	pipeline.name = fn_pipe;
	if (exists(fn_pipe + "_pipeline.star"))
		pipeline.read();

	// Check which jobs have finished
	pipeline.checkProcessCompletion();

	// Make temporary directory with all NodeNames
	pipeline.makeNodeDirectory();

    color(GUI_BACKGROUND_COLOR);
    menubar = new Fl_Menu_Bar(-3, 0, WCOL0-7, MENUHEIGHT);
    menubar->add("File/Save job settings",  FL_ALT+'s', cb_menubar_save, this);
    menubar->add("File/Display",  FL_ALT+'d', cb_display, this);
    menubar->add("File/Reactivate Run",  FL_ALT+'r', cb_menubar_reactivate_runbutton, this);
    menubar->add("File/Print all notes",  FL_ALT+'p', cb_menubar_print_notes, this);
    menubar->add("File/Show initial screen",  FL_ALT+'a', cb_show_initial_screen, this);
    menubar->add("File/About",  FL_ALT+'a', cb_menubar_about, this);
    menubar->add("File/Quit", FL_ALT+'q', cb_menubar_quit, this);
    menubar->add("Autorun/Run scheduled jobs", 0, cb_menubar_run_scheduled_jobs, this);
    menubar->add("Autorun/Stop running scheduled jobs", 0, cb_menubar_stop_run_scheduled_jobs, this);

    current_y = MENUHEIGHT + 10;

    // Add run buttons on the menubar as well
	print_CL_button = new Fl_Button(GUIWIDTH - 330, h-89, 100, 30, "Print command");
	print_CL_button->color(GUI_RUNBUTTON_COLOR);
	print_CL_button->labelsize(12);
	print_CL_button->callback( cb_print_cl, this);

	schedule_button = new Fl_Button(GUIWIDTH - 220 , h-89, 100, 30, "Schedule");
	schedule_button->color(GUI_RUNBUTTON_COLOR);
	schedule_button->labelfont(FL_ITALIC);
	schedule_button->labelsize(14);
	schedule_button->callback( cb_schedule, this);

	run_button = new Fl_Button(GUIWIDTH - 110 , h-89, 100, 30, "Run now");
	run_button->color(GUI_RUNBUTTON_COLOR);
	run_button->labelfont(FL_ITALIC);
	run_button->labelsize(14);
	run_button->callback( cb_run, this);

    // Fill browser in the right order
	browser = new Fl_Hold_Browser(10,MENUHEIGHT+10,WCOL0-20,h-MENUHEIGHT-70);
    current_job = -1;
    for (int itype = 0; itype < NR_BROWSE_TABS; itype++)
    {

        browse_grp[itype] = new Fl_Group(WCOL0, 2, 550, 600-MENUHEIGHT);
        switch (itype)
    	{
    	case PROC_IMPORT:
    	{
    		browser->add("Import");
    		job_import = new ImportJobWindow();
    		browse_grp[itype]->end();
    		break;
    	}
    	case PROC_MOTIONCORR:
    	{
    		browser->add("Motion correction");
    		job_motioncorr = new MotioncorrJobWindow();
    		break;
    	}
    	case PROC_CTFFIND:
    	{
    		browser->add("CTF estimation");
			job_ctffind = new CtffindJobWindow();
			break;
    	}
    	case PROC_MANUALPICK:
    	{
    		browser->add("Manual picking");
        	job_manualpick = new ManualpickJobWindow();
    		break;
    	}
    	case PROC_AUTOPICK:
    	{
    		browser->add("Auto-picking");
			job_autopick = new AutopickJobWindow();
			break;
    	}
    	case PROC_EXTRACT:
    	{
    		browser->add("Particle extraction");
			job_extract = new ExtractJobWindow();
			break;
    	}
    	case PROC_SORT:
    	{
    		browser->add("Particle sorting");
        	job_sort = new SortJobWindow();
        	break;
    	}
    	case PROC_CLASSSELECT:
    	{
    		browser->add("Subset selection");
			job_classselect = new ClassSelectJobWindow();
			break;
		}
    	case PROC_2DCLASS:
		{
    		browser->add("2D classification");
			job_class2d = new Class2DJobWindow();
			break;
		}
    	case PROC_3DCLASS:
    	{
    		browser->add("3D classification");
        	job_class3d = new Class3DJobWindow();
        	break;
    	}
    	case PROC_3DAUTO:
    	{
    		browser->add("3D auto-refine");
			job_auto3d = new Auto3DJobWindow();
			break;
		}
    	case PROC_POLISH:
		{
    		browser->add("Particle polishing");
			job_polish = new PolishJobWindow();
			break;
		}
    	case PROC_MASKCREATE:
    	{
    		browser->add("Mask creation");
			job_maskcreate = new MaskCreateJobWindow();
			break;
		}
    	case PROC_JOINSTAR:
    	{
    		browser->add("Join star files");
			job_joinstar = new JoinStarJobWindow();
			break;
		}
    	case PROC_SUBTRACT:
    	{
    		browser->add("Particle subtraction");
    		job_subtract = new SubtractJobWindow();
			break;
		}
    	case PROC_POST:
		{
    		browser->add("Post-processing");
			job_post = new PostJobWindow();
			break;
		}
    	case PROC_RESMAP:
		{
    		browser->add("Local resolution");
			job_resmap = new ResmapJobWindow();
			break;
		}
    	default:
    	{
           	break;
    	}
    	} // end switch
    	browse_grp[itype]->end();
    }
    browser->callback(cb_select_browsegroup);
    browser->end();
    browser->select(1); // just start from the beginning

    // Pipeline part of the GUI
#define JOBCOLWIDTH (250)
#define XJOBCOL1 (10)
#define XJOBCOL2 (JOBCOLWIDTH + 25)
#define XJOBCOL3 (2*JOBCOLWIDTH + 40)

    menubar2 = new Fl_Menu_Bar(XJOBCOL1, GUIHEIGHT_EXT_START, 100, MENUHEIGHT);
    menubar2->color(GUI_BUTTON_COLOR);
    menubar2->add("Job actions/Edit Note", 0, cb_edit_note, this);
    menubar2->add("Job actions/Alias", 0, cb_set_alias, this);
    menubar2->add("Job actions/Mark as finished", 0, cb_mark_as_finished, this);
    menubar2->add("Job actions/Clean up", 0, cb_cleanup, this);
    menubar2->add("Job actions/Delete", 0, cb_delete, this);

    // Text display with the name of the current job
    text_current_job = new Fl_Text_Display(XJOBCOL2-50 , GUIHEIGHT_EXT_START+3, 250, MENUHEIGHT-6);
    textbuff_current_job = new Fl_Text_Buffer();
    text_current_job->label("Current job: ");
    text_current_job->align(FL_ALIGN_LEFT);
    text_current_job->color(GUI_BACKGROUND_COLOR, GUI_BACKGROUND_COLOR);
    text_current_job->buffer(textbuff_current_job);

    // Left-hand side browsers for input/output nodes and processes
	display_io_node  = new Fl_Choice(XJOBCOL3, GUIHEIGHT_EXT_START+3, 250, MENUHEIGHT-6);
    display_io_node->label("Display:");
    display_io_node->color(GUI_BUTTON_COLOR);
	display_io_node->callback(cb_display_io_node, this);

	// Add browsers for finished, running and scheduled jobs
    int JOBHEIGHT = 170;
    int JOBHALFHEIGHT = (JOBHEIGHT)/2;
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
	Fl_Text_Display* textdisp3 = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_EXT_START2+JOBHALFHEIGHT+25, JOBCOLWIDTH, 25);
	Fl_Text_Display* textdisp4 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_EXT_START2, JOBCOLWIDTH, 25);
	Fl_Text_Display* textdisp5 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_EXT_START2+JOBHALFHEIGHT+25, JOBCOLWIDTH, 25);
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

    finished_job_browser  = new Fl_Select_Browser(XJOBCOL1, GUIHEIGHT_EXT_START2+25, JOBCOLWIDTH, JOBHEIGHT+25);
    running_job_browser   = new Fl_Select_Browser(XJOBCOL2, GUIHEIGHT_EXT_START2+25, JOBCOLWIDTH, JOBHALFHEIGHT);
    scheduled_job_browser = new Fl_Select_Browser(XJOBCOL2, GUIHEIGHT_EXT_START2+25+JOBHALFHEIGHT+25, JOBCOLWIDTH, JOBHALFHEIGHT);
    input_job_browser    = new Fl_Select_Browser(XJOBCOL3, GUIHEIGHT_EXT_START2+25, JOBCOLWIDTH, JOBHALFHEIGHT);
    output_job_browser   = new Fl_Select_Browser(XJOBCOL3, GUIHEIGHT_EXT_START2+25+JOBHALFHEIGHT+25, JOBCOLWIDTH, JOBHALFHEIGHT);

    // Fill the actual browsers
    fillRunningJobLists();

    // Set the callbacks
    finished_job_browser->callback(cb_select_finished_job);
    running_job_browser->callback(cb_select_running_job);
    scheduled_job_browser->callback(cb_select_scheduled_job);
    input_job_browser->callback(cb_select_input_job);
    output_job_browser->callback(cb_select_output_job);

    finished_job_browser->end();
    running_job_browser->end();
    scheduled_job_browser->end();
    input_job_browser->end();
    output_job_browser->end();

    // Display stdout and stderr of jobs
    textbuff_stdout = new Fl_Text_Buffer();
    textbuff_stderr = new Fl_Text_Buffer();
    // Disable warning message about UTF-8 transcoding
    textbuff_stdout->transcoding_warning_action=NULL;
    textbuff_stderr->transcoding_warning_action=NULL;
	disp_stdout = new Fl_Text_Display(XJOBCOL1, GUIHEIGHT_EXT_START2 + JOBHEIGHT + 60, w-20, 110);
    disp_stderr = new Fl_Text_Display(XJOBCOL1, GUIHEIGHT_EXT_START2 + JOBHEIGHT + 170, w-20, 60);
    textbuff_stdout->text("stdout will go here");
    textbuff_stderr->text("stderr will go here");
    disp_stdout->buffer(textbuff_stdout);
    disp_stderr->buffer(textbuff_stderr);
    disp_stderr->textcolor(FL_RED);
    disp_stdout->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS,0);
    disp_stderr->wrap_mode(Fl_Text_Display::WRAP_AT_BOUNDS,0);
    disp_stdout->scrollbar_width(0);
    disp_stderr->scrollbar_width(0);

    // Set and activate current selection from side-browser
    cb_select_browsegroup_i(); // make default active
    is_main_continue = false; // default is a new run

    resizable();

}

// Update the content of the finished, running and scheduled job lists
void RelionMainWindow::fillRunningJobLists()
{
    // Clear whatever was in there
	finished_job_browser->clear();
	running_job_browser->clear();
	scheduled_job_browser->clear();
	finished_processes.clear();
	running_processes.clear();
	scheduled_processes.clear();

	// Fill the Jobs browsers
    // Search backwards, so that last jobs are at the top
    for (long int i = pipeline.processList.size() -1, jobnr=0 ; i >= 0; i--, jobnr++)
    {
    	int itype = pipeline.processList[i].type;
    	switch (pipeline.processList[i].status)
    	{
    	case PROC_FINISHED:
    	{

    		finished_processes.push_back(i);
    		if (pipeline.processList[i].alias != "None")
    			finished_job_browser->add(pipeline.processList[i].alias.c_str());
    		else
    			finished_job_browser->add(pipeline.processList[i].name.c_str());
    		break;
    	}
    	case PROC_RUNNING:
    	{
    		running_processes.push_back(i);
    		if (pipeline.processList[i].alias != "None")
    			running_job_browser->add(pipeline.processList[i].alias.c_str());
    		else
    			running_job_browser->add(pipeline.processList[i].name.c_str());
    		break;
    	}
    	case PROC_SCHEDULED:
    	{
    		scheduled_processes.push_back(i);
    		if (pipeline.processList[i].alias != "None")
    		    scheduled_job_browser->add(pipeline.processList[i].alias.c_str());
    		else
    		    scheduled_job_browser->add(pipeline.processList[i].name.c_str());

    		break;
    	}
    	default:
    	{
    		REPORT_ERROR("RelionMainWindow::RelionMainWindow ERROR: unrecognised status for job" + pipeline.processList[i].name);
    	}
    	}//end switch

    }

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

	if (pipeline.processList[current_job].status == PROC_SCHEDULED)
		fn_settings = ".ScheduledJobs/" + fn_settings;

	for ( int t=0; t<NR_BROWSE_TABS; t++ )
	{
		if ( t == itype )
		{
			browse_grp[t]->show();
			browser->value(itype+1);
		}
		else
		{
			browse_grp[t]->hide();
		}
	}

	// the new job browser has reset current_job to -1....
	current_job = this_job;

	// Re-read the settings for this job
	cb_select_browsegroup_i(); // change to the corresponding jobwindow
	jobCommunicate(DONT_WRITE, DO_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

	// Update all job lists in the main GUI
	updateJobLists();

	// If an old job was loaded from the pipeline: set this to be a continuation job
    is_main_continue = true;
    cb_toggle_continue_i();

	textbuff_current_job->text(pipeline.processList[current_job].name.c_str());
	text_current_job->buffer(textbuff_current_job);

}

void RelionMainWindow::addToPipeLine(int as_status, bool do_overwrite, int this_job)
{
	int itype = (this_job > 0) ? this_job : (browser->value() - 1); // browser starts counting at 1 ...
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
		REPORT_ERROR("RelionMainWindow::addToPipeLine ERROR: unrecognised type");
	}
	}


	// Add Process to the processList of the pipeline
	Process process(oname, itype, as_status);
	long int myProcess = pipeline.addNewProcess(process, do_overwrite);

	// Add all input nodes
	for (int i=0; i < inputnodes.size(); i++)
		pipeline.addNewInputEdge(inputnodes[i], myProcess);
	// Add all output nodes
	for (int i=0; i < outputnodes.size(); i++)
		pipeline.addNewOutputEdge(myProcess, outputnodes[i]);

	// Write the pipeline to an updated STAR file
	std::vector<bool> dummy;
	pipeline.write(dummy, dummy);

}

void RelionMainWindow::jobCommunicate(bool do_write, bool do_read, bool do_toggle_continue, bool do_commandline, bool do_makedir, int this_job)
{
	int itype = (this_job > 0) ? this_job : (browser->value() - 1); // browser starts counting at 1 ....
	show_initial_screen = false;

	// always write the general settings with the (hidden) empty name
	if (do_write)
		job_import->write("");

	if (is_main_continue && do_commandline)
	{
		// For the continuation jobs: set the outputname to the one from the current job
		global_outputname = fn_settings.beforeLastOf("/") + "/";
	}

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
			job_import->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_motioncorr->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_manualpick->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_ctffind->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_autopick->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_extract->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_sort->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_classselect->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_class2d->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_class3d->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_auto3d->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_polish->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_maskcreate->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_joinstar->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_subtract->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_post->getCommands(global_outputname, commands, final_command, do_makedir);
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
			job_resmap->getCommands(global_outputname, commands, final_command, do_makedir);
		break;
	}
	} // end switch

	// set the continue button correct upon reading of old settings
	if (do_read)
	{
		// Make the choice active
		cb_toggle_continue_i();
	}

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
    jobCommunicate(DONT_WRITE, DONT_READ, DO_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

	// Update all job lists in the main GUI
	updateJobLists();

	// If an new job was loaded from the job-browser: set the GUI to be a new job
    is_main_continue = false;
    cb_toggle_continue_i(); // make default active

    textbuff_current_job->text("New job");
	text_current_job->buffer(textbuff_current_job);
	//current_job = -1;

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
		cb_fill_stdout_i();
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
    	cb_fill_stdout_i();
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
			global_manualpickjob.read(fn_job.c_str(), iscont);
		else
			REPORT_ERROR("RelionMainWindow::cb_display_io_node_i ERROR: Save a Manual picking job parameters (using the File menu) before displaying coordinate files. ");


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
				std::cout << " Only coordinates generated in the pipeline are allowed" << std::endl;
			}
	    }
	    else
	    	std::cout << " Only coordinates in .star format (generated in the pipeline) are allowed" << std::endl;
	}
	else
	{
		command = "relion_display --gui --i " + pipeline.nodeList[mynode].name + " &";
	}
	int res= system(command.c_str());

}

void RelionMainWindow::cb_display(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_display_i();
}


void RelionMainWindow::cb_display_i()
{
        std::string command = " relion_display --gui &" ;
        system(command.c_str());
}

void RelionMainWindow::cb_toggle_continue_i()
{

	if (is_main_continue)
	{
		run_button->label("Continue now");
		run_button->labelfont(FL_ITALIC);
		run_button->labelsize(13);
		schedule_button->deactivate();
	}
	else
	{
		run_button->label("Run now!");
		run_button->labelfont(FL_ITALIC);
		run_button->labelsize(16);
		schedule_button->activate();
	}

	jobCommunicate(DONT_WRITE, DONT_READ, DO_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

}

void RelionMainWindow::cb_fill_stdout_i()
{

	FileName fn_out= pipeline.processList[current_job].name + "run.out";
	FileName fn_err= pipeline.processList[current_job].name + "run.err";

	if (exists(fn_out))
	{
		// Remove annoying carriage returns
		std::string command = "awk -F\"\r\" '{if (NF>1) {print $NF} else {print}}' < " + fn_out + " > .gui_tmpout";
		int res = system(command.c_str());
		std::ifstream in(".gui_tmpout", std::ios_base::in);
		if (in.fail())
			REPORT_ERROR( (std::string) "MetaDataTable::read: File " + fn_out + " does not exists" );
		int errno = textbuff_stdout->loadfile(".gui_tmpout");
		disp_stdout->scroll(textbuff_stdout->length(), 0);
		in.close();
	}
	else
		textbuff_stdout->text("stdout will go here");

	if (exists(fn_err))
	{
		std::ifstream in(fn_err.data(), std::ios_base::in);
		if (in.fail())
			REPORT_ERROR( (std::string) "MetaDataTable::read: File " + fn_err + " does not exists" );
		int errno = textbuff_stderr->loadfile(fn_err.c_str());
		disp_stderr->scroll(textbuff_stderr->length(), 0);
		in.close();
	}
	else
		textbuff_stderr->text("stderr will go here");

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
void RelionMainWindow::cb_run(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
	// Deactivate Run button to prevent the user from accidentally submitting many jobs
	run_button->deactivate();
	// Run the job
	T->cb_run_i();
}

void RelionMainWindow::cb_run_i()
{

	// Get the command line arguments from the currently active jobwindow,
	jobCommunicate(DONT_WRITE, DONT_READ, DONT_TOGGLE_CONT, DO_GET_CL, DO_MKDIR);

	// Save temporary hidden file with this jobs settings as default for a new job
	fn_settings = "";
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

	// Also save a copy of the GUI settings with the current output name
	fn_settings = global_outputname;
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

	if (commands.size()==0)
	{
		std::cout << " Nothing to do..."<< std::endl;
		return;
	}

	// If this is a continuation job, check whether output files exist and move away!
	// This is to ensure that the continuation job goes OK
	if (is_main_continue)
	{
		bool is_refine = (pipeline.processList[current_job].type == PROC_2DCLASS ||
				pipeline.processList[current_job].type == PROC_3DCLASS ||
				pipeline.processList[current_job].type == PROC_3DAUTO);
		if (!is_refine )
		{
			for (int i = 0; i < pipeline.processList[current_job].outputNodeList.size(); i++)
			{
				int j = pipeline.processList[current_job].outputNodeList[i];
				std::string fn_node = pipeline.nodeList[j].name;
				if (exists(fn_node))
				{
					std::string mvcommand = "mv -f " + fn_node + " " + fn_node + ".old";
					int res = system(mvcommand.c_str());
				}
			}
		}
	}

	std::cout << "Executing: " << final_command << std::endl;
	int res = system(final_command.c_str());

	// Now save the job (as Running) to the PipeLine
	addToPipeLine(PROC_RUNNING, is_main_continue); // is_main_continue means: only allow to overwrite an existing process when continuing an old run...


	// Update all job lists in the main GUI
	updateJobLists();


}

// Run button call-back functions
void RelionMainWindow::cb_schedule(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_schedule_i();
}

void RelionMainWindow::cb_schedule_i()
{

	// Get the command line arguments from the currently active jobwindow,
	jobCommunicate(DONT_WRITE, DONT_READ, DONT_TOGGLE_CONT, DO_GET_CL, DONT_MKDIR);

	// Save temporary hidden file with this jobs settings as default for a new job
	fn_settings = "";
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

	// Also save a copy of the GUI settings with the current output name
	// TODO: MOVE(?) scheduled jobs to a specific directory!
	fn_settings = ".ScheduledJobs/" + global_outputname;
	FileName command = "mkdir -p " + fn_settings;
	command = command.beforeLastOf("/");
	int res = system(command.c_str());

	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

	// Now save the job (as Scheduled) to the PipeLine
	addToPipeLine(PROC_SCHEDULED, true); // true means: allow to overwrite an existing process...

	// Update all job lists in the main GUI
	updateJobLists();

}

// Run button call-back functions
void RelionMainWindow::cb_run_scheduled(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_run_scheduled_i();
}


void RelionMainWindow::cb_run_scheduled_i()
{

	int scheduled_job = current_job;

	if (pipeline.processList[scheduled_job].status != PROC_SCHEDULED)
	{
		std::cout << " You can only perform this action on scheduled jobs ... " << std::endl;
		return;
	}

    fn_settings = pipeline.processList[scheduled_job].name;

    // Get the OLD UNIQDATE of the scheduled job!!
    FileName fn_olduniqdate, fn_newuniqdate, fn_newsettings;
    findUniqueDateSubstring(fn_settings, fn_olduniqdate);

    // Set a new uniqname for this job
    fn_newsettings = fn_settings.beforeFirstOf(fn_olduniqdate) + getUniqDateString() + "/";

    is_main_continue = false;
    cb_toggle_continue_i(); // make the continue=false active
	fn_settings = fn_newsettings;
    global_outputname = fn_newsettings;
    findUniqueDateSubstring(fn_newsettings, fn_newuniqdate);

	// Run the scheduled job
	cb_run_i();

	// Store run_job id (the last one in the list) for later on..
	int run_job = pipeline.processList.size() - 1;

	// Now replace all OLD UNIQDATEs in the inputNode names of remaining scheduled job with the new UNIQDATEs
	for (size_t i = 0; i < pipeline.processList.size(); i++)
	{
		if (pipeline.processList[i].status == PROC_SCHEDULED)
		{
			// Find the OLD UNIQDATE in the inputNodeList of all Scheduled jobs
			bool found = false;
			for (size_t j = 0; j < (pipeline.processList[i]).inputNodeList.size(); j++)
			{
				int mynode = (pipeline.processList[i]).inputNodeList[j];
				int pos = (pipeline.nodeList[mynode]).name.find(fn_olduniqdate);
				if (pos != std::string::npos)
				{
					// If found: search for the name with the NEW UNIQDATE in the NodeList and re-set this entry of the inputNodeList
					std::string newname = (pipeline.nodeList[mynode]).name;
					newname.replace(pos, 13, fn_newuniqdate);
					// Find the new name in the NodeList
					for (size_t ii = 0; ii < pipeline.nodeList.size(); ii++)
					{
						if (pipeline.nodeList[ii].name == newname)
						{
							// Set this node on the inputNodeList of the scheduled job
							(pipeline.processList[i]).inputNodeList[j] = ii;
							// Set the scheduled job in the inputForProcessList of this node
							pipeline.nodeList[ii].inputForProcessList.push_back(i);
							//break;
						}
					}
					// Keep track whether the .job file needs to be changed as well...
					found = true;
				}
			}
			if (found) // only open relevant .job files
			{
				/// TODO!!! Also change the input entries in the run.job file of that job!!!
			    FileName fn_job = ".ScheduledJobs/" + pipeline.processList[i].name + "run.job";
			    std::ifstream in(fn_job.data(), std::ios_base::in);
			    FileName fn_tmp = fn_job+".tmp";
			    std::ofstream out(fn_tmp.data(), std::ios_base::out);
			    if (in.fail())
			        REPORT_ERROR( (std::string) "RelionMainWindow::cb_run_scheduled_i: File " + fn_job + " does not exists" );
			    in.seekg(0);
			    std::string line;
			    while (getline(in, line, '\n'))
			    {
					int pos = line.find(fn_olduniqdate);
					if (pos < line.length())
					{
						line.replace(pos, 13, fn_newuniqdate);
					}
					out << line << std::endl;
			    }
			    in.close();
			    out.close();
			    std::rename(fn_tmp.data(), fn_job.data());
			    std::string command = "chmod 664 " + fn_job;
			    system(command.c_str());
			}
		}
	}

	// Delete the scheduled job from the pipeline
	current_job = scheduled_job;
	cb_delete_i(false, false); // don't ask for confirmation, don't delete recursively, as output of this job may still be used later on
	current_job = run_job;

	// Update all job lists in the main GUI
	updateJobLists();

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
			if ((pipeline.processList[idel]).outputNodeList.size() > 0)
			{
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
				ask += " - " + pipeline.processList[i].name + "\n";
		}
		proceed =  fl_choice(ask.c_str(), "Don't move", "Move", NULL);
	}
	else
	{
		proceed = 1;
	}
	if (proceed)
	{

		// Delete the output directories for all selected processes from the hard disk
		for (int i = 0; i < deleteProcesses.size(); i++)
		{
			if (deleteProcesses[i])
			{
				FileName alldirs = pipeline.processList[i].name;
				alldirs = alldirs.beforeLastOf("/");
				if (pipeline.processList[i].status == PROC_SCHEDULED)
				{
					// Just remove the .ScheduledJobs entry
					alldirs = ".ScheduledJobs/" + alldirs;
					std::string command = "rm -rf " + alldirs;
					int res = system(command.c_str());
				}
				else
				{
					// Move entire output directory (with subdirectory structure) to the Trash folder
					FileName firstdirs = alldirs.beforeLastOf("/");
					std::string command = "mkdir -p Trash/" + firstdirs;
					int res = system(command.c_str());
					command= "mv -f " + alldirs + " Trash/" + firstdirs+"/.";
					res = system(command.c_str());
				}

			}
		}

		// Write new pipeline to disc and read in again
		pipeline.write(deleteNodes, deleteProcesses);
		pipeline.read();

		// Update all job lists in the main GUI
		updateJobLists();
	}

}

// Run button call-back functions
void RelionMainWindow::cb_cleanup(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_cleanup_i();
}

void RelionMainWindow::cb_cleanup_i()
{

	if (current_job < 0)
	{
		std::cout << " You can only clean up existing jobs ... " << std::endl;
		return;
	}

	std::string ask;
	ask = "Are you sure you want to delete intermediate files from " + pipeline.processList[current_job].name + "?";
	int proceed = proceed =  fl_choice(ask.c_str(), "Don't delete", "Delete", NULL);
	if (proceed)
	{
	    std::cerr << "cleanup todo" << std::endl;
	}

}


// Run button call-back functions
void RelionMainWindow::cb_set_alias(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_set_alias_i();
}

void RelionMainWindow::cb_set_alias_i()
{
	std::string alias;
	FileName before_uniqdate, uniqdate;
	before_uniqdate = pipeline.processList[current_job].name;
	size_t slashpos = findUniqueDateSubstring(before_uniqdate, uniqdate);
	before_uniqdate = before_uniqdate.beforeFirstOf(uniqdate);

	bool is_done = false;
	while (!is_done)
	{
		const char * palias;
		palias =  fl_input("Rename to (provide 'None' to use original name): ", uniqdate.c_str());
		// TODO: check for cancel
		if (palias == NULL)
			return;

		std::string al2(palias);
		alias = al2;
		bool is_unique = true;
		for (size_t i = 0; i < pipeline.processList.size(); i++)
		{
			if ( pipeline.processList[i].alias == alias && alias != "None")
			{
				is_unique = false;
				break;
			}
		}
		if (!is_unique)
			 fl_message("Alias is not unique, please provide another one");
		else
			is_done = true;
	}

	if (alias == "None")
		pipeline.processList[current_job].alias = "None";
	else
	{
		// Make sure fn_alias ends with a slash
		if (alias[alias.length()-1] != '/')
			alias += "/";
		// Don't use same alias as process name!
		if (before_uniqdate + alias == pipeline.processList[current_job].name)
			pipeline.processList[current_job].alias = "None";
		else
			pipeline.processList[current_job].alias = before_uniqdate + alias;
	}

	// Write new pipeline to disc and read in again
	std::vector<bool> dummy;
	pipeline.write(dummy, dummy);
	pipeline.read();

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
		std::cout << " You can only mark existing jobs as finished ... " << std::endl;
		return;
	}

	pipeline.processList[current_job].status = PROC_FINISHED;

	// Write new pipeline to disc and read in again
	std::vector<bool> dummy;
	pipeline.write(dummy, dummy);
	pipeline.read();

	// Update all job lists in the main GUI
	updateJobLists();

}

void RelionMainWindow::cb_edit_note(Fl_Widget*, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_edit_note_i();

}

void RelionMainWindow::cb_edit_note_i()
{
	if (current_job < 0)
	{
		std::cout << " You can only edit the note for existing jobs ... " << std::endl;
		return;
	}
	FileName fn_note = pipeline.processList[current_job].name + "note.txt";
	NoteEditorWindow* w = new NoteEditorWindow(660, 400, "Note", fn_note);
	w->show();
}


// Save button call-back function
void RelionMainWindow::cb_menubar_save(Fl_Widget* o, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_save_i();
}

void RelionMainWindow::cb_menubar_save_i()
{
	// For scheduled jobs, also save the .job file in the .Scheduled directory
	if (pipeline.processList[current_job].status == PROC_SCHEDULED)
	{
		fn_settings = ".ScheduledJobs/" + global_outputname;
		jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);
	}

	fn_settings = "";
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

}
void RelionMainWindow::cb_menubar_print_notes(Fl_Widget*, void* v)
{
  	 RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_print_notes_i();
}

void RelionMainWindow::cb_menubar_print_notes_i()
{
	std::cout << " ################################################################ " << std::endl;
	std::cout << " # Printing all note files for pipeline: " << pipeline.name << std::endl;
	for (size_t i = 0; i < pipeline.processList.size(); i++)
	{
		FileName fn_note = pipeline.processList[i].name+"note.txt";
		std::cout << " ################################################################ " << std::endl;
		std::cout << " # Job= " << pipeline.processList[i].name;
		if (pipeline.processList[i].alias != "None")
			std::cout <<" alias: " << pipeline.processList[i].alias;
		std::cout	<< std::endl;
		if (exists(fn_note))
		{
			std::ifstream in(fn_note.data(), std::ios_base::in);
			std::string line;
			if (in.fail())
        		REPORT_ERROR( (std::string) "ERROR: cannot read file " + fn_note);
    	    in.seekg(0);
    	    while (getline(in, line, '\n'))
    	    {
    	    	std::cout << line << std::endl;
    	    }
			in.close();
		}
	}

}


void RelionMainWindow::cb_menubar_reactivate_runbutton(Fl_Widget* o, void* v)
{

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_reactivate_runbutton_i();
}

void RelionMainWindow::cb_menubar_reactivate_runbutton_i()
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

void RelionMainWindow::cb_menubar_run_scheduled_jobs(Fl_Widget* o, void* v)
{

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_run_scheduled_jobs_i();
}

void RelionMainWindow::cb_menubar_run_scheduled_jobs_i()
{
	// TODO
}

void RelionMainWindow::cb_menubar_stop_run_scheduled_jobs(Fl_Widget* o, void* v)
{

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_stop_run_scheduled_jobs_i();
}

void RelionMainWindow::cb_menubar_stop_run_scheduled_jobs_i()
{
	// TODO
	std::cout <<" todo.." << std::endl;
}

void RelionMainWindow::cb_menubar_about(Fl_Widget* o, void* v)
{

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_about_i();
}

void RelionMainWindow::cb_menubar_about_i()
{
	ShowHelpText *help = new ShowHelpText("\
RELION is written by Sjors Scheres at the MRC Laboratory of Molecular Biology (scheres@mrc-lmb.cam.ac.uk).\n \
\n\
If RELION is useful in your work, please cite it in the contexts as explained under the \"Publish!\" tab, or on the RELION wiki at http://www2.mrc-lmb.cam.ac.uk/relion. \n  \
\n\
Note that RELION is completely free, open-source software. You can redistribute it and/or modify it for your own purposes, but please do make sure \
the contribution of Sjors Scheres is acknowledged appropriately. In order to maintain an overview of existing versions, he would also appreciate being \
notified of any redistribution of (modified versions of) the code. \n \
");
}


void RelionMainWindow::cb_menubar_quit(Fl_Widget* o, void* v)
{
	RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_quit_i();
}

void RelionMainWindow::cb_menubar_quit_i()
{
	exit(0);
}
