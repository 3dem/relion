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

RelionMainWindow::RelionMainWindow(int w, int h, const char* title, FileName fn_pipe):Fl_Window(w,h,title)
{

	// First setup the old part of the GUI
	h = GUIHEIGHT_OLD;

	// Read in the pipeline STAR file if it exists
	pipeline.name = fn_pipe;
	if (exists(fn_pipe + "_pipeline.star"))
		pipeline.read();

	// Check which jobs have finished
	pipeline.checkProcessCompletion();

	// Make temporary directory with all NodeNames
	pipeline.makeNodeDirectory();

	// Initialisation
	run_button = NULL;
	print_CL_button = NULL;
	cite_button = NULL;

    color(GUI_BACKGROUND_COLOR);
    menubar = new Fl_Menu_Bar(0, 0, w, MENUHEIGHT);
    menubar->add("File/Load settings",  FL_ALT+'l', cb_menubar_load, this);
    menubar->add("File/Save settings",  FL_ALT+'s', cb_menubar_save, this);
    menubar->add("File/Reactivate Run",  FL_ALT+'r', cb_menubar_reactivate_runbutton, this);
    menubar->add("File/About",  FL_ALT+'a', cb_menubar_about, this);
    menubar->add("File/Quit", FL_ALT+'q', cb_menubar_quit, this);
    current_y = MENUHEIGHT + 10;

    toggle_continue = new Fl_Toggle_Button(WCOL0, 4, 200, 22);
    toggle_continue->labelsize(14);
    toggle_continue->color(GUI_BACKGROUND_COLOR, GUI_BACKGROUND_COLOR);
    toggle_continue->callback(cb_toggle_continue, this);

    // Add run buttons on the menubar as well
	print_CL_button = new Fl_Button(GUIWIDTH - 220, h-50, 100, 30, "Print command");
	print_CL_button->color(GUI_RUNBUTTON_COLOR);
	print_CL_button->labelsize(12);
	print_CL_button->callback( cb_print_cl, this);

	run_button = new Fl_Button(GUIWIDTH - 110 , h-50, 100, 30, "Run!");
	run_button->color(GUI_RUNBUTTON_COLOR);
	run_button->labelfont(FL_ITALIC);
	run_button->labelsize(18);
	run_button->callback( cb_run, this);


    display_button = new Fl_Button(10, h-50, XCOL0-20, 30, "Display");
	display_button->color(GUI_RUNBUTTON_COLOR);
	display_button->callback( cb_display, this);

	// Browser to act as "tab selector"
    browser = new Fl_Hold_Browser(10,MENUHEIGHT+10,WCOL0-20,h-MENUHEIGHT-70);

    // Fill browser in the right order
    for (int itype = 1; itype <= NR_BROWSE_TABS; itype++)
    {

        browse_grp[itype-1] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	switch (itype)
    	{
    	case PROC_GENERAL:
    	{
            browser->add("General");
        	job_general = new GeneralJobWindow();
    		break;
    	}
    	case PROC_MANUALPICK:
    	{
    		browser->add("Micrograph inspection");
        	job_manualpick = new ManualpickJobWindow();
    		break;
    	}
    	case PROC_CTFFIND:
    	{
			browser->add("CTF estimation");
			job_ctffind = new CtffindJobWindow();
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
    	case PROC_POST:
		{
			browser->add("Post-processing");
			job_post = new PostJobWindow();
			break;
		}
    	case PROC_RESMAP:
		{
			browser->add("Local-resolution");
			job_resmap = new ResmapJobWindow();
			break;
		}
    	case PROC_PUBLISH:
		{
			browser->add("Publish!");
			job_publish = new PublishJobWindow();
			break;
		}
    	default:
    	{
           	break;
    	}
    	} // end switch
    	browse_grp[itype-1]->end();
    }

    browser->callback(cb_select_browsegroup);
    browser->end();

    // Set and activate current selection
    browser->select(1); // just start from the beginning
    cb_select_browsegroup_i(); // make default active
    toggle_continue->value(0); // 0 = new run; 1 = continue
    cb_toggle_continue_i(); // make default active


    // Pipeline part of the GUI
#define JOBCOLWIDTH (250)
#define XJOBCOL1 (10)
#define XJOBCOL2 (JOBCOLWIDTH + 25)
#define XJOBCOL3 (2*JOBCOLWIDTH + 40)


    // Add browsers for finished, running and scheduled jobs
    Fl_Text_Buffer *textbuff1 = new Fl_Text_Buffer();
    Fl_Text_Buffer *textbuff2 = new Fl_Text_Buffer();
    Fl_Text_Buffer *textbuff3 = new Fl_Text_Buffer();
    textbuff1->text("Finished jobs");
    textbuff2->text("Running jobs");
    textbuff3->text("Scheduled jobs");
    Fl_Text_Display* textdisp1 = new Fl_Text_Display(XJOBCOL1, GUIHEIGHT_OLD, JOBCOLWIDTH, 25);
	Fl_Text_Display* textdisp2 = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_OLD, JOBCOLWIDTH, 25);
	Fl_Text_Display* textdisp3 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_OLD, JOBCOLWIDTH, 25);
	textdisp1->buffer(textbuff1);
	textdisp2->buffer(textbuff2);
	textdisp3->buffer(textbuff3);
	textdisp1->color(GUI_BACKGROUND_COLOR);
	textdisp2->color(GUI_BACKGROUND_COLOR);
	textdisp3->color(GUI_BACKGROUND_COLOR);

    finished_job_browser  = new Fl_Select_Browser(XJOBCOL1, GUIHEIGHT_OLD+25, JOBCOLWIDTH, GUIHEIGHT_EXT-GUIHEIGHT_OLD-MENUHEIGHT-30);
    running_job_browser   = new Fl_Select_Browser(XJOBCOL2, GUIHEIGHT_OLD+25, JOBCOLWIDTH, GUIHEIGHT_EXT-GUIHEIGHT_OLD-MENUHEIGHT-30);
    scheduled_job_browser = new Fl_Select_Browser(XJOBCOL3, GUIHEIGHT_OLD+25, JOBCOLWIDTH, GUIHEIGHT_EXT-GUIHEIGHT_OLD-MENUHEIGHT-30);

    // Fill the Jobs browsers
    // Search backwards, so that last jobs are at the top
    for (long int i = pipeline.processList.size() -1, jobnr=0 ; i >= 0; i--, jobnr++)
    {
    	int itype = pipeline.processList[i].type;

    	//job_browse_grp[jobnr] = new Fl_Group(WCOL0, GUIHEIGHT_OLD, 550, GUIHEIGHT_EXT-GUIHEIGHT_OLD);

    	switch (pipeline.processList[i].status)
    	{
    	case PROC_FINISHED:
    	{

    		finished_processes.push_back(i);
    		finished_job_browser->add(pipeline.processList[i].name.c_str());
    		break;
    	}
    	case PROC_RUNNING:
    	{
    		running_processes.push_back(i);
    		running_job_browser->add(pipeline.processList[i].name.c_str());
    		break;
    	}
    	case PROC_SCHEDULED:
    	{
    		scheduled_processes.push_back(i);
    		scheduled_job_browser->add(pipeline.processList[i].name.c_str());
    		break;
    	}
    	default:
    	{
    		REPORT_ERROR("RelionMainWindow::RelionMainWindow ERROR: unrecognised status for job" + pipeline.processList[i].name);
    	}
    	}//end switch

		//job_browse_grp[jobnr]->end();
    }

    // Set the callbacks
    finished_job_browser->callback(cb_select_finished_job);
    running_job_browser->callback(cb_select_running_job);
    scheduled_job_browser->callback(cb_select_scheduled_job);

    finished_job_browser->end();
    running_job_browser->end();
    scheduled_job_browser->end();

    resizable();
}

void RelionMainWindow::addToPipeLine(int as_status, bool do_overwrite, int this_job)
{
	int itype = (this_job > 0) ? this_job : browser->value();
	std::vector<Node> inputnodes;
	std::vector<Node> outputnodes;
	std::string outputname;

	switch (itype)
	{
	case PROC_GENERAL:
	{
		inputnodes = job_general->pipelineInputNodes;
		outputnodes= job_general->pipelineOutputNodes;
		outputname = job_general->pipelineOutputName;
		break;
	}
	case PROC_MANUALPICK:
	{
		inputnodes = job_manualpick->pipelineInputNodes;
		outputnodes= job_manualpick->pipelineOutputNodes;
		outputname = job_manualpick->pipelineOutputName;
		break;
	}
	case PROC_CTFFIND:
	{
		inputnodes = job_ctffind->pipelineInputNodes;
		outputnodes= job_ctffind->pipelineOutputNodes;
		outputname = job_ctffind->pipelineOutputName;
		break;
	}
	case PROC_AUTOPICK:
	{
		inputnodes = job_autopick->pipelineInputNodes;
		outputnodes= job_autopick->pipelineOutputNodes;
		outputname = job_autopick->pipelineOutputName;
		break;
	}
	case PROC_EXTRACT:
	{
		inputnodes = job_extract->pipelineInputNodes;
		outputnodes= job_extract->pipelineOutputNodes;
		outputname = job_extract->pipelineOutputName;
		break;
	}
	case PROC_SORT:
	{
		inputnodes = job_sort->pipelineInputNodes;
		outputnodes= job_sort->pipelineOutputNodes;
		outputname = job_sort->pipelineOutputName;
		break;
	}
	case PROC_2DCLASS:
	{

		inputnodes = job_class2d->pipelineInputNodes;
		outputnodes= job_class2d->pipelineOutputNodes;
		outputname = job_class2d->pipelineOutputName;
		break;
	}
	case PROC_3DCLASS:
	{

		inputnodes = job_class3d->pipelineInputNodes;
		outputnodes= job_class3d->pipelineOutputNodes;
		outputname = job_class3d->pipelineOutputName;
		break;
	}
	case PROC_3DAUTO:
	{

		inputnodes = job_auto3d->pipelineInputNodes;
		outputnodes= job_auto3d->pipelineOutputNodes;
		outputname = job_auto3d->pipelineOutputName;
		break;
	}
	case PROC_POLISH:
	{

		inputnodes = job_polish->pipelineInputNodes;
		outputnodes= job_polish->pipelineOutputNodes;
		outputname = job_polish->pipelineOutputName;
		break;
	}
	case PROC_POST:
	{

		inputnodes = job_post->pipelineInputNodes;
		outputnodes= job_post->pipelineOutputNodes;
		outputname = job_post->pipelineOutputName;
		break;
	}
	case PROC_RESMAP:
	{

		inputnodes = job_resmap->pipelineInputNodes;
		outputnodes= job_resmap->pipelineOutputNodes;
		outputname = job_resmap->pipelineOutputName;
		break;
	}
	default:
	{
		REPORT_ERROR("RelionMainWindow::addToPipeLine ERROR: unrecognised type");
	}
	}


	// Add Process to the processList of the pipeline
	Process process(outputname, itype, as_status);
	long int myProcess = pipeline.addNewProcess(process, do_overwrite);
	// Add all input nodes
	for (int i=0; i < inputnodes.size(); i++)
		pipeline.addNewInputEdge(inputnodes[i], myProcess);
	// Add all output nodes
	for (int i=0; i < outputnodes.size(); i++)
		pipeline.addNewOutputEdge(myProcess, outputnodes[i]);

	// Write the pipeline to an updated STAR file
	pipeline.write();

}

void RelionMainWindow::jobCommunicate(bool do_write, bool do_read, bool do_toggle_continue, bool do_commandline, int this_job)
{
	int itype = (this_job > 0) ? this_job : browser->value();

	// always write the general settings with the (hidden) empty name
	if (do_write)
		job_general->write("");

	switch (itype)
	{
	case PROC_GENERAL:
	{
		if (do_write)
			job_general->write(fn_settings);
		if (do_read)
			job_general->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_general->toggle_new_continue(is_main_continue);
		if (do_commandline)
			job_general->getCommands(outputname, commands, final_command);
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
			job_manualpick->getCommands(outputname, commands, final_command,
					job_general->angpix.getValue(), job_general->particle_diameter.getValue());
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
			job_ctffind->getCommands(outputname, commands, final_command, job_general->angpix.getValue());
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
			job_autopick->getCommands(outputname, commands, final_command,
					job_general->angpix.getValue(), job_general->particle_diameter.getValue());
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
			job_extract->getCommands(outputname, commands, final_command,
					job_general->angpix.getValue(), job_general->particle_diameter.getValue());
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
			job_sort->getCommands(outputname, commands, final_command,
					job_general->angpix.getValue(), job_general->particle_diameter.getValue());
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
			job_class2d->getCommands(outputname, commands, final_command,
					job_general->angpix.getValue(), job_general->particle_diameter.getValue());
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
			job_class3d->getCommands(outputname, commands, final_command,
					job_general->angpix.getValue(), job_general->particle_diameter.getValue());
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
			job_auto3d->getCommands(outputname, commands, final_command,
					job_general->angpix.getValue(), job_general->particle_diameter.getValue());
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
			job_polish->getCommands(outputname, commands, final_command,
					job_general->angpix.getValue(), job_general->particle_diameter.getValue(),
					job_extract->black_dust.getValue(), job_extract->white_dust.getValue());
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
			job_post->getCommands(outputname, commands, final_command,
					job_general->angpix.getValue());
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
			job_resmap->getCommands(outputname, commands, final_command,
					job_general->angpix.getValue());
		break;
	}
	case PROC_PUBLISH:
	{
		if (do_toggle_continue)
			job_publish->toggle_new_continue(is_main_continue);
		if (do_commandline)
			job_publish->getCommands(outputname, commands, final_command);
		break;
	}
	} // end switch

	// set the continue button correct upon reading of old settings
	if (do_read)
	{
		if (is_main_continue)
			toggle_continue->value(1);
		else
			toggle_continue->value(0);
		// Make the choice active
		cb_toggle_continue_i();
	}


}

void RelionMainWindow::cb_select_browsegroup(Fl_Widget* o, void* v)
{
	RelionMainWindow* T=(RelionMainWindow*)v;
	T->cb_select_browsegroup_i();
	run_button->activate();

}

void RelionMainWindow::cb_select_browsegroup_i()
{

	// Show the 'selected' group, hide the others
    for ( int t=1; t<=NR_BROWSE_TABS; t++ )
    {
    	if ( t == (browser->value()) )
        {
        	browse_grp[t-1]->show();
        }
        else
        {
        	browse_grp[t-1]->hide();
        }
    }

    // Toggle the new tab according to the continue toggle button
    jobCommunicate(DONT_WRITE, DONT_READ, DO_TOGGLE_CONT, DONT_GET_CL);

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
		int itype = pipeline.processList[finished_processes[idx]].type;
		std::cerr << "itype= " << itype << std::endl;

		for ( int t=1; t<=NR_BROWSE_TABS; t++ )
		{
			if ( t == itype )
			{
				browse_grp[t-1]->show();
				browser->select(t);
			}
			else
			{
				browse_grp[t-1]->hide();
			}
		}

		fn_settings = pipeline.processList[finished_processes[idx]].name;

		// Re-read the settings for this job
		cb_select_browsegroup_i(); // change to the corresponding jobwindow
		jobCommunicate(DONT_WRITE, DO_READ, DONT_TOGGLE_CONT, DONT_GET_CL);
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
		int itype = pipeline.processList[running_processes[idx]].type;

		for ( int t=1; t<=NR_BROWSE_TABS; t++ )
		{
			if ( t == itype )
			{
				browse_grp[t-1]->show();
				browser->select(t);
			}
			else
			{
				browse_grp[t-1]->hide();
			}
		}
		fn_settings = pipeline.processList[running_processes[idx]].name;

		// Re-read the settings for this job
		cb_select_browsegroup_i(); // change to the corresponding jobwindow
		jobCommunicate(DONT_WRITE, DO_READ, DONT_TOGGLE_CONT, DONT_GET_CL);
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
		int itype = pipeline.processList[scheduled_processes[idx]].type;

		for ( int t=1; t<=NR_BROWSE_TABS; t++ )
		{
			if ( t == itype )
			{
				browse_grp[t-1]->show();
				browser->select(t);
			}
			else
			{
				browse_grp[t-1]->hide();
			}
		}

		fn_settings = pipeline.processList[scheduled_processes[idx]].name;
		// Re-read the settings for this job
		cb_select_browsegroup_i(); // change to the corresponding jobwindow
		jobCommunicate(DONT_WRITE, DO_READ, DONT_TOGGLE_CONT, DONT_GET_CL);
    }
}

// Display button call-back functions
void RelionMainWindow::cb_display(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_display_i();
}


void RelionMainWindow::cb_display_i()
{
	std::string command = " relion_display --gui &" ;
	int res = system(command.c_str());
}

void RelionMainWindow::cb_toggle_continue(Fl_Widget*, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_toggle_continue_i();
}

void RelionMainWindow::cb_toggle_continue_i()
{

	if (toggle_continue->value() == 1)
	{
		is_main_continue = true;
		toggle_continue->label("Continue old run");
	}
	else
	{
		is_main_continue = false;
		toggle_continue->label("Start new run");
	}

	jobCommunicate(DONT_WRITE, DONT_READ, DO_TOGGLE_CONT, DONT_GET_CL);

}

void RelionMainWindow::cb_print_cl(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_print_cl_i();
}

void RelionMainWindow::cb_print_cl_i()
{
	jobCommunicate(DONT_WRITE, DONT_READ, DONT_TOGGLE_CONT, DO_GET_CL);
    std::cout << " *** The command is:" << std::endl;
    for (int icom = 0; icom < commands.size(); icom++)
    	std::cout << commands[icom] << std::endl;

}

// Run button call-back functions
void RelionMainWindow::cb_run(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_run_i();
}

void RelionMainWindow::cb_run_i()
{

	// Get the command line arguments from the currently active jobwindow,
	jobCommunicate(DONT_WRITE, DONT_READ, DONT_TOGGLE_CONT, DO_GET_CL);

	// Save temporary hidden file with this jobs settings as default for a new job
	fn_settings = "";
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL);

	// Also save a copy of the GUI settings with the current output name
	fn_settings = outputname;
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL);

	if (commands.size()==0)
	{
		std::cout << " Nothing to do..."<< std::endl;
		return;
	}

	std::cout << "Executing: " << final_command << std::endl;
	int res = system(final_command.c_str());

	// Now save the job (as Running) to the PipeLine
	addToPipeLine(PROC_RUNNING, false);

	// Deactivate Run button to prevent the user from accidentally submitting many jobs
	run_button->deactivate();

}


// call-back functions for the menubar
void RelionMainWindow::cb_menubar_load(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_load_i();
}

void RelionMainWindow::cb_menubar_load_i()
{
    Fl_File_Chooser * G_chooser = new Fl_File_Chooser("", "*.settings", Fl_File_Chooser::SINGLE, "Choose a GUI settings file");

    G_chooser->directory(NULL);
    G_chooser->show();

    // Block until user picks something.
    //     (The other way to do this is to use a callback())
    //
    while(G_chooser->shown()) {
        Fl::wait();
    }

    // Print the results
    if ( G_chooser->value() == NULL ) {
        //fprintf(stderr, "(User hit 'Cancel')\n");
        return;
    }
    // Read in settings file with the chosen name
    fn_settings = G_chooser->value();

    // Set the jobtype to the correct one
    // deduce from the settings name
    int npos = fn_settings.length();
    	if (fn_settings.rfind(".gui_general.") < npos)
    		browser->value(1);
    	else if  (fn_settings.rfind(".gui_manualpick.") < npos)
    		browser->value(2);
    	else if  (fn_settings.rfind(".gui_ctffind.") < npos)
    		browser->value(3);
    	else if  (fn_settings.rfind(".gui_autopick.") < npos)
    		browser->value(4);
    	else if  (fn_settings.rfind(".gui_extract.") < npos)
    		browser->value(5);
    	else if  (fn_settings.rfind(".gui_sort.") < npos)
    		browser->value(6);
    	else if  (fn_settings.rfind(".gui_class2d.") < npos)
    		browser->value(7);
    	else if  (fn_settings.rfind(".gui_class3d.") < npos)
    		browser->value(8);
    	else if  (fn_settings.rfind(".gui_auto3d.") < npos)
    		browser->value(9);
    	else if  (fn_settings.rfind(".gui_polish.") < npos)
    		browser->value(10);
    	else if  (fn_settings.rfind(".gui_post.") < npos)
    		browser->value(11);
    	else if  (fn_settings.rfind(".gui_resmap.") < npos)
    		browser->value(12);
    cb_select_browsegroup_i(); // change to the corresponding jobwindow
    // Read in the settings file
    jobCommunicate(DONT_WRITE, DO_READ, DONT_TOGGLE_CONT, DONT_GET_CL);

}

// Save button call-back function
void RelionMainWindow::cb_menubar_save(Fl_Widget* o, void* v)
{
    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_menubar_save_i();
}

void RelionMainWindow::cb_menubar_save_i()
{
	fn_settings = "";
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL);

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
