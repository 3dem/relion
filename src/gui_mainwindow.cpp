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
    menubar->add("File/Save job settings",  FL_ALT+'s', cb_menubar_save, this);
    menubar->add("File/Reactivate Run",  FL_ALT+'r', cb_menubar_reactivate_runbutton, this);
    menubar->add("File/About",  FL_ALT+'a', cb_menubar_about, this);
    menubar->add("File/Quit", FL_ALT+'q', cb_menubar_quit, this);
    current_y = MENUHEIGHT + 10;

    toggle_continue = new Fl_Toggle_Button(WCOL0, 4, 200, 22);
    toggle_continue->labelsize(14);
    toggle_continue->color(GUI_BACKGROUND_COLOR, GUI_BACKGROUND_COLOR);
    toggle_continue->callback(cb_toggle_continue, this);

    // Add run buttons on the menubar as well
	print_CL_button = new Fl_Button(GUIWIDTH - 330, h-50, 100, 30, "Print command");
	print_CL_button->color(GUI_RUNBUTTON_COLOR);
	print_CL_button->labelsize(12);
	print_CL_button->callback( cb_print_cl, this);

	schedule_button = new Fl_Button(GUIWIDTH - 220 , h-50, 100, 30, "Schedule");
	schedule_button->color(GUI_RUNBUTTON_COLOR);
	schedule_button->labelfont(FL_ITALIC);
	schedule_button->labelsize(16);
	schedule_button->callback( cb_schedule, this);

	run_button = new Fl_Button(GUIWIDTH - 110 , h-50, 100, 30, "Run now");
	run_button->color(GUI_RUNBUTTON_COLOR);
	run_button->labelfont(FL_ITALIC);
	run_button->labelsize(16);
	run_button->callback( cb_run, this);


    display_button = new Fl_Button(10, h-50, XCOL0-20, 30, "Display");
	display_button->color(GUI_RUNBUTTON_COLOR);
	display_button->callback( cb_display, this);

	// Browser to act as "tab selector"
    browser = new Fl_Hold_Browser(10,MENUHEIGHT+10,WCOL0-20,h-MENUHEIGHT-70);

    // Fill browser in the right order
    current_job = -1;
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


    std::cerr << " end filling browse"<<std::endl;

    // Pipeline part of the GUI
#define JOBCOLWIDTH (250)
#define XJOBCOL1 (10)
#define XJOBCOL2 (JOBCOLWIDTH + 25)
#define XJOBCOL3 (2*JOBCOLWIDTH + 40)


    // Add browsers for finished, running and scheduled jobs
    int JOBHEIGHT = GUIHEIGHT_EXT-GUIHEIGHT_OLD-MENUHEIGHT-30-50;
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
    Fl_Text_Display* textdisp1 = new Fl_Text_Display(XJOBCOL1, GUIHEIGHT_OLD, JOBCOLWIDTH, 25);
	Fl_Text_Display* textdisp2 = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_OLD, JOBCOLWIDTH, 25);
	Fl_Text_Display* textdisp3 = new Fl_Text_Display(XJOBCOL2, GUIHEIGHT_OLD+JOBHALFHEIGHT+25, JOBCOLWIDTH, 25);
	Fl_Text_Display* textdisp4 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_OLD, JOBCOLWIDTH, 25);
	Fl_Text_Display* textdisp5 = new Fl_Text_Display(XJOBCOL3, GUIHEIGHT_OLD+JOBHALFHEIGHT+25, JOBCOLWIDTH, 25);
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

    finished_job_browser  = new Fl_Select_Browser(XJOBCOL1, GUIHEIGHT_OLD+25, JOBCOLWIDTH, JOBHEIGHT+25);
    running_job_browser   = new Fl_Select_Browser(XJOBCOL2, GUIHEIGHT_OLD+25, JOBCOLWIDTH, JOBHALFHEIGHT);
    scheduled_job_browser = new Fl_Select_Browser(XJOBCOL2, GUIHEIGHT_OLD+25+JOBHALFHEIGHT+25, JOBCOLWIDTH, JOBHALFHEIGHT);
    this_from_job_browser = new Fl_Select_Browser(XJOBCOL3, GUIHEIGHT_OLD+25, JOBCOLWIDTH, JOBHALFHEIGHT);
    this_to_job_browser   = new Fl_Select_Browser(XJOBCOL3, GUIHEIGHT_OLD+25+JOBHALFHEIGHT+25, JOBCOLWIDTH, JOBHALFHEIGHT);

    // Fill the actual browsers
    fillRunningJobLists();

    // Set the callbacks
    finished_job_browser->callback(cb_select_finished_job);
    running_job_browser->callback(cb_select_running_job);
    scheduled_job_browser->callback(cb_select_scheduled_job);
    this_from_job_browser->callback(cb_select_from_job);
    this_to_job_browser->callback(cb_select_to_job);

    finished_job_browser->end();
    running_job_browser->end();
    scheduled_job_browser->end();
    this_from_job_browser->end();
    this_to_job_browser->end();

    delete_button = new Fl_Button(XJOBCOL1 , GUIHEIGHT_EXT-50, 100, 30, "Delete");
    delete_button->color(GUI_RUNBUTTON_COLOR);
    delete_button->labelfont(FL_ITALIC);
    delete_button->labelsize(16);
    delete_button->callback( cb_delete, this);

    cleanup_button = new Fl_Button(130 , GUIHEIGHT_EXT-50, 100, 30, "Clean up");
    cleanup_button->color(GUI_RUNBUTTON_COLOR);
    cleanup_button->labelfont(FL_ITALIC);
    cleanup_button->labelsize(16);
    cleanup_button->callback( cb_cleanup, this);

    run_scheduled_button = new Fl_Button(XJOBCOL2 , GUIHEIGHT_EXT-50, 100, 30, "Run scheduled");
    run_scheduled_button->color(GUI_RUNBUTTON_COLOR);
    run_scheduled_button->labelfont(FL_ITALIC);
    run_scheduled_button->labelsize(12);
    run_scheduled_button->callback( cb_run_scheduled, this);

    // Set and activate current selection from side-browser
    browser->select(1); // just start from the beginning
    cb_select_browsegroup_i(); // make default active
    toggle_continue->value(0); // 0 = new run; 1 = continue
    cb_toggle_continue_i(); // make default active

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

    }

}

void RelionMainWindow::fillToAndFromJobLists()
{
	this_to_job_browser->clear();
	this_from_job_browser->clear();
	this_to_processes.clear();
	this_from_processes.clear();

	if (current_job >= 0)
	{
		// Where do the input nodes come from?
		for (long int inode = 0; inode < (pipeline.processList[current_job]).inputNodeList.size(); inode++)
		{
			long int mynode = (pipeline.processList[current_job]).inputNodeList[inode];
			long int myproc = (pipeline.nodeList[mynode]).outputFromProcess;
			if (myproc >= 0)
			{
				// Check if this process was already there
				bool already_there = false;
				for (long int i = 0; i < this_from_processes.size(); i++)
				{
					if (myproc == this_from_processes[i])
					{
						already_there=true;
						break;
					}
				}
				if (!already_there)
				{
					this_from_processes.push_back(myproc);
					this_from_job_browser->add(pipeline.processList[myproc].name.c_str());
				}
			}
		}
		// Where do the output nodes lead to?
		for (long int inode = 0; inode < (pipeline.processList[current_job]).outputNodeList.size(); inode++)
		{
			long int mynode = (pipeline.processList[current_job]).outputNodeList[inode];
			long int nr_outputs = (pipeline.nodeList[mynode]).inputForProcessList.size();
			for (long int iproc = 0; iproc < nr_outputs; iproc++)
			{
				long int myproc =  (pipeline.nodeList[mynode]).inputForProcessList[iproc];
				// Check if this process was already there
				bool already_there = false;
				for (long int i = 0; i < this_to_processes.size(); i++)
				{
					if (myproc == this_to_processes[i])
					{
						already_there=true;
						break;
					}
				}
				if (!already_there)
				{
					this_to_processes.push_back(myproc);
					this_to_job_browser->add(pipeline.processList[myproc].name.c_str());
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

	int itype = pipeline.processList[current_job].type;
	fn_settings = pipeline.processList[current_job].name;

	if (pipeline.processList[current_job].status == PROC_SCHEDULED)
		fn_settings = ".ScheduledJobs/" + fn_settings;
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

	// Re-read the settings for this job
	cb_select_browsegroup_i(false); // change to the corresponding jobwindow
	jobCommunicate(DONT_WRITE, DO_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

	// Update all job lists in the main GUI
	updateJobLists();

}

std::string RelionMainWindow::findUniqueDateSubstring(std::string fnt)
{
	bool found = false;
	size_t slashpos = 0;
	int i = 0;
	while (slashpos < fnt.length() && i < 5)
	{
		i++;
		slashpos = fnt.find("/", slashpos+1);
		if (std::isdigit(fnt[slashpos+1]) && std::isdigit(fnt[slashpos+2]) && std::isdigit(fnt[slashpos+3]) &&
		    std::isdigit(fnt[slashpos+4]) && std::isdigit(fnt[slashpos+5]) && std::isdigit(fnt[slashpos+6]) &&
		    fnt[slashpos+7] == '-' &&
		    std::isdigit(fnt[slashpos+8]) && std::isdigit(fnt[slashpos+9]) && std::isdigit(fnt[slashpos+10]) &&
		    std::isdigit(fnt[slashpos+11]) && std::isdigit(fnt[slashpos+12]) && std::isdigit(fnt[slashpos+13]) )
			return fnt.substr(slashpos+1,13);
	}

	REPORT_ERROR("RelionMainWindow::findUniqueDateSubstring ERROR: no uniq-date identifier found for: " + fnt);
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
	std::vector<bool> dummy;
	pipeline.write(dummy, dummy);

}

void RelionMainWindow::jobCommunicate(bool do_write, bool do_read, bool do_toggle_continue, bool do_commandline, bool do_makedir, int this_job)
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
					job_general->angpix.getValue(), job_general->particle_diameter.getValue(), do_makedir);
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
			job_ctffind->getCommands(outputname, commands, final_command, job_general->angpix.getValue(), do_makedir);
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
					job_general->angpix.getValue(), job_general->particle_diameter.getValue(), do_makedir);
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
					job_general->angpix.getValue(), job_general->particle_diameter.getValue(), do_makedir);
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
					job_general->angpix.getValue(), job_general->particle_diameter.getValue(), do_makedir);
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
					job_general->angpix.getValue(), job_general->particle_diameter.getValue(), do_makedir);
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
					job_general->angpix.getValue(), job_general->particle_diameter.getValue(), do_makedir);
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
					job_general->angpix.getValue(), job_general->particle_diameter.getValue(), do_makedir);
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
					job_extract->black_dust.getValue(), job_extract->white_dust.getValue(), do_makedir);
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
					job_general->angpix.getValue(), do_makedir);
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
					job_general->angpix.getValue(), do_makedir);
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

void RelionMainWindow::cb_select_browsegroup_i(bool change_current_job)
{

	// When filling in a new form, set the current_job to -1
	if (change_current_job)
		current_job = -1;

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
    jobCommunicate(DONT_WRITE, DONT_READ, DO_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

	// Update all job lists in the main GUI
	updateJobLists();
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

void RelionMainWindow::cb_select_from_job(Fl_Widget* o, void* v)
{
	RelionMainWindow* T=(RelionMainWindow*)v;
	T->cb_select_from_job_i();
	run_button->activate();
}

void RelionMainWindow::cb_select_from_job_i()
{
	// Show the 'selected' group, hide the others
    int idx = this_from_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
    {
    	current_job = this_from_processes[idx];
		loadJobFromPipeline();
    }
}

void RelionMainWindow::cb_select_to_job(Fl_Widget* o, void* v)
{
	RelionMainWindow* T=(RelionMainWindow*)v;
	T->cb_select_to_job_i();
	run_button->activate();
}

void RelionMainWindow::cb_select_to_job_i()
{
	// Show the 'selected' group, hide the others
    int idx = this_to_job_browser->value() - 1;
    if (idx >= 0) // only if a non-empty line was selected
    {
    	current_job = this_to_processes[idx];
		loadJobFromPipeline();
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

	jobCommunicate(DONT_WRITE, DONT_READ, DO_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

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
	fn_settings = outputname;
	std::cerr << " run fn_settings=" << fn_settings << std::endl;
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

	if (commands.size()==0)
	{
		std::cout << " Nothing to do..."<< std::endl;
		return;
	}

	std::cout << "Executing: " << final_command << std::endl;
	int res = system(final_command.c_str());

	// Now save the job (as Running) to the PipeLine
	addToPipeLine(PROC_RUNNING, true); // true means: allow to overwrite an existing process...

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
	fn_settings = ".ScheduledJobs/" + outputname;
	FileName command = "mkdir -p " + fn_settings;
	command = command.beforeLastOf("/");
	int res = system(command.c_str());

	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

	// Now save the job (as Running) to the PipeLine
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
		REPORT_ERROR("RelionMainWindow::cb_run_scheduled_i BUG: this is not a scheduled job!");

	// Get the OLD UNIQDATE of the scheduled job!!
	std::string fn_olduniqdate= findUniqueDateSubstring(fn_settings);

	// Run the scheduled job
	cb_run_i();

	// Re-set the input and output nodes
	int run_job = pipeline.processList.size() - 1;

	// Get the NEW UNIQDATE of the running job
	std::string fn_newuniqdate= findUniqueDateSubstring(fn_settings);

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
			    FileName fn_job = ".ScheduledJobs/" + pipeline.processList[i].name + ".job";
				std::cerr << "fn_job= " << fn_job << std::endl;
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
		ask = "Are you sure you want to delete the following processes? \n";
		for (size_t i = 0; i < deleteProcesses.size(); i++)
		{
			if (deleteProcesses[i])
				ask += " - " + pipeline.processList[i].name + "\n";
		}
		proceed = fl_ask(ask.c_str());
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
					std::cerr << " command= " << command << std::endl;
					int res = system(command.c_str());
				}
				else
				{
					// Move entire output directory (with subdirectory structure) to the Trash folder
					FileName firstdirs = alldirs.beforeLastOf("/");
					std::string command = "mkdir -p Trash/" + firstdirs;
					int res = system(command.c_str());
					std::cerr << " command= " << command << std::endl;
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

	std::string ask;
	ask = "Are you sure you want to delete intermediate files from " + pipeline.processList[current_job].name + "?";
	int proceed = fl_ask(ask.c_str());
	if (proceed)
	{
	    std::cerr << "cleanup todo" << std::endl;
	}

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
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL, DO_MKDIR);

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
