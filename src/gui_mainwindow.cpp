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

RelionMainWindow::RelionMainWindow(int w, int h, const char* title):Fl_Window(w,h,title)
{

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
    browser->add("General");
    browser->add("Micrograph inspection");
    browser->add("CTF estimation");
    browser->add("Auto-picking");
    browser->add("Particle extraction");
    browser->add("Particle sorting");
    browser->add("2D classification");
    browser->add("3D classification");
    browser->add("3D auto-refine");
    browser->add("Post-processing");
    browser->add("Particle polishing");
    browser->add("Local-resolution");
    browser->add("Publish!");

    // browse page
    {
        browse_grp[0] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_general = new GeneralJobWindow();
    	browse_grp[0]->end();
    }
    // browse page
    {
        browse_grp[1] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_manualpick = new ManualpickJobWindow();
       	browse_grp[1]->end();
    }
    // browse page
    {
        browse_grp[2] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_ctffind = new CtffindJobWindow();
    	browse_grp[2]->end();
    }
    // browse page
    {

        browse_grp[3] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_autopick = new AutopickJobWindow();
       	browse_grp[3]->end();
    }
    // browse page
    {

        browse_grp[4] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_extract = new ExtractJobWindow();
       	browse_grp[4]->end();
    }
    // browse page
    {

        browse_grp[5] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_sort = new SortJobWindow();
       	browse_grp[5]->end();
    }
    // browse page
    {

        browse_grp[6] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_class2d = new Class2DJobWindow();
       	browse_grp[6]->end();
    }
    // browse page
    {

        browse_grp[7] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_class3d = new Class3DJobWindow();
       	browse_grp[7]->end();
    }
    // browse page
    {

        browse_grp[8] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_auto3d = new Auto3DJobWindow();
       	browse_grp[8]->end();
    }
    // browse page
    {
        browse_grp[9] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_post = new PostJobWindow();
       	browse_grp[9]->end();
    }
    // browse page
    {

        browse_grp[10] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_polish = new PolishJobWindow();
       	browse_grp[10]->end();
    }
    // browse page
    {
        browse_grp[11] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_resmap = new ResmapJobWindow();
       	browse_grp[11]->end();
    }
    // browse page
    {
        browse_grp[12] = new Fl_Group(WCOL0, MENUHEIGHT, 550, 600-MENUHEIGHT);
    	job_publish = new PublishJobWindow();
       	browse_grp[12]->end();
    }
    browser->callback(cb_select_browsegroup);
    browser->end();

    // Set and activate current selection
    browser->select(1); // just start from the beginning
    cb_select_browsegroup_i(); // make default active
    toggle_continue->value(0); // 0 = new run; 1 = continue
    cb_toggle_continue_i(); // make default active

    resizable();
}

void RelionMainWindow::jobCommunicate(bool do_write, bool do_read, bool do_toggle_continue, bool do_commandline, int this_job)
{
	int myval = (this_job > 0) ? this_job : browser->value();

	// always write the general settings with the (hidden) empty name
	if (do_write)
		job_general->write("");

	if (myval == 1)
	{
		if (do_write)
			job_general->write(fn_settings);
		if (do_read)
			job_general->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_general->toggle_new_continue(is_main_continue);
		if (do_commandline)
			job_general->getCommands(outputname, commands, final_command);
	}
	else if (myval == 2)
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
	}
	else if (myval == 3)
	{
		if (do_write)
			job_ctffind->write(fn_settings);
		if (do_read)
			job_ctffind->read(fn_settings, is_main_continue);
		if (do_toggle_continue)
			job_ctffind->toggle_new_continue(is_main_continue);
		if (do_commandline)
			job_ctffind->getCommands(outputname, commands, final_command, job_general->angpix.getValue());
	}
	else if (myval == 4)
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
	}
	else if (myval == 5)
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
	}
	else if (myval == 6)
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
	}
	else if (myval == 7)
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
	}
	else if (myval == 8)
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
	}
	else if (myval == 9)
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
	}
	else if (myval == 10)
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
	}
	else if (myval == 11)
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
	}
	else if (myval == 12)
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
	}
	else if (myval == 13)
	{
		if (do_toggle_continue)
			job_publish->toggle_new_continue(is_main_continue);
		if (do_commandline)
			job_publish->getCommands(outputname, commands, final_command);
	}

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

// Display button call-back functions
void RelionMainWindow::cb_display(Fl_Widget* o, void* v) {

    RelionMainWindow* T=(RelionMainWindow*)v;
    T->cb_display_i();
}


void RelionMainWindow::cb_display_i()
{
	std::string command = " relion_display --gui &" ;
	system(command.c_str());
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
	// Save temporary hidden file with this jobs settings
	fn_settings = "";
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL);

	// Get the command line arguments from the currently active jobwindow,
	// save job submission script, and prepare the final command
	jobCommunicate(DONT_WRITE, DONT_READ, DONT_TOGGLE_CONT, DO_GET_CL);

	// Also save a copy of the GUI settings with the current output name
	fn_settings = outputname;
	jobCommunicate(DO_WRITE, DONT_READ, DONT_TOGGLE_CONT, DONT_GET_CL);

	if (commands.size()==0)
	{
		std::cout << " Nothing to do..."<< std::endl;
		return;
	}

	std::cout << "Executing: " << final_command << std::endl;
	system(final_command.c_str());

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
    	else if  (fn_settings.rfind(".gui_post.") < npos)
    		browser->value(10);
    	else if  (fn_settings.rfind(".gui_polish.") < npos)
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
If RELION is useful in your work, please cite us in the contexts as explained under the \"Publish!\" tab, or on the RELION wiki at http://www2.mrc-lmb.cam.ac.uk/relion. \n  \
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
