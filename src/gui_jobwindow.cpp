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
#include "src/gui_jobwindow.h"

JobWindow::JobWindow(int _x, int _y, int _w, int _h, const char* title ) : Fl_Box(x,y,w,h,title)
{
	tabs = NULL;
	tab1 = tab2 = tab3 = tab4 = tab5 = tab6 = tab7 = runtab = NULL;
	group1 = group2 = group3 = group4 = group5 = group6 = group7 = queue_group = NULL;
	x = _x; y = _y; w = _w; h = _h;
	current_y = start_y = 0;
	is_continue = false;
}

void JobWindow::clear()
{
	/* This only gives segfaults....
	if (group1 != NULL)
	{
		delete group1;
		group1 = NULL;
	}
	if (group2 != NULL)
	{
		delete group2;
		group2 = NULL;
	}
	if (group3 != NULL)
	{
		delete group3;
		group3 = NULL;
	}
	if (group4 != NULL)
	{
		delete group4;
		group4 = NULL;
	}
	if (group5 != NULL)
	{
		delete group5;
		group5 = NULL;
	}
	if (group6 != NULL)
	{
		delete group6;
		group6 = NULL;
	}
	if (group7 != NULL)
	{
		delete group7;
		group7 = NULL;
	}
	if (queue_group != NULL)
	{
		delete queue_group;
		queue_group = NULL;
	}
	if (tab1 != NULL)
	{
		delete tab1;
		tab1 = NULL;
	}
	if (tab2 != NULL)
	{
		delete tab2;
		tab2 = NULL;
	}
	if (tab3 != NULL)
	{
		delete tab3;
		tab3 = NULL;
	}
	if (tab4 != NULL)
	{
		delete tab4;
		tab4 = NULL;
	}
	if (tab5 != NULL)
	{
		delete tab5;
		tab5 = NULL;
	}
	if (tab6 != NULL)
	{
		delete tab6;
		tab6 = NULL;
	}
	if (tab7 != NULL)
	{
		delete tab7;
		tab7 = NULL;
	}
	if (runtab != NULL)
	{
		delete runtab;
		runtab = NULL;
	}
	if (tabs != NULL)
	{
		delete tabs;
		tabs = NULL;
	}
	*/
}

void JobWindow::setupTabs(int nr_tabs)
{

	current_y = y; // top of the GUI

	char * my_allow_change_dedicated = getenv ("RELION_ALLOW_CHANGE_MINIMUM_DEDICATED");
	if (my_allow_change_dedicated == NULL)
		do_allow_change_minimum_dedicated = DEFAULTMININIMUMDEDICATED;
	else
	{
		int check_allow =  textToInteger(my_allow_change_dedicated);
		do_allow_change_minimum_dedicated = (check_allow == 0) ? false : true;
	}

	// Set up tabs
    if (nr_tabs >= 1) // there is always the running tab, which is not counted on the input nr_tabs!
    {
    	tabs = new Fl_Tabs(x, current_y, w, h - MENUHEIGHT);
    	current_y += TABHEIGHT;
    	tabs->begin();
		tab1 = new Fl_Group(x, current_y , w, h - MENUHEIGHT, "");
		tab1->end();
		tab1->color(GUI_BACKGROUND_COLOR);
		tab1->selection_color(GUI_BACKGROUND_COLOR2);
		if (nr_tabs >= 2)
		{

			tab2 = new Fl_Group(x, current_y , w, h - MENUHEIGHT, "");
			tab2->end();
			tab2->color(GUI_BACKGROUND_COLOR);
			tab2->selection_color(GUI_BACKGROUND_COLOR2);
		}
		if (nr_tabs >= 3)
		{
			tab3 = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
			tab3->end();
			tab3->color(GUI_BACKGROUND_COLOR);
			tab3->selection_color(GUI_BACKGROUND_COLOR2);
		}
		if (nr_tabs >= 4)
		{
			tab4 = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
			tab4->end();
			tab4->color(GUI_BACKGROUND_COLOR);
			tab4->selection_color(GUI_BACKGROUND_COLOR2);
		}
		if (nr_tabs >= 5)
		{
			tab5 = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
			tab5->end();
			tab5->color(GUI_BACKGROUND_COLOR);
			tab5->selection_color(GUI_BACKGROUND_COLOR2);
		}
		if (nr_tabs >= 6)
		{
			tab6 = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
			tab6->end();
			tab6->color(GUI_BACKGROUND_COLOR);
			tab6->selection_color(GUI_BACKGROUND_COLOR2);
		}
		if (nr_tabs >= 7)
		{
			tab7 = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
			tab7->end();
			tab7->color(GUI_BACKGROUND_COLOR);
			tab7->selection_color(GUI_BACKGROUND_COLOR2);
		}
		if (nr_tabs >= 8)
		{
			std::cerr << "ERROR: only 7 job-specific tabs implemented..." << std::endl;
			exit(1);
		}
		current_y += 15;
	    start_y = current_y;

		runtab = new Fl_Group(x, current_y, w, h - MENUHEIGHT, "");
		runtab->label("Running");
		// Fill this in later, when we have the joboptions
		runtab->end();
		setupRunTab();
		runtab->color(GUI_BACKGROUND_COLOR);
		runtab->selection_color(GUI_BACKGROUND_COLOR2);

	    tabs->end();

    }

}

void JobWindow::setupRunTab()
{
	runtab->begin();

	resetHeight();

	bool has_parallel = false;

	if (myjob.joboptions.find("nr_mpi") != myjob.joboptions.end())
    {
		place("nr_mpi", TOGGLE_LEAVE_ACTIVE);
    	has_parallel = true;
    }

	if (myjob.joboptions.find("nr_threads") != myjob.joboptions.end())
	{
		place("nr_threads", TOGGLE_LEAVE_ACTIVE);
    	has_parallel = true;
	}

	// Add a little spacer
	if (has_parallel)
		current_y += STEPY/4;

	// Set up queue groups for running tab
    queue_group = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
    queue_group->end();

    place("do_queue", TOGGLE_LEAVE_ACTIVE, queue_group);

	queue_group->begin();

	place("queuename");

	place("qsub");

	// Two extra fields if needed
	if (myjob.joboptions.find("qsub_extra1") != myjob.joboptions.end())
		place("qsub_extra1");

	if (myjob.joboptions.find("qsub_extra2") != myjob.joboptions.end())
		place("qsub_extra2");

	place("qsubscript");

	place("min_dedicated");
	if (do_allow_change_minimum_dedicated)
		guientries["min_dedicated"].deactivate(false);
	else
		guientries["min_dedicated"].deactivate(true);

	queue_group->end();
	guientries["do_queue"].cb_menu_i(); // This is to make the default effective

    // Add a little spacer
    current_y += STEPY/4;

    place("other_args");

    runtab->end();

}

void JobWindow::place(std::string key, int deactivate_option, Fl_Group * deactivate_this_group)
{
	if (myjob.joboptions.find(key) == myjob.joboptions.end())
		std::cerr << "WARNING: cannot find " << key << " in the defined joboptions of jobtype= " << myjob.type << std::endl;

	guientries[key].place(myjob.joboptions[key], current_y, deactivate_option, deactivate_this_group, do_oldstyle);
}

void JobWindow::place2(std::string key1, std::string key2, std::string label, int deactivate_option)
{
	if (myjob.joboptions.find(key1) == myjob.joboptions.end())
		std::cerr << "WARNING: cannot find " << key1 << " in the defined joboptions of jobtype= " << myjob.type << std::endl;
	if (myjob.joboptions.find(key2) == myjob.joboptions.end())
		std::cerr << "WARNING: cannot find " << key2 << " in the defined joboptions of jobtype= " << myjob.type << std::endl;

	myjob.joboptions[key1].label_gui = label;
	myjob.joboptions[key2].label_gui = "";
	int old_y = current_y;
	guientries[key1].place(myjob.joboptions[key1], current_y, deactivate_option, NULL, do_oldstyle,
			XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	current_y = old_y;
	guientries[key2].place(myjob.joboptions[key2], current_y, deactivate_option, NULL, do_oldstyle,
			XCOL2 + (WCOL2 + COLUMN_SEPARATION) / 2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
}

void JobWindow::place3(std::string key1, std::string key2, std::string key3, std::string label, int deactivate_option)
{
	if (myjob.joboptions.find(key1) == myjob.joboptions.end())
		std::cerr << "WARNING: cannot find " << key1 << " in the defined joboptions of jobtype= " << myjob.type << std::endl;
	if (myjob.joboptions.find(key2) == myjob.joboptions.end())
		std::cerr << "WARNING: cannot find " << key2 << " in the defined joboptions of jobtype= " << myjob.type << std::endl;
	if (myjob.joboptions.find(key3) == myjob.joboptions.end())
		std::cerr << "WARNING: cannot find " << key3 << " in the defined joboptions of jobtype= "  << myjob.type<< std::endl;

	myjob.joboptions[key1].label_gui = label;
	myjob.joboptions[key2].label_gui = "";
	myjob.joboptions[key3].label_gui = "";
	int old_y = current_y;
	guientries[key1].place(myjob.joboptions[key1], current_y, deactivate_option, NULL, do_oldstyle,
			XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	current_y = old_y;
	guientries[key2].place(myjob.joboptions[key2], current_y, deactivate_option, NULL, do_oldstyle,
			XCOL2 + 1 + (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	current_y = old_y;
	guientries[key3].place(myjob.joboptions[key3], current_y, deactivate_option, NULL, do_oldstyle,
			XCOL2 + 1 + 2 * (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);

}

void JobWindow::toggle_new_continue(bool _is_continue)
{

	is_continue = _is_continue;
	myjob.is_continue = _is_continue;

	for (std::map<std::string,GuiEntry>::iterator it=guientries.begin(); it!=guientries.end(); ++it)
	{

		int my_option = (it->second).deactivate_option;
		switch (my_option)
		{
		case TOGGLE_DEACTIVATE:
		{
			(it->second).deactivate(is_continue);
			break;
		}
		case TOGGLE_REACTIVATE:
		{
			(it->second).deactivate(!is_continue);
			break;
		}
		case TOGGLE_ALWAYS_DEACTIVATE:
		{
			(it->second).deactivate(true);
			break;
		}
		case TOGGLE_LEAVE_ACTIVE:
		{
			// do nothing
			break;
		}
		default:
		{
			REPORT_ERROR("ERROR: unrecognised deactivate-option for GUI entry " + it->first);
		}
		}

	}

}

void JobWindow::resetHeight()
{
	current_y = start_y;
}

// Update all values in the Fl_Input entries from the corresponding job_options
void JobWindow::updateMyGui()
{
	for (std::map<std::string,GuiEntry>::iterator it=guientries.begin(); it!=guientries.end(); ++it)
	{
		if (myjob.joboptions.find(it->first) == myjob.joboptions.end())
			std::cerr << "WARNING: cannot find " << it->first << " in the defined joboptions!" <<std::endl;

		(it->second).setValue((myjob.joboptions[it->first]).value);
	}

}

// Update all values in the Fl_Input entries into the corresponding job_options
void JobWindow::updateMyJob()
{
	for (std::map<std::string,JobOption>::iterator it=myjob.joboptions.begin(); it!=myjob.joboptions.end(); ++it)
	{
		if (guientries.find(it->first) == guientries.end())
			std::cerr << "WARNING: cannot find " << it->first << " in the defined joboptions!" <<std::endl;

		it->second.value = std::string(((guientries[it->first]).inp)->value());
	}

}

void JobWindow::initialise(int my_job_type, bool _do_oldstyle)
{

	do_oldstyle = _do_oldstyle;
	switch (my_job_type)
	{
	case PROC_IMPORT:
	{
		myjob.initialise(my_job_type);
		initialiseImportWindow();
		break;
	}
	case PROC_MOTIONCORR:
	{
		myjob.initialise(my_job_type);
		initialiseMotioncorrWindow();
		break;
	}
	case PROC_CTFFIND:
	{
		myjob.initialise(my_job_type);
		initialiseCtffindWindow();
		break;
	}
	case PROC_MANUALPICK:
	{
		myjob.initialise(my_job_type);
		initialiseManualpickWindow();
		break;
	}
	case PROC_AUTOPICK:
	{
		myjob.initialise(my_job_type);
		initialiseAutopickWindow();
		break;
	}
	case PROC_EXTRACT:
	{
		myjob.initialise(my_job_type);
		initialiseExtractWindow();
		break;
	}
	case PROC_SORT:
	{
		myjob.initialise(my_job_type);
		initialiseSortWindow();
		break;
	}
	case PROC_CLASSSELECT:
	{
		myjob.initialise(my_job_type);
		initialiseSelectWindow();
		break;
	}
	case PROC_2DCLASS:
	{
		myjob.initialise(my_job_type);
		initialiseClass2DWindow();
		break;
	}
	case PROC_INIMODEL:
	{
		myjob.initialise(my_job_type);
		initialiseInimodelWindow();
		break;
	}
	case PROC_3DCLASS:
	{
		myjob.initialise(my_job_type);
		initialiseClass3DWindow();
		break;
	}
	case PROC_3DAUTO:
	{
		myjob.initialise(my_job_type);
		initialiseAutorefineWindow();
		break;
	}
	case PROC_MOVIEREFINE:
	{
		myjob.initialise(my_job_type);
		initialiseMovierefineWindow();
		break;
	}
	case PROC_POLISH:
	{
		myjob.initialise(my_job_type);
		initialisePolishWindow();
		break;
	}
	case PROC_MASKCREATE:
	{
		myjob.initialise(my_job_type);
		initialiseMaskcreateWindow();
		break;
	}
	case PROC_JOINSTAR:
	{
		myjob.initialise(my_job_type);
		initialiseJoinstarWindow();
		break;
	}
	case PROC_SUBTRACT:
	{
		myjob.initialise(my_job_type);
		initialiseSubtractWindow();
		break;
	}
	case PROC_POST:
	{
		myjob.initialise(my_job_type);
		initialisePostprocessWindow();
		break;
	}
	case PROC_RESMAP:
	{
		myjob.initialise(my_job_type);
		initialiseLocresWindow();
		break;
	}
	default:
	{
		REPORT_ERROR("ERROR: unrecognised job-type to add to the GUI");
	}
	}

	// read settings if hidden file exists
	myjob.read("", is_continue);

}

void JobWindow::initialiseImportWindow()
{

	setupTabs(1);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_in");

	// Add a little spacer
	current_y += STEPY/2;

	place("node_type");

	tab1->end();

	// Always deactivate the queue option
	guientries["do_queue"].deactivate_option = TOGGLE_ALWAYS_DEACTIVATE;

}

void JobWindow::initialiseMotioncorrWindow()
{

	setupTabs(4);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("input_star_mics", TOGGLE_DEACTIVATE);
	place("do_save_movies", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("first_frame_sum", TOGGLE_DEACTIVATE);
	place("last_frame_sum", TOGGLE_DEACTIVATE);
	place("angpix", TOGGLE_DEACTIVATE);

	tab1->end();

	tab2->begin();
	tab2->label("Motioncor2");
	resetHeight();

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	place("do_motioncor2", TOGGLE_DEACTIVATE, group1);

	group1->begin();

	place("fn_motioncor2_exe", TOGGLE_DEACTIVATE);
	place("fn_gain_ref", TOGGLE_DEACTIVATE);
	place("fn_defect", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place2("patch_x", "patch_y", "Number of patches X, Y", TOGGLE_DEACTIVATE);
	place("group_frames", TOGGLE_DEACTIVATE);
	place("bin_factor", TOGGLE_DEACTIVATE);
	place("bfactor", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;
	place("gpu_ids");

	// Add a little spacer
	current_y += STEPY/2;
	place("other_motioncor2_args", TOGGLE_DEACTIVATE);

	group1->end();
	guientries["do_motioncor2"].cb_menu_i(); // make default active

	tab2->end();
	tab3->begin();
	tab3->label("Unblur");
	resetHeight();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("do_unblur", TOGGLE_DEACTIVATE, group2);

	group2->begin();

	place("fn_unblur_exe", TOGGLE_DEACTIVATE);
	place("fn_summovie_exe", TOGGLE_DEACTIVATE);

	group2->end();
	guientries["do_unblur"].cb_menu_i(); // make default active

	tab3->end();
	tab4->begin();
	tab4->label("Dose-weight");
	resetHeight();

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();

	place("do_dose_weighting", TOGGLE_DEACTIVATE, group3);

	group3->begin();

	place("voltage", TOGGLE_DEACTIVATE);
	place("dose_per_frame", TOGGLE_DEACTIVATE);
	place("pre_exposure", TOGGLE_DEACTIVATE);

	group3->end();
	guientries["do_dose_weighting"].cb_menu_i(); // make default active
	tab4->end();

}

void JobWindow::initialiseCtffindWindow()
{

	setupTabs(4);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("input_star_mics", TOGGLE_DEACTIVATE);
	place("use_noDW", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("cs", TOGGLE_DEACTIVATE);
	place("kv", TOGGLE_DEACTIVATE);
	place("q0", TOGGLE_DEACTIVATE);
	place("angpix", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("dast", TOGGLE_DEACTIVATE);

	tab1->end();

	tab2->begin();
	tab2->label("Searches");
	resetHeight();

	place("box", TOGGLE_DEACTIVATE);
	place("resmin", TOGGLE_DEACTIVATE);
	place("resmax", TOGGLE_DEACTIVATE);
	place("dfmin", TOGGLE_DEACTIVATE);
	place("dfmax", TOGGLE_DEACTIVATE);
	place("dfstep", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	place("do_phaseshift", TOGGLE_DEACTIVATE, group1);
	group1->begin();

	place3("phase_min", "phase_max", "phase_step", "Phase shift - Min, Max, Step (deg)", TOGGLE_DEACTIVATE);

	group1->end();
	guientries["do_phaseshift"].cb_menu_i(); // make default active
	tab2->end();


	tab3->begin();
	tab3->label("CTFFIND-4.1");
	resetHeight();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("use_ctffind4", TOGGLE_DEACTIVATE, group2);
	group2->begin();

	place("fn_ctffind_exe", TOGGLE_DEACTIVATE);

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();

	place("do_movie_thon_rings", TOGGLE_DEACTIVATE, group3);

	group3->begin();

	place("movie_rootname", TOGGLE_DEACTIVATE);
	place("avg_movie_frames", TOGGLE_DEACTIVATE);

	group3->end();
	guientries["do_movie_thon_rings"].cb_menu_i(); // make default active

	// Add a little spacer
	current_y += STEPY/2;

	place("ctf_win", TOGGLE_DEACTIVATE);

	group2->end();
	guientries["use_ctffind4"].cb_menu_i(); // make default active

	tab3->end();

	tab4->begin();
	tab4->label("Gctf");
	resetHeight();

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();

	place("use_gctf", TOGGLE_DEACTIVATE, group4);
	group4->begin();

	place("fn_gctf_exe", TOGGLE_DEACTIVATE);
	place("do_ignore_ctffind_params", TOGGLE_DEACTIVATE);
	place("do_EPA", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;
	place("other_gctf_args", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;
	place("gpu_ids", TOGGLE_LEAVE_ACTIVE);


	group4->end();
	guientries["use_gctf"].cb_menu_i(); // make default active

	tab4->end();


}

void JobWindow::initialiseManualpickWindow()
{
	setupTabs(3);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_in", TOGGLE_DEACTIVATE);
	tab1->end();

	tab2->begin();
	tab2->label("Display");
	resetHeight();

	place("diameter");
	place("micscale");
	place("sigma_contrast");
	place("white_val");
	place("black_val");

	current_y += STEPY/2;
	place("lowpass");
	place("highpass");
	place("angpix");

	current_y += STEPY/2;
	place ("do_startend");

	current_y += STEPY/2;
	place("ctfscale");

	tab2->end();
	tab3->begin();
	tab3->label("Colors");
	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	resetHeight();
	place("do_color", TOGGLE_LEAVE_ACTIVE, group1);

	group1->begin();
	place("color_label");
	place("fn_color");
	place("blue_value");
	place("red_value");
	group1->end();
	guientries["do_color"].cb_menu_i(); // make default active

	tab3->end();

	// Always deactivate the queue option
	guientries["do_queue"].deactivate_option = TOGGLE_ALWAYS_DEACTIVATE;

}
void JobWindow::initialiseAutopickWindow()
{
	setupTabs(4);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_input_autopick", TOGGLE_DEACTIVATE);
	place("fn_refs_autopick", TOGGLE_DEACTIVATE);

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_gauss_ref", TOGGLE_DEACTIVATE, group1);

	group1->begin();

	place("gauss_max", TOGGLE_DEACTIVATE);
	group1->end();

	guientries["do_gauss_ref"].cb_menu_i();

	// Add a little spacer
	current_y += STEPY/2;
	place("angpix", TOGGLE_DEACTIVATE);
	place("particle_diameter", TOGGLE_DEACTIVATE);

	tab1->end();
	tab2->begin();
	tab2->label("References");
	resetHeight();

	//set up group
	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("lowpass", TOGGLE_DEACTIVATE);
	place("highpass", TOGGLE_DEACTIVATE);
	place("angpix_ref", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("psi_sampling_autopick", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("do_invert_refs", TOGGLE_DEACTIVATE);
	place("do_ctf_autopick", TOGGLE_DEACTIVATE, group2);

	group2->begin();

	place("do_ignore_first_ctfpeak_autopick", TOGGLE_DEACTIVATE); //(current_y, "Ignore CTFs until first peak?", false,"Set this to Yes, only if this option was also used to generate the references.");

	group2->end();
	guientries["do_ctf_autopick"].cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("autopicking");
	resetHeight();

	place("threshold_autopick");
	place("mindist_autopick");
	place("maxstddevnoise_autopick");

	current_y += STEPY/2;

	place("do_write_fom_maps");
	place("do_read_fom_maps");

	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
	place("shrink", TOGGLE_DEACTIVATE);

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
    group3->end();
    place("use_gpu", TOGGLE_LEAVE_ACTIVE, group3);

    group3->begin();
	place("gpu_ids");
    group3->end();

    guientries["use_gpu"].cb_menu_i();

	tab3->end();
	tab4->begin();
	tab4->label("Helix");
	resetHeight();

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();

	place("do_pick_helical_segments", TOGGLE_DEACTIVATE, group4);

	group4->begin();

	place("helical_tube_outer_diameter");

	current_y += STEPY/2;

	place("helical_nr_asu");
	place("helical_rise");

	current_y += STEPY/2;

	place("helical_tube_kappa_max");
	place("helical_tube_length_min");

	group4->end();

	guientries["do_pick_helical_segments"].cb_menu_i();

	tab4->end();

}
void JobWindow::initialiseExtractWindow()
{
	setupTabs(3);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

    place("star_mics", TOGGLE_DEACTIVATE);

	current_y += STEPY/2;
	place("coords_suffix", TOGGLE_DEACTIVATE);

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	place("do_reextract", TOGGLE_DEACTIVATE, group1);

	group1->begin();

	place("fndata_reextract", TOGGLE_DEACTIVATE);
	place("do_recenter", TOGGLE_DEACTIVATE);

	group1->end();
	guientries["do_reextract"].cb_menu_i();

	// Add a little spacer
	current_y += STEPY/2;

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_set_angpix", TOGGLE_DEACTIVATE, group2);

	group2->begin();
	place("angpix", TOGGLE_DEACTIVATE); //(current_y, "Pixel size (A)", 1, 0.3, 5, 0.1, "Provide the pixel size in Angstroms in the micrograph (so before any re-scaling).  If you provide input CTF parameters, then leave this value to the default of -1.");
	group2->end();
	guientries["do_set_angpix"].cb_menu_i();

	tab1->end();

	tab2->begin();
	tab2->label("extract");
	resetHeight();

	place("extract_size", TOGGLE_DEACTIVATE); //(current_y,"Particle box size (pix):", 128, 64, 512, 8, "Size of the extracted particles (in pixels). This should be an even number!");
	place("do_invert", TOGGLE_DEACTIVATE); //(current_y, "Invert contrast?", true, "If set to Yes, the contrast in the particles will be inverted.");

	// Add a little spacer
	current_y += STEPY/2;

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();
	place("do_norm", TOGGLE_DEACTIVATE, group3);

	group3->begin();


	place("bg_diameter", TOGGLE_DEACTIVATE); //(current_y, "Diameter background circle (pix): ", -1, -1, 600, 10, "Particles will be normalized to a mean value of zero and a standard-deviation of one for all pixels in the background area.\
The background area is defined as all pixels outside a circle with this given diameter in pixels (before rescaling). When specifying a negative value, a default value of 75% of the Particle box size will be used.");

	place("white_dust", TOGGLE_DEACTIVATE); //(current_y, "Stddev for white dust removal: ", -1, -1, 10, 0.1, "Remove very white pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");

	place("black_dust", TOGGLE_DEACTIVATE); //(current_y, "Stddev for black dust removal: ", -1, -1, 10, 0.1, "Remove very black pixels from the extracted particles. \
Pixels values higher than this many times the image stddev will be replaced with values from a Gaussian distribution. \n \n Use negative value to switch off dust removal.");
	group3->end();
	guientries["do_norm"].cb_menu_i();

	// Add a little spacer
	current_y += STEPY/2;

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();
	place("do_rescale", TOGGLE_DEACTIVATE, group4);
	group4->begin();
	place("rescale", TOGGLE_DEACTIVATE);
	group4->end();
	guientries["do_rescale"].cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("Helix");
	resetHeight();

	group5 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group5->end();

	place("do_extract_helix", TOGGLE_DEACTIVATE, group5);

	group5->begin();

	place("helical_tube_outer_diameter", TOGGLE_DEACTIVATE);

	current_y += STEPY/2;

	place("helical_bimodal_angular_priors", TOGGLE_DEACTIVATE);

	group6 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group6->end();

	current_y += STEPY/2;
	place("do_extract_helical_tubes", TOGGLE_DEACTIVATE, group6);

	group6->begin();

	place("do_cut_into_segments", TOGGLE_DEACTIVATE);
	place("helical_nr_asu", TOGGLE_DEACTIVATE);
	place("helical_rise", TOGGLE_DEACTIVATE);

	group6->end();

	guientries["do_extract_helical_tubes"].cb_menu_i();

	group5->end();

	guientries["do_extract_helix"].cb_menu_i();

	tab3->end();

}
void JobWindow::initialiseSortWindow()
{
	setupTabs(2);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("input_star", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	place("is_autopick", TOGGLE_DEACTIVATE, group1);

	group1->begin();

	place("autopick_refs", TOGGLE_DEACTIVATE);

	group1->end();
	guientries["is_autopick"].cb_menu_i();

	tab1->end();

	tab2->begin();
	tab2->label("References");
	resetHeight();

	place("angpix_ref", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("do_ctf", TOGGLE_DEACTIVATE, group2);

	group2->begin();
	place("do_ignore_first_ctfpeak", TOGGLE_DEACTIVATE);
	group2->end();
	guientries["do_ctf"].cb_menu_i();

	tab2->end();

}
void JobWindow::initialiseSelectWindow()
{
	setupTabs(2);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_model", TOGGLE_DEACTIVATE);
	place("fn_mic", TOGGLE_DEACTIVATE);
	place("fn_data", TOGGLE_DEACTIVATE);
	place("fn_coords", TOGGLE_DEACTIVATE);

	tab1->end();

	tab2->begin();
	tab2->label("Class options");
	resetHeight();

	place("do_recenter", TOGGLE_DEACTIVATE);
	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_regroup", TOGGLE_DEACTIVATE, group1);
	group1->begin();
	place("nr_groups", TOGGLE_DEACTIVATE);
	group1->end();
	guientries["do_regroup"].cb_menu_i();

	tab2->end();

	// Always deactivate the queue option
	guientries["do_queue"].deactivate_option = TOGGLE_ALWAYS_DEACTIVATE;


}
void JobWindow::initialiseClass2DWindow()
{
	setupTabs(6);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_img", TOGGLE_DEACTIVATE);
	place("fn_cont", TOGGLE_REACTIVATE);

	tab1->end();

	tab2->begin();
	tab2->label("CTF");

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	resetHeight();
	place("do_ctf_correction", TOGGLE_DEACTIVATE, group1);

	group1->begin();
	place("ctf_phase_flipped", TOGGLE_DEACTIVATE);
	place("ctf_intact_first_peak", TOGGLE_DEACTIVATE);
	group1->end();

	guientries["do_ctf_correction"].cb_menu_i(); // To make default effective

	tab2->end();

	tab3->begin();
	tab3->label("Optimisation");
	resetHeight();

	//set up groups
	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("nr_classes", TOGGLE_DEACTIVATE);
	place("tau_fudge");

	// Add a little spacer
	current_y += STEPY/2;

	place("nr_iter");
	place("do_subsets", TOGGLE_DEACTIVATE, group2);

	group2->begin();

	place("subset_size", TOGGLE_DEACTIVATE);
	place("max_subsets", TOGGLE_DEACTIVATE);
	group2->end();

	guientries["do_subsets"].cb_menu_i(); // to make default effective

	// Add a little spacer
	current_y += STEPY/2;

	place("particle_diameter");
	place("do_zero_mask", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("highres_limit", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	tab3->end();

	tab4->begin();
	tab4->label("Sampling");

	//set up groups
	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();

	resetHeight();

	place("dont_skip_align", TOGGLE_LEAVE_ACTIVE, group3);

	group3->begin();
	place("psi_sampling");
	place("offset_range");
	place("offset_step");
	group3->end();

	guientries["dont_skip_align"].cb_menu_i(); // to make default effective

	tab4->end();
	tab5->begin();
	tab5->label("Helix");
	resetHeight();

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();

	place("do_helix", TOGGLE_DEACTIVATE, group4);

	group4->begin();

	place("helical_tube_outer_diameter");
	place("do_bimodal_psi");
	place("range_psi");

	group4->end();
	guientries["do_helix"].cb_menu_i(); // to make default effective

	tab5->end();

	tab6->begin();
	tab6->label("Compute");
	resetHeight();

	place("do_parallel_discio");
	place("nr_pool");
	place("do_preread_images");
	place("scratch_dir");
	place("do_combine_thru_disc");

	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
    group5 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
    group5->end();

    place("use_gpu", TOGGLE_LEAVE_ACTIVE, group5);

    group5->begin();
	place("gpu_ids", TOGGLE_LEAVE_ACTIVE);
    group5->end();

    guientries["use_gpu"].cb_menu_i();

	tab6->end();



}
void JobWindow::initialiseInimodelWindow()
{
	setupTabs(5);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_img", TOGGLE_DEACTIVATE);
	place("fn_cont", TOGGLE_REACTIVATE);

	tab1->end();
	tab2->begin();
	tab2->label("CTF");

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	resetHeight();
	place("do_ctf_correction", TOGGLE_DEACTIVATE, group1);

	group1->begin();
	place("ctf_phase_flipped", TOGGLE_DEACTIVATE);
	place("ctf_intact_first_peak", TOGGLE_DEACTIVATE);
	group1->end();

	guientries["do_ctf_correction"].cb_menu_i(); // To make default effective

	tab2->end();

	tab3->begin();
	tab3->label("SGD");
	resetHeight();

	place("particle_diameter");
	place("sym_name", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("nr_iter");
	place("sgd_subset_size");
	place("sgd_write_subsets");
	place("sgd_highres_limit");
	place("sgd_sigma2fudge_halflife");

	tab3->end();
	tab4->begin();
	tab4->label("Sampling");

	resetHeight();

	place("sampling");
	place("offset_range");
	place("offset_step");

	tab4->end();

	tab5->begin();
	tab5->label("Compute");
	resetHeight();

	place("do_parallel_discio");
	place("nr_pool");
	place("do_preread_images");
	place("scratch_dir");
	place("do_combine_thru_disc");

	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
    group5 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
    group5->end();

    place("use_gpu", TOGGLE_LEAVE_ACTIVE, group5);

    group5->begin();
	place("gpu_ids", TOGGLE_LEAVE_ACTIVE);
    group5->end();

    guientries["use_gpu"].cb_menu_i();


	tab5->end();

}
void JobWindow::initialiseClass3DWindow()
{
	setupTabs(7);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_img", TOGGLE_DEACTIVATE);
	place("fn_cont", TOGGLE_REACTIVATE);
	place("fn_ref", TOGGLE_DEACTIVATE);
	place("fn_mask");

	tab1->end();
	tab2->begin();
	tab2->label("Reference");
	resetHeight();

	place("ref_correct_greyscale", TOGGLE_DEACTIVATE);
	place("ini_high", TOGGLE_DEACTIVATE);
	// Add a little spacer
	current_y += STEPY/2;

	place("sym_name", TOGGLE_DEACTIVATE);

	tab2->end();
	tab3->begin();
	tab3->label("CTF");

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	resetHeight();
	place("do_ctf_correction", TOGGLE_DEACTIVATE, group1);
	group1->begin();

	place("ctf_corrected_ref", TOGGLE_DEACTIVATE);
	place("ctf_phase_flipped", TOGGLE_DEACTIVATE);
	place("ctf_intact_first_peak", TOGGLE_DEACTIVATE);

	group1->end();

	guientries["do_ctf_correction"].cb_menu_i(); // To make default effective

	tab3->end();
	tab4->begin();
	tab4->label("Optimisation");
	resetHeight();

	//set up groups
	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("nr_classes", TOGGLE_DEACTIVATE);
	place("tau_fudge");

	// Add a little spacer
	current_y += STEPY/2;

	place("nr_iter");
	place("do_subsets", TOGGLE_DEACTIVATE, group2);

	group2->begin();
	place("subset_size", TOGGLE_DEACTIVATE);
	place("max_subsets", TOGGLE_DEACTIVATE);
	group2->end();
	guientries["do_subsets"].cb_menu_i(); // to make default effective

	// Add a little spacer
	current_y += STEPY/2;


	place("particle_diameter");
	place("do_zero_mask", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("highres_limit", TOGGLE_DEACTIVATE);

	tab4->end();

	tab5->begin();
	tab5->label("Sampling");

	//set up groups
	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();

	resetHeight();

	place("dont_skip_align", TOGGLE_LEAVE_ACTIVE, group3);
	group3->begin();

	place("sampling");
	place("offset_range");
	place("offset_step");

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();

	place("do_local_ang_searches", TOGGLE_LEAVE_ACTIVE, group4);

	group4->begin();
	place("sigma_angles");
	group4->end();
	guientries["do_local_ang_searches"].cb_menu_i(); // to make default effective

	group3->end();
	guientries["dont_skip_align"].cb_menu_i(); // to make default effective

	tab5->end();

	tab6->begin();
	tab6->label("Helix");
	resetHeight();
	group5 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group5->end();

	//helix_text", TOGGLE_DEACTIVATE); //(current_y, "Nov 21, 2015");
	place("do_helix", TOGGLE_DEACTIVATE, group5);

	group5->begin();
	place2("helical_tube_inner_diameter", "helical_tube_outer_diameter", "Tube diameter - inner, outer (A):", TOGGLE_DEACTIVATE);
	place2("range_tilt", "range_psi", "Angular search range - tilt, psi (deg):", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	group8 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group8->end();
	place("do_apply_helical_symmetry", TOGGLE_DEACTIVATE, group8);
	group8->begin();
	place("helical_nr_asu", TOGGLE_DEACTIVATE);
	place2("helical_twist_initial", "helical_rise_initial", "Initial twist (deg), rise (A):", TOGGLE_DEACTIVATE);
	place("helical_z_percentage", TOGGLE_DEACTIVATE);
	group8->end();
	guientries["do_apply_helical_symmetry"].cb_menu_i(); // to make default effective

	// Add a little spacer
	current_y += STEPY/2;

	group6 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group6->end();
	place("do_local_search_helical_symmetry", TOGGLE_DEACTIVATE, group6);
	group6->begin();
	place3("helical_twist_min","helical_twist_max", "helical_twist_inistep", "Twist search - Min, Max, Step (deg):", TOGGLE_DEACTIVATE);
	place3("helical_rise_min", "helical_rise_max", "helical_rise_inistep", "Rise search - Min, Max, Step (A):", TOGGLE_DEACTIVATE);
	group6->end();
	guientries["do_local_search_helical_symmetry"].cb_menu_i(); // to make default effective

	// Add a little spacer
	current_y += STEPY/2;

	place("helical_range_distance", TOGGLE_DEACTIVATE);
	group5->end();
	guientries["do_helix"].cb_menu_i(); // to make default effective
	tab6->end();

	tab7->begin();
	tab7->label("Compute");
	resetHeight();

	place("do_parallel_discio");
	place("nr_pool");
	place("do_preread_images");
	place("scratch_dir");
	place("do_combine_thru_disc");
	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
    group7 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
    group7->end();
	place("use_gpu", TOGGLE_LEAVE_ACTIVE, group7);
	group7->begin();
	place("gpu_ids");
    group7->end();
	guientries["use_gpu"].cb_menu_i(); // This is to make the default effective

	tab7->end();

}
void JobWindow::initialiseAutorefineWindow()
{
	setupTabs(7);
	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_img", TOGGLE_DEACTIVATE);
	place("fn_cont", TOGGLE_REACTIVATE);
	place("fn_ref", TOGGLE_DEACTIVATE);
	place("fn_mask");

	tab1->end();
	tab2->begin();
	tab2->label("Reference");
	resetHeight();

	place("ref_correct_greyscale", TOGGLE_DEACTIVATE);
	place("ini_high", TOGGLE_DEACTIVATE);
	// Add a little spacer
	current_y += STEPY/2;
	place("sym_name", TOGGLE_DEACTIVATE);

	tab2->end();
	tab3->begin();
	tab3->label("CTF");

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	resetHeight();
	place("do_ctf_correction", TOGGLE_DEACTIVATE, group1);

	group1->begin();

	place("ctf_corrected_ref", TOGGLE_DEACTIVATE);
	place("ctf_phase_flipped", TOGGLE_DEACTIVATE);
	place("ctf_intact_first_peak", TOGGLE_DEACTIVATE);

	group1->end();
	guientries["do_ctf_correction"].cb_menu_i(); // To make default effective

	tab3->end();
	tab4->begin();
	tab4->label("Optimisation");
	resetHeight();

	place("particle_diameter");
	place("do_zero_mask", TOGGLE_DEACTIVATE);
	// Add a little spacer
	current_y += STEPY/2;

	place("do_solvent_fsc");

	tab4->end();
	tab5->begin();
	tab5->label("Auto-sampling");
	resetHeight();

	place("sampling", TOGGLE_DEACTIVATE);
	place("offset_range", TOGGLE_DEACTIVATE);
	place("offset_step", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("auto_local_sampling", TOGGLE_DEACTIVATE);

	tab5->end();
	tab6->begin();
	tab6->label("Helix");
	resetHeight();
	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("do_helix", TOGGLE_DEACTIVATE, group2);
	group2->begin();
	place2("helical_tube_inner_diameter", "helical_tube_outer_diameter", "Tube diameter - inner, outer (A):",TOGGLE_DEACTIVATE);
	place2("range_tilt", "range_psi", "Angular search range - tilt, psi (deg):", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	group5 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group5->end();
	place("do_apply_helical_symmetry", TOGGLE_DEACTIVATE, group5);
	group5->begin();
	place("helical_nr_asu", TOGGLE_DEACTIVATE);
	place2("helical_twist_initial", "helical_rise_initial", "Initial twist (deg), rise (A):",TOGGLE_DEACTIVATE);
	place("helical_z_percentage", TOGGLE_DEACTIVATE);
	group5->end();
	guientries["do_apply_helical_symmetry"].cb_menu_i(); // to make default effective

	// Add a little spacer
	current_y += STEPY/2;

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();
	place("do_local_search_helical_symmetry", TOGGLE_DEACTIVATE, group3);
	group3->begin();
	place3("helical_twist_min", "helical_twist_max", "helical_twist_inistep", "Twist search - Min, Max, Step (deg):", TOGGLE_DEACTIVATE);
	place3("helical_rise_min", "helical_rise_max","helical_rise_inistep","Rise search - Min, Max, Step (A):", TOGGLE_DEACTIVATE);
	group3->end();
	guientries["do_local_search_helical_symmetry"].cb_menu_i(); // to make default effective

	// Add a little spacer
	current_y += STEPY/2;

	place("helical_range_distance", TOGGLE_DEACTIVATE);
	group2->end();
	guientries["do_helix"].cb_menu_i(); // to make default effective

	tab6->end();

	tab7->begin();
	tab7->label("Compute");
	resetHeight();

	place("do_parallel_discio");
	place("nr_pool");
	place("do_preread_images");
	place("scratch_dir");
	place("do_combine_thru_disc");

	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
    group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
    group4->end();
	place("use_gpu", TOGGLE_LEAVE_ACTIVE, group4);
	group4->begin();
	place("gpu_ids");
    group4->end();
	guientries["use_gpu"].cb_menu_i(); // This is to make the default effective

	tab7->end();

}
void JobWindow::initialiseMovierefineWindow()
{
	setupTabs(4);


	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_movie_star", TOGGLE_DEACTIVATE);
	place("movie_rootname", TOGGLE_DEACTIVATE);

    // Add a little spacer
	current_y += STEPY/2;

	place("fn_cont", TOGGLE_DEACTIVATE);

    // Add a little spacer
	current_y += STEPY/2;

	place("join_nr_mics", TOGGLE_DEACTIVATE);

	tab1->end();

	tab2->begin();
	tab2->label("extract");
	resetHeight();

	place("first_movie_frame", TOGGLE_DEACTIVATE);
	place("last_movie_frame", TOGGLE_DEACTIVATE);
	place("avg_movie_frames", TOGGLE_DEACTIVATE);
	place("max_mpi_nodes", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("extract_size", TOGGLE_DEACTIVATE);
	place("do_invert", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_rescale", TOGGLE_DEACTIVATE, group1);
	group1->begin();
	place("rescale", TOGGLE_DEACTIVATE);
	group1->end();
	guientries["do_rescale"].cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("normalise");
	resetHeight();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_norm", TOGGLE_DEACTIVATE, group2);

	group2->begin();

	place("bg_diameter", TOGGLE_DEACTIVATE);
	place("white_dust", TOGGLE_DEACTIVATE);
	place("black_dust", TOGGLE_DEACTIVATE);

	group2->end();
	guientries["do_norm"].cb_menu_i();

	tab3->end();
	tab4->begin();
	tab4->label("refine");
	resetHeight();

	place("movie_runavg_window");
	place("movie_sigma_offset");

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();

	place("do_alsorot_movies", TOGGLE_LEAVE_ACTIVE, group3);
	group3->begin();

	place("movie_sigma_angles", TOGGLE_LEAVE_ACTIVE);

	group3->end();
	guientries["do_alsorot_movies"].cb_menu_i(); // to make default effective

	tab4->end();

}
void JobWindow::initialisePolishWindow()
{
	setupTabs(5);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_in", TOGGLE_DEACTIVATE);
	place("fn_mask", TOGGLE_DEACTIVATE);

	tab1->end();
	tab2->begin();
	tab2->label("Movement");
	resetHeight();

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();


	place("do_fit_movement", TOGGLE_LEAVE_ACTIVE, group1);

	group1->begin();
	place("sigma_nb");
	group1->end();
	guientries["do_fit_movement"].cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("Damage");
	resetHeight();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("do_bfactor_weighting", TOGGLE_LEAVE_ACTIVE, group2);
	group2->begin();
	place("perframe_highres");
	place("perframe_bfac_lowres");
	place("average_frame_bfactor");

	group2->end();
	guientries["do_bfactor_weighting"].cb_menu_i();

	current_y += STEPY/2;

	place("sym_name");

	tab3->end();


	tab4->begin();
	tab4->label("Normalise");
	resetHeight();

	place("bg_diameter");
	place("white_dust");
	place("black_dust");

	tab4->end();

	tab5->begin();
	tab5->label("Helix");
	resetHeight();

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();

	place("do_helix", TOGGLE_DEACTIVATE, group4);
	group4->begin();

	place("helical_nr_asu");
	place("helical_twist");
	place("helical_rise");

	group4->end();
	guientries["do_helix"].cb_menu_i(); // to make default effective

	tab5->end();

}
void JobWindow::initialiseMaskcreateWindow()
{
	setupTabs(3);
	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_in", TOGGLE_DEACTIVATE); //(current_y, "Input 3D map:", NODE_3DREF, "", "MRC map files (*.mrc)", "Provide an input MRC map from which to start binarizing the map.");
	tab1->end();

	tab2->begin();
	tab2->label("Mask");
	resetHeight();

	place("lowpass_filter");
	place("angpix");

	// Add a little spacer
    current_y += STEPY/2;

    place("inimask_threshold");
    place("extend_inimask");
    place("width_mask_edge");

	tab2->end();

	tab3->begin();
	tab3->label("Helix");
	resetHeight();

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	place("do_helix", TOGGLE_LEAVE_ACTIVE, group1);

	group1->begin();

	place("helical_z_percentage");
	group1->end();
	guientries["do_helix"].cb_menu_i(); // to make default effective

	tab3->end();

}
void JobWindow::initialiseJoinstarWindow()
{
	setupTabs(3);
	tab1->begin();
	tab1->label("particles");
	resetHeight();

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_part", TOGGLE_DEACTIVATE, group1);
	group1->begin();
	place("fn_part1", TOGGLE_DEACTIVATE);
	place("fn_part2", TOGGLE_DEACTIVATE);
	place("fn_part3", TOGGLE_DEACTIVATE);
	place("fn_part4", TOGGLE_DEACTIVATE);
	group1->end();
	guientries["do_part"].cb_menu_i(); // make default active

	tab1->end();

	tab2->begin();
	tab2->label("micrographs");
	resetHeight();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_mic", TOGGLE_DEACTIVATE, group2);
	group2->begin();
	place("fn_mic1", TOGGLE_DEACTIVATE);
	place("fn_mic2", TOGGLE_DEACTIVATE);
	place("fn_mic3", TOGGLE_DEACTIVATE);
	place("fn_mic4", TOGGLE_DEACTIVATE);
	group2->end();
	guientries["do_mic"].cb_menu_i(); // make default active

	tab2->end();


	tab3->begin();
	tab3->label("movies");
	resetHeight();

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();
	place("do_mov", TOGGLE_DEACTIVATE, group3); //(current_y, "Combine movie STAR files?", false, "", mov_group);
	group3->begin();
	place("fn_mov1", TOGGLE_DEACTIVATE);
	place("fn_mov2", TOGGLE_DEACTIVATE);
	place("fn_mov3", TOGGLE_DEACTIVATE);
	place("fn_mov4", TOGGLE_DEACTIVATE);
	group3->end();
	guientries["do_mov"].cb_menu_i(); // make default active

	tab3->end();

}
void JobWindow::initialiseSubtractWindow()
{
	setupTabs(2);
	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_data", TOGGLE_DEACTIVATE);

    current_y += STEPY/2;
	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	place("do_subtract", TOGGLE_DEACTIVATE, group1);

	group1->begin();
	place("fn_in", TOGGLE_DEACTIVATE);
	place("fn_mask", TOGGLE_DEACTIVATE);
    group1->end();
    guientries["do_subtract"].cb_menu_i(); // make default active

    place("do_fliplabel", TOGGLE_DEACTIVATE);

	tab1->end();

	tab2->begin();
	tab2->label("CTF");
	resetHeight();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("do_ctf_correction", TOGGLE_DEACTIVATE, group2);

	group2->begin();

	place("ctf_phase_flipped", TOGGLE_DEACTIVATE);
	place("ctf_intact_first_peak", TOGGLE_DEACTIVATE);

	group2->end();
	guientries["do_ctf_correction"].cb_menu_i(); // To make default effective

	tab2->end();

}
void JobWindow::initialisePostprocessWindow()
{
	setupTabs(3);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();
	place("fn_in", TOGGLE_DEACTIVATE); //(current_y, "One of the 2 unfiltered half-maps:", NODE_HALFMAP, "", "MRC map files (*half1_class001_unfil.mrc)",  "Provide one of the two unfiltered half-reconstructions that were output upon convergence of a 3D auto-refine run.");
	place("fn_mask", TOGGLE_DEACTIVATE); //(current_y, "Solvent mask:", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a soft mask where the protein is white (1) and the solvent is black (0). Often, the softer the mask the higher resolution estimates you will get. A soft edge of 5-10 pixels is often a good edge width.");

	current_y += STEPY/2;

	place("angpix");

	tab1->end();

	tab2->begin();
	tab2->label("Sharpen");
	resetHeight();

	place("fn_mtf");

	current_y += STEPY/2;

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_auto_bfac", TOGGLE_LEAVE_ACTIVE, group1);

	group1->begin();
	place("autob_lowres");
	group1->end();
	guientries["do_auto_bfac"].cb_menu_i();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_adhoc_bfac", TOGGLE_LEAVE_ACTIVE, group2);

	group2->begin();
	place("adhoc_bfac");
	group2->end();
	guientries["do_adhoc_bfac"].cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("Filter");
	resetHeight();
	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();

	place("do_skip_fsc_weighting", TOGGLE_LEAVE_ACTIVE, group3);

	group3->begin();
	place("low_pass");
	group3->end();
	guientries["do_skip_fsc_weighting"].cb_menu_i();

	tab3->end();

}
void JobWindow::initialiseLocresWindow()
{
	setupTabs(3);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_in", TOGGLE_DEACTIVATE);

	current_y += STEPY/2;

	place("angpix", TOGGLE_DEACTIVATE);

	tab1->end();

	tab2->begin();
	tab2->label("ResMap");
	resetHeight();

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	place("do_resmap_locres", TOGGLE_DEACTIVATE, group1);

	group1->begin();

	place("fn_resmap", TOGGLE_DEACTIVATE);

	current_y += STEPY /2 ;

	place("fn_mask");

	current_y += STEPY /2 ;

	place("pval", TOGGLE_DEACTIVATE);
	place("minres", TOGGLE_DEACTIVATE);
	place("maxres", TOGGLE_DEACTIVATE);
	place("stepres", TOGGLE_DEACTIVATE);

	group1->end();
	guientries["do_resmap_locres"].cb_menu_i();

	tab2->end();

	tab3->begin();
	tab3->label("Relion");
	resetHeight();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("do_relion_locres", TOGGLE_DEACTIVATE, group2);

	group2->begin();

	//place("locres_sampling", TOGGLE_DEACTIVATE);
	//place("randomize_at", TOGGLE_DEACTIVATE);
	//current_y += STEPY /2 ;
	place("adhoc_bfac", TOGGLE_DEACTIVATE);
	place("fn_mtf", TOGGLE_DEACTIVATE);

	group2->end();
	guientries["do_relion_locres"].cb_menu_i();

	tab3->end();

}


