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
JobWindow::JobWindow(int _x, int _y, int _w, int _h, const char* title ) : Fl_Box(_x,_y,_w,_h,title)
{
	clear();
	x = _x; y = _y; w = _w; h = _h;
}

void JobWindow::clear()
{
	tabs = NULL;
	tab1 = tab2 = tab3 = tab4 = tab5 = tab6 = tab7 = runtab = NULL;
	group1 = group2 = group3 = group4 = group5 = group6 = group7 = queue_group = NULL;
	current_y = start_y = 0;
	is_continue = false;
	is_tomo = false;
	guientries.clear();
}

void JobWindow::setupTabs(int nr_tabs)
{
	current_y = y; // top of the GUI

	char * my_allow_change_dedicated = getenv ("RELION_ALLOW_CHANGE_MINIMUM_DEDICATED");
	if (my_allow_change_dedicated == NULL)
	{
		do_allow_change_minimum_dedicated = DEFAULTMININIMUMDEDICATED;
	}
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

	char * extra_count_text = getenv ("RELION_QSUB_EXTRA_COUNT");
	const char extra_count_val = (extra_count_text ? atoi(extra_count_text) : 2);
	for (int i=1; i<=extra_count_val; i++)
	{
		std::stringstream out;
		out<<i;
		const std::string i_str=out.str();
		if (myjob.joboptions.find(std::string("qsub_extra")+i_str) != myjob.joboptions.end())
		{
			place(std::string("qsub_extra")+i_str);
		}
	}

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

void JobWindow::place(std::string key, int deactivate_option, Fl_Group * deactivate_this_group, bool actually_activate)
{
	if (myjob.joboptions.find(key) == myjob.joboptions.end())
		std::cerr << "WARNING: cannot find " << key << " in the defined joboptions of jobtype= " << myjob.type << std::endl;

	guientries[key].place(myjob.joboptions[key], current_y, deactivate_option, deactivate_this_group, actually_activate);
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
	guientries[key1].place(myjob.joboptions[key1], current_y, deactivate_option, NULL, false,
	                       XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION) / 2);
	current_y = old_y;
	guientries[key2].place(myjob.joboptions[key2], current_y, deactivate_option, NULL, false,
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
	guientries[key1].place(myjob.joboptions[key1], current_y, deactivate_option, NULL, false,
	                       XCOL2, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	current_y = old_y;
	guientries[key2].place(myjob.joboptions[key2], current_y, deactivate_option, NULL, false,
	                       XCOL2 + 1 + (WCOL2 + COLUMN_SEPARATION) / 3, STEPY, (WCOL2 - COLUMN_SEPARATION * 2) / 3);
	current_y = old_y;
	guientries[key3].place(myjob.joboptions[key3], current_y, deactivate_option, NULL, false,
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
		{
			std::cerr << "ERROR: cannot find " << it->first << " in the defined joboptions!" <<std::endl;
			REPORT_ERROR("Stopping now...");
		}

		it->second.value = std::string(((guientries[it->first]).inp)->value());
	}
}

void JobWindow::initialise(int my_job_type, bool _is_tomo)
{
	is_tomo = _is_tomo;
	myjob.setTomo(_is_tomo);

	if (my_job_type == PROC_IMPORT)
	{
		myjob.initialise(my_job_type);
		initialiseImportWindow();
	}
	else if (my_job_type == PROC_MOTIONCORR)
	{
		myjob.initialise(my_job_type);
		initialiseMotioncorrWindow();
	}
	else if (my_job_type == PROC_CTFFIND)
	{
		myjob.initialise(my_job_type);
		initialiseCtffindWindow();
	}
	else if (my_job_type == PROC_MANUALPICK)
	{
		myjob.initialise(my_job_type);
		initialiseManualpickWindow();
	}
	else if (my_job_type == PROC_AUTOPICK)
	{
		myjob.initialise(my_job_type);
		initialiseAutopickWindow();
	}
	else if (my_job_type == PROC_EXTRACT)
	{
		myjob.initialise(my_job_type);
		initialiseExtractWindow();
	}
	else if (my_job_type == PROC_CLASSSELECT)
	{
		myjob.initialise(my_job_type);
		initialiseSelectWindow();
	}
	else if (my_job_type == PROC_2DCLASS)
	{
		myjob.initialise(my_job_type);
		initialiseClass2DWindow();
	}
	else if (my_job_type == PROC_INIMODEL)
	{
		myjob.initialise(my_job_type);
		initialiseInimodelWindow();
	}
	else if (my_job_type == PROC_3DCLASS)
	{
		myjob.initialise(my_job_type);
		initialiseClass3DWindow();
	}
	else if (my_job_type == PROC_3DAUTO)
	{
		myjob.initialise(my_job_type);
		initialiseAutorefineWindow();
	}
	else if (my_job_type == PROC_MULTIBODY)
	{
		myjob.initialise(my_job_type);
		initialiseMultiBodyWindow();
	}
	else if (my_job_type == PROC_MASKCREATE)
	{
		myjob.initialise(my_job_type);
		initialiseMaskcreateWindow();
	}
	else if (my_job_type == PROC_JOINSTAR)
	{
		myjob.initialise(my_job_type);
		initialiseJoinstarWindow();
	}
	else if (my_job_type == PROC_SUBTRACT)
	{
		myjob.initialise(my_job_type);
		initialiseSubtractWindow();
	}
	else if (my_job_type == PROC_POST)
	{
		myjob.initialise(my_job_type);
		initialisePostprocessWindow();
	}
	else if (my_job_type == PROC_RESMAP)
	{
		myjob.initialise(my_job_type);
		initialiseLocresWindow();
	}
	else if (my_job_type == PROC_MOTIONREFINE)
	{
		myjob.initialise(my_job_type);
		initialiseMotionrefineWindow();
	}
	else if (my_job_type == PROC_CTFREFINE)
	{
		myjob.initialise(my_job_type);
		initialiseCtfrefineWindow();
	}
	else if (my_job_type == PROC_TOMO_IMPORT)
	{
		myjob.initialise(my_job_type);
		initialiseTomoImportWindow();
	}
	else if (my_job_type == PROC_TOMO_SUBTOMO)
	{
		myjob.initialise(my_job_type);
		initialiseTomoSubtomoWindow();
	}
	else if (my_job_type == PROC_TOMO_CTFREFINE)
	{
		myjob.initialise(my_job_type);
		initialiseTomoCtfRefineWindow();
	}
	else if (my_job_type == PROC_TOMO_ALIGN)
	{
		myjob.initialise(my_job_type);
		initialiseTomoAlignWindow();
	}
	else if (my_job_type == PROC_TOMO_RECONSTRUCT)
	{
		myjob.initialise(my_job_type);
		initialiseTomoReconParWindow();
	}
	else if (my_job_type == PROC_EXTERNAL)
	{
		myjob.initialise(my_job_type);
		initialiseExternalWindow();
	}
	else
	{
		REPORT_ERROR("ERROR: unrecognised job-type to add to the GUI");
	}

	// read settings if hidden file exists
	myjob.read("", is_continue);

	// update the window
	updateMyGui();
}

void JobWindow::initialiseImportWindow()
{
	setupTabs(2);

	tab1->begin();
	tab1->label("Movies/mics");
	resetHeight();

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	place("do_raw", TOGGLE_DEACTIVATE, group1, false);

	// Add a little spacer
	current_y += STEPY/2;

	group1->begin();
	place("fn_in_raw");
	place("is_multiframe");

	// Add a little spacer
	current_y += STEPY/2;

	place("optics_group_name");
	place("fn_mtf");
	place("angpix");
	place("kV");
	place("Cs");
	place("Q0");
	place("beamtilt_x");
	place("beamtilt_y");
	group1->end();

	guientries["do_raw"].cb_menu_i(); // make default active
	tab1->end();

	tab2->begin();
	tab2->label("Others");
	resetHeight();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_other", TOGGLE_DEACTIVATE, group2, false);
	group2->begin();
	// Add a little spacer
	current_y += STEPY/2;
	place("fn_in_other");
	place("node_type");
	// Add a little spacer
	current_y += STEPY/2;
	place("optics_group_particles");
	group2->end();

	guientries["do_other"].cb_menu_i(); // make default active

	tab2->end();

	// Always deactivate the queue option
	guientries["do_queue"].deactivate_option = TOGGLE_ALWAYS_DEACTIVATE;
	myjob.joboptions["do_queue"].setString("No");
}

void JobWindow::initialiseMotioncorrWindow()
{
	setupTabs(2);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("input_star_mics", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("first_frame_sum", TOGGLE_DEACTIVATE);
	place("last_frame_sum", TOGGLE_DEACTIVATE);
	place("dose_per_frame", TOGGLE_DEACTIVATE);
	place("pre_exposure", TOGGLE_DEACTIVATE);
	place("eer_grouping", TOGGLE_DEACTIVATE);
	place("do_float16", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_dose_weighting", TOGGLE_DEACTIVATE, group1);
	group1->begin();
	place("do_save_noDW", TOGGLE_DEACTIVATE);
	group1->end();

	guientries["do_dose_weighting"].cb_menu_i(); // make default active

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_save_ps", TOGGLE_DEACTIVATE, group2);
	group2->begin();
	place("group_for_ps", TOGGLE_DEACTIVATE);
	group2->end();
	tab1->end();

	tab2->begin();
	tab2->label("Motion");
	resetHeight();

	place("bfactor", TOGGLE_DEACTIVATE);
	place2("patch_x", "patch_y", "Number of patches X, Y", TOGGLE_DEACTIVATE);
	place("group_frames", TOGGLE_DEACTIVATE);
	place("bin_factor", TOGGLE_DEACTIVATE);
	place("fn_gain_ref", TOGGLE_DEACTIVATE);
	place("gain_rot", TOGGLE_DEACTIVATE);
	place("gain_flip", TOGGLE_DEACTIVATE);
	place("fn_defect", TOGGLE_DEACTIVATE);

	current_y += STEPY/2;
	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();
	place("do_own_motioncor", TOGGLE_DEACTIVATE, group4, true);
	group4->begin();
	place("fn_motioncor2_exe", TOGGLE_DEACTIVATE);
	place("gpu_ids");
	place("other_motioncor2_args", TOGGLE_DEACTIVATE);
	group4->end();

	guientries["do_own_motioncor"].cb_menu_i(); // make default active

	tab2->end();
}

void JobWindow::initialiseCtffindWindow()
{
	setupTabs(3);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("input_star_mics", TOGGLE_DEACTIVATE);
	place("use_noDW", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	place("do_phaseshift", TOGGLE_DEACTIVATE, group1);
	group1->begin();

	place3("phase_min", "phase_max", "phase_step", "Phase shift - Min, Max, Step (deg)", TOGGLE_DEACTIVATE);

	group1->end();
	guientries["do_phaseshift"].cb_menu_i(); // make default active

	// Add a little spacer
	current_y += STEPY/2;

	place("dast", TOGGLE_DEACTIVATE);

	tab1->end();

	tab2->begin();
	tab2->label("CTFFIND-4.1");
	resetHeight();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("use_ctffind4", TOGGLE_DEACTIVATE, group2);
	group2->begin();

	place("fn_ctffind_exe", TOGGLE_DEACTIVATE);
	place("use_given_ps", TOGGLE_DEACTIVATE);
	place("slow_search", TOGGLE_DEACTIVATE);

	place("ctf_win", TOGGLE_DEACTIVATE);

	group2->end();
	guientries["use_ctffind4"].cb_menu_i(); // make default active

	// Add a little spacer
	current_y += STEPY/2;

	place("box", TOGGLE_DEACTIVATE);
	place("resmin", TOGGLE_DEACTIVATE);
	place("resmax", TOGGLE_DEACTIVATE);
	place("dfmin", TOGGLE_DEACTIVATE);
	place("dfmax", TOGGLE_DEACTIVATE);
	place("dfstep", TOGGLE_DEACTIVATE);

	tab2->end();

	tab3->begin();
	tab3->label("Gctf");
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

	tab3->end();
}

void JobWindow::initialiseManualpickWindow()
{
	setupTabs(3);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_in", TOGGLE_DEACTIVATE);

	current_y += STEPY/2;
	place ("do_startend");

	current_y += STEPY/2;

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_fom_threshold", TOGGLE_DEACTIVATE, group1);
	group1->begin();
	place("minimum_pick_fom", TOGGLE_DEACTIVATE);
	group1->end();
	guientries["do_fom_threshold"].cb_menu_i();

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

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_topaz_denoise", TOGGLE_DEACTIVATE, group2);
	group2->begin();
	place("fn_topaz_exec", TOGGLE_DEACTIVATE);
	group2->end();
	guientries["do_topaz_denoise"].cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("Colors");

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();

	resetHeight();
	place("do_color", TOGGLE_LEAVE_ACTIVE, group3);

	group3->begin();
	place("color_label");
	place("fn_color");
	place("blue_value");
	place("red_value");
	group3->end();
	guientries["do_color"].cb_menu_i(); // make default active

	tab3->end();

	// Always deactivate the queue option
	guientries["do_queue"].deactivate_option = TOGGLE_ALWAYS_DEACTIVATE;
	myjob.joboptions["do_queue"].setString("No");
}

void JobWindow::initialiseAutopickWindow()
{
	setupTabs(6);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_input_autopick", TOGGLE_DEACTIVATE);
	place("angpix", TOGGLE_DEACTIVATE);

	current_y += STEPY/2;

	place("do_refs", TOGGLE_DEACTIVATE);
	place("do_log", TOGGLE_DEACTIVATE);
	place("do_topaz", TOGGLE_DEACTIVATE);
	place("continue_manual", TOGGLE_REACTIVATE);

	tab1->end();
	tab2->begin();
	tab2->label("Laplacian");
	resetHeight();

	place("log_diam_min", TOGGLE_DEACTIVATE);
	place("log_diam_max", TOGGLE_DEACTIVATE);
	place("log_invert", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;
	place("log_maxres", TOGGLE_DEACTIVATE);
	place("log_adjust_thr");
	place("log_upper_thr");

	tab2->end();
	tab3->begin();
	tab3->label("Topaz");
	resetHeight();

	place("fn_topaz_exec");
	place("topaz_particle_diameter", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	group7 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group7->end();

	place("do_topaz_pick", TOGGLE_DEACTIVATE, group7);
	group7->end();

	group7->begin();
	place("topaz_model", TOGGLE_DEACTIVATE);
	group7->end();
	guientries["do_topaz_pick"].cb_menu_i();

	// Add a little spacer
	current_y += STEPY/2;

	group5 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group5->end();
	place("do_topaz_train", TOGGLE_DEACTIVATE, group5);
	group5->begin();

	place("topaz_nr_particles", TOGGLE_DEACTIVATE);
	place("topaz_train_picks", TOGGLE_DEACTIVATE);

	group6 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group6->end();

	place("do_topaz_train_parts", TOGGLE_DEACTIVATE, group6);

	group6->begin();
	place("topaz_train_parts", TOGGLE_DEACTIVATE);
	group6->end();
	guientries["do_topaz_train_parts"].cb_menu_i();

	group5->end();
	guientries["do_topaz_train"].cb_menu_i();

	// Add a little spacer
	current_y += STEPY/2;

	place("topaz_other_args", TOGGLE_DEACTIVATE);

	tab3->end();
	tab4->begin();
	tab4->label("References");
	resetHeight();

	place("fn_refs_autopick", TOGGLE_DEACTIVATE);

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_ref3d", TOGGLE_DEACTIVATE, group1);
	group1->begin();
	place("fn_ref3d_autopick", TOGGLE_DEACTIVATE);
	place("ref3d_symmetry", TOGGLE_DEACTIVATE);
	place("ref3d_sampling", TOGGLE_DEACTIVATE);
	group1->end();
	guientries["do_ref3d"].cb_menu_i();

	//set up group
	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	// Add a little spacer
	current_y += STEPY/2;

	place("lowpass", TOGGLE_DEACTIVATE);
	place("highpass", TOGGLE_DEACTIVATE);
	place("angpix_ref", TOGGLE_DEACTIVATE);
	place("psi_sampling_autopick", TOGGLE_DEACTIVATE);
	place("do_invert_refs", TOGGLE_DEACTIVATE);
	place("do_ctf_autopick", TOGGLE_DEACTIVATE, group2);

	group2->begin();

	place("do_ignore_first_ctfpeak_autopick", TOGGLE_DEACTIVATE); //(current_y, "Ignore CTFs until first peak?", false,"Set this to Yes, only if this option was also used to generate the references.");

	group2->end();
	guientries["do_ctf_autopick"].cb_menu_i();

	tab4->end();
	tab5->begin();
	tab5->label("autopicking");
	resetHeight();

	place("threshold_autopick");
	place("mindist_autopick");
	place("maxstddevnoise_autopick");
	place("minavgnoise_autopick");

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

	tab5->end();
	tab6->begin();
	tab6->label("Helix");
	resetHeight();

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();

	place("do_pick_helical_segments", TOGGLE_DEACTIVATE, group4);

	group4->begin();

	place("helical_tube_outer_diameter");
	place("helical_tube_length_min");
	place("helical_tube_kappa_max");

	current_y += STEPY/2;

	place("helical_nr_asu");
	place("helical_rise");

	current_y += STEPY/2;

	place("do_amyloid");

	group4->end();

	guientries["do_pick_helical_segments"].cb_menu_i();

	tab6->end();
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
	current_y += STEPY/2;

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();

	place("do_reextract", TOGGLE_DEACTIVATE, group1);

	group1->begin();

	place("fndata_reextract", TOGGLE_DEACTIVATE);
	place("do_reset_offsets", TOGGLE_DEACTIVATE);
	group7 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group7->end();
	place("do_recenter", TOGGLE_DEACTIVATE, group7);

	group7->begin();
	place3("recenter_x","recenter_y", "recenter_z", "Recenter on - X, Y, Z (pix):", TOGGLE_DEACTIVATE);
	group7->end();
	guientries["do_recenter"].cb_menu_i();

	group1->end();
	guientries["do_reextract"].cb_menu_i();

	// Add a little spacer
	current_y += STEPY/2;

	place("do_float16", TOGGLE_DEACTIVATE);

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

	// Add a little spacer
	current_y += STEPY/2;

	group7 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group7->end();
	place("do_fom_threshold", TOGGLE_DEACTIVATE, group7);
	group7->begin();
	place("minimum_pick_fom", TOGGLE_DEACTIVATE);
	group7->end();
	guientries["do_fom_threshold"].cb_menu_i();

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

void JobWindow::initialiseSelectWindow()
{
	setupTabs(4);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_model", TOGGLE_DEACTIVATE);
	place("fn_mic", TOGGLE_DEACTIVATE);
	place("fn_data", TOGGLE_DEACTIVATE);

	tab1->end();

	tab2->begin();
	tab2->label("Class options");
	resetHeight();

	group6 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group6->end();


	place("do_class_ranker", TOGGLE_DEACTIVATE, group6);
	group6->begin();
	place("rank_threshold", TOGGLE_DEACTIVATE);
	place("select_nr_parts", TOGGLE_DEACTIVATE);
	place("select_nr_classes", TOGGLE_DEACTIVATE);
	place("python_exe", TOGGLE_DEACTIVATE);

	group6->end();
	guientries["do_class_ranker"].cb_menu_i();

	current_y += STEPY/2;

	place("do_recenter", TOGGLE_DEACTIVATE);
	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_regroup", TOGGLE_DEACTIVATE, group1);
	group1->begin();
	place("nr_groups", TOGGLE_DEACTIVATE);
	group1->end();
	guientries["do_regroup"].cb_menu_i();
	tab2->end();

	tab3->begin();
	tab3->label("Subsets");
	resetHeight();

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();

	place("do_select_values", TOGGLE_DEACTIVATE, group3);
	group3->begin();
	place("select_label", TOGGLE_DEACTIVATE);
	place("select_minval", TOGGLE_DEACTIVATE);
	place("select_maxval", TOGGLE_DEACTIVATE);
	group3->end();
	guientries["do_select_values"].cb_menu_i();

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();

	// Add a little spacer
	current_y += STEPY/2;

	place("do_discard", TOGGLE_DEACTIVATE, group4);
	group4->begin();
	place("discard_label", TOGGLE_DEACTIVATE);
	place("discard_sigma", TOGGLE_DEACTIVATE);
	group4->end();
	guientries["do_discard"].cb_menu_i();

	group5 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group5->end();

	// Add a little spacer
	current_y += STEPY/2;

	place("do_split", TOGGLE_DEACTIVATE, group5);
	group5->begin();
	place("do_random", TOGGLE_DEACTIVATE);
	place("split_size", TOGGLE_DEACTIVATE);
	place("nr_split", TOGGLE_DEACTIVATE);
	group5->end();
	guientries["do_split"].cb_menu_i();

	tab3->end();

	tab4->begin();
	tab4->label("Duplicates");
	resetHeight();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_remove_duplicates", TOGGLE_DEACTIVATE, group2);
	group2->begin();
	place("duplicate_threshold", TOGGLE_DEACTIVATE);
	place("image_angpix", TOGGLE_DEACTIVATE);
	group2->end();
	guientries["do_remove_duplicates"].cb_menu_i();
	tab4->end();

	// Always deactivate the queue option
	guientries["do_queue"].deactivate_option = TOGGLE_ALWAYS_DEACTIVATE;
	myjob.joboptions["do_queue"].setString("No");
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
	place("ctf_intact_first_peak", TOGGLE_DEACTIVATE);
	group1->end();

	guientries["do_ctf_correction"].cb_menu_i(); // To make default effective

	tab2->end();

	tab3->begin();
	tab3->label("Optimisation");
	resetHeight();

	place("nr_classes", TOGGLE_DEACTIVATE);
	place("tau_fudge");

	// Add a little spacer
	current_y += STEPY/2;

	//set up groups
	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("do_em", TOGGLE_DEACTIVATE, group2);

	group2->begin();

	place("nr_iter_em");

	group2->end();

	guientries["do_em"].cb_menu_i(); // to make default effective


	//set up groups
	group5 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group5->end();

	place("do_grad", TOGGLE_DEACTIVATE, group5);

	group5->begin();

	place("nr_iter_grad");

	group5->end();

	guientries["do_grad"].cb_menu_i(); // to make default effective


	// Add a little spacer
	current_y += STEPY/2;

	place("particle_diameter");
	place("do_zero_mask", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("highres_limit", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("do_center");

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

	current_y += STEPY/2;
	place("allow_coarser");

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

	group7 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group7->end();
	place("do_restrict_xoff", TOGGLE_LEAVE_ACTIVE, group7);

	group7->begin();
	place("helical_rise", TOGGLE_LEAVE_ACTIVE);
	group7->end();

	guientries["do_restrict_xoff"].cb_menu_i();

	group4->end();
	guientries["do_helix"].cb_menu_i(); // to make default effective

	tab5->end();

	tab6->begin();
	tab6->label("Compute");
	resetHeight();

	place("do_parallel_discio");
	place("nr_pool");
	group5 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group5->end();
	place("do_preread_images", TOGGLE_LEAVE_ACTIVE, group5, true);
	group5->begin();
	place("scratch_dir");
	group5->end();
	place("do_combine_thru_disc");

	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
	group6 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group6->end();

	place("use_gpu", TOGGLE_LEAVE_ACTIVE, group6);

	group6->begin();
	place("gpu_ids", TOGGLE_LEAVE_ACTIVE);
	group6->end();

	guientries["use_gpu"].cb_menu_i();

	tab6->end();
}

void JobWindow::initialiseInimodelWindow()
{
	setupTabs(4);

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
	place("ctf_intact_first_peak", TOGGLE_DEACTIVATE);
	group1->end();

	guientries["do_ctf_correction"].cb_menu_i(); // To make default effective

	tab2->end();

	tab3->begin();
	tab3->label("Optimisation");
	resetHeight();

	place("nr_iter");
	place("tau_fudge");
	place("nr_classes", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("particle_diameter");
	place("do_solvent", TOGGLE_DEACTIVATE);
	place("sym_name", TOGGLE_DEACTIVATE);
	place("do_run_C1", TOGGLE_DEACTIVATE);

	tab3->end();

	tab4->begin();
	tab4->label("Compute");
	resetHeight();

	place("do_parallel_discio");
	place("nr_pool");
	group5 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group5->end();
	place("do_preread_images", TOGGLE_LEAVE_ACTIVE, group5, true);
	group5->begin();
	place("scratch_dir");
	group5->end();
	place("do_combine_thru_disc");

	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
	group6 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group6->end();

	place("use_gpu", TOGGLE_LEAVE_ACTIVE, group6);

	group6->begin();
	place("gpu_ids", TOGGLE_LEAVE_ACTIVE);
	group5->end();

	guientries["use_gpu"].cb_menu_i();

	tab4->end();
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
	place("do_fast_subsets", TOGGLE_DEACTIVATE);

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
	place("relax_sym");
	group4->end();
	guientries["do_local_ang_searches"].cb_menu_i(); // to make default effective

	current_y += STEPY/2;
	place("allow_coarser");

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
	place3("range_rot", "range_tilt", "range_psi", "Angular search range - rot, tilt, psi (deg):", TOGGLE_DEACTIVATE);
	place("helical_range_distance", TOGGLE_DEACTIVATE);
	place("keep_tilt_prior_fixed", TOGGLE_DEACTIVATE);

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

	group5->end();
	guientries["do_helix"].cb_menu_i(); // to make default effective
	tab6->end();

	tab7->begin();
	tab7->label("Compute");
	resetHeight();

	place("do_parallel_discio");
	place("nr_pool");
	place("do_pad1");
	group7 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group7->end();
	place("do_preread_images", TOGGLE_LEAVE_ACTIVE, group7, true);
	group7->begin();
	place("scratch_dir");
	group7->end();
	place("do_combine_thru_disc");
	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
	group8 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group8->end();
	place("use_gpu", TOGGLE_LEAVE_ACTIVE, group8);
	group8->begin();
	place("gpu_ids");
	group8->end();
	guientries["use_gpu"].cb_menu_i(); // This is to make the default effective

	tab7->end();
}

void JobWindow::initialiseAutorefineWindow()
{
	setupTabs(7);
	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	if (is_tomo)
	{
		place("in_optimisation", TOGGLE_DEACTIVATE);
		current_y += STEPY /2 ;
	}

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

	current_y += STEPY/2;
	place("auto_local_sampling", TOGGLE_DEACTIVATE);
	place("relax_sym");
	current_y += STEPY/2;
	place("auto_faster");

	tab5->end();
	tab6->begin();
	tab6->label("Helix");
	resetHeight();
	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();

	place("do_helix", TOGGLE_DEACTIVATE, group2);
	group2->begin();
	place2("helical_tube_inner_diameter", "helical_tube_outer_diameter", "Tube diameter - inner, outer (A):",TOGGLE_DEACTIVATE);
	place3("range_rot", "range_tilt", "range_psi", "Angular search range - rot, tilt, psi (deg):", TOGGLE_DEACTIVATE);
	place("helical_range_distance", TOGGLE_DEACTIVATE);
	place("keep_tilt_prior_fixed", TOGGLE_DEACTIVATE);

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

	group2->end();
	guientries["do_helix"].cb_menu_i(); // to make default effective

	tab6->end();

	tab7->begin();
	tab7->label("Compute");
	resetHeight();

	place("do_parallel_discio");
	place("nr_pool");
	place("do_pad1");
	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();
	place("do_preread_images", TOGGLE_LEAVE_ACTIVE, group4, true);
	group4->begin();
	place("scratch_dir");
	group4->end();
	place("do_combine_thru_disc");

	// Add a little spacer
	current_y += STEPY/2;

	// Set up queue groups for running tab
	group5 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group5->end();
	place("use_gpu", TOGGLE_LEAVE_ACTIVE, group5);
	group5->begin();
	place("gpu_ids");
	group5->end();
	guientries["use_gpu"].cb_menu_i(); // This is to make the default effective

	tab7->end();
}

void JobWindow::initialiseMultiBodyWindow()
{
	setupTabs(4);
	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_in", TOGGLE_DEACTIVATE);
	place("fn_cont", TOGGLE_REACTIVATE);
	place("fn_bodies", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("do_subtracted_bodies", TOGGLE_DEACTIVATE);

	tab1->end();
	tab2->begin();
	tab2->label("Auto-sampling");
	resetHeight();

	place("sampling", TOGGLE_DEACTIVATE);
	place("offset_range", TOGGLE_DEACTIVATE);
	place("offset_step", TOGGLE_DEACTIVATE);

	tab2->end();

	tab3->begin();
	tab3->label("Analyse");
	resetHeight();

	group5 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group5->end();

	place("do_analyse", TOGGLE_LEAVE_ACTIVE, group5);
	group5->begin();

	place("nr_movies");

	group6 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group6->end();

	place("do_select", TOGGLE_LEAVE_ACTIVE, group6);

	group6->begin();
	place("select_eigenval");
	place("eigenval_min");
	place("eigenval_max");
	group6->end();
	guientries["do_select"].cb_menu_i(); // This is to make the default effective

	group5->end();
	guientries["do_analyse"].cb_menu_i(); // This is to make the default effective

	tab3->end();

	tab4->begin();
	tab4->label("Compute");
	resetHeight();

	place("do_parallel_discio");
	place("nr_pool");
	place("do_pad1");
	group7 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group7->end();
	place("do_preread_images", TOGGLE_LEAVE_ACTIVE, group7, true);
	group7->begin();
	place("scratch_dir");
	group7->end();
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

	tab4->end();
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

	place("fn_opt", TOGGLE_DEACTIVATE);
	place("fn_mask", TOGGLE_DEACTIVATE);

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_data", TOGGLE_DEACTIVATE, group1);

	group1->begin();
	place("fn_data", TOGGLE_DEACTIVATE);
	group1->end();
	guientries["do_data"].cb_menu_i(); // make default active
	place("do_float16", TOGGLE_DEACTIVATE);

	current_y += STEPY/2;

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_fliplabel", TOGGLE_DEACTIVATE, group2);

	group2->begin();
	place("fn_fliplabel", TOGGLE_DEACTIVATE);
	group2->end();
	guientries["do_fliplabel"].cb_menu_i(); // make default active

	tab1->end();

	tab2->begin();
	tab2->label("Centering");
	resetHeight();

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();
	place("do_center_mask", TOGGLE_DEACTIVATE, group3, true);

	group3->begin();

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();
	place("do_center_xyz", TOGGLE_DEACTIVATE, group4);

	group4->begin();
	place3("center_x", "center_y", "center_z", "Center coordinate - X, Y, Z (pix):", TOGGLE_DEACTIVATE);
	group4->end();
	guientries["do_center_xyz"].cb_menu_i(); // To make default effective

	group3->end();
	guientries["do_center_mask"].cb_menu_i(); // To make default effective

	current_y += STEPY/2;

	place("new_box", TOGGLE_DEACTIVATE);

	tab2->end();
}

void JobWindow::initialisePostprocessWindow()
{
	setupTabs(2);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	if (is_tomo)
	{
		place("in_optimisation", TOGGLE_DEACTIVATE);
		current_y += STEPY /2 ;
	}

	place("fn_in", TOGGLE_DEACTIVATE); //(current_y, "One of the 2 unfiltered half-maps:", NODE_HALFMAP, "", "MRC map files (*half1_class001_unfil.mrc)",  "Provide one of the two unfiltered half-reconstructions that were output upon convergence of a 3D auto-refine run.");
	place("fn_mask", TOGGLE_DEACTIVATE); //(current_y, "Solvent mask:", NODE_MASK, "", "Image Files (*.{spi,vol,msk,mrc})", "Provide a soft mask where the protein is white (1) and the solvent is black (0). Often, the softer the mask the higher resolution estimates you will get. A soft edge of 5-10 pixels is often a good edge width.");

	current_y += STEPY/2;

	place("angpix");

	tab1->end();

	tab2->begin();
	tab2->label("Sharpen");
	resetHeight();

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

	current_y += STEPY/2;

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();

	place("do_skip_fsc_weighting", TOGGLE_LEAVE_ACTIVE, group3);

	group3->begin();
	place("low_pass");
	group3->end();
	guientries["do_skip_fsc_weighting"].cb_menu_i();

	current_y += STEPY/2;

	place("fn_mtf");
	place("mtf_angpix");

	tab2->end();
}

void JobWindow::initialiseLocresWindow()
{
	setupTabs(3);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("fn_in", TOGGLE_DEACTIVATE);
	place("fn_mask");

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

void JobWindow::initialiseMotionrefineWindow()
{
	setupTabs(3);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	// I/O
	place("fn_mic", TOGGLE_DEACTIVATE);
	place("fn_data", TOGGLE_DEACTIVATE);
	place("fn_post", TOGGLE_DEACTIVATE);

	current_y += STEPY /2 ;

	place("first_frame", TOGGLE_DEACTIVATE);
	place("last_frame", TOGGLE_DEACTIVATE);

	current_y += STEPY /2 ;

	place("extract_size", TOGGLE_DEACTIVATE);
	place("rescale", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("do_float16", TOGGLE_DEACTIVATE);

	tab1->end();

	tab2->begin();
	tab2->label("Train");
	resetHeight();

	// Train for optimal parameters
	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_param_optim", TOGGLE_LEAVE_ACTIVE, group2);

	group2->begin();

	place("eval_frac");
	place("optim_min_part");

	group2->end();
	guientries["do_param_optim"].cb_menu_i();

	tab2->end();

	tab3->begin();
	tab3->label("Polish");
	resetHeight();

	// Polishing
	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_polish", TOGGLE_DEACTIVATE, group1);

	current_y += STEPY /2 ;

	group1->begin();

	place("opt_params", TOGGLE_DEACTIVATE);

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();
	place("do_own_params", TOGGLE_DEACTIVATE, group4);

	group4->begin();
	place("sigma_vel", TOGGLE_DEACTIVATE);
	place("sigma_div", TOGGLE_DEACTIVATE);
	place("sigma_acc", TOGGLE_DEACTIVATE);
	group4->end();
	guientries["do_own_params"].cb_menu_i();

	current_y += STEPY /2 ;

	place("minres", TOGGLE_DEACTIVATE);
	place("maxres", TOGGLE_DEACTIVATE);

	tab3->end();
}

void JobWindow::initialiseCtfrefineWindow()
{
	setupTabs(2);

	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	// I/O
	place("fn_data", TOGGLE_DEACTIVATE);
	place("fn_post", TOGGLE_DEACTIVATE);

	tab1->end();

	tab2->begin();
	tab2->label("Fit");
	resetHeight();

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();
	place("do_aniso_mag", TOGGLE_LEAVE_ACTIVE, group3, true); //true means: activating aniso_mag will deactive higher-order aberrations

	current_y += STEPY /2 ;

	group3->begin();

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_ctf", TOGGLE_LEAVE_ACTIVE, group1);

	group1->begin();

	place("do_defocus", TOGGLE_LEAVE_ACTIVE);
	place("do_astig", TOGGLE_LEAVE_ACTIVE);
	place("do_bfactor", TOGGLE_LEAVE_ACTIVE);
	place("do_phase", TOGGLE_LEAVE_ACTIVE);

	group1->end();
	guientries["do_ctf"].cb_menu_i();

	current_y += STEPY /2 ;

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();
	place("do_tilt", TOGGLE_LEAVE_ACTIVE, group4);

	group4->begin();

	place("do_trefoil", TOGGLE_LEAVE_ACTIVE);

	group4->end();
	guientries["do_tilt"].cb_menu_i();

	current_y += STEPY /2 ;

	place("do_4thorder", TOGGLE_LEAVE_ACTIVE);

	group3->end();

	guientries["do_aniso_mag"].cb_menu_i();
	current_y += STEPY /2 ;

	place("minres", TOGGLE_DEACTIVATE);

	tab2->end();
}

void JobWindow::initialiseExternalWindow()
{
	setupTabs(2);

	tab1->begin();
	tab1->label("Input");
	resetHeight();

	// I/O
	place("fn_exe", TOGGLE_DEACTIVATE);

	current_y += STEPY /2 ;
	place("in_mov", TOGGLE_DEACTIVATE);
	place("in_mic", TOGGLE_DEACTIVATE);
	place("in_part", TOGGLE_DEACTIVATE);
	place("in_coords", TOGGLE_DEACTIVATE);
	place("in_3dref", TOGGLE_DEACTIVATE);
	place("in_mask", TOGGLE_DEACTIVATE);

	tab1->end();

	tab2->begin();
	tab2->label("Params");
	resetHeight();

	place2("param1_label", "param1_value", "Param1 label, value:", TOGGLE_LEAVE_ACTIVE);
	place2("param2_label", "param2_value", "Param2 label, value:", TOGGLE_LEAVE_ACTIVE);
	place2("param3_label", "param3_value", "Param3 label, value:", TOGGLE_LEAVE_ACTIVE);
	place2("param4_label", "param4_value", "Param4 label, value:", TOGGLE_LEAVE_ACTIVE);
	place2("param5_label", "param5_value", "Param5 label, value:", TOGGLE_LEAVE_ACTIVE);
	place2("param6_label", "param6_value", "Param6 label, value:", TOGGLE_LEAVE_ACTIVE);
	place2("param7_label", "param7_value", "Param7 label, value:", TOGGLE_LEAVE_ACTIVE);
	place2("param8_label", "param8_value", "Param8 label, value:", TOGGLE_LEAVE_ACTIVE);
	place2("param9_label", "param9_value", "Param9 label, value:", TOGGLE_LEAVE_ACTIVE);
	place2("param10_label", "param10_value", "Param10 label, value:", TOGGLE_LEAVE_ACTIVE);

	tab2->end();
}

void JobWindow::placeTomoInput(bool has_tomograms, bool has_particles,
							   bool has_trajectories, bool has_manifolds, bool has_halfmaps, bool has_postprocess)
{
	tab1->begin();
	tab1->label("I/O");
	resetHeight();

	place("in_optimisation", TOGGLE_DEACTIVATE);

	current_y += STEPY /2 ;

	if (has_particles) place("in_particles", TOGGLE_DEACTIVATE);
	if (has_tomograms) place("in_tomograms", TOGGLE_DEACTIVATE);
	if (has_trajectories) place("in_trajectories", TOGGLE_DEACTIVATE);
	if (has_manifolds) place("in_manifolds", TOGGLE_DEACTIVATE);
	if (has_halfmaps) place("in_halfmaps", TOGGLE_DEACTIVATE);
	if (has_postprocess)
	{
		place("in_refmask", TOGGLE_DEACTIVATE);
		place("in_post", TOGGLE_DEACTIVATE);
	}

	tab1->end();

}

void JobWindow::initialiseTomoImportWindow()
{
	setupTabs(3);

	tab1->begin();
	tab1->label("Tomograms");
	resetHeight();

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_tomo", TOGGLE_DEACTIVATE, group1, false);
	group1->begin();

	place("tomo_star", TOGGLE_DEACTIVATE);
	place("io_tomos", TOGGLE_DEACTIVATE);
	place("angpix", TOGGLE_DEACTIVATE);
	place("kV", TOGGLE_DEACTIVATE);
	place("Cs", TOGGLE_DEACTIVATE);
	place("Q0", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("dose", TOGGLE_DEACTIVATE);
	place("order_list", TOGGLE_DEACTIVATE);
	place("do_flipYZ", TOGGLE_DEACTIVATE);
	place("do_flipZ", TOGGLE_DEACTIVATE);
	place("hand", TOGGLE_DEACTIVATE);

	group1->end();
	guientries["do_tomo"].cb_menu_i(); // make default active

	tab1->end();

	tab2->begin();
	tab2->label("Coordinates");
	resetHeight();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_coords", TOGGLE_DEACTIVATE, group2, false);
	group2->begin();

	place("part_star", TOGGLE_DEACTIVATE);
	place("part_tomos", TOGGLE_DEACTIVATE);
	// Add a little spacer
	current_y += STEPY/2;
	place("do_coords_flipZ", TOGGLE_DEACTIVATE);

	group2->end();
	guientries["do_coords"].cb_menu_i(); // make default active

	tab2->end();

	tab3->begin();
	tab3->label("Others");
	resetHeight();

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();

	place("do_other", TOGGLE_DEACTIVATE, group3, false);
	group3->begin();

	// Add a little spacer
	current_y += STEPY/2;

	place("fn_in_other", TOGGLE_DEACTIVATE);
	place("node_type", TOGGLE_DEACTIVATE);

	// Add a little spacer
	current_y += STEPY/2;

	place("optics_group_particles", TOGGLE_DEACTIVATE);

	group3->end();
	guientries["do_other"].cb_menu_i(); // make default active

	tab3->end();
}

void JobWindow::initialiseTomoSubtomoWindow()
{
	setupTabs(2);

	placeTomoInput(true, true, true, false, false, false);

	tab2->begin();
	tab2->label("Reconstruct");
	resetHeight();

	place("box_size", TOGGLE_DEACTIVATE);
	place("crop_size", TOGGLE_DEACTIVATE);
	place("binning", TOGGLE_DEACTIVATE);

    current_y += STEPY /2 ;

    place("do_float16", TOGGLE_DEACTIVATE);

	current_y += STEPY /2 ;

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_cone_weight", TOGGLE_DEACTIVATE, group1);

	group1->begin();

	place("cone_angle", TOGGLE_DEACTIVATE);

	group1->end();
	guientries["do_cone_weight"].cb_menu_i();

	tab2->end();
}

void JobWindow::initialiseTomoCtfRefineWindow()
{
	setupTabs(3);

	placeTomoInput(true, true, true, false, true, true);

	tab2->begin();
	tab2->label("Defocus");
	resetHeight();

	place("box_size", TOGGLE_DEACTIVATE);

	current_y += STEPY /2 ;

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_defocus", TOGGLE_DEACTIVATE, group1);

	group1->begin();

	place("focus_range", TOGGLE_DEACTIVATE);

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_reg_def", TOGGLE_DEACTIVATE, group2);

	group2->begin();
	place("lambda", TOGGLE_DEACTIVATE);
	group2->end();
	guientries["do_reg_def"].cb_menu_i();

	group1->end();
	guientries["do_defocus"].cb_menu_i();


	current_y += STEPY /2 ;

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();
	place("do_scale", TOGGLE_DEACTIVATE, group3);

	group3->begin();

	place("do_frame_scale", TOGGLE_DEACTIVATE);
	place("do_tomo_scale", TOGGLE_DEACTIVATE);

	group3->end();
	guientries["do_scale"].cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("Aberrations");
	resetHeight();

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();
	place("do_odd_aberr", TOGGLE_DEACTIVATE, group3);

	group3->begin();

	place("nr_odd_aberr", TOGGLE_DEACTIVATE);

	group3->end();
	guientries["do_odd_aberr"].cb_menu_i();

	current_y += STEPY /2 ;

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();
	place("do_even_aberr", TOGGLE_DEACTIVATE, group4);

	group4->begin();

	place("nr_even_aberr", TOGGLE_DEACTIVATE);

	group4->end();
	guientries["do_even_aberr"].cb_menu_i();

	tab3->end();
}

void JobWindow::initialiseTomoAlignWindow()
{
	setupTabs(4);

	placeTomoInput(true, true, true, false, true, true);

	tab2->begin();
	tab2->label("Polish");
	resetHeight();

	place("box_size", TOGGLE_DEACTIVATE);
	place("max_error", TOGGLE_DEACTIVATE);

	current_y += STEPY /2 ;

	group3 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group3->end();
	place("do_shift_align", TOGGLE_DEACTIVATE, group3);
	group3->begin();

	place("shift_align_type", TOGGLE_DEACTIVATE);
	group3->end();
	guientries["do_shift_align"].cb_menu_i();

	tab2->end();
	tab3->begin();
	tab3->label("Motion");
	resetHeight();

	group2 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group2->end();
	place("do_motion", TOGGLE_DEACTIVATE, group2);

	current_y += STEPY /2 ;

	group2->begin();

	place("sigma_vel", TOGGLE_DEACTIVATE);
	place("sigma_div", TOGGLE_DEACTIVATE);
	place("do_sq_exp_ker", TOGGLE_DEACTIVATE);

	group2->end();
	guientries["do_motion"].cb_menu_i();

	tab3->end();
	tab4->begin();
	tab4->label("Deformations");
	resetHeight();

	group4 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group4->end();
	place("do_deform", TOGGLE_DEACTIVATE, group4);

	current_y += STEPY /2 ;

	group4->begin();

	place("def_w", TOGGLE_DEACTIVATE);
	place("def_h", TOGGLE_DEACTIVATE);
	place("def_model", TOGGLE_DEACTIVATE);
	place("lambda", TOGGLE_DEACTIVATE);

	current_y += STEPY /2 ;

	place("do_frame_def", TOGGLE_DEACTIVATE);

	group4->end();
	guientries["do_deform"].cb_menu_i();
	tab4->end();
}

void JobWindow::initialiseTomoReconParWindow()
{
	setupTabs(2);

	placeTomoInput(true, true, true, false, false, false);

	tab2->begin();
	tab2->label("Average");
	resetHeight();

	group1 = new Fl_Group(WCOL0,  MENUHEIGHT, 550, 600-MENUHEIGHT, "");
	group1->end();
	place("do_from2d", TOGGLE_DEACTIVATE, group1);

	group1->begin();
	place("box_size", TOGGLE_DEACTIVATE);
	place("crop_size", TOGGLE_DEACTIVATE);
	place("binning", TOGGLE_DEACTIVATE);
	place("snr", TOGGLE_DEACTIVATE);
	place("fn_mask", TOGGLE_DEACTIVATE);

	group1->end();
	guientries["do_from2d"].cb_menu_i();

	current_y += STEPY /2 ;

	place("sym_name", TOGGLE_DEACTIVATE);

	tab2->end();
}

