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

#ifndef SRC_GUI_JOBWINDOW_H_
#define SRC_GUI_JOBWINDOW_H_
#include "src/pipeline_jobs.h"
#include "src/gui_entries.h"


class JobWindow : public Fl_Box
{
protected:

	// Which process type is this?
	bool is_continue;

public:

	// Which job is this for?
	RelionJob myjob;

	// Is this for tomo?
	bool is_tomo;

	// All the GuiEntries of this job
	std::map<std::string, GuiEntry> guientries;

	// oldstyle GUI
	bool do_oldstyle;

	// Sizes
	int x, y, w, h;

	// Position of current cursor to place new elements
	int start_y, current_y;

	// Tabs
	Fl_Tabs *tabs;
	Fl_Group *tab1, *tab2, *tab3, *tab4, *tab5, *tab6, *tab7, *runtab;

	// Groups
	Fl_Group *group1, *group2, *group3, *group4, *group5, *group6, *group7, *group8, *queue_group;

public:
	// Constructor with x, y, w, h and a title
	JobWindow(int _x = WCOL0, int _y = 2, int _w = GUIWIDTH - WCOL0 - 10, int _h = GUIHEIGHT_OLD-65, const char* title = "");

	// Destructor
	~JobWindow() { clear(); }

	// Clear everything
	void clear();

	// set up the tabs
	void setupTabs(int nr_tabs);

	void initialise(int my_job_type, bool _is_tomo = false);

	void resetHeight();

	// Place a single entry
	void place(std::string key, int deactivate_option = TOGGLE_LEAVE_ACTIVE, Fl_Group * deactivate_this_group = NULL, bool actually_activate=false);
	// Place two entries on one line
	void place2(std::string key1, std::string key2, std::string label, int deactivate_option = TOGGLE_LEAVE_ACTIVE);
	// Place three entries on one line
	void place3(std::string key1, std::string key2, std::string key3, std::string label, int deactivate_option = TOGGLE_LEAVE_ACTIVE);

	void setupRunTab();

	// De/re-activate options upon toggling the continue button
	void toggle_new_continue(bool is_continue);

	// Write the job submission script
	void saveJobSubmissionScript(std::string newfilename, std::string outputname, std::vector<std::string> commands);

	// Initialise pipeiline stuff for each job, return outputname
	void initialisePipeline(std::string &outputname, std::string defaultname, int job_counter);

	// Prepare the final (job submission or combined (mpi) command of possibly multiple lines)
	// Returns true to go ahead, and false to cancel
	bool prepareFinalCommand(std::string &outputname, std::vector<std::string> &commands, std::string &final_command, bool do_makedir = true);

	// Update all values in the Fl_Input entries into/from the corresponding job_options
	void updateMyGui();
	void updateMyJob();

private:

	static void cb_menu_continue(Fl_Widget*, void*);
	inline void cb_menu_continue_i();


	// initialise the window for each of the jobtypes
	void initialiseImportWindow();
	void initialiseMotioncorrWindow();
	void initialiseCtffindWindow();
	void initialiseManualpickWindow();
	void initialiseAutopickWindow();
	void initialiseExtractWindow();
	void initialiseSelectWindow();
	void initialiseClass2DWindow();
	void initialiseInimodelWindow();
	void initialiseClass3DWindow();
	void initialiseAutorefineWindow();
	void initialiseMultiBodyWindow();
	void initialiseMaskcreateWindow();
	void initialiseJoinstarWindow();
	void initialiseSubtractWindow();
	void initialisePostprocessWindow();
	void initialiseLocresWindow();
	void initialiseMotionrefineWindow();
	void initialiseCtfrefineWindow();
	void initialiseExternalWindow();

	// relion-3.2: add subtomogram averaging programs by Jasenko
	void placeTomoInput(bool has_tomograms, bool has_particles,
						bool has_trajectories, bool has_manifolds, bool has_halfmaps, bool has_postprocess);
	void initialiseTomoImportWindow();
	void initialiseTomoSubtomoWindow();
	void initialiseTomoCtfRefineWindow();
	void initialiseTomoAlignWindow();
	void initialiseTomoReconParWindow();

};



#endif /* SRC_NEWGUI_JOBWINDOW_H_ */
