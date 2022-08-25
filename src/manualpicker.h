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

#ifndef MANUALPICKER_H_
#define MANUALPICKER_H_
// this define, and the undef below the FL includes, protects against another Complex definition in fltk
#define Complex tmpComplex
#include <FL/Fl.H>
#include <FL/Fl_Shared_Image.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Image.H>
#include <FL/Fl_JPEG_Image.H>
#include <FL/Fl_Box.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Float_Input.H>
#include <FL/Fl_Text_Display.H>
#undef Complex
#include "src/metadata_table.h"
#include "src/args.h"
#include "src/funcs.h"
#include "src/filename.h"
#include "src/gui_entries.h"
#include <src/jaz/single_particle/obs_model.h>

#define MWCOL1 300
#define MWCOL2 60
#define MWCOL3 60
#define MWCOL4 80
#define MXCOL0 30
#define MXCOL1 (MXCOL0 + MWCOL1 + 10)
#define MXCOL2 (MXCOL1 + MWCOL2 + 10)
#define MXCOL3 (MXCOL2 + MWCOL3 + 10)
#define MXCOL4 (MXCOL3 + MWCOL4 + 10)
#define TOTALWIDTH (MWCOL1 + MWCOL2 + MWCOL3 + MWCOL4 + MWCOL4 + 100)
#define TOTALHEIGHT 500

// The button for picking particles
void cb_viewmic(Fl_Widget* w, void* data);
// The button for viewing the CTF
void cb_viewctf(Fl_Widget* w, void* data);
// The selection button
void cb_selectmic(Fl_Widget* w, void* data);

// This class only puts scrollbars around the resizable canvas
class manualpickerGuiWindow : public Fl_Window
{
public:

	// Input, picking & output names
	FileName fn_in, fn_sel;

	// Allow saving selected micrographs?
	bool do_allow_save;

	// Save default selection immediately? (useful for always generating output files in pipeline)
	bool do_fast_save;

	// MetaDataTable of input micrographs
	MetaDataTable MDin;

	// MetaDataTable with the output coordinate files
	MetaDataTable MDcoords;

	// Observation model of input micrographs
	ObservationModel obsModel;

	// Constructor with w x h size of the window and a title
	manualpickerGuiWindow(int W, int H, const char* title=0): Fl_Window(W, H, title){}

	// Fill the window with all entries
	int fill();

private:

    static void cb_menubar_save(Fl_Widget*, void*);
    inline void cb_menubar_save_i();

    static void cb_menubar_select_all(Fl_Widget*, void*);
    inline void cb_menubar_select_all_i();

    static void cb_menubar_invert_selection(Fl_Widget*, void*);
    inline void cb_menubar_invert_selection_i();

    static void cb_menubar_quit(Fl_Widget*, void*);
    inline void cb_menubar_quit_i();

    static void cb_menubar_recount(Fl_Widget*, void*);
    inline void cb_menubar_recount_i();

    static void cb_menubar_setFOM(Fl_Widget*, void*);
    inline void cb_menubar_setFOM_i();

    void readOutputStarfile();
    void writeOutputStarfiles(bool verb = true);

};

class ManualPicker
{
public:

	// I/O Parser
	IOParser parser;

	// The input micrographs
	MetaDataTable MDin;

	// MetaDataTable with the output coordinate files
	MetaDataTable MDcoords;

	// Observation model for the input mirographs
	ObservationModel obsModel;

	// Input, picking & output names
	FileName fn_in, fn_sel;

	// Allow save selected micrographs?
	bool do_allow_save;

	// Save an output selection file immediately (with all micrographs selected)
	bool do_fast_save;


public:
	// Read command line arguments
	void read(int argc, char **argv);

	// Print usage instructions
	void usage();

	// Initialise some general stuff after reading
	void initialise();

	// General function to decide what to do
	void run();

};

#endif /* MANUALPICKER_H_ */
