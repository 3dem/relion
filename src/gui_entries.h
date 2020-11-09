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

#ifndef SRC_GUI_ENTRIES_H_
#define SRC_GUI_ENTRIES_H_

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Text_Buffer.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Float_Input.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Group.H>
#include <FL/fl_ask.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Hold_Browser.H>
#include <FL/Fl_Select_Browser.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Menu_Button.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Toggle_Button.H>
#include <FL/Fl_Widget.H>
#include <FL/Fl_Wizard.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Text_Editor.H>
#include <FL/Fl_Image.H>
#include <FL/Fl_XPM_Image.H>
#include "src/macros.h"
#include "src/strings.h"
#include "src/filename.h"
#include "src/pipeline_jobs.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <cstdio>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define GUI_MAX_RADIO_ELEMENTS 15

// forward declaration current_browse_directory, which allows CURRENT_ODIR browse buttons
extern std::string current_browse_directory;

// Create the scheduler GUI: without sliders or pull-down menus
extern bool create_scheduler_gui;


// Gui layout
//#define XCOL1 10
//#define XCOL2 260
//#define XCOL3 460
//#define XCOL4 475
//#define XCOL5 535
//Additional space in tab if more than 4 XXXextraiXXX template variables are used defined by
//environment variable RELION_QSUB_EXTRA_COUNT
#define GUIEXTRA \
	( (getenv ("RELION_QSUB_EXTRA_COUNT"))? \
	(std::max(0,(atoi(getenv ("RELION_QSUB_EXTRA_COUNT"))-4))*STEPY) : 0 )
#define MENUHEIGHT 30
#define TABHEIGHT 25
#define GUIWIDTH 800
#define GUIHEIGHT_OLD 420+GUIEXTRA
#define GUIHEIGHT_EXT_START 370+GUIEXTRA
#define GUIHEIGHT_EXT_START2 (GUIHEIGHT_EXT_START+MENUHEIGHT+10)
#define GUIHEIGHT_EXT 800+GUIEXTRA
#define XCOL0 200
#define WCOL0 200
#define XCOL1 ( (XCOL0) + 10  )
#define XCOL2 ( (XCOL0) + 280 )
#define XCOL3 ( (XCOL0) + 480 )
#define XCOL4 ( (XCOL0) + 495 )
#define XCOL5 ( (XCOL0) + 555 )
#define ENTRY_FONTSIZE 13
#define STEPY 20
#define COLUMN_SEPARATION 3
#define WCOL1 ( (XCOL2) - (XCOL1) - (COLUMN_SEPARATION) )
#define WCOL2 ( (XCOL3) - (XCOL2) - (COLUMN_SEPARATION) )
#define WCOL3 ( (XCOL4) - (XCOL3) - (COLUMN_SEPARATION) )
#define WCOL4 ( (XCOL5) - (XCOL4) - (COLUMN_SEPARATION) )
//version-1.0 #define GUI_BUTTON_COLOR (fl_rgb_color(200,255,100))
//version-1.0 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(255,155,0))
//version-1.1 #define GUI_BUTTON_COLOR (fl_rgb_color(50,150,250))
//version-1.1 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(255,50,50))
//version-1.2 #define GUI_BUTTON_COLOR (fl_rgb_color(155,150,255))
//version-1.2 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(205,53,100))
//version-1.3 #define GUI_BUTTON_COLOR (fl_rgb_color(50, 200, 50))
//version-1.3 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(255,80,80))
//version-1.4 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(205,40,150))
//version-1.4 #define GUI_BUTTON_COLOR (fl_rgb_color(60, 180, 155))
//version-1.4 #define GUI_BUTTON_DARK_COLOR (fl_rgb_color(45, 135, 120))
//version-2.0 #define GUI_BUTTON_COLOR (fl_rgb_color(0, 200, 255))
//version-2.0 #define GUI_BUTTON_DARK_COLOR (fl_rgb_color(0, 160, 200))
//version-2.0 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(70, 120, 255))
//version-2.1 #define GUI_BUTTON_COLOR (fl_rgb_color(100, 200, 50))
//version-2.1 #define GUI_BUTTON_DARK_COLOR (fl_rgb_color(70, 140, 30))
//version-2.1 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(0, 130, 0))
//version-3.0 #define GUI_BUTTON_COLOR (fl_rgb_color(255, 180, 132))
//version-3.0 #define GUI_BUTTON_DARK_COLOR (fl_rgb_color(250, 150, 124))
//version-3.0 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(235, 130, 0))
//version-3.1 #define GUI_BUTTON_COLOR (fl_rgb_color(238,130,238))
//version-3.1 #define GUI_BUTTON_DARK_COLOR (fl_rgb_color(200, 110, 200))
//version-3.1 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(170, 0, 170))
// Dont forget GUI  runbutton colour in src/displayer.h!
//version-3.2 
#define GUI_BUTTON_COLOR (fl_rgb_color(200,80,110))
#define GUI_BUTTON_DARK_COLOR (fl_rgb_color(170, 40, 70))
#define GUI_RUNBUTTON_COLOR (fl_rgb_color(160, 30, 60))
#define GUI_BACKGROUND_COLOR (fl_rgb_color(230,230,240)) // slightly blue because of blue buttons in 2.0!
#define GUI_BACKGROUND_COLOR2 (fl_rgb_color(180,180,190)) // slightly blue because of blue buttons in 2.0!
// devel-version
//#define GUI_BUTTON_COLOR (fl_rgb_color(255, 150, 150))
//#define GUI_BUTTON_DARK_COLOR (fl_rgb_color(200, 120, 120))
//#define GUI_RUNBUTTON_COLOR (fl_rgb_color(170, 0, 0))
//#define GUI_BACKGROUND_COLOR (fl_rgb_color(255,200,200)) // slightly red
//#define GUI_BACKGROUND_COLOR2 (fl_rgb_color(230,180,180)) // slightly red
//possible?#define GUI_BUTTON_COLOR (fl_rgb_color(50, 200, 255))
//devel-version
//possible #define GUI_RUNBUTTON_COLOR (fl_rgb_color(205,0,155))
#define GUI_INPUT_COLOR (fl_rgb_color(255,255,230))

#define TOGGLE_DEACTIVATE 0
#define TOGGLE_REACTIVATE 1
#define TOGGLE_ALWAYS_DEACTIVATE 2
#define TOGGLE_LEAVE_ACTIVE 3

static Fl_Menu_Item bool_options[] = {
			      {"Yes"},
			      {"No"},
			      {0} // this should be the last entry
			      };

// A text to Float converter that raises an error window.
float fltkTextToFloat(const char* str);

/** This class displays opens an additional window with (help) text
 *
 */
class ShowHelpText{

public:
    // Constructor that opens and displays the window
	ShowHelpText(const char* help = NULL);
	// Empty destructor
    ~ShowHelpText();
};


class GuiEntry{

public:

	// What to do upon toggle of continue button
	int deactivate_option;

	// Input value storage
	Fl_Input* inp;

	// JobOption
	JobOption joboption;

	// Button to show additional help text
	Fl_Button* help;

	////////////// FileName entry

	// Browse button
    Fl_Button* browse;

    ////////////// Radio entry

    // The choices
    Fl_Choice * choice;
    // The menu
    Fl_Menu_* menu;
    // Deactivate this group
    Fl_Group * my_deactivate_group;
	bool actually_activate;

    ////////////// Slider entry

    // The slider
    Fl_Slider * slider;

    /** Constructor with x,y-position from top left
	 *  wcol1, wcol2 and wcol3 are the widths of the three columns described above
	 *  title is the value displayed in the first column
	 *  defaultvalue is what will appear by default in the input value
	 *  help is the additional help text. If it is set to NULL, no help button will be displayed
	 */
    GuiEntry()
    {
    	deactivate_option = -1;
    	inp = NULL;
		help = NULL;
		browse = NULL;
		choice = NULL;
		menu = NULL;
		my_deactivate_group = NULL;
		actually_activate = false;
		slider = NULL;
    };

    /** Empty destructor
     */
	~GuiEntry() { clear(); }

    // Clear this entry
	void clear();

	/** Here really start the entry
	 */
	void initialise(int x, int y, Fl_Group * deactivate_this_group, bool actually_activate, int height, int wcol2, int wcol3);

	/** Place an entry on a window
	 */
	void place(JobOption &joboption, int &y, int _deactivate_option = TOGGLE_LEAVE_ACTIVE, Fl_Group * deactivate_this_group = NULL, bool actually_activate = false,
	           int x = XCOL2, int h = STEPY, int wcol2 = WCOL2, int wcol3 = WCOL3 );

    // Set _value in the Fl_Input on the GUI, and also in the joboptions. Also update menu/slider if necessary
    void setValue(std::string _value);

    // Deactivate this entry if the input boolean is true
    void deactivate(bool do_deactivate = true);

    /** Call-back functions for the help button
     *  The method of using two functions of static void and inline void was copied from:
     *  http://www3.telus.net/public/robark/
     */
    static void cb_help(Fl_Widget*, void*);
    void cb_help_i();

    // Call-back functions for the browse button
    static void cb_browse(Fl_Widget*, void*);
    void cb_browse_i();

    // Call-back functions for the browse button
    static void cb_browse_node(Fl_Widget*, void*);
    void cb_browse_node_i();

    // Call-back functions for the menu
    static void cb_menu(Fl_Widget*, void*);
    void cb_menu_i();

    // Call-back functions for the slider
    static void cb_slider(Fl_Widget*, void*);
    void cb_slider_i();

    static void cb_input(Fl_Widget*, void*);
    void cb_input_i();

};




#endif /* SRC_NEWGUI_ENTRIES_H_ */
