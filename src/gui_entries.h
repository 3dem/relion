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
#ifndef GUI_ENTRIES_H_
#define GUI_ENTRIES_H_
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
#include "src/macros.h"
#include "src/strings.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <cstdio>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Gui layout
//#define XCOL1 10
//#define XCOL2 260
//#define XCOL3 460
//#define XCOL4 475
//#define XCOL5 535
#define GUIWIDTH 800
#define GUIHEIGHT_OLD 400
#define GUIHEIGHT_EXT 600
#define XCOL0 200
#define WCOL0 200
#define XCOL1 ( (XCOL0) + 10  )
#define XCOL2 ( (XCOL0) + 260 )
#define XCOL3 ( (XCOL0) + 460 )
#define XCOL4 ( (XCOL0) + 475 )
#define XCOL5 ( (XCOL0) + 535 )
#define STEPY 22
#define COLUMN_SEPARATION 3
#define WCOL1 ( (XCOL2) - (XCOL1) - (COLUMN_SEPARATION) )
#define WCOL2 ( (XCOL3) - (XCOL2) - (COLUMN_SEPARATION) )
#define WCOL3 ( (XCOL4) - (XCOL3) - (COLUMN_SEPARATION) )
#define WCOL4 ( (XCOL5) - (XCOL4) - (COLUMN_SEPARATION) )
//version-1.0 #define GUI_BUTTON_COLOR (fl_rgb_color(200,255,100))
//version-1.1 #define GUI_BUTTON_COLOR (fl_rgb_color(50,150,250))
//version-1.2 #define GUI_BUTTON_COLOR (fl_rgb_color(155,150,255))
//version-1.3 #define GUI_BUTTON_COLOR (fl_rgb_color(50, 200, 50))
//version-1.4 #define GUI_BUTTON_COLOR (fl_rgb_color(60, 180, 155))
//version-1.4 #define GUI_BUTTON_DARK_COLOR (fl_rgb_color(45, 135, 120))
//devel-version
#define GUI_BUTTON_COLOR (fl_rgb_color(0, 235, 235))
#define GUI_BUTTON_DARK_COLOR (fl_rgb_color(0, 180, 180))
//possible?#define GUI_BUTTON_COLOR (fl_rgb_color(50, 200, 255))
//version-1.0 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(255,155,0))
//version-1.1 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(255,50,50))
//version-1.2 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(205,53,100))
//version-1.3 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(255,80,80))
//version-1.4 #define GUI_RUNBUTTON_COLOR (fl_rgb_color(205,40,150))
//devel-version
#define GUI_RUNBUTTON_COLOR (fl_rgb_color(70, 120, 255))
//possible #define GUI_RUNBUTTON_COLOR (fl_rgb_color(205,0,155))
#
#define GUI_BACKGROUND_COLOR (fl_rgb_color(240,240,240))
#define GUI_BACKGROUND_COLOR2 (fl_rgb_color(200,200,200))
#define GUI_INPUT_COLOR (fl_rgb_color(255,255,230))

// Replace a single instance of text in a buffer. Return true if replaced, false otherwise
bool replaceStringOnce(Fl_Text_Buffer *textbuf, std::string replacethis, std::string replaceby);
// General utility to replace strings in a text buffer.
void replaceStringAll(Fl_Text_Buffer *textbuf, std::string replacethis, std::string replaceby);
void appendLineString(Fl_Text_Buffer *textbuf, std::string copylinewiththis, int times);

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


static Fl_Menu_Item bool_options[] = {
			      {"Yes"},
			      {"No"},
			      {0} // this should be the last entry
			      };

class textOnlyEntry{

public:
	Fl_Text_Display* mydisp;
	Fl_Text_Buffer *textbuff;
	bool has_been_set;

	textOnlyEntry()
	{
		has_been_set=false;
	}
	~textOnlyEntry(){};

	void initialise(int x, int y, int width, int height, const char* text)
	{
		mydisp = new Fl_Text_Display(XCOL1, y, width, height);
		textbuff = new Fl_Text_Buffer();
		textbuff->text(text);
		mydisp->buffer(textbuff);
		mydisp->color(GUI_BACKGROUND_COLOR);
		has_been_set=true;
	}

	void place(int &y,
			const char* text,
			int width= WCOL1 + WCOL2 + WCOL3, int height = STEPY + 6, int x = XCOL1)
	{
	    // Clear if existing
		clear();

		// Add the entry to the window
		// Add 3 to step_y, otherwise the text may not fit...
		initialise(x, y, width, height, text);

	    // Update the Y-coordinate
	    y += height + 2;


	}

	void clear()
	{
		if (has_been_set)
		{
			delete mydisp;
			delete textbuff;
			has_been_set = false;
		}
	}
};

/** This is the main class to generate input entry-lines in the Gui windows.
 *  It implements three columns to be displayed:
 *  1. box with the label
 *  2. Input field with the input value
 *  3. Help button that pops up a window with additional help text
 *
 *  All specific entries (e.g. to get FileName, Boolean, etc. inherit from this class)
 *
 *
 */
class AnyEntry{

public:
    // Input value storage
	Fl_Input* inp;

	// Label
	std::string label;

    // Button to show additional help text
	Fl_Button* help;

	// The additional help text
    const char *myhelptext;

    /** Constructor with x,y-position from top left
	 *  wcol1, wcol2 and wcol3 are the widths of the three columns described above
	 *  title is the value displayed in the first column
	 *  defaultvalue is what will appear by default in the input value
	 *  help is the additional help text. If it is set to NULL, no help button will be displayed
	 */
	AnyEntry(){};

    /** Empty destructor
     */
	~AnyEntry(){};

	/** Here really start the entry
	 */
	void initialise(int x, int y, int height, int wcol2, int wcol3, const char* title, const char* defaultvalue = NULL, const char* help = NULL);

	/** Place an entry on a window
	 */
	void place(int &y,
				const char * title,
				const char* defaultvalue = NULL,
				const char* helptext = NULL,
				int x = XCOL1, int h = STEPY, int wcol2 = WCOL2, int wcol3 = WCOL3 );

	// Get the value
    std::string getValue();

    // Set the value
    void setValue(const char* inp);

    // Clear this entry
	void clear();

    // Deactivate this entry if the input boolean is true
    void deactivate(bool do_deactivate = true);

    // Save the value to a file
    void writeValue(std::ostream& out);

    // Read the value from a file
    void readValue(std::ifstream& in);

    /** Call-back functions for the help button
     *  The method of using two functions of static void and inline void was copied from:
     *  http://www3.telus.net/public/robark/
     */
    static void cb_help(Fl_Widget*, void*);
    void cb_help_i();
};


// Get a FileName value from the user (with browse button).
class FileNameEntry: public AnyEntry
{

public:
	// Browse button
    Fl_Button* browse;

    const char* pattern;

    // Constructor (with 4 column widths)
	FileNameEntry() {};

    // Destructor
	~FileNameEntry(){};

	void initialise(int x, int y, int height,
    		int wcol2, int wcol3, int wcol4,
    		const char* title,
    		const char* defaultvalue,
    		const char* _pattern = "",
    		const char* help = NULL);

	// places on one the window
	void place(int &y,
				const char * title,
				const char* defaultvalue = NULL,
				const char* pattern = "",
				const char* helptext = NULL,
				int x = XCOL1, int h = STEPY, int wcol2 = WCOL2, int wcol3 = WCOL3, int wcol4 = WCOL4 );

    // Clear this entry
	void clear();

	// Deactivate this entry if the input boolean is true
    void deactivate(bool do_deactivate = true);


private:
    // Call-back functions for the browse button
    static void cb_browse(Fl_Widget*, void*);
    void cb_browse_i();

};

class InputNodeEntry: public FileNameEntry
{

public:
    // Type of this node
    int type;

    //Output from which process?
    std::string output_from;

    // Constructor (with 4 column widths)
	InputNodeEntry() {};

    // Destructor
	~InputNodeEntry(){};

	void initialise(int x, int y, int height,
    		int wcol2, int wcol3, int wcol4,
    		const char* title,
    		int _type,
    		const char* defaultvalue,
    		const char* _pattern = "",
    		const char* help = NULL);

	// places on one the window
	void place(int &y,
				const char * title,
				int _type,
				const char* defaultvalue = NULL,
				const char* pattern = "",
				const char* helptext = NULL,
				int x = XCOL1, int h = STEPY, int wcol2 = WCOL2, int wcol3 = WCOL3, int wcol4 = WCOL4 );

    // Clear this entry
	void clear();

	// Deactivate this entry if the input boolean is true
    void deactivate(bool do_deactivate = true);


private:
    // Call-back functions for the browse button
    static void cb_browse_node(Fl_Widget*, void*);
    void cb_browse_node_i();

};



// Get an entry from a list of possible values from the user.
class RadioEntry: public AnyEntry
{
public:

    // The choices
    Fl_Choice * choice;
    // The menu
    Fl_Menu_* menu;
    // Deactivate this group
    Fl_Group * my_deactivate_group;

    // Constructor
    RadioEntry(){};

    // Destructor
    ~RadioEntry(){};

    void initialise(int x, int y, int height,
				 int wcol2, int wcol3, int wcol4,
				 const char* title,
				 Fl_Menu_Item *options,
				 Fl_Menu_Item* defaultvalue,
				 const char* help = NULL,
				 Fl_Group * deactivate_this_group = NULL);

    void place(int &y,
				const char * title,
				Fl_Menu_Item *options,
				Fl_Menu_Item* defaultvalue,
				const char* helptext = NULL,
				int x = XCOL1, int h = STEPY, int wcol2 = WCOL2, int wcol3 = WCOL3, int wcol4 = WCOL4 );

    // Clear this entry
	void clear();

	// Deactivate this entry if the input boolean is true
    void deactivate(bool do_deactivate = true);

    // Get the value
    std::string getValue();

    // Read the value from a file
    void readValue(std::ifstream& in);

	void call_menu_i()
    {
        cb_menu_i();
    }

public: // this one is public so that it can be called in mainwindow to deactivate default groups
    static void cb_menu(Fl_Widget*, void*);
    void cb_menu_i();
};

class BooleanEntry: public RadioEntry
{
public:
	// Constructor
	BooleanEntry(){};

	// Destructor
    ~BooleanEntry(){};

    void initialise(int x, int y, int height,
				 int wcol2, int wcol3, int wcol4,
				 const char* title,
				 bool defaultvalue,
				 const char* help = NULL,
				 Fl_Group * deactivate_this_group = NULL);

    void place(int &y,
				const char * title,
				bool defaultvalue = true,
				const char* help = NULL,
				Fl_Group * deactivate_this_group = NULL,
				int x = XCOL1, int h = STEPY, int wcol2 = WCOL2, int wcol3 = WCOL3, int wcol4 = WCOL4 );

    // Get the value
    bool getValue();


};

class SliderEntry:  public RadioEntry
{

public:
    // The slider
    Fl_Slider * slider;

    // Constructor
	SliderEntry(){};

	// Destructor
    ~SliderEntry(){};

    void initialise(int x, int y, int height,
				 int wcol2, int wcol3, int wcol4,
				 const char* title,
				 float defaultvalue,
                 float minvalue,
                 float maxvalue,
                 float valuestep,
				 const char* help = NULL);

    void place(int &y,
				const char* title,
				float defaultvalue,
				float minvalue,
				float maxvalue,
				float valuestep,
				const char* help,
				int x = XCOL1, int h = STEPY, int wcol2 = WCOL2, int wcol3 = WCOL3, int wcol4 = WCOL4 );

    // Clear this entry
	void clear();

	// Deactivate this entry if the input boolean is true
    void deactivate(bool do_deactivate = true);

    // Get the value
    float getValue();

    // Read the value from a file
    void readValue(std::ifstream& in);


private:
    static void cb_slider(Fl_Widget*, void*);
    void cb_slider_i();

    static void cb_input(Fl_Widget*, void*);
    void cb_input_i();


};


#endif /* GUI_ENTRIES_H_ */
