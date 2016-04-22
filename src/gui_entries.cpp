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
#include "src/gui_entries.h"

// This allows CURRENT_ODIR browse buttons
std::string current_browse_directory;

bool replaceStringOnce(Fl_Text_Buffer *textbuf, std::string findthis, std::string replaceby)
{
	const char *find = findthis.c_str();
	const char *replace = replaceby.c_str();

	// Loop through the whole string
	int pos = 0;
	if (textbuf->search_forward(pos, find, &pos))
	{
		// Found a match; update the position and replace text...
		textbuf->select(pos, pos+strlen(find));
		textbuf->remove_selection();
		textbuf->insert(pos, replace);
		pos += strlen(replace);
		return true;
	}
	else
	{
		return false;
	}
}

void replaceStringAll(Fl_Text_Buffer *textbuf, std::string findthis, std::string replaceby)
{

	bool do_search_more = true;
	while (do_search_more)
		do_search_more = replaceStringOnce(textbuf, findthis, replaceby);
}

void appendLineString(Fl_Text_Buffer *textbuf, std::string copylinewiththis, int times)
{
	const char *find = copylinewiththis.c_str();
	// Loop through the whole string
	int pos = 0;
	char * sel;

	if (textbuf->search_forward(pos, find, &pos))
	{
		// Found a match; get the entire line and append at the end
		int line0 = textbuf->line_start(pos);
		int lineF = textbuf->line_end(pos);
		textbuf->select(line0,lineF);
		sel = textbuf->selection_text();
	}
	else
	{
		std::cerr <<" appendLineString ERROR: Could not find" << copylinewiththis << " in textfile..." << std::endl;
		exit(1);
	}

	for (int n = 0; n < times; n++)
	{
		textbuf->append("\n");
		textbuf->append(sel);
		textbuf->append("\n");
	}
	textbuf->append("\n");

}

ShowHelpText::ShowHelpText(const char *help)
{
    int w=640;
    int h=480;
	Fl_Window *win = new Fl_Window(w, h);
    Fl_Text_Buffer *buff = new Fl_Text_Buffer();
    Fl_Text_Display *disp = new Fl_Text_Display(20, 20, w-40, h-40, "relion additional text.");
    disp->buffer(buff);
    disp->wrap_mode(1,79);
    win->resizable(*disp);
    win->show();
    buff->text(help);
}

ShowHelpText::~ShowHelpText(){};


// ==============================================================================
// AnyEntry =====================================================================
// ==============================================================================

void AnyEntry::initialise(int x, int y, int height,
				   int wcol2, int wcol3,
				   const char* title,
				   const char* defaultvalue,
				   const char* helptext)
{
    // Set the label
    label = title;

    // The input field
	if (defaultvalue != NULL)
	{
		inp = new Fl_Input(XCOL2, y, wcol2, height, title);

		// Set the input value
		inp->value(defaultvalue);
		inp->color(GUI_INPUT_COLOR);
		inp->textsize(ENTRY_FONTSIZE);
		inp->labelsize(ENTRY_FONTSIZE);
	}

	// Display help button if needed
    if (helptext != NULL)
    {
    	// Set the help text
		myhelptext = helptext;

		// The Help button
		help = new Fl_Button( XCOL3, y, wcol3, height, "?");
		help->callback( cb_help, this );
		help->color(GUI_BUTTON_COLOR);
		help->labelsize(ENTRY_FONTSIZE);
    }
}

void AnyEntry::place(int &y,
		const char * title,
		const char* defaultvalue,
		const char* helptext,
		int x, int h, int wcol2, int wcol3 )
{

	// Clear if existing
	clear();

	// Add the entry to the window
	initialise(x, y, h, wcol2, wcol3, title, defaultvalue, helptext);

	// Update the Y-coordinate
    y += h + 2;

}

void AnyEntry::placeOnSameYPosition(int y,
		const char * title,
		const char * title_full,
		const char* defaultvalue,
		const char* helptext,
		int x, int h, int wcol2, int wcol3 )
{

	// Clear if existing
	clear();

    // Set the label
    label = title;
    if ( (title_full) && (strlen(title_full) > 0) )
    	label_full = title_full;

    // The input field
	if (defaultvalue != NULL)
	{
		inp = new Fl_Input(x, y, wcol2, h, title);

		// Set the input value
		inp->value(defaultvalue);
		inp->color(GUI_INPUT_COLOR);
		inp->textsize(ENTRY_FONTSIZE);
		inp->labelsize(ENTRY_FONTSIZE);
	}

	// Display help button if needed
    if (helptext != NULL)
    {
    	// Set the help text
		myhelptext = helptext;

		// The Help button
		help = new Fl_Button(x + wcol2 + COLUMN_SEPARATION, y, wcol3, h, "?");
		help->callback( cb_help, this );
		help->color(GUI_BUTTON_COLOR);
		help->labelsize(ENTRY_FONTSIZE);
    }
}

std::string AnyEntry::getValue()
{
	return (std::string)inp->value();
}

void AnyEntry::setValue(const char* val)
{
	inp->value(val);
}

void AnyEntry::writeValue(std::ostream& out)
{
	// Only write entries that have been initialised
	if (label_full != "")
		out << label_full << " == " << getValue() << std::endl;
	else if (label != "")
		out << label << " == " << getValue() << std::endl;
}


void AnyEntry::readValue(std::ifstream& in)
{
	std::string label_saved;
	if (label_full != "")
		label_saved = label_full;
	else
		label_saved = label;

    if (label_saved != "")
    {
		// Start reading the ifstream at the top
		in.clear(); // reset eof if happened...
    	in.seekg(0, std::ios::beg);
		std::string line;
		while (getline(in, line, '\n'))
		{
			if (line.rfind(label_saved) == 0)
			{
				// found my label
				int equalsigns = line.rfind("==");
				std::string newval = line.substr(equalsigns + 3, line.length() - equalsigns - 3);
				inp->value(newval.c_str());
				return;
			}
		}
    }
}

void AnyEntry::clear()
{
	if (label != "")
	{
		label="";
		if (inp)
		{
			delete inp;
			inp = NULL;
		}
		if (help)
		{
			delete help;
			help = NULL;
		}
	}
}

void AnyEntry::deactivate(bool do_deactivate)
{
	if (do_deactivate)
	{
		if (inp)
			inp->deactivate();
		if (help)
			help->deactivate();
	}
	else
	{
		if (inp)
			inp->activate();
		if (help)
			help->activate();
	}

}

// Help button call-back functions
void AnyEntry::cb_help(Fl_Widget* o, void* v) {

    AnyEntry* T=(AnyEntry*)v;
    T->cb_help_i();
}

void AnyEntry::cb_help_i() {

    ShowHelpText *help = new ShowHelpText(myhelptext);

}


// ==============================================================================
// FileNameEntry ================================================================
// ==============================================================================

void FileNameEntry::initialise(int x, int y, int height,
		                     int wcol2, int wcol3, int wcol4,
		                     const char* title,
		                     const char* defaultvalue,
		                     const char* _pattern,
		                     const char* _directory,
		                     const char* helptext)
{

	AnyEntry::initialise(x,y,height,
			wcol2,wcol3,
			title,
			defaultvalue,
			helptext);

	// Store the pattern for the file chooser
	pattern = _pattern;
	directory = _directory;

    // The Browse button
    browse = new Fl_Button( XCOL4, y, WCOL4, height, "Browse");
    browse->callback( cb_browse, this );
    browse->color(GUI_BUTTON_COLOR);
    browse->labelsize(ENTRY_FONTSIZE);
}


void FileNameEntry::place(int &y,
		const char * title,
		const char* defaultvalue,
		const char* pattern,
		const char* directory,
		const char* helptext,
		int x, int h, int wcol2, int wcol3, int wcol4)
{

	// Clear if existing
	clear();

	// Add the entry to the window
	initialise(x, y, h, wcol2, wcol3, wcol4, title, defaultvalue, pattern, directory, helptext);

	// Update the Y-coordinate
    y += h + 2;

}

void FileNameEntry::clear()
{
	if (label != "")
	{
		AnyEntry::clear();
		delete browse;
	}
}

void FileNameEntry::deactivate(bool do_deactivate)
{
	AnyEntry::deactivate(do_deactivate);
	if (do_deactivate)
	{
		browse->deactivate();
	}
	else
	{
		browse->activate();
	}

}

void FileNameEntry::cb_browse(Fl_Widget* o, void* v) {

    FileNameEntry* T=(FileNameEntry*)v;
    T->cb_browse_i();
}


void FileNameEntry::cb_browse_i() {

    Fl::scheme("gtk+");
    Fl_File_Chooser * G_chooser = new Fl_File_Chooser("", pattern, Fl_File_Chooser::SINGLE, "");

    std::string test="";
    if (directory != NULL)
    {
        std::string test2(directory);
        test = test2;
    }
    if (test=="CURRENT_ODIR")
    	G_chooser->directory(current_browse_directory.c_str());
    else
    	G_chooser->directory(directory);
    G_chooser->color(GUI_BACKGROUND_COLOR);
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

    char relname[FL_PATH_MAX];
    fl_filename_relative(relname,sizeof(relname),G_chooser->value());
    inp->value(relname);
}


// ==============================================================================
// InputNodeEntry ================================================================
// ==============================================================================

void InputNodeEntry::initialise(int x, int y, int height,
		                     int wcol2, int wcol3, int wcol4,
		                     const char* title,
		                     int _type,
		                     const char* defaultvalue,
		                     const char* _pattern,
		                     const char* helptext)
{

	AnyEntry::initialise(x,y,height,
			wcol2,wcol3,
			title,
			defaultvalue,
			helptext);

	// Store the pattern for the file chooser
	pattern = _pattern;
	type = _type;
    // The Browse button
    browse = new Fl_Button( XCOL4, y, WCOL4, height, "Browse");
    browse->callback( cb_browse_node, this );
    browse->color(GUI_BUTTON_COLOR);
    browse->labelsize(ENTRY_FONTSIZE);
}


void InputNodeEntry::place(int &y,
		const char * title,
		int _type,
		const char* defaultvalue,
		const char* pattern,
		const char* helptext,
		int x, int h, int wcol2, int wcol3, int wcol4 )
{

	// Clear if existing
	clear();

	// Add the entry to the window
	initialise(x, y, h, wcol2, wcol3, wcol4, title, _type, defaultvalue, pattern, helptext);

	// Update the Y-coordinate
    y += h + 2;

}

void InputNodeEntry::clear()
{
	// TODO: add stuff to track history here
	if (label != "")
	{
		AnyEntry::clear();
		delete browse;
	}
}

void InputNodeEntry::deactivate(bool do_deactivate)
{
	// TODO: add stuff to track history here
	AnyEntry::deactivate(do_deactivate);
	if (do_deactivate)
	{
		browse->deactivate();
	}
	else
	{
		browse->activate();
	}

}

void InputNodeEntry::cb_browse_node(Fl_Widget* o, void* v) {

    InputNodeEntry* T=(InputNodeEntry*)v;
    T->cb_browse_node_i();
}


void InputNodeEntry::cb_browse_node_i() {

    Fl::scheme("gtk+");
    Fl_File_Chooser * G_chooser = new Fl_File_Chooser("", pattern, Fl_File_Chooser::SINGLE, "");

    std::string fn_dir = ".Nodes/" + integerToString(type);
    G_chooser->directory(fn_dir.c_str());
    G_chooser->color(GUI_BACKGROUND_COLOR);
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

    char relname[FL_PATH_MAX];
    fl_filename_relative(relname,sizeof(relname),G_chooser->value());

    // Get rid of the .Nodes/type/ directory-name again
    std::string replace = std::string(relname);
    std::string replace2 = (std::string::npos == replace.find(fn_dir.c_str())) ? replace : replace.substr(fn_dir.length()+1, replace.length());
    char relname2[FL_PATH_MAX];
    strcpy(relname2, replace2.c_str());

    inp->value(relname2);
}


// ==============================================================================
// RadioEntry ================================================================
// ==============================================================================
void RadioEntry::initialise(int x, int y, int height,
		                     int wcol2, int wcol3, int wcol4,
		                     const char* title,
		                     Fl_Menu_Item *options,
		                     Fl_Menu_Item* defaultvalue,
		                     const char* helptext,
		    				 Fl_Group * deactivate_this_group)
{
	AnyEntry::initialise(x,y,height,
			wcol2,wcol3,
			title,
			defaultvalue->label(),
			helptext);

    // Pull-down menu button
	//Fl_File_Chooser * G_chooser = new Fl_File_Chooser("", "", Fl_File_Chooser::SINGLE, "");

	my_deactivate_group = deactivate_this_group;
	choice = new Fl_Choice(XCOL2, y, WCOL2, height);
    choice->menu(options);
    choice->picked(defaultvalue);
    choice->callback(cb_menu, this);
    choice->textsize(ENTRY_FONTSIZE);

    menu = choice;
    //menu->color(GUI_BACKGROUND_COLOR);
    menu->color(GUI_INPUT_COLOR);
    menu->textsize(ENTRY_FONTSIZE);
}


void RadioEntry::place(int &y,
			const char * title,
			Fl_Menu_Item *options,
			Fl_Menu_Item* defaultvalue,
			const char* helptext,
			int x, int h, int wcol2, int wcol3, int wcol4 )

{
    // Clear if existing
	clear();

	// Add the entry to the window
	initialise(x, y, h, wcol2, wcol3, wcol4, title, options, defaultvalue, helptext);

    // Update the Y-coordinate
    y += h + 2;

}

void RadioEntry::clear()
{
	if (label != "")
	{
		AnyEntry::clear();
		//delete choice;
		delete menu;
	}
}

void RadioEntry::deactivate(bool do_deactivate)
{
	AnyEntry::deactivate(do_deactivate);
	if (do_deactivate)
	{
		menu->deactivate();
	}
	else
	{
		menu->activate();
	}

}

std::string RadioEntry::getValue()
{
	return (std::string)inp->value();
}

void RadioEntry::readValue(std::ifstream& in)
{
	if (label != "")
	{
		AnyEntry::readValue(in);
		const Fl_Menu_Item *p = choice->find_item(inp->value());
		if ( p )
			choice->picked(p);
		else
			std::cerr << "Error readValue: Menu item not found:" << inp->value()<< std::endl;
	}
}


void RadioEntry::cb_menu(Fl_Widget* o, void* v) {

    RadioEntry* T=(RadioEntry*)v;
    T->cb_menu_i();
}


void RadioEntry::cb_menu_i()
{

	const Fl_Menu_Item* m = menu->mvalue();
	// Set my own value
	inp->value(m->label());

	// In case this was a boolean that deactivates a group, do so:
	if (my_deactivate_group != NULL)
	if (strcmp(inp->value(), "No") == 0)
		my_deactivate_group->deactivate();
	else
		my_deactivate_group->activate();

}

// ==============================================================================
// BooleanEntry ================================================================
// ==============================================================================
void BooleanEntry::initialise(int x, int y, int height,
		                     int wcol2, int wcol3, int wcol4,
		                     const char* title,
		                     bool defaultvalue,
		                     const char* helptext,
		    				 Fl_Group * deactivate_this_group)
{

	Fl_Menu_Item* defval;

	if (defaultvalue)
		defval = &bool_options[0];
	else
		defval = &bool_options[1];
	RadioEntry::initialise(x,y,height,
			wcol2,wcol3,wcol4,
			title,
			bool_options,
			defval,
			helptext,
			deactivate_this_group);

}

void BooleanEntry::place(int &y,
		const char * title,
		bool defaultvalue,
		const char* helptext,
		Fl_Group * deactivate_this_group,
		int x, int h, int wcol2, int wcol3, int wcol4 )
{

    // Clear if existing
	clear();

	// Add the entry to the window
	initialise(x, y, h, wcol2, wcol3, wcol4, title, defaultvalue, helptext, deactivate_this_group);

    // Update the Y-coordinate
    y += h + 2;


}

bool BooleanEntry::getValue()
{
	if (strcmp(inp->value(), "Yes") == 0)
		return true;
	else
		return false;
}

// ==============================================================================
// SliderEntry ================================================================
// ==============================================================================
void SliderEntry::initialise(int x, int y, int height,
		                     int wcol2, int wcol3, int wcol4,
		                     const char* title,
		                     float defaultvalue,
		                     float minvalue,
		                     float maxvalue,
		                     float valuestep,
		                     const char* helptext)
{

	int floatwidth = 50;
	AnyEntry::initialise(x,y,height,
			floatwidth,wcol3,
			title,
			"",
			helptext);

	// Initialise label
	label = title;

	// Slider is shorter than wcol2, so that underlying input field becomes visible
	slider = new Fl_Slider(XCOL2 + floatwidth, y, wcol2 - floatwidth, height);
	slider->type(1);
	slider->callback(cb_slider, this);
	slider->minimum(minvalue);
	slider->maximum(maxvalue);
	slider->step(valuestep);
	slider->type(FL_HOR_NICE_SLIDER);
	slider->color(GUI_BACKGROUND_COLOR);
	inp->callback(cb_input, this);
	inp->when(FL_WHEN_ENTER_KEY|FL_WHEN_NOT_CHANGED);

	// Set the default in the input and the slider:
	std::string str = floatToString(defaultvalue);
	inp->value(str.c_str());
	slider->value(defaultvalue);

}

void SliderEntry::place(int &y,
		const char* title,
		float defaultvalue,
		float minvalue,
		float maxvalue,
		float valuestep,
		const char* helptext,
		int x, int h, int wcol2, int wcol3, int wcol4 )
{

    // Clear if existing
	clear();

	// Add the entry to the window
	initialise(x, y, h, wcol2, wcol3, wcol4, title, defaultvalue, minvalue, maxvalue, valuestep, helptext);

    // Update the Y-coordinate
    y += h + 2;

}

void SliderEntry::clear()
{
	if (label != "")
	{
		AnyEntry::clear();
		delete slider;
	}
}

void SliderEntry::deactivate(bool do_deactivate)
{
	AnyEntry::deactivate(do_deactivate);
	if (do_deactivate)
	{
		slider->deactivate();
	}
	else
	{
		slider->activate();
	}

}


float SliderEntry::getValue()
{
	return textToFloat(inp->value());
}

void SliderEntry::readValue(std::ifstream& in)
{
	if (label != "")
	{
		AnyEntry::readValue(in);
		// Also reset the slider
		slider->value(textToFloat(inp->value()));
	}
}

void SliderEntry::cb_slider(Fl_Widget* o, void* v) {

    SliderEntry* T=(SliderEntry*)v;
    T->cb_slider_i();
}


void SliderEntry::cb_slider_i() {

    static int recurse = 0;
    if ( recurse ) {
        return;
    } else {
        recurse = 1;
        std::string str = floatToString(slider->value());
        inp->value(str.c_str());
        slider->redraw();
        recurse = 0;
    }
}

void SliderEntry::cb_input(Fl_Widget* o, void* v) {

    SliderEntry* T=(SliderEntry*)v;
    T->cb_input_i();
}


void SliderEntry::cb_input_i() {

    static int recurse = 0;
    if ( recurse ) {
        return;
    } else {
        recurse = 1;
        //slider->value(textToFloat(my_input->value()));         // pass input's value to slider
        slider->value(textToFloat(inp->value()));         // pass input's value to slider
        //inp->value(my_input->value()); // also set the normal input for getValue!
        recurse = 0;
    }
}


