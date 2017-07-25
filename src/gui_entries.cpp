

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

float fltkTextToFloat(const char* str)
{
	float result = -999.;
    if (str == NULL)
    	fl_message("ERROR: NULL entry for TextToFloat conversion. Check your inputs!");
    else if (!sscanf(str, "%f", &result))
    	fl_message("ERROR: Invalid (non-numerical?) entry for TextToFloat conversion. Check your inputs!");

    return result;
}

// This allows CURRENT_ODIR browse buttons
std::string current_browse_directory;

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

void GuiEntry::clear()
{
	deactivate_option = -1;
	joboption.clear();
	/* This only gives segfaults....
	if (inp != NULL)
	{
		delete inp;
		inp = NULL;
	}
	if (help != NULL)
	{
		delete help;
		help = NULL;
	}
	if (browse != NULL)
	{
		delete browse;
		browse = NULL;
	}
	if (choice != NULL)
	{
		delete choice;
		choice = NULL;
	}
	if (menu != NULL)
	{
		delete menu;
		menu = NULL;
	}
	if (my_deactivate_group != NULL)
	{
		delete my_deactivate_group;
		my_deactivate_group = NULL;
	}
	if (slider != NULL)
	{
		delete slider;
		slider = NULL;
	}
	*/

}
void GuiEntry::initialise(int x, int y, Fl_Group * deactivate_this_group, int height, int wcol2, int wcol3)
{

    // The input field
	int mywidth = (joboption.joboption_type == JOBOPTION_SLIDER) ? 50 : wcol2;
	inp = new Fl_Input(x, y, mywidth, height, joboption.label_gui.c_str());
	inp->color(GUI_INPUT_COLOR);
	inp->textsize(ENTRY_FONTSIZE);
	inp->labelsize(ENTRY_FONTSIZE);
	inp->value(joboption.default_value.c_str());

	// Display help button if needed
    if (joboption.helptext != "")
    {
		// The Help button
		help = new Fl_Button( XCOL3, y, wcol3, height, "?");
		help->callback( cb_help, this );
		help->color(GUI_BUTTON_COLOR);
		help->labelsize(ENTRY_FONTSIZE);
    }

    if (joboption.joboption_type == JOBOPTION_FILENAME)
    {
        // The Browse button
        browse = new Fl_Button( XCOL4, y, WCOL4, height, "Browse");
        browse->callback( cb_browse, this );
        browse->color(GUI_BUTTON_COLOR);
        browse->labelsize(ENTRY_FONTSIZE);

    }
    else if (joboption.joboption_type == JOBOPTION_INPUTNODE)
    {

        // The Browse button
        browse = new Fl_Button( XCOL4, y, WCOL4, height, "Browse");
        browse->callback( cb_browse_node, this );
        browse->color(GUI_BUTTON_COLOR);
        browse->labelsize(ENTRY_FONTSIZE);
    }
    else if (joboption.joboption_type == JOBOPTION_RADIO || joboption.joboption_type == JOBOPTION_BOOLEAN)
    {

		choice = new Fl_Choice(XCOL2, y, WCOL2, height);
		if (joboption.joboption_type == JOBOPTION_RADIO)
		{
			if (joboption.radio_menu == RADIO_SAMPLING)
			{
				choice->menu(fl_sampling_options);
				for (int i = 0; i < 9; i++)
					if (std::string(job_sampling_options[i]) == joboption.default_value)
						choice->picked(&fl_sampling_options[i]);
			}
			else if (joboption.radio_menu == RADIO_NODETYPE)
			{
				choice->menu(fl_node_type_options);
				for (int i = 0; i < 10; i++)
					if (std::string(job_nodetype_options[i]) == joboption.default_value)
						choice->picked(&fl_node_type_options[i]);
			}
			else
				REPORT_ERROR("BUG: unrecognised radio menu type.");

		}
		else // boolean
		{
			if (deactivate_this_group != NULL)
				my_deactivate_group = deactivate_this_group;

			choice->menu(bool_options);
			if (joboption.default_value=="Yes")
				choice->picked(&bool_options[0]);
			else
				choice->picked(&bool_options[1]);
		}
		choice->callback(cb_menu, this);
		choice->textsize(ENTRY_FONTSIZE);

		menu = choice;
		//menu->color(GUI_BACKGROUND_COLOR);
		menu->color(GUI_INPUT_COLOR);
		menu->textsize(ENTRY_FONTSIZE);
    }
    else if (joboption.joboption_type == JOBOPTION_SLIDER)
    {
    	int floatwidth = 50;
    	// Slider is shorter than wcol2, so that underlying input field becomes visible
    	slider = new Fl_Slider(XCOL2 + floatwidth, y, wcol2 - floatwidth, height);
    	slider->type(1);
    	slider->callback(cb_slider, this);
    	slider->minimum(joboption.min_value);
    	slider->maximum(joboption.max_value);
    	slider->step(joboption.step_value);
    	slider->type(FL_HOR_NICE_SLIDER);
    	slider->color(GUI_BACKGROUND_COLOR);
    	inp->callback(cb_input, this);
    	inp->when(FL_WHEN_ENTER_KEY|FL_WHEN_NOT_CHANGED);

    	// Set the default in the input and the slider:
    	inp->value(joboption.default_value.c_str());
    	slider->value(textToDouble(joboption.default_value));

    }


}
void GuiEntry::place(JobOption &_joboption, int &y, int _deactivate_option, Fl_Group * deactivate_this_group, bool _do_oldstyle, int x, int h, int wcol2, int wcol3 )
{

	// Clear if existing
	clear();

	// What to do when continue is toggled
	deactivate_option = _deactivate_option;

	joboption = _joboption;

	do_oldstyle = _do_oldstyle;

	// Add the entry to the window
	initialise(x, y, deactivate_this_group, h, wcol2, wcol3);

	// Update the Y-coordinate
    y += h + 2;

}

// Set the value back from the Fl_Input into the JobOption.value
void GuiEntry::setValue(std::string _value)
{
	joboption.value = _value;
	inp->value(_value.c_str());
	// Also update menu or slider if necessary
	if (menu != NULL)
	{
		const Fl_Menu_Item *p = menu->find_item(inp->value());
		if ( p )
			menu->picked(p);
		else
			REPORT_ERROR("Error readValue: Menu item not found:" + std::string(inp->value()) + " for joboption label= " + joboption.label);
	}
	if (slider != NULL)
	{
		slider->value(fltkTextToFloat(inp->value()));
	}
}


void GuiEntry::deactivate(bool do_deactivate)
{
	if (do_deactivate)
	{
		if (inp)
			inp->deactivate();
		if (help)
			help->deactivate();
		if (browse)
			browse->deactivate();
		if (menu)
			menu->deactivate();
		if (slider)
			slider->deactivate();
	}
	else
	{
		if (inp)
			inp->activate();
		if (help)
			help->activate();
		if (browse)
			browse->activate();
		if (menu)
			menu->activate();
		if (slider)
			slider->activate();
	}

}

// Help button call-back functions
void GuiEntry::cb_help(Fl_Widget* o, void* v)
{

    GuiEntry* T=(GuiEntry*)v;
    T->cb_help_i();
}

void GuiEntry::cb_help_i()
{

    ShowHelpText *help = new ShowHelpText(joboption.helptext.c_str());

}

void GuiEntry::cb_browse(Fl_Widget* o, void* v)
{

    GuiEntry* T=(GuiEntry*)v;
    T->cb_browse_i();
}


void GuiEntry::cb_browse_i()
{

    Fl::scheme("gtk+");
    Fl_File_Chooser * G_chooser = new Fl_File_Chooser("", joboption.pattern.c_str(), Fl_File_Chooser::SINGLE, "");

    if (joboption.directory=="CURRENT_ODIR")
    	G_chooser->directory(current_browse_directory.c_str());
    else
    	G_chooser->directory(joboption.directory.c_str());
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

    FileName fn_pre, fn_jobnr, fn_post, fn_out;
    decomposePipelineSymlinkName(relname, fn_pre, fn_jobnr, fn_post);
    fn_out = fn_pre + fn_jobnr + fn_post;

    inp->value(fn_out.c_str());
}


void GuiEntry::cb_browse_node(Fl_Widget* o, void* v) {

    GuiEntry* T=(GuiEntry*)v;
    T->cb_browse_node_i();
}


void GuiEntry::cb_browse_node_i() {

    Fl::scheme("gtk+");
    Fl_File_Chooser * G_chooser = new Fl_File_Chooser("", joboption.pattern.c_str(), Fl_File_Chooser::SINGLE, "");

    std::string fn_dir = (do_oldstyle) ? "." : ".Nodes/" + integerToString(joboption.node_type);
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
    if (do_oldstyle)
    {
    	inp->value(relname);
    }
    else
    {
		std::string replace = std::string(relname);
		std::string replace2 = (std::string::npos == replace.find(fn_dir.c_str())) ? replace : replace.substr(fn_dir.length()+1, replace.length());
		char relname2[FL_PATH_MAX];
		strcpy(relname2, replace2.c_str());

		FileName fn_pre, fn_jobnr, fn_post, fn_out;
	    decomposePipelineSymlinkName(replace2, fn_pre, fn_jobnr, fn_post);
	    fn_out = fn_pre + fn_jobnr + fn_post;

	    inp->value(fn_out.c_str());
    }

}

void GuiEntry::cb_menu(Fl_Widget* o, void* v) {

    GuiEntry* T=(GuiEntry*)v;
    T->cb_menu_i();
}


void GuiEntry::cb_menu_i()
{
	const Fl_Menu_Item* m = menu->mvalue();
	// Set my own value
	inp->value(m->label());
	// In case this was a boolean that deactivates a group, do so:
	if (my_deactivate_group != NULL)
	{
		if (strcmp(inp->value(), "No") == 0)
			my_deactivate_group->deactivate();
		else
			my_deactivate_group->activate();
	}

}

void GuiEntry::cb_slider(Fl_Widget* o, void* v) {

    GuiEntry* T=(GuiEntry*)v;
    T->cb_slider_i();
}


void GuiEntry::cb_slider_i() {

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

void GuiEntry::cb_input(Fl_Widget* o, void* v) {

    GuiEntry* T=(GuiEntry*)v;
    T->cb_input_i();
}


void GuiEntry::cb_input_i() {

    static int recurse = 0;
    if ( recurse ) {
        return;
    } else {
        recurse = 1;
        slider->value(fltkTextToFloat(inp->value()));         // pass input's value to slider
        recurse = 0;
    }
}



