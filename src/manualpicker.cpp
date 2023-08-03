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

#include "src/manualpicker.h"

std::vector<int> imics;
std::vector<FileName> global_fn_mics;
std::vector<FileName> global_fn_picks;
std::vector<FileName> global_fn_ctfs;
std::vector<bool> selected;
std::vector<int> number_picked;
std::vector<Fl_Button*> viewmic_buttons;
std::vector<Fl_Button*> viewctf_buttons;
std::vector<Fl_Text_Display*> text_displays;
std::vector<Fl_Text_Display*> count_displays;
std::vector<Fl_Text_Display*> defocus_displays;
std::vector<Fl_Check_Button*> check_buttons;
int first_pick_viewed, last_pick_viewed;
int last_ctf_viewed;

bool   global_has_ctf;
bool   global_pick_startend;
RFLOAT global_angpix;
RFLOAT global_coord_scale;
RFLOAT global_lowpass;
RFLOAT global_highpass;
bool global_do_topaz_denoise;
FileName global_topaz_exe;
RFLOAT global_particle_diameter;
RFLOAT global_sigma_contrast;
RFLOAT global_black_val;
RFLOAT global_white_val;
RFLOAT global_micscale;
RFLOAT global_ctfscale;
RFLOAT global_ctfsigma;
RFLOAT global_minimum_fom;
RFLOAT global_blue_value;
RFLOAT global_red_value;
int    global_total_count;
int global_nr_simultaneous;
FileName global_fn_odir;
FileName global_pickname;
FileName global_fn_color;
FileName global_color_label;
bool global_do_color;

void cb_viewmic(Fl_Widget* w, void* data)
{
	// Get my own number back
	int *iptr = (int*)data;
	int imic = *iptr;

	const bool with_control = (Fl::event_ctrl() != 0);
	int nr_simultaneous = (with_control) ? global_nr_simultaneous : 1;

	// Update the count of the last one we picked...
	for (int mymic = first_pick_viewed; mymic <= last_pick_viewed; mymic++)
	{
		if (mymic >= 0 && mymic < count_displays.size())
		{
			MetaDataTable MDcoord;
			FileName fn_coord = global_fn_picks[mymic];
			int my_nr_picked;
			if (exists(fn_coord))
			{
				MDcoord.read(fn_coord);
				if (fabs(global_minimum_fom + 9999.) > 1e-6)
				{
					if (MDcoord.containsLabel(EMDL_PARTICLE_AUTOPICK_FOM))
					{
						my_nr_picked = 0;
						FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDcoord)
						{
							RFLOAT fom;
							MDcoord.getValue(EMDL_PARTICLE_AUTOPICK_FOM, fom);
							if (fom > global_minimum_fom) my_nr_picked++;
						}
					}
					else
					{
						my_nr_picked = MDcoord.numberOfObjects();
					}

				}
				else
				{
					my_nr_picked = MDcoord.numberOfObjects();
				}
			}
			else
			{
				my_nr_picked = 0;
			}
			Fl_Text_Buffer *textbuff2 = new Fl_Text_Buffer();
			textbuff2->text(floatToString(my_nr_picked).c_str());
			count_displays[mymic]->buffer(textbuff2);
			count_displays[mymic]->redraw();
			// Also reset the color of the button to light
			viewmic_buttons[mymic]->color(GUI_BUTTON_COLOR, GUI_BUTTON_COLOR);
			viewmic_buttons[mymic]->redraw();
		}
	}

	// Launch the picking window
	first_pick_viewed = imic;
	last_pick_viewed = XMIPP_MIN(global_fn_mics.size() - 1, imic + nr_simultaneous - 1);
	for (int mymic = first_pick_viewed; mymic <= last_pick_viewed; mymic++)
	{
		FileName fn_coord = global_fn_picks[mymic];

		int rad = ROUND(global_particle_diameter/(2. * global_angpix));
		std::string command;
		command =  "relion_display --pick  --i " + global_fn_mics[mymic];
		command += "  --coords " + fn_coord;
		command += " --scale " + floatToString(global_micscale);
		command += " --coord_scale " + floatToString(global_coord_scale);
		command += " --black "  + floatToString(global_black_val);
		command += " --white "  + floatToString(global_white_val);
		command += " --sigma_contrast "  + floatToString(global_sigma_contrast);
		command += " --particle_radius " + floatToString(rad);
		if (global_do_topaz_denoise)
		{
			command += " --topaz_denoise --topaz_exe " + global_topaz_exe;
		}
		command += " --lowpass " + floatToString(global_lowpass);
		command += " --highpass " + floatToString(global_highpass);
		command += " --angpix " + floatToString(global_angpix);
		if (global_pick_startend)
			command += " --pick_start_end ";

		if (fabs(global_minimum_fom + 9999.) > 1e-6)
		{
			command += " --minimum_pick_fom " + floatToString(global_minimum_fom);
		}
		if (global_color_label != "")
		{
			command += " --color_label " + global_color_label;
			command += " --blue " + floatToString(global_blue_value);
			command += " --red " + floatToString(global_red_value);
			if (global_fn_color != "")
				command += " --color_star " + global_fn_color;
		}

		command += " &";
		std::cerr << " command= " << command << std::endl;
		int res = system(command.c_str());
	}

	for (int i = 0; i < viewmic_buttons.size(); i++)
	{
		if (i >= first_pick_viewed && i <= last_pick_viewed)
		{
			viewmic_buttons[i]->color(GUI_BUTTON_DARK_COLOR, GUI_BUTTON_DARK_COLOR);
		}
		else
		{
			viewmic_buttons[i]->color(GUI_BUTTON_COLOR, GUI_BUTTON_COLOR);
		}
		viewmic_buttons[i]->redraw();

	}
}

void cb_viewctf(Fl_Widget* w, void* data)
{
	// Get my own number back
	int *iptr = (int*)data;
	int imic = *iptr;

	std::string command;
	command =  "relion_display --i " + global_fn_ctfs[imic];
	command += " --scale " + floatToString(global_ctfscale);
	command += " --sigma_contrast " + floatToString(global_ctfsigma);
	command += " &";
	int res = system(command.c_str());

	last_ctf_viewed = imic;
	for (int i = 0; i < viewctf_buttons.size(); i++)
	{
		if (i == last_ctf_viewed)
		{
			viewctf_buttons[i]->color(GUI_BUTTON_DARK_COLOR, GUI_BUTTON_DARK_COLOR);
		}
		else
		{
			viewctf_buttons[i]->color(GUI_BUTTON_COLOR, GUI_BUTTON_COLOR);
		}
	}
}

void cb_selectmic(Fl_Widget* w, void* data)
{
	// Get my own number back
	int *iptr = (int*)data;
	int imic = *iptr;

	Fl_Text_Buffer *textbuff2 = new Fl_Text_Buffer();
	selected[imic] = !selected[imic];
	if (selected[imic])
	{
		text_displays[imic]->color(GUI_INPUT_COLOR, GUI_INPUT_COLOR);
		text_displays[imic]->activate();
		viewmic_buttons[imic]->activate();
		count_displays[imic]->color(GUI_INPUT_COLOR, GUI_INPUT_COLOR);
		textbuff2->text(floatToString(number_picked[imic]).c_str());
		count_displays[imic]->buffer(textbuff2);
		count_displays[imic]->activate();
		if (global_has_ctf)
		{
			viewctf_buttons[imic]->activate();
			defocus_displays[imic]->color(GUI_INPUT_COLOR, GUI_INPUT_COLOR);
			defocus_displays[imic]->activate();
		}
	}
	else
	{
		text_displays[imic]->color(GUI_BACKGROUND_COLOR, GUI_BACKGROUND_COLOR);
		text_displays[imic]->deactivate();
		viewmic_buttons[imic]->deactivate();
		count_displays[imic]->color(GUI_BACKGROUND_COLOR, GUI_BACKGROUND_COLOR);
		textbuff2->text("");
		count_displays[imic]->buffer(textbuff2);
		count_displays[imic]->deactivate();
		if (global_has_ctf)
		{
			viewctf_buttons[imic]->deactivate();
			defocus_displays[imic]->color(GUI_BACKGROUND_COLOR, GUI_BACKGROUND_COLOR);
			defocus_displays[imic]->deactivate();
		}
	}
}

int manualpickerGuiWindow::fill()
{
	color(GUI_BACKGROUND_COLOR);


	Fl_Menu_Bar *menubar = new Fl_Menu_Bar(0, 0, w(), 25);
	if (do_allow_save)
	{
		menubar->add("File/Save selection",  FL_ALT+'s', cb_menubar_save, this);
		menubar->add("File/Invert selection",  FL_ALT+'i', cb_menubar_invert_selection, this);
        menubar->add("File/Select all",  FL_ALT+'a', cb_menubar_select_all, this);
	}
	menubar->add("File/Recount picked particles",  FL_ALT+'c', cb_menubar_recount, this);
	menubar->add("File/Set FOM threshold",  FL_ALT+'c', cb_menubar_setFOM, this);
	menubar->add("File/Quit", FL_ALT+'q', cb_menubar_quit, this);
	int current_y = 25;

	// Scroll bars
	Fl_Scroll scroll(0, current_y, w(), h()-current_y);
	scroll.type(Fl_Scroll::VERTICAL);

	selected.clear();
	number_picked.clear();

	global_has_ctf = MDin.containsLabel(EMDL_CTF_IMAGE);

	FileName fn_mic, fn_pick, fn_ctf;
	int ystep = 35;

	imics.clear();
	for (int ii =0; ii < MDin.numberOfObjects(); ii++)
	{
		imics.push_back(ii);
	}

	int imic =0;
	global_fn_mics.clear();
	global_fn_picks.clear();
	global_fn_ctfs.clear();
	text_displays.clear();
	viewmic_buttons.clear();
	viewctf_buttons.clear();
	number_picked.clear();
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
	{
		MDin.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
		// Display the name of the micrograph
		global_fn_mics.push_back(fn_mic);

		if (MDin.containsLabel(EMDL_MICROGRAPH_COORDINATES))
		{
			// relion 3.2+
			MDin.getValue(EMDL_MICROGRAPH_COORDINATES, fn_pick);

		}
		else
		{
			//relion 3.1-
			FileName fn_pre, fn_jobnr, fn_post;
			decomposePipelineSymlinkName(fn_mic, fn_pre, fn_jobnr, fn_post);
			fn_pick = global_fn_odir + fn_post.withoutExtension() + "_" + global_pickname + ".star";
		}
		global_fn_picks.push_back(fn_pick);

		Fl_Check_Button *mycheck = new Fl_Check_Button(4, current_y, ystep-8, ystep-8, "");
		mycheck->callback(cb_selectmic, &(imics[imic]));
		mycheck->value(1);
		if (!do_allow_save)
			mycheck->deactivate();
		selected.push_back(true);
		check_buttons.push_back(mycheck);

		Fl_Text_Buffer *textbuff = new Fl_Text_Buffer();
		textbuff->text(fn_mic.c_str());
		int ystep2 = (fn_mic.length() > MWCOL1/12) ? ystep - 5 : ystep - 10;
		Fl_Text_Display* mydisp = new Fl_Text_Display(MXCOL0, current_y, MWCOL1, ystep2);
		mydisp->scrollbar_width(5);
		mydisp->buffer(textbuff);
		mydisp->scroll(0, 9999);
		mydisp->color(GUI_INPUT_COLOR, GUI_INPUT_COLOR);
		text_displays.push_back(mydisp);

		// Button to display the micrographimage
		Fl_Button *myviewmic = new Fl_Button(MXCOL1, current_y, MWCOL2, ystep-5, "pick");
		myviewmic->color(GUI_BUTTON_COLOR);
		myviewmic->callback(cb_viewmic, &(imics[imic]));
		viewmic_buttons.push_back(myviewmic);

		// Count how many particles have been picked
		Fl_Text_Buffer *textbuff2 = new Fl_Text_Buffer();
		textbuff2->text("");
		Fl_Text_Display* mycount = new Fl_Text_Display(MXCOL2, current_y, MWCOL3, ystep-5);
		mycount->color(GUI_INPUT_COLOR, GUI_INPUT_COLOR);
		mycount->buffer(textbuff2);
		count_displays.push_back(mycount);
		number_picked.push_back(10);

		// Button to display the CTF image
		if (global_has_ctf)
		{
			MDin.getValue(EMDL_CTF_IMAGE, fn_ctf);
			global_fn_ctfs.push_back(fn_ctf);
			// Button to display the CTF image
			Fl_Button *myviewctf = new Fl_Button(MXCOL3, current_y, MWCOL4, ystep-5, "CTF");
			myviewctf->color(GUI_BUTTON_COLOR);
			myviewctf->callback(cb_viewctf, &(imics[imic]));
			viewctf_buttons.push_back(myviewctf);

			Fl_Text_Buffer *textbuffDF = new Fl_Text_Buffer();
			RFLOAT defocus;
			MDin.getValue(EMDL_CTF_DEFOCUSU, defocus);

			std::ostringstream os;
			os << defocus;
			std::string str = os.str();
			textbuffDF->text(str.c_str());

			Fl_Text_Display* myDF = new Fl_Text_Display(MXCOL4, current_y, MWCOL4, ystep-5);
			myDF->color(GUI_INPUT_COLOR, GUI_INPUT_COLOR);
			myDF->buffer(textbuffDF);
			defocus_displays.push_back(myDF);
		}

		imic++;
		current_y += ystep;
	}


	// See if the output STAR file already exists, if so apply that selection
	if (do_allow_save)
		readOutputStarfile();

	if (do_fast_save)
		cb_menubar_save_i();

	// Also count the number of particles that were already picked
	cb_menubar_recount_i();

	resizable(*this);
	show();
	return Fl::run();
}

void manualpickerGuiWindow::readOutputStarfile()
{
	FileName fn_select = global_fn_odir + "local_selection.star";
    if (exists(fn_select))
	{
        for (int imic = 0; imic < selected.size(); imic++)
			selected[imic] = true;

        MetaDataTable MDselect;
        MDselect.read(fn_select);

        // Set all unselected micrographs to false
        FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDselect)
        {

            FileName fn_mic;
            int is_selected;
            MDselect.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
            MDselect.getValue(EMDL_SELECTED, is_selected);
            if (is_selected == 0)
            {
                // find fn_mic
                for (int imic = 0; imic < selected.size(); imic++)
                {
                    if (global_fn_mics[imic] == fn_mic)
                    {
                        selected[imic] = false;
                        break;
                    }
                }
            }

        }

        // Apply the selection to the buttons
		for (int imic = 0; imic < selected.size(); imic++)
		{
			if (selected[imic])
			{
				check_buttons[imic]->value(1);
				text_displays[imic]->color(GUI_INPUT_COLOR, GUI_INPUT_COLOR);
				text_displays[imic]->activate();
				viewmic_buttons[imic]->activate();
				count_displays[imic]->color(GUI_INPUT_COLOR, GUI_INPUT_COLOR);
				count_displays[imic]->activate();
				if (global_has_ctf)
					viewctf_buttons[imic]->activate();
			}
			else
			{
				check_buttons[imic]->value(0);
				text_displays[imic]->color(GUI_BACKGROUND_COLOR, GUI_BACKGROUND_COLOR);
				text_displays[imic]->deactivate();
				viewmic_buttons[imic]->deactivate();
				count_displays[imic]->color(GUI_BACKGROUND_COLOR, GUI_BACKGROUND_COLOR);
				count_displays[imic]->deactivate();
				if (global_has_ctf)
					viewctf_buttons[imic]->deactivate();
			}
		}
	}
}

void manualpickerGuiWindow::writeOutputStarfiles(bool verb)
{
	if (!do_allow_save) return;

	MDcoords.clear();
	MetaDataTable MDmics, MDselect;
	int c = 0;
	for (int imic = 0; imic < selected.size(); imic++)
	{
		if (selected[imic])
		{
			MDmics.addObject(MDin.getObject(imic));
			FileName fn_coord = global_fn_picks[imic];
			if (exists(fn_coord))
			{
				MDcoords.addObject();
				MDcoords.setValue(EMDL_MICROGRAPH_NAME, global_fn_mics[imic]);
				MDcoords.setValue(EMDL_MICROGRAPH_COORDINATES, fn_coord);
				c++;
			}
        }

        MDselect.addObject();
        MDselect.setValue(EMDL_MICROGRAPH_NAME, global_fn_mics[imic]);
        int sel = (selected[imic]) ? 1 : 0;
        MDselect.setValue(EMDL_SELECTED, sel);

    }

	if (obsModel.opticsMdt.numberOfObjects() > 0)
	{
		obsModel.save(MDmics, fn_sel, "micrographs");
	}
	else
	{
		MDmics.write(fn_sel);
	}
	if (verb) std::cout << " Saved list of selected micrographs in: " << fn_sel << std::endl;

    // Save selection star file for manualpicker continue jobs only
    FileName fn_select = global_fn_odir + "local_selection.star";
    MDselect.write(fn_select);

	FileName fn_coords = global_fn_odir + global_pickname + ".star";
	MDcoords.setName("coordinate_files");
	MDcoords.write(fn_coords);
	if (verb) std::cout << " Saved list with " << c << " coordinate files in: " << fn_coords << std::endl;

}
void manualpickerGuiWindow::cb_menubar_save(Fl_Widget* w, void* v)
{
	manualpickerGuiWindow* T=(manualpickerGuiWindow*)v;
	T->cb_menubar_save_i();
}

void manualpickerGuiWindow::cb_menubar_save_i()
{
	writeOutputStarfiles();
	RELION_EXIT_SUCCESS;
}

void manualpickerGuiWindow::cb_menubar_select_all(Fl_Widget* w, void* v)
{
    manualpickerGuiWindow* T=(manualpickerGuiWindow*)v;
    T->cb_menubar_select_all_i();
}

void manualpickerGuiWindow::cb_menubar_select_all_i()
{
	for (int imic = 0; imic < selected.size(); imic++)
    {
        selected[imic] = true;
        check_buttons[imic]->value(1);
        text_displays[imic]->color(GUI_INPUT_COLOR, GUI_INPUT_COLOR);
        text_displays[imic]->activate();
        viewmic_buttons[imic]->activate();
        count_displays[imic]->color(GUI_INPUT_COLOR, GUI_INPUT_COLOR);
        count_displays[imic]->activate();
        if (global_has_ctf)
            viewctf_buttons[imic]->activate();
    }

}

void manualpickerGuiWindow::cb_menubar_invert_selection(Fl_Widget* w, void* v)
{
	manualpickerGuiWindow* T=(manualpickerGuiWindow*)v;
	T->cb_menubar_invert_selection_i();
}

void manualpickerGuiWindow::cb_menubar_invert_selection_i()
{
	for (int imic = 0; imic < selected.size(); imic++)
	{
		selected[imic] = !selected[imic];
		if (selected[imic])
		{
			check_buttons[imic]->value(1);
			text_displays[imic]->color(GUI_INPUT_COLOR, GUI_INPUT_COLOR);
			text_displays[imic]->activate();
			viewmic_buttons[imic]->activate();
			count_displays[imic]->color(GUI_INPUT_COLOR, GUI_INPUT_COLOR);
			count_displays[imic]->activate();
			if (global_has_ctf)
				viewctf_buttons[imic]->activate();
		}
		else
		{
			check_buttons[imic]->value(0);
			text_displays[imic]->color(GUI_BACKGROUND_COLOR, GUI_BACKGROUND_COLOR);
			text_displays[imic]->deactivate();
			viewmic_buttons[imic]->deactivate();
			count_displays[imic]->color(GUI_BACKGROUND_COLOR, GUI_BACKGROUND_COLOR);
			count_displays[imic]->deactivate();
			if (global_has_ctf)
				viewctf_buttons[imic]->deactivate();
		}
	}
}

void manualpickerGuiWindow::cb_menubar_quit(Fl_Widget* w, void* v)
{
	manualpickerGuiWindow* T=(manualpickerGuiWindow*)v;
	T->cb_menubar_quit_i();
}

void manualpickerGuiWindow::cb_menubar_quit_i()
{
	cb_menubar_save_i();
	exit(0);
}

void manualpickerGuiWindow::cb_menubar_recount(Fl_Widget* w, void* v)
{
	manualpickerGuiWindow* T=(manualpickerGuiWindow*)v;
	T->cb_menubar_recount_i();
}
void manualpickerGuiWindow::cb_menubar_recount_i()
{

	global_total_count = 0;
	int nr_sel_mic = 0;
	for (int imic = 0; imic < global_fn_mics.size(); imic++)
	{
		FileName fn_coord = global_fn_picks[imic];
		MetaDataTable MDcoord;

		int my_nr_picked;
		if (exists(fn_coord))
		{
			MDcoord.read(fn_coord);

			if (fabs(global_minimum_fom + 9999.) > 1e-6)
			{
				if (MDcoord.containsLabel(EMDL_PARTICLE_AUTOPICK_FOM))
				{
					my_nr_picked = 0;
					FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDcoord)
					{
						RFLOAT fom;
						MDcoord.getValue(EMDL_PARTICLE_AUTOPICK_FOM, fom);
						if (fom > global_minimum_fom) my_nr_picked++;
					}
				}
				else
				{
					my_nr_picked = MDcoord.numberOfObjects();
				}

			}
			else
			{
				my_nr_picked = MDcoord.numberOfObjects();
			}
		}
		else
		{
			my_nr_picked = 0;
		}

		Fl_Text_Buffer *textbuff2 = new Fl_Text_Buffer();
		if (selected[imic])
		{
			global_total_count += my_nr_picked;
			textbuff2->text(floatToString(my_nr_picked).c_str());
			count_displays[imic]->buffer(textbuff2);
			count_displays[imic]->redraw();
			nr_sel_mic++;
		}
		else
		{
			textbuff2->text("");
			count_displays[imic]->buffer(textbuff2);
		}
		number_picked[imic] = my_nr_picked;
	}
	std::cout << " Total number of picked particles: " << global_total_count << " from " << nr_sel_mic << " selected micrographs." << std::endl;
	writeOutputStarfiles();
}

void manualpickerGuiWindow::cb_menubar_setFOM(Fl_Widget* w, void* v)
{
	manualpickerGuiWindow* T=(manualpickerGuiWindow*)v;
	T->cb_menubar_setFOM_i();
	T->cb_menubar_recount_i();
}

void manualpickerGuiWindow::cb_menubar_setFOM_i()
{
	const char *pfom;
	std::string currentval = floatToString(global_minimum_fom);
	pfom =  fl_input("Minimum rlnAutopickFigureOfMerit to display: ", currentval.c_str());
	if (pfom == NULL)
		return;
	std::string newval(pfom);
	global_minimum_fom = textToFloat(newval);

}


void ManualPicker::read(int argc, char **argv)
{
	parser.setCommandLine(argc, argv);

	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "Micrograph STAR file OR filenames from which to pick particles, e.g. \"Micrographs/*.mrc\"");
	global_fn_odir = parser.getOption("--odir", "Output directory for coordinate files (default is to store next to micrographs)", "ManualPick/");
	fn_sel = parser.getOption("--selection", "STAR file with selected micrographs", "micrographs_selected.star");
	global_pickname = parser.getOption("--pickname", "Rootname for the picked coordinate files", "manualpick");
	global_angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms", "-1."));
	global_coord_scale = textToFloat(parser.getOption("--coord_scale", "Scale coordinates before display", "1.0"));
	global_particle_diameter = textToFloat(parser.getOption("--particle_diameter", "Diameter of the circles that will be drawn around each picked particle (in Angstroms)"));
	global_pick_startend = parser.checkOption("--pick_start_end", "Pick start-end coordinates of helices");
	do_allow_save = parser.checkOption("--allow_save", "Allow saving of the selected micrographs");
	do_fast_save = parser.checkOption("--fast_save", "Save a default selection of all micrographs immediately");
	global_nr_simultaneous = textToInteger(parser.getOption("--open_simultaneous", "Open this many of the next micrographs simultaneously when pressing CTRL and a Pick button", "10"));

	int mic_section = parser.addSection("Displaying options");
	global_micscale = textToFloat(parser.getOption("--scale", "Relative scale for the micrograph display", "1"));
	global_black_val = textToFloat(parser.getOption("--black", "Pixel value for black (default is auto-contrast)", "0"));
	global_white_val = textToFloat(parser.getOption("--white", "Pixel value for white (default is auto-contrast)", "0"));
	global_sigma_contrast  = textToFloat(parser.getOption("--sigma_contrast", "Set white and black pixel values this many times the image stddev from the mean (default is auto-contrast)", "0"));
	global_lowpass = textToFloat(parser.getOption("--lowpass", "Lowpass filter in Angstroms for the micrograph (0 for no filtering)","0"));
	global_highpass = textToFloat(parser.getOption("--highpass", "Highpass filter in Angstroms for the micrograph (0 for no filtering)","0"));
	global_do_topaz_denoise = parser.checkOption("--topaz_denoise", "Or instead of filtering, use Topaz denoising before picking (on GPU 0)");
	global_topaz_exe = parser.getOption("--topaz_exe", "Name of topaz executable", "topaz");
	global_ctfscale = textToFloat(parser.getOption("--ctf_scale", "Relative scale for the CTF-image display", "1"));
	global_ctfsigma = textToFloat(parser.getOption("--ctf_sigma_contrast", "Sigma-contrast for the CTF-image display", "3"));
	global_minimum_fom = textToFloat(parser.getOption("--minimum_pick_fom", "Minimum value for rlnAutopickFigureOfMerit to display picks", "-9999."));
	// coloring
	global_fn_color = parser.getOption("--color_star", "STAR file with a column for red-blue coloring (a subset of) the particles", "");
	global_color_label = parser.getOption("--color_label", "MetaDataLabel to color particles on (e.g. rlnParticleSelectZScore)", "");
	global_blue_value = textToFloat(parser.getOption("--blue", "Value of the blue color", "1."));
	global_red_value = textToFloat(parser.getOption("--red", "Value of the red color", "0."));

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}

void ManualPicker::usage()
{
	parser.writeUsage(std::cout);
}

void ManualPicker::initialise()
{
	if (fn_in.isStarFile())
	{
		// First try 2-column list of coordinate files as in relion-3.2+
		MDin.read(fn_in, "coordinate_files");
		if (MDin.numberOfObjects() > 0)
		{
			if (global_angpix < 0.)
			{
				std::cerr << " WARNING: no --angpix provided and no information about pixel size in input STAR file. Setting angpix to 1..." << std::endl;
				global_angpix = 1.;
			}

		}
		else
		{

			// Normal micrographs.star file (with optics table etc)

			ObservationModel::loadSafely(fn_in, obsModel, MDin, "micrographs");
			if (obsModel.opticsMdt.containsLabel(EMDL_MICROGRAPH_PIXEL_SIZE))
			{
				obsModel.opticsMdt.getValue(EMDL_MICROGRAPH_PIXEL_SIZE, global_angpix, 0);
				std::cout << " Setting angpix to " << global_angpix << " based on the input STAR file... " << std::endl;
			}
			else
			{
				if (global_angpix < 0.)
				{
					REPORT_ERROR("ERROR: the input STAR file does not contain the micrograph pixel size, and it is not given through --angpix.");
				}
				std::cout << " Setting angpix to " << global_angpix << " based on command-line input... " << std::endl;
				FOR_ALL_OBJECTS_IN_METADATA_TABLE(obsModel.opticsMdt)
				{
					obsModel.opticsMdt.setValue(EMDL_MICROGRAPH_PIXEL_SIZE, global_angpix);
				}
			}
		}
	}
	else
	{
		std::vector<FileName> glob_fn_mics;
		fn_in.globFiles(glob_fn_mics);
		for (int imic = 0; imic < glob_fn_mics.size(); imic++)
		{
			MDin.addObject();
			MDin.setValue(EMDL_MICROGRAPH_NAME, glob_fn_mics[imic]);
		}

		if (global_angpix < 0.)
		{
			std::cerr << " WARNING: no --angpix provided and no information about pixel size in input STAR file. Setting angpix to 1..." << std::endl;
			global_angpix = 1.;
		}
	}

	// If we down-scale the micrograph: always low-pass filter to get better displays
	if (global_micscale < 1.)
	{
		RFLOAT new_nyquist = global_angpix * 2. / global_micscale;
		if (new_nyquist > global_lowpass)
			global_lowpass = new_nyquist;
		std::cout << " Set low-pass filter to " << global_lowpass << " due to downscaling of " << global_micscale << std::endl;
	}

	std::cerr << " NOTE: in order to write the new list of coordinate STAR files, you need to re-count the particles or quite this program through the File menu. Do NOT kill the program using the operating system's window manager!" << std::endl;


}

void ManualPicker::run()
{
	Fl::scheme("gtk+");

	manualpickerGuiWindow win(TOTALWIDTH, TOTALHEIGHT, "RELION manual-picking GUI");

	// Transfer all parameters to the gui
	win.MDin = MDin;
	win.MDcoords = MDcoords;
	win.obsModel = obsModel;
	win.fn_sel = fn_sel;
	win.do_allow_save = do_allow_save;
	win.do_fast_save = do_fast_save;
	win.fill();
}
