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

#ifndef DISPLAYER_H_
#define DISPLAYER_H_

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

#include "src/image.h"
#include "src/metadata_label.h"
#include "src/metadata_table.h"
#include <src/jaz/single_particle/obs_model.h>
#include <src/matrix2d.h>
#include <src/fftw.h>
#include <src/time.h>
#include <src/args.h>


#define GUI_BACKGROUND_COLOR (fl_rgb_color(240,240,240))
#define GUI_INPUT_COLOR (fl_rgb_color(255,255,230))
#define GUI_RUNBUTTON_COLOR (fl_rgb_color(160, 30, 60))



#define SELECTED 1
#define NOTSELECTED 0
#define MULTIVIEW_WINDOW_WIDTH  720
#define MULTIVIEW_WINDOW_HEIGHT 486

#define BOX_OFFSET 4

#define MULTIVIEWER 0
#define SINGLEVIEWER 1
#define BLACK 0
#define WHITE 1

#define DEFAULTPDFVIEWER "evince"

static bool has_dragged;
static int predrag_xc;
static int predrag_yc;
static bool has_shift;
static int preshift_ipos;

static int current_selection_type;
static int colour_scheme;

class DisplayBox : public Fl_Box
{
protected:
	// Draw the actual box on the screen (this function is used by redraw())
	void draw();

public:
	int xsize_data;
	int ysize_data;
	int xoff;
	int yoff;

	// The box's selection status
	int selected;

	// The box's original position in the input MetaDataTable
	int ipos;

	// The metadata fir this image
	MetaDataTable MDimg;

	// The actual image data array
	unsigned char *img_data;
	std::string img_label;

	// For getting back close the original image values from the uchar ones...
	RFLOAT minval;
	RFLOAT maxval;
	RFLOAT scale;

	// Constructor with an image and its metadata
	DisplayBox(int X, int Y, int W, int H, const char *L=0) : Fl_Box(X,Y,W,H,L) { img_data = NULL; img_label = ""; MDimg.clear(); }

	void setData(MultidimArray<RFLOAT> &img, MetaDataContainer *MDCin, int ipos, RFLOAT minval, RFLOAT maxval,
	             RFLOAT _scale, bool do_relion_scale = false);

	// Destructor
	~DisplayBox()
	{
		MDimg.clear();
		if (img_data)
			delete [] img_data;
	};

	// Change selected status, redraw and return new status
	int toggleSelect(int set_selected);
	// Set a specific value in selected for this box
	void setSelect(int value);
	// Select, redraw and return new selected status
	int select();
	// unSelect, redraw and return new selected status
	int unSelect();
};

// This class only puts scrollbars around the resizable canvas
class basisViewerWindow : public Fl_Window
{
public:
	// Constructor with w x h size of the window and a title
	basisViewerWindow(int W, int H, const char* title=0): Fl_Window(W, H, title)
	{
		current_selection_type = 1;
	}

	int fillCanvas(int viewer_type, MetaDataTable &MDin, ObservationModel *obsModel, EMDLabel display_label, EMDLabel text_label, bool _do_read_whole_stacks, bool _do_apply_orient,
	               RFLOAT _minval, RFLOAT _maxval, RFLOAT _sigma_contrast,
	               RFLOAT _scale, RFLOAT _ori_scale, int _ncol, long int max_nr_images = -1, RFLOAT lowpass = -1.0 , RFLOAT highpass = -1.0,
	               bool do_class = false, MetaDataTable *MDdata = NULL,
	               int _nr_regroup = -1, bool do_recenter = false, bool _is_data = false, MetaDataTable *MDgroups = NULL,
	               bool do_allow_save = false, FileName fn_selected_imgs="", FileName fn_selected_parts="", int max_nr_parts_per_class = -1);
	int fillSingleViewerCanvas(MultidimArray<RFLOAT> image, RFLOAT _minval, RFLOAT _maxval, RFLOAT _sigma_contrast, RFLOAT _scale);
	int fillPickerViewerCanvas(MultidimArray<RFLOAT> image, RFLOAT _minval, RFLOAT _maxval, RFLOAT _sigma_contrast, RFLOAT _scale, RFLOAT _coord_scale,
	                           int _particle_radius, bool do_startend = false, FileName _fn_coords = "",
	                           FileName _fn_color = "", FileName _fn_mic= "", FileName _color_label = "", RFLOAT _color_blue_value = 0., RFLOAT _color_red_value = 1.,
							   RFLOAT _minimum_pick_fom = -999.);
};

class basisViewerCanvas : public Fl_Widget
{
protected:

	void draw();
	void saveImage(int ipos=0);

public:

	int ncol;
	int nrow;
	int xsize_box;
	int ysize_box;
	int xoff;
	int yoff;

	// To get positions in scrolled canvas...
	Fl_Scroll *scroll;

	// All the individual image display boxes
	std::vector<DisplayBox*> boxes;

	// Read stacks at once to speed up?
	bool do_read_whole_stacks;

	// Minimum value for rlnAutopickFigureOfMerit to display picks
	RFLOAT minimum_pick_fom;

	// Constructor with w x h size of the window and a title
	basisViewerCanvas(int X,int Y, int W, int H, const char* title=0) : Fl_Widget(X,Y,W, H, title) { }

	void SetScroll(Fl_Scroll *val) { scroll = val; }

	void fill(MetaDataTable &MDin, ObservationModel *obsModel, EMDLabel display_label, EMDLabel text_label, bool _do_apply_orient, RFLOAT _minval, RFLOAT _maxval,
	          RFLOAT _sigma_contrast, RFLOAT _scale, int _ncol, bool do_recenter = false, long int max_images = -1,
	          RFLOAT lowpass = -1.0, RFLOAT highpass = -1.0);
	void fill(MultidimArray<RFLOAT> &image, RFLOAT _minval, RFLOAT _maxval, RFLOAT _sigma_contrast, RFLOAT _scale = 1.);
	void setSelectionType();
	void setFOMThreshold();
};

class multiViewerCanvas : public basisViewerCanvas
{
protected:

	int handle(int ev);

public:

	// Flag to indicate whether this is a viewer for class averages from 2D/3D relion_refine classification runs
	bool do_class;

	// Allow saving of selection files?
	bool do_allow_save;

	// Filenames with the selected class averages and the particles from the selected classes
	FileName fn_selected_imgs, fn_selected_parts;

	// Maximum number of selected particles per class
	int max_nr_parts_per_class;

	// Flag to indicate whether this is a viewer for a data.star (to also allow regrouping)
	bool is_data;

	// Number of groups for regrouping the selected particles (for model.star)
	int nr_regroups;

	// Apply the positioning parameters in the particle metadata?
	bool do_apply_orient;

	// Save selected images recentered to their center-of-mass?
	bool do_recenter;

	// pointer to the MetaDataTable for the individually aligned particles when do_class (the data.star file)
	MetaDataTable *MDdata;

	// pointer to the MetaDataTable for the optics groups
	ObservationModel *obsModel;

	// pointer to the MetaDataTable for the groups when do_class and do_regroup (the data.star file)
	MetaDataTable *MDgroups;

	// Sjors 12mar18: read/write information-containing backup_selection for Liyi's project
	MetaDataTable MDbackup;

	// Scale for showing the original image
	RFLOAT ori_scale;

	// To know which original image to display
	EMDLabel display_label;

	// Label for text display
	EMDLabel text_label;

	// To know which contrast to apply to original image display
	RFLOAT sigma_contrast, minval, maxval;

	// Limit number of images to be shown
	long int multi_max_nr_images;

	// Name of the metadata table
	std::string metadata_table_name;

	// Constructor with w x h size of the window and a title
	multiViewerCanvas(int X,int Y, int W, int H, const char* title=0): basisViewerCanvas(X,Y,W, H, title)
	{
	}

private:

	// Functionalities for  popup menu
	void clearSelection();
	void invertSelection();
	void selectFromHereBelow(int ipos);
	void selectFromHereAbove(int ipos);
	void printMetaData(int ipos);
	void showAverage(bool selected, bool show_stddev=false);
	void showOriginalImage(int ipos);
	void showFourierAmplitudes(int ipos);
	void showFourierPhaseAngles(int ipos);
	void showHelicalLayerLineProfile(int ipos);
	void makeStarFileSelectedParticles(int save_selected, MetaDataTable &MDpart);
	void saveSelectedParticles(int save_selected);
	void showSelectedParticles(int save_selected);
	void saveTrainingSet();
	void saveSelected(int save_selected);
	void saveBackupSelection();
	// Allow re-loading of existing backup selection
public:
	void loadBackupSelection(bool do_ask = true);
};

// Generally accessible function
void regroupSelectedParticles(MetaDataTable &MDdata, MetaDataTable &MDgroups, int nr_regroups);

class popupSelectionTypeWindow : Fl_Window
{
	Fl_Choice * choice;
public:
	int result;

	// Constructor with w x h size of the window and a title
	popupSelectionTypeWindow(int W, int H, const char* title=0): Fl_Window(W, H, title){}

	int fill();

	static void cb_set(Fl_Widget* o, void* v)
	{
		popupSelectionTypeWindow* T=(popupSelectionTypeWindow*)v;
		T->cb_set_i();
	}
	inline void cb_set_i()
	{
		// void* to (small) int
		current_selection_type = static_cast<int>(reinterpret_cast<intptr_t>(choice->mvalue()->user_data()));
	}

	static void cb_close(Fl_Widget* o, void* v)
	{
		popupSelectionTypeWindow* T=(popupSelectionTypeWindow*)v;
		T->hide();
	}

};

class singleViewerCanvas : public basisViewerCanvas
{

protected:
	int handle(int ev);

public:

	// Constructor with w x h size of the window and a title
	singleViewerCanvas(int X, int Y, int W, int H, const char* title=0): basisViewerCanvas(X,Y,W, H, title) { }
private:

	// Functionalities for  popup menu
	void printMetaData();
	// void setContrast();

	// explain functionality of clicks
	void printHelp();
};

/*
class popupSetContrastWindow : Fl_Window
{
	Fl_Input *minval, *maxval, *scale, *sigma_contrast;

public:
	int result;

	// Constructor with w x h size of the window and a title
	popupSetContrastWindow(int W, int H, const char* title=0): Fl_Window(W, H, title){}

	int fill();

	static void cb_set(Fl_Widget* o, void* v)
	{
		popupSelectionTypeWindow* T=(popupSelectionTypeWindow*)v;
		T->cb_set_i();
	}
	inline void cb_set_i()
	{
		// Careful with setting the right value! Look at handle function of multiviewerCanvas
		current_minval =  textToFloat(minval->value());
		current_maxval =  textToFloat(minval->value());
		current_scale =  textToFloat(minval->value());
		current_sigma_contrast =  textToFloat(sigma_contrast->value());
	}

	static void cb_close(Fl_Widget* o, void* v)
	{
		popupSetContrastWindow* T=(popupSetContrastWindow*)v;
		T->hide();
	}

};
*/

class pickerViewerCanvas : public basisViewerCanvas
{
protected:
	int handle(int ev);
	void draw();

public:
	// MetaDataTable with all picked coordinates
	MetaDataTable MDcoords;

	int particle_radius;

	// Scale for rlnCoordinateX/Y
	RFLOAT coord_scale;

	// Filename of the picked coordinate files
	FileName fn_coords;

	// FileName of the STAR file that contains the color-based column
	FileName fn_color;

	// Label to base coloring on
	EMDLabel color_label;

	// Blue value for coloring
	RFLOAT smallest_color_value;

	// Red value for coloring
	RFLOAT biggest_color_value;

	// Red->Blue is true; blue->red is false
	bool do_blue_to_red;

	// Draw lines between start-end coordinates?
	bool do_startend;

	// Micrograph name (useful to search relevant particles in fn_color)
	FileName fn_mic;

	// Constructor with w x h size of the window and a title
	pickerViewerCanvas(int X, int Y, int W, int H, const char* title=0): basisViewerCanvas(X,Y,W, H, title) { }

	void loadCoordinates(bool ask_filename = false);

	// if a fn_zscore is given, then match the coordinates to the Zscores in the corresponding MDtable
	void findColorColumnForCoordinates();

private:

	// Functionalities for  popup menu
	void saveCoordinates(bool ask_filename = false);
	void clearCoordinates();
	void printHelp();
	void viewExtractedParticles();
};

// This class only puts scrollbars around the resizable canvas
class displayerGuiWindow : public Fl_Window
{
public:

	FileName fn_in, fn_data;

	// Some general settings for different types
	bool is_class;

	bool is_multi;

	bool is_star;

	// Allow regrouping from _data.star
	bool is_data;

	// Allow saving of the selection?
	bool do_allow_save;

	// remove duplicate-coord particles in fn_data STAR file?
	bool do_remove_duplicates;
	RFLOAT duplicate_threshold;	 // Angstrom radius

	// Number of regroups
	int nr_regroups;

	// Recenter images?
	bool do_recenter;

	// Pipeline control
	std::string pipeline_control;

	// Maximum number of images to show
	long int max_nr_images;

	// Display image in color
	int colour_scheme;

	// FileName for selected class average images and particles
	FileName fn_imgs, fn_parts;

	// Label option to display or to sort on
	std::vector<std::string> display_labels;
	std::vector<std::string> sort_labels;
	std::vector<std::string> colour_schemes;

	// Input for the display parameters
	Fl_Input *black_input, *white_input, *sigma_contrast_input, *scale_input, *lowpass_input, *highpass_input, *angpix_input;
	Fl_Input *col_input, *ori_scale_input, *max_nr_images_input, *max_parts_per_class_input;
	Fl_Check_Button *sort_button, *reverse_sort_button, *apply_orient_button, *display_label_button;
	Fl_Choice *display_choice, *sort_choice, *colour_scheme_choice;

	// Constructor with w x h size of the window and a title
	displayerGuiWindow(int W, int H, const char* title=0): Fl_Window(W, H, title),	sort_button(NULL),
			reverse_sort_button(NULL), apply_orient_button(NULL), display_label_button(NULL){}

	// Fill all except for the browser
	int fill(FileName &fn_in);

	// NUll-check value-fetch
	bool getValue(Fl_Check_Button * button)
	{
		if(button != NULL)
			return(button->value());
		else
			return(false);
	}
private:

	static void cb_display(Fl_Widget*, void*);
    inline void cb_display_i();
    void readLastSettings();
    void writeLastSettings();

};


class Displayer
{
public:

	// I/O Parser
	IOParser parser;

	// Launch the GUI for parameter input
	bool do_gui;

	// Verbosity
	int verb;

	// Which metadatalabel to display
	EMDLabel display_label, sort_label, text_label;

	// Use random sort
	bool random_sort;

	// use reverse order for sorting?
	bool reverse_sort;

	// Scale factor for displaying
	RFLOAT scale;

	// Number of rows for tiled view
	int nrow, ncol;

	// Apply orientations stored in metadatatable
	bool do_apply_orient;

	// Scale for showing the original image
	RFLOAT ori_scale;

	// Black and white values
	RFLOAT minval, maxval;

	// For setting black and white contrast to a specified times the image standard deviation from the mean
	RFLOAT sigma_contrast;

	// Particle diameter
	int particle_radius;

	// Scale for rlnCoordinateX/Y
	RFLOAT coord_scale;

	// Input & Output rootname
	FileName fn_in;

	// Ignore optics groups
	bool do_ignore_optics;

	// Filename for coordinates star file
	FileName fn_coords;

	// FileName of the STAR file that contains the color label
	FileName fn_color;

	// Which column to color on?
	FileName color_label;

	// Values for blue and red coloring
	RFLOAT color_blue_value, color_red_value;

	// Minimum value for rlnAutopickFigureOfMerit to display picks
	RFLOAT minimum_pick_fom;

	// Tablename to read from in the input STAR file
	FileName table_name;

	// Flag to pick
	bool do_pick;

	// Flag to pick start-end
	bool do_pick_startend;

	// Flag for looking at classes
	bool do_class;

	// Allow saving of selected particles or images?
	bool do_allow_save;

	// remove duplicate-coord particles in fn_data STAR file?
	bool do_remove_duplicates;
	RFLOAT duplicate_threshold;	 // Angstrom radius

	// Filenames for selected particles and selected images
	FileName fn_selected_imgs, fn_selected_parts;

	// Select maximum this number of particles from each selected classes
	int max_nr_parts_per_class;

	// Number of groups for regrouping (negative number is no regrouping)
	int nr_regroups;

	// Re-center class averages to their center-of-mass?
	bool do_recenter;

	// Flag for reading whole stacks instead of individual images
	bool do_read_whole_stacks;

	// Flag to show colour scalebar image
	bool do_colourbar;

	// data.star metadata (for do_class)
	MetaDataTable MDdata;

	// Observation model
	ObservationModel obsModel;

	// model_groups  metadata (for do_class and regrouping)
	MetaDataTable MDgroups;

	// Input metadata
	MetaDataTable MDin;

	// For the multiviewer
	std::vector<DisplayBox*> boxes;

	// Lowpass filter for picker images
	RFLOAT lowpass;

	// Highpass filter for picker images
	RFLOAT highpass;

	// Or denoise using Topaz?
	bool do_topaz_denoise;

	// Pixel size to calculate lowpass filter in Angstroms and translations in apply_orient
	RFLOAT angpix;

	// Show Fourier amplitudes?
	bool show_fourier_amplitudes;

	// Show Fourier phase angles?
	bool show_fourier_phase_angles;

	// Only show a limited number of images
	long int max_nr_images;

	// Topaz executable (for denoising of micrographs in picking mode)
	FileName fn_topaz_exe;

	// Shell for calling Topaz
	FileName fn_shell;

public:
	// Read command line arguments
	void read(int argc, char **argv);

	// Print usage instructions
	void usage();

	// Initialise some general stuff after reading
	void initialise();

	// On-the-fly topaz denoising for picking micrographs
	void topazDenoiseMap(FileName fn_in, FileName fn_odir, Image<RFLOAT> &img);

	// Decide what to do
	void run();

	// run the GUI
	int runGui();

};


#endif /* DISPLAYER_H_ */
