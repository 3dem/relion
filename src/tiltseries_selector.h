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

#ifndef TILTSERIES_SELECTOR_H_
#define TILTSERIES_SELECTOR_H_

// Protect against FLTK's Complex definition.
#define Complex tmpComplex
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Scroll.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Box.H>
#undef Complex

#include <vector>

#include "src/args.h"
#include "src/filename.h"
#include "src/jaz/tomography/tomogram_set.h"

#define TSS_CHECK_X       4
#define TSS_NAME_X       34
#define TSS_NAME_W      140
#define TSS_IMAGES_X    184
#define TSS_BUTTON_W     75
#define TSS_CTF_X       269
#define TSS_TOMOGRAM_X  354
#define TSS_TOTAL_WIDTH 443
#define TSS_TOTAL_HEIGHT 800

class TiltSeriesSelectorGuiWindow : public Fl_Window
{
public:
    TomogramSet tomogram_set;
    FileName fn_in;
    FileName fn_out;
    FileName fn_odir;

    bool do_allow_save;

    RFLOAT image_scale;
    RFLOAT image_sigma_contrast;
    RFLOAT lowpass;
    RFLOAT ctf_scale;
    RFLOAT ctf_sigma_contrast;
    RFLOAT tomogram_scale;
    RFLOAT tomogram_sigma_contrast;
    RFLOAT tomogram_lowpass;
    int nr_threads;

    TiltSeriesSelectorGuiWindow(int W, int H, const char *title = 0);

    int fill();

private:
    struct RowCallbackData
    {
        TiltSeriesSelectorGuiWindow *window;
        int index;
    };

    std::vector<bool> selected;
    std::vector<FileName> output_tiltseries_stars;
    std::vector<FileName> tomogram_files;
    std::vector<RowCallbackData> callback_data;

    std::vector<Fl_Check_Button *> check_buttons;
    std::vector<Fl_Box *> name_boxes;
    std::vector<Fl_Button *> image_buttons;
    std::vector<Fl_Button *> ctf_buttons;
    std::vector<Fl_Button *> tomogram_buttons;

    static void cb_select(Fl_Widget *, void *);
    static void cb_images(Fl_Widget *, void *);
    static void cb_ctf(Fl_Widget *, void *);
    static void cb_tomogram(Fl_Widget *, void *);

    static void cb_menubar_save(Fl_Widget *, void *);
    static void cb_menubar_select_all(Fl_Widget *, void *);
    static void cb_menubar_invert_selection(Fl_Widget *, void *);
    static void cb_menubar_quit(Fl_Widget *, void *);
    static void cb_closing(Fl_Widget *, void *);

    void selectRow(int index, bool value);
    void launchViewer(int index, bool show_ctf);
    void launchTomogramViewer(int index);
    FileName getTomogramFilename(int index) const;

    void saveSelection(bool verb = true);
    void readLocalSelection();
    void writeLocalSelection() const;

    void selectAll();
    void invertSelection();
    void quit();

    FileName makeOutputTiltSeriesName(int index) const;
    void ensureOutputDirectories() const;
    void materialiseOutputFiles();
    void refreshOutputFiles();
    void writeOuterSelection(bool verb = true);
};

class TiltSeriesSelector
{
public:
    IOParser parser;

    FileName fn_in;
    FileName fn_out;
    FileName fn_odir;

    bool do_allow_save;
    RFLOAT image_scale;
    RFLOAT image_sigma_contrast;
    RFLOAT lowpass;
    RFLOAT ctf_scale;
    RFLOAT ctf_sigma_contrast;
    RFLOAT tomogram_scale;
    RFLOAT tomogram_sigma_contrast;
    RFLOAT tomogram_lowpass;
    int nr_threads;

    TomogramSet tomogram_set;

    void read(int argc, char **argv);
    void usage();
    void initialise();
    void run();
};

#endif /* TILTSERIES_SELECTOR_H_ */
