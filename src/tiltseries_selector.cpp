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

#include "src/tiltseries_selector.h"

#include <FL/fl_ask.H>

#include <cstdlib>
#include <iostream>
#include <sstream>

#include "src/error.h"
#include "src/funcs.h"
#include "src/gui_entries.h"
#include "src/metadata_label.h"

namespace
{
std::string shellQuote(const FileName &filename)
{
    const std::string input = filename;
    std::string output("'");

    for (size_t i = 0; i < input.size(); i++)
    {
        if (input[i] == '\'')
            output += "'\\''";
        else
            output += input[i];
    }

    output += "'";
    return output;
}

void ensureTrailingSlash(FileName &directory)
{
    if (directory.length() > 0 && directory[directory.length() - 1] != '/')
        directory += "/";
}

void createParentDirectory(const FileName &filename)
{
    const FileName directory = filename.beforeLastOf("/");
    if (directory != "" && directory != filename && !exists(directory))
        mktree(directory);
}
}

TiltSeriesSelectorGuiWindow::TiltSeriesSelectorGuiWindow(
        int W, int H, const char *title)
    : Fl_Window(W, H, title),
      do_allow_save(false),
      image_scale(1.), image_sigma_contrast(0.), lowpass(-1.),
      ctf_scale(1.), ctf_sigma_contrast(3.),
      tomogram_scale(1), tomogram_sigma_contrast(3.),
      tomogram_lowpass(-1.), nr_threads(1)
{
}

FileName TiltSeriesSelectorGuiWindow::makeOutputTiltSeriesName(
        int index) const
{
    if (do_allow_save)
    {
        // Use the normal TomogramSet layout. This also avoids depending on the
        // directory layout of the input job and gives each tomogram a stable,
        // writable STAR filename below the selector output directory.
        return fn_odir + "tilt_series/" + FileName(tomogram_set.getTomogramName(index)) + ".star";
    }
    else
    {
        // Just work with the original star files
        FileName fn_star;
        tomogram_set.globalTable.getValueSafely(EMDL_TOMO_TILT_SERIES_STARFILE, fn_star, index);
        return fn_star;
    }
}

void TiltSeriesSelectorGuiWindow::ensureOutputDirectories() const
{
    if (fn_odir != "" && !exists(fn_odir))
        mktree(fn_odir);

    createParentDirectory(fn_out);

    for (size_t i = 0; i < output_tiltseries_stars.size(); i++)
        createParentDirectory(output_tiltseries_stars[i]);
}


void TiltSeriesSelectorGuiWindow::materialiseOutputFiles()
{
    if (!do_allow_save) return;

    ensureOutputDirectories();

    for (int index = 0; index < tomogram_set.size(); index++)
    {
        const FileName &working_star = output_tiltseries_stars[index];
        const std::string tomogram_name = tomogram_set.getTomogramName(index);

        if (!exists(working_star))
        {
            // First launch: create a complete writable copy for every tilt
            // series before any Images or CTF button can be pressed.
            tomogram_set.tomogramTables[index].write(working_star);
        }

        tomogram_set.globalTable.setValue(EMDL_TOMO_TILT_SERIES_STARFILE, working_star, index);

    }
}

void TiltSeriesSelectorGuiWindow::writeOuterSelection(bool verb)
{
    if (!do_allow_save) return;

    TomogramSet output_set = tomogram_set;

    for (int index = output_set.size() - 1; index >= 0; index--)
    {
        if (!selected[index])
            output_set.removeTomogram(output_set.getTomogramName(index));
    }

    output_set.globalTable.write(fn_out);
    writeLocalSelection();

    if (verb) std::cout << "Saved " << output_set.size() << " selected tilt series in: " << fn_out << std::endl;
}

FileName TiltSeriesSelectorGuiWindow::getTomogramFilename(int index) const
{
    FileName filename;

    // Prefer a denoised tomogram when both columns are present.
    if (tomogram_set.globalTable.containsLabel(
            EMDL_TOMO_DENOISED_TOMOGRAM_FILE_NAME))
    {
        tomogram_set.globalTable.getValueSafely(
                EMDL_TOMO_DENOISED_TOMOGRAM_FILE_NAME,
                filename, index);
        if (filename != "")
            return filename;
    }

    if (tomogram_set.globalTable.containsLabel(
            EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_FILE_NAME))
    {
        tomogram_set.globalTable.getValueSafely(
                EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_FILE_NAME,
                filename, index);
    }

    return filename;
}

void TiltSeriesSelectorGuiWindow::selectRow(int index, bool value)
{
    if (index < 0 || index >= selected.size())
        return;

    selected[index] = value;
    check_buttons[index]->value(value ? 1 : 0);

    if (value)
    {
        name_boxes[index]->color(GUI_INPUT_COLOR);
        name_boxes[index]->activate();
        image_buttons[index]->activate();

        if (tomogram_set.tomogramTables[index].containsLabel(EMDL_CTF_IMAGE))
            ctf_buttons[index]->activate();

        if (tomogram_files[index] != "")
            tomogram_buttons[index]->activate();
    }
    else
    {
        name_boxes[index]->color(GUI_BACKGROUND_COLOR);
        name_boxes[index]->deactivate();
        image_buttons[index]->deactivate();
        ctf_buttons[index]->deactivate();
        tomogram_buttons[index]->deactivate();
    }

    name_boxes[index]->redraw();
    image_buttons[index]->redraw();
    ctf_buttons[index]->redraw();
    tomogram_buttons[index]->redraw();
}

void TiltSeriesSelectorGuiWindow::cb_select(Fl_Widget *widget, void *data)
{
    RowCallbackData *row = static_cast<RowCallbackData *>(data);
    Fl_Check_Button *button = static_cast<Fl_Check_Button *>(widget);
    row->window->selectRow(row->index, button->value() != 0);
}

void TiltSeriesSelectorGuiWindow::launchViewer(int index, bool show_ctf)
{
    if (index < 0 || index >= output_tiltseries_stars.size())
        return;

    std::ostringstream command;
    command << "relion_display"
            << " --i " << shellQuote(output_tiltseries_stars[index])
            << " --fn_imgs " << shellQuote(output_tiltseries_stars[index])
            << " --j " << nr_threads;

    if (do_allow_save)
    {
        command << " --allow_save";
    }

    if (show_ctf)
    {
        command << " --display rlnCtfImage"
                << " --scale " << ctf_scale
                << " --sigma_contrast " << ctf_sigma_contrast;
    }
    else
    {
        command << " --scale " << image_scale
                << " --sigma_contrast " << image_sigma_contrast
                << " --lowpass " << lowpass;

    }

    command << " &";

    //std::cerr << " command= " << command.str().c_str() <<std::endl;
    const int result = std::system(command.str().c_str());
    if (result != 0)
        std::cerr << "WARNING: failed to launch: " << command.str() << std::endl;

    std::vector<Fl_Button *> &buttons = show_ctf ? ctf_buttons : image_buttons;
    for (size_t i = 0; i < buttons.size(); i++)
    {
        buttons[i]->color(i == index ? GUI_BUTTON_DARK_COLOR : GUI_BUTTON_COLOR);
        buttons[i]->redraw();
    }
}

void TiltSeriesSelectorGuiWindow::launchTomogramViewer(int index)
{
    if (index < 0 || index >= tomogram_files.size())
        return;

    const FileName &filename = tomogram_files[index];
    if (filename == "")
        return;

    std::ostringstream command;
    command << "relion_display"
            << " --i " << shellQuote(filename)
            << " --scale " << tomogram_scale
            << " --sigma_contrast " << tomogram_sigma_contrast
            << " --lowpass " << tomogram_lowpass
            << " &";

    const int result = std::system(command.str().c_str());
    if (result != 0)
        std::cerr << "WARNING: failed to launch: " << command.str() << std::endl;

    for (size_t i = 0; i < tomogram_buttons.size(); i++)
    {
        tomogram_buttons[i]->color(
                i == index ? GUI_BUTTON_DARK_COLOR : GUI_BUTTON_COLOR);
        tomogram_buttons[i]->redraw();
    }
}

void TiltSeriesSelectorGuiWindow::cb_images(Fl_Widget *, void *data)
{
    RowCallbackData *row = static_cast<RowCallbackData *>(data);
    row->window->launchViewer(row->index, false);
}

void TiltSeriesSelectorGuiWindow::cb_ctf(Fl_Widget *, void *data)
{
    RowCallbackData *row = static_cast<RowCallbackData *>(data);
    row->window->launchViewer(row->index, true);
}

void TiltSeriesSelectorGuiWindow::cb_tomogram(Fl_Widget *, void *data)
{
    RowCallbackData *row = static_cast<RowCallbackData *>(data);
    row->window->launchTomogramViewer(row->index);
}

void TiltSeriesSelectorGuiWindow::writeLocalSelection() const
{
    MetaDataTable local_selection;
    local_selection.setName("selection");

    for (int index = 0; index < tomogram_set.size(); index++)
    {
        local_selection.addObject();
        local_selection.setValue(EMDL_TOMO_NAME, tomogram_set.getTomogramName(index));
        local_selection.setValue(EMDL_SELECTED, selected[index] ? 1 : 0);
    }

    local_selection.write(fn_odir + "local_selection.star");
}

void TiltSeriesSelectorGuiWindow::readLocalSelection()
{
    const FileName filename = fn_odir + "local_selection.star";
    if (!exists(filename))
        return;

    MetaDataTable local_selection;
    local_selection.read(filename);

    for (long int row = 0; row < local_selection.numberOfObjects(); row++)
    {
        std::string name;
        int is_selected = 1;
        local_selection.getValue(EMDL_TOMO_NAME, name, row);
        local_selection.getValue(EMDL_SELECTED, is_selected, row);

        const int index = tomogram_set.getTomogramIndex(name);
        if (index >= 0)
            selectRow(index, is_selected != 0);
    }
}

void TiltSeriesSelectorGuiWindow::saveSelection(bool verb)
{
    if (!do_allow_save)
        return;

    int number_selected = 0;
    for (size_t i = 0; i < selected.size(); i++)
        if (selected[i]) number_selected++;

    if (number_selected == 0)
    {
        fl_message("No tilt series are selected.");
        return;
    }

    writeOuterSelection(verb);
}

void TiltSeriesSelectorGuiWindow::selectAll()
{
    for (int index = 0; index < selected.size(); index++)
        selectRow(index, true);
}

void TiltSeriesSelectorGuiWindow::invertSelection()
{
    for (int index = 0; index < selected.size(); index++)
        selectRow(index, !selected[index]);
}

void TiltSeriesSelectorGuiWindow::quit()
{
    hide();
}

void TiltSeriesSelectorGuiWindow::cb_menubar_save(Fl_Widget *, void *data)
{
    TiltSeriesSelectorGuiWindow *window =
            static_cast<TiltSeriesSelectorGuiWindow *>(data);
    window->saveSelection();
    RELION_EXIT_SUCCESS;
}

void TiltSeriesSelectorGuiWindow::cb_menubar_select_all(
        Fl_Widget *, void *data)
{
    static_cast<TiltSeriesSelectorGuiWindow *>(data)->selectAll();
}

void TiltSeriesSelectorGuiWindow::cb_menubar_invert_selection(
        Fl_Widget *, void *data)
{
    static_cast<TiltSeriesSelectorGuiWindow *>(data)->invertSelection();
}

void TiltSeriesSelectorGuiWindow::cb_menubar_quit(Fl_Widget *, void *data)
{
    static_cast<TiltSeriesSelectorGuiWindow *>(data)->quit();
}

void TiltSeriesSelectorGuiWindow::cb_closing(Fl_Widget *, void *data)
{
    static_cast<TiltSeriesSelectorGuiWindow *>(data)->quit();
}

int TiltSeriesSelectorGuiWindow::fill()
{
    color(GUI_BACKGROUND_COLOR);
    callback(cb_closing, this);

    ensureTrailingSlash(fn_odir);

    // The Tomogram column is always present. Individual buttons are
    // deactivated below when no denoised or reconstructed tomogram exists.
    if (w() < TSS_TOTAL_WIDTH)
        resize(x(), y(), TSS_TOTAL_WIDTH, h());

    Fl_Menu_Bar *menubar = new Fl_Menu_Bar(0, 0, w(), 25);
    if (do_allow_save)
    {
        menubar->add("File/&Save selection", 's',
                     cb_menubar_save, this);
        menubar->add("File/&Invert selection", 'i',
                     cb_menubar_invert_selection, this);
        menubar->add("File/&All selection", 'a',
                     cb_menubar_select_all, this);
    }
    menubar->add("File/&Quit", 'q', cb_menubar_quit, this);

    const int count = tomogram_set.size();
    const int row_height = 35;
    int current_y = 25;

    selected.assign(count, true);
    output_tiltseries_stars.resize(count);
    tomogram_files.assign(count, FileName(""));
    callback_data.resize(count);

    for (int index = 0; index < count; index++)
    {
        output_tiltseries_stars[index] = makeOutputTiltSeriesName(index);
        tomogram_files[index] = getTomogramFilename(index);
    }

    // Create all individual working STAR files and redirect the in-memory
    // global table before constructing buttons that launch relion_display.
    if (do_allow_save) materialiseOutputFiles();

    check_buttons.clear();
    name_boxes.clear();
    image_buttons.clear();
    ctf_buttons.clear();
    tomogram_buttons.clear();

    Fl_Scroll *scroll = new Fl_Scroll(0, current_y, w(), h() - current_y);
    scroll->type(Fl_Scroll::VERTICAL);
    scroll->begin();

    for (int index = 0; index < count; index++)
    {
        callback_data[index].window = this;
        callback_data[index].index = index;

        Fl_Check_Button *check = new Fl_Check_Button(
                TSS_CHECK_X, current_y, row_height - 8, row_height - 8, "");
        check->value(1);
        check->callback(cb_select, &callback_data[index]);
        if (!do_allow_save)
            check->deactivate();
        check_buttons.push_back(check);

        Fl_Box *name = new Fl_Box(
                TSS_NAME_X, current_y, TSS_NAME_W, row_height - 5);
        name->box(FL_DOWN_BOX);
        name->color(GUI_INPUT_COLOR);
        name->align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE | FL_ALIGN_CLIP);
        name->copy_label(tomogram_set.getTomogramName(index).c_str());
        name_boxes.push_back(name);

        Fl_Button *images = new Fl_Button(
                TSS_IMAGES_X, current_y, TSS_BUTTON_W, row_height - 5,
                "Series");
        images->color(GUI_BUTTON_COLOR);
        images->callback(cb_images, &callback_data[index]);
        image_buttons.push_back(images);

        Fl_Button *ctf = new Fl_Button(
                TSS_CTF_X, current_y, TSS_BUTTON_W, row_height - 5,
                "CTF");
        ctf->color(GUI_BUTTON_COLOR);
        ctf->callback(cb_ctf, &callback_data[index]);
        if (!tomogram_set.tomogramTables[index].containsLabel(EMDL_CTF_IMAGE))
            ctf->deactivate();
        ctf_buttons.push_back(ctf);

        Fl_Button *tomogram = new Fl_Button(
                TSS_TOMOGRAM_X, current_y,
                TSS_BUTTON_W, row_height - 5, "Tomo");
        tomogram->color(GUI_BUTTON_COLOR);
        tomogram->callback(cb_tomogram, &callback_data[index]);
        if (tomogram_files[index] == "")
            tomogram->deactivate();
        tomogram_buttons.push_back(tomogram);

        current_y += row_height;
    }

    scroll->end();

    if (do_allow_save)
        readLocalSelection();

    // Always materialise the complete output TomogramSet on launch. The
    // individual files were created above; now write the outer selection with
    // its rlnTomoTiltSeriesStarFile entries pointing into fn_odir.
    writeOuterSelection(false);

    resizable(scroll);
    show();
    return Fl::run();
}

void TiltSeriesSelector::read(int argc, char **argv)
{
    parser.setCommandLine(argc, argv);

    parser.addSection("General options");
    fn_in = parser.getOption("--i", "Input tilt_series.star containing a TomogramSet");
    fn_odir = parser.getOption("--odir","Output directory for edited per-tilt-series STAR files","ExcludeTiltImages/test/");
    ensureTrailingSlash(fn_odir);
    fn_out = parser.getOption("--selection","Output tilt_series.star containing selected tilt series",(fn_odir + "tilt_series.star").c_str());
    do_allow_save = parser.checkOption("--allow_save", "Allow saving the selected tilt series");

    parser.addSection("Display options");
    image_scale = textToFloat(parser.getOption("--scale", "Relative scale for tilt-image display", "1"));
    image_sigma_contrast = textToFloat(parser.getOption("--sigma_contrast","Sigma contrast for tilt-image display", "0"));
    lowpass = textToFloat(parser.getOption("--lowpass", "Lowpass filter to be applied to tiltseries images (a bit slower)", "-1"));
    ctf_scale = textToFloat(parser.getOption("--ctf_scale", "Relative scale for CTF-image display", "1"));
    ctf_sigma_contrast = textToFloat(parser.getOption("--ctf_sigma_contrast","Sigma contrast for CTF-image display", "3"));
    tomogram_scale = textToFloat(parser.getOption("--tomogram_scale", "Relative scale for tomogram display", "0.5"));
    tomogram_sigma_contrast = textToFloat(parser.getOption("--tomogram_sigma_contrast", "Sigma contrast for tomogram display", "3"));
    tomogram_lowpass = textToFloat(parser.getOption("--tomogram_lowpass", "Lowpass filter for tomogram display", "-1"));
    nr_threads = textToInteger(parser.getOption("--j", "Number of threads used by relion_display", "12"));

    if (parser.checkForErrors()) REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}

void TiltSeriesSelector::usage()
{
    parser.writeUsage(std::cout);
}

void TiltSeriesSelector::initialise()
{
    if (!fn_in.isStarFile())
        REPORT_ERROR("TiltSeriesSelector expects a tilt_series.star input file");

    if (!tomogram_set.read(fn_in))
        REPORT_ERROR("Could not read a TomogramSet from " + fn_in);

    ensureTrailingSlash(fn_odir);
    if (fn_odir != "" && !exists(fn_odir))
        mktree(fn_odir);
    createParentDirectory(fn_out);
}

void TiltSeriesSelector::run()
{
    Fl::scheme("gtk+");

    TiltSeriesSelectorGuiWindow window(
            TSS_TOTAL_WIDTH, TSS_TOTAL_HEIGHT,
            "RELION tilt-series selector");

    window.tomogram_set = tomogram_set;
    window.fn_in = fn_in;
    window.fn_out = fn_out;
    window.fn_odir = fn_odir;
    window.do_allow_save = do_allow_save;
    window.image_scale = image_scale;
    window.image_sigma_contrast = image_sigma_contrast;
    window.lowpass = lowpass;
    window.ctf_scale = ctf_scale;
    window.ctf_sigma_contrast = ctf_sigma_contrast;
    window.tomogram_scale = tomogram_scale;
    window.tomogram_sigma_contrast = tomogram_sigma_contrast;
    window.tomogram_lowpass = tomogram_lowpass;
    window.nr_threads = nr_threads;

    window.fill();
}
