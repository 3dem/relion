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

#include "src/displayer.h"
//#define DEBUG
//
#ifdef HAVE_PNG
#include <src/jaz/gravis/tImage.h>
#endif

const Fl_Menu_Item color_choices[] =
{
	// text, shortcut, callback, user_data, flags, type, font, size, color
	{"Red (1)",     0, (Fl_Callback*)0, (void*)1, 0, 0, 0, 0, FL_RED},
	{"Green (2)",   0, (Fl_Callback*)0, (void*)2, 0, 0, 0, 0, FL_GREEN},
	{"Blue (3)",    0, (Fl_Callback*)0, (void*)3, 0, 0, 0, 0, FL_BLUE},
	{"Cyan (4)",    0, (Fl_Callback*)0, (void*)4, 0, 0, 0, 0, FL_CYAN},
	{"Magenta (5)", 0, (Fl_Callback*)0, (void*)5, 0, 0, 0, 0, FL_MAGENTA},
	{"Yellow (6)",  0, (Fl_Callback*)0, (void*)6, 0, 0, 0, 0, FL_YELLOW},
	{0} // sentinel
};
const int NUM_COLORS = 6;

/************************************************************************/
void DisplayBox::draw()
{
	if (!img_data) return;

	short xpos = x() + xoff;
	short ypos = y() + yoff;

	/* Ensure that the full window is redrawn */
	//fl_push_clip(x(),y(),w(),h());

	/* Redraw the whole image */
	int depth = (colour_scheme) ? 3 : 1;
	fl_draw_image((const uchar *)img_data, xpos, ypos, (short)xsize_data, (short)ysize_data, depth);
	if (img_label != "")
	{
		fl_color(FL_WHITE);
		fl_draw(img_label.c_str(), xpos, ypos + fl_height());
	}
	/* Draw a red rectangle around the particle if it is selected */
	if (selected >= 1 && selected <= 6)
		fl_color(color_choices[selected - 1].labelcolor_);
	else
		fl_color(FL_BLACK);

	fl_line_style(FL_SOLID, 2);
	int x1 = xpos;
	int y1 = ypos;
	int x2 = xpos + xsize_data;
	int y2 = ypos + ysize_data;
	fl_line(x1, y1, x1, y2);
	fl_line(x1, y2, x2, y2);
	fl_line(x2, y2, x2, y1);
	fl_line(x2, y1, x1, y1);

	//fl_pop_clip();
}

unsigned char rgbToGrey(const unsigned char red, const unsigned char green, const unsigned char blue)
{
	switch (colour_scheme)
	{
	case (BLACKGREYREDSCALE):
	{
		if (red == 255) return FLOOR((RFLOAT)(255. - blue/2.));
		else return FLOOR((RFLOAT)(red/2.));
		break;
	}
	case (BLUEGREYWHITESCALE):
	{
		if (red == 0) return FLOOR((RFLOAT)(255.-blue)/2.);
		else return FLOOR((RFLOAT)(red/2. + 128.));
		break;
	}
	case (BLUEGREYREDSCALE):
	{
		unsigned char Y;
		int X;
		if (red == 0) { Y = 255-blue; X = 0; }
		else if (red == 255) { Y = 255-blue; X = 2; }
		else { Y = blue; X = 1; }
		return CEIL(85*((RFLOAT)Y/256. + X));
		break;
	}
	case (RAINBOWSCALE):
	{
		unsigned char Y;
		int X;
		if (red > 0)
		{
			if (red == 255) {Y = green; X = 0;}
			else {Y = 255-red; X = 1;}
		}
		else if (green > 0)
		{
			if (green == 255) {Y = blue; X = 2;}
			else {Y = 255 - green; X = 3;}
		}
		else {Y = 255; X = 4;}

		return 255 - CEIL(64*((RFLOAT)Y/255. + X));
		break;
	}
	case (CYANBLACKYELLOWSCALE):
	{
		if (red >0)
		{
			if (red < 255) return (unsigned char)FLOOR((RFLOAT)red / 3. + 128);
			else return (unsigned char)FLOOR((RFLOAT)green/3. + 42 + 128);
		}
		else
		{
			if (blue < 255) return (unsigned char)FLOOR((RFLOAT)-blue / 3. + 128);
			else return (unsigned char)FLOOR(-((RFLOAT)green)/3. - 42 + 128);
		}
		break;
	}
	}

	REPORT_ERROR("Logic error: should not happen");
	return 0;
}

void DisplayBox::setData(MultidimArray<RFLOAT> &img, MetaDataContainer *MDCin, int _ipos,
                         RFLOAT _minval, RFLOAT _maxval, RFLOAT _scale, bool do_relion_scale)
{
	scale = _scale;
	minval = _minval;
	maxval = _maxval;
	ipos = _ipos;
	selected = NOTSELECTED;

	// Set its own MetaDataTable
	MDimg.setIsList(true);
	MDimg.addObject(MDCin);

	// For volumes only show the central slice
	if (ZSIZE(img) > 1)
	{
		MultidimArray<RFLOAT> slice;
		img.getSlice(ZSIZE(img)/2, slice);
		img=slice;
	}

	// create array for the scaled image data
	xsize_data = CEIL(XSIZE(img) * scale);
	ysize_data = CEIL(YSIZE(img) * scale);
	xoff = (xsize_data < w() ) ? (w() - xsize_data) / 2 : 0;
	yoff = (ysize_data < h() ) ? (h() - ysize_data) / 2 : 0;
	if (colour_scheme == GREYSCALE)
	{
		img_data = new unsigned char [xsize_data * ysize_data];
	}
	else
	{
		img_data = new unsigned char [3 * xsize_data * ysize_data];
	}
	RFLOAT range = maxval - minval;
	RFLOAT step = range / 255; // 8-bit scaling range from 0 to 255
	RFLOAT* old_ptr=NULL;
	long int n;

	// For micrographs use relion-scaling to avoid bias in down-sampled positions
	// For multi-image viewers, do not use this scaling as it is slower...
	if (do_relion_scale && ABS(scale - 1.0) > 0.01)
		selfScaleToSize(img, xsize_data, ysize_data);

	// Use the same nearest-neighbor algorithm as in the copy function of Fl_Image...
	if (ABS(scale - 1.0) > 0.01 && !do_relion_scale)
	{
		int xmod   = XSIZE(img) % xsize_data;
		int xstep  = XSIZE(img) / xsize_data;
		int ymod   = YSIZE(img) % ysize_data;
		int ystep  = YSIZE(img) / ysize_data;
		int line_d = XSIZE(img);
		int dx, dy, sy, xerr, yerr;

		if (colour_scheme == GREYSCALE)
		{
			for (dy = ysize_data, sy = 0, yerr = ysize_data, n = 0; dy > 0; dy --)
			{
				for (dx = xsize_data, xerr = xsize_data, old_ptr = img.data + sy * line_d; dx > 0; dx --, n++)
				{
					img_data[n] = (char)FLOOR((*old_ptr - minval) / step);
					old_ptr += xstep;
					xerr    -= xmod;
					if (xerr <= 0)
					{
						xerr    += xsize_data;
						old_ptr += 1;
					}
				}

				sy   += ystep;
				yerr -= ymod;
				if (yerr <= 0)
				{
					yerr += ysize_data;
					sy ++;
				}
			}
		}
		else
		{
			for (dy = ysize_data, sy = 0, yerr = ysize_data, n = 0; dy > 0; dy --)
			{
				for (dx = xsize_data, xerr = xsize_data, old_ptr = img.data + sy * line_d; dx > 0; dx --, n++)
				{
					unsigned char val = FLOOR((*old_ptr - minval) / step);
					greyToRGB(colour_scheme, val, img_data[3*n], img_data[3*n+1], img_data[3*n+2]);
					old_ptr += xstep;
					xerr    -= xmod;
					if (xerr <= 0)
					{
						xerr    += xsize_data;
						old_ptr += 1;
					}
				}

				sy   += ystep;
				yerr -= ymod;
				if (yerr <= 0)
				{
					yerr += ysize_data;
					sy ++;
				}
			}
		}
	}
	else
	{
		if (colour_scheme == GREYSCALE)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(img, n, old_ptr)
			{
				img_data[n] = FLOOR((*old_ptr - minval) / step);
			}
		}
		else
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(img, n, old_ptr)
			{
				unsigned char val = FLOOR((*old_ptr - minval) / step);
				greyToRGB(colour_scheme, val, img_data[3*n], img_data[3*n+1], img_data[3*n+2]);
			}
		}
	}
}

int DisplayBox::toggleSelect(int set_selected)
{
	if (selected > 0)
		selected = 0;
	else if (selected == 0)
		selected = set_selected;
	redraw();
	return selected;
}

void DisplayBox::setSelect(int value)
{
	selected = value;
	redraw();
}

int DisplayBox::select()
{

	selected = SELECTED;
	redraw();
	return selected;
}

int DisplayBox::unSelect()
{
	selected = NOTSELECTED;
	redraw();
	return selected;
}

int basisViewerWindow::fillCanvas(int viewer_type, MetaDataTable &MDin, ObservationModel *obsModel, EMDLabel display_label, EMDLabel text_label, bool _do_read_whole_stacks, bool _do_apply_orient,
                                  RFLOAT _minval, RFLOAT _maxval, RFLOAT _sigma_contrast, RFLOAT _scale, RFLOAT _ori_scale, int _ncol, long int max_nr_images, RFLOAT lowpass, RFLOAT highpass, bool _do_class,
                                  MetaDataTable *_MDdata, int _nr_regroup, bool _do_recenter,  bool _is_data, MetaDataTable *_MDgroups,
                                  bool do_allow_save, FileName fn_selected_imgs, FileName fn_selected_parts, int max_nr_parts_per_class)
{
	// Scroll bars
	Fl_Scroll scroll(0, 0, w(), h());

	// Pre-set the canvas to the correct size
	FileName fn_img;
	Image<RFLOAT> img;
	MDin.firstObject();
	MDin.getValue(display_label, fn_img);
	img.read(fn_img, false);
	int nimgs = MDin.numberOfObjects();
	if (viewer_type == MULTIVIEWER)
	{
		int xsize_canvas = _ncol * (CEIL(XSIZE(img())*_scale) + BOX_OFFSET);
		int nrow = CEIL((RFLOAT)nimgs/_ncol);
		int ysize_canvas = nrow * (CEIL(YSIZE(img())*_scale) + BOX_OFFSET);
		multiViewerCanvas canvas(0, 0, xsize_canvas, ysize_canvas);
		canvas.multi_max_nr_images = max_nr_images;
		canvas.SetScroll(&scroll);
		canvas.do_read_whole_stacks = _do_read_whole_stacks;
		canvas.is_data = _is_data;
		canvas.ori_scale = _ori_scale;
		canvas.display_label = display_label;
		canvas.sigma_contrast = _sigma_contrast;
		canvas.minval = _minval;
		canvas.maxval = _maxval;
		canvas.do_allow_save = do_allow_save;
		canvas.fn_selected_imgs= fn_selected_imgs;
		canvas.fn_selected_parts = fn_selected_parts;
		canvas.max_nr_parts_per_class = max_nr_parts_per_class;
		canvas.fill(MDin, obsModel, display_label, text_label, _do_apply_orient, _minval, _maxval, _sigma_contrast, _scale, _ncol, _do_recenter, max_nr_images, lowpass, highpass);
		canvas.nr_regroups = _nr_regroup;
		canvas.do_recenter = _do_recenter;
		canvas.do_apply_orient = _do_apply_orient;
		canvas.obsModel = obsModel;
		canvas.text_label = text_label;
		canvas.metadata_table_name = MDin.getName();
		if (canvas.nr_regroups > 0)
			canvas.MDgroups = _MDgroups;
		if (_do_class)
		{
			canvas.do_class = true;
			canvas.MDdata = _MDdata;
		}
		else
		{
			canvas.do_class = false;
		}

		// Pre-load existing backup_selection.star file
		FileName fn_sel, fn_dir=".";
		if (fn_selected_imgs != "")
			fn_dir = fn_selected_imgs.beforeLastOf("/");
		else if (fn_selected_parts != "")
			fn_dir = fn_selected_parts.beforeLastOf("/");

		fn_dir += "/backup_selection.star";
		if (exists(fn_dir))
			canvas.loadBackupSelection(false); // false means dont ask for filename

		resizable(*this);
		show();
		return Fl::run();
	}
	else if (viewer_type == SINGLEVIEWER)
	{
		if (nimgs>1)
			REPORT_ERROR("ERROR: trying to launch a singleViewerCanvas with multiple images...");
		int xsize_canvas = CEIL(XSIZE(img())*_scale);
		int ysize_canvas = CEIL(YSIZE(img())*_scale);
		singleViewerCanvas canvas(0, 0, xsize_canvas, ysize_canvas);
		canvas.SetScroll(&scroll);
		canvas.fill(MDin, obsModel, display_label, text_label, _do_apply_orient, _minval, _maxval, _sigma_contrast, _scale, 1);
		canvas.do_read_whole_stacks = false;
		resizable(*this);
		show();
		return Fl::run();
	}

	REPORT_ERROR("Logic error: should not come here");
	return -1;
}

int basisViewerWindow::fillPickerViewerCanvas(MultidimArray<RFLOAT> image, RFLOAT _minval, RFLOAT _maxval, RFLOAT _sigma_contrast,
                                              RFLOAT _scale, RFLOAT _coord_scale, int _particle_radius, bool _do_startend, FileName _fn_coords,
                                              FileName _fn_color, FileName _fn_mic, FileName _color_label, RFLOAT _color_blue_value, RFLOAT _color_red_value,
											  RFLOAT _minimum_pick_fom)
{
	current_selection_type = 2; // Green

	// Scroll bars
	Fl_Scroll scroll(0, 0, w(), h());
	int xsize_canvas = CEIL(XSIZE(image)*_scale);
	int ysize_canvas = CEIL(YSIZE(image)*_scale);
	pickerViewerCanvas canvas(0, 0, xsize_canvas, ysize_canvas);
	canvas.particle_radius = _particle_radius;
	canvas.do_startend = _do_startend;
	canvas.coord_scale = _coord_scale;
	canvas.SetScroll(&scroll);
	canvas.fill(image, _minval, _maxval, _sigma_contrast, _scale);
	canvas.fn_coords = _fn_coords;
	canvas.fn_color = _fn_color;
	canvas.fn_mic = _fn_mic;
	canvas.color_label = EMDL::str2Label(_color_label);
	canvas.smallest_color_value = XMIPP_MIN(_color_blue_value, _color_red_value);
	canvas.biggest_color_value = XMIPP_MAX(_color_blue_value, _color_red_value);
	canvas.minimum_pick_fom = _minimum_pick_fom;
	canvas.do_blue_to_red = (_color_blue_value < _color_red_value);
	canvas.do_read_whole_stacks = false;
	if (_fn_coords != "" && exists(_fn_coords))
	{
		canvas.loadCoordinates(false);
		canvas.redraw();
	}
	resizable(*this);
	show();
	return Fl::run();
}

int basisViewerWindow::fillSingleViewerCanvas(MultidimArray<RFLOAT> image, RFLOAT _minval, RFLOAT _maxval, RFLOAT _sigma_contrast, RFLOAT _scale)
{
	// Scroll bars
	Fl_Scroll scroll(0, 0, w(), h());

	// Pre-set the canvas to the correct size
	int xsize_canvas = CEIL(XSIZE(image)*_scale);
	int ysize_canvas = CEIL(YSIZE(image)*_scale);
	singleViewerCanvas canvas(0, 0, xsize_canvas, ysize_canvas);
	canvas.SetScroll(&scroll);
	canvas.fill(image, _minval, _maxval, _sigma_contrast, _scale);
	canvas.do_read_whole_stacks = false;
	resizable(*this);
	show();
	return Fl::run();
}

void basisViewerCanvas::fill(MetaDataTable &MDin, ObservationModel *obsModel, EMDLabel display_label, EMDLabel text_label, bool _do_apply_orient, RFLOAT _minval, RFLOAT _maxval,
                            RFLOAT _sigma_contrast, RFLOAT _scale, int _ncol, bool _do_recenter, long int max_images, RFLOAT lowpass, RFLOAT highpass)
{
	ncol = _ncol;
	int nr_imgs = MDin.numberOfObjects();
	if (nr_imgs > 1)
	{
		xoff = BOX_OFFSET/2;
		yoff = BOX_OFFSET/2;
	}
	else
	{
		xoff = 0;
		yoff = 0;
	}

	int barstep;
	if (nr_imgs > 1)
	{
		std::cout << "Reading in all images..." << std::endl;
		init_progress_bar(nr_imgs);
		barstep = XMIPP_MAX(1, nr_imgs/ 60);
	}

	nrow = 0;
	long int ipos = 0;
	int irow = 0;
	int icol = 0;
	FileName fn_my_stack, fn_next_stack, fn_img, fn_tmp;
	long int my_number, my_next_number, my_stack_first_ipos = 0;
	std::vector<long int> numbers_in_stack;

	long int number_of_images = MDin.numberOfObjects();
	if (max_images > 0 && max_images < number_of_images)
		number_of_images = max_images;
	boxes.clear();
	boxes.resize(number_of_images);
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
	{
		// Read in image stacks as a whole, i.e. don't re-open and close stack for every individual image to save speed
		MDin.getValue(display_label, fn_img, ipos);
		fn_img.decompose(my_number, fn_my_stack);

		// See whether the next image has the same stackname....
		if (ipos+1 < number_of_images)
		{
			MDin.getValue(display_label, fn_tmp, ipos+1);
			fn_tmp.decompose(my_next_number, fn_next_stack);
		}
		else
			fn_next_stack = "";

		numbers_in_stack.push_back(my_number - 1); // start counting at 0!

		// If next stack is a different one, read the current stack and process all images in it
		if (fn_next_stack != fn_my_stack)
		{

			Image<RFLOAT> stack, img;
			fImageHandler hFile;
			if (do_read_whole_stacks)
				// Read the entire stack into memory
				stack.read(fn_my_stack);
			else
				// Open the stack file
				hFile.openFile(fn_my_stack);

			// 1. Process the current stack
			for (long int inum = 0; inum < numbers_in_stack.size(); inum++)
			{
				// Get the image we want from the stack
				if (do_read_whole_stacks)
					stack().getImage(numbers_in_stack[inum], img());
				else
					img.readFromOpenFile(fn_my_stack, hFile, numbers_in_stack[inum]);
				long int my_ipos = my_stack_first_ipos + inum;

				bool have_optics_group = false;
				RFLOAT angpix = 0.;

				if (_do_apply_orient || lowpass > 0. || highpass > 0.)
				{
					if (MDin.containsLabel(EMDL_IMAGE_OPTICS_GROUP))
					{
						int optics_group;
						MDin.getValue(EMDL_IMAGE_OPTICS_GROUP, optics_group, my_ipos);
						optics_group--;
						obsModel->opticsMdt.getValue(EMDL_IMAGE_PIXEL_SIZE, angpix, optics_group);
						have_optics_group = true;
					}
				}

				if (_do_apply_orient && have_optics_group)
				{
					RFLOAT psi,rot,tilt;
					Matrix1D<RFLOAT> offset(3);
					Matrix2D<RFLOAT> A;
					MDin.getValue(EMDL_ORIENT_PSI, psi, my_ipos);
					MDin.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, XX(offset), my_ipos);
					MDin.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, YY(offset), my_ipos);
					if(img().getDim()==2)
					{
						offset /= angpix;
						rotation2DMatrix(psi, A);
						MAT_ELEM(A, 0, 2) = COSD(psi) * XX(offset) - SIND(psi) * YY(offset);
						MAT_ELEM(A, 1, 2) = COSD(psi) * YY(offset) + SIND(psi) * XX(offset);
						selfApplyGeometry(img(), A, IS_NOT_INV, DONT_WRAP);
					}
					else
					{
						MDin.getValue(EMDL_ORIENT_ROT, rot, my_ipos);
						MDin.getValue(EMDL_ORIENT_TILT, tilt, my_ipos);
						MDin.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, ZZ(offset), my_ipos);
						offset /= angpix;
						Euler_rotation3DMatrix(rot,tilt,psi, A);
						MAT_ELEM(A, 0, 3) = MAT_ELEM(A, 0, 0) * XX(offset) + MAT_ELEM(A, 0, 1) * YY(offset) + MAT_ELEM(A, 0, 2) * ZZ(offset);
						MAT_ELEM(A, 1, 3) = MAT_ELEM(A, 1, 0) * XX(offset) + MAT_ELEM(A, 1, 1) * YY(offset) + MAT_ELEM(A, 1, 2) * ZZ(offset);
						MAT_ELEM(A, 2, 3) = MAT_ELEM(A, 2, 0) * XX(offset) + MAT_ELEM(A, 2, 1) * YY(offset) + MAT_ELEM(A, 2, 2) * ZZ(offset);
						selfApplyGeometry(img(), A, IS_NOT_INV, DONT_WRAP);
					}
				}
				else if(_do_apply_orient && MDin.containsLabel(EMDL_MLMODEL_IS_HELIX) && img().getDim()==3)
				{
					RFLOAT psi,rot,tilt;
					Matrix2D<RFLOAT> A;
					Euler_rotation3DMatrix(0,90,0, A);
					MAT_ELEM(A, 0, 3) = MAT_ELEM(A, 0, 0)  + MAT_ELEM(A, 0, 1)  + MAT_ELEM(A, 0, 2) ;
					MAT_ELEM(A, 1, 3) = MAT_ELEM(A, 1, 0)  + MAT_ELEM(A, 1, 1)  + MAT_ELEM(A, 1, 2) ;
					MAT_ELEM(A, 2, 3) = MAT_ELEM(A, 2, 0)  + MAT_ELEM(A, 2, 1)  + MAT_ELEM(A, 2, 2) ;
					selfApplyGeometry(img(), A, IS_NOT_INV, DONT_WRAP);
				}
				if (_do_recenter)
				{
					selfTranslateCenterOfMassToCenter(img());
				}

				if (lowpass > 0. && have_optics_group)
					lowPassFilterMap(img(), lowpass, angpix);
				if (highpass > 0. && have_optics_group)
					highPassFilterMap(img(), highpass, angpix);

				// Dont change the user-provided _minval and _maxval in the getImageContrast routine!
				RFLOAT myminval = _minval;
				RFLOAT mymaxval = _maxval;
				getImageContrast(img(), myminval, mymaxval, _sigma_contrast);

				long int my_sorted_ipos = my_ipos;
				if (MDin.containsLabel(EMDL_SORTED_IDX))
				{
					// First get the sorted index
					MDin.getValue(EMDL_SORTED_IDX, my_sorted_ipos, my_ipos);
					// Then set the original index in the sorted index, so that particles can be written out in the correct order
					MDin.setValue(EMDL_SORTED_IDX, my_ipos, my_ipos);
				}
				icol = my_sorted_ipos % ncol;
				irow = my_sorted_ipos / ncol;
				nrow = XMIPP_MAX(nrow, irow+1);
				if (my_ipos == 0)
				{
					xsize_box = CEIL(_scale * XSIZE(img())) + 2 * xoff; // 2 pixels on each side in between all images
					ysize_box = CEIL(_scale * YSIZE(img())) + 2 * yoff;
				}
				int ycoor = irow * ysize_box;
				int xcoor = icol * xsize_box;

				DisplayBox* my_box = new DisplayBox(xcoor, ycoor, xsize_box, ysize_box, "");
				my_box->setData(img(), MDin.getObject(my_ipos), my_ipos, myminval, mymaxval, _scale, false);
				if (MDin.containsLabel(text_label))
				{
					MDin.getValueToString(text_label, my_box->img_label, my_ipos);
				}
				my_box->redraw();
				boxes[my_sorted_ipos] = my_box;//boxes.push_back(my_box);
			}

			// 2. Reset numbers_in_stack and my_stack_first_ipos for next stack
			numbers_in_stack.clear();
			my_stack_first_ipos = ipos + 1;
		}

		ipos++;

		if (ipos >= number_of_images)
			break;

		if (nr_imgs > 1 && ipos % barstep == 0)
			progress_bar(ipos);
	}

	if (nr_imgs > 1)
		progress_bar(nr_imgs);
}

void basisViewerCanvas::fill(MultidimArray<RFLOAT> &image, RFLOAT _minval, RFLOAT _maxval, RFLOAT _sigma_contrast, RFLOAT _scale)
{
	xoff = yoff = 0;
	nrow = ncol = 1;
	getImageContrast(image, _minval, _maxval, _sigma_contrast);
	xsize_box = CEIL(_scale * XSIZE(image));
	ysize_box = CEIL(_scale * YSIZE(image));
	DisplayBox* my_box = new DisplayBox(0, 0, xsize_box, ysize_box, "dummy");
	MetaDataTable MDtmp;
	MDtmp.addObject();
	//FileName fn_tmp = "dummy";
	//MDtmp.setValue(EMDL_IMAGE_NAME, fn_tmp);
	my_box->setData(image, MDtmp.getObject(), 0, _minval, _maxval, _scale, true);
	my_box->redraw();
	boxes.push_back(my_box);
}

void basisViewerCanvas::draw()
{
	for (int ipos = 0 ; ipos < boxes.size(); ipos++)
		boxes[ipos]->redraw();
}

int multiViewerCanvas::handle(int ev)
{
	if (ev==FL_PUSH)
	{
		int xc = (int)Fl::event_x() - scroll->x() + scroll->hscrollbar.value();
		int yc = (int)Fl::event_y() - scroll->y() + scroll->scrollbar.value();
		int xpos = xc / xsize_box;
		int ypos = yc / ysize_box;
		int ipos = ypos * ncol + xpos;
		// Check there was no click in the area outside the boxes...
		if (xpos < ncol && ypos < nrow && ipos < boxes.size())
		{
			if (Fl::event_button() == FL_LEFT_MOUSE)
			{
				// Shift-left-click will select a whole range
				if (Fl::event_state(FL_SHIFT))
				{
					if (has_shift)
					{
						int postshift_ipos = ipos;
						int ipos0 = (postshift_ipos > preshift_ipos) ? preshift_ipos : postshift_ipos;
						int iposF = (postshift_ipos > preshift_ipos) ? postshift_ipos : preshift_ipos;
						// Select all images from ipos0 to iposF
						// TODO!!! Cannot do this here: have to define an event for the multiview window as a whole!
						// This multiview window should have all the DisplayBoxes inside it....
						for (int my_ipos = ipos0; my_ipos <= iposF; my_ipos++)
						{
							boxes[my_ipos]->select();
						}
						has_shift = false;
					}
					else
					{
						preshift_ipos = ipos;
						has_shift = true;
					}

				}
				else
				{
					boxes[ipos]->toggleSelect(current_selection_type);
				}
			}
			else  if ( Fl::event_button() == FL_RIGHT_MOUSE )
			{
				Fl_Menu_Item rclick_menu;
				if (do_class)
				{

					Fl_Menu_Item rclick_menu[] = {
						{ "Save backup selection" },
						{ "Load backup selection" },
						{ "Clear selection" },
						{ "Invert selection" },
						{ "Select all classes below" },
						{ "Select all classes above" },
						{ "Show metadata this class" },
						{ "Show original image" },
						{ "Save image as PNG" },
						{ "Show Fourier amplitudes (2x)" },
						{ "Show Fourier phase angles (2x)" },
						{ "Show helical layer line profile" },
						{ "Show particles from selected classes" },
						{ "Set selection type" },
						{ "Save selected classes" }, // idx = 14; change below when re-ordered!!
						{ "Quit" },
						{ 0 }
					};

					if (!do_allow_save)
					{
						rclick_menu[14].deactivate();
					}

				    const Fl_Menu_Item *m = rclick_menu->popup(Fl::event_x(), Fl::event_y(), 0, 0, 0);
					if ( !m )
						return 0;
					else if ( strcmp(m->label(), "Save backup selection") == 0 )
						saveBackupSelection();
					else if ( strcmp(m->label(), "Load backup selection") == 0 )
						loadBackupSelection();
					else if ( strcmp(m->label(), "Clear selection") == 0 )
						clearSelection();
					else if ( strcmp(m->label(), "Invert selection") == 0 )
						invertSelection();
					else if ( strcmp(m->label(), "Select all classes below") == 0 )
						selectFromHereBelow(ipos);
					else if ( strcmp(m->label(), "Select all classes above") == 0 )
						selectFromHereAbove(ipos);
					else if ( strcmp(m->label(), "Show metadata this class") == 0 )
						printMetaData(ipos);
					else if ( strcmp(m->label(), "Show original image") == 0 )
						showOriginalImage(ipos);
					else if ( strcmp(m->label(), "Save image as PNG") == 0 )
						saveImage(ipos);
					else if ( strcmp(m->label(), "Show Fourier amplitudes (2x)") == 0 )
						showFourierAmplitudes(ipos);
					else if ( strcmp(m->label(), "Show Fourier phase angles (2x)") == 0 )
						showFourierPhaseAngles(ipos);
					else if ( strcmp(m->label(), "Show helical layer line profile") == 0 )
						showHelicalLayerLineProfile(ipos);
					else if ( strcmp(m->label(), "Set selection type") == 0 )
						setSelectionType();
					else if ( strcmp(m->label(), "Show particles from selected classes") == 0 )
						showSelectedParticles(current_selection_type);
					else if ( strcmp(m->label(), "Save selected classes") == 0 )
					{
						saveBackupSelection();
						saveSelected(current_selection_type);
						saveSelectedParticles(current_selection_type);
						// save the exit_success file after saving already,
						// as many users close the window through the operating system's cross symbol on the window, instead of a proper exit
						RELION_EXIT_SUCCESS;
					}
					else if ( strcmp(m->label(), "Quit") == 0 )
					{
						//clean exit
						exit(RELION_EXIT_SUCCESS);
					}
				}
				else
				{
					Fl_Menu_Item rclick_menu[] = {
						{ "Save backup selection" },
						{ "Load backup selection" },
						{ "Clear selection" },
						{ "Invert selection" },
						{ "Select all below" },
						{ "Select all above" },
						{ "Show average of selection" },
						{ "Show stddev of selection" },
						{ "Show original image" },
						{ "Save image as PNG" },
						{ "Show Fourier amplitudes (2x)" },
						{ "Show Fourier phase angles (2x)" },
						{ "Show helical layer line profile" },
						{ "Set selection type" },
						{ "Show metadata" },
						{ "Save STAR with selected images" }, // idx = 15; change below when re-ordered!!
						{ "Quit" },
						{ 0 }
					};
					if (!do_allow_save)
					{
						rclick_menu[15].deactivate();
					}

					const Fl_Menu_Item *m = rclick_menu->popup(Fl::event_x(), Fl::event_y(), 0, 0, 0);
					if ( !m )
						return 0;
					else if ( strcmp(m->label(), "Save backup selection") == 0 )
						saveBackupSelection();
					else if ( strcmp(m->label(), "Load backup selection") == 0 )
						loadBackupSelection();
					else if ( strcmp(m->label(), "Clear selection") == 0 )
						clearSelection();
					else if ( strcmp(m->label(), "Invert selection") == 0 )
						invertSelection();
					else if ( strcmp(m->label(), "Select all below") == 0 )
						selectFromHereBelow(ipos);
					else if ( strcmp(m->label(), "Select all above") == 0 )
						selectFromHereAbove(ipos);
					else if ( strcmp(m->label(), "Show average of selection") == 0 )
						showAverage(SELECTED, false);
					else if ( strcmp(m->label(), "Show stddev of selection") == 0 )
						showAverage(SELECTED, true);
					else if ( strcmp(m->label(), "Show original image") == 0 )
						showOriginalImage(ipos);
					else if ( strcmp(m->label(), "Save image as PNG") == 0 )
						saveImage(ipos);
					else if ( strcmp(m->label(), "Show Fourier amplitudes (2x)") == 0 )
						showFourierAmplitudes(ipos);
					else if ( strcmp(m->label(), "Show Fourier phase angles (2x)") == 0 )
						showFourierPhaseAngles(ipos);
					else if ( strcmp(m->label(), "Show helical layer line profile") == 0 )
						showHelicalLayerLineProfile(ipos);
					else if ( strcmp(m->label(), "Set selection type") == 0 )
						setSelectionType();
					else if ( strcmp(m->label(), "Show metadata") == 0 )
						printMetaData(ipos);
					else if ( strcmp(m->label(), "Save STAR with selected images") == 0 )
					{
						saveBackupSelection();
						saveSelected(SELECTED);
						// save the exit_success file after saving already,
						// as many users close the window through the operating system's cross symbol on the window, instead of a proper exit
						RELION_EXIT_SUCCESS;
					}

					else if ( strcmp(m->label(), "Quit") == 0 )
						exit(0);

				}
				return(1);          // (tells caller we handled this event)
			}
		} // endif ipos within valid region

	}
	return 0;
}

void multiViewerCanvas::saveBackupSelection()
{
	std::vector<int> selected(boxes.size());
	for (long int ipos = 0; ipos < boxes.size(); ipos++)
	{
		long int my_sorted_ipos;
		if (boxes[ipos]->MDimg.containsLabel(EMDL_SORTED_IDX))
			boxes[ipos]->MDimg.getValue(EMDL_SORTED_IDX, my_sorted_ipos);
		else
			my_sorted_ipos = ipos;
		selected[my_sorted_ipos] = boxes[ipos]->selected;
	}

	for (long int ipos = 0; ipos < boxes.size(); ipos++)
	{
        if (MDbackup.numberOfObjects() < ipos+1)
        	MDbackup.addObject();
		// without the bool() cast, clang will interpret the formal template parameter
        // as a reference to a bit field, which is not the same as a boolean.
		MDbackup.setValue(EMDL_SELECTED, selected[ipos], ipos);
	}

	FileName fn_dir;
	if (fn_selected_imgs != "")
		fn_dir = fn_selected_imgs.beforeLastOf("/");
	else if (fn_selected_parts != "")
		fn_dir = fn_selected_parts.beforeLastOf("/");
	else
		fn_dir = ".";
	fn_dir += "/backup_selection.star";

	MDbackup.write(fn_dir);
	std::cout <<" Written out " << fn_dir << std::endl;
}

void multiViewerCanvas::loadBackupSelection(bool do_ask)
{

	FileName fn_sel, fn_dir;
	if (fn_selected_imgs != "")
		fn_dir = fn_selected_imgs.beforeLastOf("/");
	else if (fn_selected_parts != "")
		fn_dir = fn_selected_parts.beforeLastOf("/");
	else
		fn_dir = ".";
	fn_dir += "/";
	if (do_ask)
	{
		Fl_File_Chooser chooser(fn_dir.c_str(), "(backup_selection.star)",Fl_File_Chooser::SINGLE,"Choose selection file to load");      // chooser type
		chooser.show();
		// Block until user picks something.
		while(chooser.shown())
			{ Fl::wait(); }

		// User hit cancel?
		if ( chooser.value() == NULL )
			return;

		FileName fnt(chooser.value());
		fn_sel = fnt;
	}
	else
		fn_sel = fn_dir+"backup_selection.star";

	MDbackup.clear();
	MDbackup.read(fn_sel);
	if (MDbackup.numberOfObjects() != boxes.size())
	{
		std::cerr << "Warning: ignoring .relion_display_backup_selection.star with unexpected number of entries..." << std::endl;
		return;
	}

	std::vector<int> selected(boxes.size(), false);
	long int ipos = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDbackup)
	{
		MDbackup.getValue(EMDL_SELECTED, selected[ipos]);
		ipos++;
	}

	for (long int ipos = 0; ipos < boxes.size(); ipos++)
	{
		long int my_sorted_ipos;
		if (boxes[ipos]->MDimg.containsLabel(EMDL_SORTED_IDX))
			boxes[ipos]->MDimg.getValue(EMDL_SORTED_IDX, my_sorted_ipos);
		else
			my_sorted_ipos = ipos;

		boxes[ipos]->setSelect(selected[my_sorted_ipos]);
	}
}

void multiViewerCanvas::clearSelection()
{
	for (long int ipos = 0; ipos < boxes.size(); ipos++)
	{
		boxes[ipos]->unSelect();
	}
}

void multiViewerCanvas::invertSelection()
{
	for (long int ipos = 0; ipos < boxes.size(); ipos++)
	{
		boxes[ipos]->toggleSelect(current_selection_type);
	}
}

void multiViewerCanvas::selectFromHereBelow(int iposp)
{
	for (long int ipos = iposp; ipos < boxes.size(); ipos++)
	{
		boxes[ipos]->select();
	}
}

void multiViewerCanvas::selectFromHereAbove(int iposp)
{
	for (long int ipos = 0; ipos <= iposp; ipos++)
	{
		boxes[ipos]->select();
	}

}
void multiViewerCanvas::printMetaData(int main_ipos)
{
	std::ostringstream stream;

	if (do_class) {
		int myclass, iclass, nselected_classes = 0, nselected_particles = 0;
		for (long int ipos = 0; ipos < boxes.size(); ipos++)
		{
			if (boxes[ipos]->selected == SELECTED)
			{
				nselected_classes++;
				// Get class number (may not be ipos+1 if resorted!)
				boxes[ipos]->MDimg.getValue(EMDL_PARTICLE_CLASS, myclass);
				FOR_ALL_OBJECTS_IN_METADATA_TABLE(*MDdata)
				{
					MDdata->getValue(EMDL_PARTICLE_CLASS, iclass);
					if (iclass == myclass) nselected_particles++;
				}
			}
		}
		stream << "Selected " << nselected_particles << " particles in " << nselected_classes << " classes.\n";
	}
	stream << "Below is the metadata table for the last clicked class/particle.\n";

	boxes[main_ipos]->MDimg.write(stream);
	FileName str =  stream.str();

	// @ starts special symbol code in FLTK; we must escape it
	size_t pos = str.find('@', 0);
	while (pos != std::string::npos)
	{
		str.replace(pos, 1, (std::string)"@@" );
		pos = str.find('@', pos + 2);
	}
	fl_message("%s",str.c_str());
}

void multiViewerCanvas::showAverage(bool selected, bool show_stddev)
{
	int xsize = boxes[0]->xsize_data;
	int ysize = boxes[0]->ysize_data;
	MultidimArray<RFLOAT> sum(ysize, xsize);
	MultidimArray<RFLOAT> sum2(ysize, xsize);

	int nn = 0;
	for (long int ipos = 0; ipos < boxes.size(); ipos++)
	{
		if (boxes[ipos]->selected == selected)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sum)
			{
				int ival = boxes[ipos]->img_data[n];
				if (ival < 0) ival += 256;
				DIRECT_MULTIDIM_ELEM(sum, n) += ival;
				DIRECT_MULTIDIM_ELEM(sum2, n) += ival * ival;
			}
			nn++;
		}
	}
	sum  /= nn;
	sum2 /= nn;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sum)
	{
		DIRECT_MULTIDIM_ELEM(sum2, n) -=  DIRECT_MULTIDIM_ELEM(sum, n) * DIRECT_MULTIDIM_ELEM(sum, n);
	}
	sum2 *= nn / (nn - 1);

	// Show the average
	if (show_stddev)
	{
		basisViewerWindow stddev(xsize, ysize, "Stddev");
		stddev.fillSingleViewerCanvas(sum2, 0., 0., 0., 1.); // scale=1 now means: keep same scale as the one in the boxes!!!
	}
	else
	{
		basisViewerWindow avg(xsize, ysize, "Average");
		avg.fillSingleViewerCanvas(sum, 0., 0., 0., 1.); // scale=1 now means: keep same scale as the one in the boxes!!!
	}

}

void multiViewerCanvas::showOriginalImage(int ipos)
{
	// Make system call because otherwise the green drawing for distance measurements doesn't work....
	FileName fn_img;
	boxes[ipos]->MDimg.getValue(display_label, fn_img);

	std::string cl = "relion_display  --i " + fn_img + " --scale " + floatToString(ori_scale);
	cl += " --sigma_contrast " + floatToString(sigma_contrast);
	cl += " --black " + floatToString(minval);
	cl += " --white " + floatToString(maxval);

	switch (colour_scheme)
	{
	case (BLACKGREYREDSCALE): { cl += " --colour_fire"; break; }
	case (BLUEGREYWHITESCALE): { cl += " --colour_ice"; break; }
	case (BLUEGREYREDSCALE): { cl += " --colour_fire-n-ice"; break; }
	case (RAINBOWSCALE): { cl += " --colour_rainbow"; break; }
	case (CYANBLACKYELLOWSCALE): { cl += " --colour_difference"; break; }
	}
	// send job in the background
	cl += " &";

	int res = system(cl.c_str());

	/*
	FileName fn_img;
	boxes[ipos]->MDimg.getValue(display_label, fn_img);
	Image<RFLOAT> img;
	img.read(fn_img);
	basisViewerWindow win(CEIL(ori_scale*XSIZE(img())), CEIL(ori_scale*YSIZE(img())), fn_img.c_str());
	if (sigma_contrast > 0.)
	{
		win.fillSingleViewerCanvas(img(), 0., 0., sigma_contrast, ori_scale);
	}
	else
	{
		win.fillSingleViewerCanvas(img(), boxes[ipos]->minval, boxes[ipos]->maxval, 0., ori_scale);
	}
    */
}

void basisViewerCanvas::saveImage(int ipos)
{
	#ifndef HAVE_PNG
		fl_message("Cannot save an image as PNG because libPNG was not linked during compilation.");
	#else
		using namespace gravis;

		Fl_File_Chooser chooser(".",                        // directory
		                        "PNG image (*.png)\tAll Files (*)*", // filter
		                        Fl_File_Chooser::CREATE, // chooser type
		                        "Save as"); // title
		chooser.show();
		// Block until user picks something.
		while(chooser.shown())
			{ Fl::wait(); }
		// User hit cancel?
		if ( chooser.value() == NULL )
			return;

		int xsize = boxes[ipos]->xsize_data;
		int ysize = boxes[ipos]->ysize_data;
		unsigned char* img_data = boxes[ipos]->img_data;

		tImage<bRGB> pngOut(xsize, ysize);
                pngOut.fill(bRGB(0));

		for (size_t n = 0, nlim = xsize * ysize; n < nlim; n++)
		{
			if (colour_scheme == GREYSCALE)
			{
				unsigned char c = img_data[n];
				pngOut[n] = bRGB(c, c, c);
			}
			else
			{
				pngOut[n] = bRGB(img_data[3 * n], img_data[3 * n + 1], img_data[3 * n + 2]);
			}
		}

		pngOut.writePNG(chooser.value());
	#endif
}

void multiViewerCanvas::showFourierAmplitudes(int ipos)
{
	// Make system call because otherwise the green drawing for distance measurements doesn't work....
	FileName fn_img;
	Image<RFLOAT> img;
	boxes[ipos]->MDimg.getValue(display_label, fn_img);
	img.read(fn_img, false);
	if ( (ZSIZE(img()) > 1) || (NSIZE(img()) > 1) )
	{
		 fl_message("Cannot display Fourier transform of STAR files, 3D images or stacks. Please select a 2D image as input.");
		 return;
	}

	std::string cl = "relion_display  --i " + fn_img + " --scale " + floatToString(ori_scale);
	if (sigma_contrast > 0.)
    		cl += " --sigma_contrast " + floatToString(sigma_contrast);
	cl += " --show_fourier_amplitudes";
	// send job in the background
	cl += " &";

	int res = system(cl.c_str());
}

void multiViewerCanvas::showFourierPhaseAngles(int ipos)
{
	// Make system call because otherwise the green drawing for distance measurements doesn't work....
	FileName fn_img;
	Image<RFLOAT> img;
	boxes[ipos]->MDimg.getValue(display_label, fn_img);
	img.read(fn_img, false);
	if ( (ZSIZE(img()) > 1) || (NSIZE(img()) > 1) )
	{
		fl_message("Cannot display Fourier transform of STAR files, 3D images or stacks. Please select a 2D image as input.");
		return;
	}

	std::string cl = "relion_display  --i " + fn_img + " --scale " + floatToString(ori_scale);
	cl += " --show_fourier_phase_angles";
	// send job in the background
	cl += " &";

	int res = system(cl.c_str());
}

void multiViewerCanvas::showHelicalLayerLineProfile(int ipos)
{
	const char * default_pdf_viewer = getenv ("RELION_PDFVIEWER_EXECUTABLE");
	char hardcoded_pdf_viewer[]=DEFAULTPDFVIEWER;
	if (default_pdf_viewer == NULL)
	{
		default_pdf_viewer=hardcoded_pdf_viewer;
	}

	std::string mydefault = std::string(default_pdf_viewer);
	std::string command;
	FileName fn_img, fn_out;
	Image<RFLOAT> img;

	boxes[ipos]->MDimg.getValue(display_label, fn_img);
	img.read(fn_img);

	fn_out = "layerlineprofile.eps";
	if (exists(fn_out))
	{
		command = "rm -rf " + fn_out;
		int res = system(command.c_str());
	}

	helicalLayerLineProfile(img(), fn_img, fn_out);

	command = mydefault + " " + fn_out + " &";
	int res = system(command.c_str());
}

void multiViewerCanvas::makeStarFileSelectedParticles(int selected, MetaDataTable &MDpart)
{
	MDpart.clear();
	int myclass, iclass;
	for (long int ipos = 0; ipos < boxes.size(); ipos++)
	{
		if (boxes[ipos]->selected == selected)
		{
			// Get class number (may not be ipos+1 if resorted!)
			boxes[ipos]->MDimg.getValue(EMDL_PARTICLE_CLASS, myclass);
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(*MDdata)
			{
				MDdata->getValue(EMDL_PARTICLE_CLASS, iclass);
				if (iclass == myclass)
					MDpart.addObject(MDdata->getObject());
			}
		}
	}

	if (max_nr_parts_per_class > 0)
	{
		// Randomise the order, to pick random particles from each class
		// Unfortunately, this leads to random order of particles in the output file! So be it for now...
		MetaDataTable MDtmp = MDpart;
		MDpart.clear();
		MDtmp.sort(EMDL_UNDEFINED, false, false, true);
		for (long int ipos = 0; ipos < boxes.size(); ipos++)
		{
			if (boxes[ipos]->selected == selected)
			{
				int nr_selected = 0;
				boxes[ipos]->MDimg.getValue(EMDL_PARTICLE_CLASS, myclass);
				FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDtmp)
				{
					MDtmp.getValue(EMDL_PARTICLE_CLASS, iclass);
					if (iclass == myclass)
					{
						MDpart.addObject(MDtmp.getObject());
						nr_selected++;
						if (nr_selected >= max_nr_parts_per_class)
							break;
					}
				}
			}
		}
	}

	// Maintain the original image ordering
	if (MDpart.containsLabel(EMDL_SORTED_IDX))
		MDpart.sort(EMDL_SORTED_IDX);

}

void multiViewerCanvas::saveSelectedParticles(int save_selected)
{
	if (fn_selected_parts == "")
	{
		std::cout << " Not saving selected particles, as no filename was provided..." << std::endl;
		return;
	}

//#define RELION_DEVEL_ASKTRAINING
#ifdef RELION_DEVEL_ASKTRAINING
	bool do_training = false;
	std::string ask = "Is this a selection of good classes, so it can be used for Sjors' training set for automated class selection?\n \
			More info here: /public/EM/RELION/training.txt\n";
	do_training =  fl_choice("%s", "Don't use", "Use for training", NULL, ask.c_str());

	if (do_training)
		saveTrainingSet();
#endif

	MetaDataTable MDpart;
	makeStarFileSelectedParticles(save_selected, MDpart);
	if (nr_regroups > 0)
		regroupSelectedParticles(MDpart, *MDgroups, nr_regroups);
	int nparts = MDpart.numberOfObjects();
	if (nparts > 0)
	{
		obsModel->save(MDpart, fn_selected_parts, "particles");
		std::cout << "Saved " << fn_selected_parts << " with " << nparts << " selected particles." << std::endl;
	}
	else
		std::cout <<" No classes selected. Please select one or more classes..." << std::endl;
}

void regroupSelectedParticles(MetaDataTable &MDdata, MetaDataTable &MDgroups, int nr_regroups)
{
	// This function modify MDgroups, which will not be written anyway.

	if (nr_regroups <= 0)
		return;

	int max_optics_group_id = -1;

	// Find out which optics group each scale group belongs to
	// Also initialise rlnGroupNrParticles for this selection
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDgroups)
	{
		MDgroups.setValue(EMDL_IMAGE_OPTICS_GROUP, -1);
		MDgroups.setValue(EMDL_MLMODEL_GROUP_NR_PARTICLES, 0);
	}

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDdata)
	{
		long group_id, part_optics_id, group_optics_id;
		int nr_parts;
		MDdata.getValue(EMDL_MLMODEL_GROUP_NO, group_id); // 1-indexed
		MDdata.getValue(EMDL_IMAGE_OPTICS_GROUP, part_optics_id);

		MDgroups.getValue(EMDL_IMAGE_OPTICS_GROUP, group_optics_id, group_id - 1); // 0-indexed
		if (group_optics_id == -1)
		{
			MDgroups.setValue(EMDL_IMAGE_OPTICS_GROUP, part_optics_id, group_id - 1);
			if (max_optics_group_id < part_optics_id)
				max_optics_group_id = part_optics_id;
		}
		else if (group_optics_id != part_optics_id)
		{
			std::cerr << "WARNING: group_no " << group_id << " contains particles from multiple optics groups." << std::endl;
		}

		MDgroups.getValue(EMDL_MLMODEL_GROUP_NR_PARTICLES, nr_parts, group_id - 1);
		MDgroups.setValue(EMDL_MLMODEL_GROUP_NR_PARTICLES, nr_parts + 1, group_id - 1);
	}

	// First sort the MDgroups based on refined intensity scale factor
	MDgroups.sort(EMDL_MLMODEL_GROUP_SCALE_CORRECTION);

	// Store original image order
	long int nr_parts = MDdata.numberOfObjects();
	for (long int j = 0; j < nr_parts; j++)
		MDdata.setValue(EMDL_SORTED_IDX, j, j);

	// Average group size
	long average_group_size = nr_parts / nr_regroups;
	if (average_group_size < 10)
		REPORT_ERROR("Each group should have at least 10 particles");
	int fillgroupschar = (int)(floor(log(nr_regroups) / log(10))) + 1;

	std::map<long, std::string> new_group_names;
	std::map<long, std::string>::iterator it;

	// Loop through all existing, sorted groups
	long new_group_id = 0;

	// Worst case: O(old_nr_groups ^ 2) = O(mic ^ 2)
	// We can reduce this by using one more hash but this should be enough.
	for (long optics_group_id = 1; optics_group_id <= max_optics_group_id; optics_group_id++)
	{
		long nr_parts_in_new_group = 0;
		new_group_id++;

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDgroups)
		{
			long group_id, group_optics_id;
			int nr_parts;

			MDgroups.getValue(EMDL_IMAGE_OPTICS_GROUP, group_optics_id);
			if (group_optics_id != optics_group_id)
				continue;

			MDgroups.getValue(EMDL_MLMODEL_GROUP_NO, group_id);
			MDgroups.getValue(EMDL_MLMODEL_GROUP_NR_PARTICLES, nr_parts);
			nr_parts_in_new_group += nr_parts;

			if (nr_parts_in_new_group > average_group_size)
			{
				// This group is now full: start a new one
				new_group_id++;
				nr_parts_in_new_group = 0;
			}

			new_group_names[group_id] = "group_" + integerToString(new_group_id, fillgroupschar);
		}
	}

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDdata)
	{
		long group_id;
		MDdata.getValue(EMDL_MLMODEL_GROUP_NO, group_id);

		it = new_group_names.find(group_id);
		if (it != new_group_names.end())
		{
			MDdata.setValue(EMDL_MLMODEL_GROUP_NAME, new_group_names[group_id]);
		}
		else
		{
			std::cerr << "Logic error: cannot find group_id " << group_id << " during remapping." << std::endl;
			REPORT_ERROR("Failed in regrouping");
		}
	}

	MDdata.deactivateLabel(EMDL_MLMODEL_GROUP_NO); // no longer valid

	std::cout <<" Regrouped particles into " << new_group_id << " groups" << std::endl;
}

void multiViewerCanvas::showSelectedParticles(int save_selected)
{
	MetaDataTable MDpart;
	makeStarFileSelectedParticles(save_selected, MDpart);
	int nparts = MDpart.numberOfObjects();
	if (nparts > 0)
	{
		basisViewerWindow win(MULTIVIEW_WINDOW_WIDTH, MULTIVIEW_WINDOW_HEIGHT, "Particles in the selected classes");
		win.fillCanvas(MULTIVIEWER, MDpart, obsModel, EMDL_IMAGE_NAME, text_label, do_read_whole_stacks, do_apply_orient, 0., 0., 0., boxes[0]->scale, ori_scale, ncol, multi_max_nr_images);
	}
	else
		std::cout <<" No classes selected. First select one or more classes..." << std::endl;
}

void multiViewerCanvas::saveTrainingSet()
{
	FileName fn_rootdir = "/net/dstore1/teraraid3/scheres/trainingset/";

	// Make the output job directory
	char my_dir[200];
	FileName fn_projdir = std::string(getcwd(my_dir, 200));
	std::replace( fn_projdir.begin(), fn_projdir.end(), '/', '_');
	fn_projdir += "/" + (fn_selected_parts.afterFirstOf("/")).beforeLastOf("/");
	FileName fn_odir = fn_rootdir + fn_projdir;
	std::string command = "mkdir -p " + fn_odir + " ; chmod 777 " + fn_odir;
	int res = system(command.c_str());

	// Now save the selected images in a MetaData file.
	MetaDataTable MDout;
	int nsel = 0;
	for (long int ipos = 0; ipos < boxes.size(); ipos++)
	{
		MDout.addObject(boxes[ipos]->MDimg.getObject());
		if (boxes[ipos]->selected)
			MDout.setValue(EMDL_SELECTED, 1);
		else
			MDout.setValue(EMDL_SELECTED, 0);
	}

	// Maintain the original image ordering
	if (MDout.containsLabel(EMDL_SORTED_IDX))
		MDout.sort(EMDL_SORTED_IDX);

	// Copy all images
	long int nr;
	FileName fn_img, fn_new_img, fn_iroot, fn_old="";
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDout)
	{
		MDout.getValue(display_label, fn_img);
		fn_img.decompose(nr, fn_img);
		fn_new_img.compose(nr, fn_img.afterLastOf("/"));
		MDout.setValue(display_label, fn_new_img);
		if (fn_img != fn_old) // prevent multiple copies of single stack from Class2D
			copy(fn_img, fn_odir+"/"+fn_img.afterLastOf("/"));
		fn_old = fn_img;
	}
	fn_iroot = fn_img.beforeFirstOf("_class");

	// Copy rest of metadata
	fn_img = fn_iroot + "_model.star";
	copy(fn_img, fn_odir+"/"+fn_img.afterLastOf("/"));
	fn_img = fn_iroot + "_optimiser.star";
	copy(fn_img, fn_odir+"/"+fn_img.afterLastOf("/"));
	fn_img = fn_iroot + "_data.star";
	copy(fn_img, fn_odir+"/"+fn_img.afterLastOf("/"));
	fn_img = fn_iroot + "_sampling.star";
	copy(fn_img, fn_odir+"/"+fn_img.afterLastOf("/"));
	fn_iroot = fn_iroot.beforeLastOf("/");
	fn_img = fn_iroot + "/note.txt";
	copy(fn_img, fn_odir+"/"+fn_img.afterLastOf("/"));
        fn_img = fn_iroot + "/run_unmasked_classes.mrcs";
        if (exists(fn_img))
        {
            copy(fn_img, fn_odir+"/"+fn_img.afterLastOf("/"));
        }
	fn_img = fn_iroot + "/default_pipeline.star";
	copy(fn_img, fn_odir+"/"+fn_img.afterLastOf("/"));

	// Save the actual selection selection
	MDout.write(fn_odir + "/selected.star");

	// Give everyone permissions to this directory and its files
	//command = " chmod 777 -R " + fn_odir.beforeLastOf("/");
	//if (system(command.c_str()))
	//	REPORT_ERROR("ERROR in executing: " + command);

	std::cout << "Saved selection to Sjors' training directory. Thanks for helping out!" << std::endl;
}

void multiViewerCanvas::saveSelected(int save_selected)
{
	if (fn_selected_imgs == "")
		return;

	// Now save the selected images in a MetaData file.
	MetaDataTable MDout;
	int nsel = 0;
	for (long int ipos = 0; ipos < boxes.size(); ipos++)
	{
		if (boxes[ipos]->selected == save_selected)
		{
			nsel++;
			MDout.addObject(boxes[ipos]->MDimg.getObject());
		}
	}
	if (nsel > 0)
	{
		// Maintain the original image ordering
		if (MDout.containsLabel(EMDL_SORTED_IDX))
			MDout.sort(EMDL_SORTED_IDX);

		// If the images were re-centered to the center-of-mass, then output the recentered images, and change the names of the images in the MDout.
		if (do_recenter)
		{
			FileName fn_stack = fn_selected_imgs.withoutExtension()+".mrcs";
			FileName fn_img, fn_out;
			Image<RFLOAT> img;
			long int i = 0;
			long int nr_images = MDout.numberOfObjects();
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDout)
			{
				i++;
				MDout.getValue(EMDL_MLMODEL_REF_IMAGE, fn_img);
				img.read(fn_img);
				selfTranslateCenterOfMassToCenter(img());
				fn_out.compose(i, fn_stack);
				MDout.setValue(EMDL_MLMODEL_REF_IMAGE, fn_out);

                if (i == 1)
                	img.write(fn_stack, -1, (nr_images > 1), WRITE_OVERWRITE);
                else
                	img.write(fn_stack, -1, false, WRITE_APPEND);

			}
		}


		if (obsModel->opticsMdt.numberOfObjects() > 0 && !do_class)
		{
			if (metadata_table_name == "micrographs" || (!MDout.containsLabel(EMDL_IMAGE_NAME) && !MDout.containsLabel(EMDL_MICROGRAPH_MOVIE_NAME)))
			{
				obsModel->save(MDout, fn_selected_imgs, "micrographs");
				std::cout << "Saved "<< fn_selected_imgs << " with " << nsel << " selected micrographs." << std::endl;
			}
			else if (metadata_table_name == "movies" || (!MDout.containsLabel(EMDL_IMAGE_NAME) && MDout.containsLabel(EMDL_MICROGRAPH_MOVIE_NAME)))
			{
				obsModel->save(MDout, fn_selected_imgs, "movies");
				std::cout << "Saved "<< fn_selected_imgs << " with " << nsel << " selected movies." << std::endl;

			}
			else
			{
				obsModel->save(MDout, fn_selected_imgs, "particles");
				std::cout << "Saved "<< fn_selected_imgs << " with " << nsel << " selected particles." << std::endl;
			}
		}
		else
		{
			MDout.write(fn_selected_imgs);
			std::cout << "Saved "<< fn_selected_imgs << " with " << nsel << " selected images." << std::endl;
		}
	}
	else
		std::cout <<" No images to save...." << std::endl;
}

void basisViewerCanvas::setSelectionType()
{
	popupSelectionTypeWindow win(250, 50, "Set selection type");
	win.fill();
}

void basisViewerCanvas::setFOMThreshold()
{
	const char *pfom;
	std::string currentval = floatToString(minimum_pick_fom);
	pfom =  fl_input("Minimum rlnAutopickFigureOfMerit to display: ", currentval.c_str());
	if (pfom == NULL)
		return;
	std::string newval(pfom);
	minimum_pick_fom = textToFloat(newval);

}

int popupSelectionTypeWindow::fill()
{
	color(GUI_BACKGROUND_COLOR);
	choice = new Fl_Choice(50, 10, 130, 30, "type: ") ;

	choice->menu(color_choices);
	choice->color(GUI_INPUT_COLOR);

	choice->value(current_selection_type - 1);

	choice->callback(cb_set, this);

	Fl_Button * closebutton = new Fl_Button(190, 10, 50, 30, "Close");
	closebutton->color(GUI_RUNBUTTON_COLOR);
	closebutton->callback(cb_close, this);

	show();

	return Fl::run();
}

int singleViewerCanvas::handle(int ev)
{
	if (ev==FL_PUSH && Fl::event_button() == FL_LEFT_MOUSE)
	{
		int rx = (int)Fl::event_x() - scroll->x() + scroll->hscrollbar.value();
		int ry = (int)Fl::event_y() - scroll->y() + scroll->scrollbar.value();
		// Left mouse click writes value and coordinates to screen

		if (rx < boxes[0]->xsize_data && ry < boxes[0]->ysize_data && rx >= 0 && ry >=0)
		{
			unsigned char ival;
			int n = ry*boxes[0]->xsize_data + rx;
			if (colour_scheme == GREYSCALE) ival = boxes[0]->img_data[n];
			else ival = rgbToGrey(boxes[0]->img_data[3*n], boxes[0]->img_data[3*n+1], boxes[0]->img_data[3*n+2]);
			RFLOAT step = (boxes[0]->maxval - boxes[0]->minval) / 255.;
			RFLOAT dval = ival * step + boxes[0]->minval;
			int ysc = ROUND(ry/boxes[0]->scale);
			int xsc = ROUND(rx/boxes[0]->scale);
			int yscp = ysc - ROUND((boxes[0]->ysize_data/(2.* boxes[0]->scale)));
			int xscp = xsc - ROUND((boxes[0]->xsize_data/(2.* boxes[0]->scale)));
			std::cout <<" Image value at (" << xsc << "," << ysc << ") or (" << xscp << "," << yscp << ")~= " << dval <<std::endl;
		}
		return(1);
	}
	else if (ev==FL_PUSH && Fl::event_button() == FL_RIGHT_MOUSE)
	{
		Fl_Menu_Item rclick_menu[] = {
			{ "Show metadata" },
			{ "Save image as PNG" },
			{ "Help" },
			{ "Quit" },
			{ 0 }
		};
		const Fl_Menu_Item *m = rclick_menu->popup(Fl::event_x(), Fl::event_y(), 0, 0, 0);
		if ( !m )
			return 0;
		if ( strcmp(m->label(), "Show metadata") == 0 )
			printMetaData();
		else if ( strcmp(m->label(), "Save image as PNG") == 0 )
			saveImage();
		else if ( strcmp(m->label(), "Help") == 0 )
			printHelp();
		else if ( strcmp(m->label(), "Quit") == 0 )
			exit(0);
		return(1);          // (tells caller we handled this event)
	}
	else if (ev==FL_PUSH && Fl::event_button() == FL_MIDDLE_MOUSE)
	{
		// Middle-mouse dragging for measuring distances
		if (!has_dragged)
		{
			redraw();
			predrag_xc = (int)Fl::event_x();
			predrag_yc = (int)Fl::event_y();
			has_dragged = true;
			fl_color(FL_RED);
			fl_circle(predrag_xc, predrag_yc, 3);
		}
		return(1);
	}
	else if (ev==FL_DRAG  && Fl::event_button() == FL_MIDDLE_MOUSE)
	{
		fl_color(FL_RED);
		fl_circle(predrag_xc, predrag_yc, 3);
	}
	else if (ev==FL_RELEASE  && Fl::event_button() == FL_MIDDLE_MOUSE)
	{
		int postdrag_xc = (int)Fl::event_x();
		int postdrag_yc = (int)Fl::event_y();
		if (has_dragged)
		{
			fl_color(FL_RED);
			fl_circle(predrag_xc, predrag_yc, 3);
			fl_line(predrag_xc, predrag_yc, postdrag_xc, postdrag_yc);
			fl_circle(postdrag_xc, postdrag_yc, 3);
			int dx = postdrag_xc - predrag_xc;
			int dy = postdrag_yc - predrag_yc;
			RFLOAT dist = sqrt((RFLOAT)(dx*dx + dy*dy));
			std::string text =  floatToString(dist/boxes[0]->scale) + " pixels";
			fl_draw(text.c_str(), (postdrag_xc + predrag_xc)/2, (postdrag_yc + predrag_yc)/2);
			// Also write to the screen, in case the text falls outside the screen
			std::cout << "distance= " << dist/boxes[0]->scale << " pixels" << std::endl;
			has_dragged = false;
		}
		return(1);
	}

	return 0;
}

void singleViewerCanvas::printHelp()
{
	std::cout <<" + Left-mouse click: print coordinates and intensity value to screen " << std::endl;
	std::cout <<" + Middle-mouse drag: measure distances " << std::endl;
	std::cout <<" + Right-mouse click: pop-up menu" << std::endl;
}

/*
int popupSetContrastWindow::fill()
{
	color(GUI_BACKGROUND_COLOR);

	int width = 435;
	int x=150, y=15, ystep = 27, height = 25,  inputwidth = 50;
	int x2 = width - inputwidth - 50;

	// Always display these:
	scale = new Fl_Input(x, y, inputwidth, height, "Scale:");
	scale->color(GUI_INPUT_COLOR);
	scale->value("1");

	minval = new Fl_Input(x2, y, inputwidth, height, "Black value:");
	minval->value("0");
	minval->color(GUI_INPUT_COLOR);
	y += ystep;

	sigma_contrast = new Fl_Input(x, y, inputwidth, height, "Sigma contrast:");
	sigma_contrast->value("0");
	sigma_contrast->color(GUI_INPUT_COLOR);

	maxval = new Fl_Input(x2, y, inputwidth, height, "White value:");
	maxval->value("0");
	maxval->color(GUI_INPUT_COLOR);
	y += ROUND(ystep);

	Fl_Button * applybutton = new Fl_Button(width-120, y, 70, 30, "Apply!");
	applybutton->color(GUI_RUNBUTTON_COLOR);
	applybutton->callback( cb_set, this);

	Fl_Button * closebutton = new Fl_Button(width -200, y, 70, 30, "Close");
	closebutton->color(GUI_RUNBUTTON_COLOR);
	closebutton->callback( cb_close, this);

	show();

	return Fl::run();
}
*/

void pickerViewerCanvas::draw()
{
	RFLOAT scale = boxes[0]->scale;

	long int icoord = 0;
	int xcoori_start, ycoori_start;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDcoords)
	{
		icoord++;

		RFLOAT xcoor, ycoor;
		MDcoords.getValue(EMDL_IMAGE_COORD_X, xcoor);
		MDcoords.getValue(EMDL_IMAGE_COORD_Y, ycoor);

		if (MDcoords.containsLabel(EMDL_PARTICLE_AUTOPICK_FOM) && fabs(minimum_pick_fom + 9999.) > 1e-6)
		{
			RFLOAT fom;
			MDcoords.getValue(EMDL_PARTICLE_AUTOPICK_FOM, fom);
			if (fom < minimum_pick_fom) continue;
		}

		if (color_label != EMDL_UNDEFINED)
		{
			RFLOAT colval;
			if (EMDL::isInt(color_label))
			{
				int ival;
				if (!MDcoords.getValue(color_label, ival))
				{
					ival = 2; // populate as green if absent
					MDcoords.setValue(color_label, ival);
				}
				colval = (RFLOAT)ival;
				if (ival >= 1 && ival <= NUM_COLORS)
					fl_color(color_choices[ival - 1].labelcolor_);
				else
					fl_color(FL_GREEN);
			}
			else
			{
				MDcoords.getValue(color_label, colval);

				// Assume undefined values are set to -999....
				if ((colval + 999.) < XMIPP_EQUAL_ACCURACY)
				{
					fl_color(FL_GREEN);
				}
				else
				{
					colval = XMIPP_MAX(colval, smallest_color_value);
					colval = XMIPP_MIN(colval, biggest_color_value);
					unsigned char red, blue;
					if (do_blue_to_red)
					{
						red  = ROUND(255. * (colval - smallest_color_value) / (biggest_color_value - smallest_color_value));
						blue = ROUND(255. * (biggest_color_value - colval)  / (biggest_color_value - smallest_color_value));
					}
					else
					{
						blue = ROUND(255. * (colval - smallest_color_value) / (biggest_color_value - smallest_color_value));
						red  = ROUND(255. * (biggest_color_value - colval)  / (biggest_color_value - smallest_color_value));
					}
					fl_color(red, 0, blue);
				}
			}
		}
		else
		{
			fl_color(FL_GREEN);
		}

		int xcoori, ycoori;
		xcoori = ROUND(xcoor * coord_scale * scale) + scroll->x() - scroll->hscrollbar.value();
		ycoori = ROUND(ycoor * coord_scale * scale) + scroll->y() - scroll->scrollbar.value();
		fl_circle(xcoori, ycoori, particle_radius);

		if (do_startend)
		{
			if (icoord % 2 == 1)
			{
				xcoori_start = xcoori;
				ycoori_start = ycoori;
			}
			else
			{
				fl_line(xcoori_start, ycoori_start, xcoori, ycoori);
			}
		}
	}
}

int pickerViewerCanvas::handle(int ev)
{
	const int button = Fl::event_button() ;
	const bool with_shift = (Fl::event_shift() != 0);
	const bool with_control = (Fl::event_ctrl() != 0);
	const int key = Fl::event_key();
	if (ev==FL_PUSH || (ev==FL_DRAG && (button == FL_MIDDLE_MOUSE || (button == FL_LEFT_MOUSE && with_shift))))
	{
		RFLOAT scale = boxes[0]->scale;
		int xc = (int)Fl::event_x() - scroll->x() + scroll->hscrollbar.value();
		int yc = (int)Fl::event_y() - scroll->y() + scroll->scrollbar.value();
		RFLOAT xcoor = (RFLOAT)ROUND(xc / (coord_scale * scale));
		RFLOAT ycoor = (RFLOAT)ROUND(yc / (coord_scale * scale));
		RFLOAT rad2 = particle_radius * particle_radius / (coord_scale * coord_scale * scale * scale);
		if (button == FL_LEFT_MOUSE && !with_shift && !with_control)
		{
			// Left mouse for picking
			// Check the pick is not inside an existing circle
			RFLOAT xcoor_p, ycoor_p;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDcoords)
			{
				MDcoords.getValue(EMDL_IMAGE_COORD_X, xcoor_p);
				MDcoords.getValue(EMDL_IMAGE_COORD_Y, ycoor_p);
				xcoor_p -= xcoor;
				ycoor_p -= ycoor;

				if (xcoor_p*xcoor_p + ycoor_p*ycoor_p < rad2)
					return 0;
			}
			RFLOAT aux = -999., zero = 0.;
			int iaux = current_selection_type;

			// Else store new coordinate
			if (!MDcoords.isEmpty())
			{
				// If there were already entries in MDcoords, then copy the last one.
				// This will take care of re-picking in coordinate files from previous refinements
				long int last_idx = MDcoords.numberOfObjects() - 1;
				MDcoords.addObject(MDcoords.getObject(last_idx));
				RFLOAT aux2;
				if (MDcoords.getValue(EMDL_ORIENT_ROT, aux2))
					MDcoords.setValue(EMDL_ORIENT_ROT, aux);
				if (MDcoords.getValue(EMDL_ORIENT_TILT, aux2))
					MDcoords.setValue(EMDL_ORIENT_TILT, aux);
				if (MDcoords.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, aux2))
					MDcoords.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, zero);
				if (MDcoords.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, aux2))
					MDcoords.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, zero);
			}
			else
				MDcoords.addObject();
			MDcoords.setValue(EMDL_IMAGE_COORD_X, xcoor);
			MDcoords.setValue(EMDL_IMAGE_COORD_Y, ycoor);
			// No autopicking, but still always fill in the parameters for autopicking with dummy values (to prevent problems in joining autopicked and manually picked coordinates)
			MDcoords.setValue(EMDL_PARTICLE_CLASS, iaux);
			MDcoords.setValue(EMDL_ORIENT_PSI, aux);
			MDcoords.setValue(EMDL_PARTICLE_AUTOPICK_FOM, aux);

			redraw();
			return 1;
		}
		else if ((button == FL_MIDDLE_MOUSE) || (button == FL_LEFT_MOUSE && with_shift))
		{
			boxes[0]->redraw();
			// Middle mouse for deleting
			RFLOAT xcoor_p, ycoor_p;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDcoords)
			{
				MDcoords.getValue(EMDL_IMAGE_COORD_X, xcoor_p);
				MDcoords.getValue(EMDL_IMAGE_COORD_Y, ycoor_p);
				xcoor_p -= xcoor;
				ycoor_p -= ycoor;

				if (xcoor_p*xcoor_p + ycoor_p*ycoor_p < rad2)
				{
					MDcoords.removeObject();
					break;
				}
			}
			redraw();
			return 1;
		}
		else if ((button == FL_RIGHT_MOUSE) || (button == FL_LEFT_MOUSE && with_control))
		{
			redraw();
			Fl_Menu_Item rclick_menu[] = {
				{ "Save STAR with coordinates (CTRL-s)" },
//				{ "Save_as STAR with coordinates" },
				{ "Load coordinates" },
				{ "Reload coordinates" },
				{ "Clear coordinates" },
				{ "Set selection type" },
				{ "Set FOM threshold" },
				{ "Help" },
				{ "Quit (CTRL-q)" },
				{ 0 }
			};
			const Fl_Menu_Item *m = rclick_menu->popup(Fl::event_x(), Fl::event_y(), 0, 0, 0);
			if ( !m )
				return 0;
			else if ( strcmp(m->label(), "Save STAR with coordinates (CTRL-s)") == 0 )
				saveCoordinates(false);
//			else if ( strcmp(m->label(), "Save_as STAR with coordinates") == 0 )
//				saveCoordinates(true);
			else if ( strcmp(m->label(), "Load coordinates") == 0 )
				loadCoordinates(true);
			else if ( strcmp(m->label(), "Reload coordinates") == 0 )
				loadCoordinates(false);
			else if ( strcmp(m->label(), "Clear coordinates") == 0 )
				clearCoordinates();
			else if ( strcmp(m->label(), "Set selection type") == 0)
				setSelectionType();
			else if ( strcmp(m->label(), "Set FOM threshold") == 0)
			{
				setFOMThreshold();
				loadCoordinates(false);
			}
			else if ( strcmp(m->label(), "Help") == 0 )
				printHelp();
			else if ( strcmp(m->label(), "Quit (CTRL-q)") == 0 )
				exit(0);
			redraw();
			return 1; // (tells caller we handled this event)
		}
		return 0;
	}
	// Update the drawing every time something happens ....
	else if (ev==FL_RELEASE || ev==FL_LEAVE || ev==FL_ENTER || ev==FL_MOVE || ev == FL_FOCUS || ev == FL_UNFOCUS)
	{
		redraw();
		return 1;
	}
	// CTRL-s will save the coordinates in a picker window
	else if (with_control)
	{
		if (key == 's')
		{
			saveCoordinates(false);
			sleep(1); // to prevent multiple saves... dirty but don't know how to do this otherwise...
			return 1; // (tells caller we handled this event)
		}
		else if (key == 'q')
		{
			sleep(1);
			exit(0);
			return 1; // (tells caller we handled this event)
		}
		else if (key >= '1' && key <= '6')
		{
			std::cout << "debug key = " << key << std::endl;
			current_selection_type = key - '0';
			return 1;
		}
	}
	return 0;
}

void pickerViewerCanvas::saveCoordinates(bool ask_filename)
{
	FileName fn_out;
	if (ask_filename)
	{
		char *newfile;
		newfile = fl_file_chooser("Save File As?", "*.star", "");
		if (newfile == NULL)
			return;
		FileName fn_tmp(newfile);
		fn_out = fn_tmp;
	}
	else
	{
		fn_out = (fn_coords=="") ? "picked.star" : fn_coords;
	}

	FileName fn_dirs = fn_coords.beforeLastOf("/");
	if (!(exists(fn_dirs)))
	{
		std::string command = "mkdir -p " + fn_dirs;
		int res = system(command.c_str());
	}
	// Never write out columns that come from the fn_color file....
	if (fn_color != "" && color_label != EMDL_UNDEFINED)
	{
		MetaDataTable MDnew = MDcoords;
		MDnew.deactivateLabel(color_label);
		MDnew.write(fn_out); // write out a copy of the MDcoord to maintain the Z-score label active...
	}
	else
	{
		MDcoords.write(fn_out);
	}
	std::cout << "Saved "<<fn_out << " with " << MDcoords.numberOfObjects() << " selected coordinates." << std::endl;
}

void pickerViewerCanvas::loadCoordinates(bool ask_filename)
{
	clearCoordinates();
	FileName fn_coord_in;
	if (ask_filename)
	{
		char *newfile;
		newfile = fl_file_chooser("Load File?", "*.star", "");
		if (newfile==NULL)
			return;
		FileName fn_tmp(newfile);
		fn_coord_in = fn_tmp;
	}
	else
	{
		fn_coord_in = (fn_coords=="") ? "picked.star" : fn_coords;
	}
	MDcoords.read(fn_coord_in);

	if (fn_color != "")
	{
		findColorColumnForCoordinates();
	}
}

void pickerViewerCanvas::findColorColumnForCoordinates()
{
	MetaDataTable MDcolor, MDcolormic;
	MDcolor.read(fn_color);

	if (!MDcolor.containsLabel(color_label))
		REPORT_ERROR("--color STAR file does not contain the specified color_label!");

	// Pre-set all color_label column in the MDcoords to -999.
	if (EMDL::isInt(color_label))
	{
		int ival = -999;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDcoords)
		{
			MDcoords.setValue(color_label, ival);
		}
	}
	else
	{
		RFLOAT val = -999.;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDcoords)
		{
			MDcoords.setValue(color_label, val);
		}
	}

	FileName _fn_mic, _fn_img;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDcolor)
	{
		MDcolor.getValue(EMDL_MICROGRAPH_NAME, _fn_mic);
		if (fn_mic == _fn_mic)
		{
			// Get the imagename
			MDcolor.getValue(EMDL_IMAGE_NAME, _fn_img);
			long int iimg;
			std::string dum;
			_fn_img.decompose(iimg, dum);
			iimg--; // counting starts at 1 in STAR file!

			// Check that this entry in the coord file has the same xpos and ypos
			RFLOAT my_xpos, my_ypos;
			MDcoords.getValue(EMDL_IMAGE_COORD_X, my_xpos, iimg);
			MDcoords.getValue(EMDL_IMAGE_COORD_Y, my_ypos, iimg);

			RFLOAT x, y;
			MDcolor.getValue(EMDL_IMAGE_COORD_X, x);
			MDcolor.getValue(EMDL_IMAGE_COORD_Y, y);

			if ( ABS(x - my_xpos) + ABS(y - my_ypos) > XMIPP_EQUAL_ACCURACY)
			{
				std::cerr << " _fn_img= " << _fn_img << " iimg= " << iimg << " _fn_mic= " << _fn_mic << std::endl;
				std::cerr << " x= " << x << " my_xpos= " << my_xpos << std::endl;
				std::cerr << " y= " << y << " my_ypos= " << my_ypos << std::endl;
				REPORT_ERROR("The image in the --color star file does not have the same coordinates as the ones in the --coord file!");
			}
			else
			{
				if (EMDL::isInt(color_label))
				{
					int ival;
					MDcolor.getValue(color_label, ival);
					MDcoords.setValue(color_label, ival, iimg);
				}
				else
				{
					RFLOAT val;
					MDcolor.getValue(color_label, val);
					MDcoords.setValue(color_label, val, iimg);
				}
			}
		}
	}
}

void pickerViewerCanvas::clearCoordinates()
{
	boxes[0]->redraw();
	MDcoords.clear();
}

void pickerViewerCanvas::printHelp()
{
	std::cout <<" + Left-mouse   click: pick particles " << std::endl;
	std::cout <<" + Middle-mouse click: delete particles " << std::endl;
	std::cout <<" + Right-mouse click: pop-up menu" << std::endl;
	std::cout <<" + Refresh picked particles by moving out of the window and back in again ..." << std::endl;
}

void singleViewerCanvas::printMetaData()
{
	boxes[0]->MDimg.write(std::cout);
}

/*
void singleViewerCanvas::setContrast()
{
	popupSetContrastWindow win(400, 100, "Set contrast");
	win.fill();
}
*/

// Fill GUI window
int displayerGuiWindow::fill(FileName &_fn_in)
{
	// display image name at the top
	fn_in = _fn_in;
	FileName fn_short = _fn_in.removeDirectories();

	color(GUI_BACKGROUND_COLOR);

	int width = 485;
	Fl_Text_Display *mydisp = new Fl_Text_Display(15, 15, width-15, 25,"");
	Fl_Text_Buffer *textbuff = new Fl_Text_Buffer();
	textbuff->text(fn_short.c_str());
	mydisp->buffer(textbuff);
	int x=170, y=15, ystep = 27, height = 25,  inputwidth = 50, inputwidth2=30;
	int x2 = width - inputwidth - 50;
	y += ROUND(1.5*ystep);

	// General box
	Fl_Box *box1 = new Fl_Box(15, y-ROUND(0.25*ystep), width - 15, ROUND(2.5*ystep), "");
	box1->color(GUI_BACKGROUND_COLOR);
	box1->box(FL_DOWN_BOX);

	// Always display these:
	scale_input = new Fl_Input(x, y, inputwidth, height, "Scale:");
	scale_input->color(GUI_INPUT_COLOR);
	scale_input->value("1");

	black_input = new Fl_Input(x2-110, y, inputwidth, height, "Min:");
	black_input->value("0");
	black_input->color(GUI_INPUT_COLOR);
	white_input = new Fl_Input(x2, y, inputwidth, height, "Max:");
	white_input->value("0");
	white_input->color(GUI_INPUT_COLOR);
	y += ystep;

	sigma_contrast_input = new Fl_Input(x, y, inputwidth, height, "Sigma contrast:");
	sigma_contrast_input->value("0");
	sigma_contrast_input->color(GUI_INPUT_COLOR);

	colour_scheme_choice = new Fl_Choice(x2-110, y, inputwidth+110, height, "Color:");
	colour_scheme_choice->add("greyscale", 0, 0,0, FL_MENU_VALUE);
	colour_scheme_choice->add("fire", 0, 0,0, FL_MENU_VALUE);
	colour_scheme_choice->add("ice", 0, 0,0, FL_MENU_VALUE);
	colour_scheme_choice->add("fire-n-ice", 0, 0,0, FL_MENU_VALUE);
	colour_scheme_choice->add("rainbow", 0, 0,0, FL_MENU_VALUE);
	colour_scheme_choice->add("difference", 0, 0,0, FL_MENU_VALUE);
	colour_scheme_choice->picked(colour_scheme_choice->menu());
	colour_scheme_choice->color(GUI_INPUT_COLOR);

	y += ROUND(1.75*ystep);

	if (is_star)
	{

		// STAR box
		Fl_Box *box2 = new Fl_Box(15, y-ROUND(0.25*ystep), width - 15, ROUND(3.5*ystep), "");
		box2->color(GUI_BACKGROUND_COLOR);
		box2->box(FL_DOWN_BOX);

		display_choice = new Fl_Choice(x, y, width-x, height, "Display:");
		for (int i = 0; i < display_labels.size(); i++)
			display_choice->add(display_labels[i].c_str(), 0, 0,0, FL_MENU_VALUE);
		display_choice->picked(display_choice->menu());
		display_choice->color(GUI_INPUT_COLOR);
		y += ystep;

		sort_button = new Fl_Check_Button(35,y,height,height, "Sort images ");
		sort_choice = new Fl_Choice(x, y, width-x, height, "on:");
		for (int i = 0; i < sort_labels.size(); i++)
			sort_choice->add(sort_labels[i].c_str(), 0, 0,0, FL_MENU_VALUE);
		sort_choice->add("RANDOMLY", 0, 0,0, FL_MENU_VALUE);
		sort_choice->picked(sort_choice->menu());
		sort_choice->color(GUI_INPUT_COLOR);
		y += ystep;

		reverse_sort_button  = new Fl_Check_Button(35, y, inputwidth, height, "Reverse sort?");
		reverse_sort_button->color(GUI_INPUT_COLOR);
		apply_orient_button  = new Fl_Check_Button(x, y, inputwidth, height, "Apply orientations?");
		apply_orient_button->color(GUI_INPUT_COLOR);
		display_label_button = new Fl_Check_Button(x+160, y, inputwidth, height, "Display label?");
		display_label_button->color(GUI_INPUT_COLOR);
		y += ROUND(1.75*ystep);

	}

	// Optional display options
	if (is_multi)
	{
		// Multiview box
		Fl_Box *box3;
		if (do_allow_save && fn_parts != "")
			box3 = new Fl_Box(15, y-ROUND(0.25*ystep), width - 15, ROUND(2.5*ystep), "");
		else
			box3 = new Fl_Box(15, y-ROUND(0.25*ystep), width - 15, ROUND(1.5*ystep), "");

		box3->color(GUI_BACKGROUND_COLOR);
		box3->box(FL_DOWN_BOX);

		int x1p = 125;
		int x2p = x1p + 130;
		int x3p = x2p + 170;

		col_input = new Fl_Input(x1p, y, 40, height, "Nr. columns:");
		col_input->value("5");
		col_input->color(GUI_INPUT_COLOR);

		ori_scale_input = new Fl_Input(x2p, y, 40, height, "Ori scale:");
		ori_scale_input->value("1");
		ori_scale_input->color(GUI_INPUT_COLOR);

		max_nr_images_input = new Fl_Input(x3p, y, 40, height, "Max. nr. images:");
		max_nr_images_input->value("1000");
		max_nr_images_input->color(GUI_INPUT_COLOR);

		if (do_allow_save && fn_parts != "")
		{
			y += ystep;
			max_parts_per_class_input = new Fl_Input(x2p, y, 40, height, "Max nr selected parts per class:");
			max_parts_per_class_input->value("-1");
			max_parts_per_class_input->color(GUI_INPUT_COLOR);
			y += ROUND(1.75*ystep);
		}
		else
			y += ROUND(1.75*ystep);

	}
	else // is_single
	{
		// singleview box
		Fl_Box *box3 = new Fl_Box(15, y-ROUND(0.25*ystep), width - 15, ROUND(1.5*ystep), "");
		box3->color(GUI_BACKGROUND_COLOR);
		box3->box(FL_DOWN_BOX);

		lowpass_input = new Fl_Input(x, y, inputwidth2, height, "Lowpass filter (A):");
		lowpass_input->color(GUI_INPUT_COLOR);
		lowpass_input->value("0");

		highpass_input = new Fl_Input(275, y, inputwidth2, height, "Highpass:");
		highpass_input->color(GUI_INPUT_COLOR);
		highpass_input->value("0");

		angpix_input = new Fl_Input(x2+inputwidth-inputwidth2, y, inputwidth2, height, "Pixel size (A):");
		angpix_input->color(GUI_INPUT_COLOR);
		angpix_input->value("1");

		y += ROUND(1.75*ystep);
	}

	// Display button
	Fl_Button *display_button = new Fl_Button(width-100, y, 100, 30, "Display!");
	display_button->color(GUI_RUNBUTTON_COLOR);
	display_button->callback( cb_display, this);

	// Read last settings file if it is present
	readLastSettings();

	show();
	return Fl::run();
}

void displayerGuiWindow::readLastSettings()
{
	FileName fn = ".relion_display_gui_settings";
	if (!exists(fn))
		return;

	std::ifstream in(fn.c_str(), std::ios_base::in);
    if (in.fail())
    {
    	std::cerr << "Error reading last settings from file: "<< fn<<std::endl;
    	return;
    }

	in.seekg(0, std::ios::beg);
	std::string line;
	while (getline(in, line, '\n'))
	{
		int ispos = line.rfind("=");
		std::string label = line.substr(0, ispos - 1);
		std::string value = line.substr(ispos+2,line.length());
		if (label == scale_input->label())
			scale_input->value(value.c_str());
		else if (label == black_input->label())
			black_input->value(value.c_str());
		else if (label == white_input->label())
			white_input->value(value.c_str());
		else if (label == colour_scheme_choice->label())
			colour_scheme_choice->value(textToInteger(value));
		else if (label == sigma_contrast_input->label())
			sigma_contrast_input->value(value.c_str());
		else if (is_multi && label == col_input->label())
			col_input->value(value.c_str());
		else if (is_multi && label == ori_scale_input->label())
			ori_scale_input->value(value.c_str());
		else if (is_multi && label == max_nr_images_input->label())
			max_nr_images_input->value(value.c_str());
		else if (!is_multi && label == lowpass_input->label())
			lowpass_input->value(value.c_str());
		else if (!is_multi && label == highpass_input->label())
			highpass_input->value(value.c_str());
		else if (!is_multi && label == angpix_input->label())
			angpix_input->value(value.c_str());
	}

	in.close();
}

void displayerGuiWindow::writeLastSettings()
{
	std::ofstream  fh;
	FileName fn = ".relion_display_gui_settings";
	fh.open(fn.c_str(), std::ios::out);
	if (!fh)
	{
		//std::cerr << "Cannot write last settings to file: "<<fn<<std::endl;
		return;
	}

	fh << scale_input->label() << " = " << scale_input->value() << std::endl;
	fh << black_input->label() << " = " << black_input->value() << std::endl;
	fh << white_input->label() << " = " << white_input->value() << std::endl;
	fh << colour_scheme_choice->label() << " = " << colour_scheme_choice->value() << std::endl;
	fh << sigma_contrast_input->label() << " = " << sigma_contrast_input->value() << std::endl;
	if (is_multi)
	{
		fh << col_input->label() << " = " << col_input->value() << std::endl;
		fh << ori_scale_input->label() << " = " << ori_scale_input->value() << std::endl;
		fh << max_nr_images_input->label() << " = " << max_nr_images_input->value() << std::endl;
	}
	else
	{
		fh << lowpass_input->label() << " = " << lowpass_input->value() << std::endl;
		fh << highpass_input->label() << " = " << highpass_input->value() << std::endl;
		fh << angpix_input->label() << " = " << angpix_input->value() << std::endl;
	}

	fh.close();
}

// Display button call-back functions
void displayerGuiWindow::cb_display(Fl_Widget* o, void* v) {
	displayerGuiWindow* T=(displayerGuiWindow*)v;
	T->cb_display_i();
}

void displayerGuiWindow::cb_display_i()
{
	// Save last settings, so we don't need to change settings every time...
	writeLastSettings();

	// This is a rather ugly system call to the relion_display program again,
	// but I do not know how to get back to the original Displayer class from here...
	std::string cl = "relion_display ";
	cl += " --i " + fn_in;

	// Always
	cl += " --scale " + (std::string)scale_input->value();
	cl += " --black " + (std::string)black_input->value();
	cl += " --white " + (std::string)white_input->value();
	cl += " --sigma_contrast " + (std::string)sigma_contrast_input->value();

	// Get the colour scheme
	const Fl_Menu_Item* m3 = colour_scheme_choice->mvalue();
	if ((std::string)m3->label() == "fire") cl += " --colour_fire ";
	else if ((std::string)m3->label() == "ice") cl += " --colour_ice ";
	else if ((std::string)m3->label() == "fire-n-ice") cl += " --colour_fire-n-ice ";
	else if ((std::string)m3->label() == "rainbow") cl += " --colour_rainbow ";
	else if ((std::string)m3->label() == "difference") cl += " --colour_difference ";

	if (is_star)
	{

		const Fl_Menu_Item* m = display_choice->mvalue();
		cl += " --display " + (std::string)m->label();

		if (getValue(sort_button))
		{
			const Fl_Menu_Item* m2 = sort_choice->mvalue();
			if ((std::string)m2->label() == "RANDOMLY")
			{
				cl += " --random_sort ";
			}
			else
			{
				cl += " --sort " + (std::string)m2->label();

				if (getValue(reverse_sort_button))
					cl += " --reverse ";
			}
		}

		if (getValue(display_label_button))
		{
			const Fl_Menu_Item* m2 = sort_choice->mvalue();
			cl += " --text_label " + (std::string)m2->label();
		}

		if (getValue(apply_orient_button))
			cl += " --apply_orient ";
	}

	if (is_multi)
	{
		cl += " --col " + (std::string)col_input->value();
		cl += " --ori_scale " + (std::string)ori_scale_input->value();
		if (textToInteger(max_nr_images_input->value()) > 0)
		{
			if (getValue(sort_button))
				std::cerr << " WARNING: you cannot sort particles and use a maximum number of images. Ignoring the latter..." << std::endl;
			else
				cl += " --max_nr_images " + (std::string)max_nr_images_input->value();
		}
	}
	else
	{
		//check for individual images
		cl += " --lowpass " + (std::string)lowpass_input->value();
		cl += " --highpass " + (std::string)highpass_input->value();
		cl += " --angpix " + (std::string)angpix_input->value();
	}

	if (is_class)
	{
		cl += " --class ";
	}

	if (do_allow_save)
	{
		cl += " --allow_save ";
		if (fn_parts != "")
		{
			cl += " --fn_parts " + fn_parts;
			if ( textToInteger(max_parts_per_class_input->value()) > 0)
				cl += " --max_nr_parts_per_class " + (std::string)max_parts_per_class_input->value();
		}
		if (fn_imgs != "")
			cl += " --fn_imgs " + fn_imgs;
	}

	if (nr_regroups > 0)
	{
		cl += " --regroup " + integerToString(nr_regroups);
	}

	if (do_recenter)
	{
		cl += " --recenter";
	}

	if (pipeline_control != "")
	{
		cl += " --pipeline_control " + pipeline_control;
	}

	// send job in the background
	cl += " &";
	//std::cout << "Executing: " << cl << std::endl;
	int res = system(cl.c_str());
}


void Displayer::topazDenoiseMap(FileName fn_in, FileName fn_odir, Image<RFLOAT> &img)
{

    if (fn_odir[fn_odir.length()-1] != '/')
             fn_odir += "/";

	if (!(exists(fn_odir)))
	{
		std::string command = "mkdir -p " + fn_odir;
		int res = system(command.c_str());
	}

	// Now generate the topaz bash script
	FileName fn_script = fn_odir+"topaz_denoise.bash";
	FileName fn_log    = fn_odir+"topaz_denoise.out";
	FileName fn_out = fn_odir + fn_in.afterLastOf("/");
	std::ofstream  fh;
	fh.open((fn_script).c_str(), std::ios::out);
	if (!fh)
	 REPORT_ERROR( (std::string)"AutoPicker::trainTopaz cannot create file: " + fn_script);

	fh << "#!" << fn_shell  << std::endl;

	// Call Topaz to train the network
	fh << fn_topaz_exe << " denoise ";
	fh << fn_in;
	fh << " --output " << fn_odir;
	fh << " --device 0"; // pyTorch threads
	fh << " --patch-size 1024";
	fh << " --patch-padding 256";
	fh << " --normalize";
	fh << std::endl;
	fh.close();

	std::string command = fn_shell + " " + fn_script + " >> " + fn_log + " 2>&1";
	if (system(command.c_str()))
	{
		std::string warning = "WARNING: there was an error in executing: " + fn_script + "\nWARNING: skipping topaz denoising...";
		fl_message("%s", warning.c_str());
		img.read(fn_in);
	}
	else
	{

		// Now read the denoise micrograph back in
		img.read(fn_out);
		std::remove(fn_out.c_str());
		std::remove(fn_script.c_str());
		std::remove(fn_log.c_str());
	}

}


void Displayer::read(int argc, char **argv)
{
	parser.setCommandLine(argc, argv);

	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "Input STAR file, image or stack","");
	do_gui = parser.checkOption("--gui", "Use this to provide all other parameters through a GUI");
	display_label = EMDL::str2Label(parser.getOption("--display", "Metadata label to display", "rlnImageName"));
	text_label = EMDL::str2Label(parser.getOption("--text_label", "Metadata label to display text", "EMDL_UNDEFINED"));
	table_name = parser.getOption("--table", "Name of the table to read from in the input STAR file", "");
	scale = textToFloat(parser.getOption("--scale", "Relative scale", "1"));
	minval = textToFloat(parser.getOption("--black", "Pixel value for black (default is auto-contrast)", "0"));
	maxval = textToFloat(parser.getOption("--white", "Pixel value for white (default is auto-contrast)", "0"));
	sigma_contrast  = textToFloat(parser.getOption("--sigma_contrast", "Set white and black pixel values this many times the image stddev from the mean", "0"));
	do_read_whole_stacks = parser.checkOption("--read_whole_stack", "Read entire stacks at once (to speed up when many images of each stack are displayed)");
	show_fourier_amplitudes = parser.checkOption("--show_fourier_amplitudes", "Show amplitudes of 2D Fourier transform?");
	show_fourier_phase_angles = parser.checkOption("--show_fourier_phase_angles", "Show phase angles of 2D Fourier transforms?");
	if (parser.checkOption("--colour_fire", "Show images in black-grey-white-red colour scheme (highlight high signal)?")) colour_scheme = BLACKGREYREDSCALE;
	else if (parser.checkOption("--colour_ice", "Show images in blue-black-grey-white colour scheme (highlight low signal)?")) colour_scheme = BLUEGREYWHITESCALE;
	else if (parser.checkOption("--colour_fire-n-ice", "Show images in blue-grey-red colour scheme (highlight high&low signal)?")) colour_scheme = BLUEGREYREDSCALE;
	else if (parser.checkOption("--colour_rainbow", "Show images in cyan-blue-black-red-yellow colour scheme?")) colour_scheme = RAINBOWSCALE;
	else if (parser.checkOption("--colour_difference", "Show images in cyan-blue-black-red-yellow colour scheme (for difference images)?")) colour_scheme = CYANBLACKYELLOWSCALE;
	else colour_scheme = GREYSCALE;
	do_colourbar = parser.checkOption("--colour_bar", "Show colourbar image?");
	do_ignore_optics = parser.checkOption("--ignore_optics", "Ignore information about optics groups in input STAR file?");

	int disp_section  = parser.addSection("Multiviewer options");
	ncol = textToInteger(parser.getOption("--col", "Number of columns", "5"));
	do_apply_orient = parser.checkOption("--apply_orient","Apply the orientation as stored in the input STAR file angles and offsets");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in A) to calculate lowpass filter and/or translational offsets ", "-1"));
	ori_scale = textToFloat(parser.getOption("--ori_scale", "Relative scale for viewing individual images in multiviewer", "1"));
	sort_label = EMDL::str2Label(parser.getOption("--sort", "Metadata label to sort images on", "EMDL_UNDEFINED"));
	random_sort = parser.checkOption("--random_sort", "Use random order in the sorting");
	reverse_sort = parser.checkOption("--reverse", "Use reverse order (from high to low) in the sorting");
	do_class = parser.checkOption("--class", "Use this to analyse classes in input optimiser.star or model.star file");
	nr_regroups = textToInteger(parser.getOption("--regroup", "Number of groups to regroup saved particles from selected classes in (default is no regrouping)", "-1"));
	do_allow_save = parser.checkOption("--allow_save", "Allow saving of selected particles or class averages");

	fn_selected_imgs = parser.getOption("--fn_imgs", "Name of the STAR file in which to save selected images.", "");
	fn_selected_parts = parser.getOption("--fn_parts", "Name of the STAR file in which to save particles from selected classes.", "");
	max_nr_parts_per_class  = textToInteger(parser.getOption("--max_nr_parts_per_class", "Select maximum this number of particles from each selected classes.", "-1"));
	do_recenter = parser.checkOption("--recenter", "Recenter the selected images to the center-of-mass of all positive pixel values. ");
	max_nr_images = textToInteger(parser.getOption("--max_nr_images", "Only show this many images (default is show all)", "-1"));

	int pick_section  = parser.addSection("Picking options");
	do_pick = parser.checkOption("--pick", "Pick coordinates in input image");
	do_pick_startend = parser.checkOption("--pick_start_end", "Pick start-end coordinates in input image");
	fn_coords = parser.getOption("--coords", "STAR file with picked particle coordinates", "");
	coord_scale = textToFloat(parser.getOption("--coord_scale", "Scale particle coordinates before display", "1.0"));
	particle_radius = textToFloat(parser.getOption("--particle_radius", "Particle radius in pixels", "100"));
	particle_radius *= coord_scale;
	do_topaz_denoise = parser.checkOption("--topaz_denoise", "Use Topaz denoising before picking (on GPU 0)");
	fn_topaz_exe = parser.getOption("--topaz_exe", "Name of topaz executable", "topaz");
        fn_shell = parser.getOption("--bash_exe", "Name of bash shell executable", "/bin/bash");
	lowpass = textToFloat(parser.getOption("--lowpass", "Lowpass filter (in A) to filter micrograph before displaying", "0"));
	highpass = textToFloat(parser.getOption("--highpass", "Highpass filter (in A) to filter micrograph before displaying", "0"));
	minimum_pick_fom = textToFloat(parser.getOption("--minimum_pick_fom", "Minimum value for rlnAutopickFigureOfMerit to display picks", "-9999."));
	fn_color = parser.getOption("--color_star", "STAR file with a column for red-blue coloring (a subset of) the particles", "");
	color_label = parser.getOption("--color_label", "MetaDataLabel to color particles on (e.g. rlnParticleSelectZScore)", "");
	color_blue_value = textToFloat(parser.getOption("--blue", "Value of the blue color", "1."));
	color_red_value = textToFloat(parser.getOption("--red", "Value of the red color", "0."));

	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void Displayer::usage()
{
	parser.writeUsage(std::cout);
}

void Displayer::initialise()
{
	if (!do_gui && fn_in=="")
		REPORT_ERROR("Displayer::initialise ERROR: either provide --i or --gui");
	Fl::visual(FL_RGB);
	// initialise some static variables
	has_dragged = false;
	has_shift = false;

	if (do_class)
	{
		display_label = EMDL_MLMODEL_REF_IMAGE;
		table_name = "model_classes";
		FileName fn_data, fn_model = fn_in;
		if (fn_in.contains("_half1_model.star"))
			fn_data = fn_in.without("_half1_model.star") + "_data.star";
		else if (fn_in.contains("_half2_model.star"))
			fn_data = fn_in.without("_half2_model.star") + "_data.star";
		// SHWS 28nov2019: backwards compatibility
		else if (fn_in.contains("_model.star"))
			fn_data = fn_in.without("_model.star") + "_data.star";
		else if (fn_in.contains("_optimiser.star"))
		{
			fn_data = fn_in.without("_optimiser.star") + "_data.star";
			MetaDataTable MDopt;
			MDopt.read(fn_in, "optimiser_general");
			MDopt.getValue(EMDL_OPTIMISER_MODEL_STARFILE, fn_model);
		}

		if (do_ignore_optics) MDdata.read(fn_data);
		else ObservationModel::loadSafely(fn_data, obsModel, MDdata);

		// If regrouping, also read the model_groups table into memory
		if (nr_regroups > 0)
			MDgroups.read(fn_model, "model_groups");
	}

	// Also allow regrouping on data.star
	if (fn_in.contains("_data.star") && nr_regroups > 0)
	{
		FileName fn_model;
		fn_model = fn_in.without("_data.star") + "_model.star";
		bool has_model = false;
		if (exists(fn_model))
		{
			MDgroups.read(fn_model, "model_groups");
			has_model = true;
		}
		else
		{
			fn_model = fn_in.without("_data.star") + "_half1_model.star";
			if (exists(fn_model))
			{
				MDgroups.read(fn_model, "model_groups");
				has_model = true;
			}
		}
		if (!has_model)
			std::cout <<" Warning: cannot find model.star file for " << fn_in << " needed for regrouping..." << std::endl;

	}

	// Check if input STAR file contains pixel-size information
	if (!do_class && (do_apply_orient || lowpass > 0 || highpass > 0))
	{
		if (fn_in.isStarFile())
		{
			// As of v3.1 the input STAR files should always store the pixel size, no more check necessary...
			if (do_ignore_optics)
			{
				if (angpix > 0.)
				{
					obsModel.opticsMdt.setValue(EMDL_IMAGE_PIXEL_SIZE, angpix);
				}
			}
		}
		else
		{
			// if not a STAR file: always need command-line input for pixel
			obsModel.opticsMdt.addObject();
			if (angpix > 0.)
			{
				//std::cout << " Using pixel size from command-line input of " << angpix << " Angstroms" << std::endl;
				obsModel.opticsMdt.setValue(EMDL_IMAGE_PIXEL_SIZE, angpix);
			}
			else
			{
				REPORT_ERROR("Displayer::initialise ERROR: you provided a low- or highpass filter in Angstroms, so please also provide --angpix.");
			}
		}
	}


	if (show_fourier_amplitudes && show_fourier_phase_angles)
		REPORT_ERROR("Displayer::initialise ERROR: cannot display Fourier amplitudes and phase angles at the same time!");
	if (show_fourier_amplitudes || show_fourier_phase_angles)
	{
		if (do_pick || do_pick_startend)
			REPORT_ERROR("Displayer::initialise ERROR: cannot pick particles from Fourier maps!");
		if (fn_in.isStarFile())
			REPORT_ERROR("Displayer::initialise ERROR: use single 2D image files as input!");
		Image<RFLOAT> img;
		img.read(fn_in, false); // dont read data yet: only header to get size
		if ( (ZSIZE(img()) > 1) || (NSIZE(img()) > 1) )
			REPORT_ERROR("Displayer::initialise ERROR: cannot display Fourier maps for 3D images or stacks!");
	}

}

int Displayer::runGui()
{
	Fl::scheme("gtk+");

	if (fn_in == "")
	{
		// Shall I make a browser window in this GUI or in the general relion GUI?
		// Perhaps here is better..., then there will be no fn_in yet....
		// Update entire window each time the entry of the browser changes...
		Fl_File_Chooser chooser(".",                        // directory
		                        "All recognised formats (*.{star,mrc,mrcs})\tSTAR Files (*.star)\tMRC stack (*.mrcs)\tMRC image (*.mrc)\tAll Files (*)*", // filter
		                        Fl_File_Chooser::SINGLE,     // chooser type
		                        "Choose file to display");        // title
		chooser.show();
		// Block until user picks something.
		while(chooser.shown())
			{ Fl::wait(); }

		// User hit cancel?
		if ( chooser.value() == NULL )
			exit(0);
		FileName _fn_in(chooser.value());
		fn_in = _fn_in;
	}

	// make a bigger window for STAR files...
	int windowheight = fn_in.isStarFile() ? 350 : 300;

	displayerGuiWindow win(500, windowheight, "Relion display GUI");
	win.is_class = false;
	win.is_data = false;
	win.is_star = false;
	win.is_multi = false;
	win.do_allow_save = do_allow_save;
	win.nr_regroups = nr_regroups;
	win.do_recenter = do_recenter;
	win.fn_imgs = fn_selected_imgs;
	win.fn_parts = fn_selected_parts;
	win.pipeline_control = pipeline_control_outputname;

	// this GUI can never create the output itself, so disable pipeline_control
	pipeline_control_outputname = "";

	// If this is a STAR file, decide what to do
	if (fn_in.isStarFile())
	{
		MetaDataTable MD;
		win.is_star = true;
		win.is_multi = true;
		win.is_data = fn_in.contains("_data.star");
		if (fn_in.contains("_optimiser.star"))
		{
			win.fn_data = fn_in.without("_optimiser.star") + "_data.star";
			win.is_class = true;

			MetaDataTable MDopt;
			MDopt.read(fn_in, "optimiser_general");
			FileName fn_model;
			MDopt.getValue(EMDL_OPTIMISER_MODEL_STARFILE, fn_model);
			MD.read(fn_model, "model_classes");
		}
		// SHWS 28nov2019: backwards compatibility
		else if (fn_in.contains("_model.star"))
		{
			win.fn_data = fn_in.without("_model.star") + "_data.star";
			win.is_class = true;
			MD.read(fn_in, "model_classes");
		}
		else
		{
			if (do_ignore_optics) MD.read(fn_in);
			else ObservationModel::loadSafely(fn_in, obsModel, MD, "discover", 0, false); //false means dont die upon error
			// Check the MD was loaded successfully with obsModel, otherwise read as ignore_optics
			if (obsModel.opticsMdt.numberOfObjects() == 0) MD.read(fn_in);
		}

		// Get which labels are stored in this metadatatable and generate choice menus for display and sorting

		std::vector<EMDLabel> activeLabels = MD.getActiveLabels();
		for (int ilab = 0; ilab < activeLabels.size(); ilab++)
		{
			if (EMDL::isNumber(activeLabels[ilab]))
				win.sort_labels.push_back(EMDL::label2Str(activeLabels[ilab]));
		}

		// Preferred order of defaults!
		// If EMDL_IMAGE_NAME is among the labels: make that the default choice!)
		if (MD.containsLabel(EMDL_IMAGE_NAME))
			win.display_labels.push_back(EMDL::label2Str(EMDL_IMAGE_NAME));
		if (MD.containsLabel(EMDL_IMAGE_ORI_NAME))
			win.display_labels.push_back(EMDL::label2Str(EMDL_IMAGE_ORI_NAME));
		if (MD.containsLabel(EMDL_MLMODEL_REF_IMAGE))
			win.display_labels.push_back(EMDL::label2Str(EMDL_MLMODEL_REF_IMAGE));
		if (MD.containsLabel(EMDL_CTF_IMAGE))
			win.display_labels.push_back(EMDL::label2Str(EMDL_CTF_IMAGE));
		if (MD.containsLabel(EMDL_MICROGRAPH_NAME))
			win.display_labels.push_back(EMDL::label2Str(EMDL_MICROGRAPH_NAME));
		if (MD.containsLabel(EMDL_CTF_POWER_SPECTRUM))
			win.display_labels.push_back(EMDL::label2Str(EMDL_CTF_POWER_SPECTRUM));
		if (MD.containsLabel(EMDL_MICROGRAPH_MOVIE_NAME))
			win.display_labels.push_back(EMDL::label2Str(EMDL_MICROGRAPH_MOVIE_NAME));
	}
	else
	{
		// Try reading as an image/stack header
		Image<RFLOAT> img;
		img.read(fn_in, false);
		win.is_multi = (ZSIZE(img()) * NSIZE(img()) > 1);
	}

	return win.fill(fn_in);
}


void Displayer::run()
{
	if (do_gui)
	{
	}
	else if (do_colourbar)
	{
		Image<RFLOAT> img(256, 10);
		FOR_ALL_ELEMENTS_IN_ARRAY2D(img())
		{
			A2D_ELEM(img(), i, j) = (RFLOAT)j;
		}
		FileName fnt="colour scheme";
		basisViewerWindow win(CEIL(scale*XSIZE(img())), CEIL(scale*YSIZE(img())), fnt.c_str());
		win.fillSingleViewerCanvas(img(), 0., 255., 0., scale);
	}
	else if (do_pick || do_pick_startend)
	{
		Image<RFLOAT> img;

		if (do_topaz_denoise)
		{
			topazDenoiseMap(fn_in, fn_coords.beforeLastOf("/"), img);
		}
		else
		{
			img.read(fn_in); // dont read data yet: only header to get size
		}

		if (lowpass > 0.)
			lowPassFilterMap(img(), lowpass, angpix);
		if (highpass > 0.)
			highPassFilterMap(img(), highpass, angpix, 25); // use a rather soft high-pass edge of 25 pixels wide

		basisViewerWindow win(CEIL(scale*XSIZE(img())), CEIL(scale*YSIZE(img())), fn_in.c_str());
		if (fn_coords=="")
			fn_coords = fn_in.withoutExtension()+"_coords.star";
		win.fillPickerViewerCanvas(img(), minval, maxval, sigma_contrast, scale, coord_scale, ROUND(scale*particle_radius), do_pick_startend, fn_coords,
    		fn_color, fn_in, color_label, color_blue_value, color_red_value, minimum_pick_fom);
	}
	else if (fn_in.isStarFile())
	{
		if (fn_in.contains("_optimiser.star"))
		{
			MetaDataTable MDopt;
			MDopt.read(fn_in, "optimiser_general");
			FileName fn_model;
			MDopt.getValue(EMDL_OPTIMISER_MODEL_STARFILE, fn_model);
			MDin.read(fn_model, "model_classes");
		}
		// SHWS 28nov2019: Backwards compatibility
		else if (fn_in.contains("_model.star"))
		{
			MDin.read(fn_in, "model_classes");
			MetaDataTable MD2;
			MD2.read(fn_in, "model_general");
			if(MD2.containsLabel(EMDL_MLMODEL_IS_HELIX))
				MDin.addLabel(EMDL_MLMODEL_IS_HELIX);
		}
		else
		{
			if (do_ignore_optics) MDin.read(fn_in);
			else ObservationModel::loadSafely(fn_in, obsModel, MDin, "discover", 0, false); //false means dont die upon error
			// Check the MD was loaded successfully with obsModel, otherwise read as ignore_optics
			if (obsModel.opticsMdt.numberOfObjects() == 0) MDin.read(fn_in);
		}

		// Check that label to display is present in the table
		if (!MDin.containsLabel(display_label))
			REPORT_ERROR("Cannot find metadata label in input STAR file");

		// Store class number in metadata table
		if (do_class)
		{
			int iclass = 0;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
			{
				iclass++; // start counting at 1
				MDin.setValue(EMDL_PARTICLE_CLASS, iclass);
			}
		}

		if (sort_label != EMDL_UNDEFINED || random_sort)
		{
			MDin.sort(sort_label, reverse_sort, true, random_sort); // true means only store sorted_idx!
			// When sorting: never read in the whole stacks!
			do_read_whole_stacks = false;
		}


		basisViewerWindow win(MULTIVIEW_WINDOW_WIDTH, MULTIVIEW_WINDOW_HEIGHT, fn_in.c_str());
		win.fillCanvas(MULTIVIEWER, MDin, &obsModel, display_label, text_label, do_read_whole_stacks, do_apply_orient, minval, maxval, sigma_contrast, scale, ori_scale, ncol,
				max_nr_images,  lowpass, highpass, do_class, &MDdata, nr_regroups, do_recenter, fn_in.contains("_data.star"), &MDgroups,
				do_allow_save, fn_selected_imgs, fn_selected_parts, max_nr_parts_per_class);
	}
	else
	{
		// Attempt to read a single-file image
		Image<RFLOAT> img;
		img.read(fn_in, false); // dont read data yet: only header to get size

		MDin.clear();
		// display stacks
		if (NSIZE(img()) > 1)
		{
			for (int n = 0; n < NSIZE(img()); n++)
			{
				FileName fn_tmp;
				fn_tmp.compose(n+1,fn_in);
				MDin.addObject();
				MDin.setValue(EMDL_IMAGE_NAME, fn_tmp);
				MDin.setValue(EMDL_IMAGE_OPTICS_GROUP, 1);
			}
			basisViewerWindow win(MULTIVIEW_WINDOW_WIDTH, MULTIVIEW_WINDOW_HEIGHT, fn_in.c_str());
			win.fillCanvas(MULTIVIEWER, MDin, &obsModel, EMDL_IMAGE_NAME, text_label, true, false, minval, maxval, sigma_contrast, scale, ori_scale, ncol, max_nr_images, lowpass, highpass);
		}
		else if (ZSIZE(img()) > 1)
		{

			// Read volume slices from .mrc as if it were a .mrcs stack and then use normal slice viewer
			// This will not work for Spider volumes...
			if (fn_in.getFileFormat() != "mrc")
				REPORT_ERROR("Displayer::run() ERROR: only MRC maps are allowed...");

			// Use a single minval and maxval for all slice
			if (minval == maxval)
			{
				Image<RFLOAT> It;
				It.read(fn_in);
				It().computeDoubleMinMax(minval, maxval);
			}

			// Trick MD with :mrcs extension....
			for (int n = 0; n < ZSIZE(img()); n++)
			{
				FileName fn_tmp;
				fn_tmp.compose(n+1,fn_in);
				fn_tmp += ":mrcs";
				MDin.addObject();
				MDin.setValue(EMDL_IMAGE_NAME, fn_tmp);
				MDin.setValue(EMDL_IMAGE_OPTICS_GROUP, 1);
			}

			basisViewerWindow win(MULTIVIEW_WINDOW_WIDTH, MULTIVIEW_WINDOW_HEIGHT, fn_in.c_str());
			win.fillCanvas(MULTIVIEWER, MDin, &obsModel, EMDL_IMAGE_NAME, text_label, true, false, minval, maxval, sigma_contrast, scale, ori_scale, ncol, max_nr_images, lowpass, highpass);
		}
		else
		{
			img.read(fn_in); // now read image data as well (not only header)

			if (lowpass > 0.)
				lowPassFilterMap(img(), lowpass, angpix);
			if (highpass > 0.)
				highPassFilterMap(img(), highpass, angpix);

			MDin.addObject();
			MDin.setValue(EMDL_IMAGE_NAME, fn_in);
			MDin.setValue(EMDL_IMAGE_OPTICS_GROUP, 1);
			RFLOAT new_scale = scale;
			if (show_fourier_amplitudes || show_fourier_phase_angles)
				new_scale *= 2.;
			basisViewerWindow win(CEIL(new_scale*XSIZE(img())), CEIL(new_scale*YSIZE(img())), fn_in.c_str());
			if (show_fourier_amplitudes)
			{
				amplitudeOrPhaseMap(img(), img(), AMPLITUDE_MAP);
				win.fillSingleViewerCanvas(img(), minval, maxval, sigma_contrast, scale);
			}
			else if (show_fourier_phase_angles)
			{
				amplitudeOrPhaseMap(img(), img(), PHASE_MAP);
				win.fillSingleViewerCanvas(img(), -180., 180., 0., scale);
			}
			else
			{
				win.fillSingleViewerCanvas(img(), minval, maxval, sigma_contrast, scale);
			}
		}
	}
}

