/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres", "Takanori Nakane"
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

#include "src/micrograph_model.h"
#include "src/metadata_table.h"
#include "src/image.h"

/* Work in progress!
 *
 * Implement position-dependent motion model
 * Refactor relion_preprocess (windowing) & particle_polished (dose-weighted sum)
 * Use factory pattern?
 */

const RFLOAT Micrograph::NOT_OBSERVED = -9999.0;

void Micrograph::setMovie(FileName fnMovie, FileName fnGain) {
	Image<RFLOAT> Ihead;
	Ihead.read(fnMovie, false);
        
	width = XSIZE(Ihead());
	height = YSIZE(Ihead());
	nFrame = NSIZE(Ihead());

	globalShiftX.resize(nFrame, NOT_OBSERVED);
	globalShiftY.resize(nFrame, NOT_OBSERVED);

	this->fnMovie = fnMovie;
	this->fnGain = fnGain;
}

// Read from a STAR file
void Micrograph::read(FileName fn_in)
{
	// Clear current model
	clear();

	// Open input file
	std::ifstream in(fn_in.data(), std::ios_base::in);
	if (in.fail()) {
		REPORT_ERROR( (std::string) "MicrographModel::read: File " + fn_in + " cannot be read." );
	}

	MetaDataTable MDglobal;

	// Read Image metadata
	MDglobal.readStar(in, "general");

	if (!MDglobal.getValue(EMDL_IMAGE_SIZEX, width) ||
	    !MDglobal.getValue(EMDL_IMAGE_SIZEY, height) ||
	    !MDglobal.getValue(EMDL_IMAGE_SIZEZ, nFrame) ||
	    !MDglobal.getValue(EMDL_MICROGRAPH_MOVIE_NAME, fnMovie)) {
		REPORT_ERROR("MicrographModel::read: insufficient general information");
	}

	globalShiftX.resize(nFrame, NOT_OBSERVED);
	globalShiftY.resize(nFrame, NOT_OBSERVED);

	MDglobal.getValue(EMDL_MICROGRAPH_GAIN_NAME, fnGain);
	
	// Read global shifts
	int frame;
	RFLOAT shiftX, shiftY;

	MDglobal.readStar(in, "global_shift");

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDglobal)
	{
		if (!MDglobal.getValue(EMDL_MICROGRAPH_FRAME_NUMBER, frame) ||
		    !MDglobal.getValue(EMDL_MICROGRAPH_SHIFT_X, shiftX) ||
		    !MDglobal.getValue(EMDL_MICROGRAPH_SHIFT_Y, shiftY)) { 
			REPORT_ERROR("MicrographModel::read: incorrect global_shift table");
		}

		// frame is 1-indexed!
		globalShiftX[frame - 1] = shiftX;
		globalShiftY[frame - 1] = shiftY;
		std::cout << " global shift: frame #" << frame << " x " << shiftX << " Y " << shiftY << std::endl;
	}

	model = new MotionModel();	
}

// Write to a STAR file
void Micrograph::write(FileName filename) {
	std::ofstream fh;
	MetaDataTable MD;

	fh.open(filename.c_str());
	if (!fh) {
		REPORT_ERROR((std::string)"Micrograph::write: Cannot write file: " + filename);
	}

        MD.setName("general");
        MD.setIsList(true);
        MD.addObject();
	MD.setValue(EMDL_IMAGE_SIZEX, width);
        MD.setValue(EMDL_IMAGE_SIZEY, height);
        MD.setValue(EMDL_IMAGE_SIZEZ, nFrame);
        MD.setValue(EMDL_MICROGRAPH_MOVIE_NAME, fnMovie);
	if (fnGain != "") {
		MD.setValue(EMDL_MICROGRAPH_GAIN_NAME, fnGain);
	}
	MD.write(fh);

	MD.clear();
	MD.setName("global_shift");
	for (int frame = 0; frame < nFrame; frame++) {
		MD.addObject();
		MD.setValue(EMDL_MICROGRAPH_FRAME_NUMBER, frame + 1); // make 1-indexed
		MD.setValue(EMDL_MICROGRAPH_SHIFT_X, globalShiftX[frame]);
		MD.setValue(EMDL_MICROGRAPH_SHIFT_Y, globalShiftY[frame]);
	}
	MD.write(fh);

	fh.close();
}

// Get motion at frame and (x, y)
int Micrograph::getShiftAt(int frame, RFLOAT x, RFLOAT y, RFLOAT &shiftx, RFLOAT &shifty) {
	if (model != NULL) {
		model->getShiftAt(frame, x, y, shiftx, shifty);
	} else {
		shiftx = 0;
		shifty = 0;
	}

	// frame is 1-indexed!
	shiftx += globalShiftX[frame - 1];
	shifty += globalShiftY[frame - 1];
}

void Micrograph::setGlobalShift(int frame, RFLOAT shiftx, RFLOAT shifty) {
	if (frame <= 0 || frame > nFrame) {
		std::cout << "Frame: " << frame << " nFrame: " << nFrame << std::endl;
		REPORT_ERROR("Micrograph::setGlobalShift() frame out of range");
	}

	frame--;
	globalShiftX[frame] = shiftx;
	globalShiftY[frame] = shifty;
}

// TODO: Implement
int MotionModel::getShiftAt(int frame, RFLOAT x, RFLOAT y, RFLOAT &shiftx, RFLOAT &shifty) {
	shiftx = 0.0;
	shifty = 0.0;

	return 0;
}

void MotionModel::fit() {
	REPORT_ERROR("MotionModel::fit() Not implemented!");
}
