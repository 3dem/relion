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
#ifndef MICROGRAPH_MODEL_H_
#define MICROGRAPH_MODEL_H_

#include <vector>
#include "src/filename.h"

class MotionModel
{
public:
	// Fit model based on observations
	virtual void fit();

	// Get motion at frame and (x, y)
	virtual int getShiftAt(int frame, RFLOAT x, RFLOAT y, RFLOAT &shiftx, RFLOAT &shifty);
};

class Micrograph
{
public:
	static const RFLOAT NOT_OBSERVED;
	RFLOAT angpix, voltage, dose_per_frame, pre_exposure;

	// Empty Constructor is not allowed
	Micrograph(); // = delete in C++11

	// Create from a movie or a STAR file
	Micrograph(FileName filename, FileName fnGain="", RFLOAT binning=1.0) {
		model = NULL;

		clear();

		if (filename.getExtension() == "star" && fnGain == "") {
			read(filename);
		} else {
			setMovie(filename, fnGain, binning);
		}
	}

	~Micrograph()
	{
		clear();
	}

	// Initialise
	void clear()
	{
		width = 0;
		height = 0;
		nFrame = 0;
		binning = 1;

		angpix = -1;
		voltage = -1;
		dose_per_frame = -1;
		pre_exposure = -1;

		fnMovie = "";
		fnGain = "";

		globalShiftX.resize(0);
		globalShiftY.resize(0);
	
		if (model != NULL) delete model;
		model = NULL;
	}

	// Read micrograph model from a STAR file
	void read(FileName filename);

	// Write micrograph model from a STAR file
	void write(FileName filename);

	// Set target movie file
	void setMovie(FileName fnMovie, FileName fnGain="", RFLOAT binning=1.0);

	// Get gain reference file name
	FileName getGainFilename() {
		return fnGain;
	}

	// Get binning factor
	RFLOAT getBinningFactor() {
		return binning;
	}

	// Get original movie name
	FileName getMovieFilename() {
		return fnMovie;
	}

	// Get shift vector at (x, y, frame)
	// (x, y) and (shiftx, shifty) are UNBINNED pixels in the original movie
	int getShiftAt(int frame, RFLOAT x, RFLOAT y, RFLOAT &shiftx, RFLOAT &shifty);

	// Set global shift for frame
	// (shiftx, shifty) is UNBINNED pixels in the original movie
	void setGlobalShift(int frame, RFLOAT shiftx, RFLOAT shifty);

private:
	int width, height, nFrame;
	RFLOAT binning;
	FileName fnGain;
	FileName fnMovie;
	
	std::vector<RFLOAT> globalShiftX, globalShiftY;
	MotionModel *model;
};

#endif /* MICROGRAPH_MODEL_H_ */
