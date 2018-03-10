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
#include "src/matrix1d.h"

enum MotionModelVersion {
	MOTION_MODEL_NULL,
	MOTION_MODEL_THIRD_ORDER_POLYNOMIAL,
};

class MotionModel
{
public:
	// Fit model based on observations
	virtual void fit() = 0;

	virtual void read(std::ifstream &fh, std::string block_name) = 0;
	virtual void write(std::ostream &fh, std::string block_name) = 0;
	virtual int getModelVersion() const = 0;

	// Get motion at frame and (x, y); frame is 0-indexed (NOT like Micrograph::getShiftAt)
	virtual int getShiftAt(RFLOAT frame, RFLOAT x, RFLOAT y, RFLOAT &shiftx, RFLOAT &shifty) const = 0;
};

class ThirdOrderPolynomialModel: public MotionModel {
public:
	static const int NUM_COEFFS_PER_DIM;

	Matrix1D <RFLOAT> coeffX, coeffY;

	void fit() {
		REPORT_ERROR("Not implemented yet.");
	}

	void read(std::ifstream &fh, std::string block_name);
	void write(std::ostream &fh, std::string block_name);

	int getModelVersion() const {
		return MOTION_MODEL_THIRD_ORDER_POLYNOMIAL;
	}

	int getShiftAt(RFLOAT frame, RFLOAT x, RFLOAT y, RFLOAT &shiftx, RFLOAT &shifty) const;
};

class Micrograph
{
public:
	static const RFLOAT NOT_OBSERVED;
	RFLOAT angpix, voltage, dose_per_frame, pre_exposure;
	int first_frame; // First frame for local motion model. 1-indexed.
	MotionModel *model;

	// Local trajectories (not written, not read from STAR files)
	std::vector<RFLOAT> localShiftX, localShiftY, localFitX, localFitY, patchX, patchY, patchZ, patchW, patchH;
	
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
		n_frames = 0;
		first_frame = 0;
		binning = 1;

		angpix = -1;
		voltage = -1;
		dose_per_frame = -1;
		pre_exposure = -1;

		fnMovie = "";
		fnGain = "";

		globalShiftX.resize(0);
		globalShiftY.resize(0);

		localShiftX.resize(0);
		localShiftY.resize(0);
		localFitX.resize(0);
		localFitY.resize(0);
		patchX.resize(0);
		patchY.resize(0);
		patchZ.resize(0);
		patchW.resize(0);
		patchH.resize(0);
	
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
	FileName getGainFilename() const {
		return fnGain;
	}

	// Get binning factor
	RFLOAT getBinningFactor() const {
		return binning;
	}

	// Get original movie name
	FileName getMovieFilename() const {
		return fnMovie;
	}

	int getWidth() const {
		return width;
	}

	int getHeight() const {
		return height;
	}

	int getNframes() const {
		return n_frames;
	}

	// Get shift vector at (x, y, frame); frame is 1-indexed
	// (x, y) and (shiftx, shifty) are UNBINNED pixels in the original movie
	int getShiftAt(RFLOAT frame, RFLOAT x, RFLOAT y, RFLOAT &shiftx, RFLOAT &shifty) const;

	// Set global shift for frame; frame is 1-indexed
	// (shiftx, shifty) is UNBINNED pixels in the original movie
	void setGlobalShift(int frame, RFLOAT shiftx, RFLOAT shifty);

private:
	int width, height, n_frames;
	RFLOAT binning;
	FileName fnGain;
	FileName fnMovie;
	
	std::vector<RFLOAT> globalShiftX, globalShiftY;
};

#endif /* MICROGRAPH_MODEL_H_ */
