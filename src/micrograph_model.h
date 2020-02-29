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
#include "src/multidim_array.h"

enum MotionModelVersion {
	MOTION_MODEL_NULL = 0,
	MOTION_MODEL_THIRD_ORDER_POLYNOMIAL = 1,
};

class MotionModel
{
public:

	virtual ~MotionModel(){}

	// Fit model based on observations
	virtual void fit() = 0;

	virtual void read(std::ifstream &fh, std::string block_name) = 0;
	virtual void write(std::ostream &fh, std::string block_name) = 0;
	virtual int getModelVersion() const = 0;

	// Get motion at frame and (x, y);
	// NB: differences from Micrograph::getShiftAt!
	// - frame is 0-indexed 
	// - (x, y) are normalised pixels (i.e. unbinned_pixel_x / width - 0.5)
	virtual int getShiftAt(RFLOAT frame, RFLOAT x, RFLOAT y, RFLOAT &shiftx, RFLOAT &shifty) const = 0;

	virtual MotionModel* clone() const = 0;
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

	MotionModel* clone() const;
};

class Micrograph
{
public:
	// When you add a new field, don't forget to update the copy constructor!

	bool ready;
	static const RFLOAT NOT_OBSERVED;
	RFLOAT angpix, voltage, dose_per_frame, pre_exposure;
	FileName fnDefect;

	int first_frame; // First frame for local motion model. 1-indexed.
	MotionModel *model;

	// Local trajectories (not read from STAR files)
	std::vector<RFLOAT> localShiftX, localShiftY, localFitX, localFitY, patchX, patchY, patchZ, patchW, patchH;
	std::vector<int> hotpixelX, hotpixelY;

	// Default constructor
	Micrograph();

	// Copy-constructor
	Micrograph(const Micrograph& m);

	// Create from a movie or a STAR file
	Micrograph(FileName filename, FileName fnGain="", RFLOAT binning=1.0, int eer_upsampling=-1, int eer_grouping=-1);

	~Micrograph();

	Micrograph& operator = (const Micrograph& m);

	// Write micrograph model from a STAR file
	void write(FileName filename);

	// Get gain reference file name
	FileName getGainFilename() const;

	// Get binning factor
	RFLOAT getBinningFactor() const;

	// Get original movie name
	FileName getMovieFilename() const;

	int getWidth() const;
	int getHeight() const;
	int getNframes() const;

	// Get shift vector at (x, y, frame)
	// frame is 1-indexed
	// (x, y) are normalised coordinate (i.e. pixel_x / width, pixel_y / height) when normalise=false (default)
	//            unbinned pixels in the original movie when normalise=true
	// (shiftx, shifty) are always UNBINNED pixels in the original movie.
	// Returns non zero if failed (e.g. not observed).
	int getShiftAt(RFLOAT frame, RFLOAT x, RFLOAT y, RFLOAT &shiftx, RFLOAT &shifty, bool use_local=true, bool normalise=false) const;

	// Set global shift for frame; frame is 1-indexed
	// (shiftx, shifty) is UNBINNED pixels in the original movie
	void setGlobalShift(int frame, RFLOAT shiftx, RFLOAT shifty);

	// Fills a pixel mask where defect and hot pixels are true
	void fillDefectAndHotpixels(MultidimArray<bool> &mask) const;

	int getEERUpsampling() const;
	int getEERGrouping() const;

private:
	int width, height, n_frames;
	RFLOAT binning;
	FileName fnGain;
	FileName fnMovie;

	int eer_upsampling, eer_grouping;

	std::vector<RFLOAT> globalShiftX, globalShiftY;

	// Read micrograph model from a STAR file
	void read(FileName filename, bool read_hotpixels=true);

	// Set target movie file
	void setMovie(FileName fnMovie, FileName fnGain="", RFLOAT binning=1.0);

	void clear();
	void clearFields();
	void copyFieldsFrom(const Micrograph& m);

	void checkReadyFlag(std::string origin) const;
};

#endif /* MICROGRAPH_MODEL_H_ */
