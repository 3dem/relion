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


#ifndef SRC_HELIX_INIMODEL2D_H_
#define SRC_HELIX_INIMODEL2D_H_

#include "src/parallel.h"
#include "src/time.h"
#include "src/filename.h"
#include "src/metadata_table.h"
#include "src/image.h"
#include "src/euler.h"
#include "src/backprojector.h"
#include "src/transformations.h"
#include "src/fftw.h"


class HelixAlignerModel
{

public:

	// The reference images
	std::vector<MultidimArray<RFLOAT> > Aref;

	// The reconstructed xy-slice
	std::vector<MultidimArray<RFLOAT> > Arec;

	// The running sum reference images
	std::vector<MultidimArray<RFLOAT> > Asum;

	// The total sum reference images
	std::vector<MultidimArray<RFLOAT> > Asumw;

	// Total contribution to each class
	std::vector<RFLOAT> pdf;

	// Empty constructor
	HelixAlignerModel() {};

	// Destructor
	~HelixAlignerModel() {clear(); };

	// Cleaning up
	void clear();

	// To initialise model
	void initialise(int nr_classes, int ydim, int xdim);

	// To initialise the sums to zero
	void initZeroSums();

};



class HelixAligner
{

public:

	// I/O Parser
	IOParser parser;

	//Input images
	FileName fn_imgs;

	// Output rootname
	FileName fn_out;

	// Number of iterations, classes
	int nr_iter, nr_classes;

	// Pixel size of input images
	float angpix;

	// random
	int random_seed;

	// Diameter of circular mask within to extract images
	float mask_diameter, mask_radius_pix;

	// Maximum resolution to be taken into account (approximate, as adjusted to accommodate exact crossover_distance
	float maxres;

	// How many pixels away from the target resolution to search for optimal downscaled pixel size?
	int search_size;

	// Distance in Angstroms between 2 cross-overs (i.e. 180 degrees of twist)
	float crossover_distance;

	// Height in Angstroms to be taken into account
	float height;

	// How much smearing to apply to the initial reference (to start with a smoother reference along the helical axis)
	int max_smear;

	// How many pixels to search up and down?
	int max_shift;
	RFLOAT max_shift_A;

	// How many degrees to rotate?
	RFLOAT max_rotate;

	// Rotation step
	RFLOAT step_rotate;

	// The model to be refined
	HelixAlignerModel model;

	// Input micrographs
	FileName fn_mics;
	MetaDataTable MDmics;

	// STAR file with all (selected) micrographs, the suffix of the coordinates files, and the directory where the coordinate files are
	FileName fn_coord_suffix, fn_coord_dir ;

	// Width of images to be extracted
	int extract_width;

	// Filename of initial 2D reconstruction for model
	FileName fn_inimodel;

	// Only make 3d
	bool do_only_make_3d;

	// Symmetry order (Cn)
	int symmetry;

	// Number of openMP threads
	int nr_threads;

private:

	// Size of the original and downscaled images
	int ori_size, down_size;

	// Size of the rectangle
	int xrect, yrect;

	// X-Size of the images being placed inside the rectangle
	int ximg;

	// Downsized pixel size
	float down_angpix;

	// Verbosity
	int verb;

	// Pre-calculated Gaussian weight vector
	MultidimArray<RFLOAT> weight;

	// Pre-read (rotated versions of) all Xrect of the (downscaled) images into RAM
	std::vector<std::vector<MultidimArray<RFLOAT> > > Xrects;

	// Foroptimal orientation control
	std::vector<RFLOAT> psis, ori_psis, ori_yoffs;
	MetaDataTable MD;


public:

	// Empty constructor
	HelixAligner() {};

	// Destructor
	~HelixAligner() {clear(); };

	// Usage
	void usage();

	// Cleaning up
	void clear();

	void parseInitial(int argc, char **argv);

	void initialise();

	// Read in all the images
	void readImages();

	// 22 June 2017: extract helices from start-end coordinates in micrographs
	void getHelicesFromMics();

	// Initialise classes randomly
	void initialiseClasses();

	void expectationOneParticleNoFFT(long int ipart);

	void expectation();

	void maximisation();

	void reconstruct2D(int iclass);

	void writeOut(int iter);

	void reconstruct3D();

	// Run multiple iterations
	void run();



};


#endif /* SRC_HELIX_INIMODEL2D_H_ */
