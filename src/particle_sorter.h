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

#ifndef PARTICLE_SORTER_H_
#define PARTICLE_SORTER_H_

#include "src/image.h"
#include "src/multidim_array.h"
#include "src/metadata_table.h"
#include "src/projector.h"
#include "src/ctf.h"
#include <src/fftw.h>
#include <src/time.h>

#define FEATURE_DF_AVG      0
#define FEATURE_DF_SIG      1
#define FEATURE_DF_SKW      2
#define FEATURE_DF_KRT      3
#define FEATURE_DF_QUADSIG  4
#define FEATURE_DF_ROTFOURCORR 5
#define NR_FEATURES 6

#define WIDTH_FMASK_EDGEB 2


class ParticleSorter
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Input & Output rootname
	FileName fn_in, fn_ref, fn_out;

	// Input metadata
	MetaDataTable MDin;

	// Pixel size (for low-pass filter and particle diameter)
	double angpix;

	// Particle diameter (in Angstroms)
	double particle_diameter;
	int particle_radius2;

	// Low pass filetr cutoff (in Angstroms)
	double lowpass;

	// Original size of the reference images
	int particle_size;

	// Dimension of the filtered image
	int current_size;

	// Minimum Z-value to count in the sorting
	double min_z;

	// Vector with all original reference images
	std::vector<MultidimArray<double> > Mrefs;

	// FTs of the reference images for feature calculation
	std::vector<Projector > PPref;

	// Feature values for all input images
	MultidimArray<double> features;

	// Is density in micrograph inverted wrt templates?
	bool do_invert;

	// Correct the references for CTF effects?
	bool do_ctf;

	// Keep the CTFs unchanged until the first peak?
	bool intact_ctf_first_peak;

public:
	// Read command line arguments
	void read(int argc, char **argv);

	// Print usage instructions
	void usage();

	// General function to decide what to do
	void run();

	// Initialise some general stuff after reading
	void initialise();

	void calculateFeaturesOneParticle(long int ipart);

protected:

	void normaliseFeatures();

	// Write out (for now in libsvm format)
	void writeFeatures();

	void calculateStatsOneImage(MultidimArray<double> &img,
			double &mean, double &stddev, double &skew, double &kurt, double &quadrant_stddev);

	double rotationalSymmetryFourierTransform(MultidimArray<Complex > &Fimg);

};


#endif /* PARTICLE_SORTER_H_ */
