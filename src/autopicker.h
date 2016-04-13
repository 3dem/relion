/*
 * autopicker.h
 *
 *  Created on: Sep 18, 2013
 *      Author: "Sjors H.W. Scheres"
 */

#ifndef AUTOPICKER_H_
#define AUTOPICKER_H_
#include "src/image.h"
#include "src/multidim_array.h"
#include "src/metadata_table.h"
#include "src/projector.h"
#include "src/ctf.h"
#include "src/fftw.h"
#include "src/time.h"
#include "src/mask.h"
#include "src/macros.h"
#include "src/helix.h"

//#define OUTPUT_MEAN_MAP_ONLY 1
//#define OUTPUT_STDDEV_MAP_ONLY 2
//#define OUTPUT_BOTH_MEAN_AND_STDDEV_MAPS 3

class ccfPixel
{
public:
	RFLOAT x;
	RFLOAT y;
	RFLOAT fom;
	//RFLOAT psi;
	ccfPixel()
	{
		x = y = fom = (-1.);
		//psi = -1.;
	};
	//ccfPixel(int _x, int _y, RFLOAT _fom, RFLOAT _psi)
	//{
	//	x = _x; y = _y; fom = _fom; psi = _psi;
	//};
	ccfPixel(RFLOAT _x, RFLOAT _y, RFLOAT _fom)
	{
		x = _x; y = _y; fom = _fom;
	};
	bool operator<(const ccfPixel& b) const
	{
		return (fom < b.fom);
	};
};

class ccfPeak
{
public:
	int id;
	RFLOAT x;
	RFLOAT y;
	RFLOAT r;
	RFLOAT area_percentage;
	RFLOAT fom_max;
	int ref;
	RFLOAT psi;
	RFLOAT dist;
	RFLOAT fom_thres;
	int nr_peak_pixel;
	std::vector<ccfPixel> ccf_pixel_list;

	bool isValid() const
	{
		// Invalid parameters
		if ( (r < 0.) || (area_percentage < 0.) || (ccf_pixel_list.size() < 1) )
			return false;
		// TODO: check ccf values in ccf_pixel_list?
		for (int id = 0; id < ccf_pixel_list.size(); id++)
		{
			if (ccf_pixel_list[id].fom > fom_thres)
				return true;
		}
		return false;
	};
	void clear()
	{
		id = ref = nr_peak_pixel = -1;
		x = y = r = area_percentage = fom_max = psi = dist = fom_thres = (-1.);
		ccf_pixel_list.clear();
	};
	ccfPeak()
	{
		clear();
	};
	~ccfPeak()
	{
		clear();
	};
	bool operator<(const ccfPeak& b) const
	{
		if (fabs(r - b.r) < 0.01)
		{
			return (fom_max < b.fom_max);
		}
		return (r < b.r);
	};
	bool refresh()
	{
		RFLOAT x_avg, y_avg;
		int nr_valid_pixel;

		area_percentage = (-1.);

		if (ccf_pixel_list.size() < 1)
			return false;

		fom_max = (-99.e99);
		nr_valid_pixel = 0;
		x_avg = y_avg = 0.;
		for (int id = 0; id < ccf_pixel_list.size(); id++)
		{
			if (ccf_pixel_list[id].fom > fom_thres)
			{
				nr_valid_pixel++;

				if (ccf_pixel_list[id].fom > fom_max)
					fom_max = ccf_pixel_list[id].fom;

				x_avg += ccf_pixel_list[id].x;
				y_avg += ccf_pixel_list[id].y;
			}
		}
		nr_peak_pixel = nr_valid_pixel;

		if (nr_valid_pixel < 1)
			return false;

		x = x_avg / (RFLOAT)(nr_valid_pixel);
		y = y_avg / (RFLOAT)(nr_valid_pixel);
		area_percentage = (RFLOAT)(nr_valid_pixel) / ccf_pixel_list.size();

		return true;
	};
};



struct Peak
{
	int x;
	int y;
	int ref;
	RFLOAT psi;
	RFLOAT fom;
	RFLOAT relative_fom;
};

class AutoPicker
{
public:

	// For GPU-acceleration
	void* cudaPicker;
	// Available memory (in Gigabyte)
	RFLOAT available_memory;
	RFLOAT available_gpu_memory;
	RFLOAT requested_gpu_memory;

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Input & Output rootname
	FileName fn_in, fn_ref, fns_autopick, fn_odir, fn_out;

	// Pixel size (for low-pass filter and particle diameter)
	RFLOAT angpix;

	// Metadata of the micrographs
	MetaDataTable MDmic;

	// Particle diameter (in Angstroms)
	RFLOAT particle_diameter;
	int particle_radius2, decrease_radius;

	// Low pass filter cutoff (in Angstroms)
	RFLOAT lowpass;

	// High pass filter cutoff (in Angstroms)
	RFLOAT highpass;

	// Original size of the reference images
	int particle_size;

	// Dimension of the filtered image
	int current_size;

	// Vector with all original reference images
	std::vector<MultidimArray<RFLOAT> > Mrefs;

	// FTs of the reference images (either for autopicking or for feature calculation)
	std::vector<Projector > PPref;

	///// Autopicking stuff

	// Re-read precalculated best_localCCF and SPI arrays from disc
	bool do_read_fom_maps;

	// Write precalculated best_localCCF and SPI arrays to disc
	bool do_write_fom_maps;

	/// Only autopick those micrographs for which the coordinate file does not yet exist
	bool do_only_unfinished;

	// All micrographs to autopick from
	std::vector<FileName> fn_micrographs;

	// Original size of the micrographs
	int micrograph_size, micrograph_xsize, micrograph_ysize, micrograph_minxy_size;

	// decreased size micrograph
	int workSize;
	float workFrac;

	// Is density in micrograph inverted wrt templates?
	bool do_invert;

	// Correct the references for CTF effects?
	bool do_ctf;

	// use GPU hardware?
	bool do_gpu;

	// Which GPU devices to use?
	std::string gpu_ids;

	// Keep the CTFs unchanged until the first peak?
	bool intact_ctf_first_peak;

	// Are the templates 2D helical segments? If so, in-plane rotation angles (psi) are estimated for the references.
	bool autopick_helical_segments;

	RFLOAT helical_tube_curvature_factor_max;

	RFLOAT helical_tube_diameter;

	RFLOAT helical_tube_length_min;

	// Apart from keeping particle_size/2 away from the sides, should we exclude more? E.g. to get rid of Polara bar code?
	int autopick_skip_side;

	// In-plane rotational sampling (in degrees)
	RFLOAT psi_sampling;

	// Fraction of expected probability ratio to consider as peaks
	RFLOAT min_fraction_expected_Pratio;

	// Number of Angstroms any 2 particle peaks need to be apart
	RFLOAT min_particle_distance;

	// Maximum standard deviation of the noise prior to normalization to pick peaks from
	RFLOAT max_stddev_noise;

	// Removal of outlier pixel values
	RFLOAT outlier_removal_zscore;

	// Size of the downsize micrographs for autopicking
	int downsize_mic;

	// Number of non-zero pixels in the circular mask, and of its inverse (for background normalisation in do_diff2)
	int nr_pixels_circular_mask, nr_pixels_circular_invmask;

	// Array with Fourier-transform of the (circular) mask, and of its inverse
	MultidimArray<Complex > Fmsk, Finvmsk;

public:

	AutoPicker():
		available_memory(0),
		available_gpu_memory(0),
		requested_gpu_memory(0)
	{}

	// Read command line arguments
	void read(int argc, char **argv);

	// Print usage instructions
	void usage();

	// Initialise some general stuff after reading
	void initialise();

	// Set device-affinity
	int deviceInitialise();

	// General function to decide what to do
	void run();

	void pickCCFPeaks(
			const MultidimArray<RFLOAT>& Mccf,
			const MultidimArray<int>& Mclass,
			RFLOAT threshold_value,
			int peak_r_min,
			RFLOAT particle_diameter_pix,
			std::vector<ccfPeak>& ccf_peak_list,
			MultidimArray<RFLOAT>& Mccfplot,
			int micrograph_maxxy_size,
			int micrograph_minxy_size,
			int skip_side);

	void extractHelicalTubes(
			std::vector<ccfPeak>& ccf_peak_list,
			std::vector<std::vector<ccfPeak> >& tube_coord_list,
			std::vector<RFLOAT>& tube_len_list,
			std::vector<std::vector<ccfPeak> >& tube_track_list,
			RFLOAT particle_diameter_pix,
			RFLOAT curvature_factor_max,
			RFLOAT interbox_distance_pix,
			RFLOAT tube_diameter_pix);

	void exportHelicalTubes(
			const MultidimArray<RFLOAT>& Mccf,
			MultidimArray<RFLOAT>& Mccfplot,
			const MultidimArray<int>& Mclass,
			std::vector<std::vector<ccfPeak> >& tube_coord_list,
			std::vector<std::vector<ccfPeak> >& tube_track_list,
			std::vector<RFLOAT>& tube_len_list,
			FileName& fn_mic_in,
			FileName& fn_star_out,
			RFLOAT particle_diameter_pix,
			RFLOAT tube_length_min_pix,
			int micrograph_maxxy_size,
			int micrograph_xsize,
			int micrograph_ysize,
			int skip_side);

	void autoPickOneMicrograph(FileName &fn_mic);

	// Get the output coordinate filename given the micrograph filename
	FileName getOutputRootName(FileName fn_mic);
	// Uses Roseman2003 formulae to calculate stddev under the mask through FFTs
	// The FFTs of the micrograph (Fmic), micrograph-squared (Fmic2) and the mask (Fmsk) need to be provided at downsize_mic
	// The putput (Mstddev) will be at (binned) micrograph_size
	void calculateStddevAndMeanUnderMask(
			const MultidimArray<Complex > &Fmic,
			const MultidimArray<Complex > &Fmic2,
			MultidimArray<Complex > &Fmsk,
			int nr_nonzero_pixels_mask,
			MultidimArray<RFLOAT> &Mstddev,
			MultidimArray<RFLOAT> &Mmean);

	// Peak search for all pixels above a given threshold in the map
	void peakSearch(const MultidimArray<RFLOAT> &Mccf, const MultidimArray<RFLOAT> &Mpsi, const MultidimArray<RFLOAT> &Mstddev, int iref, int skip_side, std::vector<Peak> &peaks, float scale);

	// Now prune the coordinates: within min_particle_distance: all peaks are the same cluster
	// From each cluster, take the single peaks with the highest ccf
	// If then, there is another peaks at a distance of at least min_particle_distance: take that one as well, and so forth...
	void prunePeakClusters(std::vector<Peak> &peaks, int min_distance, float scale);

	// Only keep those peaks that are at the given distance apart from each other
	void removeTooCloselyNeighbouringPeaks(std::vector<Peak> &peaks, int min_distance, float scale);

};


#endif /* AUTOPICKER_H_ */
