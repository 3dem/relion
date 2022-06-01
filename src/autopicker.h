/*
 * autopicker.h
 *
 *  Created on: Sep 18, 2013
 *      Author: "Sjors H.W. Scheres"
 */

#ifndef AUTOPICKER_H_
#define AUTOPICKER_H_
#include <stdlib.h>
#include "src/image.h"
#include "src/multidim_array.h"
#include "src/metadata_table.h"
#include "src/projector.h"
#include "src/healpix_sampling.h"
#include "src/projector.h"
#include <src/jaz/single_particle/obs_model.h>
#include "src/ctf.h"
#include "src/fftw.h"
#include "src/time.h"
#include "src/mask.h"
#include "src/macros.h"
#include "src/helix.h"
#ifdef _CUDA_ENABLED
#include "src/acc/cuda/cuda_mem_utils.h"
#include "src/acc/acc_projector.h"
#include "src/acc/cuda/cuda_settings.h"
#include "src/acc/cuda/cuda_fft.h"
#include "src/acc/cuda/cuda_benchmark_utils.h"
#endif
//#define OUTPUT_MEAN_MAP_ONLY 1
//#define OUTPUT_STDDEV_MAP_ONLY 2
//#define OUTPUT_BOTH_MEAN_AND_STDDEV_MAPS 3

//#define TIMING
class ccfPixel
{
public:
	RFLOAT x, y, fom;
	//RFLOAT x, y, fom, psi;

	ccfPixel() : x(-1.), y(-1.), fom(-1.) {};
	ccfPixel(RFLOAT _x, RFLOAT _y, RFLOAT _fom) : x(_x), y(_y), fom(_fom) {};

	//ccfPixel() : x(-1.), y(-1.), fom(-1.), psi(-1.) {};
	//ccfPixel(RFLOAT _x, RFLOAT _y, RFLOAT _fom, RFLOAT _psi) : x(_x), y(_y), fom(_fom), psi(_psi) {};

	bool operator<(const ccfPixel& b) const { return (fom < b.fom); };
};

class ccfPeak
{
public:
	int id, ref, nr_peak_pixel;
	RFLOAT x, y, r, area_percentage, fom_max, psi, dist, fom_thres;
	std::vector<ccfPixel> ccf_pixel_list;

	void clear();

	ccfPeak() { clear(); };

	~ccfPeak() { clear(); };

	bool isValid() const;

	bool operator<(const ccfPeak& b) const;

	bool refresh();
};

struct Peak
{
	int x, y, ref;
	RFLOAT psi, fom, relative_fom;
};

struct AmyloidCoord
{
	RFLOAT x, y, psi, fom;
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

	// Random seed
	long int random_seed;

	// Input & Output rootname
	FileName fn_in, fn_ref, fns_autopick, fn_odir, fn_out;

	// Pixel size for the micrographs (for low-pass filter and particle diameter)
	RFLOAT angpix;

	// Pixel size for the references (for low-pass filter and particle diameter)
	RFLOAT angpix_ref;

	// Angular sampling rate for projection of 3D reference (hp=0: 60 deg, hp=1: 30 deg; hp=2: 15deg)
	int healpix_order;

	// Symmetry point group for 3D reference
	std::string symmetry;

	// Metadata of the micrographs
	MetaDataTable MDmic;

	// Optics group information
	ObservationModel obsModel;

	// Particle diameter (in Angstroms)
	RFLOAT particle_diameter;
	int particle_radius2, decrease_radius;

	// Maximum diameter for local average density calculation
	RFLOAT max_local_avg_diameter;

	// Low pass filter cutoff (in Angstroms)
	RFLOAT lowpass;

	// High pass filter cutoff (in Angstroms)
	RFLOAT highpass;

	// Original size of the reference images
	int particle_size;

	// Dimension of the filtered image
	int current_size;

	// Padding to use for Projectors
	int padding;

	// Maxmimum value in the Gaussian blob reference
	RFLOAT gauss_max_value;

	// Vector with all original reference images
	std::vector<MultidimArray<RFLOAT> > Mrefs;

	// FTs of the reference images (either for autopicking or for feature calculation)
	std::vector<Projector > PPref;

	// Use Laplacian-of-Gaussian filters instead of template-based picking
	bool do_LoG;

	// Diameter for features to be detected by the LoG filter
	RFLOAT LoG_min_diameter, LoG_max_diameter, LoG_neighbour_fudge;

	// How many times the LoG_max_diameter is searched?
	RFLOAT LoG_max_search;

	// How many sigma to adjust the FOM threshold?
	RFLOAT LoG_adjust_threshold, LoG_upper_limit;

	// Input signal is white
	bool LoG_invert, LoG_use_ctf;

	// Vector with all LoG filter FFTs
	std::vector<MultidimArray<Complex> > FT_LoGs;

	// Vector with all diameters to be sampled
	std::vector<RFLOAT> diams_LoG;

	/// Topaz wrappers
	// Use topaz train or topaz extract instead of template-based picking
	bool do_topaz_train, do_topaz_extract;

	// Topaz threshold for picking
	RFLOAT topaz_threshold;

	// Write out intermediate plots for topaz helical picking
	bool do_topaz_plot;

	// Expected number of particles per micrograph
	int topaz_nr_particles;

	// Topaz downscale factor
	int topaz_downscale;

	// Number of topaz workers for training
	int topaz_workers;

	// Topaz command executable
	FileName fn_topaz_exe;

	// sh executable
	FileName fn_shell;

	// Topaz saved model for use in extract
	FileName topaz_model;

	// Topaz particle radius for use in extract
	int topaz_radius;

	// GPU Device ID
	int device_id = -1;

	// Filename for picks or particles star file to be used for topaz training
	FileName topaz_train_picks, topaz_train_parts;

	// Metadata table for training picks
	MetaDataTable MDtrain;

	// Default ratio of picks in the test validation set
	RFLOAT topaz_test_ratio;

	// Number of topaz workers
	int nr_topaz_threads;

	// Other arguments to be passed to topaz
	FileName topaz_additional_args;

	//// Specific amyloid picker
	bool do_amyloid;

	/// Maximum psi-angle difference in subsequent amyloid segments (in degrees)
	RFLOAT amyloid_max_psidiff;

	///// Autopicking stuff

	// Re-read precalculated best_localCCF and SPI arrays from disc
	bool do_read_fom_maps;

	// Write precalculated best_localCCF and SPI arrays to disc
	bool do_write_fom_maps;
	// We impose a limit to not write insane number of images by mistake, but you can override through --fom_override
	bool no_fom_limit;

	/// Only autopick those micrographs for which the coordinate file does not yet exist
	bool do_only_unfinished;

	// Is there any work to be done?
	bool todo_anything;

	// All micrographs to autopick from
	std::vector<FileName> fn_micrographs, fn_ori_micrographs;

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

	// Extra padding around the micrographs, of this many pixels
	int extra_padding;

	// In-plane rotational sampling (in degrees)
	RFLOAT psi_sampling;

	// Fraction of expected probability ratio to consider as peaks
	RFLOAT min_fraction_expected_Pratio;

	// Number of Angstroms any 2 particle peaks need to be apart
	RFLOAT min_particle_distance;

	// Maximum standard deviation of the noise prior to normalization to pick peaks from
	RFLOAT max_stddev_noise;

	// Minimum average background density of the noise  to pick peaks from
	RFLOAT min_avg_noise;

	// Removal of outlier pixel values
	RFLOAT outlier_removal_zscore;

	// Size of the downsize micrographs for autopicking
	int downsize_mic;

	// Number of non-zero pixels in the circular mask, and of its inverse (for background normalisation in do_diff2)
	int nr_pixels_circular_mask, nr_pixels_avg_mask, nr_pixels_circular_invmask;

	// Array with Fourier-transform of the inverse of the (circular) mask
	MultidimArray<Complex > Finvmsk;

	// Array with Fourier-transform of the mask to calculate average density
	MultidimArray<Complex > Favgmsk;

	// Perform optimisation of the scale factor?
	bool do_optimise_scale;

#ifdef TIMING
    Timer timer;
	int TIMING_A0, TIMING_A1, TIMING_A2, TIMING_A3, TIMING_A4, TIMING_A5, TIMING_A6, TIMING_A7, TIMING_A8, TIMING_A9;
	int TIMING_B1, TIMING_B2, TIMING_B3, TIMING_B4, TIMING_B5, TIMING_B6, TIMING_B7, TIMING_B8, TIMING_B9;
#endif

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
	void initialise(int rank = 0);

	// Set device-affinity
	void deviceInitialise();

	// General function to decide what to do
	void run();

	// Make a PDF file with plots of numbers of particles per micrograph, average FOMs etc
	void generatePDFLogfile();

	std::vector<AmyloidCoord> findNextCandidateCoordinates(AmyloidCoord &mycoord, std::vector<AmyloidCoord> &circle,
			RFLOAT threshold_value, RFLOAT max_psidiff, int skip_side, float scale,
			MultidimArray<RFLOAT> &Mccf, MultidimArray<RFLOAT> &Mpsi);

	AmyloidCoord findNextAmyloidCoordinate(AmyloidCoord &mycoord, std::vector<AmyloidCoord> &circle, RFLOAT threshold_value,
			RFLOAT max_psidiff, RFLOAT amyloid_diameter_pix, int skip_side, float scale,
			MultidimArray<RFLOAT> &Mccf, MultidimArray<RFLOAT> &Mpsi);

	void pickAmyloids(
			MultidimArray<RFLOAT>& Mccf,
			MultidimArray<RFLOAT>& Mpsi,
			MultidimArray<RFLOAT>& Mstddev,
			MultidimArray<RFLOAT>& Mavg,
			RFLOAT threshold_value,
			RFLOAT max_psidiff,
			FileName& fn_mic_in,
			FileName& fn_star_out,
			RFLOAT amyloid_width,
			int skip_side, float scale);

	void pickCCFPeaks(
			const MultidimArray<RFLOAT>& Mccf,
			const MultidimArray<RFLOAT>& Mstddev,
			const MultidimArray<RFLOAT>& Mavg,
			const MultidimArray<int>& Mclass,
			RFLOAT threshold_value,
			int peak_r_min,
			RFLOAT particle_diameter_pix,
			std::vector<ccfPeak>& ccf_peak_list,
			MultidimArray<RFLOAT>& Mccfplot,
			int skip_side, float scale);

	void extractHelicalTubes(
			std::vector<ccfPeak>& ccf_peak_list,
			std::vector<std::vector<ccfPeak> >& tube_coord_list,
			std::vector<RFLOAT>& tube_len_list,
			std::vector<std::vector<ccfPeak> >& tube_track_list,
			RFLOAT particle_diameter_pix,
			RFLOAT curvature_factor_max,
			RFLOAT interbox_distance_pix,
			RFLOAT tube_diameter_pix, float scale);

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
			int skip_side, float scale);

	MetaDataTable getMDtrainFromParticleStar(MetaDataTable &MDparts, ObservationModel &obsModel);
	MetaDataTable readTopazCoordinates(FileName fn_coord, int _topaz_downscale = 1);

	void preprocessMicrographTopaz(FileName fn_in, FileName fn_out, int downscale);
	void trainTopaz();
	void autoPickTopazOneMicrograph(FileName &fn_mic, int rank = 0);
	void autoPickLoGOneMicrograph(FileName &fn_mic, long int imic);
	void autoPickOneMicrograph(FileName &fn_mic, long int imic);

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
	void peakSearch(const MultidimArray<RFLOAT> &Mccf, const MultidimArray<RFLOAT> &Mpsi,
			const MultidimArray<RFLOAT> &Mstddev, const MultidimArray<RFLOAT> &Mmean,
			int iref, int skip_side, std::vector<Peak> &peaks, float scale);

	// Now prune the coordinates: within min_particle_distance: all peaks are the same cluster
	// From each cluster, take the single peaks with the highest ccf
	// If then, there is another peaks at a distance of at least min_particle_distance: take that one as well, and so forth...
	void prunePeakClusters(std::vector<Peak> &peaks, int min_distance, float scale);

	// Only keep those peaks that are at the given distance apart from each other
	void removeTooCloselyNeighbouringPeaks(std::vector<Peak> &peaks, int min_distance, float scale);

#define LARGEST_ACCEPTABLE_PRIME 43

	int largestPrime(int query);

	int getGoodFourierDims(int requestedSizeRealX, int lim);
};


#endif /* AUTOPICKER_H_ */
