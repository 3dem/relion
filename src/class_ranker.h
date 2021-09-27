/***************************************************************************
 *
 * Author: "Liyi Dong and Sjors H.W. Scheres"
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

#ifndef CLASS_RANKER_H_
#define CLASS_RANKER_H_

#include <stack>
#include <unistd.h>
#include <limits.h>
#include <fstream>
#include "src/ml_optimiser.h"

static float feature_normalization_local_ps_mean=0., feature_normalization_local_ps_stddev=0.;
static float feature_normalization_local_ss_mean=0., feature_normalization_local_ss_stddev=0.;
static float feature_normalization_local_rsi_mean=0., feature_normalization_local_rsi_stddev=0.;


static std::vector<RFLOAT> feature_normalization_global_mean{
								4.352,
								3.871,
								0.001,
								0.061,
								-0.199,
								0.424,
								0.744,
								-0.380,
								0.292,
								0.024,
								0.047,
								0.116,
								0.638,
								-1.112,
								2.990
							},
			    feature_normalization_global_stddev{
								4.759,
								3.019,
								0.002,
								0.088,
								0.283,
								0.629,
								0.585,
								0.378,
								0.485,
								0.075,
								0.075,
								0.182,
								0.429,
								1.087,
								1.758
							},
			    feature_normalization_global_granulo_mean{
								837.576,
								726.097,
								602.175,
								483.469,
								328.014,
								196.059
							},
			     feature_normalization_global_granulo_stddev{
								1958.823,
								2845.330,
								3605.858,
								4001.972,
								4399.578,
								4709.751
							};

// This contains 4 moments for an image
class moments
{

public:

	RFLOAT sum, mean, stddev, skew, kurt;
    moments(): sum(0), mean(0), stddev(0), skew(0), kurt(0){}

};

// This defines all features for a single class
class classFeatures
{

public:

	// Class-wise features
    FileName name;
    int class_index;
    int is_selected;

    int protein_area, solvent_area;
    RFLOAT class_distribution, accuracy_rotation, accuracy_translation, estimated_resolution, particle_nr;
    RFLOAT class_score, relative_signal_intensity, edge_signal, scattered_signal, weighted_resolution, relative_resolution;
    RFLOAT lowpass_filtered_img_avg, lowpass_filtered_img_stddev, lowpass_filtered_img_minval, lowpass_filtered_img_maxval;
    RFLOAT CAR; // protein mask's circumference
    moments circular_mask_moments, ring_moments, inner_circle_moments, protein_moments, solvent_moments;

    // SHWS 15072020: subimages for CNN
    MultidimArray<RFLOAT> subimages;

    // Granularity features only for historical development reasons
    std::vector<RFLOAT> lbp, lbp_p, lbp_s, haralick_p, haralick_s, zernike_moments, granulo;
    RFLOAT total_entropy, protein_entropy, solvent_entropy;

    classFeatures(): name(""),
			class_index(0),
    		is_selected(0),
			particle_nr(0.),
    		class_distribution(0.),
			accuracy_rotation(999.), //? Double check
			accuracy_translation(999.), //? Double check
			estimated_resolution(999.0),
			weighted_resolution(999.0),
			relative_resolution(999.0),
			class_score(-1),
			relative_signal_intensity(0.),
			edge_signal(-1.),
			scattered_signal(-1.),
			lowpass_filtered_img_avg(0.),
			lowpass_filtered_img_stddev(0.),
			lowpass_filtered_img_minval(0.),
			lowpass_filtered_img_maxval(0.),
			total_entropy(0.),
			protein_entropy(0.),
			solvent_entropy(0.),
			protein_area(0),
			solvent_area(0),
                        CAR(0.)
    {
    }

    // Destructor needed for work with vectors
	~classFeatures() {}

	// Copy constructor needed for work with vectors
    classFeatures(classFeatures const& copy)
	{
    	name = copy.name;
    	class_index = copy.class_index;
    	is_selected = copy.is_selected;
    	class_distribution = copy.class_distribution;
    	accuracy_rotation = copy.accuracy_rotation;
    	accuracy_translation= copy.accuracy_translation;
    	estimated_resolution= copy.estimated_resolution;
    	particle_nr = copy.particle_nr;
    	class_score = copy.class_score;
    	relative_signal_intensity = copy.relative_signal_intensity;
    	edge_signal = copy.edge_signal;
    	scattered_signal = copy.scattered_signal;
    	weighted_resolution = copy.weighted_resolution;
    	relative_resolution = copy.relative_resolution;
    	lowpass_filtered_img_avg = copy.lowpass_filtered_img_avg;
    	lowpass_filtered_img_stddev = copy.lowpass_filtered_img_stddev;
    	lowpass_filtered_img_minval = copy.lowpass_filtered_img_minval;
    	lowpass_filtered_img_maxval = copy.lowpass_filtered_img_maxval;
    	circular_mask_moments = copy.circular_mask_moments;
    	ring_moments = copy.ring_moments;
    	inner_circle_moments = copy.inner_circle_moments;
    	protein_area = copy.protein_area;
    	solvent_area = copy.solvent_area;
    	protein_moments = copy.protein_moments;
    	solvent_moments = copy.solvent_moments;
    	lbp = copy.lbp;
    	lbp_p = copy.lbp_p;
    	lbp_s = copy.lbp_s;
    	haralick_p = copy.haralick_p;
    	haralick_s = copy.haralick_s;
    	zernike_moments = copy.zernike_moments;
	    granulo = copy.granulo;
	    total_entropy = copy.total_entropy;
	    protein_entropy = copy.protein_entropy;
	    solvent_entropy = copy.solvent_entropy;
	    subimages = copy.subimages;
	    CAR = copy.CAR;
	}

	// Define assignment operator in terms of the copy constructor
    classFeatures& operator=(classFeatures const& copy)
	{
    	name = copy.name;
    	class_index = copy.class_index;
    	is_selected = copy.is_selected;
    	class_distribution = copy.class_distribution;
    	accuracy_rotation = copy.accuracy_rotation;
    	accuracy_translation= copy.accuracy_translation;
    	estimated_resolution= copy.estimated_resolution;
    	particle_nr = copy.particle_nr;
    	class_score = copy.class_score;
    	relative_signal_intensity = copy.relative_signal_intensity;
    	edge_signal = copy.edge_signal;
    	scattered_signal = copy.scattered_signal;
    	weighted_resolution = copy.weighted_resolution;
    	relative_resolution = copy.relative_resolution;
    	lowpass_filtered_img_avg = copy.lowpass_filtered_img_avg;
    	lowpass_filtered_img_stddev = copy.lowpass_filtered_img_stddev;
    	lowpass_filtered_img_minval = copy.lowpass_filtered_img_minval;
    	lowpass_filtered_img_maxval = copy.lowpass_filtered_img_maxval;
    	circular_mask_moments = copy.circular_mask_moments;
    	ring_moments = copy.ring_moments;
    	inner_circle_moments = copy.inner_circle_moments;
    	protein_area = copy.protein_area;
    	solvent_area = copy.solvent_area;
    	protein_moments = copy.protein_moments;
    	solvent_moments = copy.solvent_moments;
    	lbp = copy.lbp;
    	lbp_p = copy.lbp_p;
    	lbp_s = copy.lbp_s;
    	haralick_p = copy.haralick_p;
    	haralick_s = copy.haralick_s;
    	zernike_moments = copy.zernike_moments;
	    granulo = copy.granulo;
	    total_entropy = copy.total_entropy;
	    protein_entropy = copy.protein_entropy;
	    solvent_entropy = copy.solvent_entropy;
	    subimages = copy.subimages;
	    CAR = copy.CAR;
		return *this;
	}

    std::vector<float> toNormalizedVector()
    {

		std::vector<float> out(24);

		// Order matters
		out[0] = (accuracy_rotation - feature_normalization_global_mean[0]) / feature_normalization_global_stddev[0];
		out[1] = (accuracy_translation - feature_normalization_global_mean[1]) / feature_normalization_global_stddev[1];
		out[2] = (weighted_resolution - feature_normalization_global_mean[2]) / feature_normalization_global_stddev[2];
		out[3] = (relative_resolution - feature_normalization_global_mean[3]) / feature_normalization_global_stddev[3];

		out[4] = (ring_moments.mean - feature_normalization_global_mean[4]) / feature_normalization_global_stddev[4];
		out[5] = (ring_moments.stddev - feature_normalization_global_mean[5]) / feature_normalization_global_stddev[5];
		out[6] = (protein_moments.stddev - feature_normalization_global_mean[6]) / feature_normalization_global_stddev[6];
		out[7] = (solvent_moments.mean - feature_normalization_global_mean[7]) / feature_normalization_global_stddev[7];
		out[8] = (solvent_moments.stddev - feature_normalization_global_mean[8]) / feature_normalization_global_stddev[8];
		out[9] = (scattered_signal - feature_normalization_global_mean[9]) / feature_normalization_global_stddev[9];
		out[10] = (edge_signal - feature_normalization_global_mean[10]) / feature_normalization_global_stddev[10];

		out[11] = (lowpass_filtered_img_avg - feature_normalization_global_mean[11]) / feature_normalization_global_stddev[11];
		out[12] = (lowpass_filtered_img_stddev - feature_normalization_global_mean[12]) / feature_normalization_global_stddev[12];
		out[13] = (lowpass_filtered_img_minval - feature_normalization_global_mean[13]) / feature_normalization_global_stddev[13];
		out[14] = (lowpass_filtered_img_maxval - feature_normalization_global_mean[14]) / feature_normalization_global_stddev[14];

		out[15] = (granulo[0] - feature_normalization_global_granulo_mean[0]) / feature_normalization_global_granulo_stddev[0];
		out[16] = (granulo[1] - feature_normalization_global_granulo_mean[1]) / feature_normalization_global_granulo_stddev[1];
		out[17] = (granulo[2] - feature_normalization_global_granulo_mean[2]) / feature_normalization_global_granulo_stddev[2];
		out[18] = (granulo[3] - feature_normalization_global_granulo_mean[3]) / feature_normalization_global_granulo_stddev[3];
		out[19] = (granulo[4] - feature_normalization_global_granulo_mean[4]) / feature_normalization_global_granulo_stddev[4];
		out[20] = (granulo[5] - feature_normalization_global_granulo_mean[5]) / feature_normalization_global_granulo_stddev[5];

		if (feature_normalization_local_ps_stddev < 1e-10)
		{
			out[21] = 0.;
		}
		else
		{
			out[21] = (protein_moments.sum - feature_normalization_local_ps_mean) / feature_normalization_local_ps_stddev;
		}

		if (feature_normalization_local_ss_stddev < 1e-10)
		{
			out[22] = 0.;
		}
		else
		{
			out[22] = (solvent_moments.sum - feature_normalization_local_ss_mean) / feature_normalization_local_ss_stddev;
		}

		if (feature_normalization_local_rsi_stddev < 1e-10)
		{
			out[23] = 0.;
		}
		else
		{
			out[23] = (relative_signal_intensity - feature_normalization_local_rsi_mean) / feature_normalization_local_rsi_stddev;
		}

		return out;

    }

};

class ZernikeMomentsExtractor
{
public:
	std::vector<RFLOAT> getZernikeMoments(MultidimArray<RFLOAT> img, long z_order, double radius, bool verb);

private:
	double factorial(long n);
	double zernikeR(int n, int l, double r);
	Complex zernikeZ(MultidimArray<RFLOAT> img, int n, int l, double r_max);
};

#define HARALICK_EPS 1e-6
class HaralickExtractor
{

private:

    MultidimArray<RFLOAT> matcooc; //GLCM
    MultidimArray<RFLOAT> margprobx;
    MultidimArray<RFLOAT> margproby;
    MultidimArray<RFLOAT> probsum; //sum probability
    MultidimArray<RFLOAT> probdiff; //diff probability
    RFLOAT hx, hy; //entropy of margprobx and y
    RFLOAT meanx, meany, stddevx, stddevy;
    bool initial=false; //marks if above variables are set

    RFLOAT Entropy(MultidimArray<RFLOAT> arr);
    void fast_init();
    std::vector<RFLOAT> cooc_feats();
    std::vector<RFLOAT> margprobs_feats();
    MultidimArray<RFLOAT> fast_feats(bool verbose=false);
    MultidimArray<RFLOAT> MatCooc(MultidimArray<int> img, int N, int deltax, int deltay, MultidimArray<int> *mask=NULL);

public:

    std::vector<RFLOAT> getHaralickFeatures(MultidimArray<RFLOAT> img, MultidimArray<int> *mask=NULL, bool verbose=false);

};



class ClassRanker
{

public:

	IOParser parser;
	FileName fn_out, fn_ext, fn_optimiser, fn_model, fn_data, fn_select, fn_job_score, fn_cf, fn_mask_dir;
	FileName fn_features, fn_sel_parts, fn_sel_classavgs, fn_root;

	RFLOAT minRes, job_score;
	RFLOAT radius_ratio, radius;
	RFLOAT particle_diameter, circular_mask_radius, uniform_angpix = 4.0;
	RFLOAT binary_threshold, lowpass;
    int debug, verb, start_class, end_class;
    bool do_relative_threshold;

	// Total number of particles in one jobs (always needed)
	long int total_nr_particles = 0;

    HaralickExtractor haralick_extractor;
	ZernikeMomentsExtractor zernike_extractor;

	// Also rank the classes in the input optimiser (otherwise only output feature file for network training purposes)
	bool do_ranking;
	// Perform selection of classes based on predicted scores
	bool do_select;
	// Also write out normalized feature vector
	bool do_write_normalized_features;

	// Make subimages out of class averages at standardized pixel size
	bool do_subimages, only_do_subimages;
    int nr_subimages, subimage_boxsize;

	// Automatically select classes
    RFLOAT select_min_score, select_max_score;
    int select_min_classes, select_min_parts;

    // Save some time by limiting calculations
	int only_use_this_class;
	bool do_skip_angular_errors, do_granularity_features, do_save_masks, do_save_mask_c;

	bool do_ctf_correction, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak;
	MlModel mymodel;
	Experiment mydata;
	MetaDataTable MD_optimiser, MD_select;
	std::vector<classFeatures> features_all_classes, preread_features_all_classes;

	FileName fn_pytorch_model;
	FileName fn_pytorch_script;
	FileName python_interpreter;

public:

	ClassRanker(){}

	/** Destructor
	  *
	  * Clears everything
	  *
	  * @code
	  * FourierInterpolator fourint;
	  * @endcode
	  */
	~ClassRanker() {}

	void read(int argc, char **argv, int rank = 0);

	// Print usage instructions
	void usage();

	// Initialise some stuff after reading
	void initialise();

	// Execute the program
	void run();

	/* Get path to the default pytorch model [LINUX ONLY]
	 * Check if file exists, return empty string otherwise
	 */
	static std::string get_default_pytorch_model_path();

	/* Get path to the python script for executing pytorch model [LINUX ONLY]
	 * Check if file exists, return empty string otherwise
	 */
	static std::string get_python_script_path();

private:

	int getClassIndex(FileName &name);

	MultidimArray<RFLOAT> getSubimages(MultidimArray<RFLOAT> &img, int boxsize = 16,
			int nr_images = 25, MultidimArray<int> *mask = NULL);

	moments calculateMoments(MultidimArray<RFLOAT> &img,
			RFLOAT inner_radius, RFLOAT outer_radius, MultidimArray<int> *mask = NULL);

	void calculatePvsLBP(MultidimArray<RFLOAT> I, MultidimArray<int> &p_mask, MultidimArray<int> &s_mask, classFeatures &cf);

	std::vector<RFLOAT> calculateGranulo(const MultidimArray<RFLOAT> &I);

	RFLOAT findResolution(classFeatures &cf);

	void calculateExpectedAngularErrors(int iclass, classFeatures &cf);

	RFLOAT getClassScoreFromJobScore(classFeatures &cf, RFLOAT minRes);

	void makeSolventMasks(classFeatures cf, MultidimArray<RFLOAT> img, MultidimArray<RFLOAT> &lpf, MultidimArray<int> &p_mask, MultidimArray<int> &s_mask,
				RFLOAT &scattered_signal, long &protein_area, long &solvent_area);

	void saveMasks(Image<RFLOAT> &img, MultidimArray<RFLOAT> &lpf, MultidimArray<int> &p_mask,
			MultidimArray<int> &s_mask, classFeatures &cf);

	void maskCircumference(MultidimArray<int> p_mask, RFLOAT &protein_C, classFeatures cf, bool do_save_mask_c);

	void correctCtfUntilFirstPeak(MultidimArray<RFLOAT> &in, CTF ctf);

	void localNormalisation(std::vector<classFeatures> &features_all_classes);

	void onlyGetSubimages();

	void getFeatures();

	void readFeatures();

	void writeFeatures();

	void deployTorchModel(FileName &model_path, std::vector<float> &features, std::vector<float> &subimages, std::vector<float> &score);

	void performRanking();
};
#endif /* CLASS_RANKER_H_ */
