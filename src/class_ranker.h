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
#include "src/ml_optimiser.h"

#ifdef _TORCH_ENABLED
#include <torch/script.h> // One-stop header.
#endif //_TORCH_ENABLED

// This contains 4 moments for an image
class moments
{

public:

	RFLOAT sum, mean, stddev, skew, kurt;
    moments(): sum(0), mean(0), stddev(0), skew(0), kurt(0){}

};

class NormalizedFeatures
{
public:
	RFLOAT accuracy_rotation, accuracy_translation, weighted_resolution, relative_resolution, ring_mean, ring_stddev, protein_stddev,
			solvent_mean, solvent_stddev, scattered_signal, edge_signal, relative_signal_intensity, lowpass_filtered_img_avg, lowpass_filtered_img_stddev,
			lowpass_filtered_img_minval, lowpass_filtered_img_maxval, granulo0, granulo1, granulo2, granulo3, granulo4, granulo5,
			protein_sum, solvent_sum;
	NormalizedFeatures(): accuracy_rotation(999.),
			accuracy_translation(999.),
			weighted_resolution(999.),
			relative_resolution(999.),
			ring_mean(0.),
			ring_stddev(0.),
			protein_stddev(0.),
			solvent_mean(0.),
			solvent_stddev(0.),
			scattered_signal(0.),
			edge_signal(0.),
			relative_signal_intensity(0.),
			lowpass_filtered_img_avg(0.),
			lowpass_filtered_img_stddev(0.),
			lowpass_filtered_img_minval(0.),
			lowpass_filtered_img_maxval(0.),
			protein_sum(0.),
			solvent_sum(0.),
			granulo0(0.),
			granulo1(0.),
			granulo2(0.),
			granulo3(0.),
			granulo4(0.),
			granulo5(0.)
			{}

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
    double total_entropy, protein_entropy, solvent_entropy;

    NormalizedFeatures normalized_features;

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
	    normalized_features = copy.normalized_features;

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
	    normalized_features = copy.normalized_features;

		return *this;
	}

	std::vector<float> toVector() {
		std::vector<float> out(24);

		// Order matters
		out[0] =  normalized_features.accuracy_rotation;
		out[1] =  normalized_features.accuracy_translation;
		out[2] =  normalized_features.weighted_resolution;
		out[3] =  normalized_features.relative_resolution;

		out[4] =  normalized_features.ring_mean;
		out[5] =  normalized_features.ring_stddev;
		out[6] =  normalized_features.protein_stddev;
		out[7] =  normalized_features.solvent_mean;
		out[8] =  normalized_features.solvent_stddev;
		out[9] = normalized_features.scattered_signal;
		out[10] = normalized_features.edge_signal;

		out[11] = normalized_features.lowpass_filtered_img_avg;
		out[12] = normalized_features.lowpass_filtered_img_stddev;
		out[13] = normalized_features.lowpass_filtered_img_minval;
		out[14] = normalized_features.lowpass_filtered_img_maxval;

		out[15] = normalized_features.granulo0;
		out[16] = normalized_features.granulo1;
		out[17] = normalized_features.granulo2;

		out[18] = normalized_features.granulo3;
		out[19] = normalized_features.granulo4;
		out[20] = normalized_features.granulo5;

		out[21] = normalized_features.protein_sum;
		out[22] = normalized_features.solvent_sum;
		out[23] = normalized_features.relative_signal_intensity;

		return out;
    }

};

class ZernikeMomentsExtractor
{
public:
	std::vector<double> getZernikeMoments(MultidimArray<double> img, long z_order, double radius, bool verb);

private:
	double factorial(long n);
	double zernikeR(int n, int l, double r);
	Complex zernikeZ(MultidimArray<double> img, int n, int l, double r_max);
};

#define HARALICK_EPS 1e-6
class HaralickExtractor
{

private:

	MultidimArray<double> matcooc; //GLCM
    MultidimArray<double> margprobx;
    MultidimArray<double> margproby;
    MultidimArray<double> probsum; //sum probability
    MultidimArray<double> probdiff; //diff probability
    double hx, hy; //entropy of margprobx and y
    double meanx, meany, stddevx, stddevy;
    bool initial=false; //marks if above variables are set

    double Entropy(MultidimArray<double> arr);
    void fast_init();
    std::vector<double> cooc_feats();
    std::vector<double> margprobs_feats();
    MultidimArray<double> fast_feats(bool verbose=false);
    MultidimArray<RFLOAT> MatCooc(MultidimArray<int> img, int N, int deltax, int deltay, MultidimArray<int> *mask=NULL);

public:

    std::vector<double> getHaralickFeatures(MultidimArray<RFLOAT> img, MultidimArray<int> *mask=NULL, bool verbose=false);

};



class ClassRanker
{

public:

	IOParser parser;
	FileName fn_out, fn_ext, fn_optimiser, fn_model, fn_data, fn_select, fn_job_score, fn_cf, fn_mask_dir, fn_subimages, fn_subimage_star;
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

	// Make subimages out of class averages at standardized pixel size
	bool do_subimages;
    int nr_subimages, subimage_boxsize;

	RFLOAT select_min_score, select_max_score;

    // Save some time by limiting calculations
	int only_use_this_class;
	bool do_skip_angular_errors, do_granularity_features, do_save_masks, do_save_mask_c;

	bool do_ctf_correction, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak;
	MlModel mymodel;
	Experiment mydata;
	MetaDataTable MD_optimiser, MD_select;
	std::vector<classFeatures> features_all_classes, preread_features_all_classes;
	std::vector<RFLOAT> global_mean{
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
						global_stddev{
								22.650,
								9.113,
								0.000,
								0.008,
								0.080,
								0.396,
								0.342,
								0.143,
								0.235,
								0.006,
								0.006,
								0.033,
								0.184,
								1.181,
								3.090
							},
						global_granulo_mean{
								837.576,
								726.097,
								602.175,
								483.469,
								328.014,
								196.059
							},
						global_granulo_stddev{
								3836985.599,
								8095904.406,
								13002214.522,
								16015778.201,
								19356282.669,
								22181751.608
							};
//	std::vector<std::string> features_to_global_normalize{},
//							features_to_local_normalize{};
//
	FileName fn_torch_model;

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


private:

	int getClassIndex(FileName &name);

	MultidimArray<RFLOAT> getSubimages(MultidimArray<RFLOAT> &img, int boxsize = 16,
			int nr_images = 25, MultidimArray<int> *mask = NULL);

	moments calculateMoments(MultidimArray<RFLOAT> &img,
			RFLOAT inner_radius, RFLOAT outer_radius, MultidimArray<int> *mask = NULL);

	void calculatePvsLBP(MultidimArray<RFLOAT> I, MultidimArray<int> &p_mask, MultidimArray<int> &s_mask, classFeatures &cf);

	std::vector<RFLOAT> calculateGranulo(const MultidimArray<double> &I);

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

	void getFeatures();

	void readFeatures();

	void writeFeatures();

	float deployTorchModel(FileName &model_path, std::vector<float> &features);

	void performRanking();
};
#endif /* CLASS_RANKER_H_ */
