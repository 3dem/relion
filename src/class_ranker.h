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

// This contains 4 moments for an image
class moments
{

public:

	RFLOAT mean, stddev, skew, kurt;
    moments(): mean(0), stddev(0), skew(0), kurt(0){}

};

// This defines all features for a single class
class classFeatures
{

public:

	// Class-wise features
    FileName name;
    int is_selected;
    RFLOAT class_distribution, accuracy_rotation, accuracy_translation, estimated_resolution, particle_nr;
    RFLOAT class_score, edge_signal, scattered_signal, weighted_resolution, relative_resolution;
    RFLOAT lowpass_filtered_img_avg, lowpass_filtered_img_stddev, lowpass_filtered_img_minval, lowpass_filtered_img_maxval;
    moments circular_mask_moments, ring_moments, inner_circle_moments, protein_moments, solvent_moments;

    // Granularity features only for historical development reasons
    std::vector<RFLOAT> lbp, lbp_p, lbp_s, haralick_p, haralick_s, zernike_moments, granulo;
    double total_entropy, protein_entropy, solvent_entropy;

    classFeatures(): name(""),
			is_selected(0),
			particle_nr(0),
    		class_distribution(0),
			accuracy_rotation(0),
			accuracy_translation(0),
			estimated_resolution(999.0),
			weighted_resolution(999.0),
			relative_resolution(999.0),
			class_score(-1),
			edge_signal(-1.),
			scattered_signal(-1.),
			lowpass_filtered_img_avg(0),
			lowpass_filtered_img_stddev(0),
			lowpass_filtered_img_minval(0),
			lowpass_filtered_img_maxval(0),
			total_entropy(0.),
			protein_entropy(0.),
			solvent_entropy(0.)
    {
    }

    // Destructor needed for work with vectors
	~classFeatures() {}

	// Copy constructor needed for work with vectors
    classFeatures(classFeatures const& copy)
	{
    	name = copy.name;
    	is_selected = copy.is_selected;
    	class_distribution = copy.class_distribution;
    	accuracy_rotation = copy.accuracy_rotation;
    	accuracy_translation= copy.accuracy_translation;
    	estimated_resolution= copy.estimated_resolution;
    	particle_nr = copy.particle_nr;
    	class_score = copy.class_score;
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

	}

	// Define assignment operator in terms of the copy constructor
    classFeatures& operator=(classFeatures const& copy)
	{
    	name = copy.name;
    	is_selected = copy.is_selected;
    	class_distribution = copy.class_distribution;
    	accuracy_rotation = copy.accuracy_rotation;
    	accuracy_translation= copy.accuracy_translation;
    	estimated_resolution= copy.estimated_resolution;
    	particle_nr = copy.particle_nr;
    	class_score = copy.class_score;
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

		return *this;
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
	FileName fn_out, fn_ext, fn_optimiser, fn_model, fn_select, fn_job_score, fn_cf;
	FileName fn_features, fn_sel_parts, fn_sel_classavgs;

	RFLOAT minRes, job_score;
	RFLOAT radius_ratio, radius;
	RFLOAT circular_mask_radius, uniform_angpix = 4.0;
	RFLOAT binary_threshold, lowpass;
    int debug, verb, start_class, end_class;

    HaralickExtractor haralick_extractor;
	ZernikeMomentsExtractor zernike_extractor;

	// Also rank the classes in the input optimiser (otherwise only output feature file for network training purposes)
	bool do_ranking;
	// Perform selection of classes based on predicted scores
	bool do_select;
	RFLOAT select_min_score, select_max_score;

    // Save some time by limiting calculations
	int only_use_this_class;
	bool do_skip_angular_errors, do_granularity_features;

	MlOptimiser myopt;
	MetaDataTable MD_optimiser, MD_select;
	std::vector<classFeatures> features_all_classes, preread_features_all_classes;

	bool do_save_masks, save_masks_only;

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

	long int getClassIndex(FileName &name);

	moments calculateMoments(MultidimArray<RFLOAT> &img,
			RFLOAT inner_radius, RFLOAT outer_radius, MultidimArray<int> *mask = NULL);

	void calculatePvsLBP(MultidimArray<RFLOAT> I, MultidimArray<int> &p_mask, MultidimArray<int> &s_mask, classFeatures &cf);

	std::vector<RFLOAT> calculateGranulo(const MultidimArray<double> &I);

	RFLOAT findResolution(classFeatures &cf);

	void calculateExpectedAngularErrors(int iclass, classFeatures &cf);

	RFLOAT getClassScoreFromJobScore(classFeatures &cf, RFLOAT minRes);

	void makeSolventMasks(MultidimArray<RFLOAT> img, MultidimArray<RFLOAT> &lpf, MultidimArray<int> &p_mask, MultidimArray<int> &s_mask,
				RFLOAT &scattered_signal, long &protein_area, long &solvent_area);

	void saveMasks(MultidimArray<RFLOAT> &lpf, MultidimArray<int> &p_mask,
			MultidimArray<int> &s_mask, classFeatures &cf);

	void correctCtfUntilFirstPeak(MultidimArray<RFLOAT> &in, CTF ctf);

	void getFeatures();

	void readFeatures();

	void writeFeatures();

	void performRanking();

};
#endif /* CLASS_RANKER_H_ */
