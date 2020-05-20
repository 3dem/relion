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

#ifndef _HEALPIX_SAMPLING_HH
#define _HEALPIX_SAMPLING_HH

#include "src/Healpix_2.15a/healpix_base.h"
#include "src/metadata_table.h"
#include "src/macros.h"
#include "src/multidim_array.h"
#include "src/symmetries.h"
#include "src/euler.h"
#include "src/transformations.h"
#include "src/helix.h"

// For the angular searches
#define NOPRIOR 0
#define PRIOR_ROTTILT_PSI 1

class HealpixSampling
{

public:
	/** Healpix sampling object */
    Healpix_Base healpix_base;

    /** Random perturbation */
    RFLOAT random_perturbation;

    /** Amount of random perturbation */
    RFLOAT perturbation_factor;

    /** In-plane (psi-angle) sampling rate (and original one)
     */
    RFLOAT psi_step, psi_step_ori;

    /** Healpix order (and original one) */
    int healpix_order, healpix_order_ori;

    /* Translational search range and sampling rate (as of v3.1 in Angstroms!) (and original one)
     */
    // Jun19,2015 - Shaoda, Helical refinement (all in Angstroms!)
    RFLOAT offset_range, offset_step, helical_offset_step, offset_range_ori, offset_step_ori;

    /** Flag whether this is a real 3D sampling */
    bool is_3D;

    /** Flag whether the translations are 3D (for volume refinement) */
    bool is_3d_trans;

    /** Flag whether relax symmetry */
    bool isRelax;

    /** Name of the Symmetry group */
    FileName fn_sym;
    FileName fn_sym_relax;

    /** List of symmetry operators */
    std::vector <Matrix2D<RFLOAT> > R_repository, L_repository;
    std::vector <Matrix2D<RFLOAT> > R_repository_relax, L_repository_relax;

    /** Two numbers that describe the symmetry group */
    int pgGroup;
    int pgOrder;
    int pgGroupRelaxSym;
    int pgOrderRelaxSym;

    /** Limited tilt angle range */
    RFLOAT limit_tilt;

    /** vector with the original pixel number in the healpix object */
    std::vector<int> directions_ipix;

    /** vector with sampling points described by angles */
    std::vector<RFLOAT > rot_angles, tilt_angles;

    /** vector with the psi-samples */
    std::vector<RFLOAT> psi_angles;

    /** vector with the X,Y(,Z)-translations (as of v3.1 in Angstroms!) */
    std::vector<RFLOAT> translations_x, translations_y, translations_z;


public:

    // Empty constructor
    HealpixSampling():
    	offset_range(0),
		offset_step(0),
		is_3d_trans(false),
		pgGroup(0),
		pgGroupRelaxSym(0),
		perturbation_factor(0),
		is_3D(false),
		random_perturbation(0),
		psi_step(0),
		limit_tilt(0),
		healpix_order(0),
		pgOrder(0),
		pgOrderRelaxSym(0)
    {}

    // Destructor
    ~HealpixSampling()
    {
    	directions_ipix.clear();
    	rot_angles.clear();
    	tilt_angles.clear();
    	psi_angles.clear();
    	translations_x.clear();
    	translations_y.clear();
    	translations_z.clear();
    }

    // Start from all empty vectors and meaningless parameters
    void clear();

    /** Set up the entire sampling object
     *
     * The in-plane (psi-angle) sampling is linear,
     * input_psi_sampling is modified to contain an integer number of equally-sized sampling points
     * For the 3D-case, a negative input_psi_sampling will result in a psi-sampling similar to the sqrt of a HealPix pixel area.
	 *
     * The HEALPix sampling is implemented as described by Gorski et al (2005), The Astrophysical Journal, 622:759-771
     * The order defines the number of sampling points, and thereby the angular sampling rate
     * From this paper is the following table:
     *
     * order	Npix	Theta-sampling
     * 0		12		58.6
     * 1		48		29.3
     * 2		192		14.7
     * 3		768		7.33
     * 4		3072	3.66
     * 5		12288	1.83
     * 6		49152	0.55
     * 7		196608	0.28
     * 8		786432	0.14
     * etc...
     *
     * */
    // May 6, 2015 - Shaoda & Sjors - initialise for helical translations
    void initialise(
			int ref_dim = -1,
			bool do_3d_trans = false,
			bool do_changepsi = false,
			bool do_warnpsi = false,
			bool do_local_searches_helical = false,
			bool do_helical_refine = false,
			RFLOAT rise_Angst  = 0.,
			RFLOAT twist_deg = 0.);

    // Initialize the symmetry matrices
    void initialiseSymMats(FileName fn_sym_, int & pgGroup_,
    		int & pgOrder_, std::vector <Matrix2D<RFLOAT> > & R_repository,
			std::vector <Matrix2D<RFLOAT> > & L_repository);

    // Reset the random perturbation
    void resetRandomlyPerturbedSampling();

    // Read in all information from the command line to build the sampling object
    // Jun19,2015 - Shaoda, is that available???
    //void read(IOParser &parser, int ori_size, int ref_dim);

    // Read CL options after a -continue statement.
    // Jun19,2015 - Shaoda, is that available???
    //void readContinue(int argc, char **argv, bool &directions_have_changed);

    // Read in all information from a STAR file (for restarting)
    void read(FileName fn_in);

    // Write the sampling information to a STAR file
    void write(FileName fn_out);

    /* Set the non-oversampled list of translations in Angstroms
     * For single particles, offset ranges are equal along different directions
     * For helices, x offsets (along helical axis) should be less within -+0.5 * rise
     * */
    void setTranslations(
			RFLOAT new_offset_step = -1.,
    		RFLOAT new_offset_range = -1.,
    		bool do_local_searches_helical = false,
    		bool do_helical_refine = false,
			RFLOAT new_helical_offset_step = -1.,
    		RFLOAT helical_rise_Angst = 0.,
    		RFLOAT helical_twist_deg = 0.);

    /* Add a single translation */
    void addOneTranslation(
    		RFLOAT offset_x,
    		RFLOAT offset_y,
    		RFLOAT offset_z,
    		bool do_clear = false,
    		bool do_helical_refine = false,
    		RFLOAT rot_deg = 0.,
    		RFLOAT tilt_deg = 0.,
    		RFLOAT psi_deg = 0.);

    /* Set the non-oversampled lists of directions and in-plane rotations */
    void setOrientations(int _order = -1, RFLOAT _psi_step = -1.);

    /* Add a single orientation */
    void addOneOrientation(RFLOAT rot, RFLOAT tilt, RFLOAT psi, bool do_clear = false);

    /* Write all orientations as a sphere in a bild file
     * Mainly useful for debugging */
    void writeAllOrientationsToBild(FileName fn_bild, std::string rgb = "1 0 0", RFLOAT size = 0.025);
    void writeNonZeroPriorOrientationsToBild(FileName fn_bild, RFLOAT rot_prior, RFLOAT tilt_prior,
    		std::vector<int> &pointer_dir_nonzeroprior, std::string rgb = "0 0 1", RFLOAT size = 0.025);

    /* Sjors, 9nov2015: new rot-priors for DNA-origami-bound refinements
     */
    RFLOAT calculateDeltaRot(Matrix1D<RFLOAT> my_direction, RFLOAT rot_prior);

    /* Select all orientations with zero prior probabilities
     * store all these in the vectors pointer_dir_nonzeroprior and pointer_psi_nonzeroprior
     * Also precalculate their prior probabilities and store in directions_prior and psi_prior
     */
    // Jun 04 - Shaoda & Sjors, Bimodel psi searches for helices
    void selectOrientationsWithNonZeroPriorProbability(
    		RFLOAT prior_rot, RFLOAT prior_tilt, RFLOAT prior_psi,
    		RFLOAT sigma_rot, RFLOAT sigma_tilt, RFLOAT sigma_psi,
    		std::vector<int> &pointer_dir_nonzeroprior, std::vector<RFLOAT> &directions_prior,
    		std::vector<int> &pointer_psi_nonzeroprior, std::vector<RFLOAT> &psi_prior,
			bool do_bimodal_search_psi = false,
    		RFLOAT sigma_cutoff = 3., RFLOAT sigma_tilt_from_ninety = -1., RFLOAT sigma_psi_from_zero = -1.);

    void selectOrientationsWithNonZeroPriorProbabilityFor3DHelicalReconstruction(
    		RFLOAT prior_rot, RFLOAT prior_tilt, RFLOAT prior_psi,
    		RFLOAT sigma_rot, RFLOAT sigma_tilt, RFLOAT sigma_psi,
    		std::vector<int> &pointer_dir_nonzeroprior, std::vector<RFLOAT> &directions_prior,
    		std::vector<int> &pointer_psi_nonzeroprior, std::vector<RFLOAT> &psi_prior,
    		bool do_auto_refine_local_searches,
    		RFLOAT prior_psi_flip_ratio = 0.5,
			RFLOAT prior_rot_flip_ratio = 0.5,  // KThurber
    		RFLOAT sigma_cutoff = 3.);

    // Find the symmetry mate by searching the Healpix library
    void findSymmetryMate(long int idir_, RFLOAT prior_,
    		std::vector<int> &pointer_dir_nonzeroprior,
			std::vector<RFLOAT> &directions_prior,
			std::vector<bool> &idir_flag);

    /** Get the symmetry group of this sampling object
     */
    FileName symmetryGroup();

    /* Get the original HEALPix index for this direction
     * Note that because of symmetry-equivalence removal idir no longer corresponds to the HEALPix pixel number
     *
     */
    long int getHealPixIndex(long int idir);

    /** The geometrical considerations about the symmetry below require that rot = [-180,180] and tilt [0,180]
     */
    void checkDirection(RFLOAT &rot, RFLOAT &tilt);

    /* Get the rot and tilt angles in the center of the ipix'th HEALPix sampling pixel
     * This involves calculations in the HEALPix library
     */
    void getDirectionFromHealPix(long int ipix, RFLOAT &rot, RFLOAT &tilt);

    /* Get the translational sampling step in Angstroms */
    RFLOAT getTranslationalSampling(int adaptive_oversampling = 0);

    /* Get the translational sampling step along helical axis in Angstroms */
    RFLOAT getHelicalTranslationalSampling(int adaptive_oversampling = 0);

    /* Get approximate angular sampling in degrees for any adaptive oversampling
     */
    RFLOAT getAngularSampling(int adaptive_oversampling = 0);

    /* Get the number of symmetry-unique sampling points
     * Note that because of symmetry-equivalence removal this number is not the number of original HEALPix pixels
     * In the case of orientational priors, the number of directions with non-zero prior probability is returned
     */
    long int NrDirections(int oversampling_order = 0, const std::vector<int> *pointer_dir_nonzeroprior = NULL);

    /* Get the number of in-plane (psi-angle) sampling points
     */
    long int NrPsiSamplings(int oversampling_order = 0, const std::vector<int> *pointer_psi_nonzeroprior = NULL);

    /* Get the number of in-plane translational sampling points
     */
    long int NrTranslationalSamplings(int oversampling_order = 0);

    /* Get the total number of (oversampled) sampling points, i.e. all (rot, tilt, psi, xoff, yoff) quintets
    */
    long int NrSamplingPoints(int oversampling_order = 0,
    		const std::vector<int> *pointer_dir_nonzeroprior = NULL,
    		const std::vector<int> *pointer_psi_nonzeroprior = NULL);

    /* How often is each orientation oversampled? */
    int oversamplingFactorOrientations(int oversampling_order);

    /* How often is each translation oversampled? */
    int oversamplingFactorTranslations(int oversampling_order);

    /* Get the rot and tilt angles from the precalculated sampling_points_angles vector
     * This does not involve calculations in the HEALPix library
     * Note that because of symmetry-equivalence removal idir no longer corresponds to the HEALPix pixel number
     */
    void getDirection(long int idir, RFLOAT &rot, RFLOAT &tilt);

    /* Get the value for the ipsi'th precalculated psi angle
     */
    void getPsiAngle(long int ipsi, RFLOAT &psi);

    /* Get the value for the itrans'th precalculated translations
     */
    void getTranslationInPixel(long int itrans, RFLOAT my_pixel_size, RFLOAT &trans_x, RFLOAT &trans_y, RFLOAT &trans_z);

    /* Get the position of this sampling point in the original array */
    long int getPositionSamplingPoint(int iclass, long int idir, long int ipsi, long int itrans);

    /* Get the position of this sampling point in the oversampled array */
    long int getPositionOversampledSamplingPoint(long int ipos, int oversampling_order, int iover_rot, int iover_trans);

    /* Get the vectors of (xx, yy) for a more finely (oversampled) translational sampling
     * The oversampling_order is the difference in order of the original (coarse) and the oversampled (fine) sampling
     * An oversampling_order == 0  will give rise to the same (xx, yy) pair as the original itrans.
     * An oversampling_order == 1 will give rise to 2*2 new (rot, tilt) pairs.
     * An oversampling_order == 2 will give rise to 4*4 new (rot, tilt) pairs.
     * etc.
     */
    void getTranslationsInPixel(long int itrans, int oversampling_order, RFLOAT my_pixel_size,
    		std::vector<RFLOAT > &my_translations_x,
    		std::vector<RFLOAT > &my_translations_y,
    		std::vector<RFLOAT > &my_translations_z,
			bool do_helical_refine = false);

    /* Get the vectors of (rot, tilt, psi) angle triplets for a more finely (oversampled) sampling
     * The oversampling_order is the difference in order of the original (coarse) and the oversampled (fine) sampling
     * An oversampling_order == 0  will give rise to the same (rot, tilt, psi) triplet as the original ipix.
     * An oversampling_order == 1 will give rise to 2*2*2 new (rot, tilt, psi) triplets.
     * An oversampling_order == 2 will give rise to 4*4*4 new (rot, tilt, psi) triplets.
     * etc.
     *
     * If only_nonzero_prior is true, then only the orientations with non-zero prior probabilities will be returned
     * This is for local angular searches
     */
    void getOrientations(long int idir, long int ipsi, int oversampling_order,
    		std::vector<RFLOAT > &my_rot, std::vector<RFLOAT > &my_tilt, std::vector<RFLOAT > &my_psi,
    		std::vector<int> &pointer_dir_nonzeroprior, std::vector<RFLOAT> &directions_prior,
    		std::vector<int> &pointer_psi_nonzeroprior, std::vector<RFLOAT> &psi_prior);

    /* Gets the vector of psi angles for a more finely (oversampled) sampling and
     * pushes each instance back into the oversampled_orientations vector with the given rot and tilt
     * The oversampling_order is the difference in order of the original (coarse) and the oversampled (fine) sampling
     * An oversampling_order == 0  will give rise to the same psi angle as the original ipsi.
     * An oversampling_order == 1 will give rise to 2 new psi angles
     * An oversampling_order == 2 will give rise to 4 new psi angles
     * etc.
     */
    void pushbackOversampledPsiAngles(long int ipsi, int oversampling_order,
    		RFLOAT rot, RFLOAT tilt, std::vector<RFLOAT> &oversampled_rot,
    		std::vector<RFLOAT> &oversampled_tilt, std::vector<RFLOAT> &oversampled_psi);

    /* Calculate an angular distance between two sets of Euler angles */
    RFLOAT calculateAngularDistance(RFLOAT rot1, RFLOAT tilt1, RFLOAT psi1,
    		RFLOAT rot2, RFLOAT tilt2, RFLOAT psi2);

    /* Write a BILD file describing the angular distribution
     *  R determines the radius of the sphere on which cylinders will be placed
     *  Rmax_frac determines the length of the longest cylinder (relative to R, 0.2 + +20%)
     *  width_frac determines how broad each cylinder is. frac=1 means they touch each other
     * */
    void writeBildFileOrientationalDistribution(MultidimArray<RFLOAT> &pdf_direction,
    		FileName &fn_bild, RFLOAT R, RFLOAT offset = 0.,
			const Matrix2D<RFLOAT> *Aorient = NULL, const Matrix1D<RFLOAT> *Acom = NULL,
			RFLOAT Rmax_frac = 0.3, RFLOAT width_frac = 0.5);

private:

    /* Eliminate points from the sampling_points_vector and sampling_points_angles vectors
     * that are outside the allowed tilt range.
     * Let tilt angles range from -90 to 90, then:
     * if (limit_tilt > 0) then top views (that is with ABS(tilt) > limit_tilt) are removed and side views are kept
     * if (limit_tilt < 0) then side views (that is with ABS(tilt) < limit_tilt) are removed and top views are kept
    */
    void removePointsOutsideLimitedTiltAngles();

    /* Eliminate symmetry-equivalent points from the sampling_points_vector and sampling_points_angles vectors
        This function first calls removeSymmetryEquivalentPointsGeometric,
        and then checks each point versus all others to calculate an angular distance
        If this distance is less than 0.8 times the angular sampling, the point is deleted
        This cares care of sampling points near the edge of the geometrical considerations
    */
    void removeSymmetryEquivalentPoints(RFLOAT max_ang);

    /* eliminate symmetry-related points based on simple geometrical considerations,
        symmetry group, symmetry order */
    void removeSymmetryEquivalentPointsGeometric(const int symmetry, int sym_order,
												 std::vector <Matrix1D<RFLOAT> >  &sampling_points_vector);



};
//@}
#endif
