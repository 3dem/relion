/***************************************************************************
 *
 * Author: "Shaoda He"
 * MRC Laboratory of Molecular Biology
 *
 * Kent Thurber from the NIH provided code for rot-angle priors (indicated with // KThurber comments)
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

#ifndef HELIX_H_
#define HELIX_H_

#include "src/multidim_array.h"
#include "src/macros.h"
#include "src/complex.h"
#include "src/fftw.h"
#include "src/image.h"
#include "src/transformations.h"
#include "src/euler.h"
#include "src/assembly.h"
#include <src/jaz/single_particle/obs_model.h>
#include "src/time.h"
#include <set>

#define CART_TO_HELICAL_COORDS true
#define HELICAL_TO_CART_COORDS false
#define RELION_STAR_FORMAT 0
#define XIMDISP_COORDS_FORMAT 1
#define EMAN2_FORMAT 2
#define UNIMODAL_PSI_PRIOR_FLIP_RATIO (0.)
#define BIMODAL_PSI_PRIOR_FLIP_RATIO (0.5)

class HelicalSymmetryItem
{
public:
	RFLOAT twist_deg;
	RFLOAT rise_pix;
	RFLOAT dev;

	HelicalSymmetryItem()
	{
		twist_deg = rise_pix = -1.;
		dev = (1e35);
	}
	HelicalSymmetryItem(RFLOAT _twist_deg, RFLOAT _rise_pix, RFLOAT _dev = (1e35))
	{
		twist_deg = _twist_deg;
		rise_pix = _rise_pix;
		dev = _dev;
	}
};

void makeHelicalSymmetryList(
		std::vector<HelicalSymmetryItem>& list,
		RFLOAT rise_min_pix,
		RFLOAT rise_max_pix,
		RFLOAT rise_step_pix,
		bool search_rise,
		RFLOAT twist_min_deg,
		RFLOAT twist_max_deg,
		RFLOAT twist_step_deg,
		bool search_twist);

bool calcCCofHelicalSymmetry(
		const MultidimArray<RFLOAT>& v,
		RFLOAT r_min_pix,
		RFLOAT r_max_pix,
		RFLOAT z_percentage,
		RFLOAT rise_pix,
		RFLOAT twist_deg,
		RFLOAT& cc,
		int& nr_asym_voxels);

bool localSearchHelicalSymmetry(
		const MultidimArray<RFLOAT>& v,
		RFLOAT pixel_size_A,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_inner_radius_A,
		RFLOAT cyl_outer_radius_A,
		RFLOAT z_percentage,
		RFLOAT rise_min_A,
		RFLOAT rise_max_A,
		RFLOAT rise_inistep_A,
		RFLOAT& rise_refined_A,
		RFLOAT twist_min_deg,
		RFLOAT twist_max_deg,
		RFLOAT twist_inistep_deg,
		RFLOAT& twist_refined_deg,
		std::ostream* o_ptr = NULL);

RFLOAT getHelicalSigma2Rot(
		RFLOAT helical_rise_Angst,
		RFLOAT helical_twist_deg,
		RFLOAT helical_offset_step_Angst,
		RFLOAT rot_step_deg,
		RFLOAT old_sigma2_rot);

bool checkParametersFor3DHelicalReconstruction(
		bool ignore_symmetry,
		bool do_symmetry_local_refinement,
		int nr_asu,
		RFLOAT rise_initial_A,
		RFLOAT rise_min_A,
		RFLOAT rise_max_A,
		RFLOAT twist_initial_deg,
		RFLOAT twist_min_deg,
		RFLOAT twist_max_deg,
		int box_len,
		RFLOAT pixel_size_A,
		RFLOAT z_percentage,
		RFLOAT particle_diameter_A,
		RFLOAT tube_inner_diameter_A,
		RFLOAT tube_outer_diameter_A,
		bool verboseOutput = false);

void imposeHelicalSymmetryInRealSpace(
		MultidimArray<RFLOAT>& v,
		RFLOAT pixel_size_A,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_inner_radius_A,
		RFLOAT cyl_outer_radius_A,
		RFLOAT z_percentage,
		RFLOAT rise_A,
		RFLOAT twist_deg,
		RFLOAT cosine_width_pix);

// Some functions only for specific testing
void calcRadialAverage(
		const MultidimArray<RFLOAT>& v,
		std::vector<RFLOAT>& radial_avg_val_list);

void cutZCentralPartOfSoftMask(
		MultidimArray<RFLOAT>& mask,
		RFLOAT z_percentage,
		RFLOAT cosine_width = 5.);

void createCylindricalReference(
		MultidimArray<RFLOAT>& v,
		int box_size,
		RFLOAT inner_diameter_pix,
		RFLOAT outer_diameter_pix,
		RFLOAT cosine_width = 5.);

void createCylindricalReferenceWithPolarity(
		MultidimArray<RFLOAT>& v,
		int box_size,
		RFLOAT inner_diameter_pix,
		RFLOAT outer_diameter_pix,
		RFLOAT ratio_topbottom = 0.5,
		RFLOAT cosine_width = 5.);

void transformCartesianAndHelicalCoords(
		Matrix1D<RFLOAT>& in,
		Matrix1D<RFLOAT>& out,
		RFLOAT rot_deg,
		RFLOAT tilt_deg,
		RFLOAT psi_deg,
		bool direction);

void transformCartesianAndHelicalCoords(
		RFLOAT xin,
		RFLOAT yin,
		RFLOAT zin,
		RFLOAT& xout,
		RFLOAT& yout,
		RFLOAT& zout,
		RFLOAT rot_deg,
		RFLOAT tilt_deg,
		RFLOAT psi_deg,
		int dim,
		bool direction);

// C1 Helix
// If radius_A < 0, detect the radius of original assembly and calculate r as sqrt(x^2 + y^2)
void makeSimpleHelixFromPDBParticle(
		const Assembly& ori,
		Assembly& helix,
		RFLOAT radius_A,
		RFLOAT twist_deg,
		RFLOAT rise_A,
		int nr_copy,
		bool do_center = false);

void applySoftSphericalMask(
		MultidimArray<RFLOAT>& v,
		RFLOAT sphere_diameter = -1.,
		RFLOAT cosine_width = 5.);

void extractHelicalSegmentsFromTubes_Multiple(
		FileName& suffix_in,
		FileName& suffix_out,
		int format_tag,
		int nr_asu,
		RFLOAT rise_A,
		RFLOAT pixel_size_A,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix,
		bool bimodal_angular_priors = true,
		bool cut_into_segments = true);

void convertHelicalTubeCoordsToMetaDataTable(
		FileName& fn_in,
		MetaDataTable& MD_out,
		int& total_segments,
		int& total_tubes,
		int nr_asu,
		RFLOAT rise_A,
		RFLOAT pixel_size_A,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix,
		bool bimodal_angular_priors = true,
		bool cut_into_segments = true);

void combineParticlePriorsWithKaiLocalCTF(
		FileName& fn_priors,
		FileName& fn_local_ctf,
		FileName& fn_combined);

// Files of priors: mic1_priors.star, files of local CTF: mic1_local.star
// Then suffix_priors = _priors.star, suffix_local_ctf = _local.star
void combineParticlePriorsWithKaiLocalCTF_Multiple(
		std::string& suffix_priors,
		std::string& suffix_local_ctf,
		std::string& suffix_combined);

void setNullTiltPriorsInDataStar(
		FileName& fn_in,
		FileName& fn_out);

void removeBadTiltHelicalSegmentsFromDataStar(
		FileName& fn_in,
		FileName& fn_out,
		RFLOAT max_dev_deg = 15.);

void removeBadPsiHelicalSegmentsFromDataStar(
		FileName& fn_in,
		FileName& fn_out,
		RFLOAT max_dev_deg = 15.);

void convertHelicalSegmentCoordsToStarFile_Multiple(
		FileName& suffix_coords,
		FileName& suffix_out,
		int format_tag,
		RFLOAT pixel_size_A,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT boxsize,
		bool bimodal_angular_priors = true);

void convertHelicalSegmentCoordsToMetaDataTable(
		FileName& fn_in,
		MetaDataTable& MD_out,
		int& total_segments,
		bool is_3D,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT Zdim,
		RFLOAT box_size_pix,
		bool bimodal_angular_priors = true);

void convertXimdispHelicalSegmentCoordsToMetaDataTable(
		FileName& fn_in,
		MetaDataTable& MD_out,
		int& total_segments,
		int& total_tubes,
		RFLOAT pixel_size_A,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix,
		bool bimodal_angular_priors = true);

void convertXimdispHelicalTubeCoordsToMetaDataTable(
		FileName& fn_in,
		MetaDataTable& MD_out,
		int& total_segments,
		int& total_tubes,
		int nr_asu,
		RFLOAT rise_A,
		RFLOAT pixel_size_A,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix,
		bool bimodal_angular_priors = true,
		bool cut_into_segments = true);

void convertEmanHelicalSegmentCoordsToMetaDataTable(
		FileName& fn_in,
		MetaDataTable& MD_out,
		int& total_segments,
		int& total_tubes,
		RFLOAT pixel_size_A,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix,
		bool bimodal_angular_priors = true);

void convertEmanHelicalTubeCoordsToMetaDataTable(
		FileName& fn_in,
		MetaDataTable& MD_out,
		int& total_segments,
		int& total_tubes,
		int nr_asu,
		RFLOAT rise_A,
		RFLOAT pixel_size_A,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix,
		bool bimodal_angular_priors = true,
		bool cut_into_segments = true);

void makeHelicalReference2D(
		MultidimArray<RFLOAT>& out,
		int box_size,
		RFLOAT particle_diameter_A,
		RFLOAT tube_diameter_A,
		RFLOAT pixel_size_A,
		bool is_tube_white = true);

void makeHelicalReference3DWithPolarity(
		MultidimArray<RFLOAT>& out,
		int box_size,
		RFLOAT pixel_size_A,
		RFLOAT twist_deg,
		RFLOAT rise_A,
		RFLOAT tube_diameter_A,
		RFLOAT particle_diameter_A,
		RFLOAT cyl_diameter_A,
		RFLOAT topbottom_ratio,
		int sym_Cn = 1,
		int nr_filaments_helix_with_seam = -1);

void divideStarFile(
		FileName& fn_in,
		int nr);

void mergeStarFiles(FileName& fn_in);

void sortHelicalTubeID(MetaDataTable& MD);

void simulateHelicalSegments(
		bool is_3d_tomo,
		FileName& fn_vol_in,
		FileName& fn_star_out,
		RFLOAT white_noise,
		int new_box_size,
		int nr_subunits,
		int nr_asu,
		int nr_tubes,
		bool do_bimodal_searches,
		RFLOAT cyl_outer_diameter_A,
		RFLOAT angpix,
		RFLOAT rise_A,
		RFLOAT twist_deg,
		RFLOAT sigma_psi,
		RFLOAT sigma_tilt,
		RFLOAT sigma_offset,
		int random_seed = -1);

void outputHelicalSymmetryStatus(
		int iter,
		RFLOAT rise_initial_A,
		RFLOAT rise_min_A,
		RFLOAT rise_max_A,
		RFLOAT twist_initial_deg,
		RFLOAT twist_min_deg,
		RFLOAT twist_max_deg,
		bool do_local_search_helical_symmetry,
		std::vector<RFLOAT>& rise_A,
		std::vector<RFLOAT>& twist_deg,
		RFLOAT rise_A_half1,
		RFLOAT rise_A_half2,
		RFLOAT twist_deg_half1,
		RFLOAT twist_deg_half2,
		bool do_split_random_halves,
		std::ostream& out);

void excludeLowCTFCCMicrographs(
		FileName& fn_in,
		FileName& fn_out,
		RFLOAT cc_min = -1.,
		RFLOAT EPA_lowest_res = 999999.,
		RFLOAT df_min = -999999.,
		RFLOAT df_max = 999999.);

void cutOutPartOfHelix(
		const MultidimArray<RFLOAT>& vin,
		MultidimArray<RFLOAT>& vout,
		long int new_boxdim,
		RFLOAT ang_deg,
		RFLOAT z_percentage);

// TESTING...
class HelicalSegmentPriorInfoEntry
{
public:
	std::string helical_tube_name;
	long int MDobjectID;
	RFLOAT rot_deg, psi_deg, tilt_deg;
	RFLOAT dx_A, dy_A, dz_A;
	RFLOAT track_pos_A;
	bool has_wrong_polarity;
	int subset, classID;

	RFLOAT rot_prior_deg, psi_prior_deg, tilt_prior_deg;
	RFLOAT rot_prior_deg_ori, psi_prior_deg_ori, tilt_prior_deg_ori;
	RFLOAT dx_prior_A, dy_prior_A, dz_prior_A;
	RFLOAT psi_flip_ratio;
	bool psi_prior_flip;

	void clear();

	HelicalSegmentPriorInfoEntry() { clear(); };

	~HelicalSegmentPriorInfoEntry() { clear(); };

	void checkPsiPolarity();

	bool operator<(const HelicalSegmentPriorInfoEntry &rhs) const;
};

// KThurber add this function
void flipPsiTiltForHelicalSegment(
		RFLOAT old_psi,
		RFLOAT old_tilt,
		RFLOAT& new_psi,
		RFLOAT& new_tilt);

void updatePriorsForOneHelicalTube(
		std::vector<HelicalSegmentPriorInfoEntry>& list,
		int sid,
		int eid,
		int& nr_wrong_polarity,
		bool &reverse_direction,
		RFLOAT sigma_segment_dist,
		std::vector<RFLOAT> helical_rise,
		std::vector<RFLOAT> helical_twist,
		bool is_3D_data,
		bool do_auto_refine,
		RFLOAT sigma2_rot,       // KThurber
		RFLOAT sigma2_tilt,
		RFLOAT sigma2_psi,
		RFLOAT sigma2_offset,
		RFLOAT sigma_cutoff = 3.);

void updatePriorsForHelicalReconstruction(
		MetaDataTable& MD,
		RFLOAT sigma_segment_dist,
		std::vector<RFLOAT> helical_rise,
		std::vector<RFLOAT> helical_twist,
		int helical_nstart,
		bool is_3D,
		bool do_auto_refine,
		RFLOAT sigma2_rot,
		RFLOAT sigma2_tilt,
		RFLOAT sigma2_psi,
		RFLOAT sigma2_offset,
		bool keep_tilt_prior_fixed,
		int verb);

void updateAngularPriorsForHelicalReconstructionFromLastIter(
		MetaDataTable& MD,
		bool keep_tilt_prior_fixed);

void testDataFileTransformXY(MetaDataTable& MD);

void setPsiFlipRatioInStarFile(MetaDataTable& MD, RFLOAT ratio = 0.);

void plotLatticePoints(MetaDataTable& MD,
		int x1, int y1, int x2, int y2);


void grabParticleCoordinates(
		FileName& fn_in,
		FileName& fn_out);

// Files of priors: mic1_priors.star, files of local CTF: mic1_local.star
// Then suffix_priors = _priors.star, suffix_local_ctf = _local.star
void grabParticleCoordinates_Multiple(
		std::string& suffix_fin,
		std::string& suffix_fout);

void calculateRadialAvg(MultidimArray<RFLOAT> &v, RFLOAT angpix);

void transformCartesianToHelicalCoordsForStarFiles(
		MetaDataTable& MD_in,
		MetaDataTable& MD_out);

void normaliseHelicalSegments(
		FileName& fn_in,
		FileName& fn_out_root,
		RFLOAT helical_outer_diameter_A,
		RFLOAT pixel_size_A);

// Copied online from http://paulbourke.net/miscellaneous/interpolation/
// Author: Paul Bourke, December 1999
/*
   Tension: 1 is high, 0 normal, -1 is low
   Bias: 0 is even,
         positive is towards first segment,
         negative towards the other
*/
// mu is the percentage between y1 and y2
RFLOAT HermiteInterpolate1D(
		RFLOAT y0, RFLOAT y1, RFLOAT y2, RFLOAT y3,
		RFLOAT mu,
		RFLOAT tension = 0.,
		RFLOAT bias = 0.);

void HermiteInterpolateOne3DHelicalFilament(
		MetaDataTable& MD_in,
		MetaDataTable& MD_out,
		int& total_segments,
		int nr_asu,
		RFLOAT rise_A,
		RFLOAT pixel_size_A,
		RFLOAT box_size_pix,
		int helical_tube_id,
		RFLOAT Xdim = 999999.,
		RFLOAT Ydim = 999999.,
		RFLOAT Zdim = 999999.,
		bool bimodal_angular_priors = true);

void Interpolate3DCurves(
		FileName& fn_in_root,
		FileName& fn_out_root,
		int nr_asu,
		RFLOAT rise_A,
		RFLOAT pixel_size_A,
		RFLOAT box_size_pix,
		int binning_factor = 1,
		bool bimodal_angular_priors = true);

void estimateTiltPsiPriors(
		Matrix1D<RFLOAT>& dr,
		RFLOAT& tilt_deg,
		RFLOAT& psi_deg);

void readFileHeader(
		FileName& fn_in,
		FileName& fn_out,
		int nr_bytes = 100);

void select3DsubtomoFrom2Dproj(
		MetaDataTable& MD_2d,
		MetaDataTable& MD_3d,
		MetaDataTable& MD_out);

void averageAsymmetricUnits2D(
		ObservationModel& obsModel,
		MetaDataTable &MDimgs,
		FileName fn_o_root,
		int nr_asu,
		RFLOAT rise);

#endif /* HELIX_H_ */
