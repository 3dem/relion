/***************************************************************************
 *
 * Author: "Shaoda He"
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

#define CART_TO_HELICAL_COORDS true
#define HELICAL_TO_CART_COORDS false

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
		return;
	}
	HelicalSymmetryItem(RFLOAT _twist_deg, RFLOAT _rise_pix, RFLOAT _dev = (1e35))
	{
		twist_deg = _twist_deg;
		rise_pix = _rise_pix;
		dev = _dev;
		return;
	}
};

bool IsEqualHelicalSymmetry(
		const HelicalSymmetryItem& a,
		const HelicalSymmetryItem& b);

bool IsHelicalSymmetryWithinOpenInterval(
		RFLOAT rise_pix,
		RFLOAT twist_deg,
		RFLOAT rise_ori_pix,
		RFLOAT twist_ori_deg,
		RFLOAT rise_half_range_pix,
		RFLOAT twist_half_range_deg);

void makeHelicalSymmetryList(
		std::vector<HelicalSymmetryItem>& list,
		RFLOAT rise_ori_pix,
		RFLOAT twist_ori_deg,
		RFLOAT rise_half_range_pix,
		RFLOAT twist_half_range_deg,
		RFLOAT rise_step_pix,
		RFLOAT twist_step_deg);

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
		RFLOAT r_min_A,
		RFLOAT r_max_A,
		RFLOAT z_percentage,
		RFLOAT rise_ori_A,
		RFLOAT twist_ori_deg,
		RFLOAT rise_search_max_dev_percentage,
		RFLOAT twist_search_max_dev_percentage,
		RFLOAT& rise_refined_A,
		RFLOAT& twist_refined_deg);








// Impose helical symmetry in Fourier space (It takes too much time and memory...)
/*
bool get2DZsliceIn3DVolume(
		const MultidimArray<RFLOAT>& vol_in,
		MultidimArray<RFLOAT>& img_out,
		int idz);

bool add2DZsliceInto3DVolumeSum(
		const MultidimArray<RFLOAT>& img_in,
		MultidimArray<RFLOAT>& vol_sum,
		std::vector<RFLOAT>& weight_sum,
		int idz);

void shift3DVolumeAlongZAxisInFourierSpace(
		MultidimArray<RFLOAT>& img,
		RFLOAT shift_pix);

void rotateAndSum2DZSliceInRealSpace(
		MultidimArray<RFLOAT>& img_ori,
		MultidimArray<RFLOAT>& img_sum,
		std::vector<RFLOAT>& weight_sum,
		int idz,
		RFLOAT outer_radius_pix,
		RFLOAT rot_angle_deg);

void rotate2DZSliceInFourierSpace(
		MultidimArray<RFLOAT>& img,
		RFLOAT rot_angle_deg,
		int padding_factor = 2);

// In Z axis [s, e] chunk, there must be a symmetrical segment of helix with at least one particle!
void expandZaxisInFourierSpace(
		MultidimArray<RFLOAT>& vol,
		RFLOAT outer_radius_pix,
		RFLOAT twist_deg,  // both + or -
		RFLOAT rise_pix,  // only +
		int idz_s,
		int idz_e,
		int padding_factor = 2);

void imposeHelicalSymmetryInFourierSpace(MultidimArray<RFLOAT>& vol,
		RFLOAT outer_radius_pix,
		RFLOAT twist_deg,  // both + or -
		RFLOAT rise_pix,  // only +
		int idz_s,
		int idz_e,
		int padding_factor = 2);
*/

RFLOAT getHelicalSigma2Rot(
		RFLOAT helical_rise_pix,
		RFLOAT helical_twist_deg,
		RFLOAT helical_offset_step_pix,
		RFLOAT rot_step_deg,
		RFLOAT old_sigma2_rot);

// Impose helical symmetry in real space
void imposeHelicalSymmetryInRealSpace(
		MultidimArray<RFLOAT>& vol,
		RFLOAT sphere_radius_pix,
		RFLOAT cyl_inner_radius_pix,
		RFLOAT cyl_outer_radius_pix,
		RFLOAT cosine_width,
		RFLOAT twist_deg,  // both + or -
		RFLOAT rise_pix,  // only +
		int idz_s,
		int idz_e);

// Assume all parameters are within range
RFLOAT get_lenZ_percentage_max(
		int box_len,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_outer_radius_A,
		RFLOAT pixel_size_A);

// Assume all parameters are within range
RFLOAT get_rise_A_max(
		int box_len,
		RFLOAT pixel_size_A,
		RFLOAT lenZ_percentage,
		RFLOAT nr_units_min);

bool checkHelicalParametersFor3DHelicalReference(
		int box_len,
		RFLOAT pixel_size_A,
		RFLOAT twist_deg,
		RFLOAT rise_A,
		RFLOAT lenZ_percentage,
		bool do_helical_symmetry_local_refinement,
		RFLOAT rise_max_dev_percentage,
		RFLOAT twist_max_dev_percentage,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_inner_radius_A,
		RFLOAT cyl_outer_radius_A);

void makeHelicalReferenceInRealSpace(
		MultidimArray<RFLOAT>& vol,
		RFLOAT pixel_size_A,
		RFLOAT twist_deg,
		RFLOAT rise_A,
		RFLOAT lenZ_percentage,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_inner_radius_A,
		RFLOAT cyl_outer_radius_A,
		RFLOAT cosine_width_pix);

// Search for helical symmetry (not done!)
/*
bool calcCCOfCnZSymmetry(
		const MultidimArray<RFLOAT>& v,
		RFLOAT r_min_pix,
		RFLOAT r_max_pix,
		int cn,
		RFLOAT& cc,
		int& nr_asym_voxels);

void searchCnZSymmetry(
		const MultidimArray<RFLOAT>& v,
		RFLOAT r_min_pix,
		RFLOAT r_max_pix,
		int cn_start,
		int cn_end,
		std::vector<RFLOAT>& cn_list,
		std::vector<RFLOAT>& cc_list,
		std::vector<int>& nr_asym_voxels_list,
		std::ofstream* fout_ptr = NULL);
*/

/*
bool localSearchHelicalTwist(
		const MultidimArray<RFLOAT>& v,
		RFLOAT r_min_pix,
		RFLOAT r_max_pix,
		RFLOAT z_percentage,
		RFLOAT rise_pix,
		RFLOAT twist_deg,
		RFLOAT search_half_range_deg,
		RFLOAT search_step_deg,
		RFLOAT& twist_refined_deg);

bool localSearchHelicalRise(
		const MultidimArray<RFLOAT>& v,
		RFLOAT r_min_pix,
		RFLOAT r_max_pix,
		RFLOAT z_percentage,
		RFLOAT rise_pix,
		RFLOAT twist_deg,
		RFLOAT search_half_range_pix,
		RFLOAT search_step_pix,
		RFLOAT& rise_refined_pix);

bool localSearchHelicalSymmetry(
		const MultidimArray<RFLOAT>& v,
		RFLOAT pixel_size_A,
		RFLOAT sphere_radius_A,
		RFLOAT r_min_A,
		RFLOAT r_max_A,
		RFLOAT z_percentage,
		RFLOAT rise_ori_A,
		RFLOAT twist_ori_deg,
		RFLOAT local_search_max_dev_percentage,
		RFLOAT& rise_refined_A,
		RFLOAT& twist_refined_deg);
*/

RFLOAT calcCCofPsiFor2DHelicalSegment(
		const MultidimArray<RFLOAT>& v,
		RFLOAT psi_deg,
		RFLOAT pixel_size_A,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_outer_radius_A);

RFLOAT localSearchPsiFor2DHelicalSegment(
		const MultidimArray<RFLOAT>& v,
		RFLOAT pixel_size_A,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_outer_radius_A,
		RFLOAT ori_psi_deg,
		RFLOAT search_half_range_deg,
		RFLOAT search_step_deg);

RFLOAT searchPsiFor2DHelicalSegment(
		const MultidimArray<RFLOAT>& v,
		RFLOAT pixel_size_A,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_outer_radius_A);

/*
void searchHelicalSymmetry(
		const MultidimArray<RFLOAT>& v,
		RFLOAT r_min_pix,
		RFLOAT r_max_pix,
		RFLOAT rise_pix_start,
		RFLOAT rise_pix_incr,
		int rise_pix_steps,
		RFLOAT twist_deg_start,
		RFLOAT twist_deg_incr,
		int twist_deg_steps,
		std::vector<RFLOAT>& rise_pix_list,
		std::vector<RFLOAT>& twist_deg_list,
		std::vector<RFLOAT>& cc_list,
		std::vector<int>& nr_asym_voxels_list,
		std::ofstream* fout_ptr = NULL);
*/

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

void transformCartesianAndHelicalCoords(
		Matrix1D<RFLOAT>& in,
		Matrix1D<RFLOAT>& out,
		RFLOAT psi,
		RFLOAT tilt,
		bool direction);

void transformCartesianAndHelicalCoords(
		RFLOAT xin,
		RFLOAT yin,
		RFLOAT zin,
		RFLOAT& xout,
		RFLOAT& yout,
		RFLOAT& zout,
		RFLOAT psi,
		RFLOAT tilt,
		int dim,
		bool direction);

void makeBlot(
		MultidimArray<RFLOAT>& v,
		RFLOAT y,
		RFLOAT x,
		RFLOAT r);

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

/*
void normalise2DImageSlices(
		const FileName& fn_in,
		const FileName& fn_out,
		int bg_radius,
		RFLOAT white_dust_stddev = -1.,
		RFLOAT black_dust_stddev = -1.);
*/

void enlarge3DReference(
		MultidimArray<RFLOAT>& v,
		int factor);

void applySoftSphericalMask(
		MultidimArray<RFLOAT>& v,
		RFLOAT sphere_diameter = -1.,
		RFLOAT cosine_width = 5.);

int extractCoordsForAllHelicalSegments(
		FileName& fn_in,
		FileName& fn_out,
		int nr_asu,
		RFLOAT rise_ang,
		RFLOAT pixel_size_ang,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix);

// Files of input coordinates of tops and bottoms of helices: mic1_manualpick.star
// Files of output coordinates of segments: mic1_all.star
// Then suffix_in = _manualpick.star, suffix_out = _all.star
void extractCoordsForAllHelicalSegments_Multiple(
		FileName& suffix_in,
		FileName& suffix_out,
		int nr_asu,
		RFLOAT rise_ang,
		RFLOAT pixel_size_ang,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix);

void combineParticlePriorsWithKaiLocalCTF(
		FileName& fn_priors,
		FileName& fn_local_ctf,
		FileName& fn_combined);

void addPriorsToParticleDataFile(
		FileName& fn_priors,
		FileName& fn_data,
		FileName& fn_out);

// Files of priors: mic1_priors.star, files of local CTF: mic1_local.star
// Then suffix_priors = _priors.star, suffix_local_ctf = _local.star
void combineParticlePriorsWithKaiLocalCTF_Multiple(
		std::string& suffix_priors,
		std::string& suffix_local_ctf,
		std::string& suffix_combined);

void setNullAlignmentPriorsInDataStar(
		FileName& fn_in,
		FileName& fn_out,
		bool rewrite_tilt = false,
		bool rewrite_psi = false);

void removeBadTiltHelicalSegmentsFromDataStar(
		FileName& fn_in,
		FileName& fn_out,
		RFLOAT max_dev_deg = 15.);

int transformXimdispHelicalCoordsToStarFile(
		FileName& fn_in,
		FileName& fn_out,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix);

// Files of Ximdisp coordinates: 001.coords, files of output coordinates: 001_all.star
// Then suffix_coords = .coords, suffix_out = _all.star
void transformXimdispHelicalCoordsToStarFile_Multiple(
		FileName& suffix_coords,
		FileName& suffix_out,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT boxsize);

/*
void divideHelicalSegmentsFromMultipleMicrographsIntoHalves(
		FileName& fn_in,
		FileName& fn_out);
*/

void divideHelicalSegmentsFromMultipleMicrographsIntoRandomHalves(
		FileName& fn_in,
		FileName& fn_out,
		int random_seed = -1);

void makeHelicalReference2D(
		MultidimArray<RFLOAT>& out,
		int box_size,
		RFLOAT particle_diameter_A,
		RFLOAT tube_diameter_A,
		RFLOAT pixel_size_A,
		bool is_tube_white = true);

void makeHelicalReference3D(
		MultidimArray<RFLOAT>& out,
		int box_size,
		RFLOAT pixel_size_A,
		RFLOAT twist_deg,
		RFLOAT rise_A,
		RFLOAT tube_diameter_A,
		RFLOAT particle_diameter_A,
		int sym_Cn);

void makeHelicalReconstructionStarFileFromSingle2DClassAverage(
		FileName& fn_in_class2D,
		FileName& fn_out_star,
		RFLOAT pixel_size_A,
		RFLOAT twist_deg,
		RFLOAT rise_A,
		RFLOAT tilt_deg,
		RFLOAT psi_deg,
		int nr_projections);

void divideStarFile(
		FileName& fn_in,
		int nr);

void combineStarFiles(
		FileName& fn_in);

void sortHelicalTubeID(
		MetaDataTable& MD);

void simulateHelicalSegments(
		FileName& fn_out,
		int nr_subunits,
		int nr_asu,
		int nr_tubes,
		RFLOAT twist_deg,
		RFLOAT sigma_psi = 1.,
		RFLOAT sigma_tilt = 1.);

#endif /* HELIX_H_ */
