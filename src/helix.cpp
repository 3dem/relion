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

#include "src/macros.h"
#include "src/helix.h"

// Wider search ranges for helical twist and rise
#define WIDE_HELICAL_TWIST_AND_RISE_SEARCHES

// Exclude segments close to the edges of the 2D micrographs / 3D tomograms. Please switch it on.
#define EXCLUDE_SEGMENTS_ON_THE_EDGES

void makeHelicalSymmetryList(
		std::vector<HelicalSymmetryItem>& list,
		RFLOAT rise_min_pix,
		RFLOAT rise_max_pix,
		RFLOAT rise_step_pix,
		bool search_rise,
		RFLOAT twist_min_deg,
		RFLOAT twist_max_deg,
		RFLOAT twist_step_deg,
		bool search_twist)
{
	// Assume all parameters are within range
	RFLOAT rise_pix, twist_deg;
	int rise_samplings, twist_samplings;
	std::vector<HelicalSymmetryItem> tmp_list;

	if (search_rise)
		rise_samplings = ROUND(fabs(rise_max_pix - rise_min_pix) / fabs(rise_step_pix));
	else
	{
		rise_min_pix = (rise_min_pix + rise_max_pix) / 2.;
		rise_samplings = 0;
	}
	if (search_twist)
		twist_samplings = ROUND(fabs(twist_max_deg - twist_min_deg) / fabs(twist_step_deg));
	else
	{
		twist_min_deg = (twist_min_deg + twist_max_deg) / 2.;
		twist_samplings = 0;
	}

	// Store a matrix of symmetries
	tmp_list.clear();
	for (int ii = 0; ii <= rise_samplings; ii++)
	{
		for (int jj = 0; jj <= twist_samplings; jj++)
		{
			rise_pix = rise_min_pix + rise_step_pix * RFLOAT(ii);
			twist_deg = twist_min_deg + twist_step_deg * RFLOAT(jj);
			tmp_list.push_back(HelicalSymmetryItem(twist_deg, rise_pix));
		}
	}

	// Check duplications and return this matrix
	RFLOAT twist_dev_deg, twist_avg_deg, rise_dev_pix, rise_avg_pix;
	RFLOAT err_max = (1e-10);
	bool same_twist, same_rise;
	for (int ii = 0; ii < list.size(); ii++)
	{
		for (int jj = 0; jj < tmp_list.size(); jj++)
		{
			same_twist = same_rise = false;

			twist_dev_deg = fabs(list[ii].twist_deg - tmp_list[jj].twist_deg);
			rise_dev_pix = fabs(list[ii].rise_pix - tmp_list[jj].rise_pix);
			twist_avg_deg = (fabs(list[ii].twist_deg) + fabs(tmp_list[jj].twist_deg)) / 2.;
			rise_avg_pix = (fabs(list[ii].rise_pix) + fabs(tmp_list[jj].rise_pix)) / 2.;

			if (twist_avg_deg < err_max)
			{
				if (twist_dev_deg < err_max)
					same_twist = true;
			}
			else
			{
				if ((twist_dev_deg / twist_avg_deg) < err_max)
					same_twist = true;
			}

			if (rise_avg_pix < err_max)
			{
				if (rise_dev_pix < err_max)
					same_rise = true;
			}
			else
			{
				if ( (rise_dev_pix / rise_avg_pix) < err_max)
					same_rise = true;
			}

			if (same_twist && same_rise)
				tmp_list[jj].dev = list[ii].dev;
		}
	}
	list.clear();
	list = tmp_list;
	return;
};

bool calcCCofHelicalSymmetry(
		const MultidimArray<RFLOAT>& v,
		RFLOAT r_min_pix,
		RFLOAT r_max_pix,
		RFLOAT z_percentage,
		RFLOAT rise_pix,
		RFLOAT twist_deg,
		RFLOAT& cc,
		int& nr_asym_voxels)
{
	// TODO: go through every line to find bugs!!!
	int rec_len, r_max_XY, startZ, finishZ;
	double dist_r_pix, sum_pw1, sum_pw2, sum_n, sum_chunk = 0., sum_chunk_n = 0.;
	//std::vector<RFLOAT> sin_rec, cos_rec, dev_voxel, dev_chunk;
	std::vector<RFLOAT> sin_rec, cos_rec;

	if ( (STARTINGZ(v) != FIRST_XMIPP_INDEX(ZSIZE(v))) || (STARTINGY(v) != FIRST_XMIPP_INDEX(YSIZE(v))) || (STARTINGX(v) != FIRST_XMIPP_INDEX(XSIZE(v))) )
		REPORT_ERROR("helix.cpp::calcCCofHelicalSymmetry(): The origin of input 3D MultidimArray is not at the center (use v.setXmippOrigin() before calling this function)!");

	// Check r_max
	r_max_XY = (XSIZE(v) < YSIZE(v)) ? XSIZE(v) : YSIZE(v);
	r_max_XY = (r_max_XY + 1) / 2 - 1;
	if ( r_max_pix > (((RFLOAT)(r_max_XY)) - 0.01) )  // 0.01 - avoid segmentation fault
		r_max_pix = (((RFLOAT)(r_max_XY)) - 0.01);

	// Set startZ and finishZ
	startZ = FLOOR( (-1.) * ((RFLOAT)(ZSIZE(v)) * z_percentage * 0.5) );
	finishZ = CEIL( ((RFLOAT)(ZSIZE(v))) * z_percentage * 0.5 );
	startZ = (startZ <= (STARTINGZ(v))) ? (STARTINGZ(v) + 1) : (startZ);
	finishZ = (finishZ >= (FINISHINGZ(v))) ? (FINISHINGZ(v) - 1) : (finishZ);

	// Calculate tabulated sine and cosine values
	rec_len = 2 + (CEIL((RFLOAT(ZSIZE(v)) + 2.) / rise_pix));
	sin_rec.clear();
	cos_rec.clear();
	sin_rec.resize(rec_len);
	cos_rec.resize(rec_len);
	for (int id = 0; id < rec_len; id++)
#ifdef RELION_SINGLE_PRECISION
		SINCOSF(DEG2RAD(((RFLOAT)(id)) * twist_deg), &sin_rec[id], &cos_rec[id]);
#else
		SINCOS(DEG2RAD(((RFLOAT)(id)) * twist_deg), &sin_rec[id], &cos_rec[id]);
#endif

	rise_pix = fabs(rise_pix);

	// Test a chunk of Z length = rise
	//dev_chunk.clear();
	// Iterate through all coordinates on Z, Y and then X axes
	FOR_ALL_ELEMENTS_IN_ARRAY3D(v)
	{
		RFLOAT xp, yp, zp, fx, fy, fz;

		// Test a chunk of Z length = rise
		// for(idz = startZ; (idz <= (startZ + ((int)(floor(rise_pix))))) && (idz <= finishZ); idz++)
		if ( (k < startZ) || (k > (startZ + (FLOOR(rise_pix)))) || (k > finishZ) )
			continue;

		dist_r_pix = sqrt(i * i + j * j);
		if ( (dist_r_pix < r_min_pix) || (dist_r_pix > r_max_pix) )
			continue;

		// Pick a voxel in the chunk
		//dev_voxel.clear();
		//dev_voxel.push_back(A3D_ELEM(v, k, i, j));

		// Pick other voxels according to this voxel and helical symmetry
		zp = k;
		int rot_id = 0;
		sum_pw1 = A3D_ELEM(v, k, i, j);
		sum_pw2 = A3D_ELEM(v, k, i, j)*A3D_ELEM(v, k, i, j);
		sum_n = 1.;
		while (1)
		{
			// Rise
			zp += rise_pix;
			if (zp > finishZ) // avoid segmentation fault - finishZ is always strictly smaller than FINISHINGZ(v)!
				break;

			// Twist
			rot_id++;
			xp = ((RFLOAT)(j)) * cos_rec[rot_id] - ((RFLOAT)(i)) * sin_rec[rot_id];
			yp = ((RFLOAT)(j)) * sin_rec[rot_id] + ((RFLOAT)(i)) * cos_rec[rot_id];

			// Trilinear interpolation (with physical coords)
			// Subtract STARTINGX,Y,Z to accelerate access to data
			// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
			int x0, y0, z0, x1, y1, z1;
			x0 = FLOOR(xp); fx = xp - x0; x0 -= STARTINGX(v); x1 = x0 + 1;
			y0 = FLOOR(yp); fy = yp - y0; y0 -= STARTINGY(v); y1 = y0 + 1;
			z0 = FLOOR(zp); fz = zp - z0; z0 -= STARTINGZ(v); z1 = z0 + 1;
			// DEBUG
			if ( (x0 < 0) || (y0 < 0) || (z0 < 0) || (x1 >= XSIZE(v)) || (y1 >= YSIZE(v)) || (z1 >= ZSIZE(v)) )
				std::cout << " idzidyidx= " << k << ", " << i << ", " << j << ", x0x1y0y1z0z1= " << x0 << ", " << x1 << ", " << y0 << ", " << y1 << ", " << z0 << ", " << z1 << std::endl;

			RFLOAT d000, d001, d010, d011, d100, d101, d110, d111;
			d000 = DIRECT_A3D_ELEM(v, z0, y0, x0);
			d001 = DIRECT_A3D_ELEM(v, z0, y0, x1);
			d010 = DIRECT_A3D_ELEM(v, z0, y1, x0);
			d011 = DIRECT_A3D_ELEM(v, z0, y1, x1);
			d100 = DIRECT_A3D_ELEM(v, z1, y0, x0);
			d101 = DIRECT_A3D_ELEM(v, z1, y0, x1);
			d110 = DIRECT_A3D_ELEM(v, z1, y1, x0);
			d111 = DIRECT_A3D_ELEM(v, z1, y1, x1);

			RFLOAT dx00, dx01, dx10, dx11;
			dx00 = LIN_INTERP(fx, d000, d001);
			dx01 = LIN_INTERP(fx, d100, d101);
			dx10 = LIN_INTERP(fx, d010, d011);
			dx11 = LIN_INTERP(fx, d110, d111);

			RFLOAT dxy0, dxy1, ddd;
			dxy0 = LIN_INTERP(fy, dx00, dx10);
			dxy1 = LIN_INTERP(fy, dx01, dx11);

			ddd = LIN_INTERP(fz, dxy0, dxy1);

			// Record this voxel
			sum_pw1 += ddd;
			sum_pw2 += ddd * ddd;
			sum_n += 1.;

			// dev_voxel.push_back(ddd);
		}

		sum_pw1 /= sum_n;
		sum_pw2 /= sum_n;
		//dev_chunk.push_back(sum_pw2 - sum_pw1 * sum_pw1);
		// Sum_chunk
		sum_chunk += sum_pw2 - sum_pw1 * sum_pw1;
		sum_chunk_n += 1.;

		/*
		// Calc dev of this voxel in the chunk
		if (dev_voxel.size() > 1)
		{
			sum_pw1 = sum_pw2 = 0.;
			for (int id = 0; id < dev_voxel.size(); id++)
			{
				sum_pw1 += dev_voxel[id];
				sum_pw2 += dev_voxel[id] * dev_voxel[id];
			}
			sum_pw1 /= dev_voxel.size();
			sum_pw2 /= dev_voxel.size();
			// TODO: record stddev or dev???
			dev_chunk.push_back(sum_pw2 - sum_pw1 * sum_pw1);
		}
		dev_voxel.clear();
		*/
	}

	// Calc avg of all voxels' devs in this chunk (for a specific helical symmetry)
	if (sum_chunk_n < 1)
	{
		cc = (1e10);
		nr_asym_voxels = 0;
		return false;
	}
	else
	{
		cc = (sum_chunk / sum_chunk_n);
	}
	nr_asym_voxels = sum_chunk_n;
	//dev_chunk.clear();

	return true;
};

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
		std::ostream* o_ptr)
{
	// TODO: whether iterations can exit & this function works for negative twist
	int iter, box_len, nr_asym_voxels, nr_rise_samplings, nr_twist_samplings, nr_min_samplings, nr_max_samplings, best_id, iter_not_converged;
	RFLOAT r_min_pix, r_max_pix, best_dev, err_max;
	RFLOAT rise_min_pix, rise_max_pix, rise_step_pix, rise_inistep_pix, twist_step_deg, rise_refined_pix;
	RFLOAT rise_local_min_pix, rise_local_max_pix, twist_local_min_deg, twist_local_max_deg;
	std::vector<HelicalSymmetryItem> helical_symmetry_list;
	bool out_of_range, search_rise, search_twist;

	// Check input 3D reference
	if (v.getDim() != 3)
		REPORT_ERROR("helix.cpp::localSearchHelicalSymmetry(): Input helical reference is not 3D! (v.getDim() = " + integerToString(v.getDim()) + ")");

	// Set the length of the box
	box_len = (XSIZE(v) < YSIZE(v)) ? XSIZE(v) : YSIZE(v);
	box_len = (box_len < ZSIZE(v)) ? box_len : ZSIZE(v);

	// Initialise refined helical parameters
	// Check helical parameters
	rise_refined_A = (rise_min_A + rise_max_A) / 2.;
	twist_refined_deg = (twist_min_deg + twist_max_deg) / 2.;
	checkParametersFor3DHelicalReconstruction(
			false,
			true,
			1,
			rise_refined_A,
			rise_min_A,
			rise_max_A,
			twist_refined_deg,
			twist_min_deg,
			twist_max_deg,
			box_len,
			pixel_size_A,
			z_percentage,
			sphere_radius_A * 2.,
			cyl_inner_radius_A * 2.,
			cyl_outer_radius_A * 2.);
	rise_refined_pix = rise_refined_A / pixel_size_A;

	// Initialise other parameters
	out_of_range = false;
	r_min_pix = cyl_inner_radius_A / pixel_size_A;
	r_max_pix = cyl_outer_radius_A / pixel_size_A;
	rise_inistep_pix = rise_inistep_A / pixel_size_A;
	rise_local_min_pix = rise_min_pix = rise_min_A / pixel_size_A;
	rise_local_max_pix = rise_max_pix = rise_max_A / pixel_size_A;
	twist_local_min_deg = twist_min_deg;
	twist_local_max_deg = twist_max_deg;
	if (o_ptr != NULL)
	{
		(*o_ptr) << " ### RELION helix toolbox - local searches of helical symmetry" << std::endl;
		(*o_ptr) << " --> Box size = " << ((long int)(XSIZE(v))) << ", Z(%) = " << (z_percentage * 100.)
				<< "%, pixel size = " << pixel_size_A << " Angstroms, inner diameter = " << (cyl_inner_radius_A * 2.)
				<< " Angstroms, outer diameter = " << (cyl_outer_radius_A * 2.) << " Angstroms." << std::endl;
		(*o_ptr) << " --> Searching twist from " << twist_min_deg << " to " << twist_max_deg << " degrees, rise from "
				<< rise_min_A << " to " << rise_max_A << " Angstroms." << std::endl;
	}

	// Initial searches - Iteration 1
	// Sampling of twist should be smaller than 1 degree
	// Sampling of rise should be smaller than 1%
	// And also make sure to search for at least 5*5 sampling points
	// Avoid too many searches (at most 1000*1000)
	err_max = (1e-5);
	search_rise = search_twist = true;
	nr_min_samplings = 5;
	nr_max_samplings = 1000;

	twist_inistep_deg = (twist_inistep_deg < (1e-5)) ? (1.) : (twist_inistep_deg);
	nr_twist_samplings = CEIL(fabs(twist_local_min_deg - twist_local_max_deg) / twist_inistep_deg);
	nr_twist_samplings = (nr_twist_samplings > nr_min_samplings) ? (nr_twist_samplings) : (nr_min_samplings);
	nr_twist_samplings = (nr_twist_samplings < nr_max_samplings) ? (nr_twist_samplings) : (nr_max_samplings);
	twist_step_deg = fabs(twist_local_min_deg - twist_local_max_deg) / RFLOAT(nr_twist_samplings);
	if ( fabs(twist_local_min_deg - twist_local_max_deg) < err_max)
	{
		search_twist = false;
		twist_step_deg = 0.;
		twist_min_deg = twist_max_deg = twist_local_min_deg = twist_local_max_deg = twist_refined_deg;
	}
	if (o_ptr != NULL)
	{
		if (search_twist)
			(*o_ptr) << " --> Initial searching step of twist is " << twist_step_deg << " degrees (" << nr_twist_samplings << " samplings)." << std::endl;
		else
			(*o_ptr) << " --> No need to search for twist..." << std::endl;
	}

	rise_inistep_pix = (rise_inistep_pix < (1e-5)) ? (1e30) : (rise_inistep_pix);
	rise_step_pix = 0.01 * ((fabs(rise_local_min_pix) + fabs(rise_local_max_pix)) / 2.);
	rise_step_pix = (rise_step_pix < rise_inistep_pix) ? (rise_step_pix) : (rise_inistep_pix);
	nr_rise_samplings = CEIL(fabs(rise_local_min_pix - rise_local_max_pix) / rise_step_pix);
	nr_rise_samplings = (nr_rise_samplings > nr_min_samplings) ? (nr_rise_samplings) : (nr_min_samplings);
	nr_rise_samplings = (nr_rise_samplings < nr_max_samplings) ? (nr_rise_samplings) : (nr_max_samplings);
	rise_step_pix = fabs(rise_local_min_pix - rise_local_max_pix) / RFLOAT(nr_rise_samplings);
	if ((fabs(rise_local_min_pix - rise_local_max_pix) / fabs(rise_refined_pix)) < err_max)
	{
		search_rise = false;
		rise_step_pix = 0.;
		rise_min_pix = rise_max_pix = rise_local_min_pix = rise_local_max_pix = rise_refined_pix;
	}
	if (o_ptr != NULL)
	{
		if (search_rise)
			(*o_ptr) << " --> Initial searching step of rise is " << rise_step_pix * pixel_size_A << " Angstroms (" << nr_rise_samplings << " samplings)." << std::endl;
		else
			(*o_ptr) << " --> No need to search for rise..." << std::endl;
		(*o_ptr) << " --> " << nr_twist_samplings * nr_rise_samplings << " initial samplings." << std::endl;
	}

	if ( (!search_twist) && (!search_rise) )
		return true;

	if (o_ptr != NULL)
		(*o_ptr) << std::endl << " TAG   TWIST(DEGREES)  RISE(ANGSTROMS)         DEV" << std::endl;

	// Local searches
	helical_symmetry_list.clear();
	iter_not_converged = 0;
	for (iter = 1; iter <= 100; iter++)
	{
		// TODO: please check this!!!
		// rise_step_pix and twist_step_deg should be strictly > 0 now! (if they are to be searched)
		makeHelicalSymmetryList(
				helical_symmetry_list,
				rise_local_min_pix,
				rise_local_max_pix,
				rise_step_pix,
				search_rise,
				twist_local_min_deg,
				twist_local_max_deg,
				twist_step_deg,
				search_twist);
		if (helical_symmetry_list.size() < 1)
			REPORT_ERROR("helix.cpp::localSearchHelicalSymmetry(): BUG No helical symmetries are found in the search list!");

		best_dev = (1e30);
		best_id = -1;
		for (int ii = 0; ii < helical_symmetry_list.size(); ii++)
		{
			// If this symmetry is not calculated before
			if (helical_symmetry_list[ii].dev > (1e30))
			{
				// TODO: please check this!!!
				calcCCofHelicalSymmetry(
						v,
						r_min_pix,
						r_max_pix,
						z_percentage,
						helical_symmetry_list[ii].rise_pix,
						helical_symmetry_list[ii].twist_deg,
						helical_symmetry_list[ii].dev,
						nr_asym_voxels);
				if (o_ptr != NULL)
					(*o_ptr) << " NEW" << std::flush;
			}
			else
			{
				if (o_ptr != NULL)
					(*o_ptr) << " OLD" << std::flush;
			}

			if (helical_symmetry_list[ii].dev < best_dev)
			{
				best_dev = helical_symmetry_list[ii].dev;
				best_id = ii;
			}
			if (o_ptr != NULL)
			{
				(*o_ptr) << std::setw(15) << std::setiosflags(std::ios::fixed) << helical_symmetry_list[ii].twist_deg << std::resetiosflags(std::ios::fixed)
						<< std::setw(15) << std::setiosflags(std::ios::fixed) << (helical_symmetry_list[ii].rise_pix * pixel_size_A) << std::resetiosflags(std::ios::fixed)
						<< std::setw(20) << std::setiosflags(std::ios::scientific) << helical_symmetry_list[ii].dev << std::resetiosflags(std::ios::scientific) << std::endl;
			}
		}

		// Update refined symmetry
		rise_refined_pix = helical_symmetry_list[best_id].rise_pix;
		rise_refined_A = rise_refined_pix * pixel_size_A;
		twist_refined_deg = helical_symmetry_list[best_id].twist_deg;
		if (o_ptr != NULL)
		{
			(*o_ptr) << " ################################################################################" << std::endl;
			(*o_ptr) << " ##### Refined Twist = " << twist_refined_deg << ", Rise = " << rise_refined_A << ", Dev = " << helical_symmetry_list[best_id].dev << std::endl;
			(*o_ptr) << " ################################################################################" << std::endl;
		}
		// Out of range...
		if ( (search_rise) && (rise_refined_pix < rise_min_pix) )
		{
			out_of_range = true;
			search_rise = false;
			rise_step_pix = 0.;
			rise_local_min_pix = rise_local_max_pix = rise_refined_pix = rise_min_pix;
			rise_refined_A = rise_refined_pix * pixel_size_A;
		}
		if ( (search_rise) && (rise_refined_pix > rise_max_pix) )
		{
			out_of_range = true;
			search_rise = false;
			rise_step_pix = 0.;
			rise_local_min_pix = rise_local_max_pix = rise_refined_pix = rise_max_pix;
			rise_refined_A = rise_refined_pix * pixel_size_A;
		}
		if ( (search_twist) && (twist_refined_deg < twist_min_deg) )
		{
			out_of_range = true;
			search_twist = false;
			twist_step_deg = 0.;
			twist_local_min_deg = twist_local_max_deg = twist_refined_deg = twist_min_deg;
		}
		if ( (search_twist) && (twist_refined_deg > twist_max_deg) )
		{
			out_of_range = true;
			search_twist = false;
			twist_step_deg = 0.;
			twist_local_min_deg = twist_local_max_deg = twist_refined_deg = twist_max_deg;
		}

		// Not converged in this iteration...
		if ( (iter > 1) && (!out_of_range) )
		{
			// If the symmetry does not fall into the local search range
			// Try 7*, 9*, 11* ... samplings
			bool this_iter_not_converged = false;
			if (search_rise)
			{
				if ( (rise_refined_pix < (rise_local_min_pix + rise_step_pix * 0.5))
						|| (rise_refined_pix > (rise_local_max_pix - rise_step_pix * 0.5)) )
				{
					this_iter_not_converged = true;
					rise_local_min_pix = rise_refined_pix - ((RFLOAT)(iter_not_converged) + 3.) * rise_step_pix;
					rise_local_max_pix = rise_refined_pix + ((RFLOAT)(iter_not_converged) + 3.) * rise_step_pix;
				}
			}
			if (search_twist)
			{
				if ( (twist_refined_deg < (twist_local_min_deg + twist_step_deg * 0.5))
						|| (twist_refined_deg > (twist_local_max_deg - twist_step_deg * 0.5)) )
				{
					this_iter_not_converged = true;
					twist_local_min_deg = twist_refined_deg - ((RFLOAT)(iter_not_converged) + 3.) * twist_step_deg;
					twist_local_max_deg = twist_refined_deg + ((RFLOAT)(iter_not_converged) + 3.) * twist_step_deg;
				}
			}
			if (this_iter_not_converged)
			{
				iter_not_converged++;
				if (o_ptr != NULL)
					(*o_ptr) << " !!! NR_ITERATION_NOT_CONVERGED = " << iter_not_converged << " !!!" << std::endl;
				if (iter_not_converged > 10) // Up to 25*25 samplings are allowed (original 5*5 samplings)
				{
					if (o_ptr != NULL)
						(*o_ptr) << " WARNING: Local searches of helical symmetry cannot converge. Consider a finer initial sampling of helical parameters." << std::endl;
					else
						std::cout << " WARNING: Local searches of helical symmetry cannot converge. Consider a finer initial sampling of helical parameters." << std::endl;
					return false;
				}
				continue;
			}
		}
		iter_not_converged = 0;

		// Set 5*5 finer samplings for the next iteration
		if (search_rise)
		{
			rise_local_min_pix = rise_refined_pix - rise_step_pix;
			rise_local_max_pix = rise_refined_pix + rise_step_pix;
		}
		if (search_twist)
		{
			twist_local_min_deg = twist_refined_deg - twist_step_deg;
			twist_local_max_deg = twist_refined_deg + twist_step_deg;
		}

		// When there is no need to search for either twist or rise
		if ( (search_rise) && ((rise_step_pix / fabs(rise_refined_pix)) < err_max) )
		{
			rise_local_min_pix = rise_local_max_pix = rise_refined_pix;
			search_rise = false;
		}
		if ( (search_twist) && ((twist_step_deg / fabs(twist_refined_deg)) < err_max) )
		{
			twist_local_min_deg = twist_local_max_deg = twist_refined_deg;
			search_twist = false;
		}

		// Stop searches if step sizes are too small
		if ( (!search_twist) && (!search_rise) )
			break;

		// Decrease step size
		if (search_rise)
			rise_step_pix /= 2.;
		if (search_twist)
			twist_step_deg /= 2.;
	}

	if (out_of_range)
	{
		if (o_ptr != NULL)
			(*o_ptr) << " WARNING: Refined helical symmetry is out of the search range. Check whether the initial guess of helical symmetry is reasonable. Or you may want to modify the search range." << std::endl;
		else
			std::cout << " WARNING: Refined helical symmetry is out of the search range. Check whether the initial guess of helical symmetry is reasonable. Or you may want to modify the search range." << std::endl;
		return false;
	}
	return true;
};

RFLOAT getHelicalSigma2Rot(
		RFLOAT helical_rise_Angst,
		RFLOAT helical_twist_deg,
		RFLOAT helical_offset_step_Angst,
		RFLOAT rot_step_deg,
		RFLOAT old_sigma2_rot)
{
	if ( (helical_offset_step_Angst < 0.) || (rot_step_deg < 0.) || (old_sigma2_rot < 0.) )
		REPORT_ERROR("helix.cpp::getHelicalSigma2Rot: Helical offset step, rot step or sigma2_rot cannot be negative!");

	RFLOAT nr_samplings_along_helical_axis = (fabs(helical_rise_Angst)) / helical_offset_step_Angst;
	RFLOAT rot_search_range = (fabs(helical_twist_deg)) / nr_samplings_along_helical_axis;
	RFLOAT new_rot_step = rot_search_range / 6.;
	RFLOAT factor = ceil(new_rot_step / rot_step_deg);
	RFLOAT new_sigma2_rot = old_sigma2_rot;
	//RFLOAT factor_max = 10.;
	if (factor > 1.)
	{
		// Avoid extremely big sigma_rot!!! Too expensive in time and memory!!!
		//if (factor > factor_max)
		//	factor = factor_max;
		new_sigma2_rot *= factor * factor;
	}
	return new_sigma2_rot;
};

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
		bool verboseOutput)
{
	RFLOAT nr_units_min = 2.; // Minimum nr_particles required along lenZ_max
	RFLOAT rise_range_max_percentage = 0.3334;

	// Verbose output
	if (verboseOutput)
	{
		std::cout << "##########################################################" << std::endl;
		std::cout << "   CHECKING PARAMETERS FOR 3D HELICAL RECONSTRUCTION..." << std::endl;
		std::cout << "##########################################################" << std::endl;
	}

	// Check pixel size
	if (verboseOutput)
		std::cout << " Pixel size = " << pixel_size_A << " Angstrom(s)" << std::endl;
	if (pixel_size_A < 0.001)
	{
		if (verboseOutput)
			std::cout << " ERROR! Pixel size should be larger than 0.001 Angstroms!" << std::endl;
		else
			REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Pixel size should be larger than 0.001 Angstroms!");
		return false;
	}

	// Check box size and calculate half box size
	if (verboseOutput)
		std::cout << " Box size = " << box_len << " pixels = " << (RFLOAT)(box_len) * pixel_size_A << " Angstroms" << std::endl;
	if (box_len < 10)
	{
		if (verboseOutput)
			std::cout << " ERROR! Input box size should be larger than 10!" << std::endl;
		else
			REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Input box size should be larger than 10!");
		return false;
	}
	int half_box_len = box_len / 2 - ((box_len + 1) % 2);

	// Calculate radii in pixels
	RFLOAT particle_radius_pix = particle_diameter_A * 0.5 / pixel_size_A;
	RFLOAT tube_inner_radius_pix = tube_inner_diameter_A * 0.5 / pixel_size_A;
	RFLOAT tube_outer_radius_pix = tube_outer_diameter_A * 0.5 / pixel_size_A;

	// Check particle radius
	if (verboseOutput)
	{
		std::cout << " Particle diameter = " << particle_radius_pix * 2. << " pixels = " << particle_diameter_A << " Angstroms" << std::endl;
		std::cout << " Half box size = " << half_box_len << " pixels = " << (RFLOAT)(half_box_len) * pixel_size_A << " Angstroms" << std::endl;
	}
	if ( (particle_radius_pix < 2.) || (particle_radius_pix > half_box_len) )
	{
		if (verboseOutput)
			std::cout << " ERROR! Particle radius should be > 2 pixels and < half the box size!" << std::endl;
		else
			REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Particle radius should be > 2 and < half the box size!");
		return false;
	}

	// Check inner and outer tube radii
	if (verboseOutput)
	{
		if (tube_inner_diameter_A > 0.)
			std::cout << " Inner tube diameter = " << tube_inner_radius_pix * 2. << " pixels = " << tube_inner_diameter_A << " Angstroms" << std::endl;
		std::cout << " Outer tube diameter = " << tube_outer_radius_pix * 2. << " pixels = " << tube_outer_diameter_A << " Angstroms" << std::endl;
	}
	if ( (tube_outer_radius_pix < 2.) || (tube_outer_radius_pix > half_box_len)
			//|| ( (particle_radius_pix + 0.001) < tube_outer_radius_pix ) )
			|| (particle_radius_pix < tube_outer_radius_pix) )
	{
		if (verboseOutput)
			std::cout << " ERROR! Outer tube diameter should be > 4 pixels, < particle diameter and < half the box size!" << std::endl;
		else
			REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Outer tube diameter should be > 4 pixels, < particle diameter and < half the box size");
		return false;
	}
	if ( (tube_inner_radius_pix > 0.) && ((tube_inner_radius_pix + 2.) > tube_outer_radius_pix) )
	{
		if (verboseOutput)
			std::cout << " ERROR! Inner tube diameter should be remarkably smaller (> 4 pixels) than the outer one!" << std::endl;
		else
			REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Inner tube diameter should be remarkably smaller (> 4 pixels) than the outer one!");
		return false;
	}

	// STOP CHECKING OTHER PARAMETERS IF HELICAL SYMMETRY IS IGNORED IN 3D RECONSTRUCTION!
	if (ignore_symmetry)
	{
		if (verboseOutput)
			std::cout << " You have chosen to ignore helical symmetry! Stop checking now..." << std::endl;
		return true;
	}

	// CHECKING HELICAL SYMMETRY RELATED PARAMETERS...

	// Force same helical twist and rise if local refinement is not to be performed
	if (!do_symmetry_local_refinement)
	{
		twist_min_deg = twist_max_deg = twist_initial_deg;
		rise_min_A = rise_max_A = rise_initial_A;
	}
	RFLOAT rise_avg_A = (rise_min_A + rise_max_A) / 2.;
	RFLOAT twist_avg_deg = (twist_min_deg + twist_max_deg) / 2.;

	// Check helical twist and rise
	if (verboseOutput)
	{
		if (do_symmetry_local_refinement)
		{
			std::cout << " Helical twist (min, initial, average, max) = " << twist_min_deg << ", " << twist_initial_deg << ", " << twist_avg_deg << ", " << twist_max_deg << " degrees" << std::endl;
			std::cout << " Helical rise  (min, initial, average, max) = " << rise_min_A / pixel_size_A << ", " << rise_initial_A / pixel_size_A << ", " << rise_avg_A / pixel_size_A << ", " << rise_max_A / pixel_size_A << " pixels" << std::endl;
			std::cout << " Helical rise  (min, initial, average, max) = " << rise_min_A << ", " << rise_initial_A << ", " << rise_avg_A << ", " << rise_max_A << " Angstroms" << std::endl;
		}
		else
		{
			std::cout << " Helical twist = " << twist_initial_deg << " degree(s)" << std::endl;
			std::cout << " Helical rise  = " << rise_initial_A / pixel_size_A << " pixel(s) = " << rise_initial_A << " Angstrom(s)" << std::endl;
		}
	}
	if ( (fabs(twist_min_deg) > 360.) || (fabs(twist_initial_deg) > 360.) || (fabs(twist_max_deg) > 360.) )
	{
		if (verboseOutput)
			std::cout << " ERROR! Helical twist should be > -360 and < +360 degrees!" << std::endl;
		else
			REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Helical twist should be > -360 and < +360 degrees!");
		return false;
	}
	if ( (rise_min_A < 0.001) || ((rise_min_A / pixel_size_A) < 0.001)
			|| (rise_initial_A < 0.001) || ((rise_initial_A / pixel_size_A) < 0.001)
			|| (rise_max_A < 0.001) || ((rise_max_A / pixel_size_A) < 0.001) )
	{
		if (verboseOutput)
			std::cout << " ERROR! Helical rise should be > +0.001 Angstroms and > +0.001 pixels!" << std::endl;
		else
			REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Helical rise should be > +0.001 Angstroms and > +0.001 pixels!");
		return false;
	}
	if (do_symmetry_local_refinement)
	{
		if ( (twist_min_deg > twist_max_deg) || (twist_initial_deg < twist_min_deg) || (twist_initial_deg > twist_max_deg)
				|| (rise_min_A > rise_max_A) || (rise_initial_A < rise_min_A) || (rise_initial_A > rise_max_A) )
		{
			if (verboseOutput)
				std::cout << " ERROR! The following condition must be satisfied (both for helical twist and rise): min < initial < max !" << std::endl;
			else
				REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): The following condition must be satisfied (both for helical twist and rise): min < initial < max !");
			return false;
		}
#ifndef WIDE_HELICAL_TWIST_AND_RISE_SEARCHES
		//RFLOAT rise_min_A_thres = rise_avg_A * (1. - rise_range_max_percentage);
		//RFLOAT rise_max_A_thres = rise_avg_A * (1. + rise_range_max_percentage);
		if ( (fabs(twist_avg_deg - twist_min_deg) > 180.01) || ((fabs(rise_avg_A - rise_min_A) / fabs(rise_avg_A)) > rise_range_max_percentage) )
		{
			if (verboseOutput)
				std::cout << " ERROR! Searching ranges of helical twist and rise should be < 180 degrees and < +/-33.34% from the min-max average respectively!" << std::endl;
			else
				REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Searching ranges of helical twist and rise should be < 180 degrees and < +/-33.34% from the min-max average respectively!");
			return false;
		}
#endif
	}

	// Check Z percentage
	if (verboseOutput)
		std::cout << " Z percentage = " << z_percentage << " = " << z_percentage * 100. << " %" << std::endl;
	if ( (z_percentage < 0.001) || (z_percentage > 0.999) )
	{
		if (verboseOutput)
			std::cout << " ERROR! Z percentage should at least be > 0.001 and < 0.999 (0.1% ~ 99.9%)!" << std::endl;
		else
			REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Z percentage should at least be > 0.001 and < 0.999 (0.1% ~ 99.9%)!");
		return false;
	}
	RFLOAT z_percentage_min = (nr_units_min * rise_max_A) / (pixel_size_A * (RFLOAT)(box_len));
	z_percentage_min = (z_percentage_min < 0.001) ? (0.001) : (z_percentage_min);
	RFLOAT z_percentage_max = ( (2.) * sqrt( (particle_diameter_A * particle_diameter_A / 4.) - (tube_outer_diameter_A * tube_outer_diameter_A / 4.) ) / pixel_size_A) / ((RFLOAT)(box_len));
	z_percentage_max = (z_percentage_max > 0.999) ? (0.999) : (z_percentage_max);
	if (verboseOutput)
		std::cout << " Z percentage should be > " << z_percentage_min << " and < " << z_percentage_max << " (under current settings)" << std::endl;
	if (z_percentage_min > z_percentage_max)
	{
		if (verboseOutput)
			std::cout << " ERROR! The range of Z percentage is invalid! To decrease the lower bound, make maximum rise smaller or box size larger. To increase the upper bound, make the particle diameter (along with the box size) larger and the outer tube diameter smaller!" << std::endl;
		else
			REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): The range of Z percentage is invalid! To decrease the lower bound, make maximum rise smaller or box size larger. To increase the upper bound, make the particle diameter (along with the box size) larger and the outer tube diameter smaller!");
		return false;
	}
	if ( (z_percentage < z_percentage_min) || (z_percentage > z_percentage_max) )
	{
		if (verboseOutput)
			std::cout << " ERROR! Z percentage is out of range under current settings!" << std::endl;
		else
			REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Z percentage is out of range under current settings!");
		return false;
	}

	// Check maximum rise (DO I NEED THIS???)
	RFLOAT rise_max_upper_bound_A = pixel_size_A * (RFLOAT)(box_len) * z_percentage / nr_units_min;
	if (do_symmetry_local_refinement)
	{
		if (verboseOutput)
			std::cout << " Upper bound of maximum rise = " << rise_max_upper_bound_A / pixel_size_A << " pixels = " << rise_max_upper_bound_A << " Angstroms (under current settings)" << std::endl;
		if (fabs(rise_max_A) > rise_max_upper_bound_A) // THIS CANNOT HAPPEN. ERRORS HAVE ALREADY BEEN RAISED IN Z PERCENTAGE CHECK.
		{
			if (verboseOutput)
				std::cout << " ERROR! Maximum rise exceeds its upper bound!" << std::endl;
			else
				REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Maximum rise exceeds its upper bound!");
			return false;
		}
	}

	// Check number of asymmetrical units
	RFLOAT half_nr_asu_max = 0.5 * (1. - z_percentage) * (RFLOAT(box_len)) * pixel_size_A / rise_max_A;
	int nr_asu_max = 2 * (int(floor(half_nr_asu_max))) + 1;
	nr_asu_max = (nr_asu_max < 1) ? (1) : (nr_asu_max);
	if (verboseOutput)
		std::cout << " Number of asymmetrical units = " << nr_asu << ", maximum value = " << nr_asu_max << " (under current settings)" << std::endl;
	if (nr_asu < 1)
	{
		if (verboseOutput)
			std::cout << " ERROR! Number of asymmetrical units (an integer) should be at least 1!" << std::endl;
		else
			REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Number of asymmetrical units (a positive integer) should be at least 1!");
		return false;
	}
	if ( (nr_asu > 1) && (nr_asu > nr_asu_max) )
	{
		if (verboseOutput)
			std::cout << " ERROR! Number of asymmetrical units exceeds its upper bound!" << std::endl;
		else
			REPORT_ERROR("helix.cpp::chechParametersFor3DHelicalReconstruction(): Number of asymmetrical units exceeds its upper bound!");
		return false;
	}

	// Everything seems fine :)
	return true;
}

void imposeHelicalSymmetryInRealSpace(
		MultidimArray<RFLOAT>& v,
		RFLOAT pixel_size_A,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_inner_radius_A,
		RFLOAT cyl_outer_radius_A,
		RFLOAT z_percentage,
		RFLOAT rise_A,
		RFLOAT twist_deg,
		RFLOAT cosine_width_pix)
{
	bool ignore_helical_symmetry = false;
	long int Xdim, Ydim, Zdim, Ndim, box_len;
	RFLOAT rise_pix, sphere_radius_pix, cyl_inner_radius_pix, cyl_outer_radius_pix, r_min, r_max, d_min, d_max, D_min, D_max, z_min, z_max;

	int rec_len;
	std::vector<RFLOAT> sin_rec, cos_rec;
	MultidimArray<RFLOAT> vout;

	if (v.getDim() != 3)
		REPORT_ERROR("helix.cpp::imposeHelicalSymmetryInRealSpace(): Input helical reference is not 3D! (vol.getDim() = " + integerToString(v.getDim()) + ")");
	v.getDimensions(Xdim, Ydim, Zdim, Ndim);
	box_len = (Xdim < Ydim) ? Xdim : Ydim;
	box_len = (box_len < Zdim) ? box_len : Zdim;

	// Check helical parameters
	checkParametersFor3DHelicalReconstruction(
			false,
			false,
			1,
			rise_A,
			rise_A,
			rise_A,
			twist_deg,
			twist_deg,
			twist_deg,
			box_len,
			pixel_size_A,
			z_percentage,
			sphere_radius_A * 2.,
			cyl_inner_radius_A * 2.,
			cyl_outer_radius_A * 2.);

	// Parameters of mask
	if (cosine_width_pix < 1.)
		cosine_width_pix = 1.;  // Avoid 'divided by 0' error
	rise_pix = fabs(rise_A / pixel_size_A);  // Keep helical rise as a positive number
	sphere_radius_pix = sphere_radius_A / pixel_size_A;
	cyl_inner_radius_pix = cyl_inner_radius_A / pixel_size_A;
	cyl_outer_radius_pix = cyl_outer_radius_A / pixel_size_A;
	r_min = sphere_radius_pix;
	r_max = sphere_radius_pix + cosine_width_pix;
	d_min = cyl_inner_radius_pix - cosine_width_pix;
	d_max = cyl_inner_radius_pix;
	D_min = cyl_outer_radius_pix;
	D_max = cyl_outer_radius_pix + cosine_width_pix;

	// Crop the central slices
	v.setXmippOrigin();
	z_max = ((RFLOAT)(Zdim)) * z_percentage / 2.;
	if (z_max > (((RFLOAT)(FINISHINGZ(v))) - 1.))
		z_max = (((RFLOAT)(FINISHINGZ(v))) - 1.);
	z_min = -z_max;
	if (z_min < (((RFLOAT)(STARTINGZ(v))) + 1.))
		z_min = (((RFLOAT)(STARTINGZ(v))) + 1.);

	// Init volumes
	v.setXmippOrigin();
	vout.clear();
	vout.resize(v);
	vout.setXmippOrigin();

	// Calculate tabulated sine and cosine values
	rec_len = 2 + (CEIL((RFLOAT(Zdim) + 2.) / rise_pix));
	sin_rec.clear();
	cos_rec.clear();
	sin_rec.resize(rec_len);
	cos_rec.resize(rec_len);
	for (int id = 0; id < rec_len; id++)
#ifdef RELION_SINGLE_PRECISION
		SINCOSF(DEG2RAD(((RFLOAT)(id)) * twist_deg), &sin_rec[id], &cos_rec[id]);
#else
		SINCOS(DEG2RAD(((RFLOAT)(id)) * twist_deg), &sin_rec[id], &cos_rec[id]);
#endif

	FOR_ALL_ELEMENTS_IN_ARRAY3D(v)
	{
		// Out of the mask
		RFLOAT dd = (RFLOAT)(i * i + j * j);
		RFLOAT rr = dd + (RFLOAT)(k * k);
		RFLOAT d = sqrt(dd);
		RFLOAT r = sqrt(rr);
		if ( (r > r_max) || (d < d_min) || (d > D_max) )
		{
			A3D_ELEM(v, k, i, j) = 0.;
			continue;
		}

		// How many voxels should be used to calculate the average?
		RFLOAT zi = (RFLOAT)(k);
		RFLOAT yi = (RFLOAT)(i);
		RFLOAT xi = (RFLOAT)(j);
		int rot_max = -(CEIL((zi - z_max) / rise_pix));
		int rot_min = -(FLOOR((zi - z_min) / rise_pix));
		if (rot_max < rot_min)
			REPORT_ERROR("helix.cpp::makeHelicalReferenceInRealSpace(): ERROR in imposing symmetry!");

		// Do the average
		RFLOAT pix_sum, pix_weight;
		pix_sum = pix_weight = 0.;
		for (int id = rot_min; id <= rot_max; id++)
		{
			// Get the sine and cosine value
			RFLOAT sin_val, cos_val;
			if (id >= 0)
			{
				sin_val = sin_rec[id];
				cos_val = cos_rec[id];
			}
			else
			{
				sin_val = (-1.) * sin_rec[-id];
				cos_val = cos_rec[-id];
			}

			// Get the voxel coordinates
			RFLOAT zp = zi + ((RFLOAT)(id)) * rise_pix;
			RFLOAT yp = xi * sin_val + yi * cos_val;
			RFLOAT xp = xi * cos_val - yi * sin_val;

			// Trilinear interpolation (with physical coords)
			// Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
			// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
			int x0, y0, z0, x1, y1, z1;
			RFLOAT fx, fy, fz;
			x0 = FLOOR(xp); fx = xp - x0; x0 -= STARTINGX(v); x1 = x0 + 1;
			y0 = FLOOR(yp); fy = yp - y0; y0 -= STARTINGY(v); y1 = y0 + 1;
			z0 = FLOOR(zp); fz = zp - z0; z0 -= STARTINGZ(v); z1 = z0 + 1;

			RFLOAT d000, d001, d010, d011, d100, d101, d110, d111;
			d000 = DIRECT_A3D_ELEM(v, z0, y0, x0);
			d001 = DIRECT_A3D_ELEM(v, z0, y0, x1);
			d010 = DIRECT_A3D_ELEM(v, z0, y1, x0);
			d011 = DIRECT_A3D_ELEM(v, z0, y1, x1);
			d100 = DIRECT_A3D_ELEM(v, z1, y0, x0);
			d101 = DIRECT_A3D_ELEM(v, z1, y0, x1);
			d110 = DIRECT_A3D_ELEM(v, z1, y1, x0);
			d111 = DIRECT_A3D_ELEM(v, z1, y1, x1);

			RFLOAT dx00, dx01, dx10, dx11;
			dx00 = LIN_INTERP(fx, d000, d001);
			dx01 = LIN_INTERP(fx, d100, d101);
			dx10 = LIN_INTERP(fx, d010, d011);
			dx11 = LIN_INTERP(fx, d110, d111);

			RFLOAT dxy0, dxy1;
			dxy0 = LIN_INTERP(fy, dx00, dx10);
			dxy1 = LIN_INTERP(fy, dx01, dx11);

			pix_sum += LIN_INTERP(fz, dxy0, dxy1);
			pix_weight += 1.;
		}

		if (pix_weight > 0.9)
		{
			A3D_ELEM(vout, k, i, j) = pix_sum / pix_weight;

			if ( (d > d_max) && (d < D_min) && (r < r_min) )
			{}
			else // The pixel is within cosine edge(s)
			{
				pix_weight = 1.;
				if (d < d_max)  // d_min < d < d_max : w=(0~1)
					pix_weight = 0.5 + (0.5 * cos(PI * ((d_max - d) / cosine_width_pix)));
				else if (d > D_min) // D_min < d < D_max : w=(1~0)
					pix_weight = 0.5 + (0.5 * cos(PI * ((d - D_min) / cosine_width_pix)));
				if (r > r_min) // r_min < r < r_max
				{
					pix_sum = 0.5 + (0.5 * cos(PI * ((r - r_min) / cosine_width_pix)));
					pix_weight = (pix_sum < pix_weight) ? (pix_sum) : (pix_weight);
				}
				A3D_ELEM(vout, k, i, j) *= pix_weight;
			}
		}
		else
			A3D_ELEM(vout, k, i, j) = 0.;
    }

	// Copy and exit
	v = vout;
	sin_rec.clear();
	cos_rec.clear();
	vout.clear();

	return;
};

/*
void searchCnZSymmetry(
		const MultidimArray<RFLOAT>& v,
		RFLOAT r_min_pix,
		RFLOAT r_max_pix,
		int cn_start,
		int cn_end,
		std::vector<RFLOAT>& cn_list,
		std::vector<RFLOAT>& cc_list,
		std::vector<int>& nr_asym_voxels_list,
		std::ofstream* fout_ptr)
{
	int cn, Xdim, Ydim, Zdim, Ndim, nr_asym_voxels;
	RFLOAT cc;
	bool ok_flag;

	Xdim = XSIZE(v); Ydim = YSIZE(v); Zdim = ZSIZE(v); Ndim = NSIZE(v);

	if( (Ndim != 1) || (Zdim < 5) || (Ydim < 5) || (Xdim < 5) )
		REPORT_ERROR("helix.cpp::searchCnZSymmetry(): Input 3D MultidimArray has Wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");

	if( (r_max_pix < 2.) || (r_min_pix < 0.) || ((r_max_pix - r_min_pix) < 2.)
			|| (cn_start <= 1) || (cn_end <= 1) || (cn_start > cn_end) || (cn_end > 36) )
		REPORT_ERROR("helix.cpp::searchCnZSymmetry(): Wrong parameters!");

	cn_list.clear();
	cc_list.clear();
	nr_asym_voxels_list.clear();

	for(cn = cn_start; cn <= cn_end; cn++)
	{
		ok_flag = calcCCOfCnZSymmetry(v, r_min_pix, r_max_pix, cn, cc, nr_asym_voxels);
		if(!ok_flag)
			continue;

		cn_list.push_back(cn);
		cc_list.push_back(cc);
		nr_asym_voxels_list.push_back(nr_asym_voxels);

		if(fout_ptr != NULL)
			(*fout_ptr) << "Test Cn = " << cn << ", cc = " << cc << ",                               asym voxels = " << nr_asym_voxels << std::endl;
	}
	return;
}
*/

/*
RFLOAT calcCCofPsiFor2DHelicalSegment(
		const MultidimArray<RFLOAT>& v,
		RFLOAT psi_deg,
		RFLOAT pixel_size_A,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_outer_radius_A)
{
	RFLOAT sphere_radius_pix, sphere_radius2_pix, cyl_radius_pix, r2, x, y, xp, yp, sum, sum2, nr, val;
	int x0, y0, vec_id, vec_len, half_box_len, box_len;
	std::vector<RFLOAT> sum_list, pix_list;
	Matrix2D<RFLOAT> R;

	if (pixel_size_A < 0.001)
		REPORT_ERROR("helix.cpp::calcCCofPsiFor2DHelicalSegment(): Pixel size (in Angstroms) should be larger than 0.001!");

	if (v.getDim() != 2)
		REPORT_ERROR("helix.cpp::calcCCofPsiFor2DHelicalSegment(): Input MultidimArray should be 2D!");

	if ( (YSIZE(v) < 10) || (XSIZE(v) < 10) )
		REPORT_ERROR("helix.cpp::calcCCofPsiFor2DHelicalSegment(): Input 2D MultidimArray should be larger than 10*10 pixels!");

	if ( (STARTINGY(v) != FIRST_XMIPP_INDEX(YSIZE(v))) || (STARTINGX(v) != FIRST_XMIPP_INDEX(XSIZE(v))) )
		REPORT_ERROR("helix.cpp::calcCCofPsiFor2DHelicalSegment(): The origin of input 2D MultidimArray is not at the center (use v.setXmippOrigin() before calling this function)!");

	box_len = (XSIZE(v) < YSIZE(v)) ? XSIZE(v) : YSIZE(v);
	half_box_len = box_len / 2 - ((box_len + 1) % 2);
	sphere_radius_pix = sphere_radius_A / pixel_size_A;
	cyl_radius_pix = cyl_outer_radius_A / pixel_size_A;
	if ( (sphere_radius_pix < 2.) || (sphere_radius_pix > half_box_len)
			|| (cyl_radius_pix < 2.) || (cyl_radius_pix > half_box_len)
			|| ( (sphere_radius_pix + 0.001) < cyl_radius_pix ) )
		REPORT_ERROR("helix.cpp::calcCCofPsiFor2DHelicalSegment(): Radii of spherical and/or cylindrical masks are invalid!");

	sphere_radius2_pix = sphere_radius_pix * sphere_radius_pix;
	x0 = STARTINGX(v);
	y0 = STARTINGY(v);
	vec_len = (int)((2. * cyl_radius_pix)) + 2;
	sum_list.clear();
	sum_list.resize(vec_len);
	pix_list.clear();
	pix_list.resize(vec_len);
	for (vec_id = 0; vec_id < vec_len; vec_id++)
		sum_list[vec_id] = pix_list[vec_id] = 0.;
	rotation2DMatrix(psi_deg, R, false);
	R.setSmallValuesToZero();

	FOR_ALL_ELEMENTS_IN_ARRAY2D(v)
	{
		x = j;
		y = i;
		r2 = i * i + j * j;
		if (r2 > sphere_radius2_pix)
			continue;
		//xp = x * R(0, 0) + y * R(0, 1);
		//vec_id = ROUND(xp - x0);
		yp = x * R(1, 0) + y * R(1, 1);
		vec_id = ROUND(yp - y0);
		if ( (vec_id < 0) || (vec_id >= vec_len) )
			continue;
		pix_list[vec_id]++;
		sum_list[vec_id] += A2D_ELEM(v, i, j);
	}

	nr = sum = sum2 = 0.;
	for (vec_id = 0; vec_id < vec_len; vec_id++)
	{
		if (pix_list[vec_id] > 0.5)
		{
			nr += 1.;
			val = sum_list[vec_id] / pix_list[vec_id];
			sum += val;
			sum2 += val * val;
		}
	}
	if (nr < 0.5)
		return (-1.);
	return ( (sum2 / nr) - ((sum / nr) * (sum / nr)) );
}

RFLOAT localSearchPsiFor2DHelicalSegment(
		const MultidimArray<RFLOAT>& v,
		RFLOAT pixel_size_A,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_outer_radius_A,
		RFLOAT ori_psi_deg,
		RFLOAT search_half_range_deg,
		RFLOAT search_step_deg)
{
	RFLOAT psi_deg, best_psi_deg, max_cc, cc;
	if ( (search_step_deg < 0.000001) || (search_step_deg > 2.) )
		REPORT_ERROR("helix.cpp::localSearchPsiFor2DHelicalSegment(): Search step of psi angle should be 0.000001~1.999999 degrees!");
	if ( (search_half_range_deg < 0.000001) || (search_half_range_deg < search_step_deg) )
		REPORT_ERROR("helix.cpp::localSearchPsiFor2DHelicalSegment(): Search half range of psi angle should be not be less than 0.000001 degrees or the step size!");

	max_cc = -1.;
	best_psi_deg = 0.;
	for (psi_deg = (ori_psi_deg - search_half_range_deg); psi_deg < (ori_psi_deg + search_half_range_deg); psi_deg += search_step_deg)
	{
		cc = calcCCofPsiFor2DHelicalSegment(v, psi_deg, pixel_size_A, sphere_radius_A, cyl_outer_radius_A);
		if ( (cc > 0.) && (cc > max_cc) )
		{
			max_cc = cc;
			best_psi_deg = psi_deg;
		}
		// DEBUG
		//std::cout << "psi_deg = " << psi_deg << ", cc = " << cc << std::endl;
	}
	// DEBUG
	//std::cout << "------------------------------------------------" << std::endl;
	if (max_cc < 0.)
		REPORT_ERROR("helix.cpp::localSearchPsiFor2DHelicalSegment(): Error! No cc values for any psi angles are valid! Check if the input images are blank!");
	return best_psi_deg;
}

RFLOAT searchPsiFor2DHelicalSegment(
		const MultidimArray<RFLOAT>& v,
		RFLOAT pixel_size_A,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_outer_radius_A)
{
	int nr_iter;
	RFLOAT search_half_range_deg, search_step_deg, best_psi_deg;

	search_half_range_deg = 89.999;
	search_step_deg = 0.5;
	best_psi_deg = 0.;

	for(nr_iter = 1; nr_iter <= 20; nr_iter++)
	{
		best_psi_deg = localSearchPsiFor2DHelicalSegment(v, pixel_size_A, sphere_radius_A, cyl_outer_radius_A,
				best_psi_deg, search_half_range_deg, search_step_deg);
		search_half_range_deg = 2. * search_step_deg;
		search_step_deg /= 2.;
		if (fabs(search_step_deg) < 0.00001)
			break;
	}
	return best_psi_deg;
}
*/

void calcRadialAverage(
		const MultidimArray<RFLOAT>& v,
		std::vector<RFLOAT>& radial_avg_val_list)
{
	int ii, Xdim, Ydim, Zdim, Ndim, list_size, dist;
	MultidimArray<RFLOAT> vol;
	std::vector<RFLOAT> radial_pix_counter_list;

	vol.clear();
	radial_pix_counter_list.clear();
	radial_avg_val_list.clear();

	// Check dimensions
	Xdim = XSIZE(v); Ydim = YSIZE(v); Zdim = ZSIZE(v); Ndim = NSIZE(v);

	if( (Ndim != 1) || (Zdim < 5) || (Ydim < 5) || (Xdim < 5) )
		REPORT_ERROR("helix.cpp::calcRadialAverage(): Input 3D MultidimArray has wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");

	// Resize and init vectors
	list_size = ROUND(sqrt(Xdim * Xdim + Ydim * Ydim)) + 2;
	radial_pix_counter_list.resize(list_size);
	radial_avg_val_list.resize(list_size);
	for(ii = 0; ii < list_size; ii++)
	{
		radial_pix_counter_list[ii] = 0.;
		radial_avg_val_list[ii] = 0.;
	}

	vol = v;
	vol.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
	{
		dist = ROUND(sqrt(i * i + j * j));
		if( (dist < 0) || (dist > (list_size - 1)) )
			continue;
		radial_pix_counter_list[dist] += 1.;
		radial_avg_val_list[dist] += A3D_ELEM(vol, k, i, j);
	}

	for(ii = 0; ii < list_size; ii++)
	{
		if(radial_pix_counter_list[ii] > 0.9)
			radial_avg_val_list[ii] /= radial_pix_counter_list[ii];
		else
			radial_avg_val_list[ii] = 0.;
	}

	radial_pix_counter_list.clear();
	vol.clear();

	return;
};

void cutZCentralPartOfSoftMask(
		MultidimArray<RFLOAT>& mask,
		RFLOAT z_percentage,
		RFLOAT cosine_width)
{
	int dim = mask.getDim();
	int Zdim = ZSIZE(mask);
	RFLOAT idz_s, idz_s_w, idz_e, idz_e_w, idz, val;
	if (dim != 3)
		REPORT_ERROR("helix.cpp::cutZCentralPartOfSoftMask(): Input mask should have a dimension of 3!");
	if (Zdim < 5)
		REPORT_ERROR("helix.cpp::cutZCentralPartOfSoftMask(): Z length of 3D mask is less than 5!");
	if ( (z_percentage < 0.1) || (z_percentage > 0.9) )
		REPORT_ERROR("helix.cpp::cutZCentralPartOfSoftMask(): Only 10%-90% of total Z length should be retained!");
	if (cosine_width < 0.001)
		REPORT_ERROR("helix.cpp::cutZCentralPartOfSoftMask(): Cosine width for soft edge should larger than 0!");

	idz_e = ((RFLOAT)(Zdim)) * z_percentage / 2.;
	idz_s = idz_e * (-1.);
	idz_s_w = idz_s - cosine_width;
	idz_e_w = idz_e + cosine_width;
	// DEBUG
	//std::cout << "z_len, z_percentage, cosine_width = " << Zdim << ", " << z_percentage << ", " << cosine_width << std::endl;
	//std::cout << "idz_s_w, idz_s, idz_e, idz_e_w = " << idz_s_w << ", " << idz_s << ", " << idz_e << ", " << idz_e_w << std::endl;
	mask.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
	{
		idz = ((RFLOAT)(k));
		if ( (idz > idz_s) && (idz < idz_e) )
		{}
		else if ( (idz < idz_s_w) || (idz > idz_e_w) )
			A3D_ELEM(mask, k, i, j) = 0.;
		else
		{
			val = 1.;
			if (idz < idz_s)
				val = 0.5 + 0.5 * cos(PI * (idz_s - idz) / cosine_width);
			else if (idz > idz_e)
				val = 0.5 + 0.5 * cos(PI * (idz - idz_e) / cosine_width);
			A3D_ELEM(mask, k, i, j) *= val;
		}
	}
	return;
};

void createCylindricalReference(
		MultidimArray<RFLOAT>& v,
		int box_size,
		RFLOAT inner_diameter_pix,
		RFLOAT outer_diameter_pix,
		RFLOAT cosine_width)
{
	RFLOAT r, dist, inner_radius_pix, outer_radius_pix;

	// Check dimensions
	if (box_size < 5)
		REPORT_ERROR("helix.cpp::createCylindricalReference(): Invalid box size.");

	if ( (inner_diameter_pix > outer_diameter_pix)
			|| (outer_diameter_pix < 0.) || (outer_diameter_pix > (box_size - 1))
			|| (cosine_width < 0.) )
		REPORT_ERROR("helix.cpp::createCylindricalReference(): Parameter(s) error!");

	inner_radius_pix = inner_diameter_pix / 2.;
	outer_radius_pix = outer_diameter_pix / 2.;

	v.clear();
	v.resize(box_size, box_size, box_size);
	v.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(v)
	{
		r = sqrt(i * i + j * j);
		if ( (r > inner_radius_pix) && (r < outer_radius_pix) )
		{
			A3D_ELEM(v, k, i, j) = 1.;
			continue;
		}
		dist = -9999.;
		if ( (r > outer_radius_pix) && (r < (outer_radius_pix + cosine_width)) )
			dist = r - outer_radius_pix;
		else if ( (r < inner_radius_pix) && (r > (inner_radius_pix - cosine_width)) )
			dist = inner_radius_pix - r;
		if (dist > 0.)
		{
			A3D_ELEM(v, k, i, j) = 0.5 + 0.5 * cos(PI * dist / cosine_width);
			continue;
		}
		A3D_ELEM(v, k, i, j) = 0.;
	}
	return;
}

void createCylindricalReferenceWithPolarity(
		MultidimArray<RFLOAT>& v,
		int box_size,
		RFLOAT inner_diameter_pix,
		RFLOAT outer_diameter_pix,
		RFLOAT ratio_topbottom,
		RFLOAT cosine_width)
{
	RFLOAT r, r_min, r_max, dist, top_radius_pix, bottom_radius_pix;

	// Check dimensions
	if (box_size < 5)
		REPORT_ERROR("helix.cpp::createCylindricalReferenceWithPolarity(): Invalid box size.");

	if ( (inner_diameter_pix > outer_diameter_pix)
			|| (outer_diameter_pix < 0.) || (outer_diameter_pix > (box_size - 1))
			|| (ratio_topbottom < 0.) || (ratio_topbottom > 1.) || (cosine_width < 0.) )
		REPORT_ERROR("helix.cpp::createCylindricalReferenceWithPolarity(): Parameter(s) error!");

	// Set top and bottom radii
	top_radius_pix = outer_diameter_pix / 2.;
	bottom_radius_pix = outer_diameter_pix * ratio_topbottom / 2.;
	if (inner_diameter_pix > 0.)
		bottom_radius_pix = (inner_diameter_pix / 2.) + ratio_topbottom * (outer_diameter_pix / 2. - inner_diameter_pix / 2.);

	v.clear();
	v.resize(box_size, box_size, box_size);
	v.setXmippOrigin();

	r_min = r_max = -1.;
	if (inner_diameter_pix > 0.)
		r_min = inner_diameter_pix / 2.;
    for (long int k=STARTINGZ(v); k<=FINISHINGZ(v); k++)
    {
    	r_max = top_radius_pix - (top_radius_pix - bottom_radius_pix) * ((RFLOAT)(k - STARTINGZ(v))) / ((RFLOAT)(box_size));
        for (long int i=STARTINGY(v); i<=FINISHINGY(v); i++)
        {
            for (long int j=STARTINGX(v); j<=FINISHINGX(v); j++)
            {
            	r = sqrt(i * i + j * j);
        		if ( (r > r_min) && (r < r_max) )
        		{
        			A3D_ELEM(v, k, i, j) = 1.;
        			continue;
        		}
        		dist = -9999.;
        		if ( (r > r_max) && (r < (r_max + cosine_width)) )
        			dist = r - r_max;
        		if ( (r < r_min) && (r > (r_min - cosine_width)) )
        			dist = r_min - r;
        		if (dist > 0.)
        		{
        			A3D_ELEM(v, k, i, j) = 0.5 + 0.5 * cos(PI * dist / cosine_width);
        			continue;
        		}
        		A3D_ELEM(v, k, i, j) = 0.;
            }
        }
    }
	return;
}

void transformCartesianAndHelicalCoords(
		Matrix1D<RFLOAT>& in,
		Matrix1D<RFLOAT>& out,
		RFLOAT rot_deg,
		RFLOAT tilt_deg,
		RFLOAT psi_deg,
		bool direction)
{
	int dim;
	RFLOAT x0, y0, z0;
	Matrix1D<RFLOAT> aux;
	Matrix2D<RFLOAT> A, B;

	dim = in.size();
	if( (dim != 2) && (dim != 3) )
		REPORT_ERROR("helix.cpp::transformCartesianAndHelicalCoords(): Vector of input coordinates should have 2 or 3 values!");

	aux.clear();
	aux.resize(3);
	XX(aux) = XX(in);
	YY(aux) = YY(in);
	ZZ(aux) = (dim == 3) ? (ZZ(in)) : (0.);

	if (dim == 2)
		rot_deg = tilt_deg = 0.;

	A.clear();
	A.resize(3, 3);
	// TODO: check whether rot_deg should be always set to 0 !
	// TODO: fix the --random_seed and use --perturb 0 option for testing !
	Euler_angles2matrix(rot_deg, tilt_deg, psi_deg, A, false);
	if (direction == CART_TO_HELICAL_COORDS) // Don't put minus signs before angles, use 'transpose' instead
		A = A.transpose();
	aux = A * aux;

	out.clear();
	out.resize(2);
	XX(out) = XX(aux);
	YY(out) = YY(aux);
	if (dim == 3)
	{
		out.resize(3);
		ZZ(out) = ZZ(aux);
	}
	aux.clear();
	A.clear();

	return;
}

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
		bool direction)
{
	if( (dim != 2) && (dim != 3) )
		REPORT_ERROR("helix.cpp::transformCartesianAndHelicalCoords(): Vector of input coordinates should have 2 or 3 values!");

	Matrix1D<RFLOAT> in, out;
	in.clear();
	out.clear();
	in.resize(dim);
	XX(in) = xin;
	YY(in) = yin;
	if (dim == 3)
		ZZ(in) = zin;
	transformCartesianAndHelicalCoords(in, out, rot_deg, tilt_deg, psi_deg, direction);
	xout = XX(out);
	yout = YY(out);
	if (dim == 3)
		zout = ZZ(out);
	return;
}

/*
void makeBlot(
		MultidimArray<RFLOAT>& v,
		RFLOAT y,
		RFLOAT x,
		RFLOAT r)
{
	int Xdim, Ydim, Zdim, Ndim;
	RFLOAT dist, min;
	v.getDimensions(Xdim, Ydim, Zdim, Ndim);
	if( (Ndim != 1) || (Zdim != 1) || (YXSIZE(v) <= 2) )
		return;

	min = DIRECT_A2D_ELEM(v, 0, 0);
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(v)
		min = (DIRECT_A2D_ELEM(v, i, j) < min) ? DIRECT_A2D_ELEM(v, i, j) : min;

	v.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY2D(v)
	{
		dist = (i - y) * (i - y) + (j - x) * (j - x);
		dist = sqrt(dist);
		if(dist < r)
			A2D_ELEM(v, i, j) = min;
	}
	return;
}
*/

void makeSimpleHelixFromPDBParticle(
		const Assembly& ori,
		Assembly& helix,
		RFLOAT radius_A,
		RFLOAT twist_deg,
		RFLOAT rise_A,
		int nr_copy,
		bool do_center)
{
	int nr_ori_atoms;
	Matrix1D<RFLOAT> mass_centre, shift;
	Matrix2D<RFLOAT> rotational_matrix;
	Assembly aux0, aux1;

	if (nr_copy < 3)
		REPORT_ERROR("helix.cpp::makeHelixFromPDBParticle(): More than 3 copies of original assemblies are required to form a helix!");

	nr_ori_atoms = ori.numberOfAtoms();
	std::cout << "Original assembly contains " << nr_ori_atoms << " atoms." << std::endl;

	if (nr_ori_atoms < 1)
		REPORT_ERROR("helix.cpp::makeHelixFromPDBParticle(): Original assembly contains no atoms!");

	// Calculate centre of mass of the original assembly
	mass_centre.resize(3);
	mass_centre.initZeros();
	for (int imol = 0; imol < ori.molecules.size(); imol++)
	{
		for (int ires = 0; ires < ori.molecules[imol].residues.size(); ires++)
		{
			for (int iatom = 0; iatom < ori.molecules[imol].residues[ires].atoms.size(); iatom++)
			{
				if( ((ori.molecules[imol].residues[ires].atoms[iatom]).coords).size() != 3 )
					REPORT_ERROR("helix.cpp::makeHelixFromPDBParticle(): Coordinates of atoms should have a dimension of 3!");
				mass_centre += (ori.molecules[imol].residues[ires].atoms[iatom]).coords;
			}
		}
	}
	mass_centre /= (RFLOAT)(nr_ori_atoms);

	aux0.clear();
	aux0 = ori;
	// Set the original particle on (r, 0, 0), if r > 0
	// Else just impose helical symmetry (make copies) according to the original particle
	if (do_center)
	{
		for (int imol = 0; imol < aux0.molecules.size(); imol++)
		{
			for (int ires = 0; ires < aux0.molecules[imol].residues.size(); ires++)
			{
				for (int iatom = 0; iatom < aux0.molecules[imol].residues[ires].atoms.size(); iatom++)
				{
					(aux0.molecules[imol].residues[ires].atoms[iatom]).coords -= mass_centre;
					XX((aux0.molecules[imol].residues[ires].atoms[iatom]).coords) += radius_A;
				}
			}
		}
		std::cout << "Centre of mass (in Angstroms) = " << XX(mass_centre) << ", " << YY(mass_centre) << ", " << ZZ(mass_centre) << std::endl;
		std::cout << "Bring the centre of mass to the (helical_radius, 0, 0)" << std::endl;
		std::cout << "Helical radius (for centre of mass) assigned = " << radius_A << " Angstroms" << std::endl;
	}

	// Construct the helix
	rotational_matrix.clear();
	shift.resize(3);
	shift.initZeros();
	helix.clear();
	helix.join(aux0);
	for (int ii = (((nr_copy + 1) % 2) - (nr_copy / 2)) ; ii <= (nr_copy / 2); ii++)
	{
		if (ii == 0)
			continue;

		rotation2DMatrix(((RFLOAT)(ii)) * twist_deg, rotational_matrix, true);
		ZZ(shift) = (RFLOAT)(ii) * rise_A;

		aux1.clear();
		aux1 = aux0;
		aux1.applyTransformation(rotational_matrix, shift);
		helix.join(aux1);
	}

	return;
}

/*
void normalise2DImageSlices(
		const FileName& fn_in,
		const FileName& fn_out,
		int bg_radius,
		RFLOAT white_dust_stddev,
		RFLOAT black_dust_stddev)
{
	Image<RFLOAT> stack_in, slice;
	int Xdim, Ydim, Zdim;
	long int Ndim, ii;

	if ( (fn_in.getExtension() != "mrcs") || (fn_out.getExtension() != "mrcs") )
	{
		REPORT_ERROR("helix.cpp::normalise2DImageSlices(): Input and output should be .mrcs files!");
	}
	stack_in.clear();
	stack_in.read(fn_in, false, -1, false, false); // readData = false, select_image = -1, mapData= false, is_2D = false);
	stack_in.getDimensions(Xdim, Ydim, Zdim, Ndim);
	std::cout << "File = " << fn_in.c_str() << std::endl;
	std::cout << "X, Y, Z, N dim = " << Xdim << ", " << Ydim << ", " << Zdim << ", " << Ndim << std::endl;
	std::cout << "bg_radius = " << bg_radius << ", white_dust_stddev = " << white_dust_stddev << ", black_dust_stddev = " << black_dust_stddev << std::endl;
	if( (Zdim != 1) || (Ndim < 1) || (Xdim < 3) || (Ydim < 3) )
	{
		REPORT_ERROR("helix.cpp::normalise2DImageSlices(): Invalid image dimensionality!");
	}

	for(ii = 0; ii < Ndim; ii++)
	{
		slice.read(fn_in, true, ii, false, false);
		normalise(slice, bg_radius, white_dust_stddev, black_dust_stddev, false);

		// Write this particle to the stack on disc
		// First particle: write stack in overwrite mode, from then on just append to it
		if (ii == 0)
			slice.write(fn_out, -1, (Ndim > 1), WRITE_OVERWRITE);
		else
			slice.write(fn_out, -1, false, WRITE_APPEND);
	}

	return;
}
*/

void applySoftSphericalMask(
		MultidimArray<RFLOAT>& v,
		RFLOAT sphere_diameter,
		RFLOAT cosine_width)
{
	RFLOAT r, r_max, r_max_edge;
	int dim = v.getDim();
	if (dim != 3)
		REPORT_ERROR("helix.cpp::applySoftSphericalMask(): Input image should have a dimension of 3!");
	if (cosine_width < 0.01)
		REPORT_ERROR("helix.cpp::applySoftSphericalMask(): Cosine width should be a positive value!");

	v.setXmippOrigin();
	r_max = (XSIZE(v) < YSIZE(v)) ? (XSIZE(v)) : (YSIZE(v));
	r_max = (r_max < ZSIZE(v)) ? (r_max) : (ZSIZE(v));
	// Nov11,2016 - Commented the following lines for r > 90% masks
	//if (cosine_width > 0.05 * r_max)
	//	r_max -= 2. * cosine_width;
	//r_max *= 0.45;
	if ( (sphere_diameter > 0.01) && ((sphere_diameter / 2.) < r_max) )
		r_max = sphere_diameter / 2.;
	r_max_edge = r_max + cosine_width;

	FOR_ALL_ELEMENTS_IN_ARRAY3D(v)
	{
		r = sqrt(k * k + i * i + j * j);
		if (r > r_max)
		{
			if (r < r_max_edge)
				A3D_ELEM(v, k, i, j) *= 0.5 + 0.5 * cos(PI * (r - r_max) / cosine_width);
			else
				A3D_ELEM(v, k, i, j) = 0.;
		}
	}
	return;
}

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
		bool bimodal_angular_priors,
		bool cut_into_segments)
{
	int total_segments, total_tubes, nr_segments, nr_tubes;
	FileName fns_in;
	std::vector<FileName> fn_in_list;
	MetaDataTable MD_out;

	fns_in = "*" + suffix_in;
	fns_in.globFiles(fn_in_list);
	if (fn_in_list.size() < 1)
		REPORT_ERROR("helix.cpp::extractHelicalSegmentsFromTubes_Multiple(): No input files are found!");

	total_segments = total_tubes = 0;
	for (int ii = 0; ii < fn_in_list.size(); ii++)
	{
		FileName fn_out;
		fn_out = fn_in_list[ii].beforeFirstOf(suffix_in) + suffix_out;
		if (format_tag == RELION_STAR_FORMAT)
			convertHelicalTubeCoordsToMetaDataTable(fn_in_list[ii], MD_out, nr_segments, nr_tubes, nr_asu, rise_A, pixel_size_A, Xdim, Ydim, box_size_pix, bimodal_angular_priors, cut_into_segments);
		else if (format_tag == XIMDISP_COORDS_FORMAT)
			convertXimdispHelicalTubeCoordsToMetaDataTable(fn_in_list[ii], MD_out, nr_segments, nr_tubes, nr_asu, rise_A, pixel_size_A, Xdim, Ydim, box_size_pix, bimodal_angular_priors, cut_into_segments);
		else if (format_tag == EMAN2_FORMAT)
			convertEmanHelicalTubeCoordsToMetaDataTable(fn_in_list[ii], MD_out, nr_segments, nr_tubes, nr_asu, rise_A, pixel_size_A, Xdim, Ydim, box_size_pix, bimodal_angular_priors, cut_into_segments);
		else
			REPORT_ERROR("helix.cpp::extractHelicalSegmentsFromTubes_Multiple(): BUG Invalid format tag!");
		total_segments += nr_segments;
		total_tubes += nr_tubes;
		MD_out.write(fn_out);
	}
	std::cout << " ### " << total_segments << " segments (" << total_tubes << " tubes, ~" << (total_segments * nr_asu)
			<< " subunits) are extracted from " << fn_in_list.size() << " input files. ###" << std::endl;
	return;
}

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
		bool bimodal_angular_priors,
		bool cut_into_segments)
{
	int nr_segments, MDobj_id;
	RFLOAT psi_deg, psi_rad, x1, y1, x2, y2, dx, dy, xp, yp, step_pix, half_box_size_pix, len_pix, psi_prior_flip_ratio, pitch;
	int id;
	MetaDataTable MD_in;
	std::vector<RFLOAT> x1_coord_list, y1_coord_list, x2_coord_list, y2_coord_list, pitch_list;
	std::vector<int> tube_id_list;

	// Check parameters and open files
	if ( (nr_asu < 1) || (rise_A < 0.001) || (pixel_size_A < 0.01) )
		REPORT_ERROR("helix.cpp::convertHelicalTubeCoordsToMetaDataTable(): Wrong parameters!");
	if ( (box_size_pix < 2) || (Xdim < box_size_pix) || (Ydim < box_size_pix))
		REPORT_ERROR("helix.cpp::convertHelicalTubeCoordsToMetaDataTable(): Wrong dimensions or box size!");
    if (fn_in.getExtension() != "star")
    	REPORT_ERROR("helix.cpp::convertHelicalTubeCoordsToMetaDataTable(): MetadataTable should have .star extension. Error(s) in " + fn_in);

    half_box_size_pix = box_size_pix / 2.;
    psi_prior_flip_ratio = UNIMODAL_PSI_PRIOR_FLIP_RATIO;
    if (bimodal_angular_priors)
    {
    	psi_prior_flip_ratio = BIMODAL_PSI_PRIOR_FLIP_RATIO;
    }

    // Read input STAR file
	MD_in.clear();
	MD_out.clear();
	MD_in.read(fn_in);
	if (MD_in.numberOfObjects() < 1) // Handle empty input files
		return;

    if ( (!MD_in.containsLabel(EMDL_IMAGE_COORD_X)) || (!MD_in.containsLabel(EMDL_IMAGE_COORD_Y)) )
    	REPORT_ERROR("helix.cpp::convertHelicalTubeCoordsToMetaDataTable(): Input STAR file does not contain X and Y coordinates! Error(s) in " + fn_in);
    if (MD_in.numberOfObjects() % 2)
    	REPORT_ERROR("helix.cpp::convertHelicalTubeCoordsToMetaDataTable(): Input coordinates should be in pairs! Error(s) in" + fn_in);

    // Sjors added MDin_has_id and MDin_has_pitch to allow manual calculation of different cross-over distances to be carried onto the extracted segments...
    bool MDin_has_id = MD_in.containsLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);
    bool MDin_has_pitch =  MD_in.containsLabel(EMDL_PARTICLE_HELICAL_TUBE_PITCH);
    x1_coord_list.clear();
    y1_coord_list.clear();
    x2_coord_list.clear();
    y2_coord_list.clear();
    tube_id_list.clear();
    pitch_list.clear();
    MDobj_id = 0;
    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_in)
    {
    	MDobj_id++;
    	MD_in.getValue(EMDL_IMAGE_COORD_X, xp);
    	MD_in.getValue(EMDL_IMAGE_COORD_Y, yp);
    	if (MDobj_id % 2)
    	{
    		x1_coord_list.push_back(xp);
    		y1_coord_list.push_back(yp);
    		if (MDin_has_id)
    		{
    			MD_in.getValue(EMDL_PARTICLE_HELICAL_TUBE_ID, id);
    			tube_id_list.push_back(id);
    		}
    		if (MDin_has_pitch)
    		{
    			MD_in.getValue(EMDL_PARTICLE_HELICAL_TUBE_PITCH, pitch);
    			pitch_list.push_back(pitch);
    		}
    	}
    	else
    	{
    		x2_coord_list.push_back(xp);
    		y2_coord_list.push_back(yp);
    	}
    }
    if ( (x1_coord_list.size() != x2_coord_list.size())
    		|| (x2_coord_list.size() != y1_coord_list.size())
			|| (y1_coord_list.size() != y2_coord_list.size())
			|| (y2_coord_list.size() != x1_coord_list.size()) )
    	REPORT_ERROR("helix.cpp::convertHelicalTubeCoordsToMetaDataTable(): BUG in reading input STAR file " + fn_in);
    MD_in.clear();

    // Init output STAR file
    MD_out.clear();
    MD_out.addLabel(EMDL_IMAGE_COORD_X);
    MD_out.addLabel(EMDL_IMAGE_COORD_Y);
    MD_out.addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);
    MD_out.addLabel(EMDL_ORIENT_TILT_PRIOR);
    MD_out.addLabel(EMDL_ORIENT_PSI_PRIOR);
    MD_out.addLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM);
    MD_out.addLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO);

    if (MDin_has_id)
    	MD_out.addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);
    if (MDin_has_pitch)
    	MD_out.addLabel(EMDL_PARTICLE_HELICAL_TUBE_PITCH);

	// Calculate all coordinates for helical segments
	nr_segments = 0;
    step_pix = nr_asu * rise_A / pixel_size_A;
    for (int tube_id = 0; tube_id < x1_coord_list.size(); tube_id++)
    {
    	x1 = x1_coord_list[tube_id];
    	y1 = y1_coord_list[tube_id];
    	x2 = x2_coord_list[tube_id];
    	y2 = y2_coord_list[tube_id];
    	if (MDin_has_id)
    		id = tube_id_list[tube_id];
    	if (MDin_has_pitch)
    		pitch = pitch_list[tube_id];

    	psi_rad = atan2(y2 - y1, x2 - x1);
    	psi_deg = RAD2DEG(psi_rad);
		dx = step_pix * cos(psi_rad);
		dy = step_pix * sin(psi_rad);

    	if (!cut_into_segments)
    	{
			MD_out.addObject();
	    	MD_out.setValue(EMDL_IMAGE_COORD_X, ((x1 + x2) / 2.));
	    	MD_out.setValue(EMDL_IMAGE_COORD_Y, ((y1 + y2) / 2.));
	    	MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, (tube_id + 1));
	    	MD_out.setValue(EMDL_ORIENT_TILT_PRIOR, 90.);
	    	MD_out.setValue(EMDL_ORIENT_PSI_PRIOR, -psi_deg);
	        MD_out.setValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, 0.);
	        MD_out.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, psi_prior_flip_ratio);

			if (MDin_has_id)
				MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, id);
			if (MDin_has_pitch)
				MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_PITCH, pitch);

	        nr_segments++;
    		continue;
    	}

    	xp = x1 - (dx * 0.99);
    	yp = y1 - (dy * 0.99);
    	len_pix = -step_pix;
    	while (1)
    	{
    		xp += dx;
    		yp += dy;
    		len_pix += step_pix;
    		if ( ((xp > x1) && (xp > x2)) || ((xp < x1) && (xp < x2))
    				|| ((yp > y1) && (yp > y2)) || ((yp < y1) && (yp < y2)) )
    		{
    			break;
    		}
    		else
    		{
#ifdef EXCLUDE_SEGMENTS_ON_THE_EDGES
    			// Avoid segments lying on the edges of the micrographs
    			if ( (xp < half_box_size_pix) || (xp > (Xdim - half_box_size_pix))
    					|| (yp < half_box_size_pix) || (yp > (Ydim - half_box_size_pix)) )
    			{
    				// Extract from filament start-end coordinates. It is not necessary to notify the user.
    				//std::cerr << " WARNING: Particle at (" << xp << ", " << yp << ") in coordinate file " << fn_in << " is NOT extracted because it is too close to the edge." << std::flush;
    				//std::cerr << " Box_size_pix = " << box_size_pix << ", Dimensions = " << Xdim << " * " << Ydim << std::endl;
    				continue;
    			}
#endif

    			MD_out.addObject();
    	    	MD_out.setValue(EMDL_IMAGE_COORD_X, xp);
    	    	MD_out.setValue(EMDL_IMAGE_COORD_Y, yp);
    	    	MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, (tube_id + 1));
    	    	MD_out.setValue(EMDL_ORIENT_TILT_PRIOR, 90.);
    	    	MD_out.setValue(EMDL_ORIENT_PSI_PRIOR, -psi_deg);
    	        MD_out.setValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, pixel_size_A * len_pix);
    	        MD_out.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, psi_prior_flip_ratio);

    			if (MDin_has_id)
    				MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, id);
    			if (MDin_has_pitch)
    				MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_PITCH, pitch);
   			nr_segments++;
    		}
    	}
    }

    if (nr_segments < 1)
    {
    	std::cerr << " WARNING: no segments extracted from file '" << fn_in << "'!" << std::endl;
    }
    /*
    else
    {
    	std::cout << "Input STAR file = " << fn_in << ", tubes = " << x1_coord_list.size()
    			<< ", segments = " << nr_segments << ", subunits ~ " << (nr_segments * nr_asu) << std::endl;
    }
    */
    total_segments = nr_segments;
    total_tubes = x1_coord_list.size();
}

void combineParticlePriorsWithKaiLocalCTF(
		FileName& fn_priors,
		FileName& fn_local_ctf,
		FileName& fn_combined)
{
	MetaDataTable MD_priors, MD_local_ctf;
	std::vector<RFLOAT> x, y, dU, dV, dAng, Cs, pix, mag, Q0, volt, fom, maxres, bfac, sfac, phase;
	RFLOAT _x, _y, _dU, _dV, _dAng, _Cs, _pix, _mag, _Q0, _volt, _fom, _maxres, _bfac, _sfac, _phase;
	int ii;

	if ( (fn_priors.getFileFormat() != "star") || (fn_local_ctf.getFileFormat() != "star") || (fn_combined.getFileFormat() != "star") )
		REPORT_ERROR("helix.cpp::combineParticlePriorsWithKaiLocalCTF(): MetaDataTable should have .star extension.");
//	if ( (fn_priors == fn_local_ctf) || (fn_local_ctf == fn_combined) || (fn_combined == fn_priors) )
//		REPORT_ERROR("helix.cpp::combineParticlePriorsWithKaiLocalCTF(): File names must be different.");
	if (fn_priors == fn_local_ctf)
		REPORT_ERROR("helix.cpp::combineParticlePriorsWithKaiLocalCTF(): Input file names must be different.");

	MD_priors.clear();
	MD_local_ctf.clear();
	MD_priors.read(fn_priors);
	MD_local_ctf.read(fn_local_ctf);
	if (MD_priors.numberOfObjects() != MD_local_ctf.numberOfObjects())
		REPORT_ERROR("helix.cpp::combineParticlePriorsWithKaiLocalCTF(): MetaDataTables to be combined are not of the same size.");

	if ( (!MD_priors.containsLabel(EMDL_IMAGE_COORD_X))
			|| (!MD_priors.containsLabel(EMDL_IMAGE_COORD_Y))
//			|| (!MD_local_ctf.containsLabel(EMDL_MICROGRAPH_NAME))
			|| (!MD_local_ctf.containsLabel(EMDL_IMAGE_COORD_X))
			|| (!MD_local_ctf.containsLabel(EMDL_IMAGE_COORD_Y))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_VOLTAGE))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_DEFOCUSU))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_DEFOCUSV))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_DEFOCUS_ANGLE))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_CS))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_MAGNIFICATION))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_Q0))
			)
		REPORT_ERROR("helix.cpp::combineParticlePriorsWithKaiLocalCTF(): Labels missing in MetaDataTables.");

	x.clear(); y.clear(); dU.clear(); dV.clear(); dAng.clear(); // necessary
	Cs.clear(); pix.clear(); mag.clear(); Q0.clear(); volt.clear(); // necessary
	fom.clear(); maxres.clear(); bfac.clear(); sfac.clear(); phase.clear(); // optional

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_local_ctf)
	{
		MD_local_ctf.getValue(EMDL_IMAGE_COORD_X, _x); x.push_back(_x);
		MD_local_ctf.getValue(EMDL_IMAGE_COORD_Y, _y); y.push_back(_y);
		MD_local_ctf.getValue(EMDL_CTF_DEFOCUSU, _dU); dU.push_back(_dU);
		MD_local_ctf.getValue(EMDL_CTF_DEFOCUSV, _dV); dV.push_back(_dV);
		MD_local_ctf.getValue(EMDL_CTF_DEFOCUS_ANGLE, _dAng); dAng.push_back(_dAng);
		MD_local_ctf.getValue(EMDL_CTF_CS, _Cs); Cs.push_back(_Cs);
		MD_local_ctf.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, _pix); pix.push_back(_pix);
		MD_local_ctf.getValue(EMDL_CTF_MAGNIFICATION, _mag); mag.push_back(_mag);
		MD_local_ctf.getValue(EMDL_CTF_Q0, _Q0); Q0.push_back(_Q0);
		MD_local_ctf.getValue(EMDL_CTF_VOLTAGE, _volt); volt.push_back(_volt);
		if (MD_local_ctf.containsLabel(EMDL_CTF_FOM))
			MD_local_ctf.getValue(EMDL_CTF_FOM, _fom); fom.push_back(_fom);
		if (MD_local_ctf.containsLabel(EMDL_CTF_MAXRES))
			MD_local_ctf.getValue(EMDL_CTF_MAXRES, _maxres); maxres.push_back(_maxres);
		if (MD_local_ctf.containsLabel(EMDL_CTF_BFACTOR))
			MD_local_ctf.getValue(EMDL_CTF_BFACTOR, _bfac); bfac.push_back(_bfac);
		if (MD_local_ctf.containsLabel(EMDL_CTF_SCALEFACTOR))
			MD_local_ctf.getValue(EMDL_CTF_SCALEFACTOR, _sfac); sfac.push_back(_sfac);
		if (MD_local_ctf.containsLabel(EMDL_CTF_PHASESHIFT))
			MD_local_ctf.getValue(EMDL_CTF_PHASESHIFT, _phase); phase.push_back(_phase);
	}

	if (!MD_priors.containsLabel(EMDL_CTF_DEFOCUSU))
		MD_priors.addLabel(EMDL_CTF_DEFOCUSU);
	if (!MD_priors.containsLabel(EMDL_CTF_DEFOCUSV))
		MD_priors.addLabel(EMDL_CTF_DEFOCUSV);
	if (!MD_priors.containsLabel(EMDL_CTF_DEFOCUS_ANGLE))
		MD_priors.addLabel(EMDL_CTF_DEFOCUS_ANGLE);
	if (!MD_priors.containsLabel(EMDL_CTF_CS))
		MD_priors.addLabel(EMDL_CTF_CS);
	if (!MD_priors.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
		MD_priors.addLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE);
	if (!MD_priors.containsLabel(EMDL_CTF_MAGNIFICATION))
		MD_priors.addLabel(EMDL_CTF_MAGNIFICATION);
	if (!MD_priors.containsLabel(EMDL_CTF_Q0))
		MD_priors.addLabel(EMDL_CTF_Q0);
	if (!MD_priors.containsLabel(EMDL_CTF_VOLTAGE))
		MD_priors.addLabel(EMDL_CTF_VOLTAGE);
	if (!MD_priors.containsLabel(EMDL_CTF_FOM))
		MD_priors.addLabel(EMDL_CTF_FOM);
	if (!MD_priors.containsLabel(EMDL_CTF_MAXRES))
		MD_priors.addLabel(EMDL_CTF_MAXRES);
	if (!MD_priors.containsLabel(EMDL_CTF_BFACTOR))
		MD_priors.addLabel(EMDL_CTF_BFACTOR);
	if (!MD_priors.containsLabel(EMDL_CTF_SCALEFACTOR))
		MD_priors.addLabel(EMDL_CTF_SCALEFACTOR);
	if (!MD_priors.containsLabel(EMDL_CTF_PHASESHIFT))
		MD_priors.addLabel(EMDL_CTF_PHASESHIFT);

	ii = -1;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_priors)
	{
		ii++;
		MD_priors.getValue(EMDL_IMAGE_COORD_X, _x);
		MD_priors.getValue(EMDL_IMAGE_COORD_Y, _y);
		if ( (fabs(x[ii] - _x) > 1.001) || (fabs(y[ii] - _y) > 1.001) )
			REPORT_ERROR("helix.cpp::combineParticlePriorsWithKaiLocalCTF(): Coordinates from the two MetaDataTables do not match.");

		//MD_priors.setValue(EMDL_IMAGE_COORD_X, x[ii]);
		//MD_priors.setValue(EMDL_IMAGE_COORD_Y, y[ii]);
		MD_priors.setValue(EMDL_CTF_DEFOCUSU, dU[ii]);
		MD_priors.setValue(EMDL_CTF_DEFOCUSV, dV[ii]);
		MD_priors.setValue(EMDL_CTF_DEFOCUS_ANGLE, dAng[ii]);
		MD_priors.setValue(EMDL_CTF_CS, Cs[ii]);
		MD_priors.setValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, pix[ii]);
		MD_priors.setValue(EMDL_CTF_MAGNIFICATION, mag[ii]);
		MD_priors.setValue(EMDL_CTF_Q0, Q0[ii]);
		MD_priors.setValue(EMDL_CTF_VOLTAGE, volt[ii]);
		if (MD_local_ctf.containsLabel(EMDL_CTF_FOM))
			MD_priors.setValue(EMDL_CTF_FOM, fom[ii]);
		if (MD_local_ctf.containsLabel(EMDL_CTF_MAXRES))
			MD_priors.setValue(EMDL_CTF_MAXRES, maxres[ii]);
		if (MD_local_ctf.containsLabel(EMDL_CTF_BFACTOR))
			MD_priors.setValue(EMDL_CTF_BFACTOR, bfac[ii]);
		if (MD_local_ctf.containsLabel(EMDL_CTF_SCALEFACTOR))
			MD_priors.setValue(EMDL_CTF_SCALEFACTOR, sfac[ii]);
		if (MD_local_ctf.containsLabel(EMDL_CTF_PHASESHIFT))
			MD_priors.setValue(EMDL_CTF_PHASESHIFT, phase[ii]);
	}

	MD_priors.write(fn_combined);
	return;
}

void combineParticlePriorsWithKaiLocalCTF_Multiple(
		std::string& suffix_priors,
		std::string& suffix_local_ctf,
		std::string& suffix_combined)
{
	FileName fns_priors;
	std::vector<FileName> fn_priors_list;

//	if ( (suffix_priors == suffix_local_ctf) || (suffix_priors == suffix_combined) || (suffix_combined == suffix_priors) )
//		REPORT_ERROR("helix.cpp::combineParticlePriorsWithKaiLocalCTF_Multiple(): File names error!");
	if (suffix_priors == suffix_local_ctf)
		REPORT_ERROR("helix.cpp::combineParticlePriorsWithKaiLocalCTF_Multiple(): Input file names error!");

	fns_priors = "*" + suffix_priors;
	fns_priors.globFiles(fn_priors_list);
	std::cout << "Number of input files = " << fn_priors_list.size() << std::endl;
	if (fn_priors_list.size() < 1)
		REPORT_ERROR("helix.cpp::combineParticlePriorsWithKaiLocalCTF_Multiple(): No input files are found!");

	for (int ii = 0; ii < fn_priors_list.size(); ii++)
	{
		FileName fn_local_ctf, fn_combined;
		fn_local_ctf = fn_priors_list[ii].beforeFirstOf(suffix_priors) + suffix_local_ctf;
		fn_combined = fn_priors_list[ii].beforeFirstOf(suffix_priors) + suffix_combined;
		combineParticlePriorsWithKaiLocalCTF(fn_priors_list[ii], fn_local_ctf, fn_combined);
	}
	return;
}

void setNullTiltPriorsInDataStar(
		FileName& fn_in,
		FileName& fn_out)
{
	MetaDataTable MD;
	if ( (fn_in.getFileFormat() != "star") || (fn_out.getFileFormat() != "star") )
		REPORT_ERROR("helix.cpp::addNullTiltPriorsToDataStar(): MetaDataTable should have .star extension.");
	if (fn_in == fn_out)
		REPORT_ERROR("helix.cpp::addNullTiltPriorsToDataStar(): File names must be different.");

	MD.read(fn_in);
	if (!MD.containsLabel(EMDL_ORIENT_TILT_PRIOR))
		MD.addLabel(EMDL_ORIENT_TILT_PRIOR);
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		MD.setValue(EMDL_ORIENT_TILT_PRIOR, 90.);
	}
	MD.write(fn_out);
	return;
}

void removeBadTiltHelicalSegmentsFromDataStar(
		FileName& fn_in,
		FileName& fn_out,
		RFLOAT max_dev_deg)
{
	MetaDataTable MD_in, MD_out;
	int nr_segments_old, nr_segments_new;
	RFLOAT tilt_deg;
	if ( (max_dev_deg < 0.) || (max_dev_deg > 89.) )
		REPORT_ERROR("helix.cpp::removeBadTiltParticlesFromDataStar(): Max deviations of tilt angles from 90 degree should be in the range of 0~89 degrees.");
	if ( (fn_in.getFileFormat() != "star") || (fn_out.getFileFormat() != "star") )
		REPORT_ERROR("helix.cpp::removeBadTiltParticlesFromDataStar(): MetaDataTable should have .star extension.");
	if (fn_in == fn_out)
		REPORT_ERROR("helix.cpp::removeBadTiltParticlesFromDataStar(): File names must be different.");

	MD_in.clear();
	MD_out.clear();
	MD_in.read(fn_in);
	// TODO: Use EMDL_ORIENT_TILT or EMDL_ORIENT_TILT_PRIOR ?
	if (!MD_in.containsLabel(EMDL_ORIENT_TILT))
		REPORT_ERROR("helix.cpp::removeBadTiltParticlesFromDataStar(): Input .star file contains no tilt angles.");

	nr_segments_old = nr_segments_new = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_in)
	{
		nr_segments_old++;
		MD_in.getValue(EMDL_ORIENT_TILT, tilt_deg);
		if (fabs(tilt_deg - 90.) < max_dev_deg)
		{
			nr_segments_new++;
			MD_out.addObject(MD_in.getObject());
		}
	}
	MD_out.write(fn_out);
	std::cout << " Number of segments (input / output) = " << nr_segments_old << " / " << nr_segments_new << std::endl;
	return;
}

void removeBadPsiHelicalSegmentsFromDataStar(
		FileName& fn_in,
		FileName& fn_out,
		RFLOAT max_dev_deg)
{
	MetaDataTable MD_in, MD_out;
	int nr_segments_old, nr_segments_new;
	RFLOAT psi_deg, psi_prior_deg, diff_psi;
	if ( (max_dev_deg < 0.) || (max_dev_deg > 89.) )
		REPORT_ERROR("helix.cpp::removeBadPsiParticlesFromDataStar(): Max deviations of tilt angles from 90 degree should be in the range of 0~89 degrees.");
	if ( (fn_in.getFileFormat() != "star") || (fn_out.getFileFormat() != "star") )
		REPORT_ERROR("helix.cpp::removeBadPsiParticlesFromDataStar(): MetaDataTable should have .star extension.");
	if (fn_in == fn_out)
		REPORT_ERROR("helix.cpp::removeBadPsiParticlesFromDataStar(): File names must be different.");

	MD_in.clear();
	MD_out.clear();
	MD_in.read(fn_in);
	if ( (!MD_in.containsLabel(EMDL_ORIENT_PSI)) || (!MD_in.containsLabel(EMDL_ORIENT_PSI_PRIOR)) )
		REPORT_ERROR("helix.cpp::removeBadTiltParticlesFromDataStar(): Input .star file contains no psi angles with their priors.");

	nr_segments_old = nr_segments_new = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_in)
	{
		nr_segments_old++;
		MD_in.getValue(EMDL_ORIENT_PSI, psi_deg);
		MD_in.getValue(EMDL_ORIENT_PSI_PRIOR, psi_prior_deg);
		diff_psi = ABS(psi_deg - psi_prior_deg);
		if (diff_psi > 180.)
			diff_psi = ABS(diff_psi - 360.);
		if (diff_psi < max_dev_deg)
		{
			nr_segments_new++;
			MD_out.addObject(MD_in.getObject());
		}
	}
	MD_out.write(fn_out);
	std::cout << " Number of segments (input / output) = " << nr_segments_old << " / " << nr_segments_new << std::endl;
	return;
}

void convertHelicalSegmentCoordsToStarFile_Multiple(
		FileName& suffix_coords,
		FileName& suffix_out,
		int format_tag,
		RFLOAT pixel_size_A,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT boxsize,
		bool bimodal_angular_priors)
{
	int total_segments, nr_segments, total_tubes, nr_tubes;
	FileName fns_coords;
	std::vector<FileName> fn_coords_list;
	MetaDataTable MD_out;

	fns_coords = "*" + suffix_coords;
	fns_coords.globFiles(fn_coords_list);
	if (fn_coords_list.size() < 1)
		REPORT_ERROR("helix.cpp::convertHelicalCoordsToStarFile_Multiple(): No input files are found!");

	total_segments = total_tubes = 0;
	for (int ii = 0; ii < fn_coords_list.size(); ii++)
	{
		FileName fn_out;
		fn_out = fn_coords_list[ii].beforeFirstOf(suffix_coords) + suffix_out;
		if (format_tag == XIMDISP_COORDS_FORMAT)
			convertXimdispHelicalSegmentCoordsToMetaDataTable(fn_coords_list[ii], MD_out, nr_segments, nr_tubes, pixel_size_A, Xdim, Ydim, boxsize, bimodal_angular_priors);
		else if (format_tag == EMAN2_FORMAT)
			convertEmanHelicalSegmentCoordsToMetaDataTable(fn_coords_list[ii], MD_out, nr_segments, nr_tubes, pixel_size_A, Xdim, Ydim, boxsize, bimodal_angular_priors);
		else
			REPORT_ERROR("helix.cpp::convertHelicalCoordsToStarFile_Multiple(): BUG Invalid format tag!");
		total_segments += nr_segments;
		total_tubes += nr_tubes;
		MD_out.write(fn_out);
	}
	std::cout << " ### " << total_segments << " segments (" << total_tubes << " tubes) are extracted from " << fn_coords_list.size() << " input files. ###" << std::endl;
	return;
}

void convertHelicalSegmentCoordsToMetaDataTable(
		FileName& fn_in,
		MetaDataTable& MD_out,
		int& total_segments,
		bool is_3D_data,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT Zdim,
		RFLOAT box_size_pix,
		bool bimodal_angular_priors)
{
	MetaDataTable MD_in;
	if (fn_in.getExtension() != "star")
		REPORT_ERROR("helix.cpp::convertHelicalSegmentCoordsToMetaDataTable(): Input file should have .star extension!!");
	if ( (box_size_pix < 2) || (Xdim < box_size_pix) || (Ydim < box_size_pix) || ((is_3D_data) && (Zdim < box_size_pix)) )
		REPORT_ERROR("helix.cpp::convertHelicalSegmentCoordsToMetaDataTable(): Wrong dimensions or box size!");

	RFLOAT x = 0., y = 0., z = 0.;
	RFLOAT half_box_size_pix = box_size_pix / 2.;
	RFLOAT psi_prior_flip_ratio = UNIMODAL_PSI_PRIOR_FLIP_RATIO;
	if (bimodal_angular_priors)
	{
		psi_prior_flip_ratio = BIMODAL_PSI_PRIOR_FLIP_RATIO;
	}

	MD_in.clear();
	MD_out.clear();
	MD_in.read(fn_in);
	if (MD_in.numberOfObjects() < 1) // Handle empty input files
		return;

	if ( (!MD_in.containsLabel(EMDL_IMAGE_COORD_X))
			|| (!MD_in.containsLabel(EMDL_IMAGE_COORD_Y))
			|| ( (is_3D_data) && (!MD_in.containsLabel(EMDL_IMAGE_COORD_Z)) )
			|| (!MD_in.containsLabel(EMDL_PARTICLE_HELICAL_TUBE_ID))
			|| (!MD_in.containsLabel(EMDL_ORIENT_TILT_PRIOR))
			|| (!MD_in.containsLabel(EMDL_ORIENT_PSI_PRIOR))
			|| (!MD_in.containsLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM)) )
		REPORT_ERROR("helix.cpp::convertHelicalSegmentCoordsToMetaDataTable(): Prior information of helical segments is missing in " + fn_in);

	int nr_segments = 0;
	z = 1.;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_in)
	{
		MD_in.getValue(EMDL_IMAGE_COORD_X, x);
		MD_in.getValue(EMDL_IMAGE_COORD_Y, y);
		if (is_3D_data)
			MD_in.getValue(EMDL_IMAGE_COORD_Z, z);
#ifdef EXCLUDE_SEGMENTS_ON_THE_EDGES
		// Avoid segments lying on the edges of the micrographs
		if ( (x < half_box_size_pix) || (x > (Xdim - half_box_size_pix))
				|| (y < half_box_size_pix) || (y > (Ydim - half_box_size_pix))
				|| ( (is_3D_data) && ((z < half_box_size_pix) || (z > (Zdim - half_box_size_pix))) ) )
		{
			std::cerr << " WARNING: Particle at (" << x << ", " << y << ", " << z << std::flush;
			std::cerr << ") in coordinate file " << fn_in << " is NOT extracted because it is too close to the edge." << std::flush;
			std::cerr << " Box_size_pix = " << box_size_pix << std::flush;
			std::cerr << ", Dimensions = " << Xdim << " * " << Ydim << " * " << Zdim << std::endl;
			continue;
		}
#endif
		nr_segments++;
		MD_out.addObject(MD_in.getObject());

		// TODO: check whether there is a bug...
		MD_out.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, psi_prior_flip_ratio);

	}
	total_segments = nr_segments;
	MD_in.clear();
}

void convertXimdispHelicalSegmentCoordsToMetaDataTable(
		FileName& fn_in,
		MetaDataTable& MD_out,
		int& total_segments,
		int& total_tubes,
		RFLOAT pixel_size_A,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix,
		bool bimodal_angular_priors)
{
	int nr_segments_on_edges, nr_segments, nr_tubes;
	RFLOAT x, y, x_old, y_old, psi_deg_old, psi_deg, half_box_size_pix, len_pix, psi_prior_flip_ratio;
	std::ifstream fin;
	std::string line;
	std::vector<std::string> words;

	if ( (box_size_pix < 2) || (Xdim < box_size_pix) || (Ydim < box_size_pix))
		REPORT_ERROR("helix.cpp::convertXimdispHelicalSegmentCoordsToMetaDataTable(): Wrong dimensions or box size!");

	// Header of output file
	MD_out.clear();
	MD_out.addLabel(EMDL_IMAGE_COORD_X);
	MD_out.addLabel(EMDL_IMAGE_COORD_Y);
	MD_out.addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);
	MD_out.addLabel(EMDL_ORIENT_TILT_PRIOR);
	MD_out.addLabel(EMDL_ORIENT_PSI_PRIOR);
    MD_out.addLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM);
    MD_out.addLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO);

    half_box_size_pix = box_size_pix / 2.;
    psi_prior_flip_ratio = UNIMODAL_PSI_PRIOR_FLIP_RATIO;
    if (bimodal_angular_priors)
    {
    	psi_prior_flip_ratio = BIMODAL_PSI_PRIOR_FLIP_RATIO;
    }

	fin.open(fn_in.c_str(), std::ios_base::in);
	if (fin.fail())
		REPORT_ERROR("helix.cpp::convertXimdispHelicalSegmentCoordsToMetaDataTable(): Cannot open input file " + fn_in);
	line.clear();
	getline(fin, line, '\n');
	psi_deg_old = 999999.;
	len_pix = 0.;
	nr_segments_on_edges = nr_segments = nr_tubes = 0;
	while (getline(fin, line, '\n'))
	{
		// Read in a new line of x, y, psi
		if (line.size() < 2) // End of file
			break;
		words.clear();
		tokenize(line, words);
		if (words.size() != 3)
			REPORT_ERROR("helix.cpp::convertXimdispHelicalSegmentCoordsToMetaDataTable(): Invalid input file " + fn_in);
		x = textToFloat(words[0]);
		y = textToFloat(words[1]);
		psi_deg = textToFloat(words[2]);

		// Check whether it is on a new helical tube
		if (fabs(psi_deg - psi_deg_old) > 0.1)
		{
			nr_tubes++;
			len_pix = 0.;
			x_old = x;
			y_old = y;
		}

		// Accumulate the length
		len_pix += sqrt( (x - x_old) * (x - x_old) + (y - y_old) * (y - y_old) );
		x_old = x;
		y_old = y;
		psi_deg_old = psi_deg;

#ifdef EXCLUDE_SEGMENTS_ON_THE_EDGES
		// Avoid segments lying on the edges of the micrographs
		if ( (x < half_box_size_pix) || (x > (Xdim - half_box_size_pix)) || (y < half_box_size_pix) || (y > (Ydim - half_box_size_pix)) )
		{
			nr_segments_on_edges++;

			std::cerr << " WARNING: Particle at (" << x << ", " << y << ") in coordinate file " << fn_in << " is NOT extracted because it is too close to the edge." << std::flush;
			std::cerr << " Box_size_pix = " << box_size_pix << ", Dimensions = " << Xdim << " * " << Ydim << std::endl;
			continue;
		}
#endif

		nr_segments++;
		MD_out.addObject();
		MD_out.setValue(EMDL_IMAGE_COORD_X, x);
		MD_out.setValue(EMDL_IMAGE_COORD_Y, y);
		MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, nr_tubes);
		MD_out.setValue(EMDL_ORIENT_TILT_PRIOR, 90.);
		MD_out.setValue(EMDL_ORIENT_PSI_PRIOR, -psi_deg);
	    MD_out.setValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, pixel_size_A * len_pix);
	    MD_out.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, psi_prior_flip_ratio);

	    line.clear();
	}
	fin.close();
	total_segments = nr_segments;
	total_tubes = nr_tubes;
	std::cout << "Input XIMDISP coordinates = " << fn_in.c_str() << ", micrograph size = " << Xdim << " * " << Ydim << ", box size = " << box_size_pix
			<< ", tubes = " << nr_tubes << ", " << nr_segments_on_edges << " segments excluded, " << nr_segments << " segments left." << std::endl;
}

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
		bool bimodal_angular_priors,
		bool cut_into_segments)
{
	int nr_segments, nr_tubes;
	RFLOAT xp, yp, dx, dy, x1, y1, x2, y2, psi_deg, psi_rad, half_box_size_pix, len_pix, psi_prior_flip_ratio;
	std::ifstream fin;
	std::string line;
	std::vector<std::string> words;
	std::vector<RFLOAT> x, y;

	if ( (box_size_pix < 2) || (Xdim < box_size_pix) || (Ydim < box_size_pix))
		REPORT_ERROR("helix.cpp::convertXimdispHelicalTubeCoordsToMetaDataTable(): Wrong dimensions or box size!");
	if (pixel_size_A < 0.001)
		REPORT_ERROR("helix.cpp::convertXimdispHelicalTubeCoordsToMetaDataTable(): Invalid pixel size!");
	RFLOAT step_pix = ((RFLOAT)(nr_asu)) * rise_A / pixel_size_A;
	if ( (nr_asu < 1) || (rise_A < 0.001) || (step_pix < 0.001) )
		REPORT_ERROR("helix.cpp::convertXimdispHelicalTubeCoordsToMetaDataTable(): Invalid helical rise or number of asymmetrical units!");

	// Header of output file
	MD_out.clear();
	MD_out.addLabel(EMDL_IMAGE_COORD_X);
	MD_out.addLabel(EMDL_IMAGE_COORD_Y);
	MD_out.addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);
	MD_out.addLabel(EMDL_ORIENT_TILT_PRIOR);
	MD_out.addLabel(EMDL_ORIENT_PSI_PRIOR);
    MD_out.addLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH);
    MD_out.addLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO);

	x.resize(4);
	y.resize(4);
	half_box_size_pix = box_size_pix / 2.;
    psi_prior_flip_ratio = UNIMODAL_PSI_PRIOR_FLIP_RATIO;
    if (bimodal_angular_priors)
    {
		psi_prior_flip_ratio = BIMODAL_PSI_PRIOR_FLIP_RATIO;
    }

	fin.open(fn_in.c_str(), std::ios_base::in);
	if (fin.fail())
		REPORT_ERROR("helix.cpp::convertXimdispHelicalTubeCoordsToMetaDataTable(): Cannot open input file " + fn_in);
	nr_tubes = nr_segments = 0;
	while (getline(fin, line, '\n'))
	{
		// Read in new helical tube
		if (line.size() < 2) // End of file
			break;
		words.clear();
		tokenize(line, words);
		if (words.size() != 2)
			REPORT_ERROR("helix.cpp::convertXimdispHelicalTubeCoordsToMetaDataTable(): Invalid input file " + fn_in);
		nr_tubes++;

		// Read in starting and end points for this helical tube
		for (int iline = 0; iline < 4; iline++)
		{
			line.clear();
			getline(fin, line, '\n');
			words.clear();
			tokenize(line, words);
			if (words.size() != 2)
				REPORT_ERROR("helix.cpp::convertXimdispHelicalTubeCoordsToMetaDataTable(): Invalid input file " + fn_in);
			x[iline] = textToFloat(words[0]);
			y[iline] = textToFloat(words[1]);
		}
		line.clear();
		getline(fin, line, '\n');

		// Calculate starting and end points for this helical tube
		x1 = (x[0] + x[1]) / 2.;
		y1 = (y[0] + y[1]) / 2.;
		x2 = (x[2] + x[3]) / 2.;
		y2 = (y[2] + y[3]) / 2.;
		psi_rad = atan2(y2 - y1, x2 - x1);
		psi_deg = RAD2DEG(psi_rad);
		dx = step_pix * cos(psi_rad);
		dy = step_pix * sin(psi_rad);

		if (!cut_into_segments)
		{
			MD_out.addObject();
	    	MD_out.setValue(EMDL_IMAGE_COORD_X, ((x1 + x2) / 2.));
	    	MD_out.setValue(EMDL_IMAGE_COORD_Y, ((y1 + y2) / 2.));
	    	MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, nr_tubes);
	    	MD_out.setValue(EMDL_ORIENT_TILT_PRIOR, 90.);
	    	MD_out.setValue(EMDL_ORIENT_PSI_PRIOR, -psi_deg);
	        MD_out.setValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, 0.);
	        MD_out.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, psi_prior_flip_ratio);

			continue;
		}

		// Calculate coordinates for all segments
		xp = x1 - (dx * 0.99);
		yp = y1 - (dy * 0.99);
		len_pix = -step_pix;
    	while (1)
    	{
    		xp += dx;
    		yp += dy;
    		len_pix += step_pix;
    		if ( ((xp > x1) && (xp > x2)) || ((xp < x1) && (xp < x2))
    				|| ((yp > y1) && (yp > y2)) || ((yp < y1) && (yp < y2)) )
    		{
    			break;
    		}
    		else
    		{
#ifdef EXCLUDE_SEGMENTS_ON_THE_EDGES
    			// Avoid segments lying on the edges of the micrographs
    			if ( (xp < half_box_size_pix) || (xp > (Xdim - half_box_size_pix))
    					|| (yp < half_box_size_pix) || (yp > (Ydim - half_box_size_pix)) )
    			{
    				// Extract from filament start-end coordinates. It is not necessary to notify the user.
    				//std::cerr << " WARNING: Particle at (" << xp << ", " << yp << ") in coordinate file " << fn_in << " is NOT extracted because it is too close to the edge." << std::flush;
    				//std::cerr << " Box_size_pix = " << box_size_pix << ", Dimensions = " << Xdim << " * " << Ydim << std::endl;
    				continue;
    			}
#endif

    			MD_out.addObject();
    	    	MD_out.setValue(EMDL_IMAGE_COORD_X, xp);
    	    	MD_out.setValue(EMDL_IMAGE_COORD_Y, yp);
    	    	MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, nr_tubes);
    	    	MD_out.setValue(EMDL_ORIENT_TILT_PRIOR, 90.);
    	    	MD_out.setValue(EMDL_ORIENT_PSI_PRIOR, -psi_deg);
    	        MD_out.setValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, pixel_size_A * len_pix);
    	        MD_out.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, psi_prior_flip_ratio);

    			nr_segments++;
    		}
    	}
		line.clear();
	}
	fin.close();

	total_segments = MD_out.numberOfObjects();
	total_tubes = nr_tubes;
	std::cout << "Input XIMDISP coordinates = " << fn_in.c_str() << ", micrograph size = " << Xdim << " * " << Ydim
			<< ", box size = " << box_size_pix << ", tubes = " << nr_tubes << ", segments = " << MD_out.numberOfObjects() << "." << std::endl;
}

void convertEmanHelicalSegmentCoordsToMetaDataTable(
		FileName& fn_in,
		MetaDataTable& MD_out,
		int& total_segments,
		int& total_tubes,
		RFLOAT pixel_size_A,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix,
		bool bimodal_angular_priors)
{
	int nr_segments_on_edges, nr_segments, nr_tubes;
	RFLOAT x, y, x_old, y_old, psi_deg, half_box_size_pix, len_pix, width, psi_prior_flip_ratio;
	std::ifstream fin;
	std::string line;
	std::vector<std::string> words;

	if ( (box_size_pix < 2) || (Xdim < box_size_pix) || (Ydim < box_size_pix))
		REPORT_ERROR("helix.cpp::convertEmanHelicalSegmentCoordsToMetaDataTable(): Wrong dimensions or box size!");

	// Header of output file
	MD_out.clear();
	MD_out.addLabel(EMDL_IMAGE_COORD_X);
	MD_out.addLabel(EMDL_IMAGE_COORD_Y);
	MD_out.addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);
	MD_out.addLabel(EMDL_ORIENT_TILT_PRIOR);
	MD_out.addLabel(EMDL_ORIENT_PSI_PRIOR);
    MD_out.addLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM);
    MD_out.addLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO);

    half_box_size_pix = box_size_pix / 2.;
    psi_prior_flip_ratio = UNIMODAL_PSI_PRIOR_FLIP_RATIO;
    if (bimodal_angular_priors)
	{
    	psi_prior_flip_ratio = BIMODAL_PSI_PRIOR_FLIP_RATIO;
	}

	fin.open(fn_in.c_str(), std::ios_base::in);
	if (fin.fail())
		REPORT_ERROR("helix.cpp::convertEmanHelicalSegmentCoordsToMetaDataTable(): Cannot open input file " + fn_in);
	line.clear();
	len_pix = 0.;
	nr_segments_on_edges = nr_segments = nr_tubes = 0;
	while (getline(fin, line, '\n'))
	{
		if (line.size() < 2) // End of file
			break;

		// Find lines which start with "#helix: ("
		int char_offset = 9;
		if ((line.substr(0, char_offset + 1)).find("#helix: (") != std::string::npos)
		{
			nr_tubes++;

			// Get psi angle
			RFLOAT x1, y1, x2, y2;
			char cdummy;
			std::istringstream ss(line.substr(char_offset));
			//std::cout << line.substr(char_offset) << std::endl;
			ss >> x1 >> cdummy >> y1 >> cdummy >> cdummy >> cdummy >> x2 >> cdummy >> y2 >> cdummy >> cdummy >> width;
			psi_deg = RAD2DEG(atan2(y2 - y1, x2 - x1));
			len_pix = 0.;
			x_old = y_old = (1.1e30);
			continue;
		}
		else if (line[0] == '#')
		{
			continue;
		}

		// Get x, y coordinates
		words.clear();
		tokenize(line, words);
		if (words.size() != 2)
			REPORT_ERROR("helix.cpp::convertEmanHelicalSegmentCoordsToMetaDataTable(): Invalid input file " + fn_in);
		x = textToFloat(words[0]);
		y = textToFloat(words[1]);

		// Accumulate the length
		if (x_old < (1e30))
			len_pix += sqrt( (x - x_old) * (x - x_old) + (y - y_old) * (y - y_old) );
		x_old = x;
		y_old = y;

#ifdef EXCLUDE_SEGMENTS_ON_THE_EDGES
		// Avoid segments lying on the edges of the micrographs
		if ( (x < half_box_size_pix) || (x > (Xdim - half_box_size_pix)) || (y < half_box_size_pix) || (y > (Ydim - half_box_size_pix)) )
		{
			nr_segments_on_edges++;

			std::cerr << " WARNING: Particle at (" << x << ", " << y << ") in coordinate file " << fn_in << " is NOT extracted because it is too close to the edge." << std::flush;
			std::cerr << " Box_size_pix = " << box_size_pix << ", Dimensions = " << Xdim << " * " << Ydim << std::endl;
			continue;
		}
#endif

		nr_segments++;
		MD_out.addObject();
		MD_out.setValue(EMDL_IMAGE_COORD_X, x);
		MD_out.setValue(EMDL_IMAGE_COORD_Y, y);
		MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, nr_tubes);
		MD_out.setValue(EMDL_ORIENT_TILT_PRIOR, 90.);
		MD_out.setValue(EMDL_ORIENT_PSI_PRIOR, -psi_deg);
	    MD_out.setValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, pixel_size_A * len_pix);
	    MD_out.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, psi_prior_flip_ratio);

	    line.clear();
	}

	fin.close();
	total_segments = nr_segments;
	total_tubes = nr_tubes;
	std::cout << "Input EMAN2 coordinates = " << fn_in.c_str() << ", micrograph size = " << Xdim << " * " << Ydim << ", box size = " << box_size_pix
			<< ", tubes = " << nr_tubes << ", " << nr_segments_on_edges << " segments excluded, " << nr_segments << " segments left." << std::endl;
}

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
		bool bimodal_angular_priors,
		bool cut_into_segments)
{
	int nr_segments, nr_tubes;
	RFLOAT xp, yp, dx, dy, x1, y1, x2, y2, psi_deg, psi_rad, half_box_size_pix, len_pix, psi_prior_flip_ratio;
	std::ifstream fin;
	std::string line;
	std::vector<std::string> words;

	if ( (box_size_pix < 2) || (Xdim < box_size_pix) || (Ydim < box_size_pix))
		REPORT_ERROR("helix.cpp::convertEmanHelicalTubeCoordsToMetaDataTable(): Wrong dimensions or box size!");
	if (pixel_size_A < 0.001)
		REPORT_ERROR("helix.cpp::convertEmanHelicalTubeCoordsToMetaDataTable(): Invalid pixel size!");
	RFLOAT step_pix = ((RFLOAT)(nr_asu)) * rise_A / pixel_size_A;
	if ( (nr_asu < 1) || (rise_A < 0.001) || (step_pix < 0.001) )
		REPORT_ERROR("helix.cpp::convertEmanHelicalTubeCoordsToMetaDataTable(): Invalid helical rise or number of asymmetrical units!");

	// Header of output file
	MD_out.clear();
	MD_out.addLabel(EMDL_IMAGE_COORD_X);
	MD_out.addLabel(EMDL_IMAGE_COORD_Y);
	MD_out.addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);
	MD_out.addLabel(EMDL_ORIENT_TILT_PRIOR);
	MD_out.addLabel(EMDL_ORIENT_PSI_PRIOR);
    MD_out.addLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM);
    MD_out.addLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO);

	half_box_size_pix = box_size_pix / 2.;
    psi_prior_flip_ratio = UNIMODAL_PSI_PRIOR_FLIP_RATIO;
    if (bimodal_angular_priors)
    {
    	psi_prior_flip_ratio = BIMODAL_PSI_PRIOR_FLIP_RATIO;
    }

	fin.open(fn_in.c_str(), std::ios_base::in);
	if (fin.fail())
		REPORT_ERROR("helix.cpp::convertEmanHelicalTubeCoordsToMetaDataTable(): Cannot open input file " + fn_in);
	nr_tubes = nr_segments = 0;
	line.clear();
	while (getline(fin, line, '\n'))
	{
		RFLOAT width1, width2, width3, width4;
		int tag;

		// Read in new helical tube
		if (line.size() < 2) // End of file
			break;

		// Get x1, y1
		words.clear();
		tokenize(line, words);
		if (words.size() != 5)
			REPORT_ERROR("helix.cpp::convertEmanHelicalTubeCoordsToMetaDataTable(): Invalid input file " + fn_in);
		x1 = textToFloat(words[0]);
		y1 = textToFloat(words[1]);
		width1 = textToFloat(words[2]);
		width2 = textToFloat(words[3]);
		tag = textToInteger(words[4]);
		if ( (tag != (-1)) || (fabs(width1 - width2) > 0.01) )
			REPORT_ERROR("helix.cpp::convertEmanHelicalTubeCoordsToMetaDataTable(): Invalid input file " + fn_in);
		x1 += width1 / 2.;
		y1 += width1 / 2.;

		// Get x2, y2
		line.clear();
		getline(fin, line, '\n');
		words.clear();
		tokenize(line, words);
		if (words.size() != 5)
			REPORT_ERROR("helix.cpp::convertEmanHelicalTubeCoordsToMetaDataTable(): Invalid input file " + fn_in);
		x2 = textToFloat(words[0]);
		y2 = textToFloat(words[1]);
		width3 = textToFloat(words[2]);
		width4 = textToFloat(words[3]);
		tag = textToInteger(words[4]);
		if ( (tag != (-2)) || (fabs(width3 - width4) > 0.01) || (fabs(width3 - width1) > 0.01) )
			REPORT_ERROR("helix.cpp::convertEmanHelicalTubeCoordsToMetaDataTable(): Invalid input file " + fn_in);
		x2 += width3 / 2.;
		y2 += width3 / 2.;

		nr_tubes++;

		psi_rad = atan2(y2 - y1, x2 - x1);
		psi_deg = RAD2DEG(psi_rad);
		dx = step_pix * cos(psi_rad);
		dy = step_pix * sin(psi_rad);

		// Truncate both ends of the helical tube
		RFLOAT trans_offset = 0.;
		x1 += ((width1 / 2.) - trans_offset) * cos(psi_rad);
		y1 += ((width1 / 2.) - trans_offset) * sin(psi_rad);
		x2 -= ((width1 / 2.) - trans_offset) * cos(psi_rad);
		y2 -= ((width1 / 2.) - trans_offset) * sin(psi_rad);

		if (!cut_into_segments)
		{
			MD_out.addObject();
	    	MD_out.setValue(EMDL_IMAGE_COORD_X, ((x1 + x2) / 2.));
	    	MD_out.setValue(EMDL_IMAGE_COORD_Y, ((y1 + y2) / 2.));
	    	MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, nr_tubes);
	    	MD_out.setValue(EMDL_ORIENT_TILT_PRIOR, 90.);
	    	MD_out.setValue(EMDL_ORIENT_PSI_PRIOR, -psi_deg);
	        MD_out.setValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, 0.);
	        MD_out.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, psi_prior_flip_ratio);

			continue;
		}

		// Calculate coordinates for all segments
		xp = x1 - (dx * 0.99);
		yp = y1 - (dy * 0.99);
		len_pix = -step_pix;
    	while (1)
    	{
    		xp += dx;
    		yp += dy;
    		len_pix += step_pix;
    		if ( ((xp > x1) && (xp > x2)) || ((xp < x1) && (xp < x2))
    				|| ((yp > y1) && (yp > y2)) || ((yp < y1) && (yp < y2)) )
    		{
    			break;
    		}
    		else
    		{
#ifdef EXCLUDE_SEGMENTS_ON_THE_EDGES
    			// Avoid segments lying on the edges of the micrographs
    			if ( (xp < half_box_size_pix) || (xp > (Xdim - half_box_size_pix))
    					|| (yp < half_box_size_pix) || (yp > (Ydim - half_box_size_pix)) )
    			{
    				// Extract from filament start-end coordinates. It is not necessary to notify the user.
    				//std::cerr << " WARNING: Particle at (" << xp << ", " << yp << ") in coordinate file " << fn_in << " is NOT extracted because it is too close to the edge." << std::flush;
    				//std::cerr << " Box_size_pix = " << box_size_pix << ", Dimensions = " << Xdim << " * " << Ydim << std::endl;
    				continue;
    			}
#endif

    			MD_out.addObject();
    	    	MD_out.setValue(EMDL_IMAGE_COORD_X, xp);
    	    	MD_out.setValue(EMDL_IMAGE_COORD_Y, yp);
    	    	MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, nr_tubes);
    	    	MD_out.setValue(EMDL_ORIENT_TILT_PRIOR, 90.);
    	    	MD_out.setValue(EMDL_ORIENT_PSI_PRIOR, -psi_deg);
    	        MD_out.setValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, pixel_size_A * len_pix);
    	        MD_out.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, psi_prior_flip_ratio);

    			nr_segments++;
    		}
    	}
		line.clear();
	}
	fin.close();

	total_segments = MD_out.numberOfObjects();
	total_tubes = nr_tubes;
	std::cout << "Input EMAN2 coordinates = " << fn_in.c_str() << ", micrograph size = " << Xdim << " * " << Ydim
			<< ", box size = " << box_size_pix << ", tubes = " << nr_tubes << ", segments = " << MD_out.numberOfObjects() << "." << std::endl;
}

/*
void divideHelicalSegmentsFromMultipleMicrographsIntoRandomHalves(
		FileName& fn_in,
		FileName& fn_out,
		int random_seed)
{
	RFLOAT ratio;
	std::string mic_name;
	int nr_swaps, nr_segments_subset1, nr_segments_subset2, helical_tube_id;
	bool divide_according_to_helical_tube_id;
	MetaDataTable MD;
	std::map<std::string, int> map_mics;
	std::map<std::string, int>::const_iterator ii_map;
	std::vector<std::pair<std::string, int> > vec_mics;

	if (fn_in == fn_out)
		REPORT_ERROR("helix.cpp::divideHelicalSegmentsFromMultipleMicrographsIntoRandomHalves(): File names error!");

	MD.clear();
	std::cout << " Loading input file..." << std::endl;
	MD.read(fn_in);
	init_random_generator(random_seed);

	if (!MD.containsLabel(EMDL_MICROGRAPH_NAME))
		REPORT_ERROR("helix.cpp::divideHelicalSegmentsFromMultipleMicrographsIntoRandomHalves(): Input MetadataTable should contain rlnMicrographName!");
	if (!MD.containsLabel(EMDL_PARTICLE_RANDOM_SUBSET))
		MD.addLabel(EMDL_PARTICLE_RANDOM_SUBSET);
	if (MD.containsLabel(EMDL_PARTICLE_HELICAL_TUBE_ID))
		divide_according_to_helical_tube_id = true;
	else
		divide_according_to_helical_tube_id = false;

	// Count micrograph names
	map_mics.clear();
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		mic_name.clear();
		MD.getValue(EMDL_MICROGRAPH_NAME, mic_name);
		if (divide_according_to_helical_tube_id)
		{
			MD.getValue(EMDL_PARTICLE_HELICAL_TUBE_ID, helical_tube_id);
			if (helical_tube_id < 1)
				REPORT_ERROR("helix.cpp::divideHelicalSegmentsFromMultipleMicrographsIntoRandomHalves(): Helical tube ID should be positive integer!");
			mic_name += std::string("_TUBEID_");
			mic_name += std::string(integerToString(helical_tube_id));
		}
		if ((map_mics.insert(std::make_pair(mic_name, 1))).second == false)
			map_mics[mic_name]++;
	}
	vec_mics.clear();
	for (ii_map = map_mics.begin(); ii_map != map_mics.end(); ii_map++)
		vec_mics.push_back(*ii_map);

	if (random_seed != 0)
	{
		// Randomise
		// 1. Randomise total number of swaps needed
		nr_swaps = ROUND(rnd_unif(vec_mics.size(), 2. * vec_mics.size()));
		// DEBUG
		if (divide_according_to_helical_tube_id)
			std::cout << " Helical tubes= " << vec_mics.size() << ", nr_swaps= " << nr_swaps << std::endl;
		else
			std::cout << " Micrographs= " << vec_mics.size() << ", nr_swaps= " << nr_swaps << std::endl;
		// 2. Perform swaps
		for (int ii = 0; ii < nr_swaps; ii++)
		{
			int ptr_a, ptr_b;
			std::pair<std::string, int> tmp;
			ptr_a = ROUND(rnd_unif(0, vec_mics.size()));
			ptr_b = ROUND(rnd_unif(0, vec_mics.size()));
			if ( (ptr_a == ptr_b) || (ptr_a < 0 ) || (ptr_b < 0) || (ptr_a >= vec_mics.size()) || (ptr_b >= vec_mics.size()) )
				continue;
			tmp = vec_mics[ptr_a];
			vec_mics[ptr_a] = vec_mics[ptr_b];
			vec_mics[ptr_b] = tmp;

			// DEBUG
			//std::cout << " Swap mic_id= " << ptr_a << " with mic_id= " << ptr_b << "." << std::endl;
		}
	}

	// Divide micrographs into halves
	map_mics.clear();
	nr_segments_subset1 = nr_segments_subset2 = 0;
	for (int ii = 0; ii < vec_mics.size(); ii++)
	{
		if (nr_segments_subset1 < nr_segments_subset2)
		{
			nr_segments_subset1 += vec_mics[ii].second;
			vec_mics[ii].second = 1;
		}
		else
		{
			nr_segments_subset2 += vec_mics[ii].second;
			vec_mics[ii].second = 2;
		}
		map_mics.insert(vec_mics[ii]);
	}

	if ( (nr_segments_subset1 < 1) || (nr_segments_subset2 < 1) )
		REPORT_ERROR("helix.cpp::divideHelicalSegmentsFromMultipleMicrographsIntoRandomHalves(): Number of helical segments from one of the two half sets is 0!");
	ratio = (RFLOAT(nr_segments_subset1) / RFLOAT(nr_segments_subset2));
	if ( (ratio > 1.5) || ( (1. / ratio) > 1.5) )
		REPORT_ERROR("helix.cpp::divideHelicalSegmentsFromMultipleMicrographsIntoRandomHalves(): Numbers of helical segments from two half sets are extremely unbalanced!");

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		mic_name.clear();
		MD.getValue(EMDL_MICROGRAPH_NAME, mic_name);
		if (divide_according_to_helical_tube_id)
		{
			MD.getValue(EMDL_PARTICLE_HELICAL_TUBE_ID, helical_tube_id);
			if (helical_tube_id < 1)
				REPORT_ERROR("helix.cpp::divideHelicalSegmentsFromMultipleMicrographsIntoRandomHalves(): Helical tube ID should be positive integer!");
			mic_name += std::string("_TUBEID_");
			mic_name += std::string(integerToString(helical_tube_id));
		}
		MD.setValue(EMDL_PARTICLE_RANDOM_SUBSET, map_mics[mic_name]);
	}

	// DEBUG
	std::cout << " Helical segments in two half sets = " << nr_segments_subset1 << ", " << nr_segments_subset2 << std::endl;
	std::cout << " Writing output file..." << std::endl;
	MD.write(fn_out);
	std::cout << " Done!" << std::endl;

	return;
}
*/

void makeHelicalReference2D(
		MultidimArray<RFLOAT>& out,
		int box_size,
		RFLOAT particle_diameter_A,
		RFLOAT tube_diameter_A,
		RFLOAT pixel_size_A,
		bool is_tube_white)
{
	RFLOAT particle_diameter_pix, tube_diameter_pix, p2, dist2, r2, t2;

	if (pixel_size_A < 0.0001)
		REPORT_ERROR("helix.cpp::makeHelicalReference2D(): Invalid pixel size!");

	particle_diameter_pix = particle_diameter_A / pixel_size_A;
	tube_diameter_pix = tube_diameter_A / pixel_size_A;

	if (box_size < 10)
		REPORT_ERROR("helix.cpp::makeHelicalReference2D(): Invalid box size!");
	if ( (particle_diameter_pix < 2.) || (particle_diameter_pix > box_size) )
		REPORT_ERROR("helix.cpp::makeHelicalReference2D(): Invalid particle diameter!");
	if ( (tube_diameter_pix < 1.) || (tube_diameter_pix > particle_diameter_pix) )
		REPORT_ERROR("helix.cpp::makeHelicalReference2D(): Invalid tube diameter!");

	r2 = (particle_diameter_pix / 2.) * (particle_diameter_pix / 2.);
	t2 = (tube_diameter_pix / 2.) * (tube_diameter_pix / 2.);

	out.clear();
	out.resize(box_size, box_size);
	out.initZeros();
	out.setXmippOrigin();

	FOR_ALL_ELEMENTS_IN_ARRAY2D(out)
	{
		dist2 = (RFLOAT)(i * i + j * j);
		p2 = (RFLOAT)(j * j);
		if ( (dist2 < r2) && (p2 < t2) )
		{
			if (is_tube_white)
				A2D_ELEM(out, i, j) = 1.;
			else
				A2D_ELEM(out, i, j) = (-1.);
		}
	}
	return;
};

/*
void makeHelicalReference3D(
		MultidimArray<RFLOAT>& out,
		int box_size,
		RFLOAT pixel_size_A,
		RFLOAT twist_deg,
		RFLOAT rise_A,
		RFLOAT tube_diameter_A,
		RFLOAT particle_diameter_A,
		int sym_Cn)
{
	RFLOAT rise_pix, tube_diameter_pix, particle_diameter_pix, particle_radius_pix;
	int particle_radius_max_pix;
	Matrix2D<RFLOAT> matrix1, matrix2;
	Matrix1D<RFLOAT> vec0, vec1, vec2;
	out.clear();

	if (box_size < 5)
		REPORT_ERROR("helix.cpp::makeHelicalReference3D(): Box size should be larger than 5!");
	if (pixel_size_A < 0.001)
		REPORT_ERROR("helix.cpp::makeHelicalReference3D(): Pixel size (in Angstroms) should be larger than 0.001!");
	if ( (fabs(twist_deg) < 0.01) || (fabs(twist_deg) > 179.99) || ((rise_A / pixel_size_A) < 0.001) )
		REPORT_ERROR("helix.cpp::makeHelicalReference3D(): Wrong helical twist or rise!");
	if (sym_Cn < 1)
		REPORT_ERROR("helix.cpp::makeHelicalReference3D(): Rotation symmetry Cn is invalid (n should be positive integer)!");

	rise_pix = rise_A / pixel_size_A;
	tube_diameter_pix = tube_diameter_A / pixel_size_A;
	particle_diameter_pix = particle_diameter_A / pixel_size_A;
	particle_radius_pix = particle_diameter_pix / 2.;
	particle_radius_max_pix = (CEIL(particle_diameter_pix / 2.)) + 1;

	if (particle_diameter_pix < 2.)
		REPORT_ERROR("helix.cpp::makeHelicalReference3D(): Particle diameter should be larger than 2 pixels!");
	if ( (tube_diameter_pix < 0.001) || (tube_diameter_pix > (RFLOAT)(box_size)) )
		REPORT_ERROR("helix.cpp::makeHelicalReference3D(): Tube diameter should be larger than 1 pixel and smaller than box size!");

	out.resize(box_size, box_size, box_size);
	out.initZeros();
	out.setXmippOrigin();

	RFLOAT x0, y0, z0;
	x0 = tube_diameter_pix / 2.;
	y0 = 0.;
	z0 = (RFLOAT)(FIRST_XMIPP_INDEX(box_size));
	vec0.clear();
	vec0.resize(2);
	XX(vec0) = x0;
	YY(vec0) = y0;
	vec1.clear();
	vec1.resize(2);
	vec2.clear();
	vec2.resize(2);

	for (int id = 0; ;id++)
	{
		RFLOAT rot1_deg, x1, y1, z1;
		rot1_deg = (RFLOAT)(id) * twist_deg;
		rotation2DMatrix(rot1_deg, matrix1, false);
		vec1 = matrix1 * vec0;

		x1 = XX(vec1);
		y1 = YY(vec1);
		z1 = z0 + (RFLOAT)(id) * rise_pix;
		if (z1 > LAST_XMIPP_INDEX(box_size))
			break;

		for (int Cn = 0; Cn < sym_Cn; Cn++)
		{
			RFLOAT rot2_deg, x2, y2, z2;
			rot2_deg = (360.) * (RFLOAT)(Cn) / (RFLOAT)(sym_Cn);
			rotation2DMatrix(rot2_deg, matrix2, false);
			vec2 = matrix2 * vec1;
			x2 = XX(vec2);
			y2 = YY(vec2);
			z2 = z1;

			for (int dx = -particle_radius_max_pix; dx <= particle_radius_max_pix; dx++)
			{
				for (int dy = -particle_radius_max_pix; dy <= particle_radius_max_pix; dy++)
				{
					for (int dz = -particle_radius_max_pix; dz <= particle_radius_max_pix; dz++)
					{
						RFLOAT _x, _y, _z, dist, val_old, val_new;
						int x3, y3, z3;

						x3 = ROUND(x2) + dx;
						y3 = ROUND(y2) + dy;
						z3 = ROUND(z2) + dz;

						if ( (x3 < FIRST_XMIPP_INDEX(box_size)) || (x3 > LAST_XMIPP_INDEX(box_size))
								|| (y3 < FIRST_XMIPP_INDEX(box_size)) || (y3 > LAST_XMIPP_INDEX(box_size))
								|| (z3 < FIRST_XMIPP_INDEX(box_size)) || (z3 > LAST_XMIPP_INDEX(box_size)) )
							continue;

						_x = (RFLOAT)(x3) - x2;
						_y = (RFLOAT)(y3) - y2;
						_z = (RFLOAT)(z3) - z2;

						dist = sqrt(_x * _x + _y * _y + _z * _z);
						if (dist > particle_radius_pix)
							continue;

						val_old = A3D_ELEM(out, z3, y3, x3);
						val_new = 0.5 + 0.5 * cos(PI * dist / particle_radius_pix);
						if (val_new > val_old)
							A3D_ELEM(out, z3, y3, x3) = val_new;
					}
				}
			}
		}
	}
	return;
}
*/

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
		int sym_Cn,
		int nr_filaments_helix_with_seam)
{
	RFLOAT rise_pix, tube_diameter_pix, particle_diameter_pix, particle_radius_pix, cyl_radius_pix, top_radius_pix, bottom_radius_pix;
	int particle_radius_max_pix;
	bool append_additional_densities = false;
	Matrix2D<RFLOAT> matrix1, matrix2;
	Matrix1D<RFLOAT> vec0, vec1, vec2;
	out.clear();

	if (box_size < 5)
		REPORT_ERROR("helix.cpp::makeHelicalReference3DWithPolarity(): Box size should be larger than 5!");
	if (pixel_size_A < 0.001)
		REPORT_ERROR("helix.cpp::makeHelicalReference3DWithPolarity(): Pixel size (in Angstroms) should be larger than 0.001!");
	if ( (fabs(twist_deg) > 179.99) || ((rise_A / pixel_size_A) < 0.001) )
		REPORT_ERROR("helix.cpp::makeHelicalReference3DWithPolarity(): Wrong helical twist or rise!");
	if (sym_Cn < 1)
		REPORT_ERROR("helix.cpp::makeHelicalReference3DWithPolarity(): Rotation symmetry Cn is invalid (n should be positive integer)!");
	if ( (topbottom_ratio < 0.) || (topbottom_ratio > 1.) )
		REPORT_ERROR("helix.cpp::makeHelicalReference3DWithPolarity(): Top-bottom width ratio should be 0~1!");
	if ( (nr_filaments_helix_with_seam > 1) && (sym_Cn != 1) )
		REPORT_ERROR("helix.cpp::makeHelicalReference3DWithPolarity(): Set Cn point group symmetry to 1 for a helix with seam!");
	if ( (nr_filaments_helix_with_seam > 1) && (!(topbottom_ratio > 0.9999)) )
		REPORT_ERROR("helix.cpp::makeHelicalReference3DWithPolarity(): Set top-bottom width ratio to 1 for a helix with seam!");

	rise_pix = rise_A / pixel_size_A;
	tube_diameter_pix = tube_diameter_A / pixel_size_A;
	particle_diameter_pix = particle_diameter_A / pixel_size_A;
	particle_radius_pix = particle_diameter_pix / 2.;
	particle_radius_max_pix = (CEIL(particle_diameter_pix / 2.)) + 1;
	top_radius_pix = cyl_radius_pix = 0.5 * cyl_diameter_A / pixel_size_A;
	bottom_radius_pix = top_radius_pix * topbottom_ratio;

	if (particle_diameter_pix < 2.)
		REPORT_ERROR("helix.cpp::makeHelicalReference3DWithPolarity(): Particle diameter should be larger than 2 pixels!");
	if ( (tube_diameter_pix < 0.001) || (tube_diameter_pix > (RFLOAT)(box_size)) )
		REPORT_ERROR("helix.cpp::makeHelicalReference3DWithPolarity(): Tube diameter should be larger than 1 pixel and smaller than box size!");
	if ( (cyl_radius_pix < 1.) || (cyl_radius_pix > particle_radius_pix) )
		REPORT_ERROR("helix.cpp::makeHelicalReference3DWithPolarity(): Cylindrical diameter should be > 1 pixel and < particle diameter!");

	out.resize(box_size, box_size, box_size);
	out.initZeros();
	out.setXmippOrigin();

	RFLOAT x0, y0, z0;
	// OLD
	//x0 = tube_diameter_pix / 2.;
	//y0 = 0.;
	// NEW - To generate references with Dn symmetry. TODO: Am I doing what I want?
	z0 = rise_pix * FLOOR((RFLOAT)(FIRST_XMIPP_INDEX(box_size)) / rise_pix);
	x0 = (tube_diameter_pix / 2.) * cos( (PI * twist_deg * z0) / (rise_pix * 180.) );
	y0 = (tube_diameter_pix / 2.) * sin( (PI * twist_deg * z0) / (rise_pix * 180.) );
	vec0.clear();
	vec0.resize(2);
	XX(vec0) = x0;
	YY(vec0) = y0;
	vec1.clear();
	vec1.resize(2);
	vec2.clear();
	vec2.resize(2);

	append_additional_densities = false;
	for (int id = 0; ;id++)
	{
		RFLOAT rot1_deg, x1, y1, z1;
		rot1_deg = (RFLOAT)(id) * twist_deg;
		rotation2DMatrix(rot1_deg, matrix1, false);
		vec1 = matrix1 * vec0;

		x1 = XX(vec1);
		y1 = YY(vec1);
		z1 = z0 + (RFLOAT)(id) * rise_pix;
		if (z1 > LAST_XMIPP_INDEX(box_size))
			break;

		for (int Cn = 0; Cn < sym_Cn; Cn++)
		{
			RFLOAT rot2_deg, x2, y2, z2;
			rot2_deg = (360.) * (RFLOAT)(Cn) / (RFLOAT)(sym_Cn);
			rotation2DMatrix(rot2_deg, matrix2, false);
			vec2 = matrix2 * vec1;
			x2 = XX(vec2);
			y2 = YY(vec2);
			z2 = z1;

			for (int dz = -particle_radius_max_pix; dz <= particle_radius_max_pix; dz++)
			{
				RFLOAT thres_xy = (top_radius_pix - bottom_radius_pix) * 0.5 * dz / particle_radius_pix + (top_radius_pix + bottom_radius_pix) / 2.;
				for (int dy = -particle_radius_max_pix; dy <= particle_radius_max_pix; dy++)
				{
					for (int dx = -particle_radius_max_pix; dx <= particle_radius_max_pix; dx++)
					{
						RFLOAT _x, _y, _z, dist, val_old, val_new;
						int x3, y3, z3;

						x3 = ROUND(x2) + dx;
						y3 = ROUND(y2) + dy;
						z3 = ROUND(z2) + dz;

						if ( (x3 < FIRST_XMIPP_INDEX(box_size)) || (x3 > LAST_XMIPP_INDEX(box_size))
								|| (y3 < FIRST_XMIPP_INDEX(box_size)) || (y3 > LAST_XMIPP_INDEX(box_size))
								|| (z3 < FIRST_XMIPP_INDEX(box_size)) || (z3 > LAST_XMIPP_INDEX(box_size)) )
							continue;

						_x = (RFLOAT)(x3) - x2;
						_y = (RFLOAT)(y3) - y2;
						_z = (RFLOAT)(z3) - z2;

						dist = sqrt(_x * _x + _y * _y + _z * _z);
						if (dist > particle_radius_pix)
							continue;

						val_old = A3D_ELEM(out, z3, y3, x3);
						val_new = 0.;

						// Draw the shape you want!
						if (topbottom_ratio > 0.9999) // Without polarity. Thus spheres.
						{
							val_new = 0.5 + 0.5 * cos(PI * dist / particle_radius_pix);
							if (val_new > val_old)
								A3D_ELEM(out, z3, y3, x3) = val_new;
						}
						else // With polarity
						{
							dist = sqrt(_x * _x + _y * _y);
							if (dist < thres_xy)
							{
								val_new = 0.5 + 0.5 * cos(PI * dist / thres_xy);
								val_new *= 0.5 + 0.5 * cos(PI * 0.5 * _z / particle_radius_pix); // something arbitrary
								if (val_new > val_old)
									A3D_ELEM(out, z3, y3, x3) = val_new;
							}
						}
					}
				}
			}
		}

		if (nr_filaments_helix_with_seam > 1)
		{
			if (id % nr_filaments_helix_with_seam == 0)
				append_additional_densities = (append_additional_densities == true) ? (false) : (true);

			if (append_additional_densities)
			{
				x1 *= (tube_diameter_pix + particle_diameter_pix) / tube_diameter_pix;
				y1 *= (tube_diameter_pix + particle_diameter_pix) / tube_diameter_pix;
				z1 += particle_diameter_pix / 2.;

				for (int dz = -particle_radius_max_pix / 2.; dz <= particle_radius_max_pix / 2.; dz++)
				{
					for (int dy = -particle_radius_max_pix / 2.; dy <= particle_radius_max_pix / 2.; dy++)
					{
						for (int dx = -particle_radius_max_pix / 2.; dx <= particle_radius_max_pix / 2.; dx++)
						{
							RFLOAT _x, _y, _z, dist, val_old, val_new;
							int x2, y2, z2;

							x2 = ROUND(x1 + dx);
							y2 = ROUND(y1 + dy);
							z2 = ROUND(z1 + dz);

							if ( (x2 < FIRST_XMIPP_INDEX(box_size)) || (x2 > LAST_XMIPP_INDEX(box_size))
									|| (y2 < FIRST_XMIPP_INDEX(box_size)) || (y2 > LAST_XMIPP_INDEX(box_size))
									|| (z2 < FIRST_XMIPP_INDEX(box_size)) || (z2 > LAST_XMIPP_INDEX(box_size)) )
								continue;

							_x = (RFLOAT)(x2) - x1;
							_y = (RFLOAT)(y2) - y1;
							_z = (RFLOAT)(z2) - z1;

							dist = sqrt(_x * _x + _y * _y + _z * _z);
							if (dist > (particle_radius_pix / 2.))
								continue;

							val_old = A3D_ELEM(out, z2, y2, x2);
							val_new = 0.;

							val_new = 0.5 + 0.5 * cos(2. * PI * dist / particle_radius_pix);
							if (val_new > val_old)
								A3D_ELEM(out, z2, y2, x2) = val_new;
						}
					}
				}// End of looping over x,y,z
			}
		}// End of nr_filaments_helix_with_seam > 1
	}
	return;
}

void divideStarFile(
		FileName& fn_in,
		int nr)
{
	FileName fn_out;
	MetaDataTable MD_in, MD;
	int total_lines, nr_lines, line_id, file_id;

	if (fn_in.getExtension() != "star")
		REPORT_ERROR("helix.cpp::divideStarFile: Input file should be in .star format!");
	if ( (nr < 2) || (nr > 999999) )
		REPORT_ERROR("helix.cpp::divideStarFile: The number of output files should be within range 2~999999!");

	std::cout << " Loading input file: " << fn_in << " ..." << std::endl;
	MD_in.clear();
	MD_in.read(fn_in);
	std::cout << " Input file loaded." << std::endl;

	total_lines = MD_in.numberOfObjects();
	if (total_lines < nr)
		REPORT_ERROR("helix.cpp::divideStarFile: The number of total objects is smaller than the number of output files!");
	nr_lines = total_lines / nr;
	std::cout << " Total objects = " << total_lines << ", number of output files = " << nr << std::endl;
	std::cout << " Writing output files..." << std::endl;

	line_id = 0;
	file_id = 0;
	MD.clear();
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_in)
	{
		line_id++;
		MD.addObject(MD_in.getObject());
		if ( ((line_id % nr_lines) == 0) && ((total_lines - line_id) > nr_lines) )
		{
			file_id++;
			fn_out = fn_in.withoutExtension() + "_sub" + integerToString(file_id, 6, '0');
			fn_out = fn_out.addExtension("star");
			MD.write(fn_out);
			MD.clear();
		}
	}
	if (MD.numberOfObjects() != 0)
	{
		file_id++;
		fn_out = fn_in.withoutExtension() + "_sub" + integerToString(file_id, 6, '0');
		fn_out = fn_out.addExtension("star");
		MD.write(fn_out);
		MD.clear();
	}
	std::cout << " Done!" << std::endl;
	return;
}

void mergeStarFiles(FileName& fn_in)
{
	int file_id;
	std::vector<FileName> fns_list;
	MetaDataTable MD_combined, MD_in;
	FileName fn_root, fn_out;

	fns_list.clear();
	fn_root = "*" + fn_in + "*";

	if (fn_root.globFiles(fns_list) < 2)
		REPORT_ERROR("helix.cpp::combineStarFiles: Only 0 or 1 input file! There is no need to combine!");
	for (file_id = 0; file_id < fns_list.size(); file_id++)
	{
		if (fns_list[file_id].getExtension() != "star")
			REPORT_ERROR("helix.cpp::combineStarFiles: Input files should have STAR extension!");
	}

	std::cout << " Combining STAR files: " << std::endl;
	std::cout << " BEWARE: All STAR files should contain the same header!" << std::endl;
	for (file_id = 0; file_id < fns_list.size(); file_id++)
		std::cout << "  " << fns_list[file_id] << std::endl;
	std::cout << " Loading input files..." << std::endl;

	MD_combined.clear();
	for (file_id = 0; file_id < fns_list.size(); file_id++)
	{
		MD_in.clear();
		MD_in.read(fns_list[file_id]);
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_in)
		{
			MD_combined.addObject(MD_in.getObject());
		}
		std::cout << "  " << MD_combined.numberOfObjects() << " objects loaded." << std::endl;
	}

	std::cout << " Writing the combined output file..." << std::endl;
	fn_out = fn_in + "_combined.star";
	MD_combined.write(fn_out);
	std::cout << " Done!" << std::endl;

	return;
}

void sortHelicalTubeID(MetaDataTable& MD)
{
	std::string str_particle_fullname, str_particle_name, str_comment, str_particle_id;
	int int_tube_id, nr_tubes;
	std::vector<HelicalSegmentPriorInfoEntry> list;
	std::set<std::string> tubes;
	long int MDobjectID;

	if ( (!MD.containsLabel(EMDL_IMAGE_NAME))
			|| (!MD.containsLabel(EMDL_ORIENT_TILT))
			|| (!MD.containsLabel(EMDL_ORIENT_PSI))
			|| (!MD.containsLabel(EMDL_ORIENT_TILT_PRIOR))
			|| (!MD.containsLabel(EMDL_ORIENT_PSI_PRIOR))
			|| (!MD.containsLabel(EMDL_PARTICLE_HELICAL_TUBE_ID))
			|| (!MD.containsLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM))
			|| (!MD.containsLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO)) )
		REPORT_ERROR("helix.cpp::sortHelicalTubeID: Labels of helical prior information are missing!");

	int_tube_id = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		MD.getValue(EMDL_IMAGE_NAME, str_particle_fullname);
		MD.getValue(EMDL_PARTICLE_HELICAL_TUBE_ID, int_tube_id);
		str_comment = str_particle_name + "@TUBE@" + integerToString(int_tube_id, 6);
		tubes.insert(str_comment);

		str_particle_name = str_particle_fullname.substr(str_particle_fullname.find("@") + 1);
		str_particle_id = str_particle_fullname.substr(0, str_particle_fullname.find("@"));
		str_comment = str_particle_name + "@TUBE@" + integerToString(int_tube_id, 6) + "@PARTICLE@" + str_particle_id;

		// DEBUG
		//std::cout << str_comment << std::endl;

		MD.setValue(EMDL_IMAGE_NAME, str_comment);
	}
	MD.newSort(EMDL_IMAGE_NAME);
	nr_tubes = tubes.size();
	tubes.clear();

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		MD.getValue(EMDL_IMAGE_NAME, str_comment);
		str_particle_id = str_comment.substr(str_comment.find("@PARTICLE@") + 10);
		str_particle_name = str_comment.substr(0, str_comment.find("@TUBE@"));
		str_comment = str_particle_id + "@" + str_particle_name;
		MD.setValue(EMDL_IMAGE_NAME, str_comment);
	}

	std::vector<RFLOAT> dummy;
	updatePriorsForHelicalReconstruction(
			MD, 1.,dummy, dummy, 1,
			false, false,
			0., 0., 0., 1., false, 1);

	list.clear();

	return;
}

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
		int random_seed)
{
	Image<RFLOAT> img;
	int nr_segments, tube_id;
	RFLOAT rot, psi, tilt, new_psi, new_tilt, xoff, yoff, zoff, new_xoff, new_yoff, new_zoff, step_pix, psi_flip_ratio, len_pix;
	MetaDataTable MD;
	FileName fn_mic, fn_star_out_full, fn_star_out_priors, fn_star_out_wopriors, fn_ext;
	std::ofstream fout;
	std::string command;

	if ( (fn_vol_in.contains("temporary")) || (fn_star_out.contains("temporary")) )
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Avoid 'temporary' in the input and output file names!");
	if ( (is_3d_tomo) && (fn_star_out.contains("subtomo")) )
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Avoid 'subtomo' in the input and output file names!");

	if (fn_vol_in.getExtension() != "mrc")
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Input 3D volume should be in .mrc format!");
	img.read(fn_vol_in, false); // Read the header only!
	if ( (img().getDim() != 3) || (XSIZE(img()) != YSIZE(img())) || (YSIZE(img()) != ZSIZE(img())) || (XSIZE(img()) < 10) )
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Input volume should be a 3D cube (>10*10*10)!");
	if ( (new_box_size > XSIZE(img())) || (new_box_size < 10) || (new_box_size % 2) )
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Cropped box size (an even number) should be at least 10, and smaller than the box size of the input 3D volume!");
	if (angpix < 0.001)
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Pixel size should be larger than 0.001 Angstroms!");
	if ( (rise_A < 0.001) || ((rise_A / angpix) < 0.001) || ((rise_A / angpix) > (0.3333 * new_box_size)) )
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Rise is smaller than 0.001 pixels or larger than 1/3 of the new (cropped) box size!");
	if ( (fabs(twist_deg) < 0.001) || (fabs(twist_deg) > 180.))
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Error in helical twist!");

	// TODO: raise error if nr_asu<0 or too big, n too small!
	if ( (nr_tubes < 2) || (nr_subunits < 10) || (nr_asu < 1) || (((nr_subunits / nr_asu) / nr_tubes) < 3) || ((nr_subunits / nr_asu) > 999000) || ((nr_asu * rise_A / angpix) > new_box_size) )
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Errors in the number of tubes, asymmetrical units or total subunits!");
	if ( (sigma_psi > 10.) || (sigma_tilt > 10.) )
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: sigma_psi and sigma_tilt should not be larger than 10 degrees.");
	if (sigma_offset > 50.)
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: sigma_trans should not be larger than 50 pixels.");

	if (fn_star_out.getExtension() != "star")
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Output file should be in STAR format!");
	fn_star_out_full = fn_star_out.withoutExtension() + (std::string)("_info.star");
	fn_star_out_priors = fn_star_out.withoutExtension() + (std::string)("_helical_priors.star");
	fn_star_out_wopriors = fn_star_out.withoutExtension() + (std::string)("_no_priors.star");

	if (is_3d_tomo)
		fout.open("simulate-3d-subtomos.sh", std::ios::out);
	else
		fout.open("simulate-2d-segments.sh", std::ios::out);
    if (!fout)
        REPORT_ERROR( (std::string)"helix.cpp::simulateHelicalSegments: Cannot write to .sh script!");

	sigma_tilt = (sigma_tilt > 0.) ? (sigma_tilt) : (0.);
	sigma_psi = (sigma_psi > 0.) ? (sigma_psi) : (0.);
	sigma_offset = (sigma_offset > 0.) ? (sigma_offset) : (0.);

	if (random_seed < 0)
		random_seed = time(NULL);
    init_random_generator(random_seed);

	nr_segments = nr_subunits / nr_asu;

	MD.clear();
    MD.addLabel(EMDL_ORIENT_ROT);
    MD.addLabel(EMDL_ORIENT_TILT);
    MD.addLabel(EMDL_ORIENT_PSI);
    MD.addLabel(EMDL_ORIENT_ORIGIN_X_ANGSTROM);
    MD.addLabel(EMDL_ORIENT_ORIGIN_Y_ANGSTROM);
    if (is_3d_tomo)
    	MD.addLabel(EMDL_ORIENT_ORIGIN_Z_ANGSTROM);
    MD.addLabel(EMDL_ORIENT_ROT_PRIOR);
    MD.addLabel(EMDL_ORIENT_TILT_PRIOR);
    MD.addLabel(EMDL_ORIENT_PSI_PRIOR);
    MD.addLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM);
    MD.addLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO);
    MD.addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);
    MD.addLabel(EMDL_IMAGE_NAME);
    MD.addLabel(EMDL_MICROGRAPH_NAME);
    MD.addLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE);
    MD.addLabel(EMDL_CTF_MAGNIFICATION);

    tube_id = 0;
    step_pix = nr_asu * rise_A / angpix;
    psi_flip_ratio = (do_bimodal_searches) ? (0.5) : (0.);
	for (int id = 0; id < nr_segments; id++)
	{
		if ( ( (id % (nr_segments / nr_tubes)) == 0 ) && ( (nr_segments - id) >= (nr_segments / nr_tubes) ) )
		{
			tube_id++;
			len_pix = -step_pix;

			if (is_3d_tomo)
				tilt = rnd_unif(0.01, 179.99); // If realWRAP function works well, set this to 0-180.
			else
				tilt = rnd_unif(85., 95.);
			rot = rnd_unif(0.01, 359.99);
			psi = rnd_unif(-179.99, 179.99);
		}

		len_pix += step_pix;
		rot += twist_deg * ((RFLOAT)(nr_asu));
		rot = realWRAP(rot, -180., 180.); // Does this realWRAP function work well? No...
		rot = (rot < -180.) ? (rot + 360.) : (rot);
		rot = (rot > 180.) ? (rot - 360.) : (rot);

		if (sigma_tilt < 0.0001)
			new_tilt = tilt;
		else
			new_tilt = tilt + rnd_gaus(0., sigma_tilt);
		new_tilt = (new_tilt < 0.) ? (-new_tilt) : (new_tilt); // Do NOT change the polarities
		new_tilt = (new_tilt > 180.) ? (360. - new_tilt) : (new_tilt);

		if (sigma_psi < 0.0001)
			new_psi = psi;
		else
			new_psi = psi + rnd_gaus(0., sigma_psi);
		new_psi = (new_psi < -180.) ? (new_psi + 360.) : (new_psi);
		new_psi = (new_psi > 180.) ? (new_psi - 360.) : (new_psi);

		xoff = yoff = zoff = new_xoff = new_yoff = new_zoff = 0.;
		if (sigma_offset > 0.0001)
		{
			new_xoff = (is_3d_tomo) ? (rnd_gaus(0., sigma_offset)) : (rnd_unif(-0.5 * rise_A, 0.5 * rise_A));
			new_yoff = rnd_gaus(0., sigma_offset);
			new_zoff = (is_3d_tomo) ? (rnd_unif(-0.5 * rise_A, 0.5 * rise_A)) : (0.);
			transformCartesianAndHelicalCoords(
					new_xoff, new_yoff, new_zoff,
					xoff, yoff, zoff,
					rot, new_tilt, new_psi,
					(is_3d_tomo) ? (3) : (2),
					HELICAL_TO_CART_COORDS);
		}

		MD.addObject();
    	MD.setValue(EMDL_ORIENT_ROT, rot);
    	MD.setValue(EMDL_ORIENT_TILT, new_tilt);
    	MD.setValue(EMDL_ORIENT_PSI, new_psi);
    	MD.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff);
    	MD.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff);
    	if (is_3d_tomo)
    		MD.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, zoff);
    	MD.setValue(EMDL_ORIENT_ROT_PRIOR, rot);
    	MD.setValue(EMDL_ORIENT_TILT_PRIOR, new_tilt);
    	MD.setValue(EMDL_ORIENT_PSI_PRIOR, new_psi);
    	MD.setValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, angpix * len_pix);
    	MD.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, psi_flip_ratio);
    	MD.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, tube_id);
    	if (is_3d_tomo)
    	{
    		fn_mic.compose((id + 1), "dummy");
    		fn_mic = fn_mic.beforeFirstOf("@");
    		fn_mic = (std::string)("subtomo-3d-") + (std::string)(fn_mic) + (std::string)(".mrc");
    	}
    	else
    		fn_mic.compose((id + 1), "segments-2d.mrcs");
    	MD.setValue(EMDL_IMAGE_NAME, fn_mic);
    	if (is_3d_tomo)
    		MD.setValue(EMDL_MICROGRAPH_NAME, (std::string)("tomogram-01.mrc"));
    	else
    		MD.setValue(EMDL_MICROGRAPH_NAME, (std::string)("micrograph-01.mrc"));
    	MD.setValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, 14.0);
    	MD.setValue(EMDL_CTF_MAGNIFICATION, (140000. / angpix));

    	// Generate .sh script
    	if (is_3d_tomo)
    	{
    		fout << "relion_project --i " << fn_vol_in << " --o temporary-image-3d.mrc --angpix " << angpix << " --rot "
    				<< rot << " --tilt " << new_tilt << " --psi " << new_psi << " --xoff " << xoff << " --yoff "
    				<< yoff << " --zoff " << zoff << " --3d_rot" << std::flush;
    		if (white_noise > 0.)
    			fout << " --add_noise --white_noise " << white_noise;
    		fout << std::endl;

    		fout << "relion_image_handler --i temporary-image-3d.mrc --o " << fn_mic << " --new_box " << new_box_size << std::endl;
    		//fout << "rm -rf temporary-image-3d.mrc" << std::endl;
    	}
	}
	if (is_3d_tomo)
	{
		fout << "rm -rf temporary-image-3d.mrc" << std::endl;
		fout << "relion_helix_toolbox --norm --i " << fn_star_out_full << " --o_root _norm --angpix " << angpix << " --cyl_outer_diameter " << cyl_outer_diameter_A << std::endl;
		fout << "rm -rf subtomo-3d-??????.mrc" << std::endl;
	}
	MD.write(fn_star_out_full);

	if (!is_3d_tomo)
	{
		fout << "relion_project --i " << fn_vol_in << " --o temporary-images-2d --angpix " << angpix << " --ang " << fn_star_out_full << std::flush;
		if (white_noise > 0.)
			fout << " --add_noise --white_noise " << white_noise;
		fout << std::endl;

		fout << "relion_image_handler --i temporary-images-2d.mrcs --o tt --new_box " << new_box_size << std::endl;
		fout << "rm -rf temporary-images-2d.mrcs temporary-images-2d.star" << std::endl;
		fout << "mv temporary-images-2d_tt.mrcs segments-2d.mrcs" << std::endl;
		fout << "relion_helix_toolbox --norm --i " << fn_star_out_full << " --o_root _norm --angpix " << angpix << " --cyl_outer_diameter " << cyl_outer_diameter_A << std::endl;
		fout << "rm -rf segments-2d.mrcs" << std::endl;
	}

	fout << "mv " << fn_star_out_full << " " << (fn_star_out_full.withoutExtension() + (std::string)(".txt")) << std::endl;
	fout.close();

	MD.deactivateLabel(EMDL_ORIENT_ROT_PRIOR);
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
    	MD.setValue(EMDL_ORIENT_ROT, 0.);
    	if (!is_3d_tomo)
    	{
    		MD.setValue(EMDL_ORIENT_TILT, 90.);
    		MD.setValue(EMDL_ORIENT_TILT_PRIOR, 90.);
    	}
    	MD.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, 0.);
    	MD.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, 0.);
    	if (is_3d_tomo)
    		MD.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, 0.);

		MD.getValue(EMDL_IMAGE_NAME, fn_mic);
		fn_ext = fn_mic.getExtension();
		fn_mic = fn_mic.withoutExtension() + "_norm." + fn_ext;
		MD.setValue(EMDL_IMAGE_NAME, fn_mic);
	}
	MD.write(fn_star_out_priors);

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
    	MD.setValue(EMDL_ORIENT_ROT, 0.);
    	MD.setValue(EMDL_ORIENT_TILT, 0.);
    	MD.setValue(EMDL_ORIENT_PSI, 0.);
	}
	MD.deactivateLabel(EMDL_ORIENT_TILT_PRIOR);
	MD.deactivateLabel(EMDL_ORIENT_PSI_PRIOR);
	MD.deactivateLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM);
	MD.deactivateLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO);
	MD.deactivateLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);
	MD.write(fn_star_out_wopriors);

	if (is_3d_tomo)
		command = "chmod u+x simulate-3d-subtomos.sh";
	else
		command = "chmod u+x simulate-2d-segments.sh";
	int res = system(command.c_str());

	return;
};

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
		std::ostream& out)
{
	if (iter < 1)
		REPORT_ERROR("helix.cpp::outputHelicalSymmetryStatus(): BUG iteration id cannot be less than 1!");
	if ( (do_local_search_helical_symmetry) && (iter > 1) )
	{
		out << " Local searches of helical twist from " << twist_min_deg << " to " << twist_max_deg << " degrees, rise from " << rise_min_A << " to " << rise_max_A << " Angstroms." << std::endl;
	}
	else
	{
		out << " For all classes, helical twist = " << twist_initial_deg << " degrees, rise = " << rise_initial_A << " Angstroms." << std::endl;
		return;
	}

	if (do_split_random_halves)
	{
		RFLOAT twist_avg_deg = (twist_deg_half1 + twist_deg_half2) / 2.;
		RFLOAT rise_avg_A = (rise_A_half1 + rise_A_half2) / 2.;

		// TODO: raise a warning if two sets of helical parameters are >1% apart?
		out << " (Half 1) Refined helical twist = " << twist_deg_half1 << " degrees, rise = " << rise_A_half1 << " Angstroms." << std::endl;
		out << " (Half 2) Refined helical twist = " << twist_deg_half2 << " degrees, rise = " << rise_A_half2 << " Angstroms." << std::endl;
		out << " Averaged helical twist = " << twist_avg_deg << " degrees, rise = " << rise_avg_A << " Angstroms." << std::endl;
		return;
	}
	else
	{
		if ( (rise_A.size() != twist_deg.size()) || (rise_A.size() < 1) )
			REPORT_ERROR("helix.cpp::outputHelicalSymmetryStatus(): BUG vectors rise_A and twist_deg are not of the same size!");
		for (int iclass = 0; iclass < rise_A.size(); iclass++)
		{
			out << " (Class " << (iclass + 1) << ") Refined helical twist = " << twist_deg[iclass] << " degrees, rise = " << rise_A[iclass] << " Angstroms." << std::endl;
		}
	}
}

void excludeLowCTFCCMicrographs(
		FileName& fn_in,
		FileName& fn_out,
		RFLOAT cc_min,
		RFLOAT EPA_lowest_res,
		RFLOAT df_min,
		RFLOAT df_max)
{
	bool contain_EPA_res;
	MetaDataTable MD_in, MD_out;
	int nr_mics_old, nr_mics_new;
	RFLOAT cc, EPA_res, dU, dV;
	if ( (fn_in.getFileFormat() != "star") || (fn_out.getFileFormat() != "star") )
		REPORT_ERROR("helix.cpp::excludeLowCTFCCMicrographs(): MetaDataTable should have .star extension.");
	if (fn_in == fn_out)
		REPORT_ERROR("helix.cpp::excludeLowCTFCCMicrographs(): File names must be different.");
	if (df_min > df_max)
		REPORT_ERROR("helix.cpp::excludeLowCTFCCMicrographs(): Minimum defocus threshold should be smaller the maximum.");

	MD_in.clear();
	MD_in.read(fn_in);
	if ( (!MD_in.containsLabel(EMDL_CTF_DEFOCUSU))
			|| (!MD_in.containsLabel(EMDL_CTF_DEFOCUSV))
			|| (!MD_in.containsLabel(EMDL_CTF_DEFOCUS_ANGLE))
			|| (!MD_in.containsLabel(EMDL_CTF_VOLTAGE))
			|| (!MD_in.containsLabel(EMDL_CTF_CS))
			|| (!MD_in.containsLabel(EMDL_CTF_Q0))
			|| (!MD_in.containsLabel(EMDL_CTF_MAGNIFICATION))
			|| (!MD_in.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
			|| (!MD_in.containsLabel(EMDL_CTF_FOM)) )
		REPORT_ERROR("helix.cpp::excludeLowCTFCCMicrographs(): Input STAR file should contain CTF information.");

	contain_EPA_res = MD_in.containsLabel(EMDL_CTF_MAXRES);

	nr_mics_old = nr_mics_new = 0;
	MD_out.clear();
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_in)
	{
		nr_mics_old++;
		MD_in.getValue(EMDL_CTF_FOM, cc);
		MD_in.getValue(EMDL_CTF_MAXRES, EPA_res);
		MD_in.getValue(EMDL_CTF_DEFOCUSU, dU);
		MD_in.getValue(EMDL_CTF_DEFOCUSV, dV);
		if ( (cc > cc_min) && (dU > df_min) && (dU < df_max) && (dV > df_min) && (dV < df_max) )
		{
			if ( (contain_EPA_res) && (EPA_res > EPA_lowest_res) )
			{}
			else
			{
				nr_mics_new++;
				MD_out.addObject(MD_in.getObject());
			}
		}
	}

	std::cout << " Number of micrographs (input / output) = " << nr_mics_old << " / " << nr_mics_new << std::endl;
	if (MD_out.numberOfObjects() < 1)
		std::cout << " No micrographs in output file!" << std::endl;

	MD_out.write(fn_out);
	return;
}

void cutOutPartOfHelix(
		const MultidimArray<RFLOAT>& vin,
		MultidimArray<RFLOAT>& vout,
		long int new_boxdim,
		RFLOAT ang_deg,
		RFLOAT z_percentage)
{
	long int Xdim, Ydim, Zdim, Ndim, old_boxdim;

	vout.clear();

	if (vin.getDim() != 3)
		REPORT_ERROR("helix.cpp::cutOutPartOfHelix(): Input image is not 3D!");
	if (!(z_percentage > 0.))
		REPORT_ERROR("helix.cpp::cutOutPartOfHelix(): Z length must be larger than 0!");
	if (!(ang_deg > 0.))
		REPORT_ERROR("helix.cpp::cutOutPartOfHelix(): Angular range must be larger than 0!");
	ang_deg = (ang_deg > 91.) ? (91.) : (ang_deg);

	vin.getDimensions(Xdim, Ydim, Zdim, Ndim);
	old_boxdim = (Xdim < Ydim) ? (Xdim) : (Ydim);
	old_boxdim = (Zdim < old_boxdim) ? (Zdim) : (old_boxdim);
	if ( (new_boxdim <= 0) || (new_boxdim > (old_boxdim / 2)) )
		new_boxdim = old_boxdim / 2;

	vout.initZeros(new_boxdim, new_boxdim, new_boxdim);

	// Fill in values
	long int old_ymax = YSIZE(vin) + FIRST_XMIPP_INDEX(YSIZE(vin));
	long int old_xmax = XSIZE(vin) + FIRST_XMIPP_INDEX(XSIZE(vin));
	long int old_x0 = FIRST_XMIPP_INDEX(XSIZE(vin));
	long int old_y0 = FIRST_XMIPP_INDEX(YSIZE(vin));
	long int old_z0 = FIRST_XMIPP_INDEX(ZSIZE(vin));
	long int new_z0 = FIRST_XMIPP_INDEX(ZSIZE(vout));
    for (long int zi = 0; zi < ZSIZE(vout); zi++)
    {
    	// Z subscript is out of range
    	if ( ((RFLOAT)(ABS(zi + new_z0)) / (RFLOAT)(ZSIZE(vin))) > (z_percentage / 2.) )
    		continue;

    	// Loop over X and Y
        for (long int yi = 0; yi < YSIZE(vout); yi++)
        {
        	for (long int xi = 0; xi < XSIZE(vout); xi++)
        	{
        		RFLOAT deg = (180.) * atan2((double)(yi), (double)(xi)) / PI;

        		// X or Y subscripts is out of range
        		if ( (ang_deg < 90.) && ( (deg < ((45.) - (ang_deg / 2.))) || (deg > ((45.) + (ang_deg / 2.))) ) )
        			continue;
        		if ( (yi >= old_ymax) || (xi >= old_xmax) )
        			continue;

        		// Fill in voxels
        		DIRECT_A3D_ELEM(vout, zi, yi, xi) = DIRECT_A3D_ELEM(vin, zi + new_z0 - old_z0, yi - old_y0, xi - old_x0);
        	}
        }
    }
    vout.setXmippOrigin();
}

void HelicalSegmentPriorInfoEntry::clear()
{
	helical_tube_name.clear();
	MDobjectID = -1;
	rot_deg = psi_deg = tilt_deg = 0.;
	dx_A = dy_A = dz_A = 0.;
	track_pos_A = 0.;
	has_wrong_polarity = false;
	subset = classID = 0;

	rot_prior_deg = psi_prior_deg = tilt_prior_deg = 0.;  // KThurber
	dx_prior_A = dy_prior_A = dz_prior_A = 0.;
	psi_flip_ratio = 0.;
	psi_prior_flip = false; // KThurber

};

bool HelicalSegmentPriorInfoEntry::operator<(const HelicalSegmentPriorInfoEntry &rhs) const
{
	if ( (helical_tube_name.length() == 0) || (rhs.helical_tube_name.length() == 0) )
	{
		std::cerr << "Compare # " << MDobjectID << " with # " << rhs.MDobjectID << std::endl;
		REPORT_ERROR("helix.h::HelicalSegmentPriorInfoEntry::operator<(): Name string of helical segments are empty!");
	}

	if (helical_tube_name != rhs.helical_tube_name)
		return (helical_tube_name < rhs.helical_tube_name);

	if (fabs(track_pos_A - rhs.track_pos_A) < (1e-5))
	{
		std::cerr << "Compare # " << MDobjectID << " with # " << rhs.MDobjectID << std::endl;
		REPORT_ERROR("helix.h::HelicalSegmentPriorInfoEntry::operator<(): A pair of same helical segments is found!");
	}

	return (track_pos_A < rhs.track_pos_A);
};

void HelicalSegmentPriorInfoEntry::checkPsiPolarity()
{
	RFLOAT diff_psi = ABS(psi_deg - psi_prior_deg);
	has_wrong_polarity = false;
	if (diff_psi > 180.)
		diff_psi = ABS(diff_psi - 360.);
	if (diff_psi > 90.)
		has_wrong_polarity = true;

};

// KThurber add this entire function
void flipPsiTiltForHelicalSegment(
		RFLOAT old_psi,
		RFLOAT old_tilt,
		RFLOAT& new_psi,
		RFLOAT& new_tilt)
{
	new_psi = (old_psi < 0.) ? (old_psi + 180.) : (old_psi - 180.);
	new_tilt = 180. - old_tilt;
}

//#define DEBUG_HELICAL_UPDATE_ANGULAR_PRIORS
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
		RFLOAT sigma_cutoff)
{
	RFLOAT range_rot, range_tilt, range_psi, range2_offset, psi_flip_ratio;
	std::string str_name;
	int nr_same_polarity, nr_opposite_polarity, subset, data_dim;
	bool do_avg, unimodal_angular_priors;

	// Check subscript
	if ( (list.size() < 1) || (sid < 0) || (eid >= list.size()) || (sid > eid) )
		REPORT_ERROR("helix.cpp::updatePriorsForOneHelicalTube(): Subscripts are invalid!");

	// Init
	data_dim = (is_3D_data) ? (3) : (2);
	// TODO: test: Do not do local averaging if data_dim == 3
	do_avg = (!is_3D_data) && (sigma_segment_dist > 0.01) && (list.size() > 1); // Do local average of orientations and translations or just flip tilt and psi angles?
	sigma2_rot = (sigma2_rot > 0.) ? (sigma2_rot) : (0.);  // KThurber
	sigma2_tilt = (sigma2_tilt > 0.) ? (sigma2_tilt) : (0.);
	sigma2_psi = (sigma2_psi > 0.) ? (sigma2_psi) : (0.);
	sigma2_offset = (sigma2_offset > 0.) ? (sigma2_offset) : (0.);
	range_rot = sigma_cutoff * sqrt(sigma2_rot);  // KThurber
	range_tilt = sigma_cutoff * sqrt(sigma2_tilt);
	range_psi = sigma_cutoff * sqrt(sigma2_psi);
	range2_offset = sigma_cutoff * sigma_cutoff * sigma2_offset;

	// Check helical segments and their polarity
	str_name = list[sid].helical_tube_name;
	subset = list[sid].subset;
	nr_same_polarity = nr_opposite_polarity = 1;  // Laplace smoothing
	unimodal_angular_priors = true;
	for (int id = sid; id <= eid; id++)
	{
		if (list[id].helical_tube_name != str_name)
			REPORT_ERROR("helix.cpp::updatePriorsForOneHelicalTube(): Helical segments do not come from the same tube!");
		if (list[id].subset != subset) // Do I really need this?
			REPORT_ERROR("helix.cpp::updatePriorsForOneHelicalTube(): Helical segments do not come from the same subset!");

		if (list[id].has_wrong_polarity)
		{
			flipPsiTiltForHelicalSegment(list[id].psi_deg, list[id].tilt_deg, list[id].psi_deg, list[id].tilt_deg);
			nr_opposite_polarity++;
		}
		else
		{
			nr_same_polarity++;
		}

		if ((fabs(list[id].psi_flip_ratio - UNIMODAL_PSI_PRIOR_FLIP_RATIO) > 0.01) )
			unimodal_angular_priors = false;
	}

	psi_flip_ratio = ((RFLOAT)(nr_opposite_polarity)) / (((RFLOAT)(nr_opposite_polarity)) + ((RFLOAT)(nr_same_polarity)));
	if ( (unimodal_angular_priors) && (nr_opposite_polarity <= 1) )
		psi_flip_ratio = UNIMODAL_PSI_PRIOR_FLIP_RATIO;
	nr_wrong_polarity = nr_opposite_polarity - 1;

	// Change the polarity of the entire helix if psi_flip_ratio is larger than 0.5
	if (psi_flip_ratio > 0.5)
	{
		for (int id = sid; id <= eid; id++)
		{
			flipPsiTiltForHelicalSegment(list[id].psi_prior_deg, list[id].tilt_prior_deg,
					list[id].psi_prior_deg, list[id].tilt_prior_deg);
			list[id].psi_flip_ratio = (1. - psi_flip_ratio);
		}
	}

	// Calculate new distance-averaged angular priors
	// SHWS 27042020: do two passes: one normal and one with opposite distances and find out which one is the best
	RFLOAT delta_prior_straight = 0., delta_prior_opposite = 0.;
	for (int iflip = 0; iflip < 2; iflip++)
	{

		RFLOAT delta_prior = 0.;
		for (int id = sid; id <= eid; id++)
		{
			RFLOAT this_rot, this_psi, this_tilt, center_pos, this_pos, sum_w, this_w, offset2;
			RFLOAT length_rot_vec, center_x_helix, this_x_helix;  // KThurber
			Matrix1D<RFLOAT> this_ang_vec, sum_ang_vec, this_trans_vec, center_trans_vec, sum_trans_vec;
			Matrix1D<RFLOAT> this_rot_vec, sum_rot_vec;  // KThurber

			// Init
			this_rot = this_psi = this_tilt = center_pos = this_pos = sum_w = this_w = offset2 = 0.;
			this_ang_vec.initZeros(3);
			this_rot_vec.initZeros(2);	// KThurber
			sum_ang_vec.initZeros(3);
			sum_rot_vec.initZeros(2);	// KThurber
			this_trans_vec.initZeros(data_dim);
			center_trans_vec.initZeros(data_dim);
			sum_trans_vec.initZeros(data_dim);

			// Check position
			center_pos = this_pos = list[id].track_pos_A;

			// Calculate weights
			sum_w = this_w = ((do_avg) ? (gaussian1D(this_pos, sigma_segment_dist, center_pos)) : (1.));

			// Analyze orientations
			this_psi = list[id].psi_deg; // REFRESH PSI PRIOR
			this_tilt = list[id].tilt_deg; // REFRESH TILT PRIOR
			Euler_angles2direction(this_psi, this_tilt, this_ang_vec);
			sum_ang_vec = this_ang_vec * this_w;

			// rotation angle all new KThurber
			this_rot = list[id].rot_deg;  // KThurber
			this_rot_vec(0) = cos(DEG2RAD(this_rot));
			this_rot_vec(1) = sin(DEG2RAD(this_rot));
			sum_rot_vec = this_rot_vec * this_w;
			// for adjusting rot angle by shift along helix
			center_x_helix = list[id].dx_A * cos(DEG2RAD(this_psi)) - list[id].dy_A * sin(DEG2RAD(this_psi));
			// end new KThurber

			// Analyze translations
			XX(this_trans_vec) = list[id].dx_prior_A = list[id].dx_A; // REFRESH XOFF PRIOR
			YY(this_trans_vec) = list[id].dy_prior_A = list[id].dy_A; // REFRESH YOFF PRIOR
			if (is_3D_data)
				ZZ(this_trans_vec) = list[id].dz_prior_A = list[id].dz_A; // REFRESH ZOFF PRIOR

			transformCartesianAndHelicalCoords(this_trans_vec, this_trans_vec, (is_3D_data) ? (this_rot) : (0.), (is_3D_data) ? (this_tilt) : (0.), this_psi, CART_TO_HELICAL_COORDS);
			center_trans_vec = this_trans_vec; // Record helical coordinates of the central segment
			if (!is_3D_data)
				XX(this_trans_vec) = 0.; // Do not accumulate translation along helical axis
			else
				ZZ(this_trans_vec) = 0.;
			sum_trans_vec = this_trans_vec * this_w;

			// Local averaging
			if (do_avg)
			{
				for (int idd = sid; idd <= eid; idd++)
				{
					// Find another segment
					if (id == idd)
						continue;

					// Check position
					this_pos = list[idd].track_pos_A;
					if (fabs(this_pos - center_pos) > (sigma_segment_dist * sigma_cutoff))
						continue;

					// Calculate weights
					this_w = gaussian1D(this_pos, sigma_segment_dist, center_pos);
					sum_w += this_w;

					// Analyze orientations
					// KThurber calc rot corrected for length along segment
					// This defines what the sign of pitch should be
					// KThurber unwind rotation angle in order to average
					// note should probably resolve ambiguity of rot=x or x+180 with 2d classes first
					// pitch in Angstroms, because positions are in Angstroms, pitch is 180 degree length in Angstroms
					// for adjusting rot angle by shift along helix
					RFLOAT pitch;
					if (list[idd].classID - 1 >= helical_twist.size()) REPORT_ERROR("ERROR: classID out of range...");
					if (fabs(helical_twist[list[idd].classID - 1]) > 0.)
					{
						RFLOAT pitch = helical_rise[list[idd].classID - 1] * 180. / helical_twist[list[idd].classID - 1];
						this_x_helix = list[idd].dx_A * cos(DEG2RAD(list[idd].psi_deg)) - list[idd].dy_A * sin(DEG2RAD(list[idd].psi_deg));

						// In the second pass, check the direction from large to small distances
						RFLOAT sign = (iflip == 1) ? 1. : -1.;
						this_rot = list[idd].rot_deg + sign*(180./pitch)*(this_pos - center_pos - this_x_helix + center_x_helix);
					}
					else
						this_rot = list[idd].rot_deg;

					this_rot_vec(0) = cos(DEG2RAD(this_rot));
					this_rot_vec(1) = sin(DEG2RAD(this_rot));
					sum_rot_vec += this_rot_vec * this_w;

					this_psi = list[idd].psi_deg;
					this_tilt = list[idd].tilt_deg;
					Euler_angles2direction(this_psi, this_tilt, this_ang_vec);
					sum_ang_vec += this_ang_vec * this_w;

					// Analyze translations
					XX(this_trans_vec) = list[idd].dx_A;
					YY(this_trans_vec) = list[idd].dy_A;
					if (is_3D_data)
						ZZ(this_trans_vec) = list[idd].dz_A;

					transformCartesianAndHelicalCoords(this_trans_vec, this_trans_vec, (is_3D_data) ? (this_rot) : (0.), (is_3D_data) ? (this_tilt) : (0.), this_psi, CART_TO_HELICAL_COORDS);
					if (!is_3D_data)
						XX(this_trans_vec) = 0.; // Do not accumulate translation along helical axis
					else
						ZZ(this_trans_vec) = 0.;
					sum_trans_vec += this_trans_vec * this_w;
				}

				sum_ang_vec /= sum_w;
				Euler_direction2angles(sum_ang_vec, this_psi, this_tilt);

				// KThurber added
				sum_rot_vec /= sum_w;
				length_rot_vec = sqrt(pow(sum_rot_vec(0),2) + pow(sum_rot_vec(1),2));
				if (length_rot_vec!=0)
				{
					sum_rot_vec(0) = sum_rot_vec(0) / length_rot_vec;
					sum_rot_vec(1) = sum_rot_vec(1) / length_rot_vec;
					this_rot = RAD2DEG(acos(sum_rot_vec(0)));
					if (sum_rot_vec(1) < 0.)
						this_rot = -1. * this_rot;	// if sign negative, angle is negative
				}
				else
					this_rot = list[id].rot_deg;  // don't change prior if average fails
				// KThurber end new section

				if (iflip == 0)
				{
					// Distance-averaged priors for original distances in filament
					// cannot store in rot_prior_deg, as will be needed in calculation for iflip==1 pass!
					list[id].rot_prior_deg_ori = this_rot;  // KThurber
					list[id].psi_prior_deg_ori = this_psi; // REFRESH PSI PRIOR
					list[id].tilt_prior_deg_ori = this_tilt; // REFRESH TILT PRIOR
				}
				else
				{
					// Distance-averaged priors for flipped (opposite) filament
					list[id].rot_prior_deg = this_rot;  // KThurber
					list[id].psi_prior_deg = this_psi; // REFRESH PSI PRIOR
					list[id].tilt_prior_deg = this_tilt; // REFRESH TILT PRIOR
				}

				// Keep track how different the priors are from the actual angles to determine straight/opposite
				RFLOAT diff = fabs(list[id].rot_deg - this_rot);
				if (diff > 180.) diff = fabs(diff - 360.);
				// Also count 180 degree errors (basically up-side-down flips) as OK
				// All we're after is the direction of the helix, so upside down errors should be ignored here...
				if (diff > 90.) diff = fabs(diff - 180.);
				delta_prior += diff;
				//std::cout << iflip << " " <<list[id].track_pos_A << " "<< list[id].rot_deg << " " <<this_rot << " "<< diff << " "<<length_rot_vec<< " "<< this_psi << " " <<list[id].psi_prior_deg << " " << list[id].psi_flip_ratio << std::endl;

				// Also do translational priors
				sum_trans_vec /= sum_w;
				offset2 = (YY(sum_trans_vec) - YY(center_trans_vec)) * (YY(sum_trans_vec) - YY(center_trans_vec));
				if (is_3D_data)
					offset2 += (XX(sum_trans_vec) - XX(center_trans_vec)) * (XX(sum_trans_vec) - XX(center_trans_vec));
				if (offset2 > range2_offset) // only now average translations
				{
					if (!is_3D_data)
						XX(sum_trans_vec) = XX(center_trans_vec);
					else
					{
						this_rot = list[id].rot_deg;
						ZZ(sum_trans_vec) = ZZ(center_trans_vec);
					}
					transformCartesianAndHelicalCoords(sum_trans_vec, sum_trans_vec, (is_3D_data) ? (this_rot) : (0.), (is_3D_data) ? (this_tilt) : (0.), this_psi, HELICAL_TO_CART_COORDS); // Averaged translations - use respective averaged tilt and psi
					list[id].dx_prior_A = XX(sum_trans_vec); // REFRESH XOFF PRIOR
					list[id].dy_prior_A = YY(sum_trans_vec); // REFRESH YOFF PRIOR
					if (is_3D_data)
						list[id].dz_prior_A = ZZ(sum_trans_vec); // REFRESH ZOFF PRIOR
				}

			} // end if do_avg
		} // end for id

		if (iflip == 1) delta_prior_opposite = delta_prior / (eid-sid+1);
		else delta_prior_straight = delta_prior / (eid-sid+1);

	} // end for iflip
	//std::cout << " Delta prior straight= " << delta_prior_straight << " opposite= " << delta_prior_opposite << std::endl;

	// Change the direction of the distances in the tube if the total angular difference is smaller for the opposite direction

	if (delta_prior_opposite < delta_prior_straight)
	{
		reverse_direction = true;
	}
	else
	{
		reverse_direction = false;
		for (int id = sid; id <= eid; id++)
		{
			list[id].rot_prior_deg = list[id].rot_prior_deg_ori;
			list[id].psi_prior_deg = list[id].psi_prior_deg_ori;
			list[id].tilt_prior_deg = list[id].tilt_prior_deg_ori;
		}
	}

	// for debugging
	/*
	for (int id = sid; id <= eid; id++)
	{
		std::cout << list[id].track_pos_A << " "<< list[id].rot_deg << " " <<list[id].rot_prior_deg << " "<< list[id].psi_deg << " " <<list[id].psi_prior_deg << " " << list[id].psi_flip_ratio << std::endl;
	}
	*/

}

void updatePriorsForHelicalReconstruction(
		MetaDataTable& MD,
		RFLOAT sigma_segment_dist,
		std::vector<RFLOAT> helical_rise,
		std::vector<RFLOAT> helical_twist,
		int helical_nstart,
		bool is_3D_data,
		bool do_auto_refine,
		RFLOAT sigma2_rot,
		RFLOAT sigma2_tilt,
		RFLOAT sigma2_psi,
		RFLOAT sigma2_offset,
		bool keep_tilt_prior_fixed,
		int verb)
{

	// If we're not averaging angles from neighbouring segments in the helix,
	// then just set the priors to the angles from the previous iteration
	// this effectively means that each segment is completely independent from the rest
	// the priors are only used to center the local angular searches
	if (sigma_segment_dist < 0.)
	{
		updateAngularPriorsForHelicalReconstructionFromLastIter(MD, keep_tilt_prior_fixed);
		return;
	}

	// Check labels
	if (MD.numberOfObjects() < 1)
		REPORT_ERROR("helix.cpp::updatePriorsForHelicalReconstruction: MetaDataTable is empty!");
	if (!MD.containsLabel(EMDL_IMAGE_NAME))
		REPORT_ERROR("helix.cpp::updatePriorsForHelicalReconstruction: rlnImageName is missing!");
	if ( ( (is_3D_data) && (!MD.containsLabel(EMDL_ORIENT_ROT)) )
			|| (!MD.containsLabel(EMDL_ORIENT_TILT))
			|| (!MD.containsLabel(EMDL_ORIENT_PSI))
			|| (!MD.containsLabel(EMDL_ORIENT_ORIGIN_X_ANGSTROM))
			|| (!MD.containsLabel(EMDL_ORIENT_ORIGIN_Y_ANGSTROM))
			|| ( (is_3D_data) && (!MD.containsLabel(EMDL_ORIENT_ORIGIN_Z_ANGSTROM)) )
			|| (!MD.containsLabel(EMDL_ORIENT_TILT_PRIOR))
			|| (!MD.containsLabel(EMDL_ORIENT_PSI_PRIOR))
			|| (!MD.containsLabel(EMDL_PARTICLE_HELICAL_TUBE_ID))
			|| (!MD.containsLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM))
			|| (!MD.containsLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO))
			|| ( (do_auto_refine) && (!MD.containsLabel(EMDL_PARTICLE_RANDOM_SUBSET)) ) )
		REPORT_ERROR("helix.cpp::updatePriorsForHelicalReconstruction: Labels of helical prior information are missing!");

	std::vector<HelicalSegmentPriorInfoEntry> list;
	long int MDobjectID;


	// For N-start helices, revert back to the N-start twist and rise (not the 1-start ones)
	// This is in order to reduce amplification of small deviations in twist and rise
	if (helical_nstart > 1)
	{
		// Assume same N-start for all classes
		// Shaoda's formula (which need to be inverted, as we want original N-start rise and twist back)
		// rise_1-start = rise / N
		// twist_1-start = (twist+360)/N if twist>0
		// twist_1-start = (twist-360)/N if twist<0

		for (int iclass=0; iclass < helical_rise.size(); iclass++)
		{
			helical_rise[iclass] *= helical_nstart;
			RFLOAT aux = helical_twist[iclass] * helical_nstart;
			helical_twist[iclass] = (aux > 360.) ? aux - 360. : aux + 360.;
			if (verb > 0) std::cout << " + for rotational priors go back to " << helical_nstart
					<< "-start helical twist= " << helical_twist[iclass] << " and rise= " << helical_rise[iclass] << std::endl;
		}
	}

	// Read _data.star file
	list.clear();
	MDobjectID = -1;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		HelicalSegmentPriorInfoEntry segment;
		std::string str_mic;
		int tube_id;

		segment.clear();

		MD.getValue(EMDL_MICROGRAPH_NAME, str_mic);
		MD.getValue(EMDL_PARTICLE_HELICAL_TUBE_ID, tube_id);
		segment.helical_tube_name = str_mic + integerToString(tube_id);
		MD.getValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, segment.track_pos_A);
		if (MD.containsLabel(EMDL_ORIENT_ROT)) MD.getValue(EMDL_ORIENT_ROT, segment.rot_deg);  		// KThurber
		else segment.rot_deg = 0.;
		if (MD.containsLabel(EMDL_ORIENT_ROT_PRIOR)) MD.getValue(EMDL_ORIENT_ROT_PRIOR, segment.rot_prior_deg);  	// KThurber
		//else segment.rot_prior_deg = 0.;
		else segment.rot_prior_deg = segment.rot_deg;  // SHWS, modified from KThurber!
		MD.getValue(EMDL_ORIENT_TILT, segment.tilt_deg);
		MD.getValue(EMDL_ORIENT_TILT_PRIOR, segment.tilt_prior_deg);
		MD.getValue(EMDL_ORIENT_PSI, segment.psi_deg);
		MD.getValue(EMDL_ORIENT_PSI_PRIOR, segment.psi_prior_deg);
		MD.getValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, segment.psi_flip_ratio);
		if (MD.containsLabel(EMDL_ORIENT_PSI_PRIOR_FLIP))			// KThurber2
			MD.getValue(EMDL_ORIENT_PSI_PRIOR_FLIP, segment.psi_prior_flip);
		else segment.psi_prior_flip = false;
		if (MD.containsLabel(EMDL_PARTICLE_CLASS))
			MD.getValue(EMDL_PARTICLE_CLASS, segment.classID);
		else
			segment.classID = 1;
		if (do_auto_refine)
			MD.getValue(EMDL_PARTICLE_RANDOM_SUBSET, segment.subset); // Do I really need this?

		MD.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, segment.dx_A);
		MD.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, segment.dy_A);
		if (is_3D_data)
			MD.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, segment.dz_A);

		segment.checkPsiPolarity();

		MDobjectID++;
		segment.MDobjectID = MDobjectID;
		list.push_back(segment);
	}

	// Sort the list so that segments from the same helical tube come together
	std::stable_sort(list.begin(), list.end());

	// Loop over every helical tube
	long total_opposite_polarity = 0;
	long total_opposite_rot = 0;		// KThurber
	long total_same_rot = 0;		// KThurber
	for (int sid = 0; sid < list.size(); )
	{
		// A helical tube [id_s, id_e]
		int nr_opposite_polarity = -1;
		int eid = sid; // start id (sid) and end id (eid)
		while (1)
		{
			eid++;
			if (eid >= list.size())
				break;
			if (list[eid].helical_tube_name != list[sid].helical_tube_name)
				break;
		}
		eid--;

		// Real work...
		bool reverse_direction;
		updatePriorsForOneHelicalTube(list, sid, eid, nr_opposite_polarity, reverse_direction, sigma_segment_dist, helical_rise, helical_twist,
				is_3D_data, do_auto_refine, sigma2_rot, sigma2_tilt, sigma2_psi, sigma2_offset);
		total_opposite_polarity += nr_opposite_polarity;
		if (reverse_direction) total_opposite_rot += 1;
		else total_same_rot += 1;

		// Write to _data.star file
		for (int id = sid; id <= eid; id++)
		{
			if (reverse_direction)
				MD.setValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, -1. * list[id].track_pos_A, list[id].MDobjectID);
			if (!keep_tilt_prior_fixed)
				MD.setValue(EMDL_ORIENT_TILT_PRIOR, list[id].tilt_prior_deg, list[id].MDobjectID);
			MD.setValue(EMDL_ORIENT_PSI_PRIOR, list[id].psi_prior_deg, list[id].MDobjectID);
			MD.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, list[id].psi_flip_ratio, list[id].MDobjectID);
			MD.setValue(EMDL_ORIENT_ROT_PRIOR, list[id].rot_prior_deg, list[id].MDobjectID); // KThurber
			MD.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, list[id].dx_prior_A, list[id].MDobjectID);
			MD.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, list[id].dy_prior_A, list[id].MDobjectID);
			if (is_3D_data)
				MD.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, list[id].dz_prior_A, list[id].MDobjectID);
		}

		// Next helical tube
		sid = eid + 1;
	}

	list.clear();

	if ( (verb > 0) )
	{

		long total_same_polarity = MD.numberOfObjects() - total_opposite_polarity;
		RFLOAT opposite_percentage = (100.) * ((RFLOAT)(total_opposite_polarity)) / ((RFLOAT)(MD.numberOfObjects()));
		RFLOAT opposite_percentage_rot = (100.) * ((RFLOAT)(total_opposite_rot)) / ((RFLOAT)(total_same_rot + total_opposite_rot));

		std::cout << " Number of helical segments with same / opposite polarity to their psi priors: " << total_same_polarity << " / " << total_opposite_polarity << " (" << opposite_percentage << "%)" << std::endl;
		std::cout << " Number of helices with same / reverse direction for rot priors: " << total_same_rot << " / " << total_opposite_rot << " (" << opposite_percentage_rot << "%)" << std::endl;
	}

}

void updateAngularPriorsForHelicalReconstructionFromLastIter(
		MetaDataTable& MD,
		bool keep_tilt_prior_fixed)
{
	if (MD.numberOfObjects() < 1)
		REPORT_ERROR("helix.cpp::updateAngularPriorsForHelicalReconstruction: MetaDataTable is empty!");

	bool have_tilt = MD.containsLabel(EMDL_ORIENT_TILT);
	bool have_psi = MD.containsLabel(EMDL_ORIENT_PSI);
	bool have_tilt_prior = MD.containsLabel(EMDL_ORIENT_TILT_PRIOR);
	bool have_psi_prior = MD.containsLabel(EMDL_ORIENT_PSI_PRIOR);
	bool have_rot = MD.containsLabel(EMDL_ORIENT_ROT);		// KThurber
	bool have_rot_prior = MD.containsLabel(EMDL_ORIENT_ROT_PRIOR);	// KThurber

	if ( (!have_tilt_prior) && (!have_psi_prior) && (!have_rot_prior))	// KThurber
		return;

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		RFLOAT val;
		if (have_tilt && have_tilt_prior && (!keep_tilt_prior_fixed) )
		{
			MD.getValue(EMDL_ORIENT_TILT, val);
			MD.setValue(EMDL_ORIENT_TILT_PRIOR, val);
		}
		if (have_psi && have_psi_prior)
		{
			MD.getValue(EMDL_ORIENT_PSI, val);
			MD.setValue(EMDL_ORIENT_PSI_PRIOR, val);
		}
		// KThurber add rot section
		if (have_rot && have_rot_prior)
		{
			MD.getValue(EMDL_ORIENT_ROT, val);
			MD.setValue(EMDL_ORIENT_ROT_PRIOR, val);
		}
	}
}

void setPsiFlipRatioInStarFile(MetaDataTable& MD, RFLOAT ratio)
{
	if (!MD.containsLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO))
		REPORT_ERROR("helix.cpp::setPsiFlipRatioInStarFile: Psi flip ratio is not found in this STAR file!");
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		MD.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, ratio);
	}
}

void plotLatticePoints(MetaDataTable& MD,
		int x1, int y1, int x2, int y2)
{
	MD.clear();
	MD.addLabel(EMDL_IMAGE_COORD_X);
	MD.addLabel(EMDL_IMAGE_COORD_Y);
	for (int i = -10; i <= 10; i++)
	{
		for (int j = -10; j <= 10; j++)
		{
			MD.addObject();
			MD.setValue(EMDL_IMAGE_COORD_X, RFLOAT(i * x1 + j * x2));
			MD.setValue(EMDL_IMAGE_COORD_Y, RFLOAT(i * y1 + j * y2));
		}
	}
}

void grabParticleCoordinates(
		FileName& fn_in,
		FileName& fn_out)
{
	MetaDataTable MD_in, MD_out;
	RFLOAT x, y, z;
	bool contain_z_coord = false;

	if (fn_in.getExtension() != "star")
		REPORT_ERROR("helix.cpp::grabParticleCoordinates: Input file must have STAR extension!");

	MD_in.clear();
	MD_in.read(fn_in);
	if ( (!MD_in.containsLabel(EMDL_IMAGE_COORD_X)) || (!MD_in.containsLabel(EMDL_IMAGE_COORD_Y)) )
		REPORT_ERROR("helix.cpp::grabParticleCoordinates: Input file must have X and Y coordinates!");
	contain_z_coord = MD_in.containsLabel(EMDL_IMAGE_COORD_Z);

	MD_out.clear();
	x = y = z = 0.;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_in)
	{
		MD_in.getValue(EMDL_IMAGE_COORD_X, x);
		MD_in.getValue(EMDL_IMAGE_COORD_Y, y);
		if (contain_z_coord)
			MD_in.getValue(EMDL_IMAGE_COORD_Z, z);

		MD_out.addObject();
		MD_out.setValue(EMDL_IMAGE_COORD_X, x);
		MD_out.setValue(EMDL_IMAGE_COORD_Y, y);
		if (contain_z_coord)
			MD_out.setValue(EMDL_IMAGE_COORD_Z, z);
	}
	MD_out.write(fn_out);

	return;
}

void grabParticleCoordinates_Multiple(
		std::string& suffix_fin,
		std::string& suffix_fout)
{
	FileName fns_in;
	std::vector<FileName> fn_in_list;

	if (suffix_fin == suffix_fout)
		REPORT_ERROR("helix.cpp::grabParticleCoordinates_Multiple(): File names error!");

	fns_in = "*" + suffix_fin;
	fns_in.globFiles(fn_in_list);
	std::cout << "Number of input files = " << fn_in_list.size() << std::endl;
	if (fn_in_list.size() < 1)
		REPORT_ERROR("helix.cpp::grabParticleCoordinates_Multiple(): No input files are found!");

	for (int ii = 0; ii < fn_in_list.size(); ii++)
	{
		FileName fn_out = fn_in_list[ii].beforeFirstOf(suffix_fin) + suffix_fout;
		grabParticleCoordinates(fn_in_list[ii], fn_out);
	}
	return;
}

void calculateRadialAvg(MultidimArray<RFLOAT> &v, RFLOAT angpix)
{
	std::vector<RFLOAT> rval, rcount;
	long int size, rint;

	if ( (XSIZE(v) < 5) || (YSIZE(v) < 5) || (ZSIZE(v) < 5) || (NSIZE(v) != 1) )
		REPORT_ERROR("helix.cpp::calculateRadialAvg(): Input image should be a 3D box larger than 5*5*5 !");
	if (!(angpix > 0.))
		REPORT_ERROR("helix.cpp::calculateRadialAvg(): Pixel size should be larger than 0 !");

	v.setXmippOrigin();
	size = (XSIZE(v) < YSIZE(v)) ? XSIZE(v) : YSIZE(v);
	size = size / 2 + 2;

	rval.resize(size);
	rcount.resize(size);
	for (int ii = 0; ii < rval.size(); ii++)
		rval[ii] = rcount[ii] = 0.;

	FOR_ALL_ELEMENTS_IN_ARRAY3D(v)
	{
		rint = ROUND(sqrt((RFLOAT)(i * i + j * j)));
		if (rint >= size)
			continue;

		rval[rint] += A3D_ELEM(v, k, i, j);
		rcount[rint] += 1.;
	}

	for (int ii = 0; ii < rval.size(); ii++)
	{
		if (rcount[ii] < 0.5)
			rval[ii] = 0.;
		else
			rval[ii] /= rcount[ii];
		std::cout << ii * angpix << "      " << rval[ii] << std::endl;
	}
}

void transformCartesianToHelicalCoordsForStarFiles(
		MetaDataTable& MD_in,
		MetaDataTable& MD_out)
{
	RFLOAT rot, tilt, psi, xoff, yoff, zoff;
	bool is_3d_trans = MD_in.containsLabel(EMDL_ORIENT_ORIGIN_Z_ANGSTROM);

	MD_out.clear();
	MD_out = MD_in;

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_in)
	{
		MD_in.getValue(EMDL_ORIENT_ROT, rot);
		MD_in.getValue(EMDL_ORIENT_TILT, tilt);
		MD_in.getValue(EMDL_ORIENT_PSI, psi);
		MD_in.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff);
		MD_in.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff);
		if (is_3d_trans)
			MD_in.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, zoff);

		transformCartesianAndHelicalCoords(
				xoff, yoff, zoff,
				xoff, yoff, zoff,
				rot, tilt, psi,
				((is_3d_trans) ? (3) : (2)),
				CART_TO_HELICAL_COORDS);

		MD_out.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff);
		MD_out.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff);
		if (is_3d_trans)
			MD_out.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, zoff);
		MD_out.nextObject();
	}
}

void normaliseHelicalSegments(
		FileName& fn_in,
		FileName& fn_out_root,
		RFLOAT helical_outer_diameter_A,
		RFLOAT pixel_size_A)
{
	bool is_3D_data = false, is_mrcs = false, have_tilt_prior = false, have_psi_prior = false, read_angpix_from_star = false;
	RFLOAT rot_deg = 0., tilt_deg = 0., psi_deg = 0., avg = 0., stddev = 0., val = 0., det_pixel_size = 0., mag = 0.;
	Image<RFLOAT> img0;
	MetaDataTable MD;
	FileName img_name, file_ext;

	if (fn_in.getExtension() != "star")
		REPORT_ERROR("helix.cpp::normaliseHelicalSegments(): Please provide a STAR file as input!");

	// Read STAR file
	MD.clear();
	MD.read(fn_in);
	have_tilt_prior = MD.containsLabel(EMDL_ORIENT_TILT_PRIOR);
	have_psi_prior = MD.containsLabel(EMDL_ORIENT_PSI_PRIOR);
	read_angpix_from_star = (MD.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE)) && (MD.containsLabel(EMDL_CTF_MAGNIFICATION));
	if ( (!MD.containsLabel(EMDL_IMAGE_NAME)) )
		REPORT_ERROR("helix.cpp::normaliseHelicalSegments(): MetaDataLabel _rlnImageName is missing!");
	if ( (!have_tilt_prior) && (!MD.containsLabel(EMDL_ORIENT_TILT)) )
		REPORT_ERROR("helix.cpp::normaliseHelicalSegments(): MetaDataLabel _rlnAngleTilt or _rlnAngleTiltPrior is missing!");
	if ( (!have_psi_prior) && (!MD.containsLabel(EMDL_ORIENT_PSI)) )
		REPORT_ERROR("helix.cpp::normaliseHelicalSegments(): MetaDataLabel _rlnAnglePsi or _rlnAnglePsiPrior is missing!");

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		// Read image name and angular priors
		MD.getValue(EMDL_IMAGE_NAME, img_name);
		file_ext = img_name.getExtension();
		is_mrcs = (file_ext == "mrcs");
		rot_deg = tilt_deg = psi_deg = 0.;
		if (have_tilt_prior)
			MD.getValue(EMDL_ORIENT_TILT_PRIOR, tilt_deg);
		else
			MD.getValue(EMDL_ORIENT_TILT, tilt_deg);
		if (have_psi_prior)
			MD.getValue(EMDL_ORIENT_PSI_PRIOR, psi_deg);
		else
			MD.getValue(EMDL_ORIENT_PSI, psi_deg);
		if (read_angpix_from_star)
		{
	    	MD.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, det_pixel_size);
	    	MD.getValue(EMDL_CTF_MAGNIFICATION, mag);
	    	pixel_size_A = det_pixel_size * 10000. / mag;
		}
		// DEBUG
		//std::cout << " pixel_size_A = " << pixel_size_A << std::endl;
		if (pixel_size_A < (1e-4))
			REPORT_ERROR("helix.cpp::normaliseHelicalSegments(): Invalid pixel size for image " + ((std::string)(img_name)));
		if ((helical_outer_diameter_A / pixel_size_A) < 10.)
			REPORT_ERROR("helix.cpp::normaliseHelicalSegments(): Diameter of the tubular mask should be larger than 10 pixels!");

		// Read image
		img0.clear();
		img0.read(img_name);
		is_3D_data = (ZSIZE(img0()) > 1) || (NSIZE(img0()) > 1);
		if ( (XSIZE(img0()) < (helical_outer_diameter_A / pixel_size_A)) || (YSIZE(img0()) < (helical_outer_diameter_A / pixel_size_A)) )
			REPORT_ERROR("helix.cpp::normaliseHelicalSegments(): Diameter of the tubular mask is larger than the box XY dimensions!");
		if (!is_3D_data)
			rot_deg = tilt_deg = 0.;

		// Calculate avg and stddev
		calculateBackgroundAvgStddev(
				img0,
				avg,
				stddev,
				0,
				true,
				helical_outer_diameter_A * 0.5 / pixel_size_A,
				tilt_deg,
				psi_deg);
		if (stddev < 0.0001)
		{
			std::cout << " !!! WARNING: " << img_name << "  has bg_avg = " << avg << " and bg_stddev = " << stddev << " . bg_stddev is set to 0.0001. The image cannot be properly normalised!" << std::endl;
			stddev = 0.0001;
		}
		else
			std::cout << "  Normalising " << img_name << " with bg_avg = " << avg << " and bg_stddev = " << stddev << " . " << std::endl;

		// Normalise
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img0())
		{
			val = DIRECT_MULTIDIM_ELEM(img0(),n);
			DIRECT_MULTIDIM_ELEM(img0(), n) = (val - avg) / stddev;
		}

		// Rename
		img_name = img_name.withoutExtension() + fn_out_root + "." + file_ext;

		if (is_3D_data)
		{
			img0.setSamplingRateInHeader(pixel_size_A, pixel_size_A, pixel_size_A);
			img0.setStatisticsInHeader();
		}

		// Write
		if (is_3D_data)
			img0.write(img_name);
		else
			img0.write(img_name, -1, true, WRITE_APPEND);
		img0.clear();
	}

	if (!is_3D_data)
	{
		// Read the header of .mrcs stack
		img_name = img_name.substr(img_name.find("@") + 1);
		// Set the pixel size in the file header
		img0.read(img_name);
		img0.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, pixel_size_A);
		img0.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, pixel_size_A);
		img0.write(img_name);
	}
}

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
		RFLOAT tension,
		RFLOAT bias)
{
	RFLOAT m0 = 0., m1 = 0., mu2 = 0., mu3 = 0., a0 = 0., a1 = 0., a2 = 0., a3 = 0.;

	mu2 = mu * mu;
	mu3 = mu2 * mu;
	m0  = (y1 - y0) * (1. + bias) * (1. - tension) / 2.;
	m0 += (y2 - y1) * (1. - bias) * (1. - tension) / 2.;
	m1  = (y2 - y1) * (1. + bias) * (1. - tension) / 2.;
	m1 += (y3 - y2) * (1. - bias) * (1. - tension) / 2.;
	a0  = 2. * mu3 - 3. * mu2 + 1.;
	a1  = mu3 - 2. * mu2 + mu;
	a2  = mu3 - mu2;
	a3  = (-2.) * mu3 + 3. * mu2;

	return (a0 * y1 + a1 * m0 + a2 * m1 + a3 * y2);
}

void HermiteInterpolateOne3DHelicalFilament(
		MetaDataTable& MD_in,
		MetaDataTable& MD_out,
		int& total_segments,
		int nr_asu,
		RFLOAT rise_A,
		RFLOAT pixel_size_A,
		RFLOAT box_size_pix,
		int helical_tube_id,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT Zdim,
		bool bimodal_angular_priors)
{
	RFLOAT x0, x1, x2, x3, xa, xb, y0, y1, y2, y3, ya, yb, z0, z1, z2, z3, za, zb, mu1, mu2;
	RFLOAT step_pix, chord_pix, accu_len_pix, present_len_pix, len_pix, psi_prior_flip_ratio, tilt_deg, psi_deg;
    RFLOAT half_box_size_pix = box_size_pix / 2.;
	int nr_partitions, nr_segments;
	std::vector<RFLOAT> xlist, ylist, zlist;
	Matrix1D<RFLOAT> dr;

	// DEBUG: Do not exclude particles on the edges
	// Xdim = Ydim = Zdim = 999999.;

	if (MD_in.numberOfObjects() <= 1)
		REPORT_ERROR("helix.cpp::HermiteInterpolateOne3DHelicalFilament(): MetaDataTable should have at least two points for interpolation!");
	if ( (box_size_pix < 2) || (Xdim < box_size_pix) || (Ydim < box_size_pix) || (Zdim < box_size_pix) )
		REPORT_ERROR("helix.cpp::HermiteInterpolateOne3DHelicalFilament(): Wrong dimensions or box size!");
	if (pixel_size_A < 0.001)
		REPORT_ERROR("helix.cpp::HermiteInterpolateOne3DHelicalFilament(): Invalid pixel size!");
	RFLOAT interbox_pix = ((RFLOAT)(nr_asu)) * rise_A / pixel_size_A;
	if ( (nr_asu < 1) || (rise_A < 0.001) || (interbox_pix < 0.58) )
		REPORT_ERROR("helix.cpp::HermiteInterpolateOne3DHelicalFilament(): Invalid helical rise or number of asymmetrical units!");
	if ( (!MD_in.containsLabel(EMDL_IMAGE_COORD_X))
			|| (!MD_in.containsLabel(EMDL_IMAGE_COORD_Y))
			|| (!MD_in.containsLabel(EMDL_IMAGE_COORD_Z)) )
		REPORT_ERROR("helix.cpp::HermiteInterpolateOne3DHelicalFilament(): MetaDataTable should have _rlnOriginX, _rlnOriginY and _rlnOriginZ!");

	// Header of output file
	MD_out.clear();
	MD_out.addLabel(EMDL_IMAGE_COORD_X);
	MD_out.addLabel(EMDL_IMAGE_COORD_Y);
	MD_out.addLabel(EMDL_IMAGE_COORD_Z);
	MD_out.addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);
	MD_out.addLabel(EMDL_ORIENT_TILT_PRIOR);
	MD_out.addLabel(EMDL_ORIENT_PSI_PRIOR);
    MD_out.addLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM);
    MD_out.addLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO);

	//half_box_size_pix = box_size_pix / 2.;
    psi_prior_flip_ratio = UNIMODAL_PSI_PRIOR_FLIP_RATIO;
    if (bimodal_angular_priors)
    {
     	psi_prior_flip_ratio = BIMODAL_PSI_PRIOR_FLIP_RATIO;
     }

    // Load all manually picked coordinates
    xlist.clear(); ylist.clear(); zlist.clear();
    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_in)
    {
    	MD_in.getValue(EMDL_IMAGE_COORD_X, x0);
    	MD_in.getValue(EMDL_IMAGE_COORD_Y, y0);
    	MD_in.getValue(EMDL_IMAGE_COORD_Z, z0);
    	xlist.push_back(x0);
    	ylist.push_back(y0);
    	zlist.push_back(z0);
    }

    // Interpolate
    accu_len_pix = 0.;
    present_len_pix = -1.;
    nr_segments = 0;
    dr.initZeros(3);
    for (int id = 0; id < (xlist.size() - 1); id++)
    {
        // Step size for interpolation is smaller than 1% of the inter-box distance
        // sqrt(0.57735 * 0.57735 * 0.57735 * 3) = 1.0, step size is larger than 1 pixel
    	// TODO: 1% ? Too expensive computationally? Try 10% ?
        step_pix = (interbox_pix < 57.735) ? (0.57735) : (interbox_pix / 100.); // 1%
        //step_pix = (interbox_pix < 5.7735) ? (0.57735) : (interbox_pix / 10.); // 10%

    	// Collect points 0, 1, 2, 3 for interpolations
        x0 = x1 = x2 = x3 = y0 = y1 = y2 = y3 = z0 = z1 = z2 = z3 = 0.;
    	// Point 0
    	if (id == 0)
    	{
    		x0 = 2. * xlist[id] - xlist[id + 1];
    		y0 = 2. * ylist[id] - ylist[id + 1];
    		z0 = 2. * zlist[id] - zlist[id + 1];
    	}
    	else
    	{
    		x0 = xlist[id - 1]; y0 = ylist[id - 1]; z0 = zlist[id - 1];
    	}
    	// Point 1 and Point 2
    	x1 = xlist[id]; y1 = ylist[id]; z1 = zlist[id];
    	x2 = xlist[id + 1]; y2 = ylist[id + 1]; z2 = zlist[id + 1];
    	// Point 3
    	if (id == (xlist.size() - 2))
    	{
    		x3 = 2. * xlist[id + 1] - xlist[id];
    		y3 = 2. * ylist[id + 1] - ylist[id];
    		z3 = 2. * zlist[id + 1] - zlist[id];
    	}
    	else
    	{
    		x3 = xlist[id + 2]; y3 = ylist[id + 2]; z3 = zlist[id + 2];
    	}

    	// Chord distance between point 1 and 2
    	// TODO: what will happen if the chord length is smaller than step_pix?
    	chord_pix = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
    	nr_partitions = int(CEIL(chord_pix / step_pix));
    	nr_partitions = (nr_partitions <= 0) ? (1) : (nr_partitions);

    	// Partitioning
    	for (int ip = 0; ip < nr_partitions; ip++)
    	{
    		xa = ya = za = xb = yb = zb = mu1 = mu2 = len_pix = 0.;
    		mu1 = RFLOAT( (RFLOAT(ip)) / (RFLOAT(nr_partitions)) );
    		mu2 = RFLOAT( (RFLOAT(ip) + 1.) / (RFLOAT(nr_partitions)) );
    		xa = HermiteInterpolate1D(x0, x1, x2, x3, mu1);
    		ya = HermiteInterpolate1D(y0, y1, y2, y3, mu1);
    		za = HermiteInterpolate1D(z0, z1, z2, z3, mu1);
    		xb = HermiteInterpolate1D(x0, x1, x2, x3, mu2);
    		yb = HermiteInterpolate1D(y0, y1, y2, y3, mu2);
    		zb = HermiteInterpolate1D(z0, z1, z2, z3, mu2);
    		len_pix = sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb));
    		present_len_pix += len_pix;
    		accu_len_pix += len_pix;

    		// Output one segment (xb, yb, zb)
    		if (present_len_pix > 0.)
    		{
    			present_len_pix -= interbox_pix;
#ifdef EXCLUDE_SEGMENTS_ON_THE_EDGES
    			// Check whether this segment lies on the edges of the 3D tomogram
    			if ( (xb < half_box_size_pix) || (xb > (Xdim - half_box_size_pix))
    					|| (yb < half_box_size_pix) || (yb > (Ydim - half_box_size_pix))
    					|| (zb < half_box_size_pix) || (zb > (Zdim - half_box_size_pix)) )
    			{
    				std::cout << std::resetiosflags(std::ios::fixed);
    				std::cout << " WARNING: Particle at (" << xb << ", " << yb << ", " << zb << ") is ignored because it is too close to the edge. " << std::flush;
    				std::cout << " Box_size_pix = " << box_size_pix << ", Dimensions = " << Xdim << " * " << Ydim << " * " << Zdim << " ." << std::flush;
    				std::cout << " Please choose a smaller box size OR reconstruct the 3D tomogram with a larger number of Z slices!" << std::endl;
    				continue;
    			}
#endif

    			// Add this segment to the list
    			nr_segments++;
    			MD_out.addObject();
    			MD_out.setValue(EMDL_IMAGE_COORD_X, xb);
    	    	MD_out.setValue(EMDL_IMAGE_COORD_Y, yb);
    	    	MD_out.setValue(EMDL_IMAGE_COORD_Z, zb);
    	    	MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, helical_tube_id);

        		XX(dr) = xb - xa; YY(dr) = yb - ya; ZZ(dr) = zb - za;
        		estimateTiltPsiPriors(dr, tilt_deg, psi_deg);
        		if (fabs(tilt_deg) < 0.001)
        			tilt_deg = 0.;
        		if (fabs(psi_deg) < 0.001)
        			psi_deg = 0.;
    	    	MD_out.setValue(EMDL_ORIENT_TILT_PRIOR, tilt_deg);
    	    	MD_out.setValue(EMDL_ORIENT_PSI_PRIOR, psi_deg);

    	        MD_out.setValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, accu_len_pix * pixel_size_A);
    	        MD_out.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, psi_prior_flip_ratio);
    	     }
    	}
    }
    total_segments = nr_segments;
}

void Interpolate3DCurves(
		FileName& fn_in_root,
		FileName& fn_out_root,
		int nr_asu,
		RFLOAT rise_A,
		RFLOAT pixel_size_A,
		RFLOAT box_size_pix,
		int binning_factor,
		bool bimodal_angular_priors)
{
	Image<RFLOAT> img;
	std::vector<FileName> fn_in_list;
	FileName fn_tomo, fn_in_glob, fn_in, fn_out;
	std::ifstream fin;
	std::string line;
	std::vector<std::string> words;
	std::vector<RFLOAT> xlist, ylist, zlist;
	MetaDataTable MD_in, MD_out, MD_all;
	RFLOAT x0, y0, z0, val;
	int nr_points = 0, total_segments = 0, nr_segments = 0;
	int xdim = 0, ydim = 0, zdim = 0, xdim_img = 0, ydim_img = 0, zdim_img = 0;
	long int ndim = 0, ndim_img = 0;
	char buf_word[4], buf_qword[16], tmp_char;
	bool contain_3d_points = false, flip_YZ = false;

	// General parameter checks
	if (binning_factor < 1)
		REPORT_ERROR("helix.cpp::Interpolate3DCurves(): Binning factor should be larger than 1!");
	fn_tomo = fn_in_root + ".mrc";
	if (!exists(fn_tomo))
		REPORT_ERROR("helix.cpp::Interpolate3DCurves(): Reconstructed 3D tomogram " + fn_tomo + " is not found!");
	// Read the header of 3D reconstructed tomogram
	img.clear();
	img.read(fn_tomo, false);
	img.getDimensions(xdim_img, ydim_img, zdim_img, ndim_img);
	if (zdim_img <= 1)
		REPORT_ERROR("helix.cpp::Interpolate3DCurves(): Dimension Z of reconstructed 3D tomogram " + fn_tomo + " is 1!");
	if (ndim_img != 1)
		REPORT_ERROR("helix.cpp::Interpolate3DCurves(): Dimension N of reconstructed 3D tomogram " + fn_tomo + " is not 1!");

	// Glob all files
	if (fn_in_root.length() <= 1)
		REPORT_ERROR("helix.cpp::Interpolate3DCurves(): Input rootname is an empty string!");
	//std::cout << " fn_in_root = " << fn_in_root << std::endl; // DEBUG
	//fn_in_glob = "*" + fn_in_root + "*"; // OLD
	fn_in_glob = fn_in_root + "*"; // NEW
	fn_in_glob.globFiles(fn_in_list);
	//std::cout << " fn_in_glob = " << fn_in_glob << std::endl; // DEBUG
	//std::cout << " nr_giles_globbed = " << fn_in_list.size() << std::endl;

	if (fn_in_list.size() < 1)
		REPORT_ERROR("helix.cpp::Interpolate3DCurves(): No input files are found!");

	// Check input filenames
	std::cout << " #############################################################" << std::endl;
	std::cout << " Coordinate files (.mod, .star or .txt, .coords text files with XYZ coordinates) to be processed: " << std::endl;
	for (int fid = 0; fid < fn_in_list.size(); fid++)
	{
		fn_in = fn_in_list[fid];
		std::string fn_ext = fn_in.getExtension();
		// Method 1
		//if ( (fn_ext == "") || (fn_ext == "log") || (fn_ext == "mrc") || (fn_ext == "mrcs") || (fn_ext == "ctf")
		//		|| (fn_ext == "st") || (fn_ext == "order") || (fn_ext == "tlt") || (fn_ext == "trial"))
		//{
		//	fn_in_list.erase(fn_in_list.begin() + fid);
		//	fid--;
		//	continue;
		//}
		// Method 2
		if ( (fn_ext != "mod") && (fn_ext != "star") && (fn_ext != "coords") && (fn_ext != "txt") )
		{
			fn_in_list.erase(fn_in_list.begin() + fid);
			fid--;
			continue;
		}
		if ( (fn_in.contains(fn_in_root)) && ( (fn_in.afterFirstOf(fn_in_root)).contains(fn_in_root) ) )
			REPORT_ERROR("helix.cpp::Interpolate3DCurves(): Every input filename should contain one and only one input rootname! Invalid filename: " + (std::string)(fn_in));
		std::cout << "  " << fn_in << std::endl;
	}
	std::cout << " " << fn_in_list.size() << " files (filaments) found." << std::endl;
	std::cout << " Please check whether all coordinate files are included and remove files which do not contain manually picked coordinates!" << std::endl;
	std::cout << " #############################################################" << std::endl;
	std::cout << " Are all coordinate files included? Do all of the files shown above contain manually picked XYZ coordinates? (y/n): " << std::flush;
	line.clear();
	std::cin >> line;
	if ( (line[0] != 'y') && (line[0] != 'Y') )
	{
		std::cout << std::endl << " No! Exit now..." << std::endl;
		return;
	}

	// Check endianness - for MOD files only
	int num = 1;
	bool is_little_endian = false;
	if(*(char *)&num == 1)
		is_little_endian = true;
	//std::cout << is_little_endian << std::endl;

	// Real work begins...
	MD_all.clear();
	total_segments = nr_segments = 0;
	for (int fid = 0; fid < fn_in_list.size(); fid++)
	{
		contain_3d_points = false;
		flip_YZ = false;
		xlist.clear(); ylist.clear(); zlist.clear();

		// Open an input file
		fn_in = fn_in_list[fid];
		std::cout << " ### Input filename = " << fn_in << std::endl;

		// MOD file format definition: http://bio3d.colorado.edu/imod/doc/binspec.html
		if (fn_in.getExtension() == "mod")
		{
			fin.open(fn_in.c_str(), std::ios_base::in|std::ios_base::binary);
			if (fin.fail())
				REPORT_ERROR("helix.cpp::Interpolate3DCurves(): Cannot open input file: " + (std::string)(fn_in));

			// Scheme 1 - does not work for MOD files with more than one 'OBJT's (objects)
			/*
			fin.read(reinterpret_cast<char*>(buf_qword), sizeof(buf_qword)); // Read the first line
			for (int id = 0; id < 14; id++) // Read model data structure (232 bytes)
				fin.read(reinterpret_cast<char*>(buf_qword), sizeof(buf_qword));
			for (int id = 0; id < 11; id++) // Read object data structure (first 160 out of 176 bytes)
				fin.read(reinterpret_cast<char*>(buf_qword), sizeof(buf_qword));
			// 26 lines in total omitted
			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word));
			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word));
			*/
			// Scheme 2
			fin.read(reinterpret_cast<char*>(buf_qword), sizeof(buf_qword)); // Read the first 16 bytes
			fin.read(reinterpret_cast<char*>(buf_qword), sizeof(buf_qword)); // Read the second 16 bytes
			if ( (buf_qword[0] != 'M') || (buf_qword[1] != 'o') || (buf_qword[2] != 'd') || (buf_qword[3] != 'e') || (buf_qword[4] != 'l') )
				REPORT_ERROR("helix.cpp::Interpolate3DCurves(): IMOD file header does not contain 'Model' tag!");
			for (int id = 0; id < 6; id++)
				fin.read(reinterpret_cast<char*>(buf_qword), sizeof(buf_qword)); // Read the next 96 bytes
			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word));
			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word)); // Read the next 8 bytes
			// Name of model ends
			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word)); // Xdim
			if (is_little_endian)
			{
				SWAP(buf_word[0], buf_word[3], tmp_char);
				SWAP(buf_word[1], buf_word[2], tmp_char);
			}
			xdim = *(reinterpret_cast<int*>(buf_word));
			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word)); // Ydim
			if (is_little_endian)
			{
				SWAP(buf_word[0], buf_word[3], tmp_char);
				SWAP(buf_word[1], buf_word[2], tmp_char);
			}
			ydim = *(reinterpret_cast<int*>(buf_word));
			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word)); // Zdim
			if (is_little_endian)
			{
				SWAP(buf_word[0], buf_word[3], tmp_char);
				SWAP(buf_word[1], buf_word[2], tmp_char);
			}
			zdim = *(reinterpret_cast<int*>(buf_word));
			std::cout << " Binning factor = " << binning_factor << std::endl;
			std::cout << " Dimensions XYZ   (binned, unflipped coords) = " << xdim << " * " << ydim << " * " << zdim << std::endl;
			std::cout << " Dimensions XYZ (unbinned, unflipped coords) = " << xdim * binning_factor << " * " << ydim * binning_factor << " * " << zdim * binning_factor << std::endl;
			std::cout << " Dimensions XYZ           (3D tomogram .mrc) = " << xdim_img << " * " << ydim_img << " * " << zdim_img << std::endl;
			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word)); // Number of objects
			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word)); // Flags
			if (is_little_endian)
			{
				SWAP(buf_word[0], buf_word[3], tmp_char);
				SWAP(buf_word[1], buf_word[2], tmp_char);
			}
			if ((*(reinterpret_cast<int*>(buf_word))) & 0x00010000) // Check flag #16 - flip YZ?
				flip_YZ = true;
			//std::cout << (*(reinterpret_cast<int*>(buf_word))) << std::endl;
			std::cout << " Model last viewed on Y/Z flipped or rotated image? = " << std::flush;
			if (flip_YZ)
				std::cout << "TRUE" << std::endl;
			else
				std::cout << "FALSE" << std::endl;
			contain_3d_points = false;
			while (fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word))) // Read 4-byte blocks
			{
				if ( (buf_word[0] == 'C') && (buf_word[1] == 'O') && (buf_word[2] == 'N') && (buf_word[3] == 'T') ) // Find contour section
				{
					contain_3d_points = true;
					break;
				}
			}
			if (!contain_3d_points)
				REPORT_ERROR("helix.cpp::Interpolate3DCurves(): IMOD file does not seem to contain manually picked 3D coordiantes!");

			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word)); // Number of 3D points (meshes)
			if (is_little_endian)
			{
				SWAP(buf_word[0], buf_word[3], tmp_char);
				SWAP(buf_word[1], buf_word[2], tmp_char);
			}
			nr_points = *(reinterpret_cast<int*>(buf_word));
			//std::cout << nr_points << std::endl;
			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word));
			if (nr_points <= 2)
				REPORT_ERROR("helix.cpp::Interpolate3DCurves(): Input coordinate file: " + (std::string)(fn_in) + " should contain at least 2 points!");

			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word));
			fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word));
			std::cout << " Original XYZ coordinates (unbinned):" << std::endl;
			for (int id = 0; id < nr_points; id++)
			{
				fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word));
				if (is_little_endian)
				{
					SWAP(buf_word[0], buf_word[3], tmp_char);
					SWAP(buf_word[1], buf_word[2], tmp_char);
				}
				val = ((RFLOAT)(binning_factor)) * (*(reinterpret_cast<float*>(buf_word)));
				xlist.push_back(val);

				fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word));
				if (is_little_endian)
				{
					SWAP(buf_word[0], buf_word[3], tmp_char);
					SWAP(buf_word[1], buf_word[2], tmp_char);
				}
				val = ((RFLOAT)(binning_factor)) * (*(reinterpret_cast<float*>(buf_word)));
				if (flip_YZ)
					zlist.push_back(val);
				else
					ylist.push_back(val);

				fin.read(reinterpret_cast<char*>(buf_word), sizeof(buf_word));
				if (is_little_endian)
				{
					SWAP(buf_word[0], buf_word[3], tmp_char);
					SWAP(buf_word[1], buf_word[2], tmp_char);
				}
				val = ((RFLOAT)(binning_factor)) * (*(reinterpret_cast<float*>(buf_word)));
				if (flip_YZ)
					ylist.push_back(val);
				else
					zlist.push_back(val);

				// OLD
				//std::cout << "    " << xlist[xlist.size() - 1] << "        " << ylist[ylist.size() - 1] << "        " << zlist[zlist.size() - 1] << std::endl;
				// NEW
				std::cout << std::setw(15) << std::fixed << xlist[xlist.size() - 1];
				std::cout << std::setw(15) << std::fixed << ylist[ylist.size() - 1];
				std::cout << std::setw(15) << std::fixed << zlist[zlist.size() - 1];
				std::cout << std::endl;
			}
			//std::cout << xlist.size() << std::endl;
		}
		else if (fn_in.getExtension() == "star")
		{
			MetaDataTable MD;
			RFLOAT xx = 0., yy = 0., zz = 0.;
			MD.clear();
			MD.read(fn_in);
			if ( (!MD.containsLabel(EMDL_IMAGE_COORD_X))
					|| (!MD.containsLabel(EMDL_IMAGE_COORD_Y))
					|| (!MD.containsLabel(EMDL_IMAGE_COORD_Z)) )
				REPORT_ERROR("helix.cpp::Interpolate3DCurves(): Input coordinate STAR file " + (std::string)(fn_in) + " should contain _rlnCoordinateX Y and Z!");
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
			{
				MD.getValue(EMDL_IMAGE_COORD_X, xx);
				MD.getValue(EMDL_IMAGE_COORD_Y, yy);
				MD.getValue(EMDL_IMAGE_COORD_Z, zz);
				xlist.push_back(xx);
				ylist.push_back(yy);
				zlist.push_back(zz);
			}
		}
		else
		{
			fin.open(fn_in.c_str(), std::ios_base::in);
			if (fin.fail())
				REPORT_ERROR("helix.cpp::Interpolate3DCurves(): Cannot open coordinate file: " + (std::string)(fn_in));

			// Read x, y, z coordinates into vectors and close the input file
			while (getline(fin, line, '\n'))
			{
				words.clear();
				tokenize(line, words);
				if (words.size() == 0) // Empty line.
					continue;
				if (words.size() != 3)
					REPORT_ERROR("helix.cpp::Interpolate3DCurves(): Invalid input coordinate file " + fn_in);

				xlist.push_back(textToFloat(words[0]));
				ylist.push_back(textToFloat(words[1]));
				zlist.push_back(textToFloat(words[2]));

				// Screen output
				std::cout << " Original XYZ coordinates:  " << textToFloat(words[0]) << ", "
						<< textToFloat(words[1]) << ", " << textToFloat(words[2]) << std::endl;
			}
		}
		fin.close();
		if (xlist.size() < 2)
			REPORT_ERROR("helix.cpp::Interpolate3DCurves(): Input coordinate file: " + (std::string)(fn_in) + " should contain at least 2 points!");

		// Load control points for 3D interpolations
		MD_in.clear();
		MD_in.addLabel(EMDL_IMAGE_COORD_X);
		MD_in.addLabel(EMDL_IMAGE_COORD_Y);
		MD_in.addLabel(EMDL_IMAGE_COORD_Z);
		// Mode 1 - Just use the manually picked points as control points
		//for (int id = 0; id < xlist.size(); id++)
		//{
		//	MD_in.addObject();
		//	MD_in.setValue(EMDL_IMAGE_COORD_X, xlist[id]);
	    //	MD_in.setValue(EMDL_IMAGE_COORD_Y, ylist[id]);
	    //	MD_in.setValue(EMDL_IMAGE_COORD_Z, zlist[id]);
		//}
		// Mode 2 - Manually picked points are zigzag. Choose middle points of short line segments as control points.
		// Generate smooth curves
		// However, tilt angles of the segments around start- and end- points of the filaments deviate more from ~90 degrees.
		for (int id = 0; id < xlist.size(); id++)
		{
			if (id == 0) // Start point
			{
				MD_in.addObject();
				MD_in.setValue(EMDL_IMAGE_COORD_X, xlist[0]);
		    	MD_in.setValue(EMDL_IMAGE_COORD_Y, ylist[0]);
		    	MD_in.setValue(EMDL_IMAGE_COORD_Z, zlist[0]);
		    	continue;
			}
			// Middle points of each short line segment
			MD_in.addObject();
			MD_in.setValue(EMDL_IMAGE_COORD_X, (xlist[id] + xlist[id - 1]) / 2.);
	    	MD_in.setValue(EMDL_IMAGE_COORD_Y, (ylist[id] + ylist[id - 1]) / 2.);
	    	MD_in.setValue(EMDL_IMAGE_COORD_Z, (zlist[id] + zlist[id - 1]) / 2.);
	    	if (id == (xlist.size() - 1)) // End point
	    	{
				MD_in.addObject();
				MD_in.setValue(EMDL_IMAGE_COORD_X, xlist[id]);
		    	MD_in.setValue(EMDL_IMAGE_COORD_Y, ylist[id]);
		    	MD_in.setValue(EMDL_IMAGE_COORD_Z, zlist[id]);
		    	break;
	    	}
		}
		//std::cout << MD_in.numberOfObjects() << std::endl; // DEBUG
		//fn_out = fn_in.beforeFirstOf(fn_in_root) + fn_out_root + ".star"; // DEBUG
		//MD_in.write(fn_out); // DEBUG

		// Interpolate
		HermiteInterpolateOne3DHelicalFilament(
				MD_in,
				MD_out,
				nr_segments,
				nr_asu,
				rise_A,
				pixel_size_A,
				box_size_pix,
				fid + 1,
				xdim_img, ydim_img, zdim_img,
				bimodal_angular_priors);

		// Output
	    if (MD_out.numberOfObjects() < 1)
	    	std::cout << " WARNING: No sub-tomograms have been interpolated on this helical filament!" << std::endl;
	    else
	    	MD_all.append(MD_out);
		total_segments += nr_segments;
		std::cout << " Interpolated " << nr_segments <<  " helical segments from 3D point set " << fn_in << std::endl;
		std::cout << " ========================================================= " << std::endl;
	}

	//fn_out = fn_in.beforeFirstOf(fn_in_root) + fn_out_root + ".star"; // OLD
	fn_out = fn_in_root + fn_out_root + ".star"; // NEW
	// DEBUG
	//std::cout << " fn_in = " << fn_in << std::endl;
	//std::cout << " fn_in_root = " << fn_in_root << std::endl;
	//std::cout << " fn_out_root = " << fn_out_root << std::endl;
	//std::cout << " fn_out = " << fn_out << std::endl;

	if (MD_all.numberOfObjects() < 1)
		std::cout << " ### Done! WARNING: No sub-tomograms have been interpolated! Please check whether you have done everything correctly." << std::endl;
	else
	{
		MD_all.write(fn_out);
		std::cout << " ### Done! Interpolated " << total_segments << " helical segments on " << fn_in_list.size()
				<< " filaments. Output file: " << fn_out << std::endl;
	}
}

void estimateTiltPsiPriors(
		Matrix1D<RFLOAT>& dr,
		RFLOAT& tilt_deg,
		RFLOAT& psi_deg)
{
	// euler.cpp: Euler_direction2angles: input angles = (a, b, g) then 3x3 matrix =
	//  cosg*cosb*cosa - sing*sina,  cosg*cosb*sina + sing*cosa, -cosg*sinb,
	// -sing*cosb*cosa - cosg*sina, -sing*cosb*sina + cosg*cosa,  sing*sinb,
	//       sinb*cosa,                               sinb*sina,       cosb.

	// euler.cpp: Euler_direction2angles: input angles = (0, b, g) then 3x3 matrix =
	//  cosg*cosb, sing, -cosg*sinb,
	// -sing*cosb, cosg,  sing*sinb,
	//       sinb,    0,       cosb.

	RFLOAT tilt_rad = 0., psi_rad = 0., vec_len = 0.;
	int dim = dr.size();

	if ( (dim != 2) && (dim != 3) )
		REPORT_ERROR("helix.cpp::estimateTiltPsiPriors(): Input Matrix1D should have a size of 2 or 3!");

	vec_len = XX(dr) * XX(dr) + YY(dr) * YY(dr);
	vec_len += (dim == 3) ? (ZZ(dr) * ZZ(dr)) : (0.);
	vec_len = sqrt(vec_len);
	if (vec_len < 0.0001)
		REPORT_ERROR("helix.cpp::estimateTiltPsiPriors(): Vector length is smaller than 0.0001!");

	// A * (0, 0, z) = (x', y', z')
	// x' = -z*cosg*sinb
	// y' =  z*sing*sinb
	// z' =       z*cosb
	// cosb = z' / z
	// tang = y' / (-x')
	if (dim == 3)
	{
		// Tilt (b) should be [0, +180] degrees. Psi (g) should be [-180, +180] degrees
		tilt_rad = acos(ZZ(dr) / vec_len); // 'acos' returns an angle within [0, +pi] radians for tilt
		psi_rad = atan2(YY(dr), (-1.) * XX(dr)); // 'atan2' returns an angle within [-pi, +pi] radians for rot
	}
	else
		psi_rad = (-1.) * atan2(YY(dr), XX(dr));

	if (dim == 3)
		tilt_deg = RAD2DEG(tilt_rad);
	psi_deg = RAD2DEG(psi_rad);
}

void readFileHeader(
		FileName& fn_in,
		FileName& fn_out,
		int nr_bytes)
{
	std::ifstream fin;
	std::ofstream fout;
	int nr_blocks = 0, curr_block = 0;
	char data[100];

	if (nr_bytes > 10 * 1024 * 1024)
		REPORT_ERROR("helix.cpp::readFileHeader(): Don't copy more than 10MB data!");
	fin.open(fn_in.c_str(), std::ios_base::in|std::ios_base::binary);
	if (fin.fail())
		REPORT_ERROR("helix.cpp::readFileHeader(): Cannot open input file: " + (std::string)(fn_in));
	fout.open(fn_out.c_str(), std::ios_base::out|std::ios_base::binary);
	if (fout.fail())
		REPORT_ERROR("helix.cpp::readFileHeader(): Cannot open output file: " + (std::string)(fn_out));

	nr_blocks = nr_bytes / 100;
	nr_blocks = (nr_blocks < 1) ? (1) : (nr_blocks);
	std::cout << " Copying the first " << nr_blocks * 100 << " bytes from " << fn_in << " to " << fn_out << " ..." << std::endl;

	curr_block = 0;
	while (fin.read(reinterpret_cast<char*>(data), sizeof(data)))
	{
		curr_block++;
		fout.write(reinterpret_cast<char*>(data), sizeof(data));
		if (curr_block >= nr_blocks)
			break;
	}

	fin.close();
	fout.close();
}

void select3DsubtomoFrom2Dproj(
		MetaDataTable& MD_2d,
		MetaDataTable& MD_3d,
		MetaDataTable& MD_out)
{
	//std::vector<RFLOAT> xlist, ylist, zlist, idlist;
	std::vector<std::string> mic_list, img_list;
	int id = 0;
	//RFLOAT xx = 0., yy = 0., zz = 0.;
	FileName mic_str, img_str;
	const size_t id_length = 6; // 6-digit ID
	bool auto_pixel_size = false;
	RFLOAT Dpix = -1., Mag = -1., _Dpix = -1., _Mag = -1.;

	MD_out.clear();

	if (MD_2d.numberOfObjects() < 1)
		REPORT_ERROR("helix.cpp::select3DsubtomoFrom2Dproj(): MetaDataTable 2D projections is empty!");
	if ( (!MD_2d.containsLabel(EMDL_IMAGE_COORD_X))
			|| (!MD_2d.containsLabel(EMDL_IMAGE_COORD_Y))
			|| (!MD_2d.containsLabel(EMDL_IMAGE_COORD_Z))
			|| (!MD_2d.containsLabel(EMDL_MICROGRAPH_NAME))
			|| (!MD_2d.containsLabel(EMDL_IMAGE_NAME)) )
		REPORT_ERROR("helix.cpp::select3DsubtomoFrom2Dproj(): MetaDataTable 2D projections should contain labels _rlnCoordinateXYZ, _rlnMicrographName and _rlnImageName!");
	// For particle rescaling  - Does MD_2d contain Dpix and Magnification?
	auto_pixel_size = (MD_2d.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE)) && (MD_2d.containsLabel(EMDL_CTF_MAGNIFICATION));

	if (MD_3d.numberOfObjects() < 1)
		REPORT_ERROR("helix.cpp::select3DsubtomoFrom2Dproj(): MetaDataTable 3D subtomograms is empty!");
	if ( (!MD_3d.containsLabel(EMDL_IMAGE_COORD_X))
			|| (!MD_3d.containsLabel(EMDL_IMAGE_COORD_Y))
			|| (!MD_3d.containsLabel(EMDL_IMAGE_COORD_Z))
			|| (!MD_3d.containsLabel(EMDL_MICROGRAPH_NAME))
			|| (!MD_3d.containsLabel(EMDL_IMAGE_NAME)) )
		REPORT_ERROR("helix.cpp::select3DsubtomoFrom2Dproj(): MetaDataTable 3D subtomograms should contain labels _rlnCoordinateXYZ, _rlnMicrographName and _rlnImageName!");
	if (auto_pixel_size)
	{
		if ( (!MD_3d.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE)) || (!MD_3d.containsLabel(EMDL_CTF_MAGNIFICATION)) )
			REPORT_ERROR("helix.cpp::select3DsubtomoFrom2Dproj(): If MetaDataTable 2D projections contains Dpix and Magnification, then MetaDataTable 3D subtomograms should also do!");

		// Firstly, check whether the pixel sizes of segments in MD_3d are the same
		Dpix = -1., Mag = -1., _Dpix = -1., _Mag = -1.;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_3d)
		{
			MD_3d.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, _Dpix);
			MD_3d.getValue(EMDL_CTF_MAGNIFICATION, _Mag);
			if ( (!(_Dpix > 0.)) || (!(_Mag > 0.)) )
				REPORT_ERROR("helix.cpp::select3DsubtomoFrom2Dproj(): Please ensure that all entries in MetaDataTable 3D subtomograms have valid Dpix and Magnification!");
			if ( (Dpix < 0.) && (Mag < 0.) )
			{
				Dpix = _Dpix;
				Mag = _Mag;
				continue;
			}
			if ( (fabs(_Dpix - Dpix) > 0.001) || (fabs(_Mag - Mag) > 0.001) )
				REPORT_ERROR("helix.cpp::select3DsubtomoFrom2Dproj(): Please ensure that all entries in MetaDataTable 3D subtomograms have the same pixel size!");
		}
		std::cout << " Dpix and Magnification of 3D sub-tomograms: " << Dpix << " and " << Mag << std::endl;

		// Then, check whether the pixel sizes of segments in MD_2d are the same. Store Dpix and Mag.
		Dpix = -1., Mag = -1., _Dpix = -1., _Mag = -1.;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_2d)
		{
			MD_2d.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, _Dpix);
			MD_2d.getValue(EMDL_CTF_MAGNIFICATION, _Mag);
			if ( (!(_Dpix > 0.)) || (!(_Mag > 0.)) )
				REPORT_ERROR("helix.cpp::select3DsubtomoFrom2Dproj(): Please ensure that all entries in MetaDataTable 2D projections have valid Dpix and Magnification!");
			if ( (Dpix < 0.) && (Mag < 0.) )
			{
				Dpix = _Dpix;
				Mag = _Mag;
				continue;
			}
			if ( (fabs(_Dpix - Dpix) > 0.001) || (fabs(_Mag - Mag) > 0.001) )
				REPORT_ERROR("helix.cpp::select3DsubtomoFrom2Dproj(): Please ensure that all entries in MetaDataTable 2D projections have the same pixel size!");
		}
		std::cout << " Dpix and Magnification of 2D projections:   " << Dpix << " and " << Mag << std::endl;
		std::cout << " Reset Dpix and Magnification of selected 3D sub-tomograms..." << std::endl;
		// Dpix and Mag have been stored.
	}

	// Gather all the selected subtomograms
	//xlist.clear(); ylist.clear(); zlist.clear(); idlist.clear();
	mic_list.clear();
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_2d)
	{
		//MD_2d.getValue(EMDL_IMAGE_COORD_X, xx);
		//MD_2d.getValue(EMDL_IMAGE_COORD_Y, yy);
		//MD_2d.getValue(EMDL_IMAGE_COORD_Z, zz);
		MD_2d.getValue(EMDL_MICROGRAPH_NAME, mic_str);
		MD_2d.getValue(EMDL_IMAGE_NAME, img_str);

		//id = textToInteger(img_str.beforeFirstOf("@")); // 6-digit ID
		//xlist.push_back(xx); ylist.push_back(yy); zlist.push_back(zz); idlist.push_back(id);
		mic_list.push_back(mic_str + img_str.beforeFirstOf("@"));
	}
	std::sort(mic_list.begin(), mic_list.end());

	// Scan through all the subtomograms
	MD_out.clear();
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_3d)
	{
		//MD_3d.getValue(EMDL_IMAGE_COORD_X, xx);
		//MD_3d.getValue(EMDL_IMAGE_COORD_Y, yy);
		//MD_3d.getValue(EMDL_IMAGE_COORD_Z, zz);
		MD_3d.getValue(EMDL_MICROGRAPH_NAME, mic_str);
		MD_3d.getValue(EMDL_IMAGE_NAME, img_str);

		img_str = img_str.withoutExtension();
		if (img_str.length() < id_length)
			REPORT_ERROR("helix.cpp::select3DsubtomoFrom2Dproj(): img_str.length() < " + integerToString(id_length) + " ! img_str = " + (std::string)(img_str));
		img_str = img_str.substr(img_str.length() - id_length, id_length);
		//std::cout << "  " << img_str << std::flush;
		img_str = mic_str + img_str;
		if (std::binary_search(mic_list.begin(), mic_list.end(), img_str)) // Subtomogram selected
		{
			MD_out.addObject(MD_3d.getObject());
			// For particle rescaling - reset pixel size
			if (auto_pixel_size)
			{
				MD_out.setValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, Dpix);
				MD_out.setValue(EMDL_CTF_MAGNIFICATION, Mag);
			}
		}
	}
	//std::cout << std::endl;
}

void averageAsymmetricUnits2D(
		ObservationModel& obsModel,
		MetaDataTable &MDimgs,
		FileName fn_o_root,
		int nr_asu,
		RFLOAT rise)
{

	if (nr_asu == 1)
	{
		std::cout << " WARNING: averageAsymmetricUnits2D nr_asu=1, so not doing anything ...";
		return;
	}

	int nr_asu_half = nr_asu / 2;
	RFLOAT angpix;

	FourierTransformer transformer;
	MultidimArray<Complex> Fimg, Faux, Fsum;


	long int imgno = 0;
	init_progress_bar(MDimgs.numberOfObjects());
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDimgs)
	{

		FileName fn_img;
		Matrix1D<RFLOAT> in(2), out(2);
		Image<RFLOAT> img;
		RFLOAT psi, angpix;
		int optics_group;

		MDimgs.getValue(EMDL_IMAGE_NAME, fn_img);
		MDimgs.getValue(EMDL_ORIENT_PSI, psi);
		MDimgs.getValue(EMDL_IMAGE_OPTICS_GROUP, optics_group);
		optics_group--;
		angpix = obsModel.getPixelSize(optics_group);

		img.read(fn_img);
		transformer.FourierTransform(img(), Fimg, false);
		Fsum = Fimg; // original image

		//std::cerr << " imgno= " << imgno << " fn_img= " << fn_img << " psi= " << psi << " rise= " << rise  << " angpix= " << angpix << " nr_asu= " << nr_asu << " xsize= " << XSIZE(img()) << std::endl;

		for (int i = 2; i <= nr_asu; i++)
		{
			// one way
			if (i%2 == 0)
			{
				XX(in) = rise * (i/2) / angpix;
			}
			// the other way
			else
			{
				XX(in) = -1. * rise * (i/2) / angpix;
			}
			YY(in) = 0.;

			transformCartesianAndHelicalCoords(in, out, 0., 0., psi, HELICAL_TO_CART_COORDS);
			//std::cerr << " i= " << i << " XX(in)= " << XX(in) << " YY(in)= " << YY(in) << " XX(out)= " << XX(out) << " YY(out)= " << YY(out)   << std::endl;
			shiftImageInFourierTransform(Fimg, Faux, XSIZE(img()), XX(out), YY(out));
			Fsum += Faux;
		}


		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
		{
			DIRECT_MULTIDIM_ELEM(Fimg, n) = DIRECT_MULTIDIM_ELEM(Fsum, n) / (RFLOAT)nr_asu;
		}
		transformer.inverseFourierTransform();

		// Write this particle to the stack on disc
		// First particle: write stack in overwrite mode, from then on just append to it
		MDimgs.setValue(EMDL_IMAGE_ORI_NAME, fn_img);
		fn_img.compose(imgno+1, fn_o_root + "particles.mrcs");
		if (imgno == 0)
			img.write(fn_img, -1, false, WRITE_OVERWRITE);
		else
			img.write(fn_img, -1, false, WRITE_APPEND);
		MDimgs.setValue(EMDL_IMAGE_NAME, fn_img);

		if (imgno%60==0) progress_bar(imgno);
		imgno++;
	}
	progress_bar(MDimgs.numberOfObjects());

}
