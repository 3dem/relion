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

#include "src/helix.h"

bool IsEqualHelicalSymmetry(
		const HelicalSymmetryItem& a,
		const HelicalSymmetryItem& b)
{
	if ( (fabs(a.twist_deg - b.twist_deg) < (1e-5))
			&& (fabs(a.rise_pix - b.rise_pix) < (1e-5)) )
		return true;
	return false;
};

bool IsHelicalSymmetryWithinOpenInterval(
		RFLOAT rise_pix,
		RFLOAT twist_deg,
		RFLOAT rise_ori_pix,
		RFLOAT twist_ori_deg,
		RFLOAT rise_half_range_pix,
		RFLOAT twist_half_range_deg)
{
	RFLOAT rise_min_pix, twist_min_deg, rise_max_pix, twist_max_deg;
	RFLOAT rise_error_pix, twist_error_deg;

	rise_pix = fabs(rise_pix);
	twist_deg = fabs(twist_deg);
	rise_half_range_pix = fabs(rise_half_range_pix);
	twist_half_range_deg = fabs(twist_half_range_deg);

	rise_error_pix = rise_pix * (1e-5);
	twist_error_deg = twist_deg * (1e-5);

	rise_min_pix = rise_pix - rise_half_range_pix + rise_error_pix;
	rise_max_pix = rise_pix + rise_half_range_pix - rise_error_pix;
	twist_min_deg = twist_deg - twist_half_range_deg + twist_error_deg;
	twist_max_deg = twist_deg + twist_half_range_deg - twist_error_deg;

	if ( (rise_pix < rise_min_pix)
			|| (rise_pix > rise_max_pix)
			|| (twist_deg < twist_min_deg)
			|| (twist_deg > twist_max_deg) )
	{
		return false;
	}
	return true;
};

void makeHelicalSymmetryList(
		std::vector<HelicalSymmetryItem>& list,
		RFLOAT rise_ori_pix,
		RFLOAT twist_ori_deg,
		RFLOAT rise_half_range_pix,
		RFLOAT twist_half_range_deg,
		RFLOAT rise_step_pix,
		RFLOAT twist_step_deg)
{
	// Assume all parameters are within range
	RFLOAT rise_pix, twist_deg;
	int rise_half_samplings, twist_half_samplings;
	std::vector<HelicalSymmetryItem> tmp_list;

	rise_half_samplings = ROUND(fabs(rise_half_range_pix) / fabs(rise_step_pix));
	twist_half_samplings = ROUND(fabs(twist_half_range_deg) / fabs(twist_step_deg));

	// Store a matrix of symmetries
	tmp_list.clear();
	for (int ii = -rise_half_samplings; ii <= rise_half_samplings; ii++)
	{
		for (int jj = -twist_half_samplings; jj <= twist_half_samplings; jj++)
		{
			rise_pix = rise_ori_pix + rise_step_pix * RFLOAT(ii);
			twist_deg = twist_ori_deg + twist_step_deg * RFLOAT(jj);
			tmp_list.push_back(HelicalSymmetryItem(twist_deg, rise_pix));
		}
	}

	// Check duplications and return this matrix
	for (int ii = 0; ii < list.size(); ii++)
	{
		for (int jj = 0; jj < tmp_list.size(); jj++)
		{
			if (IsEqualHelicalSymmetry(list[ii], tmp_list[jj]))
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
	int r_max_XY, startZ, finishZ;
	RFLOAT twist_rad, twist_rad_total, dist_r_pix, RFLOAT_counter, sum_pw1, sum_pw2;
	std::vector<RFLOAT> dev_voxel, dev_chunk;
	int id, idx, idy, idz, x0, y0, z0, x1, y1, z1;
	RFLOAT xp, yp, zp, fx, fy, fz, d000, d001, d010, d011, d100, d101, d110, d111, dx00, dx01, dx10, dx11, dxy0, dxy1, ddd;

	if ( (STARTINGZ(v) != FIRST_XMIPP_INDEX(ZSIZE(v)))
			|| (STARTINGY(v) != FIRST_XMIPP_INDEX(YSIZE(v)))
			|| (STARTINGX(v) != FIRST_XMIPP_INDEX(XSIZE(v))) )
	{
		REPORT_ERROR("helix.cpp::calcCCofHelicalSymmetry(): The origin of input 3D MultidimArray is not at the center (use v.setXmippOrigin() before calling this function)!");
		return false;
	}

	// Check r_max
	r_max_XY = (XSIZE(v) < YSIZE(v)) ? XSIZE(v) : YSIZE(v);
	r_max_XY = (r_max_XY + 1) / 2 - 1;
	if( r_max_pix > (((RFLOAT)(r_max_XY)) - 0.01) )  // 0.01 - avoid segmentation fault
		r_max_pix = (((RFLOAT)(r_max_XY)) - 0.01);

	// Set startZ and finishZ
	startZ = int(floor((-1.) * ((RFLOAT)(ZSIZE(v)) * z_percentage * 0.5)));
	finishZ = int(ceil(((RFLOAT)(ZSIZE(v)) * z_percentage * 0.5)));
	startZ = (startZ <= (STARTINGZ(v))) ? (STARTINGZ(v) + 1) : (startZ);
	finishZ = (finishZ >= (FINISHINGZ(v))) ? (FINISHINGZ(v) - 1) : (finishZ);

	// Get twist angle in radians
	twist_rad = twist_deg * PI / 180.;

	// Test a chunk of Z length = rise
	dev_chunk.clear();
	// Iterate through all coordinates on Z, Y and then X axes
	FOR_ALL_ELEMENTS_IN_ARRAY3D(v)
	{
		// Record X, Y, Z coordinates
		idz = k;
		idy = i;
		idx = j;

		// Test a chunk of Z length = rise
		// for(idz = startZ; (idz <= (startZ + ((int)(floor(rise_pix))))) && (idz <= finishZ); idz++)
		if ( (idz < startZ) || (idz > (startZ + ((int)(floor(rise_pix))))) || (idz > finishZ) )
			continue;

		dist_r_pix = sqrt(idy * idy + idx * idx);
		if( (dist_r_pix < r_min_pix) || (dist_r_pix > r_max_pix) )
			continue;

		// Pick a voxel in the chunk
		dev_voxel.clear();
		dev_voxel.push_back(A3D_ELEM(v, idz, idy, idx));

		// Pick other voxels according to this voxel and helical symmetry
		zp = idz;
		RFLOAT_counter = 0.;
		while(1)
		{
			// Rise
			zp += rise_pix;
			if(zp > finishZ) // avoid segmentation fault - finishZ is always strictly smaller than FINISHINGZ(v)!
				break;

			// Twist
			RFLOAT_counter += 1.;
			twist_rad_total = twist_rad * RFLOAT_counter;
			xp = ((RFLOAT)(idx)) * cos(twist_rad_total) - ((RFLOAT)(idy)) * sin(twist_rad_total);
			yp = ((RFLOAT)(idx)) * sin(twist_rad_total) + ((RFLOAT)(idy)) * cos(twist_rad_total);

			// Trilinear interpolation (with physical coords)
			// Subtract STARTINGX,Y,Z to accelerate access to data
			// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
			x0 = FLOOR(xp); fx = xp - x0; x0 -= STARTINGX(v); x1 = x0 + 1;
			y0 = FLOOR(yp); fy = yp - y0; y0 -= STARTINGY(v); y1 = y0 + 1;
			z0 = FLOOR(zp); fz = zp - z0; z0 -= STARTINGZ(v); z1 = z0 + 1;
			// DEBUG
			if( (x0 < 0) || (y0 < 0) || (z0 < 0)
					|| (x1 >= XSIZE(v)) || (y1 >= YSIZE(v)) || (z1 >= ZSIZE(v)) )
			{
				std::cout << "idzidyidx = " << idz << ", " << idy << ", " << idx
						<< ",       x0x1y0y1z0z1 = " << x0 << ", " << x1 << ", " << y0
						<< ", " << y1 << ", " << z0 << ", " << z1 << std::endl;
			}

			// Get values
			d000 = DIRECT_A3D_ELEM(v, z0, y0, x0);
			d001 = DIRECT_A3D_ELEM(v, z0, y0, x1);
			d010 = DIRECT_A3D_ELEM(v, z0, y1, x0);
			d011 = DIRECT_A3D_ELEM(v, z0, y1, x1);
			d100 = DIRECT_A3D_ELEM(v, z1, y0, x0);
			d101 = DIRECT_A3D_ELEM(v, z1, y0, x1);
			d110 = DIRECT_A3D_ELEM(v, z1, y1, x0);
			d111 = DIRECT_A3D_ELEM(v, z1, y1, x1);

			// Interpolation 3D -> 2D
			dx00 = LIN_INTERP(fx, d000, d001);
			dx01 = LIN_INTERP(fx, d100, d101);
			dx10 = LIN_INTERP(fx, d010, d011);
			dx11 = LIN_INTERP(fx, d110, d111);

			// Interpolation 2D -> 1D
			dxy0 = LIN_INTERP(fy, dx00, dx10);
			dxy1 = LIN_INTERP(fy, dx01, dx11);

			// Interpolation of 2 voxels
			ddd = LIN_INTERP(fz, dxy0, dxy1);

			// Record this voxel
			dev_voxel.push_back(ddd);
		}

		// Calc dev of this voxel in the chunk
		if(dev_voxel.size() > 1)
		{
			sum_pw1 = sum_pw2 = 0.;
			for(id = 0; id < dev_voxel.size(); id++)
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
	}

	// Calc avg of all voxels' devs in this chunk (for a specific helical symmetry)
	if(dev_chunk.size() < 1)
	{
		cc = (1e10);
		nr_asym_voxels = 0;
		return false;
	}
	else
	{
		sum_pw1 = 0.;
		for(id = 0; id < dev_chunk.size(); id++)
			sum_pw1 += dev_chunk[id];
		cc = (sum_pw1 / dev_chunk.size());
	}
	nr_asym_voxels = dev_chunk.size();
	dev_chunk.clear();

	return true;
};

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
		RFLOAT& twist_refined_deg)
{
	// TODO: whether iterations can exit & this function works for negative twist
	int iter, box_len, nr_asym_voxels, nr_half_samplings, nr_min_half_samplings, best_id;
	RFLOAT r_min_pix, r_max_pix, rise_ori_pix, best_dev;
	RFLOAT rise_half_range_pix, rise_step_pix, twist_half_range_deg, twist_step_deg;
	std::vector<HelicalSymmetryItem> helical_symmetry_list;

	if (v.getDim() != 3)
	{
		REPORT_ERROR("helix.cpp::localSearchHelicalSymmetry(): Input helical reference is not 3D! (v.getDim() = "
						+ integerToString(v.getDim()) + ")");
		return false;
	}
	box_len = (XSIZE(v) < YSIZE(v)) ? XSIZE(v) : YSIZE(v);
	box_len = (box_len < ZSIZE(v)) ? box_len : ZSIZE(v);
	if (checkHelicalParametersFor3DHelicalReference(
			box_len,
			pixel_size_A,
			twist_ori_deg,
			rise_ori_A,
			z_percentage,
			true,
			rise_search_max_dev_percentage,
			twist_search_max_dev_percentage,
			sphere_radius_A,
			r_min_A,
			r_max_A) == false)
	{
		REPORT_ERROR("helix.cpp::localSearchHelicalSymmetry(): Input Parameters error!");
		return false;
	}

	r_min_pix = r_min_A / pixel_size_A;
	r_max_pix = r_max_A / pixel_size_A;
	rise_ori_pix = rise_ori_A / pixel_size_A;

	// Initial searches - Iteration 1
	// Sampling steps should be finer than 1% the twist and the rise
	// And also make sure to search for at least 11*11 sampling points (5*2+1=11)
	nr_min_half_samplings = 5;

	rise_half_range_pix = fabs(rise_ori_pix) * rise_search_max_dev_percentage;
	nr_half_samplings = CEIL(rise_search_max_dev_percentage / 0.01);
	nr_half_samplings = (nr_half_samplings >= nr_min_half_samplings) ? (nr_half_samplings) : (nr_min_half_samplings);
	rise_step_pix = rise_half_range_pix / RFLOAT(nr_half_samplings);

	twist_half_range_deg = fabs(twist_ori_deg) * twist_search_max_dev_percentage;
	nr_half_samplings = CEIL(twist_search_max_dev_percentage / 0.01);
	nr_half_samplings = (nr_half_samplings >= nr_min_half_samplings) ? (nr_half_samplings) : (nr_min_half_samplings);
	twist_step_deg = twist_half_range_deg / RFLOAT(nr_half_samplings);

	helical_symmetry_list.clear();
	for (iter = 1; iter <= 20; iter++)
	{
		makeHelicalSymmetryList(
				helical_symmetry_list,
				rise_ori_pix,
				twist_ori_deg,
				rise_half_range_pix,
				twist_half_range_deg,
				rise_step_pix,
				twist_step_deg);

		if (helical_symmetry_list.size() < 1)
			REPORT_ERROR("helix.cpp::localSearchHelicalSymmetry(): No helical symmetries are found in the search list!");
		best_dev = (1e30);
		best_id = -1;
		for (int ii = 0; ii < helical_symmetry_list.size(); ii++)
		{
			// If this symmetry is not calculated before
			if (helical_symmetry_list[ii].dev > (1e30))
			{
				calcCCofHelicalSymmetry(
						v,
						r_min_pix,
						r_max_pix,
						z_percentage,
						helical_symmetry_list[ii].rise_pix,
						helical_symmetry_list[ii].twist_deg,
						helical_symmetry_list[ii].dev,
						nr_asym_voxels);
				// DEBUG
				//std::cout << helical_symmetry_list[ii].rise_pix << "	" << helical_symmetry_list[ii].twist_deg << "	" << helical_symmetry_list[ii].dev << std::endl;
			}
			if (helical_symmetry_list[ii].dev < best_dev)
			{
				best_dev = helical_symmetry_list[ii].dev;
				best_id = ii;
			}

			// DEBUG
			//std::cout << "  Twist = " << helical_symmetry_list[ii].twist_deg
			//		<< ", Rise = " << helical_symmetry_list[ii].rise_pix * pixel_size_A
			//		<< ", Dev = " << helical_symmetry_list[ii].dev << std::endl;
		}


		// Update refined symmetry
		rise_ori_pix = helical_symmetry_list[best_id].rise_pix;
		rise_refined_A = rise_ori_pix * pixel_size_A;
		twist_refined_deg = twist_ori_deg = helical_symmetry_list[best_id].twist_deg;

		// DEBUG
		//std::cout << " ######################################################## " << std::endl;
		//std::cout << "  ##### Refined Twist = " << twist_refined_deg << ", Rise = " << rise_refined_A
		//		<< ", Dev = " << helical_symmetry_list[best_id].dev << std::endl;
		//std::cout << " ######################################################## " << std::endl;

		if (IsHelicalSymmetryWithinOpenInterval(
				helical_symmetry_list[best_id].rise_pix,
				helical_symmetry_list[best_id].twist_deg,
				rise_ori_pix,
				twist_ori_deg,
				rise_half_range_pix,
				twist_half_range_deg) == false)
		{
			if (iter == 1)
			{
				std::cout << " WARNING: Refined helical symmetry is out of the search range. Check whether the initial helical symmetry is reasonable. Or you may want to modify the search range." << std::endl;
				return false;
			}
			else
			{
				// Finer searches, optimal symmetry is out of the current search range
				// Expand the search range and try again
				rise_half_range_pix *= (2.);
				twist_half_range_deg *= (2.);

				// Prevent extra large search range
				if ( (rise_half_range_pix / rise_step_pix) > 15)
				{
					std::cout << " WARNING: Local searches failed to converge in finding helical symmetry." << std::endl;
					return false;
				}
			}
		}
		else
		{
			// Set 5*5 finer samplings for the next iteration
			rise_half_range_pix = rise_step_pix;
			twist_half_range_deg = twist_step_deg;

			rise_step_pix /= (2.);
			twist_step_deg /= (2.);

			// Stop searches if step sizes are too small
			if ( (rise_step_pix < fabs(rise_ori_pix) * (1e-5)) || (twist_step_deg < fabs(twist_ori_deg) * (1e-5)) )
				break;
		}
	}

	return true;
};





/*
bool get2DZsliceIn3DVolume(
		const MultidimArray<RFLOAT>& vol_in,
		MultidimArray<RFLOAT>& img_out,
		int idz)
{
	int Xdim, Ydim, Zdim, Ndim;

	img_out.clear();
	Xdim = XSIZE(vol_in); Ydim = YSIZE(vol_in); Zdim = ZSIZE(vol_in); Ndim = NSIZE(vol_in);

	if( (Ndim != 1) || (Xdim < 5) || (Ydim < 5) || (Zdim < 5) )
	{
		REPORT_ERROR("helix.cpp::get2DZsliceIn3DVolume(): Input 3D MultidimArray has Wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");
		return false;
	}

	if( (idz < STARTINGZ(vol_in)) || (idz > FINISHINGZ(vol_in)) )
	{
		REPORT_ERROR("helix.cpp::get2DZsliceIn3DVolume(): Invalid Z slice! (z = [" + integerToString(STARTINGZ(vol_in)) + ", " + integerToString(FINISHINGZ(vol_in)) + "] while requesting to access idz = " + integerToString(idz) + ")");
		return false;
	}

	img_out.initZeros(Ydim, Xdim);
	FOR_ALL_ELEMENTS_IN_ARRAY2D(img_out)
	{
		A2D_ELEM(img_out, i, j)
				= A3D_ELEM(vol_in, idz, i - STARTINGY(img_out) + STARTINGY(vol_in), j - STARTINGX(img_out) + STARTINGX(vol_in));
	}
	img_out.setXmippOrigin();
	return true;
};

bool add2DZsliceInto3DVolumeSum(
		const MultidimArray<RFLOAT>& img_in,
		MultidimArray<RFLOAT>& vol_sum,
		std::vector<RFLOAT>& weight_sum,
		int idz)
{
	int Xdim_vol, Ydim_vol, Zdim_vol, Ndim_vol, Xdim_img, Ydim_img, Zdim_img, Ndim_img;

	Xdim_img = XSIZE(img_in); Ydim_img = YSIZE(img_in); Zdim_img = ZSIZE(img_in); Ndim_img = NSIZE(img_in);
	Xdim_vol = XSIZE(vol_sum); Ydim_vol = YSIZE(vol_sum); Zdim_vol = ZSIZE(vol_sum); Ndim_vol = NSIZE(vol_sum);

	if( (Ndim_img != 1) || (Zdim_img != 1) || (Ydim_img < 5) || (Xdim_img < 5)
			|| (Ndim_vol != 1) || (Zdim_vol < 5) || (Ydim_vol < 5) || (Xdim_vol < 5)
			|| (Ydim_img != Ydim_vol) || (Xdim_img != Xdim_vol)
			|| (weight_sum.size() != Zdim_vol) )
	{
		REPORT_ERROR("Wrong dimensions!");
		return false;
	}

	if( (idz < STARTINGZ(vol_sum)) || (idz > FINISHINGZ(vol_sum)) )
	{
		REPORT_ERROR("helix.cpp::add2DZsliceInto3DVolumeSum(): Invalid Z slice! (z = [" + integerToString(STARTINGZ(vol_sum)) + ", " + integerToString(FINISHINGZ(vol_sum)) + "] while requesting to access idz = " + integerToString(idz) + ")");
		return false;
	}

	for (long int i = STARTINGY(vol_sum); i <= FINISHINGY(vol_sum); i++)
	{
		for (long int j = STARTINGX(vol_sum); j <= FINISHINGX(vol_sum); j++)
		{
			A3D_ELEM(vol_sum, idz, i, j)
					+= A2D_ELEM(img_in, i - STARTINGY(vol_sum) + STARTINGY(img_in), j - STARTINGX(vol_sum) + STARTINGX(img_in));
		}
	}
	weight_sum[idz - STARTINGZ(vol_sum)] += 1.0;
	return true;
};

void shift3DVolumeAlongZAxisInFourierSpace(
		MultidimArray<RFLOAT>& img,
		RFLOAT shift_pix)
{
	const RFLOAT pi = 3.141592653589793238462643383279502884197;
	int Xdim, Ydim, Zdim, Ndim;
	RFLOAT tmp_RFLOAT, a, b, c, d, ac, bd, ab_cd;
	MultidimArray<Complex> Faux;
	FourierTransformer transformer;
	Xdim = XSIZE(img); Ydim = YSIZE(img); Zdim = ZSIZE(img); Ndim = NSIZE(img);

	if( (Ndim != 1) || (Xdim < 5) || (Ydim < 5) || (Zdim < 5) )
	{
		REPORT_ERROR("helix.cpp::shift3DVolumeAlongZAxisInFT(): Input 3D MultidimArray has wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");
		return;
	}

	Faux.clear();
	transformer.clear();

	// No need to shift
	if(fabs(shift_pix) < (1e-6))
	{
		return;
	}

	// Do FT
	CenterFFT(img, true);
	transformer.FourierTransform(img, Faux);

	// Shift
	shift_pix /= (RFLOAT)(-Zdim);
	for (long int k = 0, kp = 0; k < ZSIZE(Faux); k++, kp = (k < XSIZE(Faux)) ? k : k - ZSIZE(Faux))
	{
		tmp_RFLOAT = 2 * pi * (kp * shift_pix);
		a = cos(tmp_RFLOAT);
		b = sin(tmp_RFLOAT);
		for (long int i = 0, ip = 0 ; i < YSIZE(Faux); i++, ip = (i < XSIZE(Faux)) ? i : i - YSIZE(Faux))
		{
			for (long int j = 0, jp = 0; j < XSIZE(Faux); j++, jp = j)
			{
				c = DIRECT_A3D_ELEM(Faux, k, i, j).real;
				d = DIRECT_A3D_ELEM(Faux, k, i, j).imag;
				ac = a * c;
				bd = b * d;
				ab_cd = (a + b) * (c + d);
				DIRECT_A3D_ELEM(Faux, k, i, j) = Complex(ac - bd, ab_cd - ac - bd);
			}
		}
	}

	// Do IFT
	transformer.inverseFourierTransform(Faux, img);
	CenterFFT(img, false);
	img.setXmippOrigin();
	Faux.clear();
	transformer.clear();

	return;
};

void rotateAndSum2DZSliceInRealSpace(
		MultidimArray<RFLOAT>& img_ori,
		MultidimArray<RFLOAT>& img_sum,
		std::vector<RFLOAT>& weight_sum,
		int idz,
		RFLOAT outer_radius_pix,
		RFLOAT rot_angle_deg)
{
	const RFLOAT pi = 3.141592653589793238462643383279502884197;
	int Xdim, Ydim, Zdim, Ndim, x_int, y_int;
	RFLOAT cos_val, sin_val, r, r_max, u, v, x0, y0, x1, y1;
	//RFLOAT dd00, dd01, dd10, dd11, dylo, dyhi;

	Xdim = XSIZE(img_ori); Ydim = YSIZE(img_ori); Zdim = ZSIZE(img_ori); Ndim = NSIZE(img_ori);

	if( (Ndim != 1) || (Zdim < 5) || (Xdim != Ydim) || (Xdim < 5) )
	{
		REPORT_ERROR("helix.cpp::rotateAndSum2DZSliceInRealSpace(): Input MultidimArray has wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");
		return;
	}
	if( (img_ori.sameShape(img_sum) == false) || (weight_sum.size() != Zdim) )
	{
		REPORT_ERROR("helix.cpp::rotateAndSum2DZSliceInRealSpace(): Input and sum MultidimArrays have different dimensions!");
		return;
	}
	if( (idz < STARTINGZ(img_sum)) || (idz > FINISHINGZ(img_sum)) )
	{
		REPORT_ERROR("helix.cpp::add2DZsliceInto3DVolumeSum(): Invalid Z slice! (z = [" + integerToString(STARTINGZ(img_sum)) + ", " + integerToString(FINISHINGZ(img_sum)) + "] while requesting to access idz = " + integerToString(idz) + ")");
		return;
	}

	// Rotation angle should be in (+0 deg, +360 deg)
	u = rot_angle_deg / 360.;
	if( (u > 1.) || (u < 0.) )
	{
		rot_angle_deg = (u - floor(u)) * 360.;
	}
	// Rotation angle = 0, 360 deg, or outer_radius_pix < 0., no need to rotate
	if( (fabs(rot_angle_deg) < (1e-6)) || (fabs(rot_angle_deg - 360.) < (1e-6)) || (outer_radius_pix < 0.) )
	{
		return;
	}

	// Rotate (Trilinear interpolation)
	cos_val = cos(rot_angle_deg * pi / 180.);
	sin_val = sin(rot_angle_deg * pi / 180.);
	r_max = ((RFLOAT)((Xdim + 1) / 2 - 2)) - 0.01;
	if(r_max > outer_radius_pix)
	{
		r_max = outer_radius_pix;
	}

	img_ori.setXmippOrigin();
	img_sum.setXmippOrigin();

    for (long int i = STARTINGY(img_sum); i <= FINISHINGY(img_sum); i++)
    {
        for (long int j = STARTINGX(img_sum); j <= FINISHINGX(img_sum); j++)
        {
    		x0 = (RFLOAT)(j);
    		y0 = (RFLOAT)(i);
    		r = sqrt(x0 * x0 + y0 * y0);
    		if(r > r_max)
    		{
    			// Needed ?
    			//A3D_ELEM(img_sum, idz, i, j) = 0.;
    			continue;
    		}

    		x1 = x0 * cos_val - y0 * sin_val;
    		y1 = x0 * sin_val + y0 * cos_val;
    		x_int = (int)(floor(x1));
    		y_int = (int)(floor(y1));
    		u = x1 - ((RFLOAT)(x_int));
    		v = y1 - ((RFLOAT)(y_int));
    		A3D_ELEM(img_sum, idz, i, j) +=
    				A3D_ELEM(img_ori, idz, y_int, x_int) * (1. - u) * (1. - v) + A3D_ELEM(img_ori, idz, y_int + 1, x_int + 1) * u * v
    				+ A3D_ELEM(img_ori, idz, y_int, x_int + 1) * u * (1. - v) + A3D_ELEM(img_ori, idz, y_int + 1, x_int) * (1. - u) * v;

    		// This version creates artifacts on the surfaces of particles.
    		//dd00 = A3D_ELEM(img_ori, idz, y_int, x_int);
    		//dd01 = A3D_ELEM(img_ori, idz, y_int + 1, x_int);
    		//dd10 = A3D_ELEM(img_ori, idz, y_int, x_int + 1);
    		//dd11 = A3D_ELEM(img_ori, idz, y_int + 1, x_int + 1);
    		//dylo = LIN_INTERP(u, dd00, dd10);
    		//dyhi = LIN_INTERP(u, dd10, dd11);
    		//A3D_ELEM(img_sum, idz, i, j) += LIN_INTERP(v, dylo, dyhi);
        }
    }
    weight_sum[idz - STARTINGZ(img_sum)] += 1.0;

	return;
};


void rotate2DZSliceInFourierSpace(
		MultidimArray<RFLOAT>& img,
		RFLOAT rot_angle_deg,
		int padding_factor)
{
	const RFLOAT pi = 3.141592653589793238462643383279502884197;
	int Xdim, Ydim, Zdim, Ndim, x_int, y_int;
	RFLOAT sinc_val, sin_val, cos_val, x0, y0, x1, y1, r, r_max, u, v;
	bool is_neg_x;
	MultidimArray<RFLOAT> img_pad;
	MultidimArray<Complex> Faux, Faux_trans;
	FourierTransformer transformer;

	Xdim = XSIZE(img); Ydim = YSIZE(img); Zdim = ZSIZE(img); Ndim = NSIZE(img);

	if( (Zdim != 1) || (Ndim != 1) || (Xdim != Ydim) || (Xdim < 5) )
	{
		REPORT_ERROR("helix.cpp::rotate2DZSliceInFT(): Input 2D square MultidimArray has wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");
		return;
	}

	if( (padding_factor <= 0) || (padding_factor > 4) )
	{
		REPORT_ERROR("helix.cpp::rotate2DZSliceInFT(): Padding factor should be 1, 2, 3 or 4!");
		return;
	}

	// Rotation angle should be in (+0 deg, +360 deg)
	u = rot_angle_deg / 360.;
	if( (u > 1.) || (u < 0.) )
	{
		rot_angle_deg = (u - floor(u)) * 360.;
	}
	// Rotation angle = 0 or 360 deg, no need to rotate
	if( (fabs(rot_angle_deg) < (1e-6)) || (fabs(rot_angle_deg - 360.) < (1e-6)) )
	{
		return;
	}

	img_pad.clear();
	Faux.clear();
	Faux_trans.clear();
	transformer.clear();

	// Pad image & gridding
	img.setXmippOrigin();
	img_pad.initZeros(1, 1, Ydim * padding_factor, Xdim * padding_factor);
	img_pad.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY2D(img)
	{
		A2D_ELEM(img_pad, i, j) = A2D_ELEM(img, i, j);

		// Gridding
		r = sqrt((RFLOAT)(i * i + j * j)) / ((RFLOAT)(Xdim * padding_factor));
		if(fabs(r) > (1e-6))
		{
			sinc_val = sin(pi * r) / (pi * r);
			A2D_ELEM(img_pad, i, j) /= (sinc_val * sinc_val);
		}
	}

	// Do FT
	CenterFFT(img_pad, true);
	transformer.FourierTransform(img_pad, Faux);
	Faux_trans.resize(Faux);

	// Rotate (Trilinear interpolation)
	cos_val = cos(rot_angle_deg * pi / 180.);
	sin_val = sin(rot_angle_deg * pi / 180.);
	//r_max = ((RFLOAT)(Xdim * padding_factor / 2 - 1));  // -1 ?
	r_max = ((RFLOAT)((Xdim * padding_factor + 1) / 2 - 1));
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(Faux_trans)
	{
		x0 = RFLOAT(jp);
		y0 = RFLOAT(ip);
		r = sqrt(x0 * x0 + y0 * y0);
		if(r > r_max)
		{
			FFTW2D_ELEM(Faux_trans, ip, jp).real = FFTW2D_ELEM(Faux_trans, ip, jp).imag = 0.;
			continue;
		}

		x1 = x0 * cos_val - y0 * sin_val;
		y1 = x0 * sin_val + y0 * cos_val;
		is_neg_x = false;
		if(x1 < 0.)
		{
			x1 *= -1.;
			y1 *= -1.;
			is_neg_x = true;
		}
		x_int = (int)(floor(x1));
		y_int = (int)(floor(y1));
		u = x1 - ((RFLOAT)(x_int));
		v = y1 - ((RFLOAT)(y_int));
		FFTW2D_ELEM(Faux_trans, ip, jp) =
				FFTW2D_ELEM(Faux, y_int, x_int) * (1. - u) * (1. - v) + FFTW2D_ELEM(Faux, y_int + 1, x_int + 1) * u * v
				+ FFTW2D_ELEM(Faux, y_int, x_int + 1) * u * (1. - v) + FFTW2D_ELEM(Faux, y_int + 1, x_int) * (1. - u) * v;
		if(is_neg_x == true)
		{
			FFTW2D_ELEM(Faux_trans, ip, jp).imag *= -1.;
		}
	}
	Faux.clear();

	// Do IFT
	transformer.inverseFourierTransform(Faux_trans, img_pad);
	CenterFFT(img_pad, false);
	Faux_trans.clear();
	transformer.clear();

	// Rewindow
	r_max = ((RFLOAT)((Xdim + 1) / 2 - 1));
	img.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY2D(img)
	{
		r = sqrt(i * i + j * j);
		if(r > r_max)
		{
			A2D_ELEM(img, i, j) = 0.;
		}
		else
		{
			A2D_ELEM(img, i, j) = A2D_ELEM(img_pad, i, j);
		}
	}
	img_pad.clear();

	return;
};


// In Z axis [s, e] chunk, there must be a symmetrical segment of helix with at least one particle!
void expandZaxisInFourierSpace(
		MultidimArray<RFLOAT>& vol,
		RFLOAT outer_radius_pix,
		RFLOAT twist_deg,  // both + or -
		RFLOAT rise_pix,  // only +
		int idz_s,
		int idz_e,
		int padding_factor)
{
	int ii, jj, Xdim, Ydim, Zdim, Ndim, nr_particles, idz_s_new, idz_e_new;
	RFLOAT rise_mult_pix, twist_mult_deg;
	MultidimArray<RFLOAT> vol_ori, vol_trans, vol_sum, Zslice2D;
	std::vector<RFLOAT> weight_sum;

	vol_ori.clear();
	vol_trans.clear();
	vol_sum.clear();
	Zslice2D.clear();
	weight_sum.clear();
	Xdim = XSIZE(vol); Ydim = YSIZE(vol); Zdim = ZSIZE(vol); Ndim = NSIZE(vol);

	if( (Ndim != 1) || (Zdim < 5) || (Xdim != Ydim) || (Xdim < 5) )
	{
		REPORT_ERROR("helix.cpp::expandZaxis(): Input 3D MultidimArray has wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");
		return;
	}

	vol.setXmippOrigin();

	if( (rise_pix < 0.01) || (rise_pix > ((RFLOAT)(Zdim / 2)))
			|| (idz_s >= idz_e) || (idz_s < STARTINGZ(vol)) || (idz_e > FINISHINGZ(vol))
			|| ((idz_e - idz_s) < rise_pix) || (outer_radius_pix < 0.) )
	{
		REPORT_ERROR("helix.cpp::expandZaxis(): Wrong parameters!");
		return;
	}

	// No need to expand
	if( (idz_s <= STARTINGZ(vol)) && (idz_e >= FINISHINGZ(vol)) )
	{
		return;
	}

	// How many x = n*rise pixels are there in [idz_e, idz_e] ? (n is an integer)
	nr_particles = (int)(floor(((RFLOAT)(idz_e - idz_s)) / rise_pix));
	rise_mult_pix = nr_particles * rise_pix;
	twist_mult_deg = nr_particles * twist_deg;
	// DEBUG

	std::cout << "nr_particles = " << nr_particles
			<< ", rise_mult_pix = " << rise_mult_pix
			<< ", twist_mult_deg = " << twist_mult_deg << std::endl;
	std::cout << "Center z[" << idz_s << ", " << idz_e << "]" << std::endl;

	// Retain Z center part of volume
	weight_sum.resize(Zdim);
	for(ii = 0; ii < weight_sum.size(); ii++)
	{
		if( ((ii + STARTINGZ(vol)) < idz_s) || ((ii + STARTINGZ(vol)) > idz_e) )
		{
			weight_sum[ii] = 0.;
			continue;
		}
		weight_sum[ii] = 1.;
	}
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
	{
		if( (k < idz_s) || (k > idz_e) )
		{
			A3D_ELEM(vol, k, i, j) = 0.;
		}
	}

	// Init result
	vol_ori = vol;
	vol_sum = vol;
	idz_s_new = idz_s;
	idz_e_new = idz_e;

	// Shift upward (Z+)
	for(ii = 1; idz_e_new <= FINISHINGZ(vol_sum) ; ii++)
	{
		idz_s_new = (int)(ceil(idz_s + rise_mult_pix * ((RFLOAT)(ii))));
		idz_e_new = (int)(floor(idz_e + rise_mult_pix * ((RFLOAT)(ii))));
		// DEBUG
		std::cout << "Fill in z[" << idz_s_new << ", " << idz_e_new << "]" << std::endl;
		vol_trans = vol_ori;
		shift3DVolumeAlongZAxisInFourierSpace(vol_trans, rise_mult_pix * ((RFLOAT)(ii)) );
		for(jj = idz_s_new; (jj <= idz_e_new) && (jj >= STARTINGZ(vol_sum)) && (jj <= FINISHINGZ(vol_sum)); jj++)
		{
			if(weight_sum[jj - STARTINGZ(vol_sum)] > 0.9)
			{
				continue;
			}

			if( (padding_factor >= 2) && (padding_factor <= 4) )
			{
				get2DZsliceIn3DVolume(vol_trans, Zslice2D, jj);
				rotate2DZSliceInFourierSpace(Zslice2D, twist_mult_deg * ((RFLOAT)(ii)), padding_factor);
				add2DZsliceInto3DVolumeSum(Zslice2D, vol_sum, weight_sum, jj);
			}
			else
			{
				rotateAndSum2DZSliceInRealSpace(vol_trans, vol_sum, weight_sum, jj, outer_radius_pix, twist_mult_deg * ((RFLOAT)(ii)));
			}
		}
	}

	// Shift downward (Z-)
	for(ii = -1; idz_s_new >= STARTINGZ(vol_sum); ii--)
	{
		idz_s_new = (int)(ceil(idz_s + rise_mult_pix * ((RFLOAT)(ii))));
		idz_e_new = (int)(floor(idz_e + rise_mult_pix * ((RFLOAT)(ii))));
		// DEBUG
		std::cout << "Fill in z[" << idz_s_new << ", " << idz_e_new << "]" << std::endl;
		vol_trans = vol_ori;
		shift3DVolumeAlongZAxisInFourierSpace(vol_trans, rise_mult_pix * ((RFLOAT)(ii)) );
		for(jj = idz_e_new; (jj >= idz_s_new) && (jj >= STARTINGZ(vol_sum)) && (jj <= FINISHINGZ(vol_sum)); jj--)
		{
			if(weight_sum[jj - STARTINGZ(vol_sum)] > 0.9)
			{
				continue;
			}
			if( (padding_factor >= 2) && (padding_factor <= 4) )
			{
				get2DZsliceIn3DVolume(vol_trans, Zslice2D, jj);
				rotate2DZSliceInFourierSpace(Zslice2D, twist_mult_deg * ((RFLOAT)(ii)), padding_factor);
				add2DZsliceInto3DVolumeSum(Zslice2D, vol_sum, weight_sum, jj);
			}
			else
			{
				rotateAndSum2DZSliceInRealSpace(vol_trans, vol_sum, weight_sum, jj, outer_radius_pix, twist_mult_deg * ((RFLOAT)(ii)));
			}
		}
	}

	vol = vol_sum;

	// Destructions
	vol_ori.clear();
	vol_trans.clear();
	vol_sum.clear();
	Zslice2D.clear();
	weight_sum.clear();

	return;
};

void imposeHelicalSymmetryInFourierSpace(MultidimArray<RFLOAT>& vol,
		RFLOAT outer_radius_pix,
		RFLOAT twist_deg,  // both + or -
		RFLOAT rise_pix,  // only +
		int idz_s,
		int idz_e,
		int padding_factor)
{
	int ii, jj, Xdim, Ydim, Zdim, Ndim, idz_s_new, idz_e_new;
	RFLOAT twist_total_deg, rise_total_pix;
	MultidimArray<RFLOAT> vol_ori, vol_trans, vol_sum, Zslice2D;
	std::vector<RFLOAT> weight_sum;

	vol_ori.clear();
	vol_trans.clear();
	vol_sum.clear();
	Zslice2D.clear();
	weight_sum.clear();
	Xdim = XSIZE(vol); Ydim = YSIZE(vol); Zdim = ZSIZE(vol); Ndim = NSIZE(vol);

	if( (Ndim != 1) || (Zdim < 5) || (Xdim != Ydim) || (Xdim < 5) )
	{
		REPORT_ERROR("helix.cpp::imposeHelicalSymmetry(): Input 3D MultidimArray has wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");
		return;
	}

	vol.setXmippOrigin();

	if( (idz_s < STARTINGZ(vol)) || (idz_e > FINISHINGZ(vol))
			|| (idz_s >= idz_e) || (rise_pix < 0.01) || ((idz_e - idz_s) < rise_pix)
			|| (outer_radius_pix < 0.) )
	{
		REPORT_ERROR("helix.cpp::imposeHelicalSymmetry(): Wrong parameters!");
		return;
	}

	// Copy original volume and reset sums of weight
	vol_ori = vol;
	vol_sum = vol;
	weight_sum.resize(Zdim);
	for(ii = 0; ii < weight_sum.size(); ii++)
	{
		weight_sum[ii] = 1.0;
	}

	// Shift upward (Z+)
	for(ii = 1; ; ii++)
	{
		rise_total_pix = ((RFLOAT)(ii)) * rise_pix;
		twist_total_deg = ((RFLOAT)(ii)) * twist_deg;
		idz_s_new = ((int)(ceil((idz_s + rise_total_pix))));
		if(idz_s_new > idz_e)
		{
			break;
		}

		// DEBUG
		std::cout << "Shift upward " << ii << "..." << std::endl;

		vol_trans = vol_ori;
		shift3DVolumeAlongZAxisInFourierSpace(vol_trans, rise_total_pix);

		for(jj = idz_s_new; jj <= idz_e; jj++)
		{
			if( (padding_factor >= 2) && (padding_factor <= 4) )
			{
				get2DZsliceIn3DVolume(vol_trans, Zslice2D, jj);
				rotate2DZSliceInFourierSpace(Zslice2D, twist_total_deg, padding_factor);
				add2DZsliceInto3DVolumeSum(Zslice2D, vol_sum, weight_sum, jj);
			}
			else
			{
				rotateAndSum2DZSliceInRealSpace(vol_trans, vol_sum, weight_sum, jj, outer_radius_pix, twist_total_deg);
			}
		}
	}

	// Shift downward (Z-)
	for(ii = -1; ; ii--)
	{
		rise_total_pix = ((RFLOAT)(ii)) * rise_pix;
		twist_total_deg = ((RFLOAT)(ii)) * twist_deg;
		idz_e_new = ((int)(floor((idz_e + rise_total_pix))));
		if(idz_e_new < idz_s)
		{
			break;
		}

		// DEBUG
		std::cout << "Shift downward " << ii << "..." << std::endl;

		vol_trans = vol_ori;
		shift3DVolumeAlongZAxisInFourierSpace(vol_trans, rise_total_pix);

		for(jj = idz_e_new; jj >= idz_s; jj--)
		{
			if( (padding_factor >= 2) && (padding_factor <= 4) )
			{
				get2DZsliceIn3DVolume(vol_trans, Zslice2D, jj);
				rotate2DZSliceInFourierSpace(Zslice2D, twist_total_deg, padding_factor);
				add2DZsliceInto3DVolumeSum(Zslice2D, vol_sum, weight_sum, jj);
			}
			else
			{
				rotateAndSum2DZSliceInRealSpace(vol_trans, vol_sum, weight_sum, jj, outer_radius_pix, twist_total_deg);
			}
		}
	}

	// Sum / Weight - Only retain central part of Z - [idz_s, idz_e]
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_sum)
	{
		if( (k < idz_s) || (k > idz_e) )
		{
			A3D_ELEM(vol_sum, k, i, j) = 0.;
			continue;
		}
		A3D_ELEM(vol_sum, k, i, j) /= weight_sum[k - STARTINGZ(vol_sum)];
	}
	vol = vol_sum;

	// Destructions
	vol_ori.clear();
	vol_trans.clear();
	vol_sum.clear();
	Zslice2D.clear();
	weight_sum.clear();

	return;
};
*/

RFLOAT getHelicalSigma2Rot(
		RFLOAT helical_rise_pix,
		RFLOAT helical_twist_deg,
		RFLOAT helical_offset_step_pix,
		RFLOAT rot_step_deg,
		RFLOAT old_sigma2_rot)
{
	if ( (helical_offset_step_pix < 0.) || (rot_step_deg < 0.) || (old_sigma2_rot < 0.) )
		REPORT_ERROR("helix.cpp::getHelicalSigma2Rot: Helical offset step, rot step or sigma2_rot cannot be negative!");

	RFLOAT nr_samplings_along_helical_axis = (fabs(helical_rise_pix)) / helical_offset_step_pix;
	RFLOAT rot_search_range = (fabs(helical_twist_deg)) / nr_samplings_along_helical_axis;
	RFLOAT new_rot_step = rot_search_range / 6.;
	RFLOAT factor = CEIL(new_rot_step / rot_step_deg);
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
		int idz_e)
{
	const RFLOAT pi = 3.141592653589793238462643383279502884197;
	int ii, Zdim, x0, y0, z0, x1, y1, z1, tab_len, nn;
	RFLOAT pix_sum, pix_weight;
	RFLOAT x_ori, y_ori, z_ori, xp, yp, zp, fx, fy, fz, r, d, r_min, r_max, d_min, d_max, D_min, D_max, rot_angle_deg, cos_val, sin_val;
	RFLOAT d000, d001, d010, d011, d100, d101, d110, d111, dx00, dx01, dx10, dx11, dxy0, dxy1;
	std::vector<RFLOAT> tab_sin, tab_cos;
	//MultidimArray<RFLOAT> vol_in, vol_out;
	MultidimArray<RFLOAT> vol_out;

	tab_sin.clear();
	tab_cos.clear();
	//vol_in.clear();
	vol_out.clear();
	Zdim = ZSIZE(vol);
	//Xdim = XSIZE(vol); Ydim = YSIZE(vol); Zdim = ZSIZE(vol); Ndim = NSIZE(vol);

	/*
	if( (Ndim != 1) || (Zdim < 5) || (Xdim != Ydim) || (Xdim < 5) )
	{
		REPORT_ERROR("helix.cpp::imposeHelicalSymmetry(): Input 3D MultidimArray has wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");
		return;
	}
	*/

	vol.setXmippOrigin();

	/*
	if( (idz_s < STARTINGZ(vol)) || (idz_e > FINISHINGZ(vol))
			|| (idz_s >= idz_e) || (rise_pix < 0.01) || ((idz_e - idz_s) < (2. * rise_pix))
			|| (cyl_outer_radius_pix < 0.) )
	{
		REPORT_ERROR("helix.cpp::imposeHelicalSymmetry(): Wrong parameters!");
		return;
	}
	*/

	// Set max radius
	/*
	r_max = ((RFLOAT)((Xdim + 1) / 2 - 2)) - 0.01;
	if(r_max > cyl_outer_radius_pix)
	{
		r_max = cyl_outer_radius_pix;
	}
	*/

	// Set tab cos, sin values
	rise_pix = fabs(rise_pix); // Rise should be a positive value!
	tab_len = 10 + ROUND((Zdim + 2) / rise_pix);
	tab_sin.resize(tab_len);
	tab_cos.resize(tab_len);
	for(ii = 0; ii < tab_len; ii++)
	{
		rot_angle_deg = ((RFLOAT)(ii)) * twist_deg;
		tab_cos[ii] = cos(rot_angle_deg * pi / 180.);
		tab_sin[ii] = sin(rot_angle_deg * pi / 180.);
	}

	// Jun14,2015 - Shaoda, apply helical mask
	r_min = sphere_radius_pix;
	r_max = sphere_radius_pix + cosine_width;
	d_min = cyl_inner_radius_pix - cosine_width;
	d_max = cyl_inner_radius_pix;
	D_min = cyl_outer_radius_pix;
	D_max = cyl_outer_radius_pix + cosine_width;
	// DEBUG
	//std::cout << "r_min, r_max, d_min, d_max, D_min, D_max = "
	//		<< r_min << ", " << r_max << ", " << d_min << ", " << d_max << ", " << D_min << ", " << D_max << std::endl;

	// Copy original volume
	//vol.setXmippOrigin();
	//vol_in = vol;
	vol_out.resize(vol);
	vol_out.setXmippOrigin();

	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_out)
	{

		/*
		if( (i == STARTINGY(vol_out)) && (j == STARTINGX(vol_out)) )
		{
			std::cout << "z = " << k << std::endl;
		}
		*/


		d = ((RFLOAT)(i * i + j * j));
		r = d + ((RFLOAT)(k * k));
		d = sqrt(d);
		r = sqrt(r);
		if ( (d < d_min) || (d > D_max) || (r > r_max) )
		{
			A3D_ELEM(vol_out, k, i, j) = 0.;
			continue;
		}

		x_ori = ((RFLOAT)(j));
		y_ori = ((RFLOAT)(i));
		z_ori = ((RFLOAT)(k));
		pix_sum = pix_weight = 0.;

		if (z_ori > idz_s)
		{
			nn = FLOOR((z_ori - idz_s) / rise_pix);
			nn *= -1;
		}
		else
		{
			nn = CEIL((idz_s - z_ori) / rise_pix);
		}
		for(ii = nn; ; ii++)
		{
			zp = z_ori + ((RFLOAT)(ii)) * rise_pix;
			if(zp < idz_s)
			{
				continue;
			}
			if(zp > idz_e)
			{
				break;
			}

			if(ii >= 0)
			{
				cos_val = tab_cos[ii];
				sin_val = tab_sin[ii];
			}
			else
			{
				cos_val = tab_cos[-ii];
				sin_val = (-1.) * tab_sin[-ii];
			}

    		xp = x_ori * cos_val - y_ori * sin_val;
    		yp = x_ori * sin_val + y_ori * cos_val;

			// Trilinear interpolation (with physical coords)
			// Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
			// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
			x0 = FLOOR(xp);
			fx = xp - x0;
			x0 -= STARTINGX(vol_out);
			x1 = x0 + 1;

			y0 = FLOOR(yp);
			fy = yp - y0;
			y0 -= STARTINGY(vol_out);
			y1 = y0 + 1;

			z0 = FLOOR(zp);
			fz = zp - z0;
			z0 -= STARTINGZ(vol_out);
			z1 = z0 + 1;

			if( (z0 < 0) || (z1 >= ZSIZE(vol_out)) )
			{
				continue;
			}

			/*
			d000 = DIRECT_A3D_ELEM(vol_in, z0, y0, x0);
			d001 = DIRECT_A3D_ELEM(vol_in, z0, y0, x1);
			d010 = DIRECT_A3D_ELEM(vol_in, z0, y1, x0);
			d011 = DIRECT_A3D_ELEM(vol_in, z0, y1, x1);
			d100 = DIRECT_A3D_ELEM(vol_in, z1, y0, x0);
			d101 = DIRECT_A3D_ELEM(vol_in, z1, y0, x1);
			d110 = DIRECT_A3D_ELEM(vol_in, z1, y1, x0);
			d111 = DIRECT_A3D_ELEM(vol_in, z1, y1, x1);
			*/

			d000 = DIRECT_A3D_ELEM(vol, z0, y0, x0);
			d001 = DIRECT_A3D_ELEM(vol, z0, y0, x1);
			d010 = DIRECT_A3D_ELEM(vol, z0, y1, x0);
			d011 = DIRECT_A3D_ELEM(vol, z0, y1, x1);
			d100 = DIRECT_A3D_ELEM(vol, z1, y0, x0);
			d101 = DIRECT_A3D_ELEM(vol, z1, y0, x1);
			d110 = DIRECT_A3D_ELEM(vol, z1, y1, x0);
			d111 = DIRECT_A3D_ELEM(vol, z1, y1, x1);

			dx00 = LIN_INTERP(fx, d000, d001);
			dx01 = LIN_INTERP(fx, d100, d101);
			dx10 = LIN_INTERP(fx, d010, d011);
			dx11 = LIN_INTERP(fx, d110, d111);

			dxy0 = LIN_INTERP(fy, dx00, dx10);
			dxy1 = LIN_INTERP(fy, dx01, dx11);

			pix_sum += LIN_INTERP(fz, dxy0, dxy1);
			pix_weight += 1.;
		}
		if(pix_weight > 0.9)
		{
			A3D_ELEM(vol_out, k, i, j) = pix_sum / pix_weight;

			// Jun14,2015 - Shaoda, apply soft helical mask
			if ( (d > d_max) && (d < D_min) && (r < r_min) )
			{}
			else // The pixel is within cosine edge(s)
			{
				pix_weight = 1.;
				if (d < d_max)  // d_min < d < d_max : w=(0~1)
					pix_weight = 0.5 + (0.5 * cos(PI * ((d_max - d) / cosine_width)));
				else if (d > D_min) // D_min < d < D_max : w=(1~0)
					pix_weight = 0.5 + (0.5 * cos(PI * ((d - D_min) / cosine_width)));
				if (r > r_min) // r_min < r < r_max
				{
					pix_sum = 0.5 + (0.5 * cos(PI * ((r - r_min) / cosine_width)));
					pix_weight = (pix_sum < pix_weight) ? (pix_sum) : (pix_weight);
				}
				A3D_ELEM(vol_out, k, i, j) *= pix_weight;
				pix_weight = 0.;
			}
		}
		else
		{
			A3D_ELEM(vol_out, k, i, j) = 0.;
		}
	}

	/*
	// Jun14,2015 - Shaoda, apply soft helical mask
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_out)
	{
		d = ((RFLOAT)(i * i + j * j));
		r = d + ((RFLOAT)(k * k));
		d = sqrt(d);
		r = sqrt(r);

		if ( (d < d_min) || (d > D_max) || (r > r_max) )
		{
			A3D_ELEM(vol_out, k, i, j) = 0.;
			continue;
		}
		else if ( (d > d_max) && (d < D_min) && (r < r_min) )
		{

		}

		if ( (d > d_max) && (d < D_min) && (r < r_min) )
		{}
		else // The pixel is within cosine edge(s)
		{
			pix_weight = 1.;
			if (d < d_max)  // d_min < d < d_max : w=(0~1)
				pix_weight = 0.5 + (0.5 * cos(PI * ((d_max - d) / cosine_width)));
			else if (d > D_min) // D_min < d < D_max : w=(1~0)
				pix_weight = 0.5 + (0.5 * cos(PI * ((d - D_min) / cosine_width)));
			if (r > r_min) // r_min < r < r_max
			{
				pix_sum = 0.5 + (0.5 * cos(PI * ((r - r_min) / cosine_width)));
				pix_weight = (pix_sum < pix_weight) ? (pix_sum) : (pix_weight);
			}
			A3D_ELEM(vol_out, k, i, j) *= pix_weight;
			pix_weight = 0.;
		}
	}
	*/

	vol = vol_out;
	tab_sin.clear();
	tab_cos.clear();
	//vol_in.clear();
	vol_out.clear();

	return;
};

// Assume all parameters are within range
RFLOAT get_lenZ_percentage_max(
		int box_len,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_outer_radius_A,
		RFLOAT pixel_size_A)
{
	return (((2.) * sqrt(sphere_radius_A * sphere_radius_A - cyl_outer_radius_A * cyl_outer_radius_A) / pixel_size_A) / box_len);
};

// Assume all parameters are within range
RFLOAT get_rise_A_max(
		int box_len,
		RFLOAT pixel_size_A,
		RFLOAT lenZ_percentage,
		RFLOAT nr_units_min)
{
	return (pixel_size_A * box_len * lenZ_percentage / nr_units_min);
};

bool checkHelicalParametersFor3DHelicalReference(
		int box_len,
		RFLOAT pixel_size_A,
		RFLOAT twist_deg,
		RFLOAT rise_A,
		RFLOAT lenZ_percentage,
		bool do_helical_symmetry_local_refinement,
		RFLOAT rise_search_max_dev_percentage,
		RFLOAT twist_search_max_dev_percentage,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_inner_radius_A,
		RFLOAT cyl_outer_radius_A)
{
	long int half_box_len;
	RFLOAT nr_units_min = 2.; // Minimum nr_particles required along lenZ_max
	RFLOAT lenZ_percentage_max, rise_A_max, sphere_radius_pix, cyl_inner_radius_pix, cyl_outer_radius_pix;

	if (pixel_size_A < 0.001)
	{
		REPORT_ERROR("helix.cpp::checkHelicalParametersFor3DHelicalReference(): Pixel size (in Angstroms) should be larger than 0.001!");
		return false;
	}

	sphere_radius_pix = sphere_radius_A / pixel_size_A;
	cyl_inner_radius_pix = cyl_inner_radius_A / pixel_size_A;
	cyl_outer_radius_pix = cyl_outer_radius_A / pixel_size_A;

	if (box_len < 10)
	{
		REPORT_ERROR("helix.cpp::checkHelicalParametersFor3DHelicalReference(): Input box size should be larger than 5!");
		return false;
	}
	half_box_len = box_len / 2 - ((box_len + 1) % 2);

	if ( (fabs(twist_deg) < 0.01) || (fabs(twist_deg) > 179.99)
			|| ((rise_A / pixel_size_A) < 0.001)
			|| (lenZ_percentage < 0.001) || (lenZ_percentage > 0.999) )
	{
		REPORT_ERROR("helix.cpp::checkHelicalParametersFor3DHelicalReference(): Wrong helical twist, rise or lenZ!");
		return false;
	}
	if ( (sphere_radius_pix < 2.) || (sphere_radius_pix > half_box_len)
			|| ( (cyl_inner_radius_pix + 2.) > cyl_outer_radius_pix)
			|| (cyl_outer_radius_pix < 2.) || (cyl_outer_radius_pix > half_box_len)
			|| ( (sphere_radius_pix + 0.001) < cyl_outer_radius_pix ) )
	{
		std::cout << "sphere_radius_pix= " << sphere_radius_pix << ", half_box_len= " << half_box_len
				<< ", cyl_inner_radius_pix= " << cyl_inner_radius_pix << ", cyl_outer_radius_pix= " << cyl_outer_radius_pix << std::endl;
		REPORT_ERROR("helix.cpp::checkHelicalParametersFor3DHelicalReference(): Radii of spherical and/or cylindrical masks are invalid!");
		return false;
	}

	lenZ_percentage_max = get_lenZ_percentage_max(box_len, sphere_radius_A, cyl_outer_radius_A, pixel_size_A);
	if (lenZ_percentage > lenZ_percentage_max)
	{
		REPORT_ERROR("helix.cpp::checkHelicalParametersFor3DHelicalReference(): Central Z part is too long. (lenZ_percentage = "
				+ floatToString(lenZ_percentage) + "; lenZ_percentage < " + floatToString(lenZ_percentage_max) + ")");
		return false;
	}
	rise_A_max = get_rise_A_max(box_len, pixel_size_A, lenZ_percentage, nr_units_min);
	if (fabs(rise_A) > rise_A_max)
	{
		REPORT_ERROR("helix.cpp::checkHelicalParametersFor3DHelicalReference(): Central Z part is too short (< nr_particles_min * helical_rise_A). (rise_A = "
				+ floatToString(rise_A) + ", lenZ_percentage = " + floatToString(lenZ_percentage)
				+ ", nr_particles_min = " + floatToString(nr_units_min) + "; lenZ_percentage > "
				+ floatToString((nr_units_min * rise_A_max) / (pixel_size_A * box_len)) + ")");
		return false;
	}

	if (do_helical_symmetry_local_refinement)
	{
		if ( (rise_search_max_dev_percentage < 0.0099) || (rise_search_max_dev_percentage > 0.3301)
				|| (twist_search_max_dev_percentage < 0.0099) || (twist_search_max_dev_percentage > 0.3301) )
		{
			REPORT_ERROR("helix.cpp::checkHelicalParametersFor3DHelicalReference(): Maximum deviation of local searches should between 1% and 33%!");
			return false;
		}
	}

	return true;
};

void makeHelicalReferenceInRealSpace(
		MultidimArray<RFLOAT>& vol,
		RFLOAT pixel_size_A,
		RFLOAT twist_deg,
		RFLOAT rise_A,
		RFLOAT lenZ_percentage,
		RFLOAT sphere_radius_A,
		RFLOAT cyl_inner_radius_A,
		RFLOAT cyl_outer_radius_A,
		RFLOAT cosine_width_pix)
{
	long int Xdim, Ydim, Zdim, Ndim, box_len;
	int idz_s, idz_e;
	RFLOAT rise_pix, sphere_radius_pix, cyl_inner_radius_pix, cyl_outer_radius_pix;
	RFLOAT rr, dd, RR, DD;

	if (vol.getDim() != 3)
	{
		REPORT_ERROR("helix.cpp::makeHelicalReferenceInRealSpace(): Input helical reference is not 3D! (vol.getDim() = "
						+ integerToString(vol.getDim()) + ")");
		return;
	}
	vol.getDimensions(Xdim, Ydim, Zdim, Ndim);
	box_len = (Xdim < Ydim) ? Xdim : Ydim;
	box_len = (box_len < Zdim) ? box_len : Zdim;
	checkHelicalParametersFor3DHelicalReference(
			box_len,
			pixel_size_A,
			twist_deg,
			rise_A,
			lenZ_percentage,
			false,
			0.,
			0.,
			sphere_radius_A,
			cyl_inner_radius_A,
			cyl_outer_radius_A);

	vol.setXmippOrigin();
	idz_s = (int)(Zdim * lenZ_percentage / (-2.));
	idz_e = (int)(Zdim * lenZ_percentage / 2.);
	idz_s = (idz_s < STARTINGZ(vol)) ? (STARTINGZ(vol)) : (idz_s);
	idz_e = (idz_e > FINISHINGZ(vol)) ? (FINISHINGZ(vol)) : (idz_e);

	// DEBUG
	//std::cout << "vol, twist_deg = " << twist_deg << ", rise_pix = " << rise_pix
	//		<< ", idz_s = " << idz_s << ", idz_e = " << idz_e << std::endl;

	// Real space version
	//std::cout << "imposeHelicalSymmetry ..." << std::endl;
	rise_pix = rise_A / pixel_size_A;
	sphere_radius_pix = sphere_radius_A / pixel_size_A;
	cyl_inner_radius_pix = cyl_inner_radius_A / pixel_size_A;
	cyl_outer_radius_pix = cyl_outer_radius_A / pixel_size_A;
	imposeHelicalSymmetryInRealSpace(
			vol,
			sphere_radius_pix,
			cyl_inner_radius_pix,
			cyl_outer_radius_pix,
			cosine_width_pix,
			twist_deg,
			rise_pix,
			idz_s,
			idz_e);

	// Fourier space version
	//std::cout << "imposeHelicalSymmetry ..." << std::endl;
	//imposeHelicalSymmetryInFourierSpace(vol, helical_outer_radius_pix, helical_twist_deg, helical_rise_pix, idz_s, idz_e, padding_factor);
	//std::cout << "expandZaxis ..." << std::endl;
	//expandZaxisInFourierSpace(vol, helical_outer_radius_pix, helical_twist_deg, helical_rise_pix, idz_s, idz_e, padding_factor);

	// Spherical and cylindrical mask
	// Info around the soft edges remains intact. It will be soft-masked afterwards.
	/*
	RR = (sphere_radius_pix + cosine_width + 0.1) * (sphere_radius_pix + cosine_width + 0.1);
	DD = (cyl_outer_radius_pix + cosine_width + 0.1) * (cyl_outer_radius_pix + cosine_width + 0.1);
	vol.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
	{
		dd = ((RFLOAT)(i * i + j * j));
		rr = dd + ((RFLOAT)(k * k));
		if ( (rr > RR) || (dd > DD) )
			A3D_ELEM(vol, k, i, j) = 0.;
	}
	*/

	return;
};

/*
bool calcCCOfCnZSymmetry(
		const MultidimArray<RFLOAT>& v,
		RFLOAT r_min_pix,
		RFLOAT r_max_pix,
		int cn,
		RFLOAT& cc,
		int& nr_asym_voxels)
{
	const RFLOAT pi = 3.141592653589793238462643383279502884197;
	int id, Xdim, Ydim, Zdim, Ndim, r_max_XY;
	int x0, y0, x1, y1, z;
	RFLOAT cn_rot_rad, dist_r_pix, atan2_rad, rot_rad_total, sum_pw1, sum_pw2;
	RFLOAT xp, yp, fx, fy, d00, d01, d10, d11, dx0, dx1, dd;
	bool ok_flag;
	MultidimArray<RFLOAT> vol;
	std::vector<RFLOAT> dev_voxel, dev_chunk;

	Xdim = XSIZE(v); Ydim = YSIZE(v); Zdim = ZSIZE(v); Ndim = NSIZE(v);

	if( (Ndim != 1) || (Zdim < 5) || (Ydim < 5) || (Xdim < 5) )
	{
		REPORT_ERROR("helix.cpp::calcCCOfCnZSymmetry(): Input 3D MultidimArray has wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");
		return false;
	}

	if( (r_max_pix < 2.) || (r_min_pix < 0.) || ((r_max_pix - r_min_pix) < 2.) || (cn <= 1) || (cn > 36) )
	{
		REPORT_ERROR("helix.cpp::calcCCOfCnZSymmetry(): Wrong parameters!");
		return false;
	}

	// Check r_max
	r_max_XY = (Xdim < Ydim) ? Xdim : Ydim;
	r_max_XY = (r_max_XY + 1) / 2 - 1;
	if( r_max_pix > (((RFLOAT)(r_max_XY)) - 0.01) )  // 0.01 - avoid segmentation fault
	{
		r_max_pix = (((RFLOAT)(r_max_XY)) - 0.01);
	}

	vol.clear();
	vol = v;
	vol.setXmippOrigin();

	dev_chunk.clear();
	cn_rot_rad = (360. / ((RFLOAT)(cn))) * pi / 180.;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
	{
		dist_r_pix = sqrt(i * i + j * j);
		if( (dist_r_pix < r_min_pix) || (dist_r_pix > r_max_pix) )
		{
			continue;
		}
		atan2_rad = pi + atan2(j, i); // deg [0, 360] = rad [0, 2 * pi]
		if(atan2_rad > cn_rot_rad)
		{
			continue;
		}

		// Pick a voxel in the chunk
		dev_voxel.clear();
		dev_voxel.push_back(A3D_ELEM(vol, k, i, j));

		// Pick other voxels according to this voxel and Cn rotational symmetry
		for(id = 1; id < cn; id++)
		{
			rot_rad_total = cn_rot_rad * ((RFLOAT)(id));
			xp = ((RFLOAT)(j)) * cos(rot_rad_total) - ((RFLOAT)(i)) * sin(rot_rad_total);
			yp = ((RFLOAT)(j)) * sin(rot_rad_total) + ((RFLOAT)(i)) * cos(rot_rad_total);

			// Bilinear interpolation (with physical coords)
			// Subtract STARTINGX, STARTINGY and STARTINGZ to accelerate access to data
			// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
			x0 = FLOOR(xp); fx = xp - x0; x0 -= STARTINGX(vol); x1 = x0 + 1;
			y0 = FLOOR(yp); fy = yp - y0; y0 -= STARTINGY(vol); y1 = y0 + 1;
			z = k; z -= STARTINGZ(vol);

			// Get values
			d00 = DIRECT_A3D_ELEM(vol, z, y0, x0);
			d01 = DIRECT_A3D_ELEM(vol, z, y0, x1);
			d10 = DIRECT_A3D_ELEM(vol, z, y1, x0);
			d11 = DIRECT_A3D_ELEM(vol, z, y1, x1);

			// Interpolation 2D -> 1D
			dx0 = LIN_INTERP(fx, d00, d01);
			dx1 = LIN_INTERP(fx, d10, d11);

			// Interpolation of 2 voxels
			dd = LIN_INTERP(fy, dx0, dx1);

			// Record this voxel
			dev_voxel.push_back(dd);
		}

		// Calc dev of this voxel in the chunk
		if(dev_voxel.size() > 1)
		{
			sum_pw1 = sum_pw2 = 0.;
			for(id = 0; id < dev_voxel.size(); id++)
			{
				sum_pw1 += dev_voxel[id];
				sum_pw2 += dev_voxel[id] * dev_voxel[id];
			}
			sum_pw1 /= dev_voxel.size();
			sum_pw2 /= dev_voxel.size();
			dev_chunk.push_back(sum_pw2 - sum_pw1 * sum_pw1);
		}
		dev_voxel.clear();
	}

	// Calc avg of all voxels' devs in this chunk (for a specific helical symmetry)
	ok_flag = true;
	if(dev_chunk.size() < 1)
	{
		cc = (-1.);
		ok_flag = false;
	}
	else
	{
		sum_pw1 = 0.;
		for(id = 0; id < dev_chunk.size(); id++)
		{
			sum_pw1 += dev_chunk[id];
		}
		cc = (sum_pw1 / dev_chunk.size());
	}
	nr_asym_voxels = dev_chunk.size();

	dev_chunk.clear();
	vol.clear();

	return ok_flag;
};

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
	{
		REPORT_ERROR("helix.cpp::searchCnZSymmetry(): Input 3D MultidimArray has Wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");
		return;
	}

	if( (r_max_pix < 2.) || (r_min_pix < 0.) || ((r_max_pix - r_min_pix) < 2.)
			|| (cn_start <= 1) || (cn_end <= 1) || (cn_start > cn_end) || (cn_end > 36) )
	{
		REPORT_ERROR("helix.cpp::searchCnZSymmetry(): Wrong parameters!");
		return;
	}

	cn_list.clear();
	cc_list.clear();
	nr_asym_voxels_list.clear();

	for(cn = cn_start; cn <= cn_end; cn++)
	{
		ok_flag = calcCCOfCnZSymmetry(v, r_min_pix, r_max_pix, cn, cc, nr_asym_voxels);
		if(ok_flag == false)
		{
			continue;
		}

		cn_list.push_back(cn);
		cc_list.push_back(cc);
		nr_asym_voxels_list.push_back(nr_asym_voxels);

		if(fout_ptr != NULL)
		{
			(*fout_ptr) << "Test Cn = " << cn << ", cc = " << cc
					<< ",                               asym voxels = " << nr_asym_voxels << std::endl;
		}
	}

	return;
}
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
		RFLOAT& twist_refined_deg)
{
	int ii, nr_tests, nr_asym_voxels, id_best;
	RFLOAT nr_units_min = 2.; // Minimum nr_particles required along lenZ_max
	RFLOAT twist_test_deg, cc, min_cc, twist_best_deg, rise_max_pix;
	bool isValid;

	rise_max_pix = get_rise_A_max(ZSIZE(v), 1., z_percentage, nr_units_min);

	if ( (rise_pix < 0.001) || (rise_pix > rise_max_pix)
			|| (fabs(twist_deg) < 0.01) || (fabs(twist_deg) > 179.99)
			|| (search_half_range_deg < (1e-12)) || (search_step_deg < (1e-12)) )
	{
		REPORT_ERROR("helix.cpp::localSearchHelicalTwist(): Errors found in the helical twist and/or search parameters!");
		return false;
	}

	nr_tests = search_half_range_deg / search_step_deg;
	if (nr_tests > 5000)
		std::cout << " WARNING: More than 10000 tests are running for helical twist refinement. It is extremely time consuming..." << std::endl;

	min_cc = cc = (1e+10);
	twist_best_deg = twist_deg;
	id_best = 0;
	isValid = false;
	for(ii = -nr_tests; ii <= nr_tests; ii++)
	{
		twist_test_deg = twist_deg + ii * search_step_deg;
		if ( (fabs(twist_test_deg) < 0.01) || (fabs(twist_test_deg) > 179.99)
				|| ((twist_test_deg * twist_deg) < 0.) )
			continue;

		isValid = calcCCofHelicalSymmetry(v, r_min_pix, r_max_pix, z_percentage, rise_pix, twist_test_deg, cc, nr_asym_voxels);
		if (isValid)
		{
			if (cc < min_cc)
			{
				min_cc = cc;
				twist_best_deg = twist_test_deg;
				id_best = ii;
			}
		}
		// DEBUG
		std::cout << "Twist = " << twist_test_deg << ", Rise = " << rise_pix << ", cc = " << cc << std::endl;
	}
	twist_refined_deg = twist_best_deg;
	// DEBUG
	std::cout << "-------------------------------------------------------" << std::endl;

	if ( (min_cc > (0.99e+10)) || (ABS(id_best) == nr_tests) )
		return false;
	return true;
}

bool localSearchHelicalRise(
		const MultidimArray<RFLOAT>& v,
		RFLOAT r_min_pix,
		RFLOAT r_max_pix,
		RFLOAT z_percentage,
		RFLOAT rise_pix,
		RFLOAT twist_deg,
		RFLOAT search_half_range_pix,
		RFLOAT search_step_pix,
		RFLOAT& rise_refined_pix)
{
	int ii, nr_tests, nr_asym_voxels, id_best;
	RFLOAT nr_units_min = 2.; // Minimum nr_particles required along lenZ_max
	RFLOAT rise_test_pix, cc, min_cc, rise_best_pix, rise_max_pix;
	bool isValid;

	rise_max_pix = get_rise_A_max(ZSIZE(v), 1., z_percentage, nr_units_min);

	if ( (rise_pix < 0.001) || (rise_pix > rise_max_pix)
			|| (fabs(twist_deg) < 0.01) || (fabs(twist_deg) > 179.99)
			|| (search_half_range_pix < (1e-12)) || (search_step_pix < (1e-12)) )
	{
		REPORT_ERROR("helix.cpp::localSearchHelicalRise(): Errors found in the helical rise and/or search parameters!");
		return false;
	}

	nr_tests = search_half_range_pix / search_step_pix;
	if (nr_tests > 5000)
		std::cout << " WARNING: More than 10000 tests are running for helical rise refinement. It is extremely time consuming..." << std::endl;

	min_cc = cc = (1e+10);
	rise_best_pix = rise_pix;
	id_best = 0;
	isValid = false;
	for(ii = -nr_tests; ii <= nr_tests; ii++)
	{
		rise_test_pix = rise_pix + ii * search_step_pix;
		if ( (rise_test_pix < 0.001) || (rise_test_pix > rise_max_pix) )
			continue;

		isValid = calcCCofHelicalSymmetry(v, r_min_pix, r_max_pix, z_percentage, rise_test_pix, twist_deg, cc, nr_asym_voxels);
		if (isValid)
		{
			if (cc < min_cc)
			{
				min_cc = cc;
				rise_best_pix = rise_test_pix;
				id_best = ii;
			}
		}
		// DEBUG
		std::cout << "Twist = " << twist_deg << ", Rise = " << rise_test_pix << ", cc = " << cc << std::endl;
	}
	rise_refined_pix = rise_best_pix;
	// DEBUG
	std::cout << "-------------------------------------------------------" << std::endl;

	if ( (min_cc > (0.99e+10)) || (ABS(id_best) == nr_tests) )
		return false;
	return true;
}

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
		RFLOAT& twist_refined_deg)
{
	int nr_iter, nr_round, box_len;
	bool isParametersValid, isTwistValid, isRiseValid;
	RFLOAT r_min_pix, r_max_pix, rise_ori_pix;
	RFLOAT twist_half_range_deg, rise_half_range_pix, twist_step_deg, rise_step_pix, twist_new_deg, rise_new_pix;

	if (v.getDim() != 3)
	{
		REPORT_ERROR("helix.cpp::localSearchHelicalSymmetry(): Input helical reference is not 3D! (v.getDim() = "
						+ integerToString(v.getDim()) + ")");
		return false;
	}
	box_len = (XSIZE(v) < YSIZE(v)) ? XSIZE(v) : YSIZE(v);
	box_len = (box_len < ZSIZE(v)) ? box_len : ZSIZE(v);
	isParametersValid = checkHelicalParametersFor3DHelicalReference(
			box_len,
			pixel_size_A,
			twist_ori_deg,
			rise_ori_A,
			z_percentage,
			true,
			local_search_max_dev_percentage,
			local_search_max_dev_percentage,
			sphere_radius_A,
			r_min_A,
			r_max_A);
	if (isParametersValid == false)
	{
		REPORT_ERROR("helix.cpp::localSearchHelicalSymmetry(): Input Parameters error!");
		return false;
	}

	r_min_pix = r_min_A / pixel_size_A;
	r_max_pix = r_max_A / pixel_size_A;
	rise_ori_pix = rise_ori_A / pixel_size_A;

	// Initialization - for iteration #1
	// Step sizes of local searches should be less than 0.1% of the original parameters
	// And also at least 10 sampling points should always be guaranteed
	twist_half_range_deg = fabs(twist_ori_deg * local_search_max_dev_percentage);
	rise_half_range_pix = fabs(rise_ori_pix * local_search_max_dev_percentage);
	twist_step_deg = fabs(twist_ori_deg * 0.001);
	rise_step_pix = fabs(rise_ori_pix * 0.001);
	if (twist_half_range_deg / twist_step_deg < 5.)
		twist_step_deg = twist_half_range_deg / 5.;
	if (rise_half_range_pix / rise_step_pix < 5.)
		rise_step_pix = rise_half_range_pix / 5.;

	isTwistValid = isRiseValid = false;
	for(nr_iter = 1; nr_iter <= 20; nr_iter++)
	{
		nr_round = 0;
		while(1)
		{
			// Check convergence for a specific symmetry sampling
			nr_round++;
			if (nr_round > 20)
			{
				std::cout << " WARNING: Refined helical symmetry cannot converge." << std::endl;
				return false;
			}

			// DEBUG
			std::cout << "      ### iter = " << nr_iter
					<< "; twist half range, step (in degrees) = "
					<< twist_half_range_deg << ", " << twist_step_deg
					<< "; twist half range, step (in pixels) = "
					<< rise_half_range_pix << ", " << rise_step_pix << std::endl;

			// Perform local searches - first on twist and then rise
			isTwistValid = localSearchHelicalTwist(v, r_min_pix, r_max_pix, z_percentage,
					rise_ori_pix, twist_ori_deg, twist_half_range_deg, twist_step_deg, twist_new_deg);
			isRiseValid = localSearchHelicalRise(v, r_min_pix, r_max_pix, z_percentage,
					rise_ori_pix, twist_new_deg, rise_half_range_pix, rise_step_pix, rise_new_pix);

			// Update result
			twist_refined_deg = twist_new_deg;
			rise_refined_A = rise_new_pix * pixel_size_A;

			// DEBUG
			std::cout << "      ### iter = " << nr_iter << " ### Twist = " << twist_new_deg
					<< ", Rise = " << rise_new_pix << " ###" << std::endl;

			// Check whether the refined symmetry could be out of range
			// TODO: Only break if it is the iteration #1 ???
			if ( ((isTwistValid == false) || (isRiseValid == false)) )
			{
				std::cout << " WARNING: Refined helical symmetry is out of the search range. Check whether the initial helical symmetry is reasonable. Or you may want to modify the search range." << std::endl;
				return false;
			}

			// Refined helical symmetry converges. Search for the symmetry at a finer sampling
			if ( (fabs(twist_new_deg - twist_ori_deg) < 1.5 * twist_step_deg)
					&& (fabs(rise_new_pix - rise_ori_pix) < 1.5 * rise_step_pix) )
			{
				// Update search ranges for the next iteration
				twist_half_range_deg = 2. * twist_step_deg;
				rise_half_range_pix = 2. * rise_step_pix;

				// Update search steps for the next iteration (only try 10 sampling points)
				twist_step_deg = twist_half_range_deg / 5.;
				rise_step_pix = rise_half_range_pix / 5.;
				break;
			}

			// Update initial symmetry
			twist_ori_deg = twist_new_deg;
			rise_ori_pix = rise_new_pix;
		}

		// Stop iterations if the step size is finer than 0.0001% of the true values
		if ( (fabs(twist_step_deg / twist_new_deg) < 0.000001)
				&& (fabs(rise_step_pix / rise_new_pix) < 0.000001) )
			break;
	}

	return true;
}
*/

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
		std::ofstream* fout_ptr)
{
	int rise_pix_counter, twist_deg_counter, Xdim, Ydim, Zdim, Ndim, nr_asym_voxels;
	RFLOAT rise_pix, twist_deg, cc;
	bool ok_flag;

	Xdim = XSIZE(v); Ydim = YSIZE(v); Zdim = ZSIZE(v); Ndim = NSIZE(v);

	if( (Ndim != 1) || (Zdim < 5) || (Ydim < 5) || (Xdim < 5) )
	{
		REPORT_ERROR("helix.cpp::searchHelicalSymmetry(): Input 3D MultidimArray has wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");
		return;
	}

	if( (r_max_pix < 2.) || (r_min_pix < 0.) || ((r_max_pix - r_min_pix) < 2.)
			|| (rise_pix_start < 0.0001) || (rise_pix_incr < 0.0001) || (rise_pix_steps < 1)
			|| ( rise_pix_start > ((RFLOAT)(Zdim / 2)) )
			|| ( (rise_pix_start + ((RFLOAT)(rise_pix_steps)) * rise_pix_incr) > ((RFLOAT)(Zdim / 2)) )
			|| (twist_deg_incr < 0.0001) || (twist_deg_steps < 1) )
	{
		REPORT_ERROR("helix.cpp::searchHelicalSymmetry(): Wrong parameters!");
		return;
	}

	rise_pix_list.clear();
	twist_deg_list.clear();
	cc_list.clear();
	nr_asym_voxels_list.clear();

	for(rise_pix_counter = 0; rise_pix_counter < rise_pix_steps; rise_pix_counter++)
	{
		for(twist_deg_counter = 0; twist_deg_counter < twist_deg_steps; twist_deg_counter++)
		{
			rise_pix = rise_pix_start + ((RFLOAT)(rise_pix_counter)) * rise_pix_incr;
			twist_deg = twist_deg_start + ((RFLOAT)(twist_deg_counter)) * twist_deg_incr;

			ok_flag = calcCCofHelicalSymmetry(v, r_min_pix, r_max_pix, rise_pix, twist_deg, cc, nr_asym_voxels);

			if(ok_flag == false)
			{
				continue;
			}

			rise_pix_list.push_back(rise_pix);
			twist_deg_list.push_back(twist_deg);
			cc_list.push_back(cc);
			nr_asym_voxels_list.push_back(nr_asym_voxels);

			if(fout_ptr != NULL)
			{
				(*fout_ptr) << "Test rise = " << rise_pix << ", twist = " << twist_deg << ", cc = " << cc
						<< ",                                asym voxels = " << nr_asym_voxels << std::endl;
			}
		}
	}

	return;
};
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
	{
		REPORT_ERROR("helix.cpp::calcRadialAverage(): Input 3D MultidimArray has wrong dimensions! (Ndim = " + integerToString(Ndim) + ", Zdim = " + integerToString(Zdim) + ", Ydim = " + integerToString(Ydim) + ", Xdim = " + integerToString(Xdim) + ")");
		return;
	}

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
		{
			continue;
		}
		radial_pix_counter_list[dist] += 1.;
		radial_avg_val_list[dist] += A3D_ELEM(vol, k, i, j);
	}

	for(ii = 0; ii < list_size; ii++)
	{
		if(radial_pix_counter_list[ii] > 0.9)
		{
			radial_avg_val_list[ii] /= radial_pix_counter_list[ii];
		}
		else
		{
			radial_avg_val_list[ii] = 0.;
		}
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
		REPORT_ERROR("helix.cpp::cutZCentralPartOfSoftMask(): Z length of 3D mask os less than 5!");
	if ( (z_percentage < 0.1) || (z_percentage > 0.9) )
		REPORT_ERROR("helix.cpp::cutZCentralPartOfSoftMask(): Only 10%-90% of total Z length should be retained!");
	if (cosine_width < 0.001)
		REPORT_ERROR("helix.cpp::cutZCentralPartOfSoftMask(): Cosine width for soft edge should larger than 0!");

	idz_e = ((RFLOAT)(Zdim)) * z_percentage / 2.;
	idz_s = idz_e * (-1.);
	idz_s_w = idz_s - cosine_width;
	idz_e_w = idz_e + cosine_width;
	std::cout << "z_len, z_percentage, cosine_width = "
			<< Zdim << ", " << z_percentage << ", " << cosine_width << std::endl;
	std::cout << "idz_s_w, idz_s, idz_e, idz_e_w = "
			<< idz_s_w << ", " << idz_s << ", " << idz_e << ", " << idz_e_w << std::endl;
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
	if(box_size < 5)
		REPORT_ERROR("helix.cpp::createCylindricalReference(): Invalid box size.");

	if( (inner_diameter_pix > outer_diameter_pix)
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
		if( (r > inner_radius_pix) && (r < outer_radius_pix) )
		{
			A3D_ELEM(v, k, i, j) = 1.;
			continue;
		}
		dist = -9999.;
		if( (r > outer_radius_pix) && (r < (outer_radius_pix + cosine_width)) )
			dist = r - outer_radius_pix;
		else if( (r < inner_radius_pix) && (r > (inner_radius_pix - cosine_width)) )
			dist = inner_radius_pix - r;
		if(dist > 0.)
		{
			A3D_ELEM(v, k, i, j) = 0.5 + 0.5 * cos(PI * dist / cosine_width);
			continue;
		}
		A3D_ELEM(v, k, i, j) = 0.;
	}
	return;
}


void transformCartesianAndHelicalCoords(
		Matrix1D<RFLOAT>& in,
		Matrix1D<RFLOAT>& out,
		RFLOAT psi_deg,
		RFLOAT tilt_deg,
		bool direction)
{
	int dim;
	RFLOAT x0, y0, z0;
	Matrix1D<RFLOAT> aux;
	Matrix2D<RFLOAT> A;

	dim = in.size();
	if( (dim != 2) && (dim != 3) )
		REPORT_ERROR("helix.cpp::transformCartesianAndHelicalCoords(): Vector of input coordinates should have 2 or 3 values!");

	aux.clear();
	aux.resize(3);
	XX(aux) = XX(in);
	YY(aux) = YY(in);
	ZZ(aux) = (dim == 3) ? (ZZ(in)) : (0.);

	if(dim == 2)
	{
		tilt_deg = 0.;
	}
	if(direction == CART_TO_HELICAL_COORDS)
	{
		psi_deg *= -1.;
		tilt_deg *= -1.;
	}

	A.clear();
	A.resize(3, 3);
	Euler_angles2matrix(0., tilt_deg, psi_deg, A, false);
	aux = A * aux;

	out.clear();
	out.resize(2);
	XX(out) = XX(aux);
	YY(out) = YY(aux);
	if(dim == 3)
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
		RFLOAT psi_deg,
		RFLOAT tilt_deg,
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
	transformCartesianAndHelicalCoords(in, out, psi_deg, tilt_deg, direction);
	xout = XX(out);
	yout = YY(out);
	if (dim == 3)
		zout = ZZ(out);
	return;
}

void makeBlot(MultidimArray<RFLOAT>& v, RFLOAT y, RFLOAT x, RFLOAT r)
{
	int Xdim, Ydim, Zdim, Ndim;
	RFLOAT dist, min;
	Xdim = XSIZE(v);
	Ydim = YSIZE(v);
	Zdim = ZSIZE(v);
	Ndim = NSIZE(v);
	if( (Ndim != 1) || (Zdim != 1) || (YXSIZE(v) <= 2) )
	{
		return;
	}

	min = DIRECT_A2D_ELEM(v, 0, 0);
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(v)
	{
		min = (DIRECT_A2D_ELEM(v, i, j) < min) ? DIRECT_A2D_ELEM(v, i, j) : min;
	}

	v.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY2D(v)
	{
		dist = (i - y) * (i - y) + (j - x) * (j - x);
		dist = sqrt(dist);
		if(dist < r)
		{
			A2D_ELEM(v, i, j) = min;
		}
	}
	return;
}

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




void enlarge3DReference(
		MultidimArray<RFLOAT>& v,
		int factor)
{
	MultidimArray<RFLOAT> aux;
	int i_ori, j_ori, k_ori;
	int dim = v.getDim();
	if (dim != 3)
		REPORT_ERROR("helix.cpp::enlarge3DReference(): Input image should have a dimension of 3!");
	if (factor <= 1)
		REPORT_ERROR("helix.cpp::enlarge3DReference(): Factor should be at least 2!");
	if (factor > 4)
	{
		std::cout << " WARNING: You request for " << factor
			<< " times the size of the original volume, it might use up your memory..." << std::endl;
		std::cout << " WARNING: The enlarged volume has X, Y, Z dimensions of "
				<< XSIZE(aux) << ", " << YSIZE(aux) << ", " << ZSIZE(aux) << std::endl;
	}

	aux.clear();
	v.setXmippOrigin();
	aux.resize(factor * ZSIZE(v), factor * YSIZE(v), factor * XSIZE(v));
	aux.setXmippOrigin();

	FOR_ALL_ELEMENTS_IN_ARRAY3D(aux)
	{
		i_ori = i / factor;
		j_ori = j / factor;
		k_ori = k / factor;

		if ( (k_ori < STARTINGZ(v)) || (k_ori > FINISHINGZ(v))
				|| (i_ori < STARTINGY(v)) || (i_ori > FINISHINGY(v))
				|| (j_ori < STARTINGX(v)) || (j_ori > FINISHINGX(v)) )
			A3D_ELEM(aux, k, i, j) = 0.;
		else
			A3D_ELEM(aux, k, i, j) = A3D_ELEM(v, k_ori, i_ori, j_ori);
	}

	v.clear();
	v = aux;
	aux.clear();
	return;
}

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
	if (cosine_width > 0.05 * r_max)
		r_max -= 2. * cosine_width;
	r_max *= 0.45;
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

/*
int extractCoordsForAllHelicalSegments(
		FileName& fn_in,
		FileName& fn_out,
		int nr_asu,
		RFLOAT rise_ang,
		RFLOAT pixel_size_ang,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix)
{
	int nr_segments, nr_lines_header, tube_id;
	RFLOAT rot, tilt, psi, xoff, yoff;
	char charBuf[10000];
	RFLOAT k, b, x1, y1, x2, y2, delta_x, delta_y, now_x, now_y, step, half_box_size_pix;

	// Check parameters and open files
	if ( (nr_asu < 1) || (rise_ang < 0.001) || (pixel_size_ang < 0.01) )
		REPORT_ERROR( (std::string) "helix.cpp::extractCoordsForAllHelicalSegments(): Wrong parameters!");
	if ( (box_size_pix < 2) || (Xdim < box_size_pix) || (Ydim < box_size_pix))
		REPORT_ERROR( (std::string) "helix.cpp::transformXimdispHelicalCoordsToStarFile(): Wrong dimensions or box size!");
	std::ifstream infile(fn_in.data(), std::ios_base::in);
    if (infile.fail())
    	REPORT_ERROR( (std::string) "helix.cpp::extractCoordsForAllHelicalSegments(): File " + fn_in.c_str() + (std::string) " does not exists." );
    FileName ext = fn_in.getFileFormat();
    if (ext != "star")
    	REPORT_ERROR( (std::string) "helix.cpp::extractCoordsForAllHelicalSegments(): MetadataTable should have .star extension.");
    std::ofstream outfile(fn_out.data(), std::ios_base::out);
    if (outfile.fail())
    	REPORT_ERROR( (std::string) "helix.cpp::extractCoordsForAllHelicalSegments(): Output to file: " + fn_out.c_str() + (std::string) " failed." );

    half_box_size_pix = box_size_pix / 2.;

    // Handle input and output for headers
    nr_lines_header = 0;
    while(1)
    {
    	charBuf[0] = '\0';
    	infile.getline(charBuf, 9999);
    	if ( (strlen(charBuf) < 4) || (strstr(charBuf, "_rln") != NULL) || (strstr(charBuf, "data_") != NULL) || (strstr(charBuf, "loop_") != NULL) )
    		nr_lines_header++;
    	else
    		break;
    }
    infile.seekg(std::ios_base::beg);
    for (int ii = 0; ii < nr_lines_header; ii++)
    {
    	charBuf[0] = '\0';
    	infile.getline(charBuf, 999);
    }

    outfile << std::endl << "data_" << std::endl << std::endl << "loop_" << std::endl;
    outfile << "_rlnCoordinateX #1" << std::endl << "_rlnCoordinateY #2" << std::endl;
	outfile << "_rlnAngleRot #3" << std::endl << "_rlnAngleTilt #4" << std::endl << "_rlnAnglePsi #5" << std::endl;
	outfile << "_rlnOriginX #6" << std::endl << "_rlnOriginY #7" << std::endl;
	outfile << "_rlnHelicalTubeID #8" << std::endl;

	// Calculate all coordinates for helical segments
	nr_segments = 0;
	rot = psi = xoff = yoff = now_x = now_y = 0.;
	tilt = 90.;
    tube_id = 0;
    step = nr_asu * rise_ang / pixel_size_ang;
    tube_id = 0;
    while (infile >> x1)
    {
    	tube_id++;

    	infile >> y1 >> x2 >> y2;
    	//std::cout << "x1: " << x1 << " y1: " << y1 << " x2: " << x2 << " y2: " << y2 << std::endl;
    	if (fabs(x2 - x1) > (1e-5))
    	{
    		k = (y2 - y1) / (x2 - x1);
        	b = y1 - k * x1;
        	std::cout << "Detected Line " << tube_id << ": y = (" << k << ") x + (" << b << ")" << std::flush;

        	delta_x = step / (sqrt(k * k + 1.));
        	if (x1 > x2)
        		delta_x *= -1.;
        	delta_y = k * delta_x;
    	}
    	else
    	{
    		std::cout << "Detected Line " << tube_id << ": x = (" << x1 << ")" << std::flush;
    		delta_x = 0.;
    		delta_y = step;
    		if (y1 > y2)
    			delta_y *= -1.;
    	}
    	std::cout << "   ###   Delta x: " << delta_x << "    Delta y: " << delta_y << std::endl;

    	now_x = x1;
    	now_y = y1;
    	while (1)
    	{
    		now_x += delta_x;
    		now_y += delta_y;
    		if ( ((now_x > x1) && (now_x > x2)) || ((now_x < x1) && (now_x < x2))
    				|| ((now_y > y1) && (now_y > y2)) || ((now_y < y1) && (now_y < y2)) )
    			break;
    		else
    		{
    			// Avoid segments lying on the edges of the micrographs
    			if ( (now_x < half_box_size_pix)
    					|| (now_x > (Xdim - half_box_size_pix))
    					|| (now_y < half_box_size_pix)
    					|| (now_y > (Ydim - half_box_size_pix)) )
    			{
    				continue;
    			}

    			psi = (-180.) * atan(k) / PI;
    			outfile << std::setiosflags(std::ios::fixed) << std::setprecision(6)
						<< now_x << "	" << now_y << "	"
						<< rot << "	" << tilt << "	" << psi << "	"
						<< xoff << "	" << yoff << "	"
						<< std::resetiosflags(std::ios::fixed)
    					<< tube_id << std::endl;
    			nr_segments++;
    		}
    	}
    }

    if (nr_segments < 1)
    	std::cout << "WARNING: no segments extracted from file '"
				<< fn_in << "'!" << std::endl;
    else
    	std::cout << "Total segments / particles extracted from file '"
				<< fn_in << "' = " << nr_segments << " / " << (nr_segments * nr_asu) << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;

    infile.close();
    outfile.close();

    return nr_segments;
}
*/

int extractCoordsForAllHelicalSegments(
		FileName& fn_in,
		FileName& fn_out,
		int nr_asu,
		RFLOAT rise_ang,
		RFLOAT pixel_size_ang,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix)
{
	int nr_segments, tube_id, MDobj_id;
	RFLOAT rot, tilt, psi, xoff, yoff;
	RFLOAT k, b, x1, y1, x2, y2, delta_x, delta_y, now_x, now_y, step, half_box_size_pix;
	MetaDataTable MD_in, MD_out;
	std::vector<RFLOAT> x1_coord_list, y1_coord_list, x2_coord_list, y2_coord_list;

	// Check parameters and open files
	if ( (nr_asu < 1) || (rise_ang < 0.001) || (pixel_size_ang < 0.01) )
		REPORT_ERROR( (std::string) "helix.cpp::extractCoordsForAllHelicalSegments(): Wrong parameters!");
	if ( (box_size_pix < 2) || (Xdim < box_size_pix) || (Ydim < box_size_pix))
		REPORT_ERROR( (std::string) "helix.cpp::transformXimdispHelicalCoordsToStarFile(): Wrong dimensions or box size!");
    if ( (fn_in.getExtension() != "star") || (fn_out.getExtension() != "star") )
    	REPORT_ERROR( (std::string) "helix.cpp::extractCoordsForAllHelicalSegments(): MetadataTable should have .star extension.");

    half_box_size_pix = box_size_pix / 2.;

    // Read input STAR file
    MD_in.clear();
    MD_in.read(fn_in);
    if ( (!MD_in.containsLabel(EMDL_IMAGE_COORD_X)) || (!MD_in.containsLabel(EMDL_IMAGE_COORD_Y)) )
    	REPORT_ERROR( (std::string) "helix.cpp::extractCoordsForAllHelicalSegments(): Input STAR file does not contain X and Y coordinates!");
    if (MD_in.numberOfObjects() % 2)
    	REPORT_ERROR( (std::string) "helix.cpp::extractCoordsForAllHelicalSegments(): Input coordinates should be in pairs!");
    x1_coord_list.clear();
    y1_coord_list.clear();
    x2_coord_list.clear();
    y2_coord_list.clear();
    MDobj_id = 0;
    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_in)
    {
    	MDobj_id++;
    	MD_in.getValue(EMDL_IMAGE_COORD_X, now_x);
    	MD_in.getValue(EMDL_IMAGE_COORD_Y, now_y);
    	if (MDobj_id % 2)
    	{
    		x1_coord_list.push_back(now_x);
    		y1_coord_list.push_back(now_y);
    	}
    	else
    	{
    		x2_coord_list.push_back(now_x);
    		y2_coord_list.push_back(now_y);
    	}
    }
    if ( (x1_coord_list.size() != x2_coord_list.size())
    		|| (x2_coord_list.size() != y1_coord_list.size())
			|| (y1_coord_list.size() != y2_coord_list.size())
			|| (y2_coord_list.size() != x1_coord_list.size()) )
    	REPORT_ERROR( (std::string) "helix.cpp::extractCoordsForAllHelicalSegments(): BUG in reading input STAR file!");
    MD_in.clear();

    // Init output STAR file
    MD_out.clear();
    MD_out.addLabel(EMDL_IMAGE_COORD_X);
    MD_out.addLabel(EMDL_IMAGE_COORD_Y);
    MD_out.addLabel(EMDL_ORIENT_ROT);
    MD_out.addLabel(EMDL_ORIENT_TILT);
    MD_out.addLabel(EMDL_ORIENT_PSI);
    MD_out.addLabel(EMDL_ORIENT_ORIGIN_X);
    MD_out.addLabel(EMDL_ORIENT_ORIGIN_Y);
    MD_out.addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);

	// Calculate all coordinates for helical segments
	nr_segments = 0;
	rot = psi = xoff = yoff = now_x = now_y = 0.;
	tilt = 90.;
    step = nr_asu * rise_ang / pixel_size_ang;
    tube_id = 0;

    for (tube_id = 0; tube_id < x1_coord_list.size(); tube_id++)
    {
    	x1 = x1_coord_list[tube_id];
    	y1 = y1_coord_list[tube_id];
    	x2 = x2_coord_list[tube_id];
    	y2 = y2_coord_list[tube_id];

    	if (fabs(x2 - x1) > (1e-5))
    	{
    		k = (y2 - y1) / (x2 - x1);
        	b = y1 - k * x1;
        	std::cout << "Detected Line #" << tube_id << ": y = (" << k << ") x + (" << b << ")" << std::flush;

        	delta_x = step / (sqrt(k * k + 1.));
        	if (x1 > x2)
        		delta_x *= -1.;
        	delta_y = k * delta_x;
    	}
    	else
    	{
    		std::cout << "Detected Line #" << tube_id << ": x = (" << x1 << ")" << std::flush;
    		delta_x = 0.;
    		delta_y = step;
    		if (y1 > y2)
    			delta_y *= -1.;
    	}
    	std::cout << "   ###   Delta x: " << delta_x << "    Delta y: " << delta_y << std::endl;

    	now_x = x1;
    	now_y = y1;
    	while (1)
    	{
    		now_x += delta_x;
    		now_y += delta_y;
    		if ( ((now_x > x1) && (now_x > x2)) || ((now_x < x1) && (now_x < x2))
    				|| ((now_y > y1) && (now_y > y2)) || ((now_y < y1) && (now_y < y2)) )
    			break;
    		else
    		{
    			// Avoid segments lying on the edges of the micrographs
    			if ( (now_x < half_box_size_pix)
    					|| (now_x > (Xdim - half_box_size_pix))
    					|| (now_y < half_box_size_pix)
    					|| (now_y > (Ydim - half_box_size_pix)) )
    			{
    				continue;
    			}

    			psi = (-180.) * atan(k) / PI;

    			MD_out.addObject();
    	    	MD_out.setValue(EMDL_IMAGE_COORD_X, now_x);
    	    	MD_out.setValue(EMDL_IMAGE_COORD_Y, now_y);
    	    	MD_out.setValue(EMDL_ORIENT_ROT, rot);
    	    	MD_out.setValue(EMDL_ORIENT_TILT, tilt);
    	    	MD_out.setValue(EMDL_ORIENT_PSI, psi);
    	    	MD_out.setValue(EMDL_ORIENT_ORIGIN_X, xoff);
    	    	MD_out.setValue(EMDL_ORIENT_ORIGIN_Y, yoff);
    	    	MD_out.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, (tube_id + 1));

    			nr_segments++;
    		}
    	}
    }

    if (nr_segments < 1)
    	std::cout << "WARNING: no segments extracted from file '" << fn_in << "'!" << std::endl;
    else
    	std::cout << "Total segments / particles extracted from file '" << fn_in << "' = " << nr_segments << " / " << (nr_segments * nr_asu) << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;

    MD_out.write(fn_out);
    MD_out.clear();

    return nr_segments;
}

void extractCoordsForAllHelicalSegments_Multiple(
		FileName& suffix_in,
		FileName& suffix_out,
		int nr_asu,
		RFLOAT rise_ang,
		RFLOAT pixel_size_ang,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix)
{
	int nr_segments;
	FileName fns_in;
	std::vector<FileName> fn_in_list;

	fns_in = "*" + suffix_in;
	fns_in.globFiles(fn_in_list);
	std::cout << "Number of input files = " << fn_in_list.size() << std::endl;
	if (fn_in_list.size() < 1)
		REPORT_ERROR( (std::string) "helix.cpp::extractCoordsForAllHelicalSegments_Multiple(): No input files are found!");

	nr_segments = 0;
	for (int ii = 0; ii < fn_in_list.size(); ii++)
	{
		FileName fn_out;
		fn_out = fn_in_list[ii].beforeFirstOf(suffix_in) + suffix_out;
		nr_segments += extractCoordsForAllHelicalSegments(fn_in_list[ii], fn_out, nr_asu, rise_ang, pixel_size_ang, Xdim, Ydim, box_size_pix);
	}
	std::cout << "Total number of segments / particles from all files = "
			<< nr_segments << " / " << (nr_segments * nr_asu) << std::endl;
	return;
}

void combineParticlePriorsWithKaiLocalCTF(
		FileName& fn_priors,
		FileName& fn_local_ctf,
		FileName& fn_combined)
{
	MetaDataTable MD_priors, MD_local_ctf;
	std::vector<RFLOAT> x, y, rot, tilt, psi, xoff, yoff;
	RFLOAT _x, _y, _rot, _tilt, _psi, _xoff, _yoff;
	int ii;

	if ( (fn_priors.getFileFormat() != "star")
			|| (fn_local_ctf.getFileFormat() != "star")
			|| (fn_combined.getFileFormat() != "star") )
		REPORT_ERROR( (std::string) "helix.cpp::combineParticlePriorsWithKaiLocalCTF(): MetaDataTable should have .star extension.");
	if ( (fn_priors == fn_local_ctf)
			|| (fn_local_ctf == fn_combined)
			|| (fn_combined == fn_priors) )
		REPORT_ERROR( (std::string) "helix.cpp::combineParticlePriorsWithKaiLocalCTF(): File names must be different.");

	MD_priors.clear();
	MD_local_ctf.clear();
	MD_priors.read(fn_priors);
	MD_local_ctf.read(fn_local_ctf);
	if (MD_priors.numberOfObjects() != MD_local_ctf.numberOfObjects())
		REPORT_ERROR( (std::string) "helix.cpp::combineParticlePriorsWithKaiLocalCTF(): MetaDataTables to be combined are not of the same size.");

	if ( (!MD_priors.containsLabel(EMDL_IMAGE_COORD_X))
			|| (!MD_priors.containsLabel(EMDL_IMAGE_COORD_Y))
			|| (!MD_priors.containsLabel(EMDL_ORIENT_ROT))
			|| (!MD_priors.containsLabel(EMDL_ORIENT_TILT))
			|| (!MD_priors.containsLabel(EMDL_ORIENT_PSI))
			|| (!MD_priors.containsLabel(EMDL_ORIENT_ORIGIN_X))
			|| (!MD_priors.containsLabel(EMDL_ORIENT_ORIGIN_Y))
			|| (!MD_local_ctf.containsLabel(EMDL_MICROGRAPH_NAME))
			|| (!MD_local_ctf.containsLabel(EMDL_IMAGE_COORD_X))
			|| (!MD_local_ctf.containsLabel(EMDL_IMAGE_COORD_Y))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_VOLTAGE))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_DEFOCUSU))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_DEFOCUSV))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_DEFOCUS_ANGLE))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_CS))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_FOM))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_MAGNIFICATION))
			|| (!MD_local_ctf.containsLabel(EMDL_CTF_Q0))
			|| (!MD_local_ctf.containsLabel(EMDL_PARTICLE_AUTOPICK_FOM))
			)
		REPORT_ERROR( (std::string) "helix.cpp::combineParticlePriorsWithKaiLocalCTF(): Labels missing in MetaDataTables.");

	x.clear(); y.clear();
	rot.clear(); tilt.clear(); psi.clear();
	xoff.clear(); yoff.clear();

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_priors)
	{
		MD_priors.getValue(EMDL_IMAGE_COORD_X, _x);
		MD_priors.getValue(EMDL_IMAGE_COORD_Y, _y);
		MD_priors.getValue(EMDL_ORIENT_ROT, _rot);
		MD_priors.getValue(EMDL_ORIENT_TILT, _tilt);
		MD_priors.getValue(EMDL_ORIENT_PSI, _psi);
		MD_priors.getValue(EMDL_ORIENT_ORIGIN_X, _xoff);
		MD_priors.getValue(EMDL_ORIENT_ORIGIN_Y, _yoff);
		x.push_back(_x);
		y.push_back(_y);
		rot.push_back(_rot);
		tilt.push_back(_tilt);
		psi.push_back(_psi);
		xoff.push_back(_xoff);
		yoff.push_back(_yoff);
	}

	if (!MD_local_ctf.containsLabel(EMDL_ORIENT_ROT))
		MD_local_ctf.addLabel(EMDL_ORIENT_ROT);
	if (!MD_local_ctf.containsLabel(EMDL_ORIENT_TILT))
		MD_local_ctf.addLabel(EMDL_ORIENT_TILT);
	if (!MD_local_ctf.containsLabel(EMDL_ORIENT_PSI))
		MD_local_ctf.addLabel(EMDL_ORIENT_PSI);
	if (!MD_local_ctf.containsLabel(EMDL_ORIENT_ORIGIN_X))
		MD_local_ctf.addLabel(EMDL_ORIENT_ORIGIN_X);
	if (!MD_local_ctf.containsLabel(EMDL_ORIENT_ORIGIN_Y))
		MD_local_ctf.addLabel(EMDL_ORIENT_ORIGIN_Y);
	ii = -1;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_local_ctf)
	{
		ii++;
		MD_local_ctf.getValue(EMDL_IMAGE_COORD_X, _x);
		MD_local_ctf.getValue(EMDL_IMAGE_COORD_Y, _y);
		if ( (fabs(x[ii] - _x) > 1.001) || (fabs(y[ii] - _y) > 1.001) )
			REPORT_ERROR( (std::string) "helix.cpp::combineParticlePriorsWithKaiLocalCTF(): Coordinates from the two MetaDataTables do not match.");
		MD_local_ctf.setValue(EMDL_IMAGE_COORD_X, x[ii]);
		MD_local_ctf.setValue(EMDL_IMAGE_COORD_Y, y[ii]);
		MD_local_ctf.setValue(EMDL_ORIENT_ROT, rot[ii]);
		MD_local_ctf.setValue(EMDL_ORIENT_TILT, tilt[ii]);
		MD_local_ctf.setValue(EMDL_ORIENT_PSI, psi[ii]);
		MD_local_ctf.setValue(EMDL_ORIENT_ORIGIN_X, xoff[ii]);
		MD_local_ctf.setValue(EMDL_ORIENT_ORIGIN_Y, yoff[ii]);
	}

	MD_local_ctf.write(fn_combined);
	return;
}

void addPriorsToParticleDataFile(
		FileName& fn_priors,
		FileName& fn_data,
		FileName& fn_out)
{
	MetaDataTable MD_priors, MD_data;
	bool have_rot, have_tilt, have_dxy;
	std::string str_img_name;
	RFLOAT x1, y1, x2, y2, rot, tilt, psi, dx, dy;
	std::map<std::string, MetaDataContainer*> priors_list;
	std::map<std::string, MetaDataContainer*>::iterator prior_iterator;
	MetaDataContainer* aux;

	if ( (fn_priors.getFileFormat() != "star") || (fn_data.getFileFormat() != "star") || (fn_out.getFileFormat() != "star") )
		REPORT_ERROR("helix.cpp::combineParticlePriorsWithClass2DDataStar(): MetaDataTable should have .star extension.");
	if ( (fn_priors == fn_data) || (fn_data == fn_out) || (fn_out == fn_priors) )
		REPORT_ERROR("helix.cpp::combineParticlePriorsWithClass2DDataStar(): File names must be different.");

	std::cout << " Loading input files..." << std::endl;
	MD_priors.clear();
	MD_data.clear();
	MD_priors.read(fn_priors);
	MD_data.read(fn_data);

	if ( (!MD_priors.containsLabel(EMDL_IMAGE_COORD_X))
			|| (!MD_priors.containsLabel(EMDL_IMAGE_COORD_Y))
			|| (!MD_priors.containsLabel(EMDL_IMAGE_NAME))
			|| (!MD_priors.containsLabel(EMDL_ORIENT_PSI))

			|| (!MD_data.containsLabel(EMDL_IMAGE_COORD_X))
			|| (!MD_data.containsLabel(EMDL_IMAGE_COORD_Y))
			|| (!MD_data.containsLabel(EMDL_IMAGE_NAME))
			|| (!MD_data.containsLabel(EMDL_CTF_DEFOCUSU))
			|| (!MD_data.containsLabel(EMDL_CTF_DEFOCUSV))
			|| (!MD_data.containsLabel(EMDL_CTF_DEFOCUS_ANGLE))
			|| (!MD_data.containsLabel(EMDL_CTF_VOLTAGE))
			|| (!MD_data.containsLabel(EMDL_CTF_CS))
			|| (!MD_data.containsLabel(EMDL_CTF_Q0))
			|| (!MD_data.containsLabel(EMDL_CTF_MAGNIFICATION))
			|| (!MD_data.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE)) )
	{
		REPORT_ERROR("helix.cpp::combineParticlePriorsWithClass2DDataStar(): Labels missing in MetaDataTables.");
	}

	have_rot = have_tilt = have_dxy = true;
	if (!MD_priors.containsLabel(EMDL_ORIENT_ROT))
		have_rot = false;
	if (!MD_priors.containsLabel(EMDL_ORIENT_TILT))
		have_tilt = false;
	if ( (!MD_priors.containsLabel(EMDL_ORIENT_ORIGIN_X))
			|| (!MD_priors.containsLabel(EMDL_ORIENT_ORIGIN_Y)) )
		have_dxy = false;

	std::cout << " Number of segments in the prior / data files = "
			<< MD_priors.numberOfObjects() << " / " << MD_data.numberOfObjects() << std::endl;
	std::cout << " Adding orientational and translational priors to the data file..." << std::endl;

	priors_list.clear();
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_priors)
	{
		aux = MD_priors.getObject();
		MD_priors.getValue(EMDL_IMAGE_NAME, str_img_name);
		priors_list.insert(std::pair<std::string, MetaDataContainer*>(str_img_name, aux));
	}

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_data)
	{
		MD_data.getValue(EMDL_IMAGE_NAME, str_img_name);
		prior_iterator = priors_list.find(str_img_name);
		if (prior_iterator != priors_list.end())
		{
			aux = prior_iterator->second;

			aux->getValue(EMDL_IMAGE_COORD_X, x1);
			aux->getValue(EMDL_IMAGE_COORD_Y, y1);
			MD_data.getValue(EMDL_IMAGE_COORD_X, x2);
			MD_data.getValue(EMDL_IMAGE_COORD_Y, y2);
			if ( (fabs(x1 - x2) > 1.) || (fabs(y1 - y2) > 1.) )
				REPORT_ERROR("helix.cpp::combineParticlePriorsWithClass2DDataStar(): Some pairs of particles have the same EMDL_IMAGE_NAME values but different XY coordinates. Check if the two files come from the same set of particles.");
			else
			{
				rot = dx = dy = 0.;
				tilt = 90.;
				if (have_rot)
					aux->getValue(EMDL_ORIENT_ROT, rot);
				if (have_tilt)
					aux->getValue(EMDL_ORIENT_TILT, tilt);
				aux->getValue(EMDL_ORIENT_PSI, psi);
				if (have_dxy)
				{
					aux->getValue(EMDL_ORIENT_ORIGIN_X, dx);
					aux->getValue(EMDL_ORIENT_ORIGIN_Y, dy);
				}

				MD_data.setValue(EMDL_ORIENT_ROT, rot);
				MD_data.setValue(EMDL_ORIENT_TILT, tilt);
				MD_data.setValue(EMDL_ORIENT_PSI, psi);
				MD_data.setValue(EMDL_ORIENT_ORIGIN_X, dx);
				MD_data.setValue(EMDL_ORIENT_ORIGIN_Y, dy);
			}
		}
		else
		{
			std::cout << "  EMDL_IMAGE_NAME = " << str_img_name << std::endl;
			REPORT_ERROR("helix.cpp::combineParticlePriorsWithClass2DDataStar(): Priors are missing for some particles in the data file.");
		}
	}

	std::cout << " Writing output file..." << std::endl;
	MD_data.write(fn_out);
	std::cout << " Done!" << std::endl;

	return;
}

void combineParticlePriorsWithKaiLocalCTF_Multiple(
		std::string& suffix_priors,
		std::string& suffix_local_ctf,
		std::string& suffix_combined)
{
	FileName fns_priors;
	std::vector<FileName> fn_priors_list;

	if ( (suffix_priors == suffix_local_ctf)
			|| (suffix_priors == suffix_combined)
			|| (suffix_combined == suffix_priors) )
		REPORT_ERROR( (std::string) "helix.cpp::combineParticlePriorsWithKaiLocalCTF_Multiple(): File names error!");

	fns_priors = "*" + suffix_priors;
	fns_priors.globFiles(fn_priors_list);
	std::cout << "Number of input files = " << fn_priors_list.size() << std::endl;
	if (fn_priors_list.size() < 1)
		REPORT_ERROR( (std::string) "helix.cpp::combineParticlePriorsWithKaiLocalCTF_Multiple(): No input files are found!");

	for (int ii = 0; ii < fn_priors_list.size(); ii++)
	{
		FileName fn_local_ctf, fn_combined;
		fn_local_ctf = fn_priors_list[ii].beforeFirstOf(suffix_priors) + suffix_local_ctf;
		fn_combined = fn_priors_list[ii].beforeFirstOf(suffix_priors) + suffix_combined;
		combineParticlePriorsWithKaiLocalCTF(fn_priors_list[ii], fn_local_ctf, fn_combined);
	}
	return;
}

void setNullAlignmentPriorsInDataStar(
		FileName& fn_in,
		FileName& fn_out,
		bool rewrite_tilt,
		bool rewrite_psi)
{
	MetaDataTable MD;
	bool rot, tilt, psi, xoff, yoff;
	if ( (fn_in.getFileFormat() != "star")
			|| (fn_out.getFileFormat() != "star") )
		REPORT_ERROR( (std::string) "helix.cpp::addNullAlignmentPriorsToDataStar(): MetaDataTable should have .star extension.");
	if (fn_in == fn_out)
		REPORT_ERROR( (std::string) "helix.cpp::addNullAlignmentPriorsToDataStar(): File names must be different.");

	MD.read(fn_in);
	rot = MD.containsLabel(EMDL_ORIENT_ROT);
	tilt = MD.containsLabel(EMDL_ORIENT_TILT);
	psi = MD.containsLabel(EMDL_ORIENT_PSI);
	xoff = MD.containsLabel(EMDL_ORIENT_ORIGIN_X);
	yoff = MD.containsLabel(EMDL_ORIENT_ORIGIN_Y);
	if (!rot)
		MD.addLabel(EMDL_ORIENT_ROT);
	if (!tilt)
		MD.addLabel(EMDL_ORIENT_TILT);
	if (!psi)
		MD.addLabel(EMDL_ORIENT_PSI);
	if (!xoff)
		MD.addLabel(EMDL_ORIENT_ORIGIN_X);
	if (!yoff)
		MD.addLabel(EMDL_ORIENT_ORIGIN_Y);
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		//if (!rot)
			MD.setValue(EMDL_ORIENT_ROT, 0.);
		if ( (rewrite_tilt) || (!tilt) )
			MD.setValue(EMDL_ORIENT_TILT, 90.);
		if ( (rewrite_psi) || (!psi) )
			MD.setValue(EMDL_ORIENT_PSI, 0.);
		//if ( (rewrite_old_values) || (!xoff) )
			MD.setValue(EMDL_ORIENT_ORIGIN_X, 0.);
		//if ( (rewrite_old_values) || (!yoff) )
			MD.setValue(EMDL_ORIENT_ORIGIN_Y, 0.);
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
		REPORT_ERROR( (std::string) "helix.cpp::removeBadTiltParticlesFromDataStar(): Max deviations of tilt angles from 90 degree should be in the range of 0~89 degrees.");
	if ( (fn_in.getFileFormat() != "star")
			|| (fn_out.getFileFormat() != "star") )
		REPORT_ERROR( (std::string) "helix.cpp::removeBadTiltParticlesFromDataStar(): MetaDataTable should have .star extension.");
	if (fn_in == fn_out)
		REPORT_ERROR( (std::string) "helix.cpp::removeBadTiltParticlesFromDataStar(): File names must be different.");

	MD_in.clear();
	MD_out.clear();
	MD_in.read(fn_in);
	if (!MD_in.containsLabel(EMDL_ORIENT_TILT))
		REPORT_ERROR( (std::string) "helix.cpp::removeBadTiltParticlesFromDataStar(): Input .star file contains no tilt angles.");

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

int transformXimdispHelicalCoordsToStarFile(
		FileName& fn_in,
		FileName& fn_out,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT box_size_pix)
{
	if ( (box_size_pix < 2) || (Xdim < box_size_pix) || (Ydim < box_size_pix))
		REPORT_ERROR( (std::string) "helix.cpp::transformXimdispHelicalCoordsToStarFile(): Wrong dimensions or box size!");
	if (fn_out.isStarFile() == false)
		REPORT_ERROR( (std::string) "helix.cpp::transformXimdispHelicalCoordsToStarFile(): Output should be a star file!");

	char tmpstr[1000];
	int nr_particles_on_edges, nr_particles;
	RFLOAT x, y, psi_deg, half_box_size_pix;
	MetaDataTable MD;
	std::ifstream fin;
	std::string line;
	std::vector<std::string> words;

	half_box_size_pix = box_size_pix / 2.;
	fin.open(fn_in.c_str(), std::ios_base::in);
	if (fin.fail())
		REPORT_ERROR( (std::string) "helix.cpp::transformXimdispHelicalCoordsToStarFile(): Cannot open input file!");
	MD.clear();
	MD.addLabel(EMDL_IMAGE_COORD_X);
	MD.addLabel(EMDL_IMAGE_COORD_Y);
	MD.addLabel(EMDL_ORIENT_ROT);
	MD.addLabel(EMDL_ORIENT_TILT);
	MD.addLabel(EMDL_ORIENT_PSI);
	MD.addLabel(EMDL_ORIENT_ORIGIN_X);
	MD.addLabel(EMDL_ORIENT_ORIGIN_Y);
	nr_particles_on_edges = nr_particles = 0;
	getline(fin, line, '\n');
	while (getline(fin, line, '\n'))
	{
		words.clear();
		tokenize(line, words);
		if (words.size() != 3)
			REPORT_ERROR( (std::string) "helix.cpp::transformXimdispHelicalCoordsToStarFile(): Every line in input Ximdisp file should only contain x, y and psi angle!");
		x = textToFloat(words[0]);
		y = textToFloat(words[1]);
		psi_deg = textToFloat(words[2]);

		// Avoid segments lying on the edges of the micrographs
		if ( (x < half_box_size_pix)
				|| (x > (Xdim - half_box_size_pix))
				|| (y < half_box_size_pix)
				|| (y > (Ydim - half_box_size_pix)) )
		{
			nr_particles_on_edges++;
			continue;
		}

		nr_particles++;
		MD.addObject();
		MD.setValue(EMDL_IMAGE_COORD_X, x);
		MD.setValue(EMDL_IMAGE_COORD_Y, y);
		MD.setValue(EMDL_ORIENT_ROT, 0.);
		MD.setValue(EMDL_ORIENT_TILT, 90.);
		MD.setValue(EMDL_ORIENT_PSI, -psi_deg);
		MD.setValue(EMDL_ORIENT_ORIGIN_X, 0.);
		MD.setValue(EMDL_ORIENT_ORIGIN_Y, 0.);
	}
	MD.write(fn_out);
	std::cout << "Input Ximdisp .coords = " << fn_in.c_str()
			<< ", output STAR file = " << fn_out.c_str()
			<< ", size of micrograph = " << Xdim << " * " << Ydim
			<< ", box size = " << box_size_pix << ", "
			<< nr_particles_on_edges << " particles excluded, "
			<< nr_particles << " particles left." << std::endl;
	return nr_particles;
}

void transformXimdispHelicalCoordsToStarFile_Multiple(
		FileName& suffix_coords,
		FileName& suffix_out,
		RFLOAT Xdim,
		RFLOAT Ydim,
		RFLOAT boxsize)
{
	int nr_particles;
	FileName fns_coords;
	std::vector<FileName> fn_coords_list;

	fns_coords = "*" + suffix_coords;
	fns_coords.globFiles(fn_coords_list);
	std::cout << "Number of input files = " << fn_coords_list.size() << std::endl;
	if (fn_coords_list.size() < 1)
		REPORT_ERROR( (std::string) "helix.cpp::transformXimdispHelicalCoordsToStarFile_Multiple(): No input files are found!");

	nr_particles = 0;
	for (int ii = 0; ii < fn_coords_list.size(); ii++)
	{
		FileName fn_out;
		fn_out = fn_coords_list[ii].beforeFirstOf(suffix_coords) + suffix_out;
		nr_particles += transformXimdispHelicalCoordsToStarFile(fn_coords_list[ii], fn_out, Xdim, Ydim, boxsize);
	}
	std::cout << "Number of particles = " << nr_particles << std::endl;
	return;
}

/*
void divideHelicalSegmentsFromMultipleMicrographsIntoHalves(
		FileName& fn_in,
		FileName& fn_out)
{
	int id_mic, nr_segments_subset1, nr_segments_subset2;
	RFLOAT ratio;
	std::string mic_name;
	MetaDataTable MD;
	std::vector<std::string> mic_names;
	std::vector<int> mic_nr_segments;
	std::vector<int> mic_subset;

	if (fn_in == fn_out)
		REPORT_ERROR("helix.cpp::divideHelicalSegmentsFromMultipleMicrographsIntoHalves(): File names error!");

	MD.clear();
	MD.read(fn_in);

	if (!MD.containsLabel(EMDL_MICROGRAPH_NAME))
		REPORT_ERROR("helix.cpp::divideHelicalSegmentsFromMultipleMicrographsIntoHalves(): Input MetadataTable should contain rlnMicrographName!");
	if (!MD.containsLabel(EMDL_PARTICLE_RANDOM_SUBSET))
		MD.addLabel(EMDL_PARTICLE_RANDOM_SUBSET);

	mic_names.clear();
	mic_nr_segments.clear();
	mic_subset.clear();
	id_mic = nr_segments_subset1 = nr_segments_subset2 = 0;

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		MD.getValue(EMDL_MICROGRAPH_NAME, mic_name);
		for (id_mic = (mic_names.size() - 1); id_mic >= 0; id_mic--)
		{
			if (mic_names[id_mic] == mic_name)
			{
				mic_nr_segments[id_mic]++;
				MD.setValue(EMDL_PARTICLE_RANDOM_SUBSET, id_mic);
				break;
			}
		}
		if (id_mic < 0)
		{
			mic_names.push_back(mic_name);
			mic_nr_segments.push_back(1);
			mic_subset.push_back(1);
			id_mic = (mic_names.size() - 1);
			MD.setValue(EMDL_PARTICLE_RANDOM_SUBSET, id_mic);
		}
	}

	// Randomise

	for (id_mic = 0; id_mic < mic_names.size(); id_mic++)
	{
		if (nr_segments_subset1 < nr_segments_subset2)
		{
			mic_subset[id_mic] = 1;
			nr_segments_subset1 += mic_nr_segments[id_mic];
		}
		else
		{
			mic_subset[id_mic] = 2;
			nr_segments_subset2 += mic_nr_segments[id_mic];
		}
	}

	if ( (nr_segments_subset1 < 1) || (nr_segments_subset2 < 1) )
		REPORT_ERROR("helix.cpp::divideHelicalSegmentsFromMultipleMicrographsIntoHalves(): Number of helical segments from one of the two half sets is 0!");
	ratio = (RFLOAT(nr_segments_subset1) / RFLOAT(nr_segments_subset2));
	if ( (ratio > 1.5) || ( (1. / ratio) > 1.5) )
		REPORT_ERROR("helix.cpp::divideHelicalSegmentsFromMultipleMicrographsIntoHalves(): Numbers of helical segments from two half sets are extremely unbalanced!");

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		MD.getValue(EMDL_PARTICLE_RANDOM_SUBSET, id_mic);
		MD.setValue(EMDL_PARTICLE_RANDOM_SUBSET, mic_subset[id_mic]);
	}

	std::cout << " Helical segments in two half sets = " << nr_segments_subset1 << ", " << nr_segments_subset2 << std::endl;

	MD.write(fn_out);

	return;
}
*/

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
	particle_radius_max_pix = (int)(CEIL(particle_diameter_pix / 2.)) + 1;

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

void makeHelicalReconstructionStarFileFromSingle2DClassAverage(
		FileName& fn_in_class2D,
		FileName& fn_out_star,
		RFLOAT pixel_size_A,
		RFLOAT twist_deg,
		RFLOAT rise_A,
		RFLOAT tilt_deg,
		RFLOAT psi_deg,
		int nr_projections)
{
	RFLOAT twist_deg_acc, offset_pix_acc, rise_pix;
	MetaDataTable MD;

	if (fn_out_star.getExtension() != "star")
		REPORT_ERROR("helix.cpp::makeReconstructionStarFileFromSingle2DClassAverage: Output file should be in .star format!");
	if (pixel_size_A < 0.001)
		REPORT_ERROR("helix.cpp::makeReconstructionStarFileFromSingle2DClassAverage: Pixel size (in Angstroms) should be larger than 0.001!");
	if ((rise_A / pixel_size_A) < 0.001)
		REPORT_ERROR("helix.cpp::makeReconstructionStarFileFromSingle2DClassAverage: Wrong helical rise!");
	if (nr_projections < 1)
		REPORT_ERROR("helix.cpp::makeReconstructionStarFileFromSingle2DClassAverage: Number of projections should be larger than 1!");

	rise_pix = rise_A / pixel_size_A;

	MD.clear();
	MD.addLabel(EMDL_IMAGE_NAME);
	MD.addLabel(EMDL_ORIENT_ROT);
	MD.addLabel(EMDL_ORIENT_TILT);
	MD.addLabel(EMDL_ORIENT_PSI);
	MD.addLabel(EMDL_ORIENT_ORIGIN_X);
	MD.addLabel(EMDL_ORIENT_ORIGIN_Y);

	twist_deg_acc = -twist_deg;
	offset_pix_acc = -rise_pix;
	for (int ii = 0; ii < nr_projections; ii++)
	{
		twist_deg_acc += twist_deg;
		offset_pix_acc += rise_pix;

		MD.addObject();
		MD.setValue(EMDL_IMAGE_NAME, fn_in_class2D);
		MD.setValue(EMDL_ORIENT_ROT, twist_deg_acc);
		MD.setValue(EMDL_ORIENT_TILT, tilt_deg);
		MD.setValue(EMDL_ORIENT_PSI, psi_deg);
		MD.setValue(EMDL_ORIENT_ORIGIN_X, 0.);
		MD.setValue(EMDL_ORIENT_ORIGIN_Y, offset_pix_acc);
	}

	MD.write(fn_out_star);
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

void combineStarFiles(
		FileName& fn_in)
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

void sortHelicalTubeID(
		MetaDataTable& MD)
{
	std::string str_particle_fullname, str_particle_name, str_comment, str_particle_id, str_tube_id;
	int int_tube_id;
	bool contain_helicalTubeID;

	if ( (!MD.containsLabel(EMDL_IMAGE_NAME)) || (MD.containsLabel(EMDL_COMMENT)) )
		REPORT_ERROR("helix.cpp::sortHelicalTubeID: MetaDataTable should contain rlnImageName but not rlnComment!");
	contain_helicalTubeID = MD.containsLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);

	MD.addLabel(EMDL_COMMENT);
	int_tube_id = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		MD.getValue(EMDL_IMAGE_NAME, str_particle_fullname);
		if (contain_helicalTubeID)
			MD.getValue(EMDL_PARTICLE_HELICAL_TUBE_ID, int_tube_id);

		str_particle_name = str_particle_fullname.substr(str_particle_fullname.find("@") + 1);
		str_particle_id = str_particle_fullname.substr(0, str_particle_fullname.find("@"));
		str_comment = str_particle_name + integerToString(int_tube_id, 10) + str_particle_id;

		// DEBUG
		//std::cout << str_comment << std::endl;

		MD.setValue(EMDL_COMMENT, str_comment);
	}
	MD.newSort(EMDL_COMMENT);
	MD.deactivateLabel(EMDL_COMMENT);

	return;
}

void simulateHelicalSegments(
		FileName& fn_out,
		int nr_subunits,
		int nr_asu,
		int nr_tubes,
		RFLOAT twist_deg,
		RFLOAT sigma_psi,
		RFLOAT sigma_tilt)
{
	int nr_segments, tube_id;
	RFLOAT rot, psi, tilt;
	MetaDataTable MD;
	FileName fn_mic;

	// TODO: raise error if nr_asu<0 or too big, n too small!
	if ( (nr_tubes < 2) || (nr_subunits < 100) || (nr_asu < 1)
			|| (((nr_subunits / nr_asu) / nr_tubes) < 5) || ((nr_subunits / nr_asu) > 999999) )
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Errors in the number of tubes, asymmetrical units or total subunits!");
	if ( (sigma_psi < 0.) || (sigma_psi > 10.) || (sigma_tilt < 0.) || (sigma_tilt > 10.) )
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Errors in sigma_psi or sigma_tilt!");
	if ( (fabs(twist_deg) < 0.001) || (fabs(twist_deg) > 180.))
		REPORT_ERROR("helix.cpp::simulateHelicalSegments: Error in helical twist!");

	nr_segments = nr_subunits / nr_asu;

	MD.clear();
    MD.addLabel(EMDL_ORIENT_ROT);
    MD.addLabel(EMDL_ORIENT_TILT);
    MD.addLabel(EMDL_ORIENT_PSI);
    MD.addLabel(EMDL_ORIENT_ORIGIN_X);
    MD.addLabel(EMDL_ORIENT_ORIGIN_Y);
    MD.addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);
    MD.addLabel(EMDL_IMAGE_NAME);
    MD.addLabel(EMDL_MICROGRAPH_NAME);

    tube_id = 0;
	for (int id = 0; id < nr_segments; id++)
	{
		if ( ( (id % (nr_segments / nr_tubes)) == 0 )
				&& ( (nr_segments - id) >= (nr_segments / nr_tubes) ) )
		{
			tube_id++;
			tilt = rnd_unif(85., 95.);
			rot = rnd_unif(0.01, 359.99);
			psi = rnd_unif(-179.99, 179.99);
		}

		rot += twist_deg * ((RFLOAT)(nr_asu));
		rot = realWRAP(rot, -180., 180.);

		MD.addObject();
    	MD.setValue(EMDL_ORIENT_ROT, rot);
    	MD.setValue(EMDL_ORIENT_TILT, realWRAP(tilt + rnd_gaus(0., sigma_tilt), 0., 180.) );
    	MD.setValue(EMDL_ORIENT_PSI, realWRAP(psi + rnd_gaus(0., sigma_psi), -180., 180.) );
    	MD.setValue(EMDL_ORIENT_ORIGIN_X, 0.);
    	MD.setValue(EMDL_ORIENT_ORIGIN_Y, 0.);
    	MD.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, tube_id);
    	fn_mic.compose((id + 1), "M.mrc");
    	MD.setValue(EMDL_IMAGE_NAME, fn_mic);
    	MD.setValue(EMDL_MICROGRAPH_NAME, (std::string)("M.mrc"));
	}
	MD.write(fn_out);
	return;
};
