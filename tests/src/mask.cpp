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
#include <omp.h>
#include "src/mask.h"

// https://stackoverflow.com/questions/48273190/undefined-symbol-error-for-stdstringempty-c-standard-method-linking-error/48273604#48273604
#if defined(__APPLE__)
// explicit instantiation of std::string needed, otherwise we get a linker error on osx
// thats a bug in libc++, because of interaction with __attribute__ ((__visibility__("hidden"), __always_inline__)) in std::string
template class std::basic_string<char>;
#endif

// Workaround for compiler versions before 2018 update 2
#ifdef __INTEL_COMPILER
# if (__INTEL_COMPILER<1800)
#  pragma optimize ("", off)
# endif
# if (__INTEL_COMPILER==1800)
#  if (__INTEL_COMPILER_UPDATE<2)
#   pragma optimize ("", off)
#  endif
# endif
#endif
// Mask out corners outside sphere (replace by average value)
// Apply a soft mask (raised cosine with cosine_width pixels width)
void softMaskOutsideMap(MultidimArray<RFLOAT> &vol, RFLOAT radius, RFLOAT cosine_width, MultidimArray<RFLOAT> *Mnoise)
{

	vol.setXmippOrigin();
	RFLOAT r, radius_p, raisedcos, sum_bg = 0., sum = 0.;
	if (radius < 0)
		radius = (RFLOAT)XSIZE(vol)/2.;
	radius_p = radius + cosine_width;


	if (Mnoise == NULL)
	{
		// Calculate average background value
		FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
		{
			r = sqrt((RFLOAT)(k*k + i*i + j*j));
			if (r < radius)
				continue;
			else if (r > radius_p)
			{
				sum    += 1.;
				sum_bg += A3D_ELEM(vol, k, i, j);
			}
			else
			{
				raisedcos = 0.5 + 0.5 * cos(PI * (radius_p - r) / cosine_width );
				sum += raisedcos;
				sum_bg += raisedcos * A3D_ELEM(vol, k, i, j);
			}
		}
		sum_bg /= sum;
	}

	// Apply noisy or average background value
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
	{
		r = sqrt((RFLOAT)(k*k + i*i + j*j));
		if (r < radius)
		{
			continue;
		}
		else if (r > radius_p)
		{
			A3D_ELEM(vol, k, i, j) = (Mnoise == NULL) ? sum_bg : A3D_ELEM(*Mnoise, k, i, j);
		}
		else
		{
			raisedcos = 0.5 + 0.5 * cos(PI * (radius_p - r) / cosine_width );
			RFLOAT add = (Mnoise == NULL) ?  sum_bg : A3D_ELEM(*Mnoise, k, i, j);
			A3D_ELEM(vol, k, i, j) = (1 - raisedcos) * A3D_ELEM(vol, k, i, j) + raisedcos * add;
		}
	}

}

// May27,2015 - Shaoda, Helical refinement
void softMaskOutsideMapForHelix(
		MultidimArray<RFLOAT> &vol,
		RFLOAT psi_deg,
		RFLOAT tilt_deg,
		RFLOAT mask_sphere_radius_pix,
		RFLOAT mask_cyl_radius_pix,
		RFLOAT cosine_width,
		MultidimArray<RFLOAT> *Mnoise)
{
	Matrix1D<RFLOAT> coords;
	Matrix2D<RFLOAT> A;
	RFLOAT sum_bg, sum, R1, R2, D1, D2, r, d, noise_w, noise_w1, noise_w2, noise_val;
	int dim = vol.getDim();
	int boxsize = -1;

	// Center the box
	vol.setXmippOrigin();
	// Dimension of a particle (box) should be 2 or 3
	if ( (dim != 2) && (dim != 3) )
	{
		REPORT_ERROR("mask.cpp::softMaskOutsideMapForHelix(): Dimension of particles should be 2 or 3!");
		return;
	}
	// Check the shape of Mnoise
	if ( (Mnoise != NULL) && ((*Mnoise).sameShape(vol) == false) )
	{
		REPORT_ERROR("mask.cpp::softMaskOutsideMapForHelix(): Input particle and Mnoise should have same shape!");
		return;
	}
	// Box size is the minimum value of the 2 or 3 dimensions
	boxsize = (XSIZE(vol) < YSIZE(vol)) ? XSIZE(vol) : YSIZE(vol);
	// If it is a 2D particle, tilt angle does not apply
	if (dim == 2)
		tilt_deg = 0.;
	else
		boxsize = (boxsize < ZSIZE(vol)) ? boxsize : ZSIZE(vol);
	boxsize = boxsize / 2 - ((boxsize + 1) % 2);

	// Diameter of the cylindrical mask around the helix should not exceed the box size, otherwise noise cannot be estimated
	if ( (cosine_width < 0.)
			|| (mask_sphere_radius_pix < 1.) || (mask_sphere_radius_pix > boxsize)
			|| (mask_cyl_radius_pix < 1.) || (mask_cyl_radius_pix > boxsize)
			|| (mask_sphere_radius_pix < mask_cyl_radius_pix) )
		REPORT_ERROR("mask.cpp::softMaskOutsideMapForHelix(): Invalid radii of spherical and cylindrical masks or soft cosine widths!");

	// Spherical mask: 0 < R1 < R2
	R1 = mask_sphere_radius_pix;
	R2 = R1 + cosine_width;
	// Cylindrical mask: 0 < D1 < D2
	D1 = mask_cyl_radius_pix;
	D2 = D1 + cosine_width;

	// Init coords
	coords.clear();
	coords.resize(3);
	coords.initZeros();

	// Init rotational matrix A
	A.clear();
	A.resize(3, 3);

	// Rotate the particle (helical axes are X and Z for 2D and 3D segments respectively)
	Euler_angles2matrix(0., tilt_deg, psi_deg, A, false);
	// Don't put negative signs before tilt and psi values, use 'transpose' instead
	A = A.transpose();

	// Calculate noise weights for all voxels
	sum_bg = sum = 0.;
	if (Mnoise == NULL)
	{
		FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
		{
			// X, Y, Z coordinates
			if (dim == 3)
				ZZ(coords) = ((RFLOAT)(k));
			else
				ZZ(coords) = 0.;
			YY(coords) = ((RFLOAT)(i));
			XX(coords) = ((RFLOAT)(j));
			// Rotate
			coords = A * coords;

			// Distance from the point to helical axis (perpendicular to X axis)
			if (dim == 3)
				d = sqrt(YY(coords) * YY(coords) + XX(coords) * XX(coords));
			else
				d = ABS(YY(coords));
			if (d > D2) // Noise areas (get values for noise estimations)
			{
				sum_bg += A3D_ELEM(vol, k, i, j);
				sum += 1.;
			}
			else if (d > D1) // Edges of noise areas (get values and weights for noise estimations)
			{
				noise_w = 0.5 + 0.5 * cos(PI * (D2 - d) / cosine_width );
				sum_bg += noise_w * A3D_ELEM(vol, k, i, j);
				sum += noise_w;
			}
		}
		// Test (this should not happen)
		if (sum < 0.00001)
			REPORT_ERROR("mask.cpp::softMaskOutsideMapForHelix(): No background (noise) areas found in this particle!");
		sum_bg /= sum;
	}

	// Apply noisy or average background value
	noise_val = sum_bg;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
	{
		// X, Y, Z coordinates
		if (dim == 3)
			ZZ(coords) = ((RFLOAT)(k));
		else
			ZZ(coords) = 0.;
		YY(coords) = ((RFLOAT)(i));
		XX(coords) = ((RFLOAT)(j));

		// Rotate
		coords = A * coords;

		// Distance from the point to helical axis (perpendicular to X axis)
		if (dim == 3)
			d = sqrt(YY(coords) * YY(coords) + XX(coords) * XX(coords));
		else
			d = ABS(YY(coords));

		// Distance from the origin
		r = (RFLOAT)(i * i + j * j);
		if (dim == 3)
			r += (RFLOAT)(k * k);
		r = sqrt(r);

		// Info areas
		if ( (r < R1) && (d < D1) )
			continue;

		if (Mnoise != NULL)
			noise_val = A3D_ELEM(*Mnoise, k, i, j);

		if ( (r > R2) || (d > D2) )  // Noise areas, fill in background values
			A3D_ELEM(vol, k, i, j) = noise_val;
		else // Edges of info areas
		{
			noise_w1 = noise_w2 = 0.;
			if (r > R1)
				noise_w1 = 0.5 + 0.5 * cos(PI * (R2 - r) / cosine_width );
			if (d > D1)
				noise_w2 = 0.5 + 0.5 * cos(PI * (D2 - d) / cosine_width );
			noise_w = (noise_w1 > noise_w2) ? (noise_w1) : (noise_w2);
			A3D_ELEM(vol, k, i, j) = (1. - noise_w) * A3D_ELEM(vol, k, i, j) + noise_w * noise_val;
		}
	}
	return;
}

// Workaround for compiler versions before 2018 update 2
#ifdef __INTEL_COMPILER
# if (__INTEL_COMPILER<1800)
#  pragma optimize ("", on)
# endif
# if (__INTEL_COMPILER==1800)
#  if (__INTEL_COMPILER_UPDATE<2)
#   pragma optimize ("", on)
#  endif
# endif
#endif

void softMaskOutsideMap(MultidimArray<RFLOAT> &vol, MultidimArray<RFLOAT> &msk, bool invert_mask)
{

	if (msk.computeMax() > 1. || msk.computeMin() < 0.)
	{
		std::cerr << " msk.computeMax()= " << msk.computeMax() << " msk.computeMin()= " << msk.computeMin() << std::endl;
		REPORT_ERROR("ERROR: Values in the solvent mask should be between zero and one.");
	}
	if (!(msk.sameShape(vol)))
		REPORT_ERROR("ERROR: Solvent mask does not have the same size as the reference vol.");

	// Replace solvent by the average value in the solvent region
	RFLOAT sum = 0.;
	RFLOAT sum_bg = 0.;
	RFLOAT solv;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(msk)
	{
		solv = (invert_mask) ? DIRECT_A3D_ELEM(msk, k, i, j) : 1. - DIRECT_A3D_ELEM(msk, k, i, j);
		sum    += solv;
		sum_bg += solv * DIRECT_A3D_ELEM(vol, k, i, j);
	}
	sum_bg /= sum;

	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(msk)
	{
		solv = (invert_mask) ? DIRECT_A3D_ELEM(msk, k, i, j) : 1. - DIRECT_A3D_ELEM(msk, k, i, j);
		DIRECT_A3D_ELEM(vol, k, i, j) = ( 1. - solv) * DIRECT_A3D_ELEM(vol, k, i, j) + solv * sum_bg;
	}


}

void autoMask(MultidimArray<RFLOAT> &img_in, MultidimArray<RFLOAT> &msk_out,
		RFLOAT ini_mask_density_threshold, RFLOAT extend_ini_mask, RFLOAT width_soft_mask_edge, bool verb, int n_threads)

{
	MultidimArray<RFLOAT> msk_cp;
	int barstep, update_bar, totalbar;

	// Resize output mask
	img_in.setXmippOrigin();
	msk_out.clear();
	msk_out.resize(img_in);

	// A. Calculate initial binary mask based on density threshold
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img_in)
	{
		if (DIRECT_MULTIDIM_ELEM(img_in, n) >= ini_mask_density_threshold)
			DIRECT_MULTIDIM_ELEM(msk_out, n) = 1.;
		else
			DIRECT_MULTIDIM_ELEM(msk_out, n) = 0.;
	}

	// B. extend/shrink initial binary mask. To save memory store a temporary copy of Im in I1
	if (extend_ini_mask > 0. || extend_ini_mask < 0.)
	{
		if (verb)
		{
			if (extend_ini_mask > 0.)
				std::cout << "== Extending initial binary mask ..." << std::endl;
			else
				std::cout << "== Shrinking initial binary mask ..." << std::endl;
			init_progress_bar(MULTIDIM_SIZE(img_in) / n_threads);
			barstep = MULTIDIM_SIZE(img_in) / 120 / n_threads;
			update_bar = 0;
			totalbar =0;
		}

		int extend_size = ABS(CEIL(extend_ini_mask));
		RFLOAT extend_ini_mask2 = extend_ini_mask * extend_ini_mask;
		msk_cp = msk_out;
		if (extend_ini_mask > 0.)
		{
			#pragma omp parallel for num_threads(n_threads)
			FOR_ALL_ELEMENTS_IN_ARRAY3D(msk_cp)
			{
				// only extend zero values to 1.
				if (A3D_ELEM(msk_cp, k, i, j) < 0.001)
				{
					bool already_done = false;
					for (long int kp = k - extend_size; kp <= k + extend_size; kp++)
					{
						for (long int ip = i - extend_size; ip <= i + extend_size; ip++)
						{
							for (long int jp = j - extend_size; jp <= j + extend_size; jp++)
							{
								if ((kp >= STARTINGZ(msk_cp) && kp <= FINISHINGZ(msk_cp)) &&
									(ip >= STARTINGY(msk_cp) && ip <= FINISHINGY(msk_cp)) &&
									(jp >= STARTINGX(msk_cp) && jp <= FINISHINGX(msk_cp)))
								{
									// only check distance if neighbouring Im() is one
									if (A3D_ELEM(msk_cp, kp, ip, jp) > 0.999)
									{
										RFLOAT r2 = (RFLOAT)( (kp-k)*(kp-k) + (ip-i)*(ip-i)+ (jp-j)*(jp-j) );
										// Set original voxel to 1 if a neghouring with Im()=1 is within distance extend_ini_mask
										if (r2 < extend_ini_mask2)
										{
											A3D_ELEM(msk_out, k, i, j) = 1.;
											already_done = true;
										}
									}
								}
								if (already_done) break;
							}
							if (already_done) break;
						}
						if (already_done) break;
					}
				}
				if (verb && omp_get_thread_num() == 0)
				{
					if (update_bar > barstep)
					{
						update_bar = 0;
						progress_bar(totalbar);
					}
					update_bar++;
					totalbar++;
				}
			}
		}
		else
		{
			#pragma omp parallel for num_threads(n_threads)
			FOR_ALL_ELEMENTS_IN_ARRAY3D(msk_cp)
			{
				// only extend one values to zero.
				if (A3D_ELEM(msk_cp, k, i, j) > 0.999)
				{
					bool already_done = false;
					for (long int kp = k - extend_size; kp <= k + extend_size; kp++)
					{
						for (long int ip = i - extend_size; ip <= i + extend_size; ip++)
						{
							for (long int jp = j - extend_size; jp <= j + extend_size; jp++)
							{
								if ((kp >= STARTINGZ(msk_cp) && kp <= FINISHINGZ(msk_cp)) &&
									(ip >= STARTINGY(msk_cp) && ip <= FINISHINGY(msk_cp)) &&
									(jp >= STARTINGX(msk_cp) && jp <= FINISHINGX(msk_cp)))
								{
									// only check distance if neighbouring Im() is one
									if (A3D_ELEM(msk_cp, kp, ip, jp) < 0.001)
									{
										RFLOAT r2 = (RFLOAT)( (kp-k)*(kp-k) + (ip-i)*(ip-i)+ (jp-j)*(jp-j) );
										// Set original voxel to 1 if a neghouring with Im()=1 is within distance extend_ini_mask
										if (r2 < extend_ini_mask2)
										{
											A3D_ELEM(msk_out, k, i, j) = 0.;
											already_done = true;
										}
									}
								}
								if (already_done) break;
							}
							if (already_done) break;
						}
						if (already_done) break;
					}
				}
				if (verb && omp_get_thread_num() == 0)
				{
					if (update_bar > barstep)
					{
						update_bar = 0;
						progress_bar(totalbar);
					}
					update_bar++;
					totalbar++;
				}
			}
		}
		if (verb)
			progress_bar(MULTIDIM_SIZE(msk_out) / n_threads);
	}

	if (width_soft_mask_edge > 0.)
	{
		if (verb)
		{
			std::cout << "== Making a soft edge on the extended mask ..." << std::endl;
			init_progress_bar(MULTIDIM_SIZE(msk_out) / n_threads);
			barstep = MULTIDIM_SIZE(msk_out) / 120 / n_threads;
			update_bar = 0;
			totalbar =0;
		}
		// C. Make a soft edge to the mask
		// Note that the extended mask is now in I1, and we'll put the soft-edge mask again into Im

		msk_cp = msk_out;
		int extend_size = CEIL(width_soft_mask_edge);
		RFLOAT width_soft_mask_edge2 = width_soft_mask_edge * width_soft_mask_edge;
		#pragma omp parallel for num_threads(n_threads)
		FOR_ALL_ELEMENTS_IN_ARRAY3D(msk_cp)
		{
			// only extend zero values to values between 0 and 1.
			if (A3D_ELEM(msk_cp, k, i, j) < 0.001)
			{
				RFLOAT min_r2 = 9999.;
				for (long int kp = k - extend_size; kp <= k + extend_size; kp++)
				{
					for (long int ip = i - extend_size; ip <= i + extend_size; ip++)
					{
						for (long int jp = j - extend_size; jp <= j + extend_size; jp++)
						{
							if ((kp >= STARTINGZ(msk_cp) && kp <= FINISHINGZ(msk_cp)) &&
								(ip >= STARTINGY(msk_cp) && ip <= FINISHINGY(msk_cp)) &&
								(jp >= STARTINGX(msk_cp) && jp <= FINISHINGX(msk_cp)))
							{
								// only update distance to a neighbouring msk_cp is one
								if (A3D_ELEM(msk_cp, kp, ip, jp) > 0.999)
								{
									RFLOAT r2 = (RFLOAT)( (kp-k)*(kp-k) + (ip-i)*(ip-i)+ (jp-j)*(jp-j) );
									// Set original voxel to 1 if a neghouring with Im()=1 is within distance extend_ini_mask
									if (r2 < min_r2)
										min_r2 = r2;
								}
							}
						}
					}
				}
				if (min_r2 < width_soft_mask_edge2)
				{
					A3D_ELEM(msk_out, k, i, j) = 0.5 + 0.5 * cos( PI * sqrt(min_r2) / width_soft_mask_edge);
				}
			}
			if (verb && omp_get_thread_num() == 0)
			{
				if (update_bar > barstep)
				{
					update_bar = 0;
					progress_bar(totalbar);
				}
				update_bar++;
				totalbar++;
			}
		}
		if (verb)
			progress_bar(MULTIDIM_SIZE(msk_cp) / n_threads);
	}

}

void raisedCosineMask(MultidimArray<RFLOAT> &mask, RFLOAT radius, RFLOAT radius_p, int x, int y, int z)
{
	mask.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
	{
		// calculate distance from the origin
		RFLOAT d = sqrt((RFLOAT)((z-k)*(z-k) + (y-i)*(y-i) + (x-j)*(x-j)));
		if (d > radius_p)
			A3D_ELEM(mask, k, i, j) = 0.;
		else if (d < radius)
			A3D_ELEM(mask, k, i, j) = 1.;
		else
			A3D_ELEM(mask, k, i, j) = 0.5 - 0.5 * cos(PI * (radius_p - d) / (radius_p - radius));
	}
}

void raisedCrownMask(MultidimArray<RFLOAT> &mask, RFLOAT inner_radius, RFLOAT outer_radius, RFLOAT width, RFLOAT x, RFLOAT y, RFLOAT z)
{
	RFLOAT inner_border = inner_radius - width;
	RFLOAT outer_border = outer_radius + width;

	mask.setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY3D(mask)
	{
		RFLOAT d = sqrt((RFLOAT)((z-k)*(z-k) + (y-i)*(y-i) + (x-j)*(x-j)));
		if (d < inner_border)
			A3D_ELEM(mask, k, i, j) = 0.;
		else if (d < inner_radius)
			A3D_ELEM(mask, k, i, j) = 0.5  - 0.5 * cos(PI * (d - inner_border) / width);
		else if (d < outer_radius) 
			A3D_ELEM(mask, k, i, j) = 1.; 
		else if (d < outer_border)
			A3D_ELEM(mask, k, i, j) = 0.5 - 0.5 * cos(PI * (outer_border - d) / width);
		else
			A3D_ELEM(mask, k, i, j) = 0.;
	}
}


