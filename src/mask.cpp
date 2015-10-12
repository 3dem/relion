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
#include "src/mask.h"

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
		RFLOAT ini_mask_density_threshold, RFLOAT extend_ini_mask, RFLOAT width_soft_mask_edge, bool verb)

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
			init_progress_bar(MULTIDIM_SIZE(img_in));
			barstep = MULTIDIM_SIZE(img_in)/120;
			update_bar = 0;
			totalbar =0;
		}

		int extend_size = ABS(CEIL(extend_ini_mask));
		RFLOAT extend_ini_mask2 = extend_ini_mask * extend_ini_mask;
		msk_cp = msk_out;
		if (extend_ini_mask > 0.)
		{
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
				if (verb)
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
				if (verb)
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
			progress_bar(MULTIDIM_SIZE(msk_out));
	}

	if (width_soft_mask_edge > 0.)
	{
		if (verb)
		{
			std::cout << "== Making a soft edge on the extended mask ..." << std::endl;
			init_progress_bar(MULTIDIM_SIZE(msk_out));
			barstep = MULTIDIM_SIZE(msk_out)/120;
			update_bar = 0;
			totalbar =0;
		}
		// C. Make a soft edge to the mask
		// Note that the extended mask is now in I1, and we'll put the soft-edge mask again into Im

		msk_cp = msk_out;
		int extend_size = CEIL(width_soft_mask_edge);
		RFLOAT width_soft_mask_edge2 = width_soft_mask_edge * width_soft_mask_edge;
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
			if (verb)
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
			progress_bar(MULTIDIM_SIZE(msk_cp));
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
