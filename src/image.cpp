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
#include "src/image.h"

//#define DEBUG_REGULARISE_HELICAL_SEGMENTS

// Get size of datatype
unsigned long  gettypesize(DataType type)
{
	unsigned long	size;

	switch ( type ) {
		case UChar: case SChar:  size = sizeof(char); break;
		case UShort: case SShort: size = sizeof(short); break;
		case UInt: case Int:     size = sizeof(int); break;
		case Float:              size = sizeof(float); break;
		case Double:             size = sizeof(RFLOAT); break;
		case Boolean:            size = sizeof(bool); break;
		case Float16:            size = sizeof(short); break;
		case UHalf: REPORT_ERROR("Logic error: UHalf (4-bit) needs special consideration. Don't use this function."); break;
		default: size = 0;
	}

	return(size);
}

int datatypeString2Int(std::string s)
{
	toLower(s);
	if (!strcmp(s.c_str(),"uchar"))
	{
		return UChar;
	}
	else if (!strcmp(s.c_str(),"ushort"))
	{
		return UShort;
	}
	else if (!strcmp(s.c_str(),"short"))
	{
		return SShort;
	}
	else if (!strcmp(s.c_str(),"uint"))
	{
		return UInt;
	}
	else if (!strcmp(s.c_str(),"int"))
	{
		return Int;
	}
	else if (!strcmp(s.c_str(),"float"))
	{
		return Float;
	}
	else if (!strcmp(s.c_str(),"float16"))
	{
		return Float16;
	}
	else REPORT_ERROR("datatypeString2int; unknown datatype");
}

// Some image-specific operations
void normalise(
               Image<RFLOAT> &I,
               int bg_radius,
               RFLOAT white_dust_stddev,
               RFLOAT black_dust_stddev,
               bool do_ramp,
               bool is_helical_segment,
               RFLOAT helical_mask_tube_outer_radius_pix,
               RFLOAT tilt_deg,
               RFLOAT psi_deg)
{
	RFLOAT avg, stddev;

	if (2*bg_radius > XSIZE(I()))
		REPORT_ERROR("normalise ERROR: 2*bg_radius is larger than image size!");

	if ( (is_helical_segment) && ( (2 * (helical_mask_tube_outer_radius_pix + 1)) > XSIZE(I()) ) )
		REPORT_ERROR("normalise ERROR: Diameter of helical tube is larger than image size!");

	if (is_helical_segment)
	{
		if (I().getDim() == 2)
			tilt_deg = 0.;
	}

	if (white_dust_stddev > 0. || black_dust_stddev > 0.)
	{
		// Calculate initial avg and stddev values
		calculateBackgroundAvgStddev(I, avg, stddev, bg_radius,
		                             is_helical_segment, helical_mask_tube_outer_radius_pix, tilt_deg, psi_deg);

		// Remove white and black noise
		if (white_dust_stddev > 0.)
			removeDust(I, true, white_dust_stddev, avg, stddev);
		if (black_dust_stddev > 0.)
			removeDust(I, false, black_dust_stddev, avg, stddev);
	}

	if (do_ramp)
		subtractBackgroundRamp(I, bg_radius,
		                       is_helical_segment, helical_mask_tube_outer_radius_pix, tilt_deg, psi_deg);

	// Calculate avg and stddev (also redo if dust was removed!)
	calculateBackgroundAvgStddev(I, avg, stddev, bg_radius,
	                             is_helical_segment, helical_mask_tube_outer_radius_pix, tilt_deg, psi_deg);

	if (stddev < 1e-10)
	{
		std::cerr << " WARNING! Stddev of image " << I.name() << " is zero! Skipping normalisation..." << std::endl;
	}
	else
	{
		// Subtract avg and divide by stddev for all pixels
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I())
			DIRECT_MULTIDIM_ELEM(I(), n) = (DIRECT_MULTIDIM_ELEM(I(), n) - avg) / stddev;
	}
}

void calculateBackgroundAvgStddev(Image<RFLOAT> &I,
                                  RFLOAT &avg,
                                  RFLOAT &stddev,
                                  int bg_radius,
                                  bool is_helical_segment,
                                  RFLOAT helical_mask_tube_outer_radius_pix,
                                  RFLOAT tilt_deg,
                                  RFLOAT psi_deg)
{
	int bg_radius2 = bg_radius * bg_radius;
	RFLOAT sum, sum2, n, val, d;
	sum = sum2 = n = 0.;
	avg = stddev = 0.;

	if (is_helical_segment)
	{
		int dim = I().getDim();
		if ( (dim != 2) && (dim != 3) )
			REPORT_ERROR("image.cpp::calculateBackgroundAvgStddev(): 2D or 3D image is required!");
		if (dim == 2)
			tilt_deg = 0.;

		Matrix1D<RFLOAT> coords;
		Matrix2D<RFLOAT> A;

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
		// Refer to the code in calculateBackgroundAvgStddev() for 3D implementation

#ifdef DEBUG_REGULARISE_HELICAL_SEGMENTS
		FileName fn_test;
		Image<RFLOAT> img_test;
		int angle = ROUND(fabs(psi_deg));
		fn_test = integerToString(angle);
		if (psi_deg < 0.)
			fn_test = fn_test.addExtension("neg");
		fn_test = fn_test.addExtension("mrc");
		img_test.clear();
		img_test().resize(I());
		img_test().initZeros();
		std::cout << "FileName = " << fn_test.c_str() << std::endl;
#endif

		// Calculate avg in the background pixels
		FOR_ALL_ELEMENTS_IN_ARRAY3D(I())
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

			if (d > helical_mask_tube_outer_radius_pix)
			{
				val = A3D_ELEM(I(), k, i, j);
				sum += val;
				sum2 += val * val;
				n += 1.;

#ifdef DEBUG_REGULARISE_HELICAL_SEGMENTS
				A3D_ELEM(img_test(), k, i, j) = 1.; // Mark bg pixels as 1, others as 0
#endif
			}
		}
		if (n < 0.9)
		{
			REPORT_ERROR("image.cpp::calculateBackgroundAvgStddev(): No pixels in background are found. Radius of helical mask is too large.");
		}

		avg = sum / n;
		stddev = sqrt( (sum2 / n) - (avg * avg) );

#ifdef DEBUG_REGULARISE_HELICAL_SEGMENTS
		img_test.write(fn_test);
#endif
	}
	else
	{
		// Calculate avg in the background pixels
		FOR_ALL_ELEMENTS_IN_ARRAY3D(I())
		{
			if ( (k*k + i*i + j*j) > bg_radius2)
			{
				val = A3D_ELEM(I(), k, i, j);
				sum += val;
				sum2 += val * val;
				n += 1.;
			}
		}
		if (n < 0.9)
		{
			REPORT_ERROR("image.cpp::calculateBackgroundAvgStddev(): No pixels in background are found. Radius of circular mask is too large.");
		}

		avg = sum / n;
		stddev = sqrt( (sum2 / n) - (avg * avg) );
	}

	return;
}


void subtractBackgroundRamp(
                            Image<RFLOAT> &I,
                            int bg_radius,
                            bool is_helical_segment,
                            RFLOAT helical_mask_tube_outer_radius_pix,
                            RFLOAT tilt_deg,
                            RFLOAT psi_deg)
{

	int bg_radius2 = bg_radius * bg_radius;
	fit_point3D point;
	std::vector<fit_point3D> allpoints;
	RFLOAT pA, pB, pC, avgbg, stddevbg, minbg, maxbg;

	if (I().getDim() == 3)
		REPORT_ERROR("ERROR %% calculateBackgroundRamp is not implemented for 3D data!");

	if (is_helical_segment)  // not implemented for 3D data
	{
		Matrix1D<RFLOAT> coords;
		Matrix2D<RFLOAT> A;
		if (I().getDim() == 2)
			tilt_deg = 0.;

		// Init coords
		coords.clear();
		coords.resize(3);
		coords.initZeros();

		// Init rotational matrix A
		A.clear();
		A.resize(3, 3);

		// Rotate the particle (helical axes are X and Z for 2D and 3D segments respectively)
		// Since Z = 0, tilt_deg does not matter
		Euler_angles2matrix(0., tilt_deg, psi_deg, A, false);
		// Don't put negative signs before tilt and psi values, use 'transpose' instead
		A = A.transpose();

		FOR_ALL_ELEMENTS_IN_ARRAY2D(I()) // not implemented for 3D data
		{
			ZZ(coords) = 0.;
			YY(coords) = ((RFLOAT)(i));
			XX(coords) = ((RFLOAT)(j));
			// Rotate
			coords = A * coords;
			if (ABS(YY(coords)) > helical_mask_tube_outer_radius_pix) // not implemented for 3D data
			{
				point.x = j;
				point.y = i;
				point.z = A2D_ELEM(I(), i, j);
				point.w = 1.;
				allpoints.push_back(point);
			}
		}
		if (allpoints.size() < 5)
			REPORT_ERROR("image.cpp::subtractBackgroundRamp(): Less than 5 pixels in background are found. Radius of helical mask is too large.");
	}
	else
	{
		FOR_ALL_ELEMENTS_IN_ARRAY2D(I())
		{
			if (i*i + j*j > bg_radius2)
			{
				point.x = j;
				point.y = i;
				point.z = A2D_ELEM(I(), i, j);
				point.w = 1.;
				allpoints.push_back(point);
			}
		}
	}

	fitLeastSquaresPlane(allpoints, pA, pB, pC);

	// Substract the plane from the image
	FOR_ALL_ELEMENTS_IN_ARRAY2D(I())
	{
		A2D_ELEM(I(), i, j) -= pA * j + pB * i + pC;
	}

}


void removeDust(Image<RFLOAT> &I, bool is_white, RFLOAT thresh, RFLOAT avg, RFLOAT stddev)
{
	FOR_ALL_ELEMENTS_IN_ARRAY3D(I())
	{
		RFLOAT aux =  A3D_ELEM(I(), k, i, j);
		if (is_white && aux - avg > thresh * stddev)
			A3D_ELEM(I(), k, i, j) = rnd_gaus(avg, stddev);
		else if (!is_white && aux - avg < -thresh * stddev)
			A3D_ELEM(I(), k, i, j) = rnd_gaus(avg, stddev);
	}
}

void invert_contrast(Image<RFLOAT> &I)
{
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I())
	{
		DIRECT_MULTIDIM_ELEM(I(), n) *= -1;
	}
}

void rescale(Image<RFLOAT> &I, int mysize)
{
	int olddim = XSIZE(I());

	resizeMap(I(), mysize);

	// Also modify the scale in the MDmainheader (if present)
	RFLOAT oldscale, newscale;
	if (I.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_X, oldscale))
	{
		newscale = oldscale * (RFLOAT)olddim / (RFLOAT)mysize;
		I.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, newscale);
	}
	if (I.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_Y, oldscale))
	{
		newscale = oldscale * (RFLOAT)olddim / (RFLOAT)mysize;
		I.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, newscale);
	}
	if (I().getDim() == 3 && I.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_Z, oldscale) )
	{
		newscale = oldscale * (RFLOAT)olddim / (RFLOAT)mysize;
		I.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, newscale);
	}
}

void rewindow(Image<RFLOAT> &I, int mysize)
{
	// Check 2D or 3D dimensionality
	if (I().getDim() == 2)
	{
		I().window(FIRST_XMIPP_INDEX(mysize), FIRST_XMIPP_INDEX(mysize),
				   LAST_XMIPP_INDEX(mysize),  LAST_XMIPP_INDEX(mysize));
	}
	else if (I().getDim() == 3)
	{
		I().window(FIRST_XMIPP_INDEX(mysize), FIRST_XMIPP_INDEX(mysize), FIRST_XMIPP_INDEX(mysize),
				   LAST_XMIPP_INDEX(mysize),  LAST_XMIPP_INDEX(mysize),  LAST_XMIPP_INDEX(mysize));
	}
}

void getImageContrast(MultidimArray<RFLOAT> &image, RFLOAT &minval, RFLOAT &maxval, RFLOAT &sigma_contrast)
{
	// First check whether to apply sigma-contrast, i.e. set minval and maxval to the mean +/- sigma_contrast times the stddev
	bool redo_minmax = (sigma_contrast > 0. || minval != maxval);

	if (sigma_contrast > 0. || minval == maxval)
	{
		RFLOAT avg, stddev;
		image.computeStats(avg, stddev, minval, maxval);
		if (sigma_contrast > 0.)
		{
			minval = avg - sigma_contrast * stddev;
			maxval = avg + sigma_contrast * stddev;
			redo_minmax = true;
		}
	}

	if (redo_minmax)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(image)
		{
			RFLOAT val = DIRECT_MULTIDIM_ELEM(image, n);
			if (val > maxval)
				DIRECT_MULTIDIM_ELEM(image, n) = maxval;
			else if (val < minval)
				DIRECT_MULTIDIM_ELEM(image, n) = minval;
		}
	}
}
