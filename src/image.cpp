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


// Get size of datatype
unsigned long  gettypesize(DataType type)
{
    unsigned long   size;

    switch ( type ) {
        case UChar: case SChar:  size = sizeof(char); break;
        case UShort: case Short: size = sizeof(short); break;
        case UInt:	 case Int:   size = sizeof(int); break;
        case Float:              size = sizeof(float); break;
        case Double:             size = sizeof(double); break;
        case Bool:				  size = sizeof(bool); break;
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
    return Short;
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
  else REPORT_ERROR("datatypeString2int; unknown datatype");


}

// Some image-specific operations
void normalise(Image<double> &I, int bg_radius, double white_dust_stddev, double black_dust_stddev)
{
	int bg_radius2 = bg_radius * bg_radius;
	double avg, stddev;

	if (2*bg_radius > XSIZE(I()))
		REPORT_ERROR("normalise ERROR: 2*bg_radius is larger than image size!");

	// Calculate initial avg and stddev values
	calculateBackgroundAvgStddev(I, avg, stddev, bg_radius);

	// Remove white and black noise
	if (white_dust_stddev > 0.)
		removeDust(I, true, white_dust_stddev, avg, stddev);
	if (black_dust_stddev > 0.)
		removeDust(I, false, black_dust_stddev, avg, stddev);

	// If some dust was removed: recalculate avg and stddev
	if (white_dust_stddev > 0. || black_dust_stddev > 0.)
		calculateBackgroundAvgStddev(I, avg, stddev, bg_radius);


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

void calculateBackgroundAvgStddev(Image<double> &I, double &avg, double &stddev, int bg_radius)
{
	int bg_radius2 = bg_radius * bg_radius;
	double n = 0.;
	avg = 0.;
	stddev = 0.;

	// Calculate avg in the background pixels
	FOR_ALL_ELEMENTS_IN_ARRAY3D(I())
	{
		if (k*k + i*i + j*j > bg_radius2)
		{
			avg += A3D_ELEM(I(), k, i, j);
			n += 1.;
		}
	}
	avg /= n;

	// Calculate stddev in the background pixels
	FOR_ALL_ELEMENTS_IN_ARRAY3D(I())
	{
		if (k*k + i*i + j*j > bg_radius2)
		{
			double aux = A3D_ELEM(I(), k, i, j) - avg;
			stddev += aux * aux;
		}
	}
	stddev = sqrt(stddev/n);
}

void removeDust(Image<double> &I, bool is_white, double thresh, double avg, double stddev)
{
	FOR_ALL_ELEMENTS_IN_ARRAY3D(I())
	{
		double aux =  A3D_ELEM(I(), k, i, j);
		if (is_white && aux - avg > thresh * stddev)
			A3D_ELEM(I(), k, i, j) = rnd_gaus(avg, stddev);
		else if (!is_white && aux - avg < -thresh * stddev)
			A3D_ELEM(I(), k, i, j) = rnd_gaus(avg, stddev);
	}
}

void invert_contrast(Image<double> &I)
{
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I())
	{
		DIRECT_MULTIDIM_ELEM(I(), n) *= -1;
	}
}

void rescale(Image<double> &I, int mysize)
{
	int olddim = XSIZE(I());

	resizeMap(I(), mysize);

	// Also modify the scale in the MDmainheader (if present)
	double oldscale, newscale;
    if (I.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_X, oldscale))
    {
    	newscale = oldscale * (double)olddim / (double)mysize;
    	I.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, newscale);
    }
    if (I.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_Y, oldscale))
    {
    	newscale = oldscale * (double)olddim / (double)mysize;
    	I.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, newscale);
    }
    if (I().getDim() == 3 && I.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_Z, oldscale) )
    {
    	newscale = oldscale * (double)olddim / (double)mysize;
    	I.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, newscale);
    }

}

void rewindow(Image<double> &I, int mysize)
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

