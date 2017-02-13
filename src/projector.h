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
#ifndef __PROJECTOR_H
#define __PROJECTOR_H

#include "src/fftw.h"
#include "src/multidim_array.h"
#include "src/image.h"

#define NEAREST_NEIGHBOUR 0
#define TRILINEAR 1
#define CONVOLUTE_BLOB 2

#define FORWARD_PROJECTION 0
#define BACKWARD_PROJECTION 1

#define ACT_ON_DATA 0
#define ACT_ON_WEIGHT 1

class Projector
{
public:
	// The Fourier-space image data array
    MultidimArray<Complex > data;

    // Only points within this many pixels from the origin (in the original size) will be interpolated
    int r_max;

    // Radius of sphere within TRILINEAR interpolation will be used in NEAREST_NEIGHBOUR interpolator
    int r_min_nn;

    // Original size of the real-space map
    int ori_size;

    // Padded size of the map in Fourier-space
    int pad_size;

    // Interpolation scheme (TRILINEAR or NEAREST_NEIGHBOUR, for BackProjector also CONVOLUTE_BLOB)
    int interpolator;

    // Oversample FT by padding in real space
    float padding_factor;

    // Dimension of the reference (currently allowed 2 or 3)
    int ref_dim;

    // Dimension of the projections (2 or 3)
    int data_dim;

public:

    /** Empty constructor
     *
     * A default Projector is created.
     *
     * @code
     * Projector PPref;
     * @endcode
     */
    Projector()
    {
    	clear();
    }

   /** Constructor with parameters
     *
     * A default Projector is created.
     *
     * @code
     * Projector PPref(ori_size, NEAREST_NEIGHBOUR);
     * @endcode
     */
    Projector(int _ori_size, int _interpolator = TRILINEAR, float _padding_factor_3d = 2., int _r_min_nn = 10, int _data_dim = 2)
    {

    	// Store original dimension
    	ori_size = _ori_size;

    	// Padding factor for the map
    	padding_factor = _padding_factor_3d;

    	// Interpolation scheme
    	interpolator = _interpolator;

    	// Minimum radius for NN interpolation
    	r_min_nn = _r_min_nn;

    	// Dimension of the projections
    	data_dim = _data_dim;
    }

    /** Copy constructor
     *
     * The created Projector is a perfect copy of the input array but with a
     * different memory assignment.
     *
     * @code
     * Projector V2(V1);
     * @endcode
     */
	Projector(const Projector& op)
    {
		clear();
        *this = op;
    }

	/** Assignment.
     *
     * You can build as complex assignment expressions as you like. Multiple
     * assignment is allowed.
     */
	Projector& operator=(const Projector& op)
    {
        if (&op != this)
        {
        	data = op.data;
        	ori_size = op.ori_size;
        	pad_size = op.pad_size;
        	r_max = op.r_max;
        	r_min_nn = op.r_min_nn;
        	interpolator = op.interpolator;
        	padding_factor = op.padding_factor;
        	ref_dim = op.ref_dim;
        	data_dim  = op.data_dim;
        }
        return *this;
    }

	/** Destructor
      *
      * Clears everything
      *
      * @code
      * FourierInterpolator fourint;
      * @endcode
      */
     ~Projector()
     {
         clear();
     }

   /** Clear.
     * Initialize everything to back to default and empty arrays
     */
    void clear()
    {
    	data.clear();
    	r_max = r_min_nn = interpolator = ref_dim = data_dim = pad_size = 0;
    	padding_factor = 0.;
    }

    /*
     * Resize data array to the given size
     */
    void initialiseData(int current_size = -1);

    /*
     * Initialise data array to all zeros
     */
    void initZeros(int current_size = -1);

    /*
     *  Only get the size of the data array
     */
    long int getSize();

    /* ** Prepares a 3D map for taking slices in its 3D Fourier Transform
    *
    * This routine does the following:
    * 1. It pads the input map with zeros to obtain oversampling in the Fourier transform
    * 2. It does the Fourier transform
    * 3. It sets values beyond Nyquist for images of current_size to zero in the transform and windows the transform at max_r+1
    * Depending on whether 2D or 3D Fourier Transforms will be extracted, the map is normalized internally in a different manner
    *
    */
   void computeFourierTransformMap(MultidimArray<RFLOAT> &vol_in, MultidimArray<RFLOAT> &power_spectrum, int current_size = -1, int nr_threads = 1, bool do_gridding = true, bool do_heavy = true);

   /* This is experimental: apply a mask in Fourier-space to focus refinements on certain Fourier components
    * mask_r_min and mask_r_max are the radii of the lowest and highest frequencies (only keep crown inside)
    * mask_ang is the opening angle along z (only really useful for helices, I guess)
    */
   void applyFourierMask(int mask_r_min = 0, int mask_r_max = -1, RFLOAT mask_ang = 0.);

   /* Because we interpolate in Fourier space to make projections and/or reconstructions, we have to correct
    * the real-space maps by dividing them by the Fourier Transform of the interpolator
    * Note these corrections are made on the not-oversampled, i.e. originally sized real-space map
    */
   void griddingCorrect(MultidimArray<RFLOAT> &vol_in);

   /*
	* Get a 2D Fourier Transform from the 2D or 3D data array
	* Depending on the dimension of the map, this will be a projection or a rotation operation
	*/
	void get2DFourierTransform(MultidimArray<Complex > &img_out, Matrix2D<RFLOAT> &A, bool inv)
	{
		// Rotation of a 3D Fourier Transform
		if (data_dim == 3)
		{
			if (ref_dim != 3)
				REPORT_ERROR("Projector::get3DFourierTransform%%ERROR: Dimension of the data array should be 3");
			rotate3D(img_out, A, inv);
		}
		else
		{
			switch (ref_dim)
			{
			case 2:
			   rotate2D(img_out, A, inv);
			   break;
			case 3:
			   project(img_out, A, inv);
			   break;
			default:
			   REPORT_ERROR("Projector::get2DSlice%%ERROR: Dimension of the data array should be 2 or 3");
			}
		}
	}

	/*
	* Get a 2D slice from the 3D map (forward projection)
	*/
	void project(MultidimArray<Complex > &img_out, Matrix2D<RFLOAT> &A, bool inv);

	/*
	* Get an in-plane rotated version of the 2D map (mere interpolation)
	*/
	void rotate2D(MultidimArray<Complex > &img_out, Matrix2D<RFLOAT> &A, bool inv);

	/*
	* Get a rotated version of the 3D map (mere interpolation)
	*/
	void rotate3D(MultidimArray<Complex > &img_out, Matrix2D<RFLOAT> &A, bool inv);


};




#endif
