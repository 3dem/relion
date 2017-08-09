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

/*
 * backprojector.h
 *
 *  Created on: 24 Aug 2010
 *      Author: scheres
 */

#ifndef BACKPROJECTOR_H_
#define BACKPROJECTOR_H_

#include "src/projector.h"
#include "src/mask.h"
#include "src/tabfuncs.h"
#include "src/symmetries.h"

class BackProjector: public Projector
{
public:
	// For backward projection: sum of weights
	MultidimArray<RFLOAT> weight;

	// Tabulated blob values
	TabFtBlob tab_ftblob;

	// Symmetry object
    SymList SL;

    // Helical twist
    RFLOAT twist;

    // Helical rise
    RFLOAT rise;

    // Helical range
    int H;

    // Skip the iterative gridding part of the reconstruction
    bool skip_gridding;


public:

    /** Empty constructor
	 *
	 * A BackProjector is created.
	 *
	 * @code
	 * BackProjector BPref(orisize, 3, "d2");
	 * @endcode
	 */
	BackProjector(int _ori_size, int _ref_dim, FileName fn_sym,
			      int _interpolator = TRILINEAR, float _padding_factor_3d = 2, int _r_min_nn = 10,
			      int _blob_order = 0, RFLOAT _blob_radius = 1.9, RFLOAT _blob_alpha = 15, int _data_dim = 2, bool _skip_gridding = false)
	{
    	// Store original dimension
    	ori_size = _ori_size;

    	// Set dimensionality of the references
    	ref_dim = _ref_dim;

    	// and of the data
    	data_dim = _data_dim;

    	// Skip gridding
    	skip_gridding = _skip_gridding;

    	// Set the symmetry object
    	SL.read_sym_file(fn_sym);

    	// Padding factor for the map
    	padding_factor = _padding_factor_3d;

    	// Interpolation scheme
    	interpolator = _interpolator;

    	// Minimum radius for NN interpolation
    	r_min_nn = _r_min_nn;

    	// Precalculate tabulated ftblob values
    	//tab_ftblob.initialise(_blob_radius * padding_factor, _blob_alpha, _blob_order, 10000);
    	// Sjors 8aug2017: try to fix problems with pad1 reconstrctions
    	tab_ftblob.initialise(_blob_radius * 2., _blob_alpha, _blob_order, 10000);

	}

    /** Copy constructor
     *
     * The created BackProjector is a perfect copy of the input array but with a
     * different memory assignment.
     *
     * @code
     * BackProjector V2(V1);
     * @endcode
     */
	BackProjector(const BackProjector& op)
    {
		clear();
        *this = op;
    }

	/** Assignment.
     *
     * You can build as complex assignment expressions as you like. Multiple
     * assignment is allowed.
     */
	BackProjector& operator=(const BackProjector& op)
    {
        if (&op != this)
        {
         	// Projector stuff (is this necessary in C++?)
        	data = op.data;
        	ori_size = op.ori_size;
        	pad_size = op.pad_size;
        	r_max = op.r_max;
        	r_min_nn = op.r_min_nn;
        	interpolator = op.interpolator;
        	padding_factor = op.padding_factor;
        	ref_dim = op.ref_dim;
        	data_dim = op.data_dim;
        	skip_gridding = op.skip_gridding;
        	// BackProjector stuff
        	weight = op.weight;
        	tab_ftblob = op.tab_ftblob;
        	SL = op.SL;
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
	~BackProjector()
	{
		clear();
	}

	void clear()
	{
		skip_gridding = false;
		weight.clear();
		Projector::clear();
	}

	// Initialise data and weight arrays to the given size and set all values to zero
	void initialiseDataAndWeight(int current_size = -1);

	// Initialise data and weight arrays to the given size and set all values to zero
	void initZeros(int current_size = -1);

	/*
	* Set a 2D Fourier Transform back into the 2D or 3D data array
	* Depending on the dimension of the map, this will be a backprojection or a rotation operation
	*/
	void set2DFourierTransform(const MultidimArray<Complex > &img_in,
							   const Matrix2D<RFLOAT> &A, bool inv,
						       const MultidimArray<RFLOAT> *Mweight = NULL)
	{
		// Back-rotation of a 3D Fourier Transform
		if (img_in.getDim() == 3)
		{
			if (ref_dim != 3)
				REPORT_ERROR("Backprojector::set3DFourierTransform%%ERROR: Dimension of the data array should be 3");
			backrotate3D(img_in, A, inv, Mweight);
		}
		else
		{
			switch (ref_dim)
			{
			case 2:
				backrotate2D(img_in, A, inv, Mweight);
				break;
			case 3:
				backproject(img_in, A, inv, Mweight);
				break;
			default:
				REPORT_ERROR("Backprojector::set2DSlice%%ERROR: Dimension of the data array should be 2 or 3");
			}
		}
	}

	/*
	* Set an in-plane rotated version of the 2D map into the data array (mere interpolation)
	* If a exp_Mweight is given, rather than adding 1 to all relevant pixels in the weight array, we use exp_Mweight
	*/
	void backrotate2D(const MultidimArray<Complex > &img_in,
			          const Matrix2D<RFLOAT> &A, bool inv,
			          const MultidimArray<RFLOAT> *Mweight = NULL);

	/*
	* Set a 3D-rotated version of the 3D map into the data array (mere interpolation)
	* If a exp_Mweight is given, rather than adding 1 to all relevant pixels in the weight array, we use exp_Mweight
	*/
	void backrotate3D(const MultidimArray<Complex > &img_in,
			          const Matrix2D<RFLOAT> &A, bool inv,
			          const MultidimArray<RFLOAT> *Mweight = NULL);

	/*
	* Set a 2D slice in the 3D map (backward projection)
	* If a exp_Mweight is given, rather than adding 1 to all relevant pixels in the weight array, we use exp_Mweight
	*/
	void backproject(const MultidimArray<Complex > &img_in,
			         const Matrix2D<RFLOAT> &A, bool inv,
			         const MultidimArray<RFLOAT> *Mweight = NULL);

	/*
	 * Get only the lowest resolution components from the data and weight array
	 * (to be joined together for two independent halves in order to force convergence in the same orientation)
	 */
	void getLowResDataAndWeight(MultidimArray<Complex > &lowres_data, MultidimArray<RFLOAT> &lowres_weight,
			int lowres_r_max);

	/*
	 * Set only the lowest resolution components from the data and weight array
	 * (to be joined together for two independent halves in order to force convergence in the same orientation)
	 */
	void setLowResDataAndWeight(MultidimArray<Complex > &lowres_data, MultidimArray<RFLOAT> &lowres_weight,
			int lowres_r_max);

	/*
	 *  Get complex array at the original size as the straightforward average
	 *  padding_factor*padding_factor*padding_factor voxels
	 *  This will then be used for FSC calculation between two random halves
	 */
	void getDownsampledAverage(MultidimArray<Complex > &avg);

	/*
	 * From two of the straightforward downsampled averages, calculate an FSC curve
	 */
	void calculateDownSampledFourierShellCorrelation(MultidimArray<Complex > &avg1,
			                                         MultidimArray<Complex > &avg2,
			                                         MultidimArray<RFLOAT> &fsc);

	/* Get the 3D reconstruction
         * If do_map is true, 1 will be added to all weights
         * alpha will contain the noise-reduction spectrum
	*/
	void reconstruct(MultidimArray<RFLOAT> &vol_out,
                     int max_iter_preweight,
                     bool do_map,
                     RFLOAT tau2_fudge,
                     MultidimArray<RFLOAT> &tau2,
                     MultidimArray<RFLOAT> &sigma2,
                     MultidimArray<RFLOAT> &evidence_vs_prior,
                     MultidimArray<RFLOAT> &fourier_coverage,
                     MultidimArray<RFLOAT> fsc,
                     RFLOAT normalise = 1.,
                     bool update_tau2_with_fsc = false,
                     bool is_whole_instead_of_half = false,
                     int nr_threads = 1,
                     int minres_map = -1,
                     bool printTimes= false);


	/*  Enforce Hermitian symmetry, apply helical symmetry as well as point-group symmetry
	 */
	void symmetrise(int nr_helical_asu = 1, RFLOAT helical_twist = 0., RFLOAT helical_rise = 0.);

	/* Enforce hermitian symmetry on data and on weight (all points in the x==0 plane)
	* Because the interpolations are numerical, hermitian symmetry may be broken.
	* Repairing it here gives like a 2-fold averaging correction for interpolation errors...
    */
	void enforceHermitianSymmetry();

	/* Applies helical symmetry. Note that helical_rise is in PIXELS here, as BackProjector doesn't know angpix
	 */
	void applyHelicalSymmetry(int nr_helical_asu = 1, RFLOAT helical_twist = 0., RFLOAT helical_rise = 0.);

	/* Applies the symmetry from the SymList object to the weight and the data array
	 */
	void applyPointGroupSymmetry();


   /* Convolute in Fourier-space with the blob by multiplication in real-space
	 * Note the convlution is done on the complex array inside the transformer object!!
	 */
	void convoluteBlobRealSpace(FourierTransformer &transformer, bool do_mask = false);

	/* Calculate the inverse FFT of Fin and windows the result to ori_size
	 * Also pass the transformer, to prevent making and clearing a new one before clearing the one in reconstruct()
	 */
	void windowToOridimRealSpace(FourierTransformer &transformer, MultidimArray<RFLOAT> &Mout, int nr_threads = 1, bool printTimes = false);

   /*
	* Go from the Projector-centered fourier transform back to FFTW-uncentered one
	*/
   template <typename T>
   void decenter(MultidimArray<T> &Min, MultidimArray<T> &Mout, int my_rmax2)
   {

	   // Mout should already have the right size
	   // Initialize to zero
	   Mout.initZeros();
	   FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Mout)
	   {
		   if (kp*kp + ip*ip + jp*jp <= my_rmax2)
			   DIRECT_A3D_ELEM(Mout, k, i, j) = A3D_ELEM(Min, kp, ip, jp);
	   }
   }

#ifdef RELION_SINGLE_PRECISION
   // Fnewweight needs decentering, but has to be in double-precision for correct calculations!
   template <typename T>
   void decenter(MultidimArray<T> &Min, MultidimArray<double> &Mout, int my_rmax2)
   {

	   // Mout should already have the right size
	   // Initialize to zero
	   Mout.initZeros();
	   FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Mout)
	   {
		   if (kp*kp + ip*ip + jp*jp <= my_rmax2)
                       DIRECT_A3D_ELEM(Mout, k, i, j) = (double)A3D_ELEM(Min, kp, ip, jp);
	   }
   }
#endif

};

#endif /* BACKPROJECTOR_H_ */
