/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
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

#ifndef PAR_FFTW_H
#define PAR_FFTW_H

#include <fftw3.h>
#include "src/multidim_array.h"
#include "src/funcs.h"
#include "src/tabfuncs.h"
#include "src/complex.h"
#include "src/CPlot2D.h"

/*
  Parallelizable version of FourierTransformer from src/fftw.h:

  FFTW plans are managed globally and only one can be computed at a time.
  The original class recomputes its plan each time a new Image is passed
  to its FourierTransform() or inverseFourierTransform().
  As a consequence, when working on sets of multiple images, only one
  FT can run at a time.

  This class only recomputes the plans if the size of the image changes.
  The same plan is reused for different images of the same size.

  Otherwise, the two classes are identical.

          -- J. Zivanov, Feb. 9th 2018
*/


class ParFourierTransformer
{
public:
    /** Real array, in fact a pointer to the user array is stored. */
    MultidimArray<RFLOAT> *fReal;

     /** Complex array, in fact a pointer to the user array is stored. */
    MultidimArray<Complex > *fComplex;

    /** Fourier array  */
    MultidimArray< Complex > fFourier;

#ifdef RELION_SINGLE_PRECISION
    /* fftw Forward plan */
    fftwf_plan fPlanForward;

    /* fftw Backward plan */
    fftwf_plan fPlanBackward;
#else
    /* fftw Forward plan */
    fftw_plan fPlanForward;

    /* fftw Backward plan */
    fftw_plan fPlanBackward;
#endif

    bool plans_are_set;

// Public methods
public:
    /** Default constructor */
    ParFourierTransformer();

    /** Destructor */
    ~ParFourierTransformer();

    /** Copy constructor
     *
     * The created ParFourierTransformer is a perfect copy of the input array but with a
     * different memory assignment.
     *
     */
    ParFourierTransformer(const ParFourierTransformer& op);

    /** Compute the Fourier transform of a MultidimArray, 2D and 3D.
        If getCopy is false, an alias to the transformed data is returned.
        This is a faster option since a copy of all the data is avoided,
        but you need to be careful that an inverse Fourier transform may
        change the data.
        */
    template <typename T, typename T1>
        void FourierTransform(T& v, T1& V, bool getCopy=true)
        {
            setReal(v);
            Transform(FFTW_FORWARD);
            if (getCopy) getFourierCopy(V);
            else         getFourierAlias(V);
        }

    /** Compute the Fourier transform.
        The data is taken from the matrix with which the object was
        created. */
    void FourierTransform();

    /** Enforce Hermitian symmetry.
        If the Fourier transform risks of losing Hermitian symmetry,
        use this function to renforce it. */
    void enforceHermitianSymmetry();

    /** Compute the inverse Fourier transform.
        The result is stored in the same real data that was passed for
        the forward transform. The Fourier coefficients are taken from
        the internal Fourier coefficients */
    void inverseFourierTransform();

    /** Compute the inverse Fourier transform.
        New data is provided for the Fourier coefficients and the output
        can be any matrix1D, 2D or 3D. It is important that the output
        matrix is already resized to the right size before entering
        in this function. */
    template <typename T, typename T1>
        void inverseFourierTransform(const T& V, T1& v)
        {
            setReal(v);
            setFourier(V);
            Transform(FFTW_BACKWARD);
        }

    /** Get Fourier coefficients. */
    template <typename T>
        void getFourierAlias(T& V) {V.alias(fFourier); return;}

    /** Get Fourier coefficients. */
    MultidimArray< Complex>& getFourierReference() {return fFourier;}

    /** Get Fourier coefficients. */
    template <typename T>
        void getFourierCopy(T& V) {
            V.reshape(fFourier);
            memcpy(MULTIDIM_ARRAY(V),MULTIDIM_ARRAY(fFourier),
                MULTIDIM_SIZE(fFourier)*2*sizeof(RFLOAT));
        }

    /** Return a complete Fourier transform (two halves).
    */
    template <typename T>
        void getCompleteFourier(T& V) {
            V.reshape(*fReal);
            int ndim=3;
            if (ZSIZE(*fReal)==1)
            {
                ndim=2;
                if (YSIZE(*fReal)==1)
                    ndim=1;
            }
            switch (ndim)
            {
                case 1:
                    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(V)
                        if (i<XSIZE(fFourier))
                            DIRECT_A1D_ELEM(V,i)=DIRECT_A1D_ELEM(fFourier,i);
                        else
                            DIRECT_A1D_ELEM(V,i)=
                                conj(DIRECT_A1D_ELEM(fFourier,
                                    XSIZE(*fReal)-i));
                    break;
                case 2:
                    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(V)
                        if (j<XSIZE(fFourier))
                            DIRECT_A2D_ELEM(V,i,j)=
                                DIRECT_A2D_ELEM(fFourier,i,j);
                        else
                            DIRECT_A2D_ELEM(V,i,j)=
                                conj(DIRECT_A2D_ELEM(fFourier,
                                    (YSIZE(*fReal)-i)%YSIZE(*fReal),
                                     XSIZE(*fReal)-j));
                    break;
                case 3:
                    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(V)
                        if (j<XSIZE(fFourier))
                            DIRECT_A3D_ELEM(V,k,i,j)=
                                DIRECT_A3D_ELEM(fFourier,k,i,j);
                        else
                            DIRECT_A3D_ELEM(V,k,i,j)=
                                conj(DIRECT_A3D_ELEM(fFourier,
                                    (ZSIZE(*fReal)-k)%ZSIZE(*fReal),
                                    (YSIZE(*fReal)-i)%YSIZE(*fReal),
                                     XSIZE(*fReal)-j));
                    break;
            }
        }

    /** Set one half of the FT in fFourier from the input complete Fourier transform (two halves).
        The fReal and fFourier already should have the right sizes
    */
    template <typename T>
        void setFromCompleteFourier(T& V) {
        int ndim=3;
        if (ZSIZE(*fReal)==1)
        {
            ndim=2;
            if (YSIZE(*fReal)==1)
                ndim=1;
        }
        switch (ndim)
        {
        case 1:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fFourier)
                DIRECT_A1D_ELEM(fFourier,i)=DIRECT_A1D_ELEM(V,i);
            break;
        case 2:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(fFourier)
                DIRECT_A2D_ELEM(fFourier,i,j) = DIRECT_A2D_ELEM(V,i,j);
            break;
        case 3:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(fFourier)
                DIRECT_A3D_ELEM(fFourier,k,i,j) = DIRECT_A3D_ELEM(V,k,i,j);
            break;
        }
    }

// Internal methods
public:
    /* Pointer to the array of RFLOATs with which the plan was computed */
    RFLOAT * dataPtr;

    /* Pointer to the array of complex<RFLOAT> with which the plan was computed */
    Complex * complexDataPtr;

    /* Initialise all pointers to NULL */
    void init();

    /** Clear object */
    void clear();

    /** This calls fftw_cleanup.
    */
    void cleanup();

    /** Destroy both forward and backward fftw planes (mutex locked */
    void destroyPlans();

    /** Computes the transform, specified in Init() function
        If normalization=true the forward transform is normalized
        (no normalization is made in the inverse transform)
        If normalize=false no normalization is performed and therefore
        the image is scaled by the number of pixels.
    */
    void Transform(int sign);

    /** Get the Multidimarray that is being used as input. */
    const MultidimArray<RFLOAT> &getReal() const;
    const MultidimArray<Complex > &getComplex() const;

    /** Set a Multidimarray for input.
        The data of img will be the one of fReal. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(MultidimArray<RFLOAT> &img);

    /** Set a Multidimarray for input.
        The data of img will be the one of fComplex. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(MultidimArray<Complex > &img);

    /** Set a Multidimarray for the Fourier transform.
        The values of the input array are copied in the internal array.
        It is assumed that the container for the real image as well as
        the one for the Fourier array are already resized.
        No plan is updated. */
    void setFourier(const MultidimArray<Complex> &imgFourier);
};

#endif
