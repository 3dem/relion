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
/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef MATRIX2D_H_
#define MATRIX2D_H_

#include <string.h>
#include <iomanip>
#include "src/matrix1d.h"

/** @defgroup Matrices Matrix2D Matrices
 * @ingroup DataLibrary
 */
//@{
/** @name Matrices speed up macros */
//@{

/** Array access.
 *
 * This macro gives you access to the array (T)
 */
#define MATRIX2D_ARRAY(m) ((m).mdata)

/** For all elements in the array
 *
 * This macro is used to generate loops for the matrix in an easy way. It
 * defines internal indexes 'i' and 'j' which ranges the matrix using its
 * mathematical definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_MATRIX2D(m)
 * {
 *     std::cout << m(i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_MATRIX2D(m) \
    for (int i=0; i<(m).mdimy; i++) \
        for (int j=0; j<(m).mdimx; j++)

/** Access to a matrix element
 * v is the array, i and j define the element v_ij.
 *
  * @code
 * MAT_ELEM(m, 0, 0) = 1;
 * val = MAT_ELEM(m, 0, 0);
 * @endcode
 */
#define MAT_ELEM(m,i,j) ((m).mdata[(i)*(m).mdimx+(j)])

/** X dimension of the matrix
 */
#define MAT_XSIZE(m) ((m).mdimx)

/** Y dimension of the matrix
 */
#define MAT_YSIZE(m) ((m).mdimy)

// Forward declarations
template<typename T>
class Matrix1D;
template<typename T>
class Matrix2D;

template<typename T>
void ludcmp(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d);

template<typename T>
void lubksb(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b);

template<typename T>
void svdcmp(const Matrix2D< T >& a,
            Matrix2D< RFLOAT >& u,
            Matrix1D< RFLOAT >& w,
            Matrix2D< RFLOAT >& v);

void svbksb(Matrix2D< RFLOAT >& u,
            Matrix1D< RFLOAT >& w,
            Matrix2D< RFLOAT >& v,
            Matrix1D< RFLOAT >& b,
            Matrix1D< RFLOAT >& x);

template<typename T>
void solve(const Matrix2D<T>& A,
		   const Matrix1D<T>& b,
           Matrix1D< RFLOAT >& result,
           RFLOAT tolerance);

/** Matrix2D class */
template<typename T>
class Matrix2D
{
public:
    // The array itself
    T* mdata;

    // Destroy data
    bool destroyData;

    // Number of elements in X
    int mdimx;

    // Number of elements in Y
    int mdimy;

    // Total number of elements
    int mdim;

    /// @name Constructors
    /// @{
    /** Empty constructor
     */
    Matrix2D()
    {
        coreInit();
    }

    /** Dimension constructor
     */
    Matrix2D(int Ydim, int Xdim)
    {
        coreInit();
        resize(Ydim, Xdim);
    }

    /** Copy constructor
     */
    Matrix2D(const Matrix2D<T>& v)
    {
        coreInit();
        *this = v;
    }

    /** Destructor.
     */
    ~Matrix2D()
    {
        coreDeallocate();
    }

    /** Assignment.
     *
     * You can build as complex assignment expressions as you like. Multiple
     * assignment is allowed.
     *
     * @code
     * v1 = v2 + v3;
     * v1 = v2 = v3;
     * @endcode
     */
    Matrix2D<T>& operator=(const Matrix2D<T>& op1)
    {
        if (&op1 != this)
        {
            if (MAT_XSIZE(*this)!=MAT_XSIZE(op1) ||
                MAT_YSIZE(*this)!=MAT_YSIZE(op1))
                resize(op1);
            memcpy(mdata,op1.mdata,op1.mdim*sizeof(T));
        }

        return *this;
    }
    //@}

    /// @name Core memory operations for Matrix2D
    //@{
    /** Clear.
     */
    void clear()
    {
        coreDeallocate();
        coreInit();
    }

    /** Core init.
     * Initialize everything to 0
     */
    void coreInit()
    {
        mdimx=mdimy=mdim=0;
        mdata=NULL;
        destroyData=true;
    }

    /** Core allocate.
     */
    void coreAllocate(int _mdimy, int _mdimx)
    {
        if (_mdimy <= 0 ||_mdimx<=0)
        {
            clear();
            return;
        }

        mdimx=_mdimx;
        mdimy=_mdimy;
        mdim=_mdimx*_mdimy;
        mdata = new T [mdim];
        if (mdata == NULL)
            REPORT_ERROR("coreAllocate: No space left");
    }

    /** Core deallocate.
     * Free all mdata.
     */
    void coreDeallocate()
    {
        if (mdata != NULL && destroyData)
            delete[] mdata;
        mdata=NULL;
    }
    //@}

    /// @name Size and shape of Matrix2D
    //@{
    /** Resize to a given size
     */
    void resize(int Ydim, int Xdim)
    {

        if (Xdim == mdimx && Ydim == mdimy)
            return;

        if (Xdim <= 0 || Ydim <= 0)
        {
            clear();
            return;
        }

        T * new_mdata;
        size_t YXdim=Ydim*Xdim;

        try
        {
            new_mdata = new T [YXdim];
        }
        catch (std::bad_alloc &)
        {
            REPORT_ERROR("Allocate: No space left");
        }

        // Copy needed elements, fill with 0 if necessary
        for (int i = 0; i < Ydim; i++)
            for (int j = 0; j < Xdim; j++)
            {
                T val;
                if (i >= mdimy)
                    val = 0;
                else if (j >= mdimx)
                    val = 0;
                else
                    val = mdata[i*mdimx + j];
                new_mdata[i*Xdim+j] = val;
            }

        // deallocate old vector
        coreDeallocate();

        // assign *this vector to the newly created
        mdata = new_mdata;
        mdimx = Xdim;
        mdimy = Ydim;
        mdim = Xdim * Ydim;
    }

    /** Resize according to a pattern.
     *
     * This function resize the actual array to the same size and origin
     * as the input pattern. If the actual array is larger than the pattern
     * then the trailing values are lost, if it is smaller then 0's are
     * added at the end
     *
     * @code
     * v2.resize(v1);
     * // v2 has got now the same structure as v1
     * @endcode
     */
    template<typename T1>
    void resize(const Matrix2D<T1> &v)
    {
        if (mdimx != v.mdimx || mdimy != v.mdimy)
            resize(v.mdimy, v.mdimx);
    }

    /** Extract submatrix and assign to this object.
     */
    void submatrix(int i0, int j0, int iF, int jF)
    {
        if (i0 < 0 || j0 < 0 || iF >= MAT_YSIZE(*this) || jF >= MAT_XSIZE(*this))
            REPORT_ERROR("Submatrix indexes out of bounds");
        Matrix2D<T> result(iF - i0 + 1, jF - j0 + 1);

        FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
        MAT_ELEM(result, i, j) = MAT_ELEM(*this, i+i0, j+j0);

        *this = result;
    }

    /** Same shape.
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument
     */
    template <typename T1>
    bool sameShape(const Matrix2D<T1>& op) const
    {
        return ((mdimx == op.mdimx) && (mdimy == op.mdimy));
    }

    /** X dimension
     *
     * Returns X dimension
     */
    inline int Xdim() const
    {
        return mdimx;
    }

    /** Y dimension
     *
     * Returns Y dimension
     */
    inline int Ydim() const
    {
        return mdimy;
    }
    //@}

    /// @name Initialization of Matrix2D values
    //@{
    /** Same value in all components.
     *
     * The constant must be of a type compatible with the array type, ie,
     * you cannot  assign a RFLOAT to an integer array without a casting.
     * It is not an error if the array is empty, then nothing is done.
     *
     * @code
     * v.initConstant(3.14);
     * @endcode
     */
    void initConstant(T val)
    {
        for (int j = 0; j < mdim; j++)
            mdata[j] = val;
    }

    /** Initialize to zeros with current size.
     *
     * All values are set to 0. The current size and origin are kept. It is not
     * an error if the array is empty, then nothing is done.
     *
     * @code
     * v.initZeros();
     * @endcode
     */
    void initZeros()
    {
        memset(mdata,0,mdimx*mdimy*sizeof(T));
    }

    /** Initialize to zeros with a given size.
     */
    void initZeros(int Ydim, int Xdim)
    {
        if (mdimx!=Xdim || mdimy!=Ydim)
            resize(Ydim, Xdim);
        memset(mdata,0,mdimx*mdimy*sizeof(T));
    }

    /** Initialize to zeros following a pattern.
      *
      * All values are set to 0, and the origin and size of the pattern are
      * adopted.
      *
      * @code
      * v2.initZeros(v1);
      * @endcode
      */
    template <typename T1>
    void initZeros(const Matrix2D<T1>& op)
    {
        if (mdimx!=op.mdimx || mdimy!=op.mdimy)
            resize(op);
        memset(mdata,0,mdimx*mdimy*sizeof(T));
    }

    /** 2D Identity matrix of current size
     *
     * If actually the matrix is not squared then an identity matrix is
     * generated of size (Xdim x Xdim).
     *
     * @code
     * m.initIdentity();
     * @endcode
     */
    void initIdentity()
    {
        initIdentity(MAT_XSIZE(*this));
    }

    /** 2D Identity matrix of a given size
     *
     * A (dim x dim) identity matrix is generated.
     *
     * @code
     * m.initIdentity(3);
     * @endcode
     */
    void initIdentity(int dim)
    {
        initZeros(dim, dim);
        for (int i = 0; i < dim; i++)
            MAT_ELEM(*this,i,i) = 1;
    }
    //@}

    /// @name Operators for Matrix2D
    //@{

    /** Matrix element access
     */
	T& operator()(int i, int j)
    {
        return MAT_ELEM((*this),i,j);
    }
	
	// for constant matrices (the compiler will pick the right version)
	const T& operator()(int i, int j) const
    {
        return MAT_ELEM((*this),i,j);
    }
	
    /** Parenthesis operator for phyton
    */
    void setVal(T val,int y, int x)
    {
        MAT_ELEM((*this),y,x)=val;
    }
    /** Parenthesis operator for phyton
    */
    T getVal( int y, int x) const
    {
        return MAT_ELEM((*this),y,x);
    }

    /** v3 = v1 * k.
     */
    Matrix2D<T> operator*(T op1) const
    {
        Matrix2D<T> tmp(*this);
        for (int i=0; i < mdim; i++)
            tmp.mdata[i] = mdata[i] * op1;
        return tmp;
    }

    /** v3 = v1 / k.
     */
    Matrix2D<T> operator/(T op1) const
    {
        Matrix2D<T> tmp(*this);
        for (int i=0; i < mdim; i++)
            tmp.mdata[i] = mdata[i] / op1;
        return tmp;
    }

    /** v3 = k * v2.
     */
    friend Matrix2D<T> operator*(T op1, const Matrix2D<T>& op2)
    {
        Matrix2D<T> tmp(op2);
        for (int i=0; i < op2.mdim; i++)
            tmp.mdata[i] = op1 * op2.mdata[i];
        return tmp;
    }

    /** v3 *= k.
      */
    void operator*=(T op1)
    {
        for (int i=0; i < mdim; i++)
            mdata[i] *= op1;
    }

    /** v3 /= k.
      */
    void operator/=(T op1)
    {
        for (int i=0; i < mdim; i++)
            mdata[i] /= op1;
    }

    /** Matrix by vector multiplication
     *
     * @code
     * v2 = A*v1;
     * @endcode
     */
    Matrix1D<T> operator*(const Matrix1D<T>& op1) const
    {
        Matrix1D<T> result;

        if (mdimx != op1.size())
        {
        	std::cerr << " mdimx= " << mdimx << " opp1.size()= " << op1.size() << std::endl;
        	REPORT_ERROR("Not compatible sizes in matrix by vector");
        }

        if (!op1.isCol())
            REPORT_ERROR("Vector is not a column");

        result.initZeros(mdimy);

        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < op1.size(); j++)
                result(i) += (*this)(i, j) * op1(j);

        result.setCol();
        return result;
    }

    /** Matrix by Matrix multiplication
     *
     * @code
     * C = A*B;
     * @endcode
     */
    Matrix2D<T> operator*(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> result;
        if (mdimx != op1.mdimy)
            REPORT_ERROR("Not compatible sizes in matrix multiplication");

        result.initZeros(mdimy, op1.mdimx);
        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < op1.mdimx; j++)
                for (int k = 0; k < mdimx; k++)
                    result(i, j) += (*this)(i, k) * op1(k, j);
        return result;
    }

    /** Matrix summation
     *
     * @code
     * C = A + B;
     * @endcode
     */
    Matrix2D<T> operator+(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> result;
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR("operator+: Not same sizes in matrix summation");

        result.initZeros(mdimy, mdimx);
        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                result(i, j) = (*this)(i, j) + op1(i, j);

        return result;
    }

    /** Matrix summation
     *
     * @code
     * A += B;
     * @endcode
     */
    void operator+=(const Matrix2D<T>& op1) const
    {
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR("operator+=: Not same sizes in matrix summation");

        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                MAT_ELEM(*this,i, j) += MAT_ELEM(op1, i, j);
    }

    /** Matrix subtraction
     *
     * @code
     * C = A - B;
     * @endcode
     */
    Matrix2D<T> operator-(const Matrix2D<T>& op1) const
    {
        Matrix2D<T> result;
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR("operator-: Not same sizes in matrix summation");

        result.initZeros(mdimy, mdimx);
        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                result(i, j) = (*this)(i, j) - op1(i, j);

        return result;
    }

    /** Matrix substraction
     *
     * @code
     * A -= B;
     * @endcode
     */
    void operator-=(const Matrix2D<T>& op1) const
    {
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR("operator-=: Not same sizes in matrix summation");

        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                MAT_ELEM(*this,i, j) -= MAT_ELEM(op1, i, j);
    }

    /** Equality.
     *
     * Returns true if this object has got the same shape (origin and size)
     * than the argument and the same values (within accuracy).
     */
    bool equal(const Matrix2D<T>& op,
               RFLOAT accuracy = XMIPP_EQUAL_ACCURACY) const
    {
        if (!sameShape(op))
            return false;
        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                if (ABS( (*this)(i,j) - op(i,j) ) > accuracy)
                    return false;
        return true;
    }
    //@}

    /// @name Utilities for Matrix2D
    //@{
    /** Set very small values (ABS(val)< accuracy) equal to zero
      *
      */
    void setSmallValuesToZero(RFLOAT accuracy = XMIPP_EQUAL_ACCURACY)
    {
        for (int i = 0; i < mdimy; i++)
             for (int j = 0; j < mdimx; j++)
                 if (ABS( (*this)(i,j) ) < accuracy)
                	 (*this)(i,j) = 0.;
    }

    /// @name Utilities for Matrix2D
    //@{
    /** Maximum of the values in the array.
      *
      * The returned value is of the same type as the type of the array.
      */
    T computeMax() const
    {
        if (mdim <= 0)
            return static_cast< T >(0);

        T maxval = mdata[0];
        for (int n = 0; n < mdim; n++)
            if (mdata[n] > maxval)
                maxval = mdata[n];
        return maxval;
    }

    /** Minimum of the values in the array.
       *
       * The returned value is of the same type as the type of the array.
       */
    T computeMin() const
    {
        if (mdim <= 0)
            return static_cast< T >(0);

        T minval = mdata[0];
        for (int n = 0; n < mdim; n++)
            if (mdata[n] < minval)
                minval = mdata[n];
        return minval;
    }

    /** Produce a 2D array suitable for working with Numerical Recipes
    *
    * This function must be used only as a preparation for routines which need
    * that the first physical index is 1 and not 0 as it usually is in C. New
    * memory is needed to hold the new RFLOAT pointer array.
    */
    T** adaptForNumericalRecipes() const
    {
        T** m = NULL;
        ask_Tmatrix(m, 1, mdimy, 1, mdimx);

        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                m[i+1][j+1] = mdata[i*mdimx + j];

        return m;
    }

    /** Produce a 1D pointer suitable for working with Numerical Recipes (2)
     *
     * This function meets the same goal as the one before, however this one
     * work with 2D arrays as a single pointer. The first element of the array
     * is pointed by result[1*Xdim+1], and in general result[i*Xdim+j]
     */
    T* adaptForNumericalRecipes2() const
    {
        return mdata - 1 - mdimx;
    }

    /** Load 2D array from numerical recipes result.
     */
    void loadFromNumericalRecipes(T** m, int Ydim, int Xdim)
    {
        if (mdimx!=Xdim || mdimy!=Ydim)
            resize(Ydim, Xdim);

        for (int i = 1; i <= Ydim; i++)
            for (int j = 1; j <= Xdim; j++)
                (*this)(i - 1, j - 1) = m[i][j];
    }

    /** Kill a 2D array produced for numerical recipes
     *
     * The allocated memory is freed.
     */
    void killAdaptationForNumericalRecipes(T** m) const
    {
        free_Tmatrix(m, 1, mdimy, 1, mdimx);
    }

    /** Kill a 2D array produced for numerical recipes, 2.
     *
     * Nothing needs to be done.
     */
    void killAdaptationForNumericalRecipes2(T** m) const
        {}

    /** Write this matrix to file
      */
    void write(const FileName &fn) const
    {
        std::ofstream fhOut;
        fhOut.open(fn.c_str());
        if (!fhOut)
            REPORT_ERROR((std::string)"write: Cannot open "+fn+" for output");
        fhOut << *this;
        fhOut.close();
    }

    /** Show matrix
      */
    friend std::ostream& operator<<(std::ostream& ostrm, const Matrix2D<T>& v)
    {
        if (v.Xdim() == 0 || v.Ydim() == 0)
            ostrm << "NULL matrix\n";
        else
        {
            ostrm << std::endl;
            RFLOAT max_val = v.computeMax();
            int prec = bestPrecision(max_val, 10);

            for (int i = 0; i < v.Ydim(); i++)
            {
                for (int j = 0; j < v.Xdim(); j++)
                {
                    ostrm << std::setw(13) << floatToString((RFLOAT) v(i, j), 10, prec) << ' ';
                }
                ostrm << std::endl;
            }
        }

        return ostrm;
    }

    /** Makes a matrix from a vector
     *
     * The origin of the matrix is set such that it has one of the index origins
     * (X or Y) to the same value as the vector, and the other set to 0
     * according to the shape.
     *
     * @code
     * Matrix2D< RFLOAT > m = fromVector(v);
     * @endcode
     */
    void fromVector(const Matrix1D<T>& op1)
    {
        // Null vector => Null matrix
        if (op1.size() == 0)
        {
            clear();
            return;
        }

        // Look at shape and copy values
        if (op1.isRow())
        {
            if (mdimy!=1 || mdimx!=VEC_XSIZE(op1))
                resize(1, VEC_XSIZE(op1));

            for (int j = 0; j < VEC_XSIZE(op1); j++)
                MAT_ELEM(*this,0, j) = VEC_ELEM(op1,j);
        }
        else
        {
            if (mdimy!=1 || mdimx!=VEC_XSIZE(op1))
                resize(VEC_XSIZE(op1), 1);

            for (int i = 0; i < VEC_XSIZE(op1); i++)
                MAT_ELEM(*this,i, 0) = VEC_ELEM(op1,i);
        }
    }

    /** Makes a vector from a matrix
     *
     * An exception is thrown if the matrix is not a single row or a single
     * column. The origin of the vector is set according to the one of the
     * matrix.
     *
     * @code
     * Matrix1D< RFLOAT > v;
     * m.toVector(v);
     * @endcode
     */
    void toVector(Matrix1D<T>& op1) const
    {
        // Null matrix => Null vector
        if (mdimx == 0 || mdimy == 0)
        {
            op1.clear();
            return;
        }

        // If matrix is not a vector, produce an error
        if (!(mdimx == 1 || mdimy == 1))
            REPORT_ERROR("toVector: Matrix cannot be converted to vector");

        // Look at shape and copy values
        if (mdimy == 1)
        {
            // Row vector
            if (VEC_XSIZE(op1)!=mdimx)
                op1.resize(mdimx);

            for (int j = 0; j < mdimx; j++)
                VEC_ELEM(op1,j) = MAT_ELEM(*this,0, j);

            op1.setRow();
        }
        else
        {
            // Column vector
            if (VEC_XSIZE(op1)!=mdimy)
                op1.resize(mdimy);

            for (int i = 0; i < mdimy; i++)
                VEC_ELEM(op1,i) = MAT_ELEM(*this,i, 0);

            op1.setCol();
        }
    }

    /**Copy matrix to stl::vector
     */
    void copyToVector(std::vector<T> &v)
    {
        v.assign(mdata, mdata+mdim);
    }
    /**Copy stl::vector to matrix
      */
    void copyFromVector(std::vector<T> &v,int Xdim, int Ydim)
    {
        if (mdimx!=Xdim || mdimy!=Ydim)
            resize(Ydim, Xdim);
        copy( v.begin(), v.begin()+v.size(), mdata);
    }

    /** Get row
     *
     * This function returns a row vector corresponding to the choosen
     * row inside the nth 2D matrix, the numbering of the rows is also
     * logical not physical.
     *
     * @code
     * std::vector< RFLOAT > v;
     * m.getRow(-2, v);
     * @endcode
     */
    void getRow(int i, Matrix1D<T>& v) const
    {
        if (mdimx == 0 || mdimy == 0)
        {
            v.clear();
            return;
        }

        if (i < 0 || i >= mdimy)
            REPORT_ERROR("getRow: Matrix subscript (i) greater than matrix dimension");

        if (VEC_XSIZE(v)!=mdimx)
            v.resize(mdimx);
        for (int j = 0; j < mdimx; j++)
            VEC_ELEM(v,j) = MAT_ELEM(*this,i, j);

        v.setRow();
    }

    /** Get Column
     *
     * This function returns a column vector corresponding to the
     * choosen column.
     *
     * @code
     * std::vector< RFLOAT > v;
     * m.getCol(-1, v);
     * @endcode
     */
    void getCol(int j, Matrix1D<T>& v) const
    {
        if (mdimx == 0 || mdimy == 0)
        {
            v.clear();
            return;
        }

        if (j < 0 || j >= mdimx)
            REPORT_ERROR("getCol: Matrix subscript (j) greater than matrix dimension");

        if (VEC_XSIZE(v)!=mdimy)
            v.resize(mdimy);
        for (int i = 0; i < mdimy; i++)
            VEC_ELEM(v,i) = MAT_ELEM(*this,i, j);

        v.setCol();
    }

    /** Set Row
     *
     * This function sets a row vector corresponding to the choosen row in the 2D Matrix
     *
     * @code
     * m.setRow(-2, m.row(1)); // Copies row 1 in row -2
     * @endcode
     */
    void setRow(int i, const Matrix1D<T>& v)
    {
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR("setRow: Target matrix is empty");

        if (i < 0 || i >= mdimy)
            REPORT_ERROR("setRow: Matrix subscript (i) out of range");

        if (VEC_XSIZE(v) != mdimx)
            REPORT_ERROR("setRow: Vector dimension different from matrix one");

        if (!v.isRow())
            REPORT_ERROR("setRow: Not a row vector in assignment");

        for (int j = 0; j < mdimx; j++)
            MAT_ELEM(*this,i, j) = VEC_ELEM(v,j);
    }

    /** Set Column
     *
     * This function sets a column vector corresponding to the choosen column
     * inside matrix.
     *
     * @code
     * m.setCol(0, (m.row(1)).transpose()); // Copies row 1 in column 0
     * @endcode
     */
    void setCol(int j, const Matrix1D<T>& v)
    {
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR("setCol: Target matrix is empty");

        if (j < 0 || j>= mdimx)
            REPORT_ERROR("setCol: Matrix subscript (j) out of range");

        if (VEC_XSIZE(v) != mdimy)
            REPORT_ERROR("setCol: Vector dimension different from matrix one");

        if (!v.isCol())
            REPORT_ERROR("setCol: Not a column vector in assignment");

        for (int i = 0; i < mdimy; i++)
            MAT_ELEM(*this,i, j) = VEC_ELEM(v,i);
    }

    /** Determinant of a matrix
     *
     * An exception is thrown if the matrix is not squared or it is empty.
     *
     * @code
     * RFLOAT det = m.det();
     * @endcode
     */
    T det() const
    {
        // (see Numerical Recipes, Chapter 2 Section 5)
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR("determinant: Matrix is empty");

        if (mdimx != mdimy)
            REPORT_ERROR("determinant: Matrix is not squared");

        for (int i = 0; i < mdimy; i++)
        {
            bool all_zeros = true;
            for (int j = 0; j < mdimx; j++)
                if (ABS(MAT_ELEM((*this),i, j)) > XMIPP_EQUAL_ACCURACY)
                {
                    all_zeros = false;
                    break;
                }

            if (all_zeros)
                return 0;
        }

        // Perform decomposition
        Matrix1D< int > indx;
        T d;
        Matrix2D<T> LU;
        ludcmp(*this, LU, indx, d);

        // Calculate determinant
        for (int i = 0; i < mdimx; i++)
            d *= (T) MAT_ELEM(LU,i , i);

        return d;
    }

    /** Algebraic transpose of a Matrix
     *
     * You can use the transpose in as complex expressions as you like. The
     * origin of the vector is not changed.
     *
     * @code
     * v2 = v1.transpose();
     * @endcode
     */
    Matrix2D<T> transpose() const
    {
        Matrix2D<T> result(mdimx, mdimy);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
        MAT_ELEM(result,i,j) = MAT_ELEM((*this),j,i);
        return result;
    }

    /** Inverse of a matrix
     *
     * The matrix is inverted using a SVD decomposition. In fact the
     * pseudoinverse is returned.
     *
     * @code
     * Matrix2D< RFLOAT > m1_inv;
     * m1.inv(m1_inv);
     * @endcode
     */
    void inv(Matrix2D<T>& result) const
    {

        if (mdimx == 0 || mdimy == 0)
        {
        	REPORT_ERROR("Inverse: Matrix is empty");
        }
        // Initialise output
        result.initZeros(mdimx, mdimy);

        if (mdimx == 3 && mdimy == 3)
        {
        	MAT_ELEM(result, 0, 0) =   MAT_ELEM((*this), 2, 2)*MAT_ELEM((*this), 1, 1)-MAT_ELEM((*this), 2, 1)*MAT_ELEM((*this), 1, 2);
        	MAT_ELEM(result, 0, 1) = -(MAT_ELEM((*this), 2, 2)*MAT_ELEM((*this), 0, 1)-MAT_ELEM((*this), 2, 1)*MAT_ELEM((*this), 0, 2));
        	MAT_ELEM(result, 0, 2) =   MAT_ELEM((*this), 1, 2)*MAT_ELEM((*this), 0, 1)-MAT_ELEM((*this), 1, 1)*MAT_ELEM((*this), 0, 2);
        	MAT_ELEM(result, 1, 0) = -(MAT_ELEM((*this), 2, 2)*MAT_ELEM((*this), 1, 0)-MAT_ELEM((*this), 2, 0)*MAT_ELEM((*this), 1, 2));
        	MAT_ELEM(result, 1, 1) =   MAT_ELEM((*this), 2, 2)*MAT_ELEM((*this), 0, 0)-MAT_ELEM((*this), 2, 0)*MAT_ELEM((*this), 0, 2);
        	MAT_ELEM(result, 1, 2) = -(MAT_ELEM((*this), 1, 2)*MAT_ELEM((*this), 0, 0)-MAT_ELEM((*this), 1, 0)*MAT_ELEM((*this), 0, 2));
        	MAT_ELEM(result, 2, 0) =   MAT_ELEM((*this), 2, 1)*MAT_ELEM((*this), 1, 0)-MAT_ELEM((*this), 2, 0)*MAT_ELEM((*this), 1, 1);
        	MAT_ELEM(result, 2, 1) = -(MAT_ELEM((*this), 2, 1)*MAT_ELEM((*this), 0, 0)-MAT_ELEM((*this), 2, 0)*MAT_ELEM((*this), 0, 1));
        	MAT_ELEM(result, 2, 2) =   MAT_ELEM((*this), 1, 1)*MAT_ELEM((*this), 0, 0)-MAT_ELEM((*this), 1, 0)*MAT_ELEM((*this), 0, 1);
        	RFLOAT tmp = MAT_ELEM((*this), 0, 0) * MAT_ELEM(result, 0, 0) +
        			     MAT_ELEM((*this), 1, 0) * MAT_ELEM(result, 0, 1) +
        			     MAT_ELEM((*this), 2, 0) * MAT_ELEM(result, 0, 2);
        	result /= tmp;
        }
        else if (mdimx == 2 && mdimy == 2)
        {
        	MAT_ELEM(result, 0, 0) = MAT_ELEM((*this), 1, 1);
        	MAT_ELEM(result, 0, 1) = -MAT_ELEM((*this), 0, 1);
        	MAT_ELEM(result, 1, 0) = -MAT_ELEM((*this), 1, 0);
        	MAT_ELEM(result, 1, 1) =  MAT_ELEM((*this), 0, 0);
        	RFLOAT tmp = MAT_ELEM((*this), 0, 0) * MAT_ELEM((*this), 1, 1) -
					     MAT_ELEM((*this), 0, 1) * MAT_ELEM((*this), 1, 0);
        	result /= tmp;
        }
        else
        {

			// Perform SVD decomposition
			Matrix2D< RFLOAT > u, v;
			Matrix1D< RFLOAT > w;
			svdcmp(*this, u, w, v); // *this = U * W * V^t

			RFLOAT tol = computeMax() * XMIPP_MAX(mdimx, mdimy) * 1e-14;

			// Compute W^-1
			bool invertible = false;
			FOR_ALL_ELEMENTS_IN_MATRIX1D(w)
			{
				if (ABS(VEC_ELEM(w,i)) > tol)
				{
					VEC_ELEM(w,i) = 1.0 / VEC_ELEM(w,i);
					invertible = true;
				}
				else
					VEC_ELEM(w,i) = 0.0;
			}

			if (!invertible)
				return;

			// Compute V*W^-1
			FOR_ALL_ELEMENTS_IN_MATRIX2D(v)
			MAT_ELEM(v,i,j) *= VEC_ELEM(w,j);

			// Compute Inverse
			for (int i = 0; i < mdimx; i++)
				for (int j = 0; j < mdimy; j++)
					for (int k = 0; k < mdimx; k++)
						MAT_ELEM(result,i,j) += (T) MAT_ELEM(v,i,k) * MAT_ELEM(u,j,k);

       }

    }

    /** Inverse of a matrix
     */
    Matrix2D<T> inv() const
    {
        Matrix2D<T> result;
        inv(result);

        return result;
    }

    /** True if the matrix is identity
     *
     * @code
     * if (m.isIdentity())
     *     std::cout << "The matrix is identity\n";
     * @endcode
     */
    bool isIdentity() const
    {
        for (int i = 0; i < mdimy; i++)
            for (int j = 0; j < mdimx; j++)
                if (i != j)
                {
                    if (ABS(MAT_ELEM(*this,i,j)) > XMIPP_EQUAL_ACCURACY)
                        return false;
                }
                else
                {
                    if (ABS(MAT_ELEM(*this,i,j) - 1.) > XMIPP_EQUAL_ACCURACY )
                        return false;
                }
        return true;
    }
    //@}
};

// Implementation of the vector*matrix
// Documented in matrix1D.h
template<typename T>
Matrix1D<T> Matrix1D<T>::operator*(const Matrix2D<T>& M)
{
    Matrix1D<T> result;

    if (VEC_XSIZE(*this) != MAT_YSIZE(M))
        REPORT_ERROR("Not compatible sizes in matrix by vector");

    if (!isRow())
        REPORT_ERROR("Vector is not a row");

    result.initZeros(MAT_XSIZE(M));
    for (int j = 0; j < MAT_XSIZE(M); j++)
        for (int i = 0; i < MAT_YSIZE(M); i++)
            VEC_ELEM(result,j) += VEC_ELEM(*this,i) * MAT_ELEM(M,i, j);

    result.setRow();
    return result;
}

/**@name Matrix Related functions
 * These functions are not methods of Matrix2D
 */
//@{
/** LU Decomposition
 */
template<typename T>
void ludcmp(const Matrix2D<T>& A, Matrix2D<T>& LU, Matrix1D< int >& indx, T& d)
{
    LU = A;
    if (VEC_XSIZE(indx)!=A.mdimx)
        indx.resize(A.mdimx);
    ludcmp(LU.adaptForNumericalRecipes2(), A.mdimx,
           indx.adaptForNumericalRecipes(), &d);
}

/** LU Backsubstitution
 */
template<typename T>
void lubksb(const Matrix2D<T>& LU, Matrix1D< int >& indx, Matrix1D<T>& b)
{
    lubksb(LU.adaptForNumericalRecipes2(), indx.size(),
           indx.adaptForNumericalRecipes(),
           b.adaptForNumericalRecipes());
}

/** SVD Backsubstitution
 */
void svbksb(Matrix2D< RFLOAT >& u,
            Matrix1D< RFLOAT >& w,
            Matrix2D< RFLOAT >& v,
            Matrix1D< RFLOAT >& b,
            Matrix1D< RFLOAT >& x);

/** SVD Decomposition (through numerical recipes)
 */
template<typename T>
void svdcmp(const Matrix2D< T >& a,
            Matrix2D< RFLOAT >& u,
            Matrix1D< RFLOAT >& w,
            Matrix2D< RFLOAT >& v)
{
    // svdcmp only works with RFLOAT
    typeCast(a, u);

    // Set size of matrices
    w.initZeros(u.mdimx);
    v.initZeros(u.mdimx, u.mdimx);

    // Call to the numerical recipes routine
    svdcmp(u.mdata,
           u.mdimy, u.mdimx,
           w.vdata,
           v.mdata);
}

/** Solve system of linear equations (Ax=b) through SVD Decomposition (through numerical recipes)
 */
template<typename T>
void solve(const Matrix2D< T >& A, const Matrix1D< T >& b,
                  Matrix1D< RFLOAT >& result, RFLOAT tolerance)
{
    if (A.mdimx == 0)
        REPORT_ERROR("Solve: Matrix is empty");

    /*if (A.mdimx != A.mdimy)
        REPORT_ERROR("Solve: Matrix is not squared");*/

    if (A.mdimy != b.vdim)
        REPORT_ERROR("Solve: Different sizes of Matrix and Vector");

    /*if (b.isRow())
        REPORT_ERROR("Solve: Not correct vector shape");*/

    // First perform de single value decomposition
    // Xmipp interface that calls to svdcmp of numerical recipes
    Matrix2D< RFLOAT > u, v;
    Matrix1D< RFLOAT > w;
    svdcmp(A, u, w, v);

    // Here is checked if eigenvalues of the svd decomposition are acceptable
    // If a value is lower than tolerance, the it's zeroed, as this increases
    // the precision of the routine.
    FOR_ALL_ELEMENTS_IN_MATRIX1D(w)
    if (w(i) < tolerance)
        w(i) = 0;

    // Set size of matrices
    result.resize(b.vdim);

    // Xmipp interface that calls to svdksb of numerical recipes
    Matrix1D< RFLOAT > bd;
    typeCast(b, bd);
    svbksb(u, w, v, bd, result);
}

/** Solve system of linear equations (Ax=b), x and b being matrices through SVD Decomposition (through Gauss-Jordan numerical recipes)
 */
template<typename T>
void solve(const Matrix2D<T>& A, const Matrix2D<T>& b, Matrix2D<T>& result)
{
    if (A.mdimx == 0)
        REPORT_ERROR("Solve: Matrix is empty");

    if (A.mdimx != A.mdimy)
        REPORT_ERROR("Solve: Matrix is not squared");

    if (A.mdimy != b.mdimy)
        REPORT_ERROR("Solve: Different sizes of A and b");

    // Solve
    result = b;
    Matrix2D<T> Aux = A;
    gaussj(Aux.adaptForNumericalRecipes2(), Aux.mdimy,
           result.adaptForNumericalRecipes2(), b.mdimx);
}


/** Least-squares rigid transformation between two sets of 3D coordinates
 *
RFLOAT lsq_rigid_body_transformation(std::vector<Matrix1D<RFLOAT> > &set1, std::vector<Matrix1D<RFLOAT> > &set2,
		Matrix2D<RFLOAT> &Rot, Matrix1D<RFLOAT> &trans)
{
	Matrix2D<RFLOAT> A;
	Matrix1D<RFLOAT> avg1, avg2;

	if (set1.size() != set2.size())
		REPORT_ERROR("lsq_rigid_body_transformation ERROR: unequal set size");

	// Calculate average of set1 and set2
	avg1 = vectorR3(0., 0., 0.);
	avg2 = vectorR3(0., 0., 0.);
	for (int i = 0; i < set1.size(); i++)
	{
		if (set1[i].vdim != 3)
			REPORT_ERROR("lsq_rigid_body_transformation ERROR: not a 3-point set1");
		if (set2[i].vdim != 3)
			REPORT_ERROR("lsq_rigid_body_transformation ERROR: not a 3-point set2");
		avg1 += set1[i];
		avg2 += set2[i];
	}
	avg1 /= (RFLOAT)set1.size();
	avg2 /= (RFLOAT)set1.size();

	A.initZeros(3, 3);
	Rot.initZeros(4,4);
	for (int i = 0; i < set1.size(); i++)
	{
		// fill A
		A(0, 0) += (XX(set1[i]) - XX(avg1)) * (XX(set2[i]) - XX(avg2));
		A(0, 1) += (XX(set1[i]) - XX(avg1)) * (YY(set2[i]) - YY(avg2));
		A(0, 2) += (XX(set1[i]) - XX(avg1)) * (ZZ(set2[i]) - ZZ(avg2));
		A(1, 0) += (YY(set1[i]) - YY(avg1)) * (XX(set2[i]) - XX(avg2));
		A(1, 1) += (YY(set1[i]) - YY(avg1)) * (YY(set2[i]) - YY(avg2));
		A(1, 2) += (YY(set1[i]) - YY(avg1)) * (ZZ(set2[i]) - ZZ(avg2));
		A(2, 0) += (ZZ(set1[i]) - ZZ(avg1)) * (XX(set2[i]) - XX(avg2));
		A(2, 1) += (ZZ(set1[i]) - ZZ(avg1)) * (YY(set2[i]) - YY(avg2));
		A(2, 2) += (ZZ(set1[i]) - ZZ(avg1)) * (ZZ(set2[i]) - ZZ(avg2));
	}

	Matrix2D< RFLOAT > U, V;
	Matrix1D< RFLOAT > w;

	// TODO: check inverse, transpose etc etc!!!

	// Optimal rotation
	svdcmp(A, U, w, V);
	Rot = V.transpose() * U;

	// Optimal translation
	trans = avg1 - Rot * avg2;

	// return the squared difference term
	RFLOAT error = 0.;
	for (int i = 0; i < set1.size(); i++)
	{
		error += (Rot * set2[i] + trans - set1[i]).sum2();
	}

	return error;

}
*/

/** Conversion from one type to another.
 *
 * If we have an integer array and we need a RFLOAT one, we can use this
 * function. The conversion is done through a type casting of each element
 * If n >= 0, only the nth volumes will be converted, otherwise all NSIZE volumes
 */
template<typename T1, typename T2>
void typeCast(const Matrix2D<T1>& v1,  Matrix2D<T2>& v2)
{
    if (v1.mdim == 0)
    {
        v2.clear();
        return;
    }

    if (v1.mdimx!=v2.mdimx || v1.mdimy!=v2.mdimy)
        v2.resize(v1);
    for (unsigned long int n = 0; n < v1.mdim; n++)
        v2.mdata[n] = static_cast< T2 > (v1.mdata[n]);
}
//@}
//@}
#endif /* MATRIX2D_H_ */
