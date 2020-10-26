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

#ifndef TENSOR_3X3_H
#define TENSOR_3X3_H

#include <cmath>
#include <vector>
#include <src/error.h>
#include <src/jaz/single_particle/index_sort.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/gravis/t3Matrix.h>

extern "C"
{
        #include <src/jaz/single_particle/d3x3/dsyevh3.h>
        #include <src/jaz/single_particle/d3x3/dsyevc3.h>
}


/* Symmetric 3x3 matrix to be used as e.g. a structure tensor of a 3D volume */

template <class T>
class Tensor3x3
{
    public:

        enum Component {XX = 0, XY = 1, XZ = 2, YY = 3, YZ = 4, ZZ = 5};

        Tensor3x3(){}
        Tensor3x3(T t) : xx(t), xy(t), xz(t), yy(t), yz(t), zz(t) {}
        Tensor3x3(T xx, T xy, T xz, T yy, T yz, T zz) : xx(xx), xy(xy), xz(xz), yy(yy), yz(yz), zz(zz) {}
				
		template <class T2>
		Tensor3x3(const gravis::t3Matrix<T2>& m) : 
			xx(m(0,0)), xy(m(0,1)), xz(m(0,2)), yy(m(1,1)), yz(m(1,2)), zz(m(2,2)) {}
		
		
			T xx, xy, xz, yy, yz, zz;

			
        static Tensor3x3<T> autoDyadicProduct(const gravis::t3Vector<T>& v);
        static Tensor3x3<T> dyadicProduct(const gravis::t3Vector<T>& v0, const gravis::t3Vector<T>& v1);
		
		gravis::t3Matrix<T> toMatrix() const;
		
        void diagonalize(gravis::t3Vector<T>& eigenvalues, gravis::t3Matrix<T>& eigenvectors) const;

        T& operator[] (int idx)
        {
            switch (idx)
            {
                case XX: return xx;
                case XY: return xy;
                case XZ: return xz;
                case YY: return yy;
                case YZ: return yz;
                case ZZ: return zz;
            }
            REPORT_ERROR("Tensor3x3 operator []: invalid index");
        }

        const T& operator[] (int idx) const
        {
            switch (idx)
            {
                case XX: return xx;
                case XY: return xy;
                case XZ: return xz;
                case YY: return yy;
                case YZ: return yz;
                case ZZ: return zz;
            }
            REPORT_ERROR("Tensor3x3 operator []: invalid index");
        }

        Tensor3x3<T>& operator += (const Tensor3x3<T>& arg)
        {
            xx += arg.xx;
            xy += arg.xy;
            xz += arg.xz;
            yy += arg.yy;
            yz += arg.yz;
            zz += arg.zz;

            return *this;
        }

        Tensor3x3& operator -= (const Tensor3x3& arg)
        {
            xx -= arg.xx;
            xy -= arg.xy;
            xz -= arg.xz;
            yy -= arg.yy;
            yz -= arg.yz;
            zz -= arg.zz;

            return *this;
        }

        Tensor3x3& operator *= (T arg)
        {
            xx *= arg;
            xy *= arg;
            xz *= arg;
            yy *= arg;
            yz *= arg;
            zz *= arg;

            return *this;
		}

		Tensor3x3 operator * (float arg) const
		{
			return Tensor3x3(xx * arg, xy * arg, xz * arg, yy * arg, yz * arg, zz * arg);
		}

		Tensor3x3 operator * (double arg) const
		{
			return Tensor3x3(xx * arg, xy * arg, xz * arg, yy * arg, yz * arg, zz * arg);
		}

        Tensor3x3& operator /= (T arg)
        {
            xx /= arg;
            xy /= arg;
            xz /= arg;
            yy /= arg;
            yz /= arg;
            zz /= arg;

            return *this;
        }

        Tensor3x3 operator / (T arg) const
        {
            return Tensor3x3(xx / arg, xy / arg, xz / arg, yy / arg, yz / arg, zz / arg);
        }
};

    template <class T>
    Tensor3x3<T> Tensor3x3<T>::autoDyadicProduct(const gravis::t3Vector<T>& v)
    {
        Tensor3x3<T> out;

        out.xx = v.x * v.x;
        out.xy = v.x * v.y;
        out.xz = v.x * v.z;
        out.yy = v.y * v.y;
        out.yz = v.y * v.z;
        out.zz = v.z * v.z;

        return out;
    }

    template <class T>
    Tensor3x3<T> Tensor3x3<T>::dyadicProduct(const gravis::t3Vector<T>& v0, const gravis::t3Vector<T>& v1)
    {
        Tensor3x3<T> out;

        out.xx = v0.x * v1.x;
        out.xy = v0.x * v1.y;
        out.xz = v0.x * v1.z;
        out.yy = v0.y * v1.y;
        out.yz = v0.y * v1.z;
        out.zz = v0.z * v1.z;

        return out;
    }
	
	template <class T>
	gravis::t3Matrix<T> Tensor3x3<T>::toMatrix() const
	{
		return gravis::t3Matrix<T>(
			xx, xy, xz,
			xy, yy, yz,
			xz, yz, zz);
	}

    template <class T>
    void Tensor3x3<T>::diagonalize(gravis::t3Vector<T>& eigenvalues, gravis::t3Matrix<T>& eigenvectors) const
    {
        double A[3][3];

        A[0][0] = (double)xx;

        A[0][1] = (double)xy;
        A[0][2] = (double)xz;

        A[1][1] = (double)yy;
        A[1][2] = (double)yz;
        A[2][2] = (double)zz;

        double Q[3][3];
		std::vector<double> w(3);

        dsyevh3(A, Q, &w[0]);

        std::vector<int> inds = IndexSort<double>::sortIndices(w);

        for (int i = 0; i < 3; i++)
        {
			eigenvalues[2-i] = (T) w[inds[i]];

            for (int j = 0; j < 3; j++)
            {
				eigenvectors(j,2-i) = (T) Q[j][inds[i]];
            }
        }

    }

    template <class T>
    inline
    Tensor3x3<T> operator + (const Tensor3x3<T>& v1, const Tensor3x3<T>& v2)
    {
        return Tensor3x3<T>(
            v1.xx + v2.xx, v1.xy + v2.xy, v1.xz + v2.xz,
            v1.yy + v2.yy, v1.yz + v2.yz, v1.zz + v2.zz);
    }

    template <class T>
    inline
    Tensor3x3<T> operator - (const Tensor3x3<T>& v1)
    {
        return Tensor3x3<T>(-v1.xx, -v1.xy, -v1.xz, -v1.yy, -v1.yz, -v1.zz);
    }

    template <class T>
    inline
    Tensor3x3<T> operator - (const Tensor3x3<T>& v1, const Tensor3x3<T>& v2)
    {
        return Tensor3x3<T>(
            v1.xx - v2.xx, v1.xy - v2.xy, v1.xz - v2.xz,
            v1.yy - v2.yy, v1.yz - v2.yz, v1.zz - v2.zz);
	}

	template <class T>
	inline
	Tensor3x3<T> operator * (float f, const Tensor3x3<T>& v)
	{
		return Tensor3x3<T>(f * v.xx, f * v.xy, f * v.xz, f * v.yy, f * v.yz, f * v.zz);
	}

	template <class T>
	inline
	Tensor3x3<T> operator * (const Tensor3x3<T>& v, float f)
	{
		return Tensor3x3<T>(v.xx * f, v.xy * f, v.xz * f, v.yy * f, v.yz * f, v.zz * f);
	}

	template <class T>
	inline
	Tensor3x3<T> operator * (double f, const Tensor3x3<T>& v)
	{
		return Tensor3x3<T>(f * v.xx, f * v.xy, f * v.xz, f * v.yy, f * v.yz, f * v.zz);
	}

	template <class T>
	inline
	Tensor3x3<T> operator * (const Tensor3x3<T>& v, double f)
	{
		return Tensor3x3<T>(v.xx * f, v.xy * f, v.xz * f, v.yy * f, v.yz * f, v.zz * f);
	}

    template <class T>
    inline
    Tensor3x3<T> operator / (const Tensor3x3<T>& v, T f)
    {
        return Tensor3x3<T>(v.xx / f, v.xy / f, v.xz / f, v.yy / f, v.yz / f, v.zz / f);
    }

    template <class T>
    inline
    std::ostream& operator<< (std::ostream& os, const Tensor3x3<T>& arg)
    {
        os << std::setprecision(17) << "["
                << std::setw(8) << arg.xx << ", "
                << std::setw(8) << arg.xy << ", "
                << std::setw(8) << arg.xz << ", "
                << std::setw(8) << arg.yy << ", "
                << std::setw(8) << arg.yz << ", "
                << std::setw(8) << arg.zz << "]";
        return os;
    }


#endif
