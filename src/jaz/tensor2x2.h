#ifndef TENSOR_2X2_H
#define TENSOR_2X2_H

#include <cmath>
#include <vector>
#include <src/error.h>
#include <src/jaz/index_sort.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t2Matrix.h>

extern "C"
{
        #include <src/jaz/d3x3/dsyev2.h>
}


/* Symmetric 2x2 matrix to be used as e.g. a structure tensor of a 2D image */

template <class T>
class Tensor2x2
{
    public:

        enum Component {XX = 0, XY = 1, YY = 2};

        Tensor2x2(){}
        Tensor2x2(T t) : xx(t), xy(t), yy(t) {}
        Tensor2x2(T xx, T xy, T yy) : xx(xx), xy(xy), yy(yy) {}

        T xx, xy, yy;

        static Tensor2x2<T> autoDyadicProduct(const gravis::t2Vector<T>& v);
        static Tensor2x2<T> dyadicProduct(const gravis::t2Vector<T>& v0, const gravis::t2Vector<T>& v1);

        void diagonalize(gravis::t2Vector<T>& eigenvalues, gravis::t2Matrix<T>& eigenvectors) const;
        gravis::t2Matrix<T> toMatrix() const;

        T& operator[] (int idx)
        {
            switch (idx)
            {
                case XX: return xx;
                case XY: return xy;
                case YY: return yy;
            }
            REPORT_ERROR("Tensor2x2 operator []: invalid index");
        }

        const T& operator[] (int idx) const
        {
            switch (idx)
            {
                case XX: return xx;
                case XY: return xy;
                case YY: return yy;
            }
            REPORT_ERROR("Tensor2x2 operator []: invalid index");
        }

        Tensor2x2<T>& operator += (const Tensor2x2<T>& arg)
        {
            xx += arg.xx;
            xy += arg.xy;
            yy += arg.yy;

            return *this;
        }

        Tensor2x2& operator -= (const Tensor2x2& arg)
        {
            xx -= arg.xx;
            xy -= arg.xy;
            yy -= arg.yy;

            return *this;
        }

        Tensor2x2& operator *= (T arg)
        {
            xx *= arg;
            xy *= arg;
            yy *= arg;

            return *this;
        }

        Tensor2x2 operator * (T arg) const
        {
            return Tensor2x2(xx * arg, xy * arg, yy * arg);
        }

        Tensor2x2& operator /= (T arg)
        {
            xx /= arg;
            xy /= arg;
            yy /= arg;

            return *this;
        }

        Tensor2x2 operator / (T arg) const
        {
            return Tensor2x2(xx / arg, xy / arg, yy / arg);
        }
};

    template <class T>
    Tensor2x2<T> Tensor2x2<T>::autoDyadicProduct(const gravis::t2Vector<T>& v)
    {
        Tensor2x2<T> out;

        out.xx = v.x * v.x;
        out.xy = v.x * v.y;
        out.yy = v.y * v.y;

        return out;
    }

    template <class T>
    Tensor2x2<T> Tensor2x2<T>::dyadicProduct(const gravis::t2Vector<T>& v0, const gravis::t2Vector<T>& v1)
    {
        Tensor2x2<T> out;

        out.xx = v0.x0 * v1.x1;
        out.xy = v0.x0 * v1.y1;
        out.xz = v0.x0 * v1.z1;
        out.yy = v0.y0 * v1.y1;
        out.yz = v0.y0 * v1.z1;
        out.zz = v0.z0 * v1.z1;

        return out;
    }

    template <class T>
    void Tensor2x2<T>::diagonalize(gravis::t2Vector<T>& eigenvalues, gravis::t2Matrix<T>& eigenvectors) const
    {        
        dsyev2(xx, xy, yy, &eigenvalues[0], &eigenvalues[1], &eigenvectors(0,0), &eigenvectors(0,1));

        eigenvectors(1,0) = -eigenvectors(0,1);
        eigenvectors(1,1) = eigenvectors(0,0);
    }

    template <class T>
    gravis::t2Matrix<T> Tensor2x2<T>::toMatrix() const
    {
        return gravis::t2Matrix<T>(xx,xy,xy,yy);
    }

    template <class T>
    inline
    Tensor2x2<T> operator + (const Tensor2x2<T>& v1, const Tensor2x2<T>& v2)
    {
        return Tensor2x2<T>(
            v1.xx + v2.xx, v1.xy + v2.xy, v1.yy + v2.yy);
    }

    template <class T>
    inline
    Tensor2x2<T> operator - (const Tensor2x2<T>& v1)
    {
        return Tensor2x2<T>(-v1.xx, -v1.xy, -v1.yy);
    }

    template <class T>
    inline
    Tensor2x2<T> operator - (const Tensor2x2<T>& v1, const Tensor2x2<T>& v2)
    {
        return Tensor2x2<T>(
            v1.xx - v2.xx, v1.xy - v2.xy, v1.yy - v2.yy);
	}

	template <class T>
	inline
	Tensor2x2<T> operator * (float f, const Tensor2x2<T>& v)
	{
		return Tensor2x2<T>(f * v.xx, f * v.xy, f * v.yy);
	}

	template <class T>
	inline
	Tensor2x2<T> operator * (const Tensor2x2<T>& v, float f)
	{
		return Tensor2x2<T>(v.xx * f, v.xy * f, v.yy * f);
	}

	template <class T>
	inline
	Tensor2x2<T> operator * (double f, const Tensor2x2<T>& v)
	{
		return Tensor2x2<T>(f * v.xx, f * v.xy, f * v.yy);
	}

	template <class T>
	inline
	Tensor2x2<T> operator * (const Tensor2x2<T>& v, double f)
	{
		return Tensor2x2<T>(v.xx * f, v.xy * f, v.yy * f);
	}

    template <class T>
    inline
    Tensor2x2<T> operator / (const Tensor2x2<T>& v, T f)
    {
        return Tensor2x2<T>(v.xx / f, v.xy / f, v.yy / f);
    }

    template <class T>
    inline
    std::ostream& operator<< (std::ostream& os, const Tensor2x2<T>& arg)
    {
        os << std::setprecision(17) << "["
                << std::setw(8) << arg.xx << ", "
                << std::setw(8) << arg.xy << ", "
                << std::setw(8) << arg.yy << "]";
        return os;
    }


#endif
