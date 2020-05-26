#ifndef DUAL_CONTRAST_VOXEL_H
#define DUAL_CONTRAST_VOXEL_H

#include <src/complex.h>
#include <src/jaz/gravis/t2Matrix.h>

template <typename T>
class DualContrastVoxel
{
	public:

		DualContrastVoxel();

			tComplex<T> data_sin, data_cos;
			T weight_sin2, weight_sin_cos, weight_cos2;

		std::pair<tComplex<T>, tComplex<T>> solve(double WienerOffset);

		inline DualContrastVoxel<T>& operator += (const DualContrastVoxel<T>& rhs)
		{
			data_sin += rhs.data_sin;
			data_cos += rhs.data_cos;
			weight_sin2 += rhs.weight_sin2;
			weight_sin_cos += rhs.weight_sin_cos;
			weight_cos2 += rhs.weight_cos2;
			return *this;
		}

		inline DualContrastVoxel<T> operator + (const DualContrastVoxel<T>& rhs) const
		{
			DualContrastVoxel<T> out;

			out.data_sin = data_sin + rhs.data_sin;
			out.data_cos = data_cos + rhs.data_cos;
			out.weight_sin2 = weight_sin2 + rhs.weight_sin2;
			out.weight_sin_cos = weight_sin_cos + rhs.weight_sin_cos;
			out.weight_cos2 = weight_cos2 + rhs.weight_cos2;

			return out;
		}

		inline DualContrastVoxel<T> operator * (T scale) const
		{
			DualContrastVoxel<T> out;

			out.data_sin = scale * data_sin;
			out.data_cos = scale * data_cos;
			out.weight_sin2 = scale * weight_sin2;
			out.weight_sin_cos = scale * weight_sin_cos;
			out.weight_cos2 = scale * weight_cos2;

			return out;
		}

		DualContrastVoxel<T> conj() const
		{
			DualContrastVoxel<T> out;

			out.data_sin = data_sin.conj();
			out.data_cos = data_cos.conj();
			out.weight_sin2 = weight_sin2;
			out.weight_sin_cos = weight_sin_cos;
			out.weight_cos2 = weight_cos2;

			return out;
		}
};

template <typename T>
DualContrastVoxel<T> :: DualContrastVoxel()
	:	data_sin(0,0), data_cos(0,0),
		weight_sin2(0), weight_sin_cos(0), weight_cos2(0)
{
}

template <typename T>
std::pair<tComplex<T>, tComplex<T>> DualContrastVoxel<T> :: solve(double WienerOffset)
{
	const gravis::t2Matrix<T> A(
		weight_sin2, weight_sin_cos,
		weight_sin_cos, weight_cos2 );

	gravis::t2Matrix<T> A_inv = A;

	A_inv(0,0) += WienerOffset;
	A_inv(1,1) += WienerOffset;

	const double det = A_inv(0,0) * A_inv(1,1) - A_inv(0,1) * A_inv(1,0);

	if (std::abs(det) > 1e-16)
	{
		A_inv.invert();

		const gravis::t2Vector<T> data_real(data_sin.real, data_cos.real);
		const gravis::t2Vector<T> data_imag(data_sin.imag, data_cos.imag);

		const gravis::t2Vector<T> out_real = A_inv * data_real;
		const gravis::t2Vector<T> out_imag = A_inv * data_imag;

		return std::make_pair(
			tComplex<T>(out_real[0], out_imag[0]),
			tComplex<T>(out_real[1], out_imag[1]));
	}
	else
	{
		return std::make_pair(
			tComplex<T>(0,0),
			tComplex<T>(0,0));
	}
}


template <class T1, class T2> inline
DualContrastVoxel<T1> operator + (const DualContrastVoxel<T1>& v, const DualContrastVoxel<T2>& w)
{
	DualContrastVoxel<T1> out;

	out.data_sin =       v.data_sin       + w.data_sin;
	out.data_cos =       v.data_cos       + w.data_cos;
	out.weight_sin2 =    v.weight_sin2    + w.weight_sin2;
	out.weight_sin_cos = v.weight_sin_cos + w.weight_sin_cos;
	out.weight_cos2 =    v.weight_cos2    + w.weight_cos2;

	return out;
}

template <class T1, class T2> inline
DualContrastVoxel<T1> operator - (const DualContrastVoxel<T1>& v, const DualContrastVoxel<T2>& w)
{
	DualContrastVoxel<T1> out;

	out.data_sin =       v.data_sin       - w.data_sin;
	out.data_cos =       v.data_cos       - w.data_cos;
	out.weight_sin2 =    v.weight_sin2    - w.weight_sin2;
	out.weight_sin_cos = v.weight_sin_cos - w.weight_sin_cos;
	out.weight_cos2 =    v.weight_cos2    - w.weight_cos2;

	return out;
}

template <class T1, class T2> inline
DualContrastVoxel<T1> operator * (const DualContrastVoxel<T1>& v, const T2& scale)
{
	DualContrastVoxel<T1> out;

	out.data_sin = scale * v.data_sin;
	out.data_cos = scale * v.data_cos;
	out.weight_sin2 = scale * v.weight_sin2;
	out.weight_sin_cos = scale * v.weight_sin_cos;
	out.weight_cos2 = scale * v.weight_cos2;

	return out;
}

template <class T1, class T2> inline
DualContrastVoxel<T1> operator * (const T2& scale, const DualContrastVoxel<T1>& v)
{
	DualContrastVoxel<T1> out;

	out.data_sin = scale * v.data_sin;
	out.data_cos = scale * v.data_cos;
	out.weight_sin2 = scale * v.weight_sin2;
	out.weight_sin_cos = scale * v.weight_sin_cos;
	out.weight_cos2 = scale * v.weight_cos2;

	return out;
}

template <class T1, class T2> inline
DualContrastVoxel<T1> operator / (const DualContrastVoxel<T1>& v, const T2& scale)
{
	DualContrastVoxel<T1> out;

	out.data_sin = v.data_sin / scale;
	out.data_cos = v.data_cos / scale;
	out.weight_sin2 = v.weight_sin2 / scale;
	out.weight_sin_cos = v.weight_sin_cos / scale;
	out.weight_cos2 = v.weight_cos2 / scale;

	return out;
}

#endif
