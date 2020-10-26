#ifndef DUAL_CONTRAST_VOXEL_H
#define DUAL_CONTRAST_VOXEL_H

#include <src/complex.h>
#include <src/jaz/gravis/t2Matrix.h>

template <typename T>
class DualContrastVoxel
{
	public:

		struct Solution
		{
			tComplex<T> phase, amplitude;
			T conditionNumber;
		};


		DualContrastVoxel();

			tComplex<T> data_sin, data_cos;
			T weight_sin2, weight_sin_cos, weight_cos2;

		Solution solve(
				double WienerOffset,
				double lambda,
				bool isotropicWiener = true) const;

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
typename DualContrastVoxel<T>::Solution
	DualContrastVoxel<T>::solve(
		double WienerOffset, double lambda, bool isotropicWiener) const
{
	const double Q0 = 0.07;

	const gravis::t2Matrix<T> A0(
		weight_sin2, weight_sin_cos,
		weight_sin_cos, weight_cos2 );

	gravis::t2Matrix<T> A = A0;

	A(0,0) += WienerOffset;

	if (isotropicWiener)
	{
		A(1,1) += WienerOffset;
	}
	else
	{
		A(1,1) += WienerOffset / (Q0 * Q0);
	}

	/*
	  Encourage phase (P) and amplitude (M) structure factors to assume a ratio of  M = Q0 * P.
	  Penalise  lambda * |Q0 * P - M|^2
	  by adding  lambda * [Q0, -1]^T [Q0, -1]  to A.
	*/

	A(0,0) += lambda * Q0 * Q0;
	A(1,0) += lambda * Q0;
	A(0,1) += lambda * Q0;
	A(1,1) += lambda;

	Solution out;

	/*
	  Find eigenvalues of A:

	  det(A-tI) = 0
	  = (a00-t)*(a11-t) - a01*a10
	  = t² - (a00+a11)*t + a00*a11 - a01*a10
	   <=>
	  t = [a00+a11 +- sqrt((a00+a11)² - 4(a00*a11 - a01*a10))] / 2
	*/

	{
		const double a00 = A(0,0);
		const double a01 = A(0,1);
		const double a10 = A(1,0);
		const double a11 = A(1,1);

		const double discr = (a00+a11) * (a00+a11) - 4.0 * (a00*a11 - a01*a10);

		if (discr >= 0)
		{
			const double eig0 = (a00+a11 - sqrt(discr)) / 2.0;
			const double eig1 = (a00+a11 + sqrt(discr)) / 2.0;

			out.conditionNumber = eig0 / eig1;
		}
		else
		{
			out.conditionNumber = 0.0;
		}
	}

	const double det = A(0,0) * A(1,1) - A(0,1) * A(1,0);

	if (std::abs(det) > 1e-16)
	{
		gravis::t2Matrix<T> A_inv = A;
		A_inv.invert();

		const gravis::t2Vector<T> data_real(data_sin.real, data_cos.real);
		const gravis::t2Vector<T> data_imag(data_sin.imag, data_cos.imag);

		const gravis::t2Vector<T> out_real = A_inv * data_real;
		const gravis::t2Vector<T> out_imag = A_inv * data_imag;

		out.phase = tComplex<T>(-out_real[0], -out_imag[0]);
		out.amplitude = tComplex<T>(out_real[1], out_imag[1]);
	}
	else
	{

		out.phase = tComplex<T>(0,0);
		out.amplitude = tComplex<T>(0,0);
	}

	return out;
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
