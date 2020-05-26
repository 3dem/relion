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

		std::pair<tComplex<T>, tComplex<T>> solve();

		inline DualContrastVoxel<T>& operator+= (const DualContrastVoxel<T>& rhs)
		{
			data_sin += rhs.data_sin;
			data_cos += rhs.data_cos;
			weight_sin2 += rhs.weight_sin2;
			weight_sin_cos += rhs.weight_sin_cos;
			weight_cos2 += rhs.weight_cos2;
			return *this;
		}
};

template <typename T>
DualContrastVoxel<T> :: DualContrastVoxel()
	:	data_sin(0,0), data_cos(0,0),
		weight_sin2(0), weight_sin_cos(0), weight_cos2(0)
{
}

template <typename T>
std::pair<tComplex<T>, tComplex<T>> DualContrastVoxel<T> :: solve()
{
	const gravis::t2Matrix<T> A(
		weight_sin2, weight_sin_cos,
		weight_sin_cos, weight_cos2 );

	gravis::t2Matrix<T> A_inv = A;
	A_inv.invert();

	const gravis::t2Vector<T> data_real(data_sin.real, data_cos.real);
	const gravis::t2Vector<T> data_imag(data_sin.imag, data_cos.imag);

	const gravis::t2Vector<T> out_real = A_inv * data_real;
	const gravis::t2Vector<T> out_imag = A_inv * data_imag;

	return std::make_pair(
		tComplex<T>(out_real[0], out_imag[0]),
		tComplex<T>(out_real[1], out_imag[1]));
}

#endif
