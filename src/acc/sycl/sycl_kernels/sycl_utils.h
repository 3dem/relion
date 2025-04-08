#ifndef SYCL_UTILITIES_H
#define SYCL_UTILITIES_H

#include <cstring>
#include <cmath>
#include <limits>
#include <cstdint>
#include <sycl/sycl.hpp>

#include "src/macros.h"
#include "src/acc/sycl/sycl_settings.h"

namespace syclKernels
{

template< typename T1, typename T2 >
static inline
int floorfracf(T1 a, T2 b)
{
	return static_cast<int>(a/b);
}

template< typename T1, typename T2 >
static inline
int ceilfracf(T1 a, T2 b)
{
	return static_cast<int>(a/b + 1);
}

// 2D linear interpolation for complex data that interleaves real and
// imaginary data, rather than storing them in a separate array
__attribute__((always_inline))
static inline
void no_tex2D(XFLOAT *mdl, XFLOAT &real, XFLOAT &imag, XFLOAT xp, XFLOAT yp, int mdlX, int mdlInitY)
{
	const int x0 = sycl::floor(xp);
	const XFLOAT fx = xp - x0;

	const int y0 = sycl::floor(yp);
	const XFLOAT fy = yp - y0;

	const int offset1 = 2 * ((y0-mdlInitY)*mdlX + x0);
	const int offset2 = offset1 + 2;
	const int offset3 = offset1 + 2*mdlX;
	const int offset4 = offset3 + 2;

	//-----------------------------
	const XFLOAT d00[2] = { mdl[offset1], mdl[offset1+1] };
	const XFLOAT d01[2] = { mdl[offset2], mdl[offset2+1] };
	const XFLOAT d10[2] = { mdl[offset3], mdl[offset3+1] };
	const XFLOAT d11[2] = { mdl[offset4], mdl[offset4+1] };
	//-----------------------------
	const XFLOAT dx0[2] = { d00[0] + (d01[0] - d00[0])*fx, d00[1] + (d01[1] - d00[1])*fx };
	const XFLOAT dx1[2] = { d10[0] + (d11[0] - d10[0])*fx, d10[1] + (d11[1] - d10[1])*fx };
	//-----------------------------
	real = dx0[0] + (dx1[0] - dx0[0])*fy;
	imag = dx0[1] + (dx1[1] - dx0[1])*fy;
}

// 3D linear interpolation for complex data that interleaves real and
// imaginary data, rather than storing them in a separate array
__attribute__((always_inline))
static inline
void no_tex3D(
				XFLOAT *mdl, XFLOAT &real, XFLOAT &imag,
				XFLOAT xp, XFLOAT yp, XFLOAT zp,
				int mdlX, int mdlXY, int mdlInitY, int mdlInitZ)
{
	const int x0 = sycl::floor(xp);
	const XFLOAT fx = xp - x0;

	const int y0 = sycl::floor(yp);
	const XFLOAT fy = yp - y0;

	const int z0 = sycl::floor(zp);
	const XFLOAT fz = zp - z0;

	const int offset1 = 2*((z0-mdlInitZ)*mdlXY + (y0-mdlInitY)*mdlX + x0);
	const int offset2 = offset1 + 2;
	const int offset3 = offset1 + 2*mdlX;
	const int offset4 = offset3 + 2;
	const int offset5 = offset1 + 2*mdlXY;
	const int offset6 = offset2 + 2*mdlXY;
	const int offset7 = offset3 + 2*mdlXY;
	const int offset8 = offset4 + 2*mdlXY;

	const XFLOAT d000[2] = { mdl[offset1], mdl[offset1+1] };
	const XFLOAT d001[2] = { mdl[offset2], mdl[offset2+1] };
	const XFLOAT d010[2] = { mdl[offset3], mdl[offset3+1] };
	const XFLOAT d011[2] = { mdl[offset4], mdl[offset4+1] };
	const XFLOAT d100[2] = { mdl[offset5], mdl[offset5+1] };
	const XFLOAT d101[2] = { mdl[offset6], mdl[offset6+1] };
	const XFLOAT d110[2] = { mdl[offset7], mdl[offset7+1] };
	const XFLOAT d111[2] = { mdl[offset8], mdl[offset8+1] };

	const XFLOAT dx00[2] = { d000[0] + (d001[0] - d000[0])*fx, d000[1] + (d001[1] - d000[1])*fx };
	const XFLOAT dx10[2] = { d010[0] + (d011[0] - d010[0])*fx, d010[1] + (d011[1] - d010[1])*fx };
	const XFLOAT dx01[2] = { d100[0] + (d101[0] - d100[0])*fx, d100[1] + (d101[1] - d100[1])*fx };
	const XFLOAT dx11[2] = { d110[0] + (d111[0] - d110[0])*fx, d110[1] + (d111[1] - d110[1])*fx };

	const XFLOAT dxy0[2] = {  dx00[0] + (dx10[0] - dx00[0])*fy,
								dx00[1] + (dx10[1] - dx00[1])*fy };
	const XFLOAT dxy1[2] = {  dx01[0] + (dx11[0] - dx01[0])*fy,
								dx01[1] + (dx11[1] - dx01[1])*fy };

	real = dxy0[0] + (dxy1[0] - dxy0[0])*fz;
	imag = dxy0[1] + (dxy1[1] - dxy0[1])*fz;
}

// From https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template <typename T>
bool almostEqual(const T x, const T y, const int ulp)
{
	// the machine epsilon has to be scaled to the magnitude of the values used
	// and multiplied by the desired precision in ULPs (units in the last place)
	return std::fabs(x-y) <= std::numeric_limits<T>::epsilon() * std::fabs(x+y) * ulp
		// unless the result is subnormal
		|| std::fabs(x-y) < std::numeric_limits<T>::min();
}

// Get some ideas from https://bitbashing.io/comparing-floats.html
template <typename T>
int ULP(const T a, const T b)
{
	if (a == b) return 0;

	const auto max = std::numeric_limits<int32_t>::max();
	if ((! std::isfinite(a)) || (! std::isfinite(b)))
		return max;

	if (sizeof(T) == 4)
	{
		int32_t ia, ib;
		std::memcpy(&ia, &a, sizeof(T));
		std::memcpy(&ib, &b, sizeof(T));

		// Don't compare differently-signed floats.
		if ((ia < 0) != (ib < 0)) return max;

		// Return the absolute value of the distance in ULPs.
		return std::abs(ia - ib);
	}
	else if (sizeof(T) == 8)
	{
		int64_t ia, ib;
		std::memcpy(&ia, &a, sizeof(T));
		std::memcpy(&ib, &b, sizeof(T));

		// Don't compare differently-signed floats.
		if ((ia < 0) != (ib < 0))
			return max;

		// Return the absolute value of the distance in ULPs.
		return static_cast<int>(std::abs(ia - ib));
	}
	else
		return max;
}

template <typename T>
static bool checkFinite(const T *ptr, const size_t sz, const char *name)
{
	for (size_t i = 0; i < sz; i++)
	{
		if (! std::isfinite(ptr[i]))
		{
			std::cerr << name << " has inf or nan\n";
			return false;
		}
	}
	return true;
}

template <typename T>
static bool checkDifference(const T *ptrA, const T *ptrB, const size_t sz, const char *name, const int mULP = 16, const int maxcount = 20)
{
	int count = 0;
	for (size_t i = 0; i < sz; i++)
	{
		// Consider up to mULP(=16 by default) ULPs difference as identical
		if (! almostEqual<T>(ptrA[i], ptrB[i], mULP))
		{
			if (count < maxcount)
			{
				int ulp = ULP(ptrA[i], ptrB[i]);
				std::cerr << name << " has difference at [" << i << "]: " << ptrA[i] << " != " << ptrB[i] << "  \t[ulp=" << ULP(ptrA[i], ptrB[i]) << "]" << std::endl;
			}
			else
				return false;

			count++;
		}
	}
	return true;
}

template <typename T>
static size_t countLargerThanNumber(const T *ptrA, const size_t sz, const T val)
{
	size_t count = 0;
	for (size_t i = 0; i < sz; i++)
	{
		// Consider up to 10240 ULP(somewhat large!)difference as identical
		if (ptrA[i] >= val)
			count++;
	}
	std::cerr << count << " / " << sz << std::endl;
	return count;
}

} // end of namespace syclKernels

#endif //SYCL_UTILITIES_H
