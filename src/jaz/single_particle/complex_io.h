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

#ifndef COMPLEX_IO_H
#define COMPLEX_IO_H

#include <src/image.h>
#include <src/complex.h>
#include <string>

class ComplexIO
{
	public:

	template <typename T>
	static void write(const MultidimArray<tComplex<T> >& img, std::string fnBase, std::string fnSuffix)
	{
		Image<RFLOAT> temp(img.xdim, img.ydim, img.zdim, img.ndim);

		FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img)
		{
			DIRECT_NZYX_ELEM(temp.data, l, k, i, j) = DIRECT_NZYX_ELEM(img, l, k, i, j).real;
		}

		temp.write(fnBase + "_real" + fnSuffix);

		FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img)
		{
			DIRECT_NZYX_ELEM(temp.data, l, k, i, j) = DIRECT_NZYX_ELEM(img, l, k, i, j).imag;
		}

		temp.write(fnBase + "_imag" + fnSuffix);
	}

	template <typename T>
	static void read(Image<tComplex<T> >& img, std::string fnBase, std::string fnSuffix)
	{
		Image<RFLOAT> temp;

		temp.read(fnBase + "_real" + fnSuffix);

		img = Image<Complex>(temp.data.xdim, temp.data.ydim, temp.data.zdim, temp.data.ndim);

		FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img.data)
		{
			DIRECT_NZYX_ELEM(img.data, l, k, i, j).real = DIRECT_NZYX_ELEM(temp.data, l, k, i, j);
		}

		temp.read(fnBase + "_imag" + fnSuffix);

		FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img.data)
		{
			DIRECT_NZYX_ELEM(img.data, l, k, i, j).imag = DIRECT_NZYX_ELEM(temp.data, l, k, i, j);
		}
	}
};

#endif
