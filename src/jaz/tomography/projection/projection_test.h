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

#ifndef PROJECTION_TESTS_H
#define PROJECTION_TESTS_H

#include <src/jaz/image/raw_image.h>
#include <string>
#include <src/jaz/gravis/t3Vector.h>
	
class ProjectionTest
{
	public:
		
		template <typename T>
		static void drawGrid(
			const std::vector<gravis::d4Matrix>& proj,
			RawImage<T>& dest,
			int w, int h, int d,
			gravis::d3Vector origin = gravis::d3Vector(0.0, 0.0, 0.0),
			double spacing = 1.0);

};
			

template <typename T>
void ProjectionTest::drawGrid(
		const std::vector<gravis::d4Matrix>& proj,
		RawImage<T>& dest,
		int w, int h, int d,
		gravis::d3Vector origin,
		double spacing)
{
	const int fc = dest.zdim;
		
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < w; x++)
	{
		gravis::d4Vector pw(
			origin.x + x * spacing, 
			origin.y + y * spacing, 
			origin.z + z * spacing, 
			1.0);
	
		for (int f = 0; f < fc; f++)
		{
			gravis::d4Vector pi = proj[f] * pw;
			const int pxi = (int)round(pi.x);
			const int pyi = (int)round(pi.y);
	
			if (pxi >= 0 && pxi < dest.xdim && pyi >= 0 && pyi < dest.ydim)
			{
				dest(pxi,pyi,f) += 1;
			}
		}
	}
}

#endif

