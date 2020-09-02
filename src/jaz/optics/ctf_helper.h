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

#ifndef CTF_HELPER_H
#define CTF_HELPER_H

#include <src/ctf.h>
#include <src/jaz/image/buffered_image.h>
#include <src/image.h>
#include <vector>

class CtfHelper
{
	public:
		
		template <typename T>
		struct CTFP_CTFQ_Pair
		{
			BufferedImage<tComplex<T>> pq, qp;
		};

		template <typename T>
		static CTFP_CTFQ_Pair<T> stitchHalves(
		        const RawImage<tComplex<T>>& CTFP,
		        const RawImage<tComplex<T>>& CTFQ);

};

template <typename T>
CtfHelper::CTFP_CTFQ_Pair<T> CtfHelper::stitchHalves(
        const RawImage<tComplex<T>>& CTFP,
        const RawImage<tComplex<T>>& CTFQ)
{
	CtfHelper::CTFP_CTFQ_Pair<T> out;
	
	const int s = CTFP.ydim;
	
	out.pq = BufferedImage<tComplex<T>>(s,s);
	out.qp = BufferedImage<tComplex<T>>(s,s);
	
	const int m = s/2;
	
	for (int xo = 0; xo < s; xo++)
	for (int yo = 0; yo < s; yo++)
	{
		const int xg = xo - m;
		const int yg = yo - m;
		
		const int xn = -xg;
		const int yn = (s - yg) % s;
		
		const int xp = xg;
		const int yp = (s + yg) % s;
		
		if (xg > 0 || (xg == 0 && yg >= 0))
		{
			out.pq(xo,yo) = CTFP(xp,yp);
			out.qp(xo,yo) = CTFQ(xp,yp);
		}
		else
		{
			out.pq(xo,yo) = CTFQ(xn,yn).conj();
			out.qp(xo,yo) = CTFP(xn,yn).conj();
		}
	}
	
	return out;
}

#endif
