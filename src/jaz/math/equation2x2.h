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

#ifndef EQUATION_2X2_H
#define EQUATION_2X2_H

#include <src/jaz/gravis/t2Vector.h>

class Equation2x2
{
	public:
		
		Equation2x2();
		
		double Axx, Axy, Ayy, bx, by;
		
		Equation2x2& operator += (const Equation2x2& arg)
		{
			Axx += arg.Axx;
			Axy += arg.Axy;
			Ayy += arg.Ayy;
			bx += arg.bx;
			by += arg.by;

			return *this;
		}
};

#endif
