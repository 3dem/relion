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
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/single_particle/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <vector>

class CtfHelper
{
    public:

		
        static std::vector<CTF> loadCtffind4(std::string path, int imageCount,
                                             double voltage = 300.0, double Cs = 2.2,
                                             double Q0 = 0.1, double Bfac = 0.0,
                                             double scale = 1.0);
		
		static CTF setFromFile(
				std::stringstream& line, 
				double voltage, double Cs, double Q0, double Bfac, double scale);
};

#endif
