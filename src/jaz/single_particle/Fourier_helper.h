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

#ifndef FOURIER_HELPER_H
#define FOURIER_HELPER_H

#include <src/image.h>
#include <src/complex.h>

class FourierHelper
{
    public:

        static void FourierShift2D(MultidimArray<Complex>& img, RFLOAT xshift, RFLOAT yshift);
        static void FourierShift2D(MultidimArray<RFLOAT>& img, RFLOAT xshift, RFLOAT yshift);
};

#endif
