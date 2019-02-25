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

#ifndef SVD_HELPER_H
#define SVD_HELPER_H

#include <src/matrix2d.h>

class SvdHelper
{
    public:

        static void decompose(
            const Matrix2D<RFLOAT>& A,
            Matrix2D<RFLOAT>& U,
            Matrix1D<RFLOAT>& S,
            Matrix2D<RFLOAT>& Vt);
};

#endif
