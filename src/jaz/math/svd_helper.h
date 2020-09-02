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
            const Matrix2D<double>& A,
            Matrix2D<double>& U,
            Matrix1D<double>& S,
            Matrix2D<double>& Vt);
			
		
		/* SVD code: stolen from XMIPP and slightly altered */
		
		/*static void svdcmp(const Matrix2D< double >& a,
					Matrix1D< double >& w,
					Matrix2D< double >& v);*/
				
};

#endif
