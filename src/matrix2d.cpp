/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
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

/* Is diagonal ------------------------------------------------------------- */
#include "src/matrix2d.h"

/* Interface to numerical recipes: svbksb ---------------------------------- */
void svbksb(Matrix2D<RFLOAT> &u, Matrix1D<RFLOAT> &w, Matrix2D<RFLOAT> &v,
            Matrix1D<RFLOAT> &b, Matrix1D<RFLOAT> &x)
{
    // Call to the numerical recipes routine. Results will be stored in X
    svbksb(u.adaptForNumericalRecipes2(),
           w.adaptForNumericalRecipes(),
           v.adaptForNumericalRecipes2(),
           u.mdimy, u.mdimx,
           b.adaptForNumericalRecipes(),
           x.adaptForNumericalRecipes());
}

