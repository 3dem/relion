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

#ifndef STRUCTURE_TENSOR_H
#define STRUCTURE_TENSOR_H

#include <src/jaz/single_particle/volume.h>
#include <src/jaz/math/tensor3x3.h>
#include <src/macros.h>

class StructureTensor
{
    public:

        static void compute(const Volume<RFLOAT>& src, Volume<Tensor3x3<RFLOAT> >& dest, RFLOAT rho);

        static void computeEdgeTensor(const Volume<RFLOAT>& src, Volume<Tensor3x3<RFLOAT> >& dest);
        static void computeEigenvalues(const Volume<Tensor3x3<RFLOAT> >& J, Volume<gravis::t3Vector<RFLOAT> >& dest);

        // helper functions
        static void computeEigenvalues(const Volume<RFLOAT>& src, Volume<gravis::t3Vector<RFLOAT> >& dest, RFLOAT rho);
        static void computeSmallestEigenvalue(const Volume<RFLOAT>& src, Volume<RFLOAT>& dest, RFLOAT rho);
        static void computeMiddleEigenvalue(const Volume<RFLOAT>& src, Volume<RFLOAT>& dest, RFLOAT rho);
        static void computeEigenvalueLC(const Volume<RFLOAT>& src, Volume<RFLOAT>& dest, RFLOAT rho, RFLOAT w0, RFLOAT w1, RFLOAT w2);
};

#endif
