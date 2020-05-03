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

#include <src/jaz/single_particle/structure_tensor.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>

using namespace gravis;

void StructureTensor :: computeEdgeTensor(const Volume<RFLOAT>& src, Volume<Tensor3x3<RFLOAT> >& dest)
{
    dest.resize(src.dimx, src.dimy, src.dimz);

    FOR_ALL_VOXELS(src)
    {
        t3Vector<RFLOAT> g = FilterHelper::centralGradient(src, x, y, z);
        dest(x,y,z) = Tensor3x3<RFLOAT>::autoDyadicProduct(g);
    }
}

void StructureTensor :: compute(const Volume<RFLOAT>& src, Volume<Tensor3x3<RFLOAT> >& dest, RFLOAT rho)
{
    Volume<Tensor3x3<RFLOAT> > E(src.dimx, src.dimy, src.dimz);
    computeEdgeTensor(src, E);

    dest.resize(src.dimx, src.dimy, src.dimz);

    FilterHelper::separableGaussian(E, dest, rho);
}

void StructureTensor :: computeEigenvalues(const Volume<Tensor3x3<RFLOAT> >& J, Volume<t3Vector<RFLOAT> >& dest)
{
    dest.resize(J.dimx, J.dimy, J.dimz);

    gravis::t3Matrix<RFLOAT> Q;
    gravis::t3Vector<RFLOAT> d;

    FOR_ALL_VOXELS(J)
    {
        J(x,y,z).diagonalize(d, Q);
        dest(x,y,z) = d;
    }
}

void StructureTensor :: computeEigenvalues(const Volume<RFLOAT>& src, Volume<t3Vector<RFLOAT> >& dest, RFLOAT rho)
{
    Volume<Tensor3x3<RFLOAT> > J;
    compute(src, J, rho);

    computeEigenvalues(J, dest);
}

void StructureTensor :: computeSmallestEigenvalue(const Volume<RFLOAT>& src, Volume<RFLOAT>& dest, RFLOAT rho)
{
    Volume<Tensor3x3<RFLOAT> > J;
    compute(src, J, rho);

    dest.resize(J.dimx, J.dimy, J.dimz);

    gravis::t3Matrix<RFLOAT> Q;
    gravis::t3Vector<RFLOAT> d;

    FOR_ALL_VOXELS(J)
    {
        J(x,y,z).diagonalize(d, Q);
        dest(x,y,z) = d.z;
    }
}

void StructureTensor :: computeMiddleEigenvalue(const Volume<RFLOAT>& src, Volume<RFLOAT>& dest, RFLOAT rho)
{
    Volume<Tensor3x3<RFLOAT> > J;
    compute(src, J, rho);

    dest.resize(J.dimx, J.dimy, J.dimz);

    gravis::t3Matrix<RFLOAT> Q;
    gravis::t3Vector<RFLOAT> d;

    FOR_ALL_VOXELS(J)
    {
        J(x,y,z).diagonalize(d, Q);
        dest(x,y,z) = d.y;
    }
}

void StructureTensor :: computeEigenvalueLC(const Volume<RFLOAT>& src, Volume<RFLOAT>& dest, RFLOAT rho, RFLOAT w0, RFLOAT w1, RFLOAT w2)
{
    Volume<Tensor3x3<RFLOAT> > J;
    compute(src, J, rho);

    dest.resize(J.dimx, J.dimy, J.dimz);

    gravis::t3Matrix<RFLOAT> Q;
    gravis::t3Vector<RFLOAT> d;

    FOR_ALL_VOXELS(J)
    {
        J(x,y,z).diagonalize(d, Q);
        dest(x,y,z) = w0*d[0] + w1*d[1] + w2*d[2];
    }
}
