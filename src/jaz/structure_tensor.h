#ifndef STRUCTURE_TENSOR_H
#define STRUCTURE_TENSOR_H

#include <src/jaz/volume.h>
#include <src/jaz/tensor3x3.h>
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
