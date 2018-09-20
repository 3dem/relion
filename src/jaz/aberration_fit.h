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

#ifndef ABERRATION_FIT_H
#define ABERRATION_FIT_H

#include <src/image.h>
#include <src/metadata_table.h>

class AberrationBasis
{
    public:

        AberrationBasis(int dims);

            std::vector<double> coefficients;

        virtual void getBasisValues(double x, double y, double* dest) = 0;
        virtual void _offsetCtf(double local_Cs, double lambda,
                       double rad_azimuth, double defocus_average, double defocus_deviation,
                       double K1, double K2, double K3, MetaDataTable& mdt, int particle) = 0;

        void offsetCtf(MetaDataTable& mdt, int particle);
};

class OriginalBasis : public AberrationBasis
{
    public:

        OriginalBasis();

        void getBasisValues(double x, double y, double* dest);
        void _offsetCtf(double local_Cs, double lambda,
                   double rad_azimuth, double defocus_average, double defocus_deviation,
                   double K1, double K2, double K3, MetaDataTable& mdt, int particle);
};

class AberrationFit
{
    public:

        static OriginalBasis fitBasic(Image<RFLOAT> phase, Image<RFLOAT> weight, double angpix);
        static Image<RFLOAT> draw(AberrationBasis* fit, double angpix, int s);
};

#endif
