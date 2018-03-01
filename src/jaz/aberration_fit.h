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
};

#endif
