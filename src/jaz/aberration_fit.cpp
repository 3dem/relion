#include <src/jaz/aberration_fit.h>
#include <src/jaz/vtk_helper.h>


OriginalBasis AberrationFit::fitBasic(
    Image<RFLOAT> phase, Image<RFLOAT> weight, double angpix)
{
    Matrix2D<RFLOAT> A(5,5);
    Matrix1D<RFLOAT> b(5);

    A.initZeros();
    b.initZeros();

    const int sh = phase.data.xdim;
    const int s = phase.data.ydim;

    const double as = angpix * s;

    OriginalBasis basis;
    std::vector<double> vals(5);

    for (int yi = 0; yi < s; yi++)
    for (int xi = 0; xi < sh; xi++)
    {
        const double x = xi/as;
        const double y = yi > sh? (yi - s)/as: yi/as;

        basis.getBasisValues(x, y, &vals[0]);

        const double v = phase(yi,xi);
        const double w = weight(yi,xi);

        for (int r = 0; r < 5; r++)
        {
            b(r) += w * w * vals[r] * v;

            for (int c = 0; c < 5; c++)
            {
                A(r,c) += w * w * vals[r] * vals[c];
            }
        }
    }

    const double tol = 1e-20;
    Matrix1D<RFLOAT> sol(5);
    solve(A, b, sol, tol);

    for (int i = 0; i < 5; i++)
    {
        basis.coefficients[i] = sol(i);
    }

    return basis;
}

Image<RFLOAT> AberrationFit::draw(AberrationBasis *fit, double angpix, int s)
{
    const int sh = s/2 + 1;
    const double as = angpix * s;

    Image<RFLOAT> vis(sh,s);

    std::vector<double> vals(fit->coefficients.size(), 0.0);

    for (int yi = 0; yi < s; yi++)
    for (int xi = 0; xi < sh; xi++)
    {
        const double x = xi/as;
        const double y = yi > sh? (yi - s)/as: yi/as;

        fit->getBasisValues(x, y, &vals[0]);

        double v = 0.0;

        for (int i = 0; i < 5; i++)
        {
            v += fit->coefficients[i] * vals[i];
        }

        vis(yi,xi) = v;
    }

    return vis;
}

AberrationBasis::AberrationBasis(int dims)
    :   coefficients(dims, 0.0)
{}

void AberrationBasis::offsetCtf(MetaDataTable &mdt, int particle)
{
    // identical to CTF::read() and CTF::initialize():
    double kV, DeltafU, DeltafV, azimuthal_angle, Cs, scale, Q0, phase_shift;

    if (!mdt.getValue(EMDL_CTF_VOLTAGE, kV, particle)) kV = 200;
    if (!mdt.getValue(EMDL_CTF_DEFOCUSU, DeltafU, particle)) DeltafU = 0;
    if (!mdt.getValue(EMDL_CTF_DEFOCUSV, DeltafV, particle)) DeltafV = DeltafU;
    if (!mdt.getValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, particle)) azimuthal_angle = 0;
    if (!mdt.getValue(EMDL_CTF_CS, Cs, particle)) Cs = 0;
    //if (!mdt.getValue(EMDL_CTF_SCALEFACTOR, scale, particle)) scale = 1;
    if (!mdt.getValue(EMDL_CTF_Q0, Q0, particle)) Q0 = 0;
    //if (!mdt.getValue(EMDL_CTF_PHASESHIFT, phase_shift, particle)) phase_shift = 0;

    //std::cout << DeltafU << ", " << DeltafV << " @ " << azimuthal_angle << "°, " << Cs << ", " << Q0 << "\n";

    double local_Cs = Cs * 1e7;
    double local_kV = kV * 1e3;
    double rad_azimuth = DEG2RAD(azimuthal_angle);

    double defocus_average   = -(DeltafU + DeltafV) * 0.5;
    double defocus_deviation = -(DeltafU - DeltafV) * 0.5;
    double lambda=12.2643247 / sqrt(local_kV * (1. + local_kV * 0.978466e-6));

    double K1 = (PI / 2) * 2 * lambda;
    double K2 = (PI / 2) * local_Cs * lambda * lambda * lambda;
    double K3 = atan(Q0/sqrt(1-Q0*Q0));

    _offsetCtf(local_Cs, lambda, rad_azimuth, defocus_average, defocus_deviation,
          K1, K2, K3, mdt, particle);
}


OriginalBasis::OriginalBasis()
:   AberrationBasis(5)
{}

void OriginalBasis::getBasisValues(double x, double y, double *dest)
{
    dest[0] = 1.0; // phase shift
    dest[1] = x*x + y*y; // defocus
    dest[2] = x*x - y*y; // oblique astigmatism
    dest[3] = x*y; // vertical astigmatism
    dest[4] = (x*x + y*y)*(x*x + y*y); // primary spherical
}

void OriginalBasis::_offsetCtf(
        double local_Cs, double lambda,
        double rad_azimuth, double defocus_average, double defocus_deviation,
        double K1, double K2, double K3, MetaDataTable &mdt, int particle)
{
    /* from ctf.h:

           RFLOAT u2 = X * X + Y * Y;
           RFLOAT u4 = u2 * u2;
           RFLOAT deltaf = defocus_average + defocus_deviation*cos(2*(atan2(Y, X) - rad_azimuth))

           argument = K1 * deltaf * u2 + K2 * u4 - K5 - K3

                K1 = PI / 2 * 2 * lambda;
                K2 = PI / 2 * local_Cs * lambda * lambda * lambda;
                K3 = atan(Q0/sqrt(1-Q0*Q0));
                K5 = DEG2RAD(phase_shift);

                local_Cs = Cs * 1e7;

       astigmatism/defocus:

           K1 * deltaf * u2
           = K1 * defocus_average * u2 + defocus_deviation * K1 * cos(2*(phi - rad_azimuth)) * u2
           = K1 * defocus_average * u2 + defocus_deviation * K1 * cos(2*phi - 2*rad_azimuth) * u2
           = K1 * defocus_average * u2 + defocus_deviation * K1 * [cos(2*phi) cos(2*rad_azimuth) + sin(2*phi) sin(2*rad_azimuth)] * u2
           = K1 * defocus_average * u2 + defocus_deviation * K1 * [(cos²(phi) - sin²(phi)) cos(2*rad_azimuth) + 2 sin(phi) cos(phi) sin(2*rad_azimuth)] * u2
           = K1 * defocus_average * u2 + defocus_deviation * K1 * [(X² - Y²) cos(2*rad_azimuth) + 2 Y X sin(2*rad_azimuth)]
           = b1 (X² + Y²) + b2 (X² - Y²) + b3 (XY)

           where:  b1 = K1 * defocus_average
                   b2 = K1 * defocus_deviation * cos(2*rad_azimuth)
                   b3 = 2 * K1 * defocus_deviation * sin(2*rad_azimuth)

                   <=>

                   defocus_average = b1 / (PI lambda)
                   defocus_deviation = sqrt(b2² + b3²/4)/(PI lambda)
                   rad_azimuth = atan2(b3/2, b2) / 2                        */

    double b1 = K1 * defocus_average + coefficients[1];
    double b2 = K1 * defocus_deviation * cos(2*rad_azimuth) + coefficients[2];
    double b3 = 2 * K1 * defocus_deviation * sin(2*rad_azimuth) + coefficients[3];

    double new_defocus_average = b1 / (PI * lambda);
    double new_defocus_deviation = sqrt(b2*b2 + b3*b3/4)/(PI*lambda);
    double new_rad_azimuth = atan2(b3/2.0, b2) / 2.0;

    double azimuthal_angle = RAD2DEG(new_rad_azimuth);
    double DeltafU = -new_defocus_average - new_defocus_deviation;
    double DeltafV = new_defocus_deviation - new_defocus_average;

/*     spherical aberration:

           K2 * u4 = b4 * u4
            <=>
           PI / 2 * local_Cs * lambda * lambda * lambda = b4;
            <=>
           local_Cs = 2 * b4 / (PI lambda³)
            <=>
           Cs = 1e-7 * 2 * b4 / (PI lambda³)            */

    double b4 = PI * lambda*lambda*lambda * local_Cs / 2.0 + coefficients[4];
    double Cs = 1e-7 * 2.0 * b4 / (PI * lambda*lambda*lambda);

    /*     phase shift / amp. contrast:

               K3 = atan(Q0/sqrt(1-Q0*Q0))
                <=>
               Q0/sqrt(1-Q0²) = tan(K3)
                <=>
               Q0² = (1 - Q0²) * tan²(K3)
                <=>
               Q0²(1 + tan²K3) = tan²K3
                <=>
               Q0 = sqrt(tan²K3/(1 + tan²K3))
    */

    double b0 = K3 - coefficients[0];
    double t0 = tan(b0);
    double Q0 = sqrt(t0*t0/(1 + t0*t0));

    mdt.setValue(EMDL_CTF_DEFOCUSU, DeltafU, particle);
    mdt.setValue(EMDL_CTF_DEFOCUSV, DeltafV, particle);
    mdt.setValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, particle);
    mdt.setValue(EMDL_CTF_CS, Cs, particle);
    mdt.setValue(EMDL_CTF_Q0, Q0, particle);

    //std::cout << DeltafU << ", " << DeltafV << " @ " << azimuthal_angle << "°, " << Cs << ", " << Q0 << "\n\n";

}
