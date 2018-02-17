#include "src/jaz/obs_model.h"
#include "src/jaz/stack_helper.h"
#include "src/jaz/filter_helper.h"
#include "src/jaz/Fourier_helper.h"

ObservationModel::ObservationModel()
:   angpix(-1),
    hasTilt(false),
    anisoTilt(false)
{
}

ObservationModel::ObservationModel(double angpix)
:   angpix(angpix),
    hasTilt(false),
    anisoTilt(false)
{
}

ObservationModel::ObservationModel(double angpix, double Cs, double voltage, double beamtilt_x, double beamtilt_y)
:   angpix(angpix),
    lambda(12.2643247 / sqrt(voltage * (1.0 + voltage * 0.978466e-6))),
    Cs(Cs), beamtilt_x(beamtilt_x), beamtilt_y(beamtilt_y),
    hasTilt(true), anisoTilt(false)
{
}

Image<Complex> ObservationModel::predictObservation(
        Projector& proj, MetaDataTable& mdt, int particle,
        bool applyCtf, bool applyTilt,
        double deltaRot, double deltaTilt, double deltaPsi) const
{
    const int s = proj.ori_size;
    const int sh = s/2 + 1;

    double xoff, yoff;

    mdt.getValue(EMDL_ORIENT_ORIGIN_X, xoff, particle);
    mdt.getValue(EMDL_ORIENT_ORIGIN_Y, yoff, particle);

    double rot, tilt, psi;

    Matrix2D<RFLOAT> A3D;
    mdt.getValue(EMDL_ORIENT_ROT, rot, particle);
    mdt.getValue(EMDL_ORIENT_TILT, tilt, particle);
    mdt.getValue(EMDL_ORIENT_PSI, psi, particle);

    rot += deltaRot;
    tilt += deltaTilt;
    psi += deltaPsi;

    Euler_angles2matrix(rot, tilt, psi, A3D);

    Image<Complex> pred(sh, s);
    pred.data.initZeros();

    proj.get2DFourierTransform(pred.data, A3D, false);
    shiftImageInFourierTransform(pred(), pred(), s, s/2 - xoff, s/2 - yoff);

    if (applyCtf)
    {
        CTF ctf;
        ctf.read(mdt, mdt, particle);        

        FilterHelper::modulate(pred, ctf, angpix);
    }

    if (applyTilt)
    {
        double tx, ty;

        if (hasTilt)
        {
            tx = beamtilt_x;
            ty = beamtilt_y;
        }
        else
        {
            mdt.getValue(EMDL_IMAGE_BEAMTILT_X, tx, particle);
            mdt.getValue(EMDL_IMAGE_BEAMTILT_Y, ty, particle);
        }

        if (anisoTilt)
        {
            selfApplyBeamTilt(pred.data, -tx, -ty,
                              beamtilt_xx, beamtilt_xy, beamtilt_yy,
                              lambda, Cs, angpix, s);
        }
        else
        {
            selfApplyBeamTilt(pred.data, -tx, -ty, lambda, Cs, angpix, s);
        }
    }

    return pred;
}

std::vector<Image<Complex>> ObservationModel::predictObservations(
        Projector &proj, MetaDataTable &mdt,
        bool applyCtf, bool applyTilt, int threads) const
{
    const int pc = mdt.numberOfObjects();
    std::vector<Image<Complex>> out(pc);

    #pragma omp parallel for num_threads(threads)
    for (int p = 0; p < pc; p++)
    {
        out[p] = predictObservation(proj, mdt, p, applyCtf, applyTilt);
    }

    return out;
}

void ObservationModel::insertObservation(const Image<Complex>& img, BackProjector &bproj,
        MetaDataTable& mdt, int particle,
        bool applyCtf, bool applyTilt, double shift_x, double shift_y)
{
    const int s = img.data.ydim;
    const int sh = img.data.xdim;

    RFLOAT rot, tilt, psi;
    Matrix2D<RFLOAT> A3D;
    double tx = 0.0, ty = 0.0;

    mdt.getValue(EMDL_ORIENT_ROT, rot, particle);
    mdt.getValue(EMDL_ORIENT_TILT, tilt, particle);
    mdt.getValue(EMDL_ORIENT_PSI, psi, particle);

    Euler_angles2matrix(rot, tilt, psi, A3D);

    mdt.getValue(EMDL_ORIENT_ORIGIN_X, tx, particle);
    mdt.getValue(EMDL_ORIENT_ORIGIN_Y, ty, particle);

    tx += shift_x;
    ty += shift_y;

    MultidimArray<Complex> F2D = img.data;

    shiftImageInFourierTransform(F2D, F2D, s, tx, ty);

    MultidimArray<RFLOAT> Fctf;
    Fctf.resize(F2D);
    Fctf.initConstant(1.);

    if (applyCtf)
    {
        CTF ctf;
        ctf.read(mdt, mdt, particle);

        ctf.getFftwImage(Fctf, s, s, angpix);

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
        {
            DIRECT_MULTIDIM_ELEM(F2D, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
            DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
        }
    }

    if (applyTilt)
    {
        double my_tilt_x = hasTilt? beamtilt_x : 0.0;
        double my_tilt_y = hasTilt? beamtilt_y : 0.0;

        if (mdt.containsLabel(EMDL_IMAGE_BEAMTILT_X))
        {
            mdt.getValue(EMDL_IMAGE_BEAMTILT_X, my_tilt_x, particle);
        }

        if (mdt.containsLabel(EMDL_IMAGE_BEAMTILT_Y))
        {
            mdt.getValue(EMDL_IMAGE_BEAMTILT_Y, my_tilt_y, particle);
        }

        selfApplyBeamTilt(F2D, my_tilt_x, my_tilt_y, lambda, Cs, angpix, sh);
    }

    bproj.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
}

void ObservationModel::setAnisoTilt(double xx, double xy, double yy)
{
    beamtilt_xx = xx;
    beamtilt_xy = xy;
    beamtilt_yy = yy;
    anisoTilt = true;
}
