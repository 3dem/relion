#include "legacy_obs_model.h"
#include "stack_helper.h"
#include "img_proc/filter_helper.h"
#include "Fourier_helper.h"

#include <src/backprojector.h>

LegacyObservationModel::LegacyObservationModel()
:   angpix(-1),
    anisoTilt(false)
{
}

LegacyObservationModel::LegacyObservationModel(double angpix, double Cs, double voltage)
:   angpix(angpix),
    lambda(12.2643247 / sqrt(voltage * (1.0 + voltage * 0.978466e-6))),
    Cs(Cs),
    anisoTilt(false)
{
}

void LegacyObservationModel::predictObservation(
        Projector& proj, const MetaDataTable& mdt, int particle,
		MultidimArray<Complex>& dest,
        bool applyCtf, bool applyTilt, bool applyShift) const
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

    Euler_angles2matrix(rot, tilt, psi, A3D);

	if (dest.xdim != sh || dest.ydim != s)
	{
		dest.resize(s,sh);
	}

	dest.initZeros();

    proj.get2DFourierTransform(dest, A3D);

	if (applyShift)
	{
		shiftImageInFourierTransform(dest, dest, s, s/2 - xoff, s/2 - yoff);
	}

    if (applyCtf)
    {
        CTF ctf;
        ctf.read(mdt, mdt, particle);

        FilterHelper::modulate(dest, ctf, angpix);
    }

    if (applyTilt)
    {
        double tx = 0.0, ty = 0.0;

        mdt.getValue(EMDL_IMAGE_BEAMTILT_X, tx, particle);
        mdt.getValue(EMDL_IMAGE_BEAMTILT_Y, ty, particle);

        if (tx != 0.0 && ty != 0.0)
        {
            if (anisoTilt)
            {
                selfApplyBeamTilt(
					dest, -tx, -ty, beamtilt_xx, beamtilt_xy, beamtilt_yy, lambda, Cs, angpix, s);
            }
            else
            {
                selfApplyBeamTilt(dest, -tx, -ty, lambda, Cs, angpix, s);
            }
        }
    }
}


Image<Complex> LegacyObservationModel::predictObservation(
        Projector& proj, const MetaDataTable& mdt, int particle,
        bool applyCtf, bool applyTilt, bool applyShift) const
{
    Image<Complex> pred;

	predictObservation(proj, mdt, particle, pred.data, applyCtf, applyTilt, applyShift);
    return pred;
}

std::vector<Image<Complex>> LegacyObservationModel::predictObservations(
        Projector &proj, const MetaDataTable &mdt, int threads,
        bool applyCtf, bool applyTilt, bool applyShift) const
{
    const int pc = mdt.numberOfObjects();
    std::vector<Image<Complex>> out(pc);

    #pragma omp parallel for num_threads(threads)
    for (int p = 0; p < pc; p++)
    {
        out[p] = predictObservation(proj, mdt, p, applyCtf, applyTilt, applyShift);
    }

    return out;
}

void LegacyObservationModel::insertObservation(const Image<Complex>& img, BackProjector &bproj,
        const MetaDataTable& mdt, int particle,
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
        double my_tilt_x = 0.0;
        double my_tilt_y = 0.0;

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

    bproj.set2DFourierTransform(F2D, A3D, &Fctf);
}

void LegacyObservationModel::setAnisoTilt(double xx, double xy, double yy)
{
    beamtilt_xx = xx;
    beamtilt_xy = xy;
    beamtilt_yy = yy;
    anisoTilt = true;
}

double LegacyObservationModel::angToPix(double a, int s)
{
    return s * angpix / a;
}

double LegacyObservationModel::pixToAng(double p, int s)
{
	return s * angpix / p;
}

bool LegacyObservationModel::containsAllNeededColumns(const MetaDataTable& mdt)
{
	return (mdt.containsLabel(EMDL_ORIENT_ORIGIN_X)
         && mdt.containsLabel(EMDL_ORIENT_ORIGIN_Y)
         && mdt.containsLabel(EMDL_ORIENT_ROT)
         && mdt.containsLabel(EMDL_ORIENT_TILT)
         && mdt.containsLabel(EMDL_ORIENT_PSI)
         && mdt.containsLabel(EMDL_PARTICLE_RANDOM_SUBSET));
}
