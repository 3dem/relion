#include "src/jaz/obs_model.h"
#include "src/jaz/stack_helper.h"
#include "src/jaz/filter_helper.h"
#include "src/jaz/Fourier_helper.h"

#include <src/backprojector.h>


ObservationModel::ObservationModel()
{
}

ObservationModel::ObservationModel(const MetaDataTable &opticsMdt)
:	opticsMdt(opticsMdt),
	angpix(opticsMdt.numberOfObjects()),
	lambda(opticsMdt.numberOfObjects()),
	Cs(opticsMdt.numberOfObjects())
{
	for (int i = 0; i < opticsMdt.numberOfObjects(); i++)
	{
		RFLOAT mag, dstep;
		opticsMdt.getValue(EMDL_CTF_MAGNIFICATION, mag, i);
		opticsMdt.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, i);
		angpix[i] = 10000 * dstep / mag;
		
		double kV;		
		opticsMdt.getValue(EMDL_CTF_VOLTAGE, kV, i);		
		double V = kV * 1e3;
		lambda[i] = 12.2643247 / sqrt(V * (1.0 + V * 0.978466e-6));
		
		opticsMdt.getValue(EMDL_CTF_CS, Cs[i], i);
	}
}

void ObservationModel::predictObservation(
        Projector& proj, const MetaDataTable& mdt, long int particle,
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

    proj.get2DFourierTransform(dest, A3D, false);
	
	int opticsGroup;
	mdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	opticsGroup--;
	
	if (applyShift)
	{
		shiftImageInFourierTransform(dest, dest, s, s/2 - xoff, s/2 - yoff);
	}

    if (applyCtf)
    {
        CTF ctf;
        ctf.readByGroup(mdt, opticsMdt, particle);        

        FilterHelper::modulate(dest, ctf, angpix[opticsGroup]);
    }

    if (applyTilt)
    {
        double tx = 0.0, ty = 0.0;

        mdt.getValue(EMDL_IMAGE_BEAMTILT_X, tx, particle);
        mdt.getValue(EMDL_IMAGE_BEAMTILT_Y, ty, particle);

        if (tx != 0.0 && ty != 0.0)
        {
            selfApplyBeamTilt(dest, -tx, -ty, 
				lambda[opticsGroup], Cs[opticsGroup], angpix[opticsGroup], s);
        }
    }
}


Image<Complex> ObservationModel::predictObservation(
        Projector& proj, const MetaDataTable& mdt, long int particle,
        bool applyCtf, bool applyTilt, bool applyShift) const
{
    Image<Complex> pred;
	
	predictObservation(proj, mdt, particle, pred.data, applyCtf, applyTilt, applyShift);
    return pred;
}

std::vector<Image<Complex>> ObservationModel::predictObservations(
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

void ObservationModel::insertObservation(
		const Image<Complex>& img, BackProjector &bproj,
        const MetaDataTable& mdt, long int particle,
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
	
	int opticsGroup;
	mdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	opticsGroup--;

    if (applyCtf)
    {
        CTF ctf;
        ctf.readByGroup(mdt, opticsMdt, particle);

        ctf.getFftwImage(Fctf, s, s, angpix[opticsGroup]);

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

        selfApplyBeamTilt(F2D, my_tilt_x, my_tilt_y, 
			lambda[opticsGroup], Cs[opticsGroup], angpix[opticsGroup], sh);
    }

    bproj.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
}

double ObservationModel::angToPix(double a, int s, int opticsGroup)
{
	return s * angpix[opticsGroup] / a;
}

double ObservationModel::pixToAng(double p, int s, int opticsGroup)
{
	return s * angpix[opticsGroup] / p;
}

double ObservationModel::getPixelSize(int opticsGroup)
{
	return angpix[opticsGroup];
}

bool ObservationModel::allPixelSizesIdentical()
{
	bool out = true;
	
	for (int i = 1; i < angpix.size(); i++)
	{
		if (angpix[i] != angpix[0])
		{
			out = false;
			break;
		}
	}
	
	return out;
}

bool ObservationModel::containsAllNeededColumns(const MetaDataTable& mdt)
{
	return (mdt.containsLabel(EMDL_ORIENT_ORIGIN_X)
         && mdt.containsLabel(EMDL_ORIENT_ORIGIN_Y)
         && mdt.containsLabel(EMDL_ORIENT_ROT)
         && mdt.containsLabel(EMDL_ORIENT_TILT)
         && mdt.containsLabel(EMDL_ORIENT_PSI)
         && mdt.containsLabel(EMDL_PARTICLE_RANDOM_SUBSET));
}
