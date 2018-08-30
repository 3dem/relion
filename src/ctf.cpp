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
/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "src/ctf.h"
#include "src/args.h"
#include "src/fftw.h"
#include "src/metadata_table.h"
#include <src/jaz/obs_model.h>

/* Read -------------------------------------------------------------------- */
void CTF::readByGroup(const MetaDataTable &partMdt, ObservationModel* obs, int particle)
{
	opticsGroup = 0;
	
	if (obs != 0) 
	{
		partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	}
	
	opticsGroup--;		
	
	readValue(EMDL_CTF_VOLTAGE,       kV,              200,     particle, opticsGroup, partMdt, obs);
	readValue(EMDL_CTF_DEFOCUSU,      DeltafU,         0,       particle, opticsGroup, partMdt, obs);
	readValue(EMDL_CTF_DEFOCUSV,      DeltafV,         DeltafU, particle, opticsGroup, partMdt, obs);
	readValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, 0,       particle, opticsGroup, partMdt, obs);
	readValue(EMDL_CTF_CS,            Cs,              0,       particle, opticsGroup, partMdt, obs);
	readValue(EMDL_CTF_BFACTOR,       Bfac,            0,       particle, opticsGroup, partMdt, obs);
	readValue(EMDL_CTF_SCALEFACTOR,   scale,           1,       particle, opticsGroup, partMdt, obs);
	readValue(EMDL_CTF_Q0,            Q0,              0,       particle, opticsGroup, partMdt, obs);
	readValue(EMDL_CTF_PHASESHIFT,    phase_shift,     0,       particle, opticsGroup, partMdt, obs);
		
	initialise();
	
	obsModel = obs;
}

void CTF::readValue(EMDLabel label, RFLOAT& dest, RFLOAT defaultVal, int particle, int opticsGroup, 
					const MetaDataTable& partMdt, const ObservationModel* obs)
{
	if (!partMdt.getValue(label, dest, particle))
	{
		if (opticsGroup < 0 || !obs->opticsMdt.getValue(label, dest, opticsGroup))
		{
			dest = defaultVal;
		}
	}
}

void CTF::read(const MetaDataTable &MD1, const MetaDataTable &MD2, long int objectID)
{

	if (!MD1.getValue(EMDL_CTF_VOLTAGE, kV, objectID))
		if (!MD2.getValue(EMDL_CTF_VOLTAGE, kV, objectID))
			kV=200;

	if (!MD1.getValue(EMDL_CTF_DEFOCUSU, DeltafU, objectID))
		if (!MD2.getValue(EMDL_CTF_DEFOCUSU, DeltafU, objectID))
			DeltafU=0;

	if (!MD1.getValue(EMDL_CTF_DEFOCUSV, DeltafV, objectID))
		if (!MD2.getValue(EMDL_CTF_DEFOCUSV, DeltafV, objectID))
			DeltafV=DeltafU;

	if (!MD1.getValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, objectID))
		if (!MD2.getValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, objectID))
			azimuthal_angle=0;

	if (!MD1.getValue(EMDL_CTF_CS, Cs, objectID))
		if (!MD2.getValue(EMDL_CTF_CS, Cs, objectID))
			Cs=0;

	if (!MD1.getValue(EMDL_CTF_BFACTOR, Bfac, objectID))
		if (!MD2.getValue(EMDL_CTF_BFACTOR, Bfac, objectID))
			Bfac=0;

	if (!MD1.getValue(EMDL_CTF_SCALEFACTOR, scale, objectID))
		if (!MD2.getValue(EMDL_CTF_SCALEFACTOR, scale, objectID))
			scale=1;

	if (!MD1.getValue(EMDL_CTF_Q0, Q0, objectID))
		if (!MD2.getValue(EMDL_CTF_Q0, Q0, objectID))
			Q0=0;

	if (!MD1.getValue(EMDL_CTF_PHASESHIFT, phase_shift, objectID))
		if (!MD2.getValue(EMDL_CTF_PHASESHIFT, phase_shift, objectID))
			phase_shift=0;

	initialise();

}
void CTF::setValues(RFLOAT _defU, RFLOAT _defV, RFLOAT _defAng, RFLOAT _voltage,
		RFLOAT _Cs, RFLOAT _Q0, RFLOAT _Bfac, RFLOAT _scale, RFLOAT _phase_shift)
{
	kV              = _voltage;
	DeltafU         = _defU;
	DeltafV         = _defV;
	azimuthal_angle = _defAng;
	Cs              = _Cs;
	Bfac            = _Bfac;
	scale           = _scale;
	Q0              = _Q0;
	phase_shift     = _phase_shift;

	initialise();
}
/* Read from 1 MetaDataTable ----------------------------------------------- */
void CTF::read(const MetaDataTable &MD)
{
	MetaDataTable MDempty;
	MDempty.addObject(); // add one empty object
	read(MD, MDempty);
}

/** Write to an existing object in a MetaDataTable. */
void CTF::write(MetaDataTable &MD)
{
    MD.setValue(EMDL_CTF_VOLTAGE, kV);
    MD.setValue(EMDL_CTF_DEFOCUSU, DeltafU);
    MD.setValue(EMDL_CTF_DEFOCUSV, DeltafV);
    MD.setValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle);
    MD.setValue(EMDL_CTF_CS, Cs);
    MD.setValue(EMDL_CTF_BFACTOR, Bfac);
    MD.setValue(EMDL_CTF_SCALEFACTOR, scale);
    MD.setValue(EMDL_CTF_PHASESHIFT, phase_shift);
    MD.setValue(EMDL_CTF_Q0, Q0);
}

/* Write ------------------------------------------------------------------- */
void CTF::write(std::ostream &out)
{
    MetaDataTable MD;
	MD.addObject();
    write(MD);
    MD.write(out);
}

/* Initialise the CTF ------------------------------------------------------ */
void CTF::initialise()
{
	// Change units
    RFLOAT local_Cs = Cs * 1e7;
    RFLOAT local_kV = kV * 1e3;
    rad_azimuth = DEG2RAD(azimuthal_angle);

    // Average focus and deviation
    defocus_average   = -(DeltafU + DeltafV) * 0.5;
    defocus_deviation = -(DeltafU - DeltafV) * 0.5;

    // lambda=h/sqrt(2*m*e*kV)
    //    h: Planck constant
    //    m: electron mass
    //    e: electron charge
    // lambda=0.387832/sqrt(kV*(1.+0.000978466*kV)); // Hewz: Angstroms
    // lambda=h/sqrt(2*m*e*kV)
    lambda=12.2643247 / sqrt(local_kV * (1. + local_kV * 0.978466e-6)); // See http://en.wikipedia.org/wiki/Electron_diffraction

    // Helpful constants
    // ICE: X(u)=-PI/2*deltaf(u)*lambda*u^2+PI/2*Cs*lambda^3*u^4
    //          = K1*deltaf(u)*u^2         +K2*u^4
    K1 = PI / 2 * 2 * lambda;
    K2 = PI / 2 * local_Cs * lambda * lambda * lambda;
    K3 = atan(Q0/sqrt(1-Q0*Q0));

    K4 = -Bfac / 4.;

    // Phase shift in radian
    K5 = DEG2RAD(phase_shift);

    if (Q0 < 0. || Q0 > 1.)
    	REPORT_ERROR("CTF::initialise ERROR: AmplitudeContrast Q0 cannot be smaller than zero or larger than one!");

    if (ABS(DeltafU) < 1e-6 && ABS(DeltafV) < 1e-6 && ABS(Q0) < 1e-6 && ABS(Cs) < 1e-6)
    	REPORT_ERROR("CTF::initialise: ERROR: CTF initialises to all-zero values. Was a correct STAR file provided?");

}

double CTF::getGamma(double X, double Y)
{
    RFLOAT u2 = X * X + Y * Y;
    RFLOAT u4 = u2 * u2;

    RFLOAT deltaf = getDeltaF(X, Y);
    return K1 * deltaf * u2 + K2 * u4 - K5 - K3;
}

RFLOAT CTF::getCtfFreq(RFLOAT X, RFLOAT Y)
{
    RFLOAT u2 = X * X + Y * Y;
    RFLOAT u = sqrt(u2);

    RFLOAT deltaf = getDeltaF(X, Y);

    return 2.0 * K1 * deltaf * u + 4.0 * K2 * u * u * u;
}

/* Generate a complete CTF Image ------------------------------------------------------ */
void CTF::getFftwImage(MultidimArray<RFLOAT> &result, int orixdim, int oriydim, RFLOAT angpix,
		    		bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak, bool do_damping)
{
	RFLOAT xs = (RFLOAT)orixdim * angpix;
	RFLOAT ys = (RFLOAT)oriydim * angpix;
	
	if (obsModel != 0 && obsModel->hasEvenZernike)
	{
		const int s = result.ydim;
		const int sh = result.xdim;
		
		if (orixdim != oriydim)
		{
			REPORT_ERROR_STR("CTF::getFftwImage: symmetric aberrations are currently only "
			              << "supported for square images.\n");
		}
		
		const Image<RFLOAT>& gammaOffset = obsModel->getGammaOffset(opticsGroup, oriydim);
		
		if (   gammaOffset.data.xdim < result.xdim
		    || gammaOffset.data.ydim < result.ydim)
		{
			REPORT_ERROR_STR("CTF::getFftwImage: requested output image is larger than the original: "
				<< gammaOffset.data.xdim << "x" << gammaOffset.data.ydim << " available, "
				<< result.xdim << "x" << result.ydim << " requested\n");
		}
		
		for (int yy = 0; yy < result.ydim; yy++)
		for (int xx = 0; xx < result.xdim; xx++)
		{
			RFLOAT x = xx / xs;
			RFLOAT y = yy < sh? yy / ys : (yy - s) / ys;
			
			DIRECT_A2D_ELEM(result, yy, xx) = 
				getCTF(x, y, do_abs, do_only_flip_phases, 
					   do_intact_until_first_peak, do_damping, gammaOffset(yy,xx));
		}
	}
	else
	{
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result)
		{
			RFLOAT x = (RFLOAT)jp / xs;
			RFLOAT y = (RFLOAT)ip / ys;
			
			DIRECT_A2D_ELEM(result, i, j) = 
				getCTF(x, y, do_abs, do_only_flip_phases, 
					   do_intact_until_first_peak, do_damping);
		}
	}
}

/* Generate a complete CTFP (complex) image (with sector along angle) ------------------------------------------------------ */
void CTF::getCTFPImage(MultidimArray<Complex> &result, int orixdim, int oriydim, RFLOAT angpix,
					bool is_positive, float angle)
{

	if (angle < 0 || angle >= 360.)
		REPORT_ERROR("CTF::getCTFPImage: angle should be in [0,360>");
	// Angles larger than 180, are the inverse of the other half!
	if (angle >= 180.)
	{
		angle -= 180.;
		is_positive = !is_positive;
	}

	float anglerad = DEG2RAD(angle);

	RFLOAT xs = (RFLOAT)orixdim * angpix;
	RFLOAT ys = (RFLOAT)oriydim * angpix;
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result)
	{
		RFLOAT x = (RFLOAT)jp / xs;
		RFLOAT y = (RFLOAT)ip / ys;
		RFLOAT myangle = (x*x+y*y > 0) ? acos(y/sqrt(x*x+y*y)) : 0; // dot-product with Y-axis: (0,1)
		if (myangle >= anglerad)
			DIRECT_A2D_ELEM(result, i, j) = getCTFP(x, y, is_positive);
		else
			DIRECT_A2D_ELEM(result, i, j) = getCTFP(x, y, !is_positive);
	}

	// Special line along the vertical Y-axis, where FFTW stores both Friedel mates and Friedel symmetry needs to remain
	if (angle == 0.)
	{
		int dim = YSIZE(result);
		int hdim = dim/2;
		for (int i = hdim + 1; i < dim; i++)
			DIRECT_A2D_ELEM(result, i, 0) = conj(DIRECT_A2D_ELEM(result, dim-i, 0));
	}

}

void CTF::getCenteredImage(MultidimArray<RFLOAT> &result, RFLOAT Tm,
		    		bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak, bool do_damping)
{
	result.setXmippOrigin();
	RFLOAT xs = (RFLOAT)XSIZE(result) * Tm;
	RFLOAT ys = (RFLOAT)YSIZE(result) * Tm;

	FOR_ALL_ELEMENTS_IN_ARRAY2D(result)
	{
		RFLOAT x = (RFLOAT)j / xs;
		RFLOAT y = (RFLOAT)i / ys;
		A2D_ELEM(result, i, j) = getCTF(x, y, do_abs, do_only_flip_phases, do_intact_until_first_peak, do_damping);
	}

}
void CTF::get1DProfile(MultidimArray < RFLOAT > &result, RFLOAT angle, RFLOAT Tm,
		bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak, bool do_damping)
{

	result.setXmippOrigin();
	RFLOAT xs = (RFLOAT)XSIZE(result) * Tm; // assuming result is at the image size!

	FOR_ALL_ELEMENTS_IN_ARRAY1D(result)
	{
		RFLOAT x = (COSD(angle) * (RFLOAT)i) / xs;
		RFLOAT y = (SIND(angle) * (RFLOAT)i) / xs;
		A1D_ELEM(result, i) = getCTF(x, y, do_abs, do_only_flip_phases, do_intact_until_first_peak, do_damping);
	}

}
void CTF::applyWeightEwaldSphereCurvature(MultidimArray < RFLOAT > &result, int orixdim, int oriydim,
		RFLOAT angpix, RFLOAT particle_diameter)
{
	RFLOAT xs = (RFLOAT)orixdim * angpix;
	RFLOAT ys = (RFLOAT)oriydim * angpix;
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result)
	{
		RFLOAT x = (RFLOAT)jp / xs;
		RFLOAT y = (RFLOAT)ip / ys;
		RFLOAT deltaf = fabs(getDeltaF(x, y));
		RFLOAT inv_d = sqrt(x*x + y*y);
		RFLOAT aux = (2.*deltaf*lambda*inv_d)/(particle_diameter);
		RFLOAT A = (aux > 1.) ? 0. : (2./PI) * (acos(aux) - aux * sin(acos(aux)));
		DIRECT_A2D_ELEM(result, i, j) = 1. + A * (2.*fabs(getCTF(x, y)) - 1.);
		// Keep everything on the same scale inside RELION, where we use sin(chi), not 2sin(chi)
		DIRECT_A2D_ELEM(result, i, j) *= 0.5;

	}

}
