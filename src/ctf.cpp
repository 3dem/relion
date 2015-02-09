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

/* Read -------------------------------------------------------------------- */
void CTF::read(MetaDataTable &MD1, MetaDataTable &MD2, long int objectID)
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

	initialise();
}
void CTF::setValues(double _defU, double _defV, double _defAng, double _voltage,
		double _Cs, double _Q0, double _Bfac, double _scale)
{
	kV              = _voltage;
	DeltafU         = _defU;
	DeltafV         = _defV;
	azimuthal_angle = _defAng;
	Cs              = _Cs;
	Bfac            = _Bfac;
	scale           = _scale;
	Q0              = _Q0;

	initialise();
}
/* Read from 1 MetaDataTable ----------------------------------------------- */
void CTF::read(MetaDataTable &MD)
{
	MetaDataTable MDempty;
	MDempty.addObject(); // add one empty object
	read(MD, MDempty);
}

/* Write ------------------------------------------------------------------- */
void CTF::write(std::ostream &out)
{
    MetaDataTable MD;
    MD.addObject();
    MD.setValue(EMDL_CTF_VOLTAGE, kV);
    MD.setValue(EMDL_CTF_DEFOCUSU, DeltafU);
    MD.setValue(EMDL_CTF_DEFOCUSV, DeltafV);
    MD.setValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle);
    MD.setValue(EMDL_CTF_CS, Cs);
    MD.setValue(EMDL_CTF_BFACTOR, Bfac);
    MD.setValue(EMDL_CTF_SCALEFACTOR, scale);
    MD.setValue(EMDL_CTF_Q0, Q0);
    MD.write(out);
}

/* Default values ---------------------------------------------------------- */
void CTF::clear()
{
    kV = 200;
    DeltafU = DeltafV = azimuthal_angle = 0;
    Cs = Bfac = 0;
    Q0 = 0;
    scale = 1;
}

/* Initialise the CTF ------------------------------------------------------ */
void CTF::initialise()
{

	// Change units
    double local_Cs = Cs * 1e7;
    double local_kV = kV * 1e3;
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
    K3 = sqrt(1-Q0*Q0);
    K4 = -Bfac / 4.;

    if (Q0 < 0. || Q0 > 1.)
    	REPORT_ERROR("CTF::initialise ERROR: AmplitudeContrast Q0 cannot be smaller than zero or larger than one!");

    if (ABS(DeltafU) < 1e-6 && ABS(DeltafV) < 1e-6 && ABS(Q0) < 1e-6 && ABS(Cs) < 1e-6)
    	REPORT_ERROR("CTF::initialise: ERROR: CTF initialises to all-zero values. Was a correct STAR file provided?");

}

/* Generate a complete CTF Image ------------------------------------------------------ */
void CTF::getFftwImage(MultidimArray<double> &result, int orixdim, int oriydim, double angpix,
		    		bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak, bool do_damping)
{

	double xs = (double)orixdim * angpix;
	double ys = (double)oriydim * angpix;
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result)
	{
		double x = (double)jp / xs;
		double y = (double)ip / ys;
		DIRECT_A2D_ELEM(result, i, j) = getCTF(x, y, do_abs, do_only_flip_phases, do_intact_until_first_peak, do_damping);
	}
}

void CTF::getCenteredImage(MultidimArray<double> &result, double Tm,
		    		bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak, bool do_damping)
{
	result.setXmippOrigin();
	double xs = (double)XSIZE(result) * Tm;
	double ys = (double)YSIZE(result) * Tm;

	FOR_ALL_ELEMENTS_IN_ARRAY2D(result)
	{
		double x = (double)j / xs;
		double y = (double)i / ys;
		A2D_ELEM(result, i, j) = getCTF(x, y, do_abs, do_only_flip_phases, do_intact_until_first_peak, do_damping);
	}

}

