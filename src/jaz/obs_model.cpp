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

#include "src/jaz/obs_model.h"
#include "src/jaz/stack_helper.h"
#include "src/jaz/img_proc/filter_helper.h"
#include "src/jaz/Fourier_helper.h"
#include "src/jaz/ctf/tilt_helper.h"
#include "src/jaz/math/Zernike.h"
#include "src/jaz/vtk_helper.h"
#include "src/jaz/io/star_converter.h"

#include <src/backprojector.h>

#include <set>
#include <omp.h>

using namespace gravis;

void ObservationModel::loadSafely(
	std::string filename,
	ObservationModel& obsModel,
	MetaDataTable& particlesMdt, MetaDataTable& opticsMdt, int verb)
{
	particlesMdt.read(filename, "particles");
	opticsMdt.read(filename, "optics");

	if (particlesMdt.numberOfObjects() == 0 && particlesMdt.numberOfObjects() == 0)
	{
		if (verb > 0)
		{
			std::cerr << "WARNING: " << filename << " seems to be from a previous version of Relion. Attempting conversion...\n";
		}

		MetaDataTable oldMdt;
		oldMdt.read(filename);

		StarConverter::convert_3p0_particlesTo_3p1(oldMdt, particlesMdt, opticsMdt);
	}

	obsModel = ObservationModel(opticsMdt);

	// read pixel sizes (and make sure they are all the same)

	if (!obsModel.allPixelSizesIdentical())
	{
		REPORT_ERROR("ERROR: different pixel sizes detected. Please split your dataset by pixel size.");
	}

	// make sure all optics groups are defined

	std::vector<int> undefinedOptGroups = obsModel.findUndefinedOptGroups(particlesMdt);

	if (undefinedOptGroups.size() > 0)
	{
		std::stringstream sts;

		for (int i = 0; i < undefinedOptGroups.size(); i++)
		{
			sts << undefinedOptGroups[i];

			if (i < undefinedOptGroups.size()-1)
			{
				sts << ", ";
			}
		}

		REPORT_ERROR("ERROR: The following optics groups were not defined in "+
					 filename+": "+sts.str());
	}

	// make sure the optics groups appear in the right order (and rename them if necessary)

	if (!obsModel.opticsGroupsSorted())
	{
		if (verb > 0)
		{
			std::cerr << "   - Warning: the optics groups in " << filename
					  << " are not in the right order - renaming them now" << std::endl;
		}

		obsModel.sortOpticsGroups(particlesMdt);
	}
}

void ObservationModel::save(
		MetaDataTable &particlesMdt,
		MetaDataTable &opticsMdt,
		std::string filename)
{
	std::ofstream of(filename);

	opticsMdt.setName("optics");
	opticsMdt.write(of);

	particlesMdt.setName("particles");
	particlesMdt.write(of);
}

ObservationModel::ObservationModel()
{
}

ObservationModel::ObservationModel(const MetaDataTable &opticsMdt)
:	opticsMdt(opticsMdt),
	angpix(opticsMdt.numberOfObjects()),
	lambda(opticsMdt.numberOfObjects()),
	Cs(opticsMdt.numberOfObjects()),
	boxSizes(opticsMdt.numberOfObjects(), 0.0)
{
	if (   !opticsMdt.containsLabel(EMDL_IMAGE_PIXEL_SIZE)
	    || !opticsMdt.containsLabel(EMDL_CTF_VOLTAGE)
	    || !opticsMdt.containsLabel(EMDL_CTF_CS))
	{
		REPORT_ERROR_STR("ERROR: not all necessary variables defined in _optics.star file: "
			<< "rlnPixelSize, rlnVoltage and rlnSphericalAberration.");
	}

	// symmetrical high-order aberrations:
	hasEvenZernike = opticsMdt.containsLabel(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS);
	evenZernikeCoeffs = std::vector<std::vector<double> >(
			opticsMdt.numberOfObjects(), std::vector<double>(0));
	gammaOffset = std::vector<std::map<int,Image<RFLOAT> > >(opticsMdt.numberOfObjects());

	// antisymmetrical high-order aberrations:
	hasOddZernike = opticsMdt.containsLabel(EMDL_IMAGE_ODD_ZERNIKE_COEFFS);
	oddZernikeCoeffs = std::vector<std::vector<double> >(
			opticsMdt.numberOfObjects(), std::vector<double>(0));
	phaseCorr = std::vector<std::map<int,Image<Complex> > >(opticsMdt.numberOfObjects());

	const bool hasTilt = opticsMdt.containsLabel(EMDL_IMAGE_BEAMTILT_X)
	                  || opticsMdt.containsLabel(EMDL_IMAGE_BEAMTILT_Y);

	// anisotropic magnification:
	hasMagMatrices = opticsMdt.containsLabel(EMDL_IMAGE_MAG_MATRIX_00)
			      || opticsMdt.containsLabel(EMDL_IMAGE_MAG_MATRIX_01)
			      || opticsMdt.containsLabel(EMDL_IMAGE_MAG_MATRIX_10)
			      || opticsMdt.containsLabel(EMDL_IMAGE_MAG_MATRIX_11);

	if (hasMagMatrices) magMatrices.resize(opticsMdt.numberOfObjects());
	
	hasBoxSizes = opticsMdt.containsLabel(EMDL_IMAGE_SIZE);

	for (int i = 0; i < opticsMdt.numberOfObjects(); i++)
	{
		opticsMdt.getValue(EMDL_IMAGE_PIXEL_SIZE, angpix[i], i);
		opticsMdt.getValue(EMDL_IMAGE_SIZE, boxSizes[i], i);

		double kV;
		opticsMdt.getValue(EMDL_CTF_VOLTAGE, kV, i);
		double V = kV * 1e3;
		lambda[i] = 12.2643247 / sqrt(V * (1.0 + V * 0.978466e-6));

		opticsMdt.getValue(EMDL_CTF_CS, Cs[i], i);

		if (hasEvenZernike)
		{
			opticsMdt.getValue(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS, evenZernikeCoeffs[i], i);
		}

		if (hasOddZernike)
		{
			opticsMdt.getValue(EMDL_IMAGE_ODD_ZERNIKE_COEFFS, oddZernikeCoeffs[i], i);
		}

		if (hasTilt)
		{
			double tx(0), ty(0);

			opticsMdt.getValue(EMDL_IMAGE_BEAMTILT_X, tx, i);
			opticsMdt.getValue(EMDL_IMAGE_BEAMTILT_Y, ty, i);

			if (!hasOddZernike)
			{
				oddZernikeCoeffs[i] = std::vector<double>(6, 0.0);
			}

			TiltHelper::insertTilt(oddZernikeCoeffs[i], tx, ty, Cs[i], lambda[i]);

			hasOddZernike = true;
		}

		if (hasMagMatrices)
		{
			magMatrices[i] = Matrix2D<RFLOAT>(3,3);
			magMatrices[i].initIdentity();

			// transpose the matrix, since the transpose is used in Projector::get2DFourierTransform
			opticsMdt.getValue(EMDL_IMAGE_MAG_MATRIX_00, magMatrices[i](0,0), i);
			opticsMdt.getValue(EMDL_IMAGE_MAG_MATRIX_01, magMatrices[i](0,1), i);
			opticsMdt.getValue(EMDL_IMAGE_MAG_MATRIX_10, magMatrices[i](1,0), i);
			opticsMdt.getValue(EMDL_IMAGE_MAG_MATRIX_11, magMatrices[i](1,1), i);
		}
	}
}

void ObservationModel::predictObservation(
        Projector& proj, const MetaDataTable& partMdt, long int particle,
		MultidimArray<Complex>& dest,
        bool applyCtf, bool shiftPhases, bool applyShift)
{
    const int s = proj.ori_size;
    const int sh = s/2 + 1;

	int opticsGroup;
	partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	opticsGroup--;

    double xoff, yoff;

    partMdt.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff, particle);
    partMdt.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff, particle);

	xoff /= angpix[opticsGroup];
	yoff /= angpix[opticsGroup];

    double rot, tilt, psi;

    Matrix2D<RFLOAT> A3D;
    partMdt.getValue(EMDL_ORIENT_ROT, rot, particle);
    partMdt.getValue(EMDL_ORIENT_TILT, tilt, particle);
    partMdt.getValue(EMDL_ORIENT_PSI, psi, particle);

    Euler_angles2matrix(rot, tilt, psi, A3D);

	A3D = applyAnisoMagTransp(A3D, opticsGroup);

	if (dest.xdim != sh || dest.ydim != s)
	{
		dest.resize(s,sh);
	}

	dest.initZeros();

    proj.get2DFourierTransform(dest, A3D, false);

	if (applyShift)
	{
		shiftImageInFourierTransform(dest, dest, s, s/2 - xoff, s/2 - yoff);
	}

    if (applyCtf)
    {
        CTF ctf;
        ctf.readByGroup(partMdt, this, particle);

		Image<RFLOAT> ctfImg(sh,s);
		ctf.getFftwImage(ctfImg(), s, s, angpix[opticsGroup]);

		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			dest(y,x) *= ctfImg(y,x);
		}
    }

    if (shiftPhases && oddZernikeCoeffs.size() > opticsGroup
			&& oddZernikeCoeffs[opticsGroup].size() > 0)
    {
		const Image<Complex>& corr = getPhaseCorrection(opticsGroup, s);

		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			dest(y,x) *= corr(y,x);
		}
    }
}

Volume<t2Vector<Complex>> ObservationModel::predictComplexGradient(
		Projector &proj, const MetaDataTable &partMdt, long particle,
		bool applyCtf, bool shiftPhases, bool applyShift)
{
	if (applyCtf || applyShift)
	{
		REPORT_ERROR_STR("ObservationModel::predictComplexGradient: "
						 << "applyCtf and applyShift are currently not supported\n");
	}

	const int s = proj.ori_size;
    const int sh = s/2 + 1;

	Volume<t2Vector<Complex>> out(sh,s,1);

	int opticsGroup;
	partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	opticsGroup--;

	double xoff, yoff;

    partMdt.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff, particle);
    partMdt.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff, particle);

	xoff /= angpix[opticsGroup];
	yoff /= angpix[opticsGroup];

    double rot, tilt, psi;

    Matrix2D<RFLOAT> A3D;
    partMdt.getValue(EMDL_ORIENT_ROT, rot, particle);
    partMdt.getValue(EMDL_ORIENT_TILT, tilt, particle);
    partMdt.getValue(EMDL_ORIENT_PSI, psi, particle);

    Euler_angles2matrix(rot, tilt, psi, A3D);

	A3D = applyAnisoMagTransp(A3D, opticsGroup);

    proj.projectGradient(out, A3D);

    if (shiftPhases && oddZernikeCoeffs.size() > opticsGroup
			&& oddZernikeCoeffs[opticsGroup].size() > 0)
    {
		const Image<Complex>& corr = getPhaseCorrection(opticsGroup, s);

		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			out(x,y,0).x *= corr(y,x);
			out(x,y,0).y *= corr(y,x);
		}
    }

	return out;
}

void ObservationModel::demodulatePhase(
		const MetaDataTable& partMdt, long particle, MultidimArray<Complex>& obsImage)
{
	int opticsGroup;
	partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	opticsGroup--;

	demodulatePhase(opticsGroup, obsImage);
}

void ObservationModel::demodulatePhase(
		int opticsGroup, MultidimArray<Complex>& obsImage)
{
	const int s = obsImage.ydim;
	const int sh = obsImage.xdim;

	if (oddZernikeCoeffs.size() > opticsGroup
			&& oddZernikeCoeffs[opticsGroup].size() > 0)
    {
		const Image<Complex>& corr = getPhaseCorrection(opticsGroup, s);

		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			obsImage(y,x) *= corr(y,x).conj();
		}
    }
}

bool ObservationModel::allPixelSizesIdentical() const
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

double ObservationModel::angToPix(double a, int s, int opticsGroup) const
{
	return s * angpix[opticsGroup] / a;
}

double ObservationModel::pixToAng(double p, int s, int opticsGroup) const
{
	return s * angpix[opticsGroup] / p;
}

double ObservationModel::getPixelSize(int opticsGroup) const
{
	return angpix[opticsGroup];
}

std::vector<double> ObservationModel::getPixelSizes() const
{
	return angpix;
}

double ObservationModel::getWavelength(int opticsGroup) const
{
	return lambda[opticsGroup];
}

std::vector<double> ObservationModel::getWavelengths() const
{
	return lambda;
}

double ObservationModel::getSphericalAberration(int opticsGroup) const
{
	return Cs[opticsGroup];
}

std::vector<double> ObservationModel::getSphericalAberrations() const
{
	return Cs;
}

int ObservationModel::getBoxSize(int opticsGroup) const
{
	if (!hasBoxSizes)
	{
		REPORT_ERROR("ObservationModel::getBoxSize: box sizes not available\n");
	}
	
	return boxSizes[opticsGroup];
}

void ObservationModel::getBoxSizes(std::vector<double>& sDest, std::vector<double>& shDest) const
{
	if (!hasBoxSizes)
	{
		REPORT_ERROR("ObservationModel::getBoxSizes: box sizes not available\n");
	}
	
	sDest.resize(boxSizes.size());
	shDest.resize(boxSizes.size());
	
	for (int i = 0; i < boxSizes.size(); i++)
	{
		sDest[i] = boxSizes[i];
		shDest[i] = boxSizes[i]/2 + 1;
	}
}

Matrix2D<double> ObservationModel::getMagMatrix(int opticsGroup) const
{
	return magMatrices[opticsGroup];
}

std::vector<Matrix2D<double> > ObservationModel::getMagMatrices() const
{
	return magMatrices;
}

int ObservationModel::getOpticsGroup(const MetaDataTable &particlesMdt, int particle) const
{
	int opticsGroup;
	particlesMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	opticsGroup--;
	
	return opticsGroup;
}

int ObservationModel::numberOfOpticsGroups() const
{
	return opticsMdt.numberOfObjects();
}

bool ObservationModel::opticsGroupsSorted() const
{
	for (int i = 0; i < opticsMdt.numberOfObjects(); i++)
	{
		int og;
		opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, i);

		if (og != i+1)
		{
			return false;
		}
	}

	return true;
}

std::vector<int> ObservationModel::findUndefinedOptGroups(const MetaDataTable &partMdt) const
{
	std::set<int> definedGroups;

	for (int i = 0; i < opticsMdt.numberOfObjects(); i++)
	{
		int og;
		opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, i);

		definedGroups.insert(og);
	}

	std::vector<int> out;
	out.reserve(opticsMdt.numberOfObjects());

	for (int i = 0; i < partMdt.numberOfObjects(); i++)
	{
		int og;
		partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, i);

		if (definedGroups.find(og) == definedGroups.end())
		{
			out.push_back(og);
		}
	}

	return out;
}

void ObservationModel::sortOpticsGroups(MetaDataTable& partMdt)
{
	std::map<int,int> old2new;

	for (int i = 0; i < opticsMdt.numberOfObjects(); i++)
	{
		int og;
		opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, i);

		old2new[og] = i+1;

		opticsMdt.setValue(EMDL_IMAGE_OPTICS_GROUP, i+1, i);
	}

	for (int i = 0; i < partMdt.numberOfObjects(); i++)
	{
		int og;
		partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, i);
		partMdt.setValue(EMDL_IMAGE_OPTICS_GROUP, old2new[og], i);
	}
}

std::vector<int> ObservationModel::getOptGroupsPresent(const MetaDataTable& partMdt) const
{
	const int gc = opticsMdt.numberOfObjects();
	const int pc = partMdt.numberOfObjects();

	std::vector<bool> optGroupIsPresent(gc, false);

	for (int p = 0; p < pc; p++)
	{
		int og;
		partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, p);

		optGroupIsPresent[og-1] = true;
	}

	std::vector<int> out(0);
	out.reserve(gc);

	for (int g = 0; g < gc; g++)
	{
		if (optGroupIsPresent[g])
		{
			out.push_back(g+1);
		}
	}

	return out;
}

const Image<Complex>& ObservationModel::getPhaseCorrection(int optGroup, int s)
{
	#pragma omp critical
	{
		if (phaseCorr[optGroup].find(s) == phaseCorr[optGroup].end())
		{
			if (phaseCorr[optGroup].size() > 100)
			{
				std::cerr << "Warning: " << (phaseCorr[optGroup].size()+1)
						  << " phase shift images in cache for the same ObservationModel." << std::endl;
			}

			const int sh = s/2 + 1;
			phaseCorr[optGroup][s] = Image<Complex>(sh,s);
			Image<Complex>& img = phaseCorr[optGroup][s];

			const double as = angpix[optGroup] * s;

			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				double phase = 0.0;

				for (int i = 0; i < oddZernikeCoeffs[optGroup].size(); i++)
				{
					int m, n;
					Zernike::oddIndexToMN(i, m, n);

					const double xx = x/as;
					const double yy = y < sh? y/as : (y-s)/as;

					phase += oddZernikeCoeffs[optGroup][i] * Zernike::Z_cart(m,n,xx,yy);
				}

				img(y,x).real = cos(phase);
				img(y,x).imag = sin(phase);
			}
		}
	}

	return phaseCorr[optGroup][s];
}

const Image<RFLOAT>& ObservationModel::getGammaOffset(int optGroup, int s)
{
	#pragma omp critical
	{
		if (gammaOffset[optGroup].find(s) == gammaOffset[optGroup].end())
		{
			if (gammaOffset[optGroup].size() > 100)
			{
				std::cerr << "Warning: " << (gammaOffset[optGroup].size()+1)
						  << " gamma offset images in cache for the same ObservationModel." << std::endl;
			}

			const int sh = s/2 + 1;
			gammaOffset[optGroup][s] = Image<RFLOAT>(sh,s);
			Image<RFLOAT>& img = gammaOffset[optGroup][s];

			const double as = angpix[optGroup] * s;

			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				double phase = 0.0;

				for (int i = 0; i < evenZernikeCoeffs[optGroup].size(); i++)
				{
					int m, n;
					Zernike::evenIndexToMN(i, m, n);

					const double xx = x/as;
					const double yy = y < sh? y/as : (y-s)/as;

					phase += evenZernikeCoeffs[optGroup][i] * Zernike::Z_cart(m,n,xx,yy);
				}

				img(y,x) = phase;
			}
		}
	}

	return gammaOffset[optGroup][s];
}

Matrix2D<RFLOAT> ObservationModel::applyAnisoMagTransp(
		Matrix2D<RFLOAT> A3D_transp, int opticsGroup, double angpixDest)
{
	Matrix2D<RFLOAT> out;

	if (hasMagMatrices)
	{
		out = magMatrices[opticsGroup].transpose() * A3D_transp;
	}
	else
	{
		out = A3D_transp;
	}
	
	if (angpixDest > 0)
	{
		out *= angpixDest / angpix[opticsGroup];
	}

	return out;
}

bool ObservationModel::containsAllColumnsNeededForPrediction(const MetaDataTable& partMdt)
{
	return (partMdt.containsLabel(EMDL_ORIENT_ORIGIN_X_ANGSTROM)
         && partMdt.containsLabel(EMDL_ORIENT_ORIGIN_Y_ANGSTROM)
         && partMdt.containsLabel(EMDL_ORIENT_ROT)
         && partMdt.containsLabel(EMDL_ORIENT_TILT)
         && partMdt.containsLabel(EMDL_ORIENT_PSI)
         && partMdt.containsLabel(EMDL_PARTICLE_RANDOM_SUBSET));
}
