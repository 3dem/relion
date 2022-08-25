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

#include "obs_model.h"
#include "stack_helper.h"
#include "img_proc/filter_helper.h"
#include "Fourier_helper.h"
#include <src/jaz/single_particle/ctf/tilt_helper.h>
#include <src/jaz/math/Zernike.h>
#include <src/jaz/tomography/projection/fwd_projection.h>
#include <src/jaz/image/translation.h>
#include "vtk_helper.h"
#include "io/star_converter.h"

#include <src/backprojector.h>

#include <set>
#include <omp.h>

using namespace gravis;


ObservationModel::ObservationModel()
{
}

ObservationModel::ObservationModel(const MetaDataTable &_opticsMdt, bool do_die_upon_error)
:	opticsMdt(_opticsMdt),
	angpix(_opticsMdt.numberOfObjects()),
	lambda(_opticsMdt.numberOfObjects()),
	Cs(_opticsMdt.numberOfObjects()),
	boxSizes(_opticsMdt.numberOfObjects(), 0.0),
	CtfPremultiplied(_opticsMdt.numberOfObjects(), false)
{
	if (!(opticsMdt.containsLabel(EMDL_IMAGE_PIXEL_SIZE) ||
	      opticsMdt.containsLabel(EMDL_MICROGRAPH_PIXEL_SIZE) ||
	      opticsMdt.containsLabel(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE) ||
	      opticsMdt.containsLabel(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE))
	  || !opticsMdt.containsLabel(EMDL_CTF_VOLTAGE)
	  || !opticsMdt.containsLabel(EMDL_CTF_CS))
	{
		if (do_die_upon_error)
		{
			REPORT_ERROR_STR("ERROR: not all necessary variables defined in _optics.star file: "
			              << "rlnVoltage, rlnSphericalAberration and one of rlnImagePixelSize, rlnMicrographPixelSize or rlnOriginalMicrographPixelSize. Make sure to convert older STAR files anew in version-3.1, "
			              << "with relion_convert_star.");
		}
		else
		{
			opticsMdt.clear();
			return;
		}
	}

	// symmetrical high-order aberrations:
	hasEvenZernike = opticsMdt.containsLabel(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS);
	evenZernikeCoeffs = std::vector<std::vector<double> >(opticsMdt.numberOfObjects(), std::vector<double>(0));
	gammaOffset = std::vector<std::map<int,BufferedImage<RFLOAT> > >(opticsMdt.numberOfObjects());

	// antisymmetrical high-order aberrations:
	hasOddZernike = opticsMdt.containsLabel(EMDL_IMAGE_ODD_ZERNIKE_COEFFS);
	oddZernikeCoeffs = std::vector<std::vector<double> >(opticsMdt.numberOfObjects(), std::vector<double>(0));
	phaseCorr = std::vector<std::map<int,BufferedImage<Complex> > >(opticsMdt.numberOfObjects());

	const bool hasTilt = opticsMdt.containsLabel(EMDL_IMAGE_BEAMTILT_X)
	                  || opticsMdt.containsLabel(EMDL_IMAGE_BEAMTILT_Y);

	// anisotropic magnification:
	hasMagMatrices = opticsMdt.containsLabel(EMDL_IMAGE_MAG_MATRIX_00)
	              || opticsMdt.containsLabel(EMDL_IMAGE_MAG_MATRIX_01)
	              || opticsMdt.containsLabel(EMDL_IMAGE_MAG_MATRIX_10)
	              || opticsMdt.containsLabel(EMDL_IMAGE_MAG_MATRIX_11);

	magMatrices.resize(opticsMdt.numberOfObjects());

	hasBoxSizes = opticsMdt.containsLabel(EMDL_IMAGE_SIZE);

	if (opticsMdt.containsLabel(EMDL_IMAGE_OPTICS_GROUP_NAME))
	{
		groupNames.resize(opticsMdt.numberOfObjects());
	}

	if (opticsMdt.containsLabel(EMDL_IMAGE_MTF_FILENAME))
	{
		fnMtfs.resize(opticsMdt.numberOfObjects());
		mtfImage = std::vector<std::map<int,BufferedImage<RFLOAT> > >(opticsMdt.numberOfObjects());
	}
	if (opticsMdt.containsLabel(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE))
	{
		originalAngpix.resize(opticsMdt.numberOfObjects());
	}

	for (int i = 0; i < opticsMdt.numberOfObjects(); i++)
	{
		if (!opticsMdt.getValue(EMDL_IMAGE_PIXEL_SIZE, angpix[i], i))
		{
			if (!opticsMdt.getValue(EMDL_MICROGRAPH_PIXEL_SIZE, angpix[i], i))
			{
				if (!opticsMdt.getValue(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, angpix[i], i))
				{
					opticsMdt.getValue(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, angpix[i], i);
				}
			}
		}

		if (opticsMdt.containsLabel(EMDL_IMAGE_OPTICS_GROUP_NAME))
		{
			opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP_NAME, groupNames[i], i);
		}

		if (opticsMdt.containsLabel(EMDL_IMAGE_MTF_FILENAME))
		{
			opticsMdt.getValue(EMDL_IMAGE_MTF_FILENAME, fnMtfs[i], i);
		}

		if (opticsMdt.containsLabel(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE))
		{
			opticsMdt.getValue(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, originalAngpix[i], i);
		}

		if (opticsMdt.containsLabel(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED))
		{
			bool val;
			opticsMdt.getValue(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, val, i);
			CtfPremultiplied[i] = val;
		}

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
		}

		// always keep a set of mag matrices
		// if none are defined, keep a set of identity matrices

		magMatrices[i] = Matrix2D<RFLOAT>(2,2);
		magMatrices[i].initIdentity();

		// See if there is more than one MTF, for more rapid divideByMtf
		hasMultipleMtfs = false;
		for (int j = 1; j < fnMtfs.size(); j++)
		{
			if (fnMtfs[j] != fnMtfs[0])
			{
				hasMultipleMtfs = true;
				break;
			}
		}

		if (hasMagMatrices)
		{
			opticsMdt.getValue(EMDL_IMAGE_MAG_MATRIX_00, magMatrices[i](0,0), i);
			opticsMdt.getValue(EMDL_IMAGE_MAG_MATRIX_01, magMatrices[i](0,1), i);
			opticsMdt.getValue(EMDL_IMAGE_MAG_MATRIX_10, magMatrices[i](1,0), i);
			opticsMdt.getValue(EMDL_IMAGE_MAG_MATRIX_11, magMatrices[i](1,1), i);
		}
	}

	if (hasTilt) hasOddZernike = true;
}

void ObservationModel::predictObservation_DC(
		const RawImage<Complex>& reference_map,
		const RawImage<Complex>* amplitude_map,
		const MetaDataTable& particle_table,
		long particle_index,
		RawImage<Complex>& destination,
		double reference_pixel_size,
		bool applyCtf,
		bool shiftPhases,
		bool applyShift,
		bool applyMtf,
		bool applyCtfPadding)
{
	const int s_ref = reference_map.ydim;
	const int p = particle_index;
	const int opticsGroup = getOpticsGroup(particle_table, particle_index);

	if (!hasBoxSizes)
	{
		REPORT_ERROR_STR("ObservationModel::predictObservation: Unable to make a prediction without knowing the box size.\n");
	}

	const int s_out = boxSizes[opticsGroup];
	const int sh_out = s_out/2 + 1;
	const double ps = angpix[opticsGroup];

	const double xoff = particle_table.getDouble(EMDL_ORIENT_ORIGIN_X_ANGSTROM, p) / ps;
	const double yoff = particle_table.getDouble(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, p) / ps;

	const double rot = particle_table.getDouble(EMDL_ORIENT_ROT, p);
	const double tilt = particle_table.getDouble(EMDL_ORIENT_TILT, p);
	const double psi = particle_table.getDouble(EMDL_ORIENT_PSI, p);

	Matrix2D<RFLOAT> A3D;

	Euler_angles2matrix(rot, tilt, psi, A3D);

	A3D = applyAnisoMag(A3D, opticsGroup);
	A3D = applyScaleDifference(A3D, opticsGroup, s_ref, reference_pixel_size);

	std::vector<gravis::d4Matrix> Avec(1);
	Avec[0] = gravis::d4Matrix(
				A3D(0,0), A3D(0,1), A3D(0,2), 0.0,
				A3D(1,0), A3D(1,1), A3D(1,2), 0.0,
				A3D(2,0), A3D(2,1), A3D(2,2), 0.0,
					   0,        0,        0, 0.0 );

	if (destination.xdim != sh_out || destination.ydim != s_out)
	{
		REPORT_ERROR_STR(
			"ObservationModel::predictObservation: output image should have a size of "
			<< sh_out << "x" << s_out << " pixels, not "
			<< destination.xdim << "x" << destination.ydim);
	}

	destination.fill(Complex(0,0));

	ForwardProjection::forwardProject(reference_map, Avec, destination, 1);


	if (applyShift)
	{
		Translation::shiftInFourierSpace2D(destination, s_out/2 - xoff, s_out/2 - yoff);
	}

	if (applyCtf)
	{
		const bool dual_contrast = amplitude_map != 0;

		if (dual_contrast)
		{
			BufferedImage<Complex> ampProj(sh_out, s_out);
			ampProj.fill(Complex(0,0));

			ForwardProjection::forwardProject(*amplitude_map, Avec, ampProj, 1);

			CTF ctf;
			ctf.readByGroup(particle_table, this, particle_index);
			ctf.Q0 = 0.0;
			ctf.initialise();

			BufferedImage<RFLOAT> gamma_img(sh_out, s_out);
			ctf.drawGamma(s_out, s_out, ps, &gamma_img[0]);

			for (int y = 0; y < s_out;  y++)
			for (int x = 0; x < sh_out; x++)
			{
				const double g = gamma_img(x,y);
				destination(x,y) = -sin(g) * destination(x,y) - cos(g) * ampProj(x,y);
			}
		}
		else
		{
			CTF ctf;
			ctf.readByGroup(particle_table, this, particle_index);

			Image<RFLOAT> ctfImg_xmipp(sh_out,s_out);

			ctf.getFftwImage(ctfImg_xmipp(), s_out, s_out, angpix[opticsGroup],
					false, false, false, true, applyCtfPadding);

			RawImage<RFLOAT> ctfImg(ctfImg_xmipp);

			if (getCtfPremultiplied(opticsGroup))
			{
				ctfImg *= ctfImg;
			}

			destination *= ctfImg;
		}
	}

	if (shiftPhases && oddZernikeCoeffs.size() > opticsGroup
			&& oddZernikeCoeffs[opticsGroup].size() > 0)
	{
		const BufferedImage<Complex>& corr = getPhaseCorrection(opticsGroup, s_out);
		destination *= corr;
	}

	if (applyMtf && fnMtfs.size() > opticsGroup)
	{
		const BufferedImage<RFLOAT>& mtf = getMtfImage(opticsGroup, s_out);
		destination *= mtf;
	}
}

void ObservationModel::predictObservation(
		Projector& proj, const MetaDataTable& partMdt, long int particle,
		MultidimArray<Complex>& dest, double angpix_ref,
		bool applyCtf, bool shiftPhases, bool applyShift, bool applyMtf, bool applyCtfPadding)
{
	const int s_ref = proj.ori_size;

	int opticsGroup;
	partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	opticsGroup--;

	if (!hasBoxSizes)
	{
		REPORT_ERROR_STR("ObservationModel::predictObservation: Unable to make a prediction without knowing the box size.\n");
	}

	const int s_out = boxSizes[opticsGroup];
	const int sh_out = s_out/2 + 1;

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

	A3D = applyAnisoMag(A3D, opticsGroup);
	A3D = applyScaleDifference(A3D, opticsGroup, s_ref, angpix_ref);

	if (dest.xdim != sh_out || dest.ydim != s_out)
	{
		dest.resize(s_out,sh_out);
	}

	dest.initZeros();

	proj.get2DFourierTransform(dest, A3D);

	if (applyShift)
	{
		shiftImageInFourierTransform(dest, dest, s_out, s_out/2 - xoff, s_out/2 - yoff);
	}

	if (applyCtf)
	{
		CTF ctf;
		ctf.readByGroup(partMdt, this, particle);

		Image<RFLOAT> ctfImg(sh_out,s_out);
		ctf.getFftwImage(ctfImg(), s_out, s_out, angpix[opticsGroup],
				false, false, false, true, applyCtfPadding);

		if (getCtfPremultiplied(opticsGroup))
		{
			for (int y = 0; y < s_out;  y++)
			for (int x = 0; x < sh_out; x++)
			{
				dest(y,x) *= ctfImg(y,x) * ctfImg(y,x);
			}
		}
		else
		{
			for (int y = 0; y < s_out;  y++)
			for (int x = 0; x < sh_out; x++)
			{
				dest(y,x) *= ctfImg(y,x);
			}
		}
	}

	if (shiftPhases && oddZernikeCoeffs.size() > opticsGroup
			&& oddZernikeCoeffs[opticsGroup].size() > 0)
	{
		const BufferedImage<Complex>& corr = getPhaseCorrection(opticsGroup, s_out);

		for (int y = 0; y < s_out;  y++)
		for (int x = 0; x < sh_out; x++)
		{
			dest(y,x) *= corr(x,y);
		}
	}

	if (applyMtf && fnMtfs.size() > opticsGroup)
	{
		const BufferedImage<RFLOAT>& mtf = getMtfImage(opticsGroup, s_out);

		for (int y = 0; y < s_out;  y++)
		for (int x = 0; x < sh_out; x++)
		{
			dest(y,x) *= mtf(x,y);
		}
	}
}

Volume<t2Vector<Complex>> ObservationModel::predictComplexGradient(Projector &proj, const MetaDataTable &partMdt,
                                                                   long particle, double angpix_ref,
                                                                   bool applyCtf, bool shiftPhases, bool applyShift,
                                                                   bool applyMtf, bool applyCtfPadding)
{
	if (applyCtf || applyShift || applyCtfPadding)
	{
		REPORT_ERROR_STR("ObservationModel::predictComplexGradient: "
		              << "applyCtf and applyShift and applyCtfPadding are currently not supported\n");
	}

	const int s_ref = proj.ori_size;

	int opticsGroup;
	partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	opticsGroup--;

	const int s_out = boxSizes[opticsGroup];
	const int sh_out = s_out/2 + 1;

	Volume<t2Vector<Complex>> out(sh_out,s_out,1);

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

	A3D = applyAnisoMag(A3D, opticsGroup);
	A3D = applyScaleDifference(A3D, opticsGroup, s_ref, angpix_ref);

	proj.projectGradient(out, A3D);

	if (shiftPhases && oddZernikeCoeffs.size() > opticsGroup
			&& oddZernikeCoeffs[opticsGroup].size() > 0)
	{
		const BufferedImage<Complex>& corr = getPhaseCorrection(opticsGroup, s_out);

		for (int y = 0; y < s_out;  y++)
		for (int x = 0; x < sh_out; x++)
		{
			out(x,y,0).x *= corr(x,y);
			out(x,y,0).y *= corr(x,y);
		}
	}

	if (applyMtf && fnMtfs.size() > opticsGroup)
	{
		const BufferedImage<RFLOAT>& mtf = getMtfImage(opticsGroup, s_out);

		for (int y = 0; y < s_out;  y++)
		for (int x = 0; x < sh_out; x++)
		{
			out(x,y,0).x *= mtf(x,y);
			out(x,y,0).y *= mtf(x,y);
		}
	}

	return out;
}

void ObservationModel::divideByMtf(
		const MetaDataTable& partMdt, long particle, MultidimArray<Complex>& obsImage,
		bool do_multiply_instead, bool do_correct_average_mtf)
{
	int opticsGroup;
	partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	opticsGroup--;

	divideByMtf(opticsGroup, obsImage, do_multiply_instead, do_correct_average_mtf);
}

void ObservationModel::divideByMtf(
		int opticsGroup, MultidimArray<Complex>& obsImage,
		bool do_multiply_instead, bool do_correct_average_mtf)
{
	const int s = obsImage.ydim;
	const int sh = obsImage.xdim;

	// If there is only a single MTF and we are correcting for the average, then do nothing...
	if (do_correct_average_mtf && !hasMultipleMtfs) return;

	if (fnMtfs.size() > opticsGroup)
	{
		const BufferedImage<RFLOAT>& mtf = getMtfImage(opticsGroup, s);
		const BufferedImage<RFLOAT>& avgmtf = getAverageMtfImage(s);

		if (do_multiply_instead)
		{
			if (do_correct_average_mtf)
			{
				for (int y = 0; y < s;  y++)
				for (int x = 0; x < sh; x++)
				{
					obsImage(y,x) *= mtf(x,y);
					obsImage(y,x) /= avgmtf(x,y);
				}
			}
			else
			{
				for (int y = 0; y < s;  y++)
				for (int x = 0; x < sh; x++)
				{
					obsImage(y,x) *= mtf(x,y);
				}
			}
		}
		else
		{
			if (do_correct_average_mtf)
			{
				for (int y = 0; y < s;  y++)
				for (int x = 0; x < sh; x++)
				{
					obsImage(y,x) /= mtf(x,y);
					obsImage(y,x) *= avgmtf(x,y);
				}
			}
			else
			{
				for (int y = 0; y < s;  y++)
				for (int x = 0; x < sh; x++)
				{
					obsImage(y,x) /= mtf(x,y);
				}
			}
		}
	}
}

void ObservationModel::demodulatePhase(
		const MetaDataTable& partMdt, long particle,
		MultidimArray<Complex>& obsImage,
		bool do_modulate_instead)
{
	int opticsGroup;
	partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle);
	opticsGroup--;

	demodulatePhase(opticsGroup, obsImage, do_modulate_instead);
}

void ObservationModel::demodulatePhase(int opticsGroup, MultidimArray<Complex>& obsImage,
                                       bool do_modulate_instead)
{
	const int s = obsImage.ydim;
	const int sh = obsImage.xdim;

	if (oddZernikeCoeffs.size() > opticsGroup
			&& oddZernikeCoeffs[opticsGroup].size() > 0)
	{
		const BufferedImage<Complex>& corr = getPhaseCorrection(opticsGroup, s);

		if (do_modulate_instead)
		{
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				obsImage(y,x) *= corr(x,y);
			}
		}
		else
		{
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				obsImage(y,x) *= corr(x,y).conj();
			}
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

bool ObservationModel::allBoxSizesIdentical() const
{
	bool out = true;

	for (int i = 1; i < boxSizes.size(); i++)
	{
		if (boxSizes[i] != boxSizes[0])
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

void ObservationModel::setPixelSize(int opticsGroup, RFLOAT newPixelSize)
{
	if (opticsGroup < 0 || opticsGroup >= boxSizes.size())
	{
		REPORT_ERROR("ObservationModel::setPixelSize: wrong opticsGroup");
	}

	angpix[opticsGroup] = newPixelSize;

	phaseCorr[opticsGroup].clear();
	gammaOffset[opticsGroup].clear();

	// mtfImage can be empty
	if (mtfImage.size() > 0)
	{
		mtfImage[opticsGroup].clear();
	}
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

void ObservationModel::setBoxSize(int opticsGroup, int newBoxSize)
{
	if (opticsGroup < 0 || opticsGroup >= boxSizes.size())
	{
		REPORT_ERROR("ObservationModel::setBoxSize: wrong opticsGroup");
	}

	boxSizes[opticsGroup] = newBoxSize;

	phaseCorr[opticsGroup].clear();
	gammaOffset[opticsGroup].clear();

	// mtfImage can be empty
	if (mtfImage.size() > 0)
	{
		mtfImage[opticsGroup].clear();
	}
}

int ObservationModel::getBoxSize(int opticsGroup) const
{
	if (!hasBoxSizes)
	{
		REPORT_ERROR("ObservationModel::getBoxSize: box sizes not available. Make sure particle images are available before converting/importing STAR files from earlier versions of RELION.\n");
	}

	return boxSizes[opticsGroup];
}

void ObservationModel::getBoxSizes(std::vector<int>& sDest, std::vector<int>& shDest) const
{
	if (!hasBoxSizes)
	{
		REPORT_ERROR("ObservationModel::getBoxSizes: box sizes not available. Make sure particle images are available before converting/importing STAR files from earlier versions of RELION.\n");
	}

	sDest.resize(boxSizes.size());
	shDest.resize(boxSizes.size());

	for (int i = 0; i < boxSizes.size(); i++)
	{
		sDest[i] = boxSizes[i];
		shDest[i] = boxSizes[i]/2 + 1;
	}
}

Matrix2D<RFLOAT> ObservationModel::getMagMatrix(int opticsGroup) const
{
	return magMatrices[opticsGroup];
}

void ObservationModel::setMagMatrix(int opticsGroup, const Matrix2D<RFLOAT> &M)
{
	magMatrices[opticsGroup] = M;
}

std::vector<Matrix2D<RFLOAT> > ObservationModel::getMagMatrices() const
{
	return magMatrices;
}

int ObservationModel::getOpticsGroup(const MetaDataTable &particlesMdt, long int particle) const
{
	int opticsGroup = -1;
	if (!particlesMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, particle))
		REPORT_ERROR("ObservationModel::getOpticsGroup: Failed to get optics group for particle #" + std::to_string(particle));
	opticsGroup--;

	return opticsGroup;
}

bool ObservationModel::getCtfPremultiplied(int og) const
{
	if (og < CtfPremultiplied.size())
	{
		return CtfPremultiplied[og];
	}
	else
	{
		return false;
	}
}

void ObservationModel::setCtfPremultiplied(int og, bool val)
{
	CtfPremultiplied[og] = val;
}

std::string ObservationModel::getGroupName(int og)
{
	if (og < groupNames.size())
	{
		return groupNames[og];
	}
	else
	{
		std::stringstream sts;
		sts << (og+1);
		return sts.str();
	}
}

bool ObservationModel::allPixelAndBoxSizesIdentical(const MetaDataTable &mdt)
{
	int og0 = getOpticsGroup(mdt, 0);

	int boxSize0 = getBoxSize(og0);
	double angpix0 = getPixelSize(og0);

	bool allGood = true;

	const int pc = mdt.numberOfObjects();

	for (int p = 1; p < pc; p++)
	{
		int og = getOpticsGroup(mdt, p);

		if (og != og0)
		{
			int boxSize = getBoxSize(og);
			double angpix = getPixelSize(og);

			if (boxSize != boxSize0 || angpix != angpix0)
			{
				allGood = false;
				break;
			}
		}
	}

	return allGood;
}

bool ObservationModel::containsGroup(const MetaDataTable &mdt, int group)
{
	const int pc = mdt.numberOfObjects();

	for (int p = 0; p < pc; p++)
	{
		int og = getOpticsGroup(mdt, p);

		if (og == group)
		{
			return true;
		}
	}

	return false;
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

	for (long int i = 0; i < partMdt.numberOfObjects(); i++)
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

	for (long int i = 0; i < partMdt.numberOfObjects(); i++)
	{
		int og;
		partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, i);
		partMdt.setValue(EMDL_IMAGE_OPTICS_GROUP, old2new[og], i);
	}
}

std::vector<int> ObservationModel::getOptGroupsPresent_oneBased(const MetaDataTable& partMdt) const
{
	const int gc = opticsMdt.numberOfObjects();
	const long long int pc = partMdt.numberOfObjects();

	std::vector<bool> optGroupIsPresent(gc, false);

	for (long int p = 0; p < pc; p++)
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

std::vector<int> ObservationModel::getOptGroupsPresent_zeroBased(const MetaDataTable& partMdt) const
{
	const int gc = opticsMdt.numberOfObjects();
	const long long int pc = partMdt.numberOfObjects();

	std::vector<bool> optGroupIsPresent(gc, false);

	for (long int p = 0; p < pc; p++)
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
			out.push_back(g);
		}
	}

	return out;
}

std::vector<std::pair<int, std::vector<int>>> ObservationModel::splitParticlesByOpticsGroup(const MetaDataTable &partMdt) const
{
	std::vector<int> presentGroups = ObservationModel::getOptGroupsPresent_zeroBased(partMdt);

	const int pogc = presentGroups.size();
	const int ogc = opticsMdt.numberOfObjects();

	std::vector<int> groupToPresentGroup(ogc, -1);

	for (int pog = 0; pog < pogc; pog++)
	{
		const int og = presentGroups[pog];
		groupToPresentGroup[og] = pog;
	}

	std::vector<std::pair<int, std::vector<int>>> out(pogc);

	for (int pog = 0; pog < pogc; pog++)
	{
		out[pog] = std::make_pair(presentGroups[pog], std::vector<int>(0));
	}

	const long long int pc = partMdt.numberOfObjects();

	for (long int p = 0; p < pc; p++)
	{
		int og;
		partMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, p);
		og--;

		int pog = groupToPresentGroup[og];

		out[pog].second.push_back(p);
	}

	return out;
}

const BufferedImage<RFLOAT>& ObservationModel::getMtfImage(int optGroup, int s)
{
	#pragma omp critical(ObservationModel_getMtfImage)
	{
		if (mtfImage[optGroup].find(s) == mtfImage[optGroup].end())
		{
			if (mtfImage[optGroup].size() > 100)
			{
				std::cerr << "Warning: " << (mtfImage[optGroup].size()+1)
				          << " mtf images in cache for the same ObservationModel." << std::endl;
			}

			if (optGroup >= originalAngpix.size())
			{
				REPORT_ERROR("For MTF correction, the rlnMicrographOriginalPixelSize column is necessary in the optics table.");
			}

			MetaDataTable MDmtf;
			MDmtf.read(fnMtfs[optGroup]);
			const int mtfc = MDmtf.numberOfObjects();

			std::vector<RFLOAT> mtf_resol(mtfc), mtf_value(mtfc);

			for (int i = 0; i < mtfc; i++)
			{
				const RFLOAT resol_inv_pixel = MDmtf.getRfloat(EMDL_RESOLUTION_INVPIXEL, i);

				mtf_resol[i] = resol_inv_pixel / originalAngpix[optGroup]; // resolution needs to be given in 1/Ang
				mtf_value[i] = MDmtf.getRfloat(EMDL_POSTPROCESS_MTF_VALUE, i);

				if (mtf_value[i] < 1e-10)
				{
					std::cerr << " i= " << i <<  " mtf_value[i]= " << mtf_value[i] << std::endl;
					REPORT_ERROR("ERROR: zero or negative values encountered in MTF curve: " + fnMtfs[optGroup]);
				}

			}

			// Calculate slope of resolution (in 1/A) per element in the MTF array, in order to interpolate below
			RFLOAT res_per_elem = (mtf_resol[mtfc-1] - mtf_resol[0]) / (RFLOAT)(mtfc);

			if (res_per_elem < 1e-10)
			{
				REPORT_ERROR(" ERROR: the resolution in the MTF star file does not go up....");
			}

			const int sh = s/2 + 1;
			mtfImage[optGroup][s] = BufferedImage<RFLOAT>(sh,s);

			BufferedImage<RFLOAT>& img = mtfImage[optGroup][s];
			const double as = angpix[optGroup] * boxSizes[optGroup];

			// @TODO: replace by scaled radial average (test)

			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				const double xx = x/as;  // logical X-coordinate in 1/A
				const double yy = y < sh-1? y/as : (y-s)/as; // logical Y-coordinate in 1/A

				RFLOAT res = sqrt(xx*xx + yy*yy); // get resolution in 1/Ang
				int i_0 = FLOOR(res / res_per_elem);
				RFLOAT mtf;

				// check boundaries of the array
				if (i_0 >= mtfc-1)
				{
					mtf = mtf_value[mtfc-1];
				}
				else if (i_0 <= 0)
				{
					mtf = mtf_value[0];
				}
				else
				{
					// linear interpolation:
					RFLOAT x_0 = mtf_resol[i_0];
					RFLOAT y_0 = mtf_value[i_0];
					RFLOAT x_1 = mtf_resol[i_0 + 1];
					RFLOAT y_1 = mtf_value[i_0 + 1];

					mtf = y_0 + (y_1 - y_0)*(res - x_0)/(x_1 - x_0);
				}

				img(x,y) = mtf;
			}
		}
	}

	return mtfImage[optGroup][s];
}

const BufferedImage<RFLOAT>& ObservationModel::getAverageMtfImage(int s)
{
	#pragma omp critical(ObservationModel_getAverageMtfImage)
	{
		if (avgMtfImage.find(s) == avgMtfImage.end())
		{
			// get first mtfImage
			avgMtfImage[s] = getMtfImage(0, s);

			// Then add rest of optics groups
			for (int i = 1; i < mtfImage.size(); i++)
			{
				avgMtfImage[s] += getMtfImage(i, s);
			}

			avgMtfImage[s] /= (RFLOAT)mtfImage.size();
		}

	}

	return avgMtfImage[s];
}

const BufferedImage<Complex>& ObservationModel::getPhaseCorrection(int optGroup, int s)
{
	#pragma omp critical(ObservationModel_getPhaseCorrection)
	{
		if (phaseCorr[optGroup].find(s) == phaseCorr[optGroup].end())
		{
			if (phaseCorr[optGroup].size() > 100)
			{
				std::cerr << "Warning: " << (phaseCorr[optGroup].size()+1)
				          << " phase shift images in cache for the same ObservationModel." << std::endl;
			}

			const int sh = s/2 + 1;

			phaseCorr[optGroup][s] = BufferedImage<Complex>(sh,s);

			BufferedImage<Complex>& img = phaseCorr[optGroup][s];
			const double as = angpix[optGroup] * boxSizes[optGroup];
			const Matrix2D<RFLOAT>& M = magMatrices[optGroup];
			
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				double phase = 0.0;

				for (int i = 0; i < oddZernikeCoeffs[optGroup].size(); i++)
				{
					int m, n;
					Zernike::oddIndexToMN(i, m, n);

					const double xx0 = x/as;
					const double yy0 = y < sh-1? y/as : (y-s)/as;
					
					const double xx = M(0,0) * xx0 + M(0,1) * yy0;
					const double yy = M(1,0) * xx0 + M(1,1) * yy0;

					phase += oddZernikeCoeffs[optGroup][i] * Zernike::Z_cart(m,n,xx,yy);
				}

				img(x,y).real = cos(phase);
				img(x,y).imag = sin(phase);
			}
		}
	}

	return phaseCorr[optGroup][s];
}

const BufferedImage<RFLOAT>& ObservationModel::getGammaOffset(int optGroup, int s)
{
	#pragma omp critical(ObservationModel_getGammaOffset)
	{
		if (gammaOffset[optGroup].find(s) == gammaOffset[optGroup].end())
		{
			if (gammaOffset[optGroup].size() > 100)
			{
				std::cerr << "Warning: " << (gammaOffset[optGroup].size()+1)
				          << " gamma offset images in cache for the same ObservationModel." << std::endl;
			}

			const int sh = s/2 + 1;
			gammaOffset[optGroup][s] = BufferedImage<RFLOAT>(sh,s);
			BufferedImage<RFLOAT>& img = gammaOffset[optGroup][s];

			const double as = angpix[optGroup] * boxSizes[optGroup];
			const Matrix2D<RFLOAT>& M = magMatrices[optGroup];

			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				double phase = 0.0;

				for (int i = 0; i < evenZernikeCoeffs[optGroup].size(); i++)
				{
					int m, n;
					Zernike::evenIndexToMN(i, m, n);

					const double xx0 = x/as;
					const double yy0 = y < sh-1? y/as : (y-s)/as;
					
					const double xx = M(0,0) * xx0 + M(0,1) * yy0;
					const double yy = M(1,0) * xx0 + M(1,1) * yy0;

					phase += evenZernikeCoeffs[optGroup][i] * Zernike::Z_cart(m,n,xx,yy);
				}

				img(x,y) = phase;
			}
		}
	}

	return gammaOffset[optGroup][s];
}

Matrix2D<RFLOAT> ObservationModel::applyAnisoMag(Matrix2D<RFLOAT> A3D, int opticsGroup)
{
	Matrix2D<RFLOAT> out;

	if (hasMagMatrices)
	{
		Matrix2D<RFLOAT> mag3D(3,3);
		mag3D.initIdentity();

		mag3D(0,0) = magMatrices[opticsGroup](0,0);
		mag3D(0,1) = magMatrices[opticsGroup](0,1);
		mag3D(1,0) = magMatrices[opticsGroup](1,0);
		mag3D(1,1) = magMatrices[opticsGroup](1,1);
		out = mag3D.inv() * A3D;
	}
	else
	{
		out = A3D;
	}

	return out;
}

Matrix2D<RFLOAT> ObservationModel::applyScaleDifference(
		Matrix2D<RFLOAT> A3D, int opticsGroup, int s3D, double angpix3D)
{
	Matrix2D<RFLOAT> out = A3D;

	out *= (boxSizes[opticsGroup] * angpix[opticsGroup]) / (s3D * angpix3D);

	return out;
}

void ObservationModel::loadSafely(std::string filename, ObservationModel& obsModel,
								  MetaDataTable& particlesMdt, std::string tablename,
								  int verb, bool do_die_upon_error)
{
	MetaDataTable opticsMdt;

	std::string mytablename;

	if (tablename == "discover")
	{
		if (particlesMdt.read(filename, "particles"))
		{
			mytablename = "particles";
		}
		else if (particlesMdt.read(filename, "micrographs"))
		{
			mytablename = "micrographs";
		}
		else if (particlesMdt.read(filename, "movies"))
		{
			mytablename = "movies";
		}
	}
	else
	{
		particlesMdt.read(filename, tablename);
		mytablename = tablename;
	}
	opticsMdt.read(filename, "optics");

	if (opticsMdt.numberOfObjects() == 0)
	{
		if (verb > 0)
		{
			std::cerr << "WARNING: " << filename << " seems to be from a previous version of Relion. Attempting conversion...\n";
			std::cerr << "         You should make sure metadata in the optics group table after conversion is correct.\n";
		}

		MetaDataTable oldMdt;
		oldMdt.read(filename);

		StarConverter::convert_3p0_particlesTo_3p1(oldMdt, particlesMdt, opticsMdt, mytablename, do_die_upon_error);
		if (!do_die_upon_error && opticsMdt.numberOfObjects() == 0) return; // return an empty optics table if error was raised

		if (mytablename == "" || mytablename == "discover")
		{
			if (particlesMdt.containsLabel(EMDL_IMAGE_NAME))
				particlesMdt.setName("particles");
			else if (particlesMdt.containsLabel(EMDL_MICROGRAPH_MOVIE_NAME))
				particlesMdt.setName("movies");
			else
				particlesMdt.setName("micrographs");
		}
	}

	obsModel = ObservationModel(opticsMdt, do_die_upon_error);
	if (!do_die_upon_error && obsModel.opticsMdt.numberOfObjects() == 0) return; // return an empty optics table if error was raised

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
					 filename + ": " + sts.str());
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

	if (mytablename != "particles" && obsModel.opticsMdt.containsLabel(EMDL_IMAGE_PIXEL_SIZE))
	{
		std::cerr << "WARNING: This is not a particle STAR file but contains rlnImagePixelSize column." << std::endl;
		if (!obsModel.opticsMdt.containsLabel(EMDL_MICROGRAPH_PIXEL_SIZE))
		{
			std::cerr << "Pixel size in rlnImagePixelSize will be copied to rlnMicrographPixelSize column. Please make sure this is correct!" << std::endl;

			FOR_ALL_OBJECTS_IN_METADATA_TABLE(obsModel.opticsMdt)
			{
				RFLOAT image_angpix;
				obsModel.opticsMdt.getValue(EMDL_IMAGE_PIXEL_SIZE, image_angpix);
				obsModel.opticsMdt.setValue(EMDL_MICROGRAPH_PIXEL_SIZE, image_angpix);
			}
		}
	}
}

void ObservationModel::saveNew(
		MetaDataTable &particlesMdt,
		MetaDataTable &opticsMdt,
		std::string filename,
		std::string tablename)
{
	std::string tmpfilename = filename + ".tmp";
	std::ofstream of(tmpfilename);

	opticsMdt.setName("optics");
	opticsMdt.write(of);

	particlesMdt.setName(tablename);
	particlesMdt.write(of);

	std::rename(tmpfilename.c_str(), filename.c_str());
}

void ObservationModel::save(MetaDataTable &particlesMdt, std::string filename, std::string tablename)
{
	std::string tmpfilename = filename + ".tmp";
	std::ofstream of(tmpfilename);

	opticsMdt.setName("optics");
	opticsMdt.write(of);

	particlesMdt.setName(tablename);
	particlesMdt.write(of);

	std::rename(tmpfilename.c_str(), filename.c_str());
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
