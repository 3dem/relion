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

#include "new_reference_map.h"
#include "obs_model.h"
#include "legacy_obs_model.h"
#include "img_proc/image_op.h"
#include "refinement_helper.h"
#include "motion/motion_helper.h"
#include "new_ft.h"
#include "vtk_helper.h"
#include <src/jaz/math/fft.h>
#include <src/jaz/image/radial_avg.h>
#include <src/jaz/image/translation.h>
#include <src/jaz/image/padding.h>
#include <src/args.h>
#include <omp.h>

using namespace gravis;

NewReferenceMap::NewReferenceMap()
:	maskFn(""),
	fscFn(""),
	paddingFactor(2.0),
	hasMask(false)
{
}

void NewReferenceMap::read(IOParser& parser, int argc, char* argv[])
{
	phase_file_names[0] = parser.getOption("--m1", "Reference map, half 1", "");
	phase_file_names[1] = parser.getOption("--m2", "Reference map, half 2", "");
	amplitude_file_names[0] = parser.getOption("--a1", "Amplitude reference map, half 1", "");
	amplitude_file_names[1] = parser.getOption("--a2", "Amplitude reference map, half 2", "");
	angpix = textToDouble(parser.getOption("--angpix_ref", "Pixel size of the reference map", "-1"));
	maskFn = parser.getOption("--mask", "Reference mask", "");
	fscFn = parser.getOption("--f", "Input STAR file with the FSC of the reference (usually from PostProcess)", "");
	paddingFactor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));

	dualContrast = amplitude_file_names[0] != "";
}

void NewReferenceMap::load(int verb, bool debug)
{
	if (phase_file_names[0] == "" || phase_file_names[1] == "")
	{
		// Get half maps and masks from the PostProcess STAR file.
		FileName fn_half1, fn_half2, fn_mask;
		MetaDataTable MD;
		MD.read(fscFn, "general");
		bool star_is_valid =
				MD.getValue(EMDL_POSTPROCESS_UNFIL_HALFMAP1, fn_half1) &&
				MD.getValue(EMDL_POSTPROCESS_UNFIL_HALFMAP2, fn_half2) &&
				MD.getValue(EMDL_MASK_NAME, fn_mask);

		if (star_is_valid)
		{
			if (verb > 0)
			{
				std::cout << " + The names of the reference half maps and the mask were taken from the PostProcess STAR file.\n";
				std::cout << "   - Half map 1: " << fn_half1 << "\n";
				std::cout << "   - Half map 2: " << fn_half2 << "\n";
				std::cout << "   - Mask: " << fn_mask << std::endl;
			}
			phase_file_names[0] = fn_half1;
			phase_file_names[1] = fn_half2;
			maskFn = fn_mask;
		}
		else
		{
			REPORT_ERROR("could not get filenames for unfiltered half maps from the postprocess STAR file.");
		}
	}

	BufferedImage<RFLOAT> phase_RS[2], amp_RS[2];

	for (int half = 0; half < 2; half++)
	{
		phase_RS[half].read(phase_file_names[half]);

		if (!phase_RS[half].isCubical())
		{
			REPORT_ERROR(phase_file_names[0] + " is not cubical.\n");
		}

		presharpen(phase_RS[half]);

		if (dualContrast)
		{
			amp_RS[half].read(amplitude_file_names[half]);

			if (!amp_RS[half].isCubical())
			{
				REPORT_ERROR(amplitude_file_names[0] + " is not cubical.\n");
			}

			presharpen(amp_RS[half]);
		}
	}

	if (!phase_RS[0].hasEqualSize(phase_RS[1]))
	{
		REPORT_ERROR(phase_file_names[0] + " and " + phase_file_names[1]
				+ " are of unequal size.\n");
	}

	if (dualContrast && !amp_RS[0].hasEqualSize(amp_RS[1]))
	{
		REPORT_ERROR(amplitude_file_names[0] + " and " + amplitude_file_names[1]
				+ " are of unequal size.\n");
	}

	if (angpix < 0)
	{
		Image<RFLOAT> map0;
		map0.read(phase_file_names[0], false);
		angpix = map0.samplingRateX();

		std::cerr
			<< "WARNING: You did not specify a pixel size (--angpix_ref). Assuming "
			<< angpix << " A, as indicated in the image header of " << phase_file_names[0] << std::endl;
	}

	s = phase_RS[0].xdim;
	sh = s/2 + 1;

	if (maskFn != "")
	{
		if (verb > 0) std::cout << " + Masking references ...\n";

		mask.read(maskFn);

		for (int half = 0; half < 2; half++)
		{
			phase_RS[half] *= mask;

			if (dualContrast)
			{
				amp_RS[half] *= mask;
			}
		}

		hasMask = true;
	}

	if (paddingFactor != 1.0)
	{
		if (verb > 0) std::cout << " + Padding reference maps ...\n";

		const int border = (int)((paddingFactor * s - s) / 2.0 + 0.5);
		BufferedImage<RFLOAT> temp;

		for (int half = 0; half < 2; half++)
		{
			temp = phase_RS[half];
			phase_RS[half] = Padding::padCenter3D_full(temp, border);

			if (dualContrast)
			{
				temp = amp_RS[half];
				amp_RS[half] = Padding::padCenter3D_full(temp, border);
			}
		}
	}

	if (verb > 0) std::cout << " + Transforming references ...\n";

	for (int half = 0; half < 2; half++)
	{
		FFT::FourierTransform(phase_RS[half], phaseMap[half]);
		Translation::shiftByHalf(phaseMap[half]);

		if (dualContrast)
		{
			FFT::FourierTransform(amp_RS[half], amplitudeMap[half]);
			Translation::shiftByHalf(amplitudeMap[half]);
		}
	}

	freqWeight = BufferedImage<RFLOAT>(sh,s);

	if (fscFn != "")
	{
		drawFSC(fscFn, freqWeight1D, freqWeight);
	}
	else
	{
		freqWeight1D = std::vector<double>(sh,1.0);
		freqWeight.fill(1.0);
	}

	k_out = sh;

	for (int i = 1; i < sh; i++)
	{
		if (freqWeight1D[i] <= 0.0)
		{
			k_out = i;
			break;
		}
	}
}

Image<RFLOAT> NewReferenceMap::getHollowWeight(
		double kmin_ang, int s_out, double angpix_out)
{
	const int sh_out = s_out/2 + 1;

	Image<RFLOAT> out(sh_out, s_out);

	const double as_out = s_out * angpix_out;
	const double as_ref = s * angpix;

	for (int y = 0; y < s_out; y++)
	for (int x = 0; x < sh_out; x++)
	{
		const double x_out = x;
		const double y_out = y <= sh_out? y : y - s_out;

		const double x_ang = x_out / as_out;
		const double y_ang = y_out / as_out;

		const double x_ref = x_ang * as_ref;
		const double y_ref = y_ang * as_ref;

		const int xx_ref = (int)(x_ref + 0.5);
		const int yy_ref = y_ref >= 0.0? (int)(y_ref + 0.5) : (int)(y_ref + s + 0.5);

		double r = sqrt(x_ang * x_ang + y_ang * y_ang);

		if (r < 1.0 / kmin_ang || xx_ref >= sh || yy_ref < 0 || yy_ref >= s)
		{
			out(y,x) = 0.0;
		}
		else
		{
			out(y,x) = freqWeight(yy_ref, xx_ref);
		}
	}

	return out;
}

std::vector<Image<Complex>> NewReferenceMap::predictAll(
		const MetaDataTable& mdt,
		ObservationModel& obs,
		HalfSet hs, int threads,
		bool applyCtf, bool applyTilt, bool applyShift, bool applyMtf, bool applyCtfPadding)
{
	std::vector<Image<Complex>> out(mdt.numberOfObjects());

	const int pc = mdt.numberOfObjects();

	#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
	{
		out[p] = predict(mdt, p, obs, hs, applyCtf, applyTilt, applyShift, applyMtf, applyCtfPadding);
	}

	return out;
}

Image<Complex> NewReferenceMap::predict(
		const MetaDataTable& mdt, int p,
		ObservationModel& obs,
		HalfSet hs,
		bool applyCtf, bool applyTilt, bool applyShift, bool applyMtf, bool applyCtfPadding)
{
	const int og = obs.getOpticsGroup(mdt,p);
	const int s_og = obs.getBoxSize(og);
	const int sh_og = s/2 + 1;

	Image<Complex> pred(sh_og, s_og);

	int randSubset;
	mdt.getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
	randSubset -= 1;

	int pi = (hs == Own)? randSubset : 1 - randSubset;

	RawImage<Complex> temp(pred());

	BufferedImage<Complex>* ampMap = dualContrast? &amplitudeMap[pi] : 0;

	obs.predictObservation_DC(
		phaseMap[pi], ampMap, mdt, p, temp, angpix,
		applyCtf, applyTilt, applyShift, applyMtf, applyCtfPadding);

	return pred;
}

std::vector<Volume<gravis::t2Vector<Complex>>> NewReferenceMap::predictAllComplexGradients(
		const MetaDataTable &mdt,
		ObservationModel &obs,
		NewReferenceMap::HalfSet hs,
		int threads,
		bool applyCtf, bool applyTilt, bool applyShift, bool applyMtf, bool applyCtfPadding)
{
	std::vector<Volume<t2Vector<Complex>>> out(mdt.numberOfObjects());

	/*const int pc = mdt.numberOfObjects();

	#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
	{
		out[p] = predictComplexGradient(mdt, p, obs, hs, applyCtf, applyTilt, applyShift, applyMtf, applyCtfPadding);
	}*/

	return out;
}

Volume<t2Vector<Complex>> NewReferenceMap::predictComplexGradient(
		const MetaDataTable &mdt,
		int p, ObservationModel &obs,
		NewReferenceMap::HalfSet hs,
		bool applyCtf, bool applyTilt, bool applyShift, bool applyMtf, bool applyCtfPadding)
{
	Volume<t2Vector<Complex>> pred;

	/*int randSubset;
	mdt.getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
	randSubset -= 1;

	int pi = (hs == Own)? randSubset : 1 - randSubset;

	pred = obs.predictComplexGradient(projectors[pi], mdt, p, angpix, applyCtf, applyTilt, applyShift, applyMtf, applyCtfPadding);*/

	return pred;
}

double NewReferenceMap::angToPix(double a) const
{
	return s * angpix / a;
}

double NewReferenceMap::pixToAng(double p) const
{
	return s * angpix / p;
}

void NewReferenceMap::drawFSC(
		std::string fsc_filename,
		std::vector<double>& dest1D,
		RawImage<RFLOAT>& dest2D,
		double thresh) const
{
	MetaDataTable fscMdt;
	fscMdt.read(fsc_filename, "fsc");

	if (!fscMdt.containsLabel(EMDL_SPECTRAL_IDX))
	{
		REPORT_ERROR(fsc_filename + " does not contain a value for "
					 + EMDL::label2Str(EMDL_SPECTRAL_IDX));
	}

	if (!fscMdt.containsLabel(EMDL_POSTPROCESS_FSC_TRUE))
	{
		REPORT_ERROR(fsc_filename + " does not contain a value for "
					 + EMDL::label2Str(EMDL_POSTPROCESS_FSC_TRUE));
	}

	const int sh = fscMdt.numberOfObjects();
	const int s = 2*(sh-1);

	if (sh != dest2D.xdim)
	{
		REPORT_ERROR_STR(
			fsc_filename << " has a wrong number of FSC shells ("
			<< sh << " instead of " << dest2D.xdim << ")");
	}

	dest1D = std::vector<double>(sh);

	for (int i = 0; i < sh; i++)
	{
		int idx;
		fscMdt.getValue(EMDL_SPECTRAL_IDX, idx, i);
		fscMdt.getValue(EMDL_POSTPROCESS_FSC_TRUE, dest1D[i], i);

		if (dest1D[i] < thresh) dest1D[i] = 0.0;
	}

	for (int y = 0; y < s; y++)
	for (int x = 0; x < sh; x++)
	{
		dest2D(x,y) = RadialAvg::interpolate_FftwHalf_3D_lin(x,y,0,s,s,s,dest1D);
	}
}

void NewReferenceMap::presharpen(BufferedImage<RFLOAT>& map)
{
	const int s = map.xdim;

	const double epsilon = 0.02;

	for (int z = 0; z < s; z++)
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const double xx = x - s/2;
		const double yy = y - s/2;
		const double zz = z - s/2;

		const double r = sqrt(xx*xx + yy*yy + zz*zz);

		if (r > 0.)
		{
			const double rval = r / (s * paddingFactor);
			const double sinc = sin(PI * rval) / (PI * rval);

			map(x,y,z) /= sinc * sinc + epsilon;
		}
	}
}
