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

#include "reference_map.h"
#include <src/jaz/obs_model.h>
#include <src/jaz/legacy_obs_model.h>
#include <src/jaz/img_proc/image_op.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/motion/motion_helper.h>
#include <src/jaz/new_ft.h>
#include <src/jaz/vtk_helper.h>
#include <src/backprojector.h>
#include <src/args.h>
#include <omp.h>

using namespace gravis;

ReferenceMap::ReferenceMap()
:	reconFn0(""),
	reconFn1(""),
	maskFn(""),
	fscFn(""),
	paddingFactor(2.0),	
	hasMask(false)
{
}

void ReferenceMap::read(IOParser& parser, int argc, char* argv[])
{
	reconFn0 = parser.getOption("--m1", "Reference map, half 1");
	reconFn1 = parser.getOption("--m2", "Reference map, half 2");
	angpix = textToDouble(parser.getOption("--angpix_ref", "Pixel size of the reference map", "-1"));
	maskFn = parser.getOption("--mask", "Reference mask", "");
	fscFn = parser.getOption("--f", "Input STAR file with the FSC of the reference");
	paddingFactor = textToDouble(parser.getOption("--pad", "Padding factor", "2"));
}

void ReferenceMap::load(int verb, bool debug)
{
	Image<RFLOAT> maps[2], powSpec[2];
	
	if (debug) std::cout << "reading: " << reconFn0 << "\n";
	
	maps[0].read(reconFn0);	
	
	if ( maps[0].data.xdim != maps[0].data.ydim
	  || maps[0].data.ydim != maps[0].data.zdim)
	{
		REPORT_ERROR(reconFn0 + " is not cubical.\n");
	}
	
	if (debug) std::cout << "reading: " << reconFn1 << "\n";
	
	maps[1].read(reconFn1);
	
	if ( maps[1].data.xdim != maps[1].data.ydim
	  || maps[1].data.ydim != maps[1].data.zdim)
	{
		REPORT_ERROR(reconFn1 + " is not cubical.\n");
	}
	
	if ( maps[0].data.xdim != maps[1].data.xdim
	  || maps[0].data.ydim != maps[1].data.ydim
	  || maps[0].data.zdim != maps[1].data.zdim)
	{
		REPORT_ERROR(reconFn0 + " and " + reconFn1 + " are of unequal size.\n");
	}

	if (angpix < 0)
	{
		angpix = maps[0].samplingRateX();
		std::cerr << "WARNING: You did not specify --ref_angpix. The pixel size in the image header of " << reconFn0 << ", " << angpix << " A/px, is used." << std::endl;
	}
	s = maps[0].data.ydim;
	sh = s/2 + 1;
	
	if (maskFn != "")
	{
		if (verb > 0) std::cout << " + Masking references ...\n";
		
		Image<RFLOAT> maskedRef;
		
		mask.read(maskFn);
		
		ImageOp::multiply(mask, maps[0], maskedRef);
		maps[0] = maskedRef;
		
		ImageOp::multiply(mask, maps[1], maskedRef);
		maps[1] = maskedRef;
		
		hasMask = true;
	}
	
	if (verb > 0) std::cout << " + Transforming references ...\n";
	
	projectors[0] = Projector(s, TRILINEAR, paddingFactor, 10, 2);
	projectors[0].computeFourierTransformMap(maps[0].data, powSpec[0].data, maps[0].data.xdim);
	
	projectors[1] = Projector(s, TRILINEAR, paddingFactor, 10, 2);
	projectors[1].computeFourierTransformMap(maps[1].data, powSpec[1].data, maps[1].data.xdim);
		
	if (fscFn != "")
	{
		MetaDataTable fscMdt;
		fscMdt.read(fscFn, "fsc");
		
		if (!fscMdt.containsLabel(EMDL_SPECTRAL_IDX))
		{
			REPORT_ERROR(fscFn + " does not contain a value for "
						 + EMDL::label2Str(EMDL_SPECTRAL_IDX));
		}
		
		if (!fscMdt.containsLabel(EMDL_POSTPROCESS_FSC_TRUE))
		{
			REPORT_ERROR(fscFn + " does not contain a value for "
						 + EMDL::label2Str(EMDL_POSTPROCESS_FSC_TRUE));
		}
		
		RefinementHelper::drawFSC(&fscMdt, freqWeight1D, freqWeight);
	}
	else
	{
		freqWeight1D = std::vector<double>(sh,1.0);
		freqWeight = Image<RFLOAT>(sh,s);
		freqWeight.data.initConstant(1.0);
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

Image<RFLOAT> ReferenceMap::getHollowWeight(
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

std::vector<Image<Complex>> ReferenceMap::predictAll(
		const MetaDataTable& mdt,
		ObservationModel& obs,
		HalfSet hs, int threads,
		bool applyCtf, bool applyTilt, bool applyShift)
{
	// declare on first line to prevent copying
	std::vector<Image<Complex>> out(mdt.numberOfObjects());
	
	const int pc = mdt.numberOfObjects();
	
	#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
	{
		out[p] = predict(mdt, p, obs, hs, applyCtf, applyTilt, applyShift);
	}
	
	return out;
}

Image<Complex> ReferenceMap::predict(
		const MetaDataTable& mdt, int p,
		ObservationModel& obs,
		HalfSet hs,
		bool applyCtf, bool applyTilt, bool applyShift)
{
	Image<Complex> pred;
	
	int randSubset;
	mdt.getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
	randSubset -= 1;
	
	int pi = (hs == Own)? randSubset : 1 - randSubset;
	
	obs.predictObservation(projectors[pi], mdt, p, pred(), angpix, applyCtf, applyTilt, applyShift);
	
	return pred;
}

std::vector<Volume<gravis::t2Vector<Complex>>> ReferenceMap::predictAllComplexGradients(
		const MetaDataTable &mdt, 
		ObservationModel &obs, 
		ReferenceMap::HalfSet hs, 
		int threads, 
		bool applyCtf, bool applyTilt, bool applyShift)
{
	// declare on first line to prevent copying
	std::vector<Volume<t2Vector<Complex>>> out(mdt.numberOfObjects());
	
	const int pc = mdt.numberOfObjects();
	
	#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
	{
		out[p] = predictComplexGradient(mdt, p, obs, hs, applyCtf, applyTilt, applyShift);
	}
	
	return out;
}

Volume<t2Vector<Complex>> ReferenceMap::predictComplexGradient(
		const MetaDataTable &mdt, 
		int p, ObservationModel &obs, 
		ReferenceMap::HalfSet hs, 
		bool applyCtf, bool applyTilt, bool applyShift)
{
	Volume<t2Vector<Complex>> pred;
	
	int randSubset;
	mdt.getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
	randSubset -= 1;
	
	int pi = (hs == Own)? randSubset : 1 - randSubset;
	
	pred = obs.predictComplexGradient(projectors[pi], mdt, p, angpix, applyCtf, applyTilt, applyShift);
	
	return pred;
}

std::vector<Image<Complex>> ReferenceMap::predictAll(
		const MetaDataTable& mdt,
		const LegacyObservationModel& obs,
		HalfSet hs, int threads,
		bool applyCtf, bool applyTilt, bool applyShift)
{
	// declare on first line to prevent copying
	std::vector<Image<Complex>> out(mdt.numberOfObjects());
	
	const int pc = mdt.numberOfObjects();
	
	#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
	{
		out[p] = predict(mdt, p, obs, hs, applyCtf, applyTilt, applyShift);
	}
	
	return out;
}

Image<Complex> ReferenceMap::predict(
		const MetaDataTable& mdt, int p,
		const LegacyObservationModel& obs,
		HalfSet hs,
		bool applyCtf, bool applyTilt, bool applyShift)
{
	Image<Complex> pred;
	
	int randSubset;
	mdt.getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
	randSubset -= 1;
	
	int pi = (hs == Own)? randSubset : 1 - randSubset;
	
	pred = obs.predictObservation(projectors[pi], mdt, p, applyCtf, applyTilt, applyShift);
	
	return pred;
}

double ReferenceMap::angToPix(double a) const
{
	return s * angpix / a;
}

double ReferenceMap::pixToAng(double p) const
{
	return s * angpix / p;
}

// perhaps some other day:
/*void ReferenceMap::predictOccupancy(const MetaDataTable &particles, int threads)
{
	for (int half = 0; half < 1; half++)
	{
		occupancies[half] = Projector(s, TRILINEAR, 1.0, 10, 2);
		occupancies[half].data = MultidimArray<Complex>(1,s,s,sh);
		
		const int pc
		std::vector<CTF> ctfs(
				
		#pragma omp parallel for num_threads(threads)
		for (int z = 0; z < s; z++)
		{
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				
			}
		}
	}
}*/
