#include "reference_map.h"
#include <src/jaz/obs_model.h>
#include <src/jaz/legacy_obs_model.h>
#include <src/jaz/image_op.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/motion/motion_helper.h>
#include <src/args.h>
#include <omp.h>

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
	maskFn = parser.getOption("--mask", "Reference mask", "");
	fscFn = parser.getOption("--f", "Input STAR file with the FSC of the reference");
	paddingFactor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));
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

Image<RFLOAT> ReferenceMap::getHollowWeight(double kmin_px)
{
	Image<RFLOAT> out = freqWeight;
	
	for (int y = 0; y < s; y++)
	for (int x = 0; x < sh; x++)
	{
		double xx = x;
		double yy = y <= sh? y : y - s;
		double r = sqrt(xx*xx + yy*yy);
		
		if (r < kmin_px)
		{
			out(y,x) = 0.0;
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
	
	pred = obs.predictObservation(projectors[pi], mdt, p, applyCtf, applyTilt, applyShift);
	
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
