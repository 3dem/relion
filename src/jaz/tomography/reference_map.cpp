#include "reference_map.h"
#include <src/jaz/image/padding.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>

using namespace gravis;


TomoReferenceMap::TomoReferenceMap(){}

void TomoReferenceMap::read(IOParser &parser)
{
	mapFilenames[0] = parser.getOption("--ref1", "Reference map, half 1");
	mapFilenames[1] = parser.getOption("--ref2", "Reference map, half 2");
	maskFilename = parser.getOption("--mask", "Reference mask", "");
	fscFilename = parser.getOption("--fsc", "Star file containing the FSC of the reference", "");

	useFscThreshold = !parser.checkOption("--fsc_act", "Use the actual FSC as the frq. weight");
	fscThresholdWidth = textToDouble(parser.getOption("--fsc_thresh_width", "Width of the frq. weight flank", "5"));
}

void TomoReferenceMap::load(int boxSize)
{
	image_real.resize(2);
	image_real[0].read(mapFilenames[0]);
	image_real[1].read(mapFilenames[1]);

	if (!image_real[0].hasEqualSize(image_real[1]))
	{
		REPORT_ERROR_STR("ReferenceMap::ReferenceMap: reference maps have different sizes ("
				<< image_real[0].getSizeString() << " and "
				<< image_real[1].getSizeString() << ") ");
	}

	if (maskFilename != "")
	{
		mask.read(maskFilename);

		if (!mask.hasEqualSize(image_real[0]))
		{
			REPORT_ERROR_STR("ReferenceMap::ReferenceMap: mask and reference have different sizes ("
							 << mask.getSizeString() << " and "
							 << image_real[0].getSizeString() << ") ");
		}

		for (int i = 0; i < 2; i++)
		{
			image_real[i] *= mask;
		}
	}

	const int s = boxSize < 0? image_real[0].xdim : boxSize;
	const int sh = s/2 + 1;

	if (image_real[0].xdim < s)
	{
		const int sr = image_real[0].xdim;

		Log::print("Padding reference from "+ZIO::itoa(sr)+" to "+ZIO::itoa(boxSize));

		for (int i = 0; i < 2; i++)
		{
			image_real[i] = Padding::padCenter3D_full(image_real[i], (boxSize - sr)/2);
		}
	}


	image_FS.resize(2);

	for (int i = 0; i < 2; i++)
	{
		presharpen(image_real[i], 1.0);

		FFT::FourierTransform(image_real[i], image_FS[i], FFT::Both);
		Centering::shiftInSitu(image_FS[i]);
	}

	if (fscFilename != "")
	{
		MetaDataTable fscMdt;
		fscMdt.read(fscFilename, "fsc");

		if (!fscMdt.containsLabel(EMDL_SPECTRAL_IDX))
		{
			REPORT_ERROR("ReferenceMap::ReferenceMap: " + fscFilename + " does not contain a value for "
						 + EMDL::label2Str(EMDL_SPECTRAL_IDX));
		}

		if (!fscMdt.containsLabel(EMDL_POSTPROCESS_FSC_TRUE))
		{
			REPORT_ERROR("ReferenceMap::ReferenceMap: " + fscFilename + " does not contain a value for "
						 + EMDL::label2Str(EMDL_POSTPROCESS_FSC_TRUE));
		}

		const int sh_fsc = fscMdt.numberOfObjects();

		std::vector<double> fsc(sh_fsc);

		int firstBad = sh_fsc;

		for (int i = 0; i < sh_fsc; i++)
		{
			int idx;
			fscMdt.getValueSafely(EMDL_SPECTRAL_IDX, idx, i);
			fscMdt.getValueSafely(EMDL_POSTPROCESS_FSC_TRUE, fsc[i], i);

			if (fsc[i] < 0.143)
			{
				fsc[i] = 0.0;

				if (i < firstBad)
				{
					firstBad = i;
				}
			}
		}

		double scale = sh_fsc / (double) sh;

		freqWeight = BufferedImage<float>(sh,s);

		for (int y = 0; y < s; y++)
		for (int x = 0; x < sh; x++)
		{
			double xx = x;
			double yy = y < s/2? y : y - s;

			double r = sqrt(xx*xx + yy*yy) * sh_fsc / (double) sh;

			if (useFscThreshold)
			{
				if (r < firstBad - fscThresholdWidth/2.0)
				{
					freqWeight(x,y) = 1.f;
				}
				else if (r < firstBad + fscThresholdWidth/2.0)
				{
					double t = (r - firstBad)/fscThresholdWidth + 0.5;
					freqWeight(x,y) = 0.5 * (cos(PI*t) + 1);
				}
				else
				{
					freqWeight(x,y) = 0.f;
				}
			}
			else
			{
				int ri = (int)(scale*r+0.5);
				if (ri >= sh_fsc) ri = sh_fsc-1;

				freqWeight(x,y) = fsc[ri];
			}
		}
	}
	else
	{
		freqWeight = BufferedImage<float>(sh,s);
		freqWeight.fill(1.f);
	}
}

int TomoReferenceMap::getBoxSize() const
{
	return image_real[0].xdim;
}

void TomoReferenceMap::presharpen(
		BufferedImage<float>& map_RS,
		double padding)
{
	const int s = image_real[0].xdim;

	const double epsilon = 0.02;

	for (int z = 0; z < s; z++)
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const double xx = x - s/2;
		const double yy = y - s/2;
		const double zz = z - s/2;

		const double r = sqrt(xx*xx + yy*yy + zz*zz);
		const double rval = r / (s * padding);

		if (rval > 0.0 && rval < 1.0)
		{
			const double sinc = sin(PI * rval) / (PI * rval);

			map_RS(x,y,z) /= sinc * sinc + epsilon;
		}
		else
		{
			map_RS(x,y,z) = 0.;
		}
	}
}
