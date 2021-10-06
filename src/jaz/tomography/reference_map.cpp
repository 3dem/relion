#include "reference_map.h"
#include <src/jaz/image/padding.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/image_file_helper.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/optimisation_set.h>

using namespace gravis;


TomoReferenceMap::TomoReferenceMap() : lastShell(0) {}

void TomoReferenceMap::read(IOParser &parser)
{
	mapFilenames[0] = parser.getOption("--ref1", "Reference map, half 1");
	mapFilenames[1] = parser.getOption("--ref2", "Reference map, half 2");
	maskFilename = parser.getOption("--mask", "Reference mask", "");
	fscFilename = parser.getOption("--fsc", "Star file containing the FSC of the reference", "");

	fscThresholdWidth = textToDouble(parser.getOption("--fsc_thresh_width", "Width of the frq. weight flank", "5"));
	freqCutoff_A =  textToDouble(parser.getOption("--freq_cutoff", "Explicit cutoff frequency (in Å; negative to turn off)", "-1"));
	flatWeight = !parser.checkOption("--use_SNR_weight", "Weight each shell proportionally to its reference-map confidence");
}

void TomoReferenceMap::read(const OptimisationSet &optimisationSet)
{
	mapFilenames[0] = optimisationSet.refMap1;
	mapFilenames[1] = optimisationSet.refMap2;
	maskFilename = optimisationSet.refMask;
	fscFilename = optimisationSet.refFSC;

	fscThresholdWidth = optimisationSet.fscThresholdWidth;
	freqCutoff_A = optimisationSet.freqCutoff_A;

	flatWeight = optimisationSet.flatWeight;
}

void TomoReferenceMap::load(int boxSize, int verbosity)
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

	pixelSize = ImageFileHelper::getSamplingRate(mapFilenames[0]);

	int manual_cutoff_px = -1;

	if (freqCutoff_A > 0)
	{
		manual_cutoff_px = s * pixelSize / freqCutoff_A;
	}

	if (image_real[0].xdim < s)
	{
		const int sr = image_real[0].xdim;

		if (verbosity > 0)
		{
			Log::print("Padding reference from "+ZIO::itoa(sr)+" to "+ZIO::itoa(boxSize));
		}

		for (int i = 0; i < 2; i++)
		{
			image_real[i] = Padding::padCenter3D_full(image_real[i], (boxSize - sr)/2);
		}
	}
	else if (image_real[0].xdim > s)
	{
		REPORT_ERROR_STR("Reference map is too big: " << image_real[0].xdim
				<< " pixels instead of " << s);
	}

	const double taper_edge_width = 10;

	for (int i = 0; i < 2; i++)
	{
		Reconstruction::taper(image_real[i], taper_edge_width, true, 1);
	}

	image_FS.resize(2);

	for (int i = 0; i < 2; i++)
	{
		presharpen(image_real[i], 1.0);

		FFT::FourierTransform(image_real[i], image_FS[i], FFT::Both);
		Centering::shiftInSitu(image_FS[i]);
	}


	SNR_weight = std::vector<double>(sh, 0.0);

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

		double scale = 2 * (sh_fsc - 1) / (double) s;

		for (int r = 0; r < sh; r++)
		{
			const double r_fsc = scale * r;
			const int r0 = (int) r_fsc;
			const int r1 = r0 + 1;

			if (r1 >= sh_fsc)
			{
				SNR_weight[r] = 0.0;
			}
			else
			{
				const double f = r_fsc - r0;
				const double fsc_r = (1 - f) * fsc[r0] + f * fsc[r1];

				SNR_weight[r] = 2.0 * fsc_r / (fsc_r + 1);
			}
		}

		lastShell = (firstBad - 1) / scale;
	}
	else
	{
		std::vector<double> sigma2 = std::vector<double>(sh,0.0);
		std::vector<double> tau2 = std::vector<double>(sh,0.0);

		std::vector<int> count = std::vector<int>(sh,0);

		for (int z = 0; z < s;  z++)
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			const double xx = x;
			const double yy = y < s/2? y : y - s;
			const double zz = z < s/2? z : z - s;

			const double rd = sqrt(xx*xx + yy*yy + zz*zz);

			const int r = (int)(rd + 0.5);

			if (r < sh)
			{
				const fComplex z0 = image_FS[0](x,y,z);
				const fComplex z1 = image_FS[1](x,y,z);

				sigma2[r] += (z1 - z0).norm() / 2;
				tau2[r] += (z0.norm() + z1.norm()) / 2;

				count[r]++;
			}
		}

		for (int r = 0; r < sh; r++)
		{
			if (count[r] > 0)
			{
				tau2[r] /= count[r];
				sigma2[r] /= count[r];
			}

			SNR_weight[r] = (tau2[r] > 0.0 && tau2[r] > sigma2[r])?
						(tau2[r] - sigma2[r]) / tau2[r] : 0.0;
		}

		lastShell = -1;
		const double threshold = 0.05;

		for (int r = 0; r < sh; r++)
		{
			if (SNR_weight[r] < threshold)
			{
				lastShell = r - 1;
				break;
			}
		}

		if (lastShell < 0) lastShell = sh - 1;
	}

	if (manual_cutoff_px > 0 && manual_cutoff_px < lastShell)
	{
		lastShell = manual_cutoff_px;
	}

	if (flatWeight)
	{
		for (int r = 0; r < sh; r++)
		{
			SNR_weight[r] = 1.0;
		}
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
	const int s = map_RS.xdim;

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
