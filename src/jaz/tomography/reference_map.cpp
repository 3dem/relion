#include "reference_map.h"
#include <src/jaz/image/padding.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>

using namespace gravis;


TomoReferenceMap::TomoReferenceMap(){}

TomoReferenceMap::TomoReferenceMap(
		std::string ref1Fn, std::string ref2Fn, int boxSize, 
		std::string maskFn, std::string fscFn,
		bool rcThresh, double threshWidth)
{
	image_real.resize(2);
	image_real[0].read(ref1Fn);
	image_real[1].read(ref2Fn);
	
	if (!image_real[0].hasEqualSize(image_real[1]))
	{
		REPORT_ERROR_STR("ReferenceMap::ReferenceMap: reference maps have different sizes ("
				<< image_real[0].getSizeString() << " and " 
				<< image_real[1].getSizeString() << ") ");
	}
	
	if (maskFn != "")
	{
		mask.read(maskFn);
		
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
	
	if (image_real[0].xdim < boxSize)
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
		FFT::FourierTransform(image_real[i], image_FS[i], FFT::Both);
		Centering::shiftInSitu(image_FS[i]);
	}
	
	const int s = boxSize;
	const int sh = s/2 + 1;
	
	if (fscFn != "")
	{
		MetaDataTable fscMdt;
		fscMdt.read(fscFn, "fsc");
		
		if (!fscMdt.containsLabel(EMDL_SPECTRAL_IDX))
		{
			REPORT_ERROR("ReferenceMap::ReferenceMap: " + fscFn + " does not contain a value for "
						 + EMDL::label2Str(EMDL_SPECTRAL_IDX));
		}
		
		if (!fscMdt.containsLabel(EMDL_POSTPROCESS_FSC_TRUE))
		{
			REPORT_ERROR("ReferenceMap::ReferenceMap: " + fscFn + " does not contain a value for "
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
			
			if (rcThresh)
			{
				if (r < firstBad - threshWidth/2.0)
				{
					freqWeight(x,y) = 1.f;
				}
				else if (r < firstBad + threshWidth/2.0)
				{
					double t = (r - firstBad)/threshWidth + 0.5;
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
