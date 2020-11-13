#include "tomo_ctf_helper.h"
#include <src/jaz/optics/damage.h>
#include <src/jaz/util/zio.h>

using namespace gravis;

std::vector<CTF> TomoCtfHelper :: loadCtffind4(
		std::string path, int imageCount, 
		double voltage, double Cs, double Q0, double Bfac, double scale)
{
    /*
     example:
    # Output from CTFFind version 4.1.5, run on 2017-03-30 15:12:45
    # Input file: /beegfs/zivanov/tomograms/ts_05/frames/05_f32.mrc ; Number of micrographs: 1
    # Pixel size: 1.000 Angstroms ; acceleration voltage: 300.0 keV ; spherical aberration: 2.70 mm ; amplitude contrast: 0.07
    # Box size: 512 pixels ; min. res.: 30.0 Angstroms ; max. res.: 5.0 Angstroms ; min. def.: 5000.0 um; max. def. 50000.0 um
    # Columns: #1 - micrograph number; #2 - defocus 1 [Angstroms]; #3 - defocus 2; #4 - azimuth of astigmatism; #5 - additional phase shift [radians]; #6 - cross correlation; #7 - spacing (in Angstroms) up to which CTF rings were fit successfully
    1.000000 10295.926758 10012.275391 -38.856349 0.000000 0.030650 5.279412
    */
	
	std::vector<CTF> ctfs(imageCount);

	std::ifstream file(path);
	
	if (!file)
	{
		REPORT_ERROR_STR("CtfHelper::loadCtffind4: unable to read '" << path << "'.");
	}
		
	int currImg = 0;
	
	char text[4096];

	while (file.getline(text, 4096))
	{
		if (text[0] == '#') continue;
		
		std::stringstream line(text);

		ctfs[currImg] = setFromFile(line, voltage, Cs, Q0, Bfac, scale);
		currImg++;
		
		if (currImg >= imageCount)
		{
			break;
		}
	}
	
	if (currImg < imageCount)
	{
		REPORT_ERROR_STR("Insufficient number of CTFs found in " << path << ".\n"
						 << imageCount << " requested, " << currImg << " found.\n");
	}
    
	return ctfs;
}

std::vector<CTF> TomoCtfHelper :: loadCtfplotter(
		std::string path, int imageCount, 
		double voltage, double Cs, double Q0, double Bfac, double scale)
{
	std::vector<CTF> ctfs(0);
	ctfs.reserve(imageCount);

	std::ifstream file(path);
	
	if (!file)
	{
		REPORT_ERROR_STR("CtfHelper::loadCtfplotter: unable to read '" << path << "'.");
	}
		
	int currImg = 0;
	
	char text[4096];
	
	// get rid of the first line
	file.getline(text, 4096);

	while (file.getline(text, 4096))
	{
		if (text[0] == '#') continue;
		
		std::stringstream line(text);
		
		{
			double 
				unknown0, unknown1, unknown2, unknown3, 
				defocusU_nm, defocusV_nm, azimuth;
			
			line >> unknown0;
			line >> unknown1;
			line >> unknown2;
			line >> unknown3;
			line >> defocusU_nm;
			line >> defocusV_nm;
			line >> azimuth;
						 
			CTF ctf;
			
			ctf.setValues(10.0 * defocusU_nm, 10.0 * defocusV_nm, azimuth, 
						  voltage, Cs, Q0, Bfac, scale, 0.0);
			
			ctfs.push_back(ctf);
		}
		
		currImg++;
		
		if (currImg >= imageCount)
		{
			break;
		}
	}
	
	if (currImg < imageCount)
	{
		REPORT_ERROR_STR("Insufficient number of CTFs found in " << path << ".\n"
						 << imageCount << " requested, " << currImg << " found.\n");
	}
    
	return ctfs;
}

CTF TomoCtfHelper::setFromFile(std::stringstream& line,
						   double voltage, double Cs, double Q0, double Bfac, double scale)
{
	/*
	#1 - micrograph number;
	#2 - defocus 1 [Angstroms];
	#3 - defocus 2;
	#4 - azimuth of astigmatism;
	#5 - additional phase shift [radians];
	#6 - cross correlation;
	#7 - spacing (in Angstroms) up to which CTF rings were fit successfully
	*/

	double imgNumber, defocus1, defocus2, azimuth, phaseShift, crossCorr, bestBefore;

	line >> imgNumber;
	line >> defocus1;
	line >> defocus2;
	line >> azimuth;
	line >> phaseShift;
	line >> crossCorr;
	line >> bestBefore;

	CTF ctf;
	ctf.setValues(defocus1, defocus2, azimuth, voltage, Cs, Q0, Bfac, scale, phaseShift);
	
	return ctf;
}

void TomoCtfHelper::writeToFile(int frameIndex, const CTF& ctf, std::ofstream& file)
{
	file << frameIndex << "       " 
		 << ctf.DeltafU << "       " 
		 << ctf.DeltafV << "       " 
		 << ctf.azimuthal_angle << "       " 
		 << ctf.phase_shift << "       " 
		 << 0 << "       " 
		 << 0 << '\n';
}

void TomoCtfHelper::read(
		std::string doseFn, std::string ctfFn, 
		int s, int fc,
		double voltage, double Cs, double pixelSize,
		std::string stackFn, 
		std::vector<CTF>& centreCTFs_out, 
		BufferedImage<float>& doseWeights_out,
		std::vector<double>& cumulativeDose_out)
{
	bool do_dose_weighting = doseFn != "";
	
	centreCTFs_out = loadCtffind4(ctfFn, fc, voltage, Cs);
	
	if (fc > 0 && centreCTFs_out.size() != fc)
	{
		REPORT_ERROR_STR(ctfFn << " contains " << centreCTFs_out.size() << " frames, while " 
						 << stackFn << " contains " << fc);
	}
	
	fc = centreCTFs_out.size();
	
	if (do_dose_weighting)
	{
		if (ZIO::endsWith(doseFn, ".star"))
		{
			MetaDataTable mdt;
			mdt.read(doseFn);
			
			const int oc = mdt.numberOfObjects();
			
			if (oc != fc)
			{
				REPORT_ERROR_STR(doseFn << " contains " << oc << " frames, while " 
								 << stackFn << " contains " << fc);
			}
			
			cumulativeDose_out.resize(fc);
			
			std::vector<double> bfacs(fc), afacs(fc);
			
			for (int f = 0; f < fc; f++)
			{
				mdt.getValueSafely(EMDL_MICROGRAPH_PRE_EXPOSURE, cumulativeDose_out[f], f);
			}
			
			if (mdt.containsLabel(EMDL_CTF_BFACTOR) && mdt.containsLabel(EMDL_CTF_SCALEFACTOR))
			{
				for (int f = 0; f < fc; f++)
				{
					mdt.getValueSafely(EMDL_CTF_BFACTOR, bfacs[f], f);
					mdt.getValueSafely(EMDL_CTF_SCALEFACTOR, afacs[f], f);
				}
				
				doseWeights_out = Damage::computeBfactorWeights(bfacs, afacs, s, pixelSize);
			}
			else
			{
				doseWeights_out = Damage::weightStack_GG(cumulativeDose_out, pixelSize, s);
			}
			
		}
		else
		{
			cumulativeDose_out = ZIO::readDoubles(doseFn);
		
			if (cumulativeDose_out.size() != fc)
			{
				REPORT_ERROR_STR(doseFn << " contains " << cumulativeDose_out.size() << " frames, while " 
								 << stackFn << " contains " << fc);
			}
			
			doseWeights_out = Damage::weightStack_GG(cumulativeDose_out, pixelSize, s);
		}
		
	}
	else
	{
		doseWeights_out.fill(1.f);
	}
}
