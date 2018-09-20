#include "star_converter.h"

void StarConverter::convert_3p0_particlesTo_3p1(
		const MetaDataTable &in, MetaDataTable &outParticles, MetaDataTable &outOptics)
{
	int ver = in.getVersion();
	int curVer = MetaDataTable::getCurrentVersion();
	
	if (ver == curVer)
	{
		REPORT_ERROR_STR("StarConverter::convert_3p0_particlesTo_3p1: Star file is already at version " 
						 << curVer/10000.0);
	}
	else if (ver > curVer)
	{
		REPORT_ERROR_STR("StarConverter::convert_3p0_particlesTo_3p1: Star file is at version " 
			<< ver/10000.0 << " - this is beyond the current version of Relion (" 
			<< curVer/10000.0 << ")\n"
			<< "You are either using an outdated copy of Relion, or the file is from the future.\n");
	}
	
	const int outVer = 30001;

	const int particleCount = in.numberOfObjects();
	
	std::vector<EMDLabel> allOpticsLabels_double(0);
	
	allOpticsLabels_double.push_back(EMDL_CTF_Q0);		
	allOpticsLabels_double.push_back(EMDL_IMAGE_BEAMTILT_X);
	allOpticsLabels_double.push_back(EMDL_IMAGE_BEAMTILT_Y);
	allOpticsLabels_double.push_back(EMDL_CTF_CS);
	allOpticsLabels_double.push_back(EMDL_CTF_VOLTAGE);
	
	allOpticsLabels_double.push_back(EMDL_CTF_DETECTOR_PIXEL_SIZE);
	allOpticsLabels_double.push_back(EMDL_CTF_MAGNIFICATION);
	
	std::vector<EMDLabel> opticsLabels_double(0);
	
	for (int l = 0; l < allOpticsLabels_double.size(); l++)
	{
		if (in.labelExists(allOpticsLabels_double[l]))
		{
			opticsLabels_double.push_back(allOpticsLabels_double[l]);
		}
	}
	
	const int opticsLabelCount_double = opticsLabels_double.size();
			
	std::vector<std::vector<double>> groupValues_double(0);
	std::vector<int> opticsClasses(particleCount, -1);		
	
	for (long int p = 0; p < particleCount; p++)
	{
		int foundGroup = -1;			 
		
		std::vector<double> curVals_double(opticsLabelCount_double);
		
		for (int l = 0; l < opticsLabelCount_double; l++)
		{
			in.getValue(opticsLabels_double[l], curVals_double[l], p);
		}
		
		for (int g = 0; g < groupValues_double.size(); g++)
		{
			bool groupGood = true;
			
			for (int l = 0; l < opticsLabelCount_double; l++)
			{
				if (curVals_double[l] != groupValues_double[g][l])
				{
					groupGood = false;
					break;
				}
			}
			
			if (groupGood)
			{
				foundGroup = g;
				break;
			}
		}
		
		if (foundGroup >= 0)
		{
			opticsClasses[p] = foundGroup;
		}
		else
		{
			groupValues_double.push_back(curVals_double);				
			opticsClasses[p] = groupValues_double.size() - 1;
		}
	}
	
	outParticles = in;
	
	for (int l = 0; l < opticsLabelCount_double; l++)
	{
		outParticles.deactivateLabel(opticsLabels_double[l]);
	}
	
	outParticles.addLabel(EMDL_IMAGE_OPTICS_GROUP);
	
	for (long int p = 0; p < particleCount; p++)
	{
		outParticles.setValue(EMDL_IMAGE_OPTICS_GROUP, opticsClasses[p] + 1, p);
	}
	
	outParticles.setName("particles");
	outParticles.setVersion(outVer);
	
	
	outOptics.addLabel(EMDL_IMAGE_OPTICS_GROUP);
			
	for (int l = 0; l < opticsLabelCount_double; l++)
	{
		outOptics.addLabel(opticsLabels_double[l]);
	}
	
	for (int g = 0; g < groupValues_double.size(); g++)
	{
		outOptics.addObject();
		outOptics.setValue(EMDL_IMAGE_OPTICS_GROUP, g + 1, g);
		
		for (int l = 0; l < opticsLabelCount_double; l++)
		{
			outOptics.setValue(opticsLabels_double[l], groupValues_double[g][l], g);
		}
	}
	
	unifyPixelSize(outOptics);
	translateOffsets(outParticles, outOptics);
			
	outOptics.setName("optics");
	outOptics.setVersion(outVer);
}

void StarConverter::unifyPixelSize(MetaDataTable& outOptics)
{
	for (int i = 0; i < outOptics.numberOfObjects(); i++)
	{
		double dstep, mag;
		
		outOptics.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, i);
		outOptics.getValue(EMDL_CTF_MAGNIFICATION, mag, i);
		
		double angpix = 10000 * dstep / mag;
		
		outOptics.setValue(EMDL_MLMODEL_PIXEL_SIZE, angpix, i);
	}
	
	outOptics.deactivateLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE);
	outOptics.deactivateLabel(EMDL_CTF_MAGNIFICATION);
}

void StarConverter::translateOffsets(MetaDataTable &outParticles, const MetaDataTable &optics)
{
	for (int i = 0; i < outParticles.numberOfObjects(); i++)
	{
		int og;
		outParticles.getValue(EMDL_IMAGE_OPTICS_GROUP, og, i);
		og--;
		
		double angpix;
		optics.getValue(EMDL_MLMODEL_PIXEL_SIZE, angpix, og);
		
		double x, y, z;
		
		if (outParticles.containsLabel(EMDL_ORIENT_ORIGIN_X))
		{
			outParticles.getValue(EMDL_ORIENT_ORIGIN_X, x, i);
			outParticles.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, x*angpix, i);
		}
		
		if (outParticles.containsLabel(EMDL_ORIENT_ORIGIN_Y))
		{
			outParticles.getValue(EMDL_ORIENT_ORIGIN_Y, y, i);
			outParticles.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, y*angpix, i);
		}
		
		if (outParticles.containsLabel(EMDL_ORIENT_ORIGIN_Z))
		{
			outParticles.getValue(EMDL_ORIENT_ORIGIN_Z, z, i);
			outParticles.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, z*angpix, i);
		}
	}
	
	outParticles.deactivateLabel(EMDL_ORIENT_ORIGIN_X);
	outParticles.deactivateLabel(EMDL_ORIENT_ORIGIN_Y);
	outParticles.deactivateLabel(EMDL_ORIENT_ORIGIN_Z);
}
