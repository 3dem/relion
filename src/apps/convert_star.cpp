#include <src/metadata_table.h>

int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		std::cerr << "usage: relion_convert_star <input> <output>\n";
		std::cerr << "       will produce <output>.star and <output>_optics.star\n";
	}
	
	MetaDataTable mdt;
	mdt.read(argv[1]);
	
	int ver = mdt.getVersion();
	int curVer = MetaDataTable::getCurrentVersion();
	
	if (ver == curVer)
	{
		std::cerr << "File " << argv[1] << " is already at version " << curVer/10000.0 << "\n";
		std::exit(-1);
	}
	else if (ver > curVer)
	{
		std::cerr << "The version of " << argv[1] << " is " << ver/10000.0
			<< " - this is beyond the current version of Relion (" << curVer/10000.0 << ")\n";
		std::cerr << "You are either using an outdated copy of Relion, or the file is from the future.\n";
		std::exit(-2);
	}
	
	MetaDataTable mdtOut, optOut;
	
	if (ver < 31000)
	{
		const int particleCount = mdt.numberOfObjects();
		
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
			if (mdt.labelExists(allOpticsLabels_double[l]))
			{
				opticsLabels_double.push_back(allOpticsLabels_double[l]);
			}
		}
		
		const int opticsLabelCount_double = opticsLabels_double.size();
				
		std::vector<std::vector<double>> groupValues_double(0);
		std::vector<int> opticsClasses(particleCount, -1);		
		
		for (int p = 0; p < particleCount; p++)
		{
			int foundGroup = -1;			 
			
			std::vector<double> curVals_double(opticsLabelCount_double);
			
			for (int l = 0; l < opticsLabelCount_double; l++)
			{
				mdt.getValue(opticsLabels_double[l], curVals_double[l], p);
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
		
		mdtOut = mdt;
		
		for (int l = 0; l < opticsLabelCount_double; l++)
		{
			mdtOut.deactivateLabel(opticsLabels_double[l]);
		}
		
		mdtOut.addLabel(EMDL_IMAGE_OPTICS_GROUP);
		
		for (int p = 0; p < particleCount; p++)
		{
			mdtOut.setValue(EMDL_IMAGE_OPTICS_GROUP, opticsClasses[p] + 1, p);
		}
		
		mdtOut.write(std::string(argv[2])+".star");
		
		
		optOut.addLabel(EMDL_IMAGE_OPTICS_GROUP);
				
		for (int l = 0; l < opticsLabelCount_double; l++)
		{
			optOut.addLabel(opticsLabels_double[l]);
		}
		
		for (int g = 0; g < groupValues_double.size(); g++)
		{
			optOut.addObject();
			optOut.setValue(EMDL_IMAGE_OPTICS_GROUP, g + 1, g);
			
			for (int l = 0; l < opticsLabelCount_double; l++)
			{
				optOut.setValue(opticsLabels_double[l], groupValues_double[g][l], g);
			}
		}
		
		optOut.write(std::string(argv[2])+"_optics.star");
	}
	
	return 0;
}
