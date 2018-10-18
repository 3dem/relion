#include <src/metadata_table.h>
#include <src/jaz/io/star_converter.h>
#include <src/jaz/obs_model.h>

int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		std::cerr << "usage: relion_merge_particles <input1> <input2> ... <inputN> <output>\n";
		return -1;
	}
	
	const int srcCount = argc - 2;
	std::string destFn = argv[argc-1];
	
	std::vector<MetaDataTable> particleMdts(srcCount), opticsMdts(srcCount);
	std::vector<ObservationModel> obsModels(srcCount);
	
	std::cout << "merging: " << std::endl;
	
	for (int i = 0; i < srcCount; i++)
	{
		std::string srcFn = argv[i+1];
		
		std::cout << "    " << srcFn << std::endl;
		ObservationModel::loadSafely(srcFn, obsModels[i], particleMdts[i], opticsMdts[i]);
	}
	
	std::cout << "into: " << destFn << std::endl;
	
	
	std::vector<std::vector<int>> optGrTransl(srcCount);
	
	MetaDataTable particleOut, opticsOut;
	
	for (int i = 0; i < srcCount; i++)
	{
		const int ogc = opticsMdts[i].numberOfObjects();
		optGrTransl[i].resize(ogc);
		
		for (int g = 0; g < ogc; g++)
		{
			opticsOut.addObject(opticsMdts[i].getObject(g));
			
			const int ogNew = opticsOut.numberOfObjects() - 1;
			
			opticsOut.setValue(EMDL_IMAGE_OPTICS_GROUP, ogNew+1, ogNew);
			
			optGrTransl[i][g] = ogNew;
		}
		
		const int pc = particleMdts[i].numberOfObjects();
		
		for (int p = 0; p < pc; p++)
		{
			particleOut.addObject(particleMdts[i].getObject(p));
			
			const int pNew = particleOut.numberOfObjects() - 1;
			
			int og0; 
			particleOut.getValue(EMDL_IMAGE_OPTICS_GROUP, og0, pNew);
			og0--;
			
			int og1 = optGrTransl[i][og0];			
			particleOut.setValue(EMDL_IMAGE_OPTICS_GROUP, og1+1, pNew);
		}
	}
	
	ObservationModel::save(particleOut, opticsOut, destFn);
	
	return 0;
}
