#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/metadata_table.h>

using namespace gravis;


int main(int argc, char *argv[])
{
	IOParser parser;
	
	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	
	std::string sourceFn = parser.getOption("--i", "Input STAR file containing the source particles");
	std::string refFn = parser.getOption("--i_ref", "Input STAR file containing reference particles");
	
	const bool copyAngles = parser.checkOption("--angles", "Copy particle viewing angles from reference");
	const bool copyOffsets = parser.checkOption("--offsets", "Copy particle offsets from reference");
	
	std::string outFn = parser.getOption("--o", "Output path", "selected.star");
	
	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	MetaDataTable sourceAll, refAll;
	
	sourceAll.read(sourceFn);
	refAll.read(refFn);
	
	std::vector<MetaDataTable> sourceByMic = StackHelper::splitByMicrographName(sourceAll);
	std::vector<MetaDataTable> refByMic = StackHelper::splitByMicrographName(refAll);
	
	std::map<std::string, MetaDataTable*> micToSource;
	
	for (int m = 0; m < sourceByMic.size(); m++)
	{
		std::string micName;
		sourceByMic[m].getValue(EMDL_MICROGRAPH_NAME, micName, 0);
		
		micToSource[micName] = &sourceByMic[m];
	}
	
	MetaDataTable out;
	
	for (int m = 0; m < refByMic.size(); m++)
	{
		std::string micName;
		refByMic[m].getValue(EMDL_MICROGRAPH_NAME, micName, 0);
		
		if (micToSource.find(micName) == micToSource.end())
		{
			std::cerr << "Warning: " << micName << " not found.\n";
			continue;
		}
		
		MetaDataTable* src = micToSource[micName];
		
		const int pcRef = refByMic[m].numberOfObjects();
		const int pcSrc = src->numberOfObjects();
		
		std::vector<d2Vector> posSrc(pcSrc);
		
		for (int p = 0; p < pcSrc; p++)
		{
			src->getValue(EMDL_IMAGE_COORD_X, posSrc[p].x, p);
			src->getValue(EMDL_IMAGE_COORD_Y, posSrc[p].y, p);
		}
		
		std::vector<d2Vector> posRef(pcRef);
		
		for (int p = 0; p < pcRef; p++)
		{
			refByMic[m].getValue(EMDL_IMAGE_COORD_X, posRef[p].x, p);
			refByMic[m].getValue(EMDL_IMAGE_COORD_Y, posRef[p].y, p);
		}
		
		int missing = 0, multiple = 0;
		
		for (int p = 0; p < pcRef; p++)
		{
			int qBest = -1;
			
			for (int q = 0; q < pcSrc; q++)
			{
				double dist = (posRef[p] - posSrc[q]).length();
				
				if (dist < 1.0)
				{
					if (qBest == -1)
					{
						qBest = q;
					}
					else
					{
						qBest = -2;
					}
				}
			}
			
			if (qBest >= 0)
			{
				out.addObject(src->getObject(qBest));
				const int qNew = out.numberOfObjects() - 1;
				
				int randSubsetSrc;
				src->getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubsetSrc, qBest);
				
				int randSubsetRef;
				refByMic[m].getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubsetRef, p);
				
				if (randSubsetSrc != randSubsetRef)
				{
					if (copyAngles && copyOffsets)
					{
						out.setValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubsetRef, qNew);
					}
					else if (copyAngles != copyOffsets)
					{
						REPORT_ERROR_STR("Unable to copy only angles or only offsets, since the "
										 << "particles belong to different random subsets.\n");
					}
				}
				
				if (copyAngles)
				{
					double rot, tilt, psi;
					
					refByMic[m].getValue(EMDL_ORIENT_ROT, rot, p);
					refByMic[m].getValue(EMDL_ORIENT_TILT, tilt, p);
					refByMic[m].getValue(EMDL_ORIENT_PSI, psi, p);
					
					out.setValue(EMDL_ORIENT_ROT, rot, qNew);
					out.setValue(EMDL_ORIENT_TILT, tilt, qNew);
					out.setValue(EMDL_ORIENT_PSI, psi, qNew);
				}
				
				if (copyOffsets)
				{
					double xoff, yoff;
					
					refByMic[m].getValue(EMDL_ORIENT_ORIGIN_X, xoff, p);
					refByMic[m].getValue(EMDL_ORIENT_ORIGIN_Y, yoff, p);
					
					out.setValue(EMDL_ORIENT_ORIGIN_X, xoff, qNew);
					out.setValue(EMDL_ORIENT_ORIGIN_Y, yoff, qNew);
				}
			}			
			else if (qBest == -1)
			{
				missing++;
			}
			else // -2
			{
				multiple++;
			}
		}
		
		if (missing > 0)
		{
			std::cerr << "    Warning: " << missing << " of " << pcRef
					  << " particles missing from micrograph " << m << "\n";
		}
		
		if (multiple > 0)
		{
			std::cerr << "    Warning: " << multiple << " out of " << pcRef 
					  << " particles found multiple times in micrograph " << m << "\n"
					  << "    (all will be ignored)\n";
		}
	}
	
	out.write(outFn);
	
		
	return RELION_EXIT_SUCCESS;	
}
