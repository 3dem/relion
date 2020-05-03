#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/dynamo/catalogue.h>
#include <src/jaz/util/zio.h>
#include <src/args.h>
#include <set>

using namespace gravis;

int main(int argc, char *argv[])
{
	
	IOParser parser;
	
	std::string outFn, inFn;
	double binning;
	bool diag;
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		inFn = parser.getOption("--in", "input file");
		
		binning = textToDouble(parser.getOption("--binning", 
			"Binning level at which the coordinates were picked", "8"));
		
		outFn = parser.getOption("--o", "Output filename pattern");
		
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	std::vector<std::vector<double>> data = ZIO::readDoublesTable(inFn, 7);
	
	const int pc = data.size();
	
	std::vector<Catalogue> cats(3);
	
	for (int i = 0; i < 3; i++)
	{
		cats[i].particles.resize(pc);
	}
	
	const int tomo_col = 0;
	const int centre_col = 1;
	const int north_col = 4;
	
	const int north_index = 0;
	const int centre_index = 1;
	const int south_index = 2;
	
	std::set<int> usedTomos;
	
	for (int p = 0; p < pc; p++)
	{
		std::vector<DynamoParticle> particles(3);
		DynamoParticle& part_north = particles[north_index];
		DynamoParticle& part_centre = particles[centre_index];
		DynamoParticle& part_south = particles[south_index];
		
		int tomo = (int) (data[p][tomo_col] + 0.5);
		
		d3Vector centre(data[p][centre_col], data[p][centre_col+1], data[p][centre_col+2]);
		d3Vector north(data[p][north_col], data[p][north_col+1], data[p][north_col+2]);
		
		centre *= binning;
		north *= binning;
		
		d3Vector mid = 0.5 * (centre + north);
		
		part_north.x = north.x;
		part_north.y = north.y;
		part_north.z = north.z;
		
		part_centre.x = mid.x;
		part_centre.y = mid.y;
		part_centre.z = mid.z;
		
		part_south.x = centre.x;
		part_south.y = centre.y;
		part_south.z = centre.z;

		for (int i = 0; i < 3; i++)
		{
			DynamoParticle& pp = particles[i];
			pp.tag = p;
			pp.tomo = tomo;
			
			pp.setOrientation(north - centre, true);
			
			if (i == 0)
			{
				std::cout << pp.getAlignmentMatrixAlibi4x4(0,0,0) << "\n";
				std::cout << (north - centre).normalize() << "\n\n";
			}
			
			cats[i].particles[p] = pp;
		}
		
		usedTomos.insert(tomo);
	}
	
	cats[north_index].write(outFn + "_north.cat");
	cats[centre_index].write(outFn + "_centre.cat");
	cats[south_index].write(outFn + "_south.cat");
	
	std::ofstream tomoListFile(outFn + "_tomolist.txt");


	for (std::set<int>::iterator it = usedTomos.begin();
		 it != usedTomos.end(); it++)
	{
		int tomo = *it;
        tomoListFile << tomo << " tilt_series/R1_" << tomo << "_clear.mrc "
					 << "proj/R1_" << tomo << ".proj "
					 << "no-CTF\n";
	}
	
	return 0;
}
