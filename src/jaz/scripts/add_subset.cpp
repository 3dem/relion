#include <src/args.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/particle_set.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		std::cerr << "usage: relion_exp_add_subset <in> <out>" << std::endl;
		return -1;
	}
	
	ParticleSet ps(argv[1], "");
	
	std::cout << ps.partTable.numberOfObjects() << std::endl;
	
	ps.partTable.addLabel(EMDL_PARTICLE_RANDOM_SUBSET);
	
	for (int p = 0; p < ps.partTable.numberOfObjects(); p++)
	{
		ps.partTable.setValue(EMDL_PARTICLE_RANDOM_SUBSET, 1 + (rand()%2), p);
	}
	
	ps.write(argv[2]);
	        
	return 0;
	
}
