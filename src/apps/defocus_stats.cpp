
#include <src/args.h>
#include <src/metadata_table.h>
#include <src/jaz/single_particle/stack_helper.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	IOParser parser;
	std::string starFn;

	parser.setCommandLine(argc, argv);
	parser.addSection("General options");
	starFn = parser.getOption("--i", "Input STAR file with a list of particles");
	
	if (parser.checkForErrors()) return RELION_EXIT_FAILURE;

	MetaDataTable mdt0;
	mdt0.read(starFn);
	
	double mu(0.0), var(0.0);
	
	for (int i = 0; i < mdt0.numberOfObjects(); i++)
	{
		double u, v;
		
		mdt0.getValue(EMDL_CTF_DEFOCUSU, u, i);
		mdt0.getValue(EMDL_CTF_DEFOCUSU, v, i);
		
		double a = 0.5 * (u + v);
		
		mu += a;
	}
	
	mu /= (double) mdt0.numberOfObjects();
	
	std::cout << "mu: " << mu << "\n";
	
	for (int i = 0; i < mdt0.numberOfObjects(); i++)
	{
		double u, v;
		
		mdt0.getValue(EMDL_CTF_DEFOCUSU, u, i);
		mdt0.getValue(EMDL_CTF_DEFOCUSU, v, i);
		
		double a = 0.5 * (u + v);
		double d = a - mu;
		
		var += d*d;
	}
	
	var /= (double) (mdt0.numberOfObjects() - 1);
	
	std::cout << "sigma: " << sqrt(var) << "\n";
	
	return RELION_EXIT_SUCCESS;
}
