
#include <src/args.h>
#include <src/metadata_table.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/motion/motion_helper.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	IOParser parser;
	std::string starFn, tracksPath0, tracksPath1;

	parser.setCommandLine(argc, argv);
	parser.addSection("General options");
	starFn = parser.getOption("--i", "Input STAR file with a list of particles");
	tracksPath0 = parser.getOption("--t0", "Path to particle trajectories 0");
	tracksPath1 = parser.getOption("--t1", "Path to particle trajectories 1");

	if (parser.checkForErrors()) return 1;

	MetaDataTable mdt0;
	mdt0.read(starFn);

	std::vector<MetaDataTable> mdts = StackHelper::splitByStack(&mdt0);

	const int mgc = mdts.size();

	double mu = 0.0;
	long tpc = 0;
	int fc = 0;
	
	int bc = 10;
	std::vector<double> bins(bc+1, 0.0);

	for (int m = 0; m < mgc; m++)
	{
		std::stringstream stsg;
		stsg << m;

		if (m%100 == 0) std::cout << "micrograph " << (m+1) << " / " << mgc << "\n";

		const int pc = mdts[m].numberOfObjects();
		tpc += pc;

		std::string tag;
		mdts[m].getValue(EMDL_IMAGE_NAME, tag, 0);
		tag = tag.substr(0,tag.find_last_of('.'));
		tag = tag.substr(tag.find_first_of('@')+1);
		tag = tag.substr(tag.find_last_of('/')+1);
		
		std::string tfn0 = tracksPath0 + "/" + tag + "_tracks.star";
		std::string tfn1 = tracksPath1 + "/" + tag + "_tracks.star";
		
		std::vector<std::vector<d2Vector>> shift0;
		std::vector<std::vector<d2Vector>> shift1;

		try
		{
			shift0 = MotionHelper::readTracks(tfn0);
			shift1 = MotionHelper::readTracks(tfn1);
		}
		catch (RelionError XE)
		{
			std::cerr << "Warning: error reading tracks in " << tfn0 << ", " << tfn1 << "\n";
			continue;
		}
		
		fc = shift0[0].size();
		
		for (int p = 0; p < pc; p++)
		for (int f = 0; f < fc; f++)
		{
			double d = (shift1[p][f] - shift0[p][f]).length(); 
			
			mu += d;
			
			int b = floor((log(d) / log(10.0)) + bc);
			
			if (b > bc) b = bc;
			else if (b < 0) b = 0;
			
			bins[b] += 1.0;
		}
	}
	
	mu /= tpc;
	
	std::cout << "mean: " << mu << "\n\n";
		
	for (int b = 0; b < bc+1; b++)
	{
		// b = (log(d) / log(10.0)) + bc;
		// (b - bc)*log(10.0) = log(d)
		// exp((b-bc)*log(10)) = d
		double bv = exp((b - bc) * log(10.0));
		
		std::cout << bv << " " << bins[b]/(tpc*fc) << "\n";
	}
	
	std::cout << "\n";
	
	double cumm = 0.0;
	
	for (int b = 0; b < bc+1; b++)
	{
		double bv = exp((b - bc) * log(10.0));			
		cumm += bins[b]/(tpc*fc);	
		std::cout << bv << " " << cumm << "\n";
	}
	
	return 0;
}
