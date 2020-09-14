#include <src/args.h>
#include <src/jaz/util/zio.h>
#include <src/metadata_table.h>

int main(int argc, char *argv[])
{
	IOParser parser;
	
	std::string orderFn, tltFn, outFn;
	double fractionalDose;
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		orderFn = parser.getOption("--ol", "Frame-order list (*.csv)", "order_list.csv");
		tltFn = parser.getOption("--t", "Tilt angles list (*.tlt)");
		fractionalDose = textToDouble(parser.getOption("--fd", "Fractional dose"));
		outFn = parser.getOption("--o", "Output file");
		
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	std::vector<std::vector<double>> order = ZIO::readFixedDoublesTable(orderFn, 2, ',');
	std::vector<double> tilts = ZIO::readDoubles(tltFn);
	
	int fc = tilts.size();
	//
	
	MetaDataTable mdt;
	mdt.addLabel(EMDL_MICROGRAPH_FRAME_NUMBER);
	mdt.addLabel(EMDL_MICROGRAPH_PRE_EXPOSURE);
		
	for (int i = 0; i < fc; i++)
	{
		double minDist = 1000000.0;
		int bestJ = -1;
		
		for (int j = 0; j < order.size(); j++)
		{
			double d = order[j][1] - tilts[i];
			double dd = d*d;
			
			if (dd < minDist)
			{
				bestJ = j;
				minDist = dd;
			}
		}
		
		mdt.addObject();
		mdt.setValue(EMDL_MICROGRAPH_FRAME_NUMBER, i+1, i);
		mdt.setValue(EMDL_MICROGRAPH_PRE_EXPOSURE, bestJ * fractionalDose, i);
	}
	
	mdt.write(outFn);
	
	return 0;
	
	/*std::ofstream out(outFn);
	
	for (int i = 0; i < fc; i++)
	{
		double minDist = 1000000.0;
		int bestJ = -1;
		
		for (int j = 0; j < order.size(); j++)
		{
			double d = order[j][1] - tilts[i];
			double dd = d*d;
			
			if (dd < minDist)
			{
				bestJ = j;
				minDist = dd;
			}
		}
		
		out << (bestJ * fractionalDose) << '\n';
	}
	
	return 0;*/
}
