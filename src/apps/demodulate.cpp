#include <src/jaz/vtk_helper.h>
#include <src/jaz/obs_model.h>
#include <src/jaz/stack_helper.h>
#include <src/metadata_table.h>
#include <src/args.h>


int main(int argc, char *argv[])
{
	std::string particlesFn, outPath;
	MetaDataTable particlesMdt, opticsMdt;
	int nr_omp_threads;
	bool r31;
		
	IOParser parser;

    try
    {
        parser.setCommandLine(argc, argv);
		
        parser.addSection("General options");
		
		particlesFn = parser.getOption("--i", "Input STAR file with a list of particles");
		nr_omp_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		r31 = parser.checkOption("--r31", "Write output in Relion-3.1 format");
		outPath = parser.getOption("--out", "Output path");
		
		if (parser.checkForErrors()) return 1;		
	}
	catch (RelionError XE)
    {
        parser.writeUsage(std::cout);
        std::cerr << XE;
        exit(1);
    }
	
	if (outPath[outPath.length()-1] != '/')
	{
		outPath += "/";
	}
	
	std::string command = " mkdir -p " + outPath;
	int res = system(command.c_str());
	
	ObservationModel obsModel;
	
	ObservationModel::loadSafely(
		particlesFn, obsModel, particlesMdt, opticsMdt);
		
	
	particlesMdt.read(particlesFn);
	
	std::vector<ParFourierTransformer> fts(nr_omp_threads);
			
	std::vector<MetaDataTable> mdts = StackHelper::splitByStack(&particlesMdt);
	
	const int mc = mdts.size();
	
	for (int m = 0; m < mc; m++)
	{
		const int pc = mdts[m].size();		
		
		std::vector<Image<Complex>> obs;		
		obs = StackHelper::loadStackFS(&mdts[m], "", nr_omp_threads, &fts, false);
				
		std::string name, fullName;
		mdts[m].getValue(EMDL_IMAGE_NAME, fullName, 0);
		name = fullName.substr(fullName.find("@")+1);
		
		for (int p = 0; p < pc; p++)
		{
			obsModel.demodulatePhase(mdts[m], p, obs[p].data);
		}
		
		std::vector<Image<RFLOAT>> demodulated = StackHelper::inverseFourierTransform(obs);
		Image<RFLOAT> out = StackHelper::toSingleImage(demodulated);
		
		FileName fn_pre, fn_jobnr, fn_post;
		decomposePipelineFileName(name, fn_pre, fn_jobnr, fn_post);
		
		std::string outFn = outPath + fn_post;
		
		if (outFn.find_last_of("/") != std::string::npos)
		{
			std::string command = " mkdir -p " + outFn.substr(0, outFn.find_last_of("/"));
			int res = system(command.c_str());
		}
		
		for (int p = 0; p < pc; p++)
		{
			std::stringstream sts;
			sts << (p+1) << "@" << outFn;
			
			mdts[m].setValue(EMDL_IMAGE_NAME, sts.str(), p);
		}
		
		out.write(outFn);
	}
	
	MetaDataTable mdt1;
	
	for (int m = 0; m < mc; m++)
	{
		mdt1.append(mdts[m]);
	}
	
	if (!r31)
	{
		const int tpc = mdt1.numberOfObjects();
		
		std::vector<EMDLabel> allOpticsLabels_double(0);
		
		allOpticsLabels_double.push_back(EMDL_CTF_Q0);
		allOpticsLabels_double.push_back(EMDL_CTF_CS);
		allOpticsLabels_double.push_back(EMDL_CTF_VOLTAGE);
		allOpticsLabels_double.push_back(EMDL_CTF_DETECTOR_PIXEL_SIZE);
		allOpticsLabels_double.push_back(EMDL_CTF_MAGNIFICATION);
		
		for (int l = 0; l < allOpticsLabels_double.size(); l++)
		{
			EMDLabel lab = allOpticsLabels_double[l];
			
			mdt1.addLabel(lab);
			
			for (int p = 0; p < tpc; p++)
			{
				int opticsGroup;
				mdt1.getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, p);
				opticsGroup--;
				
				double v;
				opticsMdt.getValue(lab, v, opticsGroup);
				mdt1.setValue(lab, v, p);
			}
		}
		
		opticsMdt.deactivateLabel(EMDL_IMAGE_OPTICS_GROUP);
		
		mdt1.setVersion(30000);
	}
	else
	{
		opticsMdt.deactivateLabel(EMDL_IMAGE_BEAMTILT_X);
		opticsMdt.deactivateLabel(EMDL_IMAGE_BEAMTILT_Y);		
		opticsMdt.deactivateLabel(EMDL_IMAGE_ODD_ZERNIKE_COEFFS);
		
		opticsMdt.write(outPath+"demodulated_particles_optics.star");
	}
	
	mdt1.write(outPath+"demodulated_particles.star");
	
	std::cout << "output written into " << (outPath+"demodulated_particles.star") << "\n";
			
	return 0;
}
