#include "substack.h"
#include <src/jaz/tomography/dynamo/catalogue.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/tomography/tomolist.h>
#include <iostream>

using namespace gravis;


void SubstackProgram::run()
{
	std::string stackOutFn = outTag + "/stacks/";
	std::string projOutFn = outTag + "/matrices/";
	std::string tomoListOutFn = outTag + "/tomolist.txt";
	std::string catOutFn = outTag + "/particles.cat";
	
	int dummy;
	dummy = std::system(("mkdir -p "+stackOutFn).c_str());
	dummy = std::system(("mkdir -p "+projOutFn).c_str());
	
	Catalogue cat(particlesFn);	
	std::vector<std::vector<DynamoParticle>> particles = cat.splitByTomogram();
		
	const int tc = particles.size();
	const int s = boxSize;

    const int s0 = (int)(binning * s + 0.5);
	
	TomoList tomoList(tomoListFn);
	
	TomoList tomoListOut;
	Catalogue catOut;
	
	for (int t = 0; t < tc; t++)
	{
		const int pc = particles[t].size();
		
		if (pc == 0 || ! tomoList.isKnown(t)) continue;
						
		std::string stackFn = tomoList.getTiltSeriesFilename(t);
		std::string projFn = tomoList.getProjectionsFilename(t);
				
		BufferedImage<float> stack;
		stack.read(stackFn);
		
		const int fc = stack.zdim;
		
		int w0, h0, d0;
		std::vector<d4Matrix> projTomo = ProjectionIO::read(projFn, w0, h0, d0);
		
		if (fc != stack.zdim)
		{
			REPORT_ERROR_STR(stackFn << " contains " << stack.zdim << " frames, while " << projFn
							 << " contains " << projTomo.size());
		}
		
		std::cout << stackFn << ":" << std::endl;
		
		for (int p = 0; p < pc; p++)
		{
			if (p%10 == 0 || p < 10)
			{
				std::cout << "    " << p << "/" << pc << std::endl;
			}
			
			DynamoParticle& part = particles[t][p];
			
			const d3Vector pos = part.getPosition();
			std::vector<d4Matrix> projCut(fc);
			
			BufferedImage<float> particleStackReal(s,s,fc);
			
			TomoExtraction::extractAt3D_real(
                stack, s0, binning, projTomo, pos,
				particleStackReal, projCut, num_threads, true);
			
			std::string partName = part.getFormattedTag();
			
			std::string tsOut = stackOutFn + partName + ".mrc";
			std::string projOut = projOutFn + partName + ".proj";
			
			particleStackReal.write(tsOut);			
			ProjectionIO::write(projCut, w0, h0, d0, projOut);
			
			int tomo = tomoListOut.addTomogram(tsOut, projOut, tomoList.getCtfFilename(t));
			
			DynamoParticle p2(part);
			p2.tomo = tomo;
			catOut.particles.push_back(p2);
			
		}
	}
	
	catOut.write(catOutFn);
	tomoListOut.write(tomoListOutFn);	
}
