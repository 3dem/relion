#include <src/args.h>
#include <src/jaz/tomography/dynamo/catalogue.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/symmetry.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/util/zio.h>
#include <iostream>


using namespace gravis;

int main(int argc, char *argv[])
{
	IOParser parser;
	
	std::string particlesFn, tomoSetFn, outFn;
	int particlesToOutput;
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		particlesFn = parser.getOption("--i", "Input particle set");
		tomoSetFn = parser.getOption("--t", "Tomogram set", "tomograms.star");
		particlesToOutput = textToInteger(parser.getOption("--n", "Number of particles to output", "10"));
		outFn = parser.getOption("--o", "Output filename pattern");
		
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	TomogramSet tomogramSet(tomoSetFn);
	
	ParticleSet dataSet(particlesFn, "");
	std::vector<std::vector<ParticleIndex>> particles = dataSet.splitByTomogram(tomogramSet);
		
	const int tc = particles.size();
		
	std::ofstream output(outFn);
	
	int tpc = 0;
	
	for (int t = 0; t < tc; t++)
	{
		const int pc = particles[t].size();
		
		if (pc == 0) continue;

		Tomogram tomogram = tomogramSet.loadTomogram(t, false);
		output << "ts_" << t << ": \n";
		
		std::vector<d4Matrix> projTomo = tomogram.projectionMatrices;	
		const int fc = projTomo.size();
			
		for (int p = 0; p < pc; p++)
		{
			output << "\n\tparticle " << p << ": \n\n";
			
			const ParticleIndex part_id = particles[t][p];
			
			const d3Vector pos = dataSet.getPosition(part_id);
			
			const gravis::d4Vector pw(pos.x, pos.y, pos.z, 1.0);
			
			for (int f = 0; f < fc; f++)
			{
				const d4Vector q = projTomo[f] * pw;
				d2Vector r = gravis::d2Vector(q.x, q.y);
				
				output << "\t\t" << f << ": \t" 
					   << std::setw(8) << r.x << " \t" 
					   << std::setw(8) << r.y << '\n';
			}
			
			tpc++;
			
			if (tpc >= particlesToOutput) return 0;
		}
	}
	
	return 0;
}
