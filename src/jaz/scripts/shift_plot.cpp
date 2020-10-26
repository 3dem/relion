#include <src/ctf.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/motion/proto_alignment.h>
#include <src/jaz/tomography/projection_IO.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <omp.h>


using namespace gravis;

int main(int argc, char *argv[])
{
	std::string partFn0, partFn1, tomoSetFn0, tomoSetFn1, outDir;
	double scale;
			
	IOParser parser;
	parser.setCommandLine(argc, argv);
	
	try
	{
		int defocus_section = parser.addSection("General options");
		
		partFn0 = parser.getOption("--i0", "Initial-state particles");
		tomoSetFn0 = parser.getOption("--t0", "Initial tomogram set");
		
		partFn1 = parser.getOption("--i1", "Final-state particles");
		tomoSetFn1 = parser.getOption("--t1", "Final tomogram set");
		
		scale = textToDouble(parser.getOption("--scale", "Scale of output 3D arrows", "1"));
		
		outDir = parser.getOption("--o", "Output directory");
		
		if (parser.checkForErrors()) std::exit(-1);
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	if (outDir[outDir.length()-1] != '/')
	{
		outDir = outDir + "/";
	}
	
	int res = system(("mkdir -p "+outDir).c_str());
	
	
	ParticleSet dataSet0(partFn0, "");
	ParticleSet dataSet1(partFn1, "");

	TomogramSet tomogramSet0(tomoSetFn0);
	TomogramSet tomogramSet1(tomoSetFn1);

	std::vector<std::vector<ParticleIndex>> particles0 = dataSet0.splitByTomogram(tomogramSet0);
	std::vector<std::vector<ParticleIndex>> particles1 = dataSet1.splitByTomogram(tomogramSet1);
	
	const int tc0 = particles0.size();
	const int tc1 = particles1.size();
	
	if (tc0 != tc1)
	{
		REPORT_ERROR_STR("Unequal tomogram numbers: " << tc0 << " and " << tc1);
	}
	
	const int tc = tc0;
		
	for (int t = 0; t < tc; t++)
	{
		const int pc0 = particles0[t].size();
		const int pc1 = particles1[t].size();
		
		if (pc0 != pc1)
		{
			std::cerr << "Warning: unequal particle numbers in tomogram " << t << ": "
				  << pc0 << " and " << pc1 << std::endl;
			std::cerr << "skipping tomogram" << std::endl;
			
			continue;
		}
					 
		const int pc = pc0;
		
		if (pc == 0) continue;
		
		const std::string tag = outDir+"tomo_"+ZIO::itoa(t);
		
		Mesh out;
		
		double var1D = 0.0;
		d3Vector var3D(0.0, 0.0, 0.0);
	
		
		for (int pp = 0; pp < pc; pp++)
		{
			const ParticleIndex p = particles0[t][pp];
					
			d3Vector p0 = dataSet0.getPosition(p);
			d3Vector p1 = dataSet1.getPosition(p);
			d3Vector m = (p0 + p1) / 2.0;
			
			p0 = m + scale * (p0 - m);
			p1 = m + scale * (p1 - m);
			
			double rad = 0.1 * (p0 - p1).length();
			MeshBuilder::addCone(p0, p1, rad, 5, out);
			
			const d3Vector d = p1 - p0;
			
			var1D += d.norm2();
			
			var3D.x += d.x * d.x;
			var3D.y += d.y * d.y;
			var3D.z += d.z * d.z;
		}
		
		var1D /= pc - 1;
		var3D /= pc - 1;
		
		out.writePly(tag+"_mesh.ply");
		
		std::cout << "std. dev.: " << sqrt(var1D) << '\n';
		std::cout << "        x: " << sqrt(var3D.x) << '\n';
		std::cout << "        y: " << sqrt(var3D.y) << '\n';
		std::cout << "        z: " << sqrt(var3D.z) << std::endl;
	}
	
	return 0;
}
