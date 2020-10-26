#include <src/args.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/particle_set.h>

using namespace gravis;


void addPoint(const d3Vector& pos, const d3Matrix& A, double rad, double out_bin, Mesh& mesh)
{
	const d3Vector n(A(0,2), A(1,2), A(2,2));

	MeshBuilder::addCone(pos / out_bin, (pos + 2 * rad * n) / out_bin, rad, 12, mesh);
}

int main(int argc, char *argv[])
{
	std::string inFn, refFn, tomoSetFn, motFn, outDir;
	int t0, t1, f;
	double bin, rad, relScale;
	bool allFrames;
	

	IOParser parser;
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		inFn = parser.getOption("--i", "Input particle set");
		refFn = parser.getOption("--ref", "Other particles file to serve as reference to zoom against", "");
		tomoSetFn = parser.getOption("--t", "Tomogram set", "tomograms.star");
		
		t0 = textToInteger(parser.getOption("--t0", "First tomogram", "0"));
		t1 = textToInteger(parser.getOption("--t1", "Last tomogram", "-1"));
		
		motFn = parser.getOption("--mot", "Particle trajectories", "");
		
		f = textToInteger(parser.getOption("--f", "Frame", "0"));
		allFrames = parser.checkOption("--all_frames", "Output meshes for all frames");
		
		bin = textToDouble(parser.getOption("--bin", "Divide coordinates by this", "1.0"));
		rad = textToDouble(parser.getOption("--rad", "Size of output octahedron", "1.0"));
		relScale = textToDouble(parser.getOption("--rel_scale", "Scale the distance from reference pt. by this", "8.0"));
		
		outDir = parser.getOption("--o", "Output directory");
				
				
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	if ((f != 0 || allFrames) && (motFn == ""))
	{
		std::cerr << "A motion file (--mot) is required to select a frame (--f) or output all frames (--all_frames)." << std::endl;
		return -1;
	}
	
	if (f < 0)
	{
		std::cerr << "The frame index (--f) must not be negative." << std::endl;
		return -1;
	}
			
	
	TomogramSet tomogramSet(tomoSetFn);
	
	ParticleSet dataSet(inFn, motFn);
	std::vector<std::vector<ParticleIndex>> particles = dataSet.splitByTomogram(tomogramSet);
	

	const int tc = particles.size();
	
	if (t1 < 0) t1 = tc - 1;
		
	outDir = ZIO::makeOutputDir(outDir);
	
	
	for (int t = t0; t <= t1; t++)
	{
		const int pc = particles[t].size();
		
		if (pc == 0) continue;

		Tomogram tomogram = tomogramSet.loadTomogram(t, false);

		std::vector<int> frameSeq = IndexSort<double>::sortIndices(tomogram.cumulativeDose);
		
		const int fc = tomogram.frameCount;
		dataSet.checkTrajectoryLengths(particles[t][0], pc, fc, "backproject");
		
		
		Mesh mesh;		
		std::vector<Mesh> meshes(fc);
		
		for (int p = 0; p < pc; p++)
		{
			const ParticleIndex part_id = particles[t][p];

			const d3Matrix A = dataSet.getMatrix3x3(part_id);
			
			
			if (f == 0 && !allFrames)
			{
				const d3Vector pos = dataSet.getPosition(part_id);
				
				addPoint(pos, A, rad, bin, mesh);
			}
			else
			{
				const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(
							part_id, 0, tomogram.optics.pixelSize);
				
				if (allFrames)
				{
					const int fc = traj.size();
					
					for (int f = 0; f < fc; f++)
					{
						const d3Vector pos = traj[frameSeq[f]];
						
						addPoint(pos, A, rad, bin, meshes[f]);
					}
				}
				else
				{
					if (f >= traj.size())
					{
						if (motFn != "")
						{
							REPORT_ERROR_STR(
								"Bad frame index: only " << traj.size() << " frames available for particle "
								<< part_id.value << " in " << motFn);
						}
						else
						{
							// has been handled
						}
					}
					
					const d3Vector pos = traj[frameSeq[f]];
					
					addPoint(pos, A, rad, bin, mesh);
				}
			}
		}
		
		if (allFrames)
		{
			for (int f = 0; f < fc; f++)
			{
				meshes[f].writePly(outDir + ZIO::itoa(t) + "_frame_" + ZIO::itoa(f) + ".ply");
			}	
		}
		else
		{
			mesh.writePly(outDir + ZIO::itoa(t) + ".ply");
		}
	}

	return 0;
}
