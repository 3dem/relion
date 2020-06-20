#include <src/args.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/atomic/pdb_helper.h>


using namespace gravis;


void addPoint(d3Vector pos, d3Vector ref, double scale, double rad, double out_bin, Mesh& mesh)
{
	const d3Vector pos2 = ref + scale * (pos - ref);	
	MeshBuilder::addOctahedron(pos2 / out_bin, rad, mesh);
}

int main(int argc, char *argv[])
{
	std::string inFn, outDir;
	double psModel, psOut;
	int boxModel, boxOut;
			
	IOParser parser;
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		inFn = parser.getOption("--i", "Input PDB file");
		outDir = parser.getOption("--o", "Output directory");
		boxModel = textToInteger(parser.getOption("--box_model", "Box size of the map corresponding to the PDB file"));
		psModel = textToDouble(parser.getOption("--angpix_model", "Pixel size of the map corresponding to the PDB file"));
		boxOut = textToInteger(parser.getOption("--box_out", "Box size of the map to be compared"));
		psOut = textToDouble(parser.getOption("--angpix_out", "Pixel size of the map to be compared"));
		
				
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	outDir = ZIO::makeOutputDir(outDir);
	
	const double coordOffset = (boxOut/2)*psOut - (boxModel/2)*psModel;

	Assembly assembly;
	assembly.readPDB(inFn);

	std::map<std::string,std::vector<d3Vector>> atoms = PdbHelper::groupAtomsByElement(assembly);

	std::cout << atoms.size() << " elements found" << std::endl;
	std::cout << "    " << assembly.molecules.size() << " molecules" << std::endl;
	std::cout << "    " << assembly.looseAtoms.size() << " lose atoms" << std::endl;
	
	for (std::map<std::string,std::vector<d3Vector>>::iterator it = atoms.begin();
		 it != atoms.end(); it++)
	{
		Mesh mesh;
		
		const std::string element = it->first;
		const std::vector<d3Vector>& positions = it->second;

		std::cout << element << " (" << positions.size() << ")" << std::endl;
		
		const int pc = positions.size();
		
		for (int p = 0; p < pc; p++)
		{
			const d3Vector pos = positions[p] + coordOffset * d3Vector(1,1,1);
			MeshBuilder::addOctahedron(pos, 0.5, mesh);
		}
		
		mesh.writePly(outDir+element+".ply");
	}

	return 0;
}
