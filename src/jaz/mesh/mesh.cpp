#include "mesh.h"
#include <fstream>

using namespace gravis;

Mesh::Mesh()
	: vertices(0), vertexColors(0), triangles(0)
{}

void Mesh::writeObj(std::string filename)
{
	std::ofstream ofs(filename);
	
	for (int i = 0; i < vertices.size(); i++)
	{
		const d3Vector v = vertices[i];
		ofs << "v " << v.x << " " << v.y << " " << v.z << std::endl;
	}
	ofs << "# " << vertices.size() << " vertices." << std::endl << std::endl;
	
	for (int i = 0; i < triangles.size(); i++)
	{
		const Triangle t = triangles[i];
		ofs << "f " << (t.a + 1) << " " << (t.b + 1) << " " << (t.c + 1) << std::endl;
	}
	ofs << "# " << triangles.size() << " faces." << std::endl << std::endl;
}

void Mesh::writePly(std::string filename)
{
	std::ofstream ofs(filename);
	
	const int vc = vertices.size();
	const bool hasColors = vertexColors.size() == vc;
	
	ofs << "ply\n";
	ofs << "format ascii 1.0\n";
	ofs << "element vertex " << vertices.size() << "\n";
	
	ofs << "property float32 x\n";
	ofs << "property float32 y\n";
	ofs << "property float32 z\n";
	
	if (hasColors)
	{
		ofs << "property uint8 red\n";
		ofs << "property uint8 green\n";
		ofs << "property uint8 blue\n";
	}
	
	ofs << "element face " << triangles.size() << "\n";
	ofs << "property list uint8 int32 vertex_index\n";
	ofs << "end_header\n";
		
	for (int i = 0; i < vertices.size(); i++)
	{
		const d3Vector v = vertices[i];
		ofs << v.x << " " << v.y << " " << v.z;
		
		if (hasColors)
		{
			const dRGB c = 255.0 * vertexColors[i];
			ofs << " " << (int)(c.r) << " " << (int)(c.g) << " " << (int)(c.b);			
		}
		
		ofs << std::endl;
	}
	
	for (int i = 0; i < triangles.size(); i++)
	{
		const Triangle t = triangles[i];
		ofs << "3 " << t.a << " " << t.b << " " << t.c << std::endl;
	}
}
