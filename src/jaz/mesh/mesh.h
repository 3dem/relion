#ifndef MESH_H
#define MESH_H

#include "triangle.h"
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/gravis/tRGB.h>
#include <string>
#include <vector>

class Mesh
{
	public:
			
		Mesh();
		
			std::vector<gravis::d3Vector> vertices;
			std::vector<gravis::dRGB> vertexColors;
			std::vector<Triangle> triangles;
			
		void writeObj(std::string filename);
		void writePly(std::string filename);
		
};

#endif
