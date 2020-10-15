#ifndef MESH_BUILDER_H
#define MESH_BUILDER_H

#include "mesh.h"

class MeshBuilder
{
	public:
		
		static Mesh makeBar(
				gravis::d3Vector p0, 
				gravis::d3Vector p1, 
				double rad, 
				int segments);
		
				
		static void addBar(
						gravis::d3Vector p0, 
						gravis::d3Vector p1, 
						double rad, 
						int segments,
						Mesh& mesh);
		
		static void addColouredBar(
						gravis::d3Vector p0, 
						gravis::d3Vector p1, 
						double rad, 
						int segments,
						gravis::dRGB colour,
						Mesh& mesh);
		
		static void addPointyBar(
						gravis::d3Vector p0, 
						gravis::d3Vector p1, 
						double aspect, 
						double pointedness, 
						int segments,
						Mesh& mesh);
		
		static void addCone(
						gravis::d3Vector p0, 
						gravis::d3Vector p1, 
						double rad, 
						int segments,
						Mesh& mesh);
		
		static void addOctahedron(
						gravis::d3Vector centre, 
						double rad, 
						Mesh& mesh);
		
		static void insert(const Mesh& src, Mesh& dest);
		
		static void insert(const Mesh& src, Mesh& dest, gravis::dRGB color);
};

#endif
