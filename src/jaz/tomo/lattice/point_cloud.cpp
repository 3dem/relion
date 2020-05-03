#include "point_cloud.h"
#include <src/jaz/mesh/mesh_builder.h>
#include <src/error.h>
#include <fstream>
#include <string.h>

using namespace gravis;

PointCloud::PointCloud()
{}

void PointCloud::readCsv(std::string fn)
{
	std::ifstream ifs(fn);
	
	if (!ifs.is_open())
	{
		REPORT_ERROR("Unable to read: " + fn);
	}
	
	char text[4096];
	
	vertices.reserve(1000);

	while (ifs.good())
	{
		ifs.getline(text, 4096);
		
		const size_t tlen = strlen(text);
		
		if (tlen < 5) continue;
		
		for (int i = 0; i < tlen; i++)
		{
			if (text[i] == ',') text[i] = ' ';
		}

		std::stringstream line(text);
		
		d3Vector v;
		
		line >> v.x;
		line >> v.y;
		line >> v.z;
		
		vertices.push_back(v);
	}
}

void PointCloud::findNeighbors(double minDist, double maxDist, int expected)
{
	const int pc = vertices.size();
	neighbors.resize(pc);
	
	for (int p = 0; p < pc; p++)
	{
		neighbors[p].reserve(expected);
	}
	
	for (int p = 0; p < pc; p++)
	for (int q = p+1; q < pc; q++)
	{
		const double d = (vertices[p] - vertices[q]).length();
		
		if (d < maxDist && d > minDist)
		{
			neighbors[p].push_back(q);
			neighbors[q].push_back(p);
		}
	}
}

std::vector<d3Vector> PointCloud::getNeighborhood(int vertex) const
{
	std::vector<d3Vector> out(neighbors[vertex].size());
	
	for (int n = 0; n < neighbors[vertex].size(); n++)
	{
		out[n] = vertices[neighbors[vertex][n]] - vertices[vertex];
	}
	
	return out;
}

Mesh PointCloud::visualizeNeighbors()
{
	Mesh out;
	
	const int pc = vertices.size();
	
	out.vertices.reserve(5*pc);
	out.triangles.reserve(5*pc);
	
	for (int p = 0; p < pc; p++)
	{
		const d3Vector vp = vertices[p];
		
		for (int q = 0; q < neighbors[p].size(); q++)
		{
			const d3Vector vq = vertices[neighbors[p][q]];
			MeshBuilder::addBar(vp, vp + 0.48*(vq - vp), 0.1, 5, out);
		}
	}
	
	return out;
}

std::vector<size_t> PointCloud::distanceHist(int bins, double maxDist)
{
	const int pc = vertices.size();
	
	std::vector<size_t> out(bins, 0);
	
	for (int p = 0; p < pc; p++)
	for (int q = p+1; q < pc; q++)
	{
		const double d = (vertices[p] - vertices[q]).length();
		
		if (d < maxDist)
		{
			out[(int)(bins * d / maxDist)]++; 
		}
	}
	
	return out;
}
