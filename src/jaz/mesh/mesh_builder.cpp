#include "mesh_builder.h"
#include <src/macros.h>

using namespace gravis;

Mesh MeshBuilder::makeBar(
		d3Vector p0, 
		d3Vector p1, 
		double rad, 
		int segments)
{
	d3Vector d = p1 - p0;	
	d3Vector q;
	
	if (std::abs(d.x) < std::abs(d.y)) 
	{
		q = d.cross(d3Vector(1,0,0)).normalize();
	}
	else 
	{
		q = d.cross(d3Vector(0,1,0)).normalize();
	}
	
	d3Vector r = q.cross(d).normalize();
	
	
	Mesh bar;
		
	bar.vertices.resize(2*segments);
	
	for (int s = 0; s < segments; s++)
	{
		const double phi = s * 2.0 * PI / (double) (segments);
		const double fq = rad * cos(phi);
		const double fr = rad * sin(phi);
		
		bar.vertices[2*s] = p0 + fq * q + fr * r;
		bar.vertices[2*s + 1] = bar.vertices[2*s] + d;
	}
	
	bar.triangles.resize(2*segments);
	
	for (int s = 0; s < segments; s++)
	{
		const int a = 2*s;
		const int b = 2*s + 1;
		const int c = 2*((s+1)%segments) + 1;
		const int d = 2*((s+1)%segments);
		
		bar.triangles[2*s] = Triangle(a,b,c);
		bar.triangles[2*s+1] = Triangle(a,c,d);
	}
	
	return bar;
}

void MeshBuilder::addBar(
		d3Vector p0, 
		d3Vector p1, 
		double rad, 
		int segments, 
		Mesh& mesh)
{
	Mesh bar = makeBar(p0,p1,rad,segments);
	insert(bar, mesh);
}

void MeshBuilder::addColouredBar(
		d3Vector p0, 
		d3Vector p1, 
		double rad, 
		int segments, 
		gravis::dRGB colour,
		Mesh& mesh)
{
	Mesh bar = makeBar(p0,p1,rad,segments);
	insert(bar, mesh, colour);
}

void MeshBuilder::addPointyBar(
		d3Vector p0, 
		d3Vector p1, 
		double aspect, 
		double pointedness, 
		int segments, 
		Mesh &mesh)
{
	d3Vector d = p1 - p0;
	
	const double len = d.length();
	const double rad = len / aspect;
	
	d3Vector q;
	
	if (std::abs(d.x) < std::abs(d.y)) 
	{
		q = d.cross(d3Vector(1,0,0)).normalize();
	}
	else 
	{
		q = d.cross(d3Vector(0,1,0)).normalize();
	}
	
	d3Vector r = q.cross(d).normalize();
	
	
	Mesh bar;
		
	bar.vertices.resize(2*segments + 2);
	
	const int ip0 = 2*segments;
	const int ip1 = 2*segments + 1;
	
	const d3Vector pm0 = p0 + 0.5 * pointedness * d;
	const d3Vector pm1 = p1 - 0.5 * pointedness * d;
	
	for (int s = 0; s < segments; s++)
	{
		const double phi = s * 2.0 * PI / (double) (segments);
		const double fq = rad * cos(phi);
		const double fr = rad * sin(phi);
		
		bar.vertices[2*s]     = pm0 + fq * q + fr * r;
		bar.vertices[2*s + 1] = pm1 + fq * q + fr * r;
	}
	
	bar.vertices[ip0] = p0;
	bar.vertices[ip1] = p1;
	
	bar.triangles.resize(4*segments);
	
	/*
	    p1
	  b    c
	  a    d
	    p0
	         */
	
	for (int s = 0; s < segments; s++)
	{
		const int a = 2*s;
		const int b = 2*s + 1;
		const int c = 2*((s+1)%segments) + 1;
		const int d = 2*((s+1)%segments);
		
		bar.triangles[4*s]     = Triangle(a,b,c);
		bar.triangles[4*s + 1] = Triangle(a,c,d);
		bar.triangles[4*s + 2] = Triangle(a,d,ip0);
		bar.triangles[4*s + 3] = Triangle(c,b,ip1);
	}
	
	insert(bar, mesh);	
}

void MeshBuilder::addCone(d3Vector p0, d3Vector p1, double rad, int segments, Mesh &mesh)
{
	d3Vector d = p1 - p0;	
	d3Vector q;
	
	if (std::abs(d.x) < std::abs(d.y)) 
	{
		q = d.cross(d3Vector(1,0,0)).normalize();
	}
	else 
	{
		q = d.cross(d3Vector(0,1,0)).normalize();
	}
	
	d3Vector r = q.cross(d).normalize();
	
	
	Mesh cone;
		
	cone.vertices.resize(segments + 1);
	
	for (int s = 0; s < segments; s++)
	{
		const double phi = s * 2.0 * PI / (double) (segments);
		const double fq = rad * cos(phi);
		const double fr = rad * sin(phi);
		
		cone.vertices[s] = p0 + fq * q + fr * r;
	}
	
	cone.vertices[segments] = p1;	
	cone.triangles.resize(segments);
		
	for (int s = 0; s < segments; s++)
	{
		const int a = s;
		const int b = (s + 1) % segments;
		const int c = segments;
		
		cone.triangles[s] = Triangle(a,b,c);
	}
	
	insert(cone, mesh);	
}

void MeshBuilder::addOctahedron(d3Vector centre, double rad, Mesh &mesh)
{	
	const int v0 = mesh.vertices.size();
	
	mesh.vertices.reserve(mesh.vertices.size() + 6);
	
	mesh.vertices.push_back(centre + d3Vector(rad,0,0));
	mesh.vertices.push_back(centre - d3Vector(rad,0,0));
	mesh.vertices.push_back(centre + d3Vector(0,rad,0));
	mesh.vertices.push_back(centre - d3Vector(0,rad,0));
	mesh.vertices.push_back(centre + d3Vector(0,0,rad));
	mesh.vertices.push_back(centre - d3Vector(0,0,rad));
	
	mesh.triangles.reserve(mesh.triangles.size() + 8);
	
	mesh.triangles.push_back(Triangle(v0,   v0+2, v0+4));
	mesh.triangles.push_back(Triangle(v0+2, v0+1, v0+4));
	mesh.triangles.push_back(Triangle(v0+1, v0+3, v0+4));
	mesh.triangles.push_back(Triangle(v0+3, v0,   v0+4));
	
	mesh.triangles.push_back(Triangle(v0+2, v0,   v0+5));
	mesh.triangles.push_back(Triangle(v0+1, v0+2, v0+5));
	mesh.triangles.push_back(Triangle(v0+3, v0+1, v0+5));
	mesh.triangles.push_back(Triangle(v0,   v0+3, v0+5));
}

void MeshBuilder::insert(const Mesh& src, Mesh& dest)
{
	const int vcd = dest.vertices.size();	
	const int vcs = src.vertices.size();
	const int tcs = src.triangles.size();	
	
	for (int v = 0; v < vcs; v++)
	{
		dest.vertices.push_back(src.vertices[v]);
	}
	
	for (int t = 0; t < tcs; t++)
	{
		const Triangle t0 = src.triangles[t];
		
		dest.triangles.push_back(Triangle(t0.a + vcd, t0.b + vcd, t0.c + vcd));
	}
}

void MeshBuilder::insert(const Mesh &src, Mesh &dest, dRGB color)
{
	const int vcd = dest.vertices.size();	
	const int vcs = src.vertices.size();
	const int tcs = src.triangles.size();	
	
	for (int v = 0; v < vcs; v++)
	{
		dest.vertices.push_back(src.vertices[v]);
		dest.vertexColors.push_back(color);
	}
	
	for (int t = 0; t < tcs; t++)
	{
		const Triangle t0 = src.triangles[t];
		
		dest.triangles.push_back(Triangle(t0.a + vcd, t0.b + vcd, t0.c + vcd));
	}
}
