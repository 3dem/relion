#include "skeletonization.h"

using namespace gravis;

BufferedImage<Tensor3x3<float>> Skeletonization :: getKernel(
	const RawImage<Tensor3x3<float>>& J,
	Mode mode, float tol)
{
	BufferedImage<Tensor3x3<float>> out(J.xdim, J.ydim, J.zdim);
	
	const int w = J.xdim;
	const int h = J.ydim;
	const int d = J.zdim;
	
	const float vn = -1.f - tol;
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < w; x++)
	{
		d3Matrix Q;
		d3Vector Lvec;
		
		Tensor3x3<double> Jd(J(x,y,z).toMatrix());
		
		Jd.diagonalize(Lvec, Q);
		
		d3Matrix L;
		
		switch (mode)
		{
			case Point:
				
				L = d3Matrix(
						 1,  0,  0, 
						 0,  1,  0,
						 0,  0,  1);
			break;
				
			case Curve:
				
				L = d3Matrix(
						 1,  0,  0,
						 0,  1,  0,
						 0,  0, vn);
			break;
				
			case Surface:
				
				L = d3Matrix(
						 1,  0,  0,
						 0, vn,  0,
						 0,  0, vn);
			break;
		}
		
		d3Matrix Qt = Q;
		Qt.transpose();
		
		d3Matrix Dm = Q * L * Qt;
		
		out(x,y,z) = Tensor3x3<float>(Dm);
	}
	
	return out;
}

std::vector<d3Vector> Skeletonization::discretize(
		const RawImage<float>& skeleton,
		float minVal)
{
	std::vector<d3Vector> out(0);
	out.reserve(skeleton.xdim);
	
	const int w = skeleton.xdim;
	const int h = skeleton.ydim;
	const int d = skeleton.zdim;
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < w; x++)
	{
		if (skeleton(x,y,z) > minVal)
		{
			out.push_back(d3Vector(x,y,z));
		}
	}
	
	return out;
}
