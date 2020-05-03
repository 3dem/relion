#include "primal_dual_TV.h"
#include <src/jaz/image/gradient.h>
#include <src/jaz/image/divergence.h>
#include <src/jaz/gravis/t3Vector.h>

using namespace gravis;

BufferedImage<float> PrimalDualTV::anisotropicTV(
	RawImage<float> regionalCost, 
	RawImage<Tensor3x3<float>> diffusionTensor, 
	int maxIterations, float sigma, float tau, float nu)
{
	const int w = regionalCost.xdim;
	const int h = regionalCost.ydim;
	const int d = regionalCost.zdim;
	
	BufferedImage<float> u(w,h,d), uBar(w,h,d), divDXi(w,h,d);
	BufferedImage<f3Vector> xi(w,h,d), Dxi(w,h,d), gradUBar(w,h,d);
	
	u.fill(0.5f);
	uBar.fill(0.5f);
		
	const float eps = 0.0001f;
	
	float lastMaxDelta = -1.f;
	
	for (int it = 0; it < maxIterations; it++)
	{
		if (it%10 == 0) std::cout << "iteration " << it << " - ";
		
		Gradient::forward3D_inSitu(uBar, gradUBar);
				
		for (int z = 0; z < d; z++)
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			f3Vector newXi = xi(x,y,z);
			f3Vector grd = gradUBar(x,y,z);
			f3Matrix D = diffusionTensor(x,y,z).toMatrix();
			
			newXi += sigma * D * grd;
			
			const float xl = newXi.length();
			
			if (xl > eps)
			{
				newXi = newXi / xl;
			}
			
			xi(x,y,z) = newXi;
			Dxi(x,y,z) = D * newXi;
		}
		
		Divergence::backward3D_inSitu(Dxi, divDXi);
		
		float maxDelta = 0.f;
		
		for (int z = 0; z < d; z++)
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			const float fv = -regionalCost(x,y,z);
			
			float uv = u(x,y,z) + tau * (divDXi(x,y,z) - fv);
			
			if (uv > 1.f) uv = 1.f;
			else if (uv < 0.f) uv = 0.f;
			
			const float lastUv = u(x,y,z);
			
			u(x,y,z) = uv;
			uBar(x,y,z) = 2.f * uv - lastUv;
			
			const float delta = uv - lastUv;
			if (fabs(delta) > fabs(maxDelta)) maxDelta = delta;
		}
		
		if (it%10 == 0) std::cout << "maxDelta = " << maxDelta << "\n";
		
		lastMaxDelta = maxDelta;
	}
	
	return u;
}
