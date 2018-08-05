#include <src/jaz/vtk_helper.h>
#include <src/jaz/math/Zernike.h>

int main(int argc, char *argv[])
{
	const int N = 6;
	const int s = 400;
	
	for (int n = 1; n <= N; n++)
	{
		Image<RFLOAT> out(s,s,n+1);
		
		for (int i = 0; i <= n; i++)
		{
			const int m = 2*i - n;
			
			std::cout << n << ", " << m << " -> " << i << "\n";
			
			for (int y = 0; y < s; y++)
			for (int x = 0; x < s; x++)
			{
				double xx = 2.0*x/(double)s - 1.0;
				double yy = 2.0*y/(double)s - 1.0;
				
				double rho = sqrt(xx*xx + yy*yy);
				double phi = atan2(yy,xx);
				
				out(i,y,x) = Zernike::Z(m,n,rho,phi);
			}
		}
		
		std::stringstream sts;
		sts << "Zernike_" << n << ".vtk";
		
		VtkHelper::writeVTK(out,sts.str());
	}
	
	return 0;
}
