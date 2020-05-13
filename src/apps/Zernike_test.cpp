#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/math/Zernike.h>

int main(int argc, char *argv[])
{
	const int N = 12;
	const int s = 400;
	
	Image<RFLOAT> out(s,s,N);
	
	for (int i = 0; i < N; i++)
	{
		int m, n;
		Zernike::oddIndexToMN(i,m,n);
			
		std::cout << i << " -> " << m << ", " << n << "\n";
		
		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			double xx = 2.0*x/(double)s - 1.0;
			double yy = 2.0*y/(double)s - 1.0;
						
			out(i,y,x) = Zernike::Z_cart(m,n,xx,yy);
		}
	}

	VtkHelper::writeVTK(out, "Zernike-odd-test.vtk");
	
	return RELION_EXIT_SUCCESS;
}
