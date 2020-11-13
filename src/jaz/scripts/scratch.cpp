#include <src/jaz/math/Euler_angles_dynamo.h>
#include <src/macros.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/buffered_image.h>
#include <omp.h>


using namespace gravis;

int main(int argc, char *argv[])
{
	const int s = 32;
	BufferedImage<double> test(s,s);
	
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		test(x,y) = rand() / (double)RAND_MAX;
	}
	
	const int s2 = 512;
	BufferedImage<double> test2(s2,s2,5);
	
	const double scale = s / (double)s2;
	const double eps = 1e-6;
	
	for (int y = 0; y < s2; y++)
	for (int x = 0; x < s2; x++)
	{
		test2(x,y,0) = Interpolation::cubicXY_clip(test, scale*x, scale*y, 0);
		
		d2Vector g = Interpolation::cubicXYGrad_clip(test, scale*x, scale*y, 0);
		
		const double vxp = Interpolation::cubicXY_clip(test, scale*x + eps, scale*y, 0);
		const double vxn = Interpolation::cubicXY_clip(test, scale*x - eps, scale*y, 0);
		const double vyp = Interpolation::cubicXY_clip(test, scale*x, scale*y + eps, 0);
		const double vyn = Interpolation::cubicXY_clip(test, scale*x, scale*y - eps, 0);
		
		test2(x,y,1) = g.x;
		test2(x,y,2) = (vxp - vxn) / (2.0 * eps);
		
		test2(x,y,3) = g.y;
		test2(x,y,4) = (vyp - vyn) / (2.0 * eps);
		
	}
	
	test2.write("test2.mrc");
	
	return 0;
}
