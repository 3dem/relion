#include <src/jaz/tomography/motion/modular_alignment/spline_2D_deformation_model.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/util/zio.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	const int w = 5;
	const int h = 4;
	
	BufferedImage<d4Vector> data(w,h,2);
	
	for (int z = 0; z < 2; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	for (int i = 0; i < 4; i++)
	{
		data(x,y,z)[i] = 1.0 - 2.0 * (rand() / (double) RAND_MAX);
	}
	
	
	const int ww = 1024;
	const int hh = 768;
	
	Deformation2D::Parameters parameters;
	parameters.grid_width  = w;
	parameters.grid_height = h;

	Spline2DDeformationModel model(parameters, i2Vector(ww, hh));
	
	/*BufferedImage<double> test0(ww,hh), test1(ww,hh);
	
	for (int y = 0; y < hh; y++)
	for (int x = 0; x < ww; x++)
	{
		d2Vector def, def_x, def_y;
		
		model.computeShiftAndGradient(
			d2Vector(x,y), (double*)(data.data), def, def_x, def_y);
		
		test0(x,y) = def_y.x;
	}
	
	test0.write("spline_test_gradXy-a.mrc");
	
	for (int y = 1; y < hh-1; y++)
	for (int x = 0; x < ww; x++)
	{
		d2Vector def1, def1_x, def1_y;
		d2Vector def0, def0_x, def0_y;
		
		model.computeShiftAndGradient(
			d2Vector(x,y-1), (double*)(data.data), def0, def0_x, def0_y);
		
		model.computeShiftAndGradient(
			d2Vector(x,y+1), (double*)(data.data), def1, def1_x, def1_y);
		
		test1(x,y) = (def1.x - def0.x) / 2.0;
	}
	
	test1.write("spline_test_gradXy-n.mrc");*/
	
	
	BufferedImage<d4Vector> data0 = data;
	
	const double eps = 1e-3;
	const int dim = 0;
	
	for (int dy = 0; dy < 2; dy++)
	for (int dx = 0; dx < 2; dx++)
	for (int  i = 0;  i < 4;  i++)
	{
		data = data0;
		
		d2Vector def1, def1_x, def1_y;
		d2Vector def0, def0_x, def0_y;
		
		const d2Vector pl(1.7 * model.gridSpacing.x, 1.2 * model.gridSpacing.y);
		
		model.computeShiftAndGradient(
			pl, (double*)(data.data), def0, def0_x, def0_y);
		
		data(1+dx, 1+dy, dim)[i] += eps;
		
		model.computeShiftAndGradient(
			pl, (double*)(data.data), def1, def1_x, def1_y);
		
		const d2Vector d = def1 - def0;
		const d2Vector g0(1,2);
		
		const double m = g0.dot(d) / eps;
		
		
		BufferedImage<d4Vector> grad(w,h,2);
		grad.fill(d4Vector(0.0,0.0,0.0,0.0));
		
		model.updateDataTermGradient(pl, g0, 0, (double*)grad.data);
		
		std::cout << std::setprecision(12) << std::scientific
				  << grad(1+dx, 1+dy, dim)[i] << " \t vs \t " << m << " \t ("
				  << (grad(1+dx, 1+dy, dim)[i] - m) << ')' << std::endl;
	}

	return 0;
}
