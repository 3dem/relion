
#include <src/jaz/tomo/tomo_stack.h>
#include <src/jaz/tomo/backprojection_helper.h>
#include <src/jaz/vtk_helper.h>

using namespace gravis;

int main(int argc, char *argv[])
{			
	double angpix = 1.177;
	double scale = 1.0;
	const double bin0 = 8.0;
	const double bin = 8.0;
	
	/*TomoStack ts(
		"frames/TS_03_f*.mrc", 40, "TS_03.tlt", 
		 "TS_03.xf", "TS_03_output.txt", 
		 angpix*scale, 1.0/scale);
	
	std::cout << "loading done.\n";
	
	ts.downsample(bin);
	ts.saveImages("frames/bin8_*.mrc");
	
	std::cout << "downsampling done.\n";
	*/
	
	TomoStack ts(
			"frames/bin8_*.mrc", 40, "TS_03.tlt", 
			 "TS_03.xf", "TS_03_output.txt", 
			 angpix*bin, 1.0/bin);
		
	std::cout << "loading done.\n";
	
	
	d3Vector cc(0, 0, -250);
    //double areaSize = 200;
	
	Volume<RFLOAT> dest(200,200,80);
	Volume<RFLOAT> maskDest(200,200,80);
	
	dest.fill(0.0);
	maskDest.fill(0.0);
	
	std::cout << "filling done.\n";
	
	//BackprojectionHelper::backprojectExactWeights(ts, dest, cc, scale);
	/*, double taperX, double taperY, double taperZ, double wMin,
		int frame0, int frames*/
				
				//ts, wbp, origin, 2, 100, 100, 100);
	
	const double spacing = 3710.0/200.0;
	BackprojectionHelper::backprojectRaw(ts, dest, maskDest, cc, spacing);
	
	VtkHelper::writeVTK(
		dest, "test00.vtk", 
		cc.x/2, cc.y/2, cc.z/2, 
		spacing/2, spacing/2, spacing/2);
	
	
	return 0;
}
