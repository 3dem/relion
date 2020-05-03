
#include <src/jaz/single_particle/tomo/tomo_stack.h>
#include <src/jaz/single_particle/tomo/backprojection_helper.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/single_particle/volume_converter.h>

using namespace gravis;

void drawPoint(Volume<RFLOAT>& dest, d3Vector p, RFLOAT val);
			   
int main(int argc, char *argv[])
{			
	const double angpix = 1.177;
	const double bin = 4.0;
	const int frameCount = 40;
	const double fullWidth = 3710.0;
	
	std::string tomoFn = "frames/TS_03_f*.mrc";
	std::string tltFn = "TS_03.tlt";
	std::string xfFn = "TS_03.xf";
	std::string aliFn = "TS_03_output.txt";
	std::string fiducialsFile = "TS_03_3dmod.txt";
	std::string particlesStar = "allmotl_TS_03.star";
	
	const bool prescaled = true;
	
	TomoStack ts;
	
	if (!prescaled)
	{
		ts = TomoStack(
			tomoFn, frameCount, tltFn, 
			xfFn, aliFn, 
			angpix, 1.0);
		
		std::cout << "loading done.\n";
		
		std::stringstream sts;
		sts << bin;
		
		ts.downsample(bin);
		ts.saveImages("frames/bin"+sts.str()+"_*.mrc");
		
		std::cout << "downsampling done.\n";
	}
	else
	{
		std::stringstream sts;
		sts << bin;
		
		ts = TomoStack(
			"frames/bin"+sts.str()+"_*.mrc", 
			frameCount, tltFn, xfFn, aliFn, 
			angpix*bin, 1.0/bin);
		
		std::cout << "loading done.\n";
	}
	
	const int w = 400;
	const int h = 400;
	const int d = 200;
	
	d3Vector origin(0.0, 0.0, -400.0);
	const double spacing = fullWidth/w;
		
	Volume<RFLOAT> dest(w,h,d);
	Volume<RFLOAT> maskDest(w,h,d);
	
	dest.fill(0.0);
	maskDest.fill(0.0);
	
	std::cout << "filling done.\n";
		
	BackprojectionHelper::backprojectRaw(ts, dest, maskDest, origin, spacing);
	
	// write clean tomogram into test00.vtk
	
	VtkHelper::writeVTK(
		dest, "test00.vtk", 
		origin.x, 
		origin.y, 
		origin.z, 
		spacing, 
		spacing, 
		spacing);
	
	Image<RFLOAT> destImg;
	VolumeConverter::convert(dest, destImg);
	destImg.write("test00.mrc");
	
	// add fiducials
	
	std::ifstream fidFile(fiducialsFile);
	std::vector<d3Vector> fids;
	
	char text[4096];

	while (fidFile.getline(text, 4096))
	{
		std::stringstream line(text);
		
		d3Vector fid;
		line >> fid.x;
		line >> fid.y;
		line >> fid.z;
		
		fids.push_back(fid);
	}
	
	for (int f = 0; f < fids.size(); f++)
	{
		d3Vector fidG = (2.0 * fids[f] - origin) / spacing;
		drawPoint(dest, fidG, 1000);
	}
	
	// write tomogram with fiducials into test01.vtk
	
	VtkHelper::writeVTK(
		dest, "test01.vtk", 
		origin.x, 
		origin.y, 
		origin.z, 
		spacing, 
		spacing, 
		spacing);
	
	// add particles
	
	const double z_offset = 450.0;
	
	MetaDataTable partMdt;
	partMdt.read(particlesStar);
	
	for (int p = 0; p < partMdt.numberOfObjects(); p++)
	{
		d3Vector partCoord, partOff;
		
		partMdt.getValue(EMDL_IMAGE_COORD_X, partCoord.x, p);
		partMdt.getValue(EMDL_IMAGE_COORD_Y, partCoord.y, p);
		partMdt.getValue(EMDL_IMAGE_COORD_Z, partCoord.z, p);
		
		partMdt.getValue(EMDL_ORIENT_ORIGIN_X, partOff.x, p);
		partMdt.getValue(EMDL_ORIENT_ORIGIN_Y, partOff.y, p);
		partMdt.getValue(EMDL_ORIENT_ORIGIN_Z, partOff.z, p);
		
		d3Vector partPos = partCoord + partOff;
		partPos.z -= z_offset;
		
		d3Vector posVol = (partPos - origin) / spacing;
		
		drawPoint(dest, posVol, 1000);
	}
	
	VtkHelper::writeVTK(
		dest, "test02.vtk", 
		origin.x, 
		origin.y, 
		origin.z, 
		spacing, 
		spacing, 
		spacing);
	
	return RELION_EXIT_SUCCESS;
}

void drawPoint(Volume<RFLOAT>& dest, d3Vector p, RFLOAT val)
{
	const int w = dest.dimx;
	const int h = dest.dimy;
	const int d = dest.dimz;
	
	if (p.x > 0.0 && p.y > 0.0 && p.z > 0.0 && 
		p.x < w-1 && p.y < h-1 && p.z < d-1)
	{
		int x0 = (int) p.x;
		int y0 = (int) p.y;
		int z0 = (int) p.z;
		
		int x1 = (int) p.x + 1;
		int y1 = (int) p.y + 1;
		int z1 = (int) p.z + 1;
		
		const double ix1 = p.x - x0;
		const double iy1 = p.y - y0;
		const double iz1 = p.z - z0;
		
		const double ix0 = 1.0 - ix1;
		const double iy0 = 1.0 - iy1;
		const double iz0 = 1.0 - iz1;
		
		dest(x0,y0,z0) += val * ix0 * iy0 * iz0;
		dest(x1,y0,z0) += val * ix1 * iy0 * iz0;
		dest(x0,y1,z0) += val * ix0 * iy1 * iz0;
		dest(x1,y1,z0) += val * ix1 * iy1 * iz0;
		dest(x0,y0,z1) += val * ix0 * iy0 * iz1;
		dest(x1,y0,z1) += val * ix1 * iy0 * iz1;
		dest(x0,y1,z1) += val * ix0 * iy1 * iz1;
		dest(x1,y1,z1) += val * ix1 * iy1 * iz1;
	}
}
