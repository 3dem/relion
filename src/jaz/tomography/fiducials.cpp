#include "fiducials.h"
#include <omp.h>

using namespace gravis;


std::vector<d3Vector> Fiducials::read(
		const std::string &filename,
		double pixelSize)
{
	MetaDataTable metaDataTable;
	metaDataTable.read(filename);

	const int n = metaDataTable.numberOfObjects();

	std::vector<d3Vector> out(n);

	for (int i = 0; i < n; i++)
	{
		out[i].x = metaDataTable.getDouble(EMDL_ORIENT_ORIGIN_X_ANGSTROM, i) / pixelSize;
		out[i].y = metaDataTable.getDouble(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, i) / pixelSize;
		out[i].z = metaDataTable.getDouble(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, i) / pixelSize;
	}

	return out;
}

std::string Fiducials::write(
		const std::vector<gravis::d3Vector> &positions,
		double pixelSize,
		const std::string &tomoName,
		const std::string &path)
{
	MetaDataTable metaDataTable;

	for (int i = 0; i < positions.size(); i++)
	{
		const d3Vector posAngst = positions[i] * pixelSize;

		metaDataTable.addObject();

		metaDataTable.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, posAngst.x, i);
		metaDataTable.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, posAngst.y, i);
		metaDataTable.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, posAngst.z, i);
	}
	
	std::string filename = path + "fiducials_" + tomoName + ".star";

	metaDataTable.write(filename);
	
	return filename;
}

void Fiducials::drawMask(
		const std::vector<d3Vector> &positions,
		const d4Matrix &proj,
		double radius,
		RawImage<float> &destination,
		double value)
{
	const int w = destination.xdim;
	const int h = destination.ydim;

	const double fidRad2 = radius * radius;

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		for (int fid = 0; fid < positions.size(); fid++)
		{
			const d4Vector fid_img = proj * d4Vector(positions[fid]);

			const double dxf = fid_img.x - x;
			const double dyf = fid_img.y - y;
			const double distF2 = dxf * dxf + dyf * dyf;

			if (distF2 < fidRad2)
			{
				destination(x,y) = value;
			}
		}
	}
}

void Fiducials::erase(
        const std::vector<d3Vector>& positions, 
        double radius,
        Tomogram& tomogram, 
        int num_threads)
{
	const int w  = tomogram.stack.xdim;
	const int h  = tomogram.stack.ydim;
	const int fc = tomogram.stack.zdim;
	
	const double averaging_radius = 2 * radius;
	const double falloff = radius / 4;
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		BufferedImage<float> is_good(w,h);
		is_good.fill(1.f);
	
		drawMask(
			positions, tomogram.projectionMatrices[f],
			radius, is_good, 0.f);
		
		double frame_average = Normalization::computeMean(tomogram.stack.getSliceRef(f));
		
		for (int fid = 0; fid < positions.size(); fid++)
		{
			const d4Vector fid_img = tomogram.projectionMatrices[f] * d4Vector(positions[fid]);
			
			double sum_val = 0.0;
			double sum_wgh = 0.0;
			
			for (int dy = 0; dy < 2*averaging_radius; dy++)
			for (int dx = 0; dx < 2*averaging_radius; dx++)
			{
				const int x = (int) (fid_img.x - averaging_radius) + dx;
				const int y = (int) (fid_img.y - averaging_radius) + dy;
				
				const double dxf = x - fid_img.x;
				const double dyf = y - fid_img.y;
				
				const double r = sqrt(dxf*dxf + dyf*dyf);
				
				if (r > radius && r < averaging_radius && 
				    x >= 0 && x < w && 
				    y >= 0 && y < h)
				{
					const float m = is_good(x,y);
					
					sum_val += m * tomogram.stack(x,y,f);
					sum_wgh += m;
				}
			}
			
			const double avg_val = sum_wgh == 0.0? frame_average : sum_val / sum_wgh;
			
			for (int dy = 0; dy < 2*radius; dy++)
			for (int dx = 0; dx < 2*radius; dx++)
			{
				const int x = (int) (fid_img.x - radius) + dx;
				const int y = (int) (fid_img.y - radius) + dy;
				
				const double dxf = x - fid_img.x;
				const double dyf = y - fid_img.y;
				
				const double r = sqrt(dxf*dxf + dyf*dyf);
				
				if (r < radius && 
				    x >= 0 && x < w && 
				    y >= 0 && y < h)
				{
					double m;
					
					if (r < radius - falloff/2) 
					{
						m = 1;
					}
					else if (r < radius + falloff/2) 
					{
						m = 0.5 * (cos(PI * (r - radius + falloff/2) / falloff) + 1);
					}
					else
					{
						m = 0;
					}
					
					tomogram.stack(x,y,f) = m * avg_val + (1 - m) * tomogram.stack(x,y,f);
				}
			}
		}
	}
}
