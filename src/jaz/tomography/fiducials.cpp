#include "fiducials.h"

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
