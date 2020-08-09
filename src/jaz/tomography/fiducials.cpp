#include "fiducials.h"

using namespace gravis;


std::vector<d3Vector> Fiducials::read(
		const std::string &tomoName,
		const std::string &path)
{
	MetaDataTable metaDataTable;
	metaDataTable.read(path+"fiducials_"+tomoName+".star");

	const int n = metaDataTable.numberOfObjects();

	std::vector<d3Vector> out(n);

	for (int i = 0; i < n; i++)
	{
		out[i].x = metaDataTable.getDouble(EMDL_ORIENT_ORIGIN_X_ANGSTROM, i);
		out[i].y = metaDataTable.getDouble(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, i);
		out[i].z = metaDataTable.getDouble(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, i);
	}

	return out;
}

void Fiducials::write(
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

	metaDataTable.write(path+"fiducials_"+tomoName+".star");
}
