#include "manifold_loader.h"
#include "sphere.h"
#include "spheroid.h"
#include <memory>

std::vector<std::shared_ptr<Manifold> > ManifoldLoader::loadAll(const MetaDataTable &table)
{
	std::vector<std::shared_ptr<Manifold>> out;
	out.reserve(table.numberOfObjects());

	for (int i = 0; i < table.numberOfObjects(); i++)
	{
		const int index  = table.getInt(EMDL_TOMO_MANIFOLD_INDEX, i);
		const std::string type = table.getString(EMDL_TOMO_MANIFOLD_TYPE, i);
		const std::vector<double> parameters = table.getDoubleVector(EMDL_TOMO_MANIFOLD_PARAMETERS, i);

		if (type.length() > 0)
		{
			if (type == Spheroid::getTypeName())
			{
				out.push_back(std::shared_ptr<Spheroid>(new Spheroid(parameters, index)));
			}
			else if (type == Sphere::getTypeName())
			{
				out.push_back(std::shared_ptr<Sphere>(new Sphere(parameters, index)));
			}
			else
			{
				REPORT_ERROR("TomogramManifoldSet::constructor: unknown manifold type: " + type);
			}
		}
	}

	return out;
}
