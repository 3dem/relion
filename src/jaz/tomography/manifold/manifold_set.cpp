#include "manifold_set.h"

using namespace gravis;


TomogramManifoldSet::TomogramManifoldSet()
{
}

TomogramManifoldSet::TomogramManifoldSet(const MetaDataTable &table)
{
	for (int i = 0; i < table.numberOfObjects(); i++)
	{
		const int index  = table.getInt(EMDL_TOMO_MANIFOLD_INDEX, i);
		const std::string type = table.getString(EMDL_TOMO_MANIFOLD_TYPE, i);
		const std::vector<double> parameters = table.getDoubleVector(EMDL_TOMO_MANIFOLD_PARAMETERS, i);
		
		if (type.length() > 0)
		{
			if (type == Spheroid::getTypeName())
			{
				addSpheroid(Spheroid(parameters, index));
			}
			else
			{
				REPORT_ERROR("TomogramManifoldSet::constructor: unknown manifold type: " + type);
			}
		}
	}
}

void TomogramManifoldSet::addSpheroid(const Spheroid& spheroid)
{
	spheroids.push_back(spheroid);
}

// @TODO: clean this up when we switch to C++11 (use unique pointers) 
//                        --JZ, 25-9-2020 (!)

std::map<int, const Manifold*> TomogramManifoldSet::getMapToManifolds() const
{
	std::map<int, const Manifold*> out;
	
	for (int i = 0; i < spheroids.size(); i++)
	{
		const Spheroid* s = &spheroids[i];
		
		out[s->index] = s;
	}
	
	return out;
}

MetaDataTable TomogramManifoldSet::composeTable() const
{
	MetaDataTable out;
	
	std::map<int, const Manifold*> manifolds = getMapToManifolds();
	
	for (std::map<int, const Manifold*>::iterator it = manifolds.begin();
		 it != manifolds.end(); it++)
	{
		const Manifold* m = it->second;

		out.addObject();
		
		out.setValue(EMDL_TOMO_MANIFOLD_INDEX, m->index);
		out.setValue(EMDL_TOMO_MANIFOLD_TYPE, m->type);
		out.setValue(EMDL_TOMO_MANIFOLD_PARAMETERS, m->getParameters());
	}
	
	return out;
}



ManifoldSet::ManifoldSet()
{
}

ManifoldSet::ManifoldSet(std::string filename)
{
	std::ifstream ifs(filename);

	if (!ifs)
	{
		REPORT_ERROR("ManifoldSet::constructor: Unable to read " + filename);
	}
	else
	{
		std::vector<MetaDataTable> allTables = MetaDataTable::readAll(ifs, 10);
		const int mc = allTables.size();
		
		for (int m = 0; m < mc; m++)
		{
			const std::string tomogramName = allTables[m].getName();
			add(tomogramName, TomogramManifoldSet(allTables[m]));
		}
	}
}

void ManifoldSet::add(
		const std::string& tomogramName,
		const TomogramManifoldSet& tomogramManifoldSet)
{
	perTomogramSets[tomogramName] = tomogramManifoldSet;
}

std::map<int, const Manifold*> ManifoldSet::getManifoldsInTomogram(
		const std::string& tomogramName) const
{
	std::map<std::string, TomogramManifoldSet>::const_iterator it
			= perTomogramSets.find(tomogramName);
	
	if (it == perTomogramSets.end())
	{
		REPORT_ERROR("ManifoldSet::getManifoldsInTomogram: no manifolds found for tomogram "
					 +tomogramName);
	}
	else
	{
		const TomogramManifoldSet& tms = it->second;
		return tms.getMapToManifolds();
	}
}

void ManifoldSet::write(std::string filename) const
{
	if (filename.find_last_of('/') != std::string::npos)
	{
		std::string path = filename.substr(0, filename.find_last_of('/'));
		mktree(path);
	}

	std::ofstream ofs(filename);
	
	if (!ofs)
	{
		REPORT_ERROR("ManifoldSet::write: unable to write to "+filename);
	}
	
	for (std::map<std::string, TomogramManifoldSet>::const_iterator it = perTomogramSets.begin();
		 it != perTomogramSets.end(); it++)
	{
		const std::string tomoName = it->first;
		const TomogramManifoldSet& tms = it->second;

		MetaDataTable table = tms.composeTable();
		table.setName(tomoName);
		
		table.write(ofs);
	}
}
