#include "manifold_set.h"
#include "manifold_loader.h"

using namespace gravis;


TomogramManifoldSet::TomogramManifoldSet()
{
}

TomogramManifoldSet::TomogramManifoldSet(const MetaDataTable &table)
{
	manifolds = ManifoldLoader::loadAll(table);
}

std::map<int, const Manifold*> TomogramManifoldSet::getMapToManifolds() const
{
	std::map<int, const Manifold*> out;
	
	for (int i = 0; i < manifolds.size(); i++)
	{
		const Manifold* m = manifolds[i].get();
		
		out[m->index] = m;
	}
	
	return out;
}

MetaDataTable TomogramManifoldSet::composeTable() const
{
	MetaDataTable out;
	
	for (int i = 0; i < manifolds.size(); i++)
	{
		const Manifold* m = manifolds[i].get();

		out.addObject();
		
		out.setValue(EMDL_TOMO_MANIFOLD_INDEX, m->index);
		out.setValue(EMDL_TOMO_MANIFOLD_TYPE, m->type);
		out.setValue(EMDL_TOMO_MANIFOLD_PARAMETERS, m->getParameters());
	}
	
	return out;
}

void TomogramManifoldSet::add(Manifold* manifold)
{
	manifolds.push_back(std::shared_ptr<Manifold>(manifold));
}


// --- //


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
		return std::map<int, const Manifold*>();
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
