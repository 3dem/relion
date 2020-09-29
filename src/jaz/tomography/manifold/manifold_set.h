#ifndef MANIFOLD_SET_H
#define MANIFOLD_SET_H

#include <memory>
#include <string>
#include <src/metadata_table.h>
#include "manifold.h"

class TomogramManifoldSet
{
	public:
		
		TomogramManifoldSet();
		TomogramManifoldSet(const MetaDataTable& table);
		
			std::vector<std::shared_ptr<Manifold>> manifolds;
		
		std::map<int, const Manifold*> getMapToManifolds() const;
		MetaDataTable composeTable() const;

		void add(Manifold* manifold);
};

class ManifoldSet
{
	public:	
		
		ManifoldSet();
		ManifoldSet(std::string filename);
				
			std::map<std::string, TomogramManifoldSet> perTomogramSets;
			
			
		void add(
				const std::string& tomogramName,
				const TomogramManifoldSet& tomogramManifoldSet);

		std::map<int, const Manifold*> getManifoldsInTomogram(
				const std::string& tomogramName) const;
		
		void write(std::string filename) const;
		
};

#endif
