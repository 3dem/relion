#ifndef MANIFOLD_SET_H
#define MANIFOLD_SET_H

#include <string>
#include <src/metadata_table.h>
#include "spheroid.h"

class TomogramManifoldSet
{
	public:
		
		TomogramManifoldSet();
		TomogramManifoldSet(const MetaDataTable& table);
		
			std::vector<Spheroid> spheroids;
			
		void addSpheroid(const Spheroid& spheroid);
		
		std::map<int, const Manifold*> getMapToManifolds() const;
		
		MetaDataTable composeTable() const;
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
