#ifndef MANIFOLD_LOADER_H
#define MANIFOLD_LOADER_H

#include <memory>
#include <src/metadata_table.h>
#include "manifold.h"


class ManifoldLoader
{
	public:

		static std::vector<std::shared_ptr<Manifold>> loadAll(
				const MetaDataTable& table);

};

#endif
