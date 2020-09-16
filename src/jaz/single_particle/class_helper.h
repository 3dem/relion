#ifndef SPA_CLASS_HELPER_H
#define SPA_CLASS_HELPER_H

#include <src/metadata_table.h>
#include <vector>

class ClassHelper
{
	public:
		
		static int countClasses(const MetaDataTable& particles_table);	
		static std::vector<int> getClassSizes(const MetaDataTable& particles_table, int class_count = -1);
		static std::vector<int> sortByAscendingFrequency(const std::vector<int>& particle_counts);
};

#endif
