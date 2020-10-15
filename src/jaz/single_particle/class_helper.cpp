#include "class_helper.h"
#include <src/jaz/util/index_sort.h>


int ClassHelper::countClasses(const MetaDataTable &particles_table)
{
	int max_class = -1;

	for (long int p = 0; p < particles_table.numberOfObjects(); p++)
	{
		const int class_id = particles_table.getIntMinusOne(EMDL_PARTICLE_CLASS, p);

		if (class_id > max_class) max_class = class_id;
	}

	return max_class + 1;
}

std::vector<int> ClassHelper::getClassSizes(const MetaDataTable &particles_table, int class_count)
{
	if (class_count < 0) 
	{
		class_count = countClasses(particles_table);
	}
	
	std::vector<int> particle_count(class_count, 0);

	for (long int p = 0; p < particles_table.numberOfObjects(); p++)
	{
		const int class_id = particles_table.getIntMinusOne(
					EMDL_PARTICLE_CLASS, p);

		particle_count[class_id]++;
	}
	
	return particle_count;
}

std::vector<int> ClassHelper::sortByAscendingFrequency(const std::vector<int>& particle_counts)
{
	std::vector<int> order = IndexSort<int>::sortIndices(particle_counts);
	
	const int n = order.size();
	
	std::vector<int> reverse(n);
	
	for (int i = 0; i < n; i++)
	{
		reverse[i] = order[n - i - 1];
	}
	
	return reverse;
}
