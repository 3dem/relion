/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef INDEX_SORT_H
#define INDEX_SORT_H

#include <vector>
#include <algorithm>

template<class T>
class IndexSort
{
	public:

		static std::vector<int> sortIndices(const std::vector<T>& data)
		{
			const int s = (int) data.size();
			std::vector<int> indices(s);

			for (int i = 0; i < s; i++)
			{
				indices[i] = i;
			}

			sort(indices.begin(), indices.end(), IndexComparator(data));

			return indices;
		}

		struct IndexComparator
		{
			IndexComparator(const std::vector<T>& data) : data(data) {}

			bool operator()(const int a, const int b) const
			{
				return data[a] < data[b];
			}

			const std::vector<T>& data;
		};
};

#endif
