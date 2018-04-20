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
