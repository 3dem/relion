#ifndef CMM_LOADER_H
#define CMM_LOADER_H

#include <src/jaz/gravis/t4Vector.h>
#include <vector>


class CMM_Loader
{
	public:

		static std::vector<gravis::d4Vector> readSpheres(
				const std::string& filename,
				double binning);
};

#endif
