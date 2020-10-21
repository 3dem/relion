#ifndef ABERRATIONS_CACHE_H
#define ABERRATIONS_CACHE_H

#include <src/metadata_table.h>
#include <src/jaz/image/buffered_image.h>


class AberrationsCache
{
	public:

		AberrationsCache(const MetaDataTable& opticsTable, const int size);

			int size;
			bool hasSymmetrical, hasAntisymmetrical;
			std::vector<BufferedImage<double>> symmetrical, antisymmetrical;
			std::vector<BufferedImage<fComplex>> phaseShift;
};

#endif
