#ifndef STAR_CONVERTER_H
#define STAR_CONVERTER_H

#include <src/metadata_table.h>

class StarConverter
{
	public:
		
		static void convert_3p0_particlesTo_3p1(
						const MetaDataTable& in,
						MetaDataTable& outParticles,
						MetaDataTable& outOptics);
		
	protected:
		
		static void unifyPixelSize(MetaDataTable& outOptics);
		static void translateOffsets(MetaDataTable& outParticles, const MetaDataTable& optics);
};

#endif
