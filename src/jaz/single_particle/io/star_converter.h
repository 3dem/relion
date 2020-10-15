#ifndef STAR_CONVERTER_H
#define STAR_CONVERTER_H

#include <src/metadata_table.h>
#include <src/image.h>

class StarConverter
{
	public:

		static void convert_3p0_particlesTo_3p1(
						const MetaDataTable& in,
						MetaDataTable& outParticles,
						MetaDataTable& outOptics,
						std::string tablename = "particles",
						bool do_die_upon_error = true);

	protected:

		static void unifyPixelSize(MetaDataTable& outOptics, std::string tablename = "particles");
		static void translateOffsets(MetaDataTable& outParticles, const MetaDataTable& optics);
};

#endif
