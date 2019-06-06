#include <src/metadata_table.h>
#include <src/jaz/io/star_converter.h>

int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		std::cerr << "usage: relion_convert_star <input> <output>\n";
	}
	
	MetaDataTable mdt;
	mdt.read(argv[1]);
	
	MetaDataTable mdtOut, optOut;
	
	StarConverter::convert_3p0_particlesTo_3p1(mdt, mdtOut, optOut);
	
	std::ofstream of(argv[2]);
	
	optOut.write(of);
	mdtOut.write(of);
			
	return RELION_EXIT_SUCCESS;
}
