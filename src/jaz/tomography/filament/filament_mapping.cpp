#include "filament_mapping.h"

void FilamentMapping::write(std::string filename)
{
	index.write(filename + "_index.mrc");
	signed_dist.write(filename + "_signed_dist.mrc");
	mask.write(filename + "_mask.mrc");
}

FilamentMapping FilamentMapping::read(std::string filename)
{
	FilamentMapping out;
	
	out.index.read(filename + "_index.mrc");
	out.signed_dist.read(filename + "_signed_dist.mrc");
	out.mask.read(filename + "_mask.mrc");
	
	return out;
}
