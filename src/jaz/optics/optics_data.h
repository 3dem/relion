#ifndef DYN_OPTICS_DATA_H
#define DYN_OPTICS_DATA_H

#include <string>

class OpticsData
{
	public:
		
		OpticsData();
		OpticsData(std::string optFn, int verbosity = 1);
		
		double voltage, pixelSize, Cs;
};

#endif

