#include "optics_data.h"
#include <src/error.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>
#include <iostream>
#include <fstream>

OpticsData::OpticsData()
{
    
}

OpticsData::OpticsData(std::string optFn, int verbosity)
{
	if (verbosity > 0)
	{
		Log::beginSection("Reading "+optFn);
	}
	
	std::ifstream ifs(optFn);
	
	if (!ifs.is_open())
	{
		REPORT_ERROR("OpticsData::OpticsData: unable to read optics file: " + optFn);
	}
	
	ifs >> pixelSize;
	ifs >> voltage;
	ifs >> Cs;
	
	if (verbosity > 0)
	{
		Log::print("pixel size: "+ZIO::itoa(pixelSize)+" Å");
		Log::print("voltage:    "+ZIO::itoa(voltage)+" kV");
		Log::print("C_s:        "+ZIO::itoa(Cs)+" mm");
		
		Log::endSection();
	}
}
