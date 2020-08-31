#ifndef Z_IO_TOOLS_H
#define Z_IO_TOOLS_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <src/error.h>

class ZIO
{
	public:
		
		template <typename T>
		static void writeDat(
			const std::vector<T>& vec, 
			std::string fn,
			long int begin = 0, 
			long int length = -1);
		
		template <typename T>
		static std::vector<T> readDat(std::string fn);

		static std::vector<double> readDoubles(std::string fn);
		static std::vector<int> readInts(std::string fn);
		static std::vector<std::vector<double>> readDoublesTable(std::string fn, int cols, char delim = ' ');
		static std::vector<std::string> split(const std::string& s, const std::string &delimiter);
		
		static std::string itoa(double num);
		
		static bool beginsWith(const std::string& string, const std::string& prefix);
		static bool endsWith(const std::string& string, const std::string& suffix);
		
		static std::string makeOutputDir(const std::string& dir);
		static std::string ensureEndingSlash(const std::string& dir);
};

template <typename T>
void ZIO::writeDat(
	const std::vector<T>& vec, std::string fn,
	long int begin, long int length)
{
	std::ofstream ofs(fn);
	
	const size_t end = length >= 0? begin + length : vec.size();
	
	for (int i = begin; i < end; i++)
	{
		ofs << i << " " << vec[i] << "\n";
	}
}

template <typename T>
std::vector<T> ZIO::readDat(std::string fn)
{
	std::ifstream ifs(fn);
	
	if (!ifs.is_open())
	{
		REPORT_ERROR_STR("ZIO::readDat: unable to read " << fn);
	}
	
	std::string line;
	
	std::vector<T> out;
	out.reserve(128);
	
	while (std::getline(ifs, line))
	{
		std::stringstream sts;
		sts << line;
		
		double d;
		sts >> d;
		
		sts << line;
		
		sts >> d;
		
		out.push_back(d);
	}
	
	return out;
}


#endif
