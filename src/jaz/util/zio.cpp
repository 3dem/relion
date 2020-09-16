#include "zio.h"
#include <string.h>
#include <stdio.h>


std::vector<double> ZIO::readDoubles(std::string fn)
{
	std::vector<double> out(0);
	std::ifstream ifs(fn);	
	std::string line;
	
	if (!ifs)
	{
		REPORT_ERROR("ZIO::readDoubles: unable to read "+fn);
	}
	
	out.reserve(128);
	
	while (std::getline(ifs, line))
	{
		std::stringstream sts;
		sts << line;
		
		double d;
		sts >> d;
		
		out.push_back(d);
	}
	
	return out;
}

std::vector<int> ZIO::readInts(std::string fn)
{
	std::vector<int> out(0);
	std::ifstream ifs(fn);
	
	if (!ifs)
	{
		REPORT_ERROR("ZIO::readInts: unable to read "+fn);
	}
	
	std::string line;

	out.reserve(128);

	while (std::getline(ifs, line))
	{
		std::stringstream sts;
		sts << line;

		int i;
		sts >> i;

		out.push_back(i);
	}

	return out;
}

std::vector<std::vector<double>> ZIO::readFixedDoublesTable(std::string fn, int cols, char delim)
{
	std::vector<std::vector<double>> out(0);
	
	std::ifstream ifs(fn);

	if (!ifs)
	{
		REPORT_ERROR_STR("ZIO::readFixedDoublesTable: unable to read "+fn);
	}
	
	std::string line;
	
	out.reserve(128);
	
	int i = -1;
	
	while (std::getline(ifs, line))
	{
		out.push_back(std::vector<double>(cols, 0.0));
		i++;

		if (delim != ' ')
		{
			for (int i = 0; i < line.length(); i++)
			{
				if (line[i] == delim)
				{
					line[i] = ' ';
				}
			}
		}
		
		std::stringstream sts;
		sts << line;
		
		for (int c = 0; c < cols; c++)
		{
			double d;
			sts >> d;
			
			out[i][c] = d;
		}
	}
	
	return out;
}

std::vector<std::vector<double> > ZIO::readDoublesTable(std::string fn, char delim)
{
	std::vector<std::vector<double>> out(0);
	
	std::ifstream ifs(fn);

	if (!ifs)
	{
		REPORT_ERROR_STR("ZIO::readFixedDoublesTable: unable to read "+fn);
	}
	
	std::string line;
	
	out.reserve(128);
	
	int i = -1;
	
	while (std::getline(ifs, line))
	{
		out.push_back(std::vector<double>());
		i++;

		if (delim != ' ')
		{
			for (int i = 0; i < line.length(); i++)
			{
				if (line[i] == delim)
				{
					line[i] = ' ';
				}
			}
		}
		
		out[i].reserve(line.length()/5 + 1);
		
		std::stringstream sts;
		sts << line;
				
		double d;
		
		while (sts >> d)
		{
			out[i].push_back(d);
			std::cout << "[" << d << "]\n";
		}
	}
	
	return out;
}

std::vector<std::string> ZIO::split(const std::string &s, const std::string &delimiter)
{
	std::vector<std::string> out;
	
	char* cstr = const_cast<char*>(s.c_str());
	char* current = strtok(cstr, delimiter.c_str());
	
	while (current != NULL)
	{
		out.push_back(current);
		current = strtok(0, delimiter.c_str());
	}
	
	return out;
}

std::string ZIO::itoa(double num)
{
	std::ostringstream sts;
	sts << num;
	
	return sts.str();
}

bool ZIO::beginsWith(const std::string &string, const std::string &prefix)
{
	return string.length() >= prefix.length() 
			&& string.substr(0, prefix.length()) == prefix;
}

bool ZIO::endsWith(const std::string &string, const std::string &prefix)
{
	return string.length() >= prefix.length() 
			&& string.substr(string.length() - prefix.length()) == prefix;
}

std::string ZIO::makeOutputDir(const std::string& dir)
{
	std::string out = dir;
	
	const int len = out.length();
	
	if (len > 0)
	{
		if (out[len-1] != '/')
		{
			out = out + "/";
		}
		
		int res = system(("mkdir -p "+out).c_str());

		if (res)
		{
			REPORT_ERROR_STR("Unable to write to: " << out);
		}
	}
	
	return out;
}

std::string ZIO::ensureEndingSlash(const std::string &dir)
{
	std::string out = dir;
	
	const int len = out.length();
	
	if (len > 0)
	{
		if (out[len-1] != '/')
		{
			out = out + "/";
		}
	}
	
	return out;	
}
