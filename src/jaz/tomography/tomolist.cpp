#include "tomolist.h"
#include <fstream>
#include <sstream>
#include <src/error.h>

TomoList::TomoList()
{
	
}

TomoList::TomoList(std::string filename)
{
	std::vector<std::vector<std::string>> temp;
	std::vector<int> tempInd;
	
	std::ifstream ifs(filename);
	
	if (!ifs.is_open())
	{
		REPORT_ERROR_STR("TomoList::TomoList: unable to read " << filename << ".");
	}
	
	char line[4096];
	
	int maxInd = 0;
	
	while(ifs.getline(line, 4096))
	{
		std::istringstream iss(line);
				
		int ind;		
        iss >> ind;
		
        if (iss.fail()) continue;
		
		if (ind > maxInd) maxInd = ind;
				
		std::vector<std::string> tomoData(4);
		
		iss >> tomoData[0];
		iss >> tomoData[1];
		iss >> tomoData[2];
		iss >> tomoData[3];
				
		temp.push_back(tomoData);
		tempInd.push_back(ind);
	}
	
	table = std::vector<std::vector<std::string>>(maxInd+1, std::vector<std::string>(0));
	
	for (int i = 0; i < temp.size(); i++)
	{
		table[tempInd[i]] = temp[i];
	}
}

int TomoList::addTomogram(std::string ts, std::string proj, std::string ctf, std::string dose)
{
	table.push_back({ts, proj, ctf, dose});
	return table.size() - 1;
}

std::string TomoList::getTiltSeriesFilename(long int t) const
{
	return table[t][0];
}

std::string TomoList::getProjectionsFilename(long int t) const
{
	return table[t][1];
}

std::string TomoList::getCtfFilename(long int t) const
{
	return table[t][2];
}

std::string TomoList::getDoseFilename(long int t) const
{
	return table[t][3];
}

void TomoList::setTiltSeriesFilename(long t, std::string s)
{
	table[t][0] = s;
}

void TomoList::setProjectionsFilename(long int t, std::string s)
{
	table[t][1] = s;
}

void TomoList::setCtfFilename(long int t, std::string s)
{
	table[t][2] = s;
}

void TomoList::setDoseFilename(long int t, std::string s)
{
	table[t][3] = s;
}

bool TomoList::isKnown(long t)
{
	return t >= 0 && t < table.size() && table[t].size() > 0;
}

std::vector<long int> TomoList::getKnownIndices()
{
	std::vector<long int> out(0);
	out.reserve(table.size());
	
	for (int i = 0; i < table.size(); i++)
	{
		if (table[i].size() > 0)
		{
			out.push_back(i);
		}
	}
	
	return out;
}

void TomoList::write(std::string filename)
{
	const int tc = table.size();
	
	std::ofstream ofs(filename);
	
	for (int t = 0; t < tc; t++)
	{
		const int cs = table[t].size();
		
		ofs << t << " ";
		
		for (int c = 0; c < cs; c++)
		{
			ofs << table[t][c];
			
			if (c < cs-1)
			{
				ofs << ' ';
			}
		}
		
		ofs << '\n';
	}
}
