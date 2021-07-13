#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <src/metadata_table.h>


int main(int argc, char *argv[])
{
	std::map<std::string, std::map<int, std::string>> movies_by_tomo;

	for (int i = 1; i < argc; i++)
	{
		std::string arg = std::string(argv[i]);
		std::string arg_no_path = arg.substr(arg.find_last_of('/')+1);

		int tomo_index = -1;
		int tilt_index = -1;

		if (!sscanf(arg_no_path.c_str(), "TS_%d_%d_%*s", &tomo_index, &tilt_index))
		{
			std::cerr << "failed to extract indices from " << argv[i] << std::endl;
		}
		else
		{
			std::string tomo_name = arg_no_path.substr(0, arg_no_path.find_first_of('_', 3));

			movies_by_tomo[tomo_name][tilt_index] = arg;
		}
	}

	std::ofstream out_stream("tilt_frames.star");

	for (auto it = movies_by_tomo.begin(); it != movies_by_tomo.end(); it++)
	{
		MetaDataTable table;

		std::cout << it->first << std::endl;

		table.setName(it->first);

		for (auto itt = it->second.begin(); itt != it->second.end(); itt++)
		{
			std::cout << "    " << itt->first << "    " << itt->second << std::endl;

			table.addObject();
			table.setValue(EMDL_TOMO_TILT_MOVIE_INDEX, itt->first);
			table.setValue(EMDL_TOMO_TILT_MOVIE_FILE_NAME, itt->second);
		}

		table.write(out_stream);
	}

	out_stream.flush();

	return 0;
}
