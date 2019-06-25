/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres" "Jasenko Zivanov"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include <src/args.h>
#include <src/metadata_table.h>
#include <src/jaz/io/star_converter.h>

class star_converter
{
public:

	FileName fn_in, fn_out;
	IOParser parser;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("Options");
		fn_in = parser.getOption("--i", "Input STAR file to be converted", "None");
		fn_out = parser.getOption("--o", "Output STAR file to be written", "None");

		if (fn_in == "None" || fn_out == "None")
		{
			usage();
			REPORT_ERROR("Please specify input and output file names");
		}
	}

	void run()
	{
		MetaDataTable mdt;
		mdt.read(fn_in);
	
		MetaDataTable mdtOut, optOut;
		StarConverter::convert_3p0_particlesTo_3p1(mdt, mdtOut, optOut);
	
		std::ofstream of(fn_out);

		optOut.write(of);
		mdtOut.write(of);
		of.close();
	}
};

int main(int argc, char *argv[])
{
	star_converter app;

	try
	{
		app.read(argc, argv);
		app.run();
	}
	catch (RelionError XE)
	{
        	std::cerr << XE;
	        return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}
