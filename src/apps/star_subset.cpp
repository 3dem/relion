/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
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

#include <src/image.h>
#include <src/filename.h>
#include <src/metadata_table.h>

class star_combine_parameters
{
	public:
	FileName fn_in, fn_out;
	std::string label, tablename;
	RFLOAT minval, maxval;
	// I/O Parser
	IOParser parser;


	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");
	    fn_in = parser.getOption("--i", "Input STAR file");
	    fn_out = parser.getOption("--o", "Output STAR file ");
	    label = parser.getOption("--label", "Metadata label to base selection on (e.g. rlnCtfFigureOfMerit)");
	    minval = textToFloat(parser.getOption("--minval", "Minimum acceptable value for this label", "-99999999."));
	    maxval = textToFloat(parser.getOption("--maxval", "Maximum acceptable value for this label", "99999999."));

      	// Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line, exiting...");
	}


	void run()
	{
		MetaDataTable MDin, MDout;

		if (fn_in.contains("_model.star"))
		{
			MDin.read(fn_in, "model_classes");
		}
		else
		{
			MDin.read(fn_in);
		}

		MDout = subsetMetaDataTable(MDin, EMDL::str2Label(label), minval, maxval);

		MDout.write(fn_out);
		std::cout << " Done! Written: " << fn_out << std::endl;
	}

};


int main(int argc, char *argv[])
{
	star_combine_parameters prm;

	try
    {

		prm.read(argc, argv);

		prm.run();

    }
    catch (RelionError XE)
    {
        std::cerr << XE;
        //prm.usage();
        exit(1);
    }
    return 0;
}



