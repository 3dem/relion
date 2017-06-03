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
#include <src/metadata_table.h>
#include <src/filename.h>

class star_compare_parameters
{
	public:
   	FileName fn1, fn2, fn_both, fn_only1, fn_only2, fn_label1, fn_label2, fn_label3;
	RFLOAT eps;
   	MetaDataTable MD1, MD2, MDonly1, MDonly2, MDboth;
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
	    fn1 = parser.getOption("--i1", "1st input STAR file ");
	    fn2 = parser.getOption("--i2", "2nd input STAR file ");
	    fn_both = parser.getOption("--both", "Output STAR file with entries from both input STAR files ", "");
	    fn_only1 = parser.getOption("--only1", "Output STAR file with entries that only occur in the 1st input STAR files ", "");
	    fn_only2 = parser.getOption("--only2", "Output STAR file with entries that only occur in the 2nd input STAR files ", "");
	    fn_label1 = parser.getOption("--label1", "1st metadata label for the comparison (may be string, int or RFLOAT)", "");
	    fn_label2 = parser.getOption("--label2", "2nd metadata label for the comparison (RFLOAT only) for 2D/3D-distance)", "");
	    fn_label3 = parser.getOption("--label3", "3rd metadata label for the comparison (RFLOAT only) for 3D-distance)", "");
	    eps = textToFloat(parser.getOption("--max_dist", "Maximum distance to consider a match (for int and RFLOAT only)", "0."));

      	// Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line, exiting...");
	}


	void run()
	{
		EMDLabel label1, label2, label3;
		MD1.read(fn1);
		MD2.read(fn2);

		label1 = EMDL::str2Label(fn_label1);
		label2 = (fn_label2 == "") ? EMDL_UNDEFINED : EMDL::str2Label(fn_label2);
		label3 = (fn_label3 == "") ? EMDL_UNDEFINED : EMDL::str2Label(fn_label3);


		compareMetaDataTable(MD1, MD2, MDboth, MDonly1, MDonly2, label1, eps, label2, label3);

		std::cout << MDboth.numberOfObjects()  << " entries occur in both input STAR files." << std::endl;
		std::cout << MDonly1.numberOfObjects() << " entries occur only in the 1st input STAR file." << std::endl;
		std::cout << MDonly2.numberOfObjects() << " entries occur only in the 2nd input STAR file." << std::endl;
		if (fn_both != "")
			MDboth.write(fn_both);
		if (fn_only1 != "")
			MDonly1.write(fn_only1);
		if (fn_only2 != "")
			MDonly2.write(fn_only2);

	}

};


int main(int argc, char *argv[])
{
	star_compare_parameters prm;

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

