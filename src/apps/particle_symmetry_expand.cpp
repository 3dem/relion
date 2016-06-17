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

#include <src/args.h>
#include <src/metadata_table.h>
#include <src/symmetries.h>
#include <src/euler.h>
#include <src/time.h>

class particle_symmetry_expand_parameters
{
public:

	FileName fn_sym, fn_in, fn_out;

	// I/O Parser
	IOParser parser;


	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("Options");

		fn_in = parser.getOption("--i", "Input particle STAR file");
		fn_out = parser.getOption("--o", "Output expanded particle STAR file", "expanded.star");
                fn_sym = parser.getOption("--sym", "Symmetry point group", "C1");

       	// Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	}

	void run()
	{

		MetaDataTable DFi, DFo;
		RFLOAT rot, tilt, psi;
		RFLOAT rotp, tiltp, psip;
		Matrix2D<RFLOAT> L(3,3), R(3,3); // A matrix from the list
		SymList SL;

		SL.read_sym_file(fn_sym);
		if (SL.SymsNo() < 1)
			REPORT_ERROR("ERROR Nothing to do. Provide a point group with symmetry!");


		DFi.read(fn_in);
		int barstep = XMIPP_MAX(1, DFi.numberOfObjects()/ 60);
		init_progress_bar(DFi.numberOfObjects());
                DFo.clear();

		long int imgno = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(DFi)
		{

			// Get the original line from the STAR file
			DFo.addObject();
			DFo.setObject(DFi.getObject());

			DFi.getValue(EMDL_ORIENT_ROT, rot);
			DFi.getValue(EMDL_ORIENT_TILT, tilt);
			DFi.getValue(EMDL_ORIENT_PSI, psi);

			for (int isym = 0; isym < SL.SymsNo(); isym++)
			{

				SL.get_matrices(isym, L, R);
                                L.resize(3, 3); // Erase last row and column
                                R.resize(3, 3); // as only the relative orientation is useful and not the translation
                                Euler_apply_transf(L, R, rot, tilt, psi, rotp, tiltp, psip);
				DFo.addObject();
				DFo.setObject(DFi.getObject());
				DFo.setValue(EMDL_ORIENT_ROT, rotp);
				DFo.setValue(EMDL_ORIENT_TILT, tiltp);
				DFo.setValue(EMDL_ORIENT_PSI, psip);

			}

                        if (imgno%barstep==0) progress_bar(imgno);
                        imgno++;

		} // end loop over input MetadataTable
		progress_bar(DFi.numberOfObjects());

		DFo.write(fn_out);
		std::cout << " Done! Written: " << fn_out << " with the expanded particle set." << std::endl;

	}// end run function



};

int main(int argc, char *argv[])
{
    time_config();
    particle_symmetry_expand_parameters prm;

    try
    {
        prm.read(argc, argv);
        prm.run();
    }

    catch (RelionError XE)
    {
        prm.usage();
        std::cerr << XE;
        exit(1);
    }

    return 0;

}
