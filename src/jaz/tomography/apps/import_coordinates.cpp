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
#include <src/jaz/tomography/imod_import.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>
#include <src/ctf.h>

bool readCoordinateFile(FileName fn_coord, MetaDataTable &MD,  bool is_centered,
                        RFLOAT rescale_factor, RFLOAT add_factor, std::string tomoname="None")
{
    MD.clear();

    // STAR files are read in directly
    if (fn_coord.getExtension() == "star")
    {
        MD.read(fn_coord);

        // Set the tomoname as well, if it was not yet present
        if (!MD.containsLabel(EMDL_TOMO_NAME))
        {
            if (tomoname == "None")
                REPORT_ERROR("readCoordinateFile ERROR: File " + fn_coord + " does not contain rlnTomoName, but no tomoname was provided...");

            FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
            {
                MD.setValue(EMDL_TOMO_NAME, tomoname);
            }
        }
        else if (tomoname != "None")
        {
            // Just check the first tomoname in the input STAR file is the same as the provided one
            std::string myname;

            MD.getValue(EMDL_TOMO_NAME, myname, 0);
            if (myname != tomoname)
                REPORT_ERROR("readCoordinateFile ERROR: File " + fn_coord + " contains tomonames " + myname + " that is not compatible with provided name of " + tomoname);
        }

    }
    else // assume 3 or 6-column text file. Lines that do not start with 3 numbers are ignored
    {
        std::ifstream in(fn_coord.data(), std::ios_base::in);
        if (in.fail())
            REPORT_ERROR((std::string) "readCoordinateFile ERROR: File " + fn_coord + " does not exists");

        // Start reading the ifstream at the top
        in.seekg(0);
        std::string line;
        int n = 0;
        while (getline(in, line, '\n'))
        {
            std::vector<std::string> words;
            tokenize(line, words);

            float num1, num2, num3, num4, num5, num6;

            // Ignore lines that do not have at least two integer numbers on it (at this point I do not know dimensionality yet....)
            if (words.size() >= 3 &&
                sscanf(words[0].c_str(), "%f", &num1) &&
                sscanf(words[1].c_str(), "%f", &num2) &&
                sscanf(words[2].c_str(), "%f", &num3))
            {
                MD.addObject();
                MD.setValue(EMDL_TOMO_NAME, tomoname);

                if (is_centered)
                {

                    MD.setValue(EMDL_IMAGE_CENT_COORD_X_ANGST, num1*rescale_factor + add_factor);
                    MD.setValue(EMDL_IMAGE_CENT_COORD_Y_ANGST, num2*rescale_factor + add_factor);
                    MD.setValue(EMDL_IMAGE_CENT_COORD_Z_ANGST, num3*rescale_factor + add_factor);

                }
                else
                {

                    MD.setValue(EMDL_IMAGE_COORD_X, num1*rescale_factor + add_factor);
                    MD.setValue(EMDL_IMAGE_COORD_Y, num2*rescale_factor + add_factor);
                    MD.setValue(EMDL_IMAGE_COORD_Z, num3*rescale_factor + add_factor);

                }

                // Check if prior angles are also defined as numbers 3-5
                if (words.size() >= 6 &&
                    sscanf(words[3].c_str(), "%f", &num4) &&
                    sscanf(words[4].c_str(), "%f", &num5) &&
                    sscanf(words[5].c_str(), "%f", &num6))
                {

                    MD.setValue(EMDL_TOMO_SUBTOMOGRAM_ROT, num4);
                    MD.setValue(EMDL_TOMO_SUBTOMOGRAM_TILT, num5);
                    MD.setValue(EMDL_TOMO_SUBTOMOGRAM_PSI, num6);
                    MD.setValue(EMDL_ORIENT_ROT, 0.);
                    MD.setValue(EMDL_ORIENT_TILT, 90.);
                    MD.setValue(EMDL_ORIENT_PSI, 0.);
                    MD.setValue(EMDL_ORIENT_TILT_PRIOR, 90.);
                    MD.setValue(EMDL_ORIENT_PSI_PRIOR, 0.);
                }

            }
        }
        in.close();
    }

    return (MD.containsLabel(EMDL_TOMO_SUBTOMOGRAM_ROT) && MD.containsLabel(EMDL_TOMO_SUBTOMOGRAM_TILT) && MD.containsLabel(EMDL_TOMO_SUBTOMOGRAM_PSI));
}

int main(int argc, char *argv[])
{
	try
	{
		IOParser parser;

		FileName inStarFn, tomoFn, outDir;
		bool is_centered;
        RFLOAT scale_factor, add_factor;
        std::string remove_substring, remove_substring2;

		try
		{
			parser.setCommandLine(argc, argv);

			int gen_section = parser.addSection("General options");

			inStarFn = parser.getOption("--i", "Linux wildcard for set of coordinate files (in STAR or 3/6-column ASCII), or a STAR file with a list of coordinate files and their tomogram name");
			outDir = parser.getOption("--o", "Name of the output directory");
            remove_substring = parser.getOption("--remove_substring", "Remove this substring from the coordinate filenames to get the tomogram names", "");
            remove_substring2 = parser.getOption("--remove_substring2", "Remove also this substring from the coordinate filenames to get the tomogram names", "");
            is_centered = parser.checkOption("--centered", "Are the coordinates relative to the center of the tomogram?");
            scale_factor = textToFloat(parser.getOption("--scale_factor", "Multiply coordinates with this factor", "1.0"));
            add_factor =  textToFloat(parser.getOption("--add_factor", "Add this value to final coordinates", "0.0"));
			Log::readParams(parser);

			if (parser.checkForErrors())
			{
				exit(RELION_EXIT_FAILURE);
			}
		}
		catch (RelionError XE)
		{
			parser.writeUsage(std::cout);
			std::cerr << XE;
			exit(RELION_EXIT_FAILURE);
		}

        // Make sure outDir ends with a slash
        if (outDir[outDir.length()-1] != '/')
            outDir += "/";

		outDir = ZIO::makeOutputDir(outDir);

        MetaDataTable MDin, MDout;
        std::vector<FileName> fns_partfiles;
        // See if the input is a linux wildcard
        if (inStarFn.globFiles(fns_partfiles) == 1)
        {
            MDin.read(inStarFn);
        }
        else
        {
            MDin.clear();

            // Read in sets of filenames and get the tomogram name from the filenames before their extension
            std::vector<FileName> fns_partfiles;
            inStarFn.globFiles(fns_partfiles);
            for (size_t i = 0; i < fns_partfiles.size(); i++)
            {
                // If the partfile is a STAR file, see if it contains a rlnTomoName and use that one
                MetaDataTable MDpart;
                FileName tomoname = "None";
                if (!(MDpart.read(fns_partfiles[i]) && MDpart.containsLabel(EMDL_TOMO_NAME)))
                {
                    tomoname = (fns_partfiles[i].afterLastOf("/")).beforeLastOf(".");
                    tomoname.replaceAllSubstrings(remove_substring, "");
                    tomoname.replaceAllSubstrings(remove_substring2, "");
                }
                MDin.addObject();
                MDin.setValue(EMDL_TOMO_NAME, tomoname);
                MDin.setValue(EMDL_TOMO_IMPORT_PARTICLE_FILE, fns_partfiles[i]);
            }

        }

        if (!MDin.containsLabel(EMDL_TOMO_NAME))
            REPORT_ERROR("Input star file " + inStarFn + " does not contain tomogram names (rlnTomoName)");
        if (!MDin.containsLabel(EMDL_TOMO_IMPORT_PARTICLE_FILE))
            REPORT_ERROR("Input star file " + inStarFn + " does not contain coordinate file names (rlnTomoImportParticleFile)");

        bool some_have_angles = false;
        bool all_have_angles = true;
        FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
        {
            std::string tomoname, partfile;
            MetaDataTable MDpart;

            MDin.getValue(EMDL_TOMO_NAME, tomoname);
            MDin.getValue(EMDL_TOMO_IMPORT_PARTICLE_FILE, partfile);

            bool hasAngles = readCoordinateFile(partfile, MDpart, is_centered, scale_factor, add_factor, tomoname);
            if (hasAngles) some_have_angles = true;
            else all_have_angles = false;

            MDout.append(MDpart);
        }

        if (some_have_angles && !all_have_angles)
            REPORT_ERROR("ERROR: some, but not all coordinate files contain prior angle information...");


        MDout.setName("particles");
        MDout.write(outDir + "particles.star");
        std::cout <<" Done! Imported coordinates as output " << outDir << "particles.star" << std::endl;
        if (all_have_angles)
            std::cout << "This file contains also prior angle information" << std::endl;

	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}
