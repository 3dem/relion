
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

#include <src/jaz/single_particle/obs_model.h>
#include <src/image.h>
#include <src/metadata_table.h>

class import_parameters
{
	public:
   	FileName fn_in, fn_odir, fn_out, fn_mtf;
	bool do_write_types, do_continue, do_movies, do_micrographs, do_coordinates, do_halfmaps, do_particles, do_other;
   	FileName optics_group_name, node_type, particles_optics_group_name;
   	RFLOAT kV, Cs, Q0, beamtilt_x, beamtilt_y, pixel_size;

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
		fn_in = parser.getOption("--i", "Input (wildcard) filename");
		fn_odir = parser.getOption("--odir", "Output directory (e.g. \"Import/job001/\"");
		fn_out = parser.getOption("--ofile", "Output file name (e.g. \"movies.star\"");
		do_movies = parser.checkOption("--do_movies", "Import movies");
		do_micrographs = parser.checkOption("--do_micrographs", "Import micrographs");
		do_coordinates = parser.checkOption("--do_coordinates", "Import coordinates");
		do_halfmaps = parser.checkOption("--do_halfmaps", "Import unfiltered half maps");
		do_particles = parser.checkOption("--do_particles", "Import particle STAR files");
		particles_optics_group_name = parser.getOption("--particles_optics_group_name", "Rename optics group for all imported particles (e.g. \"opticsGroupLMBjan2019\"", "");
		do_other = parser.checkOption("--do_other", "Import anything else");

		int mic_section = parser.addSection("Specific options for movies or micrographs");
		optics_group_name = parser.getOption("--optics_group_name", "Name for this optics group", "opticsGroup1");
		fn_mtf = parser.getOption("--optics_group_mtf", "Name for this optics group's MTF", "");
		pixel_size = textToFloat(parser.getOption("--angpix", "Pixel size (Angstrom)", "1.0"));
		kV = textToFloat(parser.getOption("--kV", "Voltage (kV)", "300"));
		Cs = textToFloat(parser.getOption("--Cs", "Spherical aberration (mm)", "2.7"));
		Q0 = textToFloat(parser.getOption("--Q0", "Amplitude contrast", "0.1"));
		beamtilt_x = textToFloat(parser.getOption("--beamtilt_x", "Beam tilt (X; mrad)", "0.0"));
		beamtilt_y = textToFloat(parser.getOption("--beamtilt_y", "Beam tilt (Y; mrad)", "0.0"));
		do_continue = parser.checkOption("--continue", "Continue and old run, add more files to the same import directory");

		// Check for errors in the command-line option
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

		if (pixel_size <= 0)
			REPORT_ERROR("Pixel size must be positive!");

		if (kV <= 0)
			REPORT_ERROR("Acceleration voltage must be positive!");
	}

	void run()
	{
		std::string command;
		MetaDataTable MDout, MDopt;
		std::vector<FileName> fns_in;
		long nr_input_files = fn_in.globFiles(fns_in);
		std::ofstream  fh;

		int nr_count = 0;
		if (do_movies) nr_count++;
		if (do_micrographs) nr_count++;
		if (do_coordinates) nr_count++;
		if (do_other || do_halfmaps || do_particles) nr_count++;
		if (nr_count != 1)
		{
			REPORT_ERROR("ERROR: you can only use only one, and at least one, of the options --do_movies, --do_micrographs, --do_coordinates, --do_halfmaps or --do_other");
		}

		std::cout << " importing..." << std::endl;
		// For micrographs or movies
		if (do_movies || do_micrographs)
		{
			if (fn_in.rfind("../") != std::string::npos) // Forbid at any place
				REPORT_ERROR("Please don't import files outside the project directory.\nPlease make a symbolic link by an absolute path before importing.");

			if (fn_in.rfind("/", 0) == 0) // Forbid only at the beginning
				REPORT_ERROR("Please import files by a relative path.\nIf you want to import files outside the project directory, make a symbolic link by an absolute path and\nimport the symbolic link by a relative path.");

			std::string tablename = (do_movies) ? "movies" : "micrographs";
			bool do_new_optics_group = true;
			int old_optics_group_number, optics_group_number = 1;
			long old_nr_files = 0;

			// When continuing old jobs in the pipeliner, the old names are moved out of the way. Read it in anyway!
			FileName old_fn_out = fn_odir + fn_out;
			if (do_continue && exists(old_fn_out))
			{
				MDopt.read(old_fn_out, "optics");
				MDout.read(old_fn_out, tablename);
				old_nr_files = MDout.numberOfObjects();
				std::string old_optics_group_name;
				FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDopt)
				{
					MDopt.getValue(EMDL_IMAGE_OPTICS_GROUP_NAME, old_optics_group_name);
					if (old_optics_group_name == optics_group_name)
					{
						do_new_optics_group = false;
						MDopt.getValue(EMDL_IMAGE_OPTICS_GROUP, optics_group_number);
						break;
					}
				}
				if (do_new_optics_group)
				{
					optics_group_number = MDopt.numberOfObjects() + 1;
				}
			}

			if (do_new_optics_group)
			{
				if (!optics_group_name.validateCharactersStrict())
					REPORT_ERROR("The optics group name may contain only numbers, alphabets and hyphen(-).");

				// Generate MDopt for the optics group
				MDopt.setName("optics");
				MDopt.addObject();
				MDopt.setValue(EMDL_IMAGE_OPTICS_GROUP_NAME, optics_group_name);
				MDopt.setValue(EMDL_IMAGE_OPTICS_GROUP, optics_group_number);
				if (fn_mtf != "") MDopt.setValue(EMDL_IMAGE_MTF_FILENAME, fn_mtf);
				if (do_micrographs) MDopt.setValue(EMDL_MICROGRAPH_PIXEL_SIZE, pixel_size);
				MDopt.setValue(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, pixel_size);
				MDopt.setValue(EMDL_CTF_VOLTAGE, kV);
				MDopt.setValue(EMDL_CTF_CS, Cs);
				MDopt.setValue(EMDL_CTF_Q0, Q0);
				if (fabs(beamtilt_x) + fabs(beamtilt_y) > 0.001)
				{
					MDopt.setValue(EMDL_IMAGE_BEAMTILT_X, beamtilt_x);
					MDopt.setValue(EMDL_IMAGE_BEAMTILT_Y, beamtilt_y);
				}
			}

			// Fill in the actual data (movies/micrographs) table
			MDout.setName(tablename);
			EMDLabel mylabel = (do_movies) ? EMDL_MICROGRAPH_MOVIE_NAME : EMDL_MICROGRAPH_NAME;
			for (long i = 0; i < nr_input_files; i++)
			{
				// Check this file was not yet present in the input STAR file
				// TODO: this N^2 algorithm might get too expensive with large data sets....
				bool already_there = false;
				for (long j = 0; j < old_nr_files; j++)
				{
					FileName oldfile;
					MDout.getValue(mylabel, oldfile, j);
					if (oldfile == fns_in[i])
					{
						already_there = true;
						int old_optics_group_number;
						MDout.getValue(EMDL_IMAGE_OPTICS_GROUP, old_optics_group_number, j);
						if (old_optics_group_number != optics_group_number)
						{
							std::cerr << " fns_in[i]= " << fns_in[i] << " old_optics_group_number= " << old_optics_group_number << " optics_group_number= " << optics_group_number << std::endl;
							REPORT_ERROR("ERROR: trying to add an pre-existing image with a different optics group!");
						}
						break;
					}
				}
				if (!already_there)
				{
					MDout.addObject();
					MDout.setValue(mylabel, fns_in[i]);
					MDout.setValue(EMDL_IMAGE_OPTICS_GROUP, optics_group_number);
				}
			}

			// Write output STAR file
			fh.open((fn_odir + fn_out).c_str(), std::ios::out);
			if (!fh)
				REPORT_ERROR( (std::string)"MlModel::write: Cannot write file: " + fn_odir + fn_out);

			MDopt.write(fh);
			MDout.write(fh);
			fh.close();

			long nr_new_files = MDout.numberOfObjects();
			std::cout << " Written " << (fn_odir + fn_out) << " with " << nr_new_files << " items (" << (nr_new_files - old_nr_files) << " new items)" << std::endl;
		}
		else if (do_coordinates)
		{
			// Make the same directory structure of the coordinates
			// Copy all coordinate files into the same subdirectory in the Import directory
			// But remove directory structure from pipeline if that exists
			// Dereference symbolic links if needed
			FileName fn_dir = fn_in;
			if (fn_dir.contains("/"))
				fn_dir = fn_dir.beforeLastOf("/");
			else
				fn_dir = ".";
			FileName fn_pre, fn_jobnr, fn_post;
			if (decomposePipelineSymlinkName(fn_dir, fn_pre, fn_jobnr, fn_post))
			{
				// Make the output directory
				command = "mkdir -p " + fn_odir + fn_post;
				if (system(command.c_str())) REPORT_ERROR("ERROR: there was an error executing: " + command);
				// Copy the coordinate files one by one to prevent problems of too long command line
				for (long i = 0; i < nr_input_files; i++)
				{
					command = "cp " + fns_in[i] + " " + fn_odir + fn_post;
					if (system(command.c_str())) REPORT_ERROR("ERROR: there was an error executing: " + command);
				}
			}
			else
			{
				// Copy the coordinate files one by one to prevent problems of too long command line
				for (long i = 0; i < nr_input_files; i++)
				{
					command = "cp --parents " + fns_in[i] + " " + fn_odir;
					if (system(command.c_str())) REPORT_ERROR("ERROR: there was an error executing: " + command);
				}
			}

			// Make a suffix file, which contains the actual suffix as a suffix
			// Get the coordinate-file suffix
			FileName fn_suffix2 = fn_in.beforeLastOf("*");
			fh.open((fn_odir + fn_out).c_str(), std::ios::out);
			if (!fh)
				REPORT_ERROR( (std::string)"Import: Cannot write file: " + fn_odir + fn_out);

			fh << fn_suffix2 << "*.mrc" << std::endl;
			fh.close();
		}
		else if (do_particles)
		{
			ObservationModel obsModel;
			MetaDataTable MD;
			ObservationModel::loadSafely(fn_in, obsModel, MD);

			// Make sure rlnOpticsGroupName is set to this value
			// This is only a valid option if there was a single optics_group in the input file
			if (particles_optics_group_name != "")
			{
				if (!particles_optics_group_name.validateCharactersStrict())
					REPORT_ERROR("The optics group name may contain only numbers, alphabets and hyphen(-).");

				if (obsModel.opticsMdt.numberOfObjects() != 1)
				{
					obsModel.opticsMdt.write(std::cerr);
					REPORT_ERROR(" ERROR: cannot rename particles optics groups when multiple ones are imported!");
				}
				obsModel.opticsMdt.setValue(EMDL_IMAGE_OPTICS_GROUP_NAME, particles_optics_group_name, 0);
			}

			FileName fnt = "/" + fn_in;
			fnt = fn_odir + fnt.afterLastOf("/");
			obsModel.save(MD, fnt, "particles");
		}
		else if (do_other || do_halfmaps)
		{

			if (nr_input_files > 1)
			{
				REPORT_ERROR("ERROR: Multiple files (i.e. filename wildcards) are not allowed for the import of other types.");
			}

			// For all the rest of the imports, just copy the files in the Import/jobxxx/ directory with the same name
			FileName fnt = "/" + fn_in;
			fnt = fnt.afterLastOf("/");
			command = "cp " + fn_in + " " + fn_odir + fnt;
			if (system(command.c_str())) REPORT_ERROR("ERROR: there was an error executing: " + command);

			if (do_halfmaps)
			{

				// For unfiltered half-maps, also get the other half-map
				FileName fn_inb = fn_in;
				size_t pos = fn_inb.find("half1");
				if (pos != std::string::npos)
				{
					fn_inb.replace(pos, 5, "half2");
				}
				else
				{
					pos = fn_inb.find("half2");
					if (pos != std::string::npos)
					{
						fn_inb.replace(pos, 5, "half1");
					}
				}
				fnt = "/" + fn_inb;
				fnt = fnt.afterLastOf("/");
				command = "cp " + fn_inb + " " + fn_odir + fnt;
				if (system(command.c_str())) REPORT_ERROR("ERROR: there was an error executing: " + command);
			}
		}
		std::cout << " done!" << std::endl;
	}
};

int main(int argc, char *argv[])
{
	import_parameters prm;

	try
	{
		prm.read(argc, argv);

		prm.run();

	}
	catch (RelionError XE)
	{
			//prm.usage();
			std::cerr << XE;
			return RELION_EXIT_FAILURE;
	}
	return RELION_EXIT_SUCCESS;
}
