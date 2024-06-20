/***************************************************************************
 *
 * Author: "Dari Kimanius"
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
#include <src/error.h>
#include <src/time.h>
#include "src/exp_model.h"

// TODO: set pixel sizes in the outputs

class pseudotomo_expansion
{
	public:
   	FileName fn_part, fn_tomo, fn_motion, fn_out;
	IOParser parser;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");
        fn_part = parser.getOption("--particles", "Input particles STAR file");
        fn_tomo = parser.getOption("--tomograms", "Input tomography STAR file");
        fn_motion = parser.getOption("--motion", "Input motions STAR file");
        fn_out = parser.getOption("--o", "Output STAR file");

		// Check for errors in the command-line option
		if (parser.checkForErrors())
    			REPORT_ERROR("Errors encountered on the command line, exiting...");
	}

	void run()
    {
        Experiment data;
        MetaDataTable optics, input_part, md;

        data.read(fn_part, fn_tomo, fn_motion);

        if (!data.is_tomo)
            throw std::runtime_error("Input not Tomography dataset");

        ObservationModel obsModel;
        ObservationModel::loadSafely(fn_part, obsModel, input_part, "particles");

        md.addMissingLabels(&input_part);

        md.addLabel(EMDL_CTF_DEFOCUSU);
        md.addLabel(EMDL_CTF_DEFOCUSV);
        md.addLabel(EMDL_CTF_DEFOCUS_ANGLE);
        md.addLabel(EMDL_MICROGRAPH_PRE_EXPOSURE);

        std::vector<EMDLabel> labels = md.getActiveLabels();

        for (int part_id = 0; part_id < data.particles.size(); part_id ++)
        {
            size_t id = data.particles[part_id].id;
            MetaDataContainer* rowToCopy = input_part.getObject(id);

            RFLOAT shift3d_x, shift3d_y, shift3d_z, rot, tilt, psi;
            Matrix2D<RFLOAT> A;
            std::string img_fn;

            input_part.getValue(EMDL_IMAGE_NAME, img_fn, id);

            // Rotation
            input_part.getValue(EMDL_ORIENT_ROT, rot, id);
            input_part.getValue(EMDL_ORIENT_TILT, tilt, id);
            input_part.getValue(EMDL_ORIENT_PSI, psi, id);
            Euler_angles2matrix(rot, tilt, psi, A);

            // Translation
            input_part.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, shift3d_x, id);
            input_part.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, shift3d_y, id);
            input_part.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, shift3d_z, id);

            for (int img_id = 0; img_id < data.particles[part_id].images.size(); img_id ++)
            {
                md.addObject(rowToCopy);
                md.setValue(EMDL_IMAGE_NAME, std::to_string(img_id) + "@" + img_fn);

                // CTF
                md.setValue(EMDL_CTF_DEFOCUSU, data.particles[part_id].images[img_id].defU);
                md.setValue(EMDL_CTF_DEFOCUSV, data.particles[part_id].images[img_id].defV);
                md.setValue(EMDL_CTF_DEFOCUS_ANGLE, data.particles[part_id].images[img_id].defAngle);
                md.setValue(EMDL_MICROGRAPH_PRE_EXPOSURE, data.particles[part_id].images[img_id].dose);

                // Rotation
                Matrix2D<RFLOAT> a = data.getRotationMatrix(part_id, img_id) * A;
                Euler_matrix2angles(a, rot, tilt, psi);
                md.setValue(EMDL_ORIENT_ROT, rot);
                md.setValue(EMDL_ORIENT_TILT, tilt);
                md.setValue(EMDL_ORIENT_PSI, psi);

                // Translation
                RFLOAT shift2d_x, shift2d_y, shift2d_z;
                data.getTranslationInTiltSeries(
                        part_id, img_id,
                        shift3d_x, shift3d_y, shift3d_z,
                        shift2d_x, shift2d_y, shift2d_z
                );
                md.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, shift3d_x);
                md.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, shift3d_y);
                md.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, shift3d_z);
            }
        }

        std::cout << "Writing output to " << fn_out << std::endl;
        obsModel.save(md, fn_out, "particles");
	}
};


int main(int argc, char *argv[])
{
    pseudotomo_expansion exp;

	try
	{
        exp.read(argc, argv);
        exp.run();
	}
	catch (RelionError XE)
	{
		std::cerr << XE;
        exp.usage();
		return RELION_EXIT_FAILURE;
	}
	return RELION_EXIT_SUCCESS;
}

