
#include <src/jaz/tomography/dynamo/catalogue.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/projection/point_insertion.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/optics/damage.h>
#include <src/time.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <iostream>

using namespace gravis;


int main(int argc, char *argv[])
{
	IOParser parser;

	std::string particles_filename, output_filename;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		particles_filename = parser.getOption("--i", "Catalogue .tbl or .star file");
		output_filename = parser.getOption("--o", "Output filename pattern");

		Log::readParams(parser);

		parser.checkForErrors();

	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}




	ParticleSet data_set(particles_filename, "");
	ParticleSet data_set_copy = data_set;

	const int pc = data_set.getTotalParticleNumber();

	for (int p = 0; p < pc; p++)
	{
		double rot_particle = data_set.partTable.getDouble(EMDL_ORIENT_ROT, p);
		double tilt_particle = data_set.partTable.getDouble(EMDL_ORIENT_TILT, p);
		double psi_particle = data_set.partTable.getDouble(EMDL_ORIENT_PSI, p);

		data_set_copy.partTable.setValue(EMDL_TOMO_SUBTOMOGRAM_ROT, rot_particle, p);
		data_set_copy.partTable.setValue(EMDL_TOMO_SUBTOMOGRAM_TILT, tilt_particle, p);
		data_set_copy.partTable.setValue(EMDL_TOMO_SUBTOMOGRAM_PSI, psi_particle, p);

		data_set_copy.partTable.setValue(EMDL_ORIENT_ROT, 0.0, p);
		data_set_copy.partTable.setValue(EMDL_ORIENT_TILT, 0.0, p);
		data_set_copy.partTable.setValue(EMDL_ORIENT_PSI, 0.0, p);
	}

	data_set_copy.write(output_filename);
}
