#include <src/args.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/projection/real_backprojection.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/image/filter.h>
#include <src/jaz/image/detection.h>
#include <src/jaz/image/similarity.h>
#include <src/jaz/mesh/mesh.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/tomography/fiducials.h>
#include <src/jaz/tomography/programs/find_fiducials.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>

using namespace gravis;


int main(int argc, char *argv[])
{
	try
	{
		FindFiducialsProgram program(argc, argv);
		program.run();
	}
	catch (RelionError XE)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}
