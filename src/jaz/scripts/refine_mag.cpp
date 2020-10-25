#include <src/args.h>
#include <src/jaz/tomography/programs/mag_fit.h>


int main(int argc, char *argv[])
{
	MagFitProgram mfp(argc, argv);

	mfp.run();

	return 0;
}
