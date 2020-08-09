#include <src/args.h>
#include <src/jaz/tomography/programs/subtomo.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	SubtomoProgram sp;

	sp.readParameters(argc, argv);
	sp.run();

	return 0;
}
