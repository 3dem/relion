#include <src/args.h>
#include <src/jaz/tomography/programs/sample_manifold.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	SampleManifoldProgram sp;

	sp.readParameters(argc, argv);
	sp.run();

	return 0;
}
