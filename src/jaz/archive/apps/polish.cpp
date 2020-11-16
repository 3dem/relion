#include <src/args.h>
#include <src/jaz/tomography/programs/polish.h>


int main(int argc, char *argv[])
{
	PolishProgram pp(argc, argv);
	
	pp.run();
	
	return 0;
}
