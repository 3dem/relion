#include <src/args.h>
#include <src/jaz/tomo_programs/aberration_fit.h>


int main(int argc, char *argv[])
{
	AberrationFit af(argc, argv);
	
	af.run();
	
	return 0;
}
