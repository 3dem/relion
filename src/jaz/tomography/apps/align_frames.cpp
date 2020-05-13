#include <src/args.h>
#include <src/jaz/tomography/programs/frame_alignment_program.h>


int main(int argc, char *argv[])
{
	FrameAlignmentProgram fap(argc, argv);
	
	fap.run();
	
	return 0;
}
