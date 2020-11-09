#include <src/args.h>
#include <src/jaz/tomography/programs/template_picker.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	TemplatePickerProgram tpp;

	tpp.readParameters(argc, argv);
	tpp.run();

	return 0;
}
