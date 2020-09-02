#include "tomo_ctf.h"

#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/centering.h>
#include <src/args.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/tomography/tomo_ctf_find.h>

#include <omp.h>

using namespace gravis;


void TomoCtfProgram::run()
{
	int w0, h0, d0;

    std::vector<d4Matrix> proj = ProjectionIO::read(projFn, w0, h0, d0);

    if (thickness < 0)
    {
        thickness = d0;
        std::cout << "Using thickness from '" << projFn << "': " << thickness << std::endl;
    }
	
	BufferedImage<float> stack;
	stack.read(stackFn);
	
	const double Q0 = 0.07;
	const int s = 512;
	const bool diag = true;
	
	// add user-supplied dose numbers
	
	TomoCtfFind tcf(stack, proj, diag, pixSize, voltage, Cs, Q0, s, 50.0, 10.0);
	tcf.initialize(n_threads);
			
	//tcf.findAstigmatism();
	
	tcf.findInitialCtf(10000.0, 60000.0, 100);
	
	std::vector<double> deltaF = tcf.findDefoci(10000.0, 60000.0, 100);
	
	/*for (int f = 0; f < stack.zdim; f++)
	{
		std::cout << f << ": " << deltaF[f] << "\n";
	}*/
	
}
