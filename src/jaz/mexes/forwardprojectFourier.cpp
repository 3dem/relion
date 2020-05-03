#include <math.h>
#include <matrix.h>
#include <mex.h>

#include "include/mex_parser.h"
#include "include/dynamo_helper.h"

#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/tomo/extraction.h>
#include <src/jaz/tomo/projection/projection.h>
#include <src/jaz/tomo/reconstruction.h>

using namespace gravis;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MexParser mp(nlhs, plhs, nrhs, prhs);

	RawImage<float> map = mp.getImage<float>(      MexParser::Right, 0, "Reference map  [single, S x S x S]");
	std::vector<d3Vector> angles = mp.getR3Vector( MexParser::Right, 1, "Euler angles   [3 x #frames]");
	const int n_threads = (int) mp.getDouble(      MexParser::Right, 2, "#threads", 1);
	
	const int s = map.xdim;
	const int sh = s/2 + 1;
	const int fc = angles.size();
	
	RawImage<float> out = mp.createImage<float>(MexParser::Left,  0, "Forward-projected stack [single, S x S x #frames]", s, s, fc);
	
	if (!mp.finalize()) return;
	
	
	// prepare matrices:
	
	std::vector<d4Matrix> projections(fc);
	
	for (int f = 0; f < fc; f++)
	{
		d3Vector ang = angles[f];
		projections[f] = anglesToMatrix(ang[0], ang[1], ang[2], 0, 0, 0, 0, 0);
	}
	
	// FT map:
	
	Image<fComplex> mapFS;
	Image<float> mapCopy(map);
	FFT::FourierTransform(mapCopy, mapFS, FFT::Both);
	
	Centering::shiftInSitu(mapFS);
	
	// perform projection:
	
	Image<fComplex> dataFS(sh,s,fc), psfFS(sh,s,fc);
	Projection::forwardProject(mapFS, projections, dataFS, psfFS, 1.0, n_threads);
	
	// perform gridding correction:
	
	Reconstruction::correctStack(dataFS, psfFS, out, true, n_threads);
}
