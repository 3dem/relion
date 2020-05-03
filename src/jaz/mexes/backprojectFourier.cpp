#include <math.h>
#include <matrix.h>
#include <mex.h>

#include "include/mex_parser.h"
#include "include/dynamo_helper.h"

#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/tomo/extraction.h>
#include <src/jaz/tomo/projection/projection.h>
#include <src/jaz/tomo/reconstruction.h>

using namespace gravis;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MexParser mp(nlhs, plhs, nrhs, prhs);

	RawImage<float> stack = mp.getImage<float>(    MexParser::Right, 0, "Tilt series    [single, W x H x #frames]");
	std::vector<d2Vector> pos2D = mp.getR2Vector(  MexParser::Right, 1, "2D coordinates [2 x #frames]");
	std::vector<d3Vector> angles = mp.getR3Vector( MexParser::Right, 2, "Euler angles   [3 x #frames]");	
	const int s = (int) mp.getDouble(              MexParser::Right, 3, "Output box size 'S'");
	const double wc = mp.getDouble(                MexParser::Right, 4, "Wiener fract.", 0.1);
	const int n_threads = (int) mp.getDouble(      MexParser::Right, 5, "#threads", 1);
	
	const int w = stack.xdim;
	const int h = stack.ydim;
	const int fc = stack.zdim;
	const int sh = s/2 + 1;
	
	RawImage<float> out = mp.createImage<float>(MexParser::Left,  0, "Backprojected volume [single, S x S x S]", s, s, s);
	
	if (!mp.finalize()) return;
	
	
	// prepare matrices:
	
	std::vector<d4Matrix> projections(fc);
	
	for (int f = 0; f < fc; f++)
	{
		d3Vector ang = angles[f];
		d2Vector pos = pos2D[f];
		projections[f] = anglesToMatrix(ang[0], ang[1], ang[2], pos.x, pos.y, w, h, 0);
	}
	
	
	// extract and FFT squares:

	Image<fComplex> smallStackFS_shifted(sh,s,fc);
	std::vector<d4Matrix> projNew(fc);

	Extraction::extractShiftedSquares(
		stack, s, projections, pos2D, 
		smallStackFS_shifted, projNew, n_threads);
	
	
	// perform backprojection:	
	
	Image<fComplex> 
			dataImgFS(sh,s,s);
	
	Image<float>
			ctfImgFS(sh,s,s),
			psfImgFS(sh,s,s);
	
	dataImgFS.fill(fComplex(0.f, 0.f));
	ctfImgFS.fill(0.f);
	psfImgFS.fill(0.f);
	
	Projection::backprojectFS_bwd(
		smallStackFS_shifted, projNew, dataImgFS, psfImgFS, ctfImgFS, 1.0, n_threads);
	
	
	// perform gridding and CTF correction:	
	
	Image<float> outCopy;
	
	Reconstruction::correct3D(
		dataImgFS, psfImgFS, ctfImgFS, outCopy,
		wc * ctfImgFS(0,0,0), true, n_threads);
	
	outCopy.copyTo(out);	
}
