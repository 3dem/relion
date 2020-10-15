
#include <src/jaz/tomography/tomo_stack.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/membrane/blob_3d.h>
#include <src/jaz/membrane/blob_fit_3d.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/image/detection.h>
#include <src/jaz/image/similarity.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/single_particle/volume_converter.h>
#include <src/spherical-harmonics/SphericalHarmonics.h>

#include <omp.h>

using namespace gravis;
using namespace au::edu::anu::qm::ro;

template<class T>
void test(const RawImage<T>& img)
{
	std::cout << img.ydim << "\n";
}
		  
int main(int argc, char *argv[])
{
	BufferedImage<float> imgVec0(64,64,64);
	RawImage<float> img0(64,64,64,imgVec0.data);
	std::cout << img0.ydim << "\n";
	test(img0);
	return 0;
}
