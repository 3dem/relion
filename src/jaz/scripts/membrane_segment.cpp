
#include <src/jaz/tomography/tomo_stack.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/image/detection.h>
#include <src/jaz/image/similarity.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/tomography/dynamo/catalogue.h>
#include <src/jaz/membrane/membrane_fit.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/image/structure_tensor.h>
#include <src/jaz/segmentation/diffusion.h>
#include <src/jaz/d3x3/dsyevh3.h>
#include <src/jaz/util/new_vtk_helper.h>
#include <src/jaz/segmentation/primal_dual_TV.h>
#include <src/jaz/segmentation/diffusion_tensors.h>
#include <src/jaz/segmentation/skeletonization.h>
#include <src/jaz/membrane/membrane_segmentation.h>
#include <src/jaz/membrane/cells.h>
#include <src/jaz/tomography/tomolist.h>

#include <omp.h>

using namespace gravis;


int main(int argc, char *argv[])
{
	std::string tomoTag = "214";
			
	std::string pointsFn = "coords_"+tomoTag+"_bin1.txt";//"coords_145_bin8.txt";
	std::vector<std::vector<double>> coords = ZIO::readFixedDoublesTable(pointsFn, 6);
	
	const double binning = 16.0;
	const int pc = coords.size();
	
	const bool debug = true;
	const double rho = 1.0; // 5
	const double sigma0 = 0.0; // 0
	const int num_threads = 6;
	
	const double sigma_diff = 1;
	const double lambda_edge = 0.5;
	const double thresh_edge = 0.1;
	const double sigma_plate = 0.0002;
	
	const int iters_across = 30;
	const int iters_along = 500;
	
	
	{
		
		BufferedImage<float> membrane;
		membrane.read("debug/membrane.mrc");
		
		BufferedImage<float> softmaxDist = MembraneSegmentation::softmaxMembraneDist(membrane, 5);
		softmaxDist.write("debug/softmaxDist.mrc");
		
		BufferedImage<float> centers;
		centers.read("debug/maxCent_FS.mrc");
		
		Tapering::taper(centers, 30, 6, false);
		
		centers /= softmaxDist;
		
		centers.write("debug/maxCent_by_softmaxDist.mrc");
				
		return 0;
	}
	
	
	
	
	
	
	
	
	BufferedImage<float> tomo;
	tomo.read("bin16/R1_"+tomoTag+".mrc");
	
	Tapering::taper(tomo, 20, num_threads, true);
	
	TomoList tomoList("local_tomolist.txt");
	std::string projFn = tomoList.getProjectionsFilename(145);
	
	
	int w0, h0, d0;
	std::vector<d4Matrix> projTomo = ProjectionIO::read(projFn, w0, h0, d0);
	
	
	std::vector<d3Vector> south(pc), north(pc);
	std::vector<double> thickness(pc);
	
	for (int p = 0; p < pc ; p++)
	{
		south[p] = d3Vector(coords[p][0], coords[p][1], coords[p][2]) / binning;
		north[p] = d3Vector(coords[p][3], coords[p][4], coords[p][5]) / binning;
		thickness[p] = (north[p] - south[p]).length();
	}
	
	
	const int w = tomo.xdim;
	const int h = tomo.ydim;
	const int d = tomo.zdim;
	
	BufferedImage<float> data(w,h,d), regional(w,h,d), weight(w,h,d);
	BufferedImage<float> regional_memb(w,h,d), weight_memb(w,h,d);
	data.fill(0.f);
	weight.fill(0.f);
	regional.fill(0.f);
	weight_memb.fill(0.f);
	regional_memb.fill(0.f);
	
	for (int p = 0; p < pc ; p++)
	{
		d3Vector nd = north[p];
		
		const int r = 3;
		
		for (int dx = -r; dx <= r; dx++)
		for (int dy = -r; dy <= r; dy++)
		for (int dz = -r; dz <= r; dz++)
		{
			regional_memb(
				(int)(nd.x + 0.5) + dx,
				(int)(nd.y + 0.5) + dy,
				(int)(nd.z + 0.5) + dz) = 1.f;
			
			weight_memb(
				(int)(nd.x + 0.5) + dx,
				(int)(nd.y + 0.5) + dy,
				(int)(nd.z + 0.5) + dz) = 1.f;
		}
	}
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for	(int x = 0; x < w; x++)
	{
		const d3Vector r(x,y,z);
				
		for (int p = 0; p < pc ; p++)
		{
			const d3Vector n = north[p];
			const d3Vector s = south[p];
			const d3Vector t = 2.0 * n - s;
			
			const double dp = 0.9 * thickness[p];
			
			const double ds = (s - r).length();
			const double dt = (t - r).length();
			
			const double ra = (r - n).dot(n - s)/dp;
			
			if (dt < dp || ds < dp)
			{
				weight(x,y,z) += 1.f;
				data(x,y,z) -= ra;
				
				const double thr = 1;
				
				if (ra > thr)
				{
					regional(x,y,z) = 1.f;
				}
				else if (ra < -thr)
				{
					regional(x,y,z) = -1.f;
				}	
			}
			
			double bounary_conf = 0.2;
			double boundary = Tapering::getTaperWeight3D(x,y,z, w,h,d, 5);
			
			regional(x,y,z) += bounary_conf * (1.0 - boundary);
			weight(x,y,z) += bounary_conf * (1.0 - boundary);
			
			double conf_reg = 0.001;
			double reg_thresh = -4.0;
			
			regional(x,y,z) += conf_reg * (tomo(x,y,z) - reg_thresh);
			weight(x,y,z) += conf_reg;
		}
		
		
	}
	
	if (debug)
	{
		weight.write("debug/weight.mrc");
		regional.write("debug/regional.mrc");
		weight_memb.write("debug/weight_memb.mrc");
		regional_memb.write("debug/regional_memb.mrc");
		data.write("debug/data.mrc");
	}
	
	BufferedImage<Tensor3x3<float>> J = StructureTensor::compute3D(tomo, rho, sigma0, 20.0);
	J = StructureTensor::forwardAverage(J);
	
	/*std::vector<Image<float>> evals0 = StructureTensor::split(
				StructureTensor::computeEigenvalues3D(J));
	
	evals0[0].write("debug/Jlin_eval0.mrc");
	(evals0[0] - evals0[1]).write("debug/Jlin_eval0by1.mrc");
			
	Image<Tensor3x3<float>> J_nl = StructureTensor::computeNonLinear(
		tomo, rho, sigma0, 20.0,
		20, 20, true, sigma_diff, 0.05f, num_threads);
	
	
	std::vector<Image<float>> evals1 = StructureTensor::split(
				StructureTensor::computeEigenvalues3D(J_nl));
			
	evals1[0].write("debug/Jnonlin_eval0.mrc");
	(evals1[0] - evals1[1]).write("debug/Jnonlin_eval0by1.mrc");*/
				
	std::cout << "finding membranes..." << std::endl;
	
	BufferedImage<float> membrane = MembraneSegmentation::determineMembraniness(
				tomo, J, sigma_diff, lambda_edge, thresh_edge, 
				num_threads, iters_across, iters_along);
	
	membrane.write("debug/membrane.mrc");
	
	
	BufferedImage<float> softmaxDist = MembraneSegmentation::softmaxMembraneDist(membrane, 10);
	softmaxDist.write("debug/softmaxDist.mrc");
			
			
	BufferedImage<float> ptSymmsMembrane = Cells::findPointSymmetries(membrane, 20, tomo.xdim/24);
	ptSymmsMembrane.write("debug/membrane_ptSymms.mrc");
	
	BufferedImage<Tensor3x3<float>> surfKernel = Skeletonization::getKernel(
				J, Skeletonization::Surface, 0.1f);	
	BufferedImage<float> surf = Skeletonization::apply(membrane, surfKernel, num_threads, 1);	
	surf.write("debug/membrane_thinned.mrc");
		
	BufferedImage<float> ptSymmsMembraneThin = Cells::findPointSymmetries(surf, 20, tomo.xdim/24);
	ptSymmsMembraneThin.write("debug/thin_membrane_ptSymms.mrc");
	
	
	//std::vector<d3Vector> surfPts = Skeletonization::discretize(surf, 1.f);

	//std::cout << surfPts.size() << " surface points found." << std::endl;
	
	BufferedImage<float> maxCent(w,h,d), maxRad(w,h,d);
	
	Cells::findCenters(surf, 10, 100.f, 1.f, maxCent, maxRad, num_threads);
	
	
	maxCent.write("debug/maxCent_FS.mrc");
	maxRad.write("debug/maxRad_FS.mrc");
	
	//Image<float> ptSymms = Cells::findPointSymmetries(tomo);
	
	
	
	//NewVtkHelper::writeR3(eig0, "debug/eigenvector0.vtk");
	
	/*Image<float> diffused = Diffusion::diffuse(tomo, D, 0.05f, 500, num_threads);
	diffused.write("debug/diffused.mrc");
	
	Image<float> diffusedNoise = Diffusion::diffuse(noise, D, 0.05f, 500, num_threads);
	diffusedNoise.write("debug/diffused_noise.mrc");*/
	
	/*Image<float> segTV = Segmentation::anisotropicTV(regional, D, 100, 0.1f, 0.1f);
	segTV.write("debug/seg-TV_100.mrc");*/
	
}

