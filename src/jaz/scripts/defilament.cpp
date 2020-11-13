#include <src/jaz/math/Euler_angles_dynamo.h>
#include <src/macros.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/image/gradient.h>
#include <src/jaz/image/filter.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/local_extrema.h>
#include <src/jaz/image/structure_tensor.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/math/spline.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/tomography/filament/filament.h>
#include <src/jaz/tomography/filament/filament_fit.h>
#include <src/jaz/tomography/filament/filament_model.h>
#include <src/jaz/tomography/filament/circular_Fourier_filament_model.h>
#include <src/jaz/tomography/projection_IO.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/lbfgs.h>
#include <omp.h>
#include <stdio.h>
#include <vector>

using namespace gravis;

/*
	TODO:
  
	- consider all frames again
	- write non-circular Fourier model
	- write explicit-position fit (outside of FilamentFit?)
	- write TV segmentation with varying surface cost but constant (axis-aligned) anisotropy
	  (n.x expensive, n.y and n.z cheap)
	

	outline:
	
		variant A:
		
			- iterate similarity fit until convergence (better average, more constant-error regions)
			- perform TV segmentation
			- fit deformation to model
			- offset model parameters
			- iterate similarity fit again
			
			
		variant B:
		
			- find segmentation cost using a synthetic template
			- perform TV segmentation
			- use to intialize similarity fit
  
*/

void read_splines(
		std::string filename, 
		std::vector<std::vector<d3Vector>>& points, 
		std::vector<std::vector<double>>& radii,
		double spline_binning)
{
	std::string nextSplineKey = "<marker_set ";
	std::string nextPointKey = "<marker ";
	
	std::string formatStr = 
		"%*s %*s %d %*s %lf %*s %lf %*s %lf %*s %*f %*s %*f %*s %*f %*s %lf";
			
	
	std::ifstream ifs(filename);
	
	char buffer[1024];
	
	points.clear();
	radii.clear();
	
	int splineIndex = -1;
	
	while (ifs.getline(buffer, 1024))
	{
		std::string line(buffer);
		
		if (ZIO::beginsWith(line, nextSplineKey))
		{
			points.push_back(std::vector<d3Vector>(0));
			radii.push_back(std::vector<double>(0));
			
			splineIndex++;
		}
		else if (ZIO::beginsWith(line, nextPointKey))
		{
			for (int i = 0; i < line.length(); i++)
			{
				if (line[i] == '"') line[i] = ' ';
			}
			
			int id;
			double x, y, z, r;
			
			std::sscanf(line.c_str(), formatStr.c_str(), &id, &x, &y, &z, &r);
			
			points[splineIndex].push_back(spline_binning * d3Vector(x,y,z));
			radii[splineIndex].push_back(r);
		}
	}
}

BufferedImage<float> compute_mask(
		const BufferedImage<float>& binnedStack,
		double edge_bin,
		double sigma,
		double thresh,
		double grow_thresh,
		int iters,
		double current_binning,
		int num_threads
		)
{
	const double edge_bin_rel = edge_bin / current_binning;
	
	
	const int w  = binnedStack.xdim;
	const int h  = binnedStack.ydim;
	const int fc = binnedStack.zdim;
	
	BufferedImage<float> out(w,h,fc);
	
	BufferedImage<float> binnedStackSmooth = ImageFilter::GaussStack(binnedStack, edge_bin_rel, true);
	
	BufferedImage<float> edgeStack = Resampling::downsampleFiltStack_2D_full(
				binnedStackSmooth, edge_bin_rel, num_threads);	
	
	BufferedImage<float> edgeStackHP = ImageFilter::highpassStackGaussPadded(edgeStack, sigma);	
	
	edgeStackHP = Normalization::byNormalDistByFrame(edgeStackHP);		

	edgeStackHP = Resampling::upsampleCubic_Stack_full(
				edgeStackHP, edge_bin_rel, w, h);
			
	for (int f = 0; f < fc; f++)		
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out(x,y,f) = edgeStackHP(x,y,f) > thresh? 1.0 : 0.0;
	}
	
	for (int i = 0; i < iters; i++)
	{
		out = ImageFilter::GaussStack(out, edge_bin_rel, true);
		
		for (int f = 0; f < fc; f++)		
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			out(x,y,f) = out(x,y,f) > grow_thresh? 1.0 : 0.0;
		}
	}
	
	return out;
}

void precondition(BufferedImage<float>& binnedStack)
{
	const float mean = Normalization::computeMean(binnedStack);
	const float var = Normalization::computeVariance(binnedStack, mean);
	const float stdDev = sqrt(var);
	
	for (size_t i = 0; i < binnedStack.getSize(); i++)
	{
		binnedStack[i] = (binnedStack[i] - mean) / stdDev;
	}
}

int main(int argc, char *argv[])
{
	std::string splineFn = "spline_fit_coord/clicker_23.cmm";
	const double spline_binning = 6.0;
	
	std::string stackFn = "TS_23/TS_23_aligned_bin8.mrc";
	const double initial_binning = 8.0;
	const double max_binning = 8.0;
	
	std::string projFn = "proj/23.proj";
	
	std::string cacheDir = "cache/";
	const bool cached_filament_coordinates = true;
	
	const double radius_bin_1 = 500.0;
	const int num_iters = 100;
	
	const bool use_Fourier_cropping = false;
	const bool hi_pass = false;
	const bool low_pass = false;
	
	
	std::string outDir = "1st_max8_DS/";
	
	
	
	int dummy = std::system(("mkdir -p "+outDir).c_str());
	
	const int num_freqs = 20;
	const int num_threads = 6;
	
	
	
	std::cout << "reading projection matrices from: " << projFn << "..." << std::endl;
	
	int w3D, h3D, d3D;			
	std::vector<d4Matrix> proj0 = ProjectionIO::read(projFn, w3D, h3D, d3D);
	
	for (int f = 0; f < proj0.size(); f++)
	{
		proj0[f] /= initial_binning;
	}
	
	
	
	std::cout << "reading filament splines from: " << splineFn << "..." << std::endl;
	
	std::vector<std::vector<d3Vector>> splines;
	std::vector<std::vector<double>> radii;	
	
	read_splines(splineFn, splines, radii, spline_binning);
	
	const int filament_count = splines.size();
	
	std::vector<Filament> filaments(filament_count);
	
	for (int fil = 0; fil < filament_count; fil++)
	{
		filaments[fil] = Filament(splines[fil], radius_bin_1);		
	}
	
		
	
	std::cout << "reading image stack from: " << stackFn << "..." << std::endl;
	
	BufferedImage<float> origStack_unfilt;
	origStack_unfilt.read(stackFn);	
	
	BufferedImage<float> origStack = origStack_unfilt;
		
	if (hi_pass)
	{
		std::cout << "   filtering..." << std::endl;
		
		origStack = ImageFilter::highpassStackGaussPadded(origStack, radius_bin_1 / initial_binning);	
	}
	
	std::cout << "   preconditioning..." << std::endl;
		
	precondition(origStack);
	
	origStack.write(outDir+"initial_stack.mrc");	
	
	
	/*{
		const double rho = 7.0;
		Image<Tensor2x2<float>> J = StructureTensor::compute2D(origStack, rho, 0.0, 5.0);
		Image<f2Vector> evals = StructureTensor::computeEigenvalues2D(J);		
		std::vector<Image<float>> evalsSplit = StructureTensor::split(evals);
		
		
		Image<float> diff = evalsSplit[0] - evalsSplit[1];
		
		evalsSplit[0].write("eval_0.mrc");
		evalsSplit[1].write("eval_1.mrc");
		diff.write("evals_diff.mrc");
		
		return 0;
	}*/
	
	
	std::cout << "   computing outlier mask..." << std::endl;
	
	BufferedImage<float> marker_mask_full_size;
	
	
	{
		const double edge_bin = 16.0;
		
		const double sigma = 20.0;
		const double thresh = -3.0;
		const double grow_thresh = 0.8;
		const int iters = 3;
			
		
		marker_mask_full_size = compute_mask(
				origStack, edge_bin, sigma, thresh, grow_thresh, iters, initial_binning, num_threads);
		
		marker_mask_full_size.write(outDir+"outlier_mask.mrc");
	}	
	
	std::cout << "\ninitiating multiscale fit..." << std::endl;	
	
	
	const int scaleCount = (int) round(log(max_binning/initial_binning)/log(2.0)) + 1;
	
	std::vector<double> 
			abs_binning_levels(scaleCount), 
			rel_binning_levels(scaleCount);
	
	{
		double temp_binning = max_binning;
		
		for (int scale = 0; scale < scaleCount; scale++)
		{
			abs_binning_levels[scale] = temp_binning;
			rel_binning_levels[scale] = temp_binning / initial_binning;
			
			temp_binning /= 2.0;
		}
	}
	
	
	std::vector<std::vector<double>> current_params(filament_count);
	
	for (int fil = 0; fil < filament_count; fil++)
	{
		current_params[fil] = std::vector<double>(8 * num_freqs, 0.0);
	}
	
	for (int scale = 0; scale < scaleCount; scale++)
	{
		std::cout << "   scale #" << (scale + 1) << ", " << abs_binning_levels[scale] << std::endl;
		
		double current_binning = abs_binning_levels[scale];
		double current_rel_binning = rel_binning_levels[scale];	
		
		BufferedImage<float> binnedStack, binnedStack_unfilt, marker_mask;
				
		if (use_Fourier_cropping)
		{
			
			binnedStack = Resampling::FourierCrop_fullStack(
						origStack, current_rel_binning, num_threads, true);
			
			binnedStack_unfilt = Resampling::FourierCrop_fullStack(
							origStack_unfilt, current_rel_binning, num_threads, true);
			
			marker_mask = Resampling::FourierCrop_fullStack(
						marker_mask_full_size, current_rel_binning, num_threads, true);
			
		}
		else
		{
			binnedStack = Resampling::downsampleFiltStack_2D_full(
							origStack, current_rel_binning, num_threads);
			
			binnedStack_unfilt = Resampling::downsampleFiltStack_2D_full(
							origStack_unfilt, current_rel_binning, num_threads);
			
			marker_mask = Resampling::downsampleFiltStack_2D_full(
						marker_mask_full_size, current_rel_binning, num_threads);
		}
		
		if (low_pass) 
		{
			binnedStack = ImageFilter::GaussStack(binnedStack, 0.5, false);
		}
		
		std::string scaleTag = "bin"+ZIO::itoa(abs_binning_levels[scale]);
				
		binnedStack.write(outDir + scaleTag + "_filtered_stack.mrc");		
		binnedStack_unfilt.write(outDir + scaleTag + "_unfiltered_stack.mrc");
		

		const int w = binnedStack.xdim;
		const int h = binnedStack.ydim;
		
		marker_mask.write(outDir + scaleTag + "_mask.mrc");
			
		
		std::vector<d4Matrix> proj = proj0;
		
		for (int f = 0; f < proj.size(); f++)
		{
			proj[f] /= current_rel_binning;
		}
		
		std::vector<CircularFourierFilamentModel> models;
		std::vector<FilamentMapping> mappings(filament_count);
		std::vector<FilamentFit> initialFits;
		models.reserve(filament_count);
		
		std::cout << "      rasterising filament coordinates..." << std::endl;		
		
		BufferedImage<float> allErased = binnedStack;
		
		for (int fil = 0; fil < filament_count; fil++)
		{
			filaments[fil] = Filament(splines[fil], radius_bin_1);
			
			if (!cached_filament_coordinates)
			{
				mappings[fil] = filaments[fil].rasteriseCoordinates(
							w, h, proj, current_binning, num_threads);
				
				mappings[fil].write(
							cacheDir+scaleTag+"_filament_mapping_"+ZIO::itoa(fil));
			}
			else
			{
				mappings[fil] = FilamentMapping::read(
							cacheDir+scaleTag+"_filament_mapping_"+ZIO::itoa(fil));
			}
			
			models.push_back(CircularFourierFilamentModel(filaments[fil].arcLen, proj));
			
			initialFits.push_back(FilamentFit(filaments[fil], mappings[fil], proj, binnedStack, 
								  marker_mask, &models[fil], current_binning, num_threads));
			
			allErased = initialFits[fil].visualise(current_params[fil], true, &allErased);
		}
		
		std::cout << "      fitting filaments..." << std::endl;
		
		
		BufferedImage<float> allFitsOut = binnedStack_unfilt;
		BufferedImage<float> allErasedOut = binnedStack_unfilt;
		
		
		for (int fil = 0; fil < filament_count; fil++)
		{
			std::cout << "         " << (fil+1) << " / " << filament_count << std::endl;		
			
			BufferedImage<float> othersErased = initialFits[fil].visualise(
						current_params[fil], true, &allErased, -1);
			
			FilamentFit subtractedFit(filaments[fil], mappings[fil], proj, othersErased, 
							marker_mask, &models[fil], current_binning, num_threads);
			
			
			std::string tag = scaleTag + "_fil_" + ZIO::itoa(fil);
	
			
			
			/*current_params[fil] = LBFGS::optimize(
						current_params[fil], subtractedFit, true, num_iters, 1e-7);
			
			ZIO::writeDat(current_params[fil], outDir+"params_"+tag+".dat");*/	
			
			
			current_params[fil] = ZIO::readDat<double>(outDir+"params_"+tag+".dat");
				
				BufferedImage<float> costByOffset = subtractedFit.computeCostByOffset(
						current_params[fil], -10, 10, 42, 1.0);
				
				costByOffset.write(outDir+tag+"_costByOffset.mrc");
				
			
			
			{
				FilamentFit visFit(filaments[fil], mappings[fil], proj, binnedStack_unfilt, 
								   marker_mask, &models[fil], current_binning, num_threads);
				
				BufferedImage<float> predSub = visFit.visualise(current_params[fil], true, &allErasedOut);
				predSub.write(outDir+tag+"_sub.mrc");
				
				allFitsOut = visFit.visualise(current_params[fil], false, &allFitsOut);
				allFitsOut.write(outDir+tag+"_fit.mrc");
				
				allErasedOut = predSub;
			}
		}
	}
	
	return 0;
}
