#include "fit_blobs_3d.h"

#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/membrane/blob_fit_3d.h>
#include <src/jaz/membrane/membrane_segmentation.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/tomography/fiducials.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>


using namespace gravis;


void FitBlobs3DProgram::readParameters(int argc, char *argv[])
{
	double sphere_thickness_0;

	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		tomoSetFn = parser.getOption("--t", "Tomogram set filename", "tomograms.star");
		listFn = parser.getOption("--i", "File containing a list of tomogram-name/spheres-file pairs");

		sphere_thickness_0 = textToDouble(parser.getOption("--th", "Sphere thickness (same units as sphere centres)"));
		spheres_binning = textToDouble(parser.getOption("--sbin", "Binning factor of the sphere coordinates"));

		fiducials_radius_A = textToDouble(parser.getOption("--frad", "Fiducial marker radius [Å]", "100"));

		prior_sigma_A = textToDouble(parser.getOption("--sig", "Uncertainty std. dev. of initial position [Å]", "10"));
		max_binning = textToDouble(parser.getOption("--bin0", "Initial (maximal) binning factor", "8"));
		min_binning = textToDouble(parser.getOption("--bin1", "Final (minimal) binning factor", "2"));

		diag = parser.checkOption("--diag", "Write out diagnostic information");

		SH_bands = textToInteger(parser.getOption("--n", "Number of spherical harmonics bands", "2"));
		highpass_sigma_real_A = textToDouble(parser.getOption("--hp", "High-pass sigma [Å, real space]", "300"));
		max_iters = textToInteger(parser.getOption("--max_iters", "Maximum number of iterations", "1000"));
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));

		outPath = parser.getOption("--o", "Output filename pattern");

		Log::readParams(parser);

		if (parser.checkForErrors())
		{
			parser.writeUsage(std::cout);
			exit(1);
		}
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	sphere_thickness = sphere_thickness_0 * spheres_binning;

	outPath = ZIO::makeOutputDir(outPath);

	if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag");
	}
}

void FitBlobs3DProgram::run()
{
	TomogramSet initial_tomogram_set = TomogramSet(tomoSetFn);
	TomogramSet subtracted_tomogram_set = initial_tomogram_set;
	TomogramSet blobs_tomogram_set = initial_tomogram_set;
	
	if (!initial_tomogram_set.globalTable.labelExists(EMDL_TOMO_FIDUCIALS_STARFILE))
	{
		Log::warn("No fiducial markers present: you are advised to run relion_tomo_find_fiducials first.");
	}
	
	std::ifstream list(listFn);
	
	if (!list)
	{
		REPORT_ERROR_STR("Unable to read "+listFn);
	}
	
	std::map<std::string, std::string> tomoToSpheres;
	
	std::string line;

	BufferedImage<float> visualisation(0,0,0);
	
	while (std::getline(list, line))
	{
		std::stringstream sts;
		sts << line;
	
		std::string tomoName, spheresFn;
		sts >> tomoName;
		sts >> spheresFn;
		
		tomoToSpheres[tomoName] = spheresFn;
	}

	int tomo_index = 0;
		
	for (std::map<std::string, std::string>::iterator it = tomoToSpheres.begin();
	     it != tomoToSpheres.end(); it++)
	{
		const std::string tomoName = it->first;
		const std::string spheresFn = it->second;
		
		Log::beginSection("Tomogram " + tomoName);
		
		processTomogram(
			tomoName, spheresFn,
			initial_tomogram_set,
			subtracted_tomogram_set,
			blobs_tomogram_set,
			visualisation,
			tomo_index,
			tomoToSpheres.size());
		
		Log::endSection();

		tomo_index++;
	}

	subtracted_tomogram_set.write(outPath + "tomograms.star");
	blobs_tomogram_set.write(outPath + "blob_tomograms.star");

	if (diag)
	{
		visualisation.write(outPath + "diagnostic.mrc");
	}
}

void FitBlobs3DProgram::processTomogram(
		std::string tomoName,
		std::string spheresFn,
		TomogramSet& initial_tomogram_set,
		TomogramSet& subtracted_tomogram_set,
		TomogramSet& blobs_tomogram_set,
		BufferedImage<float>& visualisation,
		int tomo_batch_index,
		int tomo_batch_size)
{
	Log::print("Loading tilt series");
	
	spheres = readSpheresCMM(spheresFn, spheres_binning);
	
	const int tomo_index = initial_tomogram_set.getTomogramIndexSafely(tomoName);

	Tomogram tomogram0 = initial_tomogram_set.loadTomogram(tomo_index, true);
	const int frame_count = tomogram0.frameCount;

	const int w_full = tomogram0.stack.xdim;
	const int h_full = tomogram0.stack.ydim;

	const double pixel_size = tomogram0.optics.pixelSize;

	const double highpass_sigma_real = highpass_sigma_real_A / pixel_size;
	const double fiducials_radius = fiducials_radius_A / pixel_size;


	
	
	/*BufferedImage<float> original_stack = tomogram0.stack;

	tomogram0.stack = ImageFilter::highpassStackGaussPadded(
						original_stack,
						highpass_sigma_real,
						num_threads);*/
	
	        
	/*if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag/" + tomogram0.name);

		if (visualisation.xdim == 0)
		{
			visualisation.resize(w_full, h_full, 3 * tomo_batch_size);
		}

		visualisation.getSliceRef(3 * tomo_batch_index).copyFrom(
				original_stack.getConstSliceRef(frame_count / 2));
	}*/


	std::vector<d3Vector> fiducials(0);
	
	bool has_fiducials = 
	           tomogram0.fiducialsFilename.length() > 0 
	        && tomogram0.fiducialsFilename != "empty";

	//BufferedImage<float> fiducialsMask(w_full, h_full, frame_count);

	
	        
	Log::print("Filtering");
	
	
	const double segmentation_binning = 2;
	
	Tomogram tomogram_binned = tomogram0.FourierCrop(segmentation_binning, num_threads);
	
	if (has_fiducials)
	{
		Log::print("Erasing fiducial markers");
		
		if (fiducialsDir[fiducialsDir.length()-1] != '/')
		{
			fiducialsDir = fiducialsDir + "/";
		}

		fiducials = Fiducials::read(tomogram0.fiducialsFilename, tomogram0.optics.pixelSize);

		Fiducials::erase(
			fiducials, 
			fiducials_radius / segmentation_binning, 
			tomogram_binned, 
			num_threads);
	}
	
	tomogram_binned.stack.write("DEBUG_erased_binned.mrc");
	
	BufferedImage<float> preweighted_stack = RealSpaceBackprojection::preWeight(
	            tomogram_binned.stack, tomogram_binned.projectionMatrices, num_threads);
	
	Damage::applyWeight(
			preweighted_stack, 
			tomogram_binned.optics.pixelSize, 
			tomogram_binned.cumulativeDose, num_threads);
	
	
	/*BufferedImage<float> blobs_stack(w_full, h_full, frame_count);
	blobs_stack.fill(0.f);*/


	for (int blob_id = 0; blob_id < spheres.size(); blob_id++)
	{
		Log::beginSection("Blob #" + ZIO::itoa(blob_id + 1));
		
		const d4Vector sphere = spheres[blob_id];
		std::vector<double> initial = {sphere.x, sphere.y, sphere.z};

		std::vector<double> blob_coeffs = initial;

		const d3Vector sphere_position = sphere.xyz();
		const double outer_radius_full = sphere.w + sphere_thickness / 2;

		const int substack_size_binned = 2 * std::ceil(outer_radius_full / max_binning);
		const int substack_size = (int)(substack_size_binned * max_binning);

		presegmentBlob(
				sphere_position, 
				sphere.w, 
				sphere_thickness, 
		        segmentation_binning,
				preweighted_stack, 
				tomogram_binned.projectionMatrices);

		/*Tomogram tomogram_cropped = tomogram0.extractSubstack(
		            sphere_position, substack_size, substack_size);


		double current_binning = max_binning;

		while (current_binning > min_binning - 1e-6)
		{
			if (current_binning == max_binning)
			{
				Log::beginSection("Fitting at bin " + ZIO::itoa((int)current_binning));
			}
			else
			{
				Log::beginSection("Refining at bin " + ZIO::itoa((int)current_binning));
			}
			
			blob_coeffs = fitBlob(blob_id, blob_coeffs, current_binning, tomogram_cropped, fiducials);
			current_binning /= 2;
			
			Log::endSection();
		}

		SphericalHarmonics sh(SH_bands);
		Blob3D blob(blob_coeffs, spheres[blob_id].w + sphere_thickness / 2, &sh);

		#pragma omp parallel for num_threads(num_threads)
		for (int f = 0; f < frame_count; f++)
		{
			RawImage<float> stackSlice = original_stack.getSliceRef(f);
			RawImage<float> blobsSlice = blobs_stack.getSliceRef(f);

			blob.decompose(
				stackSlice, blobsSlice,
				tomogram0.projectionMatrices[f],
				fiducialsMask.getConstSliceRef(f),
				50.0);
		}*/
		
		Log::endSection();
	}

	/*const std::string tag = outPath + tomogram0.name;

	original_stack.write(tag + "_subtracted.mrc", tomogram0.optics.pixelSize);
	blobs_stack.write(tag + "_blobs.mrc", tomogram0.optics.pixelSize);

	subtracted_tomogram_set.setTiltSeriesFile(tomo_index, tag + "_subtracted.mrc");
	blobs_tomogram_set.setTiltSeriesFile(tomo_index, tag + "_blobs.mrc");

	if (diag)
	{
		visualisation.getSliceRef(3 * tomo_batch_index + 1).copyFrom(
				original_stack.getConstSliceRef(frame_count / 2));

		visualisation.getSliceRef(3 * tomo_batch_index + 2).copyFrom(
				blobs_stack.getConstSliceRef(frame_count / 2));
	}*/
}

BufferedImage<float> computeTiltSpaceMap(
        d3Vector sphere_position, 
        double mean_radius_full, 
        double radius_range, 
        double binning, 
        const RawImage<float>& preweighted_stack, 
        const std::vector<d4Matrix>& projections)
{
	const int fc = projections.size();
		
	const double perimeter_length_full = 2 * PI * mean_radius_full;
	
	const int w_map = perimeter_length_full / binning;
	const int h_map = radius_range / binning;
	
	const int w_stack = preweighted_stack.xdim;
	const int h_stack = preweighted_stack.ydim;
	
	const double min_radius_full = mean_radius_full - radius_range/2;
	
	BufferedImage<float> map(w_map, h_map, fc);
	
	for (int z = 0; z < fc; z++)
	for (int y = 0; y < h_map; y++)
	for (int x = 0; x < w_map; x++)
	{
		const int f0 = z;
		const d4Matrix& A = projections[f0];
		
		d3Vector dir_x(A(0,0), A(0,1), A(0,2));
		d3Vector dir_y(A(1,0), A(1,1), A(1,2));
		
		dir_x.normalize();
		dir_y.normalize();
		
		const double phi = 2 * PI * x / (double) w_map;
		const double r = min_radius_full + radius_range * y / (double) h_map;
		
		const d3Vector pos = sphere_position + r * (cos(phi) * dir_x + sin(phi) * dir_y);
		
		float sum = 0.f;
		float weight = 0.f;
		
		for (int f = 0; f < fc; f++)
		{
			const d4Vector pi = projections[f] * d4Vector(pos);

			if (pi.x >= 0.0 && pi.x < w_stack && pi.y >= 0.0 && pi.y < h_stack)
			{
				sum += Interpolation::linearXY_clip(preweighted_stack, pi.x, pi.y, f);
				weight += 1.0;
			}
		}
		
		if (weight > 0)
		{
			map(x,y,z) = sum / weight;
		}
		else
		{
			map(x,y,z) = 0;
		}     
	}
	
	return map;
}

std::vector<double> FitBlobs3DProgram::presegmentBlob(
        d3Vector sphere_position, 
        double mean_radius_full, 
        double radius_range, 
        double binning, 
        const RawImage<float>& preweighted_stack, 
        const std::vector<d4Matrix>& projections)
{
	BufferedImage<float> map = computeTiltSpaceMap(
			sphere_position, 
			mean_radius_full, 
			radius_range, 
			binning, 
			preweighted_stack, 
			projections);	            
	            
	map.write("DEBUG_tilt_space_map.mrc");
	
	const int fc = projections.size();
	
	const d4Matrix& A0 = projections[0];
	const d4Matrix& A1 = projections[fc - 1];
	
	const d3Vector tilt_axis_3D = 
			d3Vector(A0(2,0), A0(2,1), A0(2,2)).cross(
			d3Vector(A1(2,0), A1(2,1), A1(2,2))).normalize();
	
	std::vector<double> tilt_axis_azimuth(fc);
	
	for (int f = 0; f < fc; f++)
	{
		const d4Vector t(tilt_axis_3D.x, tilt_axis_3D.y, tilt_axis_3D.z, 0.0);
		const d4Vector tilt_axis_2D = projections[f] * t;
		            
		tilt_axis_azimuth[f] = std::atan2(tilt_axis_2D.y, tilt_axis_2D.x);
	}
	
	const int falloff = 100 / binning; // falloff
	const int width = 10 / binning; // width
	const double spacing = 40.0 / binning; // spacing
	const double ratio = 5.0; // ratio
	
	const double min_radius_full = mean_radius_full - radius_range/2;
	const double max_radius_full = mean_radius_full + radius_range/2;
		
	
	BufferedImage<float> kernel = MembraneSegmentation::constructMembraneKernel(
		map.xdim, map.ydim, map.zdim, falloff, width, spacing, ratio);   
	 
	kernel.write("DEBUG_tilt_space_kernel.mrc");
	
	BufferedImage<fComplex> map_FS, kernel_FS, correlation_FS;
	
	FFT::FourierTransform(map, map_FS);
	FFT::FourierTransform(kernel, kernel_FS);
	
	correlation_FS.resize(map_FS.xdim, map_FS.ydim, map_FS.zdim);
	
	for (int z = 0; z < map_FS.zdim; z++)
	for (int y = 0; y < map_FS.ydim; y++)
	for (int x = 0; x < map_FS.xdim; x++)
	{
		correlation_FS(x,y,z) = map_FS(x,y,z) * kernel_FS(x,y,z).conj();
	}
	
	BufferedImage<float> correlation;
	        
	FFT::inverseFourierTransform(correlation_FS, correlation);
	
	const double st = map.zdim / 4.0;
	const double s2t = 2 * st * st;
	const double s2f = 2 * falloff * falloff;
	        
	for (int f = 0; f < map.zdim; f++)
	for (int y = 0; y < map.ydim; y++)
	for (int x = 0; x < map.xdim; x++)
	{
		const double r = (min_radius_full + radius_range * y / (double) map.ydim) / max_radius_full;
		
		const double ay = map.ydim - y - 1;		        
		const double q0 = 1.0 - exp( -y*y  / s2f);
		const double q1 = 1.0 - exp(-ay*ay / s2f);
		
		const double phi = 2 * PI * x / (double) map.xdim;
		const double mz = (f - map.zdim / 2) * sin(phi - tilt_axis_azimuth[f]);	
		const double qt = exp(-mz*mz / s2t);
		
		correlation(x,y,f) *= r * q0 * q1 * qt;
	}
	
	const float corr_var = Normalization::computeVariance(correlation, 0.f);
	correlation /= sqrt(corr_var);
	
	correlation.write("DEBUG_tilt_space_correlation.mrc");
	
	
	
	/*BufferedImage<float> surfaceCost(map.xdim, map.ydim, map.zdim);
	BufferedImage<float> indicator(map.xdim, map.ydim, map.zdim);
	
	const float sig_cost = 4;
	        
	for (int f = 0; f < map.zdim; f++)
	for (int y = 0; y < map.ydim; y++)
	for (int x = 0; x < map.xdim; x++)
	{
		float c = correlation(x,y,f);
		if (c < 0) c = 0;
		
		surfaceCost(x,y,f) = exp(-c*c/(sig_cost*sig_cost));
		//indicator(x,y,f) = 1.f - y / (float) (map.ydim - 1);
		//indicator(x,y,f) = 0.5f;
		indicator(x,y,f) = y < map.ydim/2? 1 : 0;
	}
	
	surfaceCost.write("DEBUG_tilt_space_surfaceCost.mrc");
	
	const int iterations = 501;
	
	{
		BufferedImage<f3Vector> xi(map.xdim, map.ydim, map.zdim);
		xi.fill(f3Vector(0,0,0));
		
		BufferedImage<float> u(map.xdim, map.ydim, map.zdim);
		BufferedImage<float> uBar(map.xdim, map.ydim, map.zdim);
		
		u    = indicator;
		uBar = indicator;
				
		const float sigma = 0.5;
		const float tau = 0.5;
		const float nu = 1.8;
		
		
		for (int it = 0; it < iterations; it++)
		{
			for (int f = 0; f < map.zdim; f++)
			for (int y = 0; y < map.ydim; y++)
			for (int x = 0; x < map.xdim; x++)
			{
				const float u0 = uBar(x,y,f);
				const float ux = uBar((x+1)%map.xdim, y, f);
				
				float uy;
				
				if (y == map.ydim - 1)
				{
					uy = uBar(x, y-1, f);
				}
				else
				{
					uy = uBar(x, y+1, f);
				}
				
				float uz;
				
				if (f == fc-1)
				{
					uz = uBar(x, y, f-1);
				}
				else
				{
					uz = uBar(x, y, f+1);
				}
				
				
				const f3Vector grad = f3Vector(ux - u0, uy - u0, uz - u0);
				
				xi(x,y,f) += sigma * grad;
				
				if (xi(x,y,f).norm2() > 1)
				{
					xi(x,y,f) /= xi(x,y,f).length();
				}
			}
			
			for (int f = 0; f < map.zdim; f++)
			for (int y = 0; y < map.ydim; y++)
			for (int x = 0; x < map.xdim; x++)
			{
				const f3Vector g0 = xi(x,y,f);
				const f3Vector gx = xi((x+map.xdim-1)%map.xdim, y, f);
				
				f3Vector gy;
				
				if (y == 0)
				{
					gy = -xi(x, 0, f);
				}
				else
				{
					gy = xi(x, y-1, f);
				}
				
				f3Vector gz;
				
				if (f == 0)
				{
					gz = -xi(x, y, 0);
				}
				else
				{
					gz = xi(x, y, f-1);
				}
				
				const int d = 3;
				const float regional = (y < d? -1 : (y > map.ydim-d-1? 1 : 0));
				
				const float div = (g0.x - gx.x) + (g0.y - gy.y) + (g0.z - gz.z);
				
				const float lastU = u(x,y,f);
				
				const float sc = surfaceCost(x,y,f);
				
				u(x,y,f) += tau * (nu * sc * div - regional);
				
				
				if (u(x,y,f) > 1) 
				{
					u(x,y,f) = 1;
				}
				else if (u(x,y,f) < 0)
				{
					u(x,y,f) = 0;
				}
				
				uBar(x,y,f) = 2 * lastU - u(x,y,f);
			}
			
			if (it % 50 == 0)
			{
				u.write("DEBUG_tilt_space_F_"+ZIO::itoa(it)+"_u.mrc");
				uBar.write("DEBUG_tilt_space_F_"+ZIO::itoa(it)+"_uBar.mrc");
				
				BufferedImage<float> xi_y(map.xdim, map.ydim, map.zdim);
				
				for (int i = 0; i < xi.getSize(); i++)
				{
					xi_y[i] = xi[i].y;
				}
				
				xi_y.write("DEBUG_tilt_space_F_"+ZIO::itoa(it)+"_xi_y.mrc");
			}
		}
	}*/
	
	std::exit(0);
}

std::vector<double> FitBlobs3DProgram::fitBlob(
		int blob_id,
		const std::vector<double>& initial,
		double binning_factor,
		Tomogram& tomogram_cropped,
		const std::vector<d3Vector>& fiducials)
{
	const double pixelSize_full = tomogram_cropped.optics.pixelSize;
	const double pixelSize_binned = tomogram_cropped.optics.pixelSize * binning_factor;

	// 3D coordinates live in bin-1 pixel coordinates;
	const d4Vector sphere = spheres[blob_id];

	// 2D radii, thicknesses and distances are binned
	const double sphere_radius_full = sphere.w;
	const double outer_radius_full = sphere.w + sphere_thickness / 2;
	const double sphere_radius_binned = sphere_radius_full / binning_factor;
	const double outer_radius_binned = outer_radius_full / binning_factor;

	// The prior is evaluated in 2D: consider binning
	const double prior_sigma = prior_sigma_A / pixelSize_full;


	const double fiducials_radius_binned = fiducials_radius_A / pixelSize_binned;
	const double sphere_thickness_binned = sphere_thickness / binning_factor;

	const double initial_step = 2.0;

	const std::string tomoTag = tomogram_cropped.name;
	std::string blobTag = "blob_" + ZIO::itoa(blob_id);
	std::string binTag = "bin_" + ZIO::itoa((int)binning_factor);

	std::string outTag = outPath + "diag/" + tomoTag + "/" + blobTag + "/" + binTag;

	Tomogram tomogram_binned = tomogram_cropped.FourierCrop(binning_factor, num_threads);

	if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag/" + tomoTag + "/" + blobTag);
		tomogram_binned.stack.write(outTag + "_data.mrc");
	}

	std::vector<double> last_optimum = fromBin1(initial, binning_factor);
	int initial_SH_bands = initial.size() < 7? 0 : sqrt(initial.size() - 3) - 1;

	for (int current_SH_bands = initial_SH_bands; current_SH_bands <= SH_bands; current_SH_bands++)
	{
		const int SH_coeffs = (current_SH_bands + 1) * (current_SH_bands + 1);
		
		Log::print(ZIO::itoa(current_SH_bands) + " SH bands (" + ZIO::itoa(SH_coeffs + 2) + " parameters)");

		std::string tag = outTag + "_SH_" + ZIO::itoa(current_SH_bands);

		BlobFit3D bf(
			tomogram_binned, spheres[blob_id].xyz(), current_SH_bands,
			sphere_radius_binned,
			sphere_thickness_binned,
			spheres,
			fiducials, fiducials_radius_binned,
			prior_sigma, num_threads);

		std::vector<double> current_optimum(SH_coeffs + 3, 0.0);

		for (int i = 0; i < last_optimum.size(); i++)
		{
			current_optimum[i] = last_optimum[i];
		}

		std::vector<double> new_optimum;

		new_optimum = NelderMead::optimize(
				current_optimum, bf, initial_step, 0.05, max_iters, 1.0, 2.0, 0.5, 0.5, false);

		if (diag)
		{
			SphericalHarmonics sh(current_SH_bands);
			Blob3D blob2(new_optimum, outer_radius_binned, &sh);
			
			const double E0 = bf.f(current_optimum, &sh);
			const double E1 = bf.f(new_optimum, &sh);
			
			Log::print("  E: "+ZIO::itoa(E0)+" -> "+ZIO::itoa(E1));

			drawTestStack(blob2, tomogram_binned, bf.weight).write(tag+"_fit.mrc");
			drawFit(blob2, tomogram_binned, bf.weight).write(tag+"_residual.mrc");
		}

		last_optimum = new_optimum;
	}


	std::vector<double> upscaled_optimum = toBin1(last_optimum, binning_factor);

	return upscaled_optimum;
}

std::vector<d4Vector> FitBlobs3DProgram::readSpheresCMM(
		const std::string& filename,
		double binning)
{
	std::string nextPointKey = "<marker ";

	std::string formatStr =
		"<marker id=\"%d\" x=\"%lf\" y=\"%lf\" z=\"%lf\" r=\"%*f\" g=\"%*f\" b=\"%*f\" radius=\"%lf\"/>";

	std::ifstream ifs(filename);

	char buffer[1024];

	std::vector<d4Vector> spheres(0);

	while (ifs.getline(buffer, 1024))
	{
		std::string line(buffer);

		if (ZIO::beginsWith(line, nextPointKey))
		{
			int id;
			double x, y, z, rad;

			if (std::sscanf(line.c_str(), formatStr.c_str(), &id, &x, &y, &z, &rad) != 5)
			{
				REPORT_ERROR_STR("Bad syntax in " << filename << ": " << line);
			}

			spheres.push_back(binning * d4Vector(x,y,z,rad));
		}
	}

	return spheres;
}


BufferedImage<float> FitBlobs3DProgram::drawFit(
		Blob3D& blob,
		const Tomogram& tomogram,
		BufferedImage<float>& realWeight)
{
	const int fc = tomogram.frameCount;
	const int w = tomogram.stack.xdim;
	const int h = tomogram.stack.ydim;

	BufferedImage<float> out(w,h,fc);

	for (int f = 0; f < fc; f++)
	{
		std::vector<double> radAvg = blob.radialAverage(
			tomogram.stack.getConstSliceRef(f),
			tomogram.projectionMatrices[f],
			realWeight.getConstSliceRef(f));

		BufferedImage<float> diff = blob.drawError(
			tomogram.stack.getConstSliceRef(f),
			tomogram.projectionMatrices[f],
			realWeight.getConstSliceRef(f),
			radAvg);

		out.copySliceFrom(f, diff);
	}

	return out;
}


BufferedImage<float> FitBlobs3DProgram::drawTestStack(
		Blob3D& blob,
		const Tomogram& tomogram,
		BufferedImage<float>& realWeight)
{
	const int fc = tomogram.frameCount;
	const int w = tomogram.stack.xdim;
	const int h = tomogram.stack.ydim;

	BufferedImage<float> out(w,h,fc);

	for (int f = 0; f < fc; f++)
	{
		std::vector<double> radAvg = blob.radialAverage(
			tomogram.stack.getConstSliceRef(f),
			tomogram.projectionMatrices[f],
			realWeight.getConstSliceRef(f));

		BufferedImage<float> radAvgProj = blob.radialAverageProjection(
			tomogram.stack.getConstSliceRef(f),
			tomogram.projectionMatrices[f],
			radAvg);

		out.copySliceFrom(f, radAvgProj);
	}

	return out;
}

std::vector<double> FitBlobs3DProgram::toBin1(
		const std::vector<double> &parameters,
		double binning_factor)
{
	std::vector<double> upscaled_optimum = parameters;

	for (int i = 3; i < upscaled_optimum.size(); i++)
	{
		upscaled_optimum[i] *= binning_factor;
	}

	return upscaled_optimum;
}

std::vector<double> FitBlobs3DProgram::fromBin1(
		const std::vector<double> &parameters,
		double binning_factor)
{
	const int xc = parameters.size();
	std::vector<double> downscaled_optimum(xc);

	for (int i = 0; i < downscaled_optimum.size(); i++)
	{
		downscaled_optimum[i] = parameters[i];

		if (i > 2)
		{
			downscaled_optimum[i] /= binning_factor;
		}
	}

	return downscaled_optimum;
}
