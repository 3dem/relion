#include "fit_blobs_3d.h"
#include <src/jaz/membrane/tilt_space_blob_fit.h>
#include <src/jaz/membrane/blob_fit_3d.h>
#include <src/jaz/membrane/membrane_segmentation.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/tomography/fiducials.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/manifold/manifold_set.h>
#include <src/jaz/tomography/manifold/manifold_loader.h>
#include <src/jaz/tomography/manifold/sphere.h>
#include <src/jaz/tomography/manifold/spheroid.h>
#include <src/jaz/tomography/manifold/CMM_loader.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>


using namespace gravis;


void FitBlobs3DProgram::readParameters(int argc, char *argv[])
{
	IOParser parser;

	parser.setCommandLine(argc, argv);

	optimisationSet.read(
		parser,
		true,            // optimisation set
		false,  false,   // particles
		true,   true,    // tomograms
		false,  false,   // trajectories
		true,   true,    // manifolds
		false,  false);  // reference

	int gen_section = parser.addSection("General options");

	relative_radius_range = textToDouble(parser.getOption("--rr", "Sphere radius range as a multiple of radius", "1"));
	membrane_separation = textToDouble(parser.getOption("--ms", "Membrane separation [Å]", "40"));
	fiducials_radius_A = textToDouble(parser.getOption("--frad", "Fiducial marker radius [Å]", "100"));
	fit_binning = textToDouble(parser.getOption("--bin", "Binning at which to perform the fit", "2"));

	diag = parser.checkOption("--diag", "Write out diagnostic information");

	SH_bands = textToInteger(parser.getOption("--n", "Number of spherical harmonics bands", "5"));
	lowpass_sigma_real_A = textToDouble(parser.getOption("--lp", "Low-pass sigma [Å, real space]", "-1"));
	max_iters = textToInteger(parser.getOption("--max_iters", "Maximum number of iterations", "1000"));
	num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));

	outPath = parser.getOption("--o", "Output filename pattern");

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	outPath = ZIO::prepareTomoOutputDirectory(outPath, argc, argv);

	ZIO::makeDir(outPath + "Meshes");

	if (diag)
	{
		ZIO::makeDir(outPath + "diag");
	}
}

void FitBlobs3DProgram::run()
{
	TomogramSet tomogram_set(optimisationSet.tomograms);
	ManifoldSet input_manifold_set(optimisationSet.manifolds);

	ManifoldSet output_manifold_set;
	

	if (!tomogram_set.globalTable.containsLabel(EMDL_TOMO_FIDUCIALS_STARFILE))
	{
		Log::warn("No fiducial markers present: you are advised to run relion_tomo_find_fiducials first.");
	}

	const int tc = tomogram_set.size();
		
	for (int t = 0; t < tc; t++)
	{
		std::string tomogram_name = tomogram_set.getTomogramName(t);

		std::map<int, const Manifold*> manifolds_map = input_manifold_set.getManifoldsInTomogram(tomogram_name);

		if (manifolds_map.empty()) continue;

		Log::beginSection("Tomogram " + tomogram_name);

		processTomogram(
			t, tomogram_name, manifolds_map, tomogram_set, output_manifold_set);
		
		Log::endSection();
	}
	
	output_manifold_set.write(outPath + "manifolds.star");

	optimisationSet.manifolds = outPath + "manifolds.star";
	optimisationSet.write(outPath + "optimisation_set.star");
}

void FitBlobs3DProgram::processTomogram(
		int tomo_index,
		const std::string& tomogram_name,
		const std::map<int, const Manifold*>& input_manifolds_map,
		const TomogramSet& tomogram_set,
		ManifoldSet& output_manifold_set)
{
	Log::print("Loading tilt series");

	Tomogram tomogram0 = tomogram_set.loadTomogram(tomo_index, true);
	const double pixel_size = tomogram0.optics.pixelSize;
	const double fiducials_radius = fiducials_radius_A / pixel_size;

	std::vector<d3Vector> fiducials(0);
	bool has_fiducials = tomogram0.hasFiducials();

	Log::print("Filtering");
		
	Tomogram tomogram_binned = tomogram0.FourierCrop(fit_binning, num_threads);


	if (lowpass_sigma_real_A > 0.0)
	{
		tomogram_binned.stack = ImageFilter::lowpassStack(
			tomogram_binned.stack,
			lowpass_sigma_real_A / (pixel_size * fit_binning),
			30.0, false);
	}

	if (has_fiducials)
	{
		Log::print("Erasing fiducial markers");

		fiducials = Fiducials::read(tomogram0.fiducialsFilename, tomogram0.optics.pixelSize);

		Fiducials::erase(
			fiducials, 
			fiducials_radius / fit_binning,
			tomogram_binned, 
			num_threads);
	}
	
	BufferedImage<float> preweighted_stack = RealSpaceBackprojection::preWeight(
			tomogram_binned.stack, tomogram_binned.projectionMatrices, num_threads);
	
	Damage::applyWeight(
			preweighted_stack,
			tomogram_binned.optics.pixelSize,
			tomogram_binned.cumulativeDose, num_threads);
	
	std::vector<std::vector<double>> all_blob_coeffs;
	Mesh blob_meshes;
	TomogramManifoldSet output_tomogram_manifold_set;

	for (std::map<int, const Manifold*>::const_iterator it = input_manifolds_map.begin();
		 it != input_manifolds_map.end(); it++)
	{
		const int sphere_id = it->first;
		const Manifold* manifold = it->second;

		if (manifold->type != Sphere::getTypeName()) continue;

		Log::beginSection("Sphere #" + ZIO::itoa(sphere_id + 1));
		
		std::vector<double> sphere_params = manifold->getParameters();

		const d3Vector sphere_position(sphere_params[0], sphere_params[1], sphere_params[2]);
		const double sphere_radius = sphere_params[3];
		
		const std::string blob_tag = outPath + "diag/" + tomogram_name + "_blob_" + ZIO::itoa(sphere_id);
		
		std::vector<double> blob_coeffs = segmentBlob(
					sphere_position,
					sphere_radius,
					relative_radius_range * sphere_radius,
					membrane_separation,
					fit_binning,
					preweighted_stack,
					pixel_size,
					tomogram_binned.projectionMatrices,
					diag? blob_tag : "");
		
		all_blob_coeffs.push_back(blob_coeffs);
		
		Mesh blob_mesh = createMesh(blob_coeffs, pixel_size, 50, 60);

		output_tomogram_manifold_set.add(new Spheroid(blob_coeffs, sphere_id));

		MeshBuilder::insert(blob_mesh, blob_meshes);
		Log::endSection();
	}

	blob_meshes.writePly(outPath + "Meshes/" + tomogram_name + ".ply");
	
	output_manifold_set.add(tomogram_name, output_tomogram_manifold_set);
}

Mesh FitBlobs3DProgram::createMesh(
		const std::vector<double>& blob_coeffs,
		double pixel_size,
		double spacing,
		double max_tilt_deg)
{
	const double max_tilt = DEG2RAD(max_tilt_deg);
	const double rad = blob_coeffs[3];
	const int azimuth_samples = (int) std::round(2 * PI * rad / spacing);
	const int tilt_samples = (int) std::round(2 * max_tilt * rad / spacing);
	const d3Vector centre(blob_coeffs[0],blob_coeffs[1],blob_coeffs[2]);
	
	const int vertex_count = azimuth_samples * tilt_samples;
	
	Mesh out;
	out.vertices.resize(vertex_count);
	
	const int SH_params = blob_coeffs.size() - 3;
	
	const int SH_bands = (int) sqrt(SH_params - 1);
	SphericalHarmonics SH(SH_bands);
	
	std::vector<double> Y(SH_params);
	
	for (int a = 0; a < azimuth_samples; a++)
	for (int t = 0; t < tilt_samples; t++)
	{
		const double phi   = 2 * PI * a / (double) azimuth_samples;
		const double theta = -max_tilt + 2 * max_tilt * t / (double) (tilt_samples - 1);

		SH.computeY(SH_bands, sin(theta), phi, &Y[0]);

		double dist = 0.0;

		for (int b = 0; b < SH_params; b++)
		{
			dist += blob_coeffs[b+3] * Y[b];
		}

		out.vertices[t * azimuth_samples + a] =
				pixel_size * (centre + dist * d3Vector(cos(theta) * cos(phi), cos(theta) * sin(phi), sin(theta)));
	}
	
	const int triangle_count = 2 * (tilt_samples - 1) * azimuth_samples;

	out.triangles.resize(triangle_count);
	
	for (int a = 0; a < azimuth_samples; a++)
	for (int t = 0; t < tilt_samples-1; t++)
	{
		Triangle tri0;

		tri0.a =  t      * azimuth_samples +  a;
		tri0.b = (t + 1) * azimuth_samples +  a;
		tri0.c = (t + 1) * azimuth_samples + (a + 1) % azimuth_samples;
		
		Triangle tri1;
		
		tri1.a =  t      * azimuth_samples +  a;
		tri1.b = (t + 1) * azimuth_samples + (a + 1) % azimuth_samples;
		tri1.c =  t      * azimuth_samples + (a + 1) % azimuth_samples;
		
		out.triangles[2 * (t * azimuth_samples + a)    ] = tri0;
		out.triangles[2 * (t * azimuth_samples + a) + 1] = tri1;
	}
	
	return out;
}


std::vector<double> FitBlobs3DProgram::segmentBlob(
		d3Vector sphere_position,
		double mean_radius_full,
		double radius_range,
		double membrane_separation,
		double binning,
		const RawImage<float>& preweighted_stack,
		double pixel_size,
		const std::vector<d4Matrix>& projections,
		const std::string& debug_prefix)
{
	BufferedImage<float> map = TiltSpaceBlobFit::computeTiltSpaceMap(
			sphere_position, 
			mean_radius_full, 
			radius_range, 
			binning, 
			preweighted_stack, 
			projections);
	
	BufferedImage<d3Vector> directions_XZ = TiltSpaceBlobFit::computeDirectionsXZ(
				mean_radius_full, binning, projections);

	if (debug_prefix != "")
	{
		map.write(debug_prefix+"_tilt_space_map.mrc");
	}
	
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
	
	const int y_prior = 100 / binning;
	const int falloff = 100 / binning;
	const int width = 30 / binning;
	const double spacing = membrane_separation / (pixel_size * binning);
	const double ratio = 5.0;
	const double depth = 0;
	
	const double min_radius_full = mean_radius_full - radius_range/2;
	const double max_radius_full = mean_radius_full + radius_range/2;
		
	const double max_tilt = DEG2RAD(30);
	const double tilt_steps = 15;


	/*BufferedImage<float> kernel = MembraneSegmentation::constructMembraneKernel(
		map.xdim, map.ydim, map.zdim, falloff, width, spacing, ratio, depth);   
	
	if (debug_prefix != "")
	{
		kernel.write(debug_prefix+"_tilt_space_kernel.mrc");
	}
	
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
	        
	FFT::inverseFourierTransform(correlation_FS, correlation);*/




	BufferedImage<float> correlation =
		MembraneSegmentation::correlateWithMembraneMultiAngle(
			map, falloff, width, spacing, ratio, depth,
			max_tilt, tilt_steps);


	
	const double st = map.zdim / 4.0;
	const double s2t = 2 * st * st;
	const double s2f = 2 * y_prior * y_prior;
	        
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
	
	if (debug_prefix != "")
	{
		correlation.write(debug_prefix+"_tilt_space_correlation.mrc");
	}
	
	const double lambda = 0.00001;
	
	
	TiltSpaceBlobFit blob_pre_fit(0, lambda, correlation, directions_XZ);
	double h0 = blob_pre_fit.estimateInitialHeight();

	const double Y00_norm = 1.0  / (2 * sqrt(PI));

	std::vector<double> last_params = {h0 / Y00_norm};
	
	for (int current_SH_bands = 1; current_SH_bands <= SH_bands; current_SH_bands++)
	{
		TiltSpaceBlobFit blob_fit(current_SH_bands, lambda, correlation, directions_XZ);
		
		std::vector<double> params(blob_fit.getParameterCount(), 0.0);
		
		for (int i = 0; i < last_params.size() && i < params.size(); i++)
		{
			params[i] = last_params[i];
		}


		// try explicitly squishing Z to avoid local optima

		/*if (current_SH_bands == 2)
		{
			BufferedImage<float> plot0 = blob_fit.drawSolution(params, map);
			plot0.write(debug_prefix+"_squish_0.mrc");

			const int squish_samples = 20;
			const double min_squish = 0.5;
			const double max_squish = 1.0;

			double best_f = std::numeric_limits<double>::max();
			double best_squish = 1.0;

			for (int q = 0; q < squish_samples; q++)
			{
				const double squish = min_squish + q * (max_squish - min_squish)
						/ (squish_samples - 1);

				std::vector<double> squished_params = params;

				squished_params[6] = 2.0 * params[0] * (squish - 1.0) / (3.0 * sqrt(5.0));
				squished_params[0] = params[0] + (sqrt(5.0) / 2.0) * squished_params[6];

				const double f = blob_fit.f(squished_params, 0);

				std::cout << squish << ": " << f << std::endl;

				if (f < best_f)
				{
					best_f = f;
					best_squish = squish;
				}
			}

			params[6] = 2.0 * params[0] * (best_squish - 1.0) / (3.0 * sqrt(5.0));
			params[0] = params[0] + (sqrt(5.0) / 2.0) * params[6];

			BufferedImage<float> plot1 = blob_fit.drawSolution(params, map);
			plot1.write(debug_prefix+"_squish_1.mrc");
		}*/
		
		std::vector<double> final_params = LBFGS::optimize(params, blob_fit, 0, 1000, 1e-6);

		BufferedImage<float> plot = blob_fit.drawSolution(final_params, map);

		plot.write(debug_prefix+"_tilt_space_plot_SH_"+ZIO::itoa(current_SH_bands)+".mrc");
		
		/*BufferedImage<float> plot = blob_fit.drawSolution(final_params, map);
		plot.write("DEBUG_tilt_space_plot_SH_"+ZIO::itoa(SH_bands)+".mrc");*/
				
		last_params = final_params;
	}
	
	if (debug_prefix != "")
	{
		TiltSpaceBlobFit blob_fit(SH_bands, lambda, correlation, directions_XZ);
		BufferedImage<float> plot = blob_fit.drawSolution(last_params, map);
		
		plot.write(debug_prefix+"_tilt_space_plot_SH_"+ZIO::itoa(SH_bands)+".mrc");
	}
	
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
	
	std::vector<double> out(last_params.size() + 3);
	
	for (int i = 3; i < out.size(); i++)
	{
		out[i] = binning * last_params[i-3];
	}

	out[3] += min_radius_full / Y00_norm;
	
	SphericalHarmonics SH_3(1);
	std::vector<double> Y_3(4);
	SH_3.computeY(1, 0, 0, &Y_3[0]);
	
	for (int i = 0; i < 3; i++)
	{
		out[i] = sphere_position[i] - out[i+4] * Y_3[0];
		out[i+4] = 0.0;
	}
	
	if (debug_prefix != "")
	{
		BufferedImage<float> plot = TiltSpaceBlobFit::visualiseBlob(
			out, 
			mean_radius_full, 
			radius_range, 
			binning, 
			preweighted_stack, 
			projections);
		
		plot.write(debug_prefix+"_tilt_space_plot_SH_"+ZIO::itoa(SH_bands)+"_blob_space.mrc");
	}
	
	return out;
}

