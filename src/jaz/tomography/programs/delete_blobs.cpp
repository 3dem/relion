#include "delete_blobs.h"

#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/membrane/blob_fit.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>


using namespace gravis;


void DeleteBlobsProgram::readParameters(int argc, char *argv[])
{
	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		tomoSetFn = parser.getOption("--t", "Tomogram set filename", "tomograms.star");
		tomoName = parser.getOption("--tn", "Tomogram name");
		spheresFn = parser.getOption("--sn", "Spheres filename");
		spheres_binning = textToDouble(parser.getOption("--sbin", "Binning factor of the sphere coordinates"));

		diag = parser.checkOption("--diag", "Write out diagnostic information");

		SH_bands = textToInteger(parser.getOption("--n", "Number of spherical harmonics bands", "2"));
		max_iters = textToInteger(parser.getOption("--max_iters", "Maximum number of iterations", "1000"));
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));

		outer_margin = textToDouble(parser.getOption("--m_out", "Outside margin (relative to radius)", "1.5"));
		inner_margin = textToDouble(parser.getOption("--m_in",  "Inside margin (relative to radius)", "0.5"));

		outPath = parser.getOption("--o", "Output filename pattern");

		Log::readParams(parser);

		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	spheres = readSpheresCMM(spheresFn, spheres_binning);

	outPath = ZIO::makeOutputDir(outPath);
}

void DeleteBlobsProgram::run()
{
	d4Vector b = spheres[0];

	TomogramSet tomogramSet = TomogramSet(tomoSetFn);
	const int tomo_index = tomogramSet.getTomogramIndexSafely(tomoName);

	Tomogram tomogram0 = tomogramSet.loadTomogram(tomo_index, true);



	const int test_frame = 19;

	const int vesicle_id = 0;
	const d4Vector sphere = spheres[vesicle_id];
	const d3Vector sphere_position(sphere.x, sphere.y, sphere.z);
	const double sphere_radius = sphere.w;

	const double outer_radius_0 = outer_margin * sphere_radius;
	const double inner_radius_0 = inner_margin * sphere_radius;

	const bool use_masks = true;
	const double prior_sigma = 2.0;
	const double initial_step = 1.0;
	const double binning_factor = 8.0;
	const int substack_size = 1024;

	const double outer_radius = outer_radius_0 / binning_factor;
	const double inner_radius = inner_radius_0 / binning_factor;

	const int frame_count = tomogram0.frameCount;




	Tomogram tomogram1 = tomogram0.extractSubstack(sphere_position, substack_size, substack_size);
	Tomogram tomogram = tomogram1.FourierCrop(binning_factor, num_threads);

	{
		Blob blob(sphere_position, outer_radius);

		tomogram.stack.getSliceRef(test_frame).write(outPath+"original_f"+ZIO::itoa(test_frame)+".mrc");

		std::vector<double> radAvg = blob.radialAverage(tomogram, test_frame);
		BufferedImage<float> radAvgProj = blob.radialAverageProjection(tomogram, test_frame, radAvg);

		radAvgProj.write(outPath+"blob0_f"+ZIO::itoa(test_frame)+".mrc");
	}



	BlobFit bf0(
		tomogram, spheres, vesicle_id, inner_radius,
		outer_radius, 0, use_masks, prior_sigma, num_threads);

	std::vector<double> initial0 = {sphere_position.x, sphere_position.y, sphere_position.z};

	std::vector<double> opt0;

	opt0 = NelderMead::optimize(
				initial0, bf0, initial_step, 0.05, 100, 1.0, 2.0, 0.5, 0.5, diag);

	if (!diag)
	{
		std::cout << "           -> ";

		for (int i = 0; i < 3; i++)
		{
			std::cout << opt0[i] << ", ";
		}

		std::cout << std::endl;
	}
	else
	{
		Blob blob0(opt0, outer_radius, 0);
		std::vector<double> radAvg0 = blob0.radialAverage(tomogram, test_frame);
		BufferedImage<float> radAvgProj0 = blob0.radialAverageProjection(tomogram, test_frame, radAvg0);

		radAvgProj0.write("dev/ves-"+ZIO::itoa(vesicle_id)+"_rad_avg_2_L0_f"+ZIO::itoa(test_frame)+".mrc");
	}

	for (int current_SH_bands = 1; current_SH_bands <= SH_bands; current_SH_bands++)
	{
		const int SH_coeffs = (current_SH_bands + 1) * (current_SH_bands + 1);

		BlobFit bf(
			tomogram, spheres, vesicle_id, inner_radius,
			outer_radius, current_SH_bands, use_masks, prior_sigma, num_threads);

		std::vector<double> initial(SH_coeffs + 3, 0.0);

		for (int i = 0; i < opt0.size(); i++)
		{
			initial[i] = opt0[i];
		}

		std::vector<double> opt;

		opt = NelderMead::optimize(
				initial, bf, initial_step, 0.05, 200, 1.0, 2.0, 0.5, 0.5, diag);

		if (!diag)
		{
			std::cout << "           -> ";

			for (int i = 0; i < 3; i++)
			{
				std::cout << opt[i] << ", ";
			}

			std::cout << "   ";

			for (int j = 1; j <= current_SH_bands; j++)
			{
				for (int i = -j; i <= j; i++)
				{
					std::cout << opt[3 + j*j + j + i] << ", ";
				}

				std::cout << "   ";
			}

			std::cout << std::endl;
		}
		else
		{
			SphericalHarmonics sh(current_SH_bands);

			Blob blob2(opt, outer_radius, &sh);
			std::vector<double> radAvg2 = blob2.radialAverage(tomogram, test_frame);
			BufferedImage<float> radAvgProj2 = blob2.radialAverageProjection(tomogram, test_frame, radAvg2);

			std::stringstream sts;
			sts << current_SH_bands;
			radAvgProj2.write(
				"dev/ves-"+ZIO::itoa(vesicle_id)
				+"_rad_avg_2_L"+ZIO::itoa(current_SH_bands)+"_f"+ZIO::itoa(test_frame)+".mrc");
		}

		opt0 = opt;
	}


	std::vector<double> opt_upscaled = opt0;

	for (int i = 3; i < opt_upscaled.size(); i++)
	{
		opt_upscaled[i] *= binning_factor;
	}



	{
		tomogram1.stack.getSliceRef(test_frame).write(outPath+"original_f"+ZIO::itoa(test_frame)+"_bin1.mrc");


		SphericalHarmonics sh(SH_bands);
		Blob blob(opt_upscaled, outer_radius_0, &sh);

		std::vector<double> radAvg = blob.radialAverage(tomogram1, test_frame);
		BufferedImage<float> radAvgProj = blob.radialAverageProjection(tomogram1, test_frame, radAvg);

		radAvgProj.write(outPath+"blob0_f"+ZIO::itoa(test_frame)+"_bin1.mrc");
	}


	const bool refine_full = false;

	std::vector<double> opt_fin;

	if (refine_full)
	{
		BlobFit bf(
			tomogram1, spheres, vesicle_id, inner_radius_0, outer_radius_0,
			SH_bands, use_masks, prior_sigma * binning_factor, num_threads);

		opt_fin = NelderMead::optimize(
			opt_upscaled, bf, initial_step, 0.05, 200, 1.0, 2.0, 0.5, 0.5, diag);


		{
			SphericalHarmonics sh(SH_bands);
			Blob blob(opt_fin, outer_radius_0, &sh);

			std::vector<double> radAvg = blob.radialAverage(tomogram1, test_frame);
			BufferedImage<float> radAvgProj = blob.radialAverageProjection(tomogram1, test_frame, radAvg);

			radAvgProj.write(outPath+"blob0_f"+ZIO::itoa(test_frame)+"_bin1_final.mrc");
		}
	}
	else
	{
		opt_fin = opt_upscaled;
	}



	SphericalHarmonics sh(SH_bands);

	Blob blobOut(opt_fin, outer_radius_0, &sh);

	for (int f = 0; f < frame_count; f++)
	{
		blobOut.subtract(tomogram0, f, 50.0);
	}

	if (diag)
	{
		for (int f = 0; f < frame_count; f++)
		{
			blobOut.subtract(tomogram0, f, 50.0);
		}

		if (diag)
		{
			tomogram0.stack.write(outPath+"blob_"+ZIO::itoa(vesicle_id)+"_stack.mrc", tomogram0.optics.pixelSize);
		}
	}
}

std::vector<d4Vector> DeleteBlobsProgram::readSpheresCMM(
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
