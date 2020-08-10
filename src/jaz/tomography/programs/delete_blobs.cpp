#include "delete_blobs.h"

#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/membrane/blob_fit.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/tomography/fiducials.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>


using namespace gravis;


void DeleteBlobsProgram::readParameters(int argc, char *argv[])
{
	double sphere_thickness_0;

	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		tomoSetFn = parser.getOption("--t", "Tomogram set filename", "tomograms.star");
		tomoName = parser.getOption("--tn", "Tomogram name");

		spheresFn = parser.getOption("--sn", "Spheres filename");
		sphere_thickness_0 = textToDouble(parser.getOption("--th", "Sphere thickness (same units as sphere centres)"));
		spheres_binning = textToDouble(parser.getOption("--sbin", "Binning factor of the sphere coordinates"));

		fiducialsDir = parser.getOption("--fid", "Fiducial markers directory", "");
		fiducials_radius_A = textToDouble(parser.getOption("--frad", "Fiducial marker radius [Å]", "100"));

		prior_sigma_A = textToDouble(parser.getOption("--sig", "Uncertainty std. dev. of initial position [Å]", "10"));

		diag = parser.checkOption("--diag", "Write out diagnostic information");

		SH_bands = textToInteger(parser.getOption("--n", "Number of spherical harmonics bands", "2"));
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

	spheres = readSpheresCMM(spheresFn, spheres_binning);
	sphere_thickness = sphere_thickness_0 * spheres_binning;

	outPath = ZIO::makeOutputDir(outPath);
}

void DeleteBlobsProgram::run()
{
	TomogramSet tomogramSet = TomogramSet(tomoSetFn);
	const int tomo_index = tomogramSet.getTomogramIndexSafely(tomoName);

	Tomogram tomogram0 = tomogramSet.loadTomogram(tomo_index, true);

	std::vector<d3Vector> fiducials(0);
	bool has_fiducials = fiducialsDir.length() > 0;

	if (has_fiducials)
	{
		if (fiducialsDir[fiducialsDir.length()-1] != '/')
		{
			fiducialsDir = fiducialsDir + "/";
		}

		fiducials = Fiducials::read(tomogram0.name, fiducialsDir, tomogram0.optics.pixelSize);
	}

	processBlob(0, tomogram0, fiducials);
}

void DeleteBlobsProgram::processBlob(
		int blob_id,
		Tomogram& tomogram0,
		const std::vector<d3Vector>& fiducials)
{
	const double pixelSize_full = tomogram0.optics.pixelSize;
	const int frame_count = tomogram0.frameCount;

	// 3D coordinates live in bin-1 pixel coordinates
	const d4Vector sphere = spheres[blob_id];
	const d3Vector sphere_position(sphere.x, sphere.y, sphere.z);

	// 2D radii, thicknesses and distances are binned
	const double sphere_radius_full = sphere.w;
	const double outer_radius_full = sphere.w + sphere_thickness / 2;
	const double sphere_thickness_full = sphere_thickness;
	const double fiducials_radius_full = fiducials_radius_A / pixelSize_full;

	// The prior is evaluated in 2D
	const double prior_sigma = prior_sigma_A / pixelSize_full;
	const double initial_step = 1.0;

	const int test_frame = 19;


	const double binning_factor = 8.0;



	const int substack_size = 1024;
	const double pixelSize_binned = tomogram0.optics.pixelSize * binning_factor;
	const double fiducials_radius_binned = fiducials_radius_A / pixelSize_binned;

	const double sphere_radius_binned = sphere_radius_full / binning_factor;
	const double outer_radius_binned = outer_radius_full / binning_factor;

	const double sphere_thickness_binned = sphere_thickness / binning_factor;

	std::string outTag = outPath + "blob_" + ZIO::itoa(blob_id);
	std::string binTag = outPath + "bin_" + ZIO::itoa(binning_factor);
	std::string testFrameTag = "frame_" + ZIO::itoa(test_frame);


	Tomogram tomogram_cropped = tomogram0.extractSubstack(sphere_position, substack_size, substack_size);
	Tomogram tomogram_binned = tomogram_cropped.FourierCrop(binning_factor, num_threads);

	const int w_cropped = tomogram_cropped.stack.xdim;
	const int h_cropped = tomogram_cropped.stack.ydim;
	const int w_binned = tomogram_binned.stack.xdim;
	const int h_binned = tomogram_binned.stack.ydim;

	BufferedImage<float> dummyWeight_cropped(w_cropped, h_cropped);
	dummyWeight_cropped.fill(1.f);

	BufferedImage<float> dummyWeight_binned(w_binned, h_binned);
	dummyWeight_binned.fill(1.f);


	if (diag)
	{
		Blob blob(sphere_position, outer_radius_binned);

		tomogram_binned.stack.getSliceRef(test_frame)
			.write(outTag+"_000_original_"+testFrameTag+".mrc");

		drawTestFrame(blob, tomogram_binned, test_frame, dummyWeight_binned)
			.write(outTag+"_001_initial_"+testFrameTag+".mrc");
	}



	BlobFit bf0(
		tomogram_binned, spheres[blob_id].xyz(), 0,
		sphere_radius_binned,
		sphere_thickness_binned,
		spheres,
		fiducials, fiducials_radius_binned,
		prior_sigma, num_threads);

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
		Blob blob0(opt0, outer_radius_binned, 0);

		drawTestFrame(blob0, tomogram_binned, test_frame, dummyWeight_binned)
			.write(outTag+"_002_SH0_"+testFrameTag+".mrc");
	}

	for (int current_SH_bands = 1; current_SH_bands <= SH_bands; current_SH_bands++)
	{
		const int SH_coeffs = (current_SH_bands + 1) * (current_SH_bands + 1);

		/*!!BlobFit bf(
			tomogram_binned, spheres, blob_id, inner_radius,
			outer_radius_binned, current_SH_bands, use_masks, prior_sigma, num_threads);*/

		BlobFit bf(
			tomogram_binned, spheres[blob_id].xyz(), current_SH_bands,
			sphere_radius_binned,
			sphere_thickness_binned,
			spheres,
			fiducials, fiducials_radius_binned,
			prior_sigma, num_threads);

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
			Blob blob2(opt, outer_radius_binned, &sh);

			drawTestFrame(blob2, tomogram_binned, test_frame, dummyWeight_binned)
				.write(
					outPath+"ves-"+ZIO::itoa(blob_id)
					+"_rad_avg_2_L"+ZIO::itoa(current_SH_bands)
					+"_f"+ZIO::itoa(test_frame)+".mrc");
		}

		opt0 = opt;
	}


	std::vector<double> opt_upscaled = opt0;

	for (int i = 3; i < opt_upscaled.size(); i++)
	{
		opt_upscaled[i] *= binning_factor;
	}



	{
		tomogram_cropped.stack.getSliceRef(test_frame).write(outPath+"original_f"+ZIO::itoa(test_frame)+"_bin1.mrc");


		SphericalHarmonics sh(SH_bands);
		Blob blob(opt_upscaled, outer_radius_full, &sh);

		drawTestFrame(blob, tomogram_cropped, test_frame, dummyWeight_cropped)
			.write(outPath+"blob0_f"+ZIO::itoa(test_frame)+"_bin1.mrc");
	}


	const bool refine_full = false;

	std::vector<double> opt_fin;

	BufferedImage<float> full_weight;

	if (refine_full)
	{
		/*!!BlobFit bf(
			tomogram1, spheres, blob_id, inner_radius_0, outer_radius_full,
			SH_bands, use_masks, prior_sigma * binning_factor, num_threads);*/

		BlobFit bf(
			tomogram_cropped, spheres[blob_id].xyz(), SH_bands,
			sphere_radius_full,
			sphere_thickness_full,
			spheres,
			fiducials, fiducials_radius_full,
			prior_sigma, num_threads);

		full_weight = bf.weight;

		opt_fin = NelderMead::optimize(
			opt_upscaled, bf, initial_step, 0.05, 200, 1.0, 2.0, 0.5, 0.5, diag);


		{
			SphericalHarmonics sh(SH_bands);
			Blob blob(opt_fin, outer_radius_full, &sh);

			drawTestFrame(blob, tomogram_cropped, test_frame, dummyWeight_cropped)
				.write(outPath+"blob0_f"+ZIO::itoa(test_frame)+"_bin1_final.mrc");
		}
	}
	else
	{
		opt_fin = opt_upscaled;
	}




	SphericalHarmonics sh(SH_bands);

	Blob blobOut(opt_fin, outer_radius_full, &sh);

	for (int f = 0; f < frame_count; f++)
	{
		RawImage<float> slice = tomogram0.stack.getSliceRef(f);

		blobOut.subtract(
			slice,
			tomogram0.projectionMatrices[f],
			full_weight.getSliceRef(f),
			50.0);
	}

	tomogram0.stack.write(outPath+"blob_"+ZIO::itoa(blob_id)+"_stack.mrc", tomogram0.optics.pixelSize);
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

BufferedImage<float> DeleteBlobsProgram::drawTestFrame(
		Blob& blob,
		const Tomogram& tomogram,
		int test_frame,
		BufferedImage<float>& dummyWeight)
{
	std::vector<double> radAvg = blob.radialAverage(
			tomogram.stack.getConstSliceRef(test_frame),
			tomogram.projectionMatrices[test_frame],
			dummyWeight);

	BufferedImage<float> radAvgProj = blob.radialAverageProjection(
			tomogram.stack.getConstSliceRef(test_frame),
			tomogram.projectionMatrices[test_frame],
			radAvg);

	return radAvgProj;
}
