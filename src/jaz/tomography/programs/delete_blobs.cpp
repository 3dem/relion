#include "delete_blobs.h"

#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/membrane/blob_fit_3d.h>
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
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	sphere_thickness = sphere_thickness_0 * spheres_binning;

	outPath = ZIO::makeOutputDir(outPath);

	if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag");
	}
}

void DeleteBlobsProgram::run()
{
	TomogramSet initial_tomogram_set = TomogramSet(tomoSetFn);
	TomogramSet subtracted_tomogram_set = initial_tomogram_set;
	TomogramSet blobs_tomogram_set = initial_tomogram_set;
	
	if (!initial_tomogram_set.globalTable.containsLabel(EMDL_TOMO_FIDUCIALS_STARFILE))
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

void DeleteBlobsProgram::processTomogram(
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


	Log::print("Filtering");
	
	BufferedImage<float> original_stack = tomogram0.stack;

	tomogram0.stack = ImageFilter::highpassStackGaussPadded(
						original_stack,
						highpass_sigma_real,
						num_threads);
	if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag/" + tomogram0.name);

		if (visualisation.xdim == 0)
		{
			visualisation.resize(w_full, h_full, 3 * tomo_batch_size);
		}

		visualisation.getSliceRef(3 * tomo_batch_index).copyFrom(
				original_stack.getConstSliceRef(frame_count / 2));
	}


	std::vector<d3Vector> fiducials(0);
	bool has_fiducials = tomogram0.fiducialsFilename.length() > 0;

	BufferedImage<float> fiducialsMask(w_full, h_full, frame_count);

	if (has_fiducials)
	{
		Log::print("Loading fiducial markers");
		
		if (fiducialsDir[fiducialsDir.length()-1] != '/')
		{
			fiducialsDir = fiducialsDir + "/";
		}

		fiducials = Fiducials::read(tomogram0.fiducialsFilename, tomogram0.optics.pixelSize);

		#pragma omp parallel for num_threads(num_threads)
		for (int f = 0; f < frame_count; f++)
		{
			RawImage<float> maskSlice = fiducialsMask.getSliceRef(f);
			maskSlice.fill(1.f);

			Fiducials::drawMask(
					fiducials, tomogram0.projectionMatrices[f],
					fiducials_radius, maskSlice, 1e-5f);
		}
	}
	
	Log::print("Allocating blobs stack");
	
	BufferedImage<float> blobs_stack(w_full, h_full, frame_count);
	blobs_stack.fill(0.f);


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


		Tomogram tomogram_cropped = tomogram0.extractSubstack(sphere_position, substack_size, substack_size);


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
		}
		
		Log::endSection();
	}

	const std::string tag = outPath + tomogram0.name;

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
	}
}

std::vector<double> DeleteBlobsProgram::fitBlob(
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


BufferedImage<float> DeleteBlobsProgram::drawFit(
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


BufferedImage<float> DeleteBlobsProgram::drawTestStack(
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

std::vector<double> DeleteBlobsProgram::toBin1(
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

std::vector<double> DeleteBlobsProgram::fromBin1(
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
