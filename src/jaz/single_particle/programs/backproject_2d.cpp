#include "backproject_2d.h"
#include <src/args.h>
#include <src/metadata_table.h>
#include <src/ctf.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/translation.h>
#include <src/jaz/image/tapering.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/util/image_file_helper.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <omp.h>


using namespace gravis;


void Backproject2D::read(int argc, char **argv)
{
	IOParser parser;

	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		particlesFn = parser.getOption("--i", "Input file (e.g. run_it023_data.star)", "");
		reextract = parser.checkOption("--reextract", "Extract particles from the micrographs");
		do_dual_contrast = parser.checkOption("--dual_contrast", "Perform a dual-contrast reconstruction");
		SNR = textToDouble(parser.getOption("--SNR", "Assumed signal-to-noise ratio", "0.1"));
		margin = textToDouble(parser.getOption("--m", "Margin around the particle [Px]", "20"));
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		outDir = parser.getOption("--o", "Output directory");

		Log::readParams(parser);

		if (parser.checkForErrors())
		{
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
		}
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
}

void Backproject2D::run()
{
	outDir = ZIO::makeOutputDir(outDir);

	ObservationModel obs_model;
	MetaDataTable particles_table;
	ObservationModel::loadSafely(particlesFn, obs_model, particles_table);

	int max_class = -1;

	for (long int p = 0; p < particles_table.numberOfObjects(); p++)
	{
		const int class_id = particles_table.getIntMinusOne(EMDL_PARTICLE_CLASS, p);

		if (class_id > max_class) max_class = class_id;
	}

	const int class_count = max_class + 1;

	if (class_count == 1)
	{
		Log::print("1 class found");
	}
	else
	{
		Log::print(ZIO::itoa(class_count)+" classes found");
	}


	std::vector<int> class_size(class_count, 0);

	for (long int p = 0; p < particles_table.numberOfObjects(); p++)
	{
		const int class_id = particles_table.getIntMinusOne(
					EMDL_PARTICLE_CLASS, p);

		class_size[class_id]++;
	}



	const int box_size = obs_model.getBoxSize(0);
	const double pixel_size = obs_model.getPixelSize(0);



	BufferedImage<dComplex> data;
	BufferedImage<double> weight;
	BufferedImage<DualContrastVoxel<double>> data_dual_contrast;

	if (do_dual_contrast)
	{
		data_dual_contrast = BufferedImage<DualContrastVoxel<double>>
				(box_size / 2 + 1, box_size, num_threads * class_count);
	}
	else
	{
		data = BufferedImage<dComplex>(box_size / 2 + 1, box_size, num_threads * class_count);
		weight = BufferedImage<double>(box_size / 2 + 1, box_size, num_threads * class_count);

		data.fill(dComplex(0.0, 0.0));
		weight.fill(0.0);
	}

	std::vector<MetaDataTable> particles_by_micrograph = StackHelper::splitByMicrographName(particles_table);

	const int micrograph_count = particles_by_micrograph.size();

	for (int micrograph_id = 0; micrograph_id < micrograph_count; micrograph_id++)
	{
		Log::print("Micrograph " + ZIO::itoa(micrograph_id+1));

		const MetaDataTable& particles = particles_by_micrograph[micrograph_id];
		const int particle_count = particles.numberOfObjects();


		BufferedImage<float> micrograph;
		float mean_value = 0, std_dev = 1;
		double micrograph_pixel_size = pixel_size;
		double extraction_scale = 1;

		if (reextract)
		{
			const std::string micrograph_filename = particles.getString(EMDL_MICROGRAPH_NAME, 0);
			micrograph_pixel_size = ImageFileHelper::getSamplingRate(micrograph_filename);
			extraction_scale = pixel_size / micrograph_pixel_size;

			// Round extraction scale to the third decimal to avoid numbers
			// that are infinitesimally shy of the next integer
			extraction_scale = std::round(1000 * extraction_scale) / 1000.0;

			micrograph.read(micrograph_filename);

			mean_value = Normalization::computeMean(micrograph);
			std_dev = sqrt(Normalization::computeVariance(micrograph, mean_value));

		}

		const i2Vector micrograph_size(micrograph.xdim, micrograph.ydim);


		#pragma omp parallel for num_threads(num_threads)
		for (long int p = 0; p < particle_count; p++)
		{
			const int thread_id = omp_get_thread_num();
			const int class_id = particles.getIntMinusOne(EMDL_PARTICLE_CLASS, p);
			const int slice_id = thread_id * class_count + class_id;

			BufferedImage<float> particle_image_RS;
			BufferedImage<fComplex> particle_image_FS;

			const double dx_A = particles.getDouble(EMDL_ORIENT_ORIGIN_X_ANGSTROM, p);
			const double dy_A = particles.getDouble(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, p);

			d2Vector shift;

			if (reextract)
			{
				const int extraction_box_size = std::round(extraction_scale * box_size);

				const d2Vector local_shift = d2Vector(dx_A, dy_A) / pixel_size;

				const d2Vector global_position_0(
					particles.getDouble(EMDL_IMAGE_COORD_X, p),
					particles.getDouble(EMDL_IMAGE_COORD_Y, p));

				const d2Vector global_position = global_position_0 - extraction_scale * local_shift;

				i2Vector integral_position(std::round(global_position.x), std::round(global_position.y));

				for (int dim = 0; dim < 2; dim++)
				{
					shift[dim] = (integral_position[dim] - global_position[dim]) / extraction_scale;
				}

				BufferedImage<float> extraction_buffer(extraction_box_size, extraction_box_size);

				const int x0 = integral_position.x - extraction_box_size / 2;
				const int y0 = integral_position.y - extraction_box_size / 2;

				double sum = 0.0;

				for (int y = 0; y < extraction_box_size; y++)
				for (int x = 0; x < extraction_box_size; x++)
				{
					int xx = x0 + x;
					int yy = y0 + y;

					if (xx < 0) xx = 0;
					else if (xx >= micrograph_size.x) xx = micrograph_size.x - 1;

					if (yy < 0) yy = 0;
					else if (yy >= micrograph_size.y) yy = micrograph_size.y - 1;

					extraction_buffer(x,y) = -micrograph(xx, yy);

					sum += extraction_buffer(x,y);
				}

				extraction_buffer -= sum / (extraction_box_size * extraction_box_size);
				extraction_buffer /= std_dev;


				if (std::abs(micrograph_pixel_size - pixel_size) > 0.001)
				{
					particle_image_RS = Resampling::FourierCrop_fullStack(
								extraction_buffer, extraction_scale, num_threads, true);
				}
				else
				{
					particle_image_RS = extraction_buffer;
				}
			}
			else
			{
				shift = d2Vector(dx_A, dy_A) / pixel_size;

				std::string img_fn = particles.getString(EMDL_IMAGE_NAME, p);
				particle_image_RS.read(img_fn);
			}

			FFT::FourierTransform(particle_image_RS, particle_image_FS, FFT::Both);

			particle_image_FS(0,0) = 0;

			const double m = box_size / 2;
			Translation::shiftInFourierSpace2D(particle_image_FS, shift.x + m, shift.y + m);

			if (do_dual_contrast)
			{
				backrotate_particle_dual_contrast(
					particle_image_FS,
					p, particles,
					obs_model,
					data_dual_contrast.getSliceRef(slice_id));
			}
			else
			{
				backrotate_particle(
					particle_image_FS,
					p, particles,
					obs_model,
					data.getSliceRef(slice_id),
					weight.getSliceRef(slice_id));
			}
		}
	}

	for (int class_id = 0; class_id < class_count; class_id++)
	{
		const int slice_id_0 = class_id;

		for (long int t = 1; t < num_threads; t++)
		{
			const int slice_id = t * class_count + class_id;

			if (do_dual_contrast)
			{
				data_dual_contrast.getSliceRef(slice_id_0) += data_dual_contrast.getSliceRef(slice_id);
			}
			else
			{
				data.getSliceRef(slice_id_0) += data.getSliceRef(slice_id);
				weight.getSliceRef(slice_id_0) += weight.getSliceRef(slice_id);
			}
		}
	}

	if (do_dual_contrast)
	{
		BufferedImage<double> average_phase_stack(box_size, box_size, class_count);
		BufferedImage<double> average_amp_stack(box_size, box_size, class_count);

		for (int class_id = 0; class_id < class_count; class_id++)
		{
			BufferedImage<double> average_phase, average_amp;

			std::pair<BufferedImage<double>, BufferedImage<double>> average =
				reconstruct_dual_contrast(
						data_dual_contrast.getSliceRef(class_id),
						1.0/SNR);

			average_phase = average.first;
			average_amp = average.second;

			const double radius = box_size/2;

			Tapering::taperCircularly2D(average_phase, radius - margin, radius - margin + 5);
			Tapering::taperCircularly2D(average_amp, radius - margin, radius - margin + 5);

			average_phase_stack.getSliceRef(class_id).copyFrom(average_phase);
			average_amp_stack.getSliceRef(class_id).copyFrom(average_amp);
		}

		average_phase_stack.write(outDir + "class_averages_phase.mrc", pixel_size);
		average_amp_stack.write(outDir + "class_averages_amplitude.mrc", pixel_size);
	}
	else
	{
		BufferedImage<double> average_stack(box_size, box_size, class_count);

		for (int class_id = 0; class_id < class_count; class_id++)
		{
			BufferedImage<double> average;

			average = reconstruct(
						data.getSliceRef(class_id),
						weight.getSliceRef(class_id),
						1.0/SNR);

			const double radius = box_size/2;

			Tapering::taperCircularly2D(average, radius - margin, radius - margin + 5);

			average_stack.getSliceRef(class_id).copyFrom(average);
		}

		average_stack.write(outDir + "class_averages.mrc", pixel_size);
	}
}

void Backproject2D::backrotate_particle(
		const RawImage<fComplex>& image,
		long int particle_id,
		const MetaDataTable& particles_table,
		ObservationModel& obsModel,
		RawImage<dComplex> data,
		RawImage<double> weight)
{
	const int sh = data.xdim;
	const int s  = data.ydim;

	const double pixel_size = obsModel.getPixelSize(0);
	const double box_size_px = obsModel.getBoxSize(0);
	const double box_size_A = box_size_px * pixel_size;


	const double psi = DEG2RAD(particles_table.getDouble(EMDL_ORIENT_PSI, particle_id));

	const d2Matrix rot(
			 cos(psi), sin(psi),
			-sin(psi), cos(psi)	);

	CTF ctf;
	ctf.readByGroup(particles_table, &obsModel, particle_id);


	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		const d2Vector p0(x, (y < s/2? y : y - s));
		const d2Vector p1 = rot * p0;

		const fComplex z = Interpolation::linearXY_complex_FftwHalf_wrap(
					image, p1.x, p1.y);

		const double c = ctf.getCTF(p1.x / box_size_A, p1.y / box_size_A);

		data(x,y)   += c * z;
		weight(x,y) += c * c;
	}
}

void Backproject2D::backrotate_particle_dual_contrast(
		const RawImage<fComplex>& image,
		long particle_id,
		const MetaDataTable& particles_table,
		ObservationModel& obsModel,
		RawImage<DualContrastVoxel<double>> average)
{
	const int sh = average.xdim;
	const int s  = average.ydim;

	const double pixel_size = obsModel.getPixelSize(0);
	const double box_size_px = obsModel.getBoxSize(0);
	const double box_size_A = box_size_px * pixel_size;


	const double psi = DEG2RAD(particles_table.getDouble(EMDL_ORIENT_PSI, particle_id));

	const d2Matrix rot(
			 cos(psi), sin(psi),
			-sin(psi), cos(psi)	);

	CTF ctf;
	ctf.readByGroup(particles_table, &obsModel, particle_id);
	ctf.Q0 = 0.0;
	ctf.initialise();


	BufferedImage<fComplex>
		data_sin(sh,s),
		data_cos(sh,s);

	BufferedImage<float>
		weight_sin2(sh,s),
		weight_sin_cos(sh,s),
		weight_cos2(sh,s);

	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		const double xx = x / box_size_A;
		const double yy = (y < s/2? y : y - s) / box_size_A;

		const double gamma = ctf.getLowOrderGamma(xx,yy);
		const fComplex z = image(x,y);

		data_sin(x,y) = sin(gamma) * z;
		data_cos(x,y) = cos(gamma) * z;

		weight_sin2(x,y)    = sin(gamma) * sin(gamma);
		weight_sin_cos(x,y) = sin(gamma) * cos(gamma);
		weight_cos2(x,y)    = cos(gamma) * cos(gamma);
	}

	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		const d2Vector p0(x, (y < s/2? y : y - s));
		const d2Vector p1 = rot * p0;

		DualContrastVoxel<double>& t = average(x,y);

		t.data_sin += Interpolation::linearXY_complex_FftwHalf_wrap(data_sin, p1.x, p1.y);
		t.data_cos += Interpolation::linearXY_complex_FftwHalf_wrap(data_cos, p1.x, p1.y);

		t.weight_sin2    += Interpolation::linearXY_symmetric_FftwHalf_wrap(weight_sin2,    p1.x, p1.y);
		t.weight_sin_cos += Interpolation::linearXY_symmetric_FftwHalf_wrap(weight_sin_cos, p1.x, p1.y);
		t.weight_cos2    += Interpolation::linearXY_symmetric_FftwHalf_wrap(weight_cos2,    p1.x, p1.y);
	}
}


BufferedImage<double> Backproject2D::reconstruct(
		const RawImage<dComplex>& data,
		const RawImage<double>& weight,
		double Wiener_offset)
{
	const int sh = data.xdim;
	const int s  = data.ydim;

	BufferedImage<double> out_RS(s,s);
	BufferedImage<dComplex> out_FS(sh,s);

	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		const double mod = (1 - 2 * (x % 2)) * (1 - 2 * (y % 2));

		out_FS(x,y) = mod * data(x,y) / (weight(x,y) + Wiener_offset);
	}

	FFT::inverseFourierTransform(out_FS, out_RS, FFT::Both);

	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const double xx = x - s/2;
		const double yy = y - s/2;

		if (xx == 0 && yy == 0)
		{
			// sinc at 0 is 1
		}
		else
		{
			const double r = sqrt(xx*xx + yy*yy);
			const double d = r / s;
			const double sinc = sin(PI * d) / (PI * d);
			const double sinc2 = sinc * sinc;

			if (d < 0.99)
			{
				out_RS(x,y) /= sinc2;
			}
			else
			{
				out_RS(x,y) = 0.0;
			}
		}
	}

	return out_RS;
}

std::pair<BufferedImage<double>, BufferedImage<double>>
	Backproject2D::reconstruct_dual_contrast(
		const RawImage<DualContrastVoxel<double>> &data,
		double Wiener_offset)
{
	const int sh = data.xdim;
	const int s  = data.ydim;

	BufferedImage<double> out_phase_RS(s,s), out_amp_RS(s,s);
	BufferedImage<dComplex> out_phase_FS(sh,s), out_amp_FS(sh,s);

	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		const double mod = (1 - 2 * (x % 2)) * (1 - 2 * (y % 2));

		DualContrastVoxel<double>::Solution d = data(x,y).solve(Wiener_offset, 0.0);

		out_phase_FS(x,y) = mod * d.phase;
		out_amp_FS(x,y) = mod * d.amplitude;
	}

	FFT::inverseFourierTransform(out_phase_FS, out_phase_RS, FFT::Both);
	FFT::inverseFourierTransform(out_amp_FS, out_amp_RS, FFT::Both);

	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const double xx = x - s/2;
		const double yy = y - s/2;

		if (xx == 0 && yy == 0)
		{
			// sinc at 0 is 1
		}
		else
		{
			const double r = sqrt(xx*xx + yy*yy);
			const double d = r / s;
			const double sinc = sin(PI * d) / (PI * d);
			const double sinc2 = sinc * sinc;

			if (d < 0.99)
			{
				out_phase_RS(x,y) /= sinc2;
				out_amp_RS(x,y) /= sinc2;
			}
			else
			{
				out_phase_RS(x,y) = 0.0;
				out_amp_RS(x,y) = 0.0;
			}
		}
	}

	return std::make_pair(out_phase_RS, out_amp_RS);
}
