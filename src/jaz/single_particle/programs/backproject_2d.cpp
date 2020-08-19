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
#include <src/jaz/math/fft.h>


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
		SNR = textToDouble(parser.getOption("--SNR", "Assumed signal-to-noise ratio", "0.1"));
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


	std::vector<int> particle_count(class_count, 0);

	for (long int p = 0; p < particles_table.numberOfObjects(); p++)
	{
		const int class_id = particles_table.getIntMinusOne(
					EMDL_PARTICLE_CLASS, p);

		particle_count[class_id]++;
	}

	const int box_size = obs_model.getBoxSize(0);
	const double pixel_size = obs_model.getPixelSize(0);


	BufferedImage<double> average_stack(box_size, box_size, class_count);
	BufferedImage<float> particle_image_RS;
	BufferedImage<fComplex> particle_image_FS;

	for (int class_id = 0; class_id < class_count; class_id++)
	{
		BufferedImage<dComplex> data(box_size / 2 + 1, box_size);
		data.fill(dComplex(0.0, 0.0));

		BufferedImage<double> weight(box_size / 2 + 1, box_size);
		weight.fill(0.0);

		for (long int p = 0; p < particles_table.numberOfObjects(); p++)
		{
			const int class_id_p = particles_table.getIntMinusOne(
						EMDL_PARTICLE_CLASS, p);

			if (class_id_p == class_id)
			{
				const double dx_A = particles_table.getDouble(
							EMDL_ORIENT_ORIGIN_X_ANGSTROM, p);

				const double dy_A = particles_table.getDouble(
							EMDL_ORIENT_ORIGIN_Y_ANGSTROM, p);

				const d2Vector shift(dx_A / pixel_size, dy_A / pixel_size);

				std::string img_fn = particles_table.getString(EMDL_IMAGE_NAME, p);
				particle_image_RS.read(img_fn);

				FFT::FourierTransform(particle_image_RS, particle_image_FS, FFT::Both);

				const double m = box_size / 2;
				Translation::shiftInFourierSpace2D(particle_image_FS, shift.x + m, shift.y + m);

				backrotate_particle(
					particle_image_FS,
					p, particles_table,
					obs_model,
					data, weight);
			}
		}

		BufferedImage<double> average = reconstruct(data, weight, 1.0/SNR);

		Tapering::taperCircularly2D(average, 5);

		average_stack.getSliceRef(class_id).copyFrom(average);
	}

	average_stack.write(outDir + "class_averages.mrc", pixel_size);

}

void Backproject2D::backrotate_particle(
		const RawImage<fComplex> image,
		long int particle_id,
		const MetaDataTable& particles_table,
		ObservationModel& obsModel,
		RawImage<dComplex>& data,
		RawImage<double>& weight)
{
	const int sh = data.xdim;
	const int s  = data.ydim;

	const double pixel_size = obsModel.getPixelSize(0);
	const double box_size_px = obsModel.getBoxSize(0);
	const double box_size_A = box_size_px * pixel_size;


	const double psi = DEG2RAD(particles_table.getDouble(EMDL_ORIENT_PSI));

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


BufferedImage<double> Backproject2D::reconstruct(
		RawImage<dComplex>& data,
		RawImage<double>& weight,
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
