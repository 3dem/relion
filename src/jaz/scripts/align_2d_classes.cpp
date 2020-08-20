#include <string>
#include <src/metadata_table.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/single_particle/obs_model.h>
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
#include <omp.h>


using namespace gravis;


BufferedImage<float> rotate(const RawImage<float>& img, double phi);
BufferedImage<float> flip_x(const RawImage<float>& img);
BufferedImage<float> correlate(const RawImage<float>& img0, const RawImage<float>& img1);
float max_correlation(const RawImage<float>& img0, const RawImage<float>& img1);

int main(int argc, char *argv[])
{
	std::string particlesFn, class_averages_filename, outDir;
	double margin;
	int num_threads;
	
	
	IOParser parser;
	
	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		particlesFn = parser.getOption("--i", "Input file (e.g. run_it023_data.star)");
		class_averages_filename = parser.getOption("--ca", "Class averages stack");
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
		Log::warn("1 class found");
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

	std::vector<std::vector<int>> particle_by_class(class_count);

	for (long int c = 0; c < class_count; c++)
	{
		particle_by_class[c].reserve(particle_count[c]);
	}

	for (long int p = 0; p < particles_table.numberOfObjects(); p++)
	{
		const int class_id = particles_table.getIntMinusOne(
					EMDL_PARTICLE_CLASS, p);

		particle_by_class[class_id].push_back(p);
	}


	const int box_size = obs_model.getBoxSize(0);
	const double pixel_size = obs_model.getPixelSize(0);
	
	
	int biggest_class = -1;
	int biggest_class_size = 0;
	
	for (long int c = 0; c < class_count; c++)
	{
		const int n = particle_count[c];
		
		if (n > biggest_class_size)
		{
			biggest_class_size = n;
			biggest_class = c;
		}
	}
	
	BufferedImage<float> class_averages;
	class_averages.read(class_averages_filename);
	
	// align biggest class so as to maximise horizontal symmetry
	{
		const int angle_samples = 100;
		
		double best_phi = 0;
		double best_phi_CC = 0;
		
		RawImage<float> slice0 = class_averages.getSliceRef(biggest_class);
		
		
		for (int ai = 0; ai < angle_samples; ai++)
		{
			const double phi = ai * 2.0 * PI / angle_samples;
						
			BufferedImage<float> slice = rotate(slice0, phi);
			BufferedImage<float> flipped_slice = flip_x(slice0);
			
			const double max_CC = max_correlation(slice, flipped_slice);
			
			if (max_CC > best_phi_CC)
			{
				best_phi_CC = max_CC;
				best_phi = phi;
			}
		}
		
		BufferedImage<float> rotated_slice = rotate(slice0, best_phi);
		
		slice0.write("debug_class_align__original.mrc");
		rotated_slice.write("debug_class_align__rotated.mrc");
	}
	
	return 0;
}



BufferedImage<float> rotate(const RawImage<float>& img, double phi)
{
	const int s = img.xdim;
	const int m = s / 2;
	
	const d2Matrix A(
			 cos(phi), sin(phi),
	        -sin(phi), cos(phi) );
	
	BufferedImage<float> out(s,s);
	
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const d2Vector r0(x - m, y - m);
		const d2Vector r1 = A * r0;
		
		out(x,y) = Interpolation::linearXY_clip(img, m + r1.x, m + r1.y);
	}
	
	return out;
}

BufferedImage<float> flip_x(const RawImage<float>& img)
{
	const int s = img.xdim;
	
	BufferedImage<float> out(s,s);
	
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		out(x,y) = img(s-x-1, y);
	}
	
	return out;
}

BufferedImage<float> correlate(const RawImage<float>& img0, const RawImage<float>& img1)
{
	const int s  = img0.xdim;
	const int sh = s/2 + 1;
	
	BufferedImage<fComplex> img0_FS, img1_FS;
	
	FFT::FourierTransform(img0, img0_FS, FFT::Both);
	FFT::FourierTransform(img1, img1_FS, FFT::Both);
	
	BufferedImage<fComplex> product(sh,s);
	
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		product(x,y) = img0_FS(x,y) * img1_FS(x,y).conj();
	}
	
	BufferedImage<float> CC;
	
	FFT::inverseFourierTransform(product, CC, FFT::Both);
	
	return CC;
}

float max_correlation(const RawImage<float>& img0, const RawImage<float>& img1)
{
	BufferedImage<float> CC = correlate(img0, img1);	
	return Interpolation::discreteMaxXYV(CC)[2];
}
