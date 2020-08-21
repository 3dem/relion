#include <string>
#include <src/metadata_table.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/args.h>
#include <src/metadata_table.h>
#include <src/ctf.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/index_sort.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/translation.h>
#include <src/jaz/image/tapering.h>
#include <src/jaz/image/centering.h>
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
	bool flip_orientation, reorder, align_horizontally_only;
	
	
	IOParser parser;
	
	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		particlesFn = parser.getOption("--i", "Input file (e.g. run_it023_data.star)");
		class_averages_filename = parser.getOption("--ca", "Class averages stack");
		flip_orientation = parser.checkOption("--flip", "Rotate the initial alignment by 180°");
		align_horizontally_only = parser.checkOption("--horizontal_only", "Only align the biggest class horizontally, then exit");
		reorder = parser.checkOption("--sort", "Sort classes by descending cardinality");
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

	std::vector<int> class_by_size = IndexSort<int>::sortIndices(particle_count);
	const int biggest_class = class_by_size[class_count-1];


	const int angle_samples_coarse = 36;
	const int angle_samples_fine = 100;


	BufferedImage<float> class_averages;
	class_averages.read(class_averages_filename);


	std::vector<double> angle_by_class(class_count, 0.0);

	BufferedImage<float> template_slice;


	// align biggest class to maximise horizontal symmetry
	{

		
		double best_phi_coarse = 0;
		double best_phi_CC_coarse = 0;
		
		RawImage<float> slice0 = class_averages.getSliceRef(biggest_class);


		for (int ai = 0; ai < angle_samples_coarse; ai++)
		{
			const double phi = ai * 2.0 * PI / angle_samples_coarse;

			BufferedImage<float> slice = rotate(slice0, phi);
			BufferedImage<float> flipped_slice = flip_x(slice);

			const double max_CC = max_correlation(slice, flipped_slice);

			if (max_CC > best_phi_CC_coarse)
			{
				best_phi_CC_coarse = max_CC;
				best_phi_coarse = phi;
			}
		}

		double best_phi_fine = best_phi_coarse;
		double best_phi_CC_fine = best_phi_CC_coarse;

		for (int ai = 0; ai < angle_samples_fine; ai++)
		{
			const double range = 2.0 * PI / angle_samples_coarse;
			const double phi = best_phi_coarse + (2 * ai / (double) angle_samples_fine - 1) * range;

			BufferedImage<float> slice = rotate(slice0, phi);
			BufferedImage<float> flipped_slice = flip_x(slice);

			const double max_CC = max_correlation(slice, flipped_slice);

			if (max_CC > best_phi_CC_fine)
			{
				best_phi_CC_fine = max_CC;
				best_phi_fine = phi;
			}
		}

		if (flip_orientation)
		{
			best_phi_fine += PI;

			if (best_phi_fine > 2 * PI)
			{
				best_phi_fine -= 2 * PI;
			}
		}

		angle_by_class[biggest_class] = best_phi_fine;

		template_slice = rotate(slice0, best_phi_fine);

		slice0.write(outDir + "class_"+ZIO::itoa(biggest_class)+"_original.mrc");
		template_slice.write(outDir + "class_"+ZIO::itoa(biggest_class)+"_rotated.mrc");
	}



	if (align_horizontally_only) return 0;


	BufferedImage<float> class_averages_rotated = class_averages;
	class_averages_rotated.getSliceRef(biggest_class).copyFrom(template_slice);


	// align every other class against the biggest one

	for (int ci = 1; ci < class_count; ci++)
	{
		const int class_id = class_by_size[class_count - ci - 1];
		
		Log::print("Aligning class " + ZIO::itoa(class_id)); 


		double best_phi_coarse = 0;
		double best_phi_CC_coarse = 0;

		RawImage<float> slice0 = class_averages.getSliceRef(class_id);

		for (int ai = 0; ai < angle_samples_coarse; ai++)
		{
			const double phi = ai * 2.0 * PI / angle_samples_coarse;

			BufferedImage<float> slice = rotate(slice0, phi);
			const double max_CC = max_correlation(slice, template_slice);

			if (max_CC > best_phi_CC_coarse)
			{
				best_phi_CC_coarse = max_CC;
				best_phi_coarse = phi;
			}
		}

		double best_phi_fine = best_phi_coarse;
		double best_phi_CC_fine = best_phi_CC_coarse;

		for (int ai = 0; ai < angle_samples_fine; ai++)
		{
			const double range = 2.0 * PI / angle_samples_coarse;
			const double phi = best_phi_coarse + (2 * ai / (double) angle_samples_fine - 1) * range;

			BufferedImage<float> slice = rotate(slice0, phi);
			const double max_CC = max_correlation(slice, template_slice);

			if (max_CC > best_phi_CC_fine)
			{
				best_phi_CC_fine = max_CC;
				best_phi_fine = phi;
			}
		}
		
		angle_by_class[class_id] = best_phi_fine;

		BufferedImage<float> slice = rotate(slice0, best_phi_fine);
		class_averages_rotated.getSliceRef(class_id).copyFrom(slice);
	}

	class_averages_rotated.write(outDir+"class_averages_rotated.mrc");

	
	
	for (long int p = 0; p < particles_table.numberOfObjects(); p++)
	{
		const int class_id = particles_table.getIntMinusOne(EMDL_PARTICLE_CLASS, p);
		const double phi0 = particles_table.getDouble(EMDL_ORIENT_PSI, p);

		double phi1 = phi0 + RAD2DEG(angle_by_class[class_id]);
		
		if (phi1 > 360) phi1 -= 360;
		
		particles_table.setValue(EMDL_ORIENT_PSI, phi1, p);
	}
	
	if (reorder)
	{
		std::vector<int> new_class_id(class_count);
		
		for (int i = 0; i < class_count; i++)
		{
			new_class_id[class_by_size[i]] = class_count - i - 1;
		}
		
		for (long int p = 0; p < particles_table.numberOfObjects(); p++)
		{
			const int class_id = particles_table.getIntMinusOne(EMDL_PARTICLE_CLASS, p);
			
			particles_table.setValue(EMDL_PARTICLE_CLASS, new_class_id[class_id] + 1, p);
		}
	}
	
	
	
	ObservationModel::saveNew(
		particles_table, 
		obs_model.opticsMdt, 
		outDir+"rotated_particles.star");
	
	
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
