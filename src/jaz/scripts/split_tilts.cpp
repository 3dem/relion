#include <src/args.h>
#include <src/metadata_table.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/util/image_file_helper.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>

using namespace gravis;


int main(int argc, char *argv[])
{
	IOParser parser;

	double data_scale;
	int num_threads;
	bool write_stack;

	std::string in_file, gain_file, out_path;
	OptimisationSet optimisation_set;

	try
	{
		parser.setCommandLine(argc, argv);

		optimisation_set.read(
			parser,
			true,            // optimisation set
			false,  false,   // particles
			true,   true,    // tomograms
			false,  false,   // trajectories
			false,  false,   // manifolds
			false,  false);  // reference

		in_file = parser.getOption("--ti", "STAR file mapping chronological tilt indices to their movie files");
		gain_file = parser.getOption("--g", "Gain reference file", "");
		data_scale = textToDouble(parser.getOption("--ds", "Scale gain-reference multiplied pixels by this before truncating them to 16 bits", "1000"));
		write_stack = !parser.checkOption("--no_stack", "Do not write out an image stack");
		out_path = parser.getOption("--o", "Output filename");
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));

		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	out_path = ZIO::prepareTomoOutputDirectory(out_path, argc, argv);
	ZIO::makeDir(out_path+"TiltSeries");

	TomogramSet tomogram_set = TomogramSet(optimisation_set.tomograms, true);

	std::vector<MetaDataTable> tables = MetaDataTable::readAll(in_file, 100);

	const int tc = tables.size();
	const int tc2 = tomogram_set.size();

	if (tc != tc2)
	{
		REPORT_ERROR_STR("Tomogram count mismatch: " << tc << " vs. " << tc2);
	}


	BufferedImage<float> gain;
	bool has_gain = gain_file != "";

	if (has_gain)
	{
		gain.read(gain_file);
	}

	TomogramSet new_tomogram_set = tomogram_set;

	Log::beginProgress("Processing tomograms", tc);

	for (int tt = 0; tt < tc; tt++)
	{
		Log::updateProgress(tt);

		MetaDataTable& table = tables[tt];
		const std::string tomo_name = table.getName();

		const int t = tomogram_set.getTomogramIndexSafely(tomo_name);
		const Tomogram tomogram = tomogram_set.loadTomogram(t, false);

		const int fc = tomogram.frameCount;
		const int fc2 = table.numberOfObjects();

		if (fc2 != fc)
		{
			REPORT_ERROR_STR(
				"Tilt number mismatch in " << table.getName()
				<< ": " << fc2 << " vs. " << fc);
		}


		int new_frame_count = 0;
		i2Vector movie_size(-1,-1);
		std::vector<int> frames_by_tilt(fc);

		double pixel_size = 1;

		for (int ft = 0; ft < fc; ft++)
		{
			const int f = tomogram.frameSequence[ft];

			const std::string movie_file = table.getString(EMDL_TOMO_TILT_MOVIE_FILE_NAME, ft);
			t3Vector<long int> size = ImageFileHelper::getSize(movie_file);

			if (ft == 0)
			{
				pixel_size = ImageFileHelper::getSamplingRate(movie_file);
			}

			new_frame_count += size.z;

			if (movie_size.x < 0)
			{
				movie_size.x = size.x;
				movie_size.y = size.y;
			}
			else if (movie_size.x != size.x || movie_size.y != size.y)
			{
				REPORT_ERROR_STR(
					"Movie size mismatch in " << table.getName()
					<< ": " << movie_size.x << 'x' << movie_size.y
					<< " vs. " << size.x << 'x' << size.y);
			}

			frames_by_tilt[f] = size.z;
		}

		if (has_gain && (gain.xdim != movie_size.x || gain.ydim != movie_size.y))
		{
			REPORT_ERROR_STR(
				"Gain-reference size mismatch in " << table.getName()
				<< ": " << movie_size.x << 'x' << movie_size.y
				<< " vs. " << gain.xdim << 'x' << gain.ydim);
		}

		std::vector<int> tilt_to_frame(fc);
		int current_index = 0;

		for (int f = 0; f < fc; f++)
		{
			tilt_to_frame[f] = current_index;
			current_index += frames_by_tilt[f];
		}



		BufferedImage<short int> new_stack(movie_size.x, movie_size.y, new_frame_count);
		const std::string out_file = out_path + "TiltSeries/" + tomo_name + ".mrc";

		if (write_stack)
		{
			new_stack.fill(0.f);

			#pragma omp parallel for num_threads(num_threads)
			for (int ft = 0; ft < fc; ft++)
			{
				const int f = tomogram.frameSequence[ft];

				BufferedImage<float> movie;
				const std::string movie_file = table.getString(EMDL_TOMO_TILT_MOVIE_FILE_NAME, ft);
				movie.read(movie_file);

				for (int ff = 0; ff < movie.zdim; ff++)
				{
					if (has_gain)
					{
						for (size_t y = 0; y < movie.ydim; y++)
						for (size_t x = 0; x < movie.xdim; x++)
						{
							new_stack(x, y, tilt_to_frame[f] + ff) = data_scale * (gain(x,y) * movie(x,y,ff));
						}
					}
					else
					{
						for (size_t y = 0; y < movie.ydim; y++)
						for (size_t x = 0; x < movie.xdim; x++)
						{
							new_stack(x, y, tilt_to_frame[f] + ff) = movie(x,y,ff);
						}
					}
				}
			}


			new_stack.writeNoVTK(out_file, pixel_size);
		}

		new_tomogram_set.globalTable.setValue(EMDL_TOMO_TILT_SERIES_NAME, out_file, t);
		new_tomogram_set.globalTable.setValue(EMDL_TOMO_FRAME_COUNT, new_frame_count, t);

		const double fdose0 = tomogram_set.globalTable.getDouble(EMDL_TOMO_IMPORT_FRACT_DOSE, t);
		new_tomogram_set.globalTable.setValue(EMDL_TOMO_IMPORT_FRACT_DOSE, fdose0 / frames_by_tilt[0], t);


		MetaDataTable& new_table = new_tomogram_set.tomogramTables[t];
		MetaDataTable& old_table = tomogram_set.tomogramTables[t];

		new_table.clear();
		new_table.setName(old_table.getName());

		for (int f = 0; f < fc; f++)
		{
			const double dose0 = old_table.getDouble(EMDL_MICROGRAPH_PRE_EXPOSURE, f);
			const double ddose = tomogram.fractionalDose / frames_by_tilt[f];

			for (int ff = 0; ff < frames_by_tilt[f]; ff++)
			{
				new_table.addObject(old_table.getObject(f));
				new_table.setValue(EMDL_MICROGRAPH_PRE_EXPOSURE, dose0 + ff * ddose);
			}
		}
	}

	Log::endProgress();

	optimisation_set.tomograms = out_path + "tomograms.star";
	new_tomogram_set.write(optimisation_set.tomograms);
	optimisation_set.write(out_path + "optimisation_set.star");

	return 0;
}
