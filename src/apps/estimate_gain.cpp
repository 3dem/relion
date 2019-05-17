#include <iostream>
#include <omp.h>

#include <src/time.h>
#include <src/metadata_table.h>
#include <src/args.h>
#include <src/image.h>
#include <src/jaz/obs_model.h>

#define RCTIC(label) (timer.tic(label))
#define RCTOC(label) (timer.toc(label))

Timer timer;
//int TIMING_READ_MOVIE = timer.setNew("read movie");

/* A simple program to estimate the gain reference by averaging all unaligned frames and taking the inverse.
 * This program should be used only as a last resort when the gain reference was lost or badly calibrated.
 */

class estimate_gain {
	public:

	IOParser parser;

	FileName fn_movie_star, fn_out;
	int n_threads;

	void read(int argc, char **argv) {
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("Options");
		fn_movie_star = parser.getOption("--i", "Input movie STAR file");
		fn_out = parser.getOption("--o", "Output file name");
		n_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));

		if (parser.checkForErrors()) {
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
		}  
	}

	void run() {
		MetaDataTable MDin;
		ObservationModel obsModel;
		ObservationModel::loadSafely(fn_movie_star, obsModel, MDin, "movies");

		const int n_movies = MDin.numberOfObjects();
		if (n_movies == 0)
			REPORT_ERROR("No movies in the input STAR file");
		else
			std::cout << "Number of movies in the STAR file: " << n_movies << std::endl;
		
		FileName fn_img;
		if (!MDin.getValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_img, 0))
			REPORT_ERROR("The input STAR file does not contain the rlnMicrographMovieName column.");

		Image<RFLOAT> Ihead;
		std::vector<Image<RFLOAT>> Isums(n_threads);
		int ny = 0, nx = 0, total_frames = 0;

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
		{
			MDin.getValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_img);
			Ihead.read(fn_img, false, -1, false, true); // select_img -1, mmap false, is_2D true
			if (current_object == 0)
			{
				ny = YSIZE(Ihead());
				nx = XSIZE(Ihead());
				std::cout << "The size of the input: NY = " << ny << " NX = " << nx << std::endl;
				std::cout << "Checking that all movies have the same size ..." << std::endl;
				init_progress_bar(n_movies);
			}
			else if (ny != YSIZE(Ihead()) || nx != XSIZE(Ihead()))
			{
				std::cout << "Movie: " << fn_img << std::endl;
				REPORT_ERROR("The size of the movie " + fn_img + " does not much the size of the others.");
			}
			progress_bar(current_object);
		}
		progress_bar(n_movies);

		for (int i = 0; i < n_threads; i++)
		{
			Isums[i]().resize(ny, nx);
			Isums[i]().initZeros();
		}

		std::cout << "Summing frames ... " << std::endl;
		init_progress_bar(n_movies);
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
		{
			MDin.getValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_img);
			Ihead.read(fn_img, false, -1, false, true);
			const int n_frames = NSIZE(Ihead());
			total_frames += n_frames;
			#pragma omp parallel for num_threads(n_threads)
			for (int iframe = 0; iframe < n_frames; iframe++) {
				Image<RFLOAT> Iframe;
				Iframe.read(fn_img, true, iframe, false, true); // mmap false, is_2D true
				const int tid = omp_get_thread_num();
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iframe()) {
					DIRECT_MULTIDIM_ELEM(Isums[tid](), n) += DIRECT_MULTIDIM_ELEM(Iframe(), n);
				}
			}
			progress_bar(current_object);
		}
		progress_bar(n_movies);
		std::cout << "Total number of frames: " << total_frames << std::endl;

		for (int i = 1; i < n_threads; i++)
		{
			#pragma omp parallel for num_threads(n_threads)
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Isums[0]())
			{
				DIRECT_MULTIDIM_ELEM(Isums[0](), n) += DIRECT_MULTIDIM_ELEM(Isums[i](), n);
			}
		}

		#pragma omp parallel for num_threads(n_threads)
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Isums[0]())
		{
			if (DIRECT_MULTIDIM_ELEM(Isums[0](), n) == 0)
				DIRECT_MULTIDIM_ELEM(Isums[0](), n) = 1.0;
			else
				DIRECT_MULTIDIM_ELEM(Isums[0](), n) = total_frames / DIRECT_MULTIDIM_ELEM(Isums[0](), n);
		}

		Isums[0].write(fn_out);
		std::cout << "Written the estimated gain to " << fn_out << std::endl;
	}
};

int main(int argc, char **argv) {
	estimate_gain app;
	app.read(argc, argv);
	app.run();
//	timer.printTimes(false);
	return RELION_EXIT_SUCCESS;
}
