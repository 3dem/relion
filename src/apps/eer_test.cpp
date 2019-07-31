#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <omp.h>

#include <src/time.h>
#include <src/metadata_table.h>
#include <src/args.h>
#include <src/image.h>
#include <src/renderEER.h>


class eer_test {
	public:

	IOParser parser;

	FileName fn_movie, fn_out;
	int n_threads, max_frames;

	void read(int argc, char **argv) {
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("Options");
		fn_movie = parser.getOption("--i", "Input EER file");
		fn_out = parser.getOption("--o", "Output file name");
		n_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));
		max_frames = textToInteger(parser.getOption("--max_frames", "Target number of frames to average (rounded to movies; -1 means use all)", "-1"));

		if (parser.checkForErrors()) {
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
		}  
	}

	void run() {
		EERRenderer renderer(fn_movie);
		Image<float> image;

		int frames_to_read = renderer.getNFrames();
		if (max_frames > 0)
			frames_to_read = XMIPP_MIN(frames_to_read, max_frames);

		renderer.renderFrames(1, frames_to_read, image());
		image.write(fn_out);
	}
};

int main(int argc, char **argv) {
	eer_test app;
	app.read(argc, argv);
	app.run();

#ifdef TIMING
	timer.printTimes(false);
#endif

	return RELION_EXIT_SUCCESS;
}
