#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <omp.h>

#include <src/time.h>
#include <src/metadata_table.h>
#include <src/args.h>
#include <src/image.h>
#include <src/renderEER.h>
#include <src/tiff_converter.h>

#ifdef HAVE_TIFF
#include <tiffio.h>
#endif

class EERtoTIFF {
	public:

	IOParser parser;

	FileName fn_movie, fn_out;
	int n_threads, max_frames, eer_upsampling, eer_grouping, compression;
	bool use_short;

	void read(int argc, char **argv) {
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("Options");
		fn_movie = parser.getOption("--i", "Input EER file");
		fn_out = parser.getOption("--o", "Output file name");
		n_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));
		max_frames = textToInteger(parser.getOption("--max_frames", "Target number of (internal) frames to process (-1 means use all)", "-1"));
		eer_grouping = textToInteger(parser.getOption("--eer_grouping", "EER grouping", "40"));
		eer_upsampling = textToInteger(parser.getOption("--eer_upsampling", "EER upsampling (1 = 4K or 2 = 8K)", "2"));
		if (eer_upsampling != 1 && eer_upsampling != 2)
			REPORT_ERROR("eer_upsampling must be 1 or 2");
		FileName fn_compression = parser.getOption("--compression", "compresion (none, lzw, deflate=zip)", "lzw");
		use_short = parser.checkOption("--short", "use unsigned short instead of signed byte");

		if (fn_compression == "none")
                	compression = COMPRESSION_NONE;
		else if (fn_compression == "lzw")
			compression = COMPRESSION_LZW;
		else if (fn_compression == "deflate" || fn_compression == "zip")
			compression = COMPRESSION_DEFLATE;
		else
			REPORT_ERROR("Unknown --compression type");

		if (parser.checkForErrors()) {
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
		}  
	}

	template <typename T>
	void convert() {
#ifndef HAVE_TIFF
		REPORT_ERROR("To use this program, you have to re-compile RELION with libtiff.");
#else
		EERRenderer renderer;
		renderer.read(fn_movie, eer_upsampling);

		int nframes = renderer.getNFrames();
		std::cout << "Found " << nframes << " frames" << std::endl;

		TIFF *tif = TIFFOpen(fn_out.c_str(), "w");
		if (tif == NULL)
                	REPORT_ERROR("Failed to open the output TIFF file.");

		MultidimArray<T> buf;
		if (max_frames < 0)
			max_frames = nframes;
		for (int frame = 1; frame < max_frames; frame += eer_grouping)
		{
			const int frame_end = frame + eer_grouping - 1;
			if (frame_end > nframes)
				break;

			std::cout << "Writing frame " << frame << " to " << frame_end << std::endl;
			buf.initZeros(renderer.getHeight(), renderer.getWidth());
			renderer.renderFrames(frame, frame_end, buf);
			TIFFConverter::write_tiff_one_page(tif, buf, -1, compression);		
		}

		TIFFClose(tif);
	}

	void run()
	{
		if (use_short)
			convert<unsigned short>();
		else
			convert<char>();
	}
#endif
};

int main(int argc, char **argv) {
	EERtoTIFF app;
	app.read(argc, argv);
	app.run();

#ifdef TIMING
	timer.printTimes(false);
#endif

	return RELION_EXIT_SUCCESS;
}
