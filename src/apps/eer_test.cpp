#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>

#include <src/time.h>
#include <src/metadata_table.h>
#include <src/args.h>
#include <src/image.h>

#define TIMING
#ifdef TIMING
	#define RCTIC(label) (timer.tic(label))
	#define RCTOC(label) (timer.toc(label))

	Timer timer;
	int TIMING_READ_EER = timer.setNew("read EER");
	int TIMING_BUILD_INDEX = timer.setNew("build index");
	int TIMING_UNPACK_RLE = timer.setNew("unpack RLE");
#else
	#define RCTIC(label)
	#define RCTOC(label)
#endif

const char EER_FOOTER_OK[]  = "ThermoFisherECComprOK000";
const char EER_FOOTER_ERR[] = "ThermoFisherECComprERR00";
const int EER_IMAGE_PIXELS = 4096 * 4096;
const unsigned int EER_LEN_FOOTER = 24;

class eer_test {
	public:

	IOParser parser;

	FileName fn_movie, fn_out;
	int n_threads, max_frames;

	std::vector<long long> frame_starts, frame_sizes;
	unsigned char* buf;

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
		long file_size;

		RCTIC(TIMING_READ_EER);
		FILE *fh = fopen(fn_movie.c_str(), "r");
		if (fh == NULL)
			REPORT_ERROR("Failed to open the EER file.");

		fseek(fh, 0, SEEK_END);
		file_size = ftell(fh);
		fseek(fh, 0, SEEK_SET);
#ifdef DEBUG_EER
		printf("File size: %ld\n", file_size);
#endif
		buf = (unsigned char*)malloc(file_size);
		if (buf == NULL)
			REPORT_ERROR("Failed to allocate the buffer.");
		fread(buf, sizeof(char), file_size, fh);
		fclose(fh);
		RCTOC(TIMING_READ_EER);

		RCTIC(TIMING_BUILD_INDEX);
		long long pos = file_size;
		while (pos > 0)
		{
			pos -= EER_LEN_FOOTER;
			if (strncmp(EER_FOOTER_OK, (char*)buf + pos, EER_LEN_FOOTER) == 0)
			{
				pos -= 8;
				long long frame_size = *(long long*)(buf + pos) * sizeof(long long);
				pos -= frame_size;
				frame_starts.push_back(pos);
				frame_sizes.push_back(frame_size);
#ifdef DEBUG_EER
				printf("Frame: LAST-%5d Start at: %08d Frame size: %lld\n", frame_starts.size(), pos, frame_size);
#endif
			}
			else // if (strncmp(EER_FOOTER_ERR, (char*)buf + pos, EER_LEN_FOOTER) == 0)
			{
				REPORT_ERROR("Broken frame");
			}
		}
		std::reverse(frame_starts.begin(), frame_starts.end());
		std::reverse(frame_sizes.begin(), frame_sizes.end());
		RCTOC(TIMING_BUILD_INDEX);

		const long long nframes = frame_sizes.size();
		long long total_n_electron = 0;

		RCTIC(TIMING_UNPACK_RLE);
		for (int iframe = 0; iframe < nframes; iframe++)
		{
			pos = frame_starts[iframe];
			long long n_pix = 0, n_electron = 0;

			// unpack every two symbols = 12 bit * 2 = 24 bit = 3 byte
			// |bbbbBBBB|BBBBaaaa|AAAAAAAA|
			unsigned char p1, p2, s1, s2;

			// Because there is a footer, it is safe to go beyond the limit by two bytes.
			long long pos_limit = frame_starts[iframe] + frame_sizes[iframe];
			while (pos < pos_limit)
			{
				p1 = buf[pos];
				s1 = buf[pos + 1] & 0x0F; // 00001111

				p2 = (buf[pos + 1] >> 4) | (buf[pos + 2] << 4);
				s2 = buf[pos + 2] >> 4;

				// Note the order. Add p before checking the size and placing a new electron.
				n_pix += p1;
				if (n_pix == EER_IMAGE_PIXELS) break;
				if (p1 < 255)
				{
					n_electron++;
					n_pix++;
				}

				n_pix += p2;
				if (n_pix == EER_IMAGE_PIXELS) break;
				if (p2 < 255)
				{
					n_electron++;
					n_pix++;
				}

#ifdef DEBUG_EER
				printf("%d: %u %u, %u %u %d\n", pos, p1, s1, p2, s2, n_pix);
#endif
				pos += 3;
			}
			if (n_pix != EER_IMAGE_PIXELS)
				REPORT_ERROR("Number of pixels is not right.");

			total_n_electron += n_electron;
#ifdef DEBUG_EER
			printf("Decoded %lld electrons / %d pixels from frame %5d.\n", n_electron, n_pix, iframe);
#endif
		}
		RCTOC(TIMING_UNPACK_RLE);
		printf("Decoded %lld electrons in total.\n", total_n_electron);

		if (buf != NULL)
			free(buf);
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
