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

const char EER_FOOTER_OK[]  = "ThermoFisherECComprOK000";
const char EER_FOOTER_ERR[] = "ThermoFisherECComprERR00";
const unsigned int EER_LEN_FOOTER = 24;

class eer_test {
	public:

	IOParser parser;

	FileName fn_movie, fn_out;
	int n_threads, max_frames;

	std::vector<long long> frame_start;
	char* buf;

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

		FILE *fh = fopen(fn_movie.c_str(), "r");
		if (fh == NULL)
			REPORT_ERROR("Failed to open the EER file.");

		fseek(fh, 0, SEEK_END);
		file_size = ftell(fh);
		fseek(fh, 0, SEEK_SET);
#ifdef DEBUG_EER
		printf("File size: %ld\n", file_size);
#endif
		buf = (char*)malloc(file_size);
		if (buf == NULL)
			REPORT_ERROR("Failed to allocate the buffer.");
		fread(buf, sizeof(char), file_size, fh);
		fclose(fh);

		long long pos = file_size;
		while (pos > 0)
		{
			pos -= EER_LEN_FOOTER;
			if (strncmp(EER_FOOTER_OK, buf + pos, EER_LEN_FOOTER) == 0)
			{
				pos -= 8;
				long long frame_size = *(long long*)(buf + pos) * sizeof(long long);
				pos -= frame_size;
				frame_start.push_back(pos);
#ifdef DEBUG_EER
				printf("Frame: LAST-%5d Start at: %08d Frame size: %lld\n", frame_start.size(), pos, frame_size);
#endif
			}
			else // if (strncmp(EER_FOOTER_ERR, buf + pos, EER_LEN_FOOTER) == 0)
			{
				REPORT_ERROR("Broken frame");
			}
		}
		std::reverse(frame_start.begin(), frame_start.end());

		if (buf != NULL)
			free(buf);
	}
};

int main(int argc, char **argv) {
	eer_test app;
	app.read(argc, argv);
	app.run();

	return RELION_EXIT_SUCCESS;
}
