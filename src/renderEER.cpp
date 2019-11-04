#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>

#include <src/time.h>
#include <src/metadata_table.h>
#include <src/image.h>
#include <src/renderEER.h>

//#define DEBUG_EER

//#define TIMING
#ifdef TIMING
	#define RCTIC(label) (timer.tic(label))
	#define RCTOC(label) (timer.toc(label))

	Timer timer;
	int TIMING_READ_EER = timer.setNew("read EER");
	int TIMING_BUILD_INDEX = timer.setNew("build index");
	int TIMING_UNPACK_RLE = timer.setNew("unpack RLE");
	int TIMING_RENDER_ELECTRONS = timer.setNew("render electrons");
#else
	#define RCTIC(label)
	#define RCTOC(label)
#endif

int EER_grouping = 40;
int EER_upsample = 2;

const char EERRenderer::EER_FOOTER_OK[]  = "ThermoFisherECComprOK000";
const char EERRenderer::EER_FOOTER_ERR[] = "ThermoFisherECComprERR00";
const int EERRenderer::EER_IMAGE_PIXELS = 4096 * 4096;
const unsigned int EERRenderer::EER_LEN_FOOTER = 24;

template <typename T>
void EERRenderer::render8K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons)
{
	for (int i = 0; i < n_electrons; i++)
	{
		int x = ((positions[i] & 4095) << 1) | ((symbols[i] & 2) >> 1); // 4095 = 111111111111b, 3 = 00000010b
		int y = ((positions[i] >> 12) << 1) | ((symbols[i] & 8) >> 3); //  4096 = 2^12, 8 = 00001000b
			DIRECT_A2D_ELEM(image, y, x)++;
	}
}

template <typename T>
void EERRenderer::render4K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons)
{
	for (int i = 0; i < n_electrons; i++)
	{
		int x = positions[i] & 4095; // 4095 = 111111111111b
		int y = positions[i] >> 12; //  4096 = 2^12
			DIRECT_A2D_ELEM(image, y, x)++;
	}
}

EERRenderer::EERRenderer()
{
	ready = false;
}

EERRenderer::EERRenderer(FileName fn_movie)
{
	read(fn_movie);
}

void EERRenderer::read(FileName fn_movie)
{
	/* Check file size and load everything */
	RCTIC(TIMING_READ_EER);
	FILE *fh = fopen(fn_movie.c_str(), "r");
	if (fh == NULL)
		REPORT_ERROR("Failed to open the EER file: " + fn_movie);

	fseek(fh, 0, SEEK_END);
	long long file_size = ftell(fh);
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

	/* Build frame index */
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

	ready = true;	
}

EERRenderer::~EERRenderer()
{
	if (buf != NULL)
		free(buf);
}

int EERRenderer::getNFrames()
{
	if (!ready)
		REPORT_ERROR("EERRenderer::getNFrames called before ready.");

	return frame_sizes.size();
}

int EERRenderer::getWidth()
{
	if (!ready)
		REPORT_ERROR("EERRenderer::getNFrames called before ready.");

	return 4096 * EER_upsample;
}

int EERRenderer::getHeight()
{
	if (!ready)
		REPORT_ERROR("EERRenderer::getNFrames called before ready.");

	return getWidth();
}

long long EERRenderer::renderFrames(int frame_start, int frame_end, MultidimArray<float> &image)
{
	if (!ready)
		REPORT_ERROR("EERRenderer::renderNFrames called before ready.");

	if (frame_start <= 0 || frame_start > getNFrames() ||
	    frame_end < frame_start || frame_end > getNFrames())
		REPORT_ERROR("Invalid frame range was requested.");

	// Make this 0-indexed
	frame_start--;
	frame_end--;

	long long total_n_electron = 0;

	std::vector<unsigned int> positions;
	std::vector<unsigned char> symbols;
	image.initZeros(getHeight(), getWidth());

	for (int iframe = frame_start; iframe <= frame_end; iframe++)
	{
		RCTIC(TIMING_UNPACK_RLE);
		long long pos = frame_starts[iframe];
		long long n_pix = 0, n_electron = 0;
		const int max_electrons = (frame_sizes[iframe] * 8 + 11) / 12;
		if (positions.size() < max_electrons)
		{
			positions.resize(max_electrons);
			symbols.resize(max_electrons);
		}

		// unpack every two symbols = 12 bit * 2 = 24 bit = 3 byte
		//  |bbbbBBBB|BBBBaaaa|AAAAAAAA|
		// With SIMD intrinsics at the SSSE3 level, we can unpack 10 symbols (120 bits) simultaneously.
		unsigned char p1, p2, s1, s2;

		// Because there is a footer, it is safe to go beyond the limit by two bytes.
		long long pos_limit = frame_starts[iframe] + frame_sizes[iframe];
		while (pos < pos_limit)
		{
			// symbol is bit tricky. 0000YyXx; Y and X must be flipped.
			p1 = buf[pos];
			s1 = (buf[pos + 1] & 0x0F) ^ 0x0A; // 0x0F = 00001111, 0x0A = 00001010

			p2 = (buf[pos + 1] >> 4) | (buf[pos + 2] << 4);
			s2 = (buf[pos + 2] >> 4) ^ 0x0A;

			// Note the order. Add p before checking the size and placing a new electron.
			n_pix += p1;
			if (n_pix == EER_IMAGE_PIXELS) break;
			if (p1 < 255)
			{
				positions[n_electron] = n_pix;
				symbols[n_electron] = s1;
				n_electron++;
				n_pix++;
			}

			n_pix += p2;
			if (n_pix == EER_IMAGE_PIXELS) break;
			if (p2 < 255)
			{
				positions[n_electron] = n_pix;
				symbols[n_electron] = s2;
				n_electron++;
				n_pix++;
			}

#ifdef DEBUG_EER_DETAIL
			printf("%d: %u %u, %u %u %d\n", pos, p1, s1, p2, s2, n_pix);
#endif
			pos += 3;
		}

		if (n_pix != EER_IMAGE_PIXELS)
			REPORT_ERROR("Number of pixels is not right.");

		RCTOC(TIMING_UNPACK_RLE);

		RCTIC(TIMING_RENDER_ELECTRONS);
		if (EER_upsample == 2)
			render8K(image, positions, symbols, n_electron);
		else if (EER_upsample == 1)
			render4K(image, positions, symbols, n_electron);
		else
			REPORT_ERROR("Invalid EER upsamle");
		RCTOC(TIMING_RENDER_ELECTRONS);

		total_n_electron += n_electron;
#ifdef DEBUG_EER
		printf("Decoded %lld electrons / %d pixels from frame %5d.\n", n_electron, n_pix, iframe);
#endif
	}
#ifdef DEBUG_EER
	printf("Decoded %lld electrons in total.\n", total_n_electron);
#endif

#ifdef TIMING
	timer.printTimes(false);
#endif

	return total_n_electron;
}
