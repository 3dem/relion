#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>

#include <src/image.h>

#ifdef HAVE_TIFF
#include <tiffio.h>
#endif

extern int EER_grouping; // TODO: avoid global variable!
extern int EER_upsample; //       design better way; MicrographModel also needs to change
                         //       allow non-uniform slicing?

class EERRenderer {
	private:

	FileName fn_movie;
#ifdef HAVE_TIFF
	TIFF *ftiff;
#endif

	bool ready;
	bool is_legacy;
	bool read_data;

	std::vector<long long> frame_starts, frame_sizes;
	unsigned char* buf;

	static const char EER_FOOTER_OK[];
	static const char EER_FOOTER_ERR[];
	static const int EER_IMAGE_WIDTH, EER_IMAGE_HEIGHT, EER_IMAGE_PIXELS;
	static const unsigned int EER_LEN_FOOTER;
	static const uint16_t TIFF_COMPRESSION_EER;

	int nframes;
	long long file_size;
	void readLegacy(FILE *fh);
	void lazyReadFrames();

	template <typename T>
	void render8K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	template <typename T>
	void render4K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	public:

	EERRenderer();
	~EERRenderer();

	//TODO: Implement proper copy constructors. Currently, they are disabled to prevent memory corruption.
	EERRenderer(const EERRenderer&)
	{
		REPORT_ERROR("Copy constructor for EERRenderer not implemented yet.");
	}

	EERRenderer& operator=(const EERRenderer&)
	{
		REPORT_ERROR("Copy assignment operator for EERRenderer not implemented yet.");
	}

	void read(FileName _fn_movie);

	int getNFrames();
	int getWidth();
	int getHeight();

	// Frame indices are 1-indexed.
	// image is cleared.
	// This function is thread-safe (except for timing).
	// It is caller's responsibility to make sure type T does not overflow.
	template <typename T>
	long long renderFrames(int frame_start, int frame_end, MultidimArray<T> &image);

	template <typename T>
	static void upsampleEERGain(MultidimArray<T> &gain)
	{
		const long long size = EER_IMAGE_WIDTH * EER_upsample;
		RFLOAT sum = 0;

		if (EER_upsample == 2 && XSIZE(gain) == EER_IMAGE_WIDTH && YSIZE(gain) == EER_IMAGE_HEIGHT) // gain = 4K and grid = 8K
		{
			MultidimArray<T> original = gain;

			gain.resize(size, size);
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(gain)
			{
				DIRECT_A2D_ELEM(gain, i, j) = DIRECT_A2D_ELEM(original, i / 2, j / 2);
				sum += DIRECT_A2D_ELEM(gain, i, j);
			}
		}
		else if ((EER_upsample == 1 && XSIZE(gain) == EER_IMAGE_WIDTH && YSIZE(gain) == EER_IMAGE_HEIGHT) || // gain = 4K and grid = 4K
		         (EER_upsample == 2 && XSIZE(gain) == EER_IMAGE_WIDTH * 2 && YSIZE(gain) == EER_IMAGE_HEIGHT * 2)) // gain = 8K and grid = 8K
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(gain)
			{
				sum += DIRECT_MULTIDIM_ELEM(gain, n);
			}
		}	
	
		else
		{
			std::cerr << "Size of input gain: X = " << XSIZE(gain) << " Y = " << YSIZE(gain) << " Expected: X = " << (EER_IMAGE_WIDTH * 2) << " Y = " << (EER_IMAGE_HEIGHT * 2) << std::endl;
			REPORT_ERROR("Invalid gain size in EERRenderer::upsampleEERGain()");
		}

		sum /= size * size;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(gain)
		{
			if (DIRECT_MULTIDIM_ELEM(gain, n) != 0)
			{
				DIRECT_MULTIDIM_ELEM(gain, n) = sum / DIRECT_MULTIDIM_ELEM(gain, n);
			}
		}
	}

	static bool isEER(FileName fn_movie)
	{
		FileName ext = fn_movie.getExtension();
		return (ext == "eer"  || ext == "ecc");
	}
};
