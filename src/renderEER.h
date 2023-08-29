#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>

#include <src/image.h>

#include <tiffio.h>

class EERRenderer {
	private:

	// Constants
	static const char EER_FOOTER_OK[];
	static const char EER_FOOTER_ERR[];
	static const int EER_4K, EER_2K;
	static const unsigned int EER_LEN_FOOTER;
	static const uint16_t TIFF_COMPRESSION_EER8bit, TIFF_COMPRESSION_EER7bit, TIFF_COMPRESSION_EERDetailed;
	static const ttag_t TIFFTAG_EER_RLE_DEPTH, TIFFTAG_EER_SUBPIXEL_H_DEPTH, TIFFTAG_EER_SUBPIXEL_V_DEPTH;

	FileName fn_movie;

	bool ready;
	bool is_legacy; // legacy, non-TIFF container
	bool read_data;

	std::vector<long long> frame_starts, frame_sizes;
	unsigned char* buf;

	int eer_upsampling;
	int nframes, width, height;
	int preread_start, preread_end;
	uint16_t rle_bits, subpixel_bits;
	long long file_size, total_pixels;
	void readLegacy(FILE *fh);
	void lazyReadFrames();

	template <typename T>
	void render4K_to_16K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	template <typename T>
	void render4K_to_8K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	template <typename T>
	void render4K_to_4K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	template <typename T>
	void render4K_to_2K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	template <typename T>
	void render2K_to_4K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	template <typename T>
	void render2K_to_2K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	static TIFFErrorHandler prevTIFFWarningHandler;

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

	// Wrapper to the default TIFF warning handler to suppress EER private tag warnings
	static void TIFFWarningHandler(const char* module, const char* fmt, va_list ap);
	static void silenceTIFFWarnings();

	// 1-indexed
	void setFramesOfInterest(int start, int end)
	{
		if (is_legacy)
			return;

		if (read_data)
			REPORT_ERROR("Logic error in EERRenderer::setFramesOfInterest(). This must be set before rendering.");
		preread_start = start - 1;
		preread_end = end - 1;
	}

	void read(FileName _fn_movie, int eer_upsampling=1);

	// Due to a limitation in libtiff (not TIFF specification!),
	// the maximum number of frames is 65535.
	// See https://www.asmail.be/msg0055011809.html.
	int getNFrames();
	int getWidth();
	int getHeight();

	// Frame indices are 1-indexed.
	// image is cleared.
	// This function is thread-safe (except for timing).
	// It is caller's responsibility to make sure type T does not overflow.
	template <typename T>
	long long renderFrames(int frame_start, int frame_end, MultidimArray<T> &image);

	// The gain reference for EER is not multiplicative! So the inverse is taken here.
	// 0 means defect.
	// This reads the gain reference into the `gain` array but does NOT apply it to movies.
	template <typename T>
	void loadEERGain(FileName fn_gain, MultidimArray<T> &gain)
	{
		if (!ready)
			REPORT_ERROR("EERRenderer::loadEERGain called before ready.");

		const bool is_multiplicative = (fn_gain.getExtension() == "gain");
		if (is_multiplicative)
		{
			silenceTIFFWarnings();
			fn_gain += ":tif";
		}

		const int detector_width = width;
		const int detector_height = height;

		Image<T> original;
		original.read(fn_gain, true, 0, false, true); // explicitly use the first page
		const int nx_in = XSIZE(original());
		const int ny_in = YSIZE(original());
		long long size_out = getWidth();

		// Revert Y flip in TIFF reader
		if (is_multiplicative)
		{
			const int ylim = ny_in / 2;
			for (int y1 = 0; y1 < ylim; y1++)
			{
				const int y2 = ny_in - 1 - y1;
				for (int x = 0; x < nx_in; x++)
				{
					const T tmp = DIRECT_A2D_ELEM(original(), y1, x);
					DIRECT_A2D_ELEM(original(), y1, x) = DIRECT_A2D_ELEM(original(), y2, x);
					DIRECT_A2D_ELEM(original(), y2, x) = tmp;
				}
			}
		} 

		if (eer_upsampling == 2 && nx_in == detector_width && ny_in == detector_height) // (gain=det=4K, grid=8K) or (gain=det=2K, grid=4K)
		{
			gain.initZeros(size_out, size_out);
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(gain)
				DIRECT_A2D_ELEM(gain, i, j) = DIRECT_A2D_ELEM(original(), i / 2, j / 2);
		}
		else if ((eer_upsampling == 1 && nx_in == detector_width && ny_in == detector_height) || // (gain=det=grid=4K) or (gain=det=grid=2K)
		         (eer_upsampling == 2 && nx_in == size_out && ny_in == size_out)) // (det=4K, gain=grid=8K) or (det=2K, gain=grid=4K)
		{
			gain = original();
		}
		else if (eer_upsampling == 1 && nx_in == detector_width * 2 && ny_in == detector_height * 2) // (gain=8K, det=grid=4K) or (gain=4K, det=grid=2K)
		{
			gain.initZeros(size_out, size_out);
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(original())
				DIRECT_A2D_ELEM(gain, i / 2, j / 2) += DIRECT_A2D_ELEM(original(), i, j);
		}
		else
		{
			std::cerr << "Size of input gain: X = " << nx_in << " Y = " << ny_in << " Expected: X = " << size_out << " Y = " << size_out << std::endl;
			REPORT_ERROR("Unsupported gain size in EERRenderer::loadEERGain()");
		}
		
		if (!is_multiplicative)
		{
			double sum = 0;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(gain)
				sum += DIRECT_MULTIDIM_ELEM(gain, n);
			sum /= size_out * size_out;

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(gain)
			{
				if (DIRECT_MULTIDIM_ELEM(gain, n) != 0)
				{
					DIRECT_MULTIDIM_ELEM(gain, n) = sum / DIRECT_MULTIDIM_ELEM(gain, n);
				}
			}
		}
	}

	static bool isEER(FileName fn_movie)
	{
		FileName ext = fn_movie.getExtension();
		return (ext == "eer"  || ext == "ecc");
	}
};
