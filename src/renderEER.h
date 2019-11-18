#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>

#include <src/image.h>

extern int EER_grouping; // TODO: TAKANORI: Make this a command-line argument
extern int EER_upsample;

class EERRenderer {
	private:

	FileName fn_movie;
	int n_threads;

	bool ready;

	std::vector<long long> frame_starts, frame_sizes;
	unsigned char* buf;

	static const char EER_FOOTER_OK[];
	static const char EER_FOOTER_ERR[];
	static const int EER_IMAGE_PIXELS;
	static const unsigned int EER_LEN_FOOTER;

	template <typename T>
	void render8K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	template <typename T>
	void render4K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	public:

	EERRenderer();
	EERRenderer(FileName fn_movie);
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

	void read(FileName fn_movie);
	
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
		if (EER_upsample == 2)
		{
			const long long size = 4096 * EER_upsample;
			MultidimArray<T> original = gain;

			gain.resize(size, size);
			RFLOAT sum = 0;
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(gain)
			{
				DIRECT_A2D_ELEM(gain, i, j) = DIRECT_A2D_ELEM(original, i / 2, j / 2);
				sum += DIRECT_A2D_ELEM(gain, i, j);
			}
			sum /= size * size;

			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(gain)
			{
				if (DIRECT_A2D_ELEM(gain, i, j) != 0)
				{
					DIRECT_A2D_ELEM(gain, i, j) = sum / DIRECT_A2D_ELEM(gain, i, j);
				}
			}
		}
		else if (EER_upsample == 1)
			return;
		else
			REPORT_ERROR("Invalid EER_upsample");
	}
};
