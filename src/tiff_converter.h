/***************************************************************************
 *
 * Author: "Takanori Nakane"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include <cstdio>
#include <cmath>
#include <src/args.h>
#include <src/image.h>
#include <src/metadata_table.h>

#ifdef HAVE_TIFF
#include <tiffio.h>

class TIFFConverter 
{
public:
	void usage();
	void read(int argc, char **argv);
	void initialise(int _rank=0, int _total_ranks=1);
	void run();

private:
	int rank, total_ranks;

	FileName fn_in, fn_out, fn_gain, fn_compression;
	bool do_estimate, input_type, lossy, dont_die_on_error, line_by_line, only_do_unfinished;
	int deflate_level, thresh_reliable, nr_threads;
	IOParser parser;

	MetaDataTable MD;
	Image<short> defects;
	Image<float> gain;
	int nn, ny, nx, mrc_mode;

	template <typename T>
	void write_tiff_one_page(TIFF *tif, MultidimArray<T> buf, const int filter=COMPRESSION_LZW, const int level=6, const float pixel_size=-1);
	void estimate(FileName fn_movie);
	int decide_filter(int nx);

	template <typename T>
	void unnormalise(FileName fn_movie, FileName fn_tiff);

	template <typename T>
	void only_compress(FileName fn_movie, FileName fn_tiff);
	int checkMRCtype(FileName fn_movie);
	void processOneMovie(FileName fn_movie, FileName fn_tiff);
};
#endif
