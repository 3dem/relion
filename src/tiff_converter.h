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

#include <tiffio.h>

class TIFFConverter 
{
public:
	void usage();
	void read(int argc, char **argv);
	void initialise(int _rank=0, int _total_ranks=1);
	void run();

	template <typename T>
	static void write_tiff_one_page(TIFF *tif, MultidimArray<T> buf, const float pixel_size=-1, const int filter=COMPRESSION_LZW, const int level=6, const bool strip_per_line=false)
	{
		TIFFSetField(tif, TIFFTAG_SOFTWARE, "RELION");
		TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, XSIZE(buf));
		TIFFSetField(tif, TIFFTAG_IMAGELENGTH, YSIZE(buf));
		TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, strip_per_line ? 1 : YSIZE(buf));
		TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
		TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);		

		if (std::is_same<T, float>::value)
		{
			TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
			TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
		}
		else if (std::is_same<T, unsigned short>::value)
		{
			TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
			TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
		}
		else if (std::is_same<T, short>::value)
		{
			TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
			TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
		}
		else if (std::is_same<T, char>::value || std::is_same<T, unsigned char>::value )
		{
			TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
			TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
		}
		else if (std::is_same<T, signed char>::value)
		{
			TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
			TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT);
		}
		else
		{
			REPORT_ERROR("write_tiff_one_page: unknown data type");
		}

		// compression is COMPRESSION_LZW or COMPRESSION_DEFLATE or COMPRESSION_NONE
		TIFFSetField(tif, TIFFTAG_COMPRESSION, filter);
		if (filter == COMPRESSION_DEFLATE)
		{
			if (level <= 0 || level > 9)
				REPORT_ERROR("Deflate level must be 1, 2, ..., 9");
			TIFFSetField(tif, TIFFTAG_ZIPQUALITY, level);
		}

		if (pixel_size > 0)
		{
			TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_CENTIMETER); // 1 cm = 1E8 A
			TIFFSetField(tif, TIFFTAG_XRESOLUTION, 1E8 / pixel_size); // pixels / 1 cm
			TIFFSetField(tif, TIFFTAG_YRESOLUTION, 1E8 / pixel_size);
		}

		// Have to flip the Y axis
		for (int iy = 0; iy < YSIZE(buf); iy++)
			TIFFWriteScanline(tif, buf.data + (YSIZE(buf) - 1 - iy) * XSIZE(buf), iy, 0);

		TIFFWriteDirectory(tif);
	}

private:
	int rank, total_ranks;

	FileName fn_in, fn_out, fn_gain, fn_compression;
	bool do_estimate, input_type, lossy, dont_die_on_error, line_by_line, only_do_unfinished, eer_short;
	int deflate_level, thresh_reliable, nr_threads, eer_upsampling, eer_grouping;
	IOParser parser;

	MetaDataTable MD;
	Image<short> defects;
	Image<float> gain;
	int mrc_mode; // -99 for EER

	void estimate(FileName fn_movie);
	int decide_filter(int nx, bool isEER=false);

	template <typename T>
	void unnormalise(FileName fn_movie, FileName fn_tiff);

	template <typename T>
	void only_compress(FileName fn_movie, FileName fn_tiff);
	int checkMRCtype(FileName fn_movie);
	void processOneMovie(FileName fn_movie, FileName fn_tiff);
};
