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
#include <src/tiff_converter.h>
#include <src/parallel.h>

// TODO: Make less verbose
//       Lossy strategy

#ifdef HAVE_TIFF

void TIFFConverter::usage()
{
	parser.writeUsage(std::cerr);
}

void TIFFConverter::read(int argc, char **argv)
{
	parser.setCommandLine(argc, argv);

	int general_section = parser.addSection("General Options");
	fn_in = parser.getOption("--i", "Input movie to be compressed (an MRC/MRCS file or a list of movies as .star or .lst)");
	fn_out = parser.getOption("--o", "Directory for output TIFF files", "./");
	fn_gain = parser.getOption("--gain", "Estimated gain map and its reliablity map (read)", "");
	nr_threads = textToInteger(parser.getOption("--j", "Number of threads (useful only for --estimate_gain)", "1"));
	only_do_unfinished = parser.checkOption("--only_do_unfinished", "Only process non-converted movies.");
	thresh_reliable = textToInteger(parser.getOption("--thresh", "Number of success needed to consider a pixel reliable", "50"));
	do_estimate = parser.checkOption("--estimate_gain", "Estimate gain");

	int tiff_section = parser.addSection("TIFF options");
	fn_compression = parser.getOption("--compression", "compression type (none, auto, deflate (= zip), lzw)", "auto");
	deflate_level = textToInteger(parser.getOption("--deflate_level", "deflate level. 1 (fast) to 9 (slowest but best compression)", "6"));
	//lossy = parser.checkOption("--lossy", "Allow slightly lossy but better compression on defect pixels");
	dont_die_on_error = parser.checkOption("--ignore_error", "Don't die on un-expected defect pixels (can be dangerous)");
	line_by_line = parser.checkOption("--line_by_line", "Use one strip per row");

	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}

template <typename T>
void TIFFConverter::write_tiff_one_page(TIFF *tif, MultidimArray<T> buf, const int filter, const int level, const float pixel_size)
{
	TIFFSetField(tif, TIFFTAG_SOFTWARE, "relion_convert_to_tiff");
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, XSIZE(buf));
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, YSIZE(buf));
	TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, line_by_line ? 1 : YSIZE(buf));
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
		TIFFWriteScanline(tif, buf.data + (ny - 1 - iy) * XSIZE(buf), iy, 0);

	TIFFWriteDirectory(tif);
}

void TIFFConverter::estimate(FileName fn_movie)
{
	Image<float> frame;
	frame.read(fn_movie, false, -1, false, true); // select_img -1, mmap false, is_2D true
	const int nframes = NSIZE(frame());

	for (int iframe = 0; iframe < nframes; iframe++)
	{
		int error = 0, changed = 0, stable = 0, negative = 0;
		
		frame.read(fn_movie, true, iframe, false, true);

		#pragma omp parallel for num_threads(nr_threads) reduction(+:error, changed, negative)
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(frame())
		{
			const float val = DIRECT_MULTIDIM_ELEM(frame(), n);
			const float gain_here = DIRECT_MULTIDIM_ELEM(gain(), n);

			if (val == 0)
			{
				continue;
			}
			else if (val < 0)
			{
//#define DEBUG
#ifdef DEBUG
				printf(" negative: %s frame %2d pos %4d %4d obs % 8.4f gain %.4f\n", 
				       fn_movie.c_str(), iframe, n / nx, n % ny, (double)val, (double)gain_here);
#endif
				negative++;
				DIRECT_MULTIDIM_ELEM(defects(), n) = -1;
			}
			else if (gain_here > val)
			{
				DIRECT_MULTIDIM_ELEM(gain(), n) = val;
				changed++;
				DIRECT_MULTIDIM_ELEM(defects(), n) = 0;
			}
			else
			{
				const int ival = (int)round(val / gain_here);
				const float expected = gain_here * ival;
				if (fabs(expected - val) > 0.0001)
				{
#ifdef DEBUG
					printf(" mismatch: %s frame %2d pos %4d %4d obs % 8.4f expected % 8.4f gain %.4f\n",
					       fn_movie.c_str(), iframe, n / nx, n % ny, (double)val,
					       (double)expected, (double)gain_here);
#endif
					error++;
					DIRECT_MULTIDIM_ELEM(defects(), n) = -1;
				}
				else if (DIRECT_MULTIDIM_ELEM(defects(), n) >= 0)
				{
					DIRECT_MULTIDIM_ELEM(defects(), n)++;
				}
			}

		}

		#pragma omp parallel for num_threads(nr_threads) reduction(+:stable)
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(defects())
		{
			short val = DIRECT_MULTIDIM_ELEM(defects(), n);
			if (val >= thresh_reliable)
				stable++;
		}

		printf(" %s Frame %03d #Changed %10d #Mismatch %10d, #Negative %10d, #Unreliable %10d / %10d\n",
		       fn_movie.c_str(), iframe + 1, changed, error, negative, ny * nx - stable, ny * nx);
	}
}

int TIFFConverter::decide_filter(int nx)
{
	if (fn_compression == "none")
		return COMPRESSION_NONE;
	else if (fn_compression == "lzw")
		return COMPRESSION_LZW;
	else if (fn_compression == "deflate" || fn_compression == "zip")
		return COMPRESSION_DEFLATE;
	else if (fn_compression == "auto")
	{
		if (nx == 4096)
			return COMPRESSION_DEFLATE; // likely Falcon
		else
			return COMPRESSION_LZW;
	}
	else
		REPORT_ERROR("Compression type must be one of none, auto, deflate (= zip) or lzw.");

	return -1;
}

template <typename T>
void TIFFConverter::unnormalise(FileName fn_movie, FileName fn_tiff)
{
	FileName fn_tmp = fn_tiff + ".tmp";
	TIFF *tif = TIFFOpen(fn_tmp.c_str(), "w");
	if (tif == NULL)
		REPORT_ERROR("Failed to open the output TIFF file: " + fn_tiff);

	Image<float> frame;
	MultidimArray<T> buf(ny, nx);
	char msg[256];

	frame.read(fn_movie, false, -1, false, true); // select_img -1, mmap false, is_2D true
	const int nframes = NSIZE(frame());
	const float angpix = frame.samplingRateX();

	for (int iframe = 0; iframe < nframes; iframe++)
	{
		int error = 0;
		
		frame.read(fn_movie, true, iframe, false, true);

		#pragma omp parallel for num_threads(nr_threads) reduction(+:error)
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(frame())
		{
			const float val = DIRECT_MULTIDIM_ELEM(frame(), n);
			const float gain_here = DIRECT_MULTIDIM_ELEM(gain(), n);
			bool is_bad = DIRECT_MULTIDIM_ELEM(defects(), n) < thresh_reliable;
			
			if (is_bad)
			{
				// TODO: implement other strategy
				DIRECT_MULTIDIM_ELEM(buf, n) = val;
				continue;
			}

			int ival = (int)round(val / gain_here);
			const float expected = gain_here * ival;
			if (fabs(expected - val) > 0.0001)
			{
				snprintf(msg, 255, " mismatch: %s frame %2d pos %4d %4d status %5d obs % 8.4f expected % 8.4f gain %.4f\n",
					 fn_movie.c_str(), iframe, n / nx, n % ny, (double)val, DIRECT_MULTIDIM_ELEM(defects(), n),
					 (double)expected, (double)gain_here);
				std::cerr << msg << std::endl;
				if (!dont_die_on_error)
					REPORT_ERROR("Unexpected pixel value in a pixel that was considered reliable");
				error++;
			}

			if (!std::is_same<T, float>::value)
			{
				const int overflow = std::is_same<T, short>::value ? 32767: 127;
				const int underflow = std::is_same<T, short>::value ? -32768: 0;

				if (ival < underflow)
				{
					ival = underflow;
					error++;
					
					printf(" underflow: %s frame %2d pos %4d %4d obs % 8.4f expected % 8.4f gain %.4f\n",
					       fn_movie.c_str(), iframe, n / nx, n % ny, (double)val,
					       (double)expected, (double)gain_here);
				}
				else if (ival > overflow)
				{
					ival = overflow;
					error++;

					printf(" overflow: %s frame %2d pos %4d %4d obs % 8.4f expected % 8.4f gain %.4f\n",
					       fn_movie.c_str(), iframe, n / nx, n % ny, (double)val,
					       (double)expected, (double)gain_here);
				}
			}
			
			DIRECT_MULTIDIM_ELEM(buf, n) = ival;
		}

		write_tiff_one_page(tif, buf, decide_filter(nx), deflate_level, angpix);
		printf(" %s Frame %3d / %3d #Error %10d\n", fn_movie.c_str(), iframe + 1, nframes, error);
	}

	TIFFClose(tif);
	std::rename(fn_tmp.c_str(), fn_tiff.c_str());
}

template <typename T>
void TIFFConverter::only_compress(FileName fn_movie, FileName fn_tiff)
{
	FileName fn_tmp = fn_tiff + ".tmp";
	TIFF *tif = TIFFOpen(fn_tmp.c_str(), "w");
	if (tif == NULL)
		REPORT_ERROR("Failed to open the output TIFF file.");

	Image<T> frame;
	frame.read(fn_movie, false, -1, false, true); // select_img -1, mmap false, is_2D true
	const int nframes = NSIZE(frame());
	const float angpix = frame.samplingRateX();
	
	for (int iframe = 0; iframe < nframes; iframe++)
	{
		frame.read(fn_movie, true, iframe, false, true);
		write_tiff_one_page(tif, frame(), decide_filter(nx), deflate_level, angpix);
		printf(" %s Frame %3d / %3d\n", fn_movie.c_str(), iframe + 1, nframes);
	}

	TIFFClose(tif);
	std::rename(fn_tmp.c_str(), fn_tiff.c_str());
}

int TIFFConverter::checkMRCtype(FileName fn_movie)
{
	// Check data type; Unfortunately I cannot do this through Image object.
	FILE *mrcin = fopen(fn_movie.c_str(), "r");
	int headers[25];
	fread(headers, sizeof(int), 24, mrcin);
	fclose(mrcin);

	return headers[3];
}

void TIFFConverter::initialise(int _rank, int _total_ranks)
{
	rank = _rank;
	total_ranks = _total_ranks;

	if (do_estimate && total_ranks != 1)
		REPORT_ERROR("MPI parallelisation is not avaialble for --estimate_gain");

	if (fn_out.back() != '/')
		fn_out += "/";

	FileName fn_first;
	FileName fn_in_ext = fn_in.getExtension();

	if (fn_in_ext == "star")
	{
		MD.read(fn_in, "movie");

		// Support non-optics group STAR files
		if (MD.numberOfObjects() == 0)
			MD.read(fn_in, "");

		if (!MD.getValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_first, 0))
			REPORT_ERROR("The input STAR file does not contain the rlnMicrographMovieName column");

		std::cout << "The number of movies in the input: " << MD.numberOfObjects() << std::endl;
	}
	else if (fn_in_ext == "lst") // treat as a simple list
	{
		std::ifstream f;
		std::string line;
		f.open(fn_in);
		while (std::getline(f, line))
		{
			MD.addObject();
			MD.setValue(EMDL_MICROGRAPH_MOVIE_NAME, line);
		}
		f.close();

		MD.getValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_first, 0);		
	}
	else 
	{
		MD.addObject();
		MD.setValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_in);
		fn_first = fn_in;
	}

	if (fn_first.getExtension() != "mrc" && fn_first.getExtension() != "mrcs")
		REPORT_ERROR(fn_first + ": the input must be MRC or MRCS files");

	if (do_estimate)
		MD.randomiseOrder();	

	// Check type and mode of the input
	Image<RFLOAT> Ihead;
	Ihead.read(fn_first, false, -1, false, true); // select_img -1, mmap false, is_2D true
	nn = NSIZE(Ihead());
	ny = YSIZE(Ihead());
	nx = XSIZE(Ihead());
	mrc_mode = checkMRCtype(fn_first);
	if (rank == 0)
		printf("Input (NX, NY, NN) = (%d, %d, %d), MODE = %d\n\n", nx, ny, nn, mrc_mode);

	if (mrc_mode != 2 && do_estimate)
		REPORT_ERROR("The input movie is not in mode 2. Gain estimation does not make sense.");

	if (fn_gain != "")
	{
		if (mrc_mode != 2)
		{
			std::cerr << "The input movie is not in mode 2. A gain reference is irrelavant." << std::endl;
		}
		else
		{
			gain.read(fn_gain + ":mrc");
			if (rank == 0)
				std::cout << "Read " << fn_gain << std::endl;
			if (XSIZE(gain()) != nx || YSIZE(gain()) != ny)
				REPORT_ERROR("The input gain has a wrong size.");

			FileName fn_defects = fn_gain.withoutExtension() + "_reliablity." + fn_gain.getExtension();
			defects.read(fn_defects + ":mrc");
			if (rank == 0)
				std::cout << "Read " << fn_defects << "\n" << std::endl;
			if (XSIZE(defects()) != nx || YSIZE(defects()) != ny)
				REPORT_ERROR("The input reliability map has a wrong size.");
		}
	}
	else if (mrc_mode == 2)
	{
		gain().reshape(ny, nx);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(gain())
			DIRECT_MULTIDIM_ELEM(gain(), n) = 999.9;
		defects().reshape(ny, nx);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(defects())
			DIRECT_MULTIDIM_ELEM(defects(), n) = -1;

		if (!do_estimate)
			std::cerr << "WARNING: To effectively compress mode 2 MRC files, you should first estimate the gain with --estimate_gain." << std::endl;
	}

	if (fn_out.contains("/"))
		system(("mkdir -p " + fn_out.beforeLastOf("/")).c_str());

	if (!do_estimate && mrc_mode == 2)
	{
		// TODO: other strategy
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(gain())
			if (DIRECT_MULTIDIM_ELEM(defects(), n) < thresh_reliable)
				DIRECT_MULTIDIM_ELEM(gain(), n) = 1.0;

		if (rank == 0 && fn_gain != "")
		{
			gain.write(fn_out + "gain-reference.mrc");
			std::cout << "Written " + fn_out + "gain-reference.mrc. Please use this file as a gain reference when processing the converted movies.\n" << std::endl; 	
		}
	}
}

void TIFFConverter::processOneMovie(FileName fn_movie, FileName fn_tiff)
{
	if (fn_movie.getExtension() != "mrc" && fn_movie.getExtension() != "mrcs")
	{
		std::cerr << fn_movie <<  " is not MRC or MRCS file. Skipped." << std::endl;
	}

	// Check type and mode of the input
	Image<RFLOAT> Ihead;
	Ihead.read(fn_movie, false, -1, false, true); // select_img -1, mmap false, is_2D true
	if (ny != YSIZE(Ihead()) || nx != XSIZE(Ihead()) || mrc_mode != checkMRCtype(fn_movie))
		REPORT_ERROR("A movie " + fn_movie + " has a different size and/or mode from other movies.");

	if (mrc_mode == 1)
	{
		only_compress<short>(fn_movie, fn_tiff);
	}
	else if (mrc_mode == 6)
	{
		only_compress<unsigned short>(fn_movie, fn_tiff);
	}
	else if (mrc_mode == 0 || mrc_mode == 101)
	{
		only_compress<char>(fn_movie, fn_tiff);
	}
	else if (do_estimate)
	{
		estimate(fn_movie);

		// Write for each movie so that one can stop anytime
		gain.write(fn_out + "gain_estimate.bin:mrc"); // .bin to prevent people from using this by mistake
		defects.write(fn_out + "gain_estimate_reliablity.bin:mrc");

		std::cout << "\nUpdated " + fn_out + "gain_estimate.bin and " + fn_out + "gain_estimate_reliablity.bin\n" << std::endl;
	}
	else
	{
		unnormalise<float>(fn_movie, fn_tiff);
	}
}

void TIFFConverter::run()
{
	long int my_first, my_last;
	 divide_equally(MD.numberOfObjects(), total_ranks, rank, my_first, my_last); // MPI parallelization

	for (long i = my_first; i <= my_last; i++)
	{
		FileName fn_movie, fn_tiff;
		MD.getValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_movie, i);

		fn_tiff = fn_out + fn_movie.withoutExtension() + ".tif";
		if (only_do_unfinished && !do_estimate && exists(fn_tiff))
		{			
			std::cout << "Skipping already processed " << fn_movie << std::endl;
			continue;
		}

		std::cout << "Processing " << fn_movie;
		if (!do_estimate)
			std::cout  << " into " << fn_tiff;
		std::cout << std::endl;

		if (fn_tiff.contains("/"))
			system(("mkdir -p " + fn_tiff.beforeLastOf("/")).c_str());

		processOneMovie(fn_movie, fn_tiff);
	}
}
#endif
