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

// TODO: MPI parallelization
// TODO: STAR input

#ifndef HAVE_TIFF
int main(int argc, char *argv[])
{
	REPORT_ERROR("To use this program, please recompile with libtiff.");
}
#else

#include <tiffio.h>

class convert_to_tiff
{
public:
	FileName fn_in, fn_out, fn_gain, fn_defects, fn_gain_out, fn_defects_out, fn_compression;
	bool do_estimate, input_type, lossy, dont_die_on_error, line_by_line;
	int deflate_level, thresh_reliable, nr_threads;
	IOParser parser;

	Image<short> defects;
	Image<float> gain;
	int nn, ny, nx;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("Options");
		fn_in = parser.getOption("--i", "Input movie to be compressed");
		fn_out = parser.getOption("--o", "Rootname for output projections", "TIFF");
		fn_gain = parser.getOption("--gain", "Estimated gain reference to read (read)", "");
		fn_defects = parser.getOption("--defect", "Estimated unreliable pixels (read)", "");

		fn_compression = parser.getOption("--compression", "compression type (none, auto, deflate (= zip), lzw)", "auto");
		deflate_level = textToInteger(parser.getOption("--deflate_level", "deflate level. 1 (fast) to 9 (slowest but best compression)", "6"));
		lossy = parser.checkOption("--lossy", "Allow slightly lossy but better compression on defect pixels");
		dont_die_on_error = parser.checkOption("--ignore_error", "Don't die on un-expected defect pixels (can be dangerous)");
		line_by_line = parser.checkOption("--line_by_line", "Use one strip per row");

		do_estimate = parser.checkOption("--estimate_gain", "Estimate gain");
		fn_gain_out = parser.getOption("--gain_out", "Estimated gain reference (written)", "gain.mrc");
		fn_defects_out = parser.getOption("--defect_out", "Estimated unreliable pixels (written)", "defects.mrc");
		thresh_reliable = textToInteger(parser.getOption("--thresh", "Number of success needed to consider a pixel reliable", "10"));
		nr_threads = textToInteger(parser.getOption("--j", "Number of threads (More than 2 is not effective)", "1"));
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	template <typename T>
	void write_tiff_one_page(TIFF *tif, MultidimArray<T> buf, const int filter=COMPRESSION_LZW, const int level=6)
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
		else if (std::is_same<T, short>::value)
		{
			TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
			TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
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

		// Have to flip the Y axis
		for (int iy = 0; iy < YSIZE(buf); iy++)
			TIFFWriteScanline(tif, buf.data + (ny - 1 - iy) * XSIZE(buf), iy, 0);

		TIFFWriteDirectory(tif);
	}

	void estimate(FileName fn_movie)
	{
		Image<float> frame;

		for (int iframe = 0; iframe < nn; iframe++)
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
					printf(" negative: frame %2d pos %4d %4d obs % 8.4f gain %.4f\n", 
					       iframe, n / nx, n % ny, (double)val, (double)gain_here);
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
						printf(" mismatch: frame %2d pos %4d %4d obs % 8.4f expected % 8.4f gain %.4f\n",
						       iframe, n / nx, n % ny, (double)val,
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

			printf("Frame %02d #Changed %10d #Mismatch %10d, #Negative %10d, #Unreliable %10d / %10d\n",
			       iframe, changed, error, negative, ny * nx - stable, ny * nx);
		}
	}

	int decide_filter(int nx)
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
	void unnormalise()
	{
		TIFF *tif = TIFFOpen(fn_out.c_str(), "w");
		if (tif == NULL)
			REPORT_ERROR("Failed to open the output TIFF file.");

		Image<float> frame;
		MultidimArray<T> buf(ny, nx);
		char msg[256];
	
		for (int iframe = 0; iframe < nn; iframe++)
		{
			int error = 0;
			
			frame.read(fn_in, true, iframe, false, true);

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
					snprintf(msg, 255, " mismatch: frame %2d pos %4d %4d status %5d obs % 8.4f expected % 8.4f gain %.4f\n",
					       iframe, n / nx, n % ny, (double)val, DIRECT_MULTIDIM_ELEM(defects(), n),
					       (double)expected, (double)gain_here);
					std::cerr << "mismatch" << msg << std::endl;
					if (!dont_die_on_error)
						REPORT_ERROR("Unexpected pixel value in a pixel that was considered reliable");
					error++;
				}

				if (false) // TODO: for integer output
				{
					if (ival < 0)
					{
						ival = 0;
						error++;
						
						printf(" negative: frame %2d pos %4d %4d obs % 8.4f expected % 8.4f gain %.4f\n",
					               iframe, n / nx, n % ny, (double)val,
					               (double)expected, (double)gain_here);
					}
					else if (ival > 127) // TOOD: Use proper limit
					{
						ival = 127;
						error++;

						printf(" overflow: frame %2d pos %4d %4d obs % 8.4f expected % 8.4f gain %.4f\n",
					               iframe, n / nx, n % ny, (double)val,
					               (double)expected, (double)gain_here);
					}
				}
				
				DIRECT_MULTIDIM_ELEM(buf, n) = ival;
			}

			write_tiff_one_page(tif, buf, decide_filter(nx), deflate_level);
			printf("Frame %3d / %3d\n", iframe + 1, nn);
		}

		TIFFClose(tif);
	}

	template <typename T>
	void only_compress()
	{
		TIFF *tif = TIFFOpen(fn_out.c_str(), "w");
		if (tif == NULL)
			REPORT_ERROR("Failed to open the output TIFF file.");

		Image<T> frame;
		for (int iframe = 0; iframe < nn; iframe++)
		{
			frame.read(fn_in, true, iframe, false, true);
			write_tiff_one_page(tif, frame(), decide_filter(nx), deflate_level);
			printf("Frame %3d / %3d\n", iframe + 1, nn);
		}

		TIFFClose(tif);
	}

	void run()
	{
		Image<RFLOAT> Ihead;
		Ihead.read(fn_in, false, -1, false, true); // select_img -1, mmap false, is_2D true
		nn = NSIZE(Ihead());
		ny = YSIZE(Ihead());
		nx = XSIZE(Ihead());

		// Check data type; Unfortunately I cannot do this through Image object.
		FILE *mrcin = fopen(fn_in.c_str(), "r");
		int headers[25];
		fread(headers, sizeof(int), 24, mrcin);
		int input_type = headers[3];
		fclose(mrcin);

		printf("Input (NX, NY, NN) = (%d, %d, %d), MODE = %d\n", nx, ny, nn, input_type);

		if (input_type == 1 || input_type == 6)
		{
			only_compress<short>();
			return;
		}
		else if (input_type == 0 || input_type == 101)
		{
			only_compress<char>();
			return;
		}
		else if (input_type != 2)
		{
			REPORT_ERROR("Input MRC file must be in mode 0 (BYTE), 1 (SHORT), 2 (FLOAT) or 101 (4-bit)"); 
		}

		if (fn_gain != "")
		{
			gain.read(fn_gain);
			std::cout << "Read " << fn_gain << std::endl;
			if (XSIZE(gain()) != nx || YSIZE(gain()) != ny)
				REPORT_ERROR("The input gain has a wrong size.");
		}
		else
		{
			gain().reshape(ny, nx);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(gain())
				DIRECT_MULTIDIM_ELEM(gain(), n) = 999.9;
		}

		if (fn_defects != "")
		{
			defects.read(fn_defects);
			std::cout << "Read " << fn_defects << std::endl;
			if (XSIZE(defects()) != nx || YSIZE(defects()) != ny)
				REPORT_ERROR("The input defect map has a wrong size.");
		}
		else
		{
			defects().reshape(ny, nx);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(defects())
				DIRECT_MULTIDIM_ELEM(defects(), n) = -1;
		}

		if (do_estimate)
		{
			estimate(fn_in);

			defects.write(fn_defects_out);
			gain.write(fn_gain_out);
		}
		else
		{
			// TODO: other strategy
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(gain())
				if (DIRECT_MULTIDIM_ELEM(defects(), n) < thresh_reliable)
					DIRECT_MULTIDIM_ELEM(gain(), n) = 1.0;

			gain.write(fn_gain_out);

			unnormalise<float>();
		}
	}
};

int main(int argc, char *argv[])
{
	convert_to_tiff app;

	try
	{
		app.read(argc, argv);
		app.run();
	}
	catch (RelionError XE)
	{
        	std::cerr << XE;
	        return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}
#endif
