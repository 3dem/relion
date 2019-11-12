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
	FileName fn_in, fn_out, fn_gain, fn_defects, fn_gain_out, fn_defects_out;
	bool do_estimate, input_type, lossy, dont_die_on_error;
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

		lossy = parser.checkOption("--lossy", "Allow slightly lossy but better compression on defect pixels");
		dont_die_on_error = parser.checkOption("--ignore_error", "Don't die on un-expected defect pixels (can be dangerous)");

		do_estimate = parser.checkOption("--estimate_gain", "Estimate gain");
		fn_gain_out = parser.getOption("--gain_out", "Estimated gain reference (written)", "gain.mrc");
		fn_defects_out = parser.getOption("--defect_out", "Estimated unreliable pixels (written)", "defects.mrc");

		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	template <typename T>
	void write_tiff_one_page(TIFF *tif, MultidimArray<T> buf, const int compression=COMPRESSION_LZW)
	{
		TIFFSetField(tif, TIFFTAG_SOFTWARE, "relion_convert_to_tiff");
		TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, nx);
		TIFFSetField(tif, TIFFTAG_IMAGELENGTH, ny);
		TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, ny);
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
		else if (std::is_same<T, char>::value)
		{
			TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
			TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
		}
		else
		{
			REPORT_ERROR("write_tiff_one_page: unknown data type");
		}

		// compression is COMPRESSION_LZW or COMPRESSION_DEFLATE or COMPRESSION_NONE
		TIFFSetField(tif, TIFFTAG_COMPRESSION, compression);

		// Have to flip the Y axis
		for (int iy = 0; iy < ny; iy++)
			TIFFWriteScanline(tif, buf.data + (ny - 1 - iy) * nx, iy, 0);

		TIFFWriteDirectory(tif);
	}

	void estimate(FileName fn_movie)
	{
		Image<float> frame;

		for (int iframe = 0; iframe < nn; iframe++)
		{
			int error = 0, changed = 0, stable = 0, negative = 0;
			
			frame.read(fn_movie, true, iframe);

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

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(defects())
			{
				short val = DIRECT_MULTIDIM_ELEM(defects(), n);
				if (val > 10)
					stable++;
			}

			printf("Frame %02d #Changed %10d #Mismatch %10d, #Negative %10d, #Unstable %10d / %10d\n",
			       iframe, changed, error, negative, ny * nx - stable, ny * nx);
		}
	}

	template <typename T>
	void unnormalise()
	{
		TIFF *tif = TIFFOpen(fn_out.c_str(), "w");
		if (tif == NULL)
			REPORT_ERROR("Failed to open the output TIFF file.");

		Image<float> frame;
		MultidimArray<T> buf(ny, nx);

		for (int iframe = 0; iframe < nn; iframe++)
		{
			int error = 0;
			
			frame.read(fn_in, true, iframe);

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(frame())
			{
				const float val = DIRECT_MULTIDIM_ELEM(frame(), n);
				const float gain_here = DIRECT_MULTIDIM_ELEM(gain(), n);
				bool is_bad = DIRECT_MULTIDIM_ELEM(defects(), n);
				
				const int ival = (int)round(val / gain_here);
				const float expected = gain_here * ival;
				if (fabs(expected - val) > 0.0001)
				{
					printf(" mismatch: frame %2d pos %4d %4d obs % 8.4f expected % 8.4f gain %.4f\n",
					       iframe, n / nx, n % ny, (double)val,
					       (double)expected, (double)gain_here);
					is_bad = true;
					error++;
				}

				// TODO: do magic
			}

			// TODO: Write one page
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
			frame.read(fn_in, true, iframe);
			write_tiff_one_page(tif, frame, compression=COMPRESSION_DEFLATE);
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

		if (input_type == 2)
		{
			only_compress<short>();
		}
		else if (input_type == 101)
		{
			only_compress<char>();
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
			convert<float>();
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
