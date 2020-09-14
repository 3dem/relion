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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include <iostream>
#include <omp.h>

#include <src/time.h>
#include <src/metadata_table.h>
#include <src/args.h>
#include <src/image.h>

#include <src/jaz/single_particle/obs_model.h>
#include <src/renderEER.h>


/* A simple program to estimate the gain reference by averaging all unaligned frames and taking the inverse.
 * This program should be used only as a last resort when the gain reference was lost or badly calibrated.
 */

class estimate_gain
{
	public:

	IOParser parser;

	FileName fn_movie_star, fn_out;
	int n_threads, max_frames, eer_upsampling;
	bool randomise_order, dont_invert;

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("Options");
		fn_movie_star = parser.getOption("--i", "Input movie STAR file");
		fn_out = parser.getOption("--o", "Output file name");
		n_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));
		max_frames = textToInteger(parser.getOption("--max_frames", "Target number of frames to average (rounded to movies; -1 means use all)", "-1"));
		randomise_order = parser.checkOption("--random", "Randomise the order of input movies before taking subset");
		dont_invert = parser.checkOption("--dont_invert", "Don't take the inverse but simply writes the sum");		
		eer_upsampling = textToInteger(parser.getOption("--eer_upsampling", "EER upsampling (1 = 4K or 2 = 8K)", "2"));
		// --eer_upsampling 3 is only for debugging. Hidden.
		if (eer_upsampling != 1 && eer_upsampling != 2 && eer_upsampling != 3)
			REPORT_ERROR("eer_upsampling must be 1, 2 or 3");

		if (parser.checkForErrors())
		{
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
		}
	}

	void run()
	{
		MetaDataTable MDin;

		MDin.read(fn_movie_star, "movies");
		// Support non-optics group STAR files
		if (MDin.numberOfObjects() == 0)
			MDin.read(fn_movie_star, "");

		const int n_total_movies = MDin.numberOfObjects();
		if (n_total_movies == 0)
			REPORT_ERROR("No movies in the input STAR file");
		else
			std::cout << "Number of movies in the STAR file: " << n_total_movies << std::endl;
		
		FileName fn_img;
		if (!MDin.getValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_img, 0))
			REPORT_ERROR("The input STAR file does not contain the rlnMicrographMovieName column.");

		Image<RFLOAT> Ihead;
		std::vector<Image<RFLOAT>> Isums(n_threads);
		int ny = 0, nx = 0;

		if (randomise_order)
		{
			MDin.randomiseOrder();
			std::cout << "Randomised the order of input movies." << std::endl;
		}

		MDin.getValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_img, 0);
		if (EERRenderer::isEER(fn_img))
		{
			EERRenderer renderer;
			renderer.read(fn_img, eer_upsampling);
			ny = renderer.getHeight();
			nx = renderer.getWidth();

			if (!dont_invert)
				REPORT_ERROR("The input movie is EER. For EER, the gain reference is expected to be the average of counts, not the inverse as in K2/K3. Thus, you need the --dont_invert flag.");
		}
		else
		{
			Ihead.read(fn_img, false, -1, false, true); // select_img -1, mmap false, is_2D true
			ny = YSIZE(Ihead());
			nx = XSIZE(Ihead());
		}

		std::cout << "The size of the input: NY = " << ny << " NX = " << nx << std::endl;
		
		for (int i = 0; i < n_threads; i++)
		{
			Isums[i]().resize(ny, nx);
			Isums[i]().initZeros(ny, nx);
		}

		std::cout << "Summing frames ... " << std::endl;
		long n_frames_used = 0, n_movies_done = 0;

		if (max_frames > 0)
			init_progress_bar(max_frames);
		else
			init_progress_bar(MDin.numberOfObjects());

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
		{
			MDin.getValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_img);

			int n_frames;
			if (!EERRenderer::isEER(fn_img))
			{
				Ihead.read(fn_img, false, -1, false, true);
				n_frames = NSIZE(Ihead());
				
				if (ny != YSIZE(Ihead()) || nx != XSIZE(Ihead()))
				{
					std::cerr << "The size of the movie " + fn_img + " does not much the size of the others. Skipped." << std::endl;
					continue;
				}

				#pragma omp parallel for num_threads(n_threads)
				for (int iframe = 0; iframe < n_frames; iframe++)
				{
					Image<RFLOAT> Iframe;
					Iframe.read(fn_img, true, iframe, false, true); // mmap false, is_2D true
					const int tid = omp_get_thread_num();
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iframe())
					{
						DIRECT_MULTIDIM_ELEM(Isums[tid](), n) += DIRECT_MULTIDIM_ELEM(Iframe(), n);
					}
				}
			}
			else
			{
				EERRenderer renderer;
				renderer.read(fn_img, eer_upsampling);
				n_frames = renderer.getNFrames();
				
				if (ny != renderer.getHeight() || nx != renderer.getWidth())
				{
					std::cerr << "The size of the movie " + fn_img + " does not much the size of the others. Skipped." << std::endl;
					continue;
				}
			
				const int eer_grouping = (n_frames + n_threads - 1) / n_threads;
				#pragma omp parallel for num_threads(n_threads)
				for (int frame = 1; frame < n_frames; frame += eer_grouping)
				{
					int frame_end = frame + eer_grouping - 1;
					if (frame_end > n_frames)
						frame_end = n_frames;	

					const int tid = omp_get_thread_num();
//					std::cout << " Thread " << tid << ": Rendering EER (hardware) frame " << frame << " to " << frame_end << std::endl;
					MultidimArray<unsigned short> buf;
					// unfortunately this function clears the buffer
					renderer.renderFrames(frame, frame_end, buf); 
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(buf)
					{
						DIRECT_MULTIDIM_ELEM(Isums[tid](), n) += DIRECT_MULTIDIM_ELEM(buf, n);
					}
				}
			}

			n_frames_used += n_frames;
			n_movies_done++;

			if (max_frames > 0 && n_frames_used > max_frames)
				break;

			if (max_frames > 0)
				progress_bar(n_frames_used);
			else
				progress_bar(n_movies_done);
		}
		if (max_frames > 0)
			progress_bar(max_frames);
		else
			progress_bar(MDin.numberOfObjects());

		for (int i = 1; i < n_threads; i++)
		{
			#pragma omp parallel for num_threads(n_threads)
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Isums[0]())
			{
				DIRECT_MULTIDIM_ELEM(Isums[0](), n) += DIRECT_MULTIDIM_ELEM(Isums[i](), n);
			}
		}

		std::cout << "Summed " << n_frames_used << " frames from " << n_movies_done << " movies." << std::endl;

		if (dont_invert)
		{
			std::cout << "Because of --dont_invert, the sum is written as is." << std::endl;
		}
		else
		{
			double total_count = 0;
			#pragma omp parallel for num_threads(n_threads) reduction(+: total_count)
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Isums[0]())
				total_count += DIRECT_MULTIDIM_ELEM(Isums[0](), n);
			const double avg_count = total_count / ((double)nx * ny * n_frames_used);
			std::cout << "Average count per pixel per frame: " << avg_count << std::endl;

			#pragma omp parallel for num_threads(n_threads)
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Isums[0]())
			{
				if (DIRECT_MULTIDIM_ELEM(Isums[0](), n) == 0)
					DIRECT_MULTIDIM_ELEM(Isums[0](), n) = 1.0;
				else
					DIRECT_MULTIDIM_ELEM(Isums[0](), n) = n_frames_used / DIRECT_MULTIDIM_ELEM(Isums[0](), n) * avg_count;
			}
		}

		Isums[0].write(fn_out);
		std::cout << "Written the estimated gain to " << fn_out << std::endl;
	}
};

int main(int argc, char **argv)
{
	estimate_gain app;
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
