/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
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

#include "micrograph_handler.h"
#include <src/jaz/stack_helper.h>
#include <src/renderEER.h>

using namespace gravis;

MicrographHandler::MicrographHandler()
	: hasCorrMic(false),
	  nr_omp_threads(1),
	  firstFrame(0),
	  lastFrame(-1),
	  hotCutoff(-1),
	  debug(false),
	  saveMem(false),
	  ready(false),
	  last_gainFn(""),
	  corrMicFn("")
{}

void MicrographHandler::init(
		// in:
		const std::vector<MetaDataTable>& mdts,
		bool verb,
		int nr_omp_threads,
		// out:
		int& fc,
		double& dosePerFrame,
		std::string& metaFn)
{
	this->nr_omp_threads = nr_omp_threads;
	this->firstFrame = firstFrame;
	this->lastFrame = lastFrame;

	if (corrMicFn != "")
	{
		MetaDataTable corrMic;
		ObservationModel obsModel;
		// Don't die even if conversion failed. Polishing does not use obsModel from a motion correction STAR file
		ObservationModel::loadSafely(corrMicFn, obsModel, corrMic, "micrographs", verb, false);
		mic2meta.clear();

		std::string micName, metaName;

		if (!corrMic.containsLabel(EMDL_MICROGRAPH_NAME))
		{
			REPORT_ERROR(" The corrected_micrographs STAR file does not contain rlnMicrographName label.");
		}
		if (!corrMic.containsLabel(EMDL_MICROGRAPH_METADATA_NAME))
		{
			REPORT_ERROR(" The corrected_micrographs STAR file does not contain rlnMicrographMetadata label. Did you not run motion correction from the RELION-3.0 GUI?");
		}

		for (int i = 0; i < corrMic.numberOfObjects(); i++)
		{
			corrMic.getValueToString(EMDL_MICROGRAPH_NAME, micName, i);
			corrMic.getValueToString(EMDL_MICROGRAPH_METADATA_NAME, metaName, i);
			// remove the pipeline job prefix
			FileName fn_pre, fn_jobnr, fn_post;
			decomposePipelineFileName(micName, fn_pre, fn_jobnr, fn_post);

//			std::cout << fn_post << " => " << metaName << std::endl;
			mic2meta[fn_post] = metaName;
		}

		hasCorrMic = true;
	}
	else
	{
		hasCorrMic = false;
	}

	loadInitial(mdts, verb, fc, dosePerFrame, metaFn);

	ready = true;
}

std::vector<MetaDataTable> MicrographHandler::cullMissingMovies(
		const std::vector<MetaDataTable> &mdts, int verb)
{
	if (!ready)
	{
//		REPORT_ERROR("ERROR: MicrographHandler::cullMissingMovies - MicrographHandler not initialized.");
	}

	std::vector<MetaDataTable> good(0);
	std::vector<std::string> bad(0);

	const int mc = mdts.size();

	for (int m = 0; m < mc; m++)
	{
		if (isMoviePresent(mdts[m], false))
		{
			good.push_back(mdts[m]);
		}
		else
		{
			FileName fn_movie;
			mdts[m].getValueToString(EMDL_MICROGRAPH_NAME, fn_movie, 0);
			bad.push_back(fn_movie);
		}
	}

	if (verb && bad.size() > 0)
	{
		if (bad.size() == 1)
		{
			std::cerr << " - The movie for the following micrograph is missing:\n";
		}
		else
		{
			std::cerr << " - Movies for the following micrographs are missing:\n";
		}

		for (int i = 0; i < bad.size(); i++)
		{
			std::cerr << "       " << bad[i] << "\n";
		}
	}

	return good;
}

void MicrographHandler::findLowestFrameCount(
		const std::vector<MetaDataTable> &mdts, int verb)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: MicrographHandler::findLowestFrameCount - MicrographHandler not initialized.");
	}

	int fcmin = std::numeric_limits<int>::max();
	const int mc = mdts.size();

	for (int m = 0; m < mc; m++)
	{
		int fcm = determineFrameCount(mdts[m]);

		if (fcm < fcmin)
		{
			fcmin = fcm;
		}
	}

	if (lastFrame >= fcmin)
	{
		std::cout << " - Warning: some movies contain only " << fcmin
				  << " frames. Unable to load frames " << (fcmin+1)
				  << ".." << (lastFrame+1) << " ( = --last_frame).\n";
	}
	else if (verb > 0)
	{
		std::cout << " + Max. frame number available in all movies: " << fcmin << "\n";
	}

	if (lastFrame < 0 || lastFrame > fcmin-1)
	{
		lastFrame = fcmin - 1;
	}
}

std::vector<MetaDataTable> MicrographHandler::findLongEnoughMovies(
		const std::vector<MetaDataTable> &mdts, int fc, int verb)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: MicrographHandler::findLongEnoughMovies - MicrographHandler not initialized.");
	}

	std::vector<MetaDataTable> good(0);
	std::vector<std::string> bad(0);

	const int mc = mdts.size();

	for (int m = 0; m < mc; m++)
	{
		int fcm = determineFrameCount(mdts[m]);

		if (fcm < fc)
		{
			bad.push_back(getMovieFilename(mdts[m]));
		}
		else
		{
			good.push_back(mdts[m]);
		}
	}

	if (good.size() == 0)
	{
		REPORT_ERROR_STR("ERROR: Not a single movie contains the requested number of frames ("
						 << fc << ")");
	}

	if (verb && bad.size() > 0)
	{
		if (bad.size() == 1)
		{
			std::cerr << " - The following micrograph does not contain "
					  << fc << " frames. Particles in it will be ignored:\n";
		}
		else
		{
			std::cerr << " - The following micrographs do not contain "
					  << fc << " frames. Particles in them will be ignored:\n";
		}

		for (int i = 0; i < bad.size(); i++)
		{
			std::cerr << "       " << bad[i] << "\n";
		}
	}

	return good;
}

// This reads pixel sizes from a single metadata star file.
// For multi optics group scenarios, we should process only micrographs
// in the given MotionCorr STAR file. Then we can safely assume all pixel sizes are the same.
// TODO: TAKANORI: make sure in MotionCorr runner and Polish
void MicrographHandler::loadInitial(
		const std::vector<MetaDataTable>& mdts, bool verb,
		int& fc, double& dosePerFrame, std::string& metaFn)
{
	if (hasCorrMic)
	{
		std::string mgFn;
		FileName fn_pre, fn_jobnr, fn_post;

		for (int i = 0, ilim = mdts.size(); i < ilim; i++)
		{
			mdts[i].getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);

			// remove the pipeline job prefix
			decomposePipelineFileName(mgFn, fn_pre, fn_jobnr, fn_post);

			metaFn = getMetaName(fn_post, false);
			if (metaFn != "") break;
		}

		if (metaFn == "")
			REPORT_ERROR("There is no movie metadata STAR file for any micrographs!");

		if (debug)
		{
			std::cout << "first movie: " << fn_post << "\n";
			std::cout << "maps to: " << metaFn << "\n";
		}

		micrograph = Micrograph(metaFn);

		if (movie_angpix <= 0)
		{
			movie_angpix = micrograph.angpix;

			if (verb > 0)
			{
				std::cout << " + Using movie pixel size from " << metaFn << ": "
						  << movie_angpix << " A\n";
			}
		}
		else
		{
			if (verb > 0)
			{
				std::cout << " + Using movie pixel size from command line: "
						  << movie_angpix << " A\n";
			}
		}

		if (coords_angpix <= 0)
		{
			coords_angpix = micrograph.angpix * micrograph.getBinningFactor();

			if (verb > 0)
			{
				std::cout << " + Using coord. pixel size from " << metaFn << ": "
						  << coords_angpix << " A\n";
			}
		}
		else
		{
			if (verb > 0)
			{
				std::cout << " + Using coord. pixel size from command line: "
						  << coords_angpix << " A\n";
			}
		}

		dosePerFrame = micrograph.dose_per_frame;

		micrograph_size.x = micrograph.getWidth();
		micrograph_size.y = micrograph.getHeight();

		if (lastFrame >= micrograph.getNframes())
		{
			REPORT_ERROR_STR("ERROR: There are only " << micrograph.getNframes()
							 << " frames in " << metaFn << " - " << lastFrame
							 << " have been requested using the --lastFrame option.");
		}

		if (lastFrame < 0)
		{
			fc = micrograph.getNframes() - firstFrame;
		}
		else
		{
			fc = lastFrame - firstFrame + 1;
		}
	}
	else
	{
		std::string mgFn0;
		mdts[0].getValueToString(EMDL_MICROGRAPH_NAME, mgFn0, 0);
		FileName fn_pre, fn_jobnr, fn_post;
		decomposePipelineFileName(mgFn0, fn_pre, fn_jobnr, fn_post);

		Image<RFLOAT> dum;
		dum.read(fn_post, false);
		micrograph_size.x = XSIZE(dum());
		micrograph_size.y = YSIZE(dum());

		const int fc0 = dum().zdim > 1? dum().zdim : dum().ndim;

		if (lastFrame < 0)
		{
			fc = fc0 - firstFrame;
		}
		else
		{
			if (lastFrame >= fc0)
			{
				REPORT_ERROR_STR("ERROR: There are only " << micrograph.getNframes()
								 << " frames in " << metaFn << " - " << (lastFrame+1)
								 << " have been requested using the --lastFrame option.");
			}
			else
			{
				fc = lastFrame - firstFrame + 1;
			}
		}
	}
}

void MicrographHandler::validatePixelSize(RFLOAT angpix) const
{
//	std::cout << "angpix = " << angpix << " coords_angpix = " << coords_angpix << " movie_angpix = " << movie_angpix << std::endl;

	if (angpix < coords_angpix - 1e-9)
	{
		std::cerr << "WARNING: pixel size (--angpix) is smaller than the AutoPick pixel size (--coords_angpix)\n";

		if (coords_angpix < angpix + 0.01)
		{
			std::cerr << "        This is probably a rounding error. It is recommended to set --angpix ("
					  << angpix << ") to at least " << coords_angpix << "\n";

		}
	}

	if (angpix < movie_angpix - 1e-9)
	{
		std::cerr << "WARNING: pixel size (--angpix) is smaller than the movie pixel size (--movie_angpix)\n";

		if (movie_angpix < angpix + 0.01)
		{			std::cerr << "        This is probably a rounding error. It is recommended to set --angpix ("
					  << angpix << ") to at least " << movie_angpix << "\n";

		}
	}
}

std::vector<std::vector<Image<Complex>>> MicrographHandler::loadMovie(
		const MetaDataTable &mdt, int s,
		double angpix, std::vector<ParFourierTransformer>& fts,
		const std::vector<std::vector<gravis::d2Vector>>* offsets_in,
		std::vector<std::vector<gravis::d2Vector>>* offsets_out)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: MicrographHandler::loadMovie - MicrographHandler not initialized.");
	}

	std::vector<std::vector<Image<Complex>>> movie;

	const int nr_omp_threads = fts.size();

	std::string mgFn0;
	mdt.getValueToString(EMDL_MICROGRAPH_NAME, mgFn0, 0);
	FileName fn_pre, fn_jobnr, fn_post;
	decomposePipelineFileName(mgFn0, fn_pre, fn_jobnr, fn_post);

	if (hasCorrMic)
	{
		std::string metaFn = getMetaName(fn_post);
		micrograph = Micrograph(metaFn);

		FileName mgFn = micrograph.getMovieFilename();
		std::string gainFn = micrograph.getGainFilename();
		MultidimArray<bool> defectMask;

		bool hasDefect = (micrograph.fnDefect != "" || micrograph.hotpixelX.size() != 0);
		if (hasDefect)
			micrograph.fillDefectAndHotpixels(defectMask);

		if (debug)
		{
			std::cout << "loading: " << fn_post << "\n";
			std::cout << "-> meta: " << metaFn << "\n";
			std::cout << "-> data: " << mgFn << "\n";
			std::cout << "-> gain: " << gainFn << "\n";
			std::cout << "-> mask: " << micrograph.fnDefect << "\n";
			std::cout << "-> nhot: " << micrograph.hotpixelX.size() << "\n";
			std::cout << "-> hasdefect: " << (hasDefect ? 1 : 0) << std::endl;
		}

		const bool isEER = EERRenderer::isEER(mgFn);

		bool mgHasGain = false;

		if (gainFn != "")
		{
			if (gainFn != last_gainFn)
			{
				lastGainRef.read(gainFn);
				if (isEER) // TODO: Takanori: Remove this once we updated RelionCor
					EERRenderer::upsampleEERGain(lastGainRef());

				last_gainFn = gainFn;
			}

			mgHasGain = true;
		}

		if (!isEER)
		{
#define OLD_CODE
#ifdef OLD_CODE
			movie = StackHelper::extractMovieStackFS(&mdt, mgHasGain? &lastGainRef : 0, hasDefect ? &defectMask : 0,
			                                         mgFn, angpix, coords_angpix, movie_angpix, s,
			                                         nr_omp_threads, true, firstFrame, lastFrame,
			                                         hotCutoff, debug, saveMem, offsets_in, offsets_out);
#else
			// TODO: Implement gain and defect correction, and remove the old code path
			std::cout << "New code path" << std::endl;

			Image<float> mgStack;
			mgStack.read(mgFn, false);

			// lastFrame and firstFrame is 0 indexed
			const int my_lastFrame = ((mgStack.data.zdim > 1)? mgStack.data.zdim : mgStack.data.ndim) - 1;
			const int n_frames = my_lastFrame - firstFrame + 1;

			std::cout << "first = " << firstFrame << " last = " << my_lastFrame << " n_frames = " << n_frames << std::endl;

			std::vector<MultidimArray<float> > Iframes(n_frames);

			#pragma omp parallel for num_threads(nr_omp_threads)
			for (int iframe = 0; iframe < n_frames; iframe++)
			{
				Image<float> img;
				img.read(mgFn, true, iframe, false, true); 
				Iframes[iframe] = img();
			}

			movie = StackHelper::extractMovieStackFS(&mdt, Iframes, angpix, coords_angpix, movie_angpix, s,
			                                         nr_omp_threads, true,
			                                         debug, offsets_in, offsets_out);
#endif
		}
		else
		{
			EERRenderer renderer;
			renderer.read(mgFn);

			// lastFrame and firstFrame is 0 indexed
			int my_lastFrame = (lastFrame < 0) ? (renderer.getNFrames() / EER_grouping - 1) : lastFrame;
			int n_frames = my_lastFrame - firstFrame + 1;

			std::vector<MultidimArray<float> > Iframes(n_frames);

			#pragma omp parallel for num_threads(nr_omp_threads)
			for (int iframe = 0; iframe < n_frames; iframe++)
			{
				// this takes 1-indexed frame numbers
//				std::cout << "EER: iframe = " << iframe << " start = " << ((firstFrame + iframe) * EER_grouping + 1) << " end = " << ((firstFrame + iframe + 1) * EER_grouping) << std::endl;
				renderer.renderFrames((firstFrame + iframe) * EER_grouping + 1, (firstFrame + iframe + 1) * EER_grouping, Iframes[iframe]);

				if (mgHasGain)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(lastGainRef())
					{
						DIRECT_MULTIDIM_ELEM(Iframes[iframe], n) *= DIRECT_MULTIDIM_ELEM(lastGainRef(), n);
					}
				}

			}

			if (hasDefect) // TODO: TAKANORI: code duplication from RelionCor
			{
				if (XSIZE(defectMask) != XSIZE(Iframes[0]) || YSIZE(defectMask) != YSIZE(Iframes[0]))
				{
					std::cerr << "X/YSIZE of defectMask = " << XSIZE(defectMask) << " x " << YSIZE(defectMask) << std::endl;
					std::cerr << "X/YSIZE of Iframe[0] = " << XSIZE(Iframes[0]) << " x " << YSIZE(Iframes[0]) << std::endl;
					REPORT_ERROR("Invalid dfefect mask size for " + mgFn0);
				}

				MultidimArray<float> Isum;
				Isum.initZeros(Iframes[0]);
				for (int iframe = 0; iframe < n_frames; iframe++)
				{
					#pragma omp parallel for num_threads(nr_omp_threads)
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Isum)
					{
						DIRECT_MULTIDIM_ELEM(Isum, n) += DIRECT_MULTIDIM_ELEM(Iframes[iframe], n);
					}
				}
#ifdef DEBUG
				Image<float> tmp;
				tmp() = Isum;
				tmp.write("Isum.mrc");
#endif

				RFLOAT mean = 0, std = 0;
				#pragma omp parallel for reduction(+:mean) num_threads(nr_omp_threads)
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Isum) {
					mean += DIRECT_MULTIDIM_ELEM(Isum, n);
				}
				mean /= YXSIZE(Isum);
				#pragma omp parallel for reduction(+:std) num_threads(nr_omp_threads)
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Isum) {
					RFLOAT d = (DIRECT_MULTIDIM_ELEM(Isum, n) - mean);
					std += d * d;
				}
				std = std::sqrt(std / YXSIZE(Isum));

				mean /= n_frames;
				std /= n_frames;
				Isum.clear();

//				std::cout << "DEBUG: defect correction: mean = " << mean << " std = " << std << std::endl;

				// 25 neighbours; should be enough even for super-resolution images.
		                const int NUM_MIN_OK = 6;
		                const int D_MAX = 2;
		                FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(defectMask)
        		        {
		                        if (!DIRECT_A2D_ELEM(defectMask, i, j)) continue;

		                        #pragma omp parallel for num_threads(nr_omp_threads)
					for (int iframe = 0; iframe < n_frames; iframe++)
					{
		 				int n_ok = 0;
						RFLOAT val = 0;
						for (int dy= -D_MAX; dy <= D_MAX; dy++)
						{
							int y = i + dy;
							if (y < 0 || y >= YSIZE(defectMask)) continue;
							for (int dx = -D_MAX; dx <= D_MAX; dx++)
							{
								int x = j + dx;
								if (x < 0 || x >= XSIZE(defectMask)) continue;
								if (DIRECT_A2D_ELEM(defectMask, y, x)) continue;

								n_ok++;
								val += DIRECT_A2D_ELEM(Iframes[iframe], y, x);
							}
						}
//						std::cout << "n_ok = " << n_ok << " val = " << val << std::endl;
						if (n_ok > NUM_MIN_OK) DIRECT_A2D_ELEM(Iframes[iframe], i, j) = val / n_ok;
						else DIRECT_A2D_ELEM(Iframes[iframe], i, j) = rnd_gaus(mean, std);
					}
				}

#ifdef DEBUG
				Isum.initZeros(Iframes[0]);
				for (int iframe = 0; iframe < n_frames; iframe++)
				{
					#pragma omp parallel for num_threads(nr_omp_threads)
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Isum)
					{
						DIRECT_MULTIDIM_ELEM(Isum, n) += DIRECT_MULTIDIM_ELEM(Iframes[iframe], n);
					}
				}

				tmp() = Isum;
				tmp.write("Isum-fix-defect.mrc");
				exit(0);
#endif
			}

			movie = StackHelper::extractMovieStackFS(&mdt, Iframes, angpix, coords_angpix, movie_angpix, s,
			                                         nr_omp_threads, true,
			                                         debug, offsets_in, offsets_out);
		}
	}
	else
	{
		REPORT_ERROR("You can no longer use this program without micrograph metadata STAR files.");
	}

	const int pc = movie.size();

	#pragma omp parallel for num_threads(nr_omp_threads)
	for (int p = 0; p < pc; p++)
	{
		StackHelper::varianceNormalize(movie[p], false);
	}

	return movie;
}

std::vector<std::vector<Image<Complex>>> MicrographHandler::loadMovie(
		const MetaDataTable &mdt, int s, double angpix,
		std::vector<ParFourierTransformer>& fts,
		const std::vector<d2Vector>& pos,
		std::vector<std::vector<d2Vector>>& tracks,
		bool unregGlob, std::vector<d2Vector>& globComp,
		const std::vector<std::vector<gravis::d2Vector>>* offsets_in,
		std::vector<std::vector<gravis::d2Vector>>* offsets_out)
{
	std::vector<std::vector<Image<Complex>>> out = loadMovie(
				mdt, s, angpix, fts, offsets_in, offsets_out);

	if (!hasCorrMic)
	{
		tracks.resize(0);
	}
	else
	{
		const int fc0 = micrograph.getNframes();
		int fc;

		if (lastFrame >= 0)
		{
			fc = lastFrame - firstFrame + 1;
		}
		else
		{
			fc = fc0 - firstFrame;
		}

		const int pc = pos.size();

		const d2Vector inputScale(
					coords_angpix / (movie_angpix * micrograph.getWidth()),
					coords_angpix / (movie_angpix * micrograph.getHeight()));

		const double outputScale = movie_angpix / angpix;

		globComp = std::vector<d2Vector>(fc, d2Vector(0,0));

		if (unregGlob)
		{
			for (int f = 0; f < fc; f++)
			{
				RFLOAT sx, sy;
				micrograph.getShiftAt(firstFrame + f + 1, 0, 0, sx, sy, false);

				globComp[f] = -outputScale * d2Vector(sx, sy);
			}
		}

		tracks.resize(pc);

		for (int p = 0; p < pc; p++)
		{
			tracks[p] = std::vector<d2Vector>(fc);

			for (int f = 0; f < fc; f++)
			{
				d2Vector in(inputScale.x * pos[p].x - 0.5,
				            inputScale.y * pos[p].y - 0.5);

				RFLOAT sx, sy;
				micrograph.getShiftAt(firstFrame + f + 1, in.x, in.y, sx, sy, true);

				tracks[p][f] = -outputScale * d2Vector(sx,sy) - globComp[f];
			}
		}
	}

	return out;
}

std::string MicrographHandler::getMetaName(std::string micName, bool die_on_error)
{
//	std::cout << "MicrographHandler::getMetaName " << micName << std::endl;
	std::map<std::string, std::string>::iterator it = mic2meta.find(micName);

	if (it == mic2meta.end())
	{
		if (die_on_error)
			REPORT_ERROR("ERROR: MicrographHandler::getMetaName: no metadata star-file for "
			              +micName+" found in "+corrMicFn+".");
		else
			return "";
	}
	else
	{
		return it->second;
	}
}

int MicrographHandler::determineFrameCount(const MetaDataTable &mdt)
{
	int fc = 0;

	std::string mgFn;
	mdt.getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);

	FileName fn_pre, fn_jobnr, fn_post;
	decomposePipelineFileName(mgFn, fn_pre, fn_jobnr, fn_post);

	if (hasCorrMic)
	{
		std::string metaFn = getMetaName(fn_post);
		micrograph = Micrograph(metaFn);

		if (!exists(micrograph.getMovieFilename()))
		{
			return -1;
		}

		fc = micrograph.getNframes();
	}
	else
	{
		if (!exists(fn_post))
		{
			return -1;
		}

		Image<RFLOAT> dum;
		dum.read(fn_post, false);

		fc = dum().zdim > 1? dum().zdim : dum().ndim;
	}

	return fc;
}

bool MicrographHandler::isMoviePresent(const MetaDataTable &mdt, bool die_on_error)
{
	std::string mgFn;
	FileName fn_pre, fn_jobnr, fn_post;
	mdt.getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);
	decomposePipelineFileName(mgFn, fn_pre, fn_jobnr, fn_post);

	if (hasCorrMic)
	{
		std::string metaFn = getMetaName(fn_post, die_on_error);
		if (exists(metaFn))
		{
			micrograph = Micrograph(metaFn);

			return exists(micrograph.getMovieFilename());
		}
		else
		{
			return false;
		}
	}
	else
	{
		return exists(fn_post);
	}
}

std::string MicrographHandler::getMovieFilename(const MetaDataTable& mdt, bool die_on_error)
{
	std::string mgFn;
	mdt.getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);

	FileName fn_pre, fn_jobnr, fn_post;
	decomposePipelineFileName(mgFn, fn_pre, fn_jobnr, fn_post);

	if (hasCorrMic)
	{
		std::string metaFn = getMetaName(fn_post, die_on_error);

		if (exists(metaFn))
		{
			micrograph = Micrograph(metaFn);

			return micrograph.getMovieFilename();
		}
		else
		{
			return metaFn;
		}
	}
	else
	{
		return fn_post;
	}
}
