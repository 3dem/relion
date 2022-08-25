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
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/movie_loader.h>
#include <src/jaz/single_particle/spa_extraction.h>


using namespace gravis;

MicrographHandler::MicrographHandler()
	: nr_omp_threads(1),
	  firstFrame(0),
	  lastFrame(-1),
	  hotCutoff(-1),
	  debug(false),
	  saveMem(false),
	  ready(false),
	  last_gainFn(""),
	  last_movieFn(""),
	  corrMicFn(""),
	  eer_upsampling(-1),
	  eer_grouping(-1)
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

		mic2meta[fn_post] = metaName;
	}

	loadInitial(mdts, verb, fc, dosePerFrame, metaFn);

	ready = true;
}

std::vector<MetaDataTable> MicrographHandler::cullMissingMovies(
		const std::vector<MetaDataTable> &mdts, int verb)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: MicrographHandler::cullMissingMovies - MicrographHandler not initialized.");
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
	std::string mgFn;
	FileName fn_pre, fn_jobnr, fn_post;
	
	const int mc = mdts.size();

	for (int m = 0; m < mc; m++)
	{
		mdts[m].getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);

		// remove the pipeline job prefix
		decomposePipelineFileName(mgFn, fn_pre, fn_jobnr, fn_post);

		metaFn = getMetaName(fn_post, false);
		if (metaFn != "") break;
	}

	if (metaFn == "")
	{
		REPORT_ERROR("There is no movie metadata STAR file for any micrographs!");
	}

	if (debug)
	{
		std::cout << "first movie: " << fn_post << "\n";
		std::cout << "maps to: " << metaFn << "\n";
	}

	Micrograph micrograph = Micrograph(metaFn);

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

void MicrographHandler::validatePixelSize(RFLOAT angpix) const
{
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
		std::vector<std::vector<gravis::d2Vector>>* offsets_out,
		double data_angpix,
		int single_frame_relative_index)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: MicrographHandler::loadMovie - MicrographHandler not initialized.");
	}
	

	const int nr_omp_threads = fts.size();

	std::string metaFn = getMicrographMetadataFilename(mdt, true);
	Micrograph micrograph = Micrograph(metaFn);

	FileName movieFn = micrograph.getMovieFilename();
	std::string gainFn = micrograph.getGainFilename();
	MultidimArray<bool> defectMask;

	const bool mgHasGain = (gainFn != "");
	const bool hasDefect = (mgHasGain || micrograph.fnDefect != "" || micrograph.hotpixelX.size() != 0);
	
	if (hasDefect)
	{
		if (movieFn == last_movieFn)
		{
			defectMask = lastDefectMask;
		}
		else
		{
			micrograph.fillDefectAndHotpixels(defectMask);
			lastDefectMask = defectMask;
		}
	}
	last_movieFn = movieFn;

	if (debug)
	{
		std::cout << "loading: " << "\n";
		std::cout << "-> meta: " << metaFn << "\n";
		std::cout << "-> data: " << movieFn << "\n";
		std::cout << "-> gain: " << gainFn << "\n";
		std::cout << "-> mask: " << micrograph.fnDefect << "\n";
		std::cout << "-> nhot: " << micrograph.hotpixelX.size() << "\n";
		std::cout << "-> hasdefect: " << (hasDefect ? 1 : 0) << std::endl;
	}

	const bool isEER = EERRenderer::isEER(movieFn);

	if (mgHasGain)
	{
		if (gainFn != last_gainFn)
		{
			last_gainFn = gainFn;
			
			if (isEER) // TODO: Takanori: Remove this once we updated RelionCor
			{
				if (eer_upsampling < 0)
					eer_upsampling = micrograph.getEERUpsampling();
				EERRenderer::loadEERGain(gainFn, lastGainRef(), eer_upsampling);
			}
			else
			{
				lastGainRef.read(gainFn);
			}
		}

		// Mask pixels with zero gain
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(defectMask)
			if (DIRECT_MULTIDIM_ELEM(lastGainRef(), n) == 0)
				DIRECT_MULTIDIM_ELEM(defectMask, n) = true;
	}
	
	BufferedImage<float> muGraph;
	
	RawImage<RFLOAT> gainRef_new(lastGainRef);
	RawImage<bool> defectMask_new(defectMask);
	
	RawImage<RFLOAT>* gainRefToUse = mgHasGain? &gainRef_new : 0;
	RawImage<bool>* defectMaskToUse = hasDefect? &defectMask_new : 0;
	
	const bool returnSingleFrame = single_frame_relative_index >= 0;
	
	if (returnSingleFrame && offsets_in != 0 && (*offsets_in)[0].size() != 1)
	{
		REPORT_ERROR_STR("MicrographHandler::loadMovie: attempting to read one single frame "
						 << "while the initial trajectories contain more than one position");
	}

	const int frame0 = returnSingleFrame? single_frame_relative_index : firstFrame;
	const int fc = returnSingleFrame? 1 : lastFrame - firstFrame + 1;
				
	if (isEER)			
	{
		if (eer_upsampling < 0)
		{
			eer_upsampling = micrograph.getEERUpsampling();
		}
		
		if (eer_grouping < 0)
		{
			eer_grouping = micrograph.getEERGrouping();
		}

		muGraph = MovieLoader::readEER<float>(
			movieFn, gainRefToUse, defectMaskToUse,
			frame0, fc,
			eer_upsampling, eer_grouping,
			nr_omp_threads);
	}
	else
	{
		muGraph = MovieLoader::readDense<float>(
			movieFn, gainRefToUse, defectMaskToUse,
			frame0, fc,
			hotCutoff,
			nr_omp_threads);
	}
	
	std::vector<std::vector<Image<Complex>>> movie = SpaExtraction::extractMovieStackFS(
			mdt, muGraph, s,
			angpix, coords_angpix, movie_angpix, data_angpix,
			offsets_in, offsets_out, 
			nr_omp_threads);
	
	const int pc = movie.size();

	if (!returnSingleFrame)
	{
		#pragma omp parallel for num_threads(nr_omp_threads)
		for (int p = 0; p < pc; p++)
		{
			RFLOAT scale2 = StackHelper::computePower(movie[p], false);
			
			for (int f = 0; f < fc; f++)
			{
				// NOTE: sqrt(fc) is a legacy scaling factor.
				// It probably shouldn't be there.
				//                   --JZ
				RawImage<Complex>(movie[p][f]) /= s * sqrt(scale2 / fc); 
			}
		}
	}

	return movie;
}

void MicrographHandler::loadInitialTracks(
		const MetaDataTable &mdt, double angpix,
		const std::vector<d2Vector>& pos,
		std::vector<std::vector<d2Vector>>& tracks_out,
		bool unregGlob, 
		std::vector<d2Vector>& globalComponent_out)
{	
	Micrograph micrograph = loadMicrographMetadata(mdt, true);
	
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

	globalComponent_out = std::vector<d2Vector>(fc, d2Vector(0,0));

	if (unregGlob)
	{
		for (int f = 0; f < fc; f++)
		{
			RFLOAT sx, sy;
			micrograph.getShiftAt(firstFrame + f + 1, 0, 0, sx, sy, false);

			globalComponent_out[f] = -outputScale * d2Vector(sx, sy);
		}
	}

	tracks_out.resize(pc);

	for (int p = 0; p < pc; p++)
	{
		tracks_out[p] = std::vector<d2Vector>(fc);

		for (int f = 0; f < fc; f++)
		{
			d2Vector in(inputScale.x * pos[p].x - 0.5,
						inputScale.y * pos[p].y - 0.5);

			RFLOAT sx, sy;
			micrograph.getShiftAt(firstFrame + f + 1, in.x, in.y, sx, sy, true);

			tracks_out[p][f] = -outputScale * d2Vector(sx,sy) - globalComponent_out[f];
		}
	}
}

std::string MicrographHandler::getMetaName(std::string micName, bool die_on_error)
{
	std::map<std::string, std::string>::iterator it = mic2meta.find(micName);

	if (it == mic2meta.end())
	{
		if (die_on_error)
		{
			REPORT_ERROR("ERROR: MicrographHandler::getMetaName: no metadata star-file for "
			              +micName+" found in "+corrMicFn+".");
		}
		else
		{
			return "";
		}
	}
	else
	{
		return it->second;
	}
}

std::string MicrographHandler::getMicrographMetadataFilename(
		const MetaDataTable &mdt, 
		bool die_on_error)
{
	std::string mgFn;
	// all entries in the MetaDataTable mdt (particles) have the same micrograph name
	mdt.getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);
	
	FileName fn_pre, fn_jobnr, fn_post;
	decomposePipelineFileName(mgFn, fn_pre, fn_jobnr, fn_post);	
	
	return getMetaName(fn_post, die_on_error);
}

Micrograph MicrographHandler::loadMicrographMetadata(
		const MetaDataTable &mdt,
		bool die_on_error)
{
	return Micrograph(getMicrographMetadataFilename(mdt, die_on_error));
}

// TODO: TAKANORI: This needs to handle changes in EER grouping
int MicrographHandler::determineFrameCount(const MetaDataTable &mdt)
{
	Micrograph micrograph = loadMicrographMetadata(mdt, true); // TODO: TAKANORI: shouldn't read hot pixels

	if (!exists(micrograph.getMovieFilename()))
	{
		return -1;
	}
	else
	{
		return micrograph.getNframes();
	}
}

bool MicrographHandler::isMoviePresent(const MetaDataTable &mdt, bool die_on_error)
{	
	std::string metaFn = getMicrographMetadataFilename(mdt, die_on_error);
	
	if (exists(metaFn))
	{
		Micrograph micrograph = Micrograph(metaFn);
		
		return exists(micrograph.getMovieFilename());
	}
	else
	{
		return false;
	}
}

std::string MicrographHandler::getMovieFilename(const MetaDataTable& mdt, bool die_on_error)
{
	std::string metaFn = getMicrographMetadataFilename(mdt, die_on_error);

	if (exists(metaFn))
	{
		Micrograph micrograph = Micrograph(metaFn);

		return micrograph.getMovieFilename();
	}
	else
	{
		return metaFn;
	}
}
