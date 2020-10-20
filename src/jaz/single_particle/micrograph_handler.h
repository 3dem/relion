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

#ifndef MICROGRAPH_HANDLER_H
#define MICROGRAPH_HANDLER_H

#include <string>
#include <map>

#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/single_particle/parallel_ft.h>

#include <src/micrograph_model.h>
#include <src/image.h>

class MicrographHandler
{
	public:

	MicrographHandler();
	

		int nr_omp_threads, firstFrame, lastFrame;
		double movie_angpix, coords_angpix, data_angpix, hotCutoff;
		int eer_upsampling, eer_grouping;
	
		bool debug, saveMem, ready;
	
		std::string corrMicFn;
	
		gravis::t2Vector<int> micrograph_size;


	// initialise corrected/uncorrected micrograph dictionary, then
	// load first movie (or read corrected_micrographs.star) to obtain:
	//   fc, micrograph_xsize, micrograph_ysize, motionEstimator.dosePerFrame
	void init(
		// in:
		const std::vector<MetaDataTable>& mdts,
		bool verb,
		int nr_omp_threads,
		// out:
		int& fc,
		double& dosePerFrame,
		std::string& metaFn);

	void validatePixelSize(RFLOAT angpix) const;

	// remove movies from the list for which either the meta-star or the movie itself is missing
	std::vector<MetaDataTable> cullMissingMovies(
		const std::vector<MetaDataTable>& mdts, int verb);

	// find the greatest number of frames available in all micrographs
	void findLowestFrameCount(
		const std::vector<MetaDataTable>& mdts, int verb);

	// find all movies of sufficient length
	std::vector<MetaDataTable> findLongEnoughMovies(
		const std::vector<MetaDataTable>& mdts,
		int fc, int verb);

	// load a movie and extract all particles
	// returns a per-particle vector of per-frame images of size (s/2+1) x s
	std::vector<std::vector<Image<Complex>>> loadMovie(
		const MetaDataTable& mdt, int s, double angpix, 
		std::vector<ParFourierTransformer>& fts,
		const std::vector<std::vector<gravis::d2Vector>>* offsets_in = 0,
		std::vector<std::vector<gravis::d2Vector>>* offsets_out = 0,
		double data_angpix = -1,
		int single_frame_relative_index = -1); // if this isn't negative, return a single frame for all particles

	/* Write the initial tracks of particles at 'pos' into 'tracks_out' 
	   (by interpolating  the polynomial motionCor2 model).
	   If 'unregGlob' is set, also write the global component of motion into 'globalComponent_out'.*/
	void loadInitialTracks(
		const MetaDataTable& mdt, 
		double angpix, 
		const std::vector<gravis::d2Vector>& pos,
		std::vector<std::vector<gravis::d2Vector>>& tracks_out,
		bool unregGlob, 
		std::vector<gravis::d2Vector>& globalComponent_out);
	
	
	protected:

	
		Image<RFLOAT> lastGainRef;
		MultidimArray<bool> lastDefectMask;
		std::string last_gainFn, last_movieFn;
	
		std::map<std::string, std::string> mic2meta;
		

	void loadInitial(
			const std::vector<MetaDataTable>& mdts, bool verb,
			int& fc, double& dosePerFrame, std::string& metaFn);
		
	std::string getMetaName(
			std::string micName, bool die_on_error = true);
	
	std::string getMicrographMetadataFilename(
			const MetaDataTable& mdt, 
			bool die_on_error);
	
	Micrograph loadMicrographMetadata(
			const MetaDataTable& mdt, 
			bool die_on_error);
	
	int determineFrameCount(
			const MetaDataTable& mdt);
	
	bool isMoviePresent(
			const MetaDataTable& mdt, bool die_on_error = true);
	
	std::string getMovieFilename(
			const MetaDataTable& mdt, bool die_on_error = true);
};

#endif
