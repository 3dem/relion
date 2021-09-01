/***************************************************************************
 *
 * Authors: "Jasenko Zivanov & Sjors H.W. Scheres"
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

#include "motion_refiner.h"

#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/slice_helper.h>
#include <src/jaz/single_particle/spectral_helper.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/complex_io.h>
#include <src/jaz/single_particle/fftw_helper.h>
#include <src/jaz/single_particle/refinement_helper.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/damage_helper.h>
#include <src/jaz/single_particle/fsc_helper.h>
#include <src/jaz/single_particle/img_proc/image_op.h>
#include <src/jaz/single_particle/parallel_ft.h>
#include <src/renderEER.h>

#include "gp_motion_fit.h"
#include "motion_helper.h"

using namespace gravis;

MotionRefiner::MotionRefiner()
:	motionParamEstimator(),
	motionEstimator(),
	frameRecombiner()
{
}

void MotionRefiner::read(int argc, char **argv)
{
	IOParser parser;
	parser.setCommandLine(argc, argv);
	
	parser.addSection("General options");
	// TODO: fn_opt = parser.getOption("--opt", "optimiser STAR file from a previous 3D auto-refinement");
	
	starFn = parser.getOption("--i", "Input STAR file");
	outPath = parser.getOption("--o", "Output directory, e.g. MotionFit/job041/");
	
	reference.read(parser, argc, argv);
	
	micrographHandler.firstFrame = textToInteger(parser.getOption("--first_frame", "First move frame to process", "1")) - 1;
	micrographHandler.lastFrame = textToInteger(parser.getOption("--last_frame", "Last movie frame to process (default is all)", "-1")) - 1;
	only_do_unfinished = parser.checkOption("--only_do_unfinished", "Skip those steps for which output files already exist.");
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));
	
	motionEstimator.read(parser, argc, argv);
	motionParamEstimator.read(parser, argc, argv);
	frameRecombiner.read(parser, argc, argv);
	
	parser.addSection("Computational options");
	
	nr_omp_threads = textToInteger(parser.getOption("--j", "Number of (OMP) threads", "1"));
	particlesForFcc = textToInteger(parser.getOption("--B_parts", "Number of particles used for B-factor estimation (negative means all)", "-1"));
	minMG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
	maxMG = textToInteger(parser.getOption("--max_MG", "Last micrograph index (default is to process all)", "-1"));
	
	micrographHandler.saveMem = parser.checkOption("--sbs", "Load movies slice-by-slice to save memory (slower)");
	
	parser.addSection("Expert options");
	
	findShortestMovie = parser.checkOption("--find_shortest", "Load only as many frames as are present in all movies.");
	debug = parser.checkOption("--debug", "Write debugging data");
	
	micrographHandler.corrMicFn = parser.getOption("--corr_mic", "List of uncorrected micrographs (e.g. corrected_micrographs.star)");
	micrographHandler.movie_angpix = textToDouble(parser.getOption("--mps", "Pixel size of input movies (Angst/pix)", "-1"));
	micrographHandler.coords_angpix = textToDouble(parser.getOption("--cps", "Pixel size of particle coordinates in star-file (Angst/pix)", "-1"));
	micrographHandler.hotCutoff = textToDouble(parser.getOption("--hot", "Clip hot pixels to this max. value (-1 = off, TIFF only)", "-1"));
	micrographHandler.debug = parser.checkOption("--debug_mov", "Write debugging data for movie loading");
	
	movie_toReplace = parser.getOption("--mov_toReplace", "Replace this string in micrograph names...", "");
	movie_replaceBy = parser.getOption("--mov_replaceBy", "..by this one", "");

	micrographHandler.eer_upsampling = textToInteger(parser.getOption("--eer_upsampling", "EER upsampling (1 = 4K or 2 = 8K)", "-1"));
	micrographHandler.eer_grouping = textToInteger(parser.getOption("--eer_grouping", "EER grouping", "-1"));

	if (micrographHandler.eer_upsampling > 0 || micrographHandler.eer_grouping > 0)
	{
		// TODO: TAKANORI: Support changing EER upsampling and EER grouping in Polish.
		// Change in upsampling needs changes in ImageX/Y, Binning, OriginalPixelSize, ShiftX/Y, HotpixelX/Y (but not coeffs)
		// Change in grouping needs changes in ImageZ, DoseRate and interpolation of ShiftX/Y and scaling of MotionModelCoeffs
		std::cerr << "At the moment, Polishing does not support changing --eer_upsampling and --eer_grouping from the values used in MotionCorr job.\n";
		std::cerr << "You need to manually modify trajectory STAR files for it." << std::endl;
	}
	
	// Check for errors in the command-line option
	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	if (argc == 1) parser.writeUsage(std::cerr);
}

void MotionRefiner::init()
{
	if (outPath[outPath.length()-1] != '/')
	{
		outPath += "/";
	}
	
	if (verb > 0) std::cout << " + Reading " << starFn << "..." << std::endl;
	
	mdt0.read(starFn);	
	
	ObservationModel::loadSafely(starFn, obsModel, mdt0);
	
	//@CHECK
	if (!ObservationModel::containsAllColumnsNeededForPrediction(mdt0))
	{
		REPORT_ERROR_STR(starFn << " does not contain all of the required columns ("
			<< "rlnOriginX, rlnOriginY, rlnAngleRot, rlnAngleTilt, rlnAnglePsi and rlnRandomSubset)");
	}
	
	adaptMovieNames();
	
	allMdts = StackHelper::splitByMicrographName(mdt0);

	if (minMG >= allMdts.size())
	{
		std::stringstream sts0, sts1;
		sts0 << minMG;
		sts0 << allMdts.size();
		
		REPORT_ERROR("ERROR: Cannot start with micrograph "+sts0.str()
					 +" (--min_MG); only "+sts1.str()+" micrographs defined in "+starFn+".");
	}
	
	if (minMG < 0)
	{
		minMG = 0;
	}
	
	// Only work on a user-specified subset of the micrographs
	if (maxMG < 0 || maxMG >= allMdts.size())
	{
		maxMG = allMdts.size()-1;
	}
	
	if (minMG > 0 || maxMG < allMdts.size()-1)
	{
		if (verb > 0)
		{
			std::cout << "   - Will only process micrographs in range: ["
					  << minMG << "-" << maxMG << "]"  << std::endl;
		}
		
		chosenMdts.clear();
		
		for (long int g = minMG; g <= maxMG; g++)
		{
			chosenMdts.push_back(allMdts[g]);
		}
	}
	else
	{
		chosenMdts = allMdts;
	}
	
	/* There are two options on how to handle movies of varying length:
	  
	   1)  Find the lowest frame-count in the selected dataset,
		   and only extract that number of frames from all movies.
		   (findShortestMovie == true)
		   
	   2)  Remove all movies of insufficient size from the dataset.
		   If no --last_frame value is supplied, the length of the first
		   movie is used as the reference length
		   (findShortestMovie == false)
	*/
	
	{
		std::string metaFn = ""; // the first meta-star filename
		double fractDose = 0.0; // the dose in the first meta-star (@TODO: support variable dose)
		// => we don't need variable dose support as long as different motioncorr jobs are processed separately
		int fc0; // the frame count in the first movie
		
		// initialise corrected/uncorrected micrograph dictionary, then load the header
		// of the first movie (or read corrected_micrographs.star) to obtain:
		// frame count, micrograph size and the fractional dose
		micrographHandler.init(
			chosenMdts, verb, nr_omp_threads, // in
			fc0, fractDose, metaFn); // out

		chosenMdts = micrographHandler.cullMissingMovies(chosenMdts, verb);
		
		if (!findShortestMovie)
		{
			if (micrographHandler.lastFrame < 0)
			{
				fc = fc0 - micrographHandler.firstFrame;
				micrographHandler.lastFrame = fc0 - 1;
			}
			else
			{
				fc = micrographHandler.lastFrame - micrographHandler.firstFrame + 1;
			}
		}
		
		// metaFn is only needed for console output
		// verb is needed since motionEstimator has not been initialised yet
		motionEstimator.proposeDosePerFrame(fractDose, metaFn, verb);
	}
	
	if (findShortestMovie)
	{
		// make sure we don't try to load too many frames
		micrographHandler.findLowestFrameCount(chosenMdts, verb);
		fc = micrographHandler.lastFrame - micrographHandler.firstFrame + 1;
	}
	else
	{
		// remove movies of insufficient size
		chosenMdts = micrographHandler.findLongEnoughMovies(
					chosenMdts, micrographHandler.lastFrame+1, verb);
	}
	
	if (only_do_unfinished)
	{
		motionUnfinished = MotionEstimator::findUnfinishedJobs(chosenMdts, outPath);
		recombUnfinished = frameRecombiner.findUnfinishedJobs(chosenMdts, outPath);

		motionMdts = selectMicrographs(chosenMdts, motionUnfinished);
		recombMdts = selectMicrographs(chosenMdts, recombUnfinished);
		
		if (verb > 0)
		{
			if (motionMdts.size() > 0)
			{
				if (motionMdts.size() < chosenMdts.size())
				{
					std::cout << "   - Will only estimate motion for " << motionMdts.size()
							  << " unfinished micrographs" << std::endl;
				}
				else
				{
					std::cout << "   - Will estimate motion for all " << motionMdts.size()
							  << " micrographs - none are finished" << std::endl;
				}
			}
			else
			{
				std::cout << "   - Motion has already been estimated for all micrographs" << std::endl;
			}
			
			if (frameRecombiner.doingRecombination())
			{
				if (recombMdts.size() > 0)
				{
					if (recombMdts.size() < chosenMdts.size())
					{
						std::cout << "   - Will only recombine frames for " << recombMdts.size()
								  << " unfinished micrographs" << std::endl;
					}
					else
					{
						std::cout << "   - Will recombine frames for all " << recombMdts.size()
								  << " micrographs - none are finished" << std::endl;
					}
				}
				else
				{
					std::cout << "   - Frames have already been recombined for all micrographs; "
							  << "a new STAR file will be generated" << std::endl;
				}
			}

			if (debug)
			{
				std::cout << "index \t motion unfinished \t recombination unfinished:\n";

				for (int i = 0; i < chosenMdts.size(); i++)
				{
					std::cout << i << " \t " << motionUnfinished[i] << " \t " << recombUnfinished[i] << '\n';
				}

				std::cout << '\n';
			}
		}
	}
	else
	{
		motionMdts = chosenMdts;
		recombMdts = chosenMdts;

		motionUnfinished = std::vector<bool>(chosenMdts.size(), true);
		recombUnfinished = std::vector<bool>(recombMdts.size(), true);
	}
	
	estimateParams = motionParamEstimator.anythingToDo();
	estimateMotion = motionMdts.size() > 0;
	recombineFrames = frameRecombiner.doingRecombination() && 
	                  (recombMdts.size() > 0 || !exists(outPath + "shiny" + frameRecombiner.getOutputSuffix() + ".star"));
	generateStar = frameRecombiner.doingRecombination();
	
	bool doAnything = estimateParams || estimateMotion || recombineFrames;
	bool needsReference = doAnything;

	if (!doAnything) 
		exit(0); //TODO: To be replaced with RELION_EXIT_SUCCESS	
	
	if (needsReference)
	{
		if (verb > 0) std::cout << " + Reading references ..." << std::endl;
		
		reference.load(verb, debug);
	}

	micrographHandler.validatePixelSize(reference.angpix);
	
	if (estimateMotion || estimateParams)
	{
		if (verb > 0) std::cout << " + Initializing motion estimator ..." << std::endl;
		
		motionEstimator.init(
			verb, fc, nr_omp_threads, debug, outPath,
			&reference, &obsModel, &micrographHandler);
	}
	
	if (estimateParams)
	{
		//REPORT_ERROR("Parameter estimation currently not supported");
		
		if (verb > 0) std::cout << " + Initializing motion parameter estimator ..." << std::endl;
		
		motionParamEstimator.init(
			verb, nr_omp_threads, debug, outPath,
			fc, chosenMdts, &motionEstimator, &reference, &obsModel);
	}
}

void MotionRefiner::run()
{
	if (estimateParams)
	{
		motionParamEstimator.run();
		
		return;
		// @TODO: apply the optimized parameters, then continue with motion estimation
	}

	const int mgc = chosenMdts.size();
	
	const int lastTotalMgForFCC = lastTotalMicrographForFCC();

	int lastMotionMgForFCC = subtractFinishedMicrographs(lastTotalMgForFCC, motionUnfinished);
	int lastRecombMgForFCC = subtractFinishedMicrographs(lastTotalMgForFCC, recombUnfinished);

	const int firstMotionMgWithoutFCC = lastMotionMgForFCC + 1;
	const int firstTotalMgWithoutFCC = lastTotalMgForFCC + 1;

	std::vector<int> total2motion = getForwardIndices(motionUnfinished);
	std::vector<int> total2recomb = getForwardIndices(recombUnfinished);
	
	if (estimateMotion)
	{
		motionEstimator.process(motionMdts, 0, lastMotionMgForFCC, true);
	}
	
	if (recombineFrames)
	{
		double k_out_A = reference.pixToAng(reference.k_out);


		// Recombine movies that have already been aligned:
		
		frameRecombiner.init(
			allMdts, verb, reference.s, fc, k_out_A, reference.angpix,
			nr_omp_threads, outPath, debug,
			&reference, &obsModel, &micrographHandler);
		
		frameRecombiner.process(recombMdts, 0, lastRecombMgForFCC);


		// Then, align and recombine in an alternating pattern to minimize disk I/O:

		motionEstimator.setVerbosity(0);
		frameRecombiner.setVerbosity(0);

		const int left = chosenMdts.size() - firstTotalMgWithoutFCC;
		const int barstep = XMIPP_MAX(1, left/ 60);

		if (verb > 0)
		{
			std::cout << " + Aligning and combining frames for micrographs ... " << std::endl;
			init_progress_bar(mgc - firstTotalMgWithoutFCC);
		}

		for (int m = firstTotalMgWithoutFCC; m < mgc; m++)
		{
			// TODO: TAKANORI: micrograph handler can cache movie frames to avoid reading movies twice.
			//                 (if --sbs, don't cache to save memory)

			if (estimateMotion && motionUnfinished[m])
			{
				motionEstimator.process(motionMdts, total2motion[m], total2motion[m], false);
			}

			if (recombUnfinished[m])
			{
				frameRecombiner.process(recombMdts, total2recomb[m], total2recomb[m]);
			}

			const int nr_done = m - firstTotalMgWithoutFCC;

			if (verb > 0 && nr_done % barstep == 0)
			{
				progress_bar(nr_done);
			}
		}

		if (verb > 0)
		{
			progress_bar(left);
		}
	}
	else
	{
		motionEstimator.process(motionMdts, firstMotionMgWithoutFCC, motionMdts.size()-1, true);
	}
	
	if (generateStar)
	{
		combineEPSAndSTARfiles();
	}
}

int MotionRefiner::getVerbosityLevel()
{
	return verb;
}

// combine all EPS files into one logfile.pdf
void MotionRefiner::combineEPSAndSTARfiles()
{
	std::vector<FileName> fn_eps;
	
	if (verb > 0)
	{
		std::cout << " + Combining all EPS and STAR files ... " << std::endl;
	}
	
	MetaDataTable mdtAll;
	
	if (frameRecombiner.doingRecombination())
	{
		if (exists(outPath+"bfactors.eps"))
		{
			fn_eps.push_back(outPath+"bfactors.eps");
		}
		
		if (exists(outPath+"scalefactors.eps"))
		{
			fn_eps.push_back(outPath+"scalefactors.eps");
		}
	}

	std::vector<int> n_OgPresent(obsModel.numberOfOpticsGroups(), 0);
	std::vector<int> n_OgAbsent(obsModel.numberOfOpticsGroups(), 0);
	for (long g = 0; g < allMdts.size(); g++)
	{
		FileName fn_root = getOutputFileNameRoot(outPath, allMdts[g]);
		
		if (exists(fn_root+"_tracks.eps"))
		{
			fn_eps.push_back(fn_root+"_tracks.eps");
		}
		
		if (frameRecombiner.doingRecombination() && exists(fn_root+"_shiny" + frameRecombiner.getOutputSuffix() + ".star"))
		{
			MetaDataTable mdt;
			mdt.read(fn_root+"_shiny" + frameRecombiner.getOutputSuffix() + ".star");
			mdtAll.append(mdt);

			FOR_ALL_OBJECTS_IN_METADATA_TABLE(mdt)
			{
				n_OgPresent[obsModel.getOpticsGroup(mdt)]++;
			}
		}
		else
		{
			// Non-processed particles belonging to micrographs not present in the MotionCorr STAR file
			// Remove them from the output

			FOR_ALL_OBJECTS_IN_METADATA_TABLE(allMdts[g])
			{
				n_OgAbsent[obsModel.getOpticsGroup(allMdts[g])]++;
			}
		}
	}
	
	if (fn_eps.size() > 0)
	{
		joinMultipleEPSIntoSinglePDF(outPath + "logfile.pdf", fn_eps);
	}
	
	if (frameRecombiner.doingRecombination())
	{
		for (int og = 0; og < obsModel.numberOfOpticsGroups(); og++)
		{
			// If this optics group was not processed, don't change anything
			if (n_OgPresent[og] == 0)
			{
				std::cerr << "WARNING: All " << n_OgAbsent[og] << " particles in the optics group " << (og + 1) << " were removed because no particles belong to the movies in the input MotionCorr STAR file." << std::endl;

				obsModel.opticsMdt.setValue(EMDL_IMAGE_PIXEL_SIZE, -1.0, og); // mark for deletion
				continue;
			}

			if (n_OgAbsent[og] > 0)
			{
				std::cerr << "WARNING: " << n_OgAbsent[og] << " particles in the optics group " << (og + 1) << " were removed." << std::endl;
			}

			obsModel.opticsMdt.setValue(EMDL_IMAGE_PIXEL_SIZE, frameRecombiner.getOutputPixelSize(og), og);
			obsModel.opticsMdt.setValue(EMDL_IMAGE_SIZE, frameRecombiner.getOutputBoxSize(og), og);
			obsModel.opticsMdt.setValue(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, frameRecombiner.isCtfMultiplied(og), og);
			std::cout << " + Pixel size for optics group " << (og + 1) << ": " << frameRecombiner.getOutputPixelSize(og) << std::endl;
		}

		// Remove absent optics groups; After this, NOTHING should be done except for saving. obsModel's internal data structure is now corrupted!
		int og = 0;
		while (og < obsModel.opticsMdt.numberOfObjects())
		{
			RFLOAT og_angpix;
			obsModel.opticsMdt.getValue(EMDL_IMAGE_PIXEL_SIZE, og_angpix, og);
			if (og_angpix < 0)
			{
				obsModel.opticsMdt.removeObject(og);
			}
			else
			{
				og++;
			}
		}

		obsModel.save(mdtAll, outPath + "shiny" + frameRecombiner.getOutputSuffix() + ".star");
	}
	
	if (verb > 0)
	{
		std::cout << " + Done! " << std::endl;
		std::cout << " + Written logfile in " << outPath << "logfile.pdf" << std::endl;
		
		if (frameRecombiner.doingRecombination())
		{
			std::cout << " + Written new particle STAR file in "
					  << outPath << "shiny" + frameRecombiner.getOutputSuffix() + ".star" << std::endl;
		}
	}
}

// Get the output filename from the micrograph filename
FileName MotionRefiner::getOutputFileNameRoot(std::string outPath, const MetaDataTable& mdt)
{
	FileName fn_mic, fn_pre, fn_jobnr, fn_post;
	mdt.getValue(EMDL_MICROGRAPH_NAME, fn_mic, 0);
	decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
	return outPath + fn_post.withoutExtension();
}

void MotionRefiner::adaptMovieNames()
{
	if (movie_toReplace != "")
	{
		std::string name;
		
		for (int i = 0; i < mdt0.numberOfObjects(); i++)
		{
			mdt0.getValue(EMDL_MICROGRAPH_NAME, name, i);
			
			if (i == 0 && verb > 0)
			{
				std::cout << "   - Replacing: " << name << std::endl;
			}
			
			std::string::size_type pos0 = name.find(movie_toReplace);
			
			if (pos0 != std::string::npos)
			{
				std::string::size_type pos1 = pos0 + movie_toReplace.length();
				
				std::string before = name.substr(0, pos0);
				std::string after = pos1 < name.length()? name.substr(pos1) : "";
				
				name = before + movie_replaceBy + after;
			}
			
			if (i == 0 && verb > 0)
			{
				std::cout << "                -> " << name << std::endl;
			}
			
			mdt0.setValue(EMDL_MICROGRAPH_NAME, name, i);
		}
	}
}

int MotionRefiner::lastTotalMicrographForFCC()
{
	int total = chosenMdts.size() - 1;

	if (particlesForFcc > 0)
	{
		int partSoFar = 0;
		
		for (int m = 0; m < chosenMdts.size(); m++)
		{
			partSoFar += chosenMdts[m].numberOfObjects();
			
			if (partSoFar > particlesForFcc)
			{
				total = m;
				break;
			}
		}
	}

	return total;
}

int MotionRefiner::subtractFinishedMicrographs(
		int lastTotal,
		const std::vector<bool> &selection)
{
	int lastSelected = lastTotal;

	for (int m = 0; m <= lastTotal; m++)
	{
		if (!selection[m])
		{
			lastSelected--;
		}
	}

	return lastSelected;
}

std::vector<MetaDataTable> MotionRefiner::selectMicrographs(
		const std::vector<MetaDataTable> &mdts,
		const std::vector<bool> &selection) const
{
	std::vector<MetaDataTable> out(0);

	for (int i = 0; i < mdts.size(); i++)
	{
		if (selection[i])
		{
			out.push_back(mdts[i]);
		}
	}

	return out;
}

std::vector<int> MotionRefiner::getForwardIndices(const std::vector<bool>& selection) const
{
	const int mgc = selection.size();
	std::vector<int> total_to_selection(mgc, -1);

	int index = 0;

	for (int m = 0; m < mgc; m++)
	{
		if (selection[m])
		{
			total_to_selection[m] = index;
			index++;
		}
	}

	return total_to_selection;
}

std::vector<int> MotionRefiner::getBackwardIndices(const std::vector<bool>& selection) const
{
	const int mgc = selection.size();
	std::vector<int> selection_to_total(0);
	selection_to_total.reserve(mgc);

	for (int m = 0; m < mgc; m++)
	{
		if (selection[m])
		{
			selection_to_total.push_back(m);
		}
	}

	return selection_to_total;
}


