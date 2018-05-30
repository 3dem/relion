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
#include <src/jaz/image_log.h>
#include <src/jaz/slice_helper.h>
#include <src/jaz/spectral_helper.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/backprojection_helper.h>
#include <src/jaz/complex_io.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/damage_helper.h>
#include <src/jaz/fsc_helper.h>
#include <src/jaz/image_op.h>
#include <src/jaz/parallel_ft.h>

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
	minMG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
	maxMG = textToInteger(parser.getOption("--max_MG", "Last micrograph index (default is to process all)", "-1"));
	
	micrographHandler.saveMem = parser.checkOption("--sbs", "Load movies slice-by-slice to save memory (slower)");
	
	parser.addSection("Expert options");
	
	findShortestMovie = parser.checkOption("--find_shortest", "Load only as many frames as are present in all movies.");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (angst/pix) - read from STAR file by default", "-1"));
	debug = parser.checkOption("--debug", "Write debugging data");
	
	micrographHandler.corrMicFn = parser.getOption("--corr_mic", "List of uncorrected micrographs (e.g. corrected_micrographs.star)", "");
	micrographHandler.movie_angpix = textToFloat(parser.getOption("--mps", "Pixel size of input movies (Angst/pix)", "-1"));
	micrographHandler.coords_angpix = textToFloat(parser.getOption("--cps", "Pixel size of particle coordinates in star-file (Angst/pix)", "-1"));
	micrographHandler.hotCutoff = textToFloat(parser.getOption("--hot", "Clip hot pixels to this max. value (-1 = off, TIFF only)", "-1"));
	micrographHandler.debug = parser.checkOption("--debug_mov", "Write debugging data for movie loading");
	
	movie_toReplace = parser.getOption("--mov_toReplace", "Replace this string in micrograph names...", "");
	movie_replaceBy = parser.getOption("--mov_replaceBy", "..by this one", "");
	
	bfactorEstimator.read(parser, argc, argv);
			
	// Check for errors in the command-line option
	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	if (argc == 1) parser.writeUsage(std::cerr);
}

void MotionRefiner::init()
{
	if (micrographHandler.movie_angpix <= 0 && micrographHandler.corrMicFn == "")
	{
		REPORT_ERROR("ERROR: Movie pixel size (--mps) is required unless a corrected_micrographs.star (--corr_mic) is provided.");
	}
	
	if (micrographHandler.coords_angpix <= 0 && micrographHandler.corrMicFn == "")
	{
		REPORT_ERROR("ERROR: Coordinates pixel size (--cps) is required unless a corrected_micrographs.star (--corr_mic) is provided.");
	}
	
	if (outPath[outPath.length()-1] != '/')
	{
		outPath += "/";
	}
	
	if (verb > 0) std::cout << " + Reading " << starFn << "..." << std::endl;
	
	mdt0.read(starFn);	
	
	if (!ObservationModel::containsAllNeededColumns(mdt0))
	{
		REPORT_ERROR_STR(starFn << " does not contain all of the required columns ("
			<< "rlnOriginX, rlnOriginY, rlnAngleRot, rlnAngleTilt, rlnAnglePsi and rlnRandomSubset)");
	}
	
	adaptMovieNames();
	
	if (angpix <= 0.0)
	{
		RFLOAT mag, dstep;
		mdt0.getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
		mdt0.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
		angpix = 10000 * dstep / mag;
		
		if (verb > 0)
		{
			std::cout << "   - Using pixel size calculated from magnification and "
					  << "detector pixel size in the input STAR file: " << angpix << std::endl;
		}
	}
	
	allMdts = StackHelper::splitByMicrographName(&mdt0);
	
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
		int fc0; // the frame count in the first movie
		
		// initialise corrected/uncorrected micrograph dictionary, then load the header
		// of the first movie (or read corrected_micrographs.star) to obtain:
		//  frame count, micrograph size and the fractional dose
		micrographHandler.init(
					chosenMdts[0], angpix, verb, nr_omp_threads, // in
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
		motionMdts.clear();
		recombMdts.clear();
		
		motionMdts = MotionEstimator::findUnfinishedJobs(chosenMdts, outPath);
		recombMdts = FrameRecombiner::findUnfinishedJobs(chosenMdts, outPath);
		
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
		}
	}
	else
	{
		motionMdts = chosenMdts;
		recombMdts = chosenMdts;
	}
	
	estimateParams = motionParamEstimator.anythingToDo();
	estimateMotion = motionMdts.size() > 0;
	recombineFrames = frameRecombiner.doingRecombination() && (recombMdts.size() > 0);
	estimateOldBfacs = bfactorEstimator.doingAnything();
	generateStar = frameRecombiner.doingRecombination();
	
	
	bool doAnything = estimateParams || estimateMotion || recombineFrames || estimateOldBfacs;
	bool needsReference = estimateParams || estimateMotion 
			|| !frameRecombiner.outerFreqKnown() || estimateOldBfacs;
	
	if (doAnything)
	{
		obsModel = ObservationModel(angpix);
		
		// Read the first reference
		// (even if there is no motion to estimate - only to learn the image size)
		// @TODO: replace this once the data is tree-structured
		Image<RFLOAT> map0;
		map0.read(reference.reconFn0, false);
		
		// Get dimensions
		s = map0.data.xdim;
		sh = s/2 + 1;
		
		if (needsReference)
		{
			if (verb > 0) std::cout << " + Reading references ..." << std::endl;
			
			reference.load(verb, debug);
		}
		
	}
	
	if (estimateMotion || estimateParams)
	{
		if (verb > 0) std::cout << " + Initializing motion estimator ..." << std::endl;
		
		motionEstimator.init(verb, s, fc, nr_omp_threads, debug, outPath,
							 &reference, &obsModel, &micrographHandler);
	}
	
	if (estimateParams)
	{
		if (verb > 0) std::cout << " + Initializing motion parameter estimator ..." << std::endl;
		
		motionParamEstimator.init(
					verb, nr_omp_threads, debug, outPath,
					s, fc, chosenMdts, &motionEstimator, &reference, &obsModel);
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
	
	// The subsets will be used in openMPI parallelisation: instead of over g0->gc,
	// they will be over smaller subsets
	if (estimateMotion)
	{
		motionEstimator.process(motionMdts, 0, motionMdts.size()-1);
	}
	
	if (estimateOldBfacs)
	{
		bfactorEstimator.init(verb, s, fc, nr_omp_threads, outPath, debug,
							  &obsModel, &micrographHandler, &reference);
		
		bfactorEstimator.process(chosenMdts);
	}
	
	if (recombineFrames)
	{
		double k_out_A = obsModel.pixToAng(reference.k_out, s);
		
		frameRecombiner.init(
					allMdts, verb, s, fc, k_out_A, nr_omp_threads, outPath, debug,
					&obsModel, &micrographHandler);
		
		frameRecombiner.process(recombMdts, 0, recombMdts.size()-1);
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
	
	for (long g = 0; g < allMdts.size(); g++)
	{
		FileName fn_root = getOutputFileNameRoot(outPath, allMdts[g]);
		
		if (exists(fn_root+"_tracks.eps"))
		{
			fn_eps.push_back(fn_root+"_tracks.eps");
		}
		
		if (frameRecombiner.doingRecombination() && exists(fn_root+"_shiny.star"))
		{
			MetaDataTable mdt;
			mdt.read(fn_root+"_shiny.star");
			mdtAll.append(mdt);
		}
	}
	
	if (fn_eps.size() > 0)
	{
		joinMultipleEPSIntoSinglePDF(outPath + "logfile.pdf", fn_eps);
	}
	
	if (frameRecombiner.doingRecombination())
	{
		mdtAll.write(outPath+"shiny.star");
	}
	
	if (verb > 0)
	{
		std::cout << " + Done! " << std::endl;
		std::cout << " + Written logfile in " << outPath << "logfile.pdf" << std::endl;
		
		if (frameRecombiner.doingRecombination())
		{
			std::cout << " + Written new particle STAR file in "
					  << outPath << "shiny.star" << std::endl;
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


