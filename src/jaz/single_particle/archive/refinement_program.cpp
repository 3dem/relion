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

#include <src/jaz/single_particle/refinement_program.h>

#include <src/jaz/single_particle/img_proc/image_op.h>
#include <src/jaz/single_particle/refinement_helper.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/image_log.h>


RefinementProgram::RefinementProgram(bool singleReference, bool doesMovies)
:   singleReference(singleReference),
    optStar(false),
    noStar(false),
    optReference(false),
    noReference(false),
    noTilt(false),
    doesMovies(doesMovies),
    hasCorrMic(false),
    last_gainFn("")
{
}

int RefinementProgram::init(int argc, char *argv[])
{
    IOParser parser;

    try
    {
        parser.setCommandLine(argc, argv);

        parser.addSection("General options");

        starFn = parser.getOption("--i", "Input STAR file with a list of particles", optStar? "" : "NULL");

        if (!noReference)
        {
            if (singleReference)
            {
                reconFn0 = parser.getOption("--m", "Reference map", optReference? "" : "NULL");
            }
            else
            {
                reconFn0 = parser.getOption("--m1", "Reference map, half 1", optReference? "" : "NULL");
                reconFn1 = parser.getOption("--m2", "Reference map, half 2", optReference? "" : "NULL");
            }

            maskFn = parser.getOption("--mask", "Reference mask", "");
            fscFn = parser.getOption("--f", "Input STAR file with the FSC of the reference", "");
        }
        else
        {
            maskFn = "";
            fscFn = "";
        }

        outPath = parser.getOption("--out", "Output path");

        if (doesMovies)
        {            
            imgPath = parser.getOption("--mov", "Path to movies", "");
            corrMicFn = parser.getOption("--corr_mic", "List of uncorrected micrographs (e.g. corrected_micrographs.star)", "");
            preextracted = parser.checkOption("--preex", "Preextracted movie stacks");
            meta_path = parser.getOption("--meta", "Path to per-movie metadata star files", "");
            gain_path = parser.getOption("--gain_path", "Path to gain references", "");
            movie_ending = parser.getOption("--mov_end", "Ending of movie filenames", "");
            movie_toReplace = parser.getOption("--mov_toReplace", "Replace this string in micrograph names...", "");
            movie_replaceBy = parser.getOption("--mov_replaceBy", "..by this one", "");

            movie_angpix = textToFloat(parser.getOption("--mps", "Pixel size of input movies (Angst/pix)", "-1"));
            coords_angpix = textToFloat(parser.getOption("--cps", "Pixel size of particle coordinates in star-file (Angst/pix)", "-1"));

            hotCutoff = textToFloat(parser.getOption("--hot", "Clip hot pixels to this max. value (-1 = off, TIFF only)", "-1"));

            firstFrame = textToInteger(parser.getOption("--first_frame", "", "1")) - 1;
            lastFrame = textToInteger(parser.getOption("--last_frame", "", "-1")) - 1;

            saveMem = parser.checkOption("--sbs", "Load movies slice-by-slice to save memory (slower)");
        }
        else
        {
            imgPath = parser.getOption("--img", "Path to images", "");
        }

        angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (angst/pix) - read from STAR file by default", "0.0"));

        if (!noReference)
        {
            Cs = textToFloat(parser.getOption("--Cs", "Spherical aberration - read from STAR file by default", "-1"));
            kV = textToFloat(parser.getOption("--kV", "Electron energy (keV) - read from STAR file by default", "-1"));

            paddingFactor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));
        }
        else
        {
            Cs = -1;
            kV = -1;
            paddingFactor = 2;
        }

        if (noTilt)
        {
            beamtilt_x = 0.0;
            beamtilt_y = 0.0;

            applyTilt = false;

            beamtilt_xx = 1.0;
            beamtilt_xy = 0.0;
            beamtilt_yy = 1.0;
        }
        else
        {
            beamtilt_x = textToFloat(parser.getOption("--beamtilt_x", "Beamtilt in X-direction (in mrad)", "0."));
            beamtilt_y = textToFloat(parser.getOption("--beamtilt_y", "Beamtilt in Y-direction (in mrad)", "0."));
            applyTilt = ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.;

            beamtilt_xx = textToFloat(parser.getOption("--beamtilt_xx", "Anisotropic beamtilt, XX-coefficient", "1."));
            beamtilt_xy = textToFloat(parser.getOption("--beamtilt_xy", "Anisotropic beamtilt, XY-coefficient", "0."));
            beamtilt_yy = textToFloat(parser.getOption("--beamtilt_yy", "Anisotropic beamtilt, YY-coefficient", "1."));
        }

        anisoTilt = beamtilt_xx != 1.0 || beamtilt_xy != 0.0 || beamtilt_yy != 1.0;

        nr_omp_threads = textToInteger(parser.getOption("--jomp", "Number of OMP threads", "1"));
        minMG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
        maxMG = textToInteger(parser.getOption("--max_MG", "Last micrograph index", "-1"));

        debug = parser.checkOption("--debug", "Write debugging data");
        debugMov = parser.checkOption("--debug_mov", "Write debugging data for movie loading");

        int rco = readMoreOptions(parser, argc, argv);

        if (argc == 1)
        {
            parser.writeUsage(std::cerr);
            return 1;
        }

        if (parser.checkForErrors()) return 1;
        if (rco != 0) return rco;

        bool allGood = true;

        if (doesMovies && movie_angpix <= 0 && corrMicFn == "")
        {
            std::cerr << "Movie pixel size (--mps) is required unless a corrected_micrographs.star (--corr_mic) is provided.\n";
            allGood = false;
        }

        if (doesMovies && coords_angpix <= 0 && corrMicFn == "")
        {
            std::cerr << "Coordinates pixel size (--cps) is required unless a corrected_micrographs.star (--corr_mic) is provided.\n";
            allGood = false;
        }

        if (!allGood) return 12;
    }
    catch (RelionError XE)
    {
        parser.writeUsage(std::cout);
        std::cerr << XE;
        exit(1);
    }

    bool allGood = true;

    if (!noReference)
    {
        try
        {
            maps[0].read(reconFn0);
        }
        catch (RelionError XE)
        {
            std::cerr << "Unable to read map: " << reconFn0 << "\n";
            return 2;
        }

        if (!singleReference)
        {
            try
            {
                maps[1].read(reconFn1);
            }
            catch (RelionError XE)
            {
                std::cerr << "Unable to read map: " << reconFn1 << "\n";
                return 3;
            }
        }

        if (maps[0].data.xdim != maps[0].data.ydim || maps[0].data.ydim != maps[0].data.zdim)
        {
            REPORT_ERROR(reconFn0 + " is not cubical.\n");
        }

        if (!singleReference)
        {
            if (maps[1].data.xdim != maps[1].data.ydim || maps[1].data.ydim != maps[1].data.zdim)
            {
                REPORT_ERROR(reconFn1 + " is not cubical.\n");
            }

            if (   maps[0].data.xdim != maps[1].data.xdim
                || maps[0].data.ydim != maps[1].data.ydim
                || maps[0].data.zdim != maps[1].data.zdim)
            {
                REPORT_ERROR(reconFn0 + " and " + reconFn1 + " are of unequal size.\n");
            }
        }

        if (maskFn != "")
        {
            if (singleReference) std::cout << "masking reference...\n";
            else std::cout << "masking references...\n";

            Image<RFLOAT> mask, maskedRef;

            try
            {
                mask.read(maskFn);
            }
            catch (RelionError XE)
            {
                std::cout << "Unable to read mask: " << maskFn << "\n";
                return 4;
            }

            mask.read(maskFn);

            ImageOp::multiply(mask, maps[0], maskedRef);
            maps[0] = maskedRef;

            if (!singleReference)
            {
                ImageOp::multiply(mask, maps[1], maskedRef);
                maps[1] = maskedRef;
            }
        }

        s = maps[0].data.xdim;
        sh = s/2 + 1;


        if (singleReference) std::cout << "transforming reference...\n";
        else std::cout << "transforming references...\n";

        projectors[0] = Projector(s, TRILINEAR, paddingFactor, 10, 2);
        projectors[0].computeFourierTransformMap(maps[0].data, powSpec[0].data, maps[0].data.xdim);

        if (!singleReference)
        {
            projectors[1] = Projector(s, TRILINEAR, paddingFactor, 10, 2);
            projectors[1].computeFourierTransformMap(maps[1].data, powSpec[1].data, maps[1].data.xdim);
        }
    }

    useFsc = fscFn != "";
    MetaDataTable fscMdt;

    if (useFsc)
    {
        fscMdt.read(fscFn, "fsc");

        if (!fscMdt.containsLabel(EMDL_SPECTRAL_IDX))
        {
            std::cerr << fscFn << " does not contain a value for " << EMDL::label2Str(EMDL_SPECTRAL_IDX) << ".\n";
            allGood = false;
        }
        if (!fscMdt.containsLabel(EMDL_POSTPROCESS_FSC_TRUE))
        {
            std::cerr << fscFn << " does not contain a value for " << EMDL::label2Str(EMDL_POSTPROCESS_FSC_TRUE) << ".\n";
            allGood = false;
        }
    }

    if (!allGood)
    {
        return 1;
    }

    if (!noStar)
    {
        std::cout << "reading " << starFn << "...\n";

        mdt0.read(starFn);

        if (Cs < 0.0)
        {
            mdt0.getValue(EMDL_CTF_CS, Cs, 0);
            std::cout << " + Using spherical aberration from the input STAR file: " << Cs << "\n";
        }
        else
        {
            setForAll(EMDL_CTF_CS, Cs);
        }

        if (kV < 0.0)
        {
            mdt0.getValue(EMDL_CTF_VOLTAGE, kV, 0);
            std::cout << " + Using voltage from the input STAR file: " << kV << " kV\n";
        }
        else
        {
            setForAll(EMDL_CTF_VOLTAGE, kV);
        }

        if (angpix <= 0.0)
        {
            RFLOAT mag, dstep;
            mdt0.getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
            mdt0.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
            angpix = 10000 * dstep / mag;

            std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << "\n";
        }

        if (doesMovies)
        {
            if (movie_toReplace != "")
            {
                std::string name;

                for (int i = 0; i < mdt0.numberOfObjects(); i++)
                {
                    mdt0.getValue(EMDL_MICROGRAPH_NAME, name, i);

                    if (i == 0) std::cout << name << " -> ";

                    std::string::size_type pos0 = name.find(movie_toReplace);

                    if (pos0 != std::string::npos)
                    {
                        std::string::size_type pos1 = pos0 + movie_toReplace.length();

                        std::string before = name.substr(0, pos0);
                        std::string after = pos1 < name.length()? name.substr(pos1) : "";

                        name = before + movie_replaceBy + after;
                    }

                    if (i == 0) std::cout << name << "\n";

                    mdt0.setValue(EMDL_MICROGRAPH_NAME, name, i);
                }
            }
        }

        mdts = StackHelper::splitByStack(&mdt0);

        gc = maxMG >= 0? maxMG : mdts.size()-1;
        g0 = minMG;

        std::cout << "mg range: " << g0 << ".." << gc << "\n";
    }

    obsModel = LegacyObservationModel(angpix, Cs, kV * 1e3);

    if (applyTilt && anisoTilt)
	{
		obsModel.setAnisoTilt(beamtilt_xx, beamtilt_xy, beamtilt_yy);
	}

    int rc0 = _init();

    if (useFsc)
    {
        RefinementHelper::drawFSC(&fscMdt, freqWeight1D, freqWeight);
    }
    else if (!noReference)
    {
        freqWeight1D = std::vector<double>(sh, 1.0);
        freqWeight = Image<RFLOAT>(sh,s);
        freqWeight.data.initConstant(1.0);
    }

    if (doesMovies && corrMicFn != "")
    {
        MetaDataTable corrMic;
        corrMic.read(corrMicFn);

        mic2meta.clear();

        std::string micName, metaName;


        for (int i = 0; i < corrMic.numberOfObjects(); i++)
        {
            corrMic.getValueToString(EMDL_MICROGRAPH_NAME, micName, i);
            corrMic.getValueToString(EMDL_MICROGRAPH_METADATA_NAME, metaName, i);

            mic2meta[micName] = metaName;
        }

        hasCorrMic = true;
    }

    return rc0;
}

int RefinementProgram::run()
{
    return _run();
}

double RefinementProgram::angstToPixFreq(double a)
{
    return 2.0*sh*angpix/a;
}

double RefinementProgram::pixToAngstFreq(double p)
{
    return 2.0*sh*angpix/p;
}

void RefinementProgram::loadInitialMovieValues()
{
    if (preextracted)
    {
        std::string name, fullName, movieName;
        mdts[0].getValue(EMDL_IMAGE_NAME, fullName, 0);
        mdts[0].getValue(EMDL_MICROGRAPH_NAME, movieName, 0);
        name = fullName.substr(fullName.find("@")+1);

        std::string finName;

        if (imgPath == "")
        {
            finName = name;
        }
        else
        {
            finName = imgPath + "/" + movieName.substr(movieName.find_last_of("/")+1);
        }

        Image<RFLOAT> stack0;
        stack0.read(finName, false);

        const int pc0 = mdts[0].numberOfObjects();
        const bool zstack = stack0.data.zdim > 1;
        const int stackSize = zstack? stack0.data.zdim : stack0.data.ndim;

        if (lastFrame < 0) fc = stackSize / pc0 - firstFrame;
        else fc = lastFrame - firstFrame + 1;
    }
    else
    {
        if (hasCorrMic)
        {
            std::string mgFn;
            mdts[0].getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);

            std::string metaFn = mic2meta[mgFn];

            if (meta_path != "")
            {
                metaFn = meta_path + "/" + metaFn.substr(metaFn.find_last_of("/")+1);
            }

            micrograph = Micrograph(metaFn);

            if (movie_angpix <= 0)
            {
                movie_angpix = micrograph.angpix;
                std::cout << " + Using movie pixel size from " << metaFn << ": " << movie_angpix << " A\n";
            }
            else
            {
                std::cout << " + Using movie pixel size from command line: " << movie_angpix << " A\n";
            }

            if (coords_angpix <= 0)
            {
                coords_angpix = micrograph.angpix * micrograph.getBinningFactor();
                std::cout << " + Using coord. pixel size from " << metaFn << ": " << coords_angpix << " A\n";
            }
            else
            {
                std::cout << " + Using coord. pixel size from command line: " << coords_angpix << " A\n";
            }

            if (lastFrame < 0) fc = micrograph.getNframes() - firstFrame;
            else fc = lastFrame - firstFrame + 1;
        }
        else
        {
            REPORT_ERROR("You can no longer use this program without micrograph metadata STAR files.");
        }
    }
}

std::vector<std::vector<Image<Complex>>> RefinementProgram::loadMovie(
        int g, int pc, std::vector<ParFourierTransformer>& fts)
{
    std::vector<std::vector<Image<Complex>>> movie;

    if (preextracted)
    {
        movie = StackHelper::loadMovieStackFS(
            &mdts[g], imgPath, false, nr_omp_threads, &fts,
            firstFrame, lastFrame);
    }
    else
    {
        std::string mgFn;
        mdts[g].getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);

        if (hasCorrMic)
        {
            std::string metaFn = mic2meta[mgFn];

            if (meta_path != "")
            {
                metaFn = meta_path + "/" + metaFn.substr(metaFn.find_last_of("/")+1);
            }

            micrograph = Micrograph(metaFn);

            std::string mgFn = micrograph.getMovieFilename();
            std::string gainFn = micrograph.getGainFilename();

            if (movie_ending != "")
            {
                mgFn.substr(0, mgFn.find_last_of(".")+1) + movie_ending;
            }

            if (imgPath != "")
            {
                mgFn = imgPath + "/" + mgFn.substr(mgFn.find_last_of("/")+1);
            }

            bool mgHasGain = false;

            if (gainFn != "")
            {
                if (gain_path != "")
                {
                    gainFn = gain_path + "/" + gainFn.substr(gainFn.find_last_of("/")+1);
                }

                if (gainFn != last_gainFn)
                {
                    lastGainRef.read(gainFn);
                    last_gainFn = gainFn;
                }

                mgHasGain = true;
            }

            MultidimArray<bool> defectMask;

            bool hasDefect = (micrograph.fnDefect != "" || micrograph.hotpixelX.size() != 0);
            if (hasDefect)
                micrograph.fillDefectAndHotpixels(defectMask);

            movie = StackHelper::extractMovieStackFS(
                &mdts[g], mgHasGain? &lastGainRef : 0, hasDefect ? &defectMask : 0,
                mgFn, angpix, coords_angpix, movie_angpix, s,
                nr_omp_threads, true, firstFrame, lastFrame, hotCutoff, debugMov, saveMem);
        }
        else
        {
            REPORT_ERROR("You can no longer use this program without micrograph metadata STAR files.");
        }

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (int p = 0; p < pc; p++)
        {
            StackHelper::varianceNormalize(movie[p], false);
        }
    }

    if (angpix < coords_angpix)
    {
        std::cerr << "WARNING: pixel size (--angpix) is greater than the AutoPick pixel size (--coords_angpix)\n";

        if (coords_angpix < angpix + 0.01)
        {
            std::cerr << "        This is probably a rounding error. It is recommended to set --angpix ("
                      << angpix << ") to at least " << coords_angpix << "\n";

        }
    }

    if (angpix < movie_angpix)
    {
        std::cerr << "WARNING: pixel size (--angpix) is greater than the movie pixel size (--movie_angpix)\n";

        if (movie_angpix < angpix + 0.01)
        {
            std::cerr << "        This is probably a rounding error. It is recommended to set --angpix ("
                      << angpix << ") to at least " << movie_angpix << "\n";

        }
    }

    return movie;
}

void RefinementProgram::setForAll(EMDLabel label, RFLOAT value)
{
    for (int i = 0; i < mdt0.numberOfObjects(); i++)
    {
        mdt0.setValue(label, value, i);
    }
}

std::string RefinementProgram::getMicrographTag(int m)
{
    std::string tag;
    mdts[m].getValue(EMDL_IMAGE_NAME, tag, 0);
    tag = tag.substr(0,tag.find_last_of('.'));
    tag = tag.substr(tag.find_first_of('@')+1);

    return tag;
}
