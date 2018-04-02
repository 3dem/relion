#include "micrograph_handler.h"
#include <src/jaz/stack_helper.h>

using namespace gravis;

MicrographHandler::MicrographHandler()
:   hasCorrMic(false),
    nr_omp_threads(1),
    firstFrame(0),
    lastFrame(-1),
    hotCutoff(-1),
    preextracted(false),
    debug(false),
    saveMem(false),
    ready(false),
    movie_path(""),
    meta_path(""),
    movie_ending(""),
    gain_path(""),
    last_gainFn("")
{}

void MicrographHandler::init(std::string corrMicFn)
{
    MetaDataTable corrMic;
    corrMic.read(corrMicFn);

    mic2meta.clear();

    std::string micName, metaName;

    for (int i = 0; i < corrMic.numberOfObjects(); i++)
    {
        corrMic.getValueToString(EMDL_MICROGRAPH_NAME, micName, i);
        corrMic.getValueToString(EMDL_MICROGRAPH_METADATA_NAME, metaName, i);

        // remove the pipeline job prefix
        FileName fn_pre, fn_jobnr, fn_post;
        decomposePipelineFileName(micName, fn_pre, fn_jobnr, fn_post);

        mic2meta[fn_post] = metaName;
    }

    hasCorrMic = true;
    ready = true;
}

void MicrographHandler::init()
{
    hasCorrMic = false;
    ready = true;
}

void MicrographHandler::loadInitial(
        const MetaDataTable& mdt, double angpix, bool verb,
        int& fc, int& micrograph_xsize, int& micrograph_ysize,
        double& movie_angpix, double& coords_angpix,
        double& dosePerFrame)
{
    if (!ready)
    {
        REPORT_ERROR("ERROR: MicrographHandler::loadInitial - MicrographHandler not initialized.");
    }

    if (preextracted)
    {
        if (lastFrame < 0)
        {
             std::string name;
             mdt.getValue(EMDL_MICROGRAPH_NAME, name, 0);

             Image<RFLOAT> stack0;
             stack0.read(name, false);

             const int pc0 = mdt.numberOfObjects();
             const bool zstack = stack0.data.zdim > 1;
             const int stackSize = zstack? stack0.data.zdim : stack0.data.ndim;
             fc = stackSize / pc0 - firstFrame;
        }
        else
        {
            fc = lastFrame - firstFrame + 1;
        }

        FileName fn_mic;
        mdt.getValue(EMDL_MICROGRAPH_NAME, fn_mic, 0);

        Image<RFLOAT> dum;
        dum.read(fn_mic, false);
        micrograph_xsize = XSIZE(dum());
        micrograph_ysize = YSIZE(dum());
    }
    else
    {
        if (hasCorrMic)
        {
            std::string mgFn;
            mdt.getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);

            // remove the pipeline job prefix
            FileName fn_pre, fn_jobnr, fn_post;
            decomposePipelineFileName(mgFn, fn_pre, fn_jobnr, fn_post);

            //@TODO: make safe:
            std::string metaFn = mic2meta[fn_post];

            if (meta_path != "")
            {
                metaFn = meta_path + "/" + metaFn.substr(metaFn.find_last_of("/")+1);
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

            // this is safe to do - motionEstimator.read() has been called already
            if (dosePerFrame < 0)
            {
                dosePerFrame = micrograph.dose_per_frame;

                if (verb > 0)
                {
                    std::cout << " + Using dose per frame from " << metaFn << ": "
                              << dosePerFrame << " A\n";
                }
            }

            micrograph_xsize = micrograph.getWidth();
            micrograph_ysize = micrograph.getHeight();

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
            if (lastFrame < 0)
            {
                FileName mgFn;
                mdt.getValue(EMDL_MICROGRAPH_NAME, mgFn, 0);

                Image<RFLOAT> dum;
                dum.read(mgFn, false);
                micrograph_xsize = XSIZE(dum());
                micrograph_ysize = YSIZE(dum());

                fc = (dum().zdim > 1? dum().zdim : dum().ndim) - firstFrame;
            }
            else
            {
                fc = lastFrame - firstFrame + 1;
            }
        }
    }

    if (angpix < coords_angpix)
    {
        std::cerr << "WARING: pixel size (--angpix) is greater than the AutoPick pixel size (--coords_angpix)\n";

        if (coords_angpix < angpix + 0.01)
        {
            std::cerr << "        This is probably a rounding error. It is recommended to set --angpix ("
                      << angpix << ") to at least " << coords_angpix << "\n";

        }
    }

    if (angpix < movie_angpix)
    {
        std::cerr << "WARING: pixel size (--angpix) is greater than the movie pixel size (--movie_angpix)\n";

        if (movie_angpix < angpix + 0.01)
        {
            std::cerr << "        This is probably a rounding error. It is recommended to set --angpix ("
                      << angpix << ") to at least " << movie_angpix << "\n";

        }
    }
}

std::vector<std::vector<Image<Complex>>> MicrographHandler::loadMovie(
    const MetaDataTable &mdt, int s,
    double angpix, double movie_angpix, double coords_angpix,
    std::vector<ParFourierTransformer>& fts)
{
    if (!ready)
    {
        REPORT_ERROR("ERROR: MicrographHandler::loadMovie - MicrographHandler not initialized.");
    }

    std::vector<std::vector<Image<Complex>>> movie;

    const int nr_omp_threads = fts.size();

    if (preextracted)
    {
        movie = StackHelper::loadMovieStackFS(
            &mdt, "", false, nr_omp_threads, &fts,
            firstFrame, lastFrame);
    }
    else
    {
        std::string mgFn;
        mdt.getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);
        FileName fn_pre, fn_jobnr, fn_post;
        decomposePipelineFileName(mgFn, fn_pre, fn_jobnr, fn_post);

        if (hasCorrMic)
        {
            //@TODO: make safe:
            std::string metaFn = mic2meta[fn_post];

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

            if (movie_path != "")
            {
                mgFn = movie_path + "/" + mgFn.substr(mgFn.find_last_of("/")+1);
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

            movie = StackHelper::extractMovieStackFS(
                &mdt, mgHasGain? &lastGainRef : 0,
                mgFn, angpix, coords_angpix, movie_angpix, s,
                nr_omp_threads, true, firstFrame, lastFrame,
                hotCutoff, debug, saveMem);
        }
        else
        {
            // gain no longer supported without a corrected_micrographs.star:
            /*movie = StackHelper::extractMovieStackFS(
                &mdts[g], meta_path, movie_path, movie_ending,
                angpix, coords_angpix, movie_angpix, s,
                nr_omp_threads, true, firstFrame, lastFrame, hotCutoff, debugMov);*/

            movie = StackHelper::extractMovieStackFS(
                &mdt, 0,
                mgFn, angpix, coords_angpix, movie_angpix, s,
                nr_omp_threads, true, firstFrame, lastFrame,
                hotCutoff, debug, saveMem);
        }

        const int pc = movie.size();

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (int p = 0; p < pc; p++)
        {
            StackHelper::varianceNormalize(movie[p], false);
        }
    }

    return movie;
}

std::vector<std::vector<Image<Complex>>> MicrographHandler::loadMovie(
        const MetaDataTable &mdt, int s,
        double angpix, double movie_angpix, double coords_angpix,
        std::vector<ParFourierTransformer>& fts,
        const std::vector<d2Vector>& pos,
        std::vector<std::vector<d2Vector>>& tracks,
        bool unregGlob, std::vector<d2Vector>& globComp)
{
    std::vector<std::vector<Image<Complex>>> out = loadMovie(
        mdt, s, angpix, movie_angpix, coords_angpix, fts);

    if (!hasCorrMic || micrograph.model == 0)
    {
        tracks.resize(0);
    }
    else
    {
        const int fc = micrograph.getNframes();
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
                micrograph.getShiftAt(f+1, 0, 0, sx, sy, false);

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

                micrograph.getShiftAt(f+1, in.x, in.y, sx, sy, true);

                tracks[p][f] = -outputScale * d2Vector(sx,sy) - globComp[f];
            }
        }
    }

    return out;
}
