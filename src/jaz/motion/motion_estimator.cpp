#include "motion_estimator.h"
#include "motion_refiner.h"
#include "motion_helper.h"
#include "gp_motion_fit.h"

#include <src/jaz/damage_helper.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/fsc_helper.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/image_op.h>
#include <src/jaz/image_log.h>
#include <src/jaz/optimization/lbfgs.h>

using namespace gravis;

MotionEstimator::MotionEstimator(MotionRefiner& motionRefiner)
:   motionRefiner(motionRefiner)
{
}

void MotionEstimator::read(IOParser& parser, int argc, char *argv[])
{
    parser.addSection("Motion fit options (basic)");

    dosePerFrame = textToFloat(parser.getOption("--fdose", "Electron dose per frame (in e^-/A^2)", "-1"));
    sig_vel = textToFloat(parser.getOption("--s_vel", "Velocity sigma [Angst/dose]", "0.5"));
    sig_div = textToFloat(parser.getOption("--s_div", "Divergence sigma [Angst]", "1000.0"));
    sig_acc = textToFloat(parser.getOption("--s_acc", "Acceleration sigma [Angst/dose]", "-1.0"));
    k_cutoff = textToFloat(parser.getOption("--k_cut", "Freq. cutoff for parameter estimation [Pixels]", "-1.0"));
    k_cutoff_Angst = textToFloat(parser.getOption("--k_cut_A", "Freq. cutoff for parameter estimation [Angstrom]", "-1.0"));
    diag = parser.checkOption("--diag", "Write out diagnostic data");

    parser.addSection("Motion fit options (advanced)");

    dmga = textToFloat(parser.getOption("--dmg_a", "Damage model, parameter a", " 3.40"));
    dmgb = textToFloat(parser.getOption("--dmg_b", "                        b", "-1.06"));
    dmgc = textToFloat(parser.getOption("--dmg_c", "                        c", "-0.54"));

    maxIters = textToInteger(parser.getOption("--max_iters", "Maximum number of iterations", "10000"));
    optEps = textToFloat(parser.getOption("--eps", "Terminate optimization after gradient length falls below this value", "1e-4"));

    unregGlob = parser.checkOption("--unreg_glob", "Do not regularize global component of motion");
    noGlobOff = parser.checkOption("--no_glob_off", "Do not compute initial per-particle offsets");
    debugOpt = parser.checkOption("--debug_opt", "Write optimization debugging info");

    global_init = parser.checkOption("--gi", "Initialize with global trajectories instead of loading them from metadata file");
    expKer = parser.checkOption("--exp_k", "Use exponential kernel instead of sq. exponential");
    maxEDs = textToInteger(parser.getOption("--max_ed", "Maximum number of eigendeformations", "-1"));
}

void MotionEstimator::init()
{
    if (!global_init && motionRefiner.corrMicFn == "")
    {
        if (motionRefiner.verb > 0)
        {
            std::cerr << "\nWarning: in the absence of a corrected_micrographs.star file (--corr_mic), global paths are used for initialization.\n\n";
        }

        global_init = true;
    }

    if (k_cutoff_Angst > 0.0 && k_cutoff > 0.0)
    {
        REPORT_ERROR("ERROR: Cutoff frequency can only be provided in pixels (--k_cut) or Angstrom (--k_cut_A), not both.");
    }

    s = motionRefiner.s;
    sh = motionRefiner.sh;
    fc = motionRefiner.fc;

    if (k_cutoff_Angst > 0.0 && k_cutoff < 0.0)
    {
        k_cutoff = motionRefiner.angToPix(k_cutoff_Angst);
    }

    else if (k_cutoff > 0.0 && k_cutoff_Angst < 0.0)
    {
        k_cutoff_Angst = motionRefiner.pixToAng(k_cutoff);
    }

    dmgWeight = DamageHelper::damageWeights(s, motionRefiner.angpix,
            motionRefiner.firstFrame, fc, dosePerFrame, dmga, dmgb, dmgc);

    dmgWeightEval.resize(fc);

    for (int f = 0; f < fc; f++)
    {
        dmgWeight[f].data.xinit = 0;
        dmgWeight[f].data.yinit = 0;

        ImageOp::multiplyBy(dmgWeight[f], motionRefiner.freqWeight);

        // evaluation is performed beyond k_cutoff only
        // do not crop dmgWeightEval
        dmgWeightEval[f] = dmgWeight[f];

        if (k_cutoff > 0.0)
        {
            std::stringstream stsf;
            stsf << f;
            dmgWeight[f] = FilterHelper::ButterworthEnvFreq2D(
                        dmgWeight[f], k_cutoff-1, k_cutoff+1);
        }
    }
}

void MotionEstimator::process(long g_start, long g_end)
{
    int barstep;
    int my_nr_micrographs = g_end - g_start + 1;

    if (motionRefiner.verb > 0)
    {
        std::cout << " + Performing loop over all micrographs ... " << std::endl;
        init_progress_bar(my_nr_micrographs);
        barstep = XMIPP_MAX(1, my_nr_micrographs/ 60);
    }

    std::vector<ParFourierTransformer> fts(motionRefiner.nr_omp_threads);

    std::vector<Image<RFLOAT>>
            tables(motionRefiner.nr_omp_threads),
            weights0(motionRefiner.nr_omp_threads),
            weights1(motionRefiner.nr_omp_threads);

    for (int i = 0; i < motionRefiner.nr_omp_threads; i++)
    {
        FscHelper::initFscTable(sh, fc, tables[i], weights0[i], weights1[i]);
    }

    const double sig_vel_nrm = dosePerFrame * sig_vel / motionRefiner.angpix;
    const double sig_acc_nrm = dosePerFrame * sig_acc / motionRefiner.angpix;
    const double sig_div_nrm = dosePerFrame * sig_div / motionRefiner.coords_angpix;

    int pctot = 0;

    long nr_done = 0;
    FileName prevdir = "";

    for (long g = g_start; g <= g_end; g++)
    {
        const int pc = motionRefiner.mdts[g].numberOfObjects();
        if (pc == 0) continue;

        // Make sure output directory exists
        FileName newdir = motionRefiner.getOutputFileNameRoot(g);
        newdir = newdir.beforeLastOf("/");

        if (newdir != prevdir)
        {
            std::string command = " mkdir -p " + newdir;
            int res = system(command.c_str());
        }

        std::vector<std::vector<Image<Complex>>> movie;
        std::vector<std::vector<Image<RFLOAT>>> movieCC;
        std::vector<d2Vector> positions;
        std::vector<std::vector<d2Vector>> initialTracks;
        std::vector<d2Vector> globComp;

        prepMicrograph(g, fts, dmgWeight,
                movie, movieCC, positions, initialTracks, globComp);

        pctot += pc;

        std::vector<std::vector<gravis::d2Vector>> tracks = optimize(
                movieCC, initialTracks,
                sig_vel_nrm, sig_acc_nrm, sig_div_nrm,
                positions, globComp);

        updateFCC(movie, tracks, motionRefiner.mdts[g], tables, weights0, weights1);

        writeOutput(tracks, tables, weights0, weights1, positions, motionRefiner.outPath, g, 30.0);

        for (int i = 0; i < motionRefiner.nr_omp_threads; i++)
        {
            tables[i].data.initZeros();
            weights0[i].data.initZeros();
            weights1[i].data.initZeros();
        }

        nr_done++;

        if (motionRefiner.verb > 0 && nr_done % barstep == 0)
        {
            progress_bar(nr_done);
        }
    }

    if (motionRefiner.verb > 0)
    {
        progress_bar(my_nr_micrographs);
    }
}


void MotionEstimator::prepMicrograph(
        long g, std::vector<ParFourierTransformer>& fts,
        const std::vector<Image<RFLOAT>>& dmgWeight,
        std::vector<std::vector<Image<Complex>>>& movie,
        std::vector<std::vector<Image<RFLOAT>>>& movieCC,
        std::vector<d2Vector>& positions,
        std::vector<std::vector<d2Vector>>& initialTracks,
        std::vector<d2Vector>& globComp)
{
    const int pc = motionRefiner.mdts[g].numberOfObjects();

    movie = motionRefiner.loadMovie(g, pc, fts); // throws exceptions
    std::vector<double> sigma2 = StackHelper::powerSpectrum(movie);

    #pragma omp parallel for num_threads(motionRefiner.nr_omp_threads)
    for (int p = 0; p < pc; p++)
    for (int f = 0; f < fc; f++)
    {
        MotionHelper::noiseNormalize(movie[p][f], sigma2, movie[p][f]);
    }

    positions = std::vector<gravis::d2Vector>(pc);

    for (int p = 0; p < pc; p++)
    {
        motionRefiner.mdts[g].getValue(EMDL_IMAGE_COORD_X, positions[p].x, p);
        motionRefiner.mdts[g].getValue(EMDL_IMAGE_COORD_Y, positions[p].y, p);
    }

    movieCC = MotionHelper::movieCC(
            motionRefiner.projectors[0], motionRefiner.projectors[1],
            motionRefiner.obsModel, motionRefiner.mdts[g], movie,
            sigma2, dmgWeight, fts, motionRefiner.nr_omp_threads);

    initialTracks = std::vector<std::vector<d2Vector>>(pc);

    if (global_init)
    {
        std::vector<Image<RFLOAT>> ccSum = MotionHelper::addCCs(movieCC);
        std::vector<gravis::d2Vector> globTrack = MotionHelper::getGlobalTrack(ccSum);
        std::vector<gravis::d2Vector> globOffsets;

        if (noGlobOff)
        {
            globOffsets = std::vector<d2Vector>(pc, d2Vector(0,0));
        }
        else
        {
            globOffsets = MotionHelper::getGlobalOffsets(
                    movieCC, globTrack, 0.25*s, motionRefiner.nr_omp_threads);
        }

        if (diag)
        {
            ImageLog::write(ccSum, motionRefiner.getOutputFileNameRoot(g) + "_CCsum", CenterXY);
        }

        for (int p = 0; p < pc; p++)
        {
            initialTracks[p] = std::vector<d2Vector>(fc);

            for (int f = 0; f < fc; f++)
            {
                if (unregGlob)
                {
                    initialTracks[p][f] = globOffsets[p];
                }
                else
                {
                    initialTracks[p][f] = globTrack[f] + globOffsets[p];
                }
            }
        }

        globComp = unregGlob? globTrack : std::vector<d2Vector>(fc, d2Vector(0,0));
    }
    else
    {
        const d2Vector inputScale(
                motionRefiner.coords_angpix
                    / (motionRefiner.movie_angpix * motionRefiner.micrograph.getWidth()),
                motionRefiner.coords_angpix
                    / (motionRefiner.movie_angpix * motionRefiner.micrograph.getHeight()));

        const double outputScale = motionRefiner.movie_angpix / motionRefiner.angpix;

        globComp = std::vector<d2Vector>(fc, d2Vector(0,0));

        if (unregGlob)
        {
            for (int f = 0; f < fc; f++)
            {
                RFLOAT sx, sy;
                motionRefiner.micrograph.getShiftAt(f+1, 0, 0, sx, sy, false);

                globComp[f] = -outputScale * d2Vector(sx, sy);
            }
        }

        for (int p = 0; p < pc; p++)
        {
            initialTracks[p] = std::vector<d2Vector>(fc);

            for (int f = 0; f < fc; f++)
            {
                d2Vector in(inputScale.x * positions[p].x - 0.5,
                            inputScale.y * positions[p].y - 0.5);

                RFLOAT sx, sy;

                motionRefiner.micrograph.getShiftAt(f+1, in.x, in.y, sx, sy, true);

                initialTracks[p][f] = -outputScale * d2Vector(sx,sy) - globComp[f];
            }
        }
    }
}

std::vector<std::vector<d2Vector>> MotionEstimator::optimize(
    const std::vector<std::vector<Image<RFLOAT>>>& movieCC,
    const std::vector<std::vector<gravis::d2Vector>>& inTracks,
    double sig_vel_px, double sig_acc_px, double sig_div_px,
    const std::vector<gravis::d2Vector>& positions,
    const std::vector<gravis::d2Vector>& globComp)
{
    const double eps = 1e-20;

    if (sig_vel_px < eps)
    {
        sig_vel_px = eps;
    }

    if (sig_div_px < eps)
    {
        sig_div_px = eps;
    }

    const int pc = inTracks.size();

    if (pc == 0) return std::vector<std::vector<d2Vector>>(0);

    const int fc = inTracks[0].size();

    GpMotionFit gpmf(movieCC, sig_vel_px, sig_div_px, sig_acc_px,
                     maxEDs, positions,
                     globComp, motionRefiner.nr_omp_threads, expKer);

    std::vector<double> initialCoeffs;

    gpmf.posToParams(inTracks, initialCoeffs);

    std::vector<double> optCoeffs = LBFGS::optimize(
        initialCoeffs, gpmf, debugOpt, maxIters, optEps);

    std::vector<std::vector<d2Vector>> out(pc, std::vector<d2Vector>(fc));
    gpmf.paramsToPos(optCoeffs, out);

    return out;
}

void MotionEstimator::updateFCC(
        const std::vector<std::vector<Image<Complex>>>& movie,
        const std::vector<std::vector<d2Vector>>& tracks,
        const MetaDataTable& mdt,
        std::vector<Image<RFLOAT>>& tables,
        std::vector<Image<RFLOAT>>& weights0,
        std::vector<Image<RFLOAT>>& weights1)
{
    const int pc = mdt.numberOfObjects();

    #pragma omp parallel for num_threads(motionRefiner.nr_omp_threads)
    for (int p = 0; p < pc; p++)
    {
        int threadnum = omp_get_thread_num();

        Image<Complex> pred;
        std::vector<Image<Complex>> obs = movie[p];

        for (int f = 0; f < fc; f++)
        {
            shiftImageInFourierTransform(obs[f](), obs[f](), s, -tracks[p][f].x, -tracks[p][f].y);
        }

        int randSubset;
        mdt.getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
        randSubset -= 1;

        if (randSubset == 0)
        {
            pred = motionRefiner.obsModel.predictObservation(
                        motionRefiner.projectors[1], mdt, p, true, true);
        }
        else
        {
            pred = motionRefiner.obsModel.predictObservation(
                        motionRefiner.projectors[0], mdt, p, true, true);
        }

        FscHelper::updateFscTable(obs, pred, tables[threadnum],
                                  weights0[threadnum], weights1[threadnum]);
    }
}

void MotionEstimator::writeOutput(
        const std::vector<std::vector<d2Vector>>& tracks,
        const std::vector<Image<RFLOAT>>& fccData,
        const std::vector<Image<RFLOAT>>& fccWeight0,
        const std::vector<Image<RFLOAT>>& fccWeight1,
        const std::vector<d2Vector>& positions,
        std::string outPath, long g,
        double visScale)
{
    const int pc = tracks.size();

    if (pc == 0) return;

    const int fc = tracks[0].size();

    FileName fn_root = motionRefiner.getOutputFileNameRoot(g);
    MotionHelper::writeTracks(tracks, fn_root + "_tracks.star");

    Image<RFLOAT> fccDataSum(sh,fc), fccWeight0Sum(sh,fc), fccWeight1Sum(sh,fc);
    fccDataSum.data.initZeros();
    fccWeight0Sum.data.initZeros();
    fccWeight1Sum.data.initZeros();

    for (int i = 0; i < fccData.size(); i++)
    {
        for (int y = 0; y < fc; y++)
        for (int x = 0; x < sh; x++)
        {
            fccDataSum(y,x) += fccData[i](y,x);
            fccWeight0Sum(y,x) += fccWeight0[i](y,x);
            fccWeight1Sum(y,x) += fccWeight1[i](y,x);
        }
    }

    fccDataSum.write(fn_root + "_FCC_cc.mrc");
    fccWeight0Sum.write(fn_root + "_FCC_w0.mrc");
    fccWeight1Sum.write(fn_root + "_FCC_w1.mrc");

    // plot EPS graph with all observed and fitted tracks
    std::vector<std::vector<gravis::d2Vector>> visTracks(pc);

    for (int p = 0; p < pc; p++)
    {
        visTracks[p] = std::vector<gravis::d2Vector>(fc);
    }

    std::vector<gravis::d2Vector> globalTrack(fc);

    for (int f = 0; f < fc; f++)
    {
        globalTrack[f] = d2Vector(0,0);

        for (int p = 0; p < pc; p++)
        {
            globalTrack[f] += tracks[p][f];
        }

        globalTrack[f] /= pc;
        for (int p = 0; p < pc; p++)
        {
            visTracks[p][f] = positions[p] + visScale * tracks[p][f];
        }
    }

    // Make a postscript with the tracks
    FileName fn_eps = fn_root + "_tracks.eps";
    CPlot2D *plot2D=new CPlot2D(fn_eps);
    plot2D->SetXAxisSize(600);
    plot2D->SetYAxisSize(600);
    plot2D->SetDrawLegend(false);

    // Global track in the middle
    CDataSet dataSet;
    dataSet.SetDrawMarker(false);
    dataSet.SetDatasetColor(0.0,0.0,1.0);
    dataSet.SetLineWidth(1.);
    RFLOAT xshift, yshift;
    const RFLOAT xcenter =  motionRefiner.micrograph_xsize / 2.0;
    const RFLOAT ycenter =  motionRefiner.micrograph_ysize / 2.0;
    for (int f = 0; f < fc; f++)
    {
        CDataPoint point(xcenter + visScale * globalTrack[f].x, ycenter + visScale * globalTrack[f].y);
        dataSet.AddDataPoint(point);
    }
    plot2D->AddDataSet(dataSet);

    // Mark starting point global track
    CDataSet dataSetStart;
    dataSetStart.SetDrawMarker(true);
    dataSetStart.SetMarkerSize(2);
    dataSetStart.SetDatasetColor(1.0,0.0,0.0);
    CDataPoint point2(xcenter + visScale * globalTrack[0].x, ycenter + visScale * globalTrack[0].y);
    dataSetStart.AddDataPoint(point2);
    plot2D->AddDataSet(dataSetStart);

    // Now loop over all particles for local tracks
    for (int p = 0; p < pc; p++)
    {
        CDataSet fit;
        fit.SetDrawMarker(false);
        fit.SetDatasetColor(0.0,0.0,0.0);
        fit.SetLineWidth(0.5);

        for (int f = 0; f < fc; f++)
        {
            CDataPoint point(visTracks[p][f].x, visTracks[p][f].y);
            fit.AddDataPoint(point);
        }
        plot2D->AddDataSet(fit);

        // Mark start of each track
        CDataSet patch_start;
        patch_start.SetDrawMarker(true);
        patch_start.SetMarkerSize(2);
        patch_start.SetDatasetColor(1.0,0.3,0.0);
        CDataPoint point3(visTracks[p][0].x, visTracks[p][0].y);
        patch_start.AddDataPoint(point3);
        plot2D->AddDataSet(patch_start);
    }

    char title[256];
    snprintf(title, 255, "X (in pixels; trajectory scaled by %.0f)", visScale);
    plot2D->SetXAxisTitle(title);
    title[0] = 'Y';
    plot2D->SetYAxisTitle(title);

    plot2D->OutputPostScriptPlot(fn_eps);

    // Compatibility with Jasenko's diagnostic .dat files
    // NOTTODO: remove this
    // Don't! It's the only way to plot tracks on top of each other.

    if (!diag) return;

    std::ofstream rawOut(fn_root + "_tracks.dat");
    std::ofstream visOut(fn_root + "_visTracks.dat");
    std::ofstream visOut15(fn_root + "_visTracks_first15.dat");

    for (int p = 0; p < pc; p++)
    {
        rawOut << "#particle " << p << "\n";
        visOut << "#particle " << p << "\n";
        visOut15 << "#particle " << p << "\n";

        for (int f = 0; f < fc; f++)
        {
            rawOut << tracks[p][f].x << " " << tracks[p][f].y << "\n";
            visOut << visTracks[p][f].x << " " << visTracks[p][f].y << "\n";

            if (f < 15) visOut15 << visTracks[p][f].x << " " << visTracks[p][f].y << "\n";
        }

        rawOut << "\n";
        visOut << "\n";
        visOut15 << "\n";
    }

    std::ofstream glbOut(fn_root + "_globTrack.dat");

    for (int f = 0; f < fc; f++)
    {
        glbOut << globalTrack[f].x << " " << globalTrack[f].y << "\n";
    }
}

