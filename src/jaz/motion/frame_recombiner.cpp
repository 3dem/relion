#include "frame_recombiner.h"
#include "motion_refiner.h"
#include "motion_helper.h"
#include <src/jaz/stack_helper.h>
#include <src/jaz/vtk_helper.h>
#include <src/jaz/damage_helper.h>
#include <src/jaz/image_log.h>
#include <src/filename.h>

using namespace gravis;

FrameRecombiner::FrameRecombiner(MotionRefiner& motionRefiner)
:   motionRefiner(motionRefiner)
{
}

void FrameRecombiner::read(IOParser& parser, int argc, char* argv[])
{
    parser.addSection("Combine frames options");

    doCombineFrames = parser.checkOption("--combine_frames", "Combine movie frames into polished particles.");
    k0a = textToFloat(parser.getOption("--bfac_minfreq", "Min. frequency used in B-factor fit [Angst]", "20"));
    k1a = textToFloat(parser.getOption("--bfac_maxfreq", "Max. frequency used in B-factor fit [Angst]", "-1"));
    bfacFn = parser.getOption("--bfactors", "A .star file with external B/k-factors", "");
    bfac_diag = parser.checkOption("--diag_bfactor", "Write out B/k-factor diagnostic data");
}

void FrameRecombiner::init()
{
    // Split again, as a subset might have been done before for only_do_unfinished...
    motionRefiner.mdts.clear();
    motionRefiner.mdts = StackHelper::splitByMicrographName(&motionRefiner.mdt0);

    // check whether combine_frames output stack exist and if they do, then skip this micrograph
    // @TODO: turn off if not all motion had been estimated!
    if (motionRefiner.only_do_unfinished)
    {
        std::vector<MetaDataTable> unfinished_mdts;

        for (long int g = motionRefiner.minMG; g <= motionRefiner.maxMG; g++ )
        {
            bool is_done = true;

            if (!exists(motionRefiner.getOutputFileNameRoot(g)+"_shiny.mrcs"))
            {
                is_done = false;
            }

            if (!exists(motionRefiner.getOutputFileNameRoot(g)+"_shiny.star"))
            {
                is_done = false;
            }

            if (!is_done)
            {
                unfinished_mdts.push_back(motionRefiner.mdts[g]);
            }
        }

        motionRefiner.mdts = unfinished_mdts;

        if (motionRefiner.verb > 0)
        {
            if (motionRefiner.mdts.size() > 0)
            {
                std::cout << "   - Will only combine frames for " << motionRefiner.mdts.size()
                          << " unfinished micrographs" << std::endl;
            }
            else
            {
                std::cout << "   - Will not combine frames for any unfinished micrographs, "
                          << "just generate a STAR file" << std::endl;
            }
        }
    }
}

void FrameRecombiner::process(long g_start, long g_end)
{
    std::vector<Image<RFLOAT>> freqWeights;

    // Either calculate weights from FCC or from user-provided B-factors
    hasBfacs = bfacFn != "";

    if (!hasBfacs)
    {
        freqWeights = weightsFromFCC();
    }
    else
    {
        freqWeights = weightsFromBfacs();
    }

    int barstep;
    int my_nr_micrographs = g_end - g_start + 1;

    if (motionRefiner.verb > 0)
    {
        std::cout << " + Combining frames for all micrographs ... " << std::endl;
        init_progress_bar(my_nr_micrographs);
        barstep = XMIPP_MAX(1, my_nr_micrographs/ 60);
    }

    std::vector<ParFourierTransformer> fts(motionRefiner.nr_omp_threads);

    int pctot = 0;

    long nr_done = 0;
    FileName prevdir = "";

    for (long g = g_start; g <= g_end; g++)
    {
        const int pc = motionRefiner.mdts[g].numberOfObjects();
        if (pc == 0) continue;

        pctot += pc;

        std::vector<std::vector<Image<Complex>>> movie;
        movie = motionRefiner.loadMovie(g, pc, fts);

        FileName fn_root = motionRefiner.getOutputFileNameRoot(g);
        std::vector<std::vector<d2Vector>> shift;
        shift = MotionHelper::readTracks(fn_root+"_tracks.star");

        Image<RFLOAT> stack(s,s,1,pc);

        #pragma omp parallel for num_threads(motionRefiner.nr_omp_threads)
        for (int p = 0; p < pc; p++)
        {
            int threadnum = omp_get_thread_num();

            Image<Complex> sum(sh,s);
            sum.data.initZeros();

            Image<Complex> obs(sh,s);

            for (int f = 0; f < fc; f++)
            {
                shiftImageInFourierTransform(movie[p][f](), obs(), s, -shift[p][f].x, -shift[p][f].y);

                for (int y = 0; y < s; y++)
                for (int x = 0; x < sh; x++)
                {
                    sum(y,x) += freqWeights[f](y,x) * obs(y,x);
                }
            }

            Image<RFLOAT> real(s,s);

            fts[threadnum].inverseFourierTransform(sum(), real());

            for (int y = 0; y < s; y++)
            for (int x = 0; x < s; x++)
            {
                DIRECT_NZYX_ELEM(stack(), p, 0, y, x) = real(y,x);
            }
        }

        stack.write(fn_root+"_shiny.mrcs");

        if (motionRefiner.debug)
        {
            VtkHelper::writeTomoVTK(stack, fn_root+"_shiny.vtk");
        }

        for (int p = 0; p < pc; p++)
        {
            std::stringstream sts;
            sts << (p+1);
            motionRefiner.mdts[g].setValue(EMDL_IMAGE_NAME, sts.str() + "@" + fn_root+"_shiny.mrcs", p);
        }

        motionRefiner.mdts[g].write(fn_root+"_shiny.star");

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

std::vector<Image<RFLOAT>> FrameRecombiner::weightsFromFCC()
{
    if (motionRefiner.debug && motionRefiner.verb > 0)
    {
        std::cout << " + Summing up FCCs...\n";
    }

    Image<RFLOAT> fccData, fccWgh0, fccWgh1;
    Image<RFLOAT> fccDataMg, fccWgh0Mg, fccWgh1Mg;

    bool first = true;

    for (long g = 0; g < motionRefiner.mdts.size(); g++)
    {
        FileName fn_root = motionRefiner.getOutputFileNameRoot(g);

        fccDataMg.read(fn_root + "_FCC_cc.mrc");
        fccWgh0Mg.read(fn_root + "_FCC_w0.mrc");
        fccWgh1Mg.read(fn_root + "_FCC_w1.mrc");

        if (first)
        {
            sh = fccDataMg.data.xdim;
            s = 2 * (sh-1);
            fc = fccDataMg.data.ydim;

            fccData = Image<RFLOAT>(sh,fc);
            fccWgh0 = Image<RFLOAT>(sh,fc);
            fccWgh1 = Image<RFLOAT>(sh,fc);

            fccData.data.initZeros();
            fccWgh0.data.initZeros();
            fccWgh1.data.initZeros();

            first = false;
        }

        for (int y = 0; y < fc; y++)
        for (int x = 0; x < sh; x++)
        {
            fccData(y,x) += fccDataMg(y,x);
            fccWgh0(y,x) += fccWgh0Mg(y,x);
            fccWgh1(y,x) += fccWgh1Mg(y,x);
        }
    }

    Image<RFLOAT> fcc(sh,fc);

    for (int y = 0; y < fc; y++)
    for (int x = 0; x < sh; x++)
    {
        const double wgh = sqrt(fccWgh0Mg(y,x) * fccWgh1Mg(y,x));

        if (wgh > 0.0)
        {
            fcc(y,x) = fccData(y,x) / wgh;
        }
        else
        {
            fcc(y,x) = 0.0;
        }
    }

    if (motionRefiner.debug) std::cout << "done\n";

    k0 = (int) motionRefiner.angToPix(k0a);
    k1 = k1a > 0.0? (int) motionRefiner.angToPix(k1a) : sh;

    if (motionRefiner.verb > 0)
    {
        std::cout << " + Fitting B/k-factors between " << k0 << " and " << k1 << " pixels, or "
                  << k0a << " and " << k1a << " Angstrom ...\n";
    }

    std::pair<std::vector<d2Vector>,std::vector<double>> bkFacs
            = DamageHelper::fitBkFactors(fcc, k0, k1);

    std::vector<Image<RFLOAT>> freqWeights;
    freqWeights = DamageHelper::computeWeights(bkFacs.first, sh);

    const double cf = 8.0 * motionRefiner.angpix * motionRefiner.angpix * sh * sh;

    if (bfac_diag)
    {
        mktree(motionRefiner.outPath + "/bfacs");

        Image<RFLOAT> bfacFit = DamageHelper::renderBkFit(bkFacs, sh, fc);
        Image<RFLOAT> bfacFitNoScale = DamageHelper::renderBkFit(bkFacs, sh, fc, true);

        ImageLog::write(bfacFit, motionRefiner.outPath + "/bfacs/glob_Bk-fit");
        ImageLog::write(bfacFitNoScale, motionRefiner.outPath + "/bfacs/glob_Bk-fit_noScale");
        ImageLog::write(fcc, motionRefiner.outPath + "/bfacs/glob_Bk-data");
        ImageLog::write(freqWeights, motionRefiner.outPath + "/bfacs/freqWeights");

        std::ofstream bfacsDat(motionRefiner.outPath + "/bfacs/Bfac.dat");
        std::ofstream kfacsDat(motionRefiner.outPath + "/bfacs/kfac.dat");

        for (int i = 0; i < fc; i++)
        {
            double s = bkFacs.first[i].x;
            double b = -cf/(s*s);

            bfacsDat << i << " " << b << "\n";
            kfacsDat << i << " " << log(bkFacs.first[i].y) << "\n";
        }

        bfacsDat.close();
        kfacsDat.close();
    }

    MetaDataTable mdt;
    mdt.setName("perframe_bfactors");

    for (int f = 0; f < fc; f++ )
    {
        double s = bkFacs.first[f].x;
        double b = -cf/(s*s);
        double k = log(bkFacs.first[f].y);

        mdt.addObject();
        mdt.setValue(EMDL_IMAGE_FRAME_NR, f);
        mdt.setValue(EMDL_POSTPROCESS_BFACTOR, b);
        mdt.setValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, k);
    }

    mdt.write(motionRefiner.outPath + "/bfactors.star");

    // Also write out EPS plots of the B-factors and scale factors
    CPlot2D *plot2D=new CPlot2D("Polishing B-factors");
    plot2D->SetXAxisSize(600);
    plot2D->SetYAxisSize(400);
    plot2D->SetDrawLegend(false);
    plot2D->SetXAxisTitle("movie frame");
    plot2D->SetYAxisTitle("B-factor");
    mdt.addToCPlot2D(plot2D, EMDL_IMAGE_FRAME_NR, EMDL_POSTPROCESS_BFACTOR);
    plot2D->OutputPostScriptPlot(motionRefiner.outPath + "bfactors.eps");

    CPlot2D *plot2Db=new CPlot2D("Polishing scale-factors");
    plot2Db->SetXAxisSize(600);
    plot2Db->SetYAxisSize(400);
    plot2Db->SetDrawLegend(false);
    plot2Db->SetXAxisTitle("movie frame");
    plot2Db->SetYAxisTitle("Scale-factor");
    mdt.addToCPlot2D(plot2Db, EMDL_IMAGE_FRAME_NR, EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT);
    plot2Db->OutputPostScriptPlot(motionRefiner.outPath + "scalefactors.eps");

    return freqWeights;
}

std::vector<Image<RFLOAT>> FrameRecombiner::weightsFromBfacs()
{
    // initialization on the first line to avoid copying of return value
    std::vector<Image<RFLOAT>> freqWeights;

    MetaDataTable mdt;
    mdt.read(bfacFn);

    fc = mdt.numberOfObjects();

    std::vector<d2Vector> bkFacs(fc);

    const double cf = 8.0 * motionRefiner.angpix * motionRefiner.angpix * sh * sh;

    for (int f = 0; f < fc; f++)
    {
        int ff;
        mdt.getValue(EMDL_IMAGE_FRAME_NR, ff);

        double b, k;
        mdt.getValue(EMDL_POSTPROCESS_BFACTOR, b, f);
        mdt.getValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, k, f);

        bkFacs[f] = d2Vector(sqrt(-cf/b), exp(k));
    }

    freqWeights = DamageHelper::computeWeights(bkFacs, sh);

    if (bfac_diag)
    {
        mktree(motionRefiner.outPath + "bfacs");

        std::pair<std::vector<d2Vector>,std::vector<double>> bkFacs2;
        bkFacs2.first = bkFacs;
        bkFacs2.second = std::vector<double>(fc, 1.0);

        Image<RFLOAT> bfacFitNoScale = DamageHelper::renderBkFit(bkFacs2, sh, fc, true);

        ImageLog::write(bfacFitNoScale, motionRefiner.outPath + "/bfacs/glob_Bk-fit_noScale");
        ImageLog::write(freqWeights, motionRefiner.outPath + "/bfacs/freqWeights");

        std::ofstream bfacsDat(motionRefiner.outPath + "/bfacs/Bfac.dat");
        std::ofstream kfacsDat(motionRefiner.outPath + "/bfacs/kfac.dat");

        for (int i = 0; i < fc; i++)
        {
            double s = bkFacs[i].x;
            double b = -cf/(s*s);

            bfacsDat << i << " " << b << "\n";
            kfacsDat << i << " " << log(bkFacs[i].y) << "\n";
        }

        bfacsDat.close();
        kfacsDat.close();
    }

    return freqWeights;
}
