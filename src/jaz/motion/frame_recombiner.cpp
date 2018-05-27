#include "frame_recombiner.h"
#include "motion_refiner.h"
#include "motion_helper.h"

#include <src/jaz/micrograph_handler.h>
#include <src/jaz/obs_model.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/vtk_helper.h>
#include <src/jaz/damage_helper.h>
#include <src/jaz/image_log.h>
#include <src/filename.h>

using namespace gravis;

FrameRecombiner::FrameRecombiner()
{}

void FrameRecombiner::read(IOParser& parser, int argc, char* argv[])
{
    parser.addSection("Combine frames options");

    doCombineFrames = parser.checkOption("--combine_frames", "Combine movie frames into polished particles.");
    k0a = textToFloat(parser.getOption("--bfac_minfreq", "Min. frequency used in B-factor fit [Angst]", "20"));
    k1a = textToFloat(parser.getOption("--bfac_maxfreq", "Max. frequency used in B-factor fit [Angst]", "-1"));
	b_scale = textToFloat(parser.getOption("--bfac_scale", "Scale of B/k-factors", "1"));
    bfacFn = parser.getOption("--bfactors", "A .star file with external B/k-factors", "");
    bfac_diag = parser.checkOption("--diag_bfactor", "Write out B/k-factor diagnostic data");
}

void FrameRecombiner::init(
    const std::vector<MetaDataTable>& allMdts,
    int verb, int s, int fc, double maxFreq,
    int nr_omp_threads, std::string outPath, bool debug,
    ObservationModel* obsModel,
    MicrographHandler* micrographHandler)
{
    this->verb = verb;
    this->s = s;
    this->sh = s/2 + 1;
    this->fc = fc;
    this->nr_omp_threads = nr_omp_threads;
    this->outPath = outPath;
    this->debug = debug;
    this->obsModel = obsModel;
    this->micrographHandler = micrographHandler;
    this->angpix = obsModel->angpix;
	this->maxFreq = maxFreq;


    // Either calculate weights from FCC or from user-provided B-factors
    const bool hasBfacs = bfacFn != "";

    if (!hasBfacs)
    {
        freqWeights = weightsFromFCC(allMdts);
    }
    else
    {
        freqWeights = weightsFromBfacs();
    }
}

void FrameRecombiner::process(const std::vector<MetaDataTable>& mdts, long g_start, long g_end)
{
    int barstep;
    int my_nr_micrographs = g_end - g_start + 1;

    if (verb > 0)
    {
        std::cout << " + Combining frames for all micrographs ... " << std::endl;
        init_progress_bar(my_nr_micrographs);
        barstep = XMIPP_MAX(1, my_nr_micrographs/ 60);
    }

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    int pctot = 0;
    long nr_done = 0;

    for (long g = g_start; g <= g_end; g++)
    {
        const int pc = mdts[g].numberOfObjects();
        if (pc == 0) continue;

        pctot += pc;

        std::vector<std::vector<Image<Complex>>> movie;
        movie = micrographHandler->loadMovie(mdts[g], s, angpix, fts);

        FileName fn_root = MotionRefiner::getOutputFileNameRoot(outPath, mdts[g]);
        std::vector<std::vector<d2Vector>> shift;
        shift = MotionHelper::readTracks(fn_root+"_tracks.star");

        Image<RFLOAT> stack(s,s,1,pc);

        #pragma omp parallel for num_threads(nr_omp_threads)
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

        if (debug)
        {
            VtkHelper::writeTomoVTK(stack, fn_root+"_shiny.vtk");
        }

        MetaDataTable mdtOut = mdts[g];

        for (int p = 0; p < pc; p++)
        {
            std::stringstream sts;
            sts << (p+1);
            mdtOut.setValue(EMDL_IMAGE_NAME, sts.str() + "@" + fn_root+"_shiny.mrcs", p);
        }

        mdtOut.write(fn_root+"_shiny.star");

        nr_done++;

        if (verb > 0 && nr_done % barstep == 0)
        {
            progress_bar(nr_done);
        }
    }

    if (verb > 0)
    {
        progress_bar(my_nr_micrographs);
    }
}

std::vector<Image<RFLOAT>> FrameRecombiner::weightsFromFCC(
        const std::vector<MetaDataTable>& allMdts)
{
    if (debug && verb > 0)
    {
        std::cout << " + Summing up FCCs..." << std::endl;
    }

    Image<RFLOAT> fccData, fccWgh0, fccWgh1;
    Image<RFLOAT> fccDataMg, fccWgh0Mg, fccWgh1Mg;

    bool first = true;

    // Compute B/k-factors from all available FCCs (allMdts),
    // even if only a subset of micrographs (chosenMdts) is being recombined.
    for (long g = 0; g < allMdts.size(); g++)
    {
        FileName fn_root = MotionRefiner::getOutputFileNameRoot(outPath, allMdts[g]);

        if (!( exists(fn_root + "_FCC_cc.mrc")
            && exists(fn_root + "_FCC_w0.mrc")
            && exists(fn_root + "_FCC_w1.mrc")))
        {
            continue;
        }

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

    if (debug) std::cout << "done\n";

    k0 = (int) obsModel->angToPix(k0a, s);

    if (!outerFreqKnown())
    {
        k1a = maxFreq;
    }
	
	k1 = (int) obsModel->angToPix(k1a, s);

    if (verb > 0)
    {
        std::cout << " + Fitting B/k-factors between " << k0 << " and " << k1 << " pixels, or "
                  << k0a << " and " << k1a << " Angstrom ..." << std::endl;
    }

    std::pair<std::vector<d2Vector>,std::vector<double>> bkFacs
            = DamageHelper::fitBkFactors(fcc, k0, k1);
	
	for (int f = 0; f < fc; f++)
	{
		bkFacs.first[f].x = bkFacs.first[f].x / sqrt(b_scale);
		bkFacs.first[f].y = pow(bkFacs.first[f].y, b_scale);
	}

    std::vector<Image<RFLOAT>> freqWeights;
    freqWeights = DamageHelper::computeWeights(bkFacs.first, sh);

    const double cf = 8.0 * angpix * angpix * sh * sh;

    if (bfac_diag)
    {
        mktree(outPath + "/bfacs");

        Image<RFLOAT> bfacFit = DamageHelper::renderBkFit(bkFacs, sh, fc);
        Image<RFLOAT> bfacFitNoScale = DamageHelper::renderBkFit(bkFacs, sh, fc, true);

        ImageLog::write(bfacFit, outPath + "/bfacs/glob_Bk-fit");
        ImageLog::write(bfacFitNoScale, outPath + "/bfacs/glob_Bk-fit_noScale");
        ImageLog::write(fcc, outPath + "/bfacs/glob_Bk-data");
        ImageLog::write(freqWeights, outPath + "/bfacs/freqWeights");

        std::ofstream bfacsDat(outPath + "/bfacs/Bfac.dat");
        std::ofstream kfacsDat(outPath + "/bfacs/kfac.dat");

        for (int i = 0; i < fc; i++)
        {
            double s = bkFacs.first[i].x;
            double b = -cf/(s*s);

            bfacsDat << i << " " << b << std::endl;
            kfacsDat << i << " " << log(bkFacs.first[i].y) << std::endl;
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

    mdt.write(outPath + "/bfactors.star");

    // Also write out EPS plots of the B-factors and scale factors
    CPlot2D *plot2D = new CPlot2D("Polishing B-factors");
    plot2D->SetXAxisSize(600);
    plot2D->SetYAxisSize(400);
    plot2D->SetDrawLegend(false);
    plot2D->SetXAxisTitle("movie frame");
    plot2D->SetYAxisTitle("B-factor");
    mdt.addToCPlot2D(plot2D, EMDL_IMAGE_FRAME_NR, EMDL_POSTPROCESS_BFACTOR);
    plot2D->OutputPostScriptPlot(outPath + "bfactors.eps");

    CPlot2D *plot2Db=new CPlot2D("Polishing scale-factors");
    plot2Db->SetXAxisSize(600);
    plot2Db->SetYAxisSize(400);
    plot2Db->SetDrawLegend(false);
    plot2Db->SetXAxisTitle("movie frame");
    plot2Db->SetYAxisTitle("Scale-factor");
    mdt.addToCPlot2D(plot2Db, EMDL_IMAGE_FRAME_NR, EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT);
    plot2Db->OutputPostScriptPlot(outPath + "scalefactors.eps");
	
	delete plot2D;

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
	
	double bfacOff = 0.0;
	
	for (int f = 0; f < fc; f++)
    {
        double b;
        mdt.getValue(EMDL_POSTPROCESS_BFACTOR, b, f);

        if (b > bfacOff) bfacOff = b;
    }
	
    const double cf = 8.0 * angpix * angpix * sh * sh;

    for (int f = 0; f < fc; f++)
    {
        double b, k;
        mdt.getValue(EMDL_POSTPROCESS_BFACTOR, b, f);
        mdt.getValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, k, f);

		bkFacs[f] = d2Vector(sqrt(-cf/(b-bfacOff-1)), exp(k));
    }

    freqWeights = DamageHelper::computeWeights(bkFacs, sh);

    if (bfac_diag)
    {
        mktree(outPath + "bfacs");

        std::pair<std::vector<d2Vector>,std::vector<double>> bkFacs2;
        bkFacs2.first = bkFacs;
        bkFacs2.second = std::vector<double>(fc, 1.0);

        Image<RFLOAT> bfacFitNoScale = DamageHelper::renderBkFit(bkFacs2, sh, fc, true);

        ImageLog::write(bfacFitNoScale, outPath + "/bfacs/glob_Bk-fit_noScale");
        ImageLog::write(freqWeights, outPath + "/bfacs/freqWeights");

        std::ofstream bfacsDat(outPath + "/bfacs/Bfac.dat");
        std::ofstream kfacsDat(outPath + "/bfacs/kfac.dat");

        for (int i = 0; i < fc; i++)
        {
            double s = bkFacs[i].x;
            double b = -cf/(s*s);

            bfacsDat << i << " " << b << std::endl;
            kfacsDat << i << " " << log(bkFacs[i].y) << std::endl;
        }

        bfacsDat.close();
        kfacsDat.close();
    }

    return freqWeights;
}

bool FrameRecombiner::doingRecombination()
{
	return doCombineFrames;
}

bool FrameRecombiner::outerFreqKnown()
{
	return k1a > 0.0;
}

std::vector<MetaDataTable> FrameRecombiner::findUnfinishedJobs(
        const std::vector<MetaDataTable> &mdts, std::string path)
{
    std::vector<MetaDataTable> out(0);

    const int gc = mdts.size();

    for (int g = 0; g < gc; g++)
    {
        std::string fn_root = MotionRefiner::getOutputFileNameRoot(path, mdts[g]);

        if (!isJobFinished(fn_root))
        {
            out.push_back(mdts[g]);
        }
    }

    return out;
}

bool FrameRecombiner::isJobFinished(std::string filenameRoot)
{
    return exists(filenameRoot+"_shiny.mrcs")
            && exists(filenameRoot+"_shiny.star");
}
