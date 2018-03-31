
#include <unistd.h>
#include <string.h>
#include <fstream>

#include <src/args.h>
#include <src/image.h>
#include <src/fftw.h>
#include <src/complex.h>
#include <src/metadata_table.h>
#include <src/backprojector.h>
#include <src/euler.h>
#include <src/jaz/image_log.h>
#include <src/jaz/slice_helper.h>
#include <src/jaz/spectral_helper.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/backprojection_helper.h>
#include <src/jaz/volume_converter.h>
#include <src/jaz/complex_io.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/resampling_helper.h>
#include <src/jaz/ctf_helper.h>
#include <src/jaz/defocus_refinement.h>
#include <src/jaz/magnification_refinement.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/tilt_refinement.h>
#include <src/jaz/image_op.h>
#include <src/jaz/refinement_program.h>
#include <src/jaz/parallel_ft.h>

#include <omp.h>

using namespace gravis;

class TiltFit : public RefinementProgram
{
    public:

        TiltFit();

            RFLOAT kmin,
                testtilt_x, testtilt_y,
                testtilt_xx, testtilt_xy, testtilt_yy;

            bool precomputed, aniso;
            std::string precomp;

            Image<Complex> lastXY;
            Image<RFLOAT> lastW;

        int readMoreOptions(IOParser& parser, int argc, char *argv[]);
        int _init();
        int _run();
};

TiltFit :: TiltFit()
:   RefinementProgram(true)
{
    noTilt = true;
    optStar = true;
    optReference = true;
}

int main(int argc, char *argv[])
{
    TiltFit tf;

    int rc0 = tf.init(argc, argv);
    if (rc0 != 0) return rc0;

    int rc1 = tf.run();
    if (rc1 != 0) return rc1;
}

int TiltFit::readMoreOptions(IOParser& parser, int argc, char *argv[])
{
    kmin = textToFloat(parser.getOption("--kmin", "Inner freq. threshold [Angst]", "30.0"));

    precomp = parser.getOption("--precomp", "Precomputed *_xy and *_w files from previous run (optional)", "");
    aniso = parser.checkOption("--aniso", "Use anisotropic coma model");

    testtilt_x = textToFloat(parser.getOption("--testtilt_x", "Test value", "0.0"));
    testtilt_y = textToFloat(parser.getOption("--testtilt_y", "Test value", "0.0"));

    testtilt_xx = textToFloat(parser.getOption("--testtilt_xx", "Test value", "1.0"));
    testtilt_xy = textToFloat(parser.getOption("--testtilt_xy", "Test value", "0.0"));
    testtilt_yy = textToFloat(parser.getOption("--testtilt_yy", "Test value", "1.0"));

    precomputed = precomp != "";

    noReference = precomputed;
    noStar = starFn == "";

    bool allGood = true;

    if (starFn == "" && !precomputed)
    {
        std::cerr << "A .star file (--i) is required if no precomputed pixel-fit is available (--precomp).\n";
        allGood = false;
    }

    if (reconFn0 == "" && !precomputed)
    {
        std::cerr << "A reference map (--m) is required if no precomputed pixel-fit is available (--precomp).\n";
        allGood = false;
    }

    if (starFn == "" && angpix == 0)
    {
        std::cerr << "Pixel resolution (--angpix) is required if no .star file is provided (--i).\n";
        allGood = false;
    }

    if (!allGood) return 13;
    else return 0;
}

int TiltFit::_init()
{
    if (precomputed)
    {
        ComplexIO::read(lastXY, precomp+"_xy", ".mrc");
        lastW.read(precomp+"_w.mrc");

        s = lastXY.data.ydim;
        sh = s/2 + 1;
    }

    return 0;
}

int TiltFit::_run()
{
    Image<Complex> xyAccSum(sh,s);
    Image<RFLOAT> wAccSum(sh,s);

    double t0 = omp_get_wtime();

    if (!precomputed)
    {
        xyAccSum.data.initZeros();
        wAccSum.data.initZeros();

        std::vector<Image<Complex>> xyAcc(nr_omp_threads);
        std::vector<Image<RFLOAT>> wAcc(nr_omp_threads);

        for (int i = 0; i < nr_omp_threads; i++)
        {
            xyAcc[i] = Image<Complex>(sh,s);
            xyAcc[i].data.initZeros();

            wAcc[i] = Image<RFLOAT>(sh,s);
            wAcc[i].data.initZeros();
        }

        std::vector<MetaDataTable> mdts = StackHelper::splitByStack(&mdt0);

        std::vector<ParFourierTransformer> fts(nr_omp_threads);

        for (long g = minMG; g <= gc; g++)
        {
            std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

            CTF ctf0;
            ctf0.read(mdts[g], mdts[g], 0);

            const int pc = mdts[g].numberOfObjects();

            std::vector<Image<Complex>> pred(pc);
            std::vector<Image<Complex>> obsF;

            #pragma omp parallel for num_threads(nr_omp_threads)
            for (long p = 0; p < pc; p++)
            {
                pred[p] = obsModel.predictObservation(
                    projectors[0], mdts[g], p, false, false);
            }

            obsF = StackHelper::loadStackFS(&mdts[g], imgPath, nr_omp_threads, &fts);

            if (ABS(testtilt_x) > 0.0 || ABS(testtilt_y) > 0.0)
            {
                if (g == minMG)
                {
                    std::cout << "applying test tilt to images: "
                              << testtilt_x << ", " << testtilt_y << " / "
                              << testtilt_xx << ", " << testtilt_xy << ", " << testtilt_yy << "\n";
                }

                obsF = StackHelper::applyBeamTiltPar(
                            obsF, Cs, lambda, angpix,
                            testtilt_x, testtilt_y,
                            testtilt_xx, testtilt_xy, testtilt_yy,
                            nr_omp_threads);
            }

            #pragma omp parallel for num_threads(nr_omp_threads)
            for (long p = 0; p < pred.size(); p++)
            {
                int threadnum = omp_get_thread_num();

                CTF ctf(ctf0);
                TiltRefinement::updateTiltShift(
                    pred[p], obsF[p], ctf, angpix, xyAcc[threadnum], wAcc[threadnum]);
            }
        }

        for (int i = 0; i < nr_omp_threads; i++)
        {
            ImageOp::linearCombination(xyAccSum, xyAcc[i], 1.0, 1.0, xyAccSum);
            ImageOp::linearCombination(wAccSum, wAcc[i], 1.0, 1.0, wAccSum);
        }

        ComplexIO::write(xyAccSum(), outPath+"_xy", ".mrc");
        wAccSum.write(outPath+"_w.mrc");
    }
    else
    {
        xyAccSum = lastXY;
        wAccSum = lastW;
    }

    Image<RFLOAT> wgh, phase, fit, phaseFull, fitFull;

    FilterHelper::getPhase(xyAccSum, phase);

    Image<Complex> xyNrm(sh,s);

    if (useFsc)
    {
        FilterHelper::multiply(wAccSum, freqWeight, wgh);
    }
    else
    {
        wgh = wAccSum;
    }

    //Image<RFLOAT> debug(sh,s);

    for (int y = 0; y < s; y++)
    for (int x = 0; x < sh; x++)
    {
        double xx = x;
        double yy = y <= sh? y : y - s;
        double r = sqrt(xx*xx + yy*yy);

        if (r == 0 || 2.0*sh*angpix/r > kmin)
        {
            wgh(y,x) = 0.0;
        }

        //debug(y,x) = r == 0? 0.0 : sh/(2.0*angpix*r) - kmin;

        xyNrm(y,x) = wAccSum(y,x) > 0.0? xyAccSum(y,x)/wAccSum(y,x) : Complex(0.0, 0.0);
    }

    //ImageLog::write(debug, outPath+"_debug");

    Image<RFLOAT> wghFull;
    FftwHelper::decenterDouble2D(wgh(), wghFull());
    ImageLog::write(wghFull, outPath+"_weight_full");


    FftwHelper::decenterUnflip2D(phase.data, phaseFull.data);

    ImageLog::write(phaseFull, outPath+"_delta_phase");
    wgh.write(outPath+"_weight.mrc");

    double shift_x, shift_y, tilt_x, tilt_y;

    TiltRefinement::fitTiltShift(phase, wgh, Cs, lambda, angpix,
                                 &shift_x, &shift_y, &tilt_x, &tilt_y, &fit);



    FftwHelper::decenterUnflip2D(fit.data, fitFull.data);
    ImageLog::write(fitFull, outPath+"_delta_phase_fit");

    std::ofstream os(outPath+"_beam_tilt_0.txt");
    os << "beamtilt_x = " << tilt_x << "\n";
    os << "beamtilt_y = " << tilt_y << "\n";
    os.close();

    double tilt_xx, tilt_xy, tilt_yy;

    if (aniso)
    {
        TiltRefinement::optimizeAnisoTilt(
            xyNrm, wgh, Cs, lambda, angpix, false,
            shift_x, shift_y, tilt_x, tilt_y,
            &shift_x, &shift_y, &tilt_x, &tilt_y,
            &tilt_xx, &tilt_xy, &tilt_yy, &fit);
    }
    else
    {
        TiltRefinement::optimizeTilt(
            xyNrm, wgh, Cs, lambda, angpix, false,
            shift_x, shift_y, tilt_x, tilt_y,
            &shift_x, &shift_y, &tilt_x, &tilt_y, &fit);
    }


    FftwHelper::decenterUnflip2D(fit.data, fitFull.data);
    ImageLog::write(fitFull, outPath+"_delta_phase_iter_fit");

    std::ofstream os2(outPath+"_beam_tilt_1.txt");
    os2 << "beamtilt_x = " << tilt_x << "\n";
    os2 << "beamtilt_y = " << tilt_y << "\n";
    os2.close();

    setForAll(EMDL_IMAGE_BEAMTILT_X, tilt_x);
    setForAll(EMDL_IMAGE_BEAMTILT_Y, tilt_y);

    mdt0.write(outPath+"_particles.star");

    double t1 = omp_get_wtime();

    std::cout << "elapsed: " << (t1 - t0) << " s\n";

    return 0;
}
