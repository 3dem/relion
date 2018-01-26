
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
#include <src/jaz/vtk_helper.h>
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
#include <src/jaz/motion_refinement.h>
#include <src/jaz/image_op.h>
#include <src/jaz/refinement_program.h>

#include <omp.h>

using namespace gravis;

class TiltFit : public RefinementProgram
{
    public:

        TiltFit();

            RFLOAT kmin;
            bool precomputed;
            std::string precomp;

            Image<Complex> lastXY;
            Image<RFLOAT> lastW;

        void readMoreOptions(IOParser& parser, int argc, char *argv[]);
        int _init();
        int _run();
};

TiltFit :: TiltFit()
:   RefinementProgram(true)
{
    noTilt = true;
}

int main(int argc, char *argv[])
{
    TiltFit tf;

    int rc0 = tf.init(argc, argv);
    if (rc0 != 0) return rc0;

    int rc1 = tf.run();
    if (rc1 != 0) return rc1;
}

void TiltFit::readMoreOptions(IOParser& parser, int argc, char *argv[])
{
    kmin = textToFloat(parser.getOption("--kmin", "Inner freq. threshold [Angst]", "30.0"));

    precomp = parser.getOption("--precomp", "Precomputed *_xy and *_w files from previous run (optional)", "");
    precomputed = precomp != "";

    noReference = precomputed;
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

        std::vector<FourierTransformer> fts(nr_omp_threads);

        const long gc = maxMG >= 0? maxMG+1 : mdts.size();

        for (long g = minMG; g < gc; g++)
        {
            std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

            CTF ctf0;
            ctf0.read(mdts[g], mdts[g], 0);

            std::vector<Image<Complex> > pred;
            std::vector<Image<Complex> > obsF;

            if (nr_omp_threads > 1)
            {
                pred = StackHelper::projectStackPar(&projectors[0], &mdts[g], nr_omp_threads);
            }
            else
            {
                pred = StackHelper::projectStack(&projectors[0], &mdts[g]);
            }

            obsF = StackHelper::loadStackFS(&mdts[g], imgPath, nr_omp_threads, &fts);

            #pragma omp parallel for num_threads(nr_omp_threads)
            for (long p = 0; p < pred.size(); p++)
            {
                int threadnum = omp_get_thread_num();

                CTF ctf(ctf0);
                TiltRefinement::updateTiltShift(pred[p], obsF[p], ctf, angpix, xyAcc[threadnum], wAcc[threadnum]);
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

    if (useFsc)
    {
        FilterHelper::multiply(wAccSum, freqWeight, wgh);
    }
    else
    {
        wgh = wAccSum;
    }

    Image<RFLOAT> debug(sh,s);

    for (int y = 0; y < s; y++)
    for (int x = 0; x < sh; x++)
    {
        double xx = x;
        double yy = y <= sh? y : y - s;
        double r = sqrt(xx*xx + yy*yy);

        if (r == 0 || sh/(2.0*angpix*r) > kmin)
        {
            wgh(y,x) = 0.0;
        }

        debug(y,x) = r == 0? 0.0 : sh/(2.0*angpix*r) - kmin;
    }

    VtkHelper::writeVTK(debug, outPath+"_debug.vtk");

    Image<RFLOAT> wghFull;
    FftwHelper::decenterDouble2D(wgh(), wghFull());
    VtkHelper::writeVTK(wghFull, outPath+"_weight.vtk");


    FftwHelper::decenterUnflip2D(phase.data, phaseFull.data);

    VtkHelper::writeVTK(phaseFull, outPath+"_delta_phase.vtk");
    phaseFull.write(outPath+"_delta_phase.mrc");
    wgh.write(outPath+"_weight.mrc");


    RFLOAT shift_x, shift_y, tilt_x, tilt_y;

    TiltRefinement::fitTiltShift(phase, wgh, Cs, lambda, angpix,
                                 &shift_x, &shift_y, &tilt_x, &tilt_y, &fit);

    FftwHelper::decenterUnflip2D(fit.data, fitFull.data);
    VtkHelper::writeVTK(fitFull, outPath+"_delta_phase_fit.vtk");
    fitFull.write(outPath+"_delta_phase_fit.mrc");

    std::ofstream os(outPath+"_beam_tilt.txt");
    os << "beamtilt_x = " << tilt_x << "\n";
    os << "beamtilt_y = " << tilt_y << "\n";
    os.close();

    double t1 = omp_get_wtime();

    std::cout << "elapsed: " << (t1 - t0) << " s\n";

    return 0;
}
