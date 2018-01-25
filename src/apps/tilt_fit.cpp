
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

#include <omp.h>

using namespace gravis;

#define BILLION 1000000000L

int main(int argc, char *argv[])
{
    std::string starFn, reconFn, tiltFn, inPath, fscFn, maskFn, precomp;
    long maxMG = -1;
    int nr_omp_threads;
    bool useFsc, precomputed;
    RFLOAT angpix, pad, kmin;

    Image<RFLOAT> map, dummy, lastW;
    Image<Complex> lastXY;
    MetaDataTable fscMdt;

    IOParser parser;

    try
    {
        parser.setCommandLine(argc, argv);

        parser.addSection("General options");

        starFn = parser.getOption("--i", "Input STAR file with the projection images and their orientations", "");
        reconFn = parser.getOption("--m", "Input MRC file of a initial reconstruction", "");
        fscFn = parser.getOption("--f", "Input STAR file with the FSC of the initial reconstruction (optional)", "");
        maskFn = parser.getOption("--mask", "Mask for the initial reconstruction", "");

        precomp = parser.getOption("--precomp", "Precomputed *_xy and *_w files from previous run (optional)", "");
        precomputed = precomp != "";

        useFsc = fscFn != "";

        tiltFn = parser.getOption("--out", "Output filename prefix", "");
        inPath = parser.getOption("--img", "Path to images", "");

        angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (angst/pix)", "0."));
        pad = textToFloat(parser.getOption("--pad", "Padding factor", "2.0"));
        kmin = textToFloat(parser.getOption("--kmin", "Inner freq. threshold [Angst]", "30.0"));

        nr_omp_threads = textToInteger(parser.getOption("--jomp", "Number of OMP threads", "1"));

        maxMG = textToInteger(getParameter(argc, argv, "--max_MG", "-1"));

        if (!precomputed)
        {
            try
            {
                map.read(reconFn);
            }
            catch (RelionError XE)
            {
                std::cout << "Unable to read map: " << reconFn << "\n";
                exit(1);
            }
        }

        if (useFsc)
        {
            fscMdt.read(fscFn, "fsc");

            bool allGood = true;

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

        if (precomputed)
        {
            ComplexIO::read(lastXY, precomp+"_xy", ".mrc");
            lastW.read(precomp+"_w.mrc");
        }
    }
    catch (RelionError XE)
    {
        parser.writeUsage(std::cout);
        std::cerr << XE;
        exit(1);
    }

    if (!precomputed && (map.data.xdim != map.data.ydim || map.data.ydim != map.data.zdim))
    {
        REPORT_ERROR(reconFn + " is not cubical.\n");
    }

    const int s = precomputed? lastXY.data.ydim : map.data.xdim;
    const int sh = s/2 + 1;

    if (!precomputed && maskFn != "")
    {
        Image<RFLOAT> mask, maskedRef;
        mask.read(maskFn);

        ImageOp::multiply(mask, map, maskedRef);
        map = maskedRef;
    }

    Image<Complex> xyAccSum(sh,s);
    Image<RFLOAT> wAccSum(sh,s);

    Image<RFLOAT> imgSnr;

    if (useFsc)
    {
        RefinementHelper::computeSNR(&fscMdt, imgSnr);
    }

    double t0 = omp_get_wtime();


    MetaDataTable mdt0;
    mdt0.read(starFn);

    RFLOAT Cs, lambda, kV;

    mdt0.getValue(EMDL_CTF_CS, Cs, 0);
    mdt0.getValue(EMDL_CTF_VOLTAGE, kV, 0);

    RFLOAT V = kV * 1e3;
    lambda = 12.2643247 / sqrt(V * (1.0 + V * 0.978466e-6));


    if (!precomputed)
    {
        xyAccSum.data.initZeros();
        wAccSum.data.initZeros();

        Projector projector(s, TRILINEAR, pad, 10, 2);
        projector.computeFourierTransformMap(map.data, dummy.data, 2*s);

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

        if (angpix <= 0.0)
        {
            RFLOAT mag, dstep;
            mdts[0].getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
            mdts[0].getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
            angpix = 10000 * dstep / mag;
        }

        std::vector<FourierTransformer> fts(nr_omp_threads);

        const long gc = maxMG >= 0? maxMG : mdts.size();

        for (long g = 0; g < gc; g++)
        {
            std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

            CTF ctf0;
            ctf0.read(mdts[g], mdts[g], 0);

            std::vector<Image<Complex> > pred;
            std::vector<Image<Complex> > obsF;

            if (nr_omp_threads > 1)
            {
                pred = StackHelper::projectStackPar(&projector, &mdts[g], nr_omp_threads);
            }
            else
            {
                pred = StackHelper::projectStack(&projector, &mdts[g]);
            }

            obsF = StackHelper::loadStackFS(&mdts[g], inPath, nr_omp_threads, &fts);

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

        ComplexIO::write(xyAccSum(), tiltFn+"_xy", ".mrc");
        wAccSum.write(tiltFn+"_w.mrc");
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
        FilterHelper::multiply(wAccSum, imgSnr, wgh);
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

    VtkHelper::writeVTK(debug, tiltFn+"_debug.vtk");

    Image<RFLOAT> wghFull;
    FftwHelper::decenterDouble2D(wgh(), wghFull());
    VtkHelper::writeVTK(wghFull, tiltFn+"_weight.vtk");


    FftwHelper::decenterUnflip2D(phase.data, phaseFull.data);

    VtkHelper::writeVTK(phaseFull, tiltFn+"_delta_phase.vtk");
    phaseFull.write(tiltFn+"_delta_phase.mrc");
    wgh.write(tiltFn+"_weight.mrc");


    RFLOAT shift_x, shift_y, tilt_x, tilt_y;

    TiltRefinement::fitTiltShift(phase, wgh, Cs, lambda, angpix,
                                 &shift_x, &shift_y, &tilt_x, &tilt_y, &fit);

    FftwHelper::decenterUnflip2D(fit.data, fitFull.data);
    VtkHelper::writeVTK(fitFull, tiltFn+"_delta_phase_fit.vtk");
    fitFull.write(tiltFn+"_delta_phase_fit.mrc");

    std::ofstream os(tiltFn+"_beam_tilt.txt");
    os << "beamtilt_x = " << tilt_x << "\n";
    os << "beamtilt_y = " << tilt_y << "\n";
    os.close();

    double t1 = omp_get_wtime();

    std::cout << "elapsed: " << (t1 - t0) << " s\n";

    return 0;
}
