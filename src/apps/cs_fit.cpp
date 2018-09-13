
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
#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/backprojection_helper.h>
#include <src/jaz/volume_converter.h>
#include <src/jaz/complex_io.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/resampling_helper.h>
#include <src/jaz/ctf_helper.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/img_proc/image_op.h>
#include <src/jaz/refinement_program.h>
#include <src/jaz/parallel_ft.h>

#include <src/jaz/ctf/defocus_helper.h>

#include <omp.h>

using namespace gravis;

class CsFit : public RefinementProgram
{
    public:

        CsFit();

        bool aniso;
        double range, baseCs;
        int sectors, samples;
        RFLOAT defocusRange;

        int readMoreOptions(IOParser& parser, int argc, char *argv[]);
        int _init();
        int _run();
};

int main(int argc, char *argv[])
{
    CsFit mt;

    int rc0 = mt.init(argc, argv);
    if (rc0 != 0) return rc0;

    int rc1 = mt.run();
    if (rc1 != 0) return rc1;
}

CsFit::CsFit()
: RefinementProgram(true)
{}

int CsFit::readMoreOptions(IOParser& parser, int argc, char *argv[])
{
    aniso = parser.checkOption("--aniso", "Anisotropic C_s");
    range = textToFloat(parser.getOption("--range", "C_s scan range", "0.5"));
    baseCs = textToFloat(parser.getOption("--base", "C_s base value", "2.7"));
    sectors = textToInteger(parser.getOption("--sec", "Number of radial sectors", "8"));
    samples = textToInteger(parser.getOption("--n", "Number of C_s samples", "31"));
    defocusRange = textToFloat(parser.getOption("--range", "Defocus scan range (in A)", "500."));

    return 0;
}

int CsFit::_init()
{
    return 0;
}

int CsFit::_run()
{
    double t0 = omp_get_wtime();

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    std::vector<double> csVals(samples);

    for (int i = 0; i < samples; i++)
    {
        csVals[i] = baseCs + (i / (double)(samples-1) - 0.5) * range;
    }

    Image<RFLOAT> sectorIndex(sh,s);

    for (int y = 0; y < s; y++)
    for (int x = 0; x < sh; x++)
    {
        double xx = x;
        double yy = y < sh? y : y - s;

        double phi;

        if (y == 0 && x == 0)
        {
            phi = 0.0;
        }
        else
        {
            phi = atan2(yy,xx);
        }

        int sec = (int)(2*sectors * (phi + PI/2) / (2.0 * PI) + 0.5 + 2*sectors) % (2*sectors);

        sectorIndex(y,x) = (RFLOAT)sec;
    }

    ImageLog::write(sectorIndex, "debug/sectors");

    std::vector<std::vector<double>> totalCostIso(nr_omp_threads);

    for (int i = 0; i < nr_omp_threads; i++)
    {
        totalCostIso[i] = std::vector<double>(samples, 0.0);
    }

    for (long g = minMG; g <= gc; g++)
    {
        std::stringstream stsg;
        stsg << g;

        std::cout << "micrograph " << g << " / " << gc <<"\n";

        const int pc = mdts[g].numberOfObjects();

        std::vector<Image<Complex> > obsF;

        obsF = StackHelper::loadStackFS(&mdts[g], imgPath, nr_omp_threads, &fts);

        if (applyTilt)
        {
            obsF = StackHelper::applyBeamTiltPar(
                        obsF, Cs, lambda, angpix,
                        beamtilt_x, beamtilt_y,
                        beamtilt_xx, beamtilt_xy, beamtilt_yy,
                        nr_omp_threads);
        }

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (long p = 0; p < pc; p++)
        {
            int threadnum = omp_get_thread_num();

            std::stringstream stsp;
            stsp << p;

            CTF ctf0;
            ctf0.read(mdts[g], mdts[g], p);

            Image<Complex> pred;

            pred = obsModel.predictObservation(
                projectors[0], mdts[g], p, false, true);

            if (aniso)
            {

            }
            else
            {
                for (int i = 0; i < samples; i++)
                {
                    ctf0.Cs = csVals[i];
                    ctf0.initialise();

                    double u, v;

                    double err = DefocusHelper::findDefocus1D(
                        pred, obsF[p], freqWeight, ctf0, angpix, &u, &v, defocusRange, 11, 1);

                    totalCostIso[threadnum][i] += err;
                }
            }
        }
    }

    if (aniso)
    {

    }
    else
    {
        for (int i = 1; i < nr_omp_threads; i++)
        {
            for (int j = 0; j < samples; j++)
            {
                totalCostIso[0][j] += totalCostIso[i][j];
            }
        }

        std::ofstream out(outPath + "_Cs_iso.dat");

        for (int j = 0; j < samples; j++)
        {
            out << csVals[j] << " " << std::setprecision(12) << totalCostIso[0][j] << "\n";
        }
    }

    double t1 = omp_get_wtime();

    std::cout << "elapsed: " << (t1 - t0) << "s \n";

    return 0;
}
