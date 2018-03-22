
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
#include <src/jaz/motion_refinement.h>
#include <src/jaz/image_op.h>
#include <src/jaz/Fourier_helper.h>
#include <src/jaz/fsc_helper.h>
#include <src/jaz/damage_helper.h>
#include <src/jaz/interpolation.h>
#include <src/jaz/distribution_helper.h>
#include <src/jaz/noise_helper.h>
#include <src/jaz/convolution_helper.h>
#include <src/jaz/motion_em.h>
#include <src/jaz/local_motion_fit.h>
#include <src/jaz/gradient_descent.h>
#include <src/jaz/refinement_program.h>
#include <src/jaz/parallel_ft.h>

#include <omp.h>

using namespace gravis;

class FrameRecomb : public RefinementProgram
{
    public:

        FrameRecomb();

            int k0, k1;
            double k0a, k1a;
            std::string trackFn;

        int readMoreOptions(IOParser& parser, int argc, char *argv[]);
        int _init();
        int _run();
};

int main(int argc, char *argv[])
{
    FrameRecomb mt;

    int rc0 = mt.init(argc, argv);
    if (rc0 != 0) return rc0;

    int rc1 = mt.run();
    if (rc1 != 0) return rc1;
}

FrameRecomb::FrameRecomb()
: RefinementProgram(true, true)
{
    noReference = true;
    noTilt = true;
}

int FrameRecomb::readMoreOptions(IOParser& parser, int argc, char *argv[])
{
    trackFn = parser.getOption("--t", "Directory with input tracks and FCCs");

    k0a = textToFloat(parser.getOption("--k0", "Min. frequency used in B-factor fit (Angst)", "20"));
    k1a = textToFloat(parser.getOption("--k1", "Max. frequency used in B-factor fit (Angst)", "-1"));

    return 0;
}

int FrameRecomb::_init()
{
    return 0;
}

int FrameRecomb::_run()
{
    if (debug) std::cout << "summing up FCCs...\n";

    Image<RFLOAT> fccData, fccWgh0, fccWgh1;
    Image<RFLOAT> fccDataMg, fccWgh0Mg, fccWgh1Mg;

    bool first = true;

    for (long g = 0; g < mdts.size(); g++)
    {
        std::string tag = getMicrographTag(g);
        std::string tfn = trackFn + "/" + tag;

        try
        {
            fccDataMg.read(tfn + "_FCC_cc.mrc");
            fccWgh0Mg.read(tfn + "_FCC_w0.mrc");
            fccWgh1Mg.read(tfn + "_FCC_w1.mrc");

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
        }
        catch (RelionError e)
        {
            std::cerr << "    Warning: unable to read FCCs in " << tfn << "_FCC_*.mrc";
            continue;
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


    k0 = (int) angstToPixFreq(k0a);
    k1 = k1a > 0.0? (int) angstToPixFreq(k1a) : sh;

    std::cout << "fitting B/k-factors between " << k0 << " and " << k1 << " pixels...\n";

    std::pair<std::vector<d2Vector>,std::vector<double>> bkFacs
            = DamageHelper::fitBkFactors(fcc, k0, k1);

    std::vector<Image<RFLOAT>> freqWeights;
    freqWeights = DamageHelper::computeWeights(bkFacs.first, sh);

    if (debug)
    {
        Image<RFLOAT> bfacFit = DamageHelper::renderBkFit(bkFacs, sh, fc);
        Image<RFLOAT> bfacFitNoScale = DamageHelper::renderBkFit(bkFacs, sh, fc, true);

        ImageLog::write(bfacFit, "bfacs/glob_Bk-fit");
        ImageLog::write(bfacFitNoScale, "bfacs/glob_Bk-fit_noScale");
        ImageLog::write(fcc, "bfacs/glob_Bk-data");
        ImageLog::write(freqWeights, "bfacs/freqWeights");

        std::ofstream bfacsDat("bfacs/Bfac.dat");
        std::ofstream kfacsDat("bfacs/kfac.dat");

        const double cf = -8.0 * angpix*angpix * sh*sh;

        for (int i = 0; i < fc; i++)
        {
            double s = bkFacs.first[i].x;
            double b = cf/(s*s);

            bfacsDat << i << " " << b << "\n";
            kfacsDat << i << " " << log(bkFacs.first[i].y) << "\n";
        }

        bfacsDat.close();
        kfacsDat.close();
    }

    loadInitialMovieValues();

    std::cout << "mg range: " << g0 << ".." << gc << "\n";

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    double t0 = omp_get_wtime();

    for (long g = g0; g <= gc; g++)
    {
        std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

        std::stringstream stsg;
        stsg << g;

        const int pc = mdts[g].numberOfObjects();

        std::vector<std::vector<Image<Complex>>> movie;

        try
        {
            movie = loadMovie(g, pc, fts);
        }
        catch (RelionError XE)
        {
            std::cerr << "warning: unable to load micrograph #" << (g+1) << "\n";
            continue;
        }

        const int sh = movie[0][0]().xdim;
        const int s = movie[0][0]().ydim;

        std::string tag = getMicrographTag(g);
        std::vector<std::vector<d2Vector>> shift;

        std::string tfn = trackFn + "/" + tag + "_tracks.star";

        try
        {
            shift = MotionRefinement::readTracks(tfn);
        }
        catch (RelionError XE)
        {
            std::cerr << "Warning: error reading tracks in " << tfn << "\n";
            continue;
        }

        std::string imgName;
        mdts[g].getValue(EMDL_IMAGE_NAME, imgName, 0);
        imgName = outPath + "/" + imgName.substr(imgName.find_last_of("/")+1);

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

        stack.write(imgName);

        if (debug)
        {
            VtkHelper::writeTomoVTK(stack, imgName+".vtk");
        }

    } // micrographs

    double t1 = omp_get_wtime();
    double diff = t1 - t0;
    std::cout << "elapsed (total): " << diff << " sec\n";
}
