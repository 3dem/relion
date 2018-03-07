
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
            std::string trackFn, fccFn;

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
    trackFn = parser.getOption("--t", "Input tracks");
    fccFn = parser.getOption("--cc", "Input MRC file of per-micrograph Fourier cross-correlations");

    // @TODO: change these to Angstrom
    k0 = textToInteger(parser.getOption("--k0", "Min. frequency used in B-factor fit", "20"));
    k1 = textToInteger(parser.getOption("--k1", "Max. frequency used in B-factor fit", "100"));

    return 0;
}

int FrameRecomb::_init()
{
    return 0;
}

int FrameRecomb::_run()
{
    Image<RFLOAT> fcc;

    std::vector<Image<RFLOAT>> freqWeights;

    fcc.read(fccFn);

    sh = fcc.data.xdim;
    s = 2 * (sh-1);
    fc = fcc.data.ydim;

    std::pair<std::vector<d2Vector>,std::vector<double>> bkFacs
            = DamageHelper::fitBkFactors(fcc, k0, k1);

    freqWeights = DamageHelper::computeWeights(bkFacs.first, sh);

    if (debug)
    {
        Image<RFLOAT> bfacFit = DamageHelper::renderBkFit(bkFacs, sh, fc);
        Image<RFLOAT> bfacFitNoScale = DamageHelper::renderBkFit(bkFacs, sh, fc, true);

        VtkHelper::writeVTK(bfacFit, "bfacs/glob_Bk-fit.vtk");
        VtkHelper::writeVTK(bfacFitNoScale, "bfacs/glob_Bk-fit_noScale.vtk");
        VtkHelper::writeVTK(fcc, "bfacs/glob_Bk-data.vtk");
        VtkHelper::write(freqWeights, "bfacs/freqWeights.vtk");
    }

    std::string name, fullName, movieName;
    mdts[0].getValue(EMDL_IMAGE_NAME, fullName, 0);
    mdts[0].getValue(EMDL_MICROGRAPH_NAME, movieName, 0);
    name = fullName.substr(fullName.find("@")+1);

    if (preextracted)
    {
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

        fc = stackSize / pc0;

        s = stack0.data.xdim;
        sh = s/2 + 1;
    }
    else
    {
        std::vector<std::vector<Image<Complex>>> movie = StackHelper::extractMovieStackFS(
            &mdts[0], meta_path, imgPath, angpix, angpix, movie_angpix, s,
            nr_omp_threads, false, hotCutoff, debug);

        fc = movie[0].size();
    }

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
            if (preextracted)
            {
                movie = StackHelper::loadMovieStackFS(
                    &mdts[g], imgPath, false, nr_omp_threads, &fts);
            }
            else
            {
                movie = StackHelper::extractMovieStackFS(
                    &mdts[g], meta_path, imgPath,
                    angpix, angpix, movie_angpix, s,
                    nr_omp_threads, true, hotCutoff, debug);

                #pragma omp parallel for num_threads(nr_omp_threads)
                for (int p = 0; p < pc; p++)
                {
                    StackHelper::varianceNormalize(movie[p], false);
                }
            }
        }
        catch (RelionError XE)
        {
            std::cerr << "warning: unable to load micrograph #" << (g+1) << "\n";
            continue;
        }

        const int sh = movie[0][0]().xdim;
        const int s = movie[0][0]().ydim;

        std::string tfn = trackFn + "_mg" + stsg.str() + "_tracks.dat";
        std::ifstream trackIn(tfn);

        std::vector<std::vector<d2Vector>> shift(pc);

        for (int p = 0; p < pc; p++)
        {
            if (!trackIn.good())
            {
                std::cerr << "relion_recombine: error reading tracks in " << tfn << "\n";
                return 5;
            }

            shift[p] = std::vector<d2Vector>(fc);

            char dummy[4069];
            trackIn.getline(dummy, 4069);

            for (int f = 0; f < fc; f++)
            {
                char dummy[4069];
                trackIn.getline(dummy, 4069);

                std::istringstream sts(dummy);

                sts >> shift[p][f].x;
                sts >> shift[p][f].y;
            }

            trackIn.getline(dummy, 4069);
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
