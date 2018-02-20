
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
#include <src/jaz/parallel_ft.h>

#include <omp.h>

using namespace gravis;

class AberrationPlot : public RefinementProgram
{
    public:

        AberrationPlot();

        int readMoreOptions(IOParser& parser, int argc, char *argv[]);
        int _init();
        int _run();
};

AberrationPlot :: AberrationPlot()
:   RefinementProgram(true)
{
}

int main(int argc, char *argv[])
{
    AberrationPlot tf;

    int rc0 = tf.init(argc, argv);
    if (rc0 != 0) return rc0;

    int rc1 = tf.run();
    if (rc1 != 0) return rc1;
}

int AberrationPlot::readMoreOptions(IOParser& parser, int argc, char *argv[])
{
    return 0;
}

int AberrationPlot::_init()
{
    return 0;
}

int AberrationPlot::_run()
{
    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    double t0 = omp_get_wtime();

    const bool differential = false;

    if (differential)
    {
        std::vector<Image<double>>
                A(nr_omp_threads, Image<double>(sh,s)),
                b(nr_omp_threads, Image<double>(sh,s));

        const double as = (double)s * angpix;

        for (long g = minMG; g <= gc; g++)
        {
            std::stringstream stsg;
            stsg << g;

            std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

            const int pc = mdts[g].numberOfObjects();

            std::vector<Image<Complex> > pred;
            std::vector<Image<Complex> > obsF;

            pred = obsModel.predictObservations(projectors[0], mdts[g], false, true, nr_omp_threads);
            obsF = StackHelper::loadStackFS(&mdts[g], imgPath, nr_omp_threads, &fts);

            #pragma omp parallel for num_threads(nr_omp_threads)
            for (long p = 0; p < pc; p++)
            {
                int t = omp_get_thread_num();

                CTF ctf0;
                ctf0.read(mdts[g], mdts[g], p);
                //ctf0.Cs = 0.0;
                ctf0.initialise();

                for (int y = 0; y < s;  y++)
                for (int x = 0; x < sh; x++)
                {
                    const double xf = x;
                    const double yf = y < sh? y : y - s;
                    const double gamma_i = ctf0.getGamma(xf/as, yf/as);
                    const double cg = cos(gamma_i);
                    const double sg = sin(gamma_i);

                    Complex zobs = obsF[p](y,x);
                    Complex zprd = pred[p](y,x);

                    double zz = zobs.real*zprd.real + zobs.imag*zprd.imag;
                    double nr = zprd.norm();

                    A[t](y,x) += nr*cg*cg;
                    b[t](y,x) += cg*(sg*nr+zz);
                }
            }
        }

        for (int t = 1; t < nr_omp_threads; t++)
        {
            for (int y = 0; y < s;  y++)
            for (int x = 0; x < sh; x++)
            {
                A[0](y,x) += A[t](y,x);
                b[0](y,x) += b[t](y,x);
            }
        }

        Image<RFLOAT> dgamma(sh,s);

        for (int y = 0; y < s;  y++)
        for (int x = 0; x < sh; x++)
        {
            if (A[0](y,x) != 0.0)
            {
                dgamma(y,x) = b[0](y,x)/A[0](y,x);
            }
        }

        VtkHelper::writeVTK(dgamma, outPath+"_dgamma.vtk");
    }
    else
    {
        std::vector<Image<double>>
            Axx(nr_omp_threads, Image<double>(sh,s)),
            Axy(nr_omp_threads, Image<double>(sh,s)),
            Ayy(nr_omp_threads, Image<double>(sh,s)),
            bx(nr_omp_threads, Image<double>(sh,s)),
            by(nr_omp_threads, Image<double>(sh,s));

        const double as = (double)s * angpix;

        for (long g = minMG; g <= gc; g++)
        {
            std::stringstream stsg;
            stsg << g;

            std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

            const int pc = mdts[g].numberOfObjects();

            std::vector<Image<Complex> > pred;
            std::vector<Image<Complex> > obsF;

            pred = obsModel.predictObservations(projectors[0], mdts[g], false, true, nr_omp_threads);
            obsF = StackHelper::loadStackFS(&mdts[g], imgPath, nr_omp_threads, &fts);

            #pragma omp parallel for num_threads(nr_omp_threads)
            for (long p = 0; p < pc; p++)
            {
                int t = omp_get_thread_num();

                CTF ctf0;
                ctf0.read(mdts[g], mdts[g], p);
                //ctf0.Cs = 0.0;
                ctf0.initialise();

                for (int y = 0; y < s;  y++)
                for (int x = 0; x < sh; x++)
                {
                    const double xf = x;
                    const double yf = y < sh? y : y - s;
                    const double gamma_i = ctf0.getGamma(xf/as, yf/as);
                    const double cg = cos(gamma_i);
                    const double sg = sin(gamma_i);

                    Complex zobs = obsF[p](y,x);
                    Complex zprd = pred[p](y,x);

                    double zz = zobs.real*zprd.real + zobs.imag*zprd.imag;
                    double nr = zprd.norm();

                    Axx[t](y,x) += nr*sg*sg;
                    Axy[t](y,x) += nr*cg*sg;
                    Ayy[t](y,x) += nr*cg*cg;

                    bx[t](y,x) += zz*sg;
                    by[t](y,x) += zz*cg;
                }
            }
        }

        for (int t = 1; t < nr_omp_threads; t++)
        {
            for (int y = 0; y < s;  y++)
            for (int x = 0; x < sh; x++)
            {
                Axx[0](y,x) += Axx[t](y,x);
                Axy[0](y,x) += Axy[t](y,x);
                Ayy[0](y,x) += Ayy[t](y,x);

                bx[0](y,x) += bx[t](y,x);
                by[0](y,x) += by[t](y,x);
            }
        }

        Image<RFLOAT> xx(sh,s), xy(sh,s), phase(sh,s);

        for (int y = 0; y < s;  y++)
        for (int x = 0; x < sh; x++)
        {
            d2Matrix A(
                Axx[0](y,x), Axy[0](y,x),
                Axy[0](y,x), Ayy[0](y,x));

            d2Vector b(bx[0](y,x), by[0](y,x));

            double det = A(0,0)*A(1,1) - A(1,0)*A(0,1);

            if (det != 0.0)
            {
                d2Matrix Ai = A;
                Ai.invert();

                d2Vector opt = Ai*b;

                xx(y,x) = opt.x;
                xy(y,x) = opt.y;

                phase(y,x) = std::abs(opt.x) > 0.0? atan2(opt.y, opt.x) : 0.0;
            }
        }

        VtkHelper::writeVTK(xx, outPath+"_cos.vtk");
        VtkHelper::writeVTK(xy, outPath+"_sin.vtk");
        VtkHelper::writeVTK(phase, outPath+"_phase.vtk");
    }

    return 0;
}
