
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
#include <src/jaz/legacy/vtk_helper.h>
#include <src/jaz/legacy/slice_helper.h>
#include <src/jaz/legacy/spectral_helper.h>
#include <src/jaz/legacy/img_proc/filter_helper.h>
#include <src/jaz/legacy/volume_converter.h>
#include <src/jaz/legacy/complex_io.h>
#include <src/jaz/legacy/fftw_helper.h>
#include <src/jaz/legacy/resampling_helper.h>
#include <src/jaz/legacy/ctf_helper.h>
#include <src/jaz/legacy/refinement_helper.h>
#include <src/jaz/legacy/stack_helper.h>
#include <src/jaz/legacy/img_proc/image_op.h>
#include <src/jaz/legacy/Fourier_helper.h>
#include <src/jaz/legacy/fsc_helper.h>
#include <src/jaz/legacy/damage_helper.h>
#include <src/jaz/legacy/interpolation.h>
#include <src/jaz/legacy/distribution_helper.h>
#include <src/jaz/legacy/noise_helper.h>
#include <src/jaz/legacy/convolution_helper.h>
#include <src/jaz/legacy/local_motion_fit.h>

#include <omp.h>

using namespace gravis;

int main(int argc, char *argv[])
{
    std::string starFn;
    int s, threads;

    IOParser parser;

    try
    {
        parser.setCommandLine(argc, argv);

        parser.addSection("General options");

        starFn = parser.getOption("--i", "Input *.star file");
        s = textToInteger(parser.getOption("--s", "Image size"));
        threads = textToInteger(parser.getOption("--j", "Number of threads"));

        parser.checkForErrors();
    }
    catch (RelionError XE)
    {
        parser.writeUsage(std::cout);
        std::cerr << XE;
        return RELION_EXIT_FAILURE;
    }

    MetaDataTable mdt;
    mdt.read(starFn);

    RFLOAT mag, dstep;
    mdt.getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
    mdt.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
    double angpix = 10000 * dstep / mag;

    int sh = s/2 + 1;

    std::vector<int> percentages = {50, 60, 70, 80, 90, 100, 110};
    const int perCnt = percentages.size();

    std::vector<std::vector<Image<RFLOAT>>> countN(perCnt);

    for (int c = 0; c < perCnt; c++)
    {
        countN[c] = std::vector<Image<RFLOAT>>(threads);

        for (int t = 0; t < threads; t++)
        {
            countN[c][t] = Image<RFLOAT>(sh,s);
        }
    }

    const int pc = mdt.numberOfObjects();


    #pragma omp parallel for num_threads(threads)
    for (int p = 0; p < pc; p++)
    {
        if (p%1000 == 0) std::cout << p << " / " << pc << "\n";

        int th = omp_get_thread_num();

        CTF ctf;
        ctf.read(mdt, mdt, p);

        RFLOAT as = (RFLOAT)s * angpix;

        for (long int i = 0; i < s;  i++) \
        for (long int j = 0; j < sh; j++)
        {
            const int x = j;
            const int y = i < sh? i : i - s;

            RFLOAT cf = ctf.getCtfFreq(x/as, y/as) / (as * PI);
            double cfa = std::abs(cf);

            for (int c = 0; c < perCnt; c++)
            {
                if (cfa < percentages[c]/100.0)
                {
                    countN[c][th](i,j) += 1.0/pc;
                }
            }
        }
    }
	
	{
		std::string command = " mkdir -p ctf_test";
		int ret = system(command.c_str());
	}

    for (int c = 0; c < perCnt; c++)
    {
        for (int t = 1; t < threads; t++)
        {
            for (long int i = 0; i < s;  i++) \
            for (long int j = 0; j < sh; j++)
            {
                countN[c][0](i,j) += countN[c][t](i,j);
            }
        }

        std::stringstream stsc;
        stsc << percentages[c];

        VtkHelper::writeVTK(countN[c][0], "ctf_test/below_"+stsc.str()+"%_nyq.vtk");
    }
}
