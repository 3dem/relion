
#include <src/args.h>
#include <src/metadata_table.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/motion_refinement.h>

using namespace gravis;

int main(int argc, char *argv[])
{
    IOParser parser;
    std::string starFn, tracksPath;
    int fc;
    double angpix;

    parser.setCommandLine(argc, argv);
    parser.addSection("General options");
    starFn = parser.getOption("--i", "Input STAR file with a list of particles");
    tracksPath = parser.getOption("--t", "Path to particle trajectories");
    angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (Angst/pix)"));
    fc = textToInteger(parser.getOption("--fc", "Frame count"));

    if (parser.checkForErrors()) return 1;

    MetaDataTable mdt0;
    mdt0.read(starFn);

    std::vector<MetaDataTable> mdts = StackHelper::splitByStack(&mdt0);

    const int mgc = mdts.size();

    std::vector<double> varPerFrame(fc, 0.0);
    double var = 0.0;
    int tpc;

    for (int m = 0; m < mgc; m++)
    {
        std::stringstream stsg;
        stsg << m;

        if (m%100 == 0) std::cout << "micrograph " << (m+1) << " / " << mgc << "\n";

        const int pc = mdts[m].numberOfObjects();
        tpc += pc;

        std::string tfn = tracksPath + "_mg" + stsg.str() + "_tracks.dat";
        std::vector<std::vector<d2Vector>> shift;

        try
        {
            shift = MotionRefinement::readTrack(tfn, pc, fc);
        }
        catch (RelionError XE)
        {
            std::cerr << "warning: error reading tracks in " << tfn << "\n";
            continue;
        }

        for (int p = 0; p < pc; p++)
        for (int f = 0; f < fc-1; f++)
        {
            d2Vector v = angpix * (shift[p][f+1] - shift[p][f]);

            var += v.norm2();
            varPerFrame[f] += v.norm2();
        }
    }

    for (int f = 0; f < fc-1; f++)
    {
        std::cout << f << ": " << sqrt(varPerFrame[f]/tpc) << " A\n";
    }

    std::cout << "\ntotal: " << sqrt(var/(tpc*fc)) << " A\n";
}
