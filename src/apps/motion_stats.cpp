
#include <src/args.h>
#include <src/metadata_table.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/motion/motion_helper.h>

using namespace gravis;

int main(int argc, char *argv[])
{
    IOParser parser;
    std::string starFn, tracksPath;
    int fc;
    double angpix, fdose;

    parser.setCommandLine(argc, argv);
    parser.addSection("General options");
    starFn = parser.getOption("--i", "Input STAR file with a list of particles");
    tracksPath = parser.getOption("--t", "Path to particle trajectories");
    fdose = textToFloat(parser.getOption("--fdose", "Electron dose per frame (e^-/A^2)"));
    angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (Angst/pix)"));

    fc = textToInteger(parser.getOption("--fc", "Frame count"));

    if (parser.checkForErrors()) return 1;

    MetaDataTable mdt0;
    mdt0.read(starFn);

    std::vector<MetaDataTable> mdts = StackHelper::splitByStack(&mdt0);

    const int mgc = mdts.size();

    std::vector<double> varPerFrame(fc, 0.0);
    double var = 0.0;
    long tpc = 0;

    for (int m = 0; m < mgc; m++)
    {
        std::stringstream stsg;
        stsg << m;

        if (m%100 == 0) std::cout << "micrograph " << (m+1) << " / " << mgc << "\n";

        const int pc = mdts[m].numberOfObjects();
        tpc += pc;

        std::string tag;
        mdts[m].getValue(EMDL_IMAGE_NAME, tag, 0);
        tag = tag.substr(0,tag.find_last_of('.'));
        tag = tag.substr(tag.find_first_of('@')+1);

        std::string tfn = tracksPath + "/" + tag + "_tracks.star";

        std::vector<std::vector<d2Vector>> shift;

        try
        {
            shift = MotionHelper::readTracks(tfn);
        }
        catch (RelionError XE)
        {
            std::cerr << "Warning: error reading tracks in " << tfn << "\n";
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
        std::cout << fdose*f << " " << sqrt(varPerFrame[f]/tpc) << "\n";
    }

    std::cout << "\ntotal: " << sqrt(var/(tpc*fc)) << " A\n";
}
