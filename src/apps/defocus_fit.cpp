
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

class DefocusFit : public RefinementProgram
{
    public:

        RFLOAT defocusRange, Cs;
        bool fitAstigmatism, diag, ownCs;

        void readMoreOptions(IOParser& parser, int argc, char *argv[]);
        int _init();
        int _run();
};

int main(int argc, char *argv[])
{
    DefocusFit mt;

    int rc0 = mt.init(argc, argv);
    if (rc0 != 0) return rc0;

    int rc1 = mt.run();
    if (rc1 != 0) return rc1;
}

void DefocusFit::readMoreOptions(IOParser& parser, int argc, char *argv[])
{
    fitAstigmatism = parser.checkOption("--astig", "Estimate independent astigmatism for each particle");
    diag = parser.checkOption("--diag", "Write out defocus errors");
    Cs = textToFloat(parser.getOption("--Cs", "Spherical aberration", "-1"));
    ownCs = Cs != -1;
    defocusRange = textToFloat(parser.getOption("--range", "Defocus scan range (in A)", "2000."));
}

int DefocusFit::_init()
{
    return 0;
}

int DefocusFit::_run()
{
    if (diag) std::cout << "running in diagnostic mode\n";
    if (ownCs) std::cout << "assuming Cs = " << Cs << "\n";

    double t0 = omp_get_wtime();

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    MetaDataTable mdtAll;
    mdtAll.reserve(mdt0.numberOfObjects());

    const long gc = maxMG >= 0? maxMG+1 : mdts.size();

    for (long g = minMG; g < gc; g++)
    {
        std::stringstream stsg;
        stsg << g;

        std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

        const int pc = mdts[g].numberOfObjects();

        std::vector<Image<Complex> > obsF;

        obsF = StackHelper::loadStackFS(&mdts[g], imgPath, nr_omp_threads, &fts);

        if (applyTilt)
        {
            if (nr_omp_threads > 1)
            {
                obsF = StackHelper::applyBeamTiltPar(obsF, Cs, lambda, angpix, beamtilt_x, beamtilt_y, nr_omp_threads);
            }
            else
            {
                obsF = StackHelper::applyBeamTilt(obsF, Cs, lambda, angpix, beamtilt_x, beamtilt_y);
            }
        }

        if (diag)
        {
            std::ofstream ofst(outPath+"_diag_m"+stsg.str()+".dat");
            std::ofstream ofsto(outPath+"_diag_m"+stsg.str()+"_opt.dat");

            for (long p = 0; p < pc; p++)
            {
                std::cout << "    " << p << " / " << pc << "\n";

                CTF ctf0;
                ctf0.read(mdts[g], mdts[g], p);

                if (ownCs)
                {
                    ctf0.Cs = Cs;
                    ctf0.initialise();
                }

                int randSubset;
                mdts[g].getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
                randSubset -= 1;

                Image<Complex> pred;

                pred = obsModel.predictObservation(
                    projectors[randSubset], mdts[g], p, false, true);

                std::vector<d2Vector> cost = DefocusRefinement::diagnoseDefocus(
                    pred, obsF[p], freqWeight,
                    ctf0, angpix, defocusRange, 100, nr_omp_threads);

                double cMin = cost[0][1];
                double dOpt = cost[0][0];

                for (int i = 0; i < cost.size(); i++)
                {
                    ofst << cost[i][0] << " " << cost[i][1] << "\n";

                    if (cost[i][1] < cMin)
                    {
                        cMin = cost[i][1];
                        dOpt = cost[i][0];
                    }
                }

                ofsto << dOpt << " " << cMin << "\n";

                ofst << "\n";
            }
        }
        else
        {
            #pragma omp parallel for num_threads(nr_omp_threads)
            for (long p = 0; p < pc; p++)
            {
                std::stringstream stsp;
                stsp << p;

                CTF ctf0;
                ctf0.read(mdts[g], mdts[g], p);

                if (ownCs)
                {
                    ctf0.Cs = Cs;
                    ctf0.initialise();
                }

                int randSubset;
                mdts[g].getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
                randSubset -= 1;

                Image<Complex> pred;

                pred = obsModel.predictObservation(
                    projectors[randSubset], mdts[g], p, false, true);

                CTF ctf(ctf0);

                if (fitAstigmatism)
                {
                    double u, v, phi;
                    DefocusRefinement::findAstigmatismNM(pred, obsF[p], freqWeight, ctf0, angpix, &u, &v, &phi);

                    mdts[g].setValue(EMDL_CTF_DEFOCUSU, u, p);
                    mdts[g].setValue(EMDL_CTF_DEFOCUSV, v, p);
                    mdts[g].setValue(EMDL_CTF_DEFOCUS_ANGLE, phi, p);

                    ctf.DeltafU = u;
                    ctf.DeltafV = v;
                    ctf.azimuthal_angle = phi;
                    ctf.initialise();
                }
                else
                {
                    double u, v;
                    DefocusRefinement::findDefocus1D(pred, obsF[p], freqWeight, ctf0, angpix, &u, &v, defocusRange);

                    mdts[g].setValue(EMDL_CTF_DEFOCUSU, u, p);
                    mdts[g].setValue(EMDL_CTF_DEFOCUSV, v, p);

                    ctf.DeltafU = u;
                    ctf.DeltafV = v;
                    ctf.initialise();
                }

                if (ownCs)
                {
                    mdts[g].setValue(EMDL_CTF_CS, Cs, p);
                }
            }

            mdtAll.append(mdts[g]);
        }
    }

    double t1 = omp_get_wtime();

    std::cout << "elapsed: " << (t1 - t0) << "s \n";

    if (!diag)
    {
        mdtAll.write(outPath);
    }

    return 0;
}
