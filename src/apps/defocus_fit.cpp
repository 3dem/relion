
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

        RFLOAT defocusRange, kmin;
        bool fitAstigmatism, noGlobAstig, diag, fitPhase, fitCs, globOnly;

        int readMoreOptions(IOParser& parser, int argc, char *argv[]);
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

int DefocusFit::readMoreOptions(IOParser& parser, int argc, char *argv[])
{
    fitAstigmatism = parser.checkOption("--astig", "Estimate independent astigmatism for each particle");
    noGlobAstig = parser.checkOption("--no_glob_astig", "Skip per-micrograph astigmatism estimation");
    fitPhase = !parser.checkOption("--no_phase", "Skip per-micrograph phase shift (amplitude contrast) estimation");
    diag = parser.checkOption("--diag", "Write out defocus errors");
    globOnly = parser.checkOption("--glob", "Only perform per-micrograph fit");
    defocusRange = textToFloat(parser.getOption("--range", "Defocus scan range (in A)", "2000."));
    kmin = textToFloat(parser.getOption("--kmin", "Inner freq. threshold [Angst]", "20.0"));
    fitCs = parser.checkOption("--fit_Cs", "Fit spherical aberration (per micrograph)");

    return 0;
}

int DefocusFit::_init()
{
    return 0;
}

int DefocusFit::_run()
{
    if (diag) std::cout << "running in diagnostic mode\n";

    double t0 = omp_get_wtime();

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    MetaDataTable mdtAll;
    mdtAll.reserve(mdt0.numberOfObjects());

    for (long g = minMG; g <= gc; g++)
    {
        std::stringstream stsg;
        stsg << g;

        std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

        const int pc = mdts[g].numberOfObjects();

        std::vector<Image<Complex> > obsF;

        obsF = StackHelper::loadStackFS(&mdts[g], imgPath, nr_omp_threads, &fts);

        std::vector<Image<Complex>> preds(pc);

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (long p = 0; p < pc; p++)
        {
            int randSubset;
            mdts[g].getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
            randSubset -= 1;

            preds[p] = obsModel.predictObservation(
                projectors[randSubset], mdts[g], p, false, true);
        }

        for (int y = 0; y < s; y++)
        for (int x = 0; x < sh; x++)
        {
            double xx = x;
            double yy = y <= sh? y : y - s;
            double r = sqrt(xx*xx + yy*yy);

            if (r == 0 || 2.0*sh*angpix/r > kmin)
            {
                freqWeight(y,x) = 0.0;
            }
        }

        if (!noGlobAstig)
        {
            CTF ctf0;
            ctf0.read(mdts[g], mdts[g], 0);

            Image<RFLOAT> dataVis;

            if (diag)
            {
                Image<RFLOAT> dotp0(sh,s), cc(sh,s), wgh0(sh,s), wgh1(sh,s), dotp0_full(s,s);
                dotp0.data.initZeros();
                wgh0.data.initZeros();
                wgh1.data.initZeros();

                for (long p = 0; p < pc; p++)
                {
                    for (long y = 0; y < s; y++)
                    for (long x = 0; x < sh; x++)
                    {
                        Complex vx = DIRECT_A2D_ELEM(preds[p].data, y, x);
                        const Complex vy = DIRECT_A2D_ELEM(obsF[p].data, y, x);

                        dotp0(y,x) += vy.real*vx.real + vy.imag*vx.imag;
                        wgh0(y,x) += vx.norm();
                        wgh1(y,x) += vy.norm();
                    }
                }

                for (long y = 0; y < s; y++)
                for (long x = 0; x < sh; x++)
                {
                    double nrm = sqrt(wgh0(y,x) * wgh1(y,x));
                    cc(y,x) = nrm > 0.0? 10.0 * dotp0(y,x) / nrm : 0.0;
                }

                FftwHelper::decenterDouble2D(cc.data, dotp0_full.data);

                VtkHelper::writeVTK(dotp0_full, outPath+"_astig_data_m"+stsg.str()+".vtk");
                dataVis = FilterHelper::polarBlur(dotp0_full, 10.0);
                VtkHelper::writeVTK(dataVis, outPath+"_astig_data_m"+stsg.str()+"_blurred.vtk");

                Image<RFLOAT> ctfFit(s,s);
                ctf0.getCenteredImage(ctfFit.data, angpix, false, false, false, false);
                VtkHelper::writeVTK(ctfFit, outPath+"_astig0_m"+stsg.str()+".vtk");

                ctfFit.data.xinit = 0;
                ctfFit.data.yinit = 0;

                Image<RFLOAT> vis = FilterHelper::sectorBlend(dataVis, ctfFit, 12);
                VtkHelper::writeVTK(vis, outPath+"_astig_data_m"+stsg.str()+"_vis0.vtk");
            }

            RFLOAT u, v, phi, phase, newCs;

            mdts[g].getValue(EMDL_CTF_PHASESHIFT, phase, 0);
            mdts[g].getValue(EMDL_CTF_CS, newCs, 0);

            if (fitCs)
            {
                std::cout << "initial phi and Cs: " << phase << ", " << newCs << "\n";

                DefocusRefinement::findAstigmatismPhaseAndCsNM(
                        preds, obsF, freqWeight, ctf0, angpix, &u, &v, &phi, &phase, &newCs);

                std::cout << "final phi and Cs: " << phase << ", " << newCs << "\n";
            }
            else if (fitPhase)
            {
                std::cout << "initial phase shift: " << phase << "\n";

                DefocusRefinement::findAstigmatismAndPhaseNM(
                        preds, obsF, freqWeight, ctf0, angpix, &u, &v, &phi, &phase);

                std::cout << "final phase shift: " << phase << "\n";
            }
            else
            {
                DefocusRefinement::findAstigmatismNM(preds, obsF, freqWeight, ctf0, angpix, &u, &v, &phi);
            }

            for (long p = 0; p < pc; p++)
            {
                mdts[g].setValue(EMDL_CTF_DEFOCUSU, u, p);
                mdts[g].setValue(EMDL_CTF_DEFOCUSV, v, p);
                mdts[g].setValue(EMDL_CTF_DEFOCUS_ANGLE, phi, p);

                if (fitPhase) mdts[g].setValue(EMDL_CTF_PHASESHIFT, phase, p);
                if (fitCs) mdts[g].setValue(EMDL_CTF_CS, newCs, p);
            }

            if (diag)
            {
                CTF ctf1;
                ctf1.read(mdts[g], mdts[g], 0);

                Image<RFLOAT> ctfFit(s,s);
                ctf1.getCenteredImage(ctfFit.data, angpix, false, false, false, false);
                VtkHelper::writeVTK(ctfFit, outPath+"_astig1_m"+stsg.str()+".vtk");
                ctfFit.data.xinit = 0;
                ctfFit.data.yinit = 0;

                Image<RFLOAT> vis = FilterHelper::sectorBlend(dataVis, ctfFit, 12);
                VtkHelper::writeVTK(vis, outPath+"_astig_data_m"+stsg.str()+"_vis1.vtk");
            }
        }

        if (globOnly) continue;

        if (diag)
        {
            std::ofstream ofst(outPath+"_diag_m"+stsg.str()+".dat");
            std::ofstream ofsto(outPath+"_diag_m"+stsg.str()+"_opt.dat");

            for (long p = 0; p < pc; p++)
            {
                std::cout << "    " << p << " / " << pc << "\n";

                CTF ctf0;
                ctf0.read(mdts[g], mdts[g], p);

                std::vector<d2Vector> cost = DefocusRefinement::diagnoseDefocus(
                    preds[p], obsF[p], freqWeight,
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

                CTF ctf(ctf0);

                if (fitAstigmatism)
                {
                    double u, v, phi;
                    DefocusRefinement::findAstigmatismNM(
                                preds[p], obsF[p], freqWeight, ctf0,
                                angpix, &u, &v, &phi);

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
                    DefocusRefinement::findDefocus1D(
                                preds[p], obsF[p], freqWeight, ctf0,
                                angpix, &u, &v, defocusRange);

                    mdts[g].setValue(EMDL_CTF_DEFOCUSU, u, p);
                    mdts[g].setValue(EMDL_CTF_DEFOCUSV, v, p);

                    ctf.DeltafU = u;
                    ctf.DeltafV = v;
                    ctf.initialise();
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
