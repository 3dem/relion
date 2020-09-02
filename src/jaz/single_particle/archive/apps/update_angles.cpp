
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
#include <src/matrix2d.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/slice_helper.h>
#include <src/jaz/single_particle/spectral_helper.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/backprojection_helper.h>
#include <src/jaz/single_particle/volume_converter.h>
#include <src/jaz/single_particle/complex_io.h>
#include <src/jaz/single_particle/fftw_helper.h>
#include <src/jaz/single_particle/resampling_helper.h>
#include <src/jaz/single_particle/ctf_helper.h>
#include <src/jaz/single_particle/defocus_refinement.h>
#include <src/jaz/single_particle/magnification_refinement.h>
#include <src/jaz/single_particle/refinement_helper.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/tilt_refinement.h>
#include <src/jaz/single_particle/img_proc/image_op.h>
#include <src/jaz/single_particle/Fourier_helper.h>
#include <src/jaz/single_particle/fsc_helper.h>
#include <src/jaz/single_particle/damage_helper.h>
#include <src/jaz/single_particle/interpolation.h>
#include <src/jaz/single_particle/distribution_helper.h>
#include <src/jaz/single_particle/noise_helper.h>
#include <src/jaz/single_particle/convolution_helper.h>
#include <src/jaz/single_particle/motion_em.h>
#include <src/jaz/single_particle/local_motion_fit.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/single_particle/parallel_ft.h>

#include <omp.h>

using namespace gravis;

int main(int argc, char *argv[])
{
    std::string starFn, reconFn0, reconFn1, maskFn, outPath, inPath, fscFn;

    bool debug, applyTilt, useFsc;

    long maxMG = -1, minMG = 0;
    RFLOAT angpix, paddingFactor, beamtilt_x, beamtilt_y, deltaAngle;
    int nr_omp_threads, kmax;

    IOParser parser;

    try
    {
        parser.setCommandLine(argc, argv);

        parser.addSection("General options");

        starFn = parser.getOption("--i", "Input STAR file", "");
        reconFn0 = parser.getOption("--m0", "Reference, half 1", "");
        reconFn1 = parser.getOption("--m1", "Reference, half 2", "");
        maskFn = parser.getOption("--mask", "Reference mask", "");
        fscFn = parser.getOption("--f", "Input STAR file with the FSC of the reference", "");
        outPath = parser.getOption("--out", "Output path", "");
        inPath = parser.getOption("--img", "Path to images", "");

        deltaAngle = textToFloat(parser.getOption("--delta", "Initial angle shift (in degrees)", "1.0"));
        angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (angst/pix)", "0.0"));
        paddingFactor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));

        beamtilt_x = textToFloat(parser.getOption("--beamtilt_x", "Beamtilt in the X-direction (in mrad)", "0."));
        beamtilt_y = textToFloat(parser.getOption("--beamtilt_y", "Beamtilt in the Y-direction (in mrad)", "0."));
        applyTilt = (ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.);

        kmax = textToInteger(parser.getOption("--kmax", "Max. frequency used for alignment", "-1"));

        nr_omp_threads = textToInteger(parser.getOption("--jomp", "Number of OMP threads", "1"));
        maxMG = textToInteger(parser.getOption("--max_MG", "first micrograph index", "-1"));
        minMG = textToInteger(parser.getOption("--min_MG", "last micrograph index", "0"));

        debug = parser.checkOption("--debug", "TBD");

        if (reconFn0 == "" || reconFn1 == "")
        {
            std::cout << "An initial reconstruction for per-micrograph B-factors (--m) is required.\n";
            return 666;
        }

    }
    catch (RelionError XE)
    {
        parser.writeUsage(std::cout);
        std::cerr << XE;
        exit(1);
    }
    bool allGood = true;

    useFsc = fscFn != "";
    MetaDataTable fscMdt;

    if (useFsc)
    {
        fscMdt.read(fscFn, "fsc");

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

    if (!allGood)
    {
        return 1;
    }

    Image<RFLOAT> map0, map1, dummy;
    Projector projector0, projector1;

    try
    {
        map0.read(reconFn0);
    }
    catch (RelionError XE)
    {
        std::cout << "Unable to read map: " << reconFn0 << "\n";
        exit(1);
    }
    try
    {
        map1.read(reconFn1);
    }
    catch (RelionError XE)
    {
        std::cout << "Unable to read map: " << reconFn1 << "\n";
        exit(1);
    }

    if (map0.data.xdim != map0.data.ydim || map0.data.ydim != map0.data.zdim)
    {
        REPORT_ERROR(reconFn0 + " is not cubical.\n");
    }

    if (map1.data.xdim != map1.data.ydim || map1.data.ydim != map1.data.zdim)
    {
        REPORT_ERROR(reconFn1 + " is not cubical.\n");
    }

    if (   map0.data.xdim != map1.data.xdim
        || map0.data.ydim != map1.data.ydim
        || map0.data.zdim != map1.data.zdim)
    {
        REPORT_ERROR(reconFn0 + " and " + reconFn1 + " are of unequal size.\n");
    }

    if (maskFn != "")
    {
        std::cout << "masking references...\n";
        Image<RFLOAT> mask, maskedRef;

        try
        {
            mask.read(maskFn);
        }
        catch (RelionError XE)
        {
            std::cout << "Unable to read mask: " << maskFn << "\n";
            exit(1);
        }

        mask.read(maskFn);

        ImageOp::multiply(mask, map0, maskedRef);
        map0 = maskedRef;

        ImageOp::multiply(mask, map1, maskedRef);
        map1 = maskedRef;
    }

    const int s = map0.data.xdim;
    const int sh = s/2 + 1;

    Image<RFLOAT> imgSnr;

    if (useFsc)
    {
        RefinementHelper::computeSNR(&fscMdt, imgSnr);
    }
    else
    {
        imgSnr = Image<RFLOAT>(sh,s);
        imgSnr.data.initConstant(1.0);
    }

    std::cout << "transforming references...\n";

    projector0 = Projector(s, TRILINEAR, paddingFactor, 10, 2);
    projector0.computeFourierTransformMap(map0.data, dummy.data, map0.data.xdim);

    projector1 = Projector(s, TRILINEAR, paddingFactor, 10, 2);
    projector1.computeFourierTransformMap(map1.data, dummy.data, map1.data.xdim);

    MetaDataTable mdt0;
    mdt0.read(starFn);

    std::vector<MetaDataTable> mdts = StackHelper::splitByStack(&mdt0);

    RFLOAT Cs, lambda, kV;

    mdt0.getValue(EMDL_CTF_CS, Cs, 0);
    mdt0.getValue(EMDL_CTF_VOLTAGE, kV, 0);

    RFLOAT V = kV * 1e3;
    lambda = 12.2643247 / sqrt(V * (1.0 + V * 0.978466e-6));

    if (angpix <= 0.0)
    {
        RFLOAT mag, dstep;
        mdts[0].getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
        mdts[0].getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
        angpix = 10000 * dstep / mag;
    }

    ObservationModel obsModel(angpix);

    if (applyTilt)
    {
        obsModel = ObservationModel(angpix, Cs, kV * 1e3, beamtilt_x, beamtilt_y);
    }

    const long gc = maxMG >= 0? maxMG : mdts.size()-1;
    const long g0 = minMG;

    std::cout << "mg range: " << g0 << ".." << gc << "\n";

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    double t0 = omp_get_wtime();

    const bool quadratic = true;

    MetaDataTable mdtAll;
    mdtAll.reserve(mdt0.numberOfObjects());

    for (long g = g0; g <= gc; g++)
    {
        std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

        std::stringstream stsg;
        stsg << g;

        const int pc = mdts[g].numberOfObjects();

        std::vector<Image<Complex>> obsF
                = StackHelper::loadStackFS(&mdts[g], inPath, nr_omp_threads, &fts);

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (long p = 0; p < pc; p++)
        {
            int randSubset;
            mdts[g].getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
            randSubset -= 1;

            if (quadratic)
            {
                Matrix2D<RFLOAT> A(27,10);
                Matrix1D<RFLOAT> b(27);

                for (int rot = -1; rot <= 1; rot++)
                for (int tilt = -1; tilt <= 1; tilt++)
                for (int psi = -1; psi <= 1; psi++)
                {
                    Image<Complex> pred;

                    if (randSubset == 0)
                    {
                        pred = obsModel.predictObservation(
                            projector0, mdts[g], p, true, true,
                            rot*deltaAngle, tilt*deltaAngle, psi*deltaAngle);
                    }
                    else
                    {
                        pred = obsModel.predictObservation(
                            projector1, mdts[g], p, true, true,
                            rot*deltaAngle, tilt*deltaAngle, psi*deltaAngle);
                    }

                    const double index = 9*(rot+1) + 3*(tilt+1) + (psi+1);

                    b(index) = 0.0;

                    for (int y = 0; y < s; y++)
                    for (int x = 0; x < sh; x++)
                    {
                        double yy = y < sh? y : y - s;
                        double r = sqrt(x*x + yy*yy);
                        if (r > kmax) continue;

                        b(index) += imgSnr(y,x) * (pred(y,x) - obsF[p](y,x)).norm();
                    }

                    A(index, 0) = rot*rot;
                    A(index, 1) = 2.0*rot*tilt;
                    A(index, 2) = 2.0*rot*psi;
                    A(index, 3) = 2.0*rot;

                    A(index, 4) = tilt*tilt;
                    A(index, 5) = 2.0*tilt*psi;
                    A(index, 6) = 2.0*tilt;

                    A(index, 7) = psi*psi;
                    A(index, 8) = 2.0*psi;

                    A(index, 9) = 1.0;
                }

                const double tol = 1e-20;
                Matrix1D<RFLOAT> x(10);
                solve(A, b, x, tol);

                d3Matrix C3(x(0), x(1), x(2),
                            x(1), x(4), x(5),
                            x(2), x(5), x(7));

                d3Vector d(x(3), x(6), x(8));

                d3Matrix C3i = C3;
                C3i.invert();

                d3Vector min = -C3i * d;

                if (debug) std::cout << p << ": " << min*deltaAngle << "\n";

                if (min.length() > 1.0) min /= min.length();

                double rot, tilt, psi;

                mdts[g].getValue(EMDL_ORIENT_ROT, rot, p);
                mdts[g].getValue(EMDL_ORIENT_TILT, tilt, p);
                mdts[g].getValue(EMDL_ORIENT_PSI, psi, p);

                rot += min[0]*deltaAngle;
                tilt += min[1]*deltaAngle;
                psi += min[2]*deltaAngle;

                mdts[g].setValue(EMDL_ORIENT_ROT, rot, p);
                mdts[g].setValue(EMDL_ORIENT_TILT, tilt, p);
                mdts[g].setValue(EMDL_ORIENT_PSI, psi, p);

            }
            else
            {

            }
        }

        mdtAll.append(mdts[g]);
    }

    mdtAll.write(outPath);

    double t1 = omp_get_wtime();
    double diff = t1 - t0;
    std::cout << "elapsed (total): " << diff << " sec\n";
}
