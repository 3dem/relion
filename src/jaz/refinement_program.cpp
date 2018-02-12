#include <src/jaz/refinement_program.h>

#include <src/jaz/image_op.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/vtk_helper.h>

RefinementProgram::RefinementProgram(bool singleReference, bool doesMovies)
:   singleReference(singleReference),
    optStar(false),
    noStar(false),
    optReference(false),
    noReference(false),
    noTilt(false),
    doesMovies(doesMovies)
{
}

int RefinementProgram::init(int argc, char *argv[])
{
    IOParser parser;

    try
    {
        parser.setCommandLine(argc, argv);

        parser.addSection("General options");

        starFn = parser.getOption("--i", "Input STAR file", optStar? "" : "NULL");

        if (singleReference)
        {
            reconFn0 = parser.getOption("--m", "Reference map", optReference? "" : "NULL");
        }
        else
        {
            reconFn0 = parser.getOption("--m1", "Reference map, half 1", optReference? "" : "NULL");
            reconFn1 = parser.getOption("--m2", "Reference map, half 2", optReference? "" : "NULL");
        }

        maskFn = parser.getOption("--mask", "Reference mask", "");
        fscFn = parser.getOption("--f", "Input STAR file with the FSC of the reference", "");
        outPath = parser.getOption("--out", "Output path");

        if (doesMovies)
        {
            imgPath = parser.getOption("--mov", "Path to movies", "");
            preextracted = parser.checkOption("--preex", "Preextracted movie stacks");
            meta_path = parser.getOption("--meta", "Path to per-movie metadata star files", "");
            nogain = parser.checkOption("--nogain", "Ignore gain reference");

            bin = textToInteger(parser.getOption("--bin", "Binning level for optimization and output (e.g. 2 for 2x2)", "1"));
            coords_bin = textToInteger(parser.getOption("--cbin", "Binning level of input coordinates", "1"));
            movie_bin = textToInteger(parser.getOption("--mbin", "Binning level of input movies", "1"));
            bin_type_str = parser.getOption("--bintype", "Binning method (box, gauss or fourier)", "fourier");
        }
        else
        {
            imgPath = parser.getOption("--img", "Path to images", "");
        }

        angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (angst/pix) - read from STAR file by default", "0.0"));
        Cs = textToFloat(parser.getOption("--Cs", "Spherical aberration - read from STAR file by default", "-1"));
        kV = textToFloat(parser.getOption("--kV", "Electron energy (keV) - read from STAR file by default", "-1"));

        paddingFactor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));

        if (noTilt)
        {
            applyTilt = false;
        }
        else
        {
            beamtilt_x = textToFloat(parser.getOption("--beamtilt_x", "Beamtilt in X-direction (in mrad)", "0."));
            beamtilt_y = textToFloat(parser.getOption("--beamtilt_y", "Beamtilt in Y-direction (in mrad)", "0."));
            applyTilt = ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.;
        }

        beamtilt_xx = textToFloat(parser.getOption("--beamtilt_xx", "Anisotropic beamtilt, XX-coefficient", "1."));
        beamtilt_xy = textToFloat(parser.getOption("--beamtilt_xy", "Anisotropic beamtilt, XY-coefficient", "0."));
        beamtilt_yy = textToFloat(parser.getOption("--beamtilt_yy", "Anisotropic beamtilt, YY-coefficient", "1."));

        anisoTilt = beamtilt_xx != 1.0 || beamtilt_xy != 0.0 || beamtilt_yy != 1.0;

        nr_omp_threads = textToInteger(parser.getOption("--jomp", "Number of OMP threads", "1"));
        minMG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
        maxMG = textToInteger(parser.getOption("--max_MG", "Last micrograph index", "-1"));

        debug = parser.checkOption("--debug", "Write debugging data");

        int rco = readMoreOptions(parser, argc, argv);

        if (argc == 1)
        {
            parser.writeUsage(std::cerr);
            return 1;
        }

        if (parser.checkForErrors()) return 1;
        if (rco != 0) return rco;

        if (doesMovies)
        {
            if (bin_type_str == "box") binType = StackHelper::BoxBin;
            else if (bin_type_str == "gauss") binType = StackHelper::GaussBin;
            else if (bin_type_str == "fourier") binType = StackHelper::FourierCrop;
            else
            {
                REPORT_ERROR("Illegal binning type: " + bin_type_str + " (supported: box, gauss or fourier)");
            }
        }
    }
    catch (RelionError XE)
    {
        parser.writeUsage(std::cout);
        std::cerr << XE;
        exit(1);
    }

    bool allGood = true;

    if (!noReference)
    {
        try
        {
            maps[0].read(reconFn0);
        }
        catch (RelionError XE)
        {
            std::cout << "Unable to read map: " << reconFn0 << "\n";
            return 2;
        }

        if (!singleReference)
        {
            try
            {
                maps[1].read(reconFn1);
            }
            catch (RelionError XE)
            {
                std::cout << "Unable to read map: " << reconFn1 << "\n";
                return 3;
            }
        }

        if (maps[0].data.xdim != maps[0].data.ydim || maps[0].data.ydim != maps[0].data.zdim)
        {
            REPORT_ERROR(reconFn0 + " is not cubical.\n");
        }

        if (!singleReference)
        {
            if (maps[1].data.xdim != maps[1].data.ydim || maps[1].data.ydim != maps[1].data.zdim)
            {
                REPORT_ERROR(reconFn1 + " is not cubical.\n");
            }

            if (   maps[0].data.xdim != maps[1].data.xdim
                || maps[0].data.ydim != maps[1].data.ydim
                || maps[0].data.zdim != maps[1].data.zdim)
            {
                REPORT_ERROR(reconFn0 + " and " + reconFn1 + " are of unequal size.\n");
            }
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
                return 4;
            }

            mask.read(maskFn);

            ImageOp::multiply(mask, maps[0], maskedRef);
            maps[0] = maskedRef;

            if (!singleReference)
            {
                ImageOp::multiply(mask, maps[1], maskedRef);
                maps[1] = maskedRef;
            }
        }

        s = maps[0].data.xdim;
        sh = s/2 + 1;

        std::cout << "transforming references...\n";

        projectors[0] = Projector(s, TRILINEAR, paddingFactor, 10, 2);
        projectors[0].computeFourierTransformMap(maps[0].data, powSpec[0].data, maps[0].data.xdim);

        if (!singleReference)
        {
            projectors[1] = Projector(s, TRILINEAR, paddingFactor, 10, 2);
            projectors[1].computeFourierTransformMap(maps[1].data, powSpec[1].data, maps[1].data.xdim);
        }
    }

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

    if (!noStar)
    {
        std::cout << "reading " << starFn << "...\n";

        mdt0.read(starFn);

        if (Cs < 0.0)
        {
            mdt0.getValue(EMDL_CTF_CS, Cs, 0);
            std::cout << " + Using spherical aberration from the input STAR file: " << Cs << "\n";
        }
        else
        {
            for (int i = 0; i < mdt0.numberOfObjects(); i++)
            {
                mdt0.setValue(EMDL_CTF_CS, Cs, i);
            }
        }

        if (kV < 0.0)
        {
            mdt0.getValue(EMDL_CTF_VOLTAGE, kV, 0);
            std::cout << " + Using voltage from the input STAR file: " << kV << " kV\n";
        }
        else
        {
            for (int i = 0; i < mdt0.numberOfObjects(); i++)
            {
                mdt0.setValue(EMDL_CTF_VOLTAGE, kV, i);
            }
        }

        if (angpix <= 0.0)
        {
            RFLOAT mag, dstep;
            mdt0.getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
            mdt0.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
            angpix = 10000 * dstep / mag;

            std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << "\n";
        }

        mdts = StackHelper::splitByStack(&mdt0);

        gc = maxMG >= 0? maxMG : mdts.size()-1;
        g0 = minMG;

        std::cout << "mg range: " << g0 << ".." << gc << "\n";
    }

    RFLOAT V = kV * 1e3;
    lambda = 12.2643247 / sqrt(V * (1.0 + V * 0.978466e-6));

    obsModel = ObservationModel(angpix);

    if (applyTilt)
    {
        obsModel = ObservationModel(angpix, Cs, kV * 1e3, beamtilt_x, beamtilt_y);

        if (anisoTilt)
        {
            obsModel.setAnisoTilt(beamtilt_xx, beamtilt_xy, beamtilt_yy);
        }
    }

    int rc0 = _init();

    if (useFsc)
    {
        RefinementHelper::drawFSC(&fscMdt, freqWeight);
    }
    else
    {
        freqWeight = Image<RFLOAT>(sh,s);
        freqWeight.data.initConstant(1.0);
    }

    return rc0;
}

int RefinementProgram::run()
{
    return _run();
}
