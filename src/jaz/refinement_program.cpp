#include <src/jaz/refinement_program.h>

#include <src/jaz/image_op.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/vtk_helper.h>

RefinementProgram::RefinementProgram(bool singleReference)
:   singleReference(singleReference)
{
}

int RefinementProgram::init(int argc, char *argv[])
{
    IOParser parser;

    try
    {
        parser.setCommandLine(argc, argv);

        parser.addSection("General options");

        starFn = parser.getOption("--i", "Input STAR file", "");

        if (singleReference)
        {
            reconFn0 = parser.getOption("--m", "Reference map", "");
        }
        else
        {
            reconFn0 = parser.getOption("--m0", "Reference map, half 1", "");
            reconFn1 = parser.getOption("--m1", "Reference map, half 2", "");
        }

        maskFn = parser.getOption("--mask", "Reference mask", "");
        fscFn = parser.getOption("--f", "Input STAR file with the FSC of the reference", "");
        outPath = parser.getOption("--out", "Output path", "");
        imgPath = parser.getOption("--img", "Path to images", "");

        angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (angst/pix)", "0.0"));
        paddingFactor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));

        beamtilt_x = textToFloat(parser.getOption("--beamtilt_x", "Beamtilt in the X-direction (in mrad)", "0."));
        beamtilt_y = textToFloat(parser.getOption("--beamtilt_y", "Beamtilt in the Y-direction (in mrad)", "0."));
        applyTilt = (ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.);

        nr_omp_threads = textToInteger(parser.getOption("--jomp", "Number of OMP threads", "1"));
        maxMG = textToInteger(parser.getOption("--max_MG", "First micrograph index", "-1"));
        minMG = textToInteger(parser.getOption("--min_MG", "Last micrograph index", "0"));

        debug = parser.checkOption("--debug", "Write debugging data");

        readMoreOptions(parser, argc, argv);

        if (parser.checkForErrors()) return 1;
    }
    catch (RelionError XE)
    {
        parser.writeUsage(std::cout);
        std::cerr << XE;
        exit(1);
    }

    bool allGood = true;



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

    if (useFsc)
    {
        RefinementHelper::drawFSC(&fscMdt, freqWeight);
    }
    else
    {
        freqWeight = Image<RFLOAT>(sh,s);
        freqWeight.data.initConstant(1.0);
    }

    if (!allGood)
    {
        return 1;
    }

    std::cout << "reading " << starFn << "...\n";

    mdt0.read(starFn);
    mdts = StackHelper::splitByStack(&mdt0);

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

    obsModel = ObservationModel(angpix);

    if (applyTilt)
    {
        obsModel = ObservationModel(angpix, Cs, kV * 1e3, beamtilt_x, beamtilt_y);
    }

    gc = maxMG >= 0? maxMG : mdts.size()-1;
    g0 = minMG;

    std::cout << "mg range: " << g0 << ".." << gc << "\n";

    return _init();
}

int RefinementProgram::run()
{
    return _run();
}
