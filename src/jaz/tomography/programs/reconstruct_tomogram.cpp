#include "reconstruct_tomogram.h"
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/util/log.h>
#include <src/args.h>
#include <src/backprojector.h>
#include <src/parallel.h>


#include <omp.h>

using namespace gravis;

double ReconstructSnrOptimisation::gradAndValue(const std::vector<double> &x, std::vector<double> &gradDest) const
{

    int size = x.size();
    if (size != ctf2.size()) REPORT_ERROR("ReconstructSnrOptimisation ERROR: input x has incorrect size!");
    if (size != gradDest.size()) REPORT_ERROR("ReconstructSnrOptimisation ERROR: input gradDest has incorrect size!");

    // Convert input vector to Matrix1D
    Matrix1D<double> xx(size);
    for (int i = 0; i < size; i++)
        VEC_ELEM(xx, i) = x[i];

    // Calculate the current value of the target function
    Matrix1D<double> Dy = D * xx;
    Matrix1D<double> diff = ctf2 * xx - snr;
    double value = diff.sum2() + (lambda / 2.) * Dy.sum2();

    // Calculate the gradient
    Matrix1D<double> grad = ctf2 * diff + lambda * D.transpose() * D * xx;

    // Convert output Matrix1D grad to vector
    for (int i = 0; i < size; i++)
        gradDest[i] = VEC_ELEM(grad, i);

    return value;
}

void TomoBackprojectProgram::readParameters(int argc, char *argv[])
{
	IOParser parser;

	parser.setCommandLine(argc, argv);

	optimisationSet.read(
		parser,
		true,             // optimisation set
		false,   false,   // particles
		true,    true,    // tomograms
		false,   false,   // trajectories
		false,   false,   // manifolds
		false,   false);  // reference

	int gen_section = parser.addSection("General options");

	tomoName = parser.getOption("--tn", "Tomogram name", "*");
	outFn = parser.getOption("--o", "Output filename (or output directory in case of reconstructing multiple tomograms)");
 	do_even_odd_tomograms = parser.checkOption("--generate_split_tomograms", "Reconstruct tomograms from even/odd movie frames or tilt image index for denoising");

    w = textToInteger(parser.getOption("--w", "Width"));
	h = textToInteger(parser.getOption("--h", "Height" ));
	d = textToInteger(parser.getOption("--d", "Thickness"));

    fourierInversion = parser.checkOption("--fourier", "Use a Fourier-inversion reconstruction algorithm");
    lambda = textToDouble(parser.getOption("--lambda", "Regularisation constant for CTF-correction of the SNRs in the Fourier-inversion algorithm", "10"));
	ctf_intact_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");
    applyWeight = !parser.checkOption("--no_weight", "Do not perform weighting in Fourier space using a Wiener filter");
	applyPreWeight = parser.checkOption("--pre_weight", "Pre-weight the 2D slices prior to backprojection");
    FourierCrop = parser.checkOption("--Fc", "Downsample the 2D images by Fourier cropping");
    do_only_unfinished = parser.checkOption("--only_do_unfinished", "Only reconstruct those tomograms that haven't finished yet");
    SNR = textToDouble(parser.getOption("--SNR", "SNR assumed by the Wiener filter", "10"));

	applyCtf = parser.checkOption("--ctf", "Perform CTF correction");
    doWiener = !parser.checkOption("--skip_wiener", "Do multiply images with CTF, but don't divide by CTF^2 in Wiener filter");

    if (!doWiener) applyCtf = true;

    zeroDC = !parser.checkOption("--keep_mean", "Do not zero the DC component of each frame");

	taperDist = textToDouble(parser.getOption("--td", "Tapering distance", "0.0"));
	taperFalloff = textToDouble(parser.getOption("--tf", "Tapering falloff", "0.0"));

    // SHWS & Aburt 19Jul2022: use zero-origins from relion-4.1 onwards....
    x0 = textToDouble(parser.getOption("--x0", "X origin", "0.0"));
    y0 = textToDouble(parser.getOption("--y0", "Y origin", "0.0"));
    z0 = textToDouble(parser.getOption("--z0", "Z origin", "0.0"));

	spacing = textToDouble(parser.getOption("--bin", "Binning", "1.0"));
    angpix_spacing = textToDouble(parser.getOption("--binned_angpix", "OR: desired pixel size after binning", "-1"));

    tiltAngleOffset = textToDouble(parser.getOption("--tiltangle_offset", "Offset applied to all tilt angles (in deg)", "0"));
    BfactorPerElectronDose = textToDouble(parser.getOption("--bfactor_per_edose", "B-factor dose-weighting per electron/A^2 dose (default is use Niko's model)", "0"));
    n_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));

    do_2dproj = parser.checkOption("--do_proj", "Use this to skip calculation of 2D projection of the tomogram along the Z-axis");
    centre_2dproj = textToInteger(parser.getOption("--centre_proj", "Central Z-slice for 2D projection (in tomogram pixels from the middle)", "0"));
    thickness_2dproj = textToInteger(parser.getOption("--thickness_proj", "Thickness of the 2D projection (in tomogram pixels)", "10"));

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	if (applyPreWeight)
	{
		applyWeight = false;
	}

	ZIO::ensureParentDir(outFn);
}
void TomoBackprojectProgram::initialise(bool verbose)
{
    if (!tomogramSet.read(optimisationSet.tomograms))
        REPORT_ERROR("ERROR: there was a problem reading the tomogram set");

    tomoIndexTodo.clear();

    if (tomoName == "*")
    {
        do_multiple = true;

        for (int idx = 0; idx < tomogramSet.size(); idx++)
        {
            if (do_only_unfinished && exists(getOutputFileName(idx,false,false)))
                continue;
            tomoIndexTodo.push_back(idx);
        }

        if (outFn[outFn.size()-1] != '/') outFn += '/';
        FileName fn_dir = outFn + "tomograms/";

    }
    else
    {
        int myidx = tomogramSet.getTomogramIndex(tomoName);
        if (myidx < 0) REPORT_ERROR("ERROR: cannot find specific tomogram name \"" + tomoName + "\" in the input star file");
        tomoIndexTodo.push_back(myidx);
        do_multiple = false;
    }

    if (verbose)
    {
        std::cout << " + Reconstructing " << tomoIndexTodo.size() << " tomograms: " << std::endl;
        for (int idx = 0; idx < tomoIndexTodo.size(); idx++)
        {
            std::cout << "  - " << tomogramSet.getTomogramName(tomoIndexTodo[idx]) << std::endl;
        }
        if (fabs(tiltAngleOffset) > 0.)
        {
            std::cout << " + Applying a tilt angle offset of " << tiltAngleOffset << " degrees" << std::endl;
        }

        if (do_2dproj)
        {
            std::cout << " + Making 2D projections " << std::endl;
            std::cout << "    - centered at " << centre_2dproj << " tomogram pixels from the centre of the tomogram" << std::endl;
            std::cout << "    - and a thickness of " << thickness_2dproj << " tomogram pixels" << std::endl;
        }
    }

}

void TomoBackprojectProgram::run(int rank, int size)
{
    long my_first_idx, my_last_idx;
    divide_equally(tomoIndexTodo.size(), size, rank , my_first_idx, my_last_idx);

    int barstep, nr_todo = my_last_idx-my_first_idx+1;
    if (rank == 0)
    {
        std::cout << " + Reconstructing ... " << std::endl;
        init_progress_bar(nr_todo);
        barstep = XMIPP_MAX(1, nr_todo / 60);
    }
    for (long idx = my_first_idx; idx <= my_last_idx; idx++)
    {
        // Abort through the pipeline_control system
        if (pipeline_control_check_abort_job())
            exit(RELION_EXIT_ABORTED);

        if (fabs(tiltAngleOffset) > 0.)
        {
            tomogramSet.applyTiltAngleOffset(idx, tiltAngleOffset);
        }

        if (fourierInversion)
        {
            reconstructOneTomogramFourier(tomoIndexTodo[idx]);
        }
        else if (do_even_odd_tomograms)
        {
            reconstructOneTomogram(tomoIndexTodo[idx],true,false); // true/false indicates to reconstruct tomogram from even frames
            reconstructOneTomogram(tomoIndexTodo[idx],false,true); // false/true indicates from odd frames
        }
        else
        {
            reconstructOneTomogram(tomoIndexTodo[idx],false,false);
        }

        if (rank == 0 && idx % barstep == 0)
            progress_bar(idx);
    }

    if (rank == 0) progress_bar(nr_todo);

}


void TomoBackprojectProgram::writeOutput(bool do_all_metadata)
{
    // If we were doing multiple tomograms, then also write the updated tomograms.star.
    if (do_multiple)
    {
        // for MPI program
        if (do_all_metadata) setMetaDataAllTomograms();

        tomogramSet.write(outFn + "tomograms.star");

        std::cout << " Written out: " << outFn << "tomograms.star" << std::endl;

    }
    else if (fabs(tiltAngleOffset) > 0.)
    {
        // Also write out the modified metadata file
        int idx=tomoIndexTodo[0];
        FileName fn_star;
        tomogramSet.globalTable.getValue(EMDL_TOMO_TILT_SERIES_STARFILE, fn_star, idx);
        FileName fn_newstar = getOutputFileWithNewUniqueDate(fn_star, outFn);
        tomogramSet.tomogramTables[idx].write(fn_newstar);

    }

}

void TomoBackprojectProgram::initialiseCtfScaleFactors(int tomoIndex, Tomogram &tomogram)
{
    // Skip initialisation of scale factor if it is already present in the tilt series star file
    if (tomogramSet.tomogramTables[tomoIndex].containsLabel(EMDL_CTF_SCALEFACTOR))
        return;

    const int fc = tomogram.frameCount;
    for (int f = 0; f < fc; f++)
    {
        RFLOAT ytilt;
        tomogramSet.tomogramTables[tomoIndex].getValueSafely(EMDL_TOMO_YTILT, ytilt, f);
        RFLOAT scale = cos(DEG2RAD(ytilt));
        tomogramSet.tomogramTables[tomoIndex].setValue(EMDL_CTF_SCALEFACTOR, scale, f);
        tomogram.centralCTFs[f].scale = scale;
    }
}

MultidimArray<RFLOAT>  TomoBackprojectProgram::getCtfCorrectedSNR(const MultidimArray<RFLOAT> &FSC, const MultidimArray<RFLOAT> &Fctf, double lambda, int verb)
{

    MultidimArray<RFLOAT> myCTF, oriSNR, corrSNR;
    myCTF.resize(FSC);
    oriSNR.resize(FSC);
    corrSNR.resize(FSC);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FSC)
    {
        // Avoid negative FSCs!
        RFLOAT myFSC = XMIPP_MAX(0.0001, DIRECT_MULTIDIM_ELEM(FSC, n) );
        // This line is because SNR is for sum of two halves
        RFLOAT FSCp = sqrt(2. * myFSC / (1. + myFSC));
        // Convert FSC to (CTF^2-affected) SNR
        if (FSCp > 0.999) DIRECT_MULTIDIM_ELEM(oriSNR, n) = 999.;
        else DIRECT_MULTIDIM_ELEM(oriSNR, n) = FSCp / (1. - FSCp);
        DIRECT_MULTIDIM_ELEM(myCTF, n) = DIRECT_A2D_ELEM(Fctf, 0, n);
    }

    // Now correct the SNR for the CTF^2 factor, but cannot divide by zeros, so do a regularised optimisation that imposes smoothness on SNRs
    ReconstructSnrOptimisation problem(myCTF, oriSNR, lambda);

    // Ignore the spatial frequencies until the first peak in the CTF: just divide SNR by CTF^2 for those in the block below this one
    // Because SNRs can be very high here, they can mess up the optimisation below, just set all values before first peak to the SNR at first_peak
    int first_peak = problem.getFirstPeak();

    std::vector<double> initial(XSIZE(oriSNR));
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(oriSNR)
    {
        if (first_peak < XSIZE(oriSNR) &&  n < first_peak)
            initial[n] = DIRECT_MULTIDIM_ELEM(oriSNR, first_peak);
        else
            initial[n] = DIRECT_MULTIDIM_ELEM(oriSNR, n);
    }
    std::vector<double> newSNR= LBFGS::optimize(initial, problem, false, 100, 1e-7, 1e-6);


    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(corrSNR)
    {
        if (n < first_peak)
        {
            if (DIRECT_MULTIDIM_ELEM(myCTF, n) < 0.00001)
            {
                REPORT_ERROR("TomoBackprojectProgram::getCtfCorrectedSNR ERROR: zero CTF before reaching first peak. Did you set Amplitude contrast to a non-zero value? Are the images perhaps over-focus?");
            }
            DIRECT_MULTIDIM_ELEM(corrSNR, n) = DIRECT_MULTIDIM_ELEM(oriSNR, n) / (DIRECT_MULTIDIM_ELEM(myCTF, n) * DIRECT_MULTIDIM_ELEM(myCTF, n));
        }
        else
        {
            DIRECT_MULTIDIM_ELEM(corrSNR, n) = newSNR[n];
        }
    }

    if (verb > 0)
    {
        std::cout << " first_peak= " << first_peak << std::endl;
        std::cout << "# shell corrSNR newSNR oriSNR corrSNR*CTF^2 oriSNR*CTF^2 CTF^2 FSC" << std::endl;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(corrSNR)
        {
            std::cout << n
                      << " " << DIRECT_MULTIDIM_ELEM(corrSNR, n)
                      << " " << newSNR[n]
                      << " " << DIRECT_MULTIDIM_ELEM(oriSNR, n)
                      << " " << DIRECT_MULTIDIM_ELEM(corrSNR, n) * DIRECT_MULTIDIM_ELEM(myCTF, n) * DIRECT_MULTIDIM_ELEM(myCTF, n)
                      << " " << DIRECT_MULTIDIM_ELEM(oriSNR, n) * DIRECT_MULTIDIM_ELEM(myCTF, n) * DIRECT_MULTIDIM_ELEM(myCTF, n)
                      << " " << DIRECT_MULTIDIM_ELEM(myCTF, n) * DIRECT_MULTIDIM_ELEM(myCTF, n)
                      << " " << DIRECT_MULTIDIM_ELEM(FSC, n)
                      << std::endl;
        }
    }

    return corrSNR;
}

void TomoBackprojectProgram::reconstructOneTomogramFourier(int tomoIndex)
{

    Tomogram tomogram1, tomogram2;

    tomogram1 = tomogramSet.loadTomogram(tomoIndex, true, true, false, w, h, d);
    tomogram2 = tomogramSet.loadTomogram(tomoIndex, true, false, true, w, h, d);

    MetaDataTable& m = tomogramSet.tomogramTables[tomoIndex];

    if (!tomogram1.hasMatrices) REPORT_ERROR("ERROR; tomograms do not have tilt series alignment parameters to calculate projectionMatrices!");

	const int fc = tomogram1.frameCount;

    BufferedImage<float> stackAct;

    // Make sure to make all images squared before using FFTs
    int square_box = (tomogram1.stack.xdim == tomogram1.stack.ydim) ? tomogram1.stack.xdim : XMIPP_MAX(tomogram1.stack.xdim, tomogram1.stack.ydim);
    // Make sqrt(2) bigger to account for the empty corners that otherwise appear with the spherical mask....
    square_box *= sqrt(2.);

    double pixelSizeAct = tomogramSet.getTiltSeriesPixelSize(tomoIndex);
    int new_box = square_box;
    if (angpix_spacing > 0.)
    {
        // downscale to new pixel size
        new_box = ROUND(square_box * (pixelSizeAct / angpix_spacing));
        new_box -= new_box % 2; //make even in case it is not already

        RFLOAT real_angpix_spacing = square_box * pixelSizeAct / new_box;
        if (fabs(real_angpix_spacing - angpix_spacing) / angpix_spacing > 0.001)
            std::cerr << "WARNING: Although the requested pixel size is " << angpix_spacing << " A/px, the actual pixel size will be " << real_angpix_spacing << " A/px due to rounding of the box size to an even number. The latter value is set to the image header. You can overwrite the header pixel size by --force_header_angpix." << std::endl;

        angpix_spacing = real_angpix_spacing;
        spacing = angpix_spacing / pixelSizeAct;
    }
    tomogramSet.globalTable.setValue(EMDL_TOMO_TOMOGRAM_BINNING, spacing, tomoIndex);

    float padding_factor = 1;
    bool skip_gridding = true;
    BackProjector BP = BackProjector(new_box, 3,
                                     "C1", TRILINEAR, padding_factor,
                                     10, 0, 1.9, 15, 2, skip_gridding);
	BP.initZeros();

    #pragma omp parallel for num_threads(n_threads)
    for (int f = 0; f < fc; f++)
    {
        //std::cerr << " f= " << f << std::endl;

        CTF ctf = tomogram1.centralCTFs[f];
        // Don't use CTF scale factors, as we will measure SNRs using the FSC!
        ctf.scale = 1.0;
        // Skip any frames that are over-focused, as first-peak calculations below will be invalid. These frames are probably bad anyway....
        if (ctf.DeltafU < 100. || ctf.DeltafV < 100.)
            continue;

        // Get frame from Jasenko's stack
        MultidimArray<RFLOAT> frame1(tomogram1.stack.ydim, tomogram1.stack.xdim);
        MultidimArray<RFLOAT> frame2(tomogram1.stack.ydim, tomogram1.stack.xdim);
        for (long int y=0; y<tomogram1.stack.ydim; y++)
            for (long int x=0; x<tomogram1.stack.xdim; x++)
            {
                DIRECT_A2D_ELEM(frame1, y, x) = tomogram1.stack(x, y, f);
                DIRECT_A2D_ELEM(frame2, y, x) = tomogram2.stack(x, y, f);
            }

        // Make square (plus factor 1.4 padding)
        frame1.setXmippOrigin();
        frame2.setXmippOrigin();
        frame1.window(FIRST_XMIPP_INDEX(square_box), FIRST_XMIPP_INDEX(square_box),
                   LAST_XMIPP_INDEX(square_box), LAST_XMIPP_INDEX(square_box));
        frame2.window(FIRST_XMIPP_INDEX(square_box), FIRST_XMIPP_INDEX(square_box),
                      LAST_XMIPP_INDEX(square_box), LAST_XMIPP_INDEX(square_box));

        // Mirror the image back out into the padding area to prevent low-resolution artifacts
        int first_x = FIRST_XMIPP_INDEX(tomogram1.stack.xdim);
        int last_x = LAST_XMIPP_INDEX(tomogram1.stack.xdim);
        int first_y = FIRST_XMIPP_INDEX(tomogram1.stack.ydim);
        int last_y = LAST_XMIPP_INDEX(tomogram1.stack.ydim);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(frame1)
        {
            int jp = j, ip = i;
            bool do_change = false;
            if (j < first_x)      {jp = 2 * first_x  - j; do_change = true;}
            else if (j > last_x)  {jp = 2 * last_x - j; do_change = true;}
            if (i < first_y)      {ip = 2 * first_y  - i; do_change = true;}
            else if (i > last_y)  {ip = 2 * last_y - i; do_change = true;}
            if (do_change)
            {
                A2D_ELEM(frame1, i, j) = A2D_ELEM(frame1, ip, jp);
                A2D_ELEM(frame2, i, j) = A2D_ELEM(frame2, ip, jp);
            }

        }

        // Downscale
        if (new_box != square_box)
        {
            resizeMap(frame1, new_box);
            resizeMap(frame2, new_box);
        }

        // Get the transformation matrix
        const Matrix2D<RFLOAT> A(3,3);
        for (int row= 0; row < 3; row++)
            for (int col = 0; col < 3; col++)
                MAT_ELEM(A, row, col) = tomogram1.projectionMatrices[f](row, col);

        RFLOAT xshift, yshift;
        m.getValueSafely(EMDL_TOMO_XSHIFT_ANGST, xshift, f);
        m.getValueSafely(EMDL_TOMO_YSHIFT_ANGST, yshift, f);

        // FT and get SNRs
        MultidimArray<Complex> FT1, FT2;
        // Not entirely sure this is necessary, but FFTW transformers have been troublesome with threads in the past...
        #pragma omp critical
        {
            FourierTransformer transformer;
            transformer.FourierTransform(frame1, FT1);
            transformer.FourierTransform(frame2, FT2);
        };
        MultidimArray<RFLOAT> FSC;
        getFSC(FT1, FT2, FSC);

        // Now that we have the FSC, sum the two halves together
        FT1 += FT2;
        // Center and shift
        CenterFFTbySign(FT1);
        shiftImageInFourierTransform(FT1, FT1, XSIZE(frame1), -xshift/angpix_spacing, -yshift/angpix_spacing);

        // Get CTF
        MultidimArray<RFLOAT>  Fctf, Ftmp;
        Fctf.resize(YSIZE(FT1), XSIZE(FT1));
        ctf.getFftwImage(Fctf, new_box, new_box, angpix_spacing, false, false, ctf_intact_first_peak, false);

        // Calculate the CTF-corrected SNR from the FSC
        MultidimArray<RFLOAT> SNR = getCtfCorrectedSNR(FSC, Fctf, lambda, 0);

        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(FT1)
        {
            long int idx = XMIPP_MIN(ROUND(sqrt(ip*ip + jp*jp)), XSIZE(SNR)-1);
            // Wiener filter = Sum(CTF*SNR*X) / (Sum(CTF^2*SNR) + 1.)
            DIRECT_A2D_ELEM(FT1, i, j)  *= DIRECT_MULTIDIM_ELEM(SNR, idx) * DIRECT_A2D_ELEM(Fctf, i, j);
            DIRECT_A2D_ELEM(Fctf, i, j) *= DIRECT_MULTIDIM_ELEM(SNR, idx) * DIRECT_A2D_ELEM(Fctf, i, j);
        }

        #pragma omp critical
        {
            BP.set2DFourierTransform(FT1, A, &Fctf);
        };

    }

    Image<RFLOAT> vol;
    vol().resize(new_box, new_box, new_box);
    vol().setXmippOrigin();
    bool do_map = true;
    MultidimArray<RFLOAT> tau2(new_box);
    tau2.initConstant(1.);
    BP.reconstruct(vol(), 0, do_map, tau2);
    if (!do_multiple) Log::print("Writing output");

    vol().window(FIRST_XMIPP_INDEX(d/spacing), FIRST_XMIPP_INDEX(h/spacing), FIRST_XMIPP_INDEX(w/spacing),
                 LAST_XMIPP_INDEX(d/spacing), LAST_XMIPP_INDEX(h/spacing), LAST_XMIPP_INDEX(w/spacing), 0.);

    vol.setSamplingRateInHeader(angpix_spacing);
    FileName fn_vol = getOutputFileName(tomoIndex, false, false, false);
    vol.write(fn_vol);

    // Also add the tomogram sizes and name to the tomogramSet
    tomogramSet.globalTable.setValue(EMDL_TOMO_SIZE_X, w, tomoIndex);
    tomogramSet.globalTable.setValue(EMDL_TOMO_SIZE_Y, h, tomoIndex);
    tomogramSet.globalTable.setValue(EMDL_TOMO_SIZE_Z, d, tomoIndex);
    tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_FILE_NAME, fn_vol, tomoIndex);

    if (do_2dproj)
    {
        const int w1 = w / spacing + 0.5;
        const int h1 = h / spacing + 0.5;
        BufferedImage<float> proj(w1, h1);
        proj.fill(0.f);
        int minz = vol().zdim/2 + centre_2dproj - thickness_2dproj/2;
        int maxz = vol().zdim/2 + centre_2dproj + thickness_2dproj/2;
        for (int z = 0; z < vol().zdim; z++)
        {
            if (z >= minz && z <= maxz)
            {
                for (int y = 0; y < vol().ydim; y++)
                    for (int x = 0; x < vol().xdim; x++)
                        proj(x, y) += DIRECT_A3D_ELEM(vol(), z, y, x);
            }
        }
        proj.write(getOutputFileName(tomoIndex, false, false, true), angpix_spacing);

        tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_PROJ2D_FILE_NAME, getOutputFileName(tomoIndex, false, false, true), tomoIndex);

    }

}

void TomoBackprojectProgram::reconstructOneTomogram(int tomoIndex, bool doEven, bool doOdd)
{

    Tomogram tomogram;

    if (doEven)
    {
    	tomogram = tomogramSet.loadTomogram(tomoIndex, true, true, false, w, h, d);
    }
    else if (doOdd)
    {
    	tomogram = tomogramSet.loadTomogram(tomoIndex, true, false, true, w, h, d);
    }
    else
    {
        tomogram = tomogramSet.loadTomogram(tomoIndex, true, false, false, w, h, d);
    }

    // Initialise CTF scale factors to cosine(tilt) if they're not present yet
    initialiseCtfScaleFactors(tomoIndex, tomogram);

	if (zeroDC) Normalization::zeroDC_stack(tomogram.stack);
	
	const int fc = tomogram.frameCount;

	BufferedImage<float> stackAct;
	std::vector<d4Matrix> projAct(fc);

	double pixelSizeAct = tomogramSet.getTiltSeriesPixelSize(tomoIndex);

    if (angpix_spacing > 0.)
    {
        spacing = angpix_spacing / pixelSizeAct;
    }

    if (!tomogram.hasMatrices) REPORT_ERROR("ERROR; tomograms do not have tilt series alignment parameters to calculate projectionMatrices!");

    const int w1 = w / spacing + 0.5;
	const int h1 = h / spacing + 0.5;
	const int t1 = d / spacing;

	if (std::abs(spacing - 1.0) < 1e-2)
	{
		projAct = tomogram.projectionMatrices;
		stackAct = tomogram.stack;
	}
	else
	{
		for (int f = 0; f < fc; f++)
		{
			projAct[f] = tomogram.projectionMatrices[f] / spacing;
			projAct[f](3,3) = 1.0;
		}
		
		if (std::abs(spacing - 1.0) > 1e-2)
		{
			if (!do_multiple) Log::print("Resampling image stack");
			
			if (FourierCrop)
			{
				stackAct = Resampling::FourierCrop_fullStack(
						tomogram.stack, spacing, n_threads, true);
			}
			else
			{
				stackAct = Resampling::downsampleFiltStack_2D_full(
						tomogram.stack, spacing, n_threads);
			}
			
			pixelSizeAct *= spacing;

            tomogramSet.globalTable.setValue(EMDL_TOMO_TOMOGRAM_BINNING, spacing, tomoIndex);
		}
		else
		{
			stackAct = tomogram.stack;
		}
	}
	
	const int w_stackAct = stackAct.xdim;
	const int h_stackAct = stackAct.ydim;
	const int wh_stackAct = w_stackAct/2 + 1;
	
	
	d3Vector orig(x0, y0, z0);
	BufferedImage<float> out(w1, h1, t1);
	out.fill(0.f);
	
	BufferedImage<float> psfStack;

	if (applyCtf)
	{
		// modulate stackAct with CTF (mind the spacing)
		
		psfStack.resize(w_stackAct, h_stackAct, fc);
		BufferedImage<fComplex> debug(wh_stackAct, h_stackAct, fc);
		
		#pragma omp parallel for num_threads(n_threads)
		for (int f = 0; f < fc; f++)
		{
			BufferedImage<float> frame = stackAct.getSliceRef(f);
			
			BufferedImage<fComplex> frameFS;
			FFT::FourierTransform(frame, frameFS, FFT::Both);
			
			CTF ctf = tomogram.centralCTFs[f];
			
			
			BufferedImage<fComplex> ctf2ImageFS(wh_stackAct, h_stackAct);
			
			const double box_size_x = pixelSizeAct * w_stackAct;
			const double box_size_y = pixelSizeAct * h_stackAct;
			
			for (int y = 0; y < h_stackAct;  y++)
			for (int x = 0; x < wh_stackAct; x++)
			{
				const double xA = x / box_size_x;
				const double yA = (y < h_stackAct/2? y : y - h_stackAct) / box_size_y;
				
                const float c = ctf.getCTF(xA, yA, false, false,
                                           true, false, 0.0, false);

				
				ctf2ImageFS(x,y) = fComplex(c*c,0);
				frameFS(x,y) *= c;
			}
			
			FFT::inverseFourierTransform(frameFS, frame, FFT::Both);
			stackAct.getSliceRef(f).copyFrom(frame);
			
			FFT::inverseFourierTransform(ctf2ImageFS, frame, FFT::Both);
			psfStack.getSliceRef(f).copyFrom(frame);
			
			debug.getSliceRef(f).copyFrom(ctf2ImageFS);
		}
	}	
	
	if (applyPreWeight)
	{
		stackAct = RealSpaceBackprojection::preWeight(stackAct, projAct, n_threads);
	}

    if (!do_multiple) Log::print("Backprojecting");
	
	RealSpaceBackprojection::backproject(
		stackAct, projAct, out, n_threads,
		orig, spacing, RealSpaceBackprojection::Linear, taperFalloff, taperDist);
	
	
	if ((applyWeight || applyCtf) && doWiener)
	{
		BufferedImage<float> psf(w1, h1, t1);
		psf.fill(0.f);
		
		if (applyCtf)
		{
			RealSpaceBackprojection::backproject(
					psfStack, projAct, psf, n_threads, 
					orig, spacing, RealSpaceBackprojection::Linear, taperFalloff, taperDist);
		}
		else
		{
			RealSpaceBackprojection::backprojectPsf(
					stackAct, projAct, psf, n_threads, orig, spacing);
		}
		
		Reconstruction::correct3D_RS(out, psf, out, 1.0 / SNR, n_threads);
	}

    if (!do_multiple) Log::print("Writing output");

    const double samplingRate = tomogramSet.getTiltSeriesPixelSize(tomoIndex) * spacing;

    if (doEven)
    	out.write(getOutputFileName(tomoIndex, true, false), samplingRate);
    else if (doOdd)
    	out.write(getOutputFileName(tomoIndex, false, true), samplingRate);
    else 
        out.write(getOutputFileName(tomoIndex, false, false), samplingRate);

    // Also add the tomogram sizes and name to the tomogramSet
    tomogramSet.globalTable.setValue(EMDL_TOMO_SIZE_X, w, tomoIndex);
    tomogramSet.globalTable.setValue(EMDL_TOMO_SIZE_Y, h, tomoIndex);
    tomogramSet.globalTable.setValue(EMDL_TOMO_SIZE_Z, d, tomoIndex);

    if (doEven)
    	tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_HALF1_FILE_NAME, getOutputFileName(tomoIndex, true, false), tomoIndex);
    else if (doOdd)
        tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_HALF2_FILE_NAME, getOutputFileName(tomoIndex, false, true), tomoIndex);
    else  
        tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_FILE_NAME, getOutputFileName(tomoIndex, false, false), tomoIndex);

    if (do_2dproj)
    {
        BufferedImage<float> proj(w1, h1);
        proj.fill(0.f);
        int minz = out.zdim/2 + centre_2dproj - thickness_2dproj/2;
        int maxz = out.zdim/2 + centre_2dproj + thickness_2dproj/2;
        for (int z = 0; z < out.zdim; z++)
        {
            if (z >= minz && z <= maxz)
            {
                for (int y = 0; y < out.ydim; y++)
                    for (int x = 0; x < out.xdim; x++)
                        proj(x, y) += out(x, y, z);
            }
        }
        if (doEven)
            proj.write(getOutputFileName(tomoIndex, true, false, true), samplingRate);
        else if (doOdd)
            proj.write(getOutputFileName(tomoIndex, false, true, true), samplingRate);
        else
            proj.write(getOutputFileName(tomoIndex, false, false, true), samplingRate);

        if (doEven)
            tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_PROJ2D_HALF1_FILE_NAME, getOutputFileName(tomoIndex, true, false, true), tomoIndex);
        else if (doOdd)
            tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_PROJ2D_HALF2_FILE_NAME, getOutputFileName(tomoIndex, false, true, true), tomoIndex);
        else
            tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_PROJ2D_FILE_NAME, getOutputFileName(tomoIndex, false, false, true), tomoIndex);

    }


}

void TomoBackprojectProgram::setMetaDataAllTomograms()
{

    for (int tomoIndex = 0; tomoIndex < tomogramSet.size(); tomoIndex++)
    {

        // SHWS 19apr2023: need to do this again for all tomograms: after completion of MPI job, leader does not know about the tomograms of the followers.
        Tomogram tomogram;
        tomogram = tomogramSet.loadTomogram(tomoIndex, false);
        if (!tomogram.hasMatrices) REPORT_ERROR("ERROR: tomograms do not have tilt series alignment parameters to calculate projectionMatrices");

        double pixelSizeAct = tomogramSet.getTiltSeriesPixelSize(tomoIndex);
        if (angpix_spacing > 0.) spacing = angpix_spacing / pixelSizeAct;
        if (std::abs(spacing - 1.0) > 1e-2)
            tomogramSet.globalTable.setValue(EMDL_TOMO_TOMOGRAM_BINNING, spacing, tomoIndex);

        // Also add the tomogram sizes and name to the tomogramSet
        tomogramSet.globalTable.setValue(EMDL_TOMO_SIZE_X, w, tomoIndex);
        tomogramSet.globalTable.setValue(EMDL_TOMO_SIZE_Y, h, tomoIndex);
        tomogramSet.globalTable.setValue(EMDL_TOMO_SIZE_Z, d, tomoIndex);

        // And the Bfactor per e/A^2 dose, if provided
        if (BfactorPerElectronDose > 0.)
            tomogramSet.globalTable.setValue(EMDL_CTF_BFACTOR_PERELECTRONDOSE, BfactorPerElectronDose, tomoIndex);

        if (do_even_odd_tomograms)
        {
            tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_HALF1_FILE_NAME,
                                             getOutputFileName(tomoIndex, true, false), tomoIndex);
            tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_HALF2_FILE_NAME,
                                             getOutputFileName(tomoIndex, false, true), tomoIndex);
        }
        else
        {
            tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_FILE_NAME,
                                             getOutputFileName(tomoIndex, false, false), tomoIndex);
        }

        if (do_2dproj)
        {
            if (do_even_odd_tomograms)
            {
                tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_PROJ2D_HALF1_FILE_NAME, getOutputFileName(tomoIndex, true, false, true), tomoIndex);
                tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_PROJ2D_HALF2_FILE_NAME, getOutputFileName(tomoIndex, false, true, true), tomoIndex);
            }
            else
            {
                tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_PROJ2D_FILE_NAME, getOutputFileName(tomoIndex, false, false, true), tomoIndex);
            }
        }
    }

}

FileName TomoBackprojectProgram::getOutputFileName(int index, bool nameEven, bool nameOdd, bool is_2dproj)
{
    // If we're reconstructing many tomograms, or the output filename is a directory: use standardized output filenames
    FileName fn_result = outFn;

    std::string dirname = (is_2dproj) ? "projections/" : "tomograms/";
    if (do_even_odd_tomograms)
    {
		if (nameEven)
		{
			fn_result += dirname + "rec_" + tomogramSet.getTomogramName(index)+"_half1.mrc";
		}
		else if (nameOdd)
		{
			fn_result += dirname + "rec_" + tomogramSet.getTomogramName(index)+"_half2.mrc";
		}
    }
	else
	{
	    fn_result += dirname + "rec_" + tomogramSet.getTomogramName(index)+".mrc";
	}

    if (!exists(fn_result.beforeLastOf("/"))) mktree(fn_result.beforeLastOf("/"));

    return fn_result;

}
