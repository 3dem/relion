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
#include <src/parallel.h>

#include <omp.h>

using namespace gravis;


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

	applyWeight = !parser.checkOption("--no_weight", "Do not perform weighting in Fourier space using a Wiener filter");
	applyPreWeight = parser.checkOption("--pre_weight", "Pre-weight the 2D slices prior to backprojection");
    FourierCrop = parser.checkOption("--Fc", "Downsample the 2D images by Fourier cropping");
    do_only_unfinished = parser.checkOption("--only_do_unfinished", "Only reconstruct those tomograms that haven't finished yet");
    SNR = textToDouble(parser.getOption("--SNR", "SNR assumed by the Wiener filter", "10"));

	applyCtf = !parser.checkOption("--noctf", "Ignore the CTF");
    doWiener = !parser.checkOption("--skip_wiener", "Do multiply images with CTF, but don't divide by CTF^2 in Wiener filter");

    if (!doWiener) applyCtf = true;

    zeroDC = !parser.checkOption("--keep_mean", "Do not zero the DC component of each frame");

	taperDist = textToDouble(parser.getOption("--td", "Tapering distance", "0.0"));
	taperFalloff = textToDouble(parser.getOption("--tf", "Tapering falloff", "0.0"));

    // SHWS & Aburt 19Jul2022: use zero-origins from relion-4.1 onwards....
	//x0 = textToDouble(parser.getOption("--x0", "X origin", "1.0"));
	//y0 = textToDouble(parser.getOption("--y0", "Y origin", "1.0"));
	//z0 = textToDouble(parser.getOption("--z0", "Z origin", "1.0"));

    x0 = textToDouble(parser.getOption("--x0", "X origin", "0.0"));
    y0 = textToDouble(parser.getOption("--y0", "Y origin", "0.0"));
    z0 = textToDouble(parser.getOption("--z0", "Z origin", "0.0"));


	spacing = textToDouble(parser.getOption("--bin", "Binning", "1.0"));
    angpix_spacing = textToDouble(parser.getOption("--binned_angpix", "OR: desired pixel size after binning", "-1"));

	n_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));


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

    if (do_even_odd_tomograms)
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

        if (do_all_metadata) setMetaDataAllTomograms();

        if (do_multiple) tomogramSet.write(outFn + "tomograms.star");

        std::cout << " Written out: " << outFn << "tomograms.star" << std::endl;

    }

}


void TomoBackprojectProgram::getProjectMatrices(Tomogram &tomogram, MetaDataTable &tomogramTable)
{
/* From Alister Burt
 *
 * tilt_image_center = tilt_image_dimensions / 2
 * specimen_center = tomogram_dimensions / 2
 *
 * # Transformations, defined in order of application
 * s0 = S(-specimen_center)  # put specimen center-of-rotation at the origin
 * r0 = Rx(euler_angles['rlnTomoXTilt'])  # rotate specimen around X-axis
 * r1 = Ry(euler_angles['rlnTomoYTilt'])  # rotate specimen around Y-axis
 * r2 = Rz(euler_angles['rlnTomoZRot'])  # rotate specimen around Z-axis
 * s1 = S(specimen_shifts)  # shift projected specimen in xy (camera) plane
 * s2 = S(tilt_image_center)  # move specimen back into tilt-image coordinate system
 *
 * # compose matrices
 * transformations = s2 @ s1 @ r2 @ r1 @ r0 @ s0
 *
 */

    if (!(tomogramTable.containsLabel(EMDL_TOMO_XTILT) &&
          tomogramTable.containsLabel(EMDL_TOMO_YTILT) &&
          tomogramTable.containsLabel(EMDL_TOMO_ZROT) &&
          tomogramTable.containsLabel(EMDL_TOMO_XSHIFT_ANGST) &&
          tomogramTable.containsLabel(EMDL_TOMO_YSHIFT_ANGST)))

        REPORT_ERROR("ERROR: at least one of the input tilt series star file(s) does not contain projections matrices, NOR rlnTomoXTilt, rlnTomoYtilt, rlnTomoZRot, rlnTomoXShiftAngst or rlnTomoYShiftAng.");


    double pixelSizeAct = tomogram.optics.pixelSize;
    const int fc = tomogram.frameCount;
    for (int f = 0; f < fc; f++)
    {
        d4Matrix s0, s1, s2, r0, r1, r2, out;

        // Get specimen center
        t3Vector<double> specimen_center((double)int(w/2), (double)int(h/2), (double)int(d/2) );
        s0 = s0. translation(-specimen_center);

        // Get specimen shifts (in pixels)
        double xshift, yshift;
        tomogramTable.getValueSafely(EMDL_TOMO_XSHIFT_ANGST, xshift, f);
        tomogramTable.getValueSafely(EMDL_TOMO_YSHIFT_ANGST, yshift, f);
        t3Vector<double> specimen_shifts(xshift/pixelSizeAct, yshift/pixelSizeAct, 0.);
        s1 = s1.translation(specimen_shifts);

        // Get tilt image center
        std::vector<long int> tilt_image_center_int = tomogram.stack.getSizeVector();
        t3Vector<double> tilt_image_center((double)int(tilt_image_center_int[0]/2), (double)int(tilt_image_center_int[1]/2), 0.);
        s2 = s2.translation(tilt_image_center);

        // Get rotation matrices
        t3Vector<double> xaxis(1., 0., 0.), yaxis(0., 1., 0.), zaxis(0., 0., 1.);
        double xtilt, ytilt, zrot;
        tomogramTable.getValueSafely(EMDL_TOMO_XTILT, xtilt, f);
        tomogramTable.getValueSafely(EMDL_TOMO_YTILT, ytilt, f);
        tomogramTable.getValueSafely(EMDL_TOMO_ZROT, zrot, f);

        r0 = r0.rotation(xaxis, xtilt);
        r1 = r1.rotation(yaxis, ytilt);
        r2 = r2.rotation(zaxis, zrot);

        tomogram.projectionMatrices[f] = s2 * s1 * r2 * r1 * r0 * s0;

        // Set the four rows of the projectionMatrix back into the metaDataTable
        std::vector<EMDLabel> rows({
            EMDL_TOMO_PROJECTION_X,
            EMDL_TOMO_PROJECTION_Y,
            EMDL_TOMO_PROJECTION_Z,
            EMDL_TOMO_PROJECTION_W });

        for (int i = 0; i < 4; i++)
        {
            std::vector<double> vals(4);
            for (int j = 0; j < 4; j++)
            {
                vals[j] = tomogram.projectionMatrices[f](i,j);
            }
            tomogramTable.setValue(rows[i], vals, f);
        }

    }
}

void TomoBackprojectProgram::reconstructOneTomogram(int tomoIndex, bool doEven, bool doOdd)
{
    Tomogram tomogram;

    if (doEven)
    {
    	tomogram = tomogramSet.loadTomogram(tomoIndex, true, true, false);
    }
    else if (doOdd)
    {
    	tomogram = tomogramSet.loadTomogram(tomoIndex, true, false, true);
    }
    else
    {
        tomogram = tomogramSet.loadTomogram(tomoIndex, true);
    }
    
	if (zeroDC) Normalization::zeroDC_stack(tomogram.stack);
	
	const int fc = tomogram.frameCount;

	BufferedImage<float> stackAct;
	std::vector<d4Matrix> projAct(fc);

	double pixelSizeAct = tomogramSet.getTiltSeriesPixelSize(tomoIndex);

    if (angpix_spacing > 0.)
    {
        spacing = angpix_spacing / pixelSizeAct;
    }

    if (!tomogram.hasMatrices) getProjectMatrices(tomogram, tomogramSet.tomogramTables[tomoIndex]);

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
				
				const float c = ctf.getCTF(xA, yA);
				
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
    {
    tomogramSet.globalTable.setValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_FILE_NAME, getOutputFileName(tomoIndex, false, false), tomoIndex);
    }  
}

void TomoBackprojectProgram::setMetaDataAllTomograms()
{

    for (int idx = 0; idx < tomoIndexTodo.size(); idx++)
    {

        int tomoIndex = tomoIndexTodo[idx];

        double pixelSizeAct = tomogramSet.getTiltSeriesPixelSize(tomoIndex);
        if (angpix_spacing > 0.) spacing = angpix_spacing / pixelSizeAct;
        if (std::abs(spacing - 1.0) > 1e-2)
            tomogramSet.globalTable.setValue(EMDL_TOMO_TOMOGRAM_BINNING, spacing, tomoIndex);

        // Also add the tomogram sizes and name to the tomogramSet
        tomogramSet.globalTable.setValue(EMDL_TOMO_SIZE_X, w, tomoIndex);
        tomogramSet.globalTable.setValue(EMDL_TOMO_SIZE_Y, h, tomoIndex);
        tomogramSet.globalTable.setValue(EMDL_TOMO_SIZE_Z, d, tomoIndex);

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
    }

}

FileName TomoBackprojectProgram::getOutputFileName(int index, bool nameEven, bool nameOdd)
{
    // If we're reconstructing many tomograms, or the output filename is a directory: use standardized output filenames
    FileName fn_result = outFn;

    if (do_even_odd_tomograms)
    {
		if (nameEven)
		{
			fn_result += "tomograms/rec_" + tomogramSet.getTomogramName(index)+"_half1.mrc";
		}
		else if (nameOdd)
		{
			fn_result += "tomograms/rec_" + tomogramSet.getTomogramName(index)+"_half2.mrc";
		}
    }
	else
	{
	    fn_result += "tomograms/rec_" + tomogramSet.getTomogramName(index)+".mrc";
	}

    if (!exists(fn_result.beforeLastOf("/"))) mktree(fn_result.beforeLastOf("/"));

    return fn_result;

}
