/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
#include "src/amyloid_finder.h"

#include <array>
//#define DEBUG_BOUNDS

void AmyloidFinder::read(int argc, char **argv, int rank)
{
    parser.setCommandLine(argc, argv);


    int general_section = parser.addSection("General options");
    fn_in = parser.getOption("--i", "Input image (.mrc) or STAR file with micrographs");
    fn_out = parser.getOption("--pickname", "Rootname for coordinate STAR files", "autopick");
    fn_odir = parser.getOption("--odir", "Output directory for coordinate files (default is to store next to micrographs)", "AutoPick/");
    do_only_unfinished = parser.checkOption("--only_do_unfinished", "Only estimate CTFs for those tomograms for which there is not yet a logfile with Final values.");
    nr_threads = textToInteger(parser.getOption("--j", "Number of threads to us in parallel", "1"));

    int search_section = parser.addSection("Filament searching options ");
    psi_step = textToFloat(parser.getOption("--psi_step", "Angular sampling rate (in degrees)", "5."));
    shift_step = textToInteger(parser.getOption("--shift_step", "Step in shifts to search (in downscaled pixels)", "5"));
    search_filament_width = textToFloat(parser.getOption("--search_filament_width", "Width of searching image (in A)", "50"));
    search_filament_length = textToFloat(parser.getOption("--search_filament_length", "Length of searching image (in A)", "250"));
    do_skip_fom = parser.checkOption("--skip_fom", "Skip FOM calculation.");

    int pick_section = parser.addSection("Filament tracing options ");
    threshold = textToFloat(parser.getOption("--threshold", "Threshold in Z-scores for coordinate picking", "0.5"));
    trace_filament_width = textToFloat(parser.getOption("--trace_filament_width", "Minimum width occupied by a traced filaments (in A)", "100"));
    trace_filament_length = textToFloat(parser.getOption("--trace_filament_length", "Minimum length of traced filaments (in A)", "300"));
    psi_jump_threshold = textToFloat(parser.getOption("--psi_jump_threshold", "Maximum difference in PSI values between consecutive elements of a skeletonised branch (in degrees)", "45"));
    do_plot = parser.checkOption("--plot", "Display images with intermediate tracing results for each micrograph");
    fn_exe =  parser.getOption("--exe", "Name of python script for filament tracing", "relion_python_trace_amyloids");
    fn_other_args = parser.getOption("--other_args", "Other arguments for the python script", "");
    fn_model_path = parser.getOption("--model_path", "Name of the model to execute for filament tracing","amytracer-v2.0");
    do_carbon = parser.checkOption("--detect_carbon", "Detect carbon and ignore filaments on there.");
    fn_carbon_model_path = parser.getOption("--carbon_model_path", "Name of the model to execute for carbon detection","carbonpicker-v1.0");
    carbon_threshold = textToFloat(parser.getOption("--carbon_threshold", "Threshold for carbon detection", "0.9"));
	do_skip_tracing = parser.checkOption("--skip_tracing", "Skip tracing.");
    do_gpu = parser.checkOption("--gpu", "Use GPU acceleration when availiable");
    gpu_ids = parser.getOption("--gpu", "Device ids for each MPI-thread","default");

    int expert_section = parser.addSection("Expert options (typically no need to change)");
    signal_minres = textToFloat(parser.getOption("--signal_minres", "Minimum resolution value for signal (in A)", "4.85"));
    signal_maxres = textToFloat(parser.getOption("--signal_maxres", "Maximum resolution value for signal (in A)", "4.65"));
    nonsignal_minres = textToFloat(parser.getOption("--nonsignal_minres", "Minimum resolution value for non-signal (in A)", "4.4"));
    nonsignal_maxres = textToFloat(parser.getOption("--nonsignal_maxres", "Maximum resolution value for non-signal (in A)", "4.2"));
    down_angpix = textToFloat(parser.getOption("--down_angpix", "Pixel size for downscaled images (needs to include signal frequency!)", "2.1"));
    angpix = textToFloat(parser.getOption("--force_angpix", "Force this pixel size, regardless of what is in the image header", "-1"));
    verb =textToInteger(parser.getOption("--verb", "Verbosity", "1"));

    // Check for errors in the command-line option
    if (parser.checkForErrors())
        REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}

void AmyloidFinder::usage()
{
    parser.writeUsage(std::cout);
}

void AmyloidFinder::initialise(bool is_leader)
{
    // Make sure fn_odir ends with a slash
    if (fn_odir[fn_odir.length()-1] != '/')
        fn_odir += "/";

    fn_ori_micrographs.clear();
    fn_ori_micrographs_fom.clear();
    fn_ori_micrographs_psi.clear();

    todo_micrographs_fom.clear();
    todo_micrographs_tracing.clear();
    if (fn_in.isStarFile())
    {
        MetaDataTable MDin;
        MDin.read(fn_in, "micrographs");

        if (do_skip_fom && !(MDin.containsLabel(EMDL_MICROGRAPH_AUTOPICK_FOM) && MDin.containsLabel(EMDL_MICROGRAPH_AUTOPICK_PSI)))
        {
            REPORT_ERROR("ERROR: You are skipping FOM calculation, but the input STAR file does not contain the FOM/PSI images.");
        }
        FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
        {
            FileName fn_mic;
            MDin.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
            fn_ori_micrographs.push_back(fn_mic);
            if (MDin.containsLabel(EMDL_MICROGRAPH_AUTOPICK_FOM) && MDin.containsLabel(EMDL_MICROGRAPH_AUTOPICK_PSI))
            {
                MDin.getValue(EMDL_MICROGRAPH_AUTOPICK_FOM, fn_mic);
                fn_ori_micrographs_fom.push_back(fn_mic);
                MDin.getValue(EMDL_MICROGRAPH_AUTOPICK_PSI, fn_mic);
                fn_ori_micrographs_psi.push_back(fn_mic);
            }
        }
    }
    else
    {
        // Read a single micrograph
        fn_ori_micrographs.push_back(fn_in);
    }

    if (!do_only_unfinished)
    {
        if (!do_skip_fom) todo_micrographs_fom = fn_ori_micrographs;
        if (!do_skip_tracing)
        {
            todo_micrographs_tracing = fn_ori_micrographs;
            idx_todo_micrographs_tracing.resize(todo_micrographs_tracing.size());
            std::iota(idx_todo_micrographs_tracing.begin(), idx_todo_micrographs_tracing.end(), 0);
        }
    }
    else
    {
        // If we're continuing an old run, see which micrographs have not been finished yet...
        // A. for FOM/PSI calculation
        if (!do_skip_fom)
        {
            if (verb > 0)
            {
                std::cout << " + Skipping those micrographs for which FOM and PSI images already exist" << std::endl;
            }
            todo_micrographs_fom.clear();
            for (long int imic = 0; imic < fn_ori_micrographs.size(); imic++)
            {
                FileName fn_fom = getOutputRootName(fn_ori_micrographs[imic]) + "_" + fn_out + "_fom.mrc";
                FileName fn_psi = getOutputRootName(fn_ori_micrographs[imic]) + "_" + fn_out + "_psi.mrc";
                if (!exists(fn_fom) || !exists(fn_psi))
                    todo_micrographs_fom.push_back(fn_ori_micrographs[imic]);
            }
        }
        // B. For tracing
        if (!do_skip_tracing)
        {
            if (do_skip_fom && (fn_ori_micrographs_fom.size() == 0 || fn_ori_micrographs_psi.size() == 0))
            {
                REPORT_ERROR("ERROR: you cannot skip FOM calculation without providing autopick FOM and PSI images in the input STAR file!");
            }

            if (verb > 0)
            {
                std::cout << " + Skipping those micrographs for which coordinate file already exists" << std::endl;
            }
            for (long int imic = 0; imic < fn_ori_micrographs.size(); imic++)
            {
                FileName fn_tmp = getOutputRootName(fn_ori_micrographs[imic]) + "_" + fn_out + ".star";
                if (!exists(fn_tmp))
                {
                    todo_micrographs_tracing.push_back(fn_ori_micrographs[imic]);
                    idx_todo_micrographs_tracing.push_back(imic);
                }
            }

        }
    }

    if (verb > 0)
    {
        std::cout << " + Calculating FOM images for " << todo_micrographs_fom.size() << " micrographs... " << std::endl;
        std::cout << " + Tracing filaments for " << todo_micrographs_tracing.size() << " micrographs... " << std::endl;
    }

    // Read in header of first image
    Image<RFLOAT> Iin;
    Iin.read(fn_ori_micrographs[0], false);
    ori_xsize = XSIZE(Iin());
    ori_ysize = YSIZE(Iin());
    Iin().setXmippOrigin();
    if (angpix < 0.)
    {
        angpix = Iin.samplingRateX();
        if (verb > 0) std::cout << " - Using pixel size from the header of : " << fn_in << " = " << angpix << std::endl;
    }
    if (nonsignal_maxres < 2*down_angpix) REPORT_ERROR("ERROR: the down_angpix is not enough to support the maximum resolution of the signal!");
    if (angpix > down_angpix) REPORT_ERROR("ERROR: this program requires input images with a pixel size of at least down_angpix (" + floatToString(down_angpix) + ")!");

    if (todo_micrographs_fom.size() > 0)
    {


        // Width and length in the downscaled pixels
        iwidthmax = ROUND(search_filament_width / down_angpix );
        ilengthmax = CEIL(search_filament_length / down_angpix );

        down_xsize = FLOOR( (ori_xsize * angpix) / down_angpix );
        down_ysize = FLOOR( (ori_ysize * angpix) / down_angpix );
        if (ilengthmax %2 != 0) ilengthmax++;
        nr_psi = ROUND(180./psi_step);
        psi_step = 180./nr_psi;

        // Calculate Fourier shells for amyloid signal
        imin_signal = FLOOR(ilengthmax*down_angpix/signal_minres);
        imax_signal = CEIL(ilengthmax*down_angpix/signal_maxres);
        imin_nonsignal = FLOOR(ilengthmax*down_angpix/nonsignal_minres);
        imax_nonsignal = CEIL(ilengthmax*down_angpix/nonsignal_maxres);

        // Box size, orginal and cropped: set size of rectangular image to largest dimension
        large_box = sqrt(2.)*XMIPP_MAX(ori_xsize, ori_ysize);
        large_box += ROUND(XMIPP_MAX(search_filament_width, search_filament_length) / angpix);
        if (large_box%2 != 0) large_box++;
        // Also calculate size of cropped box:
        crop_box = large_box * angpix/down_angpix;
        if (crop_box%2 != 0) crop_box++;

        // Output some information to the user
        if (verb > 0)
        {
            std::cout << " + Number of 1D rows for filament width (in downscaled pixels): " << iwidthmax << std::endl;
            std::cout << " + Length of 1D rows for filament (in downscaled pixels): " << ilengthmax << std::endl;
            std::cout << " + Number of in-plane rotations to sample: " << nr_psi << " with step of " << psi_step << " degrees" << std::endl;
            std::cout << " + Original size of the input micrographs: " << ori_xsize << " x " << ori_ysize << " pixels" << std::endl;
            std::cout << " + Size of image to sample (in downscaled pixels): " <<  down_xsize << " x " << down_ysize << std::endl;
            std::cout << " + Fourier shells for the amyloid signal (in downscaled pixels): " << imin_signal  << " - " << imax_signal << std::endl;
            std::cout << " + Fourier shells for the non-signal control (in downscaled pixels): " << imin_nonsignal  << " - " << imax_nonsignal << std::endl;
            std::cout << "  ========================== " << std::endl;
        }
    }

}


#if defined _CUDA_ENABLED
void AmyloidFinder::deviceInitialise()
{
	int devCount;
	accGPUGetDeviceCount(&devCount);

	std::vector < std::vector < std::string > > allThreadIDs;
	untangleDeviceIDs(gpu_ids, allThreadIDs);

	// Sequential initialisation of GPUs on all ranks
	if (!std::isdigit(*gpu_ids.begin()))
		device_id = 0;
	else
		device_id = textToInteger((allThreadIDs[0][0]).c_str());

	if (verb>0)
	{
		std::cout << " + Using GPU device " << device_id << std::endl;
	}
}
#endif

FileName AmyloidFinder::getOutputRootName(FileName fn_mic)
{
	FileName fn_pre, fn_jobnr, fn_post;
	decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
	return fn_odir + fn_post.withoutExtension();
}

RFLOAT AmyloidFinder::getPsiAngle(int ipsi)
{
    return 2.3 + ipsi * psi_step;
}

MultidimArray<RFLOAT> AmyloidFinder::growNonSignalMask(MultidimArray<RFLOAT> &inmask, int extend_size)
{

    MultidimArray<RFLOAT> Mresult = inmask;
    RFLOAT extend_ini_mask2 = extend_size * extend_size;

#pragma omp parallel for num_threads(nr_threads)
    for (long int i=STARTINGY(inmask)+extend_size; i<=FINISHINGY(inmask)-extend_size; i++)
    {
        for (long int j=STARTINGX(inmask)+extend_size; j<=FINISHINGX(inmask)-extend_size; j++)
        {
            // only extend from 1 values
            if (A2D_ELEM(inmask, i, j) > 0.99)
            {
                for (long int ip = i - extend_size; ip <= i + extend_size; ip++)
                {
                    for (long int jp = j - extend_size; jp <= j + extend_size; jp++)
                    {
                        // only check distance if neighbouring pixel is zero
                        if (A2D_ELEM(inmask, ip, jp) < 0.01)
                        {
                            RFLOAT r2 = (RFLOAT)( (ip-i)*(ip-i)+ (jp-j)*(jp-j) );
                            // Set original voxel to 1 if a neghouring with Im()=1 is within distance extend_ini_mask
                            if (r2 < extend_ini_mask2)
                            {
                                A2D_ELEM(Mresult, ip, jp) = 1.;
                            }
                        }
                    }
                }
            }
        }
    }

    return Mresult;

}

namespace
{
struct AmyloidPsiScore
{
    MultidimArray<RFLOAT> signal;
    MultidimArray<RFLOAT> nonsignal;
};

void calculateAmyloidPsiScore(
        AmyloidFinder &finder,
        const MultidimArray<RFLOAT> &Mbig,
        int ipsi,
        std::vector<FourierTransformer> &rotation_transformers,
        std::vector<FourierTransformer> &line_transformers,
        AmyloidPsiScore &out)
{
    const int half_nr_psi = finder.nr_psi / 2;
    const int source_ipsi = ipsi % half_nr_psi;
    const RFLOAT source_psi = finder.getPsiAngle(source_ipsi);

    MultidimArray<RFLOAT> Mrot;
    Mrot.setXmippOrigin();
    Mrot.initZeros(finder.large_box, finder.large_box);
    rotate(Mbig, Mrot, source_psi, 'Z', true);

    MultidimArray<Complex > FT, FT2;
    rotation_transformers[0].FourierTransform(Mrot, FT, false);
    windowFourierTransform(FT, FT2, finder.crop_box);
    rotation_transformers[0].clear();
    Mrot.reshape(finder.crop_box, finder.crop_box);
    rotation_transformers[0].inverseFourierTransform(FT2, Mrot);
    Mrot.setXmippOrigin();

    if (ipsi >= half_nr_psi)
    {
        MultidimArray<RFLOAT> Mrot90;
        Mrot90.initZeros(finder.crop_box, finder.crop_box);
        Mrot90.setXmippOrigin();
        for (long int i = STARTINGY(Mrot90) + 1; i <= FINISHINGY(Mrot90) - 1; i++)
        {
            for (long int j = STARTINGX(Mrot90) + 1; j <= FINISHINGX(Mrot90) - 1; j++)
            {
                A2D_ELEM(Mrot90, i, j) = A2D_ELEM(Mrot, -j, i);
            }
        }
        Mrot = Mrot90;
    }

    MultidimArray<RFLOAT> scores_perline, nonscores_perline;
    scores_perline.initZeros(finder.crop_box, finder.crop_box);
    scores_perline.setXmippOrigin();
    nonscores_perline.initZeros(finder.crop_box, finder.crop_box);
    nonscores_perline.setXmippOrigin();

    out.signal.initZeros(finder.crop_box / finder.shift_step, finder.crop_box / finder.shift_step);
    out.signal.setXmippOrigin();
    out.nonsignal.initZeros(finder.crop_box / finder.shift_step, finder.crop_box / finder.shift_step);
    out.nonsignal.setXmippOrigin();

    const int my_skip_side_length = finder.ilengthmax / 2;
#pragma omp parallel for num_threads(finder.nr_threads)
    for (int ypos = my_skip_side_length; ypos < YSIZE(Mrot) - my_skip_side_length; ypos += 1)
    {
        const int cen_ypos = ypos - YSIZE(Mrot) / 2;
        const int tid = omp_get_thread_num();
        MultidimArray<RFLOAT> oneline(finder.ilengthmax);
        MultidimArray<Complex> FTline(finder.ilengthmax / 2 + 1);

        for (int xpos = my_skip_side_length; xpos < XSIZE(Mrot) - my_skip_side_length; xpos += 1)
        {
            const int cen_xpos = xpos - XSIZE(Mrot) / 2;

            for (int iline = 0; iline < finder.ilengthmax; iline++)
                DIRECT_A1D_ELEM(oneline, iline) = A2D_ELEM(Mrot, cen_ypos, cen_xpos + iline - finder.ilengthmax / 2);

            line_transformers[tid].FourierTransform(oneline, FTline, false);

            for (int isig = finder.imin_signal; isig <= finder.imax_signal; isig++)
                A2D_ELEM(scores_perline, cen_ypos, cen_xpos) += norm(DIRECT_A1D_ELEM(FTline, isig));

            for (int isig = finder.imin_nonsignal; isig <= finder.imax_nonsignal; isig++)
                A2D_ELEM(nonscores_perline, cen_ypos, cen_xpos) += norm(DIRECT_A1D_ELEM(FTline, isig));
        }
    }

    const int my_skip_side_width = finder.iwidthmax / 2;
#pragma omp parallel for num_threads(finder.nr_threads)
    for (int ypos = 0; ypos < YSIZE(Mrot) / 2 - my_skip_side_width; ypos += finder.shift_step)
    {
        for (int ipassy = 0; ipassy < 2; ipassy++)
        {
            const int cen_ypos = (ipassy == 0) ? ypos : -ypos;
            if (ypos == 0 && ipassy == 1) continue;

            for (int xpos = 0; xpos < XSIZE(Mrot) / 2 - my_skip_side_width; xpos += finder.shift_step)
            {
                for (int ipass = 0; ipass < 2; ipass++)
                {
                    const int cen_xpos = (ipass == 0) ? xpos : -xpos;
                    if (xpos == 0 && ipass == 1) continue;

                    for (int iwidth = 0; iwidth < finder.iwidthmax; iwidth++)
                    {
                        A2D_ELEM(out.signal, cen_ypos / finder.shift_step, cen_xpos / finder.shift_step) +=
                                A2D_ELEM(scores_perline, cen_ypos + iwidth - finder.iwidthmax / 2, cen_xpos);
                        A2D_ELEM(out.nonsignal, cen_ypos / finder.shift_step, cen_xpos / finder.shift_step) +=
                                A2D_ELEM(nonscores_perline, cen_ypos + iwidth - finder.iwidthmax / 2, cen_xpos);
                    }
                }
            }
        }
    }

    const RFLOAT psi = finder.getPsiAngle(ipsi);
    selfRotate(out.signal, -psi);
    selfRotate(out.nonsignal, -psi);
}
}

void AmyloidFinder::getScoreForOneMicrograph(MultidimArray<RFLOAT> &image, MultidimArray<RFLOAT> &Mscore,
                                             MultidimArray<RFLOAT> &Mangle, RFLOAT &skew, RFLOAT &kurt, bool myverb)
{

    MultidimArray<RFLOAT> Mbig(large_box, large_box);
    Mbig.setXmippOrigin();
    for (long int i=STARTINGY(Mbig); i<=FINISHINGY(Mbig); i++)
    {
        long int ip = i;
        //if (i < STARTINGY(image)) ip += YSIZE(image);
        //else if (i > FINISHINGY(image)) ip -= YSIZE(image);
        if (i < STARTINGY(image)) ip = 2*STARTINGY(image) - i;
        else if (i > FINISHINGY(image)) ip = 2*FINISHINGY(image) - i;

        for (long int j=STARTINGX(Mbig); j<=FINISHINGX(Mbig); j++)
        {
            long int jp = j;
            if (j < STARTINGX(image)) jp = 2*STARTINGX(image) - j;
            else if (j > FINISHINGX(image)) jp = 2*FINISHINGX(image) - j;

            A2D_ELEM(Mbig, i, j) = A2D_ELEM(image, ip, jp);
        }
    }

    std::vector<FourierTransformer> rotation_transformers(1);
    std::vector<FourierTransformer> line_transformers(nr_threads);
    MultidimArray<RFLOAT> oneline_tmp(ilengthmax);
    for (int i = 0; i < nr_threads; i++)
        line_transformers[i].setReal(oneline_tmp);

    Mangle.initZeros(down_ysize/shift_step, down_xsize/shift_step);
    Mangle.setXmippOrigin();
    Mscore.initZeros(Mangle);
    MultidimArray<RFLOAT> Msum, Mnonsum, Mnonscore, Mneighbour, Mneighbour2;
    Msum.initZeros(Mangle);
    Mnonsum.initZeros(Mangle);
    Mnonscore.initZeros(Mangle);
    Mneighbour.initZeros(Mangle);
    Mneighbour2.initZeros(Mangle);

    if (myverb)
    {
        std::cout << " - Searching over all orientations and coordinates ..." << std::endl;
        init_progress_bar(nr_psi);
    }

    std::array<AmyloidPsiScore, 3> psi_scores;
    int prev_slot = 0;
    int curr_slot = 1;
    int next_slot = 2;

    calculateAmyloidPsiScore(*this, Mbig, nr_psi - 1, rotation_transformers, line_transformers, psi_scores[prev_slot]);
    calculateAmyloidPsiScore(*this, Mbig, 0, rotation_transformers, line_transformers, psi_scores[curr_slot]);

    const int xsize = XSIZE(Mscore);
    const int ysize = YSIZE(Mscore);
    for (int ipsi = 0; ipsi < nr_psi; ipsi++)
    {
        calculateAmyloidPsiScore(*this, Mbig, (ipsi + 1) % nr_psi, rotation_transformers, line_transformers, psi_scores[next_slot]);

        const RFLOAT mypsi = getPsiAngle(ipsi);
        const AmyloidPsiScore &previous = psi_scores[prev_slot];
        const AmyloidPsiScore &current = psi_scores[curr_slot];
        const AmyloidPsiScore &next = psi_scores[next_slot];

        for (int ypos = 0; ypos < ysize; ypos ++)
        {
            int cen_ypos = ypos - ysize/2;
            for (int xpos = 0; xpos < xsize; xpos ++)
            {
                int cen_xpos = xpos - xsize/2;

                RFLOAT myscore = A2D_ELEM(current.signal, cen_ypos, cen_xpos);
                RFLOAT mynonscore = A2D_ELEM(current.nonsignal, cen_ypos, cen_xpos);
                A2D_ELEM(Msum, cen_ypos, cen_xpos) += myscore;
                A2D_ELEM(Mnonsum, cen_ypos, cen_xpos) += mynonscore;

                if (myscore > A2D_ELEM(Mscore, cen_ypos, cen_xpos))
                {
                    A2D_ELEM(Mscore, cen_ypos, cen_xpos) = myscore;
                    A2D_ELEM(Mangle, cen_ypos, cen_xpos) = mypsi;
                    A2D_ELEM(Mneighbour, cen_ypos, cen_xpos) =
                            A2D_ELEM(previous.signal, cen_ypos, cen_xpos) +
                            A2D_ELEM(next.signal, cen_ypos, cen_xpos);
                }

                if (mynonscore > A2D_ELEM(Mnonscore, cen_ypos, cen_xpos))
                {
                    A2D_ELEM(Mnonscore, cen_ypos, cen_xpos) = mynonscore;
                    A2D_ELEM(Mneighbour2, cen_ypos, cen_xpos) =
                            A2D_ELEM(previous.nonsignal, cen_ypos, cen_xpos) +
                            A2D_ELEM(next.nonsignal, cen_ypos, cen_xpos);
                }
            }
        }

        if (myverb) progress_bar(ipsi);

        const int old_prev_slot = prev_slot;
        prev_slot = curr_slot;
        curr_slot = next_slot;
        next_slot = old_prev_slot;
    }
    if (myverb) progress_bar(nr_psi);

//#define DEBUG_FOM
#ifdef DEBUG_FOM
    Image<RFLOAT> It, It2;
    It()= Mscore;
    FileName fnt="Mscore.spi";
    It.write(fnt);
    std::cerr <<" written: "<<fnt << std::endl;
    It2()= Mnonscore;
    fnt="Mnonscore.spi";
    It2.write(fnt);
    std::cerr <<" written: "<<fnt << std::endl;
#endif

    // Now need to subtract the max Mscore for the best psi (plus sum of its two neighbouring ipsi), as the mean and stddev should be calculated for the non-signal!
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Msum)
    {
        DIRECT_MULTIDIM_ELEM(Msum, n) -= DIRECT_MULTIDIM_ELEM(Mscore, n) + DIRECT_MULTIDIM_ELEM(Mneighbour, n);
        DIRECT_MULTIDIM_ELEM(Msum, n) /= (RFLOAT)(nr_psi-3);
        // calculate as normalised score:
        // (max_psi - adjusted_mean_psi) / adjusted_mean_psi
        // where adjusted_mean is the average of the score over all psi-values, except the maximum
        RFLOAT Zscore_signal = 0.;
        if (DIRECT_MULTIDIM_ELEM(Msum, n) > 0.)
            Zscore_signal = (DIRECT_MULTIDIM_ELEM(Mscore, n) - DIRECT_MULTIDIM_ELEM(Msum, n)) / DIRECT_MULTIDIM_ELEM(Msum, n);
        DIRECT_MULTIDIM_ELEM(Mscore, n) = Zscore_signal;

#ifdef DEBUG_FOM
       DIRECT_MULTIDIM_ELEM(It(), n) = Zscore_signal;
#endif

        // Also for non-signal
        DIRECT_MULTIDIM_ELEM(Mnonsum, n) -= DIRECT_MULTIDIM_ELEM(Mnonscore, n) + DIRECT_MULTIDIM_ELEM(Mneighbour2, n);
        DIRECT_MULTIDIM_ELEM(Mnonsum, n) /= (RFLOAT)(nr_psi-3);
        RFLOAT Zscore_nonsignal = 0.;
        if (DIRECT_MULTIDIM_ELEM(Mnonsum, n) > 0.)
        {
            Zscore_nonsignal = (DIRECT_MULTIDIM_ELEM(Mnonscore, n) - DIRECT_MULTIDIM_ELEM(Mnonsum, n)) / DIRECT_MULTIDIM_ELEM(Mnonsum, n);

#ifdef DEBUG_FOM
            DIRECT_MULTIDIM_ELEM(It2(), n) = Zscore_nonsignal;
#endif
            // binarize to generate a non-signal mask
            Zscore_nonsignal = (Zscore_nonsignal < 0.7) ? 0 : 1;
        }
        DIRECT_MULTIDIM_ELEM(Mnonscore, n) = Zscore_nonsignal;

    }

    // Grow the nonsignal mask a bit, as ice crystals give artefacts near their borders
    Mnonscore = growNonSignalMask(Mnonscore, iwidthmax);

#ifdef DEBUG_FOM
    fnt="Zscore_signal.spi";
    It.write(fnt);
    fnt="Zscore_nonsignal.spi";
    It2.write(fnt);
    It()=Mnonscore;
    It.write("grownmask.spi");
#endif

    // Apply inverse non-signal mask to the Mscore to calculate final FOM image
    RFLOAT sum=0., sum2=0.;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mscore)
    {
        DIRECT_MULTIDIM_ELEM(Mscore, n) *= (1. - DIRECT_MULTIDIM_ELEM(Mnonscore, n));
        // Also calculate mean and stddev of final combined score over the whole micrograph, to later calculate skewness and kurtosis for signal detection
        sum  += DIRECT_MULTIDIM_ELEM(Mscore, n);
        sum2 += DIRECT_MULTIDIM_ELEM(Mscore, n) * DIRECT_MULTIDIM_ELEM(Mscore, n);

    }

#ifdef DEBUG_FOM
    It()=Mscore;
    fnt="fom.spi";
    It.write(fnt);
#endif

    // Output skewness and kurtosis of Mscore distribution to detect which micrographs have filaments
    RFLOAT n = NZYXSIZE(Msum);
    sum /= n;
    sum2 /= n;
    sum2 = sqrt(sum2-sum*sum);
    skew = 0.;
    kurt = 0.;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Msum)
    {
        RFLOAT aux = (DIRECT_MULTIDIM_ELEM(Mscore, n) - sum)/sum2;
        skew += aux*aux*aux;
        kurt += aux*aux*aux*aux;
    }
    kurt *= n*(n+1)/((n-1)*(n-2)*(n-3));
    skew *= n/((n-1)*(n-2));

    if (myverb) progress_bar(nr_psi);

}


void AmyloidFinder::calculateFOMOneMicrograph(FileName fn_mic, bool myverb)
{

    FileName fn_root = getOutputRootName(fn_mic);
    FileName fn_fom = fn_root + "_" + fn_out + "_fom.mrc";
    FileName fn_psi = fn_root + "_" + fn_out + "_psi.mrc";
    FileName fn_skew = fn_root + "_" + fn_out + "_skew.txt";
    MultidimArray<RFLOAT> Mscore, Mangle;

    if (!exists(fn_fom) || !exists(fn_psi))
    {

        Image<RFLOAT> Iin;
        Iin.read(fn_mic);
        Iin().setXmippOrigin();
        if (XSIZE(Iin()) != ori_xsize || YSIZE(Iin()) != ori_ysize || fabs(angpix - Iin.samplingRateX()) > 0.001)
            REPORT_ERROR("ERROR: incorrect size or pixel size for image " + fn_mic);

        RFLOAT skew, kurt;
        getScoreForOneMicrograph(Iin(), Mscore, Mangle, skew, kurt, myverb);

        Image<RFLOAT> Ipsi, Izscore;
        Ipsi.setSamplingRateInHeader(down_angpix*shift_step);
        Izscore.setSamplingRateInHeader(down_angpix*shift_step);
        Ipsi()=Mangle;
        Izscore()=Mscore;
        // Set the skewness and kurtosis in the header of the FOM image
        Izscore.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, skew);
        Izscore.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, kurt);
        Ipsi.write(fn_psi);
        Izscore.write(fn_fom);
    }

    if (myverb) std::cout << "done!" << std::endl;



}


void AmyloidFinder::runFOMBatch(long int my_first, long int my_last)
{

    long int my_nr = my_last - my_first + 1;
    if (my_nr <= 0) return;

    int barstep;
    if (verb > 0)
    {
        std::cout << " Calculating FOMs ..." << std::endl;
        init_progress_bar(my_nr);
        barstep = XMIPP_MAX(1, my_nr / 60);
    }

    FileName fn_olddir="";
    for (long int imic = my_first; imic <= my_last; imic++)
    {

        // Abort through the pipeline_control system
        if (pipeline_control_check_abort_job())
            exit(RELION_EXIT_ABORTED);

        // Check new-style outputdirectory exists and make it if not!
        FileName fn_oroot = getOutputRootName(todo_micrographs_fom[imic]);
        FileName fn_dir = fn_oroot.beforeLastOf("/");
        if (fn_dir != fn_olddir)
        {
            // Make a Particles directory
            mktree(fn_dir);
            fn_olddir = fn_dir;
        }
#ifdef TIMING
        timer.tic(TIMING_A5);
#endif
        calculateFOMOneMicrograph(todo_micrographs_fom[imic], todo_micrographs_fom.size() == 1);
        if (verb > 0 && (imic-my_first+1)%barstep == 0) progress_bar(imic - my_first + 1);
#ifdef TIMING
        timer.toc(TIMING_A5);
#endif
    }

    if (verb > 0) progress_bar(my_nr);

}

void AmyloidFinder::runTracingBatch(long int my_first, long int my_last, int my_rank)
{
    long int my_nr = my_last - my_first + 1;
    if (my_nr <= 0) return;

    if (verb > 0) std::cout << " - Tracing filaments ..." << std::endl;

    // TODO!!! Make a temp STAR file with mic, fom and psi names to pass to Jenny's program!!!!
    // takes from my_first to my_last!
    FileName fn_tracing_star = fn_odir + "input_trace_rank" + integerToString(my_rank) + ".star";
    MetaDataTable MDtrace;
    FileName fn_fom, fn_psi;
    for (long int imic = my_first; imic <= my_last; imic++)
    {
        MDtrace.addObject();
        FileName fn_root = getOutputRootName(todo_micrographs_tracing[imic]);
        if (do_skip_fom)
        {
            long int imic_ori = idx_todo_micrographs_tracing[imic];
            fn_fom = fn_ori_micrographs_fom[imic_ori];
            fn_psi = fn_ori_micrographs_psi[imic_ori];
        }
        else
        {
            fn_fom = fn_root + "_" + fn_out + "_fom.mrc";
            fn_psi = fn_root + "_" + fn_out + "_psi.mrc";
        }
        // remove leading job number from fn_mic filename
        FileName fn_pick = fn_root + "_" + fn_out + ".star";
        MDtrace.setValue(EMDL_MICROGRAPH_NAME, todo_micrographs_tracing[imic]);
        MDtrace.setValue(EMDL_MICROGRAPH_AUTOPICK_FOM, fn_fom);
        MDtrace.setValue(EMDL_MICROGRAPH_AUTOPICK_PSI, fn_psi);
        MDtrace.setValue(EMDL_MICROGRAPH_COORDINATES, fn_pick);
    }
    MDtrace.write(fn_tracing_star);

    Image<RFLOAT> It;
    It.read(fn_fom, false);
    down_angpix = It.samplingRateX();

    // hardcoded python script for now...
    FileName command = fn_exe;

    command += " -i " + fn_tracing_star;
    command += " -m " + fn_model_path;
    if (do_gpu)
    {
        command += " -d cuda:" + integerToString(device_id);
    }
    else
    {
        command += " -d cpu";
        command += " -j " + integerToString(nr_threads);
    }
    command += " -t " + floatToString(threshold);
    command += " -r " + floatToString(trace_filament_width/2);
    command += " -l " + floatToString(trace_filament_length);
    command += " -p " + floatToString(psi_jump_threshold);
    command += " -s " + floatToString(down_angpix/angpix);
    command += " -a " + pipeline_control_outputname+RELION_JOB_ABORT_NOW;
    command += " -v " + integerToString(verb);
    if (do_carbon)
    {
        command += " -c ";
        command += " -cm " + fn_carbon_model_path;
        command +=  " --carbon_threshold " + floatToString(carbon_threshold);
    }
    if (do_plot)
        command += " --plot ";

    command += " " + fn_other_args;

    std::cerr << command << std::endl;
    int res = system(command.c_str());

    if (pipeline_control_check_abort_job())
        exit(RELION_EXIT_ABORTED);
    else if (res != 0) exit(RELION_EXIT_FAILURE);

}


void AmyloidFinder::run()
{

    runFOMBatch(0, todo_micrographs_fom.size() - 1);
    runTracingBatch(0, todo_micrographs_tracing.size() - 1);

}


void AmyloidFinder::finalise()
{

    long int barstep = XMIPP_MAX(1, fn_ori_micrographs.size() / 60);
	if (verb > 0)
	{
		std::cout << " Generating  output list of coordinate files ... " << std::endl;
		init_progress_bar(fn_ori_micrographs.size());
	}

    MetaDataTable MDin;
    ObservationModel obsModel;
    ObservationModel::loadSafely(fn_in, obsModel, MDin, "micrographs", verb);

    MetaDataTable MDcoords;
    MDcoords.setName("coordinate_files");
	long total_nr_picked = 0;
	int nr_coord_files = 0;
	for (long int imic = 0; imic < fn_ori_micrographs.size(); imic++)
	{

        FileName fn_root = getOutputRootName(fn_ori_micrographs[imic]);
        FileName fn_fom = (do_skip_fom) ? fn_ori_micrographs_fom[imic] : fn_root + "_" + fn_out + "_fom.mrc";

        if (!do_skip_fom)
		{
            FileName fn_psi = fn_root + "_" + fn_out + "_psi.mrc";

            Image<RFLOAT> Ifom;
            Ifom.read(fn_fom, false);
            RFLOAT kurt = 0., skew = 0.;
            Ifom.MDMainHeader.getValue(EMDL_IMAGE_STATS_AVG, skew);
            Ifom.MDMainHeader.getValue(EMDL_IMAGE_STATS_STDDEV, kurt);
            if (isnan(skew)) skew = 0.;
            if (isnan(kurt)) kurt = 0.;

            MDin.setValue(EMDL_MICROGRAPH_SCORE_KURTOSIS, kurt, imic);
            MDin.setValue(EMDL_MICROGRAPH_SCORE_SKEWNESS, skew, imic);
            MDin.setValue(EMDL_MICROGRAPH_AUTOPICK_FOM, fn_fom, imic);
            MDin.setValue(EMDL_MICROGRAPH_AUTOPICK_PSI, fn_psi, imic);
        }

        if (!do_skip_tracing)
        {

            FileName fn_pick = fn_root + "_" + fn_out + ".star";

            MetaDataTable MD;
            MD.read(fn_pick);
			long nr_pick = MD.numberOfObjects();
			total_nr_picked += nr_pick;
            MDin.setValue(EMDL_MLMODEL_GROUP_NR_PARTICLES, nr_pick, imic);

            MDcoords.addObject();
            MDcoords.setValue(EMDL_MICROGRAPH_NAME, fn_ori_micrographs[imic]);
            MDcoords.setValue(EMDL_MICROGRAPH_COORDINATES, fn_pick);
            MDcoords.setValue(EMDL_MICROGRAPH_AUTOPICK_FOM, fn_fom);
            nr_coord_files++;

        }

		if (verb > 0 && imic % 60 == 0) progress_bar(imic);

	}

    if (verb > 0) progress_bar(fn_ori_micrographs.size());

    // Make histograms of skewness and kurtosis of FOM values for all micrographs
    FileName fn_eps;
    std::vector<FileName> all_fn_eps;
    std::vector<RFLOAT> histX, histY;

    if (!do_skip_fom)
    {

        FileName fn_mics = fn_odir + "micrographs_" + fn_out + ".star";
        obsModel.save(MDin, fn_mics, "micrographs");
        if (verb > 0) std::cout << " Saved output micrograph STAR file with FOM images in: " << fn_mics << std::endl;

        CPlot2D *plot2De=new CPlot2D("Skewness of FOM for all micrographs");
        MDin.addToCPlot2D(plot2De, EMDL_UNDEFINED, EMDL_MICROGRAPH_SCORE_SKEWNESS, 1.);
        plot2De->SetDrawLegend(false);
        fn_eps = fn_odir + "all_FOM_skew.eps";
        plot2De->OutputPostScriptPlot(fn_eps);
        all_fn_eps.push_back(fn_eps);
        delete plot2De;
        if (MDin.numberOfObjects() > 3)
        {
            CPlot2D *plot2Df=new CPlot2D("");
            MDin.columnHistogram(EMDL_MICROGRAPH_SCORE_SKEWNESS,histX,histY,0, plot2Df);
            fn_eps = fn_odir + "histogram_FOM_skew.eps";
            plot2Df->SetTitle("Histogram of FOM skewness per micrograph");
            plot2Df->OutputPostScriptPlot(fn_eps);
            all_fn_eps.push_back(fn_eps);
            delete plot2Df;
        }

        CPlot2D *plot2Dg=new CPlot2D("Kurtosis of FOM for all micrographs");
        MDin.addToCPlot2D(plot2Dg, EMDL_UNDEFINED, EMDL_MICROGRAPH_SCORE_KURTOSIS, 1.);
        plot2Dg->SetDrawLegend(false);
        fn_eps = fn_odir + "all_FOM_kurt.eps";
        plot2Dg->OutputPostScriptPlot(fn_eps);
        all_fn_eps.push_back(fn_eps);
        delete plot2Dg;
        if (MDin.numberOfObjects() > 3)
        {
            CPlot2D *plot2Dh=new CPlot2D("");
            MDin.columnHistogram(EMDL_MICROGRAPH_SCORE_KURTOSIS,histX,histY,0, plot2Dh);
            fn_eps = fn_odir + "histogram_FOM_kurt.eps";
            plot2Dh->SetTitle("Histogram of FOM kurtosis per micrograph");
            plot2Dh->OutputPostScriptPlot(fn_eps);
            all_fn_eps.push_back(fn_eps);
            delete plot2Dh;
        }

    }

    if (!do_skip_tracing)
    {

        CPlot2D *plot2De=new CPlot2D("Nr of picked particles for all micrographs");
        MDin.addToCPlot2D(plot2De, EMDL_UNDEFINED, EMDL_MLMODEL_GROUP_NR_PARTICLES, 1.);
        plot2De->SetDrawLegend(false);
        fn_eps = fn_odir + "all_nrparts.eps";
        plot2De->OutputPostScriptPlot(fn_eps);
        all_fn_eps.push_back(fn_eps);
        delete plot2De;
        if (MDin.numberOfObjects() > 3)
        {
            CPlot2D *plot2Df=new CPlot2D("");
            MDin.columnHistogram(EMDL_MLMODEL_GROUP_NR_PARTICLES,histX,histY,0, plot2Df);
            fn_eps = fn_odir + "histogram_nrparts.eps";
            plot2Df->SetTitle("Histogram of nr of picked particles per micrograph");
            plot2Df->OutputPostScriptPlot(fn_eps);
            all_fn_eps.push_back(fn_eps);
            delete plot2Df;
        }

        FileName fn_coords = fn_odir + fn_out + ".star";

        MetaDataTable MDhead;
        MDhead.setIsList(true);
        MDhead.setName("general");
        MDhead.addObject();
        std::string picktype = "lines";
        MDhead.setValue(EMDL_MICROGRAPH_PICKTYPE, picktype);

        std::vector<MetaDataTable> MDins;
        MDins.push_back(MDhead);
        MDins.push_back(MDcoords);
        writeMultipleTablesToStar(MDins, fn_coords);

    }

    joinMultipleEPSIntoSinglePDF(fn_odir + "logfile.pdf", all_fn_eps);

}
