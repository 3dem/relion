#include "subtomo.h"
#include <src/jaz/tomography/dynamo/catalogue.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/projection/point_insertion.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/optics/aberrations_cache.h>
#include <src/time.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/math/Euler_angles_relion.h>
#include <mpi.h>
#include <iostream>

using namespace gravis;


SubtomoProgram::SubtomoProgram()
: do_not_write_any(false),
  run_from_MPI(false)
{

}

void SubtomoProgram::readBasicParameters(IOParser& parser)
{
	optimisationSet.read(
		parser,
		true,           // optimisation set
		true,   true,   // particles
		true,   true,   // tomograms
		true,   false,  // trajectories
		false,  false,  // manifolds
		false,  false); // reference

	int gen_section = parser.addSection("Reconstruction options");

	boxSize = textToInteger(parser.getOption("--b", "Binned projection box size"));
	cropSize = textToInteger(parser.getOption("--crop", "Output box size", "-1"));
	binning = textToDouble(parser.getOption("--bin", "Binning factor", "1"));
    do_stack2d = parser.checkOption("--stack2d", "Write out 2D stacks of cropped images for each particle, instead of pseudo-subtomograms");
	write_multiplicity = parser.checkOption("--multi", "Write out multiplicity volumes");
	SNR = textToDouble(parser.getOption("--SNR", "Assumed signal-to-noise ratio (negative means use a heuristic)", "-1"));
    min_frames = textToInteger(parser.getOption("--min_frames", "Minimum number of lowest-dose tilt series frames that needs to be inside the box", "1"));
    maxDose = textToDouble(parser.getOption("--max_dose", "Maximum dose (in e/A^2) of tilt series frames to be included (negative means use all frames)", "-1"));

	do_cone_weight = parser.checkOption("--cone_weight", "Weight down a double cone along Z");
	const double alpha = 0.5 * textToDouble(parser.getOption("--cone_angle", "Opening angle of the cone in degrees", "10"));
	cone_slope = sin(DEG2RAD(alpha));
	cone_sig0 = textToDouble(parser.getOption("--cone_sig0", "Cone width at Z = 0", "2"));

	do_circle_crop = !parser.checkOption("--no_circle_crop", "Do not crop 2D images to a circle");
	do_circle_precrop = parser.checkOption("--circle_precrop", "Crop 2D images to the large circle (--b) prior to CTF modulation");
	do_narrow_circle_crop = true;
	do_gridding_precorrection = parser.checkOption("--grid_precorr", "Perform gridding pre-correction on 2D images");

	taper = textToDouble(parser.getOption("--taper", "Taper against the sphere by this number of pixels", "5"));
	env_sigma = textToDouble(parser.getOption("--env", "Sigma of a Gaussian envelope applied before cropping", "-1"));

	do_whiten = parser.checkOption("--whiten", "Whiten the noise by flattening the power spectrum");
	do_center = !parser.checkOption("--no_center", "Do not subtract the mean from the voxel values");

	flip_value = !parser.checkOption("--no_ic", "Do not invert contrast (keep particles dark)");
    do_ctf = !parser.checkOption("--no_ctf", "Do not apply CTFs");
	write_combined = !parser.checkOption("--no_comb", "Do not write the concatenated CTF-multiplicity image");
	write_ctf = parser.checkOption("--ctf", "Write 3D CTFs");
	write_divided = parser.checkOption("--div", "Write CTF-corrected subtomograms");
	write_normalised = parser.checkOption("--nrm", "Write multiplicity-normalised subtomograms");

    apply_offsets = !parser.checkOption("--dont_apply_offsets", "By default, rlnOriginX/Y/ZAngst are combined with rlnCoordinateX/Y/Z to construct the particles in their refined translations. Use this argument to skip that.");
    apply_orientations = parser.checkOption("--apply_orientations", "rlnAngle<Rot/Tilt/Psi> are combined with rlnTomoSubtomogram<Rot/Tilt/Psi> to construct the particles in their refined orientations. This will also apply translations!");
    if (apply_orientations) apply_offsets = true;

	only_do_unfinished = parser.checkOption("--only_do_unfinished", "Only process undone subtomograms");

	write_float16  = parser.checkOption("--float16", "Write in half-precision 16 bit floating point numbers (MRC mode 12), instead of 32 bit (MRC mode 0).");


	diag = parser.checkOption("--diag", "Write out diagnostic information");

	num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));

	freqCutoffFract = textToDouble(parser.getOption("--cutoff_fract", "Ignore shells for which the dose weight falls below this value", "0.01"));

	outDir = parser.getOption("--o", "Output filename pattern");

	run_from_GUI = is_under_pipeline_control();

    do_real_subtomo = parser.checkOption("--real_subtomo", "Extract true subtomograms and write out projections of those out as 2D stacks");

}

void SubtomoProgram::readParameters(int argc, char *argv[])
{
	IOParser parser;

	parser.setCommandLine(argc, argv);

	readBasicParameters(parser);

	do_sum_all = parser.checkOption("--sum", "Sum up all subtomograms (for debugging)");
	do_not_write_any = parser.checkOption("--no_writing", "Do not write out any files, only a sum");

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	if (do_gridding_precorrection)
	{
		do_narrow_circle_crop = true;
	}

	outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);
}

void SubtomoProgram::run()
{
	TomogramSet tomogramSet(optimisationSet.tomograms, true);

	ParticleSet particleSet(optimisationSet.particles, optimisationSet.trajectories, true, &tomogramSet);
	std::vector<std::vector<ParticleIndex> > particles = particleSet.splitByTomogram(tomogramSet, true);

    if (cropSize < 0) cropSize = boxSize;
	
	const long int s2D = boxSize;

	const long int s3D = cropSize;
	const long int sh3D = s3D / 2 + 1;

	const long int s02D = (int)(binning * s2D + 0.5);

	const double relative_box_scale = cropSize / (double) boxSize;
	const double binned_pixel_size = binning * particleSet.getTiltSeriesPixelSize(0);

    initialise(particleSet, particles, tomogramSet);

    BufferedImage<float> sum_data, sum_weights;

    if (do_stack2d) do_sum_all = false;

	if (do_sum_all)
	{
		sum_data.resize(s3D,s3D,s3D);
		sum_data.fill(0.0);

		sum_weights.resize(sh3D,s3D,s3D);
		sum_weights.fill(0.0);
	}

	AberrationsCache aberrationsCache(particleSet.optTable, s2D, binned_pixel_size);


	std::vector<int> tomoIndices = ParticleSet::enumerate(particles);

	processTomograms(
		tomoIndices,
		tomogramSet,
		particleSet,
		particles,
		aberrationsCache,
		s02D,
		s2D,
		s3D,
		relative_box_scale,
		1,
		sum_data,
		sum_weights);


	if (do_sum_all)
	{
		const double pixel_size = binning * tomogramSet.getTiltSeriesPixelSize(0);

		sum_data.write(outDir + "sum_data.mrc", pixel_size);
		Centering::fftwHalfToHumanFull(sum_weights).write(outDir + "sum_weight.mrc", pixel_size);

		BufferedImage<float> dataImgDivRS(s3D,s3D,s3D);
		dataImgDivRS.fill(0.0);

		if (SNR > 0.0)
		{
			Reconstruction::ctfCorrect3D_Wiener(
				sum_data, sum_weights, dataImgDivRS,
				1.0 / SNR, num_threads);
		}
		else
		{
			Reconstruction::ctfCorrect3D_heuristic(
				sum_data, sum_weights, dataImgDivRS,
				0.001, num_threads);
		}

		dataImgDivRS.write(outDir + "sum_div.mrc", pixel_size);
	}
}

void SubtomoProgram::initialise(
		ParticleSet& particleSet,
		const std::vector<std::vector<ParticleIndex>>& particles,
		const TomogramSet& tomogramSet,
        bool verbose)
{
	const int tc = tomogramSet.size();

	int firstGoodTomo = 0;

	for (int t = 0; t < tc; t++)
	{
		if (particles[t].size() > 0)
		{
			firstGoodTomo = t;
			break;
		}
	}

	const std::string firstName = particleSet.getName(particles[firstGoodTomo][0]);

	directoriesPerTomogram = firstName.find_first_of('/') == std::string::npos;

    if (verbose)
    {
        if (directoriesPerTomogram)
        {
            Log::print("No slashes found in first particle name: creating subdirectories for each tomogram");
        }
        else
        {
            Log::print("Slash found in first particle name: not creating subdirectories for each tomogram");
        }
    }

	for (int t = 0; t < tc; t++)
	{
		if (particles[t].size() > 0)
		{
			ZIO::ensureParentDir(getOutputFilename(
				particles[t][0], t, particleSet, tomogramSet));
		}
	}

    if (verbose) writeParticleSet(particleSet, particles, tomogramSet);
}

std::string SubtomoProgram::getOutputFilename(
		ParticleIndex p,
		int tomogramIndex,
		const ParticleSet& particleSet,
		const TomogramSet& tomogramSet)
{
	if (directoriesPerTomogram)
	{
		return outDir + "Subtomograms/"
				+ tomogramSet.getTomogramName(tomogramIndex) + "/"
				+ particleSet.getName(p);
	}
	else
	{
		return outDir + "Subtomograms/" + particleSet.getName(p);
	}
}

void SubtomoProgram::writeParticleSet(
		const ParticleSet& particleSet,
		const std::vector<std::vector<ParticleIndex>>& particles,
		const TomogramSet& tomogramSet)
{
	const int tc = particles.size();

	ParticleSet copy2d, copy = particleSet;
	copy.clearParticles();
    copy.is_stack2d = do_stack2d;

    if (do_real_subtomo)
    {
        copy2d = particleSet;
        copy2d.clearParticles();
        copy2d.is_stack2d = false;
    }

	int particles_removed = 0;

    RFLOAT real_subtomo_binning = 0.;
    int cropSize_tomogram;
	for (int t = 0; t < tc; t++)
	{
		const int pc = particles[t].size();

        int mintilt_idx;

        if (pc == 0) continue;

        Tomogram tomogram = tomogramSet.loadTomogram(t, false);

        if (do_real_subtomo)
        {
            RFLOAT tomogram_binning;
            tomogramSet.globalTable.getValue(EMDL_TOMO_TOMOGRAM_BINNING, tomogram_binning, t);
            if (t == 0) real_subtomo_binning = tomogram_binning;
            else if (real_subtomo_binning != tomogram_binning) REPORT_ERROR("ERROR: not all tomograms have the same binning; can't do real subtomos");
            cropSize_tomogram = ROUND(cropSize * binning / tomogram_binning);
            if (cropSize_tomogram%2 != 0)
            {
                if ( (cropSize * binning / tomogram_binning) - cropSize_tomogram >= 0. ) cropSize_tomogram++;
                else cropSize_tomogram--;
            }

        }
		for (int p = 0; p < pc; p++)
		{
			const ParticleIndex part_id = particles[t][p];

			const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
						part_id, tomogram.frameCount, tomogram.centre, tomogram.optics.pixelSize, !apply_offsets);

            //SHWS 24jan2024 only for debugging: // int my_min_frames = (min_frames < 0 ) ? tomogram.frameCount / 2 : min_frames;
			//SHWS 24jan2024 only for debugging: // if (tomogram.isVisibleFirstFrames(traj, boxSize / 2.0, my_min_frames))
			std::vector<bool> isVisible;
            if (tomogram.getVisibilityMinFramesMaxDose(traj, binning * cropSize / 2.0, maxDose, min_frames, isVisible))
            {

                // Get all the particle metadata
                const ParticleIndex new_id = copy.addParticle(particleSet, part_id);

                // Also set isVisible in the output particle STAR file
                if (do_stack2d)
                {
                    std::vector<int> isVisibleInt(isVisible.size(), 0);
                    for (int f = 0; f < tomogram.frameCount; f++)
                        if (isVisible[f]) isVisibleInt[f] = 1;
                    copy.partTable.setValue(EMDL_TOMO_VISIBLE_FRAMES, isVisibleInt);
                }

				const int opticsGroup = particleSet.getOpticsGroup(part_id);
				const double tiltSeriesPixelSize = particleSet.getTiltSeriesPixelSize(opticsGroup);

				const std::string filenameRoot = getOutputFilename(
					part_id, t, particleSet, tomogramSet);

                std::string outData = (do_stack2d) ? filenameRoot + "_stack2d.mrcs" : filenameRoot + "_data.mrc";
                std::string outWeight = (do_stack2d) ? "" : filenameRoot + "_weights.mrc";

                copy.setImageFileNames(outData, outWeight, new_id);

                if (apply_offsets)
                {
                    const d3Vector pos = particleSet.getPosition(part_id, tomogram.centre, true);
                    copy.setParticleOffset(new_id, d3Vector(0,0,0));
                    copy.setParticleCoordDecenteredPixel(new_id, pos, tomogram.centre, tiltSeriesPixelSize);
                }

				if (apply_orientations)
				{
                    d3Matrix A = particleSet.getMatrix3x3(part_id);
					const gravis::d3Vector ang = Euler::matrixToAngles(A);

                    copy.partTable.setValue(EMDL_TOMO_SUBTOMOGRAM_ROT, RAD2DEG(ang[0]), new_id.value);
                    copy.partTable.setValue(EMDL_TOMO_SUBTOMOGRAM_TILT, RAD2DEG(ang[1]), new_id.value);
                    copy.partTable.setValue(EMDL_TOMO_SUBTOMOGRAM_PSI, RAD2DEG(ang[2]), new_id.value);

                    copy.partTable.setValue(EMDL_ORIENT_ROT, 0.0, new_id.value);
                    copy.partTable.setValue(EMDL_ORIENT_TILT, 0.0, new_id.value);
                    copy.partTable.setValue(EMDL_ORIENT_PSI, 0.0, new_id.value);
				}

                if (do_real_subtomo)
                {
                    // Also make a particle star file for 2D classification
                    mintilt_idx = tomogramSet.getImageIndexWithSmallestVisibleTiltAngle(t, isVisible);

                    copy2d.partTable.addObject( copy.partTable.getObject(new_id.value) );

                    std::string outData2d = integerToString(mintilt_idx+1) + "@" + outData;
                    copy2d.partTable.setValue(EMDL_IMAGE_NAME, outData2d);
                }
			}
			else
			{
				particles_removed++;
			}
		}
	}

	if (particles_removed == 1)
	{
		Log::warn("One particle was removed because it was too close to the edge in all images.");
	}
	else if (particles_removed > 1)
	{
		Log::warn(ZIO::itoa(particles_removed)+" particles were removed because they were too close to the edge in all images.");
	}

	for (int og = 0; og < copy.numberOfOpticsGroups(); og++)
	{

        bool is_premultiplied = (do_stack2d) ? do_ctf : true;
        if (do_real_subtomo)
        {
            is_premultiplied = false;
            copy.optTable.setValue(EMDL_OPTIMISER_DATA_ARE_CTF_CORRECTED, true, og);
        }
		copy.optTable.setValue(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, is_premultiplied, og);
		int datadim = (do_stack2d || do_real_subtomo) ? 2 : 3;
        copy.optTable.setValue(EMDL_IMAGE_DIMENSIONALITY, datadim, og);

        RFLOAT mybinning = (do_real_subtomo) ? real_subtomo_binning : binning;
		copy.optTable.setValue(EMDL_TOMO_SUBTOMOGRAM_BINNING, mybinning, og);
        const double ps_img = copy.optTable.getDouble(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, og);
        const double ps_out = mybinning * ps_img;
		copy.optTable.setValue(EMDL_IMAGE_PIXEL_SIZE, ps_out, og);
		int mysize = (do_real_subtomo) ? cropSize_tomogram : cropSize;
        copy.optTable.setValue(EMDL_IMAGE_SIZE, mysize, og);


	}

    if (copy.partTable.containsLabel(EMDL_IMAGE_COORD_X)) copy.partTable.deactivateLabel(EMDL_IMAGE_COORD_X);
    if (copy.partTable.containsLabel(EMDL_IMAGE_COORD_Y)) copy.partTable.deactivateLabel(EMDL_IMAGE_COORD_Y);
    if (copy.partTable.containsLabel(EMDL_IMAGE_COORD_Z)) copy.partTable.deactivateLabel(EMDL_IMAGE_COORD_Z);
    // remove CTF image name that could be there from old pseudo-subtomo jobs if we're writing 2D stacks
    if (do_stack2d && copy.partTable.containsLabel(EMDL_CTF_IMAGE))  copy.partTable.deactivateLabel(EMDL_CTF_IMAGE);

    copy.write(outDir + "particles.star");

    if (do_real_subtomo)
    {

        if (copy2d.partTable.containsLabel(EMDL_IMAGE_COORD_X)) copy2d.partTable.deactivateLabel(EMDL_IMAGE_COORD_X);
        if (copy2d.partTable.containsLabel(EMDL_IMAGE_COORD_Y)) copy2d.partTable.deactivateLabel(EMDL_IMAGE_COORD_Y);
        if (copy2d.partTable.containsLabel(EMDL_IMAGE_COORD_Z)) copy2d.partTable.deactivateLabel(EMDL_IMAGE_COORD_Z);
        if (copy2d.partTable.containsLabel(EMDL_CTF_IMAGE))  copy2d.partTable.deactivateLabel(EMDL_CTF_IMAGE);

        copy2d.optTable.clear();
        copy2d.optTable = copy.optTable;
        copy2d.write(outDir + "particles_for_class2d.star");
    }


	if (copy.hasMotion && particles_removed > 0)
	{
		copy.writeTrajectories(outDir + "motion.star");
		optimisationSet.trajectories = outDir + "motion.star";
	}

	optimisationSet.particles = outDir + "particles.star";
	optimisationSet.write(outDir + "optimisation_set.star");
}

BufferedImage<float> SubtomoProgram::extractSubtomogramsAndReProject(
        ParticleIndex part_id, MultidimArray<RFLOAT> &recTomo,
        const Tomogram& tomogram, const ParticleSet &particleSet,
        const std::vector<bool> &isVisible, RFLOAT tomogram_angpix)
{

    // get decentered coordinates of the particle (in tilt series pixels)
    d3Vector pos = particleSet.getPosition(part_id, tomogram.centre, apply_offsets);
    // convert from tilt series pixel to tomogram pixel
    pos *= tomogram.optics.pixelSize / tomogram_angpix;

    const float sign = flip_value ? -1.f : 1.f;

    // Center the subtomo in a box that is as wide as the cropSize in Angstroms, but only fill the center sphere. Leave the rest at zero
    int cropSize_tomogram = ROUND(cropSize * binning * tomogram.optics.pixelSize / tomogram_angpix);
    if (cropSize_tomogram%2 != 0)
    {
        if ( (cropSize * binning * tomogram.optics.pixelSize / tomogram_angpix) - cropSize_tomogram >= 0. ) cropSize_tomogram++;
        else cropSize_tomogram--;
    }
    //std::cerr << " tomogram_angpix= " << tomogram_angpix << std::endl;
    //std::cerr << " cropSize_tomogram= " << cropSize_tomogram << " cropSize  * tomogram.optics.pixelSize / tomogram_angpix = " << cropSize * tomogram.optics.pixelSize / tomogram_angpix << std::endl;

    MultidimArray<RFLOAT> subtom(cropSize_tomogram, cropSize_tomogram, cropSize_tomogram);
    subtom.setXmippOrigin();
    // Only extract central sphere with diameter cropSize_tomogram!
    RFLOAT edge_width = 2.;
    long int r2_max_smaller = ROUND((cropSize_tomogram - 2*edge_width) * (cropSize_tomogram - 2*edge_width) / 4.);
    long int r2_max = ROUND(cropSize_tomogram * cropSize_tomogram / 4.);
    RFLOAT sum= 0., sum2 = 0., nn= 0.;
    for (long int k= FIRST_XMIPP_INDEX(cropSize_tomogram); k<= LAST_XMIPP_INDEX(cropSize_tomogram); k++)
        for (long int i= FIRST_XMIPP_INDEX(cropSize_tomogram); i<= LAST_XMIPP_INDEX(cropSize_tomogram); i++)
            for (long int j= FIRST_XMIPP_INDEX(cropSize_tomogram); j<= LAST_XMIPP_INDEX(cropSize_tomogram); j++)
            {
                int kp = k + pos.z;
                int ip = i + pos.y;
                int jp = j + pos.x;
                long int r2 = k*k + i*i + j*j;
                if (r2 <= r2_max &&
                        jp > 0  && jp < XSIZE(recTomo)  &&
                        ip > 0  && ip < YSIZE(recTomo)  &&
                        kp > 0  && kp < ZSIZE(recTomo))
                {
                    A3D_ELEM(subtom, k, i, j) = sign * DIRECT_A3D_ELEM(recTomo, kp, ip, jp);
                    if (r2 >= r2_max_smaller)
                    {
                        sum += A3D_ELEM(subtom, k, i, j);
                        sum2 += A3D_ELEM(subtom, k, i, j) * A3D_ELEM(subtom, k, i, j);
                        nn += 1.;
                    }
                }
                else
                    A3D_ELEM(subtom, k, i, j) = 0.;
            }

    // Calculate the background mean in a shell of edge_width around the particle
    RFLOAT bg_mean = sum / nn;
    RFLOAT stddev = sqrt( (sum2 / nn - bg_mean * bg_mean));
    if (stddev < 0.00000001)
        REPORT_ERROR("ERROR: all-zero subtomogram should not happen here!");

    // Apply background subtraction and a slightly softer edge
    for (long int k= FIRST_XMIPP_INDEX(cropSize_tomogram); k<= LAST_XMIPP_INDEX(cropSize_tomogram); k++)
        for (long int i= FIRST_XMIPP_INDEX(cropSize_tomogram); i<= LAST_XMIPP_INDEX(cropSize_tomogram); i++)
            for (long int j= FIRST_XMIPP_INDEX(cropSize_tomogram); j<= LAST_XMIPP_INDEX(cropSize_tomogram); j++)
            {
                int kp = k + pos.z;
                int ip = i + pos.y;
                int jp = j + pos.x;
                long int r2 = k*k + i*i + j*j;
                if (r2 <= r2_max &&
                        jp > 0  && jp < XSIZE(recTomo)  &&
                        ip > 0  && ip < YSIZE(recTomo)  &&
                        kp > 0  && kp < ZSIZE(recTomo))
                {
                    A3D_ELEM(subtom, k, i, j) -= bg_mean;
                    //A3D_ELEM(subtom, k, i, j) /= stddev;
                    // Also apply a slightly softer edge
                    if (r2 >= r2_max_smaller)
                    {
                        RFLOAT frac = (cropSize_tomogram/2. - sqrt(r2))/edge_width;
                        A3D_ELEM(subtom, k, i, j) *= frac;
                    }
                }
            }

    //Image<RFLOAT> It;
    //It()= subtom;
    //It.write("subtom.spi");
    //exit(0);

    MultidimArray<RFLOAT> dummy;
    Projector projector(cropSize_tomogram);
    projector.computeFourierTransformMap(subtom, dummy);

    MultidimArray<RFLOAT> img(cropSize_tomogram, cropSize_tomogram);
    MultidimArray<Complex> F2D;
    FourierTransformer transformer;
    transformer.setReal(img);
    transformer.getFourierAlias(F2D);

    BufferedImage<float> resultImg(cropSize_tomogram, cropSize_tomogram, tomogram.frameCount);
    for (int f = 0; f < tomogram.frameCount; f++)
    {

        if (!isVisible[f]) continue;

        d4Matrix Aproj = tomogram.projectionMatrices[f];
        Matrix2D<RFLOAT> A(3,3);
        for (int row= 0; row < 3; row++)
            for (int col = 0; col < 3; col++)
                MAT_ELEM(A, row, col) = Aproj(row, col);

        // Get the 2D slice out of the 3D Fourier transform
        F2D.initZeros();
        projector.get2DFourierTransform(F2D, A);
        shiftImageInFourierTransform(F2D, F2D, cropSize_tomogram, cropSize_tomogram/2, cropSize_tomogram/2);
        transformer.inverseFourierTransform();

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(img)
        {
            resultImg(j, i, f) = DIRECT_A2D_ELEM(img, i, j);
        }

    }

    return resultImg;

}

void SubtomoProgram::processTomograms(
		const std::vector<int>& tomoIndices,
		const TomogramSet& tomogramSet,
		const ParticleSet& particleSet,
		const std::vector<std::vector<ParticleIndex>>& particles,
		const AberrationsCache& aberrationsCache,
		long int s02D,
		long int s2D,
		long int s3D,
		double relative_box_scale,
		int verbosity,
		BufferedImage<float>& sum_data,
		BufferedImage<float>& sum_weights )
{
	const int tc = tomoIndices.size();
	const int sh2D = s2D / 2 + 1;
	const int sh3D = s3D / 2 + 1;

	for (int tt = 0; tt < tc; tt++)
	{
		const int t = tomoIndices[tt];

		const int pc = particles[t].size();
		if (pc == 0) continue;

		if (run_from_GUI && pipeline_control_check_abort_job())
		{
			if (run_from_MPI)
			{
				MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_ABORTED);
				exit(RELION_EXIT_ABORTED);
			}
			else
			{
				exit(RELION_EXIT_ABORTED);
			}
		}

		if (verbosity > 0)
		{
			Log::beginSection("Tomogram " + ZIO::itoa(tt+1) + " / " + ZIO::itoa(tc));
			Log::print("Loading");
		}

		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		tomogram.validateParticleOptics(particles[t], particleSet);

        // If using the real_subtomo approach, then need to read in the reconstructed tomogram volume
        Image<RFLOAT> recTomo;
        RFLOAT tomogram_angpix, tomogram_binning;
        if (do_real_subtomo)
        {
            FileName fn_tomo, fn_tomo2="";
            if (tomogramSet.globalTable.containsLabel(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_FILE_NAME))
            {
                tomogramSet.globalTable.getValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_FILE_NAME, fn_tomo, t);
            }
            else if (tomogramSet.globalTable.containsLabel(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_HALF1_FILE_NAME))
            {
                tomogramSet.globalTable.getValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_HALF1_FILE_NAME, fn_tomo, t);
                tomogramSet.globalTable.getValue(EMDL_TOMO_RECONSTRUCTED_TOMOGRAM_HALF2_FILE_NAME, fn_tomo2, t);
            }
            else
                REPORT_ERROR("ERROR: cannot find rlnTomoReconstructedTomogram or rlnTomoReconstructedTomogramHalf1 in tomogram star file");

            tomogramSet.globalTable.getValue(EMDL_TOMO_TOMOGRAM_BINNING, tomogram_binning, t);
            recTomo.read(fn_tomo);
            if (fn_tomo2 != "")
            {
                Image<RFLOAT> recTomo2;
                recTomo2.read(fn_tomo2);
                recTomo() += recTomo2();
            }
            tomogram_angpix =  tomogram.optics.pixelSize * tomogram_binning;
        }

		const int fc = tomogram.frameCount;

		particleSet.checkTrajectoryLengths(particles[t], fc, "subtomo");

		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s2D, binning);
		BufferedImage<float> noiseWeights;

		if (do_whiten)
		{
			noiseWeights = tomogram.computeNoiseWeight(s2D, binning);
		}

		BufferedImage<int> xRanges = tomogram.findDoseXRanges(doseWeights, freqCutoffFract);

		const int inner_thread_num = 1;
		const int outer_thread_num = num_threads / inner_thread_num;

		// @TODO: define input and output pixel sizes!

		const double binnedPixelSize = tomogram.optics.pixelSize * binning;
 		if (verbosity > 0)
		{
            Log::beginProgress(
				"Extracting particles",
				(int)ceil(pc/(double)outer_thread_num));
		}

		omp_lock_t writelock;
		if (do_sum_all) omp_init_lock(&writelock);

		#pragma omp parallel for num_threads(outer_thread_num)
		for (int p = 0; p < pc; p++) {
            const int th = omp_get_thread_num();

            if (verbosity > 0 && th == 0) {
                Log::updateProgress(p);
            }

            const ParticleIndex part_id = particles[t][p];

            const std::string filenameRoot = getOutputFilename(
                    part_id, t, particleSet, tomogramSet);

            std::string outData = (do_stack2d) ? filenameRoot + "_stack2d.mrcs" : filenameRoot + "_data.mrc";
            std::string outWeight = (do_stack2d) ? "" : filenameRoot + "_weights.mrc";
            std::string outCTF = filenameRoot + "_CTF2.mrc";
            std::string outDiv = filenameRoot + "_div.mrc";
            std::string outMulti = filenameRoot + "_multi.mrc";
            std::string outNrm = filenameRoot + "_data_nrm.mrc";
            std::string outWeightNrm = filenameRoot + "_CTF2_nrm.mrc";

            if (only_do_unfinished && ZIO::fileExists(outData)) {
                continue;
            }

            const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
                    part_id, fc, tomogram.centre, tomogram.optics.pixelSize, !apply_offsets);

            std::vector<bool> isVisible;
            if (!tomogram.getVisibilityMinFramesMaxDose(traj, binning * cropSize / 2.0, maxDose, min_frames, isVisible))
                continue;

            std::vector<d4Matrix> projCut(fc), projPart(fc);

            BufferedImage<fComplex> particleStack = BufferedImage<fComplex>(sh2D, s2D, fc);
            BufferedImage<float> weightStack(sh2D, s2D, fc);

            if (do_real_subtomo)
            {

                // This will extract the true subtomograms and calculate their FTs in the directions of the tilt series
                BufferedImage<float> subtomo_reprojs = extractSubtomogramsAndReProject(part_id, recTomo(),
                                                                tomogram, particleSet, isVisible, tomogram_angpix);
                BufferedImage<float> visible_reprojs = NewStackHelper::getVisibleSlices(subtomo_reprojs, isVisible);
                visible_reprojs.write(outData, tomogram_angpix, write_float16);

            }
            else
            {

                TomoExtraction::extractAt3D_Fourier(
                        tomogram.stack, s02D, binning, tomogram, traj, isVisible,
                        particleStack, projCut, inner_thread_num, do_circle_precrop);

                if (!do_ctf) weightStack.fill(1.f);

                const int og = particleSet.getOpticsGroup(part_id);

                const BufferedImage<double> *gammaOffset =
                        aberrationsCache.hasSymmetrical ? &aberrationsCache.symmetrical[og] : 0;

                const float sign = flip_value ? -1.f : 1.f;
                for (int f = 0; f < fc; f++)
                {
                    if (!isVisible[f]) continue;

                    d3Matrix A;

                    if (apply_orientations) {
                        A = particleSet.getMatrix3x3(part_id);
                    } else {
                        A = particleSet.getSubtomogramMatrix(part_id);
                    }

                    projPart[f] = projCut[f] * d4Matrix(A);

                    if (do_ctf) {
                        const d3Vector pos = particleSet.getPosition(part_id, tomogram.centre, apply_offsets);

                        CTF ctf = tomogram.getCtf(f, pos);
                        BufferedImage<float> ctfImg(sh2D, s2D);
                        ctf.draw(s2D, s2D, binnedPixelSize, gammaOffset, &ctfImg(0, 0, 0));

                        // Apply doseWeigths until Nyquist frequency! Otherwise, convolution artefacts when do_circle_crop invFFT/FFT below
                        for (int y = 0; y < s2D; y++) {
                            for (int x = 0; x < sh2D; x++) {
                                const double c = ctfImg(x, y) * doseWeights(x, y, f);

                                particleStack(x, y, f) *= sign * c;
                                weightStack(x, y, f) = c * c;
                            }
                        }
                    } // end if do_ctf
                } // end for f

                // If we're not doing CTF premultiplication, we may still want to invert the contrast
                if (!do_ctf) particleStack *= sign;

                aberrationsCache.correctObservations(particleStack, og);

                if (do_whiten) {
                    particleStack *= noiseWeights;
                    weightStack *= noiseWeights;
                }

                // Make sure output greyscale of 2D stacks does not depend on binning
                if (do_stack2d) particleStack/= (float)binning;

                const int boundary = (boxSize - cropSize) / 2;

                BufferedImage<float> particlesRS;
                if (do_gridding_precorrection || do_circle_crop || do_stack2d) {

                    particlesRS = NewStackHelper::inverseFourierTransformStack(particleStack);

                }

                if (do_circle_crop) {
                    const double crop_boundary = do_narrow_circle_crop ? boundary : 0.0;
                    TomoExtraction::cropCircle(particlesRS, crop_boundary, 5, num_threads);
                }

                if (do_gridding_precorrection) {
                    TomoExtraction::griddingPreCorrect(particlesRS, boundary, num_threads);
                }

                if (do_stack2d) {

                    BufferedImage<float> cropParticlesRS = Padding::unpadCenter2D_full(particlesRS, boundary);
                    BufferedImage<float> cropParticlesRS2 = NewStackHelper::getVisibleSlices(cropParticlesRS, isVisible);
                    cropParticlesRS2.write(outData, binnedPixelSize, write_float16);

                } else {

                    if (do_gridding_precorrection || do_circle_crop)
                    {
                        particleStack = NewStackHelper::FourierTransformStack(particlesRS);
                    }

                    // Write 3D subtomograms
                    BufferedImage<fComplex> dataImgFS(sh3D, s3D, s3D);
                    dataImgFS.fill(fComplex(0.0, 0.0));

                    BufferedImage<float> ctfImgFS(sh3D, s3D, s3D),
                            dataImgRS(s3D, s3D, s3D), dataImgDivRS(s3D, s3D, s3D),
                            multiImageFS(sh3D, s3D, s3D);

                    ctfImgFS.fill(0.0);
                    dataImgRS.fill(0.0);
                    dataImgDivRS.fill(0.0);

                    for (int f = 0; f < fc; f++)
                    {
                        if (isVisible[f])
                        {
                            FourierBackprojection::backprojectSlice_forward_with_multiplicity(
                                    &xRanges(0, f),
                                    particleStack.getSliceRef(f),
                                    weightStack.getSliceRef(f),
                                    projPart[f] * relative_box_scale,
                                    dataImgFS,
                                    ctfImgFS,
                                    multiImageFS);
                        }
                    }

                    Centering::shiftInSitu(dataImgFS);

                    // correct FT scale after the implicit cropping:
                    float normfft = (float) s2D / ((float) s3D * sqrt((float)(binning)) );
                    if (fabs(normfft - 1.) > 0.0001)
                        dataImgFS *= normfft;

                    FFT::inverseFourierTransform(dataImgFS, dataImgRS, FFT::Both);

                    if (do_cone_weight) {
                        FFT::FourierTransform(dataImgRS, dataImgFS);

                        d3Matrix R = particleSet.getMatrix3x3(part_id);

                        for (int z = 0; z < s3D; z++)
                            for (int y = 0; y < s3D; y++)
                                for (int x = 0; x < sh3D; x++) {
                                    const d3Vector p0(
                                            x,
                                            y < s3D / 2 ? y : y - s3D,
                                            z < s3D / 2 ? z : z - s3D);

                                    const d3Vector p = R * p0;

                                    const double rho = sqrt(p.x * p.x + p.y * p.y);
                                    const double t = rho / (std::abs(p.z) * cone_slope + cone_sig0);

                                    const double m = 1.0 - exp(-0.5 * t * t);

                                    dataImgFS(x, y, z) *= m;
                                    ctfImgFS(x, y, z) *= m;
                                    multiImageFS(x, y, z) *= m; // apply to both multiplicity and weight?
                                }

                        FFT::inverseFourierTransform(dataImgFS, dataImgRS);
                    }

                    // What if we didn't? The 2D image is already tapered.
                    //Reconstruction::taper(dataImgRS, taper, do_center, inner_thread_num);

                    if (do_sum_all) {
                        omp_set_lock(&writelock);

                        sum_data += dataImgRS;
                        sum_weights += ctfImgFS;

                        omp_unset_lock(&writelock);
                    }

                    if (do_not_write_any) continue;


                    dataImgRS.write(outData, binnedPixelSize, write_float16);

                    if (write_combined) {
                        BufferedImage<float> ctfAndMultiplicity(sh3D, s3D, 2 * s3D);
                        ctfAndMultiplicity.getSlabRef(0, s3D).copyFrom(ctfImgFS);
                        ctfAndMultiplicity.getSlabRef(s3D, s3D).copyFrom(multiImageFS);

                        ctfAndMultiplicity.write(outWeight, 1.0 / binnedPixelSize, write_float16);
                    }

                    if (write_ctf) {
                        Centering::fftwHalfToHumanFull(ctfImgFS).write(outCTF, 1.0 / binnedPixelSize, write_float16);
                    }

                    if (write_multiplicity) {
                        Centering::fftwHalfToHumanFull(multiImageFS).write(outMulti, 1.0 / binnedPixelSize, write_float16);
                    }

                    if (write_normalised) {
                        BufferedImage<float> ctfImgFSnrm = ctfImgFS;
                        BufferedImage<fComplex> dataImgCorrFS;

                        FFT::FourierTransform(dataImgRS, dataImgCorrFS, FFT::Both);

                        for (long int i = 0; i < ctfImgFSnrm.getSize(); i++) {
                            const float n = multiImageFS[i];
                            ctfImgFSnrm[i] = n > 0.f ? ctfImgFS[i] / n : 0.f;
                            dataImgCorrFS[i] = n > 0.f ? dataImgCorrFS[i] / n : fComplex(0.f, 0.f);
                        }

                        FFT::inverseFourierTransform(dataImgCorrFS, dataImgDivRS, FFT::Both);

                        dataImgDivRS.write(outNrm, binnedPixelSize, write_float16);
                        Centering::fftwHalfToHumanFull(ctfImgFSnrm).write(outWeightNrm, 1.0 / binnedPixelSize,
                                                                          write_float16);
                    }

                    if (write_divided) {
                        if (SNR > 0.0) {
                            Reconstruction::ctfCorrect3D_Wiener(
                                    dataImgRS, ctfImgFS, dataImgDivRS,
                                    1.0 / SNR, inner_thread_num);
                        } else {
                            Reconstruction::ctfCorrect3D_heuristic(
                                    dataImgRS, ctfImgFS, dataImgDivRS,
                                    0.001, inner_thread_num);
                        }

                        Reconstruction::taper(dataImgDivRS, taper, do_center, inner_thread_num);
                        dataImgDivRS.write(outDiv, binnedPixelSize, write_float16);
                    }
                } // end if do_stack2d
            } // end if do_real_subtomo
        } // end loop particles p

		if (verbosity > 0)
		{
			Log::endProgress();
			Log::endSection(); // tomogram
		}
	}
}

BufferedImage<float> SubtomoProgram::cropAndTaper(const BufferedImage<float>& imgFS, int boundary, int num_threads) const
{
	BufferedImage<fComplex> ctfImgFS_complex = FFT::toComplex(imgFS);
	
	BufferedImage<float> ctfImgRS;
	FFT::inverseFourierTransform(ctfImgFS_complex, ctfImgRS);
	
	ctfImgRS = Centering::fftwFullToHumanFull(ctfImgRS);
	ctfImgRS = Padding::unpadCenter3D_full(ctfImgRS, boundary);

    Reconstruction::GaussEnvelope(ctfImgRS, env_sigma, do_center, num_threads);
	Reconstruction::taper(ctfImgRS, taper, do_center, num_threads);
	
	BufferedImage<float> ctfImgRS_cent = Centering::humanFullToFftwFull(ctfImgRS);
	
	FFT::FourierTransform(ctfImgRS_cent, ctfImgFS_complex);
	
	return FFT::toReal(ctfImgFS_complex);
}
