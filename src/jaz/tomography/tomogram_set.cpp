#include "tomogram_set.h"
#include "motion/Fourier_2D_deformation.h"
#include "motion/spline_2D_deformation.h"
#include "motion/linear_2D_deformation.h"
#include <fstream>
#include <sstream>
#include <src/error.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/image_file_helper.h>

using namespace gravis;


TomogramSet::TomogramSet()
{
	globalTable.setName("global");
}

TomogramSet::TomogramSet(FileName filename, bool verbose)
{
    bool success = read(filename, verbose);
}

bool TomogramSet::read(FileName filename, bool verbose)
{

    globalTable.read(filename, "global");
    globalTable.setName("global");

    const int tc = globalTable.numberOfObjects();
    tomogramTables.resize(tc);
    if (tc == 0) return false;

    if (!globalTable.containsLabel(EMDL_TOMO_NAME))
    {
        REPORT_ERROR("ERROR: input starfile for TomogramSet " + filename + " does not contain rlnTomoName label ");
    }

    if (!globalTable.containsLabel(EMDL_TOMO_TILT_SERIES_STARFILE))
    {
        if (verbose) std::cerr << "Warning: " << filename
                  << " does not have rlnTomoTiltSeriesStarFile labels. It may be written in an old format. If so, will try to convert ..."
                  << std::endl;

        // This may be a tomograms.star file in the old, original relion-4.0 format. Try to convert
        FileName mydir = filename.beforeLastOf("/");

        std::ifstream ifs(filename);
        if (!ifs)
        {
            REPORT_ERROR_STR("TomogramSet::TomogramSet: Unable to read " << filename);
        }
        else
        {
            std::vector<MetaDataTable> allTables = MetaDataTable::readAll(ifs, tc+1);
            for (int t = 0; t < tc; t++)
            {
                FileName expectedNewName = globalTable.getString(EMDL_TOMO_NAME, t);
                FileName name = allTables[t+1].getName();


                if (name != expectedNewName)
                {
                    REPORT_ERROR_STR("TomogramSet::TomogramSet: file is corrupted " << filename);
                }

                tomogramTables[t] = allTables[t+1];
                tomogramTables[t].setName(expectedNewName);

                // Also set the name for the titlseries STAR file
                FileName fn_star = filename.beforeLastOf("/") + "/tilt_series/" + expectedNewName + ".star";
                globalTable.setValue(EMDL_TOMO_TILT_SERIES_STARFILE, fn_star, t);

            }
        }
    }
    else
    {
        // The new way of reading in the tomogram STAR files
        for (int t = 0; t < tc; t++)
        {
            FileName name = globalTable.getString(EMDL_TOMO_NAME, t);
            std::string fn_star = globalTable.getString(EMDL_TOMO_TILT_SERIES_STARFILE, t);
            tomogramTables[t].read(fn_star, name);

            if ((tomogramTables[t]).numberOfObjects() == 0)
            {
                REPORT_ERROR("ERROR: could not read data from " + fn_star + ". Does the table have the correct name: data_" + name + "?");
            }

            if (!(tomogramTables[t]).containsLabel(EMDL_MICROGRAPH_PRE_EXPOSURE))
            {
                REPORT_ERROR("ERROR: tomogramTable " + fn_star + " does not contain compulsory rlnMicrographPreExposure label");
            }
        }
    }
    return true;

}

void TomogramSet::write(FileName filename)
{
    FileName fn_outdir = filename.beforeLastOf("/") + "/";

    const int tc = tomogramTables.size();

    // Change all the filenames in tomograms.star
    for (int t = 0; t < tc; t++)
    {
        FileName fn_star;
        globalTable.getValue(EMDL_TOMO_TILT_SERIES_STARFILE, fn_star, t);
        FileName fn_newstar = getOutputFileWithNewUniqueDate(fn_star, fn_outdir);
        globalTable.setValue(EMDL_TOMO_TILT_SERIES_STARFILE, fn_newstar, t);

        // Create output directory if necessary
        FileName newdir = fn_newstar.beforeLastOf("/");
        if (!exists(newdir)) mktree(newdir);

        // Write the individual tomogram starfile
        if (tomogramTables[t].containsLabel(EMDL_TOMO_PROJECTION_X)) tomogramTables[t].deactivateLabel(EMDL_TOMO_PROJECTION_X);
        if (tomogramTables[t].containsLabel(EMDL_TOMO_PROJECTION_Y)) tomogramTables[t].deactivateLabel(EMDL_TOMO_PROJECTION_Y);
        if (tomogramTables[t].containsLabel(EMDL_TOMO_PROJECTION_Z)) tomogramTables[t].deactivateLabel(EMDL_TOMO_PROJECTION_Z);
        if (tomogramTables[t].containsLabel(EMDL_TOMO_PROJECTION_W)) tomogramTables[t].deactivateLabel(EMDL_TOMO_PROJECTION_W);

        tomogramTables[t].write(fn_newstar);
    }

    // Also write the (now modified with fn_newstars) tilt_series.star file in the root directory
    globalTable.write(filename);

}

Tomogram TomogramSet::loadTomogram(int index, bool loadImageData, bool loadEvenFramesOnly, bool loadOddFramesOnly,
                                   int _w0, int _h0, int _d0) const //Set loadEven/OddFramesOnly to True to loadImageData from rlnTomoMicrographNameEven/Odd rather than rlnMicrographName
{
	Tomogram out;

    std::string tomoName;
    i3Vector stackSize;

    globalTable.getValueSafely(EMDL_TOMO_NAME, tomoName, index);
    const MetaDataTable& m = tomogramTables[index];

    if (_w0 != -999 && _h0 != -999 && _d0 != -999)
    {
        out.w0 = _w0;
        out.h0 = _h0;
        out.d0 = _d0;
    }
    else if (globalTable.containsLabel(EMDL_TOMO_SIZE_X) &&
        globalTable.containsLabel(EMDL_TOMO_SIZE_Y) &&
        globalTable.containsLabel(EMDL_TOMO_SIZE_Z))
    {
        globalTable.getValueSafely(EMDL_TOMO_SIZE_X, out.w0, index);
        globalTable.getValueSafely(EMDL_TOMO_SIZE_Y, out.h0, index);
        globalTable.getValueSafely(EMDL_TOMO_SIZE_Z, out.d0, index);
    }
    else
    {
        out.w0 = -999;
        out.h0 = -999;
        out.d0 = -999;
    }

    // Select only a subset of the tilt series images with the lowest dose
    out.frameCount = tomogramTables[index].numberOfObjects();

    if (globalTable.containsLabel(EMDL_CTF_BFACTOR_PERELECTRONDOSE))
        globalTable.getValue(EMDL_CTF_BFACTOR_PERELECTRONDOSE, out.BfactorPerElectronDose, index);
    else
        out.BfactorPerElectronDose = 0.;

    if (globalTable.containsLabel(EMDL_TOMO_TILT_SERIES_NAME))
    {
        // option A: Kino's original IMOD import functionality
        globalTable.getValueSafely(EMDL_TOMO_TILT_SERIES_NAME, out.tiltSeriesFilename, index);
        globalTable.getValueSafely(EMDL_TOMO_FRAME_COUNT, out.frameCount, index);

        if (loadImageData)
        {
            out.stack.read(out.tiltSeriesFilename);

            out.hasImage = true;

            stackSize.x = out.stack.xdim;
            stackSize.y = out.stack.ydim;
            stackSize.z = out.stack.zdim;
       }
        else
        {
            out.hasImage = false;

            t3Vector<long int> isl = ImageFileHelper::getSize(out.tiltSeriesFilename);

            stackSize.x = isl.x;
            stackSize.y = isl.y;
            stackSize.z = isl.z;

            // Still set image size in header, as this is used for example in TomoBackprojectProgram::getProjectMatrices
            out.stack.xdim = stackSize.x;
            out.stack.ydim = stackSize.y;
            out.stack.zdim = stackSize.z;
        }

    }
    else
    {
        // option B: new functionality to work directly with images from RELION's motioncorr runner

        out.tiltSeriesFilename = "";
        out.frameCount = m.numberOfObjects();

        // Get image size from the first image in the tomogramTable
        if (!m.containsLabel(EMDL_MICROGRAPH_NAME) && !loadEvenFramesOnly && !loadOddFramesOnly)
        {
            REPORT_ERROR("ERROR: tomogramTable for " + tomoName + " does not contain a rlnMicrographName label, yet the globalTable also does not contain a rlnTomoTiltSeriesName label!");
        }
	
        if (!m.containsLabel(EMDL_MICROGRAPH_ODD) && loadOddFramesOnly)
        {
            REPORT_ERROR("ERROR: tomogramTable for " + tomoName + " does not contain a rlnTomoMicrographNameOdd label");
        }

        if (!m.containsLabel(EMDL_MICROGRAPH_EVEN) && loadEvenFramesOnly)
        {
            REPORT_ERROR("ERROR: tomogramTable for " + tomoName + " does not contain a rlnTomoMicrographNameEven label");
        }
	
        std::string fn_img;
        if (loadEvenFramesOnly)
        {
            m.getValueSafely(EMDL_MICROGRAPH_EVEN, fn_img, 0);
        }
        else if (loadOddFramesOnly)
        {
            m.getValueSafely(EMDL_MICROGRAPH_ODD, fn_img, 0);
        }
        else
        {
        	m.getValueSafely(EMDL_MICROGRAPH_NAME, fn_img, 0);
        }

        Image<RFLOAT> I;
        I.read(fn_img,false); // false means don't read the actual image data, only the header

        stackSize.x = XSIZE(I());
        stackSize.y = YSIZE(I());
        stackSize.z = out.frameCount;

        if (loadImageData)
        {
            Image<RFLOAT> myStack(stackSize.x, stackSize.y, stackSize.z);
            for (int f = 0; f < out.frameCount; f++)
            {

                if (loadEvenFramesOnly)
                {
                    m.getValueSafely(EMDL_MICROGRAPH_EVEN, fn_img, f);
                }
                else if (loadOddFramesOnly)
                {
                    m.getValueSafely(EMDL_MICROGRAPH_ODD, fn_img, f);
                }
                else
                {
                    m.getValueSafely(EMDL_MICROGRAPH_NAME, fn_img, f);
                }

                Image<RFLOAT> I2;
                I2.read(fn_img);

                if (XSIZE(I2()) != stackSize.x || YSIZE(I2()) != stackSize.y)
                {
                    REPORT_ERROR("ERROR: unequal image dimensions in the individual tilt series images of tomogram: " + tomoName);
                }
                myStack().setSlice(f, I2());

            }

            out.stack.copyDataAndSizeFrom(myStack);
            out.hasImage = true;
        }
        else
        {
            // Still set image size in header, as this is used for example in TomoBackprojectProgram::getProjectMatrices
            out.stack.xdim = stackSize.x;
            out.stack.ydim = stackSize.y;
            out.stack.zdim = stackSize.z;
            out.hasImage = false;
        }

    }

	out.imageSize = stackSize.xy();
    out.centre = d3Vector(out.w0/2.0, out.h0/2.0, out.d0/2.0);

	globalTable.getValueSafely(EMDL_TOMO_HANDEDNESS, out.handedness, index);

    // Now that we do notioncorrection on tomogram_sets, the micrograph pixel size may not have been set yet...
    if (globalTable.containsLabel(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE))
        globalTable.getValue(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, out.optics.pixelSize, index);
    else
        out.optics.pixelSize = -999.;

    globalTable.getValueSafely(EMDL_CTF_VOLTAGE, out.optics.voltage, index);
    globalTable.getValueSafely(EMDL_CTF_CS, out.optics.Cs, index);
    globalTable.getValueSafely(EMDL_CTF_Q0, out.optics.Q0, index);

	out.hasDeformations = (globalTable.containsLabel(EMDL_TOMO_DEFORMATION_GRID_SIZE_X) &&
                           globalTable.containsLabel(EMDL_TOMO_DEFORMATION_GRID_SIZE_Y) );

	i2Vector deformationGridSize;
	std::string deformationType = "";

	if (out.hasDeformations)
	{
		deformationGridSize.x = globalTable.getInt(EMDL_TOMO_DEFORMATION_GRID_SIZE_X, index);
		deformationGridSize.y = globalTable.getInt(EMDL_TOMO_DEFORMATION_GRID_SIZE_Y, index);
		deformationType = globalTable.getString(EMDL_TOMO_DEFORMATION_TYPE, index);

		out.imageDeformations.resize(out.frameCount);

		if (deformationType != "spline" && deformationType != "Fourier" && deformationType != "linear")
		{
			REPORT_ERROR_STR(
				"TomogramSet::loadTomogram: illegal deformation type '"
				<< deformationType << "'");
		}
	}

	out.cumulativeDose.resize(out.frameCount);
	out.centralCTFs.resize(out.frameCount);
	out.projectionMatrices.resize(out.frameCount);
    out.hasMatrices = (m.containsLabel(EMDL_TOMO_YTILT) &&
                       m.containsLabel(EMDL_TOMO_ZROT) &&
                       m.containsLabel(EMDL_TOMO_XSHIFT_ANGST) &&
                       m.containsLabel(EMDL_TOMO_YSHIFT_ANGST));
    bool has_old_matrix = (m.containsLabel(EMDL_TOMO_PROJECTION_X) &&
                       m.containsLabel(EMDL_TOMO_PROJECTION_Y) &&
                       m.containsLabel(EMDL_TOMO_PROJECTION_Z) &&
                       m.containsLabel(EMDL_TOMO_PROJECTION_W));
    if (has_old_matrix)
    {
        std::cerr << " WARNING: tomogram " << tomoName << " has relion-4 definition of projection matrices; converting them now... " << std::endl;
        out.hasMatrices = true;
    }

	for (int f = 0; f < out.frameCount; f++)
	{

        if (has_old_matrix)
        {
            d4Matrix& P = out.projectionMatrices[f];
            std::vector<EMDLabel> rows({
                                               EMDL_TOMO_PROJECTION_X,
                                               EMDL_TOMO_PROJECTION_Y,
                                               EMDL_TOMO_PROJECTION_Z,
                                               EMDL_TOMO_PROJECTION_W });

            for (int i = 0; i < 4; i++)
            {
                std::vector<double> vals;
                m.getValueSafely(rows[i], vals, f);

                for (int j = 0; j < 4; j++)
                {
                    P(i,j) = vals[j];
                }
            }

        }
        else if (out.hasMatrices)
        {

            RFLOAT xtilt, ytilt, zrot, xshift_angst, yshift_angst;

            if (m.containsLabel(EMDL_TOMO_XTILT))
                m.getValue(EMDL_TOMO_XTILT, xtilt, f);
            else
                xtilt = 0.;

            m.getValueSafely(EMDL_TOMO_YTILT, ytilt, f);
            m.getValueSafely(EMDL_TOMO_ZROT, zrot, f);
            m.getValueSafely(EMDL_TOMO_XSHIFT_ANGST, xshift_angst, f);
            m.getValueSafely(EMDL_TOMO_YSHIFT_ANGST, yshift_angst, f);

            out.setProjectionMatrix(f, xtilt, ytilt, zrot, xshift_angst, yshift_angst);

        }

		CTF& ctf = out.centralCTFs[f];

		m.getValueSafely(EMDL_CTF_DEFOCUSU, ctf.DeltafU, f);
		m.getValueSafely(EMDL_CTF_DEFOCUSV, ctf.DeltafV, f);
		m.getValueSafely(EMDL_CTF_DEFOCUS_ANGLE, ctf.azimuthal_angle, f);

		ctf.Q0 = out.optics.Q0;
		ctf.Cs = out.optics.Cs;
		ctf.kV = out.optics.voltage;

		if (m.containsLabel(EMDL_CTF_SCALEFACTOR))
		{
			ctf.scale = m.getDouble(EMDL_CTF_SCALEFACTOR, f);
		}
        else
        {
            ctf.scale = 1.;
        }

        if (m.containsLabel(EMDL_CTF_BFACTOR))
        {
            ctf.Bfac = m.getDouble(EMDL_CTF_BFACTOR, f);
        }
        else
        {
            ctf.Bfac = 0.;
        }

        if (m.containsLabel(EMDL_CTF_PHASESHIFT))
        {
            ctf.phase_shift = m.getDouble(EMDL_CTF_PHASESHIFT, f);
        }
        else
        {
            ctf.phase_shift = 0.;
        }

		ctf.initialise();

		m.getValueSafely(EMDL_MICROGRAPH_PRE_EXPOSURE, out.cumulativeDose[f], f);

		if (out.hasDeformations && m.containsLabel(EMDL_TOMO_DEFORMATION_COEFFICIENTS))
		{
			const std::vector<double> coeffs = m.getDoubleVector(
					EMDL_TOMO_DEFORMATION_COEFFICIENTS, f);

			if (deformationType == "spline")
			{
				out.imageDeformations[f] = std::make_shared<Spline2DDeformation>(
							stackSize.xy(), deformationGridSize, &coeffs[0]);
			}
			else if (deformationType == "Fourier")
			{
				out.imageDeformations[f] = std::make_shared<Fourier2DDeformation>(
							stackSize.xy(), deformationGridSize, &coeffs[0]);
			}
			else if (deformationType == "linear")
			{
				out.imageDeformations[f] = std::make_shared<Linear2DDeformation>(
							stackSize.xy(), &coeffs[0]);
			}
		}
	}

	out.frameSequence = IndexSort<double>::sortIndices(out.cumulativeDose);

	if (globalTable.containsLabel(EMDL_TOMO_IMPORT_FRACT_DOSE))
	{
		out.fractionalDose = globalTable.getDouble(EMDL_TOMO_IMPORT_FRACT_DOSE, index);
	}
	else
	{
		out.fractionalDose = out.cumulativeDose[out.frameSequence[1]] - out.cumulativeDose[out.frameSequence[0]];
	}

	out.name = tomoName;

	if (globalTable.containsLabel(EMDL_TOMO_FIDUCIALS_STARFILE))
	{
		 globalTable.getValue(EMDL_TOMO_FIDUCIALS_STARFILE, out.fiducialsFilename, index);
	}
	else
	{
		out.fiducialsFilename = "";
	}

	if (globalTable.containsLabel(EMDL_TOMO_DEFOCUS_SLOPE))
	{
		 globalTable.getValue(EMDL_TOMO_DEFOCUS_SLOPE, out.defocusSlope, index);
	}
	else
	{
		out.defocusSlope = 1.0;
	}

	return out;
}

int TomogramSet::size() const
{
	return tomogramTables.size();
}

void TomogramSet::setProjectionAngles(int tomogramIndex, int frame, RFLOAT xtilt, RFLOAT ytilt, RFLOAT zrot, RFLOAT xshift_angst, RFLOAT yshift_angst)
{
	MetaDataTable& m = tomogramTables[tomogramIndex];

    m.setValue(EMDL_TOMO_XTILT, xtilt, frame);
    m.setValue(EMDL_TOMO_YTILT, ytilt, frame);
    m.setValue(EMDL_TOMO_ZROT,  zrot, frame);
    m.setValue(EMDL_TOMO_XSHIFT_ANGST, xshift_angst, frame);
    m.setValue(EMDL_TOMO_YSHIFT_ANGST, yshift_angst, frame);
}

void TomogramSet::setCtf(int tomogramIndex, int frame, const CTF& ctf)
{
	MetaDataTable& m = tomogramTables[tomogramIndex];
	
	m.setValue(EMDL_CTF_DEFOCUSU, ctf.DeltafU, frame);
	m.setValue(EMDL_CTF_DEFOCUSV, ctf.DeltafV, frame);
	m.setValue(EMDL_CTF_DEFOCUS_ANGLE, ctf.azimuthal_angle, frame);
	m.setValue(EMDL_CTF_SCALEFACTOR, ctf.scale, frame);
}

void TomogramSet::setDose(int tomogramIndex, int frame, double dose)
{
	MetaDataTable& m = tomogramTables[tomogramIndex];
	
	m.setValue(EMDL_MICROGRAPH_PRE_EXPOSURE, dose, frame);
}

void TomogramSet::setTiltSeriesFile(int tomogramIndex, const std::string &filename)
{
	globalTable.setValue(EMDL_TOMO_TILT_SERIES_NAME, filename, tomogramIndex);
}

void TomogramSet::setFiducialsFile(int tomogramIndex, const std::string &filename)
{
	globalTable.setValue(EMDL_TOMO_FIDUCIALS_STARFILE, filename, tomogramIndex);
}

void TomogramSet::setDefocusSlope(int tomogramIndex, double slope)
{
	globalTable.setValue(EMDL_TOMO_DEFOCUS_SLOPE, slope, tomogramIndex);
}
void TomogramSet::applyTiltAngleOffset(int tomogramIndex, double offset)
{
    MetaDataTable& m = tomogramTables[tomogramIndex];

    FOR_ALL_OBJECTS_IN_METADATA_TABLE(m) {
        double tilt;
        m.getValueSafely(EMDL_TOMO_YTILT, tilt);
        tilt += offset;
        m.setValue(EMDL_TOMO_YTILT, tilt);
    }
}

void TomogramSet::setDeformation(
	int tomogramIndex,
	gravis::i2Vector gridSize,
	const std::string& deformationType,
	const std::vector<std::vector<double>>& coeffs)
{
	globalTable.setValue(EMDL_TOMO_DEFORMATION_GRID_SIZE_X, gridSize.x, tomogramIndex);
	globalTable.setValue(EMDL_TOMO_DEFORMATION_GRID_SIZE_Y, gridSize.y, tomogramIndex);
	globalTable.setValue(EMDL_TOMO_DEFORMATION_TYPE, deformationType, tomogramIndex);

	MetaDataTable& mdt = tomogramTables[tomogramIndex];

	const int fc = coeffs.size();

	for (int f = 0; f < fc; f++)
	{
		mdt.setValue(EMDL_TOMO_DEFORMATION_COEFFICIENTS, coeffs[f], f);
	}
}

void TomogramSet::clearDeformation()
{
	globalTable.deactivateLabel(EMDL_TOMO_DEFORMATION_GRID_SIZE_X);
	globalTable.deactivateLabel(EMDL_TOMO_DEFORMATION_GRID_SIZE_Y);

	for (int t = 0; t < tomogramTables.size(); t++)
	{
		tomogramTables[t].deactivateLabel(EMDL_TOMO_DEFORMATION_COEFFICIENTS);
	}
}

int TomogramSet::getTomogramIndex(std::string tomogramName) const
{
	const int tc = globalTable.numberOfObjects();

	for (int t = 0; t < tc; t++)
	{
		std::string name_t;
		globalTable.getValueSafely(EMDL_TOMO_NAME, name_t, t);

		if (name_t == tomogramName)
		{
			return t;
		}
	}

	return -1;
}

std::string TomogramSet::getTomogramName(int index) const
{
	std::string name;
	globalTable.getValueSafely(EMDL_TOMO_NAME, name, index);

	return name;
}


int TomogramSet::getTomogramIndexSafely(std::string tomogramName) const
{
	int t = getTomogramIndex(tomogramName);

	if (t < 0)
	{
		REPORT_ERROR_STR("No tomogram named '" << tomogramName << "' found in the set");
	}
	else
	{
		return t;
	}
}

int TomogramSet::getFrameCount(int index) const
{
	return tomogramTables[index].numberOfObjects();
}

int TomogramSet::getMaxFrameCount() const
{
	int max_val = 0;

	for (int t = 0; t < tomogramTables.size(); t++)
	{
		const int fc = tomogramTables[t].numberOfObjects();

		if (fc > max_val)
		{
			max_val = fc;
		}
	}

	return max_val;
}

double TomogramSet::getOriginalPixelSize(int index) const
{
    if (!globalTable.containsLabel(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE))
        REPORT_ERROR("ERROR: cannot find rlnMicrographOriginalPixelSize label in star file.");
	return globalTable.getDouble(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, index);
}

double TomogramSet::getTiltSeriesPixelSize(int index) const
{
	if (!globalTable.containsLabel(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE))
        REPORT_ERROR("ERROR: cannot find rlnTomoTiltSeriesPixelSize label in star file.");
    return globalTable.getDouble(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, index);
}

std::string TomogramSet::getOpticsGroupName(int index) const
{
	if (!globalTable.containsLabel(EMDL_IMAGE_OPTICS_GROUP_NAME))
	{
		return "opticsGroup1";
	}
	else
	{
		return globalTable.getString(EMDL_IMAGE_OPTICS_GROUP_NAME, index);
	}
}

void TomogramSet::generateSingleMetaDataTable(MetaDataTable &MDout, ObservationModel &obsModel)
{
    MDout.clear();
    for (long int t = 0; t < tomogramTables.size(); t++)
    {
        // Store all the necessary optics stuff in an opticsGroup per tomogram
        RFLOAT moviePixelSize, voltage, Cs, Q0;
        std::string tomo_name = getTomogramName(t);
        globalTable.getValueSafely(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, moviePixelSize, t);
        globalTable.getValueSafely(EMDL_CTF_VOLTAGE, voltage, t);
        globalTable.getValueSafely(EMDL_CTF_CS, Cs, t);
        globalTable.getValueSafely(EMDL_CTF_Q0, Q0, t);
        obsModel.opticsMdt.addObject();
        obsModel.opticsMdt.setValue(EMDL_IMAGE_OPTICS_GROUP_NAME, tomo_name);
        obsModel.opticsMdt.setValue(EMDL_IMAGE_OPTICS_GROUP, t+1);
        obsModel.opticsMdt.setValue(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, moviePixelSize);
        obsModel.opticsMdt.setValue(EMDL_CTF_VOLTAGE, voltage);
        obsModel.opticsMdt.setValue(EMDL_CTF_CS, Cs);
        obsModel.opticsMdt.setValue(EMDL_CTF_Q0, Q0);

        if (globalTable.containsLabel(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE))
        {
            RFLOAT pixSize;
            globalTable.getValueSafely(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, pixSize, t);
            obsModel.opticsMdt.setValue(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, pixSize);
        }

        FOR_ALL_OBJECTS_IN_METADATA_TABLE(tomogramTables[t])
        {
            tomogramTables[t].setValue(EMDL_IMAGE_OPTICS_GROUP, t+1);
        }

        MDout.append(tomogramTables[t]);
    }
}

void TomogramSet::convertBackFromSingleMetaDataTable(MetaDataTable &MDin)
{
    if (!MDin.containsLabel(EMDL_MICROGRAPH_PRE_EXPOSURE))
    {
        REPORT_ERROR("BUG: MDin should contain a rlnMicrographPreExposure label");
    }

    for (long int t = 0; t < tomogramTables.size(); t++)
    {
        MetaDataTable MDsub = subsetMetaDataTable(MDin, EMDL_IMAGE_OPTICS_GROUP, t+1, t+1);

        // Make sure no one unsorted the tilt images in each serie...
        MDsub.newSort(EMDL_MICROGRAPH_PRE_EXPOSURE);

        if (MDsub.numberOfObjects() != tomogramTables[t].numberOfObjects())
        {
            REPORT_ERROR("ERROR: unequal number of Objects in tiltserie starfiles");
        }

        FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsub)
        {
            tomogramTables[t].setObject(MDsub.getObject(), current_object);
        }

       tomogramTables[t].deactivateLabel(EMDL_IMAGE_OPTICS_GROUP);

    }
}
