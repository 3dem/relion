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

TomogramSet::TomogramSet(std::string filename, bool verbose)
{
    if (!read(filename, verbose))
    {
        // This may be a tomograms.star file in the old, original relion-4.0 format. Try to convert
        // TODO: make conversion from old tomograms.star to new tilt_series.star!

        std::ifstream ifs(filename);
        if (!ifs)
        {
            REPORT_ERROR_STR("TomogramSet::TomogramSet: Unable to read " << filename);
        }
        else
        {
            globalTable.readStar(ifs, "global");

            if (!globalTable.containsLabel(EMDL_TOMO_NAME))
            {
                REPORT_ERROR("ERROR: input starfile for TomogramSet " + filename + " does not contain rlnTomoName label ");
            }

            // remove information from optics groups and move into globalObsModel
            ObservationModel globalObsModel(globalTable);
            globalTable.deactivateLabel(EMDL_IMAGE_OPTICS_GROUP_NAME);
            globalTable.deactivateLabel(EMDL_IMAGE_OPTICS_GROUP);
            globalTable.deactivateLabel(EMDL_CTF_CS);
            globalTable.deactivateLabel(EMDL_CTF_VOLTAGE);
            globalTable.deactivateLabel(EMDL_CTF_Q0);

            const int tc = globalTable.numberOfObjects();

            tomogramTables.resize(tc);
            tomogramObsModels.resize(tc);
            tomogramNames.resize(tc);

            std::vector<MetaDataTable> allTables = MetaDataTable::readAll(ifs, tc+1);

            for (int t = 0; t < tc; t++)
            {
                const std::string expectedNewName = globalTable.getString(EMDL_TOMO_NAME, t);
                const std::string name = allTables[t+1].getName();

                if (name != expectedNewName)
                {
                    REPORT_ERROR_STR("TomogramSet::TomogramSet: file is corrupted " << filename);
                }

                tomogramTables[t] = allTables[t+1];
                tomogramTables[t].setName("tilt_images");
                tomogramObsModels[t] = globalObsModel;
                tomogramNames[t] = expectedNewName;
            }
        }

    }

    globalTable.setName("global");

}

bool TomogramSet::read(std::string filename, bool verbose)
{

    globalTable.read(filename, "global");

    const int tc = globalTable.numberOfObjects();

    if (tc == 0) return false;

    if (!globalTable.containsLabel(EMDL_TOMO_TILT_SERIES_STARFILE))
    {
        std::cerr << "Warning: " << filename
                  << " does not have rlnTomoTiltSeriesStarFile labels. It may be written in an old format. If so, will try to convert ..."
                  << std::endl;
        return false;
    }

    if (!globalTable.containsLabel(EMDL_TOMO_NAME))
    {
        REPORT_ERROR("ERROR: input starfile for TomogramSet " + filename + " does not contain rlnTomoName label ");
    }

    tomogramTables.resize(tc);
    tomogramObsModels.resize(tc);
    tomogramNames.resize(tc);

    for (int t = 0; t < tc; t++)
    {
        tomogramNames[t] = globalTable.getString(EMDL_TOMO_NAME, t);
        std::string fn_star = globalTable.getString(EMDL_TOMO_TILT_SERIES_STARFILE, t);
        ObservationModel::loadSafely(fn_star, tomogramObsModels[t], tomogramTables[t], "tilt_images", verbose);
        // Make sure tilt images are sorted on their index (used to convert back from single large metadatatable)
        tomogramTables[t].newSort(EMDL_TOMO_TILT_MOVIE_INDEX);
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
        tomogramObsModels[t].save(tomogramTables[t], fn_newstar, "tilt_images");
    }

    // Also write the (now modified with fn_newstars) tilt_series.star file in the root directory
    globalTable.write(filename);

}

Tomogram TomogramSet::loadTomogram(int index, bool loadImageData) const
{
	Tomogram out;

	std::string tomoName, stackFn;

	globalTable.getValueSafely(EMDL_TOMO_NAME, tomoName, index);
	globalTable.getValueSafely(EMDL_TOMO_TILT_SERIES_NAME, stackFn, index);
	globalTable.getValueSafely(EMDL_TOMO_FRAME_COUNT, out.frameCount, index);

	i3Vector stackSize;

	if (loadImageData)
	{
		out.stack.read(stackFn);
		out.hasImage = true;

		stackSize.x = out.stack.xdim;
		stackSize.y = out.stack.ydim;
		stackSize.z = out.stack.zdim;
	}
	else
	{
		out.hasImage = false;

		t3Vector<long int> isl = ImageFileHelper::getSize(stackFn);

		stackSize.x = isl.x;
		stackSize.y = isl.y;
		stackSize.z = isl.z;
	}

	out.imageSize = stackSize.xy();

	out.tiltSeriesFilename = stackFn;

	globalTable.getValueSafely(EMDL_TOMO_SIZE_X, out.w0, index);
	globalTable.getValueSafely(EMDL_TOMO_SIZE_Y, out.h0, index);
	globalTable.getValueSafely(EMDL_TOMO_SIZE_Z, out.d0, index);

	out.centre = d3Vector(out.w0/2.0, out.h0/2.0, out.d0/2.0);

	globalTable.getValueSafely(EMDL_TOMO_HANDEDNESS, out.handedness, index);

	double Q0;

	globalTable.getValueSafely(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, out.optics.pixelSize, index);
	globalTable.getValueSafely(EMDL_CTF_VOLTAGE, out.optics.voltage, index);
	globalTable.getValueSafely(EMDL_CTF_CS, out.optics.Cs, index);
	globalTable.getValueSafely(EMDL_CTF_Q0, Q0, index);

	out.hasDeformations = (
		globalTable.containsLabel(EMDL_TOMO_DEFORMATION_GRID_SIZE_X) &&
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

	const MetaDataTable& m = tomogramTables[index];

	out.cumulativeDose.resize(out.frameCount);
	out.centralCTFs.resize(out.frameCount);
	out.projectionMatrices.resize(out.frameCount);

	for (int f = 0; f < out.frameCount; f++)
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

		CTF& ctf = out.centralCTFs[f];

		m.getValueSafely(EMDL_CTF_DEFOCUSU, ctf.DeltafU, f);
		m.getValueSafely(EMDL_CTF_DEFOCUSV, ctf.DeltafV, f);
		m.getValueSafely(EMDL_CTF_DEFOCUS_ANGLE, ctf.azimuthal_angle, f);

		ctf.Q0 = Q0;
		ctf.Cs = out.optics.Cs;
		ctf.kV = out.optics.voltage;

		if (m.containsLabel(EMDL_CTF_SCALEFACTOR))
		{
			ctf.scale = m.getDouble(EMDL_CTF_SCALEFACTOR, f);
		}

		ctf.initialise();

		m.getValueSafely(EMDL_MICROGRAPH_PRE_EXPOSURE, out.cumulativeDose[f], f);

		if (out.hasDeformations &&
			m.containsLabel(EMDL_TOMO_DEFORMATION_COEFFICIENTS))
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

	if (globalTable.containsLabel(EMDL_IMAGE_OPTICS_GROUP_NAME))
	{
		out.opticsGroupName = globalTable.getString(EMDL_IMAGE_OPTICS_GROUP_NAME, index);
	}
	else
	{
		out.opticsGroupName = "opticsGroup1";
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

void TomogramSet::addTomogram(
		std::string tomoName, std::string stackFilename,
		const std::vector<gravis::d4Matrix>& projections, 
		int w, int h, int d, 
		const std::vector<double>& dose,
		double fractionalDose,
		const std::vector<CTF>& ctfs, 
		double handedness, 
		double pixelSize,
		const std::string& opticsGroupName)
{
	const int index = globalTable.numberOfObjects();
	const int fc = projections.size();
	
	globalTable.addObject();
	
	globalTable.setValue(EMDL_TOMO_NAME, tomoName, index);
	globalTable.setValue(EMDL_TOMO_TILT_SERIES_NAME, stackFilename, index);
	globalTable.setValue(EMDL_TOMO_FRAME_COUNT, fc, index);
	
	globalTable.setValue(EMDL_TOMO_SIZE_X, w, index);
	globalTable.setValue(EMDL_TOMO_SIZE_Y, h, index);
	globalTable.setValue(EMDL_TOMO_SIZE_Z, d, index);	
	globalTable.setValue(EMDL_TOMO_HANDEDNESS, handedness, index);
	
	const CTF& ctf0 = ctfs[0];
	
	globalTable.setValue(EMDL_IMAGE_OPTICS_GROUP_NAME, opticsGroupName, index);
	globalTable.setValue(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, pixelSize, index);
	globalTable.setValue(EMDL_CTF_VOLTAGE, ctf0.kV, index);
	globalTable.setValue(EMDL_CTF_CS, ctf0.Cs, index);
	globalTable.setValue(EMDL_CTF_Q0, ctf0.Q0, index);
	globalTable.setValue(EMDL_TOMO_IMPORT_FRACT_DOSE, fractionalDose, index);
	
	if (tomogramTables.size() != index)
	{
		REPORT_ERROR_STR("TomogramSet::add: corrupted tomogram set: tomogramTables.size() = "
						 << tomogramTables.size() << ", globalTable.numberOfObjects() = " 
						 << globalTable.numberOfObjects());
	}
	
	tomogramTables.push_back(MetaDataTable());
	MetaDataTable& m = tomogramTables[index];
	m.setName(tomoName);
		
	for (int f = 0; f < fc; f++)
	{
		m.addObject();
		
		setProjection(index, f, projections[f]);
		setCtf(index, f, ctfs[f]);	
		setDose(index, f, dose[f]);
	}
}

int TomogramSet::size() const
{
	return tomogramTables.size();
}

void TomogramSet::setProjections(int tomogramIndex, const std::vector<d4Matrix>& proj)
{
	MetaDataTable& m = tomogramTables[tomogramIndex];
	
	const int fc = proj.size();
	
	if (fc != m.numberOfObjects())
	{
		REPORT_ERROR_STR("TomogramSet::setProjections: frame count mismatch");
	}
	
	for (int f = 0; f < fc; f++)
	{
		setProjection(tomogramIndex, f, proj[f]);
	}
}

void TomogramSet::setProjection(int tomogramIndex, int frame, const d4Matrix& P)
{
	MetaDataTable& m = tomogramTables[tomogramIndex];
	
	m.setValue(EMDL_TOMO_PROJECTION_X, std::vector<double>{P(0,0), P(0,1), P(0,2), P(0,3)}, frame);
	m.setValue(EMDL_TOMO_PROJECTION_Y, std::vector<double>{P(1,0), P(1,1), P(1,2), P(1,3)}, frame);
	m.setValue(EMDL_TOMO_PROJECTION_Z, std::vector<double>{P(2,0), P(2,1), P(2,2), P(2,3)}, frame);
	m.setValue(EMDL_TOMO_PROJECTION_W, std::vector<double>{P(3,0), P(3,1), P(3,2), P(3,3)}, frame);
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

double TomogramSet::getPixelSize(int index) const
{
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
    std::vector<FileName> fn_stars;
    for (long int t = 0; t < globalTable.numberOfObjects(); t++)
    {
        FileName fn_star = globalTable.getString(EMDL_TOMO_TILT_SERIES_STARFILE, t);
        fn_stars.push_back(fn_star);
    }

    combineStarfiles(fn_stars, MDout, obsModel);

}

void TomogramSet::convertBackFromSingleMetaDataTable(MetaDataTable &MDin, ObservationModel &obsModel)
{
    if (!MDin.containsLabel(EMDL_TOMO_TILT_MOVIE_INDEX))
    {
        REPORT_ERROR("BUG: the MDin that is passed to TomogramSet::convertBackFromSingleMetaDataTable should contain a rlnTomoTiltMovieIndex label");
    }

    for (long int t = 0; t < globalTable.numberOfObjects(); t++)
    {
        tomogramObsModels[t] = obsModel;
        MetaDataTable MDjoin = subsetMetaDataTable(MDin, EMDL_TOMO_NAME, tomogramNames[t], false);

        MDjoin.newSort(EMDL_TOMO_TILT_MOVIE_INDEX);

        if (MDjoin.numberOfObjects() != tomogramTables[t].numberOfObjects())
        {
            REPORT_ERROR("ERROR: unequal number of Objects in tiltserie starfiles");
        }

        FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDjoin)
        {
            tomogramTables[t].setObject(MDjoin.getObject(), current_object);
        }

        tomogramObsModels[t].removeUnusedOpticsGroups(tomogramTables[t]);

    }
}
