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
	std::ifstream ifs(filename);

	bool namesAreOld = false;

	if (!ifs)
	{
		REPORT_ERROR_STR("TomogramSet::TomogramSet: Unable to read " << filename);
	}
	else
	{
		globalTable.readStar(ifs, "global");
		
		const int tc = globalTable.numberOfObjects();
		
		tomogramTables.resize(tc);
		
		std::vector<MetaDataTable> allTables = MetaDataTable::readAll(ifs, tc+1);
		
		for (int t = 0; t < tc; t++)
		{
			const std::string expectedOldName = "tomo_" + ZIO::itoa(t);
			const std::string expectedNewName = globalTable.getString(EMDL_TOMO_NAME, t);
			const std::string name = allTables[t+1].getName();

			if (name == expectedOldName)
			{
				namesAreOld = true;
			}
			else if (name != expectedNewName)
			{
				REPORT_ERROR_STR("TomogramSet::TomogramSet: file is corrupted " << filename);
			}

			tomogramTables[t] = allTables[t+1];
			tomogramTables[t].setName(expectedNewName);
		}	
	}
	
	globalTable.setName("global");

	if (verbose && namesAreOld)
	{
		Log::warn("Tomogram set " + filename + " is out of date. You are recommended to run relion_exp_update_tomogram_set on it.");
	}
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

void TomogramSet::write(std::string filename) const
{
	const int tc = tomogramTables.size();

	if (filename.find_last_of('/') != std::string::npos)
	{
		std::string path = filename.substr(0, filename.find_last_of('/'));
		mktree(path);
	}

	std::ofstream ofs(filename);
	
	if (!ofs)
	{
		REPORT_ERROR("TomogramSet::write: unable to write to "+filename);
	}
	
	globalTable.write(ofs);

	for (int t = 0; t < tc; t++)
	{
		tomogramTables[t].write(ofs);
	}
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
