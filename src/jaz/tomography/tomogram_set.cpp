#include "tomogram_set.h"
#include <fstream>
#include <sstream>
#include <src/error.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/util/zio.h>

using namespace gravis;


TomogramSet::TomogramSet()
{
	globalTable.setName("global");
}

TomogramSet::TomogramSet(std::string filename)
{
	std::ifstream ifs(filename);

	if (!ifs.fail())
	{
		globalTable.readStar(ifs, "global");
		
		const int tc = globalTable.numberOfObjects();
		
		tomogramTables.resize(tc);
		
		std::vector<MetaDataTable> allTables = MetaDataTable::readAll(ifs, tc+1);
		
		for (int t = 0; t < tc; t++)
		{
			tomogramTables[t] = allTables[t+1];
			tomogramTables[t].setName("tomo_" + ZIO::itoa(t));
		}	
	}
	else
	{
		REPORT_ERROR_STR("TomogramSet::TomogramSet: Unable to read " << filename);
	}
	
	globalTable.setName("global");
}

int TomogramSet::addTomogram(
		std::string tomoName, std::string stackFilename,
		const std::vector<gravis::d4Matrix>& projections, 
		int w, int h, int d, 
		const std::vector<double>& dose, 
		const std::vector<CTF>& ctfs, 
		double handedness, 
		double pixelSize)
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
	
	globalTable.setValue(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, pixelSize, index);
	globalTable.setValue(EMDL_CTF_VOLTAGE, ctf0.kV, index);
	globalTable.setValue(EMDL_CTF_CS, ctf0.Cs, index);
	globalTable.setValue(EMDL_CTF_Q0, ctf0.Q0, index);
	
	
	if (tomogramTables.size() != index)
	{
		REPORT_ERROR_STR("TomogramSet::add: corrupted tomogram set: tomogramTables.size() = "
						 << tomogramTables.size() << ", globalTable.numberOfObjects() = " 
						 << globalTable.numberOfObjects());
	}
	
	tomogramTables.push_back(MetaDataTable());
	MetaDataTable& m = tomogramTables[index];
	m.setName("tomo_" + ZIO::itoa(index));
		
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

Tomogram TomogramSet::loadTomogram(int index, bool loadImageData) const
{
	Tomogram out;
	
	std::string tomoName, stackFn;
	
	globalTable.getValueSafely(EMDL_TOMO_NAME, tomoName, index);
	globalTable.getValueSafely(EMDL_TOMO_TILT_SERIES_NAME, stackFn, index);
	globalTable.getValueSafely(EMDL_TOMO_FRAME_COUNT, out.frameCount, index);
		
	if (loadImageData)
	{
		out.stack.read(stackFn);
		out.hasImage = true;
	}
	else
	{
		out.hasImage = false;
	}
		
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
		
		ctf.initialise();
		
		m.getValueSafely(EMDL_MICROGRAPH_PRE_EXPOSURE, out.cumulativeDose[f], f);
	}
	
	out.frameSequence = IndexSort<double>::sortIndices(out.cumulativeDose);
	out.name = tomoName;
	
	if (globalTable.labelExists(EMDL_TOMO_FIDUCIALS_STARFILE))
	{
		 globalTable.getValue(EMDL_TOMO_FIDUCIALS_STARFILE, out.fiducialsFilename, index);
	}
	else
	{
		out.fiducialsFilename;
	}
	
	return out;
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
