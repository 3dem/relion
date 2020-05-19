#include "relion_data_set.h"
#include <src/euler.h>

using namespace gravis;


RelionDataSet::RelionDataSet()
: DataSet(Relion)
{}

RelionDataSet::RelionDataSet(std::string filename)
: DataSet(Relion)
{
	optTable.read(filename, "optics");
	partTable.read(filename, "particles");
	
	// partTable.newSort(EMDL_MICROGRAPH_NAME);
	
	/*  ISSUE: 
	  
		- Tomograms have to be already grouped inside the particles table.
		- Sorting them here might rearrange them.
		- Their order has to be the same as in the tomolist (since there is no tomo index field yet)
	*/
	
	if (!optTable.labelExists(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE))
	{
		REPORT_ERROR("RelionDataSet::RelionDataSet: "
					 + EMDL::label2Str(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE)
					 + " missing from optics MetaDataTable.\n");
	}
	
	const int groupCount = optTable.numberOfObjects();
	
	originalPixelSizes.resize(groupCount);
	binnedPixelSizes.resize(groupCount);
	
	for (int i = 0; i < groupCount; i++)
	{
		optTable.getValueSafely(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, originalPixelSizes[i], i);
		optTable.getValueSafely(EMDL_IMAGE_PIXEL_SIZE, binnedPixelSizes[i], i);
	}
}

std::vector<std::vector<int> > RelionDataSet::splitByTomogram() const
{
	std::vector<std::vector<int>> out(0);

	if (!partTable.labelExists(EMDL_MICROGRAPH_NAME))
	{
		REPORT_ERROR("RelionDataSet::splitByTomogram: "
					 + EMDL::label2Str(EMDL_MICROGRAPH_NAME)
					 + " missing from MetaDataTable.\n");
	}

	const long lc = partTable.numberOfObjects();
	std::string lastName = "", curName;
	long curInd = -1;

	for (int i = 0; i < lc; i++)
	{
		partTable.getValueSafely(EMDL_MICROGRAPH_NAME, curName, i);

		if (curName != lastName)
		{
			lastName = curName;
			curInd++;
			out.push_back(std::vector<int>(0));
		}

		out[curInd].push_back(i);
	}

	return out;
}

int RelionDataSet::getTotalParticleNumber() const
{
	return partTable.numberOfObjects();
}

d3Vector RelionDataSet::getPosition(long int particle_id) const
{
	d3Vector pos, off;
	
	partTable.getValueSafely(EMDL_IMAGE_COORD_X, pos.x, particle_id);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Y, pos.y, particle_id);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Z, pos.z, particle_id);
	
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_X_ANGSTROM, off.x, particle_id);
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, off.y, particle_id);
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, off.z, particle_id);
		
	const int og = getOpticsGroup(particle_id);
	d3Vector out = (binnedPixelSizes[og] * pos - off) / originalPixelSizes[og];
	
	out.x += 1.0;
	out.y += 1.0;
	out.z += 1.0;
	
	return out;
}

d3Matrix RelionDataSet::getMatrix3x3(long int particle_id) const
{
	Matrix2D<RFLOAT> A;
	
	double phi, theta, psi;
	
	partTable.getValueSafely(EMDL_ORIENT_ROT,  phi,   particle_id);
	partTable.getValueSafely(EMDL_ORIENT_TILT, theta, particle_id);
	partTable.getValueSafely(EMDL_ORIENT_PSI,  psi,   particle_id);
	
	Euler_angles2matrix(phi, theta, psi, A, false);
	
	return d3Matrix(
		A(0,0), A(0,1), A(0,2), 
		A(1,0), A(1,1), A(1,2), 
		A(2,0), A(2,1), A(2,2) );
}

d4Matrix RelionDataSet::getMatrix4x4(long int particle_id, double w, double h, double d) const
{
	Matrix2D<RFLOAT> A;
	
	double phi, theta, psi;
	
	partTable.getValueSafely(EMDL_ORIENT_ROT,  phi,   particle_id);
	partTable.getValueSafely(EMDL_ORIENT_TILT, theta, particle_id);
	partTable.getValueSafely(EMDL_ORIENT_PSI,  psi,   particle_id);
	
	Euler_angles2matrix(phi, theta, psi, A, false);
	
	d3Vector pos = getPosition(particle_id);
	
	int cx = ((int)w) / 2;
	int cy = ((int)h) / 2;
	int cz = ((int)d) / 2;
	
	gravis::d4Matrix Tc(
		1, 0, 0, -cx,
		0, 1, 0, -cy,
		0, 0, 1, -cz,
		0, 0, 0, 1);
	
	d4Matrix R(
		A(0,0), A(0,1), A(0,2), pos.x, 
		A(1,0), A(1,1), A(1,2), pos.y, 
		A(2,0), A(2,1), A(2,2), pos.z, 
		0.0,    0.0,    0.0,    1.0   );
	
	d4Matrix Ts(
		1, 0, 0, pos.x,
		0, 1, 0, pos.y,
		0, 0, 1, pos.z,
		0, 0, 0, 1);
	
	return Ts * R * Tc;
}

std::string RelionDataSet::getName(long int particle_id) const
{
	std::stringstream sts;
	sts << particle_id;
	
	return sts.str();
}

int RelionDataSet::getHalfSet(long int particle_id) const
{
	int s;
	partTable.getValueSafely(EMDL_PARTICLE_RANDOM_SUBSET, s, particle_id);
	return s - 1;
}


void RelionDataSet::moveParticleTo(long int particle_id, gravis::d3Vector pos)
{
	partTable.setValue(EMDL_IMAGE_COORD_X, pos.x - 1.0, particle_id);
	partTable.setValue(EMDL_IMAGE_COORD_Y, pos.y - 1.0, particle_id);
	partTable.setValue(EMDL_IMAGE_COORD_Z, pos.z - 1.0, particle_id);
	
	partTable.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, 0.0, particle_id);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, 0.0, particle_id);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, 0.0, particle_id);
}

void RelionDataSet::shiftParticleBy(long int particle_id, gravis::d3Vector shift)
{
	double x, y, z;
	
	partTable.getValueSafely(EMDL_IMAGE_COORD_X, x, particle_id);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Y, y, particle_id);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Z, z, particle_id);
	
	partTable.setValue(EMDL_IMAGE_COORD_X, x + shift.x, particle_id);
	partTable.setValue(EMDL_IMAGE_COORD_Y, y + shift.y, particle_id);
	partTable.setValue(EMDL_IMAGE_COORD_Z, z + shift.z, particle_id);
}

void RelionDataSet::write(std::string fn) const
{
	std::ofstream ofs(fn);
	
	optTable.write(ofs);
	partTable.write(ofs);
}

void RelionDataSet::setImageFileNames(std::string data, std::string weight, long int particle_id)
{
	partTable.setValue(EMDL_IMAGE_NAME, data, particle_id);
	partTable.setValue(EMDL_CTF_IMAGE, weight, particle_id);
}

void RelionDataSet::getParticleOffset(long particle_id, double& x, double& y, double& z) const
{
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_X_ANGSTROM, x, particle_id);
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, y, particle_id);
	partTable.getValueSafely(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, z, particle_id);
}

void RelionDataSet::setParticleOffset(long particle_id, double x, double y, double z)
{
	partTable.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, x, particle_id);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, y, particle_id);
	partTable.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, z, particle_id);
}

void RelionDataSet::getParticleCoord(long particle_id, double& x, double& y, double& z) const
{
	partTable.getValueSafely(EMDL_IMAGE_COORD_X, x, particle_id);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Y, y, particle_id);
	partTable.getValueSafely(EMDL_IMAGE_COORD_Z, z, particle_id);
}

void RelionDataSet::setParticleCoord(long particle_id, double x, double y, double z)
{
	partTable.setValue(EMDL_IMAGE_COORD_X, x, particle_id);
	partTable.setValue(EMDL_IMAGE_COORD_Y, y, particle_id);
	partTable.setValue(EMDL_IMAGE_COORD_Z, z, particle_id);
}

int RelionDataSet::getOpticsGroup(long particle_id) const
{
	if (!partTable.labelExists(EMDL_IMAGE_OPTICS_GROUP))
	{
		REPORT_ERROR("RelionDataSet::getPixelSize: pixel size (rlnImagePixelSize) missing from optics table");
	}
	
	int out;
	partTable.getValueSafely(EMDL_IMAGE_OPTICS_GROUP, out, particle_id);
	return out - 1;
}

int RelionDataSet::numberOfOpticsGroups() const
{
	return optTable.numberOfObjects();
}

double RelionDataSet::getBinnedPixelSize(int opticsGroup) const
{
	if (!optTable.labelExists(EMDL_IMAGE_PIXEL_SIZE))
	{
		REPORT_ERROR("RelionDataSet::getBinnedPixelSize: pixel size (rlnImagePixelSize) missing from optics table");
	}
	
	double out;
	optTable.getValueSafely(EMDL_IMAGE_PIXEL_SIZE, out, opticsGroup);
	return out;
}

double RelionDataSet::getOriginalPixelSize(int opticsGroup) const
{
	if (!optTable.labelExists(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE))
	{
		REPORT_ERROR("RelionDataSet::getOriginalPixelSize: pixel size (rlnMicrographOriginalPixelSize) missing from optics table");
	}
	
	double out;
	optTable.getValueSafely(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, out, opticsGroup);
	return out;
}
