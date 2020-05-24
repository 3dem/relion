#include "dynamo_data_set.h"
#include <fstream>

using namespace gravis;


DynamoDataSet::DynamoDataSet()
{}

DynamoDataSet::DynamoDataSet(std::string filename)
{
	std::ifstream ifs(filename);
	
	const int maxlen = 4096;
	char buffer[maxlen];
	
	while (ifs.getline(buffer, maxlen))
	{
		particles.push_back(DynamoParticle(buffer));
	}
}

std::vector<std::vector<int>> DynamoDataSet::splitByTomogram() const
{
	int maxTomo = 0;
	
	for (int p = 0; p < particles.size(); p++)
	{
		const int t = particles[p].tomo;
		
		if (t > maxTomo) maxTomo = t;
	}
	
	std::vector<std::vector<int>> out(maxTomo+1);
	
	for (int p = 0; p < particles.size(); p++)
	{
		out[particles[p].tomo].push_back(p);
	}
	
	return out;
}

int DynamoDataSet::getTotalParticleNumber() const
{
	return particles.size();
}

d3Vector DynamoDataSet::getPosition(long int particle_id) const
{
	return particles[particle_id].getPosition();
}

d3Matrix DynamoDataSet::getMatrix3x3(long int particle_id) const
{
	d4Matrix A = particles[particle_id].getAlignmentMatrixAlibi4x4(0,0,0);
	return d3Matrix::extract(A);
}

d4Matrix DynamoDataSet::getMatrix4x4(long int particle_id, double w, double h, double d) const
{
	return particles[particle_id].getAlignmentMatrixAlibi4x4(w,h,d);
}

std::string DynamoDataSet::getName(long int particle_id) const
{
	return particles[particle_id].getFormattedTag();
}

int DynamoDataSet::getHalfSet(long int particle_id) const
{
	return particle_id % 2;
}

void DynamoDataSet::moveParticleTo(long int particle_id, gravis::d3Vector pos)
{
	particles[particle_id].x = pos.x;
	particles[particle_id].y = pos.y;
	particles[particle_id].z = pos.z;
}

void DynamoDataSet::shiftParticleBy(long int particle_id, gravis::d3Vector shift)
{
	particles[particle_id].x += shift.x;
	particles[particle_id].y += shift.y;
	particles[particle_id].z += shift.z;
}


void DynamoDataSet::write(std::string fn) const
{
	std::ofstream ofs(fn);

    for (int p = 0; p < particles.size(); p++)
    {
        ofs << particles[p] << '\n';
    }
}
