#include "data_set.h"
#include "dynamo_data_set.h"
#include "relion_data_set.h"
#include "motion/trajectory_set.h"

#include <src/macros.h>

using namespace gravis;


DataSet* DataSet::load(std::string filename, std::string motionFilename)
{
	size_t dot = filename.find_last_of(".");
	
	if (dot == std::string::npos)
	{
		REPORT_ERROR_STR("DataSet::load: file name has no ending: " << filename);
	}
	
	std::string ending = filename.substr(dot+1);
	
	DataSet* out;
	
	if (ending == "star")
	{
		out = new RelionDataSet(filename);
	}
	else
	{
		out = new DynamoDataSet(filename);
	}
	
	out->hasMotion = motionFilename != "";
			
	if (out->hasMotion)
	{
		out->motionTrajectories = Trajectory::read(motionFilename);
	}
	
	return out;
}

std::vector<d3Vector> DataSet::getTrajectoryInPix(long particle_id, int fc, double pixelSize) const
{
	const d3Vector p0 = getPosition(particle_id);
	
	if (hasMotion)
	{
		return motionTrajectories[particle_id].getShiftsInPix(p0, pixelSize);
	}
	else
	{
		return std::vector<d3Vector>(fc, p0);
	}
}

void DataSet::checkTrajectoryLengths(int p0, int np, int fc, std::string caller) const
{
	if (hasMotion)
	{
		for (int p = p0; p < p0 + np; p++)
		{
			if (motionTrajectories[p].shifts_Ang.size() != fc)
			{
				REPORT_ERROR_STR(caller << ": bad trajectory lengths; expected " << fc << " frames, found "
								 << motionTrajectories[p].shifts_Ang.size());
			}
		}
	}
}
