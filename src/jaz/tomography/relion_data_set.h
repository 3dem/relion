#ifndef RELION_TOMO_DATA_SET_H
#define RELION_TOMO_DATA_SET_H

#include "data_set.h"
#include <src/metadata_table.h>


class RelionDataSet : public DataSet
{
	public:
		
		RelionDataSet();
		RelionDataSet(std::string filename);
		
			MetaDataTable partTable, optTable;
			std::vector<double> binnedPixelSizes, originalPixelSizes;
		
		std::vector<std::vector<int>> splitByTomogram() const;
		int getTotalParticleNumber() const;
		gravis::d3Vector getPosition(long int particle_id) const;
		gravis::d3Matrix getMatrix3x3(long int particle_id) const;
		gravis::d4Matrix getMatrix4x4(long int particle_id, double w, double h, double d) const;	
		std::string getName(long int particle_id) const;
		int getHalfSet(long int particle_id) const;
		
		void moveParticleTo(long int particle_id, gravis::d3Vector pos);
		void shiftParticleBy(long int particle_id, gravis::d3Vector shift);
				
		void write(std::string fn) const;
		
		void setImageFileNames(std::string data, std::string weight, long int particle_id);
		
		void getParticleOffset(long int particle_id, double& x, double& y, double& z) const;
		void setParticleOffset(long int particle_id, double x, double y, double z);
		
		void getParticleCoord(long int particle_id, double& x, double& y, double& z) const;
		void setParticleCoord(long int particle_id, double x, double y, double z);
		
		int getOpticsGroup(long int particle_id) const;
		
		double getBinnedPixelSize(int opticsGroup) const;
		double getOriginalPixelSize(int opticsGroup) const;
};

#endif
