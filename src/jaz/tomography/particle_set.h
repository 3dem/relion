#ifndef RELION_TOMO_PARTICLE_SET_H
#define RELION_TOMO_PARTICLE_SET_H

#include <src/metadata_table.h>
#include <src/jaz/gravis/t3Matrix.h>
#include <src/jaz/gravis/t4Matrix.h>
#include "motion/trajectory.h"
#include "motion/trajectory_set.h"

class TomogramSet;

class ParticleIndex
{
	public:

		ParticleIndex(int value) : value(value) {}

		int value;
};

class ParticleSet
{
	public:

		ParticleSet();
		ParticleSet(std::string filename, std::string motionFilename = "");
		
			MetaDataTable partTable, optTable;

			bool hasMotion;
			std::vector<Trajectory> motionTrajectories;

		
		std::vector<std::vector<int>> splitByTomogram(const TomogramSet& tomogramSet) const;
		int getTotalParticleNumber() const;
		
		gravis::d3Vector getPosition(long int particle_id) const;
		
		gravis::d3Matrix getSubtomogramMatrix(long int particle_id) const;
		gravis::d3Matrix getParticleMatrix(long int particle_id) const;
		gravis::d3Matrix getMatrix3x3(long int particle_id) const;
		gravis::d4Matrix getMatrix4x4(long int particle_id, double w, double h, double d) const;
		
		std::string getName(long int particle_id) const;
		int getHalfSet(long int particle_id) const;
		
		void moveParticleTo(long int particle_id, gravis::d3Vector pos);
		void shiftParticleBy(long int particle_id, gravis::d3Vector shift);
				
		void write(std::string fn) const;
		
		void setImageFileNames(std::string data, std::string weight, long int particle_id);
		
		gravis::d3Vector getParticleOffset(long int particle_id) const;
		void setParticleOffset(long int particle_id, const gravis::d3Vector& v);
		
		gravis::d3Vector getParticleCoord(long int particle_id) const;
		void setParticleCoord(long int particle_id, const gravis::d3Vector& v);

		int getOpticsGroup(long int particle_id) const;
		int numberOfOpticsGroups() const;
		
		double getBinnedPixelSize(int opticsGroup) const;
		double getOriginalPixelSize(int opticsGroup) const;


		std::vector<gravis::d3Vector> getTrajectoryInPixels(long int particle_id, int fc, double pixelSize) const;
		void checkTrajectoryLengths(int p0, int np, int fc, std::string caller) const;
		
		static gravis::d3Matrix convert(const Matrix2D<double>& A);
};

#endif
