#ifndef RELION_TOMO_PARTICLE_SET_H
#define RELION_TOMO_PARTICLE_SET_H

#include <src/metadata_table.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t3Matrix.h>
#include <src/jaz/gravis/t4Matrix.h>
#include "motion/trajectory.h"
#include "motion/trajectory_set.h"

class TomogramSet;

class ParticleIndex
{
	public:

		ParticleIndex(){}
		explicit ParticleIndex(long int value) : value(value) {}

		long int value;
};

class ParticleSet
{
	public:

		ParticleSet();
		ParticleSet(std::string filename, std::string motionFilename = "", bool verbose = true);
		
			MetaDataTable partTable, optTable;

			bool hasMotion;
			std::vector<Trajectory> motionTrajectories;


		void reserve(int particleNumber);
		ParticleIndex addParticle(const ParticleSet& particleSet, ParticleIndex index);
		void clearParticles();
		
		std::vector<std::vector<ParticleIndex>> splitByTomogram(const TomogramSet& tomogramSet, bool verbose = true) const;
		int getTotalParticleNumber() const;
		
		gravis::d3Vector getPosition(ParticleIndex particle_id) const;
		
		gravis::d3Matrix getSubtomogramMatrix(ParticleIndex particle_id) const;
		gravis::d3Matrix getParticleMatrix(ParticleIndex particle_id) const;
		gravis::d3Matrix getMatrix3x3(ParticleIndex particle_id) const;
		gravis::d4Matrix getMatrix4x4(ParticleIndex particle_id, double w, double h, double d) const;

		gravis::t4Vector<gravis::d3Matrix> getMatrixDerivativesOverParticleAngles(ParticleIndex particle_id) const;
		
		std::string getName(ParticleIndex particle_id) const;
		int getHalfSet(ParticleIndex particle_id) const;
		
		void moveParticleTo(ParticleIndex particle_id, gravis::d3Vector pos);
		void shiftParticleBy(ParticleIndex particle_id, gravis::d3Vector shift);
				
		void write(const std::string& filename) const;
		void writeTrajectories(const std::string& filename) const;
		
		void setImageFileNames(std::string data, std::string weight, ParticleIndex particle_id);
		
		gravis::d3Vector getParticleOffset(ParticleIndex particle_id) const;
		void setParticleOffset(ParticleIndex particle_id, const gravis::d3Vector& v);
		
		gravis::d3Vector getParticleCoord(ParticleIndex particle_id) const;
		void setParticleCoord(ParticleIndex particle_id, const gravis::d3Vector& v);

		int getOpticsGroup(ParticleIndex particle_id) const;
		void setOpticsGroup(ParticleIndex particle_id, int zeroBasedId);
		int numberOfOpticsGroups() const;
		
		double getBinnedPixelSize(int opticsGroup) const;
		double getOriginalPixelSize(int opticsGroup) const;


		std::vector<gravis::d3Vector> getTrajectoryInPixels(ParticleIndex particle_id, int fc, double pixelSize, bool from_original_coordinate = false) const;
		void checkTrajectoryLengths(const std::vector<ParticleIndex> &tomogram_particles, int fc, std::string caller) const;

		// Split tomograms into segments of similar total particle count to facilitate load balancing.
		// Specifically, minimise the number of particles in the segment containing the most particles.
		static std::vector<std::vector<int>> splitEvenly(
				const std::vector<std::vector<ParticleIndex>>& particlesByTomogram,
				int segment_count);

		// Simplified version of above for single-MPI-node versions.
		static std::vector<int> enumerate(
				const std::vector<std::vector<ParticleIndex>>& particlesByTomogram);

		static std::vector<int> enumerateNonEmpty(
				const std::vector<std::vector<ParticleIndex>>& particlesByTomogram);
};

#endif
