#ifndef DATA_SET_H
#define DATA_SET_H

#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/gravis/t3Matrix.h>
#include <src/jaz/gravis/t4Matrix.h>
#include "motion/trajectory.h"

class DataSet
{
	public:
		
		typedef enum {Relion, Dynamo} Type;
		
		static DataSet* load(std::string filename, std::string motionFilename);
		
		
		DataSet(Type type) : type(type) {}
		
			const Type type;
			bool hasMotion;
			std::vector<Trajectory> motionTrajectories;
		
		virtual std::vector<std::vector<int>> splitByTomogram() const = 0;
		virtual int getTotalParticleNumber() const = 0;
		virtual gravis::d3Vector getPosition(long int particle_id) const = 0;
		virtual gravis::d3Matrix getMatrix3x3(long int particle_id) const = 0;
		virtual gravis::d4Matrix getMatrix4x4(long int particle_id, double w, double h, double d) const = 0;
		virtual std::string getName(long int particle_id) const = 0;
		virtual int getHalfSet(long int particle_id) const = 0;
		virtual void moveParticleTo(long int particle_id, gravis::d3Vector pos) = 0;
		virtual void shiftParticleBy(long int particle_id, gravis::d3Vector shift) = 0;
		virtual void write(std::string fn) const = 0;

			std::vector<gravis::d3Vector> getTrajectoryInPix(long int particle_id, int fc, double pixelSize) const;
			void checkTrajectoryLengths(int p0, int np, int fc, std::string caller) const;
};

#endif
