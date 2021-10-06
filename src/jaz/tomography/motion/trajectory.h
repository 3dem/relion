#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <src/jaz/gravis/t3Vector.h>
#include <vector>

class ParticleIndex;
class ParticleSet;

class Trajectory
{
	public:
		
		static std::vector<Trajectory> read(
				std::string filename,
				ParticleSet& particleSet);

		static void write(
				const std::vector<std::vector<Trajectory>>& shifts,
				const ParticleSet& particleSet,
				const std::vector<std::vector<ParticleIndex>>& particles,
				std::string filename);
		
		
		Trajectory();
		Trajectory(int fc);
		
			std::vector<gravis::d3Vector> shifts_Ang;
			
			
		std::vector<gravis::d3Vector> getShiftsInPix(gravis::d3Vector origin, double pixelSize) const;
		
		Trajectory& operator += (const Trajectory& t)
		{
			const int fc = shifts_Ang.size();
			const int fct = t.shifts_Ang.size();
			
			if (fc != fct)
			{
				shifts_Ang = t.shifts_Ang;
			}
			else
			{
				for (int f = 0; f < fc; f++)
				{
					shifts_Ang[f] += t.shifts_Ang[f];
				}
			}

			return *this;	
		}
		
};

Trajectory operator + (const Trajectory& t1, const Trajectory& t2);


#endif
