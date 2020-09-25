#ifndef RIGID_ALIGNMENT
#define RIGID_ALIGNMENT

#include <src/jaz/gravis/t3Vector.h>

class RigidAlignment
{
	public:
		
		RigidAlignment(double x, double y, double z, double rot, double tilt, double psi);
		
			gravis::d3Vector position;
			double rot, tilt, psi;
};

#endif
