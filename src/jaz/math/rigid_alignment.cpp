#include "rigid_alignment.h"

RigidAlignment::RigidAlignment(
        double x, double y, double z, 
        double rot, double tilt, double psi)
:	
	position(x,y,z),
	rot(rot), tilt(tilt), psi(psi)
{	
}
