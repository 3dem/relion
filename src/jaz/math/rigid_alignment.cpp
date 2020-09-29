#include "rigid_alignment.h"

using namespace gravis;

RigidAlignment::RigidAlignment(gravis::d3Vector position, double rot, double tilt, double psi)
	:
	position(position),
	rot(rot), tilt(tilt), psi(psi)
{
}
