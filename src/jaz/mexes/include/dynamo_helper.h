#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/image/buffered_image.h>

gravis::d4Matrix anglesToMatrix(
	double tdrot, double tilt, double narot, 
	double cx, double cy,
	double w, double h, double d	)
{
	const double phi = DEG2RAD(tdrot);
	const double theta = DEG2RAD(tilt);
	const double chi = DEG2RAD(narot);
	
	const double sinPhi = sin(phi), cosPhi = cos(phi);
	const double sinTheta = sin(theta), cosTheta = cos(theta);
	const double sinChi = sin(chi), cosChi = cos(chi);
	
	// Center at size/2
	gravis::d4Matrix Tc(
		1, 0, 0, -w/2,
		0, 1, 0, -h/2,
		0, 0, 1, -d/2,
		0, 0, 0, 1);
	
	// The particle rotates clockwise about its z axis.
	gravis::d4Matrix R0(
		 cosPhi, sinPhi, 0, 0,
		-sinPhi, cosPhi, 0, 0,
			  0,      0, 1, 0,
			  0,      0, 0, 1);
	
	// the particle rotates clockwise about its new x axis.
	gravis::d4Matrix R1(
		 1,         0,        0, 0,
		 0,  cosTheta, sinTheta, 0,
		 0, -sinTheta, cosTheta, 0,
		 0,         0,        0, 1);
	
	// the particle rotates clockwise about its new z axis.
	gravis::d4Matrix R2(
		 cosChi, sinChi, 0, 0,
		-sinChi, cosChi, 0, 0,
			  0,      0, 1, 0,
			  0,      0, 0, 1);
	
	// Shift to 3D position of center
	gravis::d4Matrix Ts(
		1, 0, 0, cx,
		0, 1, 0, cy,
		0, 0, 1, 0,
		0, 0, 0, 1);
	
	return Ts * R2 * R1 * R0 * Tc;
}
