#ifndef DYN_EULER_ANGLES_H
#define DYN_EULER_ANGLES_H

#include <stdexcept>
#include <limits>

#include <src/jaz/gravis/t4Matrix.h>

class EulerDynamo
{
	public:

		inline static gravis::d4Matrix anglesToMatrix4(double phi, double theta, double chi);
		inline static gravis::d3Matrix anglesToMatrix3(double phi, double theta, double chi);

		inline static gravis::t4Vector<gravis::d3Matrix>
			anglesToMatrixAndDerivatives(double phi, double theta, double chi);

		inline static gravis::d3Vector matrixToAngles(const gravis::d4Matrix& A);
		inline static gravis::d3Vector matrixToAngles(const gravis::d3Matrix& A);
};

inline gravis::d4Matrix EulerDynamo::anglesToMatrix4(double phi, double theta, double chi)
{
	const double sp = sin(phi),   cp = cos(phi);
	const double st = sin(theta), ct = cos(theta);
	const double sc = sin(chi),   cc = cos(chi);
	
	return gravis::d4Matrix (
		 cp * cc - sp * ct * sc,   sp * cc + cp * ct * sc,  st * sc,  0,
		-cp * sc - sp * ct * cc,  -sp * sc + cp * ct * cc,  st * cc,  0,
		 sp * st,                 -cp * st,                 ct,       0,
		 0,                        0,                        0,       1  );
}

inline gravis::d3Matrix EulerDynamo::anglesToMatrix3(double phi, double theta, double chi)
{
	const double sp = sin(phi),   cp = cos(phi);
	const double st = sin(theta), ct = cos(theta);
	const double sc = sin(chi),   cc = cos(chi);
	
	return gravis::d3Matrix (
		 cp * cc - sp * ct * sc,   sp * cc + cp * ct * sc,  st * sc,
		-cp * sc - sp * ct * cc,  -sp * sc + cp * ct * cc,  st * cc,
		 sp * st,                 -cp * st,                 ct       );
}

inline gravis::t4Vector<gravis::d3Matrix> 
	EulerDynamo::anglesToMatrixAndDerivatives(double phi, double theta, double chi)
{
	const double sp = sin(phi),   cp = cos(phi);
	const double st = sin(theta), ct = cos(theta);
	const double sc = sin(chi),   cc = cos(chi);
	
	const double dsp = cp, dcp = -sp;
	const double dst = ct, dct = -st;
	const double dsc = cc, dcc = -sc;
	
	return gravis::t4Vector<gravis::d3Matrix>(
				
		// d_A / d_phi:
		gravis::d3Matrix (
			 dcp * cc - dsp * ct * sc,  dsp * cc + dcp * ct * sc,  0,
			-dcp * sc - dsp * ct * cc, -dsp * sc + dcp * ct * cc,  0,
			 dsp * st,                 -dcp * st,                  0  ),
				
		// d_A / d_theta:
		gravis::d3Matrix (
			-sp * dct * sc,  cp * dct * sc,  dst * sc,
			-sp * dct * cc,  cp * dct * cc,  dst * cc,
			 sp * dst,      -cp * dst,       dct       ),
				
		// d_A / d_chi:
		gravis::d3Matrix (
			 cp * dcc - sp * ct * dsc,   sp * dcc + cp * ct * dsc,  st * dsc,
			-cp * dsc - sp * ct * dcc,  -sp * dsc + cp * ct * dcc,  st * dcc,
			 0,                          0,                         0  ),
				
		// A:
		gravis::d3Matrix (
			 cp * cc - sp * ct * sc,   sp * cc + cp * ct * sc,  st * sc,
			-cp * sc - sp * ct * cc,  -sp * sc + cp * ct * cc,  st * cc,
			 sp * st,                 -cp * st,                 ct       )
	);
		
}

inline gravis::d3Vector EulerDynamo::matrixToAngles(const gravis::d4Matrix& A)
{
	const double st2 = A(2,0) * A(2,0) + A(2,1) * A(2,1);
	const double theta = atan2(sqrt(st2), A(2,2));
	
	double phi(0.0), chi(0.0);

	const double eps = std::numeric_limits<double>::min();
	
	if (st2 > eps) // can distinguish phi and chi
	{
		phi = atan2(A(2,0), -A(2,1));
		chi = atan2(A(0,2),  A(1,2));
	}
	else // gimbal lock!
	{
		phi = atan2(A(0,1), A(0,0));
		chi = 0.0;
	}
	
	return gravis::d3Vector(phi, theta, chi);
}
 
inline gravis::d3Vector EulerDynamo::matrixToAngles(const gravis::d3Matrix& A)
{
	const double st2 = A(2,0) * A(2,0) + A(2,1) * A(2,1);
	const double theta = atan2(sqrt(st2), A(2,2));
	
	double phi(0.0), chi(0.0);

	const double eps = std::numeric_limits<double>::min();
	
	if (st2 > eps) // can distinguish phi and chi
	{
		phi = atan2(A(2,0), -A(2,1));
		chi = atan2(A(0,2),  A(1,2));
	}
	else // gimbal lock!
	{
		phi = atan2(A(0,1), A(0,0));
		chi = 0.0;
	}
	
	return gravis::d3Vector(phi, theta, chi);
}

#endif
