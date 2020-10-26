#ifndef DYN_TAIT_BRYAN_ANGLES_H
#define DYN_TAIT_BRYAN_ANGLES_H

#include <src/jaz/gravis/t4Matrix.h>
#include <src/error.h>


class TaitBryan
{
	public:
		
		inline static gravis::d4Matrix anglesToMatrix4(double phi, double theta, double chi);
		inline static gravis::d3Matrix anglesToMatrix3(double phi, double theta, double chi);
		
		inline static gravis::t4Vector<gravis::d3Matrix> 
			anglesToMatrixAndDerivatives(double phi, double theta, double chi);
		
		inline static gravis::d3Vector matrixToAngles(const gravis::d4Matrix& A);
		inline static gravis::d3Vector matrixToAngles(const gravis::d3Matrix& A);
};

inline gravis::d4Matrix TaitBryan::anglesToMatrix4(double phi, double theta, double chi)
{
	return gravis::d4Matrix(anglesToMatrix3(phi, theta, chi));
}

inline gravis::d3Matrix TaitBryan::anglesToMatrix3(double phi, double theta, double chi)
{
	const double sp = sin(phi),   cp = cos(phi);
	const double st = sin(theta), ct = cos(theta);
	const double sc = sin(chi),   cc = cos(chi);
	
	gravis::d3Matrix Rx(
		  1,   0,   0,
		  0,  cp, -sp,
		  0,  sp,  cp  );
	
	gravis::d3Matrix Ry(
		 ct,   0,  st,
		  0,   1,   0,
		-st,   0,  ct  );
	
	gravis::d3Matrix Rz(
		 cc, -sc,   0,
		 sc,  cc,   0,
		  0,   0,   1  );
	
	return Rz * Ry * Rx;
}

inline gravis::t4Vector<gravis::d3Matrix> 
	TaitBryan::anglesToMatrixAndDerivatives(double phi, double theta, double chi)
{
	const double sp = sin(phi),   cp = cos(phi);
	const double st = sin(theta), ct = cos(theta);
	const double sc = sin(chi),   cc = cos(chi);
	
	const double dsp = cp, dcp = -sp;
	const double dst = ct, dct = -st;
	const double dsc = cc, dcc = -sc;
	
	gravis::d3Matrix Rx(
		  1,   0,   0,
		  0,  cp, -sp,
		  0,  sp,  cp  );
	
	gravis::d3Matrix Ry(
		 ct,   0,  st,
		  0,   1,   0,
		-st,   0,  ct  );
	
	gravis::d3Matrix Rz(
		 cc, -sc,   0,
		 sc,  cc,   0,
		  0,   0,   1  );
	
	gravis::d3Matrix dRx(
		  0,    0,    0,
		  0,  dcp, -dsp,
		  0,  dsp,  dcp  );
	
	gravis::d3Matrix dRy(
		 dct,   0,  dst,
		  0,    0,   0,
		-dst,   0,  dct  );
	
	gravis::d3Matrix dRz(
		 dcc, -dsc,   0,
		 dsc,  dcc,   0,
		  0,   0,     0  );
	
	return gravis::t4Vector<gravis::d3Matrix>(
		// d_A / d_phi:
		Rz * Ry * dRx,
		// d_A / d_theta:
		Rz * dRy * Rx,
		// d_A / d_chi:
		dRz * Ry * Rx,
		// A:
		Rz * Ry * Rx
	);		
}

inline gravis::d3Vector TaitBryan::matrixToAngles(const gravis::d4Matrix& A)
{
	REPORT_ERROR("TaitBryan::matrixToAngles: not implemented");
}
 
inline gravis::d3Vector TaitBryan::matrixToAngles(const gravis::d3Matrix& A)
{
	REPORT_ERROR("TaitBryan::matrixToAngles: not implemented");
}

#endif
