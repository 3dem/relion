#ifndef RELION_JAZ_EULER_ANGLES_H
#define RELION_JAZ_EULER_ANGLES_H

#include <src/jaz/gravis/t4Matrix.h>
#include <src/euler.h>
#include <iostream>

class Euler
{
	public:

		inline static gravis::d4Matrix anglesToMatrix4(double phi, double theta, double chi);
		inline static gravis::d3Matrix anglesToMatrix3(double phi, double theta, double chi);

		inline static gravis::t4Vector<gravis::d3Matrix>
			anglesToMatrixAndDerivatives(double phi, double theta, double chi);

		inline static gravis::d3Vector matrixToAngles(const gravis::d4Matrix& A);
		inline static gravis::d3Vector matrixToAngles(const gravis::d3Matrix& A);


		inline static void test();
};

inline gravis::d4Matrix Euler::anglesToMatrix4(double phi, double theta, double chi)
{
	const double sp = sin(phi),   cp = cos(phi);
	const double st = sin(theta), ct = cos(theta);
	const double sc = sin(chi),   cc = cos(chi);
	
	return gravis::d4Matrix (
		 cc * ct * cp - sc * sp,   cc * ct * sp + sc * cp,  -cc * st,   0,
		-sc * ct * cp - cc * sp,  -sc * ct * sp + cc * cp,   sc * st,   0,
		 st * cp,                 -st * sp,                  ct,        0,
		 0,                        0,                        0,         1  );
}

inline gravis::d3Matrix Euler::anglesToMatrix3(double phi, double theta, double chi)
{
	const double sp = sin(phi),   cp = cos(phi);
	const double st = sin(theta), ct = cos(theta);
	const double sc = sin(chi),   cc = cos(chi);
	
	return gravis::d3Matrix (
		 cc * ct * cp - sc * sp,   cc * ct * sp + sc * cp,  -cc * st,
		-sc * ct * cp - cc * sp,  -sc * ct * sp + cc * cp,   sc * st,
		 st * cp,                  st * sp,                  ct       );
}

inline gravis::t4Vector<gravis::d3Matrix> 
	Euler::anglesToMatrixAndDerivatives(double phi, double theta, double chi)
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
			 cc * ct * dcp - sc * dsp,   cc * ct * dsp + sc * dcp,   0,
			-sc * ct * dcp - cc * dsp,  -sc * ct * dsp + cc * dcp,   0,
			 st * dcp,                   st * dsp,                   0 ),
				
		// d_A / d_theta:
		gravis::d3Matrix (
			 cc * dct * cp,   cc * dct * sp,  -cc * dst,
			-sc * dct * cp,  -sc * dct * sp,   sc * dst,
			 dst * cp,        dst * sp,        dct       ),

		// d_A / d_chi:
		gravis::d3Matrix (
			 dcc * ct * cp - dsc * sp,   dcc * ct * sp + dsc * cp,  -dcc * st,
			-dsc * ct * cp - dcc * sp,  -dsc * ct * sp + dcc * cp,   dsc * st,
			 0,                          0,                          0         ),

		gravis::d3Matrix (
			 cc * ct * cp - sc * sp,   cc * ct * sp + sc * cp,  -cc * st,
			-sc * ct * cp - cc * sp,  -sc * ct * sp + cc * cp,   sc * st,
			 st * cp,                  st * sp,                  ct       )
	);
		
}

inline gravis::d3Vector Euler::matrixToAngles(const gravis::d4Matrix& A)
{
	const double st2 = A(0,2) * A(0,2) + A(1,2) * A(1,2);
	const double theta = atan2(sqrt(st2), A(2,2));
	
	double phi(0.0), chi(0.0);

	const double eps = std::numeric_limits<double>::min();
	
	if (st2 > eps) // can distinguish phi and chi
	{
		phi = atan2(A(2,0), -A(2,1));
		chi = atan2(A(2,1),  A(2,0));
	}
	else // gimbal lock!
	{
		if (A(2,2) > 0)
		{
			phi = 0;
			chi = atan2(-A(1,0), A(0,0));
		}
		else
		{
			phi = 0;
			chi = atan2(A(1,0), -A(0,0));
		}
	}
	
	return gravis::d3Vector(phi, theta, chi);
}
 
inline gravis::d3Vector Euler::matrixToAngles(const gravis::d3Matrix& A)
{
	const double st2 = A(0,2) * A(0,2) + A(1,2) * A(1,2);
	const double theta = atan2(sqrt(st2), A(2,2));

	double phi(0.0), chi(0.0);

	const double eps = std::numeric_limits<double>::min();

	if (st2 > eps) // can distinguish phi and chi
	{
		phi = atan2(A(2,1),  A(2,0));
		chi = atan2(A(1,2), -A(0,2));
	}
	else // gimbal lock!
	{
		if (A(2,2) > 0)
		{
			phi = 0;
			chi = atan2(-A(1,0), A(0,0));
		}
		else
		{
			phi = 0;
			chi = atan2(A(1,0), -A(0,0));
		}
	}
	
	return gravis::d3Vector(phi, theta, chi);
}

void Euler::test()
{
	for (int i = 0; i < 5; i++)
	{
		const double phi   = 2 * PI * rand() / (double) RAND_MAX - PI;
		const double theta = 2 * PI * rand() / (double) RAND_MAX - PI;
		const double chi   = 2 * PI * rand() / (double) RAND_MAX - PI;

		std::cout << phi << ", " << theta   << ", " << chi << ":\n" << std::endl;

		const gravis::d3Matrix A = anglesToMatrix3(phi, theta, chi);

		Matrix2D<RFLOAT> A0(3,3);
		Euler_angles2matrix(RAD2DEG(phi), RAD2DEG(theta), RAD2DEG(chi), A0);

		for (int r = 0; r < 3; r++)
		{
			for (int c = 0; c < 3; c++)
			{
				std::cout << A(r,c) << " / " << A0(r,c) << " / " << (A(r,c) - A0(r,c)) << " \t ";
			}

			std::cout << std::endl;
		}

		const gravis::d3Vector ang = matrixToAngles(A);

		RFLOAT phi0_deg, theta0_deg, chi0_deg;

		Euler_matrix2angles(A0, phi0_deg, theta0_deg, chi0_deg);

		const double phi0 =   DEG2RAD(phi0_deg);
		const double theta0 = DEG2RAD(theta0_deg);
		const double chi0 =   DEG2RAD(chi0_deg);

		std::cout << std::endl;

		std::cout << ang[0] << " / " << phi0   << " / " << (ang[0] - phi0)   << std::endl;
		std::cout << ang[1] << " / " << theta0 << " / " << (ang[1] - theta0) << std::endl;
		std::cout << ang[2] << " / " << chi0   << " / " << (ang[2] - chi0)   << std::endl;

		std::cout << "\n\n";

		const double delta = 1e-6;

		gravis::t4Vector<gravis::d3Matrix> grad = Euler::anglesToMatrixAndDerivatives(phi, theta, chi);

		const gravis::d3Matrix A_phi   = anglesToMatrix3(phi + delta, theta, chi);
		const gravis::d3Matrix A_theta = anglesToMatrix3(phi, theta + delta, chi);
		const gravis::d3Matrix A_chi   = anglesToMatrix3(phi, theta, chi + delta);

		const gravis::d3Matrix num_grad_phi   = (A_phi - A) / delta;
		const gravis::d3Matrix num_grad_theta = (A_theta - A) / delta;
		const gravis::d3Matrix num_grad_chi   = (A_chi - A) / delta;

		std::cout << "phi: \n";
		std::cout << grad[0] << '\n';
		std::cout << num_grad_phi << '\n';
		std::cout << grad[0] - num_grad_phi << "\n\n";

		std::cout << "theta: \n";
		std::cout << grad[1] << '\n';
		std::cout << num_grad_theta << '\n';
		std::cout << grad[1] - num_grad_theta << "\n\n";

		std::cout << "chi: \n";
		std::cout << grad[2] << '\n';
		std::cout << num_grad_chi << '\n';
		std::cout << grad[2] - num_grad_chi << "\n\n";


		std::cout << "\n\n\n\n";

		std::cout << std::endl << std::endl;

	}

}

#endif
