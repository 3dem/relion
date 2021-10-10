#include "symmetry.h"

std::vector<gravis::d4Matrix> Symmetry::getPointGroupMatrices(
		std::string symmGroup)
{
	SymList SL;
	SL.read_sym_file(symmGroup);
	
	const int sc = SL.SymsNo();
	
	std::vector<gravis::d4Matrix> R(sc);
			
	for (int isym = 0; isym < sc; isym++)
	{
		Matrix2D<RFLOAT> L0, R0;
		SL.get_matrices(isym, L0, R0);
		
		R[isym] = gravis::d4Matrix(
					R0(0,0), R0(0,1), R0(0,2), 0,
					R0(1,0), R0(1,1), R0(1,2), 0,
					R0(2,0), R0(2,1), R0(2,2), 0,
					R0(3,0), R0(3,1), R0(3,2), 0 );
	}
	
	return R;
}

std::vector<gravis::d4Matrix> Symmetry::getHelicalSymmetryMatrices(
		int units, double twist, double rise)
{
	const int h_min = -units/2;
	const int h_max = -h_min + units%2;
	const int num = h_max - h_min - 1;

	std::vector<gravis::d4Matrix> R(num);

	Matrix2D<RFLOAT> R0(4, 4);

	for (int hh = h_min; hh < h_max; hh++)
	{
		const int index = hh < 0? hh - h_min : hh - h_min - 1;

		if (hh != 0)
		{
			const double rot_ang = hh * (-twist);
			rotation3DMatrix(rot_ang, 'Z', R0);

			R[index] = gravis::d4Matrix(
						R0(0,0), R0(0,1), R0(0,2), 0,
						R0(1,0), R0(1,1), R0(1,2), 0,
						R0(2,0), R0(2,1), R0(2,2), -hh * rise,
						R0(3,0), R0(3,1), R0(3,2), 0 );

		}
	}

	return R;
}
