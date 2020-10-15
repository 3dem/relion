#include "symmetry.h"

std::vector<gravis::d4Matrix> Symmetry::getMatrices(std::string symmGroup)
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
					R0(0,0), R0(0,1), R0(0,2), R0(0,3),
					R0(1,0), R0(1,1), R0(1,2), R0(1,3),
					R0(2,0), R0(2,1), R0(2,2), R0(2,3),
					R0(3,0), R0(3,1), R0(3,2), R0(3,3) );
	}
	
	return R;
}
