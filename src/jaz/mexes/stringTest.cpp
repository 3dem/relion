#include "include/mex_parser.h"
#include "include/dynamo_helper.h"

#include <iostream>


using namespace gravis;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MexParser mp(nlhs, plhs, nrhs, prhs);

	std::string s = mp.getString(MexParser::Right, 0, "test string");
	
	if (!mp.finalize()) return;
	
	const int l = s.length();
	
	for (int i = 0; i < l; i++)
	{
		std::cout << s[l-i-1];
	}
	
	std::cout << std::endl;
}
