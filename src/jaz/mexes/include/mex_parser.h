#ifndef DYN_MEX_PARSER_H
#define DYN_MEX_PARSER_H

#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <string>
#include <map>
#include <sstream>

#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/image/buffered_image.h>

// ----- ----- ----- -----  ----- ----- ----- ----- - ----- ----- ----- -----  ----- ----- ----- -----


// compile with "mex -R2018a ..."

#define DYN_MATLAB_R2018a 1


// ----- ----- ----- -----  ----- ----- ----- ----- - ----- ----- ----- -----  ----- ----- ----- -----


class MexParser
{
	public:
				
		typedef enum
		{
			Left = 0, 
			Right = 1
		}
		Side;

		
		MexParser(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
	
		
			double getDouble(Side side, int id, std::string desc);	
			
			double getDouble(Side side, int id, std::string desc, double defaultVal);	
			
			template<typename T>
			RawImage<T> getImage(Side side, int id, std::string desc);	
			
			template<typename T>
			RawImage<T> createImage(Side side, int id, std::string desc, int w, int h, int d);	
			
			std::vector<gravis::d2Vector> getR2Vector(Side side, int id, std::string desc);	
			
			std::vector<gravis::d3Vector> getR3Vector(Side side, int id, std::string desc);
			
			std::string getString(Side side, int id, std::string desc);
			
			bool checkOption(Side side, int id, std::string desc, bool def = false);
			
			
			bool finalize();
	
	
	protected:
		
			int nlhs, nrhs;
			mxArray **plhs;
			const mxArray **prhs;
		
			int argNums[2];
			bool anyErrors;
			std::stringstream errors;
			std::map<int, std::string> usage[2];
			
			std::string sideNames[2];
		
		std::string formatUsage();
			
		template <typename T>
		bool checkArg(Side side, int id, std::string desc, bool toBeAllocated = false);
		
		template <typename T>
		bool checkType(Side side, int id);
		
		void wrongType(Side side, int id, std::string expected, std::string found);
		
		template<typename T>
		T* getPointer(const mxArray* arr);
		
		template<typename T>
		mxArray* allocateArray(int id, int w, int h, int d);
				
		const mxArray* getArgument(Side side, int id);
	
};


// template specializations have to come first

template<>
double* MexParser::getPointer<double>(const mxArray* arr)
{
	#if DYN_MATLAB_R2018a
	
		return mxGetDoubles(arr);
	
	#else
	
		cannot compile without defining a standard!
	
	#endif
}

template<>
float* MexParser::getPointer<float>(const mxArray* arr)
{	
	#if DYN_MATLAB_R2018a
	
		return mxGetSingles(arr);
	
	#else
	
		cannot compile without defining a standard!
	
	#endif
}

// 

MexParser::MexParser(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
:	nlhs(nlhs),
	nrhs(nrhs),
	plhs(plhs),
	prhs(prhs),
	anyErrors(false)
{
	argNums[Left] = nlhs;
	argNums[Right] = nrhs;
	sideNames[Left] = "left";
	sideNames[Right] = "right";
}

double MexParser::getDouble(Side side, int id, std::string desc)
{
	if (!checkArg<double>(side, id, desc))
	{
		return 0.0;
	}
	
	return mxGetDoubles(prhs[id])[0];
}

double MexParser::getDouble(Side side, int id, std::string desc, double defaultVal)
{
	std::stringstream sts;
	sts << desc << " (default: " << defaultVal << ")";
	
	if (argNums[side] <= id) 
	{
		usage[side][id] = sts.str();
		return defaultVal;	
	}
	else 
	{
		return getDouble(side, id, sts.str());
	}
}

template<typename T>
RawImage<T> MexParser::getImage(Side side, int id, std::string desc)
{	
	if (!checkArg<T>(side,id,desc))
	{
		return RawImage<T>(0,0,0,0);
	}
	else
	{
		const mxArray* arr = getArgument(side, id);
		
		mwSize dim_num = mxGetNumberOfDimensions(arr);	
		const mwSize* dims = mxGetDimensions(arr);
		
		const mwSize w = dims[0];
		const mwSize h = dims[1];
		const mwSize d = dim_num < 3? 1 : dims[2];
		
		T* data = getPointer<T>(arr);
		
		return RawImage<T>(w, h, d, data);
	}
}

template<typename T>
RawImage<T> MexParser::createImage(Side side, int id, std::string desc, int w, int h, int d)
{	
	if (side == Right)
	{
		errors << "* Cannot allocate images as right-hand side arguments (argument #" 
			   << (id+1) << ", " << desc << ").\n";
		
		anyErrors = true;
		
		return RawImage<T>(0,0,0,0);
	}
	
	if (!checkArg<T>(Left,id,desc,true))
	{
		return RawImage<T>(0,0,0,0);
	}
	else
	{
		mxArray* arr = allocateArray<T>(id, w,h,d);
		
		return RawImage<T>(w, h, d, getPointer<T>(arr));
	}
}

std::vector<gravis::d2Vector> MexParser::getR2Vector(Side side, int id, std::string desc)
{
	if (!checkArg<double>(side,id,desc))
	{
		return std::vector<gravis::d2Vector>(0);
	}
	
	const mxArray* arr = getArgument(side, id);		
	const mwSize* dims = mxGetDimensions(arr);
	
	if (dims[0] != 2)
	{
		errors << "* " << sideNames[side] << "-hand side argument #" << (id + 1) 
			<< " (" << desc << ") must have 2 rows.\n";
	 
		anyErrors = true;
		
		return std::vector<gravis::d2Vector>(0);
	}
	
	const mwSize num = dims[1];
	double* data = getPointer<double>(arr);
	
	std::vector<gravis::d2Vector> out(num);
	
	for (int i = 0; i < num; i++)
	{
		out[i] = gravis::d2Vector(data[2*i], data[2*i + 1]);
	}
	
	return out;
}

std::vector<gravis::d3Vector> MexParser::getR3Vector(Side side, int id, std::string desc)
{
	if (!checkArg<double>(side,id,desc))
	{
		return std::vector<gravis::d3Vector>(0);
	}
	
	const mxArray* arr = getArgument(side, id);		
	const mwSize* dims = mxGetDimensions(arr);
	
	if (dims[0] != 3)
	{
		errors << "* " << sideNames[side] << "-hand side argument #" << (id + 1) 
			<< " (" << desc << ") must have 3 rows.\n";
	 
		anyErrors = true;
		
		return std::vector<gravis::d3Vector>(0);
	}
	
	const mwSize num = dims[1];
	double* data = getPointer<double>(arr);
	
	std::vector<gravis::d3Vector> out(num);
	
	for (int i = 0; i < num; i++)
	{
		out[i] = gravis::d3Vector(data[3*i], data[3*i + 1], data[3*i + 2]);
	}
	
	return out;
}

std::string MexParser::getString(MexParser::Side side, int id, std::string desc)
{
	if (!checkArg<std::string>(side, id, desc))
	{
		return "";
	}
	
	const mxArray* arg = getArgument(side, id);
	
	int size = (mxGetM(arg) * mxGetN(arg)) + 1;
	char* buf = (char*) mxCalloc(size, sizeof(char));
	
	mxGetString(arg, buf, size);
	  
	std::string out(buf);
	
	mxFree(buf);
	
	return out;
}

bool MexParser::checkOption(Side side, int id, std::string desc, bool def)
{
	double d = getDouble(side, id, desc, def? 1.0 : 0.0);
	
	if (d == 0.0) 
	{
		return false;
	}
	else if (d == 1.0)
	{
		return true;
	}
	else
	{
		errors << "* Value of " << sideNames[side] << "-hand side argument #" << (id + 1) 
				<< " has to be either 0 or 1\n";
		
		anyErrors = true;
		
		return true;
	}
}

bool MexParser::finalize()
{
	if (anyErrors)
	{
		std::string text = formatUsage();
		
		mexPrintf(text.c_str());
		mexErrMsgTxt(text.c_str());
		
		return false;
	}
	else
	{
		return true;
	}
}


// ----- - ----- ----- - ----- ----- - ----- ----- - ----- - - - ----- - ----- ----- - ----- ----- - ----- ----- - ----- //

std::string MexParser::formatUsage()
{
	std::stringstream sts;
	
	if (nlhs != 0 || nrhs != 0)
	{
		sts << "\nerrors encountered:\n";	
		sts << errors.str() << "\n";
	}
	
	sts << "Expected arguments:\n\n";
	
	for (int s = 0; s < 2; s++)
	{
		if (usage[s].size() > 0)
		{
			sts << sideNames[s] << "-hand side:\n";
			
			for (std::map<int,std::string>::iterator it = usage[s].begin();
				 it != usage[s].end(); it++)
			{
				sts << (it->first + 1) << ": \t" << it->second << "\n";
			}
			
			sts << "\n";
		}
	}
	
	return sts.str();
}

template <typename T>
bool MexParser::checkArg(Side side, int id, std::string desc, bool toBeAllocated)
{
	usage[side][id] = desc;
	
	if (id >= argNums[side])
	{
		errors << "* Missing " << sideNames[side] << "-hand side argument #" << (id + 1) 
		       << " (only " << argNums[side] << " provided).\n";
		
		anyErrors = true;
		
		return false;
	}
	else
	{
		return toBeAllocated || checkType<T>(side, id);
	}
}

template <>
bool MexParser::checkType<double>(Side side, int id)
{
	const mxArray* arr = getArgument(side, id);
	
	if (!mxIsDouble(arr))
	{
		wrongType(side, id, "double", mxGetClassName(arr));
		
		return false;
	}
	else 
	{
		return true;
	}
}

template <>
bool MexParser::checkType<float>(Side side, int id)
{
	const mxArray* arr = getArgument(side, id);
	
	if (!mxIsSingle(arr))
	{
		wrongType(side, id, "single", mxGetClassName(arr));
		
		return false;
	}
	else 
	{
		return true;
	}
}

template <>
bool MexParser::checkType<std::string>(Side side, int id)
{
	const mxArray* arr = getArgument(side, id);
	
	if (!mxIsChar(arr))
	{
		wrongType(side, id, "char", mxGetClassName(arr));
		
		return false;
	}
	else 
	{
		return true;
	}
}

void MexParser::wrongType(Side side, int id, std::string expected, std::string found)
{
	errors << "* Type of " << sideNames[side] << "-hand side argument #" << (id + 1) 
			<< " should be " << expected << " (" << found << " found)\n";
	
	anyErrors = true;
}

template<>
mxArray* MexParser::allocateArray<double>(int id, int w, int h, int d)
{
	const mwSize dims[] = {(mwSize)w, (mwSize)h, (mwSize)d};	
	plhs[id] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	return plhs[id];
}

template<>
mxArray* MexParser::allocateArray<float>(int id, int w, int h, int d)
{
	const mwSize dims[] = {(mwSize)w, (mwSize)h, (mwSize)d};	
	plhs[id] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
	return plhs[id];
}

const mxArray* MexParser::getArgument(Side side, int id)
{
	switch (side)
	{
		case Left: return plhs[id];
		case Right: return prhs[id];
	}
}

#endif
