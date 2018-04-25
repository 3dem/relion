#ifndef EQUATION_2X2_h
#define EQUATION_2X2_h

#include <src/jaz/gravis/t2Vector.h>

class Equation2x2
{
	public:
		
		Equation2x2();
		
		double Axx, Axy, Ayy, bx, by;
		
		Equation2x2& operator += (const Equation2x2& arg)
        {
            Axx += arg.Axx;
			Axy += arg.Axy;
			Ayy += arg.Ayy;
			bx += arg.bx;
			by += arg.by;

            return *this;
        }
};

#endif
