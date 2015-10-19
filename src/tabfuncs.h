/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef TABFUNCS_H_
#define TABFUNCS_H_

#include "src/multidim_array.h"
#include "src/funcs.h"

// Class to tabulate some functions
class TabFunction
{

protected:
	MultidimArray<RFLOAT> tabulatedValues;
	RFLOAT  sampling;
public:
	// Empty constructor
	TabFunction() {}

	// Destructor
    virtual ~TabFunction()
    {
    	tabulatedValues.clear();
    }

    /** Copy constructor
     *
     * The created TabFunction is a perfect copy of the input array but with a
     * different memory assignment.
     */
    TabFunction(const TabFunction& op)
    {
    	tabulatedValues.clear();
    	*this = op;
    }

	/** Assignment.
     *
     * You can build as complex assignment expressions as you like. Multiple
     * assignment is allowed.
     */
    TabFunction& operator=(const TabFunction& op)
    {
        if (&op != this)
        {
         	// Projector stuff (is this necessary in C++?)
        	tabulatedValues = op.tabulatedValues;
        	sampling = op.sampling;
        }
        return *this;
    }


};

class TabSine : public TabFunction
{
public:
	// Empty constructor
	TabSine() {}

	// Constructor (with parameters)
	void initialise(const int _nr_elem = 5000);

	//Pre-calculate table values
	void fillTable(const int _nr_elem = 5000);

	// Value access
	RFLOAT operator()(RFLOAT val) const;

};

class TabCosine : public TabFunction
{
public:
	// Empty constructor
	TabCosine() {}

	void initialise(const int _nr_elem = 5000);

	//Pre-calculate table values
	void fillTable(const int _nr_elem = 5000);

	// Value access
	RFLOAT operator()(RFLOAT val) const;

};

class TabBlob : public TabFunction
{

private:
	RFLOAT radius;
	RFLOAT alpha;
	int order;

public:
	// Empty constructor
	TabBlob() {}

	// Constructor (with parameters)
	void initialise(RFLOAT _radius, RFLOAT _alpha, int _order, const int _nr_elem = 10000);

	//Pre-calculate table values
	void fillTable(const int _nr_elem = 5000);

	// Value access
	RFLOAT operator()(RFLOAT val) const;

};

class TabFtBlob : public TabFunction
{

private:
	RFLOAT radius;
	RFLOAT alpha;
	int order;

public:
	// Empty constructor
	TabFtBlob() {}

	 // Constructor (with parameters)
	void initialise(RFLOAT _radius, RFLOAT _alpha, int _order, const int _nr_elem = 10000);

	//Pre-calculate table values
	void fillTable(const int _nr_elem = 5000);

	// Value access
	RFLOAT operator()(RFLOAT val) const;

};


#endif /* TABFUNCS_H_ */
