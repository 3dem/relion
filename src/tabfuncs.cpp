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
#include "src/tabfuncs.h"

void TabSine::initialise(const int _nr_elem)
{
	sampling = 2 * PI / (RFLOAT) _nr_elem;
	TabSine::fillTable(_nr_elem);
}
//Pre-calculate table values
void TabSine::fillTable(const int _nr_elem)
{
	tabulatedValues.resize(_nr_elem);
	for (int i = 0; i < _nr_elem; i++)
	{
		RFLOAT xx = (RFLOAT) i * sampling;
		tabulatedValues(i) = sin(xx);
	}
}
// Value access
RFLOAT TabSine::operator()(RFLOAT val) const
{
	int idx = (int)( ABS(val) / sampling);
	RFLOAT retval = DIRECT_A1D_ELEM(tabulatedValues, idx % XSIZE(tabulatedValues));
	return (val < 0 ) ? -retval : retval;
}

void TabCosine::initialise(const int _nr_elem)
{
	sampling = 2 * PI / (RFLOAT) _nr_elem;
	TabCosine::fillTable(_nr_elem);
}
//Pre-calculate table values
void TabCosine::fillTable(const int _nr_elem)
{
	tabulatedValues.resize(_nr_elem);
	for (int i = 0; i < _nr_elem; i++)
	{
		RFLOAT xx = (RFLOAT) i * sampling;
		tabulatedValues(i) = cos(xx);
	}
}
// Value access
RFLOAT TabCosine::operator()(RFLOAT val) const
{
	int idx = (int)( ABS(val) / sampling);
	return DIRECT_A1D_ELEM(tabulatedValues, idx % XSIZE(tabulatedValues));
}

void TabBlob::initialise(RFLOAT _radius, RFLOAT _alpha, int _order, const int _nr_elem)
{
	radius = _radius;
	alpha = _alpha;
	order = _order;
	sampling = radius / _nr_elem;
	TabBlob::fillTable(_nr_elem);
}
//Pre-calculate table values
void TabBlob::fillTable(const int _nr_elem)
{
	tabulatedValues.resize(_nr_elem);
	for (int i = 0; i < _nr_elem; i++)
	{
		RFLOAT xx = (RFLOAT) i * sampling;
		tabulatedValues(i) = kaiser_value(xx, radius, alpha, order);
	}
}
// Value access
RFLOAT TabBlob::operator()(RFLOAT val) const
{
	int idx = (int)( ABS(val) / sampling);
	if (idx >= XSIZE(tabulatedValues))
		return 0.;
	else
		return DIRECT_A1D_ELEM(tabulatedValues, idx);
}

void TabFtBlob::initialise(RFLOAT _radius, RFLOAT _alpha, int _order, const int _nr_elem)
{
	radius = _radius;
	alpha = _alpha;
	order = _order;
	sampling = 0.5 / (RFLOAT)_nr_elem;
	TabFtBlob::fillTable(_nr_elem);
}
//Pre-calculate table values
void TabFtBlob::fillTable(const int _nr_elem)
{
	tabulatedValues.resize(_nr_elem);
	for (int i = 0; i < _nr_elem; i++)
	{
		RFLOAT xx = (RFLOAT) i * sampling;
		tabulatedValues(i) = kaiser_Fourier_value(xx, radius, alpha, order);
	}
}
// Value access
RFLOAT TabFtBlob::operator()(RFLOAT val) const
{
	int idx = (int)( ABS(val) / sampling);
	if (idx >= XSIZE(tabulatedValues))
		return 0.;
	else
		return DIRECT_A1D_ELEM(tabulatedValues, idx);
}

