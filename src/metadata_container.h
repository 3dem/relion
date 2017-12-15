/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
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

#ifndef METADATA_CONTAINER_H
#define METADATA_CONTAINER_H

#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "src/funcs.h"
#include "src/metadata_label.h"

class MetaDataTable;

class MetaDataContainer
{    
    public:

        MetaDataTable* table;

        std::vector<RFLOAT> RFLOATs;
        std::vector<int> ints;
        std::vector<long> longs;
        std::vector<bool> bools;
        std::vector<std::string> strings;

        MetaDataContainer();
        MetaDataContainer(MetaDataTable* table, long RFLOATCount, long intCount,
                          long longCount, long boolCount, long stringCount);
        MetaDataContainer(MetaDataTable* table, MetaDataContainer* mdc);

        void getValue(long offset, RFLOAT& dest) const;
        void getValue(long offset, int& dest) const;
        void getValue(long offset, long& dest) const;
        void getValue(long offset, bool& dest) const;
        void getValue(long offset, std::string& dest) const;

	// Even when RELION_SINGLE_PRECISION, this must be double
	// as floating point literals are double
        void setValue(long offset, const double& src);
        void setValue(long offset, const int& src);
        void setValue(long offset, const long& src);
        void setValue(long offset, const bool& src);
        void setValue(long offset, const std::string& src);
};

#endif
