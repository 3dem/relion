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
#include "src/metadata_container.h"

MetaDataContainer::MetaDataContainer()
    :   RFLOATs(0), ints(0), longs(0), bools(0), strings(0)
{}


MetaDataContainer::MetaDataContainer(
        MetaDataTable *table, long RFLOATCount, long intCount,
        long longCount, long boolCount, long stringCount)
: table(table),
  RFLOATs(RFLOATCount, 0),
  ints(intCount, 0),
  longs(longCount, 0),
  bools(boolCount, false),
  strings(stringCount, "")
{}

MetaDataContainer::MetaDataContainer(
        MetaDataTable *table, MetaDataContainer* mdc)
: table(table),
  RFLOATs(mdc->RFLOATs),
  ints(mdc->ints),
  longs(mdc->longs),
  bools(mdc->bools),
  strings(mdc->strings)
{}

void MetaDataContainer::getValue(long offset, RFLOAT& dest) const
{
    dest = RFLOATs[offset];
}

void MetaDataContainer::getValue(long offset, int& dest) const
{
    dest = ints[offset];
}

void MetaDataContainer::getValue(long offset, long& dest) const
{
    dest = longs[offset];
}

void MetaDataContainer::getValue(long offset, bool& dest) const
{
    dest = bools[offset];
}

void MetaDataContainer::getValue(long offset, std::string& dest) const
{
    dest = strings[offset];
}


void MetaDataContainer::setValue(long offset, const RFLOAT& src)
{
    RFLOATs[offset] = src;
}

void MetaDataContainer::setValue(long offset, const int& src)
{
    ints[offset] = src;
}

void MetaDataContainer::setValue(long offset, const long& src)
{
    longs[offset] = src;
}

void MetaDataContainer::setValue(long offset, const bool& src)
{
    bools[offset] = src;
}

void MetaDataContainer::setValue(long offset, const std::string& src)
{
    strings[offset] = src;
}
