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
/***************************************************************************
 *
 * Authors:     J.R. Bilbao-Castro (jrbcast@ace.ual.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
#include "src/metadata_container.h"

void MetaDataContainer::copy(const MetaDataContainer &MDc)
{

    clear();
	if (this != &MDc)
    {
		labels = MDc.labels;

        for (std::map<EMDLabel, RFLOAT *>::const_iterator It = MDc.RFLOATs.begin(); It != MDc.RFLOATs.end(); It++)
            addValue(It->first, *It->second);
        for (std::map<EMDLabel, int *>::const_iterator It = MDc.ints.begin(); It != MDc.ints.end(); It++)
            addValue(It->first, *It->second);
        for (std::map<EMDLabel, long *>::const_iterator It = MDc.longs.begin(); It != MDc.longs.end(); It++)
            addValue(It->first, *It->second);
        for (std::map<EMDLabel, bool *>::const_iterator It = MDc.bools.begin(); It != MDc.bools.end(); It++)
            addValue(It->first, *It->second);
        for (std::map<EMDLabel, std::string *>::const_iterator It = MDc.strings.begin(); It != MDc.strings.end(); It++)
            addValue(It->first, *It->second);
    }
}

void MetaDataContainer::copy_select(const MetaDataContainer &MDc, std::vector<EMDLabel> only_labels)
{

    clear();
	if (this != &MDc)
    {
		labels = only_labels;

        for (std::map<EMDLabel, RFLOAT *>::const_iterator It = MDc.RFLOATs.begin(); It != MDc.RFLOATs.end(); It++)
        {
        	if (vectorContainsLabel(only_labels, It->first))
        		addValue(It->first, *It->second);
        }
        for (std::map<EMDLabel, int *>::const_iterator It = MDc.ints.begin(); It != MDc.ints.end(); It++)
        {
        	if (vectorContainsLabel(only_labels, It->first))
        		addValue(It->first, *It->second);
        }
        for (std::map<EMDLabel, long *>::const_iterator It = MDc.longs.begin(); It != MDc.longs.end(); It++)
        {
        	if (vectorContainsLabel(only_labels, It->first))
        		addValue(It->first, *It->second);
        }
        for (std::map<EMDLabel, bool *>::const_iterator It = MDc.bools.begin(); It != MDc.bools.end(); It++)
        {
        	if (vectorContainsLabel(only_labels, It->first))
        		addValue(It->first, *It->second);
        }
        for (std::map<EMDLabel, std::string *>::const_iterator It = MDc.strings.begin(); It != MDc.strings.end(); It++)
        {
        	if (vectorContainsLabel(only_labels, It->first))
        		addValue(It->first, *It->second);
        }
    }
}

MetaDataContainer& MetaDataContainer::operator =(const MetaDataContainer &MDc)
{
    copy(MDc);
    return *this;
}

MetaDataContainer::MetaDataContainer(const MetaDataContainer &MDc)
{
    copy(MDc);
}


void MetaDataContainer::addValueFromString(const EMDLabel &lCode, const std::string &value)
{

	if (EMDL::isString(lCode))
	{
		addValue(lCode, value);
	}
	else
	{
		std::istringstream i(value);
		// Look for a RFLOAT value
		if (EMDL::isDouble(lCode))
		{
			RFLOAT RFLOATValue;
			i >> RFLOATValue;
			addValue(lCode, RFLOATValue);
		}
		else if (EMDL::isInt(lCode))
		{
			int intValue;
			i >> intValue;
			addValue(lCode, intValue);
		}
		else if (EMDL::isLong(lCode))
		{
			long int longValue;
			i >> longValue;
			addValue(lCode, longValue);
		}
		else if (EMDL::isBool(lCode))
		{
			bool boolValue;
			i >> boolValue;
			addValue(lCode, boolValue);
		}
	}
}

#ifdef RELION_SINGLE_PRECISION
void MetaDataContainer::addValue(EMDLabel name, const double &value)
{
    addValue(name, (RFLOAT)value);
}
#endif

void MetaDataContainer::addValue(EMDLabel name, const RFLOAT &value)
{
	if (! valueExists(name))
		labels.push_back(name);

    if (EMDL::isDouble(name))
    {
		if (RFLOATs[name])
			delete RFLOATs[name];
		RFLOATs[name] = new RFLOAT(value);
    }
    else
    	REPORT_ERROR("addValue for RFLOAT: label " + EMDL::label2Str(name) + " is not of type RFLOAT!");
}

void MetaDataContainer::addValue(EMDLabel name, const int &value)
{
	if (! valueExists(name))
		labels.push_back(name);

    if (EMDL::isInt(name))
    {
		if (ints[name])
			delete ints[name];
		ints[name] = new int(value);
    }
    else
    	REPORT_ERROR("addValue for int: label " + EMDL::label2Str(name) + " is not of type int!");
}

void MetaDataContainer::addValue(EMDLabel name, const long int &value)
{
	if (! valueExists(name))
		labels.push_back(name);

    if (EMDL::isLong(name))
    {
		if (longs[name])
			delete longs[name];
		longs[name] = new long(value);
    }
    else
    	REPORT_ERROR("addValue for long: label " + EMDL::label2Str(name) + " is not of type long!");
}

void MetaDataContainer::addValue(EMDLabel name, const bool value)
{
	if (! valueExists(name))
		labels.push_back(name);

    if (EMDL::isBool(name))
    {
		if (bools[name])
			delete bools[name];
		bools[name] = new bool(value);
    }
    else
    	REPORT_ERROR("addValue for bool: label " + EMDL::label2Str(name) + " is not of type bool!");
}

void MetaDataContainer::addValue(EMDLabel name, const std::string &value)
{
	if (! valueExists(name))
		labels.push_back(name);

	if (EMDL::isString(name))
    {
		if (strings[name])
			delete strings[name];
		strings[name] = new std::string(value);
    }
    else
    	REPORT_ERROR("addValue for string: label " + EMDL::label2Str(name) + " is not of type string!");
}

/** Creates a new label-value pair, with the default value for the corresponding type */
void MetaDataContainer::addDefaultValue(EMDLabel name)
{
	if (! valueExists(name))
		labels.push_back(name);

	if (EMDL::isDouble(name))
    {
		if (RFLOATs[name])
			delete RFLOATs[name];
		RFLOATs[name] = new RFLOAT(0.);
    }
    else if (EMDL::isInt(name))
    {
		if (ints[name])
			delete ints[name];
		ints[name] = new int(0.);
    }
    else if (EMDL::isLong(name))
    {
		if (longs[name])
			delete longs[name];
		longs[name] = new long(0.);
    }
    else if (EMDL::isBool(name))
    {
		if (bools[name])
			delete bools[name];
		bools[name] = new bool(false);
    }
    else if (EMDL::isString(name))
    {
		if (strings[name])
			delete strings[name];
		strings[name] = new std::string("");
    }
    else
    	REPORT_ERROR("MetaDataContainer::addDefaultValu: unrecognised data type for label " + EMDL::label2Str(name));
}

bool MetaDataContainer::getValue( const EMDLabel name, RFLOAT &value)
{
    std::map<EMDLabel, RFLOAT *>::iterator element = RFLOATs.find(name);
    if (element == RFLOATs.end())
        return false;
    else
    	value = *element->second;
    return true;
}

bool MetaDataContainer::getValue( const EMDLabel name, int &value)
{
    std::map<EMDLabel, int *>::iterator element = ints.find(name);
    if (element == ints.end())
        return false;
    else
    	value = *element->second;
    return true;
}

bool MetaDataContainer::getValue( const EMDLabel name, long &value)
{
    std::map<EMDLabel, long *>::iterator element = longs.find(name);
    if (element == longs.end())
        return false;
    else
    	value = *element->second;
    return true;
}

bool MetaDataContainer::getValue( const EMDLabel name, bool &value)
{
    std::map<EMDLabel, bool *>::iterator element = bools.find(name);
    if (element == bools.end())
        return false;
    else
    	value = *element->second;
    return true;
}

bool MetaDataContainer::getValue( const EMDLabel name, std::string &value)
{
    std::map<EMDLabel, std::string *>::iterator element = strings.find(name);
    if (element == strings.end())
        return false;
    else
    	value = *element->second;
    return true;
}

bool MetaDataContainer::valueExists(EMDLabel name)
{
	if (! EMDL::isValidLabel(name))
		REPORT_ERROR("Unrecognised label type in MetaDataContainer valueExists");

	if (find(labels.begin(), labels.end(), name) != labels.end())
		return true;
    return false;
}

std::vector<EMDLabel> MetaDataContainer::getLabels()
{
	return labels;
}

bool MetaDataContainer::writeValueToStream(std::ostream &outstream, EMDLabel inputLabel)
{
	if (valueExists(inputLabel))
    {
		if (EMDL::isDouble(inputLabel))
        {
            RFLOAT d;
            if (! getValue(inputLabel, d))
        		REPORT_ERROR("Double value not found in writeValueToStream.");
            if ((ABS(d) > 0. && ABS(d) < 0.001) || ABS(d) > 100000.)
                outstream << std::setw(12) << std::scientific;
            else
                outstream << std::setw(12) << std::fixed;
            outstream << d;
        }
        else if (EMDL::isString(inputLabel))
        {
        	std::string s;
        	if (! getValue(inputLabel, s))
        		REPORT_ERROR("String value not found in writeValueToStream.");
            outstream << s;
        }
        else if (EMDL::isInt(inputLabel))
        {
        	int i;
        	if (! getValue(inputLabel, i))
    			REPORT_ERROR("Integer value not found in writeValueToStream.");
        	outstream << std::setw(12) << std::fixed;
            outstream << i;
        }
        else if (EMDL::isLong(inputLabel))
        {
        	long l;
        	if (! getValue(inputLabel, l))
        		REPORT_ERROR("Long value not found in writeValueToStream.");
        	outstream << std::setw(12) << std::fixed;
            outstream << l;
        }
        else if (EMDL::isBool(inputLabel))
        {
        	bool b;
        	if (! getValue(inputLabel, b))
        		REPORT_ERROR("Boolean value not found in writeValueToStream.");
        	outstream << std::setw(12) << std::fixed;
        	outstream << b;
        }
		return true;
    }
    else
    {
        return false;
    }
}

bool MetaDataContainer::writeValueToString(std::string &outString, EMDLabel inLabel)
{
    std::ostringstream oss;
    bool result = writeValueToStream(oss, inLabel);
    outString = result ? oss.str() : std::string("");

    return result;
}

void MetaDataContainer::keepOnlyLabels(std::vector<EMDLabel> only_labels)
{
	MetaDataContainer MDc = *this;
	copy_select(MDc, only_labels);
}
