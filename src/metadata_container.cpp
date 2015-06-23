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

void MetaDataContainer::insertVoidPtr(EMDLabel name, void * value)
{
//     If values[name] already existed, first free that memory
    if (values[name])
    {
    	if (EMDL::isDouble(name))
            delete (double*)values[name];
        else if (EMDL::isInt(name))
            delete (int*)values[name];
        else if (EMDL::isLong(name))
            delete (long int*)values[name];
        else if (EMDL::isBool(name))
            delete (bool*)values[name];
        else if (EMDL::isString(name))
        	delete (std::string*)values[name];
        else
            REPORT_ERROR("Unrecognised label type in MetaDataContainer clear");
    }
    values[name] = value;
}

void * MetaDataContainer::getVoidPtr(EMDLabel name)
{
    std::map<EMDLabel, void *>::iterator element;
    element = values.find(name);
    if (element == values.end())
    {
        REPORT_ERROR((std::string) "Label " + EMDL::label2Str(name) + " not found on getVoidPtr()\n" );
    }
    else
    {
        return element->second;
    }
}

void MetaDataContainer::copy(const MetaDataContainer &MDc)
{

    clear();
	if (this != &MDc)
    {
        void * aux;
        EMDLabel lCode;
        std::map<EMDLabel, void *>::const_iterator It;
        for (It = (MDc.values).begin(); It != (MDc.values).end(); It++)
        {
            aux = It->second;
            lCode = It->first;

            if (EMDL::isDouble(lCode))
            {
                addValue(lCode, *((double *) aux));
            }
            else if (EMDL::isString(lCode))
            {
                addValue(lCode, *((std::string *) aux));
            }
            else if (EMDL::isInt(lCode))
            {
                addValue(lCode, *((int *) aux));
            }
            else if (EMDL::isLong(lCode))
            {
                addValue(lCode, *((long int *) aux));
            }
            else if (EMDL::isBool(lCode))
            {
                addValue(lCode, *((bool *) aux));
            }

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
		// Look for a double value
		if (EMDL::isDouble(lCode))
		{
			double doubleValue;
			i >> doubleValue;
			addValue(lCode, doubleValue);
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

void MetaDataContainer::addValue(EMDLabel name, const double &value)
{
    if (EMDL::isDouble(name))
    {
    	void * newValue = (void *) (new double(value));
    	insertVoidPtr(name, newValue);
    }
    else
    	REPORT_ERROR("addValue for double: label " + EMDL::label2Str(name) + " is not of type double!");
}

void MetaDataContainer::addValue(EMDLabel name, const int &value)
{
    if (EMDL::isInt(name))
    {
    	void * newValue = (void *) (new int(value));
    	insertVoidPtr(name, newValue);
    }
    else
    	REPORT_ERROR("addValue for int: label " + EMDL::label2Str(name) + " is not of type int!");
}

void MetaDataContainer::addValue(EMDLabel name, const long int &value)
{
    if (EMDL::isLong(name))
    {
    	void * newValue = (void *) (new long int(value));
    	insertVoidPtr(name, newValue);
    }
    else
    	REPORT_ERROR("addValue for long: label " + EMDL::label2Str(name) + " is not of type long!");
}

void MetaDataContainer::addValue(EMDLabel name, const bool &value)
{
    if (EMDL::isBool(name))
    {
    	void * newValue = (void *) (new bool(value));
    	insertVoidPtr(name, newValue);
    }
    else
    	REPORT_ERROR("addValue for bool: label " + EMDL::label2Str(name) + " is not of type bool!");
}

void MetaDataContainer::addValue(EMDLabel name, const std::string &value)
{
    if (EMDL::isString(name))
    {
    	void * newValue = (void *) (new std::string(value));
    	insertVoidPtr(name, newValue);
    }
    else
    	REPORT_ERROR("addValue for string: label " + EMDL::label2Str(name) + " is not of type string!");
}

/** Creates a new label-value pair, with the default value for the corresponding type */
void MetaDataContainer::addDefaultValue(EMDLabel name)
{
	void * newValue;
	if (EMDL::isDouble(name))
    {
    	newValue = (void *) (new double(0.));
    }
    else if (EMDL::isInt(name) || EMDL::isLong(name))
    {
    	newValue = (void *) (new int(0));
    }
    else if (EMDL::isBool(name))
    {
    	newValue = (void *) (new bool(false));
    }
    else if (EMDL::isString(name))
    {
    	newValue = (void *) (new std::string(""));
    }
    else
    	REPORT_ERROR("MetaDataContainer::addDefaultValu: unrecognised data type for label " + EMDL::label2Str(name));

    insertVoidPtr(name, newValue);

}

bool MetaDataContainer::getValue( const EMDLabel name, double &value)
{
    std::map<EMDLabel, void *>::iterator element;

    element = values.find(name);

    if (element == values.end())
    {
        return false;
    }
    else
    {
    	if (EMDL::isDouble(element->first))
    		value = *((double *) element->second);
    	else
    		REPORT_ERROR("getValue for double: label " + EMDL::label2Str(element->first) + " is not of type double!");

    	return true;
    }
}

bool MetaDataContainer::getValue( const EMDLabel name, int &value)
{
    std::map<EMDLabel, void *>::iterator element;

    element = values.find(name);

    if (element == values.end())
    {
        return false;
    }
    else
    {
    	if (EMDL::isInt(element->first))
    		value = *((int *) element->second);
    	else
    		REPORT_ERROR("getValue for int: label " + EMDL::label2Str(element->first) + " is not of type int!");

    	return true;
    }
}

bool MetaDataContainer::getValue( const EMDLabel name, long int &value)
{
    std::map<EMDLabel, void *>::iterator element;

    element = values.find(name);

    if (element == values.end())
    {
        return false;
    }
    else
    {
    	if (EMDL::isLong(element->first))
    		value = *((long int *) element->second);
    	else
    		REPORT_ERROR("getValue for long int: label " + EMDL::label2Str(element->first) + " is not of type long int!");

    	return true;
    }
}

bool MetaDataContainer::getValue( const EMDLabel name, bool &value)
{
    std::map<EMDLabel, void *>::iterator element;

    element = values.find(name);

    if (element == values.end())
    {
        return false;
    }
    else
    {
    	if (EMDL::isBool(element->first))
    		value = *((bool *) element->second);
    	else
    		REPORT_ERROR("getValue for bool: label " + EMDL::label2Str(element->first) + " is not of type bool!");

    	return true;
    }
}

bool MetaDataContainer::getValue( const EMDLabel name, std::string &value)
{
    std::map<EMDLabel, void *>::iterator element;

    element = values.find(name);

    if (element == values.end())
    {
        return false;
    }
    else
    {
    	if (EMDL::isString(element->first))
    		value = *((std::string *) element->second);
    	else
    		REPORT_ERROR("getValue for string: label " + EMDL::label2Str(element->first) + " is not of type string!");

    	return true;
    }
}

bool MetaDataContainer::valueExists(EMDLabel name)
{
    if (values.find(name) == values.end())
    {
        return false;
    }
    else
    {
        return true;
    }
}
//A template exists for pairexists different from string
bool MetaDataContainer::pairExists(EMDLabel name, const std::string &value)
{
    // Traverse all the structure looking for objects
    // that satisfy search criteria
    std::map<EMDLabel, void *>::iterator It;

    It = values.find(name);

    if (It != values.end())
    {
        if (*((std::string *) (It->second)) == value)
        {
            return true;
        }
    }

    return false;
}

std::vector<EMDLabel> MetaDataContainer::getLabels()
{
	std::vector<EMDLabel> result;
	std::map<EMDLabel, void *>::iterator It;
	for (It = values.begin(); It != values.end(); It++)
		result.push_back(It->first);

	return result;
}

bool MetaDataContainer::writeValueToStream(std::ostream &outstream,
        EMDLabel inputLabel)
{
	if (valueExists(inputLabel))
    {
#ifdef DEBUG_MDC
		std::cerr << " EMDL::label2Str(inputLabel)= " << EMDL::label2Str(inputLabel) << std::endl;
#endif
		if (EMDL::isDouble(inputLabel))
        {
            double d;
            d = *((double*) (getVoidPtr(inputLabel)));
            if ((ABS(d) > 0. && ABS(d) < 0.001) || ABS(d) > 100000.)
                outstream << std::setw(12) << std::scientific;
            else
                outstream << std::setw(12) << std::fixed;
            outstream << d;
#ifdef DEBUG_MDC
            std::cerr << " d= " << d << std::endl;
#endif
        }
        else if (EMDL::isString(inputLabel))
        {
            outstream << *((std::string*) (getVoidPtr(inputLabel)));
#ifdef DEBUG_MDC
            std::cerr << " *((std::string*) (getVoidPtr(inputLabel)))= " << *((std::string*) (getVoidPtr(inputLabel))) << std::endl;
#endif
        }
        else if (EMDL::isInt(inputLabel))
        {
        	outstream << std::setw(12) << std::fixed;
            outstream << *((int*) (getVoidPtr(inputLabel)));
#ifdef DEBUG_MDC
            std::cerr << " *((int*) (getVoidPtr(inputLabel)))= " << *((int*) (getVoidPtr(inputLabel))) << std::endl;
#endif
        }
        else if (EMDL::isLong(inputLabel))
        {
        	outstream << std::setw(12) << std::fixed;
            outstream << *((long int*) (getVoidPtr(inputLabel)));
#ifdef DEBUG_MDC
            std::cerr << " *((long int*) (getVoidPtr(inputLabel)))= " << *((long int*) (getVoidPtr(inputLabel))) << std::endl;
#endif
        }
        else if (EMDL::isBool(inputLabel))
        {
        	outstream << std::setw(12) << std::fixed;
        	outstream << *((bool*) (getVoidPtr(inputLabel)));
#ifdef DEBUG_MDC
        	std::cerr << " *((bool*) (getVoidPtr(inputLabel)))= " << *((bool*) (getVoidPtr(inputLabel))) << std::endl;
#endif
        }
		return true;
    }
    else
    {
        return false;
    }
}

bool MetaDataContainer::writeValueToString(std::string &outString,
        EMDLabel inLabel)
{
    std::ostringstream oss;
    bool result = writeValueToStream(oss, inLabel);
    outString = result ? oss.str() : std::string("");

    return result;
}

