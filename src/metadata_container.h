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

#ifndef METADATA_CONTAINER_H
#define METADATA_CONTAINER_H

#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "src/funcs.h"
//#include "src/xmipp/strings.h"
#include "src/metadata_label.h"

//useful to init values to zero
//static double zeroD=0.;
//static double    oneD=1.;
//static bool  falseb=false;

class MetaDataContainer
{
    /** Container for attribute-value pairs.
     * Note that void * allows to use mixed types */
    std::map<EMDLabel, void *> values;

    void insertVoidPtr(EMDLabel name, void * value);
    void * getVoidPtr(EMDLabel name);
    void copy(const MetaDataContainer &MDc);

public:

    /**Assignment operator
     *
     */
    MetaDataContainer& operator =(const MetaDataContainer &MDc);

    /** Empty Constructor */
    MetaDataContainer() {};

    /** Copy constructor */
    MetaDataContainer(const MetaDataContainer &MDc);

    /** Destructor */
    ~MetaDataContainer()
    {
    	clear();
    }

    /** Create a new attribute-value pair of string type */
    void addValueFromString(const EMDLabel &label, const std::string &value);

    /** Creates a new label-value pair, and checks the type of the label is the same as that of value */
    void addValue(EMDLabel name, const double &value);
    void addValue(EMDLabel name, const int &value);
    void addValue(EMDLabel name, const long int &value);
    void addValue(EMDLabel name, const bool &value);
    void addValue(EMDLabel name, const std::string &value);

    /** Creates a new label-value pair, with the default value for the corresponding type */
    void addDefaultValue(EMDLabel name);

    /** clean metadatacontainer
     *
     */
    void clear(void)
    {
        // Manually delete the map of pointers!
        for (std::map<EMDLabel, void *>::iterator it = values.begin();
             it != values.end(); ++it)
        {
        	if (EMDL::isDouble(it->first))
                delete (double*)it->second;
            else if (EMDL::isInt(it->first))
                delete (int*)it->second;
            else if (EMDL::isLong(it->first))
                delete (long int*)it->second;
            else if (EMDL::isBool(it->first))
                delete (bool*)it->second;
            else if (EMDL::isString(it->first))
                delete (std::string*)it->second;
            else
                REPORT_ERROR("Unrecognised label type in MetaDataContainer clear");
        }
    	values.clear();
    }

    /** Get a value for a given name.
     *  If the metadata container contains the name,
     *  the function will check the type of value with the type of the label and report an error if they do not match,
     *  If they do match, the value will be get and the function returns true
     *  If the name does not exist in the container, the function returns false
     *
     */
    bool getValue( const EMDLabel name, double &value);
    bool getValue( const EMDLabel name, int &value);
    bool getValue( const EMDLabel name, long int &value);
    bool getValue( const EMDLabel name, bool &value);
    bool getValue( const EMDLabel name, std::string &value);

    // Check whether this label is present in the container
    bool valueExists(EMDLabel name);

    //string is not part of the template because - is not defined for string
    bool pairExists(EMDLabel name, const std::string &value);

    template<class T>
    bool pairExists(EMDLabel name, const T& value)
    {
        // Traverse all the structure looking for objects
        // that satisfy search criteria
        std::map<EMDLabel, void *>::iterator It;

        It = values.find(name);

        if (It != values.end())
        {
            if (ABS( *((T *)(It->second)) - value )
                    < XMIPP_EQUAL_ACCURACY)
            {
                return true;
            }
        }

        return false;
    }

    // Get all labels in this container
    std::vector<EMDLabel> getLabels();

    bool writeValueToStream(std::ostream &outstream, EMDLabel inputLabel);
    bool writeValueToString(std::string &outString, EMDLabel inputLabel);
};

#endif
