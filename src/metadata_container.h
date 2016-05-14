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
//static RFLOAT zeroD=0.;
//static RFLOAT    oneD=1.;
//static bool  falseb=false;

class MetaDataContainer
{
	//Labels is just to keep the order of things
	std::vector<EMDLabel> labels;

    std::map<EMDLabel, RFLOAT *> RFLOATs;
    std::map<EMDLabel, int *> ints;
    std::map<EMDLabel, long *> longs;
    std::map<EMDLabel, bool *> bools;
    std::map<EMDLabel, std::string *> strings;

    void copy(const MetaDataContainer &MDc);
    // Copy only the provided labels
    void copy_select(const MetaDataContainer &MDc, std::vector<EMDLabel> only_labels);

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
#ifdef RELION_SINGLE_PRECISION
    void addValue(EMDLabel name, const double &value);
#endif
    void addValue(EMDLabel name, const RFLOAT &value);
    void addValue(EMDLabel name, const int &value);
    void addValue(EMDLabel name, const long int &value);
    void addValue(EMDLabel name, const bool value);
    void addValue(EMDLabel name, const std::string &value);

    /** Creates a new label-value pair, with the default value for the corresponding type */
    void addDefaultValue(EMDLabel name);

    /** clean metadatacontainer
     *
     */
    void clear(void)
    {
    	for (std::map<EMDLabel, RFLOAT *>::iterator it = RFLOATs.begin(); it != RFLOATs.end(); ++it)
    		delete it->second;
    	RFLOATs.clear();
    	for (std::map<EMDLabel, int *>::iterator it = ints.begin(); it != ints.end(); ++it)
    		delete it->second;
    	ints.clear();
    	for (std::map<EMDLabel, long *>::iterator it = longs.begin(); it != longs.end(); ++it)
    		delete it->second;
    	longs.clear();
    	for (std::map<EMDLabel, bool *>::iterator it = bools.begin(); it != bools.end(); ++it)
    		delete it->second;
    	bools.clear();
    	for (std::map<EMDLabel, std::string *>::iterator it = strings.begin(); it != strings.end(); ++it)
    		delete it->second;
    	strings.clear();
    	labels.clear();
    }

    /** Get a value for a given name.
     *  If the metadata container contains the name,
     *  the function will check the type of value with the type of the label and report an error if they do not match,
     *  If they do match, the value will be get and the function returns true
     *  If the name does not exist in the container, the function returns false
     *
     */
    bool getValue( const EMDLabel name, RFLOAT &value);
    bool getValue( const EMDLabel name, int &value);
    bool getValue( const EMDLabel name, long int &value);
    bool getValue( const EMDLabel name, bool &value);
    bool getValue( const EMDLabel name, std::string &value);

    // Check whether this label is present in the container
    bool valueExists(EMDLabel name);

    bool pairExists(EMDLabel name, const RFLOAT value)
    {
        std::map<EMDLabel, RFLOAT *>::iterator It = RFLOATs.find(name);
        if (It != RFLOATs.end())
            if (ABS( *It->second - value ) < XMIPP_EQUAL_ACCURACY)
                return true;
        return false;
    }

    bool pairExists(EMDLabel name, const int value)
    {
        std::map<EMDLabel, int *>::iterator It = ints.find(name);
        if (It != ints.end())
            if (*It->second == value)
                return true;
        return false;
    }

    bool pairExists(EMDLabel name, const long value)
    {
        std::map<EMDLabel, long *>::iterator It = longs.find(name);
        if (It != longs.end())
            if (*It->second == value)
                return true;
        return false;
    }

    bool pairExists(EMDLabel name, const bool value)
    {
        std::map<EMDLabel, bool *>::iterator It = bools.find(name);
        if (It != bools.end())
            if (*It->second == value)
                return true;
        return false;
    }

    bool pairExists(EMDLabel name, const std::string &value)
    {
        std::map<EMDLabel, std::string *>::iterator It = strings.find(name);
        if (It != strings.end())
            if (*It->second == value)
                return true;
        return false;
    }

    // Return a container with only the selected containers
    void keepOnlyLabels(std::vector<EMDLabel> only_labels);

    // Get all labels in this container
    std::vector<EMDLabel> getLabels();

    bool writeValueToStream(std::ostream &outstream, EMDLabel inputLabel);
    bool writeValueToString(std::string &outString, EMDLabel inputLabel);
};

#endif
