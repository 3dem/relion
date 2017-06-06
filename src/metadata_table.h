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

#ifndef METADATA_TABLE_H
#define METADATA_TABLE_H

#include <map>
#include <vector>
#include <iostream>
#include <iterator>
#include <sstream>
#include <time.h>
#include <stdio.h>
#include <sstream>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include "src/funcs.h"
#include "src/args.h"
#include "src/CPlot2D.h"
#include "src/metadata_container.h"

/** For all objects.
 @code
 FOR_ALL_OBJECTS_IN_METADATA(metadata) {
   RFLOAT rot;
   DF.getValue( EMDL_ANGLEROT, rot);
 }
 @endcode
 */
#define FOR_ALL_OBJECTS_IN_METADATA_TABLE(kkkk_metadata) \
        for(long int current_object = (kkkk_metadata).firstObject(); \
             current_object != MetaDataTable::NO_MORE_OBJECTS && current_object!= MetaDataTable::NO_OBJECTS_STORED; \
             current_object=(kkkk_metadata).nextObject())

/** MetaDataTable Manager.
 *
 */
class MetaDataTable
{
    // Effectively stores all metadata
    std::vector<MetaDataContainer *> objects;

    // Current object id
    long int current_objectID;

    // Is this a 2D table or a 1D list?
    bool isList;

    // Name of the metadata table
    std::string name;

    // A comment for the metadata table
    std::string comment;

public:

    /** What labels have been read from a docfile/metadata file
     *   and/or will be stored on a new metadata file when "save" is
     *   called
     **/
    std::vector<EMDLabel> activeLabels;

    /** When reading a column formated file, if a label is found that
     *   does not exists as a EMDLabel, it is ignored. For further
     *   file processing, such columns must be ignored and this structure
     *   allows to do that
     **/
    std::vector<unsigned int> ignoreLabels;

    /** Empty Constructor.
     *
     * The MetaDataTable is created with no data stored on it.
     */
    MetaDataTable();

    /** Copy constructor
     *
     * Created a new MetaDataTable by copying all data from an existing MetaDataTable object.
     */
    MetaDataTable(const MetaDataTable & c);


    /** Assignment operator
     * @ingroup MetaDataConstructors
     *
     * Copies MetaDataTable from an existing MetaDataTable object.
     */
    MetaDataTable& operator =(const MetaDataTable &MD);

    /** Destructor
     * @ingroup MetaDataTableConstructors
     *
     * Frees all used memory and destroys object.
     */
    ~MetaDataTable();

    /** Set to false if this is a MetaDataTable
     *  set to true  if this is a simple list of parameters
     *
     */
    void setIsList(bool is_list);

    // returns true of the metadatatable contains no objects
    bool isEmpty() const;

    // Returns the number of elements in the table
    long int numberOfObjects() const;

    // Removes all objects from the metadatatable
    void clear();

    /**Set Comment
     */
    void setComment(const std::string Comment = "");

    /**Get Comment
     */
    std::string getComment() const;

    /** Does this table contain a comment?
     */
    bool containsComment() const;

    /**Set Name
     */
    void setName(const std::string Name = "");

    /**Get Name
     */
    std::string getName() const;

    size_t size(void)
    {
        return objects.size();
    }

    /*  Get value for any label.
     *  If the label does not exist, false is returned
     *  Otherwise, value is set correctly and true is returned
     */

    template<class T>
    bool getValue(EMDLabel name, T &value,
                  long int objectID = -1) const
    {
        if (isEmpty())
        	return false;
        MetaDataContainer * aux = getObject(objectID);

		// Inside getValue of the container there will be a check of the correct type
        bool result = aux->getValue(name, value);
        return result;
    }

    // Read/set a new pair/value for an specified object. If no objectID is given, that
    // pointed by the class iterator is used
    bool setValueFromString(const EMDLabel &label, const std::string &value, long int objectID = -1);

    template<class T>
    bool setValue(EMDLabel name, const T &value, long int objectID=-1)
    {
        if (!isEmpty() && EMDL::isValidLabel(name))
        {

            long int auxID = (objectID == -1) ? current_objectID : objectID;
            MetaDataContainer * aux = objects[auxID];

            // Check whether label is correct (belongs to the enum in the metadata_container header
            // and whether it is present in the activeLabels vector. If not, add it to all the other
            // objects with default values
            std::vector<EMDLabel>::iterator location;
            std::map<long int, MetaDataContainer *>::iterator It;
            location = std::find(activeLabels.begin(), activeLabels.end(), name);
            if (location == activeLabels.end())
            {
                activeLabels.push_back(name);
                // Add this label to the rest of the objects in this class
                for (long int idx = 0; idx < objects.size(); idx++)
                {
                    if (objects[idx] != aux)
                    {
                    	objects[idx]->addValue(name, T());
                    }
                }
            }
            aux->addValue(name, value);
            return true;
        }
        else
        {
            std::cerr << " %%WARNING%%: unsuccessful setValue" << std::endl;
            std::cerr << " isEmpty()= " << isEmpty() << std::endl;
        	std::cerr << " EMDL::isValidLabel(name)= " << EMDL::isValidLabel(name) << std::endl;
            std::cerr << " name= " << EMDL::label2Str(name) << " value= " << value <<std::endl;
            return false;
        }
    }

    // No copying of entire MetaDataTable involved!
    void newSort(const EMDLabel name, bool do_reverse = false, bool do_sort_after_at = false);

    // Sort the order of the elements based on the values in the input label (only numbers, no strings/bools!)
    // Have to pass dummy T parameter to get the template thing running?
    void sort(EMDLabel name, bool do_reverse = false, bool only_set_index = false, bool do_random = false)
    {

    	if (do_random)
    		srand (time(NULL));    		  /* initialize random seed: */
    	else if (!(EMDL::isInt(name) || EMDL::isLong(name) || EMDL::isDouble(name)) )
    		REPORT_ERROR("MetadataTable::sort%% ERROR: can only sorted numbers");

    	std::vector<std::pair<RFLOAT,long int> > vp;
    	vp.reserve(objects.size());
    	long int i = 0;
    	FOR_ALL_OBJECTS_IN_METADATA_TABLE(*this)
    	{
    		RFLOAT dval;
    		if (do_random)
    		{
    			dval = (RFLOAT)rand();
    		}
    		else if (EMDL::isInt(name))
    		{
    			int val;
    			getValue(name, val);
    			dval = (RFLOAT) val;
    		}
    		else if (EMDL::isLong(name))
    		{
    			long int val;
    			getValue(name, val);
    			dval = (RFLOAT) val;
    		}
    		else
    		{
    			//  EMDL::isDouble(name)
    			getValue(name, dval);
    		}
    		vp.push_back(std::make_pair(dval, i));
    		i++;
    	}

    	std::sort(vp.begin(), vp.end());
    	if (do_reverse && !do_random)
    		std::reverse(vp.begin(), vp.end());

    	if (only_set_index)
    	{
    		// Add an extra column with the sorted position of each entry
			for (long int j = 0; j < vp.size(); j++)
			{
				(*this).setValue(EMDL_SORTED_IDX, j, vp[j].second);
			}
    	}
    	else
    	{
			// Change the actual order in the MetaDataTable
    		MetaDataTable MDaux;
			for (long int j = 0; j < vp.size(); j++)
			{
				MDaux.addObject();
				MDaux.setObject((*this).getObject(vp[j].second));
			}
			*this = MDaux;
    	}
    	// return pointer to the beginning of the table
    	firstObject();
    }

    bool valueExists(EMDLabel name)
    {
        return (objects[current_objectID])->valueExists(name);
    }

    /** Check whether a label is contained in metadata.
     */
    bool containsLabel(const EMDLabel label) const;

    /** Get all active labels */
    std::vector<EMDLabel> getActiveLabels() const;

    /** Deactivate a column from a table, so that it is no longer written out
      */
    void deactivateLabel(EMDLabel label);

      /** Append the app table to this one
     */
    void append(MetaDataTable &app);

    /** Add a new label to the metadata.
     */
    bool addLabel(const EMDLabel label);

    /** Adds a new object to the objects map.
     *   If objectID == -1 the new ID will be that for the last object inserted + 1, else
     *   the given objectID is used.
     *   If there is already an object whose objectID == input objectID, the old one will be replaced by the new one
     *   If data !=NULL, the object will be set to contain these data
     **/
    long int addObject(MetaDataContainer * data = NULL, long int objectID = -1);

    /** Remove an object from the table. If objectID is not given, the current object will be removed
     * This function resets the current pointer to the last entry and returns the lastObject in the table */
    long int removeObject(long int objectID = -1);

    /* Get metadatacontainer for objectID (is current_objectID when -1)
     */
    MetaDataContainer * getObject(long int objectID = -1) const;

    /* Set metadatacontainer for current metadata object
     * This function assumes there already exists an object with objectID
     * If objectID==-1,, then the (assumedly existing) currrent_objectID will be set
     *
     * Use addObject() to set an object that does not yet exist
     */
    void setObject(MetaDataContainer * data, long int objectID = -1);

    // Possible error codes for the map
    enum errors
    {
        NO_OBJECTS_STORED = -1, // NOTE: Do not change this value (-1)
        NO_MORE_OBJECTS = -2,
        NO_OBJECT_FOUND = -3
    };

    long int firstObject();
    long int nextObject();
    long int lastObject();
    long int goToObject(long int objectID);

    /* Read a STAR loop structure
      */
    long int readStarLoop(std::ifstream& in, std::vector<EMDLabel> *labelsVector = NULL, std::string grep_pattern = "", bool do_only_count = false);

    /* Read a STAR list
     * The function returns true if the list is followed by a loop, false otherwise
     */
    bool readStarList(std::ifstream& in, std::vector<EMDLabel> *labelsVector = NULL);

    /* Read a MetaDataTable from a STAR-format data block
	 *
     * If the data block contains a list and a table, the function will return 2,
     *  the first time it is called and the list is read into the MetaDataTable
     *  in that case the function needs to be called another time. The second time
     *  it will read the _loop structure into the MetaDataTable and 1 will be returned
     *
     * If the data block contains only a list or a table, it is read in the MetaDataTable and the function will return 1
     *
     * If no data block is found the function will return 0 and the MetaDataTable remains empty
     *
     */
    long int readStar(std::ifstream& in, const std::string &name = "", std::vector<EMDLabel> *labelsVector = NULL, std::string grep_pattern = "", bool do_only_count = false);

    // Read a MetaDataTable (get fileformat from extension)
    long int read(const FileName &filename, const std::string &name = "", std::vector<EMDLabel> *labelsVector = NULL, std::string grep_pattern = "", bool do_only_count = false);

    // Write a MetaDataTable in STAR format
    void write(std::ostream& out = std::cout);

    // Write to a single file
    void write(const FileName & fn_out);


    void writeValueToString(std::string & result,
                            const std::string & inputLabel);


    void addToCPlot2D(CPlot2D *plot2D, EMDLabel xaxis, EMDLabel yaxis,
    		double red=0., double green=0., double blue=0., double linewidth = 1.0, std::string marker="");


};

void compareMetaDataTable(MetaDataTable &MD1, MetaDataTable &MD2,
		MetaDataTable &MDboth, MetaDataTable &MDonly1, MetaDataTable &MDonly2,
		EMDLabel label1, RFLOAT eps = 0., EMDLabel label2 = EMDL_UNDEFINED, EMDLabel label3 = EMDL_UNDEFINED);

// Join 2 metadata tables. Only include labels that are present in both of them.
MetaDataTable combineMetaDataTables(std::vector<MetaDataTable> &MDin);

// Feb14,2017 - Shaoda, Check whether the two MetaDataTables contain the same set of activeLabels
bool compareLabels(const MetaDataTable &MD1, const MetaDataTable &MD2);

#endif
