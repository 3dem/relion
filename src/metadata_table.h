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
 * Redesigned by: J. Zivanov in June 2017
 * MRC Laboratory of Molecular Biology
 *
 * Original author:     J.R. Bilbao-Castro (jrbcast@ace.ual.es)
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
#include "src/metadata_label.h"

/** For all objects.
 @code
 FOR_ALL_OBJECTS_IN_METADATA(metadata) {
   RFLOAT rot;
   DF.getValue( EMDL_ANGLEROT, rot);
 }
 @endcode

 @TODO: remove "&& current_object >= 0"
        and make "nextObject()" return "current_object++"
        after "enum errors" has been removed (see below)
 */
#define FOR_ALL_OBJECTS_IN_METADATA_TABLE(mdt_arg) \
        for(long int current_object = (mdt_arg).firstObject(); \
             current_object < (mdt_arg).numberOfObjects() \
                && current_object >= 0; \
             current_object = (mdt_arg).nextObject())

/*  class MetaDataTable:
 *
 *  - stores a table of values for an arbitrary subset of predefined EMDLabels
 *  - each column corresponds to a label
 *  - each row represents a data point
 *  - the rows are stored in per-type contiguous blocks of memory
 *  -
 */
class MetaDataTable
{
    // Effectively stores all metadata
    std::vector<MetaDataContainer*> objects;

    // Maps labels to corresponding indices in the vectors in MetaDataContainer.
    // The length of label2offset is always equal to the number of defined labels (~320)
    // e.g.:
    // the value of "defocus-U" for row r is stored in:
    //   objects[r]->doubles[label2offset[EMDL_CTF_DEFOCUSU]]
    // the value of "image name" is stored in:
    //   objects[r]->strings[label2offset[EMDL_IMAGE_NAME]]
    std::vector<long> label2offset;

    // Current object id
    long current_objectID;

    // Number of labels of each type
    long doubleLabels, intLabels, boolLabels, stringLabels;

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



    MetaDataTable();

    // Copy constructor and assignment operator:
    // Fill the new table with *copies* of all objects
    MetaDataTable(const MetaDataTable & c);
    MetaDataTable& operator = (const MetaDataTable &MD);

    ~MetaDataTable();



    void setIsList(bool is_list);

    bool isEmpty() const;
    size_t numberOfObjects() const;
    size_t size(void); // @TODO: redundant
    void clear();

    void setComment(const std::string Comment);
    std::string getComment() const;
    bool containsComment() const;

    void setName(const std::string Name);
    std::string getName() const;


    //  getValue: returns true if the label exists
    template<class T>
    bool getValue(EMDLabel label, T& value, long objectID = -1) const;

    bool getValueToString(EMDLabel label, std::string &value, long int objectID = -1) const;

    // Set the value of label for a specified object.
    // If no objectID is given, the internal iterator 'current_objectID' is used
    template<class T>
    bool setValue(EMDLabel name, const T &value, long int objectID = -1);

    bool setValueFromString(EMDLabel label, const std::string &value, long int objectID = -1);


    // Sort the order of the elements based on the values in the input label
    // (only numbers, no strings/bools)
    void sort(EMDLabel name, bool do_reverse = false, bool only_set_index = false, bool do_random = false);
    void newSort(const EMDLabel name, bool do_reverse = false, bool do_sort_after_at = false, bool do_sort_before_at = false);

    // Check whether a label is defined in the table. (@TODO: change name)
    bool labelExists(EMDLabel name) const;

    // Check whether a label is contained in activeLabels. (@TODO: change name)
    bool containsLabel(const EMDLabel label) const;

    std::vector<EMDLabel> getActiveLabels() const; // @TODO: redundant; activeLabels is public

    // Deactivate a column from a table, so that it is no longer written out
    void deactivateLabel(EMDLabel label);

    // add a new label and update all objects
    void addLabel(EMDLabel label);

    // add missing labels that are present in 'app'
    void addMissingLabels(const MetaDataTable* app);

    // add all rows from app to the end of the table and
    // insert all missing labels
    void append(const MetaDataTable& app);

    // Get metadatacontainer for objectID (current_objectID if objectID < 0)
    MetaDataContainer* getObject(long objectID = -1) const;

    /* setObject(data, objectID)
     *   copies values from 'data' to object 'objectID'.
     *   The target object is assumed to exist.
     *   If objectID < 0, then current_objectID is set.
     *   Undefined labels are inserted.
     *
     *   Use addObject() to set an object that does not yet exist */
    void setObject(MetaDataContainer* data, long objectID = -1);

    /*  setValuesOfDefinedLabels(data, objectID)
     *   copies values from 'data' to object 'objectID'.
     *   The target object is assumed to exist.
     *   If objectID < 0, then current_objectID is set.
     *   Only already defined labels are considered.
     *
     *   Use addValuesOfDefinedLabels() to add an object that does not yet exist */
    void setValuesOfDefinedLabels(MetaDataContainer* data, long objectID = -1);

    // reserve memory for this many lines
    void reserve(size_t capacity);

    /*  addObject()
     *    Adds a new object and initializes the defined labels with default values.
     *    Afterwards, 'current_objectID' points to the newly added object.*/
    void addObject();

    /*  addObject(data)
     *    Adds a new object and sets its values to those from 'data'.
     *    The set of labels for the table is extended as necessary.
     *    Afterwards, 'current_objectID' points to the newly added object.*/
    void addObject(MetaDataContainer* data);

    /*  addValuesOfDefinedLabels(data)
     *    Adds a new object and sets the already defined values to those from 'data'.
     *    Labels from 'data' that are not already defined are ignored.
     *    Afterwards, 'current_objectID' points to the newly added object.*/
    void addValuesOfDefinedLabels(MetaDataContainer* data);

    /*  removeObject(objectID)
     *    If objectID is not given, 'current_objectID' will be removed.
     *    'current_objectID' is set to the last object in the list. */
    void removeObject(long objectID = -1);

    long firstObject();
    long nextObject();
    // @TODO: remove nextObject() after removing calls in:
    // - "particle_reposition.cpp"
    // - "helix.cpp"
    // - "preprocessing.cpp"

    long goToObject(long objectID);

    // Read a STAR loop structure
    long int readStarLoop(std::ifstream& in, std::vector<EMDLabel> *labelsVector = NULL, std::string grep_pattern = "", bool do_only_count = false);

    /* Read a STAR list
     * The function returns true if the list is followed by a loop, false otherwise */
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
     */
    long int readStar(std::ifstream& in, const std::string &name = "", std::vector<EMDLabel> *labelsVector = NULL, std::string grep_pattern = "", bool do_only_count = false);

    // Read a MetaDataTable (get file format from extension)
    long int read(const FileName &filename, const std::string &name = "", std::vector<EMDLabel> *labelsVector = NULL, std::string grep_pattern = "", bool do_only_count = false);

    // Write a MetaDataTable in STAR format
    void write(std::ostream& out = std::cout) const;

    // Write to a single file
    void write(const FileName & fn_out) const;

    void addToCPlot2D(CPlot2D *plot2D, EMDLabel xaxis, EMDLabel yaxis,
    		double red=0., double green=0., double blue=0., double linewidth = 1.0, std::string marker="");

    void printLabels(std::ostream& ost);

    // legacy error codes:
    // @TODO: remove after changing:
    //   - particle_reposition.cpp, line ~127
    //   - preprocessing.cpp, line ~299
    enum errors
    {
        NO_OBJECTS_STORED = -1,
        NO_MORE_OBJECTS = -2,
        NO_OBJECT_FOUND = -3
    };

private:

    // Check if 'id' corresponds to an actual object.
    // Crash if it does not.
    void checkObjectID(long id, std::string caller) const;

    /*  setObjectUnsafe(data)
     *    Same as setObject, but assumes that all labels are present. */
    void setObjectUnsafe(MetaDataContainer* data, long objId);

};

void compareMetaDataTable(MetaDataTable &MD1, MetaDataTable &MD2,
		MetaDataTable &MDboth, MetaDataTable &MDonly1, MetaDataTable &MDonly2,
		EMDLabel label1, double eps = 0., EMDLabel label2 = EMDL_UNDEFINED, EMDLabel label3 = EMDL_UNDEFINED);

// Join 2 metadata tables. Only include labels that are present in both of them.
MetaDataTable combineMetaDataTables(std::vector<MetaDataTable> &MDin);

// Feb14,2017 - Shaoda, Check whether the two MetaDataTables contain the same set of activeLabels
bool compareLabels(const MetaDataTable &MD1, const MetaDataTable &MD2);



template<class T>
bool MetaDataTable::getValue(EMDLabel label, T& value, long objectID) const
{
    if (label < 0 || label >= EMDL_LAST_LABEL) return false;

    const long off = label2offset[label];
    if (off > -1)
    {
        if (objectID < 0)
        {
            objectID = current_objectID;
        }
        else
        	checkObjectID(objectID,  "MetaDataTable::getValue");

        objects[objectID]->getValue(off, value);
        return true;
    }
    else
    {
        return false;
    }
}

template<class T>
bool MetaDataTable::setValue(EMDLabel name, const T &value, long int objectID)
{
    if (name < 0 || name >= EMDL_LAST_LABEL) return false;

    long off = label2offset[name];

    if (off < 0)
    {
        addLabel(name);
        off = label2offset[name];
    }

    if (objectID < 0)
    {
        objectID = current_objectID;
    }
    else checkObjectID(objectID,  "MetaDataTable::setValue");

    if (off > -1)
    {
        objects[objectID]->setValue(off, value);
        return true;
    }
    else
    {
        return false;
    }
}

#endif
