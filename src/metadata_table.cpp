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
 * Authors:      J.R. Bilbao-Castro (jrbcast@ace.ual.es)
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

#include "src/metadata_table.h"
#include "src/metadata_label.h"


MetaDataTable::MetaDataTable()
:   objects(0),
    label2offset(EMDL_LAST_LABEL, -1),
    current_objectID(0),
    doubleLabels(0),
    intLabels(0),
    boolLabels(0),
    stringLabels(0),
    isList(false),
    name(""),
    comment(""),
    activeLabels(0),
    ignoreLabels(0)
{
}

MetaDataTable::MetaDataTable(const MetaDataTable &MD)
:   objects(MD.objects.size()),
    label2offset(MD.label2offset),
    current_objectID(0),
    doubleLabels(MD.doubleLabels),
    intLabels(MD.intLabels),
    boolLabels(MD.boolLabels),
    stringLabels(MD.stringLabels),
    isList(MD.isList),
    name(MD.name),
    comment(MD.comment),
    activeLabels(MD.activeLabels),
    ignoreLabels(MD.ignoreLabels)
{
    for (size_t idx = 0; idx < MD.objects.size(); idx++)
    {
        objects[idx] = new MetaDataContainer(*(MD.objects[idx]));
        objects[idx]->table = this;
    }
}

MetaDataTable& MetaDataTable::operator = (const MetaDataTable &MD)
{
    if (this != &MD)
    {
        clear();

        objects.resize(MD.objects.size());
        label2offset = MD.label2offset;
        current_objectID = 0;
        doubleLabels = MD.doubleLabels;
        intLabels = MD.intLabels;
        boolLabels = MD.boolLabels;
        stringLabels = MD.stringLabels;
        isList = MD.isList;
        name = MD.name;
        comment = MD.comment;

        activeLabels = MD.activeLabels;
        ignoreLabels = MD.ignoreLabels;

        for (long int idx = 0; idx < MD.objects.size(); idx++)
        {
            objects[idx] = new MetaDataContainer(this, MD.objects[idx]);
            objects[idx]->table = this;
        }
    }

    return *this;
}

void MetaDataTable::setIsList(bool is_list)
{
    isList = is_list;
}

MetaDataTable::~MetaDataTable()
{
    for (long i = 0; i < objects.size(); i++)
    {
        delete objects[i];
    }
}

bool MetaDataTable::isEmpty() const
{
    return (objects.size()==0);
}

size_t MetaDataTable::numberOfObjects() const
{
	return objects.size();
}

void MetaDataTable::clear()
{
    for (long i = 0; i < objects.size(); i++)
    {
        delete objects[i];
    }
    objects.clear();

    label2offset = std::vector<long>(EMDL_LAST_LABEL, -1);
    current_objectID = 0;

    doubleLabels = 0;
    intLabels = 0;
    boolLabels = 0;
    stringLabels = 0;

    isList = false;
    name = "";
    comment = "";

    activeLabels.clear();
    ignoreLabels.clear();
}

void MetaDataTable::setComment(const std::string newComment)
{
	comment = newComment;
}

std::string MetaDataTable::getComment() const
{
    return comment;
}

bool MetaDataTable::containsComment() const
{
	return (comment != std::string(""));
}

void MetaDataTable::setName(const std::string newName)
{
	name = newName;
}

std::string MetaDataTable::getName() const
{
    return name;
}

bool MetaDataTable::getValueToString(EMDLabel label, std::string &value, long objectID) const
{
    if (EMDL::isString(label))
    {
        getValue(label, value, objectID);
    }
    else
    {
        std::stringstream sts(value);

        if (EMDL::isDouble(label))
        {
            double v;
            getValue(label, v, objectID);
            if ((ABS(v) > 0. && ABS(v) < 0.001) || ABS(v) > 100000.)
                 sts << std::setw(12) << std::scientific;
             else
                 sts << std::setw(12) << std::fixed;
              sts << v;
        }
        else if (EMDL::isInt(label))
        {
            long v;
            getValue(label, v, objectID);
            sts << std::setw(12) << std::fixed;
            sts << v;
        }
        else if (EMDL::isBool(label))
        {
            bool v;
            getValue(label, v, objectID);
            sts << std::setw(12) << std::fixed;
            sts << v;
        }

        value = sts.str();
    }

}

size_t MetaDataTable::size()
{
    return objects.size();
}

bool MetaDataTable::setValueFromString(EMDLabel label, const std::string &value,
                        long int objectID)
{
    if (EMDL::isString(label))
    {
        return setValue(label, value, objectID);
    }
    else
    {
        std::istringstream i(value);

        if (EMDL::isDouble(label))
        {
            double v;
            i >> v;
            return setValue(label, v, objectID);
        }
        else if (EMDL::isInt(label))
        {
            long v;
            i >> v;
            return setValue(label, v, objectID);
        }
        else if (EMDL::isBool(label))
        {
            bool v;
            i >> v;
            return setValue(label, v, objectID);
        }
    }
}

    // comparators used for sorting

    struct MdDoubleComparator
    {
        MdDoubleComparator(long index) : index(index) {}

        bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const
        {
            return lh->doubles[index] < rh->doubles[index];
        }

        long index;
    };

    struct MdIntComparator
    {
        MdIntComparator(long index) : index(index) {}

        bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const
        {
            return lh->ints[index] < rh->ints[index];
        }

        long index;
    };

    struct MdStringComparator
    {
        MdStringComparator(long index) : index(index) {}

        bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const
        {
            return lh->strings[index] < rh->strings[index];
        }

        long index;
    };

    struct MdStringAfterAtComparator
    {
        MdStringAfterAtComparator(long index) : index(index) {}

        bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const
        {
            std::string slh = lh->strings[index];
            std::string srh = rh->strings[index];
            slh = slh.substr(slh.find("@")+1);
            srh = srh.substr(srh.find("@")+1);
            return slh < srh;
        }

        long index;
    };

    struct MdStringBeforeAtComparator
    {
        MdStringBeforeAtComparator(long index) : index(index) {}

        bool operator()(MetaDataContainer *lh, MetaDataContainer *rh) const
        {
            std::string slh = lh->strings[index];
            std::string srh = rh->strings[index];
            slh = slh.substr(0, slh.find("@"));
            srh = srh.substr(0, srh.find("@"));
            std::stringstream stslh, stsrh;
            stslh << slh;
            stsrh << srh;
            long ilh, irh;
            stslh >> ilh;
            stsrh >> irh;

            return ilh < irh;
        }

        long index;
    };

void MetaDataTable::sort(EMDLabel name, bool do_reverse, bool only_set_index, bool do_random)
{

    if (do_random)
        srand (time(NULL));    		  /* initialize random seed: */
    else if (!(EMDL::isInt(name) || EMDL::isDouble(name)) )
        REPORT_ERROR("MetadataTable::sort%% ERROR: can only sorted numbers");

    std::vector<std::pair<double,long int> > vp;
    vp.reserve(objects.size());
    long int i = 0;
    FOR_ALL_OBJECTS_IN_METADATA_TABLE(*this)
    {
        double dval;
        if (do_random)
        {
            dval = (double)rand();
        }
        else if (EMDL::isInt(name))
        {
            long val;
            getValue(name, val);
            dval = (double) val;
        }
        else // EMDL::isDouble(name)
        {
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
        for (long j = 0; j < vp.size(); j++)
        {
            (*this).setValue(EMDL_SORTED_IDX, j, vp[j].second);
        }
    }
    else
    {
        // Change the actual order in the MetaDataTable
        std::vector<MetaDataContainer*> objs(objects.size());

        for (long j = 0; j < vp.size(); j++)
        {
            objs[j] = objects[vp[j].second];
        }

        objects = objs;
    }
    // reset pointer to the beginning of the table
    firstObject();
}

void MetaDataTable::newSort(const EMDLabel label, bool do_reverse, bool do_sort_after_at, bool do_sort_before_at)
{
    if (EMDL::isString(label))
    {
        if (do_sort_after_at)
        {
            std::stable_sort(objects.begin(), objects.end(),
                             MdStringAfterAtComparator(label2offset[label]));
        }
        else if (do_sort_before_at)
        {
            std::stable_sort(objects.begin(), objects.end(),
                             MdStringBeforeAtComparator(label2offset[label]));
        }
        else
        {
            std::stable_sort(objects.begin(), objects.end(), MdStringComparator(label2offset[label]));
        }
    }
    else if (EMDL::isDouble(label))
    {
        std::stable_sort(objects.begin(), objects.end(), MdDoubleComparator(label2offset[label]));
    }
    else if (EMDL::isInt(label))
    {
        std::stable_sort(objects.begin(), objects.end(), MdIntComparator(label2offset[label]));
    }
    else
    {
        REPORT_ERROR("Cannot sort this label: " + EMDL::label2Str(label));
    }

    if (do_reverse)
    {
        std::reverse(objects.begin(), objects.end());
    }
}

bool MetaDataTable::labelExists(EMDLabel name) const
{
    if (name < 0 || name >= EMDL_LAST_LABEL) return false;
    else return (label2offset[name] >= 0);
}

bool MetaDataTable::containsLabel(const EMDLabel label) const
{
    return vectorContainsLabel(activeLabels, label);
}

std::vector<EMDLabel> MetaDataTable::getActiveLabels() const
{
    return activeLabels;
}

void MetaDataTable::deactivateLabel(EMDLabel label)
{
    std::vector<EMDLabel>::iterator location
            = std::find(activeLabels.begin(), activeLabels.end(), label);

    if (location != activeLabels.end())
    {
        activeLabels.erase(location);
    }
}

void MetaDataTable::addLabel(EMDLabel label)
{
    if (label < 0 || label >= EMDL_LAST_LABEL)
    {
        std::stringstream sts;
        sts << label;
        REPORT_ERROR("MetaDataTable::addLabel: unrecognised label: " + sts.str());
    }

    if (label2offset[label] < 0)
    {
        activeLabels.push_back(label);
        long id;

        if (EMDL::isDouble(label))
        {
            id = doubleLabels;

            for (long i = 0; i < objects.size(); i++)
            {
                objects[i]->doubles.push_back(0);
            }

            doubleLabels++;
        }
        else if (EMDL::isInt(label))
        {
            id = intLabels;

            for (long i = 0; i < objects.size(); i++)
            {
                objects[i]->ints.push_back(0);
            }

            intLabels++;
        }
        else if (EMDL::isBool(label))
        {
            id = boolLabels;

            for (long i = 0; i < objects.size(); i++)
            {
                objects[i]->bools.push_back(false);
            }

            boolLabels++;
        }
        else if (EMDL::isString(label))
        {
            id = stringLabels;

            for (long i = 0; i < objects.size(); i++)
            {
                objects[i]->strings.push_back("empty");
            }

            stringLabels++;
        }

        label2offset[label] = id;
    }
}

void MetaDataTable::addMissingLabels(const MetaDataTable* mdt)
{
    for (long i = 0; i < mdt->activeLabels.size(); i++)
    {
        EMDLabel l = mdt->activeLabels[i];

        if (label2offset[l] < 0)
        {
            addLabel(l);
        }
    }
}

void MetaDataTable::append(const MetaDataTable& mdt)
{

	if (activeLabels.size() == 0)
	{
		// If the current one is empty, add missing labels and append the new one:
		addMissingLabels(&mdt);
	}
	else
	{
		// If the current one is not-empty, check all labels are the same before appending. Otherwise, raise error
		if (!compareLabels(*this, mdt))
			REPORT_ERROR("ERROR in appending metadata tables with not the same columns!");
	}

	// Now append
    objects.reserve(objects.size() + mdt.numberOfObjects());
    for (long i = 0; i < mdt.objects.size(); i++)
    {
        objects.push_back(new MetaDataContainer(
            this, doubleLabels, intLabels, boolLabels, stringLabels));

        setObjectUnsafe(mdt.getObject(i), objects.size() - 1);
    }

    // reset pointer to the beginning of the table
    firstObject();
}


MetaDataContainer* MetaDataTable::getObject(long objectID) const
{
    if (objectID < 0) objectID = current_objectID;

    checkObjectID(objectID,  "MetaDataTable::getObject");

    return objects[objectID];
}

void MetaDataTable::setObject(MetaDataContainer* data, long objectID)
{
    if (objectID < 0) objectID = current_objectID;

    checkObjectID(objectID,  "MetaDataTable::setObject");
    addMissingLabels(data->table);

    setObjectUnsafe(data, objectID);
}

void MetaDataTable::setValuesOfDefinedLabels(MetaDataContainer* data, long objectID)
{
    if (objectID < 0) objectID = current_objectID;

    checkObjectID(objectID,  "MetaDataTable::setValuesOfDefinedLabels");

    setObjectUnsafe(data, objectID);
}

void MetaDataTable::reserve(size_t capacity)
{
    objects.reserve(capacity);
}

void MetaDataTable::setObjectUnsafe(MetaDataContainer* data, long objectID)
{
    MetaDataContainer* obj = objects[objectID];

    for (long i = 0; i < data->table->activeLabels.size(); i++)
    {
        EMDLabel label = data->table->activeLabels[i];

        long myOff = label2offset[label];
        long srcOff = data->table->label2offset[label];

        if (myOff < 0) continue;

        if (EMDL::isDouble(label))
        {
            obj->doubles[myOff] = data->doubles[srcOff];
        }
        else if (EMDL::isInt(label))
        {
            obj->ints[myOff] = data->ints[srcOff];
        }
        else if (EMDL::isBool(label))
        {
            obj->bools[myOff] = data->bools[srcOff];
        }
        else if (EMDL::isString(label))
        {
            obj->strings[myOff] = data->strings[srcOff];
        }
    }
}

void MetaDataTable::addObject()
{
    objects.push_back(new MetaDataContainer(
        this, doubleLabels, intLabels, boolLabels, stringLabels));

    current_objectID = objects.size()-1;
}

void MetaDataTable::addObject(MetaDataContainer* data)
{
    objects.push_back(new MetaDataContainer(
        this, doubleLabels, intLabels, boolLabels, stringLabels));

    setObject(data, objects.size()-1);
    current_objectID = objects.size()-1;
}

void MetaDataTable::addValuesOfDefinedLabels(MetaDataContainer* data)
{
    objects.push_back(new MetaDataContainer(
        this, doubleLabels, intLabels, boolLabels, stringLabels));

    setValuesOfDefinedLabels(data, objects.size()-1);
    current_objectID = objects.size()-1;
}

void MetaDataTable::removeObject(long objectID)
{
    long i = (objectID < 0) ? current_objectID : objectID;

    checkObjectID(i,  "MetaDataTable::removeObject");

    delete objects[i];
    objects.erase(objects.begin() + i);

    current_objectID = objects.size() - 1;
}

long int MetaDataTable::firstObject()
{
    current_objectID = 0;
    return 0;
}

long int MetaDataTable::nextObject()
{
    current_objectID++;

    if (current_objectID >= objects.size())
    {
        return NO_MORE_OBJECTS;
    }
    else
    {
        return current_objectID;
    }
}

long int MetaDataTable::goToObject(long int objectID)
{
    checkObjectID(objectID, "MetaDataTable::goToObject");

    current_objectID = objectID;
    return current_objectID;
}

long int MetaDataTable::readStarLoop(std::ifstream& in, std::vector<EMDLabel> *desiredLabels, std::string grep_pattern, bool do_only_count)
{
	setIsList(false);

	//Read column labels
    int labelPosition = 0;
    std::string line, token;

    // First read all the column labels
    while (getline(in, line, '\n'))
    {
		line = simplify(line);
		// TODO: handle comments...
		if (line[0] == '#' || line[0] == '\0' || line[0] == ';')
			continue;

		if (line[0] == '_') // label definition line
		{
			//Only take string from "_" until "#"
            size_t pos0 = line.find("_");
            size_t pos1 = line.find("#");

            token = line.substr(pos0 + 1, pos1 - pos0 - 2);

            EMDLabel label = EMDL::str2Label(token);

			//std::cerr << " label= XX" << label << "XX token= XX" << token<<"XX" << std::endl;
			if (desiredLabels != NULL && !vectorContainsLabel(*desiredLabels, label))
				label = EMDL_UNDEFINED; //ignore if not present in desiredLabels

			if (label == EMDL_UNDEFINED)
			{
				//std::cerr << "Warning: ignoring the following (undefined) label:" <<token << std::endl;
				REPORT_ERROR("ERROR: Unrecognised metadata label: " + token);
				ignoreLabels.push_back(labelPosition);
			}
			else
            {
                addLabel(label);
            }

			labelPosition++;
		}
		else // found first data line
		{
			break;
		}
    }

    // Then fill the table (dont read another line until the one from above has been handled)
    bool is_first = true;
    long int nr_objects = 0;

    while (is_first || getline(in, line, '\n'))
    {
        is_first = false;

        if (grep_pattern == "" || line.find(grep_pattern) != std::string::npos)
    	{
    		line = simplify(line);
			// Stop at empty line
			if (line[0] == '\0')
				break;

			nr_objects++;
			if (!do_only_count)
			{
				// Add a new line to the table
				addObject();

				// Parse data values
				std::stringstream os2(line);
				std::string value;
				labelPosition = 0;
				int counterIgnored = 0;
				while (os2 >> value)
				{
					// TODO: handle comments here...
					if (std::find(ignoreLabels.begin(), ignoreLabels.end(), labelPosition) != ignoreLabels.end())
					{
						// Ignore this column
						counterIgnored++;
						labelPosition++;
						continue;
					}
                    setValueFromString(activeLabels[labelPosition - counterIgnored], value);
					labelPosition++;
				}
			}
    	} // end if grep_pattern
    }

    return nr_objects;
}

bool MetaDataTable::readStarList(std::ifstream& in, std::vector<EMDLabel> *desiredLabels)
{
	setIsList(true);
    addObject();
    long int objectID = objects.size() - 1;

    std::string line, firstword, value;
    std::vector<std::string> words;

    bool also_has_loop = false;

    // Read data and fill structures accordingly
    while (getline(in, line, '\n'))
    {
    	 tokenize(line, words);

    	 // Ignore empty lines
    	 if (words.size() == 0)
    		 continue;
    	 else
    		 firstword = words[0];

    	 // Get label-value pairs
    	 if (firstword[0] == '_')
    	 {
             EMDLabel label = EMDL::str2Label(firstword.substr(1)); // get rid of leading underscore
        	 if (words.size() != 2)
        		 REPORT_ERROR("MetaDataTable::readStarList: did not encounter a single word after "+firstword);
    		 value = words[1];

             if (desiredLabels != NULL && !vectorContainsLabel(*desiredLabels, label))
             {
                 label = EMDL_UNDEFINED; //ignore if not present in desiredLabels
             }
    		 if (label != EMDL_UNDEFINED)
             {
                 addLabel(label);
				 setValueFromString(label, value, objectID);
			 }
    	 }
    	 // Check whether there is a comment or an empty line
    	 else if (firstword[0] == '#' || firstword[0] == ';')
    	 {
    		 // TODO: handle comments?
    		 continue;
    	 }
    	 // Check whether a loop structure comes after this list
    	 else if (firstword.find("loop_") == 0)
    	 {
    		 also_has_loop = true;
    		 return also_has_loop;
    	 }
    	 // Check whether this data blocks ends (because a next one is there)
    	 else if (firstword.find("data_") == 0)
    	 {
    		 // Should I reverse the pointer one line?
    		 return also_has_loop;
    	 }
     }
     // Reached the end of the file
     return also_has_loop;
}

long int MetaDataTable::readStar(std::ifstream& in, const std::string &name, std::vector<EMDLabel> *desiredLabels, std::string grep_pattern, bool do_only_count)
{
    std::stringstream ss;
    std::string line, token, value;
    clear();
    bool also_has_loop;

    // Start reading the ifstream at the top
    in.seekg(0);

    // Proceed until the next data_ or _loop statement
    // The loop statement may be necessary for data blocks that have a list AND a table inside them
    while (getline(in, line, '\n'))
    {
    	// Find data_ lines
    	if (line.find("data_") != std::string::npos)
    	{
    		token = line.substr(line.find("data_") + 5);
    		// If a name has been given, only read data_thatname
    		// Otherwise, just read the first data_ block
    		if (name == "" || name == token)
    		{
    			setName(token);
    			// Get the next item that starts with "_somelabel" or with "loop_"
    			int current_pos = in.tellg();
    			while (getline(in, line, '\n'))
    			{
    				trim(line);
    				if (line.find("loop_") != std::string::npos)
    				{
    					return readStarLoop(in, desiredLabels, grep_pattern, do_only_count);
    				}
    				else if (line[0] == '_')
    				{
    					// go back one line in the ifstream
    					in.seekg(current_pos);
    					also_has_loop = readStarList(in, desiredLabels);
    					return (also_has_loop) ? 0 : 1;
    				}
    			}
    		}
    	}
    }

    // Clear the eofbit so we can perform more actions on the stream.
    in.clear();

    return 0;
}

long int MetaDataTable::read(const FileName &filename, const std::string &name, std::vector<EMDLabel> *desiredLabels, std::string grep_pattern, bool do_only_count)
{

    // Clear current table
    clear();

    // Check for an :star extension
    FileName fn_read = filename.removeFileFormat();

    std::ifstream in(fn_read.data(), std::ios_base::in);
    if (in.fail())
        REPORT_ERROR( (std::string) "MetaDataTable::read: File " + fn_read + " does not exist" );

    FileName ext = filename.getFileFormat();
    if (ext =="star")
    {
        //REPORT_ERROR("readSTAR not implemented yet...");
        return readStar(in, name, desiredLabels, grep_pattern, do_only_count);
    }
    else
    {
        REPORT_ERROR("MetaDataTable::read ERROR: metadata table should have a .star extension");
    }

    in.close();

    // Go to the first object
    firstObject();

}

void MetaDataTable::write(std::ostream& out) const
{
    // Only write tables that have something in them
    if (isEmpty())
        return;

    out << "\n";
    out << "data_" << getName() <<"\n";
    if (containsComment())
    	out << "# "<< comment << "\n";
    out << "\n";

    if (!isList)
    {
        // Write loop header structure
        out << "loop_ \n";
        for (long i = 0; i < activeLabels.size(); i++)
        {
            EMDLabel l = activeLabels[i];

            if (l != EMDL_COMMENT && l != EMDL_SORTED_IDX) // EMDL_SORTED_IDX is only for internal use, never write it out!
            {
                out << "_" << EMDL::label2Str(l) << " #" << i+1 << " \n";
            }
        }

        // Write actual data block
        for (long int idx = 0; idx < objects.size(); idx++)
        {
            std::string entryComment = "";

            for (long i = 0; i < activeLabels.size(); i++)
            {
                EMDLabel l = activeLabels[i];

                if (l != EMDL_COMMENT && l != EMDL_SORTED_IDX)
                {
                    out.width(10);
                    std::string val;
                    getValueToString(l, val, idx);
                    out << val << " ";
                }
                if (l == EMDL_COMMENT)
                {
                    getValue(EMDL_COMMENT, entryComment, idx);
                }
            }
            if (entryComment != std::string(""))
            {
            	out << "# " << entryComment;
            }
            out << "\n";
        }
        // Finish table with a white-line
        out << " \n";

    }
    else
    {
        // Get first object. In this case (row format) there is a single object
        std::string entryComment = "";
        int maxWidth=10;

        for (long i = 0; i < activeLabels.size(); i++)
        {
            EMDLabel l = activeLabels[i];

            if (l != EMDL_COMMENT)
            {
                int w = EMDL::label2Str(l).length();

                if (w > maxWidth) maxWidth = w;
            }
            else
            {
                getValue(EMDL_COMMENT, entryComment, 0);
            }
        }

        for (long i = 0; i < activeLabels.size(); i++)
        {
            EMDLabel l = activeLabels[i];

            if (l != EMDL_COMMENT)
            {
                int w = EMDL::label2Str(l).length();
                out << "_" << EMDL::label2Str(l) << std::setw(12 + maxWidth - w) << " ";

                std::string val;
                getValueToString(l, val, 0);
                out << val << "\n";
            }
        }
        if (entryComment != std::string(""))
        {
        	out << "# " << entryComment << "\n";
        }

        // End a data block with a white line
        out << " \n";
    }
}

void MetaDataTable::write(const FileName &fn_out) const
{
    std::ofstream  fh;
    fh.open((fn_out).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR( (std::string)"MetaDataTable::write: cannot write to file: " + fn_out);
    write(fh);
    fh.close();

}

void MetaDataTable::addToCPlot2D(CPlot2D *plot2D, EMDLabel xaxis, EMDLabel yaxis,
		double red, double green, double blue, double linewidth, std::string marker)
{
	CDataSet dataSet;
	if (marker=="")
	{
		dataSet.SetDrawMarker(false);
	}
	else
	{
		dataSet.SetDrawMarker(true);
		dataSet.SetMarkerSymbol(marker);
	}
	dataSet.SetLineWidth(linewidth);
	dataSet.SetDatasetColor(red, green, blue);
	dataSet.SetDatasetTitle(EMDL::label2Str(yaxis));

    double mydbl;
    long int myint;
	double xval, yval;
	for (long int idx = 0; idx < objects.size(); idx++)
    {
		const long offx = label2offset[xaxis];
		if (offx < 0)
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: cannot find x-axis label");


		if (EMDL::isDouble(xaxis))
    	{
    		objects[idx]->getValue(offx, mydbl);
    		xval = mydbl;
    	}
    	else if (EMDL::isInt(xaxis))
    	{
    		objects[idx]->getValue(offx, myint);
    		xval = myint;
    	}
    	else
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: can only plot x-axis double, int or long int");

		const long offy = label2offset[yaxis];
		if (offy < 0)
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: cannot find y-axis label");

		if (EMDL::isDouble(yaxis))
    	{
    		objects[idx]->getValue(offy, mydbl);
    		yval = mydbl;
    	}
    	else if (EMDL::isInt(yaxis))
    	{
    		objects[idx]->getValue(offy, myint);
    		yval = myint;
    	}
    	else
            REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: can only plot y-axis double, int or long int");

    	CDataPoint point(xval, yval);
		dataSet.AddDataPoint(point);

	}

    plot2D->AddDataSet(dataSet);

}

void MetaDataTable::printLabels(std::ostream &ost)
{
    for (int i = 0; i < activeLabels.size(); i++)
    {
        ost << EMDL::label2Str(activeLabels[i]) << "\n";
    }
}

void MetaDataTable::checkObjectID(long id, std::string caller) const
{
    if (id >= objects.size() || id < 0)
    {
        std::stringstream sts0, sts1;
        sts0 << id;
        sts1 << objects.size();
        REPORT_ERROR(caller+": object " + sts0.str()
                     + " out of bounds! (" + sts1.str() + " objects present)");
    }
}

void compareMetaDataTable(MetaDataTable &MD1, MetaDataTable &MD2,
		MetaDataTable &MDboth, MetaDataTable &MDonly1, MetaDataTable &MDonly2,
		EMDLabel label1, double eps, EMDLabel label2, EMDLabel label3)
{
	if (!MD1.containsLabel(label1))
		REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD1 does not contain the specified label1.");
	if (!MD2.containsLabel(label1))
		REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD2 does not contain the specified label1.");

	if (label2 != EMDL_UNDEFINED)
	{
		if (!EMDL::isDouble(label1) || !EMDL::isDouble(label2))
			REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR 2D or 3D distances are only allowed for doubles.");
		if (!MD1.containsLabel(label2))
			REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD1 does not contain the specified label2.");
		if (!MD2.containsLabel(label2))
			REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD2 does not contain the specified label2.");
	}

	if (label3 != EMDL_UNDEFINED)
	{
		if (!EMDL::isDouble(label3))
			REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR 3D distances are only allowed for doubles.");
		if (!MD1.containsLabel(label3))
			REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD1 does not contain the specified label3.");
		if (!MD2.containsLabel(label3))
			REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD2 does not contain the specified label3.");
	}

	MDboth.clear();
	MDonly1.clear();
	MDonly2.clear();

	std::string mystr1, mystr2;
	long int myint1, myint2;
	double myd1, myd2, mydy1 = 0., mydy2 = 0., mydz1 = 0., mydz2 = 0.;


	// loop over MD1
	std::vector<long int> to_remove_from_only2;
	for (long int current_object1 = MD1.firstObject();
	              current_object1 != MetaDataTable::NO_MORE_OBJECTS && current_object1 != MetaDataTable::NO_OBJECTS_STORED;
	              current_object1 = MD1.nextObject())
	{
		if (EMDL::isString(label1))
			MD1.getValue(label1, mystr1);
		else if (EMDL::isInt(label1))
			MD1.getValue(label1, myint1);
		else if (EMDL::isDouble(label1))
		{
			MD1.getValue(label1, myd1);
			if (label2 != EMDL_UNDEFINED)
				MD1.getValue(label2, mydy1);
			if (label3 != EMDL_UNDEFINED)
				MD1.getValue(label3, mydz1);
		}
		else
			REPORT_ERROR("compareMetaDataTableEqualLabel ERROR: only implemented for strings, integers or doubles");

		// loop over MD2
		bool have_in_2 = false;
		for (long int current_object2 = MD2.firstObject();
		              current_object2 != MetaDataTable::NO_MORE_OBJECTS && current_object2 != MetaDataTable::NO_OBJECTS_STORED;
		              current_object2 = MD2.nextObject())
		{

			if (EMDL::isString(label1))
			{
				MD2.getValue(label1, mystr2);
				if (strcmp(mystr1.c_str(), mystr2.c_str()) == 0)
				{
					have_in_2 = true;
					to_remove_from_only2.push_back(current_object2);
					MDboth.addObject(MD1.getObject());
					break;
				}
			}
			else if (EMDL::isInt(label1))
			{
				MD2.getValue(label1, myint2);
				if ( ABS(myint2 - myint1) <= ROUND(eps) )
				{
					have_in_2 = true;
					to_remove_from_only2.push_back(current_object2);
					MDboth.addObject(MD1.getObject());
					break;
				}
			}
			else if (EMDL::isDouble(label1))
			{
				MD2.getValue(label1, myd2);
				if (label2 != EMDL_UNDEFINED)
					MD2.getValue(label2, mydy2);
				if (label3 != EMDL_UNDEFINED)
					MD2.getValue(label3, mydz2);

				double dist = sqrt( (myd1 - myd2) * (myd1 - myd2) +
						            (mydy1 - mydy2) * (mydy1 - mydy2) +
						            (mydz1 - mydz2) * (mydz1 - mydz2) );
				if ( ABS(dist) <= eps )
				{
					have_in_2 = true;
					to_remove_from_only2.push_back(current_object2);
					//std::cerr << " current_object1= " << current_object1 << std::endl;
					//std::cerr << " myd1= " << myd1 << " myd2= " << myd2 << " mydy1= " << mydy1 << " mydy2= " << mydy2 << " dist= "<<dist<<std::endl;
					//std::cerr << " to be removed current_object2= " << current_object2 << std::endl;
					MDboth.addObject(MD1.getObject());
					break;
				}
			}
		}

		if (!have_in_2)
		{
			MDonly1.addObject(MD1.getObject());
		}
	}



	for (long int current_object2 = MD2.firstObject();
				current_object2 != MetaDataTable::NO_MORE_OBJECTS && current_object2 != MetaDataTable::NO_OBJECTS_STORED;
				current_object2 = MD2.nextObject())
	{

		bool to_be_removed = false;
		for (long int i = 0; i < to_remove_from_only2.size(); i++)
		{
			if (to_remove_from_only2[i] == current_object2)
			{
				to_be_removed = true;
				break;
			}
		}
		if (!to_be_removed)
		{
			//std::cerr << " doNOT remove current_object2= " << current_object2 << std::endl;
			MDonly2.addObject(MD2.getObject(current_object2));
		}
	}


}

MetaDataTable combineMetaDataTables(std::vector<MetaDataTable> &MDin)
{
	MetaDataTable MDc;

	if (MDin.size() == 0)
    {
		REPORT_ERROR("combineMetaDataTables ERROR: No input STAR files selected!");
    }
	else if (MDin.size() == 1 )
    {
		MDc = MDin[0];
    }
	else
	{
		bool some_labels_missing = false;
		// Find which labels occur in all input tables
		std::vector<EMDLabel> labelsc;
		std::vector<EMDLabel> labels1 = MDin[0].getActiveLabels();

		// Loop over all labels
		for (size_t i = 0; i < labels1.size(); i++)
		{
			// Check their presence in each of the input files
			bool is_present = true;

            for (size_t j = 1; j < MDin.size(); j++)
			{
				is_present = vectorContainsLabel(MDin[j].getActiveLabels(), labels1[i]);

				if (!is_present)
				{
					some_labels_missing = true;
					break;
				}
			}

			if (is_present)
			{
				labelsc.push_back(labels1[i]);
			}
		}

		if (!some_labels_missing)
		{
			// Just append entire tables
			for (size_t j = 0; j < MDin.size(); j++)
			{
				MDc.append(MDin[j]);
			}
		}
		else
		{
            // Select only the labels in common, do this per line!

            for (size_t i = 0; i < labelsc.size(); i++)
            {
                MDc.addLabel(labelsc[i]);
            }

            long totalLines = 0;

            for (size_t i = 0; i < MDin.size(); i++)
			{
                totalLines += MDin[i].numberOfObjects();
            }

            MDc.reserve(totalLines);

            for (size_t i = 0; i < MDin.size(); i++)
            {
                for (size_t j = 0; j < MDin[i].numberOfObjects(); j++)
                {
                    MDc.addValuesOfDefinedLabels(MDin[j].getObject(j));
                }
			}
		}
	}

	return MDc;
}

bool compareLabels(const MetaDataTable &MD1, const MetaDataTable &MD2)
{
	std::vector<EMDLabel> labels1, labels2;

	labels1 = MD1.getActiveLabels();
	labels2 = MD2.getActiveLabels();

	if (labels1.size() != labels2.size())
		return false;

	std::stable_sort(labels1.begin(), labels1.end());
	std::stable_sort(labels2.begin(), labels2.end());

	for (size_t id = 0; id < labels1.size(); id++)
	{
		if (labels1[id] != labels2[id])
			return false;
	}
	return true;
}
