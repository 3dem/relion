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

// Label for sorting
static EMDLabel sortlabel;
bool compareString(MetaDataContainer *lh, MetaDataContainer *rh)
{
	std::string lhstr, rhstr;
	lh->getValue(sortlabel, lhstr);
	rh->getValue(sortlabel, rhstr);
	return lhstr < rhstr;
}
bool compareStringNoAt(MetaDataContainer *lh, MetaDataContainer *rh)
{
	std::string lhstr, rhstr;
	lh->getValue(sortlabel, lhstr);
	rh->getValue(sortlabel, rhstr);
	lhstr = lhstr.substr(lhstr.find("@")+1);
	rhstr = rhstr.substr(rhstr.find("@")+1);
	return lhstr < rhstr;
}
bool compareDouble(MetaDataContainer *lh, MetaDataContainer *rh)
{
	RFLOAT lhstr, rhstr;
	lh->getValue(sortlabel, lhstr);
	rh->getValue(sortlabel, rhstr);
	return lhstr < rhstr;
}
bool compareInteger(MetaDataContainer *lh, MetaDataContainer *rh)
{
	int lhstr, rhstr;
	lh->getValue(sortlabel, lhstr);
	rh->getValue(sortlabel, rhstr);
	return lhstr < rhstr;
}
bool compareLong(MetaDataContainer *lh, MetaDataContainer *rh)
{
	long lhstr, rhstr;
	lh->getValue(sortlabel, lhstr);
	rh->getValue(sortlabel, rhstr);
	return lhstr < rhstr;
}

void MetaDataTable::newSort(const EMDLabel label, bool do_reverse, bool do_sort_after_at)
{

	sortlabel = label;
	if (EMDL::isString(label))
	{
		if (do_sort_after_at)
			std::stable_sort(objects.begin(), objects.end(), compareStringNoAt);
		else
			std::stable_sort(objects.begin(), objects.end(), compareString);
	}
	else if (EMDL::isDouble(label))
		std::stable_sort(objects.begin(), objects.end(), compareDouble);
	else if (EMDL::isInt(label))
		std::stable_sort(objects.begin(), objects.end(), compareInteger);
	else if (EMDL::isLong(label))
		std::stable_sort(objects.begin(), objects.end(), compareLong);
	else
		REPORT_ERROR("Cannot sort this label: " + EMDL::label2Str(label));

	if (do_reverse)
		std::reverse(objects.begin(), objects.end());

}


MetaDataTable::MetaDataTable()
{
    clear();
}

MetaDataTable::MetaDataTable(const MetaDataTable &MD)
{
    clear();
    this->setComment(MD.getComment());
    this->setName(MD.getName());
    this->isList = MD.isList;
    this->activeLabels = MD.activeLabels;
    this->objects.clear();
    this->objects.resize(MD.objects.size());
    for (unsigned long int idx = 0; idx < MD.objects.size(); idx++)
    {
    	 //long int idx = this->addObject();
    	this->objects[idx] = new MetaDataContainer(*(MD.objects[idx]));
    }
	current_objectID = 0;

}

MetaDataTable& MetaDataTable::operator =(const MetaDataTable &MD)
{
    if (this != &MD)
    {
        clear();
        this->setComment(MD.getComment());
        this->setName(MD.getName());
        this->isList = MD.isList;
        this->activeLabels = MD.activeLabels;
        this->objects.resize(MD.objects.size());
        for (long int idx = 0; idx < MD.objects.size(); idx++)
        {
        	//long int idx = this->addObject();
        	this->objects[idx] = new MetaDataContainer(*(MD.objects[idx]));
        }
    	current_objectID = 0;
    }
    return *this;
}

void MetaDataTable::setIsList(bool is_list)
{
    isList = is_list;
}

MetaDataTable::~MetaDataTable()
{
    clear();
}

bool MetaDataTable::isEmpty() const
{
    return (objects.size()==0);
}

long int MetaDataTable::numberOfObjects() const
{
	return objects.size();
}

void MetaDataTable::clear()
{
    for (long int i = 0; i < objects.size(); i++)
    {
    	objects[i]->clear();
    	delete objects[i];
    }

    objects.clear();
    comment.clear();
    name.clear();

    current_objectID = -1;
    activeLabels.clear();
    ignoreLabels.clear();

    isList = false;

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

bool MetaDataTable::setValueFromString(const EMDLabel &label, const std::string &value,
                        long int objectID)
{

	if (objectID == -1)
		objectID = current_objectID;

	if (objectID >= objects.size())
		REPORT_ERROR("MetaDataTable::setValueFromString: objectID >= objects.size()");

	objects[objectID]->addValueFromString(label, value);

	return true;

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
    if (containsLabel(label))
    {
    	std::vector<EMDLabel>::iterator location;
    	location = std::find(activeLabels.begin(), activeLabels.end(), label);
    	activeLabels.erase(location);
    }
}

void MetaDataTable::append(MetaDataTable &app)
{

	bool some_labels_missing = false;
	// Find which labels occur in all input tables
	std::vector<EMDLabel> mylabels = (*this).getActiveLabels();
	std::vector<EMDLabel> labelsin = app.getActiveLabels();

	// Loop over all labels to see whether they are the same in the two tables
	bool have_all_labels = false;
	if (mylabels.size() == 0)
	{
		// empty curent one, just append
		have_all_labels = true;
	}
	else if (mylabels.size() == labelsin.size())
	{
		// non-empty current one: need same number of labels and all labels overlap
		have_all_labels = true;
		for (size_t i = 0; i < mylabels.size(); i++)
		{
			if (!vectorContainsLabel(labelsin, mylabels[i]))
			{
				have_all_labels = false;
				break;
			}
		}
	}

	if (!have_all_labels)
		REPORT_ERROR("ERROR in appending metadata tables with not the same columns!");

	// Go to the end of the table
	current_objectID = objects.size();
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(app)
	{
		addObject(app.getObject());
	}
	// Reset pointer to the beginning of the table
	current_objectID = 0;

}

bool MetaDataTable::addLabel(const EMDLabel label)
{
    if (containsLabel(label))
        return false;
    activeLabels.push_back(label);
    return true;
}

long int MetaDataTable::addObject(MetaDataContainer * data, long int objectID)
{
    long int result;

    if (objectID == -1)
    {
    	result = objects.size();
        objects.resize(result+1);
    }
    else
    {
        if (objectID >= objects.size())
        	REPORT_ERROR("MetaDataTable::addObject: objectID >= objects.size()");

        result = objectID;
        // First free memory of old container if it existed
        if (objects[result])
        {
        	objects[result]->clear();
        	delete objects[result];
        }
    }

    if (data == NULL)
        objects[result] = new MetaDataContainer();
    else
        objects[result] = new MetaDataContainer(*(data));

    // Set iterator pointing to the newly added object
    current_objectID = result;

    // Set default values for the existing labels
    if (data == NULL)
    {
    	/* Sjors 21nov2014: I don't see why this would be necessary...
    	std::vector<EMDLabel>::iterator It;
		for (It = activeLabels.begin(); It != activeLabels.end(); It++)
		{
			(objects[result])->addDefaultValue(*It);
		}
		*/
    }
    else
    {
    	// Set all the labels from the data MDC as active
    	std::vector<EMDLabel> newlabels;
    	newlabels = data->getLabels();
    	for (int i = 0; i < newlabels.size(); i++)
    	{
    		std::vector<EMDLabel>::iterator location;
    		EMDLabel label = newlabels[i];
    		location = std::find(activeLabels.begin(), activeLabels.end(), label);
    		if (location == activeLabels.end())
    		{
    			activeLabels.push_back(label);

                // Add this label to the rest of the objects in this class
                for (long int idx = 0; idx < objects.size(); idx++ )
                {
                    if (idx != result)
                    	objects[idx]->addDefaultValue(label);
                }
    		}
    	}
    }

    return result;
}

long int MetaDataTable::removeObject(long int objectID)
{
	long int i = (objectID == -1) ? current_objectID : objectID;

    if (objects[i])
    {
    	objects[i]->clear();
        delete objects[i];
    }
    objects.erase(objects.begin() + i);

    return lastObject();
}

MetaDataContainer * MetaDataTable::getObject(const long int objectID) const
{
    if (isEmpty())
    {
        // The objects map is empty, error
        REPORT_ERROR("Requested objectID not found (no objects stored). Exiting... ");
    }

    MetaDataContainer * aux;
    if (objectID == -1)
        aux = objects[current_objectID];
    else
    {
#ifdef DEBUG_CHECKSIZES
		if (objectID >= objects.size())
		{
			std::cerr<< "objectID= "<<objectID<<" objects.size()= "<< objects.size() <<std::endl;
			REPORT_ERROR("MetaDataTable::getObject: objectID >= objects.size()");
		}
#endif
    	aux = objects[objectID];
    }

    if (aux == NULL)
    {
        // This objectID does not exist, finish execution
        REPORT_ERROR("Requested objectID not found. Exiting... ");
    }

    return aux;
}

void MetaDataTable::setObject(MetaDataContainer * data, long int objectID)
{

	long int idx = (objectID == -1) ? current_objectID : objectID;

#ifdef DEBUG_CHECKSIZES
	if (idx >= objects.size())
		REPORT_ERROR("MetaDataTable::setObject: idx >= objects.size()");
#endif

        // First delete old container if it exists
        if (this->objects[idx])
        {
        	this->objects[idx]->clear();
        	delete this->objects[idx];
        }
	this->objects[idx] = new MetaDataContainer(*data);

	// Set all the labels from the data MDC as active
	std::vector<EMDLabel>::iterator location;
	std::vector<EMDLabel> newlabels;
	newlabels = data->getLabels();
	for (int i = 0; i < newlabels.size(); i++)
	{
		EMDLabel label = newlabels[i];
		location = std::find(activeLabels.begin(), activeLabels.end(), label);
		if (location == activeLabels.end())
		{
			activeLabels.push_back(label);

            // Add this label with default values to the rest of the objects in this class
            for (long int idx2 = 0; idx2 < objects.size(); idx2++ )
            {
               if (idx2 != idx)
              	objects[idx2]->addDefaultValue(label);
            }
		}
	}


}

long int MetaDataTable::firstObject()
{
    long int result = 0;

    if (!isEmpty())
    {
    	current_objectID = 0;
        result = 0;
    }
    else
    {
        result = NO_OBJECTS_STORED; // Map is empty
    }

    return result;
}

long int MetaDataTable::nextObject()
{
    long int result = 0;

    if (!isEmpty())
    {
    	current_objectID++;

        if (current_objectID < objects.size())
        {
            result = current_objectID;
        }
        else
        {
            result = NO_MORE_OBJECTS;
            current_objectID = lastObject();
        }
    }
    else
    {
        result = NO_OBJECTS_STORED;
        current_objectID = -1;
    }

    return result;
}


long int MetaDataTable::lastObject()
{
    long int result = 0;

    if (!isEmpty())
    {
        result = objects.size() - 1;
        current_objectID = result;
    }
    else
    {
        result = NO_OBJECTS_STORED;
        current_objectID = -1;
    }

    return result;
}


long int MetaDataTable::goToObject(long int objectID)
{
	if (objectID < objects.size())
	{
		current_objectID = objectID;
		return current_objectID;
	}
	else
	{
		REPORT_ERROR("MetaDataTable::goToObject: objectID >= objects.size()");
	}
}

long int MetaDataTable::readStarLoop(std::ifstream& in, std::vector<EMDLabel> *desiredLabels, std::string grep_pattern, bool do_only_count)
{
	setIsList(false);

	//Read column labels
    int labelPosition = 0;
    EMDLabel label;
    std::string line, token, value;

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
			token = line.substr(line.find("_") + 1, line.find("#") - 2);
			label = EMDL::str2Label(token);
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
				activeLabels.push_back(label);

			labelPosition++;
		}
		else // found first data line
		{
			break;
		}
    }

    // Then fill the table (dont read another line until the one from above has been handled)
    bool is_first= true;
    long int nr_objects = 0;
    while (is_first || getline(in, line, '\n'))
    {
		is_first=false;

    	if (grep_pattern == "" || line.find(grep_pattern) !=  std::string::npos )
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
    long int objectID = addObject();
    EMDLabel label;
    std::string line, firstword, value;
    std::vector<std::string> words;
    std::stringstream ss;
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
    		 label = EMDL::str2Label(firstword.substr(1)); // get rid of leading underscore
        	 if (words.size() != 2)
        		 REPORT_ERROR("MetaDataTable::readStarList: did not encounter a single word after "+firstword);
    		 value = words[1];

    		 if (desiredLabels != NULL && !vectorContainsLabel(*desiredLabels, label))
				label = EMDL_UNDEFINED; //ignore if not present in desiredLabels
    		 if (label != EMDL_UNDEFINED)
			 {
				 activeLabels.push_back(label);
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
        REPORT_ERROR( (std::string) "MetaDataTable::read: File " + fn_read + " does not exists" );

    FileName ext = filename.getFileFormat();
    if (ext =="star")
    {
        //REPORT_ERROR("readSTAR not implemented yet...");
        return readStar(in, name, desiredLabels, grep_pattern, do_only_count);
    }
    else
    {
        REPORT_ERROR("MetaDataTable::read ERROR: metadatatable should have .star extension");
    }

    in.close();

}

void MetaDataTable::write(std::ostream& out)
{

    // Only write tables that have something in them
    if (isEmpty())
        return;

    std::vector<EMDLabel>::iterator strIt;
    std::string entryComment;

    out << "\n";
    out << "data_" << getName() <<"\n";
    if (containsComment())
    	out << "# "<< comment << "\n";
    out << "\n";

    if (!isList)
    {
        // Write loop header structure
    	out << "loop_ \n";
        int ii = 0;
    	for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
        {
            ii++;
            if (*strIt != EMDL_COMMENT && *strIt != EMDL_SORTED_IDX) // EMDL_SORTED_IDX is only for internal use, never write it out!
    		//if (*strIt != EMDL_COMMENT)
            {
                out << "_" << EMDL::label2Str(*strIt) << " #" << ii << " \n";
            }
        }

        // Write actual data block
        for (long int idx = 0; idx < objects.size(); idx++)
        {
        	entryComment = "";
        	for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
            {
        		if (*strIt != EMDL_COMMENT && *strIt != EMDL_SORTED_IDX)
            	//if (*strIt != EMDL_COMMENT)
                {
                    out.width(10);
                    objects[idx]->writeValueToStream(out, *strIt);
                    out << " ";
                }
            	if (*strIt == EMDL_COMMENT)
                {
                	objects[idx]->getValue(EMDL_COMMENT, entryComment);
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
        MetaDataContainer * object = getObject();

        entryComment = "";
        int maxWidth=10;
        for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
        {
            if (*strIt != EMDL_COMMENT)
            {
                int w=EMDL::label2Str(*strIt).length();
                if (w>maxWidth)
                    maxWidth=w;
            }
            else
            	object->getValue(EMDL_COMMENT, entryComment);
        }

        for (strIt = activeLabels.begin(); strIt != activeLabels.end(); strIt++)
        {
            if (*strIt != EMDL_COMMENT)
            {
            	int w = EMDL::label2Str(*strIt).length();
            	out << "_" << EMDL::label2Str(*strIt) << std::setw(12 + maxWidth - w) << " ";
                object->writeValueToStream(out, *strIt);
                out << "\n";
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

void MetaDataTable::write(const FileName &fn_out)
{
    std::ofstream  fh;
    fh.open((fn_out).c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR( (std::string)"MetaDataTable::write Cannot write to file: " + fn_out);
    write(fh);
    fh.close();

}

void MetaDataTable::writeValueToString(std::string & result,
                                  const std::string &inputLabel)
{
    MetaDataContainer * aux = getObject();
    aux->writeValueToString(result, EMDL::str2Label(inputLabel));
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

    RFLOAT mydbl;
    int myint;
    long int mylong;
	double xval, yval;
	for (long int idx = 0; idx < objects.size(); idx++)
    {
    	if (EMDL::isDouble(xaxis))
    	{
    		objects[idx]->getValue(xaxis, mydbl);
    		xval = mydbl;
    	}
    	else if (EMDL::isInt(xaxis))
    	{
    		objects[idx]->getValue(xaxis, myint);
    		xval = myint;
    	}
    	else if (EMDL::isLong(xaxis))
    	{
    		objects[idx]->getValue(xaxis, mylong);
    		xval = mylong;
    	}
    	else
    		REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: can only plot xaxis double, int or long int");

    	if (EMDL::isDouble(yaxis))
    	{
    		objects[idx]->getValue(yaxis, mydbl);
    		yval = mydbl;
    	}
    	else if (EMDL::isInt(yaxis))
    	{
    		objects[idx]->getValue(yaxis, myint);
    		yval = myint;
    	}
    	else if (EMDL::isLong(yaxis))
    	{
    		objects[idx]->getValue(yaxis, mylong);
    		yval = mylong;
    	}
    	else
    		REPORT_ERROR("MetaDataTable::addToCPlot2D ERROR: can only plot yaxis double, int or long int");

    	CDataPoint point(xval, yval);
		dataSet.AddDataPoint(point);

	}

    plot2D->AddDataSet(dataSet);

}

void compareMetaDataTable(MetaDataTable &MD1, MetaDataTable &MD2,
		MetaDataTable &MDboth, MetaDataTable &MDonly1, MetaDataTable &MDonly2,
		EMDLabel label1, RFLOAT eps, EMDLabel label2, EMDLabel label3)
{
	if (!MD1.containsLabel(label1))
		REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD1 does not contain the specified label1.");
	if (!MD2.containsLabel(label1))
		REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD2 does not contain the specified label1.");

	if (label2 != EMDL_UNDEFINED)
	{
		if (!EMDL::isDouble(label1) || !EMDL::isDouble(label2))
			REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR 2D or 3D distances are only allowed for RFLOATs.");
		if (!MD1.containsLabel(label2))
			REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD1 does not contain the specified label2.");
		if (!MD2.containsLabel(label2))
			REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD2 does not contain the specified label2.");
	}

	if (label3 != EMDL_UNDEFINED)
	{
		if (!EMDL::isDouble(label3))
			REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR 3D distances are only allowed for RFLOATs.");
		if (!MD1.containsLabel(label3))
			REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD1 does not contain the specified label3.");
		if (!MD2.containsLabel(label3))
			REPORT_ERROR("compareMetaDataTableEqualLabel::ERROR MD2 does not contain the specified label3.");
	}

	MDboth.clear();
	MDonly1.clear();
	MDonly2.clear();

	std::string mystr1, mystr2;
	int myint1, myint2;
	RFLOAT myd1, myd2, mydy1 = 0., mydy2 = 0., mydz1 = 0., mydz2 = 0.;


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
			REPORT_ERROR("compareMetaDataTableEqualLabel ERROR: only implemented for strings, integers or RFLOATs");

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

				RFLOAT dist = sqrt( (myd1 - myd2) * (myd1 - myd2) +
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
		REPORT_ERROR("combineMetaDataTables ERROR: No input STAR files selected!");
	else if (MDin.size() == 1 )
		MDc = MDin[0];
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
			for (size_t j = 0; j < MDin.size(); j++)
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
			for (size_t j = 0; j < MDin.size(); j++)
			{
				FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin[j])
				{
					MetaDataContainer *aux =  MDin[j].getObject();
					aux->keepOnlyLabels(labelsc);
					MDc.addObject(aux);
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

	// DEBUG
	//std::cout << " label1.size() & label2.size() = " << labels1.size() << ", " << labels2.size() << std::endl;
	if (labels1.size() != labels2.size())
		return false;

	// DEBUG
	//std::cout << " before_stable_sort: labels1[]: " << std::flush;
	//for (size_t id = 0; id < labels1.size(); id++)
	//	std::cout << labels1[id] << ", " << std::flush;
	//std::cout << std::endl;
	//std::cout << " before_stable_sort: labels2[]: " << std::flush;
	//for (size_t id = 0; id < labels2.size(); id++)
	//	std::cout << labels2[id] << ", " << std::flush;
	//std::cout << std::endl;

	std::stable_sort(labels1.begin(), labels1.end());
	std::stable_sort(labels2.begin(), labels2.end());

	// DEBUG
	//std::cout << "  after_stable_sort: labels1[]: " << std::flush;
	//for (size_t id = 0; id < labels1.size(); id++)
	//	std::cout << labels1[id] << ", " << std::flush;
	//std::cout << std::endl;
	//std::cout << "  after_stable_sort: labels2[]: " << std::flush;
	//for (size_t id = 0; id < labels2.size(); id++)
	//	std::cout << labels2[id] << ", " << std::flush;
	//std::cout << std::endl;

	for (size_t id = 0; id < labels1.size(); id++)
	{
		if (labels1[id] != labels2[id])
			return false;
	}
	return true;
};
