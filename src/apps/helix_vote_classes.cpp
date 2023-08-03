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

#include <src/image.h>
#include <src/funcs.h>
#include <src/memory.h>
#include <src/euler.h>
#include <src/time.h>
#include <src/metadata_table.h>
#include <src/args.h>
#include "src/macros.h"

struct Helix
{
	// On which micrograph is this helix
	FileName fn_mic;

	// What is the helical tube id
	int tube_id;

	// What are the xy-coordinates of the start and end points
	RFLOAT start_x, start_y, end_x, end_y;

	// Length of the helix
	RFLOAT length;

	// Count how often each class is observed in this helix
	MultidimArray<RFLOAT> count_classes;

	// My group
	int my_group;

};

class helix_analyse_classes_parameters
{
	public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	FileName fn_in, fn_coord_suffix, fn_coord_list, fn_out, fn_groups, fn_group_names;

	int nr_classes, nr_groups;

	bool do_norm;

	// Minimum score to assign group based on voting
	float voting_threshold, consistency_check_ratio;
    int min_nr_picks, min_picks_group;

	// All helices in the data set
	std::vector<Helix> helices;

	std::vector<std::string> group_names;

	// Map to go from micrograph name to coordinate file name
	std::map<FileName, FileName> micname2coordname;

	std::map<int, int> class2group;

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int gen_section = parser.addSection("General options");
		fn_in = parser.getOption("--i", "The _data.star file with the classes to be analysed");
		nr_classes = textToInteger(parser.getOption("--nr_classes", "Number of classes in the input star file"));
		fn_coord_suffix = parser.getOption("--coord_suffix", "The suffix for the coordinate files, e.g. \"_picked.star\" or \".box\"","");
		fn_coord_list = parser.getOption("--pick", "Alternative to coord_suffix: a 2-column STAR file with micrographs and coordinate files","");
		fn_out = parser.getOption("--o", "Output directory name", "HelixAnalyseClasses/");
		fn_groups = parser.getOption("--groups", "A string with grouping of comma-separated class numbers (with ':' for separation of groups, e.g. 1,4,5:6,2)", "");
		fn_group_names = parser.getOption("--group_names", "A string with :-separated names for all groups (e.g. phf:sf)", "");
		voting_threshold = textToFloat(parser.getOption("--voting_threshold", "Minimum fraction to assign a helix to a group", "0.0"));
		consistency_check_ratio = textToFloat(parser.getOption("--consistency_check", "Check this fraction of particles to be on line of start-end coordinates", "0.05"));
		min_nr_picks = textToInteger(parser.getOption("--min_nr_picks", "Select filaments that have at least this many picks for the group indicated below", "-1"));
        min_picks_group = textToInteger(parser.getOption("--min_picks_group", "Number of the group (first one is 1) to select based on minimum number of particles", "-1"));
        do_norm = parser.checkOption("--norm", "Perform normalisation before voting?");
		verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

	    if (parser.checkForErrors(verb))
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	}

	bool is_point_on_line(RFLOAT px, RFLOAT py, RFLOAT start_x, RFLOAT start_y, RFLOAT end_x, RFLOAT end_y)
	{
	    RFLOAT dx = end_x - start_x;
	    RFLOAT dy = end_y - start_y;
	    RFLOAT dpx = px - start_x;
	    RFLOAT dpy = py - start_y;
	    RFLOAT dd = fabs( (dx*dpy) - (dpx*dy) ) / sqrt( (dx*dx) + (dy*dy) );
	    RFLOAT dotp = dpx*dx + dpy*dy;
	    return (dd < 2 && 0 <= dotp && dotp  <= dx*dx + dy*dy);

	}


	void initialise()
	{

        if (verb>0) std::cout << " Initialising..." << std::endl;

		init_random_generator();

		if (fn_out[fn_out.length()-1] != '/') fn_out += "/";

        // This code is copied from preprocessing.cpp...
		// Either get coordinate filenames from coord_list, or from the fn_coord_suffix
		if (fn_coord_list != "")
		{
			MetaDataTable MDcoords;
			MDcoords.read(fn_coord_list);
			if (!MDcoords.containsLabel(EMDL_MICROGRAPH_NAME))
				REPORT_ERROR("ERROR: coordinates list star file does not contain rlnMicrographName label");
			// Make sure picks are sorted on micrographname
			MDcoords.newSort(EMDL_MICROGRAPH_NAME);
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDcoords)
			{
				FileName fn_mic, fn_coord;
				MDcoords.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
				MDcoords.getValue(EMDL_MICROGRAPH_COORDINATES, fn_coord);
				micname2coordname[fn_mic] = fn_coord;
			}

			if (fn_out == "") fn_out = fn_coord_list.beforeLastOf("/") + "/";
		}
		else if (fn_coord_suffix != "")
		{

			FileName fn_mic;
			std::ifstream fh;
			fh.open((fn_coord_suffix).c_str(), std::ios_base::in);
			fh >> fn_mic;
			fh.close();

			MetaDataTable MDmics;
			MDmics.read(fn_mic,"micrographs");
			if (!MDmics.containsLabel(EMDL_MICROGRAPH_NAME))
				REPORT_ERROR("ERROR: input star file with micrographs from coord_suffix does not contain rlnMicrographName label");

			// Make sure micrographs are sorted on micrographname
			MDmics.newSort(EMDL_MICROGRAPH_NAME);

			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDmics)
			{
				FileName fn_mic, fn_pre, fn_jobnr, fn_post;
				MDmics.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
				decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
				FileName fn_coord = fn_coord_suffix.beforeLastOf("/") + "/" + fn_post.withoutExtension() + fn_coord_suffix.afterFirstOf("coords_suffix");
				micname2coordname[fn_mic] = fn_coord;
			}
			if (fn_out == "") fn_out = fn_coord_suffix.beforeLastOf("/") + "/";

		}
		else
			REPORT_ERROR("ERROR: provide --coord_suffix OR --pick");

		if (verb > 0)
		{
			std::cout << " + Reading in all start-end coordinate files ..." << std::endl;
			init_progress_bar(micname2coordname.size());
		}

		// Read in all coordinate files to fill helices vector
		std::map<FileName, FileName>::iterator it;
        int i = 0;
		for ( it = micname2coordname.begin(); it != micname2coordname.end(); it++ )
        {

        	FileName fn_mic = it->first;
        	FileName fn_coord = it->second;
        	MetaDataTable MDcoord;
        	if (exists(fn_coord))
        	{
        		MDcoord.read(fn_coord);
				bool is_start = true;
				int tube_id = 1;
				Helix myhelix;
				FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDcoord)
				{
					RFLOAT xcoord;
					RFLOAT ycoord;
					MDcoord.getValue(EMDL_IMAGE_COORD_X, xcoord);
					MDcoord.getValue(EMDL_IMAGE_COORD_Y, ycoord);
					if (is_start)
					{
						myhelix.fn_mic = fn_mic;
						myhelix.tube_id = tube_id;
						myhelix.start_x = xcoord;
						myhelix.start_y = ycoord;
						myhelix.my_group = -1;
						is_start = false;
						tube_id++;
					}
					else
					{
						myhelix.end_x = xcoord;
						myhelix.end_y = ycoord;
						myhelix.length = sqrt( (myhelix.end_x-myhelix.start_x)*(myhelix.end_x-myhelix.start_x) +
											   (myhelix.end_y-myhelix.start_y)*(myhelix.end_y-myhelix.start_y) );
						helices.push_back(myhelix);
						is_start = true;
					}
				}
        	}
        	i++;
        	if (verb > 0)
        	{
        		progress_bar(i);
        	}
        }

       	// Get mapping from class number to group number
        if (fn_groups != "")
        {
        	std::vector<std::string> groups;
        	int numfound = splitString(fn_groups, ":", groups);
        	nr_groups = groups.size();
        	for (int igroup=0; igroup<groups.size(); igroup++)
        	{
        		std::vector<std::string> classes;
        		int numfound = splitString(groups[igroup], ",", classes);

            	for (int iclass=0; iclass<classes.size(); iclass++)
            	{
            		int myclass = textToInteger(classes[iclass]);
            		class2group[myclass] = igroup+1;
            	}
        	}

        	if (fn_group_names != "")
        	{
        		int numfound = splitString(fn_group_names, ":", group_names);
        		if (group_names.size() != nr_groups) REPORT_ERROR("ERROR: incorrect number of group names");
        	}
        	else
        	{
        		group_names.clear();
        		for (int igroup=0; igroup < groups.size(); igroup++)
        		{
        			std::string myname = "group_" + integerToString(igroup+1);
        			group_names.push_back(myname);
        		}
        	}

            bool has_unclear = false;
        	for (int i=1; i <= nr_classes; i++)
        	{
            	if (class2group.find(i) == class2group.end())
        		{
        			class2group[i] = nr_groups+1;
        			has_unclear = true;
        		}
        	}
        	if (has_unclear)
        	{
        		nr_groups++;
        		group_names.push_back("notgrouped");
        	}

    		if (verb>1)
    		{
				std::cout << " + Mapping different classes to groups: " << std::endl;
    			std::map<int, int>::iterator it;
				for ( it = class2group.begin(); it != class2group.end(); it++ )
				{
					if (fn_group_names != "")
					{
						std::cout << "   - class " << it->first << " -> group " << group_names[it->second-1] << std::endl;
					}
					else
					{
						std::cout << "   - class " << it->first << " -> group " << it->second << std::endl;
					}
			}
    		}

        }
        else
        {
        	group_names.clear();
        	for (int i=1; i <= nr_classes; i++)
        	{
        		class2group[i]=i;
    			std::string myname = "class_" + integerToString(i);
    			group_names.push_back(myname);
        	}
        	nr_groups = nr_classes;
        }

        // Initialise count vectors in all helices
        MultidimArray<RFLOAT> zeros;
        zeros.initZeros(nr_groups);
        for (long int idx = 0; idx < helices.size(); idx++)
        {
        	helices[idx].count_classes = zeros;
        }

        if (verb > 0) std::cout << " + Reading input STAR file ... " << std::endl;

        MetaDataTable MDin;
		MDin.read(fn_in,"particles");
		if (!MDin.containsLabel(EMDL_PARTICLE_HELICAL_TUBE_ID))
			REPORT_ERROR("ERROR: input star file does not contain label for helical tube ID");
		if (!MDin.containsLabel(EMDL_MICROGRAPH_NAME))
			REPORT_ERROR("ERROR: input star file does not contain rlnMicrographName label");
		if (!MDin.containsLabel(EMDL_PARTICLE_CLASS))
			REPORT_ERROR("ERROR: input star file does not contain label for class number");

		// Make sure particles are sorted on micrographname
		MDin.newSort(EMDL_MICROGRAPH_NAME);

		int barstep;
		if (verb > 0)
		{
			init_progress_bar(MDin.numberOfObjects());
			barstep = MDin.numberOfObjects()/100;
		}

		// Now count instances of all groups in all helices
        long int idx = 0, prev_idx = 0, cc = 0;
		FileName fn_mic_old = "";
        FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
        {
        	FileName fn_mic;
        	int tube_id, iclass, mygroup;
        	MDin.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
        	MDin.getValue(EMDL_PARTICLE_HELICAL_TUBE_ID, tube_id);
        	MDin.getValue(EMDL_PARTICLE_CLASS, iclass);
        	mygroup = class2group[iclass];
        	if (mygroup > nr_groups) REPORT_ERROR("ERROR: mygroup>nr_groups");
        	if (fn_mic != fn_mic_old)
        	{
        		prev_idx = idx;
        		fn_mic_old = fn_mic;
        	}
        	bool found = false;
        	for (idx = 0; idx < helices.size(); idx++)
        	{
        		if (helices[idx].fn_mic == fn_mic && helices[idx].tube_id == tube_id)
        		{
        			DIRECT_A1D_ELEM(helices[idx].count_classes, mygroup-1) += 1.;
					found = true;
        			if (rnd_unif() < consistency_check_ratio)
        			{
        				RFLOAT x, y;
        				MDin.getValue(EMDL_IMAGE_COORD_X, x);
        				MDin.getValue(EMDL_IMAGE_COORD_Y, y);
        				if (!is_point_on_line(x, y, helices[idx].start_x, helices[idx].start_y,
        						helices[idx].end_x, helices[idx].end_y))
        				{
        					std::cerr << " fn_mic= " << fn_mic << " tube_id= " << tube_id << std::endl;
        					std::cerr << " x= " << x << " y= " << y << std::endl;
        					std::cerr << " helices[idx].start_x= " << helices[idx].start_x << " helices[idx].start_y= " << helices[idx].start_y << std::endl;
        					std::cerr << " helices[idx].end_x= " << helices[idx].end_x << " helices[idx].end_y= " << helices[idx].end_y << std::endl;
        					REPORT_ERROR("ERROR: point is not on line between start-end coordinates");
        				}

        			}

        			break;
        		}
        	}

        	if (!found)
        	{
        		std::cerr << " fn_mic= " << fn_mic << " tube_id= " << tube_id << " iclass= " << iclass << "helices.size()= "<<helices.size() << std::endl;
        		REPORT_ERROR("ERROR: cannot find helix.. ");
        	}

        	cc++;
    		if (verb > 0 && cc%barstep==0)
    		{
    			progress_bar(cc);
    		}

        }
        progress_bar(MDin.numberOfObjects());


        if (do_norm)
        {
			// Normalize all count_classes vectors
			for (long int idx = 0; idx < helices.size(); idx++)
			{

				RFLOAT sum = helices[idx].count_classes.sum();
				if (sum > 0.)
				{
					helices[idx].count_classes /= sum;
					if (verb > 3)
					{
						std::cout << helices[idx].fn_mic << " id= " << helices[idx].tube_id << " l= " << helices[idx].length << " " << helices[idx].count_classes << std::endl;
					}
				}
			}
        }

        if (verb>0) std::cout << " + Done initialising!" << std::endl;


	}


	void run()
	{
        if (min_nr_picks > 0 && min_picks_group > 0)
        {

            if (verb > 0) std::cout << " + Selecting filaments with at least " << min_nr_picks << " particles for group " << min_picks_group << " ..." << std::endl;

            // Selecting
            for (long int idx = 0; idx < helices.size(); idx++)
            {
                if (helices[idx].count_classes(min_picks_group - 1) >= min_nr_picks)
                {
                    helices[idx].my_group = min_picks_group;
                    if (verb > 2) {
                        std::cout << "group= " << group_names[min_picks_group - 1] << helices[idx].fn_mic << " id= "
                                  << helices[idx].tube_id << " l= " << helices[idx].length << " "
                                  << helices[idx].count_classes << std::endl;
                    }
                }
                else
                {
                    helices[idx].my_group = nr_groups;
                }
            }

            // Write out voting results
            write();
            if (verb > 0) std::cout << " + Done selecting!" << std::endl;

        }
        else
        {

            if (verb > 0) std::cout << " + Voting ..." << std::endl;

            // Voting
            for (long int idx = 0; idx < helices.size(); idx++) {
                long int igr;
                if (helices[idx].count_classes.maxIndex(igr) >= voting_threshold && group_names[igr] != "notgrouped" &&
                    helices[idx].count_classes.sum() > 0.) {
                    if (verb > 2) {
                        std::cout << "group= " << group_names[igr] << helices[idx].fn_mic << " id= "
                                  << helices[idx].tube_id << " l= " << helices[idx].length << " "
                                  << helices[idx].count_classes << std::endl;
                    }
                    helices[idx].my_group = igr + 1;
                }
            }

            // Write out voting results
            write();
            if (verb > 0) std::cout << " + Done voting!" << std::endl;
        }

        /*
        // Make 2D coincidence matrix

        MultidimArray<RFLOAT> CIM(group_names.size(), group_names.size());
        for (long int idx = 0; idx < helices.size(); idx++)
        {

        	for (int i = 0; i < XSIZE(helices[idx].count_classes); i++)
        	{
        		RFLOAT count_i = DIRECT_A1D_ELEM(helices[idx].count_classes, i);
        		for (int j = 0; j < XSIZE(helices[idx].count_classes); j++)
            	{
            		RFLOAT count_j = DIRECT_A1D_ELEM(helices[idx].count_classes, j);
                   	DIRECT_A2D_ELEM(CIM, i, j) += (count_i+count_j);
            		DIRECT_A2D_ELEM(CIM, j, i) += (count_i+count_j);
            	}
        	}
        }

    	for (int i = 0; i < XSIZE(helices[0].count_classes); i++)
    	{

        	for (int j = 0; j < XSIZE(helices[0].count_classes); j++)
        	{
        		std::cout << DIRECT_A2D_ELEM(CIM, j, i) << " ";
        	}
        	std::cout << std::endl;
    	}
    	*/

	}

	void write()
	{

		// tmp
		std::vector<std::string> outgroup_names = group_names;
		//Add 1 to size() to include rest group
		for (int ogr = 0; ogr < outgroup_names.size() + 1; ogr++)
		{
	        MetaDataTable MDcoord, MDpick;
			FileName fn_pick, mygroupname;
	        MDpick.setName("coordinate_files");
	        int my_outgroup;

	        if (ogr == outgroup_names.size())
	        {
	        	mygroupname = "rest";
	        	my_outgroup = -1;
	        }
	        else
	        {
	        	mygroupname = outgroup_names[ogr];
	        	my_outgroup = ogr + 1;
	        }
        	fn_pick = fn_out + mygroupname+ ".star";

			int count_helices = 0;
        	for (long int idx = 0; idx < helices.size(); idx++)
	        {

				if (helices[idx].my_group == my_outgroup)
				{

					MDcoord.addObject();
					MDcoord.setValue(EMDL_IMAGE_COORD_X, helices[idx].start_x);
					MDcoord.setValue(EMDL_IMAGE_COORD_Y, helices[idx].start_y);
					MDcoord.setValue(EMDL_PARTICLE_AUTOPICK_FOM, 0.);
					MDcoord.setValue(EMDL_PARTICLE_CLASS, 0); // Dummy values to avoid problems in JoinStar
					MDcoord.setValue(EMDL_ORIENT_PSI, 0.0);
					MDcoord.addObject();
					MDcoord.setValue(EMDL_IMAGE_COORD_X, helices[idx].end_x);
					MDcoord.setValue(EMDL_IMAGE_COORD_Y, helices[idx].end_y);
					MDcoord.setValue(EMDL_PARTICLE_AUTOPICK_FOM, 0.);
					MDcoord.setValue(EMDL_PARTICLE_CLASS, 0); // Dummy values to avoid problems in JoinStar
					MDcoord.setValue(EMDL_ORIENT_PSI, 0.0);

					count_helices++;
				}

				if (MDcoord.numberOfObjects() > 0 && (idx + 1 == helices.size() || helices[idx+1].fn_mic != helices[idx].fn_mic) )
				{

					// Write out MDcoords and add line to MDpick
					FileName fn_pre, fn_jobnr, fn_post, fn_coord;
					decomposePipelineFileName(helices[idx].fn_mic, fn_pre, fn_jobnr, fn_post);
					fn_coord = fn_out + fn_post.withoutExtension() + "_" + mygroupname + ".star";
					int res = mktree(fn_coord.beforeLastOf("/"));
					MDcoord.write(fn_coord);
					MDpick.addObject();
					MDpick.setValue(EMDL_MICROGRAPH_NAME, helices[idx].fn_mic);
					MDpick.setValue(EMDL_MICROGRAPH_COORDINATES, fn_coord);

					// Re-use MDcoords for next micrograph
					MDcoord.clear();
				}

	        }

			if (MDpick.numberOfObjects() > 0)
			{
				MDpick.write(fn_pick);
				if (verb > 0) std::cout << "  - written out: " << fn_pick << " with " << count_helices << " filaments" << std::endl;
			}
		}


	}


};


int main(int argc, char *argv[])
{
	helix_analyse_classes_parameters prm;

	try
	{
		prm.read(argc, argv);

		prm.initialise();

		prm.run();
	}
	catch (RelionError XE)
	{
		std::cerr << XE;
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}
