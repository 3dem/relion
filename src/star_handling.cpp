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
#include "src/star_handling.h"

void combineStarfiles(std::vector<FileName> &fns_in,
                      MetaDataTable &MDout, ObservationModel &obsModel,
                      FileName fn_check, bool do_ignore_optics, std::string tablename_in, int verb)
{
    MetaDataTable MDin0;
    std::vector<MetaDataTable> MDsin, MDoptics;
    std::vector<ObservationModel> obsModels;
    ObservationModel myobsModel0;

    // Read the first table into the global obsModel
    if (do_ignore_optics) MDin0.read(fns_in[0], tablename_in);
    else ObservationModel::loadSafely(fns_in[0], obsModel, MDin0, "discover", 1);
    MDsin.push_back(MDin0);

    // Read all the rest of the tables into local obsModels
    for (int i = 1; i < fns_in.size(); i++)
    {
        ObservationModel myobsModel;
        MetaDataTable MDin; // define again, as reading from previous one may linger here...
        if (do_ignore_optics) MDin.read(fns_in[i], tablename_in);
        else ObservationModel::loadSafely(fns_in[i], myobsModel, MDin, "discover", 1);
        MDsin.push_back(MDin);
        obsModels.push_back(myobsModel);
    }

    // Combine optics groups with the same EMDL_IMAGE_OPTICS_GROUP_NAME, make new ones for those with a different name
    if (!do_ignore_optics)
    {
        std::vector<std::string> optics_group_uniq_names;

        // Initialise optics_group_uniq_names with the first table
        FOR_ALL_OBJECTS_IN_METADATA_TABLE(obsModel.opticsMdt)
        {
            std::string myname;
            obsModel.opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP_NAME, myname);
            optics_group_uniq_names.push_back(myname);
        }

        // Now check uniqueness of the other tables
        for (int MDs_id = 1; MDs_id < fns_in.size(); MDs_id++)
        {
            const int obs_id = MDs_id - 1;

            std::vector<int> new_optics_groups;
            FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsin[MDs_id])
            {
                int tmp;
                MDsin[MDs_id].getValue(EMDL_IMAGE_OPTICS_GROUP, tmp);
                new_optics_groups.push_back(tmp);
            }

            MetaDataTable unique_opticsMdt;
            unique_opticsMdt.addMissingLabels(&obsModels[obs_id].opticsMdt);

            FOR_ALL_OBJECTS_IN_METADATA_TABLE(obsModels[obs_id].opticsMdt)
            {
                std::string myname;
                int my_optics_group;
                obsModels[obs_id].opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP_NAME, myname);
                obsModels[obs_id].opticsMdt.getValue(EMDL_IMAGE_OPTICS_GROUP, my_optics_group);

                // Check whether this name is unique
                bool is_uniq = true;
                int new_group;
                for (new_group = 0; new_group < optics_group_uniq_names.size(); new_group++)
                {
                    if (optics_group_uniq_names[new_group] == myname)
                    {
                        is_uniq = false;
                        break;
                    }
                }
                new_group ++; // start counting of groups at 1, not 0!

                if (is_uniq)
                {
                    if (verb > 0) std::cout << " + Adding new optics_group with name: " << myname << std::endl;

                    optics_group_uniq_names.push_back(myname);
                    // Add the line to the global obsModel
                    obsModels[obs_id].opticsMdt.setValue(EMDL_IMAGE_OPTICS_GROUP, new_group);

                    unique_opticsMdt.addObject();
                    unique_opticsMdt.setObject(obsModels[obs_id].opticsMdt.getObject());
                }
                else
                {
                    if (verb > 0) std::cout << " + Joining optics_groups with the same name: " << myname << std::endl;
                    if (verb > 0) std::cerr << " + WARNING: if these are different data sets, you might want to rename optics groups instead of joining them!" << std::endl;
                    if (verb > 0) std::cerr << " + WARNING: if so, manually edit the rlnOpticsGroupName column in the optics_groups table of your input STAR files." << std::endl;
                }

                if (my_optics_group != new_group)
                {
                    if (verb > 0) std::cout << " + Renumbering group " << myname << " from " << my_optics_group << " to " << new_group << std::endl;
                }

                // Update the optics_group entry for all particles in the MDsin
                for (long int current_object2 = MDsin[MDs_id].firstObject();
                     current_object2 < MDsin[MDs_id].numberOfObjects() && current_object2 >= 0;
                     current_object2 = MDsin[MDs_id].nextObject())
                {
                    int old_optics_group;
                    MDsin[MDs_id].getValue(EMDL_IMAGE_OPTICS_GROUP, old_optics_group, current_object2);
                    if (old_optics_group == my_optics_group)
                        new_optics_groups[current_object2] = new_group;
                }
            }

            obsModels[obs_id].opticsMdt = unique_opticsMdt;

            FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsin[MDs_id])
            {
                MDsin[MDs_id].setValue(EMDL_IMAGE_OPTICS_GROUP, new_optics_groups[current_object]);

                // Also rename the rlnGroupName to not have groups overlapping from different optics groups
                std::string name;
                if (MDsin[MDs_id].getValue(EMDL_MLMODEL_GROUP_NAME, name))
                {
                    name = "optics"+integerToString(new_optics_groups[current_object])+"_"+name;
                    MDsin[MDs_id].setValue(EMDL_MLMODEL_GROUP_NAME, name);
                }
            }
        }

        // Make one vector for combination of the optics tables
        MDoptics.push_back(obsModel.opticsMdt);
        for (int i = 1; i < fns_in.size(); i++)
        {
            MDoptics.push_back(obsModels[i - 1].opticsMdt);
        }

        // Check if anisotropic magnification and/or beam_tilt are present in some optics groups, but not in others.
        // If so, initialise the others correctly
        bool has_beamtilt = false, has_not_beamtilt = false;
        bool has_anisomag = false, has_not_anisomag = false;
        bool has_odd_zernike = false, has_not_odd_zernike = false;
        bool has_even_zernike = false, has_not_even_zernike = false;
        bool has_ctf_premultiplied = false, has_not_ctf_premultiplied = false;
        for (int i = 0; i < fns_in.size(); i++)
        {
            if (MDoptics[i].containsLabel(EMDL_IMAGE_BEAMTILT_X) ||
                MDoptics[i].containsLabel(EMDL_IMAGE_BEAMTILT_Y))
            {
                has_beamtilt = true;
            }
            else
            {
                has_not_beamtilt = true;
            }
            if (MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_00) &&
                MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_01) &&
                MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_10) &&
                MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_11))
            {
                has_anisomag = true;
            }
            else
            {
                has_not_anisomag = true;
            }
            if (MDoptics[i].containsLabel(EMDL_IMAGE_ODD_ZERNIKE_COEFFS))
            {
                has_odd_zernike = true;
            }
            else
            {
                has_not_odd_zernike = true;
            }
            if (MDoptics[i].containsLabel(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS))
            {
                has_even_zernike = true;
            }
            else
            {
                has_not_even_zernike = true;
            }
            if (MDoptics[i].containsLabel(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED))
            {
                has_ctf_premultiplied = true;
            }
            else
            {
                has_not_ctf_premultiplied = true;
            }
        }
#ifdef DEBUG
        printf("has_beamtilt = %d, has_not_beamtilt = %d, has_anisomag = %d, has_not_anisomag = %d, has_odd_zernike = %d, has_not_odd_zernike = %d, has_even_zernike = %d, has_not_even_zernike = %d, has_ctf_premultiplied = %d, has_not_ctf_premultiplied = %d\n", has_beamtilt, has_not_beamtilt, has_anisomag, has_not_anisomag, has_odd_zernike, has_not_odd_zernike, has_even_zernike, has_not_even_zernike, has_ctf_premultiplied, has_not_ctf_premultiplied);
#endif

        for (int i = 0; i < fns_in.size(); i++)
        {
            if (has_beamtilt && has_not_beamtilt)
            {
                if (!MDoptics[i].containsLabel(EMDL_IMAGE_BEAMTILT_X))
                {
                    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDoptics[i])
                    {
                        MDoptics[i].setValue(EMDL_IMAGE_BEAMTILT_X, 0.);
                    }
                }
                if (!MDoptics[i].containsLabel(EMDL_IMAGE_BEAMTILT_Y))
                {
                    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDoptics[i])
                    {
                        MDoptics[i].setValue(EMDL_IMAGE_BEAMTILT_Y, 0.);
                    }
                }
            }

            if (has_anisomag && has_not_anisomag)
            {
                if (!(MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_00) &&
                      MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_01) &&
                      MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_10) &&
                      MDoptics[i].containsLabel(EMDL_IMAGE_MAG_MATRIX_11)) )
                {
                    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDoptics[i])
                    {
                        MDoptics[i].setValue(EMDL_IMAGE_MAG_MATRIX_00, 1.);
                        MDoptics[i].setValue(EMDL_IMAGE_MAG_MATRIX_01, 0.);
                        MDoptics[i].setValue(EMDL_IMAGE_MAG_MATRIX_10, 0.);
                        MDoptics[i].setValue(EMDL_IMAGE_MAG_MATRIX_11, 1.);
                    }
                }
            }

            if (has_odd_zernike && has_not_odd_zernike)
            {
                std::vector<RFLOAT> six_zeros(6, 0);
                if (!MDoptics[i].containsLabel(EMDL_IMAGE_ODD_ZERNIKE_COEFFS))
                {
                    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDoptics[i])
                    {
                        MDoptics[i].setValue(EMDL_IMAGE_ODD_ZERNIKE_COEFFS, six_zeros);
                    }
                }
            }

            if (has_even_zernike && has_not_even_zernike)
            {
                std::vector<RFLOAT> nine_zeros(9, 0);
                if (!MDoptics[i].containsLabel(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS))
                {
                    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDoptics[i])
                    {
                        MDoptics[i].setValue(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS, nine_zeros);
                    }
                }
            }

            if (has_ctf_premultiplied && has_not_ctf_premultiplied)
            {
                if (!MDoptics[i].containsLabel(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED))
                {
                    FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDoptics[i])
                    {
                        MDoptics[i].setValue(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, false);
                    }
                }
            }
        }

        // Now combine all optics tables into one
        obsModel.opticsMdt = MetaDataTable::combineMetaDataTables(MDoptics);
    }

    // Combine the particles tables
    MDout = MetaDataTable::combineMetaDataTables(MDsin);

    //Deactivate the group_name column
    MDout.deactivateLabel(EMDL_MLMODEL_GROUP_NO);

    if (fn_check != "")
    {
        EMDLabel label = EMDL::str2Label(fn_check);
        if (!MDout.containsLabel(label))
            REPORT_ERROR("ERROR: the output file does not contain the label to check for duplicates. Is it present in all input files?");

        /// Don't want to mess up original order, so make a MDsort with only that label...
        FileName fn_this, fn_prev = "";
        MetaDataTable MDsort;
        FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDout)
        {
            MDout.getValue(label, fn_this);
            MDsort.addObject();
            MDsort.setValue(label, fn_this);
        }
        // sort on the label
        MDsort.newSort(label);
        long int nr_duplicates = 0;
        FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsort)
        {
            MDsort.getValue(label, fn_this);
            if (fn_this == fn_prev)
            {
                nr_duplicates++;
                std::cerr << " WARNING: duplicate entry: " << fn_this << std::endl;
            }
            fn_prev = fn_this;
        }

        if (nr_duplicates > 0)
            std::cerr << " WARNING: Total number of duplicate "<< fn_check << " entries: " << nr_duplicates << std::endl;
    }

    MDout.setName(MDin0.getName());

}

