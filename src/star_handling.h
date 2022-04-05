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


#ifndef RELION_STAR_HANDLING_H
#define RELION_STAR_HANDLING_H

#include <src/metadata_table.h>
#include <src/filename.h>
#include <src/jaz/single_particle/obs_model.h>

// Read in all the STAR files in the vector; if do_ignore_optics is true, use tablename_in, otherwise fill obsModel
// fn_check is a string for duplication checking, leave empty for no checking
void combineStarfiles(std::vector<FileName> &fns_in,
                      MetaDataTable &MDout, ObservationModel &obsModel,
                      FileName fn_check = "", bool do_ignore_optics = false, std::string tablename_in = "", int verb = 0);


#endif //RELION_STAR_HANDLING_H
