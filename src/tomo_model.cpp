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

#include "src/tomo_model.h"


RFLOAT ExpTiltSerie::getOpticsPixelSize(int optics_group)
{
    return obsModel.getPixelSize(optics_group);
}

int ExpTiltSerie::getOpticsImageSize(int optics_group)
{
    return obsModel.getBoxSize(optics_group);
}



bool TomographyExperiment::read(FileName fn_in, int verb)
{

    clear();

    if (!fn_in.isStarFile())
	{
        if (verb > 0) std::cerr << " ERROR: input filename for TomographyExperiment is not a STAR file" << std::endl;
        return false;
    }

    MDtiltseries.read(fn_in, "tilt_series");

    if (MDtiltseries.numberOfObjects() == 0)
    {
         if (verb > 0) std::cerr << " ERROR: input starfile for TomographyExperiment " << fn_in << "  does not contain any entries in tilt_series table" << std::endl;
         return false;
    }

    for (long int ts_id = 0; ts_id < MDtiltseries.numberOfObjects(); ts_id++)
    {
        ExpTiltSerie mytiltserie;
        mytiltserie.id = ts_id;

        FileName fn_star;
        MDtiltseries.getValue(EMDL_TOMO_TILT_SERIES_STARFILE, fn_star, ts_id);
        MDtiltseries.getValue(EMDL_TOMO_NAME, mytiltserie.name, ts_id);

        ObservationModel::loadSafely(fn_star, mytiltserie.obsModel, mytiltserie.MDtiltimages, "movies", verb);
        if (mytiltserie.obsModel.opticsMdt.numberOfObjects() == 0)
        {
           if (verb > 0) std::cerr << " ERROR: input starfile for tilt serie " << fn_star << " does not contain any optics groups" << std::endl;
            return false;
        }

        // Get all the necessary information about this tiltseries
        for (long int timg_id = 0; timg_id < mytiltserie.MDtiltimages.numberOfObjects(); timg_id++)
        {
            ExpTiltImage mytiltimage;
            mytiltimage.id = timg_id;
            mytiltserie.MDtiltimages.getValue(EMDL_IMAGE_OPTICS_GROUP, mytiltimage.optics_group, timg_id);
            mytiltimage.tiltseries_id = ts_id;
            mytiltserie.tiltimages.push_back(mytiltimage);
        }

        tiltseries.push_back(mytiltserie);

    }

    return true;

}


void TomographyExperiment::write(FileName fn_outdir)
{
    // Make user fn_outdir ends in a slah: it should be a directory
    if (fn_outdir[fn_outdir.size()-1]!='/') fn_outdir+='/';

    // Write all the star files with the individual starfiles for the tilt images in the ouput directory
    for (long int ts_id = 0; ts_id < MDtiltseries.numberOfObjects(); ts_id++)
    {

        FileName fn_star;
        MDtiltseries.getValue(EMDL_TOMO_TILT_SERIES_STARFILE, fn_star, ts_id);
        FileName fn_newstar = getOutputFileWithNewUniqueDate(fn_star, fn_outdir);
        MDtiltseries.setValue(EMDL_TOMO_TILT_SERIES_STARFILE, fn_newstar, ts_id);

        // Create output directory if neccesary
        FileName newdir = fn_newstar.beforeLastOf("/");
        if (!exists(newdir))
        {
            mktree(newdir);
        }

        tiltseries[ts_id].obsModel.save(tiltseries[ts_id].MDtiltimages, fn_newstar, "tilt_images");

    }

    // Also write the (now modified with fn_newstars) tilt_series.star file in the root directory
    FileName fn_out = fn_outdir + "tilt_series.star";
    MDtiltseries.setName("tilt_series");
    MDtiltseries.write(fn_out);

}

long int TomographyExperiment::numberOfTiltImages()
{

    long int result = 0;
    for (long int i = 0; i < tiltseries.size(); i++)
    {
        result += tiltseries[i].tiltimages.size();
    }
    return result;

}

void TomographyExperiment::generateSingleMetaDataTable(MetaDataTable &MDout, ObservationModel &obsModel)
{
    std::vector<FileName> fn_stars;
    for (long int ts_id = 0; ts_id < MDtiltseries.numberOfObjects(); ts_id++)
    {
        FileName fn_star;
        MDtiltseries.getValue(EMDL_TOMO_TILT_SERIES_STARFILE, fn_star, ts_id);
        fn_stars.push_back(fn_star);
    }

    combineStarfiles(fn_stars, MDout, obsModel);

}

void TomographyExperiment::convertBackFromSingleMetaDataTable(MetaDataTable &MDin, ObservationModel &obsModel)
{
     for (long int ts_id = 0; ts_id < MDtiltseries.numberOfObjects(); ts_id++)
     {
         tiltseries[ts_id].obsModel = obsModel;
         tiltseries[ts_id].MDtiltimages = subsetMetaDataTable(MDin, EMDL_TOMO_NAME, tiltseries[ts_id].name, false);
         tiltseries[ts_id].obsModel.removeUnusedOpticsGroups(tiltseries[ts_id].MDtiltimages);
     }
}
