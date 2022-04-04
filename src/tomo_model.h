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

#ifndef RELION_TOMO_MODEL_H
#define RELION_TOMO_MODEL_H

#include "src/metadata_table.h"
#include <src/jaz/single_particle/obs_model.h>

class TomographyExperiment
{
public:
    // All groups in the experiment
    std::vector<Tomogram> tomograms;

    // Observation model holding the data for all optics groups
    ObservationModel obsModel;

    // Empty Constructor
    TomographyExperiment()
    {
        clear();
    }

    ~TomographyExperiment()
    {
        clear();
    }

    void clear()
    {
        tomograms.clear();
    }

    // Calculate the total number of tomograms in this tomography experiment
    long int numberOfTomograms()
    {
        return tomograms.size();
    }

    // Calculate the total number of tilt series images in this tomography experiment
    long int numberOfTiltSeriesImages();

    // Calculate the total number of optics groups in this experiment
    int numberOfOpticsGroups()
    {
        return obsModel.numberOfOpticsGroups();
    }

    // Get the pixel size for this optics group
    RFLOAT getOpticsPixelSize(int optics_group);

    // Get the original image size for this optics group
    int getOpticsImageSize(int optics_group);

    // Get the random_subset for this particle
    int getRandomSubset(long int part_id);

    // Get the group_id for the N'th image for this particle
    long int getGroupId(long int part_id, int img_id);

    // Get the optics group to which the N'th image for this particle belongs
    int getOpticsGroup(long int part_id, int img_id);

    // Get the original position in the input STAR file for the N'th image for this particle
    int getOriginalImageId(long int part_id, int img_id);

    // Get the pixel size for the N-th image of this particle
    RFLOAT getImagePixelSize(long int part_id, int img_id);

    // Get the vector of number of images per group_id
    void getNumberOfImagesPerGroup(std::vector<long int> &nr_particles_per_group, int random_subset = 0);

    // Get the vector of number of images per group_id
    void getNumberOfImagesPerOpticsGroup(std::vector<long int> &nr_particles_per_group, int random_subset = 0);

    // Get the metadata-row for this image in a separate MetaDataTable
    MetaDataTable getMetaDataImage(long int part_id, int img_id);

    // Which micrograph (or tomogram) doe this particle image comes from?
    FileName getMicrographName(long int ori_img_id);
    FileName getMicrographName(long int part_id, int img_id);

    // Add a particle
    long int addParticle(std::string part_name, int random_subset = 0);

    // Add an image to the given particle
    int addImageToParticle(long int part_id, std::string img_name, long int ori_img_id, long int group_id,
                           int optics_group, bool unique);

    // Add a group
    long int addGroup(std::string mic_name, int optics_group);

    // Read from file
    void read(FileName fn_in, int verb = 0);

    // Write to file
    void write(FileName fn_root);

};

#endif //RELION_TOMO_MODEL_H
