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
#include "src/star_handling.h"
#include <src/jaz/single_particle/obs_model.h>

class ExpTiltImage
{
public:

    // ID of this tiltimage inside its own tiltserie
	long int id;

    // ID of the tiltserie that this tiltimage came from
	long int tiltseries_id;

    // Optics group of the image
    int optics_group;

    // Empty Constructor
	ExpTiltImage() {}

	// Destructor needed for work with vectors
	~ExpTiltImage() {}

	// Copy constructor needed for work with vectors
	ExpTiltImage(ExpTiltImage const& copy)
	{
		id = copy.id;
        tiltseries_id = copy.tiltseries_id;
        optics_group = copy.optics_group;
	}

	// Define assignment operator in terms of the copy constructor
	ExpTiltImage& operator=(ExpTiltImage const& copy)
	{
		id = copy.id;
		tiltseries_id = copy.tiltseries_id;
        optics_group = copy.optics_group;
		return *this;
	}

};

class ExpTiltSerie
{
public:
    // All tiltimages in the tiltserie
    std::vector<ExpTiltImage> tiltimages;

	// ID of the tiltserie
	long int id;

    // Name of the tiltserie (by this name it will be recognised upon reading)
	std::string name;

    // Observation model holding the data for all optics groups
    ObservationModel obsModel;

    // table with all the metadata
    MetaDataTable MDtiltimages;

	// Empty Constructor
	ExpTiltSerie() {}

	// Destructor needed for work with vectors
	~ExpTiltSerie() {}

	// Copy constructor needed for work with vectors
	ExpTiltSerie(ExpTiltSerie const& copy)
	{
		tiltimages = copy.tiltimages;
        id = copy.id;
		name = copy.name;
        obsModel = copy.obsModel;
        MDtiltimages = copy.MDtiltimages;
	}

	// Define assignment operator in terms of the copy constructor
	ExpTiltSerie& operator=(ExpTiltSerie const& copy)
	{
        tiltimages = copy.tiltimages;
        id = copy.id;
		name = copy.name;
        obsModel = copy.obsModel;
        MDtiltimages = copy.MDtiltimages;
		return *this;
	}
    long int numberOfTiltImages()
    {
        return MDtiltimages.numberOfObjects();
    }

    // Calculate the total number of optics groups in this tiltserie
    int numberOfOpticsGroups()
    {
        return obsModel.numberOfOpticsGroups();
    }

    // Get the pixel size for this optics group
    RFLOAT getOpticsPixelSize(int optics_group);

    // Get the original image size for this optics group
    int getOpticsImageSize(int optics_group);


};


class TomographyExperiment
{
public:
    // All tiltseries in the experiment
    std::vector<ExpTiltSerie> tiltseries;

    // Metadata table with information about all tiltseries
    MetaDataTable MDtiltseries;

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
        tiltseries.clear();
        MDtiltseries.clear();
    }

    // Read from file, return false if this is not a TomographyExperiment
    bool read(FileName fn_in, int verb = 0);

    // Write to file
    void write(FileName fn_root);

    // Calculate the total number of tiltseries in this tomography experiment
    long int numberOfTiltseries()
    {
        return tiltseries.size();
    }

    // Calculate the total number of tilt series images in this tomography experiment
    long int numberOfTiltImages();

    // Make one big metadatatable with all movies/micrographs (to be used for motioncorr and ctffind runners)
    void generateSingleMetaDataTable(MetaDataTable &MDout, ObservationModel &obsModel);

    // Convert back from one big metadatatable into separate STAR files for each tilt serie
    void convertBackFromSingleMetaDataTable(MetaDataTable &MDin, ObservationModel &obsModel);
};

#endif //RELION_TOMO_MODEL_H
