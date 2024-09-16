#ifndef TOMOGRAM_SET_H
#define TOMOGRAM_SET_H

#include <string>
#include <vector>
#include <src/metadata_table.h>
#include <src/filename.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/ctf.h>
#include <src/transformations.h>

class TomogramSet
{
	public:

        MetaDataTable globalTable;
        std::vector<MetaDataTable> tomogramTables;

        TomogramSet();
        TomogramSet(FileName filename, bool verbose = true);

        // return false if this is not a TomogramSet
        bool read(FileName filename, bool verbose = true);
        void write(FileName filename);

        void removeTomogram(std::string tomogramName);

        // If max_dose is positive, then only images with cumulativeDose less than or equal to max_dose will be loaded.
		Tomogram loadTomogram(int index, bool loadImageData, bool loadEvenFrames = false, bool loadOddFrames = false, int w0 = -999, int h0 =-999, int d0 = -999 ) const;

		int size() const;
        void setProjectionAngles(int tomogramIndex, int frame, RFLOAT xtilt, RFLOAT ytilt, RFLOAT zrot, RFLOAT xshift_angst, RFLOAT yshift_angst);

		void setCtf(int tomogramIndex, int frame, const CTF& ctf);
		void setDose(int tomogramIndex, int frame, double dose);
        void setTiltSeriesFile(int tomogramIndex, const std::string& filename);
		void setFiducialsFile(int tomogramIndex, const std::string& filename);
		void setDefocusSlope(int tomogramIndex, double slope);
        void applyTiltAngleOffset(int tomogramIndex, double offset);

		void setDeformation(
				int tomogramIndex,
				gravis::i2Vector gridSize,
				const std::string& deformationType,
				const std::vector<std::vector<double>>& coeffs);

		void clearDeformation();


		int getTomogramIndex(std::string tomogramName) const;
        std::string getTomogramName(int index) const;
		int getTomogramIndexSafely(std::string tomogramName) const;
		int getFrameCount(int index) const;
		int getMaxFrameCount() const;
		double getOriginalPixelSize(int index) const;
		double getTiltSeriesPixelSize(int index) const;
		std::string getOpticsGroupName(int index) const;

        std::vector<int> getFrameDoseOrder(int tomogramIndex) const;
        std::vector<int> getFrameDoseOrderIndex(int tomogramIndex) const;
        std::vector<int> getFrameTiltOrder(int tomogramIndex) const;
        std::vector<int> getFrameTiltOrderIndex(int tomogramIndex) const;

        int getImageIndexWithSmallestVisibleTiltAngle(int index, std::vector<bool> isVisible) const;

        // SHWS 6Apr2022: Make one big metadatatable with all movies/micrographs (to be used for motioncorr and ctffind runners)
        void generateSingleMetaDataTable(MetaDataTable &MDout, ObservationModel &obsModel);

        // SHWS 6Apr2022: Convert back from one big metadatatable into separate STAR files for each tilt serie
        void convertBackFromSingleMetaDataTable(MetaDataTable &MDin);

};

#endif
