#ifndef OBS_MODEL_H
#define OBS_MODEL_H

#include <src/image.h>
#include <src/fftw.h>
#include <src/complex.h>
#include <src/metadata_table.h>
#include <src/projector.h>

class BackProjector;

class ObservationModel
{
    public:

        ObservationModel();
        ObservationModel(const MetaDataTable &opticsMdt);
		
		MetaDataTable opticsMdt;
		std::vector<double> angpix, lambda, Cs;

			
		void predictObservation(
				Projector &proj, const MetaDataTable &mdt, long int particle,
				MultidimArray<Complex>& dest,
				bool applyCtf = true, bool applyTilt = true, bool applyShift = true) const;
		
		Image<Complex> predictObservation(
				Projector &proj, const MetaDataTable &mdt, long int particle,
				bool applyCtf = true, bool applyTilt = true, bool applyShift = true) const;

        std::vector<Image<Complex>> predictObservations(
                Projector &proj, const MetaDataTable &mdt, int threads,
                bool applyCtf = true, bool applyTilt = true, bool applyShift = true) const;

        void insertObservation(
                const Image<Complex>& img, BackProjector &bproj,
                const MetaDataTable& mdt, long int particle,
                bool applyCtf, bool applyTilt,
                double shift_x = 0.0, double shift_y = 0.0);

		
        double angToPix(double a, int s, int opticsGroup);
        double pixToAng(double p, int s, int opticsGroup);
		
		double getPixelSize(int opticsGroup = 0);
		bool allPixelSizesIdentical();
		
		
		/* duh */
		int numberOfOpticsGroups();
		
		/* Check whether the optics groups appear in the correct order.
		   This makes it possible to access a group g through:
		   
		       opticsMdt.getValue(label, dest, g-1);
		*/
		bool opticsGroupsSorted();
		
		/* Find all optics groups used in particles table partMdt
		   that are not defined in opticsMdt (should return an empty vector) */
		std::vector<int> findUndefinedOptGroups(const MetaDataTable& partMdt);
				
		/* Rename optics groups to enforce the correct order
		   and translate the indices in particle table partMdt.
		   (Merely changing the order in opticsMdt would fail if groups were missing.) */
		void sortOpticsGroups(MetaDataTable& partMdt);
		
		/* Return the set of optics groups present in partMdt */
		std::vector<int> getOptGroupsPresent(const MetaDataTable& partMdt);
		
		static bool containsAllNeededColumns(const MetaDataTable& mdt);

};

#endif
