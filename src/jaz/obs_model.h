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
		
		static void loadSafely(
				std::string particlesFn, std::string opticsFn, 
				ObservationModel& obsModel, 
				MetaDataTable& particlesMdt, MetaDataTable& opticsMdt);

        ObservationModel();
        ObservationModel(const MetaDataTable &opticsMdt);
		
		MetaDataTable opticsMdt;
		std::vector<double> angpix, lambda, Cs;
		std::vector<std::vector<double>> oddZernikeCoeffs;
		std::vector<std::map<int,Image<Complex>>> phaseCorr;

		
		// prediction
			
		void predictObservation(
				Projector &proj, const MetaDataTable &partMdt, long int particle,
				MultidimArray<Complex>& dest,
				bool applyCtf = true, bool shiftPhases = true, bool applyShift = true);
		
		Image<Complex> predictObservation(
				Projector &proj, const MetaDataTable &partMdt, long int particle,
				bool applyCtf = true, bool shiftPhases = true, bool applyShift = true);

        std::vector<Image<Complex>> predictObservations(
                Projector &proj, const MetaDataTable &partMdt, int threads,
                bool applyCtf = true, bool shiftPhases = true, bool applyShift = true);

		
		// correction
		
		void demodulatePhase(
				const MetaDataTable &partMdt, 
				long int particle,
				MultidimArray<Complex>& obsImage);
		
		/* for now, the programs assume that all optics groups have the same pixel size */
		bool allPixelSizesIdentical() const;
		
        double angToPix(double a, int s, int opticsGroup = 0) const;
        double pixToAng(double p, int s, int opticsGroup = 0) const;
		
		double getPixelSize(int opticsGroup = 0) const;
		
		
		/* duh */
		int numberOfOpticsGroups() const;
		
		/* Check whether the optics groups appear in the correct order.
		   This makes it possible to access a group g through:
		   
		       opticsMdt.getValue(label, dest, g-1);
		*/
		bool opticsGroupsSorted() const;
		
		/* Find all optics groups used in particles table partMdt
		   that are not defined in opticsMdt (should return an empty vector) */
		std::vector<int> findUndefinedOptGroups(const MetaDataTable& partMdt) const;
				
		/* Rename optics groups to enforce the correct order
		   and translate the indices in particle table partMdt.
		   (Merely changing the order in opticsMdt would fail if groups were missing.) */
		void sortOpticsGroups(MetaDataTable& partMdt);
		
		/* Return the set of optics groups present in partMdt */
		std::vector<int> getOptGroupsPresent(const MetaDataTable& partMdt) const;
		
		
		const Image<Complex>& getPhaseCorrection(int optGroup, int s);
		
		static bool containsAllNeededColumns(const MetaDataTable& partMdt);

};

#endif
