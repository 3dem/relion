/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
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

#ifndef OBS_MODEL_H
#define OBS_MODEL_H

#include <src/image.h>
#include <src/fftw.h>
#include <src/complex.h>
#include <src/metadata_table.h>
#include <src/projector.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t2Matrix.h>


class BackProjector;


class ObservationModel
{

    public:

		ObservationModel();
		ObservationModel(const MetaDataTable &opticsMdt, bool do_die_upon_error = true);


			MetaDataTable opticsMdt;

			bool
				hasEvenZernike, hasOddZernike, hasMagMatrices,
				hasBoxSizes, hasMultipleMtfs;



	protected:

		// cached values - protected to prevent users from accidentally changing them,
		// expecting the changes to propagate into the optics star-file
		std::vector<double> angpix, originalAngpix, lambda, Cs;
		std::vector<int> boxSizes;
		std::vector<bool> CtfPremultiplied;
		std::vector<std::vector<double> > evenZernikeCoeffs, oddZernikeCoeffs;
		std::vector<Matrix2D<RFLOAT> > magMatrices;
		std::vector<std::string> fnMtfs, groupNames;

		// cached aberration effects for a set of given image sizes
		// e.g.: phaseCorr[opt. group][img. height](x,y)
		std::vector<std::map<int,BufferedImage<Complex> > > phaseCorr;
		std::vector<std::map<int,BufferedImage<RFLOAT> > > gammaOffset, mtfImage;
		std::map<int,BufferedImage<RFLOAT> > avgMtfImage;


	public:

		// Prediction //

		void predictObservation_DC(
				const RawImage<Complex>& reference_map,
				const RawImage<Complex>* amplitude_map,
				const MetaDataTable& particle_table, long int particle_index,
				RawImage<Complex>& destination,
				double reference_pixel_size,
				bool applyCtf = true, bool shiftPhases = true, bool applyShift = true,
				bool applyMtf = true, bool applyCtfPadding = false);

		Volume<gravis::t2Vector<Complex> > predictComplexGradient(
				Projector &proj, const MetaDataTable &partMdt,
				long int particle, double angpix_ref,
				bool applyCtf = true, bool shiftPhases = true, bool applyShift = true,
				bool applyMtf = true, bool applyCtfPadding = false);

		void predictObservation(
				Projector &proj, const MetaDataTable &partMdt, long int particle,
				MultidimArray<Complex>& dest, double angpix_ref,
				bool applyCtf = true, bool shiftPhases = true, bool applyShift = true,
				bool applyMtf = true, bool applyCtfPadding = false);


		// Correction //

		// divide by MTF of detector (using cache)
		void divideByMtf(
				const MetaDataTable& partMdt, long particle, MultidimArray<Complex>& obsImage,
				bool do_multiply_instead = false, bool do_correct_average_mtf = true);

		void divideByMtf(
				int opticsGroup, MultidimArray<Complex>& obsImage,
				bool do_multiply_instead = false, bool do_correct_average_mtf = true);

		// 2D image with the MTF (cached)
		// Nyquist X is positive, Y is negative (non-FFTW!!)
		const BufferedImage<RFLOAT>& getMtfImage(int optGroup, int s);

		// 2D image with the average MTF (cached)
		const BufferedImage<RFLOAT>& getAverageMtfImage(int s);


		// apply effect of antisymmetric aberration (using cache)
		void demodulatePhase(
				int optGroup, MultidimArray<Complex>& obsImage,
				bool do_modulate_instead = false);

		// syntactic sugar
		void demodulatePhase(
				const MetaDataTable &partMdt,
				long int particle, MultidimArray<Complex>& obsImage,
				bool do_modulate_instead = false);

		// effect of antisymmetric aberration (cached)
		// Nyquist X is positive, Y is negative (non-FFTW!!)
		const BufferedImage<Complex>& getPhaseCorrection(int optGroup, int s);

		// effect of symmetric aberration (cached)
		// Nyquist X is positive, Y is negative (non-FFTW!!)
		const BufferedImage<RFLOAT>& getGammaOffset(int optGroup, int s);

		Matrix2D<RFLOAT> applyAnisoMag(Matrix2D<RFLOAT> A3D, int opticsGroup);

		Matrix2D<RFLOAT> applyScaleDifference(
				Matrix2D<RFLOAT> A3D, int opticsGroup, int s3D, double angpix3D);


		// Saving and loading

		// tablename can be "particles", "micrographs" or "movies".
		// If tablename is "discover", the function will try to read
		// the data table with all three names (in that order).

		static void loadSafely(
				std::string filename, ObservationModel& obsModel,
				MetaDataTable& particlesMdt, std::string tablename = "particles",
				int verb = 0, bool do_die_upon_error = true);

		static void saveNew(
				MetaDataTable& particlesMdt, MetaDataTable& opticsMdt,
				std::string filename, std::string _tablename = "particles");

		void save(
				MetaDataTable& particlesMdt, std::string filename,
				std::string _tablename = "particles");


		// Bureaucracy

		static bool containsAllColumnsNeededForPrediction(
				const MetaDataTable& partMdt);

		bool allPixelSizesIdentical() const;
		bool allBoxSizesIdentical() const;

		double angToPix(double a, int s, int opticsGroup) const;
		double pixToAng(double p, int s, int opticsGroup) const;

		double getPixelSize(int opticsGroup) const;
		std::vector<double> getPixelSizes() const;

		double getWavelength(int opticsGroup) const;
		std::vector<double> getWavelengths() const;

		double getSphericalAberration(int opticsGroup) const;
		std::vector<double> getSphericalAberrations() const;

		int getBoxSize(int opticsGroup) const;
		void getBoxSizes(std::vector<int>& sDest, std::vector<int>& shDest) const;

		// These do NOT update the metadata table!
		// These are only to change prediction etc.
		void setBoxSize(int opticsGroup, int newBoxSize);
		void setPixelSize(int opticsGroup, RFLOAT newPixelSize);

		Matrix2D<RFLOAT> getMagMatrix(int opticsGroup) const;
		std::vector<Matrix2D<RFLOAT> > getMagMatrices() const;
		void setMagMatrix(int opticsGroup, const Matrix2D<RFLOAT>& M);

		// returns a zero-indexed value (it exists to ensure this)
		int getOpticsGroup(const MetaDataTable& particlesMdt, long int particle = -1) const;

		bool getCtfPremultiplied(int og) const;
		void setCtfPremultiplied(int og, bool val);

		std::string getGroupName(int og);

		bool allPixelAndBoxSizesIdentical(const MetaDataTable& mdt);
		bool containsGroup(const MetaDataTable& mdt, int group);

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
		std::vector<int> getOptGroupsPresent_oneBased(const MetaDataTable& partMdt) const;

		/* Return the set of optics groups present in partMdt */
		std::vector<int> getOptGroupsPresent_zeroBased(const MetaDataTable& partMdt) const;

		std::vector<std::pair<int, std::vector<int> > > splitParticlesByOpticsGroup(
				const MetaDataTable& partMdt) const;
};

#endif
