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

#ifndef NEW_SPA_REFERENCE_MAP_H
#define NEW_SPA_REFERENCE_MAP_H

#include <src/image.h>
#include <src/projector.h>
#include <src/jaz/single_particle/volume.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/image/buffered_image.h>

#include <vector>

class ObservationModel;
class LegacyObservationModel;

class NewReferenceMap
{
	public:

		typedef enum {Own, Opposite} HalfSet;

		NewReferenceMap();

			// input parameters:
			std::string
				phase_file_names[2], amplitude_file_names[2],
				maskFn, fscFn;

			double paddingFactor;

			// data:
			BufferedImage<Complex> phaseMap[2], amplitudeMap[2];
			BufferedImage<RFLOAT> freqWeight, mask;
			std::vector<double> freqWeight1D;
			int k_out, s, sh;
			bool dualContrast, hasMask;
			double angpix;


		void read(IOParser& parser, int argc, char *argv[]);
		void load(int verb, bool debug);

		Image<RFLOAT> getHollowWeight(double kmin_ang, int s_out, double angpix_out);

		std::vector<Image<Complex>> predictAll(
				const MetaDataTable& mdt,
				ObservationModel& obs,
				HalfSet hs, int threads,
				bool applyCtf = true,
				bool applyTilt = true,
				bool applyShift = true,
				bool applyMtf = true,
				bool applyCtfPadding = false);

		Image<Complex> predict(
				const MetaDataTable& mdt, int p,
				ObservationModel& obs,
				HalfSet hs,
				bool applyCtf = true,
				bool applyTilt = true,
				bool applyShift = true,
				bool applyMtf = true,
				bool applyCtfPadding = false);

		std::vector<Volume<gravis::t2Vector<Complex>>> predictAllComplexGradients(
				const MetaDataTable& mdt,
				ObservationModel& obs,
				HalfSet hs, int threads,
				bool applyCtf = true,
				bool applyTilt = true,
				bool applyShift = true,
				bool applyMtf = true,
				bool applyCtfPadding = false);

		Volume<gravis::t2Vector<Complex>> predictComplexGradient(
				const MetaDataTable& mdt, int p,
				ObservationModel& obs,
				HalfSet hs,
				bool applyCtf = true,
				bool applyTilt = true,
				bool applyShift = true,
				bool applyMtf = true,
				bool applyCtfPadding = false);

		double angToPix(double a) const;
		double pixToAng(double p) const;

		void drawFSC(
				std::string fsc_filename,
				std::vector<double>& dest1D,
				RawImage<RFLOAT>& dest2D,
				double thresh = 0.143) const;

	protected:

		void presharpen(BufferedImage<RFLOAT>& map);
};

#endif
