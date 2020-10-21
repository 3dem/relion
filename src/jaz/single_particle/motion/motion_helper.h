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

#ifndef MOTION_HELPER_H
#define MOTION_HELPER_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/projector.h>
#include <src/complex.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/single_particle/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/parallel_ft.h>
#include <vector>


class MotionHelper
{
    public:

        static std::vector<std::vector<Image<RFLOAT>>> movieCC(
                const std::vector<std::vector<Image<Complex>>>& movie,
                const std::vector<Image<Complex>>& preds,
                const std::vector<Image<RFLOAT>>& damageWeights,
                double pad, int threads);

        // deprecated: use the one above!
        /*static std::vector<std::vector<Image<RFLOAT>>> movieCC(
                Projector& projector0,
                Projector& projector1,
                const ObservationModel& obsModel,
                const MetaDataTable& viewParams,
                const std::vector<std::vector<Image<Complex>>>& movie,
                const std::vector<double>& sigma2,
                const std::vector<Image<RFLOAT>>& damageWeights,
                std::vector<ParFourierTransformer>& fts, int threads);*/

        /*static std::vector<gravis::d2Vector> getGlobalTrack(
                const std::vector<std::vector<Image<RFLOAT>>>& movieCC, double cc_pad);*/

        static std::vector<Image<RFLOAT>> addCCs(
                const std::vector<std::vector<Image<RFLOAT>>>& movieCC);

        static std::vector<gravis::d2Vector> getGlobalTrack(
                const std::vector<Image<RFLOAT>>& movieCcSum, double cc_pad);

        static std::vector<gravis::d2Vector> getGlobalOffsets(
                const std::vector<std::vector<Image<RFLOAT>>>& movieCC,
                const std::vector<std::vector<gravis::d2Vector>>& initialTracks, 
				double cc_pad, double sigma, int wMax, int hMax, int threads);

        static void noiseNormalize(
                const Image<Complex>& img,
                const std::vector<double> &sigma2,
                Image<Complex>& dest);

        static void writeTracks(
                const std::vector<std::vector<gravis::d2Vector>>& tracksInPix,
                std::string fn, double angpix);

	// both stack number and frame number are 0-indexed in the array and STAR file
        static std::vector<std::vector<gravis::d2Vector>> readTracksInPix(
                std::string fn, double angpix);

};

#endif
