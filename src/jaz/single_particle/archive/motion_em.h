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

#ifndef MOTION_EM_H
#define MOTION_EM_H

#include <src/projector.h>
#include <src/metadata_table.h>
#include <src/image.h>
#include <src/complex.h>
#include <src/jaz/single_particle/legacy_obs_model.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/single_particle/parallel_ft.h>
#include <vector>

class MotionEM
{
    public:

        MotionEM(
                Projector& projector0,
                Projector& projector1,
                const ObservationModel& obsModel,
                MetaDataTable& viewParams,
                const std::vector<std::vector<Image<Complex>>>& movie,
                const std::vector<gravis::d2Vector>& globalPositions,
                const std::vector<double>& sigma2,
                const std::vector<Image<RFLOAT>>& damageWeights,
                double sig_pos,
                const std::vector<double>& sig_vel_initial,
                const std::vector<double>& sig_div_initial,
                double sig_cutoff,
                int threads);

            Projector& projector0;
            Projector& projector1;
            const ObservationModel& obsModel;
            MetaDataTable& viewParams;
            const std::vector<std::vector<Image<Complex>>>& movie;
            const std::vector<gravis::d2Vector>& globalPositions;
            const std::vector<double>& sigma2;
            const std::vector<Image<RFLOAT>>& damageWeights;
            double sig_pos, sig_cutoff;
            std::vector<double> sig_vel, sig_div;

            int threads;
            std::vector<ParFourierTransformer> fts_full, fts_pos, fts_vel;

            int pc, fc,
                s_full, sh_full,
                s_pos, sh_pos;

            std::vector<int> s_vel, sh_vel, sig_vel_class;

            bool initialized;

            std::vector<std::vector<Image<RFLOAT>>> posProb, velProb, initialCC;
            std::vector<Image<Complex>> pred;
            std::vector<Image<RFLOAT>> e_sum;

        void estimate(int iterations);
        void computeInitial();
        void iterate();

        void updateVelocities();
        void consolidateVelocities(int maxPc = -1);
        void smoothVelocities();

        Image<RFLOAT> blurVelocity(const Image<Complex>& velProbFs, double sigma, int f, int threadnum);
        Image<RFLOAT> adaptSize(const Image<RFLOAT>& img, int s);

        void updatePositions(bool backward, int maxPc = -1);

        std::vector<gravis::d2Vector> getTrack(int particle);
        std::vector<gravis::d2Vector> getGlobalTrack();



};

#endif
