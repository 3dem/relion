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

#ifndef MOTION_REFINEMENT_H
#define MOTION_REFINEMENT_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/projector.h>
#include <src/complex.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/single_particle/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/single_particle/legacy_obs_model.h>
#include <src/jaz/single_particle/parallel_ft.h>
#include <vector>

class ParticleMotionFit : public Optimization
{
    public:

        ParticleMotionFit(
                const std::vector<Image<float> >& correlation,
                RFLOAT lambda_vel, RFLOAT lambda_acc);

        double f(const std::vector<double>& x, void* tempStorage) const;


    private:

        const std::vector<Image<float> >& correlation;
        RFLOAT lambda_vel, lambda_acc;
};

class MotionFit : public DifferentiableOptimization
{
    public:

        MotionFit(
                const std::vector<std::vector<Image<float> > >& correlation,
                const std::vector<std::vector<RFLOAT> >& distWeights,
                RFLOAT lambda, RFLOAT mu);

        double f(const std::vector<double>& x, void* tempStorage) const;
        double f_data(const std::vector<double>& x) const;

        void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;


    private:

        const std::vector<std::vector<Image<float> > >& correlation;
        const std::vector<std::vector<RFLOAT> >& distWeights;
        RFLOAT lambda, mu;
};

class MotionRefinement
{
    public:

        static Image<RFLOAT> recompose(const std::vector<Image<RFLOAT> >& obs,
                                       const std::vector<double>& pos);

        static Image<RFLOAT> recompose(const std::vector<Image<Complex> >& obs,
                                       const std::vector<double>& pos);

        static Image<RFLOAT> averageStack(const std::vector<Image<RFLOAT> >& obs);

        static Image<RFLOAT> averageStack(const std::vector<Image<Complex> >& obs);

        static std::vector<std::vector<Image<RFLOAT>>> movieCC(
                Projector& projector0,
                Projector& projector1,
                const ObservationModel& obsModel,
                MetaDataTable& viewParams,
                const std::vector<std::vector<Image<Complex>>>& movie,
                const std::vector<double>& sigma2,
                const std::vector<Image<RFLOAT>>& damageWeights,
                std::vector<ParFourierTransformer>& fts, int threads);

        static std::vector<gravis::d2Vector> getGlobalTrack(
                const std::vector<std::vector<Image<RFLOAT>>>& movieCC);

        static std::vector<Image<RFLOAT>> addCCs(
                const std::vector<std::vector<Image<RFLOAT>>>& movieCC);

        static std::vector<gravis::d2Vector> getGlobalTrack(
                const std::vector<Image<RFLOAT>>& movieCcSum);

        static std::vector<gravis::d2Vector> getGlobalOffsets(
                const std::vector<std::vector<Image<RFLOAT>>>& movieCC,
                const std::vector<gravis::d2Vector>& globTrack,
                double sigma, int threads);

        static Image<float> crossCorrelation2D(const Image<Complex>& obs,
                const Image<Complex>& predConj,
                const Image<RFLOAT>& wgh, const std::vector<double> &sigma2);

        static Image<float> crossCorrelation2D(
                const Image<Complex>& obs,
                const Image<Complex>& predConj,
                const std::vector<double> &sigma2,
                bool probability = false,
                bool normalize = false);

        static void noiseNormalize(
                const Image<Complex>& img,
                const std::vector<double> &sigma2,
                Image<Complex>& dest);

        static std::vector<std::vector<gravis::d2Vector>> readTrack(
                std::string fn,
                int pc, int fc);

        static void writeTracks(
                const std::vector<std::vector<gravis::d2Vector>>& tracks,
                std::string fn);

        static std::vector<std::vector<gravis::d2Vector>> readTracks(
                std::string fn);

        static gravis::d3Vector measureValueScaleReal(
                const Image<Complex>& obs,
                const Image<Complex>& ref);

        static gravis::d3Vector measureValueScale(
                const Image<Complex>& obs,
                const Image<Complex>& ref);

        static void testCC(
                const Image<Complex>& obs,
                const Image<Complex>& predConj,
                const std::vector<double> &sigma2);

        static Image<RFLOAT> zeroPad(
                const Image<RFLOAT>& img, RFLOAT ratio = 2.0, RFLOAT taper = 0.1);

        static std::vector<Image<float> > collectiveMotion(
                const std::vector<std::vector<Image<float> > >& correlation);        

        static std::vector<std::vector<Image<float>>> blockMotion(
                const std::vector<std::vector<Image<float> > >& correlation,
                std::vector<gravis::d2Vector> positions, int parts, int micrographWidth,
                std::vector<int>& numbers);

        static std::vector<gravis::d2Vector> findMaxima(std::vector<Image<float> > &corrSum);

        static std::vector<std::vector<gravis::d2Vector> > computeInitialPositions(
                const std::vector<std::vector<Image<float> > >& correlation);

        static std::vector<std::vector<gravis::d2Vector> > optimize(
                const std::vector<std::vector<Image<float> > >& correlation,
                const std::vector<gravis::d2Vector>& positions,
                const std::vector<std::vector<gravis::d2Vector> >& initial,
                double lambda, double mu, double sigma);




        static std::vector<std::vector<Image<RFLOAT> > > visualize(
                const std::vector<std::vector<gravis::d2Vector> >& positions,
                int pc, int fc, int w, int h);

        static std::vector<Image<RFLOAT> > collapsePaths(
                const std::vector<std::vector<Image<RFLOAT> > >& paths);

        static std::vector<std::vector<gravis::d2Vector> > unpack(
                const std::vector<double>& pos, int pc, int fc);

        static std::vector<double> pack(
                const std::vector<std::vector<gravis::d2Vector> >& pos);

        static std::vector<std::vector<gravis::d2Vector> > readCollectivePaths(std::string filename);

        static void writeCollectivePaths(
                const std::vector<std::vector<gravis::d2Vector> >& data, std::string filename);

        static std::vector<std::pair<int, std::vector<gravis::d2Vector> > >
                readPaths(std::string fn, int imgNum, int blockNum, int frameNum);

        static std::vector<std::vector<gravis::d2Vector>> centerBlocks(
                std::string filenameFull, std::string filenameBlocks,
                int imgCount, int blockCount, int frameCount,
                std::vector<int> &particleNumbers,
                std::vector<int> &totalParticleNumbers);

};

#endif
