#ifndef MOTION_HELPER_H
#define MOTION_HELPER_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/projector.h>
#include <src/complex.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/obs_model.h>
#include <src/jaz/parallel_ft.h>
#include <vector>


class MotionHelper
{
    public:

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

        static void noiseNormalize(
                const Image<Complex>& img,
                const std::vector<double> &sigma2,
                Image<Complex>& dest);

        static void writeTracks(
                const std::vector<std::vector<gravis::d2Vector>>& tracks,
                std::string fn);

        static std::vector<std::vector<gravis::d2Vector>> readTracks(
                std::string fn);

};

#endif
