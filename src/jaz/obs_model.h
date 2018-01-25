#ifndef OBS_MODEL_H
#define OBS_MODEL_H

#include <src/image.h>
#include <src/fftw.h>
#include <src/complex.h>
#include <src/metadata_table.h>
#include <src/projector.h>
#include <src/backprojector.h>


class ObservationModel
{
    public:

        ObservationModel();
        ObservationModel(double angpix);
        ObservationModel(double angpix, double Cs, double voltage, double beamtilt_x, double beamtilt_y);

            double angpix, lambda, Cs;
            double beamtilt_x, beamtilt_y;
            bool hasTilt;

        Image<Complex> predictObservation(
                Projector &proj, MetaDataTable &mdt, int particle,
                bool applyCtf, bool applyTilt,
                double deltaRot = 0.0,
                double deltaTilt = 0.0,
                double deltaPsi = 0.0) const;

        void insertObservation(
                const Image<Complex>& img, BackProjector &bproj,
                MetaDataTable& mdt, int particle,
                bool applyCtf, bool applyTilt,
                double shift_x = 0.0, double shift_y = 0.0);

};

#endif
