#ifndef MICROGRAPH_HANDLER_H
#define MICROGRAPH_HANDLER_H

#include <string>
#include <map>

#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/parallel_ft.h>

#include <src/micrograph_model.h>
#include <src/image.h>

class MicrographHandler
{
    public:

        MicrographHandler();

            int nr_omp_threads, firstFrame, lastFrame;
            double movie_angpix, coords_angpix, hotCutoff;

            bool preextracted, debug, saveMem, ready;

            std::string corrMicFn;
            std::string movie_path, meta_path, movie_ending;
            std::string gain_path, last_gainFn;

            gravis::t2Vector<int> micrograph_size;


        // initialise corrected/uncorrected micrograph dictionary, then
        // load first movie (or read corrected_micrographs.star) to obtain:
        //  fc, micrograph_xsize, micrograph_ysize, motionEstimator.dosePerFrame
        void init(
            // in:
                const MetaDataTable& mdt,
                double angpix, bool verb,
                int nr_omp_threads,
            // out:
                int& fc,
                double& dosePerFrame,
                std::string& metaFn);

        // load a movie and extract all particles
        // returns a per-particle vector of per-frame images of size (s/2+1) x s
        std::vector<std::vector<Image<Complex>>> loadMovie(const MetaDataTable& mdt, int s,
                double angpix, std::vector<ParFourierTransformer>& fts);

        /* Load a movie as above and also write tracks of particles at 'pos' into 'tracks'.
           If 'unregGlob' is set, also write the global component of motion into 'globComp'.*/
        std::vector<std::vector<Image<Complex>>> loadMovie(
                const MetaDataTable& mdt, int s, double angpix,
                std::vector<ParFourierTransformer>& fts,
                const std::vector<gravis::d2Vector>& pos,
                std::vector<std::vector<gravis::d2Vector>>& tracks,
                bool unregGlob, std::vector<gravis::d2Vector>& globComp);

    protected:

            Micrograph micrograph;
            Image<RFLOAT> lastGainRef;

            bool hasCorrMic;
            std::map<std::string, std::string> mic2meta;

        void loadInitial(
                const MetaDataTable& mdt, double angpix, bool verb,
                int& fc, double& dosePerFrame, std::string& metaFn);

        std::string getMetaName(std::string micName);

};

#endif
